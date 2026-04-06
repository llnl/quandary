import os
import pytest
from quandary import Quandary

# Mark all tests in this file as regression tests
pytestmark = pytest.mark.regression

BASE_DATADIR = "QUANDARY_BASE_DATADIR"


def quandary_simulate(datadir, mpi_exec):
    return Quandary(
        Ne=[2],
        Ng=[0],
        freq01=[4.0],
        selfkerr=[0.2],
        T=1.0,
        nsteps=10,
        maxiter=1,
        spline_order=0,
    ).simulate(
        datadir=datadir,
        mpi_exec=mpi_exec,
        maxcores=2
    )


def quandary_optimize(datadir, mpi_exec):
    return Quandary(
        Ne=[2],
        targetstate=[0.0, 1.0],
        initialcondition="basis",
        tol_infidelity=1e-2,
        nsteps=1,
        maxiter=1,
        spline_order=0
    ).optimize(
        datadir=datadir,
        mpi_exec=mpi_exec,
        maxcores=2
    )


test_cases = [
    quandary_optimize,
    quandary_simulate,
]


@pytest.mark.parametrize("quandary", test_cases)
def test_relative_output_path_without_env_var(quandary, request, cd_tmp_path, clean_env_var, mpi_exec):
    datadir_name = request.node.name
    datadir_path = os.path.join(os.getcwd(), datadir_name)

    quandary(datadir=datadir_name, mpi_exec=mpi_exec)

    assert_output_files(datadir_path)


@pytest.mark.parametrize("quandary", test_cases)
def test_absolute_output_path_without_env_var(quandary, request, tmp_path, clean_env_var, mpi_exec):
    datadir_name = request.node.name
    datadir_path = os.path.join(tmp_path, datadir_name)

    quandary(datadir=datadir_path, mpi_exec=mpi_exec)

    assert_output_files(datadir_path)


@pytest.mark.parametrize("quandary", test_cases)
def test_relative_output_path_with_env_var(quandary, request, tmp_path, clean_env_var, mpi_exec):
    base_dir = str(tmp_path)
    os.environ[BASE_DATADIR] = base_dir
    datadir_name = request.node.name
    datadir_path = os.path.join(base_dir, datadir_name)

    quandary(datadir=datadir_name, mpi_exec=mpi_exec)

    assert_output_files(datadir_path)


@pytest.mark.parametrize("quandary", test_cases)
def test_absolute_output_path_with_env_var(quandary, request, tmp_path, clean_env_var, mpi_exec):
    os.environ[BASE_DATADIR] = "should_not_use_this/path"
    datadir_name = request.node.name
    datadir_path = os.path.join(tmp_path, datadir_name)

    quandary(datadir=datadir_path, mpi_exec=mpi_exec)

    assert_output_files(datadir_path)
    assert not os.path.exists(os.environ[BASE_DATADIR])


@pytest.mark.parametrize("quandary", test_cases)
def test_nonexistent_base_directory(quandary, request, tmp_path, clean_env_var, mpi_exec):
    nonexistent_path = os.path.join(tmp_path, "nonexistent_directory")
    os.environ[BASE_DATADIR] = nonexistent_path
    datadir_name = "some_output_dir"

    with pytest.raises(ValueError) as excinfo:
        quandary(datadir=datadir_name, mpi_exec=mpi_exec)

    assert "non-existent path" in str(excinfo.value)
    assert nonexistent_path in str(excinfo.value)


@pytest.mark.parametrize("quandary", test_cases)
def test_file_as_base_directory(quandary, request, tmp_path, clean_env_var, mpi_exec):
    file_path = os.path.join(tmp_path, "this_is_a_file.txt")
    with open(file_path, 'w') as f:
        f.write("This is a file, not a directory")

    os.environ[BASE_DATADIR] = file_path
    datadir_name = "some_output_dir"

    with pytest.raises(ValueError) as excinfo:
        quandary(datadir=datadir_name, mpi_exec=mpi_exec)

    assert "not a directory" in str(excinfo.value)
    assert file_path in str(excinfo.value)


@pytest.fixture
def cd_tmp_path(tmp_path):
    """Change to a temporary directory for the test and return afterward."""
    original_cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield tmp_path
    finally:
        os.chdir(original_cwd)


@pytest.fixture
def clean_env_var():
    """Fixture to ensure env var is restored to previous state after the tests"""
    orig_value = os.environ.get(BASE_DATADIR)

    if BASE_DATADIR in os.environ:
        del os.environ[BASE_DATADIR]

    yield

    if orig_value is not None:
        os.environ[BASE_DATADIR] = orig_value
    elif BASE_DATADIR in os.environ:
        del os.environ[BASE_DATADIR]


def assert_output_files(datadir):
    expected_output_files = [
        "config.toml",
        "optim_history.dat",
        "params.dat",
        "control0.dat"
    ]

    assert os.path.exists(datadir), f"directory {datadir} does not exist"
    for file in expected_output_files:
        assert os.path.exists(os.path.join(datadir, file)), f"file {file} does not exist"
