"""Tests for QUANDARY_BASE_DATADIR environment variable support in the new interface."""

import os
import pytest
from quandary.new import setup_quandary, optimize, simulate

# Mark all tests in this file as regression tests
pytestmark = pytest.mark.regression

BASE_DATADIR = "QUANDARY_BASE_DATADIR"


def run_optimize(datadir):
    """Run a minimal optimization with the given output directory."""
    setup = setup_quandary(
        nessential=[2],
        transition_frequency=[4.0],
        selfkerr=[0.2],
        final_time=1.0,
        ntime=10,
        spline_order=0,
        output_directory=datadir,
        verbose=False,
    )
    setup.optim_maxiter = 1
    return optimize(
        setup,
        targetstate=[0.0, 1.0],
        quiet=True,
    )


def run_simulate(datadir):
    """Run a minimal simulation with the given output directory."""
    setup = setup_quandary(
        nessential=[2],
        transition_frequency=[4.0],
        selfkerr=[0.2],
        final_time=1.0,
        ntime=10,
        spline_order=0,
        output_directory=datadir,
        verbose=False,
    )
    return simulate(
        setup,
        quiet=True,
    )


test_cases = [
    run_optimize,
    run_simulate,
]


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
    """Ensure QUANDARY_BASE_DATADIR is restored to its previous state after the test."""
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
        "optim_history.dat",
        "params.dat",
        "control0.dat",
    ]

    assert os.path.exists(datadir), f"directory {datadir} does not exist"
    for f in expected_output_files:
        assert os.path.exists(os.path.join(datadir, f)), f"file {f} does not exist in {datadir}"


@pytest.mark.parametrize("run_quandary", test_cases)
def test_relative_path_without_env_var(run_quandary, request, cd_tmp_path, clean_env_var):
    datadir_name = request.node.name
    datadir_path = os.path.join(os.getcwd(), datadir_name)

    run_quandary(datadir=datadir_name)

    assert_output_files(datadir_path)


@pytest.mark.parametrize("run_quandary", test_cases)
def test_absolute_path_without_env_var(run_quandary, request, tmp_path, clean_env_var):
    datadir_name = request.node.name
    datadir_path = os.path.join(str(tmp_path), datadir_name)

    run_quandary(datadir=datadir_path)

    assert_output_files(datadir_path)


@pytest.mark.parametrize("run_quandary", test_cases)
def test_relative_path_with_env_var(run_quandary, request, tmp_path, clean_env_var):
    base_dir = str(tmp_path)
    os.environ[BASE_DATADIR] = base_dir
    datadir_name = request.node.name
    datadir_path = os.path.join(base_dir, datadir_name)

    run_quandary(datadir=datadir_name)

    assert_output_files(datadir_path)


@pytest.mark.parametrize("run_quandary", test_cases)
def test_absolute_path_with_env_var(run_quandary, request, tmp_path, clean_env_var):
    os.environ[BASE_DATADIR] = "should_not_use_this/path"
    datadir_name = request.node.name
    datadir_path = os.path.join(str(tmp_path), datadir_name)

    run_quandary(datadir=datadir_path)

    assert_output_files(datadir_path)
    assert not os.path.exists(os.environ[BASE_DATADIR])


@pytest.mark.parametrize("run_quandary", test_cases)
def test_nonexistent_base_directory(run_quandary, request, tmp_path, clean_env_var):
    nonexistent_path = os.path.join(str(tmp_path), "nonexistent_directory")
    os.environ[BASE_DATADIR] = nonexistent_path
    datadir_name = "some_output_dir"

    with pytest.raises(ValueError) as excinfo:
        run_quandary(datadir=datadir_name)

    assert "non-existent path" in str(excinfo.value)
    assert nonexistent_path in str(excinfo.value)


@pytest.mark.parametrize("run_quandary", test_cases)
def test_file_as_base_directory(run_quandary, request, tmp_path, clean_env_var):
    file_path = os.path.join(str(tmp_path), "this_is_a_file.txt")
    with open(file_path, 'w') as f:
        f.write("This is a file, not a directory")

    os.environ[BASE_DATADIR] = file_path
    datadir_name = "some_output_dir"

    with pytest.raises(ValueError) as excinfo:
        run_quandary(datadir=datadir_name)

    assert "not a directory" in str(excinfo.value)
    assert file_path in str(excinfo.value)
