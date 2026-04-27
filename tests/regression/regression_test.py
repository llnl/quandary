import glob
import os
import re
import subprocess
import pandas as pd
import pytest
from pydantic import BaseModel, TypeAdapter
from typing import List

from tests.utils.common import build_mpi_command

# Mark all tests in this file as regression tests
pytestmark = pytest.mark.regression

REL_TOL = 1.0e-7
ABS_TOL = 1.0e-15

# GPU runs have non-deterministic FP reduction order that perturbs optimization
# trajectories, so we compare only final converged results with looser tolerance.
GPU_REL_TOL = 5.0e-2
GPU_ABS_TOL = 1.0e-4
GPU_OPTION_TOKENS = ("kokkos", "cuda", "cusparse", "rocm", "hip", "viennacl")
OPTIM_HISTORY_FILE = "optim_history.dat"
OPTIM_FAVG_HEADER = "F_avg"

BASE_DIR = "base"
DATA_OUT_DIR = "data_out"

TEST_PATH = os.path.dirname(os.path.realpath(__file__))
TEST_CASES_PATH = os.path.join(TEST_PATH, "test_cases.json")
QUANDARY_PATH = os.path.join(TEST_PATH, "..", "..", "quandary")


class Case(BaseModel):
    simulation_name: str
    files_to_compare: List[str]
    number_of_processes: List[int]


def load_test_cases():
    with open(TEST_CASES_PATH) as test_cases_file:
        ta = TypeAdapter(List[Case])
        test_cases = ta.validate_json(test_cases_file.read())
        return test_cases


TEST_CASES = load_test_cases()


def is_gpu_mode(petsc_options: str) -> bool:
    """Detect GPU execution from tokens in PETSc runtime options (e.g. aijkokkos, aijcusparse)."""
    lower = petsc_options.lower()
    return any(token in lower for token in GPU_OPTION_TOKENS)


def find_column_index(file_path: str, column_name: str) -> int:
    """Find a column's positional index by parsing quoted tokens from the header line.

    pandas with sep=\"\\s+\" splits header names like \"LS step\" into separate tokens,
    so DataFrame column names misalign with data. Pulling \"...\" tokens directly from
    the file gives the true ordering.
    """
    with open(file_path, 'r') as f:
        header = f.readline()
    names = re.findall(r'"([^"]*)"', header)
    return names.index(column_name)


@pytest.mark.parametrize("test_case", TEST_CASES, ids=lambda x: x.simulation_name)
def test_eval(test_case: Case, request):
    exact = request.config.getoption("--exact")
    mpi_exec = request.config.getoption("--mpi-exec")
    mpi_opt = request.config.getoption("--mpi-opt")
    config_format = request.config.getoption("--config-format")
    petsc_options = request.config.getoption("--petsc-options")
    gpu_mode = is_gpu_mode(petsc_options)

    simulation_name = test_case.simulation_name
    files_to_compare = test_case.files_to_compare
    number_of_processes_list = test_case.number_of_processes

    simulation_dir = os.path.join(TEST_PATH, simulation_name)
    config_file = os.path.join(simulation_dir, f"{simulation_name}.{config_format}")

    # Verify config file exists
    if not os.path.exists(config_file):
        pytest.skip(f"Config file {config_file} not found for format '{config_format}'")

    for number_of_processes in number_of_processes_list:
        run_test(simulation_dir, number_of_processes, config_file, files_to_compare, exact, mpi_exec, mpi_opt, petsc_options, gpu_mode)


def run_test(simulation_dir, number_of_processes, config_file, files_to_compare, exact, mpi_exec, mpi_opt, petsc_options, gpu_mode):
    os.chdir(simulation_dir)

    command = build_mpi_command(
        mpi_exec=mpi_exec,
        num_processes=number_of_processes,
        mpi_opt=mpi_opt,
        quandary_path=QUANDARY_PATH,
        config_file=config_file,
        petsc_options=petsc_options)
    print(f"Running command: \"{' '.join(command)}\"")
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    print("STDOUT:\n", result.stdout)
    print("STDERR:\n", result.stderr)
    assert result.returncode == 0

    # GPU FP non-determinism makes L-BFGS land on different iterations and pulse
    # parameters, so all derived files (params, grad, rho, populations, expected)
    # diverge from CPU baselines. Only the final fidelity is a meaningful check.
    if gpu_mode and OPTIM_HISTORY_FILE in files_to_compare:
        files_to_compare = [OPTIM_HISTORY_FILE]

    matching_files = [file for pattern in files_to_compare
                      for file in glob.glob(os.path.join(simulation_dir, BASE_DIR, pattern))]
    for expected in matching_files:
        file_name = os.path.basename(expected)
        output = os.path.join(simulation_dir, DATA_OUT_DIR, file_name)
        compare_files(file_name, output, expected, exact, gpu_mode)


def compare_files(file_name, output, expected, exact, gpu_mode):
    df_output = pd.read_csv(output, sep="\\s+", header=get_header(output))
    df_expected = pd.read_csv(expected, sep="\\s+", header=get_header(expected))

    if gpu_mode:
        rtol, atol = GPU_REL_TOL, GPU_ABS_TOL
        # Compare only final-row F_avg: trajectory and iter count differ under GPU FP noise.
        if file_name == OPTIM_HISTORY_FILE:
            favg_col = find_column_index(expected, OPTIM_FAVG_HEADER)
            df_output = df_output.iloc[[-1], [favg_col]].reset_index(drop=True)
            df_expected = df_expected.iloc[[-1], [favg_col]].reset_index(drop=True)
    else:
        rtol, atol = REL_TOL, ABS_TOL

    pd.testing.assert_frame_equal(df_output, df_expected, rtol=rtol, atol=atol, obj=file_name, check_exact=exact)


def get_header(path):
    with open(path, 'r') as file:
        first_line = file.readline().strip()
        return 0 if first_line.startswith('#') else None
