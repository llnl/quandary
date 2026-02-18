import glob
import os
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


@pytest.mark.parametrize("test_case", TEST_CASES, ids=lambda x: x.simulation_name)
def test_eval(test_case: Case, request):
    exact = request.config.getoption("--exact")
    mpi_exec = request.config.getoption("--mpi-exec")
    mpi_opt = request.config.getoption("--mpi-opt")
    config_format = request.config.getoption("--config-format")

    simulation_name = test_case.simulation_name
    files_to_compare = test_case.files_to_compare
    number_of_processes_list = test_case.number_of_processes

    simulation_dir = os.path.join(TEST_PATH, simulation_name)
    config_file = os.path.join(simulation_dir, f"{simulation_name}.{config_format}")

    # Verify config file exists
    if not os.path.exists(config_file):
        pytest.skip(f"Config file {config_file} not found for format '{config_format}'")

    for number_of_processes in number_of_processes_list:
        run_test(simulation_dir, number_of_processes, config_file, files_to_compare, exact, mpi_exec, mpi_opt)


def run_test(simulation_dir, number_of_processes, config_file, files_to_compare, exact, mpi_exec, mpi_opt):
    os.chdir(simulation_dir)

    command = build_mpi_command(
        mpi_exec=mpi_exec,
        num_processes=number_of_processes,
        mpi_opt=mpi_opt,
        quandary_path=QUANDARY_PATH,
        config_file=config_file)
    print(f"Running command: \"{' '.join(command)}\"")
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    print("STDOUT:\n", result.stdout)
    print("STDERR:\n", result.stderr)
    assert result.returncode == 0

    matching_files = [file for pattern in files_to_compare
                      for file in glob.glob(os.path.join(simulation_dir, BASE_DIR, pattern))]
    for expected in matching_files:
        file_name = os.path.basename(expected)
        output = os.path.join(simulation_dir, DATA_OUT_DIR, file_name)
        compare_files(file_name, output, expected, exact)


def compare_files(file_name, output, expected, exact):
    df_output = pd.read_csv(output, sep="\\s+", header=get_header(output))
    df_expected = pd.read_csv(expected, sep="\\s+", header=get_header(expected))
    pd.testing.assert_frame_equal(df_output, df_expected, rtol=REL_TOL, atol=ABS_TOL, obj=file_name, check_exact=exact)


def get_header(path):
    with open(path, 'r') as file:
        first_line = file.readline().strip()
        return 0 if first_line.startswith('#') else None
