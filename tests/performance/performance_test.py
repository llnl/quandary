import os
import re
import subprocess
from typing import List

import pytest
from pydantic import BaseModel, TypeAdapter
from tests.utils.common import build_mpi_command

# Mark all tests in this file as performance tests
pytestmark = pytest.mark.performance

TEST_PATH = os.path.dirname(os.path.realpath(__file__))
TEST_CASES_PATH = os.path.join(TEST_PATH, "test_cases.json")
TEST_CONFIG_PATH = os.path.join(TEST_PATH, "configs")
QUANDARY_PATH = os.path.join(TEST_PATH, "..", "..", "quandary")


def pytest_benchmark_configure(config):
    config.option.benchmark_disable_gc = True
    config.option.benchmark_min_rounds = 1


class Case(BaseModel):
    simulation_name: str
    number_of_processes: List[int]
    repetitions: int


def load_test_cases():
    with open(TEST_CASES_PATH) as test_cases_file:
        ta = TypeAdapter(List[Case])
        test_cases = ta.validate_json(test_cases_file.read())
        return test_cases


TEST_CASES = load_test_cases()


def get_parametrize_values():
    """Generate a list of (simulation_name, nproc) tuples for parameterization."""
    params = []
    ids = []
    for test_case in TEST_CASES:
        for nproc in test_case.number_of_processes:
            params.append((test_case.simulation_name, nproc))
            ids.append(f"{test_case.simulation_name}_nproc_{nproc}")
    return params, ids


params, ids = get_parametrize_values()


def run_quandary(command):
    """Helper function to run quandary command and extract metrics"""
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    assert result.returncode == 0

    memory_match = re.search(r"Global Memory:\s+([\d.]+) MB", result.stdout)
    memory = float(memory_match.group(1)) if memory_match else 0
    return memory


@pytest.mark.parametrize(
    "simulation_name,number_of_processes",
    params,
    ids=ids
)
def test_eval(benchmark, simulation_name, number_of_processes, request):
    mpi_exec = request.config.getoption("--mpi-exec")
    mpi_opt = request.config.getoption("--mpi-opt")

    repetitions = next(
        (tc.repetitions for tc in TEST_CASES if tc.simulation_name == simulation_name),
        1
    )

    config_file = os.path.join(TEST_CONFIG_PATH, simulation_name + ".toml")

    command = build_mpi_command(
        mpi_exec=mpi_exec,
        num_processes=number_of_processes,
        mpi_opt=mpi_opt,
        quandary_path=QUANDARY_PATH,
        config_file=config_file)
    print(f"Running command: \"{' '.join(command)}\"")

    memory = benchmark.pedantic(
        run_quandary,
        args=(command,),
        rounds=repetitions,
        iterations=1
    )
    print(f"{memory=} MB")

    benchmark.extra_info['number_of_processors'] = number_of_processes
    benchmark.extra_info['memory_mb'] = memory
