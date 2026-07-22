"""Shared pytest fixtures for Python tests."""

import os

import pytest


@pytest.fixture(scope="session", autouse=True)
def quandary_on_path():
    """Make the repo-root quandary binary discoverable without manual PATH setup."""
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    os.environ["PATH"] = repo_root + os.pathsep + os.environ.get("PATH", "")


@pytest.fixture
def mpi_exec(request):
    """Get MPI executor and options from pytest options."""
    executor = request.config.getoption("--mpi-exec")
    options = request.config.getoption("--mpi-opt")
    if executor != "mpirun":
        return f"{executor} {options} -n"
    return f"mpirun {options} -np"
