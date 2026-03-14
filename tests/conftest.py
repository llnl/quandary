"""Common pytest configuration."""


def pytest_addoption(parser):
    """Add common command line options to pytest."""
    parser.addoption(
        "--mpi-exec",
        action="store",
        default="mpirun",
        help="Path to the MPI executable (e.g., mpirun or srun)"
    )

    parser.addoption(
        "--mpi-opt",
        action="store",
        default="",
        help="Extra options to pass to mpi exec command)"
    )

    parser.addoption(
        "--petsc-options",
        action="store",
        default="",
        help="PETSc options to pass to quandary (e.g., '-vec_type kokkos -mat_type aijkokkos')"
    )
