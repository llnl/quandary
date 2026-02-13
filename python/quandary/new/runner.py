"""Functional API for running Quandary simulations and optimizations.

MPI Lifecycle Management:
    MPI initialization and finalization are managed by mpi4py, which is imported
    in __init__.py. This allows multiple run() calls in the same Python process
    (e.g., Jupyter notebooks, parameter sweeps) without MPI re-initialization errors.

    The C++ runQuandary function is called with finalize=false to prevent it from
    calling MPI_Finalize, allowing mpi4py to manage finalization at process exit.
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys
from typing import TYPE_CHECKING, Optional

from mpi4py import MPI

from .. import _quandary_impl
from .._quandary_impl import Config, Setup
from .results import get_results as _get_results, Results

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


def _compute_optimal_core_distribution(maxcores: int, ninit: int) -> int:
    """Compute optimal MPI core distribution for Quandary.

    Quandary parallelizes over initial conditions (most efficient) and can also
    use PETSc parallelism within each initial condition. This function finds the
    best factorization: n_procs = ncores_init * ncores_petsc.

    Args:
        maxcores: Max number of MPI processes requested.
        ninit: Number of initial conditions.

    Returns:
        ncores: Actual total cores used (= ncores_init * ncores_petsc)

    """
    # Set default number of cores to the number of initial conditions, unless otherwise specified. Make sure ncores is an integer divisible of ninit.
    ncores_init = ninit
    if maxcores > -1:
        ncores_init = min(ninit, maxcores)
    for i in range(ninit, 0, -1):
        if ninit % i == 0:  # i is a factor of ninit
            if i <= ncores_init:
                ncores_init = i
                break
    # Set remaining number of cores for petsc
    ncores_petsc = 1
    if maxcores > ncores_init and maxcores % ncores_init == 0:
        ncores_petsc = int(maxcores / ncores_init)
    ncores = ncores_init * ncores_petsc

    # Log the distribution
    if ncores != maxcores:
        logger.info(
            f"Adjusted n_procs from {maxcores} to {ncores} for optimal parallelization across {ninit} initial conditions "
        )
    logger.info(
        f"MPI distribution: {ncores_init} cores for initial conditions × "
        f"{ncores_petsc} cores for PETSc = {ncores} total cores"
    )

    return ncores


def validate(setup: Setup, quiet: bool = False) -> Config:
    """Validate a Setup and return an immutable Config with all defaults applied.

    This allows inspecting computed defaults, checking for errors, and
    serializing to TOML without running a simulation.

    Args:
        setup: A Setup object with all required fields set.
        quiet: If True, suppress console output during validation.

    Returns:
        An immutable Config with all defaults computed and validation applied.

    Raises:
        RuntimeError: If the configuration is invalid.
    """
    return Config(setup, quiet)


def run(
    setup: Setup,
    max_n_procs: Optional[int] = None,
    quiet: bool = False,
    mpi_exec: str = "mpirun",
    nproc_flag: str = "-np",
    python_exec: Optional[str] = None,
    working_dir: str = ".",
) -> Results:
    """Run a Quandary simulation or optimization.

    This function validates the setup, runs the simulation/optimization,
    and returns results with the validated configuration attached.

    If n_procs is specified and not already in an MPI context, spawns a subprocess
    with the MPI launcher (useful for Jupyter notebooks). If already running in an
    MPI context (launched with mpirun), MPI will be used automatically. Otherwise
    runs serially.

    Args:
        setup: A Setup object.
        max_n_procs: Max number of MPI processes to use. If specified and not in MPI context,
            spawns subprocess. The actual number of cores will be determined based on the number of initial conditions.
            If None, runs directly (serial or existing MPI context).
        quiet: If True, suppress console output.
        mpi_exec: MPI launcher command (e.g., "mpirun", "srun", "flux run"). Default: "mpirun".
            Only used when spawning subprocess.
        nproc_flag: Flag for specifying number of processes (e.g., "-np", "-n"). Default: "-np".
            Only used when spawning subprocess.
        python_exec: Path to Python executable. Defaults to sys.executable.
            Only used when spawning subprocess.
        working_dir: Working directory for subprocess. Defaults to current directory.
            Only used when spawning subprocess.

    Returns:
        Results containing output data and the validated configuration.

    Raises:
        RuntimeError: If the configuration is invalid or execution fails.
        subprocess.CalledProcessError: If spawned subprocess fails.

    Examples:
        >>> setup = Setup()
        >>> setup.nlevels = [2]
        >>> setup.ntime = 1000
        >>> setup.dt = 0.01
        >>> setup.transition_frequency = [4.1]
        >>> setup.runtype = RunType.SIMULATION
        >>> results = run(setup)
        >>> print(f"Infidelity: {results.infidelity}")

        >>> # Spawn MPI subprocess (e.g., in Jupyter)
        >>> results = run(setup, n_procs=4)
    """
    # Check MPI context
    comm = MPI.COMM_WORLD
    size = comm.Get_size()

    # If max_n_procs specified, handle subprocess spawning
    if max_n_procs is not None:
        if size > 1:
            # Already in MPI context - can't spawn subprocess
            raise RuntimeError(
                f"Cannot spawn subprocess with max_n_procs={max_n_procs} - already running in MPI context with {size} processes.\n"
                f"Either:\n"
                f"  1. Remove max_n_procs parameter to use existing MPI context\n"
                f"  2. Run without mpirun and let run() spawn the subprocess"
            )
        # Not in MPI context (size == 1), spawn subprocess
        logger.info(f"Spawning subprocess with {max_n_procs} processes using {mpi_exec}")
        return _run_subprocess(
            setup=setup,
            max_n_procs=max_n_procs,
            quiet=quiet,
            mpi_exec=mpi_exec,
            nproc_flag=nproc_flag,
            python_exec=python_exec,
            working_dir=working_dir,
        )

    # Otherwise run directly (serial or using existing MPI context)
    if size > 1:
        logger.info(f"Running in existing MPI context with {size} processes")

    # Validate configuration
    validated_config = Config(setup, quiet)

    # Run simulation/optimization
    return_code = _quandary_impl.run(validated_config, quiet)

    if return_code != 0:
        raise RuntimeError(f"Quandary execution failed with return code {return_code}")

    # Load results with validated config
    results = _get_results(validated_config)

    return results


def _run_subprocess(
    setup: Setup,
    max_n_procs: int,
    quiet: bool = False,
    mpi_exec: str = "mpirun",
    nproc_flag: str = "-np",
    python_exec: Optional[str] = None,
    working_dir: str = ".",
) -> Results:
    """Internal: Spawn MPI subprocess to run Quandary.

    Called by run() when max_n_procs is specified. Writes config to TOML,
    spawns subprocess with MPI launcher, and returns results.

    The number of processes is automatically optimized based on the number of
    initial conditions for efficient parallelization.
    """
    # Validate configuration
    validated_config = Config(setup, quiet)

    # Compute optimal core distribution based on number of initial conditions
    n_init = validated_config.n_initial_conditions
    total_cores = _compute_optimal_core_distribution(max_n_procs, n_init)

    # Use current Python interpreter if not specified
    if python_exec is None:
        python_exec = sys.executable

    # Ensure working_dir is a string
    if not isinstance(working_dir, str):
        raise ValueError(f"working_dir must be a string, got {type(working_dir)}")

    # Create the output directory if it doesn't exist
    os.makedirs(validated_config.output_directory, exist_ok=True)

    # Write the TOML config file to the output directory (use absolute path for subprocess)
    config_file = os.path.abspath(os.path.join(validated_config.output_directory, "config.toml"))
    toml_content = validated_config.to_toml()
    with open(config_file, "w") as f:
        f.write(toml_content)

    # Python code to run Quandary from the TOML file
    python_code = f'from quandary.new import run_from_file; run_from_file("{config_file}", quiet={quiet})'

    # Build the command with optimized core count
    cmd = [mpi_exec, nproc_flag, str(total_cores), python_exec, "-c", python_code]

    # Run the subprocess
    result = subprocess.run(
        cmd,
        cwd=working_dir,
        capture_output=quiet,
        text=True,
    )

    if result.returncode != 0:
        if quiet and result.stderr:
            logger.error(f"Quandary failed: {result.stderr}")
        result.check_returncode()

    # Load results with validated config
    results = _get_results(validated_config)

    return results
