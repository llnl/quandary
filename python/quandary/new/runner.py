"""Functional API for running Quandary simulations and optimizations.

MPI Lifecycle Management:
    MPI initialization and finalization are managed by mpi4py, which is imported
    in __init__.py. This allows multiple optimize()/simulate() calls in the same Python process
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
from .setup_helpers import (
    _setup_optimization,
    _setup_simulation,
    _setup_eval_controls,
)

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


def _is_interactive():
    """Detect if running in an interactive environment (e.g., Jupyter notebook)."""
    import __main__ as main
    return not hasattr(main, '__file__')


def _compute_optimal_core_distribution(maxcores: int, ninit: int) -> int:
    """Compute optimal MPI core distribution for Quandary.

    Parallelization over initial conditions is the most efficient strategy.
    This function finds the largest divisor of ninit that fits within maxcores.
    PETSc parallelism within each initial condition is not used here — for HPC
    runs requiring it, launch with mpirun directly.

    Args:
        maxcores: Max number of MPI processes requested.
        ninit: Number of initial conditions.

    Returns:
        ncores: Actual total cores used (<= min(maxcores, ninit))

    """
    ncores = min(ninit, maxcores) if maxcores > -1 else ninit
    # Round down to nearest divisor of ninit
    for i in range(ncores, 0, -1):
        if ninit % i == 0:
            ncores = i
            break

    if ncores != maxcores:
        logger.info(
            f"Adjusted n_procs from {maxcores} to {ncores} for optimal parallelization across {ninit} initial conditions"
        )

    return ncores



def _run(
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

    if _is_interactive():
        # In interactive environment without MPI, spawn subprocess with MPI launcher for better performance
        if max_n_procs is None:
            max_n_procs = 4  # Default to 4 processes if not specified
        return _run_subprocess(
            setup=setup,
            max_n_procs=max_n_procs,
            quiet=quiet,
            mpi_exec=mpi_exec,
            nproc_flag=nproc_flag,
            python_exec=python_exec,
            working_dir=working_dir,
        )

    if max_n_procs is not None:
        logger.warning(
            f"Already running in MPI context with {size} processes. Ignoring max_n_procs={max_n_procs} and using existing MPI context."
        )

    logger.info(f"Running directly in existing MPI context with {size} processes")
    return _run_directly(setup, quiet)


def _run_directly(setup: Setup, quiet: bool = False) -> Results:
    """Internal: Run Quandary directly without spawning subprocess.

    This is used when max_n_procs is not specified, or when already in an MPI context.
    """
    # Validate configuration
    validated_config = Config(setup, quiet)
    logger.debug("Configuration:\n%s", validated_config)

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
    logger.debug("Configuration:\n%s", validated_config)

    # Compute optimal core distribution based on number of initial conditions
    n_init = validated_config.n_initial_conditions
    total_cores = _compute_optimal_core_distribution(max_n_procs, n_init)

    # Use current Python interpreter if not specified
    if python_exec is None:
        python_exec = sys.executable

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
    logger.info(f"Spawning subprocess with {total_cores} processes using {mpi_exec}")
    logger.debug(f"Subprocess command: {' '.join(cmd)}")
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


def optimize(
    setup: Setup,
    targetgate=None,
    targetstate=None,
    target_levels=None,
    gate_rot_freq=None,
    pcof=None,
    randomize_init_ctrl: bool = False,
    init_amplitude_ghz=None,
    dry_run: bool = False,
    max_n_procs: Optional[int] = None,
    quiet: bool = False,
    mpi_exec: str = "mpirun",
    nproc_flag: str = "-np",
    python_exec: Optional[str] = None,
    working_dir: str = ".",
) -> Results:
    """Run an optimization.

    Copies the setup internally; the original is not modified.

    Parameters:
    ----------
    setup : Setup
        Physics setup from setup_physics(). Not modified.
    targetgate : array-like, optional
        Target unitary gate.
    targetstate : array-like, optional
        Target state vector.
    target_levels : List[int], optional
        Target product state, e.g. [0, 0, 1].
    gate_rot_freq : List[float], optional
        Gate rotation frequencies [GHz].
    pcof : array-like, optional
        B-spline coefficients for warm-start.
    randomize_init_ctrl : bool
        Initialize controls randomly. Default: False.
    init_amplitude_ghz : float, optional
        Initial control amplitude [GHz].
    dry_run : bool
        If True, validate and return Results with config populated but do not
        run. Use print(results.config) to inspect the full configuration.
        Default: False.
    max_n_procs : int, optional
        Max MPI processes. Spawns subprocess if set.
    quiet : bool
        Suppress output. Default: False.
    mpi_exec : str
        MPI launcher. Default: "mpirun".
    nproc_flag : str
        Flag for process count. Default: "-np".
    python_exec : str, optional
        Python executable path. Default: sys.executable.
    working_dir : str
        Working directory. Default: ".".

    Returns:
    -------
    Results
        If dry_run=True, only results.config is populated.
    """
    configured = _setup_optimization(
        setup,
        targetgate=targetgate,
        targetstate=targetstate,
        target_levels=target_levels,
        gate_rot_freq=gate_rot_freq,
        pcof=pcof,
        randomize_init_ctrl=randomize_init_ctrl,
        init_amplitude_ghz=init_amplitude_ghz,
    )
    if dry_run:
        return Results(config=Config(configured, quiet))
    return _run(
        configured,
        max_n_procs=max_n_procs,
        quiet=quiet,
        mpi_exec=mpi_exec,
        nproc_flag=nproc_flag,
        python_exec=python_exec,
        working_dir=working_dir,
    )


def simulate(
    setup: Setup,
    pcof=None,
    pt0=None,
    qt0=None,
    dry_run: bool = False,
    max_n_procs: Optional[int] = None,
    quiet: bool = False,
    mpi_exec: str = "mpirun",
    nproc_flag: str = "-np",
    python_exec: Optional[str] = None,
    working_dir: str = ".",
) -> Results:
    """Run a simulation.

    Copies the setup internally; the original is not modified.

    Parameters:
    ----------
    setup : Setup
        Physics setup from setup_physics(). Not modified.
    pcof : array-like, optional
        B-spline control coefficients.
    pt0 : list of ndarray, optional
        Real part of control pulses [MHz] per oscillator.
        Auto-downsampled to B-splines. Must be paired with qt0.
    qt0 : list of ndarray, optional
        Imaginary part of control pulses [MHz] per oscillator.
    dry_run : bool
        If True, validate and return Results with config populated but do not
        run. Use print(results.config) to inspect the full configuration.
        Default: False.
    max_n_procs : int, optional
        Max MPI processes. Spawns subprocess if set.
    quiet : bool
        Suppress output. Default: False.
    mpi_exec : str
        MPI launcher. Default: "mpirun".
    nproc_flag : str
        Flag for process count. Default: "-np".
    python_exec : str, optional
        Python executable path. Default: sys.executable.
    working_dir : str
        Working directory. Default: ".".

    Returns:
    -------
    Results
        If dry_run=True, only results.config is populated.
    """
    configured = _setup_simulation(setup, pcof=pcof, pt0=pt0, qt0=qt0)
    if dry_run:
        return Results(config=Config(configured, quiet))
    return _run(
        configured,
        max_n_procs=max_n_procs,
        quiet=quiet,
        mpi_exec=mpi_exec,
        nproc_flag=nproc_flag,
        python_exec=python_exec,
        working_dir=working_dir,
    )


def evaluate_controls(
    setup: Setup,
    pcof,
    points_per_ns: float = 1.0,
    dry_run: bool = False,
    max_n_procs: Optional[int] = None,
    quiet: bool = False,
    mpi_exec: str = "mpirun",
    nproc_flag: str = "-np",
    python_exec: Optional[str] = None,
    working_dir: str = ".",
) -> Results:
    """Evaluate control pulses at a specific sample rate.

    Copies the setup internally; the original is not modified.

    Parameters:
    ----------
    setup : Setup
        Physics setup from setup_physics(). Not modified.
    pcof : array-like
        B-spline control coefficients.
    points_per_ns : float
        Sample rate [points per ns]. Default: 1.0.
    dry_run : bool
        If True, validate and return Results with config populated but do not
        run. Use print(results.config) to inspect the full configuration.
        Default: False.
    max_n_procs : int, optional
        Max MPI processes. Spawns subprocess if set.
    quiet : bool
        Suppress output. Default: False.
    mpi_exec : str
        MPI launcher. Default: "mpirun".
    nproc_flag : str
        Flag for process count. Default: "-np".
    python_exec : str, optional
        Python executable path. Default: sys.executable.
    working_dir : str
        Working directory. Default: ".".

    Returns:
    -------
    Results
        If dry_run=True, only results.config is populated.
    """
    configured = _setup_eval_controls(setup, pcof=pcof, points_per_ns=points_per_ns)
    if dry_run:
        return Results(config=Config(configured, quiet))
    return _run(
        configured,
        max_n_procs=max_n_procs,
        quiet=quiet,
        mpi_exec=mpi_exec,
        nproc_flag=nproc_flag,
        python_exec=python_exec,
        working_dir=working_dir,
    )
