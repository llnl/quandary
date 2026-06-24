"""Functional API for running Quandary simulations and optimizations.

MPI Lifecycle Management:
    MPI initialization and finalization are managed by mpi4py, which is imported
    in __init__.py. This allows multiple optimize()/simulate() calls in the same Python process
    (e.g., Jupyter notebooks, parameter sweeps) without MPI re-initialization errors.

    The C++ runQuandary function detects that MPI is already initialized (by mpi4py)
    and skips MPI_Init/MPI_Finalize, letting mpi4py manage the MPI lifecycle.
"""
from __future__ import annotations
import logging
import os
import shlex
import subprocess
import sys
from typing import Optional
from mpi4py import MPI
import numpy as np
from .. import _quandary_impl
from .._quandary_impl import Config, ConfigInput, RunType
from .results import get_results, Results
from .config import resolve_output_dir, set_controls, set_target, set_initial_condition

logger = logging.getLogger(__name__)

def simulate(
    config_input: ConfigInput,
    spline_coefficients=None,
    p_samples=None,
    q_samples=None,
    control_randomize: bool = True,
    control_amplitude: Optional[float] = None,
    initial_condition=None,
    target=None,
    gate_rot_freq=None,
    dry_run: bool = False,
    max_n_procs: Optional[int] = None,
    quiet: bool = True,
    mpi_exec: str = "mpirun",
    nproc_flag: str = "-np",
    python_exec: Optional[str] = None,
    working_dir: str = ".",
) -> Results:
    """Run a simulation.

    Copies the config_input internally; the original is not modified.

    Parameters
    ----------
    config_input : ConfigInput
        Physics config_input from create_config().
    spline_coefficients : array-like, optional
        B-spline control coefficients.
    p_samples : sequence of ndarray, optional
        Real part of control pulses [MHz] per oscillator.
        Fitted to B-splines coefficients. Must be paired with q_samples.
    q_samples : sequence of ndarray, optional
        Imaginary part of control pulses [MHz] per oscillator.
    control_randomize : bool
        Initialize controls randomly. Default: True.
    control_amplitude : float, optional
        Initial control amplitude [GHz]. When omitted, uses
        config_input.control_initializations if set, otherwise defaults from C++ code (zero controls)
    initial_condition : sequence of complex or InitialConditionSettings, optional
        Either a state vector (arbitrary superposition), or direct struct specification (advanced). 
        Default: All basis states in the essential dimensions.
    target : array-like, optional 
        Optimization target. Either 2D (unitary gate) or 1D (state vector)
    gate_rot_freq : sequence of float, optional
        Gate rotation frequencies [GHz].
    dry_run : bool
        If True, validate and return Results with config populated but do not
        run. Use ``print(results.config)`` to inspect the full configuration. Default: False.
    max_n_procs : int, optional
        Max MPI processes. Spawns subprocess if set.
    quiet : bool
        Suppress output. Default: True.
    mpi_exec : str
        MPI launcher. Default: "mpirun".
    nproc_flag : str
        Flag for process count. Default: "-np".
    python_exec : str, optional
        Python executable path. Default: sys.executable.
    working_dir : str
        Working directory. Default: ".".

    Returns
    -------
    Results
        If dry_run=True, only results.config is populated.
    """

    config_input = config_input.copy()
    config_input.output_directory = resolve_output_dir(config_input.output_directory)
    config_input.runtype = RunType.SIMULATION

    set_controls(config_input, spline_coefficients=spline_coefficients, p_samples=p_samples, q_samples=q_samples, control_randomize=control_randomize, control_amplitude=control_amplitude)
    set_target(config_input, target, gate_rot_freq=gate_rot_freq)
    set_initial_condition(config_input, initial_condition=initial_condition)

    if dry_run:
        return Results(config=Config(config_input, quiet))
    return _run(
        config_input,
        max_n_procs=max_n_procs,
        quiet=quiet,
        mpi_exec=mpi_exec,
        nproc_flag=nproc_flag,
        python_exec=python_exec,
        working_dir=working_dir,
    )

def optimize(
    config_input: ConfigInput,
    spline_coefficients=None,
    p_samples=None,
    q_samples=None,
    control_randomize: bool = True,
    control_amplitude: Optional[float] = None,
    initial_condition=None,
    target=None,
    gate_rot_freq=None,
    dry_run: bool = False,
    max_n_procs: Optional[int] = None,
    quiet: bool = True,
    mpi_exec: str = "mpirun",
    nproc_flag: str = "-np",
    python_exec: Optional[str] = None,
    working_dir: str = ".",
) -> Results:
    """Run an optimization.

    Copies the config_input internally; the original is not modified.

    Parameters
    ----------
    config_input : ConfigInput
        Physics config_input from create_config().
    spline_coefficients : array-like, optional
        For Warm-start: Initial B-spline coefficients. ConfigInput must contain the same spline parameterization and carrier frequencies.
    p_samples : sequence of ndarray, optional
        For warm-start: Real part of control pulses [MHz] per oscillator.
        Will be fitted to B-splines coefficients. Must be paired with q_samples.
    q_samples : sequence of ndarray, optional
        For warm-start: Imaginary part of control pulses [MHz] per oscillator.
    control_randomize : bool
        Initialize controls randomly. Default: True.
    control_amplitude : float, optional
        Initial control amplitude [GHz]. When omitted, uses
        config_input.control_initializations if set, otherwise defaults from C++ code (zero controls)
    initial_condition : sequence of complex or InitialConditionSettings, optional
        Either a state vector (arbitrary superposition), or direct struct specification (advanced). 
        Default: All basis states in the essential dimensions.
    target : array-like, optional 
        Optimization target. Either 2D (unitary gate) or 1D (state vector)
    gate_rot_freq : sequence of float, optional
        Gate rotation frequencies [GHz].
    dry_run : bool
        If True, validate and return Results with config populated but do not
        run. Use ``print(results.config)`` to inspect the full configuration. Default: False.
    max_n_procs : int, optional
        Max MPI processes. Spawns subprocess if set.
    quiet : bool
        Suppress output. Default: True.
    mpi_exec : str
        MPI launcher. Default: "mpirun".
    nproc_flag : str
        Flag for process count. Default: "-np".
    python_exec : str, optional
        Python executable path. Default: sys.executable.
    working_dir : str
        Working directory. Default: ".".

    Returns
    -------
    Results
        If dry_run=True, only results.config is populated.
    """

    config_input = config_input.copy()
    config_input.output_directory = resolve_output_dir(config_input.output_directory)
    config_input.runtype = RunType.OPTIMIZATION
    set_controls(config_input, spline_coefficients=spline_coefficients, p_samples=p_samples, q_samples=q_samples, control_randomize=control_randomize, control_amplitude=control_amplitude)
    set_target(config_input, target, gate_rot_freq=gate_rot_freq)
    set_initial_condition(config_input, initial_condition=initial_condition)

    if dry_run:
        return Results(config=Config(config_input, quiet))
    return _run(
        config_input,
        max_n_procs=max_n_procs,
        quiet=quiet,
        mpi_exec=mpi_exec,
        nproc_flag=nproc_flag,
        python_exec=python_exec,
        working_dir=working_dir,
    )

def evaluate_controls(
    config_input: ConfigInput,
    spline_coefficients=None,
    p_samples=None,
    q_samples=None,
    control_randomize: bool = True,
    control_amplitude: Optional[float] = None,
    points_per_ns: float = 1.0,
    dry_run: bool = False,
    max_n_procs: Optional[int] = None,
    quiet: bool = True,
    mpi_exec: str = "mpirun",
    nproc_flag: str = "-np",
    python_exec: Optional[str] = None,
    working_dir: str = ".",
) -> Results:
    """Evaluate control pulses at a specific sample rate.

    Copies the config_input internally; the original is not modified.

    Parameters
    ----------
    config_input : ConfigInput
        Physics config_input from create_config().
    spline_coefficients : array-like
        B-spline control coefficients.
    p_samples : sequence of ndarray, optional
        Real part of control pulses [MHz] per oscillator.
        Will be fitted to B-splines coefficients. Must be paired with q_samples.
    q_samples : sequence of ndarray, optional
        Imaginary part of control pulses [MHz] per oscillator.
    control_randomize : bool
        Initialize controls randomly. Default: True.
    control_amplitude : float, optional
        Initial control amplitude [GHz]. When omitted, uses
        config_input.control_initializations if set, otherwise defaults from C++ code (zero controls)
    points_per_ns : float
        Sample rate [points per ns]. Default: 1.0.
    dry_run : bool
        If True, validate and return Results with config populated but do not
        run. Use ``print(results.config)`` to inspect the full configuration.
        Default: False.
    max_n_procs : int, optional
        Max MPI processes. Spawns subprocess if set.
    quiet : bool
        Suppress output. Default: True.
    mpi_exec : str
        MPI launcher. Default: "mpirun".
    nproc_flag : str
        Flag for process count. Default: "-np".
    python_exec : str, optional
        Python executable path. Default: sys.executable.
    working_dir : str
        Working directory. Default: ".".

    Returns
    -------
    Results
        If dry_run=True, only results.config is populated.
    """

    config_input = config_input.copy()
    config_input.output_directory = resolve_output_dir(config_input.output_directory)
    config_input.runtype = RunType.EVALCONTROLS

    set_controls(config_input, spline_coefficients=spline_coefficients, p_samples=p_samples, q_samples=q_samples, control_randomize=control_randomize, control_amplitude=control_amplitude)

    # Recalculate time grid: keep total time, change resolution
    total_time = config_input.total_time
    nsteps = int(np.floor(total_time * points_per_ns))
    nsteps = max(nsteps, 1)  # Ensure at least one step
    config_input.dt = total_time / nsteps

    # Run or dry run, return Results struct
    if dry_run:
        return Results(config=Config(config_input, quiet))
    return _run(
        config_input,
        max_n_procs=max_n_procs,
        quiet=quiet,
        mpi_exec=mpi_exec,
        nproc_flag=nproc_flag,
        python_exec=python_exec,
        working_dir=working_dir,
    )


# ---------------------------------
# Private run helpers
# ---------------------------------

def _is_interactive():
    """Detect if running in an interactive Python session (Jupyter, IPython, or plain REPL)."""
    # Check for IPython/Jupyter
    try:
        from IPython import get_ipython
        shell = get_ipython()
        if shell is not None and shell.__class__.__name__ in (
            'ZMQInteractiveShell',       # Jupyter notebook/lab
            'TerminalInteractiveShell',  # IPython terminal
        ):
            return True
    except ImportError:
        pass
    # Check for plain Python REPL
    import __main__
    return not hasattr(__main__, '__file__')


def _compute_optimal_core_distribution(maxcores: int, ninit: int) -> int:
    """Compute optimal MPI core distribution for Quandary.

    Parallelization over initial conditions is the most efficient strategy.
    This function finds the largest divisor of ninit that fits within maxcores.
    PETSc parallelism within each initial condition is not used here -- for HPC
    runs requiring it, launch with mpirun directly.

    Parameters
    ----------
    maxcores : int
        Max number of MPI processes requested.
    ninit : int
        Number of initial conditions.

    Returns
    -------
    ncores : int
        Actual total cores used (<= min(maxcores, ninit)).
    """
    ncores = min(ninit, maxcores) if maxcores > -1 else ninit
    # Round down to nearest divisor of ninit
    for i in range(ncores, 0, -1):
        if ninit % i == 0:
            ncores = i
            break

    if ncores != maxcores:
        logger.info(
            f"Adjusted n_procs from {maxcores} to {ncores} for optimal parallelization "
            f"across {ninit} initial conditions"
        )

    return ncores

def _run(
    config_input: ConfigInput,
    max_n_procs: Optional[int] = None,
    quiet: bool = True,
    mpi_exec: str = "mpirun",
    nproc_flag: str = "-np",
    python_exec: Optional[str] = None,
    working_dir: str = ".",
) -> Results:
    """Run a Quandary simulation or optimization.

    This function validates the config_input, runs the simulation/optimization,
    and returns results with the validated configuration attached.

    If max_n_procs is specified and not already in an MPI context, spawns a subprocess
    with the MPI launcher (useful for Jupyter notebooks). If already running in an
    MPI context (launched with mpirun), MPI will be used automatically. Otherwise
    runs serially.

    Parameters
    ----------
    config_input : ConfigInput
        A ConfigInput object.
    max_n_procs : int, optional
        Max number of MPI processes to use. If specified and not in MPI
        context, spawns subprocess. The actual number of cores will be
        determined based on the number of initial conditions. If None,
        runs directly (serial or existing MPI context).
    quiet : bool
        Suppress console output. Default: True.
    mpi_exec : str
        MPI launcher command (e.g., "mpirun", "srun", "flux run").
        Default: "mpirun". Only used when spawning subprocess.
    nproc_flag : str
        Flag for specifying number of processes (e.g., "-np", "-n").
        Default: "-np". Only used when spawning subprocess.
    python_exec : str, optional
        Path to Python executable. Default: sys.executable.
        Only used when spawning subprocess.
    working_dir : str
        Working directory for subprocess. Default: ".".
        Only used when spawning subprocess.

    Returns
    -------
    Results
        Output data and the validated configuration.

    Raises
    ------
    RuntimeError
        If the configuration is invalid or execution fails.
    subprocess.CalledProcessError
        If spawned subprocess fails.
    """
    # Check MPI context
    comm = MPI.COMM_WORLD
    size = comm.Get_size()

    # if _is_interactive():
    if True:
        # In interactive environment without MPI, spawn subprocess with MPI launcher for better performance
        if max_n_procs is None:
            #max_n_procs = 4  # Default to 4 processes if not specified
            max_n_procs = 2  # Default to 4 processes if not specified
        return _run_subprocess(
            config_input=config_input,
            max_n_procs=max_n_procs,
            quiet=quiet,
            mpi_exec=mpi_exec,
            nproc_flag=nproc_flag,
            python_exec=python_exec,
            working_dir=working_dir,
        )

    if max_n_procs is not None:
        logger.warning(
            f"Already running in MPI context with {size} processes. "
            f"Ignoring max_n_procs={max_n_procs} and using existing MPI context."
        )

    logger.info(f"Running directly in existing MPI context with {size} processes")
    return _run_directly(config_input, quiet)


def _run_directly(config_input: ConfigInput, quiet: bool = True) -> Results:
    """Internal: Run Quandary directly without spawning subprocess.

    This is used when max_n_procs is not specified, or when already in an MPI context.
    """
    # Validate configuration
    validated_config = Config(config_input, quiet)
    logger.debug("Configuration:\n%s", validated_config)

    # Run simulation/optimization
    return_code = _quandary_impl.run(validated_config, quiet)

    if return_code != 0:
        raise RuntimeError(f"Quandary execution failed with return code {return_code}")

    # Load results with validated config
    results = get_results(validated_config.output_directory)

    return results


def _run_subprocess(
    config_input: ConfigInput,
    max_n_procs: int,
    quiet: bool = True,
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
    validated_config = Config(config_input, quiet)
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

    # Python code to run Quandary from the TOML file.
    # Pass dynamic values through argv to avoid launcher/shell quoting issues.
    python_code = "import sys; from quandary.new import run_from_file; run_from_file(sys.argv[1], quiet=bool(int(sys.argv[2])))"
    quiet_flag = "1" if quiet else "0"

    # Build command: bypass MPI launcher for single-process execution.
    # Some remote environments restrict mpirun even for -np 1.
    if total_cores <= 1:
        cmd = [python_exec, "-c", python_code, config_file, quiet_flag]
    else:
        cmd = [
            mpi_exec,
            nproc_flag,
            str(total_cores),
            python_exec,
            "-c",
            python_code,
            config_file,
            quiet_flag,
        ]

    # Run the subprocess
    if total_cores <= 1:
        logger.info("Spawning single-process subprocess without MPI launcher")
    else:
        logger.info(f"Spawning subprocess with {total_cores} processes using {mpi_exec}")
    logger.debug(f"Subprocess command: {shlex.join(cmd)}")
    capture_output = quiet
    result = subprocess.run(
        cmd,
        cwd=working_dir,
        stdout=subprocess.PIPE if capture_output else None,
        stderr=subprocess.PIPE if capture_output else None,
        text=True,
    )

    if result.returncode != 0:
        message_lines = [
            f"Quandary failed with return code {result.returncode}.",
            f"Command: {shlex.join(cmd)}",
        ]
        if capture_output and result.stderr and result.stderr.strip():
            message_lines.append(f"stderr:\n{result.stderr.strip()}")
        if capture_output and result.stdout and result.stdout.strip():
            message_lines.append(f"stdout:\n{result.stdout.strip()}")
        if not capture_output:
            message_lines.append(
                "Subprocess output was streamed directly to terminal (quiet=False); see logs above for details."
            )
        raise RuntimeError("\n\n".join(message_lines))

    # Load results with validated config
    results = get_results(validated_config.output_directory)

    return results


