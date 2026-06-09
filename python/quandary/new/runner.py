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
import subprocess
import sys
from typing import TYPE_CHECKING, Optional
from collections.abc import Sequence

from mpi4py import MPI
import numpy as np

from .. import _quandary_impl
from .._quandary_impl import Config, Setup, RunType, ControlType, ControlInitializationType
from .results import get_results as _get_results, Results
from .setup_helpers import set_target, set_initial_condition, resolve_output_dir, _get_output_dir
from .utils import fit_bspline0, fit_bspline2nd
from ._structs import (
    ControlParameterizationSettings,
    ControlInitializationSettings,
)

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


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

def _configure_run(
    setup: Setup,
    runtype: RunType,
    target: Optional[np.ndarray] = None,
    gate_rot_freq: Optional[Sequence[float]] = None,
    initial_condition=None,
    pcof=None,
    pt0=None,
    qt0=None,
    control_randomize: bool = True,
    control_amplitude: Optional[float] = None,
) -> Setup:
    """Return a copy of setup configured for the specified run type."""
    
    setup = setup.copy()
    setup.output_directory = resolve_output_dir(setup.output_directory)
    setup.runtype = runtype

    if target is not None:
        set_target(setup, target, gate_rot_freq=gate_rot_freq)
    if initial_condition is not None:
        set_initial_condition(setup, initial_condition=initial_condition)

    # Set up control initialization.
    # Priority: pt/qt > pcof > explicit amplitude > existing setup.control_initializations > default from C++ code.

    # If pt/qt are provide, fit pulses to bspline coefficients. 
    if pt0 is not None or qt0 is not None:
        if pcof is not None:
            raise ValueError("Cannot specify both pcof and pt0/qt0")
        if pt0 is None:
            pt0 = np.zeros_like(qt0)
        if qt0 is None:
            qt0 = np.zeros_like(pt0)
        
        # If control parameterization was not specified, choose Bspline 2nd order for better fitting quality.
        existing_control_params = setup.control_parameterizations or []
        if len(existing_control_params) == 0:
            control_type = ControlType.BSPLINE
        else:
            control_type = existing_control_params[0].control_type

        # Fit control parameters to either Bspline 0-th order or Bspline 2nd order. 
        if control_type == ControlType.BSPLINE0:
            nsplines = [max(2, setup.ntime + 1) for _ in range(len(setup.nessential))]
            pcof = fit_bspline0(
                pt0=pt0, qt0=qt0,
                nsplines=nsplines[0],
                spline_knot_spacing=setup.dt,
                ntime=setup.ntime, 
                dt=setup.dt,
                nessential=setup.nessential,
            )
            # Zero out carrier frequencies (pulses already include carrier)
            setup.carrier_frequencies = [[0.0] for _ in range(len(setup.nessential))]
        elif control_type == ControlType.BSPLINE:
            n_osc = len(setup.nessential)
            if len(existing_control_params) == 0:
                # Default to spline_knot_spacing of 3ns.
                spline_knot_spacing = 3.0
                computed_nspline = int(np.max([np.ceil(setup.ntime * setup.dt / spline_knot_spacing + 2), 5]))
                nsplines = [computed_nspline for _ in range(n_osc)]
            else:
                nsplines = [param.nspline for param in existing_control_params]
            pcof = fit_bspline2nd(0.0, 
                                  setup.ntime*setup.dt, 
                                  pt0, qt0, 
                                  nsplines, 
                                  carrier_frequencies=setup.carrier_frequencies,
                                  inputs_in_mhz=True )

        # Set control parameterization
        control_params = []
        for i in range(len(setup.nessential)):
            param = ControlParameterizationSettings()
            param.control_type = control_type
            param.nspline = nsplines[i]
            control_params.append(param)
        setup.control_parameterizations = control_params

    # If pcof is provided (either user input or set above by fitting splines to pt/qt), write coefficients to file and set control initialization to load from that file.
    if pcof is not None and len(pcof) > 0:
        output_dir = _get_output_dir(setup)
        os.makedirs(output_dir, exist_ok=True)
        pcof_file = os.path.join(output_dir, "pcof_init.dat")
        np.savetxt(pcof_file, pcof, fmt='%20.13e')

        control_inits = []
        for _ in range(len(setup.nessential)):
            init = ControlInitializationSettings()
            init.init_type = ControlInitializationType.FILE
            init.filename = pcof_file
            control_inits.append(init)
        setup.control_initializations = control_inits

    elif control_amplitude is not None:
        # Explicit amplitude — create uniform per-oscillator inits
        control_inits = []
        init_type = (
            ControlInitializationType.RANDOM if control_randomize
            else ControlInitializationType.CONSTANT
        )
        for _ in range(len(setup.nessential)):
            init = ControlInitializationSettings()
            init.init_type = init_type
            init.amplitude = control_amplitude
            control_inits.append(init)

        setup.control_initializations = control_inits

    return setup

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

    If max_n_procs is specified and not already in an MPI context, spawns a subprocess
    with the MPI launcher (useful for Jupyter notebooks). If already running in an
    MPI context (launched with mpirun), MPI will be used automatically. Otherwise
    runs serially.

    Parameters
    ----------
    setup : Setup
        A Setup object.
    max_n_procs : int, optional
        Max number of MPI processes to use. If specified and not in MPI
        context, spawns subprocess. The actual number of cores will be
        determined based on the number of initial conditions. If None,
        runs directly (serial or existing MPI context).
    quiet : bool
        Suppress console output. Default: False.
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
            f"Already running in MPI context with {size} processes. "
            f"Ignoring max_n_procs={max_n_procs} and using existing MPI context."
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
        stdout=subprocess.PIPE if quiet else None,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.returncode != 0:
        if result.stderr:
            raise RuntimeError(f"Quandary failed:\n{result.stderr.strip()}")
        raise RuntimeError("Quandary failed (see output above)")

    # Load results with validated config
    results = _get_results(validated_config)

    return results


def optimize(
    setup: Setup,
    pcof=None,
    pt0=None,
    qt0=None,
    control_randomize: bool = True,
    control_amplitude: Optional[float] = None,
    initial_condition=None,
    target=None,
    gate_rot_freq=None,
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

    Parameters
    ----------
    setup : Setup
        Physics setup from setup_quandary().
    pcof : array-like, optional
        For Warm-start: Initial B-spline coefficients. Setup must contain the same spline parameterization and carrier frequencies.
    pt0 : sequence of ndarray, optional
        For warm-start: Real part of control pulses [MHz] per oscillator.
        Will be fitted to B-splines coefficients. Must be paired with qt0.
    qt0 : sequence of ndarray, optional
        For warm-start: Imaginary part of control pulses [MHz] per oscillator.
    control_randomize : bool
        Initialize controls randomly. Default: True.
    control_amplitude : float, optional
        Initial control amplitude [GHz]. When omitted, uses
        setup.control_initializations if set, otherwise defaults from C++ code (zero controls)
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
        Suppress output. Default: False.
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
    configured = _configure_run(setup,
        runtype = RunType.OPTIMIZATION,
        pcof=pcof,
        pt0=pt0,
        qt0=qt0,
        target=target,
        gate_rot_freq=gate_rot_freq,
        initial_condition=initial_condition,
        control_randomize=control_randomize,
        control_amplitude=control_amplitude,
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
    control_randomize: bool = True,
    control_amplitude: Optional[float] = None,
    initial_condition=None,
    target=None,
    gate_rot_freq=None,
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

    Parameters
    ----------
    setup : Setup
        Physics setup from setup_quandary().
    pcof : array-like, optional
        B-spline control coefficients.
    pt0 : sequence of ndarray, optional
        Real part of control pulses [MHz] per oscillator.
        Fitted to B-splines coefficients. Must be paired with qt0.
    qt0 : sequence of ndarray, optional
        Imaginary part of control pulses [MHz] per oscillator.
    control_randomize : bool
        Initialize controls randomly. Default: True.
    control_amplitude : float, optional
        Initial control amplitude [GHz]. When omitted, uses
        setup.control_initializations if set, otherwise defaults from C++ code (zero controls)
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
        Suppress output. Default: False.
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
    configured = _configure_run(setup, 
        runtype=RunType.SIMULATION, 
        pcof=pcof, 
        pt0=pt0, 
        qt0=qt0, 
        target=target, 
        gate_rot_freq=gate_rot_freq,
        initial_condition=initial_condition, 
        control_randomize=control_randomize, 
        control_amplitude=control_amplitude, 
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


def evaluate_controls(
    setup: Setup,
    pcof=None,
    pt0=None,
    qt0=None,
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

    Parameters
    ----------
    setup : Setup
        Physics setup from setup_quandary().
    pcof : array-like
        B-spline control coefficients.
    pt0 : sequence of ndarray, optional
        Real part of control pulses [MHz] per oscillator.
        Will be fitted to B-splines coefficients. Must be paired with qt0.
    qt0 : sequence of ndarray, optional
        Imaginary part of control pulses [MHz] per oscillator.
    points_per_ns : float
        Sample rate [points per ns]. Default: 1.0.
    dry_run : bool
        If True, validate and return Results with config populated but do not
        run. Use ``print(results.config)`` to inspect the full configuration.
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

    Returns
    -------
    Results
        If dry_run=True, only results.config is populated.
    """

    # Either pcof0 or pt0/qt0 must be provided to evaluate controls
    if pcof is None and (pt0 is None or qt0 is None):
        raise ValueError("Must provide either pcof or both pt0 and qt0 to evaluate control pulses.")

    # Configure the run with provided setup
    configured = _configure_run(setup, 
        pcof=pcof, 
        pt0=pt0, 
        qt0=qt0,
        runtype = RunType.EVALCONTROLS
    )

    # Recalculate time grid: keep total time, change resolution
    total_time = configured.ntime * configured.dt
    configured.ntime = int(np.floor(total_time * points_per_ns))
    configured.dt = total_time / configured.ntime

    # Run or dry run, return Results struct
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
