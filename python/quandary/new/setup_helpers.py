"""Helper functions to create and configure Setup objects for simulation and optimization."""

import logging
import os
from typing import List, Optional

import numpy as np

from .._quandary_impl import (
    RunType,
    ControlType,
    InitialConditionType,
    TargetType,
    GateType,
    ControlInitializationType,
    OutputType,
)
from ._structs import Setup
from .._quandary_impl import (
    InitialConditionSettings,
    OptimTargetSettings,
    ControlParameterizationSettings,
    ControlInitializationSettings,
)
from .quantum_operators import hamiltonians, get_resonances
from .time_estimation import estimate_timesteps
from .utils import downsample_pulses

logger = logging.getLogger(__name__)

_DEFAULT_OUTPUT_DIR = "./run_dir"


def resolve_output_dir(datadir: str) -> str:
    """Resolve the output directory using the QUANDARY_BASE_DATADIR environment variable.

    Parameters
    ----------
    datadir : str
        Output directory path. If QUANDARY_BASE_DATADIR is set, relative paths
        will be resolved against it.

    Returns
    -------
    str
        Resolved absolute or relative path for the data directory.

    Raises
    ------
    ValueError
        If QUANDARY_BASE_DATADIR is set but doesn't exist or isn't a directory.
    """
    if os.path.isabs(datadir):
        return datadir

    base_dir = os.environ.get("QUANDARY_BASE_DATADIR")
    if base_dir:
        if not os.path.exists(base_dir):
            raise ValueError(f"Environment variable QUANDARY_BASE_DATADIR points to non-existent path: {base_dir}")
        if not os.path.isdir(base_dir):
            raise ValueError(f"Environment variable QUANDARY_BASE_DATADIR is not a directory: {base_dir}")

        datadir = os.path.join(base_dir, datadir)

    return os.path.normpath(datadir)


def _get_output_dir(setup):
    """Get output directory from setup, falling back to default."""
    return setup.output_directory or _DEFAULT_OUTPUT_DIR


def setup_quandary(
    nessential: List[int],
    transition_frequency: List[float],
    final_time: float,
    ntime: Optional[int] = None,
    dt: Optional[float] = None,
    selfkerr: Optional[List[float]] = None,
    nguard: Optional[List[int]] = None,
    rotation_frequency: Optional[List[float]] = None,
    crosskerr_coupling: Optional[List[float]] = None,
    dipole_coupling: Optional[List[float]] = None,
    carrier_frequency: Optional[List[List[float]]] = None,
    Pmin: int = 150,
    control_amplitude_bounds: Optional[List[float]] = None,
    nspline: Optional[int] = None,
    spline_knot_spacing: Optional[float] = None,
    spline_order: Optional[int] = None,
    control_zero_boundary_condition: Optional[bool] = None,
    initial_condition: Optional[InitialConditionSettings] = None,
    initial_levels: Optional[List[int]] = None,
    initial_state: Optional[List[complex]] = None,
    output_directory: str = _DEFAULT_OUTPUT_DIR,
    verbose: bool = True,
) -> Setup:
    """Create a Setup with physics parameters configured.

    Automatically computes Hamiltonians, timesteps, and carrier frequencies.
    Does NOT set the runtype -- use optimize() or simulate() to configure for
    a specific run type.

    **Time discretization:** specify at most 2 of (final_time, ntime, dt). If
    all 3 are given they must satisfy ``final_time == ntime * dt``. When ntime
    and dt are both omitted, ntime is auto-computed from the Hamiltonian
    eigenvalues.

    **Spline configuration:** specify at most one of (nspline, spline_knot_spacing).
    If neither is given the C++ default (10 splines) is used.

    Parameters
    ----------
    nessential : list of int
        Number of essential energy levels per qubit.
    transition_frequency : list of float
        01-transition frequencies [GHz] per qubit.
    final_time : float
        Pulse duration [ns].
    ntime : int, optional
        Number of timesteps. Auto-computed from the Hamiltonian if omitted.
    dt : float, optional
        Timestep size [ns]. Auto-computed if omitted.
    selfkerr : list of float, optional
        Anharmonicities [GHz] per qubit. Default: zeros.
    nguard : list of int, optional
        Number of guard levels per qubit. Default: zeros.
    rotation_frequency : list of float, optional
        Rotating frame frequencies [GHz] per qubit.
        Default: transition_frequency.
    crosskerr_coupling : list of float, optional
        ZZ coupling strengths [GHz]. Format: [g01, g02, ..., g12, g13, ...].
        Default: no coupling.
    dipole_coupling : list of float, optional
        Dipole-dipole coupling strengths [GHz]. Format: [J01, J02, ..., J12, ...].
        Default: no coupling.
    carrier_frequency : list of list of float, optional
        Carrier frequencies [GHz] per oscillator. When provided, the
        automatic eigenvalue-based computation (get_resonances) is skipped.
        Useful when the Hamiltonian has degenerate eigenvalues that cause
        the automatic computation to fail.
    Pmin : int
        Minimum time steps per period of the fastest oscillation (used for
        auto-computing ntime). Default: 150.
    control_amplitude_bounds : list of float, optional
        Maximum control amplitudes [GHz] per qubit. Used for timestep
        estimation and as optimization bounds. Default: 0.01 GHz (= 10 MHz).
    nspline : int, optional
        Number of B-spline basis functions per oscillator. Cannot be combined
        with spline_knot_spacing.
    spline_knot_spacing : float, optional
        Spacing between B-spline knots [ns]. Derives nspline from final_time
        and this spacing. Cannot be combined with nspline.
    spline_order : int, optional
        B-spline order: 2 (quadratic, default) or 0 (piecewise constant).
        Affects the knot spacing formula and carrier frequency defaults.
    control_zero_boundary_condition : bool, optional
        Force control pulses to start and end at zero.
    initial_condition : InitialConditionSettings, optional
        Direct struct specification (advanced). For convenience, prefer
        initial_levels or initial_state.
    initial_levels : list of int, optional
        Product-state initial condition, e.g. [0, 0, 1] for |001>.
    initial_state : list of complex, optional
        Arbitrary superposition initial state. Written to file automatically.
    output_directory : str
        Output directory for generated files (initial state, etc.).
        Default: "./run_dir".
    verbose : bool
        Print setup information. Default: True.

    Returns
    -------
    Setup
        Setup with physics parameters configured (no runtype set).

    Examples
    --------
    >>> # Default: BASIS (all basis states)
    >>> setup = setup_quandary(nessential=[3], transition_frequency=[4.1], final_time=100)
    >>> # Product state |001>:
    >>> setup = setup_quandary(nessential=[2, 2, 2], transition_frequency=[4.1, 4.5, 4.9],
    ...                        final_time=100, initial_levels=[0, 0, 1])
    >>> # Superposition (|0> + |1>)/sqrt(2):
    >>> setup = setup_quandary(nessential=[2], transition_frequency=[4.1],
    ...                        final_time=100, initial_state=[1/np.sqrt(2), 1/np.sqrt(2)])
    """
    # Set defaults
    nqubits = len(nessential)
    if selfkerr is None:
        selfkerr = [0.0] * nqubits
    if nguard is None:
        nguard = [0] * nqubits
    if rotation_frequency is None:
        rotation_frequency = transition_frequency[:]
    if crosskerr_coupling is None:
        crosskerr_coupling = []
    if dipole_coupling is None:
        dipole_coupling = []
    # Use 0.01 GHz (= 10 MHz) as default amplitude for timestep estimation only.
    # When user doesn't specify bounds, optimization bounds are left unset (C++ default: 1e12).
    bounds_for_estimation = control_amplitude_bounds if control_amplitude_bounds is not None \
        else [0.01] * nqubits

    # Resolve output directory against QUANDARY_BASE_DATADIR if set
    output_directory = resolve_output_dir(output_directory)

    # Handle initial condition (only one of initial_condition, initial_levels, initial_state)
    num_init_specs = sum([initial_condition is not None, initial_levels is not None, initial_state is not None])
    if num_init_specs > 1:
        raise ValueError("Can only specify one of: initial_condition, initial_levels, initial_state")

    if initial_levels is not None:
        # Product state like |001>
        if len(initial_levels) != nqubits:
            raise ValueError(f"initial_levels must have length {nqubits}, got {len(initial_levels)}")
        initial_condition = InitialConditionSettings()
        initial_condition.condition_type = InitialConditionType.PRODUCT_STATE
        initial_condition.levels = initial_levels
    elif initial_state is not None:
        # Arbitrary superposition, write to file
        dim_ess = int(np.prod(nessential))
        initial_state_array = np.array(initial_state, dtype=complex)
        if len(initial_state_array) != dim_ess:
            raise ValueError(f"initial_state must have length {dim_ess} (product of nessential), "
                             f"got {len(initial_state_array)}")
        os.makedirs(output_directory, exist_ok=True)
        init_state_file = os.path.join(output_directory, "initial_state.dat")
        state_vec = np.concatenate((initial_state_array.real, initial_state_array.imag))
        np.savetxt(init_state_file, state_vec, fmt='%20.13e')
        initial_condition = InitialConditionSettings()
        initial_condition.condition_type = InitialConditionType.FROMFILE
        initial_condition.filename = init_state_file
    elif initial_condition is None:
        # Default: BASIS
        initial_condition = InitialConditionSettings()
        initial_condition.condition_type = InitialConditionType.BASIS

    # Build Hamiltonians
    nlevels = [nessential[i] + nguard[i] for i in range(nqubits)]
    Hsys, Hc_re, Hc_im = hamiltonians(
        N=nlevels,
        transition_frequency=transition_frequency,
        selfkerr=selfkerr,
        crosskerr_coupling=crosskerr_coupling,
        dipole_coupling=dipole_coupling,
        rotation_frequency=rotation_frequency,
        verbose=verbose,
    )

    # Handle time discretization: final_time = ntime * dt
    # User should specify at most 2 of the 3 parameters, unless they are consistent.
    if ntime is None and dt is None:
        # Neither specified - auto-compute ntime from Hamiltonian
        ntime = estimate_timesteps(
            final_time=final_time,
            Hsys=Hsys,
            Hc_re=Hc_re,
            Hc_im=Hc_im,
            control_amplitude_bounds=bounds_for_estimation,
            Pmin=Pmin,
        )
        dt = final_time / ntime
    elif ntime is None and dt is not None:
        # Only dt specified - compute ntime
        ntime = int(np.ceil(final_time / dt))
    elif ntime is not None and dt is None:
        # Only ntime specified - compute dt
        dt = final_time / ntime
    else:
        # Both specified - check consistency
        assert ntime is not None and dt is not None
        computed_final_time = ntime * dt
        if abs(computed_final_time - final_time) > 1e-10:
            raise ValueError(
                f"Inconsistent time parameters: final_time={final_time}, "
                f"ntime={ntime}, dt={dt}. Must satisfy final_time = ntime * dt. "
                f"Got ntime * dt = {computed_final_time}."
            )

    # Compute carrier frequencies (skip if user provided them)
    if carrier_frequency is None:
        carrier_frequency, _ = get_resonances(
            nessential=nessential,
            nguard=nguard,
            Hsys=Hsys,
            Hc_re=Hc_re,
            Hc_im=Hc_im,
            rotation_frequency=rotation_frequency,
            verbose=verbose,
        )

    # Validate spline_order
    if spline_order is not None and spline_order not in (0, 2):
        raise ValueError(f"spline_order must be 0 or 2, got {spline_order}")

    # Compute nspline from knot spacing if given
    if nspline is not None and spline_knot_spacing is not None:
        raise ValueError("Cannot specify both nspline and spline_knot_spacing")
    if spline_knot_spacing is not None:
        order = spline_order if spline_order is not None else 2
        if order == 0:
            nspline = int(np.max([np.rint(final_time / spline_knot_spacing + 1), 2]))
        else:
            enforce_bc = control_zero_boundary_condition if control_zero_boundary_condition is not None else True
            minspline = 5 if enforce_bc else 3
            nspline = int(np.max([np.ceil(final_time / spline_knot_spacing + 2), minspline]))

    if verbose:
        logger.info("Configuration computed:")
        logger.info(f"  Total time: {ntime * dt:.6f} ns")
        logger.info(f"  Time steps: {ntime}")
        logger.info(f"  dt: {dt:.6f} ns")
        logger.info(f"  Carrier frequencies: {carrier_frequency}")
        if nspline is not None:
            logger.info(f"  B-spline basis functions: {nspline}")

    # Create Setup with common fields
    setup = Setup()
    setup.nlevels = nlevels
    setup.nessential = nessential
    setup.ntime = ntime
    setup.dt = dt
    setup.transition_frequency = transition_frequency
    setup.rotation_frequency = rotation_frequency
    setup.selfkerr = selfkerr
    if len(crosskerr_coupling) > 0:
        setup.crosskerr_coupling = crosskerr_coupling
    if len(dipole_coupling) > 0:
        setup.dipole_coupling = dipole_coupling
    if control_amplitude_bounds is not None:
        setup.control_amplitude_bounds = control_amplitude_bounds

    setup.carrier_frequencies = carrier_frequency
    setup.initial_condition = initial_condition
    setup.output_directory = output_directory
    if control_zero_boundary_condition is not None:
        setup.control_zero_boundary_condition = control_zero_boundary_condition

    # Set control parameterizations if nspline or spline_order was specified
    if nspline is not None or spline_order is not None:
        order = spline_order if spline_order is not None else 2
        control_type = ControlType.BSPLINE0 if order == 0 else ControlType.BSPLINE
        control_params = []
        for _ in range(nqubits):
            param = ControlParameterizationSettings()
            param.control_type = control_type
            if nspline is not None:
                param.nspline = nspline
            control_params.append(param)
        setup.control_parameterizations = control_params
        # Order 0 uses zero carrier frequencies by default
        if order == 0:
            setup.carrier_frequencies = [[0.0] for _ in range(nqubits)]

    # Set default output observables
    setup.output_observables = [OutputType.POPULATION, OutputType.EXPECTED_ENERGY, OutputType.FULLSTATE]

    return setup


def _setup_optimization(
    setup: Setup,
    targetgate=None,
    targetstate=None,
    target_levels: Optional[List[int]] = None,
    gate_rot_freq: Optional[List[float]] = None,
    pcof=None,
    randomize_initial_control: bool = False,
    control_initialization_amplitude: Optional[float] = None,
) -> Setup:
    """Return a copy of setup configured for optimization."""
    setup = setup.copy()
    setup.output_directory = resolve_output_dir(setup.output_directory)
    setup.runtype = RunType.OPTIMIZATION

    # Validate: only one target type
    num_targets = sum([targetgate is not None, targetstate is not None, target_levels is not None])
    if num_targets > 1:
        raise ValueError("Can only specify one of: targetgate, targetstate, target_levels")

    output_dir = _get_output_dir(setup)
    os.makedirs(output_dir, exist_ok=True)

    if targetgate is not None:
        # Gate optimization: write gate to file
        gate_array = np.array(targetgate, dtype=complex)
        gate_file = os.path.join(output_dir, "targetgate.dat")
        gate_vec = np.concatenate((
            gate_array.real.ravel(order='F'),
            gate_array.imag.ravel(order='F')
        ))
        np.savetxt(gate_file, gate_vec, fmt='%20.13e')

        target = OptimTargetSettings()
        target.target_type = TargetType.GATE
        target.gate_type = GateType.FILE
        target.filename = gate_file
        if gate_rot_freq is not None:
            target.gate_rot_freq = gate_rot_freq
        setup.optim_target = target

    elif targetstate is not None:
        # State-to-state: write target state to file
        target_array = np.array(targetstate, dtype=complex)
        state_file = os.path.join(output_dir, "targetstate.dat")
        state_vec = np.concatenate((target_array.real.ravel(order='F'), target_array.imag.ravel(order='F')))
        np.savetxt(state_file, state_vec, fmt='%20.13e')

        target = OptimTargetSettings()
        target.target_type = TargetType.STATE
        target.filename = state_file
        setup.optim_target = target

    elif target_levels is not None:
        # State-to-state: product state target like |001>
        target = OptimTargetSettings()
        target.target_type = TargetType.STATE
        target.levels = target_levels
        setup.optim_target = target

    # Set up control initialization (only if user explicitly provides parameters)
    if pcof is not None and len(pcof) > 0:
        # Warm-start from provided coefficients
        pcof_file = os.path.join(output_dir, "pcof_init.dat")
        np.savetxt(pcof_file, pcof, fmt='%20.13e')

        control_inits = []
        for _ in range(len(setup.nessential)):
            init = ControlInitializationSettings()
            init.init_type = ControlInitializationType.FILE
            init.filename = pcof_file
            control_inits.append(init)
        setup.control_initializations = control_inits
    elif control_initialization_amplitude is not None:
        # Fresh optimization with explicit amplitude (use randomize_initial_control to choose type)
        control_inits = []
        init_type = (
            ControlInitializationType.RANDOM if randomize_initial_control
            else ControlInitializationType.CONSTANT
        )

        for _ in range(len(setup.nessential)):
            init = ControlInitializationSettings()
            init.init_type = init_type
            init.amplitude = control_initialization_amplitude
            control_inits.append(init)

        setup.control_initializations = control_inits

    return setup


def _setup_simulation(
    setup: Setup,
    pcof=None,
    pt0=None,
    qt0=None,
) -> Setup:
    """Return a copy of setup configured for simulation."""
    setup = setup.copy()
    setup.output_directory = resolve_output_dir(setup.output_directory)
    setup.runtype = RunType.SIMULATION

    if pt0 is not None and qt0 is not None:
        if pcof is not None:
            raise ValueError("Cannot specify both pcof and pt0/qt0")

        ntime = setup.ntime
        dt = setup.dt
        nsplines = max(2, ntime + 1)
        spline_knot_spacing = dt

        # Set control parameterization to piecewise constant (BSPLINE0)
        control_params = []
        for _ in range(len(setup.nessential)):
            param = ControlParameterizationSettings()
            param.control_type = ControlType.BSPLINE0
            param.nspline = nsplines
            control_params.append(param)
        setup.control_parameterizations = control_params

        # Zero out carrier frequencies (pulses already include carrier)
        setup.carrier_frequencies = [[0.0] for _ in range(len(setup.nessential))]

        # Downsample time-domain pulses to B-spline coefficients
        pcof = downsample_pulses(
            pt0=pt0, qt0=qt0,
            nsplines=nsplines,
            spline_knot_spacing=spline_knot_spacing,
            ntime=ntime, dt=dt,
            nessential=setup.nessential,
        )

    if pcof is not None and len(pcof) > 0:
        output_dir = _get_output_dir(setup)
        os.makedirs(output_dir, exist_ok=True)
        pcof_file = os.path.join(output_dir, "pcof_init.dat")
        np.savetxt(pcof_file, pcof, fmt='%20.13e')

        init = ControlInitializationSettings()
        init.init_type = ControlInitializationType.FILE
        init.filename = pcof_file
        setup.control_initializations = [init]

    return setup


def _setup_eval_controls(
    setup: Setup,
    pcof,
    points_per_ns: float = 1.0,
) -> Setup:
    """Return a copy of setup configured for control evaluation at a given sample rate."""
    setup = setup.copy()
    setup.output_directory = resolve_output_dir(setup.output_directory)
    # Recalculate time grid: keep total time, change resolution
    total_time = setup.ntime * setup.dt
    setup.ntime = int(np.floor(total_time * points_per_ns))
    setup.dt = total_time / setup.ntime

    setup.runtype = RunType.EVALCONTROLS

    # Write pcof to file and set control initializations
    if pcof is not None and len(pcof) > 0:
        output_dir = _get_output_dir(setup)
        os.makedirs(output_dir, exist_ok=True)
        pcof_file = os.path.join(output_dir, "pcof_init.dat")
        np.savetxt(pcof_file, pcof, fmt='%20.13e')

        init = ControlInitializationSettings()
        init.init_type = ControlInitializationType.FILE
        init.filename = pcof_file
        setup.control_initializations = [init]

    return setup
