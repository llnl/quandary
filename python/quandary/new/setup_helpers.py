"""Helper functions to create and configure Setup objects for simulation and optimization."""

import logging
import os
from typing import List, Optional

import numpy as np

from .._quandary_impl import (
    Setup,
    RunType,
    ControlType,
    InitialConditionType,
    InitialConditionSettings,
    OptimTargetSettings,
    TargetType,
    GateType,
    ControlParameterizationSettings,
    ControlInitializationSettings,
    ControlInitializationType,
    OutputType,
)
from .quantum_operators import hamiltonians, get_resonances
from .time_estimation import estimate_timesteps
from .utils import downsample_pulses

logger = logging.getLogger(__name__)


def setup_physics(
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
    Pmin: int = 150,
    amplitude_bound: Optional[List[float]] = None,
    nspline: Optional[int] = None,
    spline_knot_spacing: Optional[float] = None,
    spline_order: Optional[int] = None,
    control_zero_boundary_condition: Optional[bool] = None,
    initialcondition: Optional[InitialConditionSettings] = None,
    verbose: bool = True,
) -> Setup:
    """Create a Setup with physics parameters configured.

    Automatically computes Hamiltonians, timesteps, and carrier frequencies.
    Does NOT set the runtype - use setup_optimization() or setup_simulation()
    to configure for a specific run type.

    Time discretization: You can specify at most 2 of (final_time, ntime, dt).
    If you specify all 3, they must be consistent (final_time = ntime * dt).
    If ntime and dt are not provided, ntime will be auto-computed from Hamiltonian.

    Spline configuration: You can specify at most one of (nspline, spline_knot_spacing).
    If neither is given, the C++ default (10 splines) is used.

    Required Parameters:
    -------------------
    nessential : List[int]
        Number of essential energy levels per qubit
    transition_frequency : List[float]
        01-transition frequencies [GHz] per qubit
    final_time : float
        Pulse duration [ns]

    Optional Parameters:
    -------------------
    ntime : int
        Number of timesteps. If not provided, will be auto-computed from Hamiltonian.
    dt : float
        Timestep size [ns]. If not provided, will be auto-computed.
    selfkerr : List[float]
        Anharmonicities [GHz] per qubit. Default: zeros
    nguard : List[int]
        Number of guard levels per qubit. Default: zeros
    rotation_frequency : List[float]
        Frequency of rotations for computational frame [GHz] per qubit.
        Default: transition_frequency
    crosskerr_coupling : List[float]
        ZZ coupling strength [GHz]. Format: [g01, g02, ..., g12, g13, ...]
    dipole_coupling : List[float]
        Dipole-dipole coupling strength [GHz]. Format: [J01, J02, ..., J12, J13, ...]
    Pmin : int
        Minimum points to resolve shortest period (determines timesteps). Default: 150
    amplitude_bound : List[float]
        Maximum control amplitudes [GHz] per qubit. Used for timestep estimation
        and as optimization bounds. Default: 0.01 GHz (= 10 MHz)
    nspline : int
        Number of B-spline basis functions per oscillator. Cannot be combined
        with spline_knot_spacing.
    spline_knot_spacing : float
        Spacing between B-spline knots [ns]. Computes nspline from
        final_time and spacing. Cannot be combined with nspline.
    spline_order : int
        B-spline order: 2 (quadratic, default) or 0 (piecewise constant).
        Affects the knot spacing formula and carrier frequency defaults.
    control_zero_boundary_condition : bool
        Force control pulses to start and end at zero.
        Affects minimum nspline when using spline_knot_spacing.
    initialcondition : InitialConditionSettings
        Initial state specification. Default: basis states
    verbose : bool
        Print setup information. Default: True

    Returns:
    -------
    Setup
        Setup with physics parameters configured (no runtype set).

    Example:
    -------
    >>> setup = setup_physics(nessential=[3], transition_frequency=[4.1], final_time=100)
    >>> setup_optimization(setup, targetgate=[[0,1],[1,0]])
    >>> results = run(setup)
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
    if amplitude_bound is None:
        amplitude_bound = [0.01] * nqubits
    if initialcondition is None:
        initialcondition = InitialConditionSettings()
        initialcondition.condition_type = InitialConditionType.BASIS

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
            amplitude_bound=amplitude_bound,
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

    # Compute carrier frequencies
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
    setup.control_amplitude_bounds = amplitude_bound

    setup.carrier_frequencies = carrier_frequency
    setup.initial_condition = initialcondition
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


def setup_optimization(
    setup: Setup,
    targetgate=None,
    pcof=None,
    randomize_init_ctrl: bool = False,
    init_amplitude_ghz: Optional[float] = None,
):
    """Configure a Setup for optimization.

    Modifies setup in-place to set runtype to OPTIMIZATION and configure
    the target gate and control initialization.

    Parameters:
    ----------
    setup : Setup
        Setup object (e.g. from setup_physics()).
    targetgate : array-like, optional
        Target unitary gate. Can be a plain Python list like [[0,1],[1,0]]
        or numpy array. Will be converted to complex automatically.
    pcof : array-like, optional
        B-spline control coefficients to use as initial guess (warm-start).
        If provided, writes to file and initializes from file, overriding
        randomize_init_ctrl and init_amplitude_ghz. Useful for continuing
        optimization from a previous result. Typically from results.pcof.
    randomize_init_ctrl : bool
        Initialize controls randomly (True) or at constant value (False).
        Only used if init_amplitude_ghz is provided. Default: False.
    init_amplitude_ghz : float, optional
        Initial control amplitude [GHz] per oscillator. Default: None.

    Notes:
    ------
    Control initialization behavior:
    - If pcof provided: initialize from file (warm-start)
    - If init_amplitude_ghz provided: initialize with specified amplitude
      (random or constant based on randomize_init_ctrl)
    - If neither provided: don't set control_initializations, use C++ defaults

    Example:
    -------
    >>> # Use C++ defaults for control initialization
    >>> setup = setup_physics(nessential=[3], transition_frequency=[4.1], final_time=100)
    >>> setup_optimization(setup, targetgate=[[0,0,1],[0,1,0],[1,0,0]])
    >>> results = run(setup)

    >>> # With explicit initial amplitude (in GHz, matching C++ interface):
    >>> setup_optimization(setup, targetgate=gate, init_amplitude_ghz=0.01)

    >>> # Constant initialization instead of random:
    >>> setup_optimization(setup, targetgate=gate, randomize_init_ctrl=False, init_amplitude_ghz=0.01)

    >>> # Warm-start optimization from previous result:
    >>> setup_optimization(setup, targetgate=gate, pcof=previous_results.pcof)
    """
    setup.runtype = RunType.OPTIMIZATION

    if targetgate is not None:
        output_dir = setup.output_directory if setup.output_directory else "./run_dir"
        os.makedirs(output_dir, exist_ok=True)

        # Convert to numpy array with complex dtype
        gate_array = np.array(targetgate, dtype=complex)

        # Write gate to file (column-major order)
        gate_file = os.path.join(output_dir, "targetgate.dat")
        gate_vec = np.concatenate((
            gate_array.real.ravel(order='F'),
            gate_array.imag.ravel(order='F')
        ))
        np.savetxt(gate_file, gate_vec, fmt='%20.13e')

        # Set up target settings
        target = OptimTargetSettings()
        target.target_type = TargetType.GATE
        target.gate_type = GateType.FILE
        target.filename = gate_file
        setup.optim_target = target

    # Set up control initialization (only if user explicitly provides parameters)
    if pcof is not None and len(pcof) > 0:
        # Warm-start from provided coefficients
        output_dir = setup.output_directory if setup.output_directory else "./run_dir"
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
    elif init_amplitude_ghz is not None:
        # Fresh optimization with explicit amplitude (use randomize_init_ctrl to choose type)
        control_inits = []
        init_type = ControlInitializationType.RANDOM if randomize_init_ctrl else ControlInitializationType.CONSTANT

        for _ in range(len(setup.nessential)):
            init = ControlInitializationSettings()
            init.init_type = init_type
            init.amplitude = init_amplitude_ghz
            control_inits.append(init)

        setup.control_initializations = control_inits
    # else: Don't set control_initializations, let C++ use its defaults


def setup_simulation(
    setup: Setup,
    pcof=None,
    pt0=None,
    qt0=None,
):
    """Configure a Setup for simulation.

    Modifies setup in-place to set runtype to SIMULATION and optionally
    configure control parameters from B-spline coefficients or time-domain
    pulses.

    Parameters:
    ----------
    setup : Setup
        Setup object (e.g. from setup_physics()).
    pcof : array-like, optional
        B-spline control coefficients. If provided, writes to file and
        configures control initialization to load from file.
    pt0 : list of ndarray, optional
        Real part of control pulses [MHz] per oscillator. Must be provided
        together with qt0. Automatically downsampled to B-spline coefficients
        using piecewise constant (order 0) representation.
    qt0 : list of ndarray, optional
        Imaginary part of control pulses [MHz] per oscillator. Must be
        provided together with pt0.

    Example:
    -------
    >>> setup = setup_physics(nessential=[3], transition_frequency=[4.1], final_time=100)
    >>> setup_simulation(setup, pcof=results_opt.pcof)
    >>> results = run(setup)
    """
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
        output_dir = setup.output_directory if setup.output_directory else "./run_dir"
        os.makedirs(output_dir, exist_ok=True)
        pcof_file = os.path.join(output_dir, "pcof_init.dat")
        np.savetxt(pcof_file, pcof, fmt='%20.13e')

        init = ControlInitializationSettings()
        init.init_type = ControlInitializationType.FILE
        init.filename = pcof_file
        setup.control_initializations = [init]


def setup_eval_controls(
    setup: Setup,
    pcof,
    points_per_ns: float = 1.0,
):
    """Configure a Setup for evaluating control pulses at a specific sample rate.

    Modifies setup in-place to set runtype to EVALCONTROLS with a new time grid
    based on points_per_ns. Total pulse duration is preserved.

    Parameters:
    ----------
    setup : Setup
        Setup object (e.g. from setup_physics()). Typically a .copy() of the
        optimization setup.
    pcof : array-like
        B-spline control coefficients. Typically from results.pcof.
    points_per_ns : float
        Sample rate for output [points per ns]. Default: 1.0 (1ns spacing).
        Examples: 1.0 = 1ns, 10.0 = 0.1ns, 0.5 = 2ns spacing.

    Example:
    -------
    >>> setup_eval = setup.copy()
    >>> setup_eval_controls(setup_eval, pcof=results_opt.pcof, points_per_ns=64)
    >>> eval_results = run(setup_eval)
    """
    # Recalculate time grid: keep total time, change resolution
    total_time = setup.ntime * setup.dt
    setup.ntime = int(np.floor(total_time * points_per_ns))
    setup.dt = total_time / setup.ntime

    setup.runtype = RunType.EVALCONTROLS

    # Write pcof to file and set control initializations
    if pcof is not None and len(pcof) > 0:
        output_dir = setup.output_directory if setup.output_directory else "./run_dir"
        os.makedirs(output_dir, exist_ok=True)
        pcof_file = os.path.join(output_dir, "pcof_init.dat")
        np.savetxt(pcof_file, pcof, fmt='%20.13e')

        init = ControlInitializationSettings()
        init.init_type = ControlInitializationType.FILE
        init.filename = pcof_file
        setup.control_initializations = [init]
