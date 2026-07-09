"""Helper functions to create and configure ConfigInput objects for simulation and optimization."""

import logging
import os
from collections.abc import Sequence
from typing import Optional
import numpy as np
from .._quandary_impl import ControlType, ControlInitializationType, InitialConditionType, TargetType, GateType, DecoherenceType, OutputType, inputFromFile
from .types import Config, InitialConditionSettings, OptimTargetSettings, ControlParameterizationSettings, ControlInitializationSettings
from .physics import hamiltonians, get_resonances, fit_bspline0, fit_bspline2nd, estimate_timestep_size

logger = logging.getLogger(__name__)

_DEFAULT_OUTPUT_DIR = "./data_out"

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

def load_config(
    filename: str,
    quiet: bool = False,
) -> Config:
    """Set up a Configuration by parsing an existing TOML configuration file.

    Parameters
    ----------
    filename : str
        Path to the TOML configuration file.
    quiet : bool, optional
        Whether to suppress logging output during parsing. Default: False.

    Returns
    -------
    Config
        Config object populated from the TOML file (without validation/finalization).
    """
    # Make sure file exists:
    if not os.path.isfile(filename):
        raise ValueError(f"Configuration file not found: {filename}")

    # inputFromFile returns the nanobind C++ base type.
    return Config(inputFromFile(filename, quiet=quiet))


def validate(config: Config, quiet: bool = False) -> Config:
    """Validate and return a new finalized Config instance.

    This helper mirrors ``Config.validate()`` while providing a module-level
    function for ergonomic imports like ``from quandary.new import validate``.
    """
    return config.validate(quiet=quiet)


def create_config(
    nessential: Sequence[int],
    total_time: float,
    dt: Optional[float] = None,
    transition_frequency: Optional[Sequence[float]] = None,
    selfkerr: Optional[Sequence[float]] = None,
    nguard: Optional[Sequence[int]] = None,
    rotation_frequency: Optional[Sequence[float]] = None,
    crosskerr_coupling: Optional[Sequence[float]] = None,
    dipole_coupling: Optional[Sequence[float]] = None,
    decay_time: Optional[Sequence[float]] = None,
    dephase_time: Optional[Sequence[float]] = None,
    carrier_frequency: Optional[Sequence[Sequence[float]]] = None,
    target: Optional[Sequence[complex]] = None,
    gate_rot_freq: Optional[Sequence[float]] = None,
    Pmin: int = 150,
    control_amplitude_bound: Optional[Sequence[float]] = None,
    nspline: Optional[int] = None,
    spline_knot_spacing: Optional[float] = None,
    spline_order: Optional[int] = None,
    control_zero_boundary_condition: Optional[bool] = None,
    control_randomize: bool = True,
    control_amplitude: Optional[float] = None,
    hamiltonian_Hsys: Optional[np.ndarray] = None,
    hamiltonian_Hc: Optional[Sequence[np.ndarray]] = None,
    initial_condition: Optional[Sequence[complex]] = None,
    cw_amp_thres: Optional[float] = 1e-7,
    output_directory: str = _DEFAULT_OUTPUT_DIR,
) -> Config:
    """Create a Config with physics parameters configured.

    Automatically computes Hamiltonians, timesteps, and carrier frequencies.

    **Time discretization:** total_time is required, timestep size (dt) is optional, 
    and will be computed from Hamiltonian eigenvalues if omitted. 

    **Spline configuration:** specify at most one of (nspline, spline_knot_spacing).
    If neither is given the C++ default (10 splines) is used.

    Parameters
    ----------
    nessential : sequence of int
        Number of essential energy levels per qubit.
    total_time : float
        Pulse duration [ns].
    dt : float, optional
        Timestep size [ns]. Auto-computed if omitted.
    transition_frequency : sequence of float, optional
        01-transition frequencies [GHz] per qubit. Default: zeros (suitable
        for custom Hamiltonians).
    selfkerr : sequence of float, optional
        Anharmonicities [GHz] per qubit. Default: zeros.
    nguard : sequence of int, optional
        Number of guard levels per qubit. Default: zeros.
    rotation_frequency : sequence of float, optional
        Rotating frame frequencies [GHz] per qubit.
        Default: transition_frequency.
    crosskerr_coupling : sequence of float, optional
        ZZ coupling strengths [GHz]. Format: [g01, g02, ..., g12, g13, ...].
        Default: no coupling.
    dipole_coupling : sequence of float, optional
        Dipole-dipole coupling strengths [GHz]. Format: [J01, J02, ..., J12, ...].
        Default: no coupling.
    decay_time : sequence of float, optional
        T1 relaxation times [ns] per qubit. 
        Default: no decay.
    dephase_time : sequence of float, optional
        T2 dephasing times [ns] per qubit. 
        Default: no dephasing.
    carrier_frequency : sequence of sequence of float, optional
        Carrier frequencies [GHz] per oscillator. When provided, the
        automatic eigenvalue-based computation (get_resonances) is skipped.
        Useful when the Hamiltonian has degenerate eigenvalues that cause
        the automatic computation to fail.
    target : array-like, optional 
        Target for optimization and fidelity computation. Either a 2D (unitary gate) or 1D (state vector)
    gate_rot_freq : sequence of float, optional
        Gate rotation frequencies [GHz] per qubit.
    Pmin : int
        Minimum time steps per period of the fastest oscillation (used for
        auto-computing dt). Default: 150.
    control_amplitude_bound : sequence of float, optional
        Maximum control amplitudes [GHz] per qubit. Sets optimization bounds.
        Default: unbounded.
    nspline : int, optional
        Number of B-spline basis functions per oscillator. Cannot be combined
        with spline_knot_spacing.
    spline_knot_spacing : float, optional
        Spacing between B-spline knots [ns]. Derives nspline from total_time
        and this spacing. Cannot be combined with nspline.
    spline_order : int, optional
        B-spline order: 2 (quadratic, default) or 0 (piecewise constant).
        Affects the knot spacing formula and carrier frequency defaults.
    control_zero_boundary_condition : bool, optional
        Force control pulses to start and end at zero.
    control_randomize : bool
        Initialize controls randomly. Default: True.
    control_amplitude : float, optional
        Initial control amplitude [GHz]. When omitted, uses
        config.control_initializations if set, otherwise defaults from C++ code (zero controls)
    hamiltonian_Hsys : ndarray, optional
        Custom system Hamiltonian matrix (complex, in rad/ns). When provided,
        the standard pulse-driven superconducting-qubit Hamiltonian model is
        bypassed.
    hamiltonian_Hc : sequence of ndarray, optional
        Custom control Hamiltonian matrices (complex), one per oscillator.
        Can be used with or without hamiltonian_Hsys.
    initial_condition : sequence of complex or InitialConditionSettings, optional
        Either a state vector (arbitrary superposition), or direct struct specification (advanced). 
        Default: All basis states in the essential dimensions.
    cw_amp_thres : float, optional
        Minimum threshold on growth rate for each carrier. Default: 1e-7.
    output_directory : str
        Output directory for generated files (initial state, etc.).
        Default: "./data_out".

    Returns
    -------
    Config
        Config with physics parameters configured (no runtype set).

    Examples
    --------
    >>> config = create_config(nessential=[3], transition_frequency=[4.1], total_time=100.0)
    >>> config = create_config(nessential=[2, 3], transition_frequency=[4.0, 5.0], selfkerr=[0.2, 0.3], total_time=50.0, control_amplitude_bound=[0.5, 0.3], nspline=20, control_zero_boundary_condition=True)
    """

    # Set defaults
    nqubits = len(nessential)
    if transition_frequency is None:
        transition_frequency = [0.0] * nqubits
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
    bounds_for_estimation = control_amplitude_bound if control_amplitude_bound is not None \
        else [0.01] * nqubits

    # Resolve output directory against QUANDARY_BASE_DATADIR if set
    output_directory = resolve_output_dir(output_directory)

    # Build Hamiltonians
    if hamiltonian_Hsys is not None:
        # Custom Hsys provided — use it directly
        Hsys = np.asarray(hamiltonian_Hsys, dtype=complex)
        if hamiltonian_Hc is not None:
            # Custom Hc provided — split into real/imag for timestep/carrier estimation
            Hc_re = [np.asarray(h, dtype=complex).real for h in hamiltonian_Hc]
            Hc_im = [np.asarray(h, dtype=complex).imag for h in hamiltonian_Hc]
        else:
            # No custom Hc — build standard a+aT operators for estimation
            _, Hc_re, Hc_im = hamiltonians(
                nessential=nessential,
                nguard=nguard,
                transition_frequency=transition_frequency,
                selfkerr=selfkerr,
                crosskerr_coupling=crosskerr_coupling,
                dipole_coupling=dipole_coupling,
                rotation_frequency=rotation_frequency,
            )
    else:
        # Standard model
        Hsys, Hc_re, Hc_im = hamiltonians(
            nessential=nessential,
            nguard=nguard,
            transition_frequency=transition_frequency,
            selfkerr=selfkerr,
            crosskerr_coupling=crosskerr_coupling,
            dipole_coupling=dipole_coupling,
            rotation_frequency=rotation_frequency,
        )
        if hamiltonian_Hc is not None:
            # Custom Hc with standard Hsys — still need Hc_re/Hc_im for estimation
            Hc_re = [np.asarray(h, dtype=complex).real for h in hamiltonian_Hc]
            Hc_im = [np.asarray(h, dtype=complex).imag for h in hamiltonian_Hc]

    # Handle time discretization
    if dt is None:
        # Auto-compute dt from Hamiltonian eigenvalues
        dt_est = estimate_timestep_size(
            Hsys=Hsys,
            Hc_re=Hc_re,
            Hc_im=Hc_im,
            control_amplitude_bound=bounds_for_estimation,
            Pmin=Pmin,
        )
        # Make sure that total_time is an integer multiple of dt to avoid issues with the last timestep.
        k = np.ceil(total_time / dt_est)
        dt = total_time / k

    # Compute carrier frequencies (skip if user provided them)
    if carrier_frequency is None:
        carrier_frequency, _ = get_resonances(
            nessential=nessential,
            nguard=nguard,
            Hsys=Hsys,
            Hc_re=Hc_re,
            Hc_im=Hc_im,
            rotation_frequency=rotation_frequency,
            cw_amp_thres=cw_amp_thres,
        )
    
    # Validate spline_order
    if spline_order is not None and spline_order not in (0, 2):
        raise ValueError(f"spline_order must be 0 or 2, got {spline_order}")

    # Compute nspline from knot spacing if given
    if spline_knot_spacing is not None:
        order = spline_order if spline_order is not None else 2
        if order == 0:
            computed_nspline = int(np.max([np.rint(total_time/ spline_knot_spacing + 1), 2]))
        else:
            enforce_bc = control_zero_boundary_condition if control_zero_boundary_condition is not None else True
            minspline = 5 if enforce_bc else 3
            computed_nspline = int(np.max([np.ceil(total_time/ spline_knot_spacing + 2), minspline]))
        if nspline is not None and nspline != computed_nspline:
            raise ValueError(
                f"Inconsistent spline parameters: nspline={nspline} but "
                f"spline_knot_spacing={spline_knot_spacing} implies nspline={computed_nspline}"
            )
        nspline = computed_nspline

    logger.info("Configuration computed:")
    logger.info(f"  Total time: {total_time:.6f} ns")
    logger.info(f"  dt: {dt:.6f} ns")
    logger.info(f"  Carrier frequencies: {carrier_frequency}")
    if nspline is not None:
        logger.info(f"  B-spline basis functions: {nspline}")

    # Create ConfigInput with common fields
    config = Config()
    config.nlevels = [nessential[i] + nguard[i] for i in range(nqubits)]
    config.nessential = nessential
    config.total_time = total_time
    config.dt = dt
    config.transition_frequency = transition_frequency
    config.rotation_frequency = rotation_frequency
    config.selfkerr = selfkerr
    if len(crosskerr_coupling) > 0:
        config.crosskerr_coupling = crosskerr_coupling
    if len(dipole_coupling) > 0:
        config.dipole_coupling = dipole_coupling
    if control_amplitude_bound is not None:
        config.control_amplitude_bound = control_amplitude_bound
    config.carrier_frequencies = carrier_frequency
    config.output_directory = output_directory
    if control_zero_boundary_condition is not None:
        config.control_zero_boundary_condition = control_zero_boundary_condition
    
    # Set decoherence if provided
    set_decoherence(config, decay_time=decay_time, dephase_time=dephase_time)

    # Write target and initial conditions to file, if provided, and store in Config.
    set_target(config, target, gate_rot_freq=gate_rot_freq)
    set_initial_condition(config, initial_condition=initial_condition)

    # Write custom Hamiltonian files and set paths on Config
    if hamiltonian_Hsys is not None or hamiltonian_Hc is not None:
        hsys_path, hc_path = _write_hamiltonian_files(
            output_directory,
            Hsys=hamiltonian_Hsys,
            Hc=hamiltonian_Hc,
        )
        if hsys_path is not None:
            config.hamiltonian_file_Hsys = hsys_path
        if hc_path is not None:
            config.hamiltonian_file_Hc = hc_path

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
        config.control_parameterizations = control_params
        # Order 0 uses zero carrier frequencies by default
        if order == 0:
            config.carrier_frequencies = [[0.0] for _ in range(nqubits)]

    # Set initial control pulse, if provided
    set_controls(config, control_randomize=control_randomize, control_amplitude=control_amplitude)

    # Set default output observables
    config.output_observables = [OutputType.POPULATION, OutputType.EXPECTED_ENERGY, OutputType.FULLSTATE]

    return config


def set_target(
    config: Config,
    target: np.ndarray,
    gate_rot_freq: Optional[Sequence[float]] = None,
) -> None:
    """Set the optimization target on a Config (in-place).

    The type is inferred from the dimensionality of ``target``: a 2-D array
    is treated as a gate, a 1-D array as a state.  Writes the data to a file
    in the output directory and sets ``config.optim_target``.

    Parameters
    ----------
    config : Config
        Config to modify in-place.
    target : array-like
        Target for optimization and fidelity computation. Either a 2D (unitary gate) or 1D (state vector).
    gate_rot_freq : sequence of float, optional
        Rotation frequencies for the gate. Default: None.

    Returns
    -------
    None
        Modifies config in-place to set the optimization target.
    """
    if target is None:
        return

    target_array = np.asarray(target, dtype=complex)
    dim_ess = int(np.prod(config.nessential))

    output_dir = _get_output_dir(config)
    os.makedirs(output_dir, exist_ok=True)

    if target_array.ndim == 2:
        # Gate optimization: verify dimensions and write gate to file
        if target_array.shape != (dim_ess, dim_ess):
            raise ValueError(f"Target gate must have shape ({dim_ess}, {dim_ess}), got {target_array.shape}")
        gate_file = os.path.join(output_dir, "target_gate.dat")
        gate_vec = np.concatenate((
            target_array.real.ravel(order='F'),
            target_array.imag.ravel(order='F')
        ))
        np.savetxt(gate_file, gate_vec, fmt='%20.13e')

        optim_target = OptimTargetSettings()
        optim_target.target_type = TargetType.GATE
        optim_target.gate_type = GateType.FILE
        optim_target.filename = gate_file
        if gate_rot_freq is not None:
            optim_target.gate_rot_freq = gate_rot_freq
        config.optim_target = optim_target

    else:
        # State-to-state: verify dimensions and write target state to file
        # dimension must match to nessential product
        if target_array.shape != (dim_ess,):
            raise ValueError(f"Target state must have shape ({dim_ess},), got {target_array.shape}")
        state_file = os.path.join(output_dir, "target_state.dat")
        state_vec = np.concatenate((target_array.real.ravel(order='F'), target_array.imag.ravel(order='F')))
        np.savetxt(state_file, state_vec, fmt='%20.13e')

        optim_target = OptimTargetSettings()
        optim_target.target_type = TargetType.STATE
        optim_target.filename = state_file
        config.optim_target = optim_target

def set_initial_condition(
    config: Config,
    initial_condition = None,
) -> None:
    """Set the initial condition on a Config (in-place).

    The initial condition can be either as a state vector (arbitrary superposition) or an InitialConditionSettings struct for advanced use cases. If it's a state vector, it is written to a file in the output directory and the struct is configured to load from that file. If it's already an InitialConditionSettings struct, it is used directly. If None, the C++ default is used. 

    Parameters
    ----------
    config : Config
        Config to modify in-place.
    initial_condition : array-like or InitialConditionSettings, optional
        Either a state vector (arbitrary superposition) or an InitialConditionSettings struct. If None, the C++ default is used (all basis states in the essential dimensions).

    Returns
    -------
    None
        Modifies config in-place to set the initial condition.
    """
    if initial_condition is None:
        return  

    if isinstance(initial_condition, InitialConditionSettings):
        # Already an InitialConditionSettings struct, use it directly
        config.initial_condition = initial_condition

    else:
        # Initial condition is a state vector (arbitrary superposition). Check dimensions and write to file
        dim_ess = int(np.prod(config.nessential))
        initial_state_array = np.array(initial_condition, dtype=complex)
        if len(initial_state_array) != dim_ess:
            raise ValueError(f"initial_condition state must have length {dim_ess} (product of nessential), "
                             f"got {len(initial_state_array)}")
        output_dir = _get_output_dir(config)
        os.makedirs(output_dir, exist_ok=True)

        init_state_file = os.path.join(output_dir, "initial_state.dat")
        state_vec = np.concatenate((initial_state_array.real, initial_state_array.imag))
        np.savetxt(init_state_file, state_vec, fmt='%20.13e')
        initial_condition = InitialConditionSettings()
        initial_condition.condition_type = InitialConditionType.FROMFILE
        initial_condition.filename = init_state_file
        config.initial_condition = initial_condition


def set_controls(
    config: Config,
    spline_coefficients=None,
    p_samples=None,
    q_samples=None,
    control_amplitude: Optional[float] = None,
    control_randomize: bool = True,
) -> None:
    """Set the control parameterization and initialization on a Config (in-place).

    The control pulse is defined either by providing the B-spline coefficients (spline_coefficients), or by lists of control pulses at each time point (p_samples, q_samples), or by an explicit amplitude (control_amplitude) for uniform or random initialization.    

    Priority: p_samples/q_samples > spline_coefficients > explicit amplitude > existing cfg.control_initializations > default from C++ code.

    Parameters
    ----------
    config : Config
        Config to modify in-place.
    spline_coefficients : array-like, optional
        B-spline coefficients for the control pulses. If provided, they are written to a file and the control initialization is set to load from that file.
    p_samples : array-like, optional
        Control pulse samples for the p quadrature for each qubit and each timestep. If provided, they are fitted to B-splines and the coefficients are written to a file. Requires q_samples to be provided as well. Overrides spline_coefficients if both are provided.
    q_samples : array-like, optional
        Control pulse samples for the q quadrature for each qubit and each timestep. If provided, they are fitted to B-splines and the coefficients are written to a file. Requires p_samples to be provided as well. Overrides spline_coefficients if both are provided.
    control_amplitude : float, optional
        Amplitude for uniform or random initialization. Overrides existing control_initializations if provided.
    control_randomize : bool, optional
        Whether to randomize the control pulses if an explicit amplitude is provided. Default is True.

    Returns
    -------
    None
        Modifies cfg in-place to set the control parameterization and initialization.
    """

    # If p_samples/q_samples are provide, fit pulses to bspline coefficients. 
    if p_samples is not None or q_samples is not None:
        if spline_coefficients is not None:
            raise ValueError("Cannot specify both spline_coefficients and p_samples/q_samples")
        if p_samples is None:
            p_samples = np.zeros_like(q_samples)
        if q_samples is None:
            q_samples = np.zeros_like(p_samples)
        
        # If control parameterization was not specified, choose Bspline 2nd order for better fitting quality.
        existing_control_params = config.control_parameterizations or []
        if len(existing_control_params) == 0:
            control_type = ControlType.BSPLINE
        else:
            control_type = existing_control_params[0].control_type

        # Fit control parameters to either Bspline 0-th order or Bspline 2nd order. 
        if control_type == ControlType.BSPLINE0:
            nsteps = p_samples.shape[1]
            nsplines = [max(2, nsteps + 1) for _ in range(len(config.nessential))]
            spline_coefficients = fit_bspline0(
                p_samples=p_samples, q_samples=q_samples,
                nsplines=nsplines[0],
                spline_knot_spacing=config.dt,
                dt=config.dt,
                nessential=config.nessential,
            )
            # Zero out carrier frequencies (pulses already include carrier)
            config.carrier_frequencies = [[0.0] for _ in range(len(config.nessential))]
        elif control_type == ControlType.BSPLINE:
            n_osc = len(config.nessential)
            if len(existing_control_params) == 0:
                # Default to spline_knot_spacing of 3ns.
                spline_knot_spacing = 3.0
                computed_nspline = int(np.max([np.ceil(config.total_time / spline_knot_spacing + 2), 5]))
                nsplines = [computed_nspline for _ in range(n_osc)]
            else:
                nsplines = [param.nspline for param in existing_control_params]
            spline_coefficients = fit_bspline2nd(0.0, 
                                  config.total_time, 
                                  p_samples, q_samples, 
                                  nsplines, 
                                  carrier_frequencies=config.carrier_frequencies,
                                  inputs_in_mhz=True )

        # Set control parameterization
        control_params = []
        for i in range(len(config.nessential)):
            param = ControlParameterizationSettings()
            param.control_type = control_type
            param.nspline = nsplines[i]
            control_params.append(param)
        config.control_parameterizations = control_params

    # If spline_coefficients is provided (either function input or set above by fitting splines to p_samples/q_samples), write coefficients to file and set control initialization to load from that file.
    if spline_coefficients is not None and len(spline_coefficients) > 0:
        output_dir = _get_output_dir(config)
        os.makedirs(output_dir, exist_ok=True)
        spline_coefficients_file = os.path.join(output_dir, "spline_coefficients_init.dat")
        np.savetxt(spline_coefficients_file, spline_coefficients, fmt='%20.13e')

        control_inits = []
        for _ in range(len(config.nessential)):
            init = ControlInitializationSettings()
            init.init_type = ControlInitializationType.FILE
            init.filename = spline_coefficients_file
            control_inits.append(init)
        config.control_initializations = control_inits

    elif control_amplitude is not None:
        # Explicit amplitude — create uniform per-oscillator inits
        control_inits = []
        init_type = (
            ControlInitializationType.RANDOM if control_randomize
            else ControlInitializationType.CONSTANT
        )
        for _ in range(len(config.nessential)):
            init = ControlInitializationSettings()
            init.init_type = init_type
            init.amplitude = control_amplitude
            control_inits.append(init)

        config.control_initializations = control_inits

def set_decoherence(config: Config, decay_time=None, dephase_time=None) -> None:
    """Set the decoherence parameters on a Config (in-place).

    Configures the decay and dephasing times, and sets the decoherence type accordingly. If decay or dephasing times are provided, Lindblad's master equation is solved. If neither is provided, decoherence_type is left unset (C++ default: no decoherence), solving Schroedinger's equation.

    Parameters
    ----------
    config : Config
        Config to modify in-place.
    decay_time : sequence of float, optional
        T1 relaxation times [ns] per qubit. If provided, sets decoherence_type to DECAY or BOTH.
    dephase_time : sequence of float, optional
        T2 dephasing times [ns] per qubit. If provided, sets decoherence_type to DEPHASE or BOTH.

    Returns
    -------
    None
        Modifies config in-place to set decoherence parameters.
    """

    # Set decoherence type if provided
    decoherence_type = None 
    nqubits = len(config.nessential)
    if decay_time is not None or dephase_time is not None:
        if decay_time is not None and len(decay_time) != nqubits:
            raise ValueError(f"decay_time must have length {nqubits}, got {len(decay_time)}")
        if dephase_time is not None and len(dephase_time) != nqubits:
            raise ValueError(f"dephase_time must have length {nqubits}, got {len(dephase_time)}")
        if decay_time is not None and dephase_time is not None:
            decoherence_type = DecoherenceType.BOTH
        elif decay_time is not None:
            decoherence_type = DecoherenceType.DECAY
            dephase_time = None
        elif dephase_time is not None:
            decoherence_type = DecoherenceType.DEPHASE
            decay_time = None

    if decoherence_type is not None:
        config.decoherence_type = decoherence_type
    if decay_time is not None:
        config.decay_time = decay_time
    if dephase_time is not None:
        config.dephase_time = dephase_time


# ---------------------------------
# Private config I/O helpers
# ---------------------------------

def _get_output_dir(config):
    """Get output directory from config, falling back to default."""
    return config.output_directory or _DEFAULT_OUTPUT_DIR


def _write_hamiltonian_files(output_directory, Hsys=None, Hc=None):
    """Write custom Hamiltonian matrices to sparse COO format files.

    Parameters
    ----------
    output_directory : str
        Directory to write the files into (created if needed).
    Hsys : ndarray, optional
        System Hamiltonian matrix (complex). Written as sparse COO with
        columns: row col real imag.
    Hc : sequence of ndarray, optional
        Control Hamiltonian matrices (complex), one per oscillator. Written
        as sparse COO with columns: oscillator row col real imag.

    Returns
    -------
    hsys_path : str or None
        Path to the Hsys file, or None if Hsys was not provided.
    hc_path : str or None
        Path to the Hc file, or None if Hc was not provided.
    """
    os.makedirs(output_directory, exist_ok=True)
    hsys_path = None
    hc_path = None

    if Hsys is not None:
        hsys_path = os.path.join(output_directory, "hamiltonian_Hsys.dat")
        H = np.asarray(Hsys, dtype=complex)
        with open(hsys_path, "w", newline='\n') as f:
            f.write("# row col Hsys_real Hsys_imag\n")
            nz = np.nonzero(H)
            for i, j in zip(*nz):
                v = H[i, j]
                f.write(f"{i} {j} {v.real:.13e} {v.imag:.13e}\n")

    if Hc is not None and len(Hc) > 0:
        hc_path = os.path.join(output_directory, "hamiltonian_Hc.dat")
        with open(hc_path, "w", newline='\n') as f:
            for iosc, Hc_osc in enumerate(Hc):
                Hc_mat = np.asarray(Hc_osc, dtype=complex)
                nz = np.nonzero(Hc_mat)
                f.write("# oscillator row col Hc_real Hc_imag\n")
                for i, j in zip(*nz):
                    v = Hc_mat[i, j]
                    f.write(f"{iosc} {i} {j} {v.real:.13e} {v.imag:.13e}\n")

    paths = [p for p in (hsys_path, hc_path) if p is not None]
    if paths:
        logger.info("Hamiltonian operators written to %s", ", ".join(paths))

    return hsys_path, hc_path
