"""Factory functions to create QuandaryConfig with automatic setup from physics parameters."""

import logging
import os
from typing import List, Optional

import numpy as np

from .._quandary_impl import (
    QuandaryConfig,
    RunType,
    InitialConditionType,
    InitialConditionSettings,
    OptimTargetSettings,
    TargetType,
    GateType,
    ControlInitializationSettings,
    ControlInitializationType,
    OutputType,
)
from .quantum_operators import hamiltonians, get_resonances
from .time_estimation import estimate_timesteps

logger = logging.getLogger(__name__)


def _setup_physics(
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
    initialcondition: Optional[InitialConditionSettings] = None,
    verbose: bool = True,
) -> QuandaryConfig:
    """Internal: Build base config with common physics parameters.

    Computes Hamiltonians, timesteps, and carrier frequencies, then creates
    a QuandaryConfig with all common fields set. Runtype is not set here.

    Time discretization: You can specify at most 2 of (final_time, ntime, dt).
    If you specify all 3, they must be consistent (final_time ≈ ntime * dt).
    If ntime and dt are not provided, ntime will be auto-computed from Hamiltonian.

    Returns:
        QuandaryConfig with physics parameters configured.
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
        amplitude_bound = [10.0] * nqubits
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
        # Note: actual final time will be ntime * dt (may differ slightly from requested final_time)
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

    if verbose:
        logger.info("Configuration computed:")
        logger.info(f"  Total time: {ntime * dt:.6f} ns")
        logger.info(f"  Time steps: {ntime}")
        logger.info(f"  dt: {dt:.6f} ns")
        logger.info(f"  Carrier frequencies: {carrier_frequency}")

    # Create config with common fields
    config = QuandaryConfig()
    config.nlevels = nlevels
    config.nessential = nessential
    config.ntime = ntime
    config.dt = dt
    config.transition_frequency = transition_frequency
    config.rotation_frequency = rotation_frequency
    config.selfkerr = selfkerr
    if len(crosskerr_coupling) > 0:
        config.crosskerr_coupling = crosskerr_coupling
    if len(dipole_coupling) > 0:
        config.dipole_coupling = dipole_coupling
    config.control_amplitude_bounds = amplitude_bound

    config.carrier_frequencies = carrier_frequency
    config.initial_condition = initialcondition

    # Set default output observables
    config.output_observables = [OutputType.POPULATION, OutputType.EXPECTED_ENERGY, OutputType.FULLSTATE]

    return config


def create_simulation_config(
    nessential: List[int],
    transition_frequency: List[float],
    final_time: float,
    ntime: Optional[int] = None,
    dt: Optional[float] = None,
    pcof = None,
    selfkerr: Optional[List[float]] = None,
    nguard: Optional[List[int]] = None,
    rotation_frequency: Optional[List[float]] = None,
    crosskerr_coupling: Optional[List[float]] = None,
    dipole_coupling: Optional[List[float]] = None,
    Pmin: int = 150,
    amplitude_bound: Optional[List[float]] = None,
    initialcondition: Optional[InitialConditionSettings] = None,
    verbose: bool = True,
) -> QuandaryConfig:
    """
    Create a simulation config with automatic Hamiltonian and timestep computation.

    This factory function automatically:
    - Builds Hamiltonians using hamiltonians()
    - Estimates timesteps using estimate_timesteps()
    - Computes carrier frequencies using get_resonances()

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
        Note: Specify at most 2 of (final_time, ntime, dt). They must satisfy final_time = ntime * dt.
    selfkerr : List[float]
        Anharmonicities [GHz] per qubit. Default: zeros
    nguard : List[int]
        Number of guard levels per qubit. Default: zeros
    rotation_frequency : List[float]
        Frequency of rotations for computational frame [GHz] per qubit. Default: transition_frequency
    crosskerr_coupling : List[float]
        ZZ coupling strength [GHz]. Format: [g01, g02, ..., g12, g13, ...]
    dipole_coupling : List[float]
        Dipole-dipole coupling strength [GHz]. Format: [J01, J02, ..., J12, J13, ...]
    Pmin : int
        Minimum points to resolve shortest period (determines timesteps). Default: 150
    amplitude_bound : List[float]
        Maximum control amplitudes [MHz] per qubit for timestep estimation
    initialcondition : InitialConditionSettings
        Initial state specification. Default: basis states
    verbose : bool
        Print setup information. Default: True

    Returns:
    -------
    QuandaryConfig
        Configured config ready for simulation. User can modify any fields before
        passing to run().

    Example:
    -------
    >>> # Simple simulation with auto-computed timesteps
    >>> config = create_simulation_config(nessential=[2], transition_frequency=[4.1], final_time=50.0)
    >>> # Or specify exact number of timesteps
    >>> config = create_simulation_config(nessential=[2], transition_frequency=[4.1], final_time=50.0, ntime=1000)
    >>> results = run(config)
    """
    config = _setup_physics(
        nessential=nessential,
        transition_frequency=transition_frequency,
        final_time=final_time,
        ntime=ntime,
        dt=dt,
        selfkerr=selfkerr,
        nguard=nguard,
        rotation_frequency=rotation_frequency,
        crosskerr_coupling=crosskerr_coupling,
        dipole_coupling=dipole_coupling,
        Pmin=Pmin,
        amplitude_bound=amplitude_bound,
        initialcondition=initialcondition,
        verbose=verbose,
    )

    config.runtype = RunType.SIMULATION

    # If pcof is provided, write to file and set up control initialization
    if pcof is not None:
        # Write vector to temp file in output directory
        output_dir = config.output_directory if hasattr(config, 'output_directory') and config.output_directory else "./run_dir"
        os.makedirs(output_dir, exist_ok=True)
        pcof_file = os.path.join(output_dir, "pcof_init.dat")
        np.savetxt(pcof_file, pcof, fmt='%20.13e')

        # Set up control initialization to load from file
        control_inits = []
        for iosc in range(len(nessential)):
            init = ControlInitializationSettings()
            init.init_type = ControlInitializationType.FILE
            init.filename = pcof_file
            control_inits.append(init)
        config.control_initializations = control_inits

    return config


def create_optimization_config(
    nessential: List[int],
    transition_frequency: List[float],
    final_time: float,
    ntime: Optional[int] = None,
    dt: Optional[float] = None,
    targetgate = None,
    output_directory: str = "./run_dir",
    selfkerr: Optional[List[float]] = None,
    nguard: Optional[List[int]] = None,
    rotation_frequency: Optional[List[float]] = None,
    crosskerr_coupling: Optional[List[float]] = None,
    dipole_coupling: Optional[List[float]] = None,
    Pmin: int = 150,
    amplitude_bound: Optional[List[float]] = None,
    initialcondition: Optional[InitialConditionSettings] = None,
    verbose: bool = True,
) -> QuandaryConfig:
    """
    Create an optimization config with automatic Hamiltonian and timestep computation.

    Similar to create_simulation_config() but sets runtype to OPTIMIZATION.

    Required Parameters:
    -------------------
    nessential : List[int]
        Number of essential energy levels per qubit
    transition_frequency : List[float]
        01-transition frequencies [GHz] per qubit
    final_time : float
        Pulse duration [ns]
    targetgate : array-like (list, numpy array, etc.)
        Target unitary gate. Can be plain Python list like [[0,1],[1,0]]
        or numpy array. Will be converted to complex automatically.

    Optional Parameters:
    -------------------
    ntime : int
        Number of timesteps. If not provided, will be auto-computed from Hamiltonian.
    dt : float
        Timestep size [ns]. If not provided, will be auto-computed.
        Note: Specify at most 2 of (final_time, ntime, dt). They must satisfy final_time = ntime * dt.
    output_directory : str
        Output directory for results and gate files. Default: "./run_dir"
    selfkerr : List[float]
        Anharmonicities [GHz] per qubit. Default: zeros
    nguard : List[int]
        Number of guard levels per qubit. Default: zeros
    rotation_frequency : List[float]
        Frequency of rotations for computational frame [GHz] per qubit. Default: transition_frequency
    crosskerr_coupling : List[float]
        ZZ coupling strength [GHz]
    dipole_coupling : List[float]
        Dipole-dipole coupling strength [GHz]
    Pmin : int
        Minimum points to resolve shortest period. Default: 150
    amplitude_bound : List[float]
        Maximum control amplitudes [MHz] per qubit
    initialcondition : InitialConditionSettings
        Initial state specification. Default: basis states
    verbose : bool
        Print setup information. Default: True

    Returns:
    -------
    QuandaryConfig
        Configured config ready for optimization.

    Example:
    -------
    >>> from quandary.new import create_optimization_config, run
    >>>
    >>> config = create_optimization_config(
    ...     nessential=[2], transition_frequency=[4.1], final_time=50.0,
    ...     targetgate=[[0, 1], [1, 0]]  # X gate
    ... )
    >>> # Or with explicit timesteps:
    >>> config = create_optimization_config(
    ...     nessential=[2], transition_frequency=[4.1], final_time=50.0, ntime=1000,
    ...     targetgate=[[0, 1], [1, 0]]
    ... )
    >>> results = run(config)
    """
    config = _setup_physics(
        nessential=nessential,
        transition_frequency=transition_frequency,
        final_time=final_time,
        ntime=ntime,
        dt=dt,
        selfkerr=selfkerr,
        nguard=nguard,
        rotation_frequency=rotation_frequency,
        crosskerr_coupling=crosskerr_coupling,
        dipole_coupling=dipole_coupling,
        Pmin=Pmin,
        amplitude_bound=amplitude_bound,
        initialcondition=initialcondition,
        verbose=verbose,
    )

    config.runtype = RunType.OPTIMIZATION
    config.output_directory = output_directory

    # Set up target gate if provided
    if targetgate is not None:
        os.makedirs(output_directory, exist_ok=True)

        # Convert to numpy array with complex dtype
        gate_array = np.array(targetgate, dtype=complex)

        # Write gate to file (column-major order)
        gate_file = os.path.join(output_directory, "targetgate.dat")
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
        config.optim_target = target

    return config
