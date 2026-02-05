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
)
from .quantum_operators import hamiltonians, get_resonances
from .time_estimation import estimate_timesteps

logger = logging.getLogger(__name__)


def _setup_physics(
    Ne: List[int],
    freq01: List[float],
    T: float,
    selfkerr: Optional[List[float]] = None,
    Ng: Optional[List[int]] = None,
    rotfreq: Optional[List[float]] = None,
    crosskerr: Optional[List[float]] = None,
    Jkl: Optional[List[float]] = None,
    Pmin: int = 150,
    maxctrl_MHz: Optional[List[float]] = None,
    initialcondition: Optional[InitialConditionSettings] = None,
    verbose: bool = True,
) -> QuandaryConfig:
    """Internal: Build base config with common physics parameters.

    Computes Hamiltonians, timesteps, and carrier frequencies, then creates
    a QuandaryConfig with all common fields set. Runtype is not set here.

    Returns:
        QuandaryConfig with physics parameters configured.
    """
    # Set defaults
    nqubits = len(Ne)
    if selfkerr is None:
        selfkerr = [0.0] * nqubits
    if Ng is None:
        Ng = [0] * nqubits
    if rotfreq is None:
        rotfreq = freq01[:]
    if crosskerr is None:
        crosskerr = []
    if Jkl is None:
        Jkl = []
    if maxctrl_MHz is None:
        maxctrl_MHz = [10.0] * nqubits
    if initialcondition is None:
        initialcondition = InitialConditionSettings()
        initialcondition.condition_type = InitialConditionType.BASIS

    # Build Hamiltonians
    Ntot = [Ne[i] + Ng[i] for i in range(nqubits)]
    Hsys, Hc_re, Hc_im = hamiltonians(
        N=Ntot,
        freq01=freq01,
        selfkerr=selfkerr,
        crosskerr=crosskerr,
        Jkl=Jkl,
        rotfreq=rotfreq,
        verbose=verbose,
    )

    # Estimate timesteps
    nsteps = estimate_timesteps(
        T=T,
        Hsys=Hsys,
        Hc_re=Hc_re,
        Hc_im=Hc_im,
        maxctrl_MHz=maxctrl_MHz,
        Pmin=Pmin,
    )

    # Compute carrier frequencies
    carrier_frequency, _ = get_resonances(
        Ne=Ne,
        Ng=Ng,
        Hsys=Hsys,
        Hc_re=Hc_re,
        Hc_im=Hc_im,
        rotfreq=rotfreq,
        verbose=verbose,
    )

    if verbose:
        logger.info("Configuration computed:")
        logger.info(f"  Total time: {T} ns")
        logger.info(f"  Time steps: {nsteps}")
        logger.info(f"  dt: {T/nsteps:.6f} ns")
        logger.info(f"  Carrier frequencies: {carrier_frequency}")

    # Create config with common fields
    config = QuandaryConfig()
    config.nlevels = Ntot
    config.nessential = Ne
    config.ntime = nsteps
    config.dt = T / nsteps
    config.transfreq = freq01
    config.rotfreq = rotfreq
    config.selfkerr = selfkerr
    if len(crosskerr) > 0:
        config.crosskerr = crosskerr
    if len(Jkl) > 0:
        config.Jkl = Jkl

    config.carrier_frequencies = carrier_frequency
    config.initial_condition = initialcondition

    return config


def create_simulation_config(
    Ne: List[int],
    freq01: List[float],
    T: float,
    pcof = None,
    selfkerr: Optional[List[float]] = None,
    Ng: Optional[List[int]] = None,
    rotfreq: Optional[List[float]] = None,
    crosskerr: Optional[List[float]] = None,
    Jkl: Optional[List[float]] = None,
    Pmin: int = 150,
    maxctrl_MHz: Optional[List[float]] = None,
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
    Ne : List[int]
        Number of essential energy levels per qubit
    freq01 : List[float]
        01-transition frequencies [GHz] per qubit
    T : float
        Pulse duration [ns]

    Optional Parameters:
    -------------------
    selfkerr : List[float]
        Anharmonicities [GHz] per qubit. Default: zeros
    Ng : List[int]
        Number of guard levels per qubit. Default: zeros
    rotfreq : List[float]
        Frequency of rotations for computational frame [GHz] per qubit. Default: freq01
    crosskerr : List[float]
        ZZ coupling strength [GHz]. Format: [g01, g02, ..., g12, g13, ...]
    Jkl : List[float]
        Dipole-dipole coupling strength [GHz]. Format: [J01, J02, ..., J12, J13, ...]
    Pmin : int
        Minimum points to resolve shortest period (determines timesteps). Default: 150
    maxctrl_MHz : List[float]
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
    >>> # Simple simulation with random/default controls
    >>> config = create_simulation_config(Ne=[2], freq01=[4.1], T=50.0)
    >>> results = run(config)
    """
    config = _setup_physics(
        Ne=Ne,
        freq01=freq01,
        T=T,
        selfkerr=selfkerr,
        Ng=Ng,
        rotfreq=rotfreq,
        crosskerr=crosskerr,
        Jkl=Jkl,
        Pmin=Pmin,
        maxctrl_MHz=maxctrl_MHz,
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
        for iosc in range(len(Ne)):
            init = ControlInitializationSettings()
            init.init_type = ControlInitializationType.FILE
            init.filename = pcof_file
            control_inits.append(init)
        config.control_initializations = control_inits

    return config


def create_optimization_config(
    Ne: List[int],
    freq01: List[float],
    T: float,
    targetgate = None,
    output_directory: str = "./run_dir",
    selfkerr: Optional[List[float]] = None,
    Ng: Optional[List[int]] = None,
    rotfreq: Optional[List[float]] = None,
    crosskerr: Optional[List[float]] = None,
    Jkl: Optional[List[float]] = None,
    Pmin: int = 150,
    maxctrl_MHz: Optional[List[float]] = None,
    initialcondition: Optional[InitialConditionSettings] = None,
    verbose: bool = True,
) -> QuandaryConfig:
    """
    Create an optimization config with automatic Hamiltonian and timestep computation.

    Similar to create_simulation_config() but sets runtype to OPTIMIZATION.

    Required Parameters:
    -------------------
    Ne : List[int]
        Number of essential energy levels per qubit
    freq01 : List[float]
        01-transition frequencies [GHz] per qubit
    T : float
        Pulse duration [ns]
    targetgate : array-like (list, numpy array, etc.)
        Target unitary gate. Can be plain Python list like [[0,1],[1,0]]
        or numpy array. Will be converted to complex automatically.

    Optional Parameters:
    -------------------
    output_directory : str
        Output directory for results and gate files. Default: "./run_dir"
    selfkerr : List[float]
        Anharmonicities [GHz] per qubit. Default: zeros
    Ng : List[int]
        Number of guard levels per qubit. Default: zeros
    rotfreq : List[float]
        Frequency of rotations for computational frame [GHz] per qubit. Default: freq01
    crosskerr : List[float]
        ZZ coupling strength [GHz]
    Jkl : List[float]
        Dipole-dipole coupling strength [GHz]
    Pmin : int
        Minimum points to resolve shortest period. Default: 150
    maxctrl_MHz : List[float]
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
    ...     Ne=[2], freq01=[4.1], T=50.0,
    ...     targetgate=[[0, 1], [1, 0]]  # X gate
    ... )
    >>> results = run(config)
    """
    config = _setup_physics(
        Ne=Ne,
        freq01=freq01,
        T=T,
        selfkerr=selfkerr,
        Ng=Ng,
        rotfreq=rotfreq,
        crosskerr=crosskerr,
        Jkl=Jkl,
        Pmin=Pmin,
        maxctrl_MHz=maxctrl_MHz,
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
