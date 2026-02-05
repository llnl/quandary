"""Factory functions to create QuandaryConfig with automatic setup from physics parameters."""

import logging
from typing import List, Optional

from .._quandary_impl import QuandaryConfig, RunType
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
    initialcondition: str = "basis",
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
    config.initial_condition.type = initialcondition

    return config


def create_simulation_config(
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
    targetgate: Optional[List[List[complex]]] = None,
    initialcondition: str = "basis",
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
    targetgate : List[List[complex]]
        Target gate for computing infidelity
    initialcondition : str
        Initial state specification. Default: "basis"
    verbose : bool
        Print setup information. Default: True

    Returns:
    -------
    QuandaryConfig
        Configured config ready for simulation. User can modify any fields before
        passing to run().

    Example:
    -------
    >>> config = create_simulation_config(Ne=[3, 3], freq01=[4.1, 4.2], T=100.0)
    >>> # Modify advanced settings with full autocomplete:
    >>> config.usematfree = False
    >>> config.output_frequency = 10
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

    return config


def create_optimization_config(
    Ne: List[int],
    freq01: List[float],
    T: float,
    targetgate: Optional[List[List[complex]]] = None,
    targetstate: Optional[List[complex]] = None,
    selfkerr: Optional[List[float]] = None,
    Ng: Optional[List[int]] = None,
    rotfreq: Optional[List[float]] = None,
    crosskerr: Optional[List[float]] = None,
    Jkl: Optional[List[float]] = None,
    Pmin: int = 150,
    maxctrl_MHz: Optional[List[float]] = None,
    initialcondition: str = "basis",
    verbose: bool = True,
) -> QuandaryConfig:
    """
    Create an optimization config with automatic Hamiltonian and timestep computation.

    Similar to create_simulation_config() but sets runtype to OPTIMIZATION and
    includes optimization-specific parameters.

    Required Parameters:
    -------------------
    Ne : List[int]
        Number of essential energy levels per qubit
    freq01 : List[float]
        01-transition frequencies [GHz] per qubit
    T : float
        Pulse duration [ns]

    At least one target must be specified:
    targetgate : List[List[complex]]
        Target unitary gate
    targetstate : List[complex]
        Target state vector

    Optional Parameters:
    -------------------
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
    initialcondition : str
        Initial state specification. Default: "basis"
    verbose : bool
        Print setup information. Default: True

    Returns:
    -------
    QuandaryConfig
        Configured config ready for optimization.

    Example:
    -------
    >>> import numpy as np
    >>> X_gate = np.array([[0, 1], [1, 0]], dtype=complex)
    >>> config = create_optimization_config(Ne=[2], freq01=[4.1], T=50.0, targetgate=X_gate)
    >>> config.tol_infidelity = 1e-6
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
    config.target_gate = targetgate
    config.target_state = targetstate

    return config
