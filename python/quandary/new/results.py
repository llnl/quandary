"""Quandary results parsing and data structures."""

from __future__ import annotations

import glob
import logging
import os
from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np

from .._quandary_impl import Config, DecoherenceType

logger = logging.getLogger(__name__)


@dataclass
class Results:
    """Results from a Quandary simulation or optimization.

    Attributes:
        config: Validated configuration with all defaults applied. This is always present
            and represents the exact configuration that was used to generate these results.
        time: Array of time points [ns].
        pt: Control pulses p(t) for each oscillator [MHz]. Access: pt[oscillator][time_index].
        qt: Control pulses q(t) for each oscillator [MHz]. Access: qt[oscillator][time_index].
        ft: Lab-frame control pulses f(t) for each oscillator [MHz]. Access: ft[oscillator][time_index].
        uT: Evolved states at final time T. This is the (unitary) solution operator if
            the initial conditions span the full basis. Access: uT[:, initial_condition].
        pcof: Control parameters (B-spline coefficients).
        infidelity: Final infidelity (1 - fidelity). Only meaningful for
            optimization runs; defaults to 1.0 for simulations.
        optim_hist: Optimization history with keys: 'iter', 'objective', 'gradient',
            'ls_step', 'fidelity', 'cost', 'tikhonov', 'penalty', 'state_variation',
            'energy', 'control_variation'. Empty dict for simulation runs.
        expected_energy: Expected energy evolution for each oscillator and initial condition.
            Access: expected_energy[oscillator][initial_condition][time_index].
        population: Population evolution for each oscillator and initial condition.
            Access: population[oscillator][initial_condition][level, time_index].
    """

    config: Config
    time: np.ndarray = field(default_factory=lambda: np.array([]))
    pt: List[np.ndarray] = field(default_factory=list)
    qt: List[np.ndarray] = field(default_factory=list)
    ft: List[np.ndarray] = field(default_factory=list)
    uT: np.ndarray = field(default_factory=lambda: np.array([]))
    pcof: np.ndarray = field(default_factory=lambda: np.array([]))
    infidelity: float = 1.0
    optim_hist: Dict[str, np.ndarray] = field(default_factory=dict)
    expected_energy: List[List[np.ndarray]] = field(default_factory=list)
    population: List[List[np.ndarray]] = field(default_factory=list)


def get_results(config: Config) -> Results:
    """Load results from Quandary output files.

    This function parses output files from a Quandary run and returns them
    in a structured format. The output directory is read from the config.

    Args:
        config: Validated configuration object from the run.

    Returns:
        Results containing all parsed output data and config.

    Example:
        >>> # After running
        >>> validated_config = Config(my_config, quiet=False)
        >>> results = get_results(validated_config)
        >>> print(f"Infidelity: {results.infidelity}")
        >>> print(results.config)
    """
    # Get output directory from config
    datadir = config.output_directory

    # Get parameters from config
    lindblad = config.decoherence_type != DecoherenceType.NONE
    n_init = config.n_initial_conditions

    # Create results object with the provided config
    results = Results(config=config)

    # Detect from files
    control_files = sorted(glob.glob(os.path.join(datadir, "control*.dat")))
    rho_files = sorted(glob.glob(os.path.join(datadir, "rho_Re.iinit*.dat")))
    n_osc = len(control_files)

    # For Lindblad, we only want diagonal initial conditions for some outputs
    n_init_diag = n_init if not lindblad else int(np.sqrt(n_init))

    # Read control parameters (params.dat)
    params_file = os.path.join(datadir, "params.dat")
    if os.path.exists(params_file):
        try:
            results.pcof = np.loadtxt(params_file)
        except (OSError, ValueError) as e:
            logger.warning(f"Failed to read control parameters from {params_file}: {e}")

    # Read optimization history (optim_history.dat)
    # Column names are read from the file header to avoid hardcoding indices.
    # Expected columns (from output.cpp): iter, Objective, ||Pr(grad)||, LS step,
    # F_avg, Terminal cost, Tikhonov-regul, Penalty-term, State variation,
    # Energy-term, Control variation
    _OPTIM_HIST_KEYS = [
        "iter", "objective", "gradient", "ls_step", "fidelity",
        "cost", "tikhonov", "penalty", "state_variation", "energy",
        "control_variation",
    ]
    optim_file = os.path.join(datadir, "optim_history.dat")
    if os.path.exists(optim_file):
        try:
            data = np.loadtxt(optim_file)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            results.optim_hist = {
                key: data[:, i]
                for i, key in enumerate(_OPTIM_HIST_KEYS)
                if i < data.shape[1]
            }
            # Infidelity from last iteration
            fidelity_col = _OPTIM_HIST_KEYS.index("fidelity")
            results.infidelity = 1.0 - data[-1, fidelity_col]
        except (OSError, ValueError, IndexError) as e:
            logger.warning(
                f"Failed to read optimization history from {optim_file}: {e}"
            )

    # Read control pulses for each oscillator
    # Convert from GHz (file) to MHz (output)
    ghz_to_mhz = 1e3
    for iosc in range(n_osc):
        ctrl_file = os.path.join(datadir, f"control{iosc}.dat")
        if os.path.exists(ctrl_file):
            try:
                data = np.loadtxt(ctrl_file)
                if iosc == 0:
                    results.time = data[:, 0]
                results.pt.append(data[:, 1] * ghz_to_mhz)
                results.qt.append(data[:, 2] * ghz_to_mhz)
                results.ft.append(data[:, 3] * ghz_to_mhz)
            except (OSError, ValueError, IndexError) as e:
                logger.warning(f"Failed to read control pulses from {ctrl_file}: {e}")

    # Read expected energy for each oscillator and initial condition
    for iosc in range(n_osc):
        osc_expected = []
        for iinit in range(n_init_diag):
            # For Lindblad, only read diagonal initial conditions
            iid = iinit if not lindblad else iinit * n_init_diag + iinit
            filename = os.path.join(datadir, f"expected{iosc}.iinit{iid:04d}.dat")
            if os.path.exists(filename):
                try:
                    data = np.loadtxt(filename)
                    osc_expected.append(data[:, 1])
                except (OSError, ValueError, IndexError) as e:
                    logger.warning(
                        f"Failed to read expected energy from {filename}: {e}"
                    )
        if osc_expected:
            results.expected_energy.append(osc_expected)

    # Read population for each oscillator and initial condition
    for iosc in range(n_osc):
        osc_population = []
        for iinit in range(n_init_diag):
            iid = iinit if not lindblad else iinit * n_init_diag + iinit
            filename = os.path.join(datadir, f"population{iosc}.iinit{iid:04d}.dat")
            if os.path.exists(filename):
                try:
                    data = np.loadtxt(filename)
                    # Population data: first column is time, rest are level populations
                    osc_population.append(data[:, 1:].T)
                except (OSError, ValueError, IndexError) as e:
                    logger.warning(f"Failed to read population from {filename}: {e}")
        if osc_population:
            results.population.append(osc_population)

    # Read final state/density matrix for each initial condition
    # First, determine dimensions from the first rho file
    if rho_files and n_init > 0:
        try:
            # Read first file to get dimensions
            first_rho = np.loadtxt(rho_files[0], skiprows=1)
            ndim = first_rho.shape[1] - 1  # First column is time

            # Initialize uT array
            results.uT = np.zeros((ndim, n_init), dtype=complex)

            # Read each initial condition
            for iinit in range(n_init):
                file_index = f"{iinit:04d}"
                re_file = os.path.join(datadir, f"rho_Re.iinit{file_index}.dat")
                im_file = os.path.join(datadir, f"rho_Im.iinit{file_index}.dat")

                if os.path.exists(re_file):
                    try:
                        re_data = np.loadtxt(re_file, skiprows=1)
                        # Take last time step, skip time column
                        results.uT[:, iinit] = re_data[-1, 1:]
                    except (OSError, ValueError, IndexError) as e:
                        logger.warning(
                            f"Failed to read real part of rho from {re_file}: {e}"
                        )

                if os.path.exists(im_file):
                    try:
                        im_data = np.loadtxt(im_file, skiprows=1)
                        results.uT[:, iinit] += 1j * im_data[-1, 1:]
                    except (OSError, ValueError, IndexError) as e:
                        logger.warning(
                            f"Failed to read imaginary part of rho from {im_file}: {e}"
                        )
        except (OSError, ValueError, IndexError) as e:
            logger.warning(f"Failed to determine dimensions from {rho_files[0]}: {e}")

    return results
