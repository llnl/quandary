"""Quandary results parsing and data structures."""

from __future__ import annotations
import glob
import logging
import os
from dataclasses import dataclass, field
from typing import Dict, List
import matplotlib.pyplot as plt
import numpy as np
from .._quandary_impl import Config, DecoherenceType

logger = logging.getLogger(__name__)

@dataclass
class Results:
    """Results from a Quandary simulation or optimization.

    Attributes
    ----------
    config : Config
        Validated configuration with all defaults applied. Always present;
        represents the exact configuration used to generate these results.
    time : ndarray
        Time points [ns].
    p_samples : list of ndarray
        Control pulses p(t) [MHz] per oscillator.
        Access: ``p_samples[oscillator][time_index]``.
    q_samples : list of ndarray
        Control pulses q(t) [MHz] per oscillator.
        Access: ``q_samples[oscillator][time_index]``.
    ft : list of ndarray
        Lab-frame control pulses f(t) [MHz] per oscillator.
        Access: ``ft[oscillator][time_index]``.
    uT : ndarray
        Evolved states at final time T. This is the (unitary) solution
        operator if the initial conditions span the full basis.
        Access: ``uT[:, initial_condition]``.
    spline_coefficients : ndarray
        Control parameters (B-spline coefficients).
    infidelity : float
        Final infidelity (1 - fidelity). Only meaningful if a target is set in the configuration. Defaults to 1.0
    optim_hist : dict of str to ndarray
        Optimization history. Keys: ``'iter'``, ``'objective'``,
        ``'gradient'``, ``'ls_step'``, ``'fidelity'``, ``'cost'``,
        ``'tikhonov'``, ``'penalty'``, ``'state_variation'``,
        ``'energy'``, ``'control_variation'``. Empty dict for simulations.
    expected_energy : list of list of ndarray
        Expected energy evolution per oscillator and initial condition.
        Access: ``expected_energy[oscillator][initial_condition][time_index]``.
    population : list of list of ndarray
        Population evolution per oscillator and initial condition.
        Access: ``population[oscillator][initial_condition][level, time_index]``.
    """
    config: Config
    time: np.ndarray = field(default_factory=lambda: np.array([]))
    p_samples: List[np.ndarray] = field(default_factory=list)
    q_samples: List[np.ndarray] = field(default_factory=list)
    ft: List[np.ndarray] = field(default_factory=list)
    uT: np.ndarray = field(default_factory=lambda: np.array([]))
    spline_coefficients: np.ndarray = field(default_factory=lambda: np.array([]))
    infidelity: float = 1.0
    optim_hist: Dict[str, np.ndarray] = field(default_factory=dict)
    expected_energy: List[List[np.ndarray]] = field(default_factory=list)
    population: List[List[np.ndarray]] = field(default_factory=list)


def get_results(config: Config) -> Results:
    """Load results from Quandary output files.

    Parses output files from a Quandary run and returns them in a
    structured format. The output directory is read from the config.

    Parameters
    ----------
    config : Config
        Validated configuration object from the run.

    Returns
    -------
    Results
        All parsed output data and config.

    Examples
    --------
    >>> config = Config.from_file("data_out/config_log.toml", quiet=True)
    >>> results = get_results(config)
    >>> print(f"Infidelity: {results.infidelity}")
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
            results.spline_coefficients = np.loadtxt(params_file)
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
                results.p_samples.append(data[:, 1] * ghz_to_mhz)
                results.q_samples.append(data[:, 2] * ghz_to_mhz)
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


def plot_pulse(results):
    """Plot the control pulse for all qubits.

    Parameters
    ----------
    results : Results
        Results from optimize(), simulate(), or evaluate_controls().
    """
    Ne = results.config.nessential
    time = results.time
    p_samples = results.p_samples
    q_samples = results.q_samples

    plt.figure()
    nrows = len(Ne)
    ncols = 1
    for iosc in range(len(Ne)):
        plt.subplot(nrows, ncols, iosc + 1)
        plt.plot(time, p_samples[iosc], "r", label="p(t)")
        plt.plot(time, q_samples[iosc], "b", label="q(t)")
        plt.xlabel('time (ns)')
        plt.ylabel('Drive strength [MHz]')
        maxp = max(np.abs(p_samples[iosc]))
        maxq = max(np.abs(q_samples[iosc]))
        plt.title('Qubit ' + str(iosc) + '\n max. drive ' + str(round(maxp, 1)) + ", " +
                  str(round(maxq, 1)) + " MHz")
        plt.legend(loc='lower right')
        plt.xlim([0.0, time[-1]])
    plt.subplots_adjust(hspace=0.6)
    plt.draw()
    plt.show()


def plot_expectedEnergy(results):
    """Plot evolution of expected energy levels.

    Parameters
    ----------
    results : Results
        Results from optimize(), simulate(), or evaluate_controls().
    """
    Ne = results.config.nessential
    time = results.time
    expectedEnergy = results.expected_energy

    ninit = len(expectedEnergy[0])
    nplots = ninit  # one plot for each initial state
    ncols = 2 if nplots >= 4 else 1  # 2 rows if more than 3 plots
    nrows = int(np.ceil(nplots / ncols))
    figsizex = 6.4 * nrows * 0.75
    figsizey = 4.8 * nrows * 0.75
    plt.figure(figsize=(figsizex, figsizey))
    for iplot in range(nplots):
        iinit = iplot
        plt.subplot(nrows, ncols, iplot + 1)
        emax = 1.0
        for iosc in range(len(Ne)):
            label = 'Qubit ' + str(iosc) if len(Ne) > 1 else ''
            plt.plot(time, expectedEnergy[iosc][iinit], label=label)
            emax_iosc = np.max(expectedEnergy[iosc][iinit])
            emax = max(emax, emax_iosc)  # keep track of max energy level for setting ylim
        plt.xlabel('time (ns)')
        plt.ylabel('expected energy')

        plt.ylim([0.0 - 1e-2, emax + 1e-2])
        plt.xlim([0.0, time[-1]])
        binary_ID = iplot if len(Ne) == 1 else bin(iplot).replace("0b", "").zfill(len(Ne))
        plt.title("from |" + str(binary_ID) + ">")
        plt.legend(loc='lower right')
    plt.subplots_adjust(hspace=0.5)
    plt.subplots_adjust(wspace=0.5)
    plt.draw()
    plt.show()


def plot_population(results):
    """Plot evolution of population.

    Parameters
    ----------
    results : Results
        Results from optimize(), simulate(), or evaluate_controls().
    """
    Ne = results.config.nessential
    time = results.time
    population = results.population

    ninit = len(population[0])
    nplots = ninit  # one plot for each initial state
    ncols = 2 if nplots >= 4 else 1  # 2 rows if more than 3 plots
    nrows = int(np.ceil(nplots / ncols))
    figsizex = 6.4 * nrows * 0.75
    figsizey = 4.8 * nrows * 0.75
    plt.figure(figsize=(figsizex, figsizey))

    # Iterate over initial conditions (one plot for each)
    for iplot in range(nplots):
        iinit = iplot
        plt.subplot(nrows, ncols, iplot + 1)
        for iosc in range(len(Ne)):
            for istate in range(Ne[iosc]):
                label = 'Qubit ' + str(iosc) if len(Ne) > 1 else ''
                label = label + " |" + str(istate) + ">"
                plt.plot(time, population[iosc][iinit][istate], label=label)
        plt.xlabel('time (ns)')
        plt.ylabel('population')
        plt.ylim([0.0 - 1e-4, 1.0 + 1e-2])
        plt.xlim([0.0, time[-1]])
        binary_ID = iplot if len(Ne) == 1 else bin(iplot).replace("0b", "").zfill(len(Ne))
        plt.title("from |" + str(binary_ID) + ">")
        plt.legend(loc='lower right')
    plt.subplots_adjust(hspace=0.5)
    plt.subplots_adjust(wspace=0.5)
    plt.draw()
    plt.show()


def plot_results_1osc(results, oscillator=0):
    """Plot all results of one oscillator (pulses, FFT, populations, expected energy).

    Parameters
    ----------
    results : Results
        Results from optimize(), simulate(), or evaluate_controls().
    oscillator : int
        Index of oscillator to plot. Default: 0.
    """
    t = results.time
    p = results.p_samples[oscillator]
    q = results.q_samples[oscillator]
    expectedEnergy = results.expected_energy[oscillator]
    population = results.population[oscillator]
    Ne = results.config.nessential[oscillator]

    fig, ax = plt.subplots(2, 3, figsize=(20, 8))
    fig.subplots_adjust(hspace=0.3)

    # Plot pulses
    ax[0, 0].plot(t, p, label='I')  # y label: MHz
    ax[0, 0].plot(t, q, label='Q')  # y label: MHz
    ax[0, 0].set_ylabel('Pulse amplitude (MHz)')
    ax[0, 0].set_xlabel('Time (ns)')
    ax[0, 0].legend()
    ax[0, 0].grid()

    # Compute and plot FFT
    zlist = np.array(p) * 1e-3 + 1j * np.array(q) * 1e-3
    fft = np.fft.fft(zlist)
    dt = results.config.dt
    fftfr = np.fft.fftfreq(len(zlist), d=dt)

    ax[0, 1].scatter(fftfr * 1e3, np.abs(fft)**2)
    ax[0, 1].set_ylabel('FFT')
    ax[0, 1].set_xlabel('Frequency (MHz)')
    ax[0, 1].grid()
    ax[0, 1].set_title('FFT')
    ax[0, 1].set_yscale('log')
    ax[0, 1].set_xlim(-500, 500)
    ax[0, 1].set_ylim(1e-8, 1e5)

    # Plot Populations for each initial condition
    for iinit in range(len(population)):  # for each of the initial states
        row = 1
        col = iinit

        for istate in range(Ne):  # for each essential level
            label = "|" + str(istate) + ">"
            ax[row, col].plot(t, population[iinit][istate], label=label)

        ax[row, col].set_xlabel('Time (ns)')
        ax[row, col].set_ylabel('Population')
        ax[row, col].legend()
        ax[row, col].set_title('Populations from |%d>' % iinit)
        ax[row, col].grid()

    # Plot expected Energy
    row, col = 0, 2
    for iinit in range(len(expectedEnergy)):
        label = 'from |' + str(iinit) + '>'
        ax[row, col].plot(t, expectedEnergy[iinit], label=label)
    ax[row, col].set_xlabel('Time (ns)')
    ax[row, col].set_ylabel('Expected Energy Level')
    ax[row, col].legend()
    ax[row, col].set_title('Expected Energy Level')
    ax[row, col].grid()

    plt.draw()
    plt.show()
