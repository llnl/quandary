"""Visualization utilities for Quandary results."""

import numpy as np
import matplotlib.pyplot as plt


def plot_pulse(results):
    """Plot the control pulse for all qubits.

    Args:
        results: Results object from run()
    """
    Ne = results.config.nessential
    time = results.time
    pt = results.pt
    qt = results.qt

    plt.figure()
    nrows = len(Ne)
    ncols = 1
    for iosc in range(len(Ne)):
        plt.subplot(nrows, ncols, iosc + 1)
        plt.plot(time, pt[iosc], "r", label="p(t)")
        plt.plot(time, qt[iosc], "b", label="q(t)")
        plt.xlabel('time (ns)')
        plt.ylabel('Drive strength [MHz]')
        maxp = max(np.abs(pt[iosc]))
        maxq = max(np.abs(qt[iosc]))
        plt.title('Qubit ' + str(iosc) + '\n max. drive ' + str(round(maxp, 1)) + ", " +
                  str(round(maxq, 1)) + " MHz")
        plt.legend(loc='lower right')
        plt.xlim([0.0, time[-1]])
    plt.subplots_adjust(hspace=0.6)
    plt.draw()
    plt.show()


def plot_expectedEnergy(results):
    """Plot evolution of expected energy levels.

    Args:
        results: Results object from run()
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

    Args:
        results: Results object from run()
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

    Args:
        results: Results object from run()
        oscillator: Index of oscillator to plot (default: 0)
    """
    t = results.time
    p = results.pt[oscillator]
    q = results.qt[oscillator]
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
