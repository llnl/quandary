"""Quantum operator construction and analysis utilities."""

import logging
import numpy as np

logger = logging.getLogger(__name__)


def lowering(n):
    """Lowering operator of dimension n."""
    return np.diag(np.sqrt(np.arange(1, n)), k=1)


def number(n):
    """Number operator of dimension n."""
    return np.diag(np.arange(n))


def map_to_oscillators(id, nessential, nguard):
    """Return the local energy level of each oscillator for a given global index id."""
    # len(nessential) = number of subsystems
    nlevels = [nessential[i] + nguard[i] for i in range(len(nessential))]
    localIDs = []

    index = int(id)
    for iosc in range(len(nessential)):
        postdim = np.prod(nlevels[iosc + 1:])
        localIDs.append(int(index / postdim))
        index = index % postdim

    return localIDs


def _eigen_and_reorder(H0, verbose=False):
    """Internal function that computes eigen decomposition.

    Re-orders it to make the eigenvector matrix as close to the identity as possible.
    """

    # Get eigenvalues and vectors and sort them in ascending order
    Ntot = H0.shape[0]
    evals, evects = np.linalg.eig(H0)
    reord = np.argsort(evals)
    evals = evals[reord]
    evects = evects[:, reord]

    # Find the column index corresponding to the largest element in each row of evects
    max_col = np.zeros(Ntot, dtype=np.int32)
    for row in range(Ntot):
        max_col[row] = np.argmax(np.abs(evects[row, :]))

    # test the error detection
    # max_col[1] = max_col[0]

    # loop over all columns and check max_col for duplicates
    Ndup_col = 0
    for row in range(Ntot - 1):
        for k in range(row + 1, Ntot):
            if max_col[row] == max_col[k]:
                Ndup_col += 1
                logger.error(f"Error: detected identical max_col = {max_col[row]} for rows {row} and {k}")

    if Ndup_col > 0:
        logger.error(f"Found {Ndup_col} duplicate column indices in max_col array")
        raise ValueError('Permutation of eigen-vector matrix failed')

    evects = evects[:, max_col]
    evals = evals[max_col]

    # Make sure all diagonal elements are positive
    for j in range(Ntot):
        if evects[j, j] < 0.0:
            evects[:, j] = - evects[:, j]

    return evals, evects


def get_resonances(*, nessential, nguard, Hsys, Hc_re=[], Hc_im=[], rotation_frequency=[], cw_amp_thres=1e-7,
                   cw_prox_thres=1e-2, verbose=True, stdmodel=True):
    """
    Computes system resonances, to be used as carrier wave frequencie.
    Returns resonance frequencies in GHz and corresponding growth rates.
    """

    if verbose:
        logger.info(f"\nComputing carrier frequencies, ignoring growth rate slower than: {cw_amp_thres} "
                    f"and frequencies closer than: {cw_prox_thres} [GHz]")

    nqubits = len(nessential)
    n = Hsys.shape[0]

    # Get eigenvalues of system Hamiltonian (GHz)
    Hsys_evals, Utrans = _eigen_and_reorder(Hsys, verbose)
    Hsys_evals = Hsys_evals.real  # Eigenvalues may have a small imaginary part due to numerical imprecision
    Hsys_evals = Hsys_evals / (2 * np.pi)

    # Look for resonances in the symmetric and anti-symmetric control Hamiltonians for each qubit
    resonances = []
    speed = []
    for q in range(nqubits):

        # Transform symmetric and anti-symmetric control Hamiltonians using eigenvectors (reordered)
        Hsym_trans = Utrans.conj().T @ Hc_re[q] @ Utrans
        Hanti_trans = Utrans.conj().T @ Hc_im[q] @ Utrans

        resonances_a = []
        speed_a = []
        if verbose:
            logger.info(f"  Resonances in oscillator # {q}")

        for Hc_trans in (Hsym_trans, Hanti_trans):

            # Iterate over non-zero elements in transformed control
            for i in range(n):
                # Only consider transitions from lower to higher levels
                for j in range(i):

                    # Look for non-zero elements (skip otherwise)
                    if abs(Hc_trans[i, j]) < 1e-14:
                        continue

                    # Get the resonance frequency
                    delta_f = Hsys_evals[i] - Hsys_evals[j]
                    if abs(delta_f) < 1e-10:
                        delta_f = 0.0

                    # Get involved oscillator levels
                    ids_i = map_to_oscillators(i, nessential, nguard)
                    ids_j = map_to_oscillators(j, nessential, nguard)

                    # make sure both indices correspond to essential energy levels
                    is_ess_i = all(ids_i[k] < nessential[k] for k in range(len(nessential)))
                    is_ess_j = all(ids_j[k] < nessential[k] for k in range(len(nessential)))

                    if (is_ess_i and is_ess_j):
                        # Ignore resonances that are too close by comparing to all previous resonances
                        if any(abs(delta_f - f) < cw_prox_thres for f in resonances_a):
                            if verbose:
                                logger.info(f"    Ignoring resonance from {ids_j} to {ids_i}, freq {delta_f}, "
                                            f"growth rate= {abs(Hc_trans[i, j])} "
                                            f"being too close to one that already exists.")
                        # Ignore resonances with growth rate smaller than user-defined threshold
                        elif abs(Hc_trans[i, j]) < cw_amp_thres:
                            if verbose:
                                logger.info(f"    Ignoring resonance from {ids_j} to {ids_i}, freq {delta_f}, "
                                            f"growth rate= {abs(Hc_trans[i, j])} growth rate is too slow.")
                        # Otherwise, add resonance to the list
                        else:
                            resonances_a.append(delta_f)
                            speed_a.append(abs(Hc_trans[i, j]))
                            if verbose:
                                logger.info(f"    Resonance from {ids_j} to {ids_i}, freq {delta_f}, "
                                            f"growth rate= {abs(Hc_trans[i, j])}")

        # Append resonances for this qubit to overall list
        resonances.append(resonances_a)
        speed.append(speed_a)

    # Prepare output for carrier frequencies (om) and growth_rate
    Nfreq = np.zeros(nqubits, dtype=int)
    om = [[0.0] for _ in range(nqubits)]
    growth_rate = [[] for _ in range(nqubits)]
    for q in range(len(resonances)):
        Nfreq[q] = max(1, len(resonances[q]))  # at least one being 0.0
        om[q] = np.zeros(Nfreq[q])
        if len(resonances[q]) > 0:
            om[q] = np.array(resonances[q])
        growth_rate[q] = np.ones(Nfreq[q])
        if len(speed[q]) > 0:
            growth_rate[q] = np.array(speed[q])

    return om, growth_rate


def hamiltonians(*, N, transition_frequency, selfkerr, crosskerr_coupling=[], dipole_coupling=[],
                 rotation_frequency=[], verbose=True):
    """
    Create standard Hamiltonian operators to model pulse-driven superconducting qubits.

    Parameters:
    -----------
    N                    :  Total levels per oscillator (essential plus guard levels for each qubit)
    transition_frequency :  Groundstate frequency for each qubit [GHz]
    selfkerr             :  Self-kerr coefficient for each qubit [GHz]

    Optional Parameters:
    ---------------------
    crosskerr_coupling   : Cross-Kerr coupling strength [GHz]
    dipole_coupling      : dipole-dipole coupling strength [GHz]
    rotation_frequency   : Rotational frequencies for each qubit
    verbose              : Switch to turn on more output. Default: True

    Returns:
    --------
    Hsys   : System Hamiltonian (time-independent), units rad/ns
    Hc_re  : Real parts of control Hamiltonian operators for each qubit. Unit-less.
    Hc_im  : Imag parts of control Hamiltonian operators for each qubit. Unit-less.
    """

    if len(rotation_frequency) == 0:
        rotation_frequency = np.zeros(len(N))

    nqubits = len(N)
    assert len(selfkerr) == nqubits
    assert len(transition_frequency) == nqubits
    assert len(rotation_frequency) == nqubits

    n = np.prod(N)  # System size

    # Set up lowering operators in full dimension
    Amat = []
    for i in range(len(N)):
        ai = lowering(N[i])
        for j in range(i):
            ai = np.kron(np.identity(N[j]), ai)
        for j in range(i + 1, len(N)):
            ai = np.kron(ai, np.identity(N[j]))
        Amat.append(ai)

    # Set up system Hamiltonian: Duffing oscillators
    Hsys = np.zeros((n, n))
    for q in range(nqubits):
        domega_radns = 2.0 * np.pi * (transition_frequency[q] - rotation_frequency[q])
        selfkerr_radns = 2.0 * np.pi * selfkerr[q]
        Hsys += domega_radns * Amat[q].T @ Amat[q]
        Hsys -= selfkerr_radns / 2.0 * Amat[q].T @ Amat[q].T @ Amat[q] @ Amat[q]

    # Add cross cerr coupling, if given
    if len(crosskerr_coupling) > 0:
        idkl = 0
        for q in range(nqubits):
            for p in range(q + 1, nqubits):
                if abs(crosskerr_coupling[idkl]) > 1e-14:
                    crosskerr_radns = 2.0 * np.pi * crosskerr_coupling[idkl]
                    Hsys -= crosskerr_radns * Amat[q].T @ Amat[q] @ Amat[p].T @ Amat[p]
                idkl += 1

    # Add dipole_coupling coupling term.
    # Note that if the rotating frame frequencies are different amongst oscillators, then this contributes
    # to a *time-dependent* system Hamiltonian. Here, we treat this as time-independent, because this
    # Hamiltonian here is *ONLY* used to compute the time-step size and resonances, and it is NOT passed
    # to the quandary code. Quandary sets up the standard model with a time-dependent system Hamiltonian
    # if the frequencies of rotation differ amongst oscillators.
    if len(dipole_coupling) > 0:
        idkl = 0
        for q in range(nqubits):
            for p in range(q + 1, nqubits):
                if abs(dipole_coupling[idkl]) > 1e-14:
                    dipole_coupling_radns = 2.0 * np.pi * dipole_coupling[idkl]
                    Hsys += dipole_coupling_radns * (Amat[q].T @ Amat[p] + Amat[q] @ Amat[p].T)
                idkl += 1

    # Set up control Hamiltonians
    Hc_re = [Amat[q] + Amat[q].T for q in range(nqubits)]
    Hc_im = [Amat[q] - Amat[q].T for q in range(nqubits)]

    if verbose:
        logger.info(f"*** {nqubits} coupled quantum systems setup ***")
        logger.info(f"System Hamiltonian frequencies [GHz]: f01 = {transition_frequency}, "
                    f"rot. freq = {rotation_frequency}")
        logger.info(f"Selfkerr= {selfkerr}")
        logger.info(f"Coupling: X-Kerr= {crosskerr_coupling}, J-C= {dipole_coupling}")

    return Hsys, Hc_re, Hc_im
