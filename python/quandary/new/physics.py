"""Quantum operator construction and analysis utilities."""

import logging
import numpy as np
from collections.abc import Sequence

logger = logging.getLogger(__name__)


# -----------------------------------------------
# Operators, Hamiltonians, model construction 
# -----------------------------------------------

def lowering(n):
    """Lowering operator of dimension n."""
    return np.diag(np.sqrt(np.arange(1, n)), k=1)


def number(n):
    """Number operator of dimension n."""
    return np.diag(np.arange(n))


def get_resonances(*, nessential, nguard, Hsys, Hc_re=None, Hc_im=None, rotation_frequency=None, cw_amp_thres=1e-7,
                   cw_prox_thres=1e-2, stdmodel=True):
    """Compute system resonances from Hamiltonian eigenvalues. Use to compute carrier wave frequencies.

    Parameters
    ----------
    nessential : sequence of int
        Number of essential energy levels per oscillator.
    nguard : sequence of int
        Number of guard levels per oscillator.
    Hsys : ndarray
        System Hamiltonian [rad/ns].
    Hc_re : sequence of ndarray
        Real parts of control Hamiltonian operators for each oscillator.
    Hc_im : sequence of ndarray
        Imaginary parts of control Hamiltonian operators for each oscillator.
    rotation_frequency : sequence of float
        Rotating frame frequencies [GHz] per oscillator. Default: zeros.
    cw_amp_thres : float
        Minimum growth rate to include a resonance. Resonances with a
        smaller coupling matrix element are discarded. Default: 1e-7.
    cw_prox_thres : float
        Minimum frequency separation [GHz] between retained resonances.
        Resonances closer than this to an existing one are discarded.
        Default: 0.01.
    stdmodel : bool
        Reserved for future use. Default: True.

    Returns
    -------
    om : sequence of ndarray
        Carrier wave frequencies [GHz] for each oscillator. At least one
        frequency (0.0) is returned per oscillator.
    growth_rate : sequence of ndarray
        Coupling strengths (growth rates) corresponding to each carrier
        frequency, for each oscillator.
    """

    if Hc_re is None:
        Hc_re = []
    if Hc_im is None:
        Hc_im = []
    if rotation_frequency is None:
        rotation_frequency = []

    logger.info(f"\nComputing carrier frequencies, ignoring growth rate slower than: {cw_amp_thres} "
                f"and frequencies closer than: {cw_prox_thres} [GHz]")

    nqubits = len(nessential)
    n = Hsys.shape[0]

    # Get eigenvalues of system Hamiltonian (GHz)
    Hsys_evals, Utrans = _eigen_and_reorder(Hsys)
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
                    ids_i = _map_to_oscillators(i, nessential, nguard)
                    ids_j = _map_to_oscillators(j, nessential, nguard)

                    # make sure both indices correspond to essential energy levels
                    is_ess_i = all(ids_i[k] < nessential[k] for k in range(len(nessential)))
                    is_ess_j = all(ids_j[k] < nessential[k] for k in range(len(nessential)))

                    if (is_ess_i and is_ess_j):
                        # Ignore resonances that are too close by comparing to all previous resonances
                        if any(abs(delta_f - f) < cw_prox_thres for f in resonances_a):
                            logger.info(f"    Ignoring resonance from {ids_j} to {ids_i}, freq {delta_f}, "
                                        f"growth rate= {abs(Hc_trans[i, j])} "
                                        f"being too close to one that already exists.")
                        # Ignore resonances with growth rate smaller than user-defined threshold
                        elif abs(Hc_trans[i, j]) < cw_amp_thres:
                            logger.info(f"    Ignoring resonance from {ids_j} to {ids_i}, freq {delta_f}, "
                                        f"growth rate= {abs(Hc_trans[i, j])} growth rate is too slow.")
                        # Otherwise, add resonance to the list
                        else:
                            resonances_a.append(delta_f)
                            speed_a.append(abs(Hc_trans[i, j]))
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


def hamiltonians(*, 
    nessential: Sequence[int], 
    transition_frequency: Sequence[float], 
    nguard: Sequence[int]=None, 
    selfkerr: Sequence[float]=None, 
    crosskerr_coupling: Sequence[float]=None, 
    dipole_coupling: Sequence[float]=None,
    rotation_frequency: Sequence[float]=None
) -> tuple:
    """Create standard Hamiltonian operators for pulse-driven superconducting qubits.

    Parameters
    ----------
    nessential : sequence of int, required
        Number of essential levels per oscillator.
    transition_frequency : sequence of float, required
        Ground-state transition frequency for each qubit [GHz].
    nguard : sequence of int, optional
        Number of guard levels per oscillator. Default: zeros.
    selfkerr : sequence of float, optional
        Self-Kerr (anharmonicity) coefficient for each qubit [GHz].
    crosskerr_coupling : sequence of float, optional
        Cross-Kerr coupling strengths [GHz]. Format: [g01, g02, ..., g12, ...].
        Default: no coupling.
    dipole_coupling : sequence of float, optional
        Dipole-dipole (Jaynes-Cummings) coupling strengths [GHz].
        Format: [J01, J02, ..., J12, ...]. Default: no coupling.
    rotation_frequency : sequence of float, optional
        Rotating frame frequencies for each qubit [GHz]. Default: zeros
        (lab frame).

    Returns
    -------
    Hsys : ndarray
        System Hamiltonian (time-independent) [rad/ns].
    Hc_re : sequence of ndarray
        Real (symmetric) control Hamiltonian operators for each qubit.
        Dimensionless; scaled by the control amplitude at runtime.
    Hc_im : sequence of ndarray
        Imaginary (anti-symmetric) control Hamiltonian operators for each
        qubit. Dimensionless; scaled by the control amplitude at runtime.
    """
    if nguard is None:
        nguard = [0 for _ in nessential]
    if selfkerr is None:
        selfkerr = []
    if crosskerr_coupling is None:
        crosskerr_coupling = []
    if dipole_coupling is None:
        dipole_coupling = []
    if rotation_frequency is None:
        rotation_frequency = []

    if len(rotation_frequency) == 0:
        rotation_frequency = np.zeros(len(nessential))
    if len(selfkerr) == 0:
        selfkerr = np.zeros(len(nessential))

    nqubits = len(nessential)
    assert len(selfkerr) == nqubits
    assert len(transition_frequency) == nqubits
    assert len(rotation_frequency) == nqubits

    Ntot = np.array(nessential) + np.array(nguard)
    n = np.prod(Ntot)  # System size

    # Set up lowering operators in full dimension
    Amat = []
    for i in range(len(Ntot)):
        ai = lowering(Ntot[i])
        for j in range(i):
            ai = np.kron(np.identity(Ntot[j]), ai)
        for j in range(i + 1, len(Ntot)):
            ai = np.kron(ai, np.identity(Ntot[j]))
        Amat.append(ai)

    # Set up system Hamiltonian: Duffing oscillators
    Hsys = np.zeros((n, n))
    for q in range(nqubits):
        domega_radns = 2.0 * np.pi * (transition_frequency[q] - rotation_frequency[q])
        selfkerr_radns = 2.0 * np.pi * selfkerr[q]
        Hsys += domega_radns * Amat[q].T @ Amat[q]
        Hsys -= selfkerr_radns / 2.0 * Amat[q].T @ Amat[q].T @ Amat[q] @ Amat[q]

    # Add cross-Kerr coupling, if given
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

    logger.info(f"*** {nqubits} coupled quantum systems ***")
    logger.info(f"System Hamiltonian frequencies [GHz]: f01 = {transition_frequency}, "
                f"rot. freq = {rotation_frequency}")
    logger.info(f"Selfkerr= {selfkerr}")
    logger.info(f"Coupling: X-Kerr= {crosskerr_coupling}, J-C= {dipole_coupling}")

    return Hsys, Hc_re, Hc_im


# -----------------------------------------------
# Fidelity metrics
# -----------------------------------------------

def gate_infidelity(A, B):
    """Calculate infidelity between two unitary operators or propagators.

    Parameters
    ----------
    A, B : ndarray of shape (n, n) or (m, n) with m >= n
        Unitary operators or propagators.  Rectangular matrices (e.g.
        nlevels x nessential when guard levels are present) are supported;
        the essential-level square block is extracted automatically.

    Returns
    -------
    float
        Infidelity ``1 - |Tr(A^H B)|^2 / d^2`` where *d* is the number
        of columns (essential-level dimension).
    """
    A = np.asarray(A)
    B = np.asarray(B)
    d = min(A.shape[0], A.shape[1])
    A = A[:d, :]
    B = B[:d, :]
    return 1.0 - np.abs(np.trace(A.conj().T @ B))**2 / d**2


def state_infidelity(psi, phi):
    """Calculate infidelity between two pure quantum states.

    Parameters
    ----------
    psi, phi : array-like, 1-D
        State vectors of the same dimension.

    Returns
    -------
    float
        Infidelity ``1 - |<psi|phi>|^2``.
    """
    psi = np.asarray(psi)
    phi = np.asarray(phi)
    return 1.0 - np.abs(np.vdot(psi, phi))**2


# -----------------------------------------------
# Control parameterization helpers 
# -----------------------------------------------

def fit_bspline0(*, p_samples=None, q_samples=None, nsplines, spline_knot_spacing, dt, nessential):
    """Fit 0-th order Bspline (piecewise constant) to given control pulses. 

    Parameters
    ----------
    p_samples : sequence of ndarray
        Real part of control pulses [MHz] per oscillator.
        Shape: ``(ntime+1,)`` per oscillator.
    q_samples : sequence of ndarray
        Imaginary part of control pulses [MHz] per oscillator.
        Shape: ``(ntime+1,)`` per oscillator.
    nsplines : int
        Number of B-spline basis functions.
    spline_knot_spacing : float
        Spacing between B-spline knots [ns].
    dt : float
        Time step size [ns].
    nessential : sequence of int
        Number of essential levels per oscillator.

    Returns
    -------
    spline_coefficients : ndarray
        Control parameter vector (B-spline 0-th order coefficients) in rad/ns.
    """
    if p_samples is None:
        p_samples = []
    if q_samples is None:
        q_samples = []

    Nsys = len(nessential)
    if len(p_samples) == Nsys and len(q_samples) == Nsys:
        sizes_ok = True
        for iosc in range(Nsys):
            if sizes_ok and len(p_samples[iosc]) >= 2 and len(p_samples[iosc]) == len(q_samples[iosc]):
                sizes_ok = True
            else:
                sizes_ok = False
        if sizes_ok:
            # do the downsampling and construct spline_coefficients
            spline_coefficients = np.zeros(0)  # to hold the downsampled numpy array for the control vector
            fact = 2e-3 * np.pi  # conversion factor from MHz to rad/ns

            for iosc in range(Nsys):
                Nelem = np.size(p_samples[iosc])
                p_seg = p_samples[iosc]
                q_seg = q_samples[iosc]

                seg_re = np.zeros(nsplines)  # to hold downsampled amplitudes
                seg_im = np.zeros(nsplines)
                # downsample p_seg, q_seg
                for i_spl in range(nsplines):
                    # the B-spline0 coefficients correspond to the time levels
                    t_spl = (i_spl) * spline_knot_spacing
                    # i = max(0, np.rint(t_spl/dt).astype(int))# given t_spl, find the closest time step index
                    i = int(np.rint(t_spl / dt))
                    i = min(i, Nelem - 1)  # make sure i is in range
                    seg_re[i_spl] = fact * p_seg[i]
                    seg_im[i_spl] = fact * q_seg[i]

                spline_coefficients = np.append(spline_coefficients, seg_re)  # append segment to the global control vector
                spline_coefficients = np.append(spline_coefficients, seg_im)
            return spline_coefficients
        else:
            raise ValueError(
                "fit_bspline0: size mismatch in p_samples/q_samples -- each oscillator must have "
                "matching arrays of length >= 2."
            )
    elif len(p_samples) > 0 or len(q_samples) > 0:
        raise ValueError(
            f"fit_bspline0: p_samples and q_samples must each have one entry per oscillator "
            f"(expected {Nsys}, got p_samples={len(p_samples)}, q_samples={len(q_samples)})."
        )
    else:
        return np.array([])


def fit_bspline2nd(
    t0,
    T,
    p_data,
    q_data,
    nsplines,
    carrier_frequencies=None,
    inputs_in_mhz=True,
    residual_rtol=1e-3,
    raise_on_bad_residual=True,
):
    """
    Fit 2nd order Bspline to given pulses p(t) and q(t). This matches Oscillator::evalControl for ControlType::BSPLINE:

        p(t) = sum_f [ cos(w_f t) * Bl1_f(t) - sin(w_f t) * Bl2_f(t) ]
        q(t) = sum_f [ sin(w_f t) * Bl1_f(t) + cos(w_f t) * Bl2_f(t) ]

    where Bl1_f(t) = sum_l alpha1[f,l] * B_l(t) and
          Bl2_f(t) = sum_l alpha2[f,l] * B_l(t).

    Coefficient layout in the returned vector matches Quandary's parameter
    layout for BSpline2nd with one control window (skip=0):

        [f0 alpha1(0..L-1), f0 alpha2(0..L-1),
         f1 alpha1(0..L-1), f1 alpha2(0..L-1), ...]

    Parameters
    ----------
    t0, T : float
        Start and stop times for the control window [ns].
    p_data, q_data : array-like
        Rotating-frame controls. Accepted shapes:
        - Single oscillator: (ntime,)
        - Multiple oscillators: (n_osc, ntime) or list of 1-D arrays
    nsplines : int or sequence of int
        Number of 2nd-order B-spline basis functions. If scalar, applied to
        all oscillators; if sequence, must have one value per oscillator.
    carrier_frequencies : sequence, optional
        Carrier frequencies [GHz]. Accepted forms:
        - Single oscillator: sequence of float
        - Multiple oscillators: sequence of sequence of float
        - Multiple oscillators with shared carriers: single sequence of float
          (applied to all oscillators)
        If None, defaults to [0.0] per oscillator.
    inputs_in_mhz : bool, optional
        If True (default), p_data/q_data are interpreted as MHz and converted to
        internal units (rad/ns) before fitting.
    residual_rtol : float, optional
        Relative residual tolerance for the least-squares fit. Default is 1e-3.
    raise_on_bad_residual : bool, optional
        If True (default), raise RuntimeError when the residual exceeds
        residual_rtol.
    """

    # Helper function to fit one oscillator's data to Bspline2nd coefficients
    def _fit_one_oscillator(t_vec, p_vec, q_vec, nspl, carriers):
        p_vec = np.asarray(p_vec, dtype=float).reshape(-1)
        q_vec = np.asarray(q_vec, dtype=float).reshape(-1)
        carriers = np.asarray(carriers, dtype=float).reshape(-1)

        if carriers.size == 0:
            raise ValueError("fit_bspline2nd: carrier_frequencies must be non-empty")
        if t_vec.size != p_vec.size or t_vec.size != q_vec.size:
            raise ValueError(
                "fit_bspline2nd: time grid, p_data, and q_data must have the same length"
            )
        if nspl < 3:
            raise ValueError("fit_bspline2nd: nsplines must be >= 3 for BSpline2nd")

        # Convert input p and q to rad/ns
        if inputs_in_mhz:
            mhz_to_radns = 2e-3 * np.pi
            p_vec = mhz_to_radns * p_vec
            q_vec = mhz_to_radns * q_vec

        # Convert carrier frequencies from GHz to rad/ns
        carriers = 2.0 * np.pi * carriers

        t0_local = float(t0)
        T_local = float(T)
        if T_local <= t0_local:
            raise ValueError("fit_bspline2nd: expected T > t0")

        # Basis functions from C++ Bspline2nd implementation
        dtknot = (T_local - t0_local) / (nspl - 2)
        width = 3.0 * dtknot
        tcenter = [t0_local + dtknot * ((i + 1) - 1.5) for i in range(nspl)]
        def basisfunction(l, t):
            tau = (t - tcenter[l]) / width
            if tau < -0.5 or tau >= 0.5:
                return 0.0
            if tau < -1 / 6:
                return 9 / 8 + 9 / 2 * tau + 9 / 2 * tau**2
            if tau < 1 / 6:
                return 3 / 4 - 9 * tau**2
            return 9 / 8 - 9 / 2 * tau + 9 / 2 * tau**2

        # Collocation matrix for spline basis values B_l(t_i), shape (N, L)
        B = np.array([[basisfunction(l, t) for l in range(nspl)] for t in t_vec])

        ntime = t_vec.size
        ncar = carriers.size
        nvars_per_carrier = 2 * nspl
        nvars = ncar * nvars_per_carrier

        # Build full linear system A x = y with y = [p; q], shape (2N,)
        A = np.zeros((2 * ntime, nvars), dtype=float)
        y = np.concatenate([p_vec, q_vec])

        for f_idx, w in enumerate(carriers):
            coswt = np.cos(w * t_vec)
            sinwt = np.sin(w * t_vec)

            base = f_idx * nvars_per_carrier
            a1_slice = slice(base, base + nspl)
            a2_slice = slice(base + nspl, base + 2 * nspl)

            # p rows
            A[:ntime, a1_slice] = coswt[:, None] * B
            A[:ntime, a2_slice] = -sinwt[:, None] * B
            # q rows
            A[ntime:, a1_slice] = sinwt[:, None] * B
            A[ntime:, a2_slice] = coswt[:, None] * B

        x, _, _, _ = np.linalg.lstsq(A, y, rcond=None)

        # Check relative residual quality of the fit.
        y_fit = A @ x
        residual_norm = np.linalg.norm(y_fit - y)
        y_norm = max(np.linalg.norm(y), np.finfo(float).eps)
        rel_residual = residual_norm / y_norm
        if rel_residual > residual_rtol:
            msg = (
                "fit_bspline2nd: poor least-squares fit quality "
                f"(relative residual={rel_residual:.3e} > tolerance={residual_rtol:.3e})"
            )
            if raise_on_bad_residual:
                raise RuntimeError(msg)
            logger.warning(msg)

        return x
    # end helper function

    t0 = float(t0)
    T = float(T)
    if T <= t0:
        raise ValueError("fit_bspline2nd: expected T > t0")

    p_arr = np.asarray(p_data, dtype=float)
    q_arr = np.asarray(q_data, dtype=float)

    if p_arr.shape != q_arr.shape:
        raise ValueError("fit_bspline2nd: p_data and q_data must have the same shape")

    # Single-oscillator case
    if p_arr.ndim == 1:
        t_data = np.linspace(t0, T, p_arr.size)
        carriers = [0.0] if carrier_frequencies is None else carrier_frequencies
        return _fit_one_oscillator(t_data, p_arr, q_arr, int(nsplines), carriers)

    # Multi-oscillator case: expect shape (n_osc, ntime)
    if p_arr.ndim != 2:
        raise ValueError(
            "fit_bspline2nd: p_data and q_data must be 1-D or 2-D (n_osc, ntime)"
        )

    t_data = np.linspace(t0, T, p_arr.shape[1])
    n_osc = p_arr.shape[0]

    if np.isscalar(nsplines):
        nspline_list = [int(nsplines)] * n_osc
    else:
        nspline_list = [int(v) for v in nsplines]
        if len(nspline_list) != n_osc:
            raise ValueError(
                f"fit_bspline2nd: nsplines must have length {n_osc}, got {len(nspline_list)}"
            )

    if carrier_frequencies is None:
        carrier_list = [[0.0] for _ in range(n_osc)]
    else:
        cf_arr = np.asarray(carrier_frequencies, dtype=object)
        if cf_arr.ndim == 1 and (len(cf_arr) == 0 or np.isscalar(cf_arr[0])):
            carrier_list = [carrier_frequencies for _ in range(n_osc)]
        else:
            carrier_list = list(carrier_frequencies)
            if len(carrier_list) != n_osc:
                raise ValueError(
                    f"fit_bspline2nd: carrier_frequencies must have length {n_osc}, got {len(carrier_list)}"
                )

    spline_coefficients_blocks = []
    for iosc in range(n_osc):
        spline_coefficients_blocks.append(
            _fit_one_oscillator(
                t_data,
                p_arr[iosc],
                q_arr[iosc],
                nspline_list[iosc],
                carrier_list[iosc],
            )
        )

    return np.concatenate(spline_coefficients_blocks)

# -----------------------------------------------
# Time step estimation and refinement
# -----------------------------------------------

def estimate_timestep_size(*, Hsys=None, Hc_re=None, Hc_im=None, control_amplitude_bound=None, Pmin=40):
    """Estimate the time steps size based on largest eigenvalue of the Hamiltonians.

    The estimate does not account for quickly varying control signals or a large number
    of splines. Verify that at least 2-3 time steps per spline are present to
    resolve the control function.

    Parameters
    ----------
    Hsys : ndarray
        System Hamiltonian [rad/ns].
    Hc_re : sequence of ndarray, optional
        Real parts of control Hamiltonian operators for each oscillator.
    Hc_im : sequence of ndarray, optional
        Imaginary parts of control Hamiltonian operators for each oscillator.
    control_amplitude_bound : sequence of float, optional
        Estimated max control amplitudes [GHz] per oscillator. Used to scale
        the control Hamiltonians when computing the largest eigenvalue.
        Default: 0.01 GHz per oscillator.
    Pmin : int
        Minimum number of time steps per period of the fastest oscillation.
        Default: 40.

    Returns
    -------
    suggested_dt : float
        Estimated timestep size [ns].
    """
    if Hsys is None:
        Hsys = []
    if Hc_re is None:
        Hc_re = []
    if Hc_im is None:
        Hc_im = []
    if control_amplitude_bound is None:
        control_amplitude_bound = []

    # Get estimated control pulse amplitude [GHz]
    est_control_amplitude_bound = control_amplitude_bound[:]
    if len(control_amplitude_bound) == 0:
        est_control_amplitude_bound = [0.01 for _ in range(max(len(Hc_re), len(Hc_im)))]

    # Set up Hsys +  maxctrl*Hcontrol
    K1 = np.copy(Hsys)

    for i in range(len(Hc_re)):
        est_radns = est_control_amplitude_bound[i] * 2.0 * np.pi
        if len(Hc_re[i]) > 0:
            K1 += est_radns * Hc_re[i]
    for i in range(len(Hc_im)):
        est_radns = est_control_amplitude_bound[i] * 2.0 * np.pi
        if len(Hc_im[i]) > 0:
            K1 = K1 + 1j * est_radns * Hc_im[i]  # can't use += due to type!

    # Estimate time step
    eigenvalues = np.linalg.eigvals(K1)
    maxeig = np.max(np.abs(eigenvalues))
    ctrl_fac = 1.0
    samplerate = ctrl_fac * maxeig * Pmin / (2 * np.pi)
    suggested_dt = 1.0 / samplerate

    return suggested_dt


def timestep_richardson_est(config_input, spline_coefficients, tol=1e-8, order=2, **kwargs):
    """Decrease timestep size until Richardson error estimate meets threshold.

    Parameters
    ----------
    config_input : ConfigInput
        Quandary configuration input. A copy is made internally; the caller's
        config_input is not modified.
    spline_coefficients : array-like
        B-spline control coefficients to evaluate at each refinement level.
    tol : float
        Richardson error tolerance on the infidelity. Default: 1e-8.
    order : int
        Time-stepping order for Richardson extrapolation. Default: 2.
    **kwargs
        Additional keyword arguments passed to simulate() (e.g. quiet, max_n_procs).

    Returns
    -------
    errs_J : sequence of float
        Richardson error estimates on the infidelity at each refinement step.
    errs_u : sequence of float
        Richardson error estimates on the propagator norm at each refinement step.
    dts : sequence of float
        Timestep sizes [ns] used at each refinement step.
    """
    from .run import simulate

    # Factor by which ntime is multiplied (dt halved) each step
    m = 2

    config_input = config_input.copy()
    kwargs.setdefault("quiet", True)

    results = simulate(config_input, spline_coefficients=spline_coefficients, **kwargs)
    Jcurr = results.infidelity
    uT = results.uT.copy()

    errs_J = []
    errs_u = []
    dts = []
    for i in range(10):
        dt_org = config_input.dt
        config_input.dt = config_input.dt / m

        results = simulate(config_input, spline_coefficients=spline_coefficients, **kwargs)

        err_J = np.abs(Jcurr - results.infidelity) / (m**order - 1.0)
        err_u = np.linalg.norm(np.subtract(uT, results.uT)) / (m**order - 1.0)
        errs_J.append(err_J)
        errs_u.append(err_u)
        dts.append(dt_org)

        logger.info(f"  i={i}, dt={dt_org:.6f}: err_J={err_J:.3e}, err_u={err_u:.3e}")

        if err_J < tol:
            logger.info(f"  Tolerance reached for dt={config_input.dt:.6f}")
            break

        Jcurr = results.infidelity
        uT = results.uT.copy()

    return errs_J, errs_u, dts


# --------------------------------
# Internal helper functions
# --------------------------------

def _map_to_oscillators(id, nessential, nguard):
    """Return the local level index of each oscillator for a given global state index.

    Parameters
    ----------
    id : int
        Global state index in the tensor-product Hilbert space.
    nessential : sequence of int
        Number of essential levels per oscillator.
    nguard : sequence of int
        Number of guard levels per oscillator.

    Returns
    -------
    localIDs : sequence of int
        Level index for each oscillator corresponding to the global index.
    """
    # len(nessential) = number of subsystems
    nlevels = [nessential[i] + nguard[i] for i in range(len(nessential))]
    localIDs = []

    index = int(id)
    for iosc in range(len(nessential)):
        postdim = np.prod(nlevels[iosc + 1:])
        localIDs.append(int(index / postdim))
        index = index % postdim

    return localIDs


def _eigen_and_reorder(H0):
    """Compute eigen decomposition, re-ordered so the eigenvector matrix
    is as close to the identity as possible.
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
