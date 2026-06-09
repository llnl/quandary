"""General utility functions for Quandary."""

import logging

import numpy as np

logger = logging.getLogger(__name__)


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


def fit_bspline0(*, pt0=None, qt0=None, nsplines, spline_knot_spacing, ntime, dt, nessential):
    """Fit 0-th order Bspline (piecewise constant) to given control pulses. 

    Parameters
    ----------
    pt0 : sequence of ndarray
        Real part of control pulses [MHz] per oscillator.
        Shape: ``(ntime+1,)`` per oscillator.
    qt0 : sequence of ndarray
        Imaginary part of control pulses [MHz] per oscillator.
        Shape: ``(ntime+1,)`` per oscillator.
    nsplines : int
        Number of B-spline basis functions.
    spline_knot_spacing : float
        Spacing between B-spline knots [ns].
    ntime : int
        Number of time steps.
    dt : float
        Time step size [ns].
    nessential : sequence of int
        Number of essential levels per oscillator.

    Returns
    -------
    pcof0 : ndarray
        Control parameter vector (B-spline 0-th order coefficients) in rad/ns.
    """
    if pt0 is None:
        pt0 = []
    if qt0 is None:
        qt0 = []

    Nsys = len(nessential)
    if len(pt0) == Nsys and len(qt0) == Nsys:
        sizes_ok = True
        for iosc in range(Nsys):
            if sizes_ok and len(pt0[iosc]) >= 2 and len(pt0[iosc]) == len(qt0[iosc]):
                sizes_ok = True
            else:
                sizes_ok = False
        if sizes_ok:
            # do the downsampling and construct pcof0
            pcof0 = np.zeros(0)  # to hold the downsampled numpy array for the control vector
            fact = 2e-3 * np.pi  # conversion factor from MHz to rad/ns

            for iosc in range(Nsys):
                Nelem = np.size(pt0[iosc])
                p_seg = pt0[iosc]
                q_seg = qt0[iosc]

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

                pcof0 = np.append(pcof0, seg_re)  # append segment to the global control vector
                pcof0 = np.append(pcof0, seg_im)
            return pcof0
        else:
            raise ValueError(
                "fit_bspline0: size mismatch in pt0/qt0 -- each oscillator must have "
                "matching arrays of length >= 2."
            )
    elif len(pt0) > 0 or len(qt0) > 0:
        raise ValueError(
            f"fit_bspline0: pt0 and qt0 must each have one entry per oscillator "
            f"(expected {Nsys}, got pt0={len(pt0)}, qt0={len(qt0)})."
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
    Fit 2nd order Bspline to given pulses p(t) and q(t).

    The model matches Oscillator::evalControl for ControlType::BSPLINE:

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
        - Single oscillator: ``(ntime,)``
        - Multiple oscillators: ``(n_osc, ntime)`` or list of 1-D arrays
    nsplines : int or sequence of int
        Number of 2nd-order B-spline basis functions. If scalar, applied to
        all oscillators; if sequence, must have one value per oscillator.
    carrier_frequencies : sequence, optional
        Carrier frequencies [GHz]. Accepted forms:
        - Single oscillator: sequence of float
        - Multiple oscillators: sequence of sequence of float
        - Multiple oscillators with shared carriers: single sequence of float
          (applied to all oscillators)
        If None, defaults to ``[0.0]`` per oscillator.
    inputs_in_mhz : bool, optional
        If True (default), p_data/q_data are interpreted as MHz and converted to
        internal units (rad/ns) before fitting.
    residual_rtol : float, optional
        Relative residual tolerance for the least-squares fit, measured as
        ``||A x - y|| / ||y||``. Default is ``5e-2``.
    raise_on_bad_residual : bool, optional
        If True, raise ``RuntimeError`` when the residual exceeds
        ``residual_rtol``. If False (default), emit a warning.

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

    pcof_blocks = []
    for iosc in range(n_osc):
        pcof_blocks.append(
            _fit_one_oscillator(
                t_data,
                p_arr[iosc],
                q_arr[iosc],
                nspline_list[iosc],
                carrier_list[iosc],
            )
        )

    return np.concatenate(pcof_blocks)



def estimate_timesteps(*, final_time=1.0, Hsys=None, Hc_re=None, Hc_im=None, control_amplitude_bounds=None, Pmin=40):
    """Estimate the number of time steps based on eigenvalues of Hamiltonians.

    The estimate does not account for quickly varying signals or a large number
    of splines. Verify that at least 2-3 time steps per spline are present to
    resolve the control function.

    Parameters
    ----------
    final_time : float
        Total simulation time [ns]. Default: 1.0.
    Hsys : ndarray
        System Hamiltonian [rad/ns].
    Hc_re : sequence of ndarray, optional
        Real parts of control Hamiltonian operators for each oscillator.
    Hc_im : sequence of ndarray, optional
        Imaginary parts of control Hamiltonian operators for each oscillator.
    control_amplitude_bounds : sequence of float, optional
        Estimated max control amplitudes [GHz] per oscillator. Used to scale
        the control Hamiltonians when computing the largest eigenvalue.
        Default: 0.01 GHz per oscillator.
    Pmin : int
        Minimum number of time steps per period of the fastest oscillation.
        Default: 40.

    Returns
    -------
    ntime : int
        Estimated number of time steps.
    """
    if Hsys is None:
        Hsys = []
    if Hc_re is None:
        Hc_re = []
    if Hc_im is None:
        Hc_im = []
    if control_amplitude_bounds is None:
        control_amplitude_bounds = []

    # Get estimated control pulse amplitude [GHz]
    est_control_amplitude_bounds = control_amplitude_bounds[:]
    if len(control_amplitude_bounds) == 0:
        est_control_amplitude_bounds = [0.01 for _ in range(max(len(Hc_re), len(Hc_im)))]

    # Set up Hsys +  maxctrl*Hcontrol
    K1 = np.copy(Hsys)

    for i in range(len(Hc_re)):
        est_radns = est_control_amplitude_bounds[i] * 2.0 * np.pi
        if len(Hc_re[i]) > 0:
            K1 += est_radns * Hc_re[i]
    for i in range(len(Hc_im)):
        est_radns = est_control_amplitude_bounds[i] * 2.0 * np.pi
        if len(Hc_im[i]) > 0:
            K1 = K1 + 1j * est_radns * Hc_im[i]  # can't use += due to type!

    # Estimate time step
    eigenvalues = np.linalg.eigvals(K1)
    maxeig = np.max(np.abs(eigenvalues))
    ctrl_fac = 1.0
    samplerate = ctrl_fac * maxeig * Pmin / (2 * np.pi)
    ntime = int(np.ceil(final_time * samplerate))

    return ntime


def timestep_richardson_est(setup, pcof, tol=1e-8, order=2, **kwargs):
    """Decrease timestep size until Richardson error estimate meets threshold.

    Parameters
    ----------
    setup : Setup
        Quandary setup configuration. A copy is made internally; the caller's
        setup is not modified.
    pcof : array-like
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
    from .runner import simulate

    # Factor by which ntime is multiplied (dt halved) each step
    m = 2

    setup = setup.copy()
    kwargs.setdefault("quiet", True)

    results = simulate(setup, pcof=pcof, **kwargs)
    Jcurr = results.infidelity
    uT = results.uT.copy()

    errs_J = []
    errs_u = []
    dts = []
    for i in range(10):
        dt_org = setup.dt
        setup.ntime = setup.ntime * m
        setup.dt = setup.dt / m

        results = simulate(setup, pcof=pcof, **kwargs)

        err_J = np.abs(Jcurr - results.infidelity) / (m**order - 1.0)
        err_u = np.linalg.norm(np.subtract(uT, results.uT)) / (m**order - 1.0)
        errs_J.append(err_J)
        errs_u.append(err_u)
        dts.append(dt_org)

        logger.info(f"  i={i}, dt={dt_org:.6f}: err_J={err_J:.3e}, err_u={err_u:.3e}")

        if err_J < tol:
            logger.info(f"  Tolerance reached: ntime={setup.ntime}, dt={setup.dt:.6f}")
            break

        Jcurr = results.infidelity
        uT = results.uT.copy()

    return errs_J, errs_u, dts
