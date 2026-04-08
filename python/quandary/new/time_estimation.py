"""Time step estimation utilities for Quandary simulations."""

import logging
import numpy as np

logger = logging.getLogger(__name__)


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
