"""Time step estimation utilities for Quandary simulations."""

import numpy as np


def estimate_timesteps(*, T=1.0, Hsys=[], Hc_re=[], Hc_im=[], maxctrl_MHz=[], Pmin=40):
    """
    Helper function to estimate the number of time steps based on eigenvalues of the system Hamiltonian and maximum control Hamiltonians. Note: The estimate does not account for quickly varying signals or a large number of splines. Double check that at least 2-3 points per spline are present to resolve control function. #TODO: Automate this.
    """

    # Get estimated control pulse amplitude
    est_ctrl_MHz = maxctrl_MHz[:]
    if len(maxctrl_MHz) == 0:
        est_ctrl_MHz = [10.0 for _ in range(max(len(Hc_re), len(Hc_im)))]

    # Set up Hsys +  maxctrl*Hcontrol
    K1 = np.copy(Hsys)

    for i in range(len(Hc_re)):
        est_radns = est_ctrl_MHz[i]*2.0*np.pi/1e+3
        if len(Hc_re[i])>0:
            K1 += est_radns * Hc_re[i]
    for i in range(len(Hc_im)):
        est_radns = est_ctrl_MHz[i]*2.0*np.pi/1e+3
        if len(Hc_im[i])>0:
            K1 = K1 + 1j * est_radns * Hc_im[i] # can't use += due to type!

    # Estimate time step
    eigenvalues = np.linalg.eigvals(K1)
    maxeig = np.max(np.abs(eigenvalues))
    # ctrl_fac = 1.2  # Heuristic, assuming that the total Hamiltonian is dominated by the system part.
    ctrl_fac = 1.0
    samplerate = ctrl_fac * maxeig * Pmin / (2 * np.pi)
#     print(f"{samplerate=}")
    nsteps = int(np.ceil(T * samplerate))

    return nsteps


def timestep_richardson_est(quandary, tol=1e-8, order=2, quandary_exec=""):
    """Decrease timestep size until Richardson error estimate meets threshold."""

    # Factor by which timestep size is decreased
    m = 2

    # Initialize
    quandary.verbose=False
    t, pt, qt, infidelity, _, _ = quandary.simulate(quandary_exec=quandary_exec, datadir="TS_test")

    Jcurr = infidelity
    uT = quandary.uT.copy()

    # Loop
    errs_J = []
    errs_u = []
    dts = []
    for i in range(10):

        # Update configuration number of timesteps. Note: dt will be set in dump()
        dt_org = quandary.T / quandary.nsteps
        quandary.nsteps= quandary.nsteps* m
        quandary.dT = quandary.T / quandary.nsteps

        # Get u(dt/m)
        quandary.verbose=False
        t, pt, qt, infidelity, _, _ = quandary.simulate(quandary_exec=quandary_exec, datadir="TS_test")

        # Richardson error estimate
        err_J = np.abs(Jcurr - infidelity) / (m**order-1.0)
        # err_u = np.abs(uT[1,1]- quandary.uT[1,1]) / (m**order - 1.0)
        err_u = np.linalg.norm(np.subtract(uT, quandary.uT)) / (m**order - 1.0)
        errs_J.append(err_J)
        errs_u.append(err_u)
        dts.append(dt_org)

        # Output
        print(" -> Error at i=", i, ", dt = ", dt_org, ": err_J = ", err_J, " err_u=", err_u)


        # Stop if tolerance is reached
        if err_J < tol:
            print("\n -> Tolerance reached. N=", quandary.nsteps, ", dt=",dt_org)
            break

        # Update
        Jcurr = infidelity
        uT = np.copy(quandary.uT)

    return errs_J, errs_u, dts
