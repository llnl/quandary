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


def downsample_pulses(*, pt0=None, qt0=None, nsplines, spline_knot_spacing, ntime, dt, nessential):
    """Downsample control pulses from high-resolution (pt0, qt0) to B-spline coefficients.


    This function is used to convert pulses defined at every time step to B-spline
    coefficients for piecewise constant (spline order 0) representation.

    Parameters:
    -----------
    pt0 : Sequence[ndarray]
        Real part of control pulses [MHz] for each oscillator. Size: (ntime+1,) per oscillator.
    qt0 : Sequence[ndarray]
        Imaginary part of control pulses [MHz] for each oscillator. Size: (ntime+1,) per oscillator.
    nsplines : int
        Number of B-spline basis functions.
    spline_knot_spacing : float
        Spacing between B-spline knots [ns].
    ntime : int
        Number of time steps.
    dt : float
        Time step size [ns].
    nessential : Sequence[int]
        Number of essential levels per oscillator.

    Returns:
    --------
    pcof0 : ndarray
        Control parameter vector (B-spline coefficients) in rad/ns.
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
                "downsample_pulses: size mismatch in pt0/qt0 -- each oscillator must have "
                "matching arrays of length >= 2."
            )
    elif len(pt0) > 0 or len(qt0) > 0:
        raise ValueError(
            f"downsample_pulses: pt0 and qt0 must each have one entry per oscillator "
            f"(expected {Nsys}, got pt0={len(pt0)}, qt0={len(qt0)})."
        )
    else:
        return np.array([])
