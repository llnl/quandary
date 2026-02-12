"""General utility functions for Quandary."""

import logging

import numpy as np

logger = logging.getLogger(__name__)


def infidelity_(A, B):
    """Calculate infidelity between quantum states."""
    dim = int(np.sqrt(A.size))
    return 1.0 - np.abs(np.trace(A.conj().transpose() @ B))**2 / dim**2


def downsample_pulses(*, pt0=[], qt0=[], nsplines, spline_knot_spacing, ntime, dt, nessential):
    """
    Downsample control pulses from high-resolution (pt0, qt0) to B-spline coefficients.

    This function is used to convert pulses defined at every time step to B-spline
    coefficients for piecewise constant (spline order 0) representation.

    Parameters:
    -----------
    pt0 : List[ndarray]
        Real part of control pulses [MHz] for each oscillator. Size: (ntime+1,) per oscillator.
    qt0 : List[ndarray]
        Imaginary part of control pulses [MHz] for each oscillator. Size: (ntime+1,) per oscillator.
    nsplines : int
        Number of B-spline basis functions.
    spline_knot_spacing : float
        Spacing between B-spline knots [ns].
    ntime : int
        Number of time steps.
    dt : float
        Time step size [ns].
    nessential : List[int]
        Number of essential levels per oscillator.

    Returns:
    --------
    pcof0 : ndarray
        Control parameter vector (B-spline coefficients) in rad/ns.
    """
    Nsys = len(nessential)
    if len(pt0) == Nsys and len(qt0) == Nsys:
        sizes_ok = True
        for iosc in range(Nsys):
            if sizes_ok and len(pt0[iosc]) >= 2 and len(pt0[iosc]) == len(qt0[iosc]):
                sizes_ok = True
            else:
                sizes_ok = False
        # print("simulate(): sizes_ok = ", sizes_ok)
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
            # print("simulation(): downsampling of (pt0, qt0) completed")
            return pcof0
        else:
            print("simulation(): detected a mismatch in the sizes of (pt0, qt0)")
            return np.array([])
    elif len(pt0) > 0 or len(qt0) > 0:
        print("simulation(): the length of pt0 or qt0 != Nsys = ", Nsys)
        return np.array([])
    else:
        return np.array([])
