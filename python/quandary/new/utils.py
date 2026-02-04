"""General utility functions for Quandary."""

import logging
import os
import shutil
import tempfile

import numpy as np

from .._quandary_impl import Config, run
from .results import get_results

logger = logging.getLogger(__name__)


def eval_controls(config, pcof, points_per_ns=1.0, output_directory=None, quiet=True):
    """Evaluate B-spline control parameters at a specific sample rate.

    This function evaluates B-spline coefficients (pcof) at a different time
    resolution than the optimization/simulation used. Useful for deploying
    optimized controls to hardware with different sample rates.

    Parameters
    ----------
    config : Config
        Validated configuration from a previous run. Contains all system and
        control settings needed for B-spline evaluation. The original config
        is not modified - a copy is made internally.
    pcof : ndarray
        B-spline coefficients to evaluate [rad/ns]. Typically from results.pcof.
    points_per_ns : float
        Sample rate for output [points per ns]. Default: 1.0 (1ns spacing).
        Examples: 1.0 = 1ns, 10.0 = 0.1ns, 0.5 = 2ns spacing.
    output_directory : str, optional
        Directory for output files. If None, creates a temporary directory.
    quiet : bool
        Suppress console output. Default: True.

    Returns
    -------
    time : ndarray
        Time points [ns] at the specified sample rate.
    pt : list of ndarray
        Real part of control pulses [MHz] for each oscillator.
    qt : list of ndarray
        Imaginary part of control pulses [MHz] for each oscillator.

    Example
    -------
    >>> # Optimize at fine timestep
    >>> results = run(opt_config)
    >>> # Evaluate at hardware sample rate (1ns)
    >>> time, pt, qt = eval_controls(results.config, results.pcof, points_per_ns=1.0)
    >>> # Deploy pt, qt to hardware AWG
    """
    # Set up output directory
    cleanup_temp = False
    if output_directory is None:
        output_directory = tempfile.mkdtemp(prefix='quandary_eval_')
        cleanup_temp = True
    else:
        os.makedirs(output_directory, exist_ok=True)

    # Write pcof to file for C++ to read
    pcof_file = os.path.join(output_directory, 'init_params.dat')
    np.savetxt(pcof_file, pcof, fmt='%.14e')

    # Create a copy of the config, then modify the copy
    eval_config = Config(config)
    eval_config.setup_for_eval_controls(points_per_ns, pcof_file, output_directory)

    # Run in EVALCONTROLS mode (just evaluates and writes controls)
    return_code = run(eval_config, quiet=quiet)

    if return_code != 0:
        raise RuntimeError(f"eval_controls failed with return code {return_code}")

    # Load results
    results = get_results(datadir=output_directory)

    # Clean up temporary directory if we created one
    if cleanup_temp:
        try:
            shutil.rmtree(output_directory)
        except Exception as e:
            logger.warning(f"Could not remove temporary directory {output_directory}: {e}")

    return results.time, results.pt, results.qt


def infidelity_(A, B):
    """Calculate infidelity between quantum states."""
    dim = int(np.sqrt(A.size))
    return 1.0 - np.abs(np.trace(A.conj().transpose() @ B))**2 / dim**2


def downsample_pulses(*, pt0=[], qt0=[], nsplines, spline_knot_spacing, nsteps, dT, Ne):
    """
    Downsample control pulses from high-resolution (pt0, qt0) to B-spline coefficients.

    This function is used to convert pulses defined at every time step to B-spline
    coefficients for piecewise constant (spline order 0) representation.

    Parameters:
    -----------
    pt0 : List[ndarray]
        Real part of control pulses [MHz] for each oscillator. Size: (nsteps+1,) per oscillator.
    qt0 : List[ndarray]
        Imaginary part of control pulses [MHz] for each oscillator. Size: (nsteps+1,) per oscillator.
    nsplines : int
        Number of B-spline basis functions.
    spline_knot_spacing : float
        Spacing between B-spline knots [ns].
    nsteps : int
        Number of time steps.
    dT : float
        Time step size [ns].
    Ne : List[int]
        Number of essential levels per oscillator.

    Returns:
    --------
    pcof0 : ndarray
        Control parameter vector (B-spline coefficients) in rad/ns.
    """
    Nsys = len(Ne)
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
                dt = (nsteps * dT) / (Nelem - 1)  # time step corresponding to (pt0, qt0)
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
