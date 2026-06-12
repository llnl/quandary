import numpy as np
import os
from pytest import approx
from quandary.new import create_config, simulate, evaluate_controls


def test_fit_bspline2nd(tmp_path, request):
    """
    Test 2nd-order Bsplines are fitted properly to existing controls and parameterization
    """
    datadir_path = os.path.join(tmp_path, request.node.name)

    T = 5.0
    setup = create_config(
        nessential=[2,2],
        transition_frequency=[4.0, 4.4],
        total_time=T,
        dt = 0.005,
        carrier_frequency=[[0.01, 0.02, 0.03], [0.22, 0.33, 0.44]],
        spline_order=2,
        spline_knot_spacing=0.1, # ~10 splines
        control_amplitude = 0.01,
        control_randomize = True,
        output_directory=datadir_path,
    )

    # simulate so to get spline parameters. 
    results_orig = simulate(setup, quiet=True)
    p_orig = results_orig.p_samples
    q_orig = results_orig.q_samples
    time_orig = results_orig.time

    # Fit 2nd-order Bspline to the original control parameters
    results_fitted = simulate(setup, p_samples=p_orig, q_samples=q_orig, quiet=True)

    # Compare original and fitted p/q_samples, time, and fidelity.
    assert np.allclose(results_orig.time, results_fitted.time, atol=1e-13), f"Expected time ~ {results_orig.time}, got {results_fitted.time}"

    # Loose tolerance for p/q_samples, they may not match exactly since the spline fitting is not exact. spline_coefficients can be different, not compared here.
    assert np.allclose(results_orig.p_samples, results_fitted.p_samples, atol=1e-2), f"Expected p_samples ~ {results_orig.p_samples}, got {results_fitted.p_samples}"
    assert np.allclose(results_orig.q_samples, results_fitted.q_samples, atol=1e-2), f"Expected q_samples ~ {results_orig.q_samples}, got {results_fitted.q_samples}"

    # Compare expected energy and population, which should be close but not exact. 
    assert np.allclose(results_orig.expected_energy, results_fitted.expected_energy, atol=1e-5), f"Expected expected_energy ~ {results_orig.expected_energy}, got {results_fitted.expected_energy}"
    assert np.allclose(results_orig.population, results_fitted.population, atol=1e-5), f"Expected populations ~ {results_orig.population}, got {results_fitted.population}"
