import numpy as np
import os
from pytest import approx
from quandary.new import setup_quandary, evaluate_controls


def test_evaluate_controls_updates_timestep(tmp_path, request):
    """
    Test that evaluate_controls properly updates timestep (dt) to match sampling rate.
    """
    datadir_path = os.path.join(tmp_path, request.node.name)

    T = 5.0
    setup = setup_quandary(
        nessential=[2],
        transition_frequency=[4.0],
        final_time=T,
        output_directory=datadir_path,
        verbose=False,
    )

    original_dt = setup.dt
    original_ntime = setup.ntime

    # Test evaluate_controls with different sampling rate
    points_per_ns = 2
    results = evaluate_controls(setup, pcof=np.array([]), points_per_ns=points_per_ns, quiet=True)

    expected_nsteps = int(np.floor(T * points_per_ns))
    expected_dT = T / expected_nsteps

    assert results.time[0] == approx(0.0)
    assert results.time[-1] == approx(T)
    assert results.time[1] - results.time[0] == approx(expected_dT)

    # Verify original setup is unchanged (new interface copies setup internally)
    assert setup.dt == approx(original_dt)
    assert setup.ntime == original_ntime
