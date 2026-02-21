import os
import pytest
import numpy as np
from quandary.new import setup_quandary, optimize
from utils import assert_results_equal

# Mark all tests in this file as regression tests
pytestmark = pytest.mark.regression

EXPECTED_LENGTH = 1652
EXPECTED_INFIDELITY = 8.695620910992297e-06

EXPECTED_PT = [
    [
        -0.351017944518694, -0.9111853674706331, 1.56337057296288, 1.32102312178454, -2.17680392763229,
        -1.91465552269025, -2.40030147272628, 2.44449966044797, 2.09663827839874, -1.0145362400273399
    ],
]

EXPECTED_QT = [
    [
        -0.35489683713910297, -2.8226401868623903, -2.64838539616493, -1.78300097541017, -2.75117332619166,
        -2.7605889552578198, -2.5975353420470504, -2.74733785636068, -2.76505669888879, -2.09765274564999
    ],
]

EXPECTED_ENERGY = [
    [
        [
            0.0, 0.00692816781265531, 0.0328301908913508, 0.0695204206507133, 0.115813024750246,
            0.172087630256885, 0.251613716534252, 0.336745657239306, 0.426441526468596, 0.502139791825633
        ],
    ],
]

EXPECTED_POPULATION = [
    [
        [
            1.0, 0.993078206283365, 0.967173742450001, 0.930495435459176, 0.884259384596979,
            0.828002837889172, 0.748574947700958, 0.663432310188568, 0.573644729851791, 0.497863551562967
        ],
    ],
]

# Compare output to expected result for 10 points
NUM_SAMPLES = 10
SAMPLE_INDICES = [int(i * (EXPECTED_LENGTH - 1) / (NUM_SAMPLES - 1)) for i in range(NUM_SAMPLES)]


def test_example_state_to_state(tmp_path, request):
    """Test state-to-state preparation using new Python interface."""
    datadir_path = os.path.join(tmp_path, request.node.name)

    Ne = [2]
    Ng = [1]
    freq01 = [4.10595]
    selfkerr = [0.2198]
    T = 50.0
    maxctrl_MHz = 4.0
    initialcondition = [1.0, 0.0]
    targetstate = [1.0/np.sqrt(2), 1.0/np.sqrt(2)]
    n_osc = 1
    n_levels = 1

    setup = setup_quandary(
        nessential=Ne,
        nguard=Ng,
        transition_frequency=freq01,
        selfkerr=selfkerr,
        control_amplitude_bounds=[maxctrl_MHz/1000],
        initial_state=initialcondition,
        final_time=T,
        spline_knot_spacing=3.0,              # Old default: new default is 10 splines
        control_zero_boundary_condition=False, # Old default: new default is True
        output_directory=datadir_path,
        verbose=False,
    )

    # Match old interface defaults
    setup.rand_seed = 4321
    setup.optim_tol_final_cost = 1e-5  # tol_infidelity=1e-5 in old test
    setup.optim_penalty_leakage = 0.1
    setup.optim_penalty_energy = 0.1
    setup.optim_penalty_dpdm = 0.01
    setup.output_optimization_stride = 1
    setup.linearsolver_maxiter = 20

    # Match old initialization amplitude
    num_carrier_freqs = len(setup.carrier_frequencies[0])
    init_amplitude = 10.0 / 1000.0 / np.sqrt(2) / num_carrier_freqs

    results = optimize(
        setup,
        targetstate=targetstate,
        control_initialization_amplitude=init_amplitude,
        quiet=True,
    )

    assert_results_equal(
        t=results.time,
        pt=results.pt,
        qt=results.qt,
        infidelity=results.infidelity,
        energy=results.expected_energy,
        population=results.population,
        T=T,
        n_osc=n_osc,
        n_levels=n_levels,
        expected_length=EXPECTED_LENGTH,
        expected_infidelity=EXPECTED_INFIDELITY,
        expected_pt=EXPECTED_PT,
        expected_qt=EXPECTED_QT,
        expected_energy=EXPECTED_ENERGY,
        expected_population=EXPECTED_POPULATION,
        sample_indices=SAMPLE_INDICES,
    )
