import os
import pytest
import numpy as np
from quandary.new import setup_quandary, optimize
from utils import assert_results_equal


# Mark all tests in this file as regression tests
pytestmark = pytest.mark.regression

EXPECTED_LENGTH = 1222
EXPECTED_INFIDELITY = 9.967903512697607e-05

EXPECTED_PT = [
    [
        -1.1455758286673898, 7.30666099687124, -1.56254374996937, -2.36557441557185, 1.45983158306447,
        1.1348196056112798, 0.17027119825244802, 6.66935674855051, 3.83949098703669, -1.29770951519605
    ],
    [
        -0.0857364353164705, 0.26682652763546205, -2.84259629787855, 5.8730580454732495, -7.71845217474982,
        4.108702767526379, 4.31045393001853, -8.29308842920559, 2.51709250120214, 0.100284940848674
    ],
]

EXPECTED_QT = [
    [
        -0.102518227951347, 8.03543317276026, -4.12672901503846, -1.84670226508142, -2.3095502994813,
        2.02088269375347, -3.66417041158405, -0.41116620626056, 0.40459736675858904, 0.226771921733744
    ],
    [
        0.0852976327317638, -0.30191635506298403, 1.15730296070651, 1.98832750472986, 1.29849602408214,
        -7.59955466164568, 7.1300223152945605, -5.95519648260564, 2.36244366579429, -0.3462938582071
    ],
]

EXPECTED_ENERGY = [
    [
        [
            0.0, 0.00370220612622186, 0.0801332542644481, 0.15688463267758, 0.548903485801343,
            0.605435971064497, 0.588636450251953, 0.600128445340578, 0.171512211227719, 4.78498051761215e-05
        ],
        [
            0.0, 0.0249069977367134, 0.094112746268439, 0.210177591590367, 0.613751525451786,
            0.551293437285607, 0.627396931858309, 0.755186728733611, 0.201908774636426, 7.46873419305587e-05
        ],
        [
            1.0, 0.972984763890452, 0.887143848728789, 0.794289503333496, 0.443574516832106,
            0.416077263717592, 0.326537124922318, 0.257799553542498, 0.82670168506467, 0.999943362518778
        ],
        [
            1.0, 0.998406032199003, 0.938610150832428, 0.838648272670286, 0.393770472430779,
            0.427193328442698, 0.457429493477111, 0.386885273563193, 0.799877329825436, 0.999934100640022
        ],
    ],
    [
        [
            0.0, 0.00210990545133481, 0.00704882244254717, 0.244548352176713, 0.553119340084485,
            0.321745694464255, 0.477150970394771, 0.561363914356036, 0.0615657334786886, 8.22413835454486e-05
        ],
        [
            1.0, 0.974117273615089, 0.954893484094794, 0.747556963093925, 0.451943034138268,
            0.672297699308123, 0.517590897628196, 0.43531927925783, 0.922195291055279, 0.999899924632047
        ],
        [
            0.0, 0.072458882749064, 0.187064188322364, 0.329675230692441, 0.540287018259676,
            0.648822763889229, 0.379007502887368, 0.597276056588351, 0.980989706070835, 0.99992765780649
        ],
        [
            1.0, 0.951313938373316, 0.850993505339821, 0.678219454514544, 0.45465060808428,
            0.357133842978274, 0.626250630042886, 0.4060407503267, 0.0352492695595698, 9.01763846159662e-05
        ],
    ],
]

EXPECTED_POPULATION = [
    [
        [
            1.0, 0.996297793906866, 0.919866745769127, 0.843115367467209, 0.451096514349396,
            0.394564029095284, 0.41136354991753, 0.399871554840764, 0.828487788975188, 0.999952150419256
        ],
        [
            1.0, 0.975093002264432, 0.905887253741386, 0.789822408421996, 0.386248474568972,
            0.448706562744717, 0.372603068181163, 0.244813271324819, 0.798091225434512, 0.9999253127422
        ],
        [
            0.0, 0.027015236110399, 0.112856151311767, 0.205710496725285, 0.556425483232299,
            0.583922736355513, 0.673462875159617, 0.742200446549557, 0.173298315048093, 5.66376102744832e-05
        ],
        [
            0.0, 0.00159396784041352, 0.0613898492118107, 0.161351727426722, 0.6062295276719,
            0.572806671669023, 0.542570506643452, 0.613114726569719, 0.200122670317697, 6.58995172308447e-05
        ],
    ],
    [
        [
            1.0, 0.997890094581753, 0.992951177591028, 0.755451647968075, 0.446880660066254,
            0.678254305695527, 0.522849029774711, 0.438636085825306, 0.938434266724218, 0.999917758840887
        ],
        [
            0.0, 0.0258827263860563, 0.0451065159150301, 0.252443036918438, 0.548056965882489,
            0.327702300722201, 0.482409102411276, 0.5646807208006, 0.0778047090156592, 0.000100075452083218
        ],
        [
            1.0, 0.927541117251786, 0.812935811718192, 0.670324769366341, 0.459712981804729,
            0.351177236183876, 0.620992497194567, 0.402723943503705, 0.0190102940419282, 7.23423225620632e-05
        ],
        [
            0.0, 0.0486860616661003, 0.149006494704418, 0.321780545582464, 0.545349392018398,
            0.642866157133447, 0.373749370077677, 0.593959249806213, 0.964750730583563, 0.999909823772636
        ],
    ],
]

# Compare output to expected result for 10 points
NUM_SAMPLES = 10
SAMPLE_INDICES = [int(i * (EXPECTED_LENGTH - 1) / (NUM_SAMPLES - 1)) for i in range(NUM_SAMPLES)]


def test_example_piecewise_constant_controls(tmp_path, request):
    """Test CNOT gate optimization with piecewise constant controls using new Python interface."""
    datadir_path = os.path.join(tmp_path, request.node.name)

    freq01 = [4.80595, 4.8601]
    Jkl = [0.005]
    favg = sum(freq01)/len(freq01)
    rotfreq = favg*np.ones(len(freq01))
    T = 200.0

    unitary = np.identity(4)
    unitary[2, 2] = 0.0
    unitary[3, 3] = 0.0
    unitary[2, 3] = 1.0
    unitary[3, 2] = 1.0

    spline_order = 0
    nsplines = 1000
    gamma_variation = 1.0
    control_enforce_BC = True

    n_osc = 2
    n_levels = 4

    setup = setup_quandary(
        nessential=[2, 2],
        transition_frequency=freq01,
        final_time=T,
        dipole_coupling=Jkl,
        rotation_frequency=list(rotfreq),
        spline_order=spline_order,
        nspline=nsplines,
        control_zero_boundary_condition=control_enforce_BC,
        output_directory=datadir_path,
    )

    # Match old interface defaults
    setup.rand_seed = 1234
    setup.optim_tol_final_cost = 1e-4
    setup.optim_penalty_leakage = 0.1
    setup.optim_penalty_energy = 0.1
    setup.optim_penalty_dpdm = 0.01
    setup.output_optimization_stride = 1
    setup.linearsolver_maxiter = 20
    setup.optim_penalty_variation = gamma_variation

    # Match old initialization amplitude
    # With spline_order=0, carrier_frequencies are [[0.0], ...] so N_carriers=1
    num_carrier_freqs = len(setup.carrier_frequencies[0])
    init_amplitude = 10.0 / 1000.0 / np.sqrt(2) / num_carrier_freqs

    results = optimize(
        setup,
        targetgate=unitary,
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
