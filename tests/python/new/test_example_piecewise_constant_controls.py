import os
import pytest
import numpy as np
from quandary.new import create_config, optimize
from utils import assert_results_equal


# Mark all tests in this file as regression tests
pytestmark = pytest.mark.regression

EXPECTED_LENGTH = 1222
EXPECTED_INFIDELITY = 0.05503281532406801

EXPECTED_PT = [
    [
        -3.05520454819145, 7.1275586302858995, -2.8593579534932, -2.6353722407657396, 1.4783302843688102,
        0.15115726768055798, -0.49071058357207703, 6.65478489289128, 3.21749917073204, -3.57039637534921
    ],
    [
        -0.226419880828817, 1.93778097874674, -3.96960926020117, 4.16712590759117, -5.00228666422289,
        1.8946343277331599, 3.73373656244173, -4.20036548543164, 3.19345567048068, -0.0549681037826391
    ],
]

EXPECTED_QT = [
    [
        1.18665614617633, 5.64790168523117, -2.28586553136449, -0.058002377214981304, -3.40875353226556,
        2.5677053669047103, -2.37757220951372, 1.07262234264789, 2.48792697043541, -0.950748358024412
    ],
    [
        0.268560926841457, 1.52726026756229, -1.13869506369099, 0.8368563366082821, -0.0737437721088765,
        -3.60991869848623, 5.68093124784212, -3.4260279749075897, 2.00730218119136, -1.13684849732594
    ],
]

EXPECTED_ENERGY = [
    [
        [
            0.0, 0.0169544739630093, 0.130068513567157, 0.180525004337624, 0.46625933682294,
            0.520240055493922, 0.476350886984444, 0.516968889924946, 0.209710844111083, 0.011188117656886
        ],
        [
            0.0, 0.0300763330919059, 0.159396397311751, 0.264465210006719, 0.481364851338765,
            0.487782632542616, 0.578299892718885, 0.731101958176481, 0.273767180967537, 0.016874253653599
        ],
        [
            1.0, 0.965748558360874, 0.831459400619983, 0.785054184998802, 0.511461944705542,
            0.398761975292377, 0.340961187835901, 0.29859800358922, 0.771500781699245, 0.986505844897817
        ],
        [
            1.0, 0.987220634678893, 0.879075688744442, 0.769955600865248, 0.540913867591047,
            0.593215337222911, 0.604388033104609, 0.453331149341416, 0.745021194410056, 0.985431784018733
        ],
    ],
    [
        [
            0.0, 0.026874444211076, 0.155144671497274, 0.404460281916067, 0.416998823647702,
            0.266253274377566, 0.436600824910961, 0.494421856573668, 0.169357250680674, 0.0512551773973343
        ],
        [
            1.0, 0.952111059097199, 0.81469071925784, 0.601060422932871, 0.582388524786683,
            0.710166043916773, 0.545138082963778, 0.490099448208278, 0.8219990379843, 0.945480331044952
        ],
        [
            0.0, 0.0215374645927611, 0.0569732829486814, 0.168636956812742, 0.372846364306243,
            0.499894979014515, 0.47157954902154, 0.688774693413366, 0.934129633549857, 0.970024337071052
        ],
        [
            1.0, 0.999477031967644, 0.973191326172259, 0.825842338499298, 0.627766287508744,
            0.523685702898274, 0.546681543283012, 0.326704002317612, 0.0745140782391182, 0.0332401548170014
        ],
    ],
]

EXPECTED_POPULATION = [
    [
        [
            1.0, 0.983045526069858, 0.869931486472239, 0.819474995721435, 0.533740663238269,
            0.479759944571415, 0.523649113083584, 0.483031110146872, 0.790289155965344, 0.988811882427472
        ],
        [
            1.0, 0.96992366691134, 0.840603602725134, 0.735534790049687, 0.518635148720171,
            0.512217367521015, 0.421700107347784, 0.268898041895806, 0.72623281910872, 0.983125746425983
        ],
        [
            0.0, 0.0342514416442723, 0.168540599399332, 0.214945815042795, 0.488538055337846,
            0.601238024755056, 0.659038812214143, 0.701401996464143, 0.228499218357771, 0.0134941551665292
        ],
        [
            0.0, 0.0127793653221087, 0.120924311293583, 0.23004439919337, 0.459086132470242,
            0.40678466284274, 0.395611966963867, 0.54666885073082, 0.254978805665867, 0.014568216060655
        ],
    ],
    [
        [
            1.0, 0.973125555821791, 0.844855328542122, 0.595539718142993, 0.583001176413507,
            0.73374672568777, 0.563399175157067, 0.505578143498149, 0.830642749395753, 0.948744822687024
        ],
        [
            0.0, 0.0478889409060468, 0.185309280779045, 0.398939577123534, 0.417611475272253,
            0.289833956146858, 0.45486191710289, 0.509900551864008, 0.178000962091957, 0.0545196690346298
        ],
        [
            1.0, 0.978462535412385, 0.943026717070634, 0.831363043228855, 0.627153635737145,
            0.500105021032919, 0.528420451028504, 0.311225306639997, 0.0658703665071596, 0.0299756629932941
        ],
        [
            0.0, 0.000522968033357478, 0.0268086738657661, 0.174157661559321, 0.372233712552544,
            0.476314297167377, 0.453318456785464, 0.673295997754624, 0.925485921836805, 0.966759845262387
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

    setup = create_config(
        nessential=[2, 2],
        transition_frequency=freq01,
        total_time=T,
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
    setup.optim_maxiter = 20

    # Match old initialization amplitude
    # With spline_order=0, carrier_frequencies are [[0.0], ...] so N_carriers=1
    num_carrier_freqs = len(setup.carrier_frequencies[0])
    init_amplitude = 10.0 / 1000.0 / np.sqrt(2) / num_carrier_freqs

    results = optimize(
        setup,
        target=unitary,
        control_amplitude=init_amplitude,
        quiet=True,
    )

    assert_results_equal(
        t=results.time,
        p_samples=results.p_samples,
        q_samples=results.q_samples,
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
