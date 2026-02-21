import os
import pytest
import numpy as np
from quandary.new import setup_quandary, simulate
from utils import assert_results_equal

# Mark all tests in this file as regression tests
pytestmark = pytest.mark.regression

EXPECTED_LENGTH = 1001
EXPECTED_INFIDELITY = 1.0

EXPECTED_PT = np.zeros((8, 10))

EXPECTED_QT = np.zeros((8, 10))

EXPECTED_ENERGY = [
    [
        [
            1.0, 0.836078135681189, 0.677080386733117, 0.670628110768053, 0.92300799094264,
            0.837975552236461, 0.805823129197361, 0.853504682689832, 0.866549990398791, 0.806457068240172
        ],
    ],
    [
        [
            1.0, 0.636388342297137, 0.487092110002376, 0.54413722703352, 0.771657513370998,
            0.666408072335883, 0.612024569159776, 0.824676325668549, 0.419674820100051, 0.702460959050306
        ],
    ],
    [
        [
            1.0, 0.66748361234842, 0.46347452987769, 0.744184045857701, 0.802417472822974,
            0.770651675593258, 0.63234540776416, 0.383254206449323, 0.785911746247381, 0.708761215708732
        ],
    ],
    [
        [
            1.0, 0.66009154910948, 0.492257221761028, 0.619014865091762, 0.513571520786539,
            0.653392374584252, 0.629386909048527, 0.660832003911264, 0.682292516326034, 0.707144820937568
        ],
    ],
    [
        [
            0.0, 0.500108453041719, 0.534431835905109, 0.47752894868167, 0.242450496754744,
            0.293746749349346, 0.518698522063187, 0.637304278577435, 0.659394914028851, 0.310680732457375
        ],
    ],
    [
        [
            0.0, 0.269069449733598, 0.243034570485257, 0.253350482543049, 0.30230026530336,
            0.14570551680724, 0.267274075325126, 0.112004940016193, 0.274423918272602, 0.237851814417754
        ],
    ],
    [
        [
            0.0, 0.243117466829247, 0.555330150338723, 0.391204137042619, 0.197412587591839,
            0.179316513251891, 0.195590432087981, 0.371427834282674, 0.151967524907041, 0.118392848324015
        ],
    ],
    [
        [
            0.0, 0.187662990975368, 0.547299194929018, 0.299952183030103, 0.247182152491542,
            0.452803545922466, 0.338856955450839, 0.156995728517847, 0.159784569848526, 0.408250541009659
        ],
    ],
]

EXPECTED_POPULATION = [
    [
        [
            0.0, 0.16392186432285, 0.322919613274963, 0.329371889244066, 0.0769920090735191,
            0.162024447783738, 0.194176870826878, 0.146495317338447, 0.133450009633529, 0.193542931796223
        ],
    ],
    [
        [
            0.0, 0.363611657706903, 0.512907890005703, 0.455862772978599, 0.228342486645161,
            0.333591927684317, 0.387975430864463, 0.17532367435973, 0.580325179932268, 0.297539040986089
        ],
    ],
    [
        [
            0.0, 0.33251638765562, 0.536525470130389, 0.255815954154418, 0.197582527193185,
            0.229348324426941, 0.367654592260079, 0.616745793578956, 0.214088253784938, 0.291238784327663
        ],
    ],
    [
        [
            0.0, 0.339908450894559, 0.507742778247051, 0.380985134920358, 0.48642847922962,
            0.346607625435947, 0.370613090975712, 0.339167996117016, 0.317707483706285, 0.292855179098828
        ],
    ],
    [
        [
            1.0, 0.499891546962321, 0.465568164102971, 0.522471051330449, 0.757549503261415,
            0.706253250670854, 0.481301477961052, 0.362695721450844, 0.340605086003469, 0.689319267579021
        ],
    ],
    [
        [
            1.0, 0.730930550270442, 0.756965429522822, 0.74664951746907, 0.697699734712799,
            0.854294483212959, 0.732725924699113, 0.887995060012086, 0.725576081759717, 0.762148185618641
        ],
    ],
    [
        [
            1.0, 0.756882533174792, 0.444669849669357, 0.6087958629695, 0.80258741242432,
            0.820683486768308, 0.804409567936258, 0.628572165745605, 0.848032475125278, 0.88160715171238
        ],
    ],
    [
        [
            1.0, 0.812337009028671, 0.452700805079061, 0.700047816982017, 0.752817847524617,
            0.547196454097733, 0.6611430445734, 0.843004271510433, 0.840215430183794, 0.591749459026737
        ],
    ],
]

# Compare output to expected result for 10 points
NUM_SAMPLES = 10
SAMPLE_INDICES = [int(i * (EXPECTED_LENGTH - 1) / (NUM_SAMPLES - 1)) for i in range(NUM_SAMPLES)]


def mapCoeffs_SpinChainToQuandary(N: int, h: list, U: list, J: list):
    """
    Map spin chain coefficients J, U and h onto Quandary coefficients
    """
    # 01 transition frequencies [GHz] per site (omega_q)
    freq01 = np.zeros(N)
    for i in range(1, N-1):
        freq01[i] = (-2*h[i] - 2*U[i] - 2*U[i-1]) / (2*np.pi)
    freq01[0] = (-2*h[0] - 2*U[0]) / (2*np.pi)
    freq01[N-1] = (-2*h[N-1] - 2*U[N-2]) / (2*np.pi)

    # Jkl and Xi term [GHz]
    Jkl = []
    crosskerr = []
    couplingID = 0
    for i in range(N):
        for j in range(i+1, N):
            if j == i+1:  # linear chain coupling
                valJ = - 2*J[couplingID] / (2*np.pi)
                valC = - 4*U[couplingID] / (2*np.pi)
            else:
                valJ = 0.0
                valC = 0.0
            Jkl.append(valJ)
            crosskerr.append(valC)
        couplingID += 1
    return freq01, crosskerr, Jkl


def test_example_spinchain(tmp_path, request):
    """Test spin chain simulation using new Python interface."""
    datadir_path = os.path.join(tmp_path, request.node.name)

    N = 8
    U_amp = 1.0
    J_amp = 1.0
    np.random.seed(9001)  # Set seed for reproducibility
    h = np.random.uniform(-U_amp, U_amp, N)
    U = np.zeros(N)
    J = J_amp * np.ones(N)

    # Specify the initial state (domain wall |111...000>)
    initstate = np.zeros(N, dtype=int)
    for i in range(int(N/2)):
        initstate[i] = 1

    # Set the simulation duration and step size
    T = 10.0
    dT = 0.01

    freq01, crosskerr, Jkl = mapCoeffs_SpinChainToQuandary(N, h, U, J)

    n_osc = N
    n_levels = 1

    setup = setup_quandary(
        nessential=[2]*N,
        nguard=[0]*N,
        transition_frequency=list(freq01),
        rotation_frequency=list(np.zeros(N)),
        crosskerr_coupling=crosskerr,
        dipole_coupling=Jkl,
        carrier_frequency=[[0.0] for _ in range(N)],
        initial_levels=[int(x) for x in initstate],
        final_time=T,
        dt=dT,
        output_directory=datadir_path,
        verbose=False,
    )

    # Match old interface defaults (relevant for simulation)
    setup.output_optimization_stride = 1
    setup.linearsolver_maxiter = 20

    # Run forward simulation with explicit zero controls.
    # Passing pt0/qt0 ensures control_initializations is set so the C++ outputs
    # control*.dat files (needed by assert_results_equal to check len(pt) == n_osc).
    ntime = setup.ntime
    pt0 = [np.zeros(ntime + 1) for _ in range(N)]
    qt0 = [np.zeros(ntime + 1) for _ in range(N)]
    results = simulate(setup, pt0=pt0, qt0=qt0, quiet=True)

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
