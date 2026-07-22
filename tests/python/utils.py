import numpy as np
from pytest import approx

REL_TOL = 1.0e-3
ABS_TOL = 1.0e-10


def assert_results_equal(
        t, p_samples, q_samples, infidelity, energy, population, T, n_osc, n_levels, sample_indices,
        expected_length, expected_infidelity, expected_pt, expected_qt, expected_energy, expected_population):
    """
    Utility function to assert that the results of quantum run match the expected values within given tolerances.
    """
    assert t[0] == 0.0 and t[-1] == T
    assert len(t) == expected_length
    assert infidelity == approx(expected_infidelity, rel=REL_TOL, abs=ABS_TOL)

    assert len(p_samples) == n_osc
    assert len(q_samples) == n_osc

    for i in range(n_osc):
        pt_samples = [p_samples[i][idx] for idx in sample_indices]
        qt_samples = [q_samples[i][idx] for idx in sample_indices]
        np.testing.assert_allclose(pt_samples, expected_pt[i], rtol=REL_TOL, atol=ABS_TOL)
        np.testing.assert_allclose(qt_samples, expected_qt[i], rtol=REL_TOL, atol=ABS_TOL)

    assert len(energy) == n_osc
    assert len(population) == n_osc
    assert len(energy[0]) == n_levels
    assert len(population[0]) == n_levels

    for i in range(n_osc):
        for j in range(n_levels):
            energy_data = energy[i][j]
            energy_samples = [energy_data[idx] for idx in sample_indices]
            np.testing.assert_allclose(energy_samples, expected_energy[i][j], rtol=REL_TOL, atol=ABS_TOL)

            # Note: only comparing population for the first level (index 0)
            pop_data = population[i][j]
            pop_samples = [pop_data[0, idx] for idx in sample_indices]
            np.testing.assert_allclose(pop_samples, expected_population[i][j], rtol=REL_TOL, atol=ABS_TOL)


def print_expected_values(infidelity, p_samples, q_samples, energy, population, sample_indices, n_osc):
    """
    Utility function to print actual values in the format needed for EXPECTED arrays.
    Call this function with actual test results to get copy-pasteable expected values.
    """
    print()
    print(f"EXPECTED_LENGTH = {len(p_samples[0])}")
    print(f"EXPECTED_INFIDELITY = {infidelity}")
    print()

    print("EXPECTED_PT = [")
    for i in range(n_osc):
        pt_samples = [p_samples[i][idx] for idx in sample_indices]
        print("    [")
        # Format with max 5 values per line for readability
        for j in range(0, len(pt_samples), 5):
            chunk = pt_samples[j:j+5]
            if j + 5 >= len(pt_samples):
                print(f"        {', '.join(map(str, chunk))}")
            else:
                print(f"        {', '.join(map(str, chunk))},")
        print("    ],")
    print("]")
    print()

    print("EXPECTED_QT = [")
    for i in range(n_osc):
        qt_samples = [q_samples[i][idx] for idx in sample_indices]
        print("    [")
        for j in range(0, len(qt_samples), 5):
            chunk = qt_samples[j:j+5]
            if j + 5 >= len(qt_samples):
                print(f"        {', '.join(map(str, chunk))}")
            else:
                print(f"        {', '.join(map(str, chunk))},")
        print("    ],")
    print("]")
    print()

    print("EXPECTED_ENERGY = [")
    for i in range(n_osc):
        print("    [")
        for j in range(len(energy[i])):
            energy_data = energy[i][j]
            energy_samples = [energy_data[idx] for idx in sample_indices]
            print("        [")
            for k in range(0, len(energy_samples), 5):
                chunk = energy_samples[k:k+5]
                if k + 5 >= len(energy_samples):
                    print(f"            {', '.join(map(str, chunk))}")
                else:
                    print(f"            {', '.join(map(str, chunk))},")
            print("        ],")
        print("    ],")
    print("]")
    print()

    print("EXPECTED_POPULATION = [")
    for i in range(n_osc):
        print("    [")
        for j in range(len(population[i])):
            pop_data = population[i][j]
            pop_samples = [pop_data[0, idx] for idx in sample_indices]
            print("        [")
            for k in range(0, len(pop_samples), 5):
                chunk = pop_samples[k:k+5]
                if k + 5 >= len(pop_samples):
                    print(f"            {', '.join(map(str, chunk))}")
                else:
                    print(f"            {', '.join(map(str, chunk))},")
            print("        ],")
        print("    ],")
    print("]")
