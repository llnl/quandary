"""Test custom Hamiltonian support in setup_quandary.

Builds a standard 2-level qubit Hamiltonian manually as numpy arrays, passes them
via hamiltonian_Hsys/hamiltonian_Hc to setup_quandary, and verifies the simulation
produces the same result as the standard model path.
"""

import os
import pytest
from quandary.new import setup_quandary, simulate, hamiltonians


# System parameters for a single 2-level qubit
NESSENTIAL = [2]
TRANSITION_FREQUENCY = [4.10595]
SELFKERR = [0.2198]
ROTATION_FREQUENCY = [4.10595]
FINAL_TIME = 10.0
NTIME = 100


@pytest.fixture
def standard_setup(tmp_path):
    """Create a setup using the standard model (no custom Hamiltonian)."""
    return setup_quandary(
        nessential=NESSENTIAL,
        transition_frequency=TRANSITION_FREQUENCY,
        selfkerr=SELFKERR,
        rotation_frequency=ROTATION_FREQUENCY,
        final_time=FINAL_TIME,
        ntime=NTIME,
        output_directory=os.path.join(tmp_path, "std"),
    )


@pytest.fixture
def standard_hamiltonians():
    """Build the standard Hamiltonian matrices as numpy arrays."""
    nlevels = NESSENTIAL  # no guard levels
    Hsys, Hc_re, Hc_im = hamiltonians(
        N=nlevels,
        transition_frequency=TRANSITION_FREQUENCY,
        selfkerr=SELFKERR,
        rotation_frequency=ROTATION_FREQUENCY,
    )
    Hc = [Hc_re[i] + 1j * Hc_im[i] for i in range(len(Hc_re))]
    return Hsys, Hc


def test_custom_hsys_and_hc(tmp_path, standard_setup, standard_hamiltonians):
    """Custom Hsys + Hc should produce the same result as the standard model."""
    Hsys, Hc = standard_hamiltonians

    setup = setup_quandary(
        nessential=NESSENTIAL,
        transition_frequency=TRANSITION_FREQUENCY,
        selfkerr=SELFKERR,
        rotation_frequency=ROTATION_FREQUENCY,
        final_time=FINAL_TIME,
        ntime=NTIME,
        hamiltonian_Hsys=Hsys,
        hamiltonian_Hc=Hc,
        output_directory=os.path.join(tmp_path, "custom_both"),
    )

    assert setup.hamiltonian_file_Hsys is not None
    assert setup.hamiltonian_file_Hc is not None

    results_std = simulate(standard_setup, quiet=True)
    results_custom = simulate(setup, quiet=True)

    assert abs(results_std.infidelity - results_custom.infidelity) < 1e-6


def test_custom_hsys_only(tmp_path, standard_setup, standard_hamiltonians):
    """Custom Hsys without custom Hc should also work (standard a+aT operators used for Hc)."""
    Hsys, _ = standard_hamiltonians

    setup = setup_quandary(
        nessential=NESSENTIAL,
        transition_frequency=TRANSITION_FREQUENCY,
        selfkerr=SELFKERR,
        rotation_frequency=ROTATION_FREQUENCY,
        final_time=FINAL_TIME,
        ntime=NTIME,
        hamiltonian_Hsys=Hsys,
        output_directory=os.path.join(tmp_path, "custom_hsys"),
    )

    assert setup.hamiltonian_file_Hsys is not None
    assert setup.hamiltonian_file_Hc is None

    results_std = simulate(standard_setup, quiet=True)
    results_custom = simulate(setup, quiet=True)

    assert abs(results_std.infidelity - results_custom.infidelity) < 1e-6


def test_custom_hc_only(tmp_path, standard_setup, standard_hamiltonians):
    """Custom Hc without custom Hsys should work (standard Hsys model used)."""
    _, Hc = standard_hamiltonians

    setup = setup_quandary(
        nessential=NESSENTIAL,
        transition_frequency=TRANSITION_FREQUENCY,
        selfkerr=SELFKERR,
        rotation_frequency=ROTATION_FREQUENCY,
        final_time=FINAL_TIME,
        ntime=NTIME,
        hamiltonian_Hc=Hc,
        output_directory=os.path.join(tmp_path, "custom_hc"),
    )

    assert setup.hamiltonian_file_Hsys is None
    assert setup.hamiltonian_file_Hc is not None

    results_std = simulate(standard_setup, quiet=True)
    results_custom = simulate(setup, quiet=True)

    assert abs(results_std.infidelity - results_custom.infidelity) < 1e-6


def test_hamiltonian_files_written(tmp_path, standard_hamiltonians):
    """Verify that Hamiltonian files are actually written to disk."""
    Hsys, Hc = standard_hamiltonians
    outdir = os.path.join(tmp_path, "file_check")

    setup_quandary(
        nessential=NESSENTIAL,
        transition_frequency=TRANSITION_FREQUENCY,
        selfkerr=SELFKERR,
        rotation_frequency=ROTATION_FREQUENCY,
        final_time=FINAL_TIME,
        ntime=NTIME,
        hamiltonian_Hsys=Hsys,
        hamiltonian_Hc=Hc,
        output_directory=outdir,
    )

    assert os.path.isfile(os.path.join(outdir, "hamiltonian_Hsys.dat"))
    assert os.path.isfile(os.path.join(outdir, "hamiltonian_Hc.dat"))

    # Verify Hsys file format: header + sparse entries
    with open(os.path.join(outdir, "hamiltonian_Hsys.dat")) as f:
        lines = f.readlines()
    assert lines[0].startswith("#")
    # Each data line should have 4 fields: row col real imag
    for line in lines[1:]:
        parts = line.strip().split()
        assert len(parts) == 4

    # Verify Hc file format: header + sparse entries with oscillator index
    with open(os.path.join(outdir, "hamiltonian_Hc.dat")) as f:
        lines = f.readlines()
    assert lines[0].startswith("#")
    for line in lines[1:]:
        if line.startswith("#"):
            continue
        parts = line.strip().split()
        assert len(parts) == 5  # oscillator row col real imag
