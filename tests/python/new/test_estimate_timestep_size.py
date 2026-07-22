"""Tests for estimate_timestep_size function."""

import numpy as np
import pytest
from quandary.new import estimate_timestep_size, hamiltonians

# System parameters for a two 3-level qubit
NESSENTIAL = [2,2]
NGUARD = [1,1]
TRANSITION_FREQUENCY = [5.18, 5.12] 
SELFKERR = [0.22, 0.23]
ROTFREQ = 5.15*np.ones(len(TRANSITION_FREQUENCY))
JKL = [0.01]# Expected output 

# Expected timestep size 
EXPECTED_DT = 0.03535913101225683 

def test_estimate_timestep_size_2qubit():
    """Test timestep estimation for two-qubit system."""

    Hsys, Hc_re, Hc_im = hamiltonians(
        nessential=NESSENTIAL,
        nguard=NGUARD,
        transition_frequency=TRANSITION_FREQUENCY,
        selfkerr=SELFKERR,
        rotation_frequency=ROTFREQ,
        dipole_coupling=JKL,
    )
    
    dt = estimate_timestep_size(
        Hsys=Hsys,
        Hc_re=Hc_re,
        Hc_im=Hc_im,
        control_amplitude_bound=[0.1, 0.1],
        Pmin=40,
    )
    
    assert isinstance(dt, float)
    assert dt > 0
    # assert that dt matches to expected dt up to machine precision
    assert np.isclose(dt, EXPECTED_DT, atol=1e-13), f"Expected dt ~ {EXPECTED_DT}, got {dt}"


def test_estimate_timestep_size_pmin_effect():
    """Test that smaller Pmin gives larger timestep."""
    Hsys, Hc_re, Hc_im = hamiltonians(
        nessential=NESSENTIAL,
        nguard=NGUARD,
        transition_frequency=TRANSITION_FREQUENCY,
        selfkerr=SELFKERR,
        rotation_frequency=ROTFREQ,
        dipole_coupling=JKL,
    )
    
    dt_pmin40 = estimate_timestep_size(
        Hsys=Hsys,
        Hc_re=Hc_re,
        Hc_im=Hc_im,
        control_amplitude_bound=[0.1, 0.1],
        Pmin=40,
    )
    
    dt_pmin20 = estimate_timestep_size(
        Hsys=Hsys,
        Hc_re=Hc_re,
        Hc_im=Hc_im,
        control_amplitude_bound=[0.1, 0.1],
        Pmin=20,
    )
    
    assert dt_pmin20 > dt_pmin40, "Smaller Pmin should give larger timestep"

