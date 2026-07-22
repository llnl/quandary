"""Tests for load_config_input function."""

import os
import pytest
import numpy as np
from quandary.new import load_config_input, ConfigInput


def _create_sample_toml(tmp_path):
    """Create a minimal TOML config file for testing."""
    toml_content = """
[system]
nlevels = [2, 2]
nessential = [2, 2]
total_time = 10.0
ntime = 17
dt = 0.588235294117647
transition_frequency = [5.18, 5.12]
selfkerr = [0.21, 0.21]
crosskerr_coupling = 0.0
dipole_coupling = 0.01
rotation_frequency = [5.15, 5.15]
decoherence = {
  type = "none",
  decay_time = [0.0, 0.0],
  dephase_time = [0.0, 0.0]
}
initial_condition = {type = "diagonal", subsystem = [0, 1]}
[control]
parameterization = {type = "spline", num = 6}
carrier_frequency = {value = [-0.0316227766016841, 0.0316227766016832]}
initialization = {type = "constant", amplitude = 0.01}
amplitude_bound = 0.025
zero_boundary_condition = false
[optimization]
target = {type = "gate", gate_type = "file", filename = "data_newinterface/target_gate.dat", gate_rot_freq = [0.0, 0.0]}
objective = "jtrace"
tolerance = { grad_abs = 0.0001, grad_rel = 0.0001, final_cost = 1e-08, infidelity = 1e-05 }
maxiter = 10
tikhonov = { coeff = 0.0001, use_x0 = false }
penalty = { leakage = 0.1, energy = 0.1, dpdm = 0.01, variation = 0, weightedcost = 0, weightedcost_width = 0 }
[output]
directory = "data_newinterface_sim"
observables = ["population", "expectedenergy", "fullstate"]
timestep_stride = 1
optimization_stride = 1
[solver]
runtype = "simulation"
usematfree = true
linearsolver = { type = "gmres", maxiter = 10 }
timestepper = "imr"
rand_seed = 1
"""
    
    config_file = os.path.join(tmp_path, "test_config.toml")
    os.makedirs(tmp_path, exist_ok=True)
    with open(config_file, "w") as f:
        f.write(toml_content)
    return str(config_file)

def test_load_config_input_system_parameters(tmp_path, request):
    """Test that config inputs are loaded correctly."""

    datadir_path = os.path.join(tmp_path, request.node.name)

    config_file = _create_sample_toml(datadir_path)
    config = load_config_input(config_file, quiet=True)
    
    assert isinstance(config, ConfigInput)
    assert np.allclose(config.nlevels, [2, 2])
    assert np.allclose(config.nessential, [2, 2])
    assert np.isclose(config.dt, 0.588235294117647)
    assert np.allclose(config.transition_frequency, [5.18, 5.12])
    assert np.allclose(config.selfkerr, [0.21, 0.21])
    assert np.allclose(config.rotation_frequency, [5.15, 5.15])
    assert np.allclose(config.carrier_frequencies, [[-0.0316227766016841, 0.0316227766016832]])

