"""Tests for get_results function."""

import numpy as np
import os
import pytest
from quandary.new import create_config, get_results, Results, simulate


def test_get_results(tmp_path, request):
    datadir_path = os.path.join(tmp_path, request.node.name)

    Ne = [2,2]
    Ng = [0,0]
    freq01 = [5.18, 5.12] 
    selfkerr = [0.21, 0.21]
    favg = sum(freq01)/len(freq01) 
    rotfreq = favg*np.ones(len(freq01))
    Jkl = [0.01]
    T = 10.0
    Pmin = 20   
    maxctrl_MHz = 25.0           
    bound = [maxctrl_MHz*1e-3 for _ in range(len(Ne))] 
    spline_knot_spacing =  3.0   # Bspline spacing [ns]
    control_boundary_zero = False 
    initctrl_amp = 10.0*1e-3  # Scale to GHz
    initctrl_randomize = True

    def get_QFT_gate(dim):
        gate_Hd =  np.zeros((dim, dim), dtype=complex)
        om_d = np.exp(1j*2*np.pi/dim)
        for j in range(dim):
            for k in range(dim):
                gate_Hd[j,k] = om_d**(j*k) / np.sqrt(dim)
        return gate_Hd

    target = get_QFT_gate(np.prod(Ne)) # Gate optimization

    cfg = create_config(
        total_time = T,
        nessential = Ne,
        nguard = Ng,
        transition_frequency = freq01,
        rotation_frequency = rotfreq,
        selfkerr = selfkerr,
        dipole_coupling = Jkl,
        Pmin = Pmin,
        spline_knot_spacing = spline_knot_spacing,
        control_amplitude_bound = bound,
        control_zero_boundary_condition=control_boundary_zero,
        target = target,
        control_amplitude=initctrl_amp, 
        control_randomize=initctrl_randomize,
        output_directory = datadir_path,
    )
    result_original = simulate(cfg, quiet=True)

    # Load the results from the output directory 
    datadir = cfg.output_directory
    result_test = get_results(datadir)
    
    """ Test that get_results returns a Results and config."""
    assert isinstance(result_test, Results)
    assert result_test.config is not None
    assert np.allclose(result_test.config.nlevels, [Ne[i] + Ng[i] for i in range(len(Ne))], atol=1e-13), f"Expected nlevels ~ {[Ne[i] + Ng[i] for i in range(len(Ne))]}, got {result_test.config.nlevels}"
    assert np.allclose(result_test.config.nessential, Ne, atol=1e-13), f"Expected nessential ~ {Ne}, got {result_test.config.nessential}"
    assert np.allclose(result_test.config.transition_frequency, freq01, atol=1e-13), f"Expected transition_frequency ~ {freq01}, got {result_test.config.transition_frequency}"
    assert np.allclose(result_test.config.selfkerr, selfkerr, atol=1e-13), f"Expected selfkerr ~ {selfkerr}, got {result_test.config.selfkerr}"
    assert np.allclose(result_test.config.rotation_frequency, rotfreq, atol=1e-13), f"Expected rotation_frequency ~ {rotfreq}, got {result_test.config.rotation_frequency}"
    assert np.allclose(result_test.config.dipole_coupling, Jkl, atol=1e-13), f"Expected dipole_coupling ~ {Jkl}, got {result_test.config.dipole_coupling}"
    assert result_test.config.total_time == T, f"Expected total_time = {T}, got {result_test.config.total_time}"

    """Test that control pulses are loaded with correct dimensions"""
    assert len(result_test.p_samples) > 0
    assert len(result_test.q_samples) > 0
    assert isinstance(result_test.p_samples[0], np.ndarray)
    assert isinstance(result_test.q_samples[0], np.ndarray)
    assert(len(result_test.p_samples) == len(result_original.p_samples))
    assert(len(result_test.q_samples) == len(result_original.q_samples))

    """Test that control pulses are the same in the original result and the loaded result."""
    for i in range(len(result_original.p_samples)):
        assert np.allclose(result_test.p_samples[i], result_original.p_samples[i], atol=1e-13), f"p_samples mismatch for control {i}"
        assert np.allclose(result_test.q_samples[i], result_original.q_samples[i], atol=1e-13), f"q_samples mismatch for control {i}"

    """Test that time vector aligns with pulse samples and matches the original."""
    assert len(result_test.time) > 0
    assert len(result_test.time) == len(result_test.p_samples[0])
    assert len(result_test.time) == len(result_test.q_samples[0])
    assert np.allclose(result_test.time, result_original.time, atol=1e-13), "Time vector mismatch"

    """Test that spline coefficients are loaded and match the original."""
    assert len(result_test.spline_coefficients) > 0
    assert isinstance(result_test.spline_coefficients, np.ndarray)
    assert len(result_test.spline_coefficients) == len(result_original.spline_coefficients)
    assert np.allclose(result_test.spline_coefficients, result_original.spline_coefficients, atol=1e-13), f"Spline coefficients mismatch"

    """Test that optimization history is loaded and matches the original."""
    assert len(result_test.optim_hist) > 0
    assert 'objective' in result_test.optim_hist
    assert 'fidelity' in result_test.optim_hist
    assert len(result_test.optim_hist['objective']) > 0
    assert len(result_test.optim_hist['fidelity']) > 0
    assert np.allclose(result_test.optim_hist['objective'], result_original.optim_hist['objective'], atol=1e-13), "Optimization history objective mismatch"
    assert np.allclose(result_test.optim_hist['fidelity'], result_original.optim_hist['fidelity'], atol=1e-13), "Optimization history fidelity mismatch"
    assert np.allclose(result_test.optim_hist['gradient'], result_original.optim_hist['gradient'], atol=1e-13), "Optimization history gradient mismatch"

    """Test that population data is loaded and matches the original."""
    assert len(result_test.population) > 0
    assert len(result_test.population[0]) > 0  # At least one initial condition
    assert isinstance(result_test.population[0][0], np.ndarray)
    for i in range(len(result_original.population)):
        for j in range(len(result_original.population[i])):
            assert np.allclose(result_test.population[i][j], result_original.population[i][j], atol=1e-13), f"Population data mismatch for oscillator {i}, initial condition {j}"

    """Test that expected energy data is loaded and matches original."""
    assert len(result_test.expected_energy) > 0
    assert len(result_test.expected_energy[0]) > 0  # At least one initial condition
    assert isinstance(result_test.expected_energy[0][0], np.ndarray)
    for i in range(len(result_test.expected_energy)):
        for j in range(len(result_test.expected_energy[i])):
            assert np.allclose(result_test.expected_energy[i][j], result_original.expected_energy[i][j], atol=1e-13), f"Expected energy data mismatch for oscillator {i}, initial condition {j}"

