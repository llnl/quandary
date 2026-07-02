"""Tests for internal machinery: dry_run, core distribution, and subprocess execution."""
import numpy as np
from quandary.new import ConfigInput, create_config, optimize, simulate, evaluate_controls
from quandary.new.run import (
    _compute_optimal_core_distribution,
)

import quandary.new.run as patch_run


def _make_setup(output_directory="./run_dir"):
    """Create a minimal setup for testing."""
    return create_config(
        nessential=[2],
        transition_frequency=[4.0],
        selfkerr=[0.2],
        total_time=1.0,
        dt=0.1,
        spline_order=0,
        output_directory=output_directory,
    )


def test_config_input_copy_preserves_wrapper_type():
    """ConfigInput.copy() keeps Python wrapper behavior while copying values."""
    setup = ConfigInput()
    setup.nlevels = np.array([2, 3])
    setup.dt = 0.1

    copied = setup.copy()

    assert type(copied) is ConfigInput
    assert copied is not setup
    assert copied.nlevels == [2, 3]
    assert copied.dt == 0.1
    assert repr(copied).startswith("ConfigInput(")

    setup.nlevels = [4, 5]
    assert copied.nlevels == [2, 3]

    copied.nlevels = np.array([6, 7])
    assert copied.nlevels == [6, 7]


class TestDryRun:
    """Validate that dry_run=True returns Results with config but no solver output."""

    def test_optimize_dry_run_has_config(self, tmp_path):
        """optimize(dry_run=True) returns Results with a populated config."""
        setup = _make_setup(str(tmp_path / "opt"))
        results = optimize(
            setup,
            target=[0.0, 1.0],
            dry_run=True,
            quiet=True,
        )
        assert results.config is not None
        assert results.config.output_directory == str(tmp_path / "opt")
        # No solver output
        assert len(results.time) == 0
        assert len(results.p_samples) == 0

    def test_simulate_dry_run_has_config(self, tmp_path):
        """simulate(dry_run=True) returns Results with a populated config."""
        setup = _make_setup(str(tmp_path / "sim"))
        results = simulate(
            setup,
            dry_run=True,
            quiet=True,
        )
        assert results.config is not None
        assert results.config.output_directory == str(tmp_path / "sim")
        assert len(results.time) == 0

    def test_evaluate_controls_dry_run_has_config(self, tmp_path):
        """evaluate_controls(dry_run=True) returns Results with a populated config."""
        setup = _make_setup(str(tmp_path / "eval"))
        spline_coefficients = np.zeros(10)
        results = evaluate_controls(
            setup,
            spline_coefficients=spline_coefficients,
            dry_run=True,
            quiet=True,
        )
        assert results.config is not None
        assert len(results.time) == 0

    def test_dry_run_does_not_mutate_setup(self, tmp_path):
        """The original setup object is not modified by a dry_run call."""
        setup = _make_setup(str(tmp_path / "orig"))
        original_dt = setup.dt
        original_output_dir = setup.output_directory

        optimize(
            setup,
            target=[0.0, 1.0],
            dry_run=True,
            quiet=True,
        )

        assert setup.dt == original_dt
        assert setup.output_directory == original_output_dir
        # Runtype should not have been set on the original
        assert setup.runtype is None

    def test_dry_run_output_directory_flows_through(self, tmp_path):
        """Output directory set in create_config flows to the dry_run config."""
        custom_dir = str(tmp_path / "custom_output")
        setup = _make_setup(custom_dir)
        results = optimize(
            setup,
            target=[0.0, 1.0],
            dry_run=True,
            quiet=True,
        )
        assert results.config.output_directory == custom_dir


class TestInteractiveSubprocess:
    """Compare an interactive run via a subprocess to a direct run through python."""

    def test_optimize_uses_subprocess_when_interactive(self, monkeypatch, tmp_path):
        """When _is_interactive() is true, optimize() spawns the subprocess path."""

        # setup = _make_setup(str(tmp_path / "interactive"))
        dir = "./tmp_interactive"
        setup = _make_setup(dir)

        # Run directly through python
        result_direct = optimize(setup, target=[0.0, 1.0], quiet=True)

        # Run through subprocess by monkeypatching _is_interactive to return True
        monkeypatch.setattr(patch_run, "_is_interactive", lambda: True)
        result_subprocess = optimize(setup, target=[0.0, 1.0], quiet=True)

        assert np.allclose(result_direct.time, result_subprocess.time)
        assert np.allclose(result_direct.p_samples, result_subprocess.p_samples)
        assert np.allclose(result_direct.q_samples, result_subprocess.q_samples)
        assert np.isclose(result_direct.infidelity, result_subprocess.infidelity)
        assert np.allclose(result_direct.spline_coefficients, result_subprocess.spline_coefficients)
        assert np.allclose(result_direct.optim_hist["gradient"], result_subprocess.optim_hist["gradient"])


class TestCoreDistribution:
    """Unit tests for _compute_optimal_core_distribution."""

    def test_exact_divisor(self):
        """maxcores is an exact divisor of ninit."""
        assert _compute_optimal_core_distribution(4, 8) == 4

    def test_rounds_down_to_divisor(self):
        """maxcores is not a divisor -- rounds down to nearest one."""
        assert _compute_optimal_core_distribution(3, 8) == 2

    def test_maxcores_greater_than_ninit(self):
        """maxcores > ninit -- capped at ninit."""
        assert _compute_optimal_core_distribution(16, 4) == 4

    def test_single_core(self):
        """Single core always works."""
        assert _compute_optimal_core_distribution(1, 7) == 1

    def test_prime_ninit(self):
        """Prime ninit with maxcores < ninit -- can only use 1 core."""
        assert _compute_optimal_core_distribution(4, 7) == 1

    def test_prime_ninit_exact(self):
        """maxcores == prime ninit -- use all."""
        assert _compute_optimal_core_distribution(7, 7) == 7

    def test_unlimited_cores(self):
        """maxcores == -1 means unlimited -- use all ninit."""
        assert _compute_optimal_core_distribution(-1, 6) == 6

    def test_ninit_one(self):
        """Only one initial condition -- always 1 core."""
        assert _compute_optimal_core_distribution(8, 1) == 1
