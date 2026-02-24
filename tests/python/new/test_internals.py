"""Tests for internal machinery: dry_run, core distribution, and subprocess execution."""
import numpy as np
from quandary.new import setup_quandary, optimize, simulate, evaluate_controls
from quandary.new.runner import (
    _compute_optimal_core_distribution,
)


def _make_setup(output_directory="./run_dir"):
    """Create a minimal setup for testing."""
    return setup_quandary(
        nessential=[2],
        transition_frequency=[4.0],
        selfkerr=[0.2],
        final_time=1.0,
        ntime=10,
        spline_order=0,
        output_directory=output_directory,
        verbose=False,
    )


class TestDryRun:
    """Validate that dry_run=True returns Results with config but no solver output."""

    def test_optimize_dry_run_has_config(self, tmp_path):
        """optimize(dry_run=True) returns Results with a populated config."""
        setup = _make_setup(str(tmp_path / "opt"))
        results = optimize(
            setup,
            targetstate=[0.0, 1.0],
            dry_run=True,
            quiet=True,
        )
        assert results.config is not None
        assert results.config.output_directory == str(tmp_path / "opt")
        # No solver output
        assert len(results.time) == 0
        assert len(results.pt) == 0

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
        pcof = np.zeros(10)
        results = evaluate_controls(
            setup,
            pcof=pcof,
            dry_run=True,
            quiet=True,
        )
        assert results.config is not None
        assert len(results.time) == 0

    def test_dry_run_does_not_mutate_setup(self, tmp_path):
        """The original setup object is not modified by a dry_run call."""
        setup = _make_setup(str(tmp_path / "orig"))
        original_ntime = setup.ntime
        original_output_dir = setup.output_directory

        optimize(
            setup,
            targetstate=[0.0, 1.0],
            dry_run=True,
            quiet=True,
        )

        assert setup.ntime == original_ntime
        assert setup.output_directory == original_output_dir
        # Runtype should not have been set on the original
        assert setup.runtype is None

    def test_dry_run_output_directory_flows_through(self, tmp_path):
        """Output directory set in setup_quandary flows to the dry_run config."""
        custom_dir = str(tmp_path / "custom_output")
        setup = _make_setup(custom_dir)
        results = optimize(
            setup,
            targetstate=[0.0, 1.0],
            dry_run=True,
            quiet=True,
        )
        assert results.config.output_directory == custom_dir


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
