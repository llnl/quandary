"""Tests that wrong-sized vector settings produce user-friendly ValidationError messages."""
import pytest
from quandary.new import setup_quandary, Config
from quandary._quandary_impl import ValidationError


def _make_setup(**overrides):
    """Create a minimal 1-qubit setup, applying overrides."""
    defaults = dict(
        nessential=[2],
        transition_frequency=[4.0],
        selfkerr=[0.2],
        final_time=1.0,
        ntime=10,
        spline_order=0,
        verbose=False,
    )
    defaults.update(overrides)
    return setup_quandary(**defaults)


def _make_2q_setup(**overrides):
    """Create a minimal 2-qubit setup, applying overrides."""
    defaults = dict(
        nessential=[2, 2],
        transition_frequency=[4.0, 4.5],
        selfkerr=[0.2, 0.2],
        final_time=1.0,
        ntime=10,
        spline_order=0,
        verbose=False,
    )
    defaults.update(overrides)
    return setup_quandary(**defaults)


class TestVectorSizeValidation:
    """Wrong-sized vectors should raise ValidationError with a clear message."""

    def test_selfkerr_wrong_size(self):
        setup = _make_setup()
        setup.selfkerr = [0.2, 0.3]  # 2 values for 1 qubit
        with pytest.raises(ValidationError, match="selfkerr.*must have exactly 1 elements.*got 2"):
            Config(setup, True)

    def test_rotation_frequency_wrong_size(self):
        setup = _make_setup()
        setup.rotation_frequency = [4.0, 5.0]  # 2 values for 1 qubit
        with pytest.raises(ValidationError, match="rotation_frequency.*must have exactly 1 elements.*got 2"):
            Config(setup, True)

    def test_decay_time_wrong_size(self):
        setup = _make_setup()
        setup.decay_time = [100.0, 200.0]  # 2 values for 1 qubit
        with pytest.raises(ValidationError, match="decay_time.*must have exactly 1 elements.*got 2"):
            Config(setup, True)

    def test_dephase_time_wrong_size(self):
        setup = _make_setup()
        setup.dephase_time = [100.0, 200.0]  # 2 values for 1 qubit
        with pytest.raises(ValidationError, match="dephase_time.*must have exactly 1 elements.*got 2"):
            Config(setup, True)

    def test_control_amplitude_bounds_wrong_size(self):
        setup = _make_setup(control_amplitude_bounds=[0.01, 0.02])  # 2 values for 1 qubit
        with pytest.raises(ValidationError, match="control_amplitude_bounds.*must have exactly 1 elements.*got 2"):
            Config(setup, True)

    def test_carrier_frequencies_wrong_size(self):
        setup = _make_setup()
        setup.carrier_frequencies = [[1.0], [2.0]]  # 2 oscillators for 1 qubit
        with pytest.raises(ValidationError, match="carrier_frequencies.*must have exactly 1 elements.*got 2"):
            Config(setup, True)

    def test_transition_frequency_wrong_size(self):
        setup = _make_setup()
        setup.transition_frequency = [4.0, 5.0]  # 2 values for 1 qubit
        with pytest.raises(ValidationError, match="transition_frequency.*must have exactly 1 elements.*got 2"):
            Config(setup, True)

    def test_crosskerr_coupling_wrong_size(self):
        setup = _make_2q_setup()
        setup.crosskerr_coupling = [0.01, 0.02]  # expects 1 pair
        with pytest.raises(ValidationError, match="crosskerr_coupling.*must have exactly 1 elements.*got 2"):
            Config(setup, True)

    def test_dipole_coupling_wrong_size(self):
        setup = _make_2q_setup()
        setup.dipole_coupling = [0.01, 0.02, 0.03]  # expects 1 pair
        with pytest.raises(ValidationError, match="dipole_coupling.*must have exactly 1 elements.*got 3"):
            Config(setup, True)

    def test_nessential_wrong_size(self):
        setup = _make_setup()
        setup.nessential = [2, 3]  # 2 values for 1 qubit (nlevels has 1 entry)
        with pytest.raises(ValidationError, match="nessential.*must have exactly 1 elements.*got 2"):
            Config(setup, True)


class TestCorrectSizeAccepted:
    """Correctly-sized vectors should not raise."""

    def test_1q_setup_validates(self):
        setup = _make_setup()
        config = Config(setup, True)
        assert config is not None

    def test_2q_setup_validates(self):
        setup = _make_2q_setup()
        config = Config(setup, True)
        assert config is not None

    def test_omitted_optionals_use_defaults(self):
        """When optional vector fields are omitted, defaults are used without error."""
        setup = setup_quandary(
            nessential=[2],
            transition_frequency=[4.0],
            final_time=1.0,
            ntime=10,
            verbose=False,
        )
        config = Config(setup, True)
        assert config is not None
