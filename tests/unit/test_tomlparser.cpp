#include <gtest/gtest.h>
#include <mpi.h>

#include "config.hpp"
#include "defs.hpp"

class TomlParserTest : public ::testing::Test {};

TEST_F(TomlParserTest, ParseBasicSettings) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2]
        transition_frequency = [4.1]
        rotation_frequency = [0.0]
        ntime = 500
        dt = 0.05
        decoherence = {type = "none"}
        initial_condition = {type = "basis"}
      )",
      false);

  EXPECT_EQ(config.getNTime(), 500);
  EXPECT_DOUBLE_EQ(config.getDt(), 0.05);
  EXPECT_EQ(config.getDecoherenceType(), DecoherenceType::NONE);
}

TEST_F(TomlParserTest, ParseVectorSettings) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 3]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.8]
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "basis"}
      )",
      false);

  auto nlevels = config.getNLevels();
  EXPECT_EQ(nlevels.size(), 2);
  EXPECT_EQ(nlevels[0], 2);
  EXPECT_EQ(nlevels[1], 3);

  auto transfreq = config.getTransitionFrequency();
  EXPECT_EQ(transfreq.size(), 2);
  EXPECT_DOUBLE_EQ(transfreq[0], 4.1);
  EXPECT_DOUBLE_EQ(transfreq[1], 4.8);
}

TEST_F(TomlParserTest, ParseJklSettingsAllToAll) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2, 2, 2]
        transition_frequency = [4.1, 4.8, 5.2, 5.5]
        ntime = 1000
        dt = 0.1
        initial_condition = {type = "basis"}
        dipole_coupling = 0.5
      )",
      false);

  auto Jkl = config.getDipoleCoupling();
  size_t num_osc = config.getNumOsc();
  size_t num_pairs_osc = (num_osc - 1) * num_osc / 2;
  EXPECT_EQ(Jkl.size(), num_pairs_osc);
  for (size_t i = 0; i < num_pairs_osc; ++i) {
    EXPECT_DOUBLE_EQ(Jkl[i], 0.5);
  }
}

TEST_F(TomlParserTest, ParseJklSettingsOneCoupling) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2, 2, 2]
        transition_frequency = [4.1, 4.8, 5.2, 5.5]
        ntime = 1000
        dt = 0.1
        initial_condition = {type = "basis"}
        dipole_coupling = [
        { subsystem=[1,2], value=0.4 },
        ]
      )",
      false);

  auto Jkl = config.getDipoleCoupling();
  size_t num_osc = config.getNumOsc();
  size_t num_pairs_osc = (num_osc - 1) * num_osc / 2;
  EXPECT_EQ(Jkl.size(), num_pairs_osc);
  EXPECT_DOUBLE_EQ(Jkl[0], 0.0);
  EXPECT_DOUBLE_EQ(Jkl[1], 0.0);
  EXPECT_DOUBLE_EQ(Jkl[2], 0.0);
  EXPECT_DOUBLE_EQ(Jkl[3], 0.4);
  EXPECT_DOUBLE_EQ(Jkl[4], 0.0);
  EXPECT_DOUBLE_EQ(Jkl[5], 0.0);
}

TEST_F(TomlParserTest, ParseJklSettingsPerPair) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2, 2, 2]
        transition_frequency = [4.1, 4.8, 5.2, 5.5]
        ntime = 1000
        dt = 0.1
        initial_condition = {type = "basis"}
        dipole_coupling = [
        { subsystem=[0,1], value=0.1 },
        { subsystem=[0,2], value=0.2 },
        { subsystem=[0,3], value=0.3 },
        { subsystem=[1,2], value=0.4 },
        { subsystem=[1,3], value=0.5 },
        { subsystem=[2,3], value=0.6 }
        ]
      )",
      false);

  auto Jkl = config.getDipoleCoupling();
  size_t num_osc = config.getNumOsc();
  size_t num_pairs_osc = (num_osc - 1) * num_osc / 2;
  EXPECT_EQ(Jkl.size(), num_pairs_osc);
  EXPECT_DOUBLE_EQ(Jkl[0], 0.1);
  EXPECT_DOUBLE_EQ(Jkl[1], 0.2);
  EXPECT_DOUBLE_EQ(Jkl[2], 0.3);
  EXPECT_DOUBLE_EQ(Jkl[3], 0.4);
  EXPECT_DOUBLE_EQ(Jkl[4], 0.5);
  EXPECT_DOUBLE_EQ(Jkl[5], 0.6);
}

TEST_F(TomlParserTest, ParseCrosskerrSettingsAllToAll) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2, 2, 2]
        transition_frequency = [4.1, 4.8, 5.2, 5.5]
        ntime = 1000
        dt = 0.1
        initial_condition = {type = "basis"}
        crosskerr_coupling = 0.5
      )",
      false);

  auto crosskerr_coupling = config.getCrossKerrCoupling();
  size_t num_osc = config.getNumOsc();
  size_t num_pairs_osc = (num_osc - 1) * num_osc / 2;
  EXPECT_EQ(crosskerr_coupling.size(), num_pairs_osc);
  for (size_t i = 0; i < num_pairs_osc; ++i) {
    EXPECT_DOUBLE_EQ(crosskerr_coupling[i], 0.5);
  }
}

TEST_F(TomlParserTest, ParseCrosskerrSettingsOneCoupling) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2, 2, 2]
        transition_frequency = [4.1, 4.8, 5.2, 5.5]
        ntime = 1000
        dt = 0.1
        initial_condition = {type = "basis"}
        crosskerr_coupling = [
        { subsystem=[1,2], value=0.4 },
        ]
      )",
      false);

  auto crosskerr_coupling = config.getCrossKerrCoupling();
  size_t num_osc = config.getNumOsc();
  size_t num_pairs_osc = (num_osc - 1) * num_osc / 2;
  EXPECT_EQ(crosskerr_coupling.size(), num_pairs_osc);
  EXPECT_DOUBLE_EQ(crosskerr_coupling[0], 0.0);
  EXPECT_DOUBLE_EQ(crosskerr_coupling[1], 0.0);
  EXPECT_DOUBLE_EQ(crosskerr_coupling[2], 0.0);
  EXPECT_DOUBLE_EQ(crosskerr_coupling[3], 0.4);
  EXPECT_DOUBLE_EQ(crosskerr_coupling[4], 0.0);
  EXPECT_DOUBLE_EQ(crosskerr_coupling[5], 0.0);
}

TEST_F(TomlParserTest, ParseCrosskerrSettingsPerPair) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2, 2, 2]
        transition_frequency = [4.1, 4.8, 5.2, 5.5]
        ntime = 1000
        dt = 0.1
        initial_condition = {type = "basis"}
        crosskerr_coupling = [
        { subsystem=[0,1], value=0.1 },
        { subsystem=[0,2], value=0.2 },
        { subsystem=[0,3], value=0.3 },
        { subsystem=[1,2], value=0.4 },
        { subsystem=[1,3], value=0.5 },
        { subsystem=[2,3], value=0.6 }
        ]
      )",
      false);

  auto crosskerr_coupling = config.getCrossKerrCoupling();
  size_t num_osc = config.getNumOsc();
  size_t num_pairs_osc = (num_osc - 1) * num_osc / 2;
  EXPECT_EQ(crosskerr_coupling.size(), num_pairs_osc);
  EXPECT_DOUBLE_EQ(crosskerr_coupling[0], 0.1);
  EXPECT_DOUBLE_EQ(crosskerr_coupling[1], 0.2);
  EXPECT_DOUBLE_EQ(crosskerr_coupling[2], 0.3);
  EXPECT_DOUBLE_EQ(crosskerr_coupling[3], 0.4);
  EXPECT_DOUBLE_EQ(crosskerr_coupling[4], 0.5);
  EXPECT_DOUBLE_EQ(crosskerr_coupling[5], 0.6);
}

TEST_F(TomlParserTest, ParseOutputSettings_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2]
        transition_frequency = [4.1, 4.8]
        ntime = 1000
        dt = 0.1
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [output]
        observables = ["population", "fullstate"]
      )",
      false);

  // Verify output settings
  auto output = config.getOutputObservables();
  EXPECT_EQ(output.size(), 2); // Population
  EXPECT_EQ(output[0], OutputType::POPULATION);
  EXPECT_EQ(output[1], OutputType::FULLSTATE);
}

TEST_F(TomlParserTest, ParseStructSettings) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2]
        transition_frequency = [4.1]
        ntime = 1000
        dt = 0.1
        rotation_frequency = [0.0]
        initial_condition = {type = "diagonal", subsystem = [0]}

        [optimization]
        target = {type = "gate", gate_type = "cnot"}
      )",
      false);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::GATE);
  EXPECT_EQ(target.gate_type.value(), GateType::CNOT);

  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::DIAGONAL);
  EXPECT_EQ(initcond.subsystem.value(), std::vector<size_t>{0});
}

TEST_F(TomlParserTest, ApplyDefaults) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1]
        rotation_frequency = [0.0]
        initial_condition = {type = "basis"}
      )",
      false);

  // Check defaults were applied
  EXPECT_EQ(config.getNTime(), 1000); // Default ntime
  EXPECT_DOUBLE_EQ(config.getDt(), 0.1); // Default dt
  EXPECT_EQ(config.getDecoherenceType(), DecoherenceType::NONE); // Default
}


TEST_F(TomlParserTest, Decoherence_Decay) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [3, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.8]
        rotation_frequency = [0.0, 0.0]
        decoherence = {type = "decay", decay_time = [30.0, 40.0], dephase_time = [20.0, 25.0]}
        initial_condition = {type = "basis"}
      )",
      false);
  EXPECT_EQ(config.getDecoherenceType(), DecoherenceType::DECAY);
  EXPECT_DOUBLE_EQ(config.getDecayTime()[0], 30.0);
  EXPECT_DOUBLE_EQ(config.getDecayTime()[1], 40.0);
  EXPECT_DOUBLE_EQ(config.getDephaseTime()[0], 0.0); // Overwritten to 0
  EXPECT_DOUBLE_EQ(config.getDephaseTime()[1], 0.0); // Overwritten to 0
}

TEST_F(TomlParserTest, Decoherence_Dephase) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [3, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.8]
        rotation_frequency = [0.0, 0.0]
        decoherence = {type = "dephase", decay_time = [30.0, 40.0], dephase_time = [20.0, 25.0]}
        initial_condition = {type = "basis"}
      )",
      false);
  EXPECT_EQ(config.getDecoherenceType(), DecoherenceType::DEPHASE);
  EXPECT_DOUBLE_EQ(config.getDecayTime()[0], 0.0); // Overwritten to 0
  EXPECT_DOUBLE_EQ(config.getDecayTime()[1], 0.0); // Overwritten to 0
  EXPECT_DOUBLE_EQ(config.getDephaseTime()[0], 20.0);
  EXPECT_DOUBLE_EQ(config.getDephaseTime()[1], 25.0);
}


TEST_F(TomlParserTest, InitialCondition_FromFile) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1]
        rotation_frequency = [0.0]
        initial_condition = {type = "file", filename = "test.dat"}
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::FROMFILE);
  EXPECT_EQ(initcond.filename.value(), "test.dat");
  EXPECT_EQ(config.getNInitialConditions(), 1);
}

TEST_F(TomlParserTest, InitialCondition_ProductState) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [3, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.8]
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "state", levels = [1, 0]}
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::PRODUCT_STATE);
  EXPECT_EQ(initcond.levels.value(), std::vector<size_t>({1, 0}));
  EXPECT_EQ(config.getNInitialConditions(), 1);
}

TEST_F(TomlParserTest, InitialCondition_Performance) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1]
        rotation_frequency = [0.0]
        initial_condition = {type = "performance"}
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::PERFORMANCE);
  EXPECT_EQ(config.getNInitialConditions(), 1);
}

TEST_F(TomlParserTest, InitialCondition_Ensemble) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [3, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.8]
        rotation_frequency = [0.0, 0.0]
        decoherence = {type = "decay"}
        initial_condition = {type = "ensemble", subsystem = [0, 1]}
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::ENSEMBLE);
  EXPECT_EQ(initcond.subsystem.value(), std::vector<size_t>({0, 1}));
  EXPECT_EQ(config.getNInitialConditions(), 1);
}

TEST_F(TomlParserTest, InitialCondition_ThreeStates) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [3]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1]
        rotation_frequency = [0.0]
        decoherence = {type = "decay"}
        initial_condition = {type = "3states"}
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::THREESTATES);
  EXPECT_EQ(config.getNInitialConditions(), 3);
}

TEST_F(TomlParserTest, InitialCondition_NPlusOne_SingleOscillator) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [3]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1]
        rotation_frequency = [0.0]
        decoherence = {type = "decay"}
        initial_condition = {type = "nplus1"}
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::NPLUSONE);
  // For nlevels = [3], system dimension N = 3, so n_initial_conditions = N + 1 = 4
  EXPECT_EQ(config.getNInitialConditions(), 4);
}

TEST_F(TomlParserTest, InitialCondition_NPlusOne_MultipleOscillators) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 3]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.8]
        rotation_frequency = [0.0, 0.0]
        decoherence = {type = "decay"}
        initial_condition = {type = "nplus1"}
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::NPLUSONE);
  // For nlevels = [2, 3], system dimension N = 2 * 3 = 6, so n_initial_conditions = N + 1 = 7
  EXPECT_EQ(config.getNInitialConditions(), 7);
}

TEST_F(TomlParserTest, InitialCondition_Diagonal_Schrodinger) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [3, 2]
        ntime = 1000
        dt = 0.1
        nessential = [3, 2]
        transition_frequency = [4.1, 4.8]
        rotation_frequency = [0.0, 0.0]
        decoherence = {type = "none"}
        initial_condition = {type = "diagonal", subsystem = [1]}
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::DIAGONAL);
  EXPECT_EQ(initcond.subsystem.value(), std::vector<size_t>({1}));
  // For Schrodinger solver (decoherence_type = none), n_initial_conditions = nessential[1] = 2
  EXPECT_EQ(config.getNInitialConditions(), 2);
}

TEST_F(TomlParserTest, InitialCondition_Basis_Schrodinger) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [3, 2]
        ntime = 1000
        dt = 0.1
        nessential = [3, 2]
        transition_frequency = [4.1, 4.8]
        rotation_frequency = [0.0, 0.0]
        decoherence = {type = "none"}
        initial_condition = {type = "basis", subsystem = [1]}
      )",
      false);
  // For Schrodinger solver, BASIS is converted to DIAGONAL, so n_initial_conditions = nessential[1] = 2
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::DIAGONAL);
  EXPECT_EQ(initcond.subsystem.value(), std::vector<size_t>({1}));
  EXPECT_EQ(config.getNInitialConditions(), 2);
}

TEST_F(TomlParserTest, InitialCondition_Basis_Lindblad) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [3, 2]
        ntime = 1000
        dt = 0.1
        nessential = [3, 2]
        transition_frequency = [4.1, 4.8]
        rotation_frequency = [0.0, 0.0]
        decoherence = {type = "decay"}
        initial_condition = {type = "basis", subsystem = [1]}
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::BASIS);
  EXPECT_EQ(initcond.subsystem.value(), std::vector<size_t>({1}));
  // For Lindblad solver, n_initial_conditions = nessential[1]^2 = 2^2 = 4
  EXPECT_EQ(config.getNInitialConditions(), 4);
}

TEST_F(TomlParserTest, ControlParameterizations_Spline0) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1]
        rotation_frequency = [0.0]
        initial_condition = {type = "basis"}

        [control]
        parameterization = { type = "spline0", num = 150, tstart = 0.0, tstop = 1.0 }
      )",
      false);

  const auto& control_parameterizations = config.getControlParameterizations(0);
  EXPECT_EQ(control_parameterizations.type, ControlType::BSPLINE0);
  EXPECT_EQ(control_parameterizations.nspline.value(), 150);
  EXPECT_DOUBLE_EQ(control_parameterizations.tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_parameterizations.tstop.value(), 1.0);
}

TEST_F(TomlParserTest, ControlParameterizations_Spline) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.1]
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [control]
        parameterization = [
          { subsystem = 0, type = "spline", num = 10 },
          { subsystem = 1, type = "spline", num = 20, tstart = 0.0, tstop = 1.0 }
        ]
      )",
      false);

  // Check first oscillator with one parameterization
  const auto& control_seg0 = config.getControlParameterizations(0);
  EXPECT_EQ(control_seg0.type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg0.nspline.value(), 10);

  // Check second oscillator
  const auto& control_seg1 = config.getControlParameterizations(1);
  EXPECT_EQ(control_seg1.type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg1.nspline.value(), 20);
  EXPECT_DOUBLE_EQ(control_seg1.tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg1.tstop.value(), 1.0);
}

TEST_F(TomlParserTest, ControlParameterizations_Defaults) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.8, 5.2]
        rotation_frequency = 0.0
        initial_condition = {type = "basis"}

        [control]
        parameterization = [
          { subsystem = 1, type = "spline0", num = 150, tstart = 0.0, tstop = 1.0 }
        ]
      )",
      false);

  // Check first oscillator has default settings
  const auto& control_seg0 = config.getControlParameterizations(0);
  EXPECT_EQ(control_seg0.type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg0.nspline.value(), 10);

  // Check second oscillator has given settings
  const auto& control_seg1 = config.getControlParameterizations(1);
  EXPECT_EQ(control_seg1.type, ControlType::BSPLINE0);
  EXPECT_EQ(control_seg1.nspline.value(), 150);
  EXPECT_DOUBLE_EQ(control_seg1.tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg1.tstop.value(), 1.0);

  // Check third oscillator has default settings
  const auto& control_seg2 = config.getControlParameterizations(2);
  EXPECT_EQ(control_seg2.type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg2.nspline.value(), 10);

  // Check
}

TEST_F(TomlParserTest, ControlParameterizations_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.8]
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [control]
        parameterization = { type = "spline0", num = 150, tstart = 0.0, tstop = 1.0 }
      )",
      false);

  // Check first oscillator has given settings
  const auto& control_seg0 = config.getControlParameterizations(0);
  EXPECT_EQ(control_seg0.type, ControlType::BSPLINE0);
  EXPECT_EQ(control_seg0.nspline.value(), 150);
  EXPECT_DOUBLE_EQ(control_seg0.tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg0.tstop.value(), 1.0);

  // Check second oscillator has given settings
  const auto& control_seg1 = config.getControlParameterizations(1);
  EXPECT_EQ(control_seg1.type, ControlType::BSPLINE0);
  EXPECT_EQ(control_seg1.nspline.value(), 150);
  EXPECT_DOUBLE_EQ(control_seg1.tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg1.tstop.value(), 1.0);
}

TEST_F(TomlParserTest, ControlInitialization_Defaults) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.1, 4.1]
        rotation_frequency = [0.0, 0.0, 0.0]
        initial_condition = {type = "basis"}

        [control]
        initialization = [
           { subsystem = 1, type = "random", amplitude = 2.0 }
        ]
      )",
      false);

  // Check first oscillator has default settings
  const auto& control_init0 = config.getControlInitializations(0);
  EXPECT_EQ(control_init0.type, ControlInitializationType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init0.amplitude.value(), 0.0);

  // Check second oscillator has given settings
  const auto& control_init1 = config.getControlInitializations(1);
  EXPECT_EQ(control_init1.type, ControlInitializationType::RANDOM);
  EXPECT_DOUBLE_EQ(control_init1.amplitude.value(), 2.0);

  // Check third oscillator has default settings
  const auto& control_init2 = config.getControlInitializations(2);
  EXPECT_EQ(control_init2.type, ControlInitializationType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init2.amplitude.value(), 0.0);
}

TEST_F(TomlParserTest, ControlInitializationSettings) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2, 2, 2, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.1, 4.1, 4.1, 4.1]
        rotation_frequency = [0.0, 0.0, 0.0, 0.0, 0.0]
        initial_condition = {type = "basis"}

        [control]
        initialization = [
          { subsystem = 0, type = "constant", amplitude = 1.0, phase = 1.1 },
          { subsystem = 1, type = "constant", amplitude = 2.0 },
          { subsystem = 2, type = "random", amplitude = 3.0, phase = 3.1 },
          { subsystem = 3, type = "random", amplitude = 4.0 },
          { subsystem = 4, type = "random", amplitude = 5.0, phase = 5.1 }
        ]

        parameterization = [
          { subsystem = 0, type = "spline_amplitude", num = 10, scaling = 1.0 },
          { subsystem = 2, type = "spline_amplitude", num = 10, scaling = 1.0 },
          { subsystem = 4, type = "spline_amplitude", num = 10, scaling = 1.0 }
        ]
      )",
      false);

  // Check first oscillator
  const auto& control_init0 = config.getControlInitializations(0);
  EXPECT_EQ(control_init0.type, ControlInitializationType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init0.amplitude.value(), 1.0);
  EXPECT_DOUBLE_EQ(control_init0.phase.value(), 1.1);

  // Check second oscillator
  const auto& control_init1 = config.getControlInitializations(1);
  EXPECT_EQ(control_init1.type, ControlInitializationType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init1.amplitude.value(), 2.0);

  // Check third oscillator
  const auto& control_init2 = config.getControlInitializations(2);
  EXPECT_EQ(control_init2.type, ControlInitializationType::RANDOM);
  EXPECT_DOUBLE_EQ(control_init2.amplitude.value(), 3.0);
  EXPECT_DOUBLE_EQ(control_init2.phase.value(), 3.1);

  // Check fourth oscillator
  const auto& control_init3 = config.getControlInitializations(3);
  EXPECT_EQ(control_init3.type, ControlInitializationType::RANDOM);
  EXPECT_DOUBLE_EQ(control_init3.amplitude.value(), 4.0);

  // Check fifth oscillator with two parameterizations (init is copied to match parameterization count)
  const auto& control_init4 = config.getControlInitializations(4);
  EXPECT_EQ(control_init4.type, ControlInitializationType::RANDOM);
  EXPECT_DOUBLE_EQ(control_init4.amplitude.value(), 5.0);
  EXPECT_DOUBLE_EQ(control_init4.phase.value(), 5.1);
}

TEST_F(TomlParserTest, ControlInitialization_File) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1]
        rotation_frequency = [0.0]
        initial_condition = {type = "basis"}

        [control]
        initialization = { type = "file", filename = "myparams.dat" }
      )",
      false);

  const auto& control_init0 = config.getControlInitializations(0);
  EXPECT_EQ(control_init0.type, ControlInitializationType::FILE);
  EXPECT_EQ(control_init0.filename.value(), "myparams.dat");
}

TEST_F(TomlParserTest, ControlInitialization_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.1]
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [control]
        initialization = { type = "constant", amplitude = 1.0 }
      )",
      false);

  // Check first oscillator
  const auto& control_init0 = config.getControlInitializations(0);
  EXPECT_EQ(control_init0.type, ControlInitializationType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init0.amplitude.value(), 1.0);

  // Check second oscillator
  const auto& control_init1 = config.getControlInitializations(1);
  EXPECT_EQ(control_init1.type, ControlInitializationType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init1.amplitude.value(), 1.0);
}

TEST_F(TomlParserTest, ControlInitialization_DefaultWithOverrides) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [3, 3, 3]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.1, 4.1]
        rotation_frequency = [0.0, 0.0, 0.0]
        initial_condition = {type = "basis"}

        [control]
        initialization = [
          { subsystem = 1, type = "random", amplitude = 0.05 }
        ]
      )",
      false);

  // Check first oscillator has default settings
  const auto& control_init0 = config.getControlInitializations(0);
  EXPECT_EQ(control_init0.type, ControlInitializationType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init0.amplitude.value(), 0.0);

  // Check second oscillator has override settings
  const auto& control_init1 = config.getControlInitializations(1);
  EXPECT_EQ(control_init1.type, ControlInitializationType::RANDOM);
  EXPECT_DOUBLE_EQ(control_init1.amplitude.value(), 0.05);

  // Check third oscillator has default settings
  const auto& control_init2 = config.getControlInitializations(2);
  EXPECT_EQ(control_init2.type, ControlInitializationType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init2.amplitude.value(), 0.0);
}

TEST_F(TomlParserTest, ControlBounds) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 3.3]
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [control]
        parameterization = { type = "spline", num = 10, tstart = 0.0, tstop = 1.0 }
        amplitude_bound = 1.5
      )",
      false);

  // Check control bound
  const double control_bound0 = config.getControlAmplitudeBound(0);
  const double control_bound1 = config.getControlAmplitudeBound(1);
  EXPECT_DOUBLE_EQ(control_bound0, 1.5);
  EXPECT_DOUBLE_EQ(control_bound1, 1.5);
}

TEST_F(TomlParserTest, ControlBounds_Defaults) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 3.0]
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "basis"}
      )",
      false);

  // Check control bound
  const double control_bound0 = config.getControlAmplitudeBound(0);
  const double control_bound1 = config.getControlAmplitudeBound(1);
  EXPECT_DOUBLE_EQ(control_bound0, 1e12);
  EXPECT_DOUBLE_EQ(control_bound1, 1e12);
}

TEST_F(TomlParserTest, ControlBounds_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 3.3]
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [control]
        amplitude_bound = [1.0, 2.0]
      )",
      false);

  // Check control bound
  const double control_bound0 = config.getControlAmplitudeBound(0);
  const double control_bound1 = config.getControlAmplitudeBound(1);
  EXPECT_DOUBLE_EQ(control_bound0, 1.0);
  EXPECT_DOUBLE_EQ(control_bound1, 2.0);
}

TEST_F(TomlParserTest, CarrierFrequenciesDefaults) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2,2]
        ntime = 1000
        dt = 0.1
        transition_frequency = 4.1
        rotation_frequency = 0.0
        initial_condition = {type = "basis"}

        [control]
        carrier_frequency = [
          {subsystem = 1, value = [4.0, 5.0]}
        ]
      )",
      false);

  const auto& carrier_freq0 = config.getCarrierFrequencies(0);
  EXPECT_EQ(carrier_freq0.size(), 1);
  EXPECT_EQ(carrier_freq0[0], 0.0);
  const auto& carrier_freq1 = config.getCarrierFrequencies(1);
  EXPECT_EQ(carrier_freq1.size(), 2);
  EXPECT_EQ(carrier_freq1[0], 4.0);
  EXPECT_EQ(carrier_freq1[1], 5.0);
}

TEST_F(TomlParserTest, CarrierFrequenciesPerSubsystem) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2,2]
        ntime = 1000
        dt = 0.1
        transition_frequency = 4.1
        initial_condition = {type = "basis"}

        [control]
        carrier_frequency = [
          {subsystem = 0, value = [1.0, 2.0, 3.0]},
          {subsystem = 1, value = [4.0, 5.0]}
        ]
      )",
      false);

  const auto& carrier_freq0 = config.getCarrierFrequencies(0);
  EXPECT_EQ(carrier_freq0.size(), 3);
  EXPECT_EQ(carrier_freq0[0], 1.0);
  EXPECT_EQ(carrier_freq0[1], 2.0);
  EXPECT_EQ(carrier_freq0[2], 3.0);
  const auto& carrier_freq1 = config.getCarrierFrequencies(1);
  EXPECT_EQ(carrier_freq1.size(), 2);
  EXPECT_EQ(carrier_freq1[0], 4.0);
  EXPECT_EQ(carrier_freq1[1], 5.0);
}

TEST_F(TomlParserTest, CarrierFrequencies_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.1]
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [control]
        carrier_frequency = { value = [1.0, 2.0] }
      )",
      false);

  const auto& carrier_freq0 = config.getCarrierFrequencies(0);
  EXPECT_EQ(carrier_freq0.size(), 2);
  EXPECT_EQ(carrier_freq0[0], 1.0);
  EXPECT_EQ(carrier_freq0[1], 2.0);

  const auto& carrier_freq1 = config.getCarrierFrequencies(1);
  EXPECT_EQ(carrier_freq1.size(), 2);
  EXPECT_EQ(carrier_freq1[0], 1.0);
  EXPECT_EQ(carrier_freq1[1], 2.0);
}
TEST_F(TomlParserTest, CarrierFrequencies_AllOscillatorsShorthand) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.1]
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [control]
        carrier_frequency = [1.0, 2.0]
      )",
      false);

  const auto& carrier_freq0 = config.getCarrierFrequencies(0);
  EXPECT_EQ(carrier_freq0.size(), 2);
  EXPECT_EQ(carrier_freq0[0], 1.0);
  EXPECT_EQ(carrier_freq0[1], 2.0);

  const auto& carrier_freq1 = config.getCarrierFrequencies(1);
  EXPECT_EQ(carrier_freq1.size(), 2);
  EXPECT_EQ(carrier_freq1[0], 1.0);
  EXPECT_EQ(carrier_freq1[1], 2.0);
}

TEST_F(TomlParserTest, CarrierFrequency_InvalidID) {
  ASSERT_THROW(
      Config::fromTomlString(
          R"(
            [system]
            nlevels = [2]
            transition_frequency = [4.1]
            ntime = 1000
            dt = 0.1
            rotation_frequency = [0.0]
            initial_condition = {type = "basis"}

            [control]
            carrier_frequency = [{subsystem = 5, value = [4.0]}]
          )",
          false),
      validators::ValidationError);
}



TEST_F(TomlParserTest, OptimTarget_GateType) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2,2,2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 3.4, 5.6]
        rotation_frequency = [0.0, 0.0, 0.0]
        initial_condition = {type = "basis"}

        [optimization]
        target = {type = "gate", gate_type = "cnot"}
      )",
      false);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::GATE);
  EXPECT_EQ(target.gate_type.value(), GateType::CNOT);
  EXPECT_FALSE(target.filename.has_value());
  EXPECT_FALSE(target.levels.has_value());
  EXPECT_TRUE(target.gate_rot_freq.has_value());
  const auto& gate_rot_freq = target.gate_rot_freq.value();
  EXPECT_EQ(gate_rot_freq.size(), 3);
  EXPECT_DOUBLE_EQ(gate_rot_freq[0], 0.0);
  EXPECT_DOUBLE_EQ(gate_rot_freq[1], 0.0);
  EXPECT_DOUBLE_EQ(gate_rot_freq[2], 0.0);
}

TEST_F(TomlParserTest, OptimTarget_GateRotFreq) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2,2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 3.4]
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [optimization]
        target = {type = "gate", gate_type = "cnot", gate_rot_freq = [1.0, 2.0]}
      )",
      false);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::GATE);
  EXPECT_FALSE(target.filename.has_value());
  EXPECT_TRUE(target.gate_rot_freq.has_value());
  const auto& gate_rot_freq = target.gate_rot_freq.value();
  EXPECT_EQ(gate_rot_freq.size(), 2);
  EXPECT_DOUBLE_EQ(gate_rot_freq[0], 1.0);
  EXPECT_DOUBLE_EQ(gate_rot_freq[1], 2.0);
}

TEST_F(TomlParserTest, OptimTarget_GateFromFile) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1]
        rotation_frequency = [0.0]
        initial_condition = {type = "basis"}

        [optimization]
        target = {type = "gate", filename = "/path/to/gate.dat"}
      )",
      false);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::GATE);
  EXPECT_EQ(target.gate_type.value(), GateType::FILE);
  EXPECT_EQ(target.filename.value(), "/path/to/gate.dat");
}

TEST_F(TomlParserTest, OptimTarget_ProductState) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [3, 3, 3]
        ntime = 1000
        dt = 0.1
        transition_frequency = 4.1
        rotation_frequency = 0.0
        initial_condition = {type = "basis"}

        [optimization]
        target = {type = "state", levels = [0,1,2]}
      )",
      false);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::STATE);
  EXPECT_TRUE(target.levels.has_value());
  EXPECT_FALSE(target.filename.has_value());
  const auto& levels = target.levels.value();
  EXPECT_EQ(levels.size(), 3);
  EXPECT_EQ(levels[0], 0);
  EXPECT_EQ(levels[1], 1);
  EXPECT_EQ(levels[2], 2);
}

TEST_F(TomlParserTest, OptimTarget_FromFile) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1]
        rotation_frequency = [0.0]
        initial_condition = {type = "basis"}

        [optimization]
        target = {type = "state", filename = "/path/to/target.dat"}
      )",
      false);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::STATE);
  EXPECT_TRUE(target.filename.has_value());
  EXPECT_EQ(target.filename.value(), "/path/to/target.dat");
}

TEST_F(TomlParserTest, OptimTarget_DefaultNone) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1]
        rotation_frequency = [0.0]
        initial_condition = {type = "basis"}
      )",
      false);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::NONE);
  EXPECT_FALSE(target.gate_type.has_value());
  EXPECT_FALSE(target.levels.has_value());
  EXPECT_FALSE(target.filename.has_value());
}

TEST_F(TomlParserTest, OptimWeightsVectorNormalization) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 2]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.1]
        rotation_frequency = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [optimization]
        weights = [0.25, 0.50, 0.25, 0.0]
      )",
      false);

  const auto& weights = config.getOptimWeights();
  EXPECT_EQ(weights.size(), 4);
  EXPECT_DOUBLE_EQ(weights[0], 0.25);
  EXPECT_DOUBLE_EQ(weights[1], 0.5);
  EXPECT_DOUBLE_EQ(weights[2], 0.25);
  EXPECT_DOUBLE_EQ(weights[3], 0.0);
}

TEST_F(TomlParserTest, OptimWeightsScalarNotAllowed) {
  // Single value format is no longer allowed - must use array notation
  ASSERT_THROW(
      Config::fromTomlString(
          R"(
            [system]
            nlevels = [2, 2]
            ntime = 1000
            dt = 0.1
            transition_frequency = [4.1, 4.1]
            rotation_frequency = [0.0, 0.0]
            initial_condition = {type = "basis"}

            [optimization]
            weights = 1.0
          )",
          false),
      validators::ValidationError);
}

TEST_F(TomlParserTest, OptimWeightsDefault) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [4,3]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 3.4]
        initial_condition = {type = "diagonal", subsystem = [0]}
      )",
      false);

  const auto& weights = config.getOptimWeights();
  EXPECT_EQ(weights.size(), 4);
  EXPECT_DOUBLE_EQ(weights[0], 0.25); // normalized
  EXPECT_DOUBLE_EQ(weights[1], 0.25);
  EXPECT_DOUBLE_EQ(weights[2], 0.25);
  EXPECT_DOUBLE_EQ(weights[3], 0.25);
}

TEST_F(TomlParserTest, ComprehensiveNonDefaultSettings_PrintConfigValidation) {
  // Create a comprehensive configuration with non-default settings
  // Uses multi-oscillator system with per-subsystem settings to test output formatting
  std::string input_toml = R"(
[system]
nlevels = [3, 4]
nessential = [2, 3]
ntime = 2500
dt = 0.02
transition_frequency = [5.5, 6.2]
selfkerr = [-0.2, -0.25]
rotation_frequency = [2.1, 2.5]
crosskerr_coupling = [
  {subsystem = [0, 1], value = 0.05}
]
dipole_coupling = [
  {subsystem = [0, 1], value = 0.01}
]
decoherence = {type = "decay", decay_time = [25.0, 30.0], dephase_time = [0.0, 0.0]}
initial_condition = {type = "diagonal", subsystem = [0, 1]}
hamiltonian_file_Hsys = "/path/to/system_hamiltonian.dat"
hamiltonian_file_Hc = "/path/to/control_hamiltonian.dat"

[control]
parameterization = [
  {subsystem = 0, type = "spline", num = 25, tstart = 0.5, tstop = 1.5},
  {subsystem = 1, type = "spline0", num = 30}
]
carrier_frequency = [
  {subsystem = 0, value = [2.5, 3.0]},
  {subsystem = 1, value = [1.5, 2.0, 2.5]}
]
initialization = [
  {subsystem = 0, type = "random", amplitude = 1.5, phase = 0.5},
  {subsystem = 1, type = "constant", amplitude = 0.1}
]
amplitude_bound = [5.0, 6.0]
zero_boundary_condition = false

[optimization]
target = {type = "gate", gate_type = "cnot", gate_rot_freq = [1.5, 2.0]}
objective = "jtrace"
weights = [0.5, 0.3, 0.1, 0.05, 0.03, 0.02]
tolerance = { grad_abs = 1e-06, grad_rel = 0.001, final_cost = 1e-07, infidelity = 0.0001 }
maxiter = 150
tikhonov = { coeff = 0.002, use_x0 = true }
penalty = { leakage = 0.1, weightedcost = 0.8, weightedcost_width = 0.5, dpdm = 0.05, energy = 0.02, variation = 0.03 }

[output]
directory = "/custom/output/path"
observables = ["population", "expectedenergy"]
timestep_stride = 5
optimization_stride = 25

[solver]
runtype = "optimization"
usematfree = true
linearsolver = { type = "gmres", maxiter = 20 }
timestepper = "imr"
rand_seed = 12345
)";

  Config config = Config::fromTomlString(input_toml, false);

  // Test that printed config is valid TOML that can be parsed
  std::stringstream printed_output;
  config.printConfig(printed_output);
  std::string output = printed_output.str();

  // Verify output is valid TOML by parsing it
  ASSERT_NO_THROW({
    Config::fromTomlString(output, false);
  });
}

TEST_F(TomlParserTest, Transfreq_ScalarValue) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 3, 4]
        ntime = 1000
        dt = 0.1
        transition_frequency = 5.5
        initial_condition = {type = "basis"}
      )",
      false);

  auto transfreq = config.getTransitionFrequency();
  EXPECT_EQ(transfreq.size(), 3);
  EXPECT_DOUBLE_EQ(transfreq[0], 5.5);
  EXPECT_DOUBLE_EQ(transfreq[1], 5.5);
  EXPECT_DOUBLE_EQ(transfreq[2], 5.5);
}

TEST_F(TomlParserTest, Transfreq_ArrayValue) {
  Config config = Config::fromTomlString(
      R"(
        [system]
        nlevels = [2, 3, 4]
        ntime = 1000
        dt = 0.1
        transition_frequency = [4.1, 4.8, 5.2]
        initial_condition = {type = "basis"}
      )",
      false);

  auto transfreq = config.getTransitionFrequency();
  EXPECT_EQ(transfreq.size(), 3);
  EXPECT_DOUBLE_EQ(transfreq[0], 4.1);
  EXPECT_DOUBLE_EQ(transfreq[1], 4.8);
  EXPECT_DOUBLE_EQ(transfreq[2], 5.2);
}

TEST_F(TomlParserTest, Transfreq_WrongSizeError) {
  ASSERT_THROW(
      Config::fromTomlString(
          R"(
            [system]
            nlevels = [2, 3, 4]
            ntime = 1000
            dt = 0.1
            transition_frequency = [4.1, 4.8]
            initial_condition = {type = "basis"}
          )",
          false),
      validators::ValidationError);
}

// Test that hasLength validator correctly rejects arrays with wrong number of elements.
// The subsystem field for coupling parameters (dipole_coupling, crosskerr_coupling) must have exactly 2 elements.
TEST_F(TomlParserTest, CouplingSubsystem_HasLength_TooManyElements) {
  ASSERT_THROW(
      Config::fromTomlString(
          R"(
            [system]
            nlevels = [2, 2, 2, 2]
            transition_frequency = [4.1, 4.8, 5.2, 5.5]
            ntime = 1000
            dt = 0.1
            initial_condition = {type = "basis"}
            dipole_coupling = [
              { subsystem = [0, 1, 2], value = 0.5 }
            ]
          )",
          false),
      validators::ValidationError);
}

TEST_F(TomlParserTest, CouplingSubsystem_HasLength_TooFewElements) {
  ASSERT_THROW(
      Config::fromTomlString(
          R"(
            [system]
            nlevels = [2, 2, 2, 2]
            transition_frequency = [4.1, 4.8, 5.2, 5.5]
            ntime = 1000
            dt = 0.1
            initial_condition = {type = "basis"}
            crosskerr_coupling = [
              { subsystem = [1], value = 0.5 }
            ]
          )",
          false),
      validators::ValidationError);
}
