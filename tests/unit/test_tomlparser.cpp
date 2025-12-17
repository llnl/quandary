#include <gtest/gtest.h>
#include <mpi.h>

#include "config.hpp"
#include "defs.hpp"

class TomlParserTest : public ::testing::Test {
 protected:
  MPILogger logger = MPILogger(0, false);
};

TEST_F(TomlParserTest, ParseBasicSettings) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        ntime = 500
        dt = 0.05
        collapse_type = "none"
        initial_condition = {type = "basis"}
      )",
      logger);

  EXPECT_EQ(config.getNTime(), 500);
  EXPECT_DOUBLE_EQ(config.getDt(), 0.05);
  EXPECT_EQ(config.getCollapseType(), LindbladType::NONE);
}

TEST_F(TomlParserTest, ParseVectorSettings) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 3]
        transfreq = [4.1, 4.8, 5.2]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}
      )",
      logger);

  auto nlevels = config.getNLevels();
  EXPECT_EQ(nlevels.size(), 2);
  EXPECT_EQ(nlevels[0], 2);
  EXPECT_EQ(nlevels[1], 3);

  auto transfreq = config.getTransFreq();
  EXPECT_EQ(transfreq.size(), 3);
  EXPECT_DOUBLE_EQ(transfreq[0], 4.1);
  EXPECT_DOUBLE_EQ(transfreq[1], 4.8);
  EXPECT_DOUBLE_EQ(transfreq[2], 5.2);
}

TEST_F(TomlParserTest, ParseOutputSettings) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}
        [[write]]
        oscID = 0
        type = ["population"]
        [[write]]
        oscID = 1
        type = ["population", "expectedEnergy"]
      )",
      logger);

  // Verify output settings
  auto output = config.getOutput();
  EXPECT_EQ(output.size(), 2); // 2 oscillators
  EXPECT_EQ(output[0].size(), 1); // 1 output
  EXPECT_EQ(output[0][0], OutputType::POPULATION);
  EXPECT_EQ(output[1].size(), 2); // 2 outputs
  EXPECT_EQ(output[1][0], OutputType::POPULATION);
  EXPECT_EQ(output[1][1], OutputType::EXPECTED_ENERGY);
}

TEST_F(TomlParserTest, ParseOutputSettings_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}
        [[write]]
        type = ["population"]
      )",
      logger);

  // Verify output settings
  auto output = config.getOutput();
  EXPECT_EQ(output.size(), 2); // 2 oscillators
  EXPECT_EQ(output[0].size(), 1); // 1 output
  EXPECT_EQ(output[0][0], OutputType::POPULATION);
  EXPECT_EQ(output[1].size(), 1); // 1 output
  EXPECT_EQ(output[1][0], OutputType::POPULATION);
}

TEST_F(TomlParserTest, ParseStructSettings) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "diagonal", oscIDs = [0]}
        optim_target = {target_type = "gate", gate_type = "cnot"}
      )",
      logger);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::GATE);
  EXPECT_EQ(target.gate_type.value(), GateType::CNOT);

  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::DIAGONAL);
  EXPECT_EQ(initcond.osc_IDs.value(), std::vector<size_t>{0});
}

TEST_F(TomlParserTest, ApplyDefaults) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}
      )",
      logger);

  // Check defaults were applied
  EXPECT_EQ(config.getNTime(), 1000); // Default ntime
  EXPECT_DOUBLE_EQ(config.getDt(), 0.1); // Default dt
  EXPECT_EQ(config.getCollapseType(), LindbladType::NONE); // Default
}

TEST_F(TomlParserTest, InitialCondition_FromFile) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "file", filename = "test.dat"}
      )",
      logger);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::FROMFILE);
  EXPECT_EQ(initcond.filename.value(), "test.dat");
  EXPECT_EQ(config.getNInitialConditions(), 1);
}

TEST_F(TomlParserTest, InitialCondition_Pure) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [3, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "pure", levels = [1, 0]}
      )",
      logger);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::PURE);
  EXPECT_EQ(initcond.levels.value(), std::vector<size_t>({1, 0}));
  EXPECT_EQ(config.getNInitialConditions(), 1);
}

TEST_F(TomlParserTest, InitialCondition_Performance) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "performance"}
      )",
      logger);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::PERFORMANCE);
  EXPECT_EQ(config.getNInitialConditions(), 1);
}

TEST_F(TomlParserTest, InitialCondition_Ensemble) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [3, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        collapse_type = "decay"
        initial_condition = {type = "ensemble", oscIDs = [0, 1]}
      )",
      logger);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::ENSEMBLE);
  EXPECT_EQ(initcond.osc_IDs.value(), std::vector<size_t>({0, 1}));
  EXPECT_EQ(config.getNInitialConditions(), 1);
}

TEST_F(TomlParserTest, InitialCondition_ThreeStates) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [3]
        transfreq = [4.1]
        rotfreq = [0.0]
        collapse_type = "decay"
        initial_condition = {type = "3states"}
      )",
      logger);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::THREESTATES);
  EXPECT_EQ(config.getNInitialConditions(), 3);
}

TEST_F(TomlParserTest, InitialCondition_NPlusOne_SingleOscillator) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [3]
        transfreq = [4.1]
        rotfreq = [0.0]
        collapse_type = "decay"
        initial_condition = {type = "nplus1"}
      )",
      logger);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::NPLUSONE);
  // For nlevels = [3], system dimension N = 3, so n_initial_conditions = N + 1 = 4
  EXPECT_EQ(config.getNInitialConditions(), 4);
}

TEST_F(TomlParserTest, InitialCondition_NPlusOne_MultipleOscillators) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 3]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        collapse_type = "decay"
        initial_condition = {type = "nplus1"}
      )",
      logger);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::NPLUSONE);
  // For nlevels = [2, 3], system dimension N = 2 * 3 = 6, so n_initial_conditions = N + 1 = 7
  EXPECT_EQ(config.getNInitialConditions(), 7);
}

TEST_F(TomlParserTest, InitialCondition_Diagonal_Schrodinger) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [3, 2]
        nessential = [3, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        collapse_type = "none"
        initial_condition = {type = "diagonal", oscIDs = [1]}
      )",
      logger);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::DIAGONAL);
  EXPECT_EQ(initcond.osc_IDs.value(), std::vector<size_t>({1}));
  // For Schrodinger solver (collapse_type = none), n_initial_conditions = nessential[1] = 2
  EXPECT_EQ(config.getNInitialConditions(), 2);
}

TEST_F(TomlParserTest, InitialCondition_Basis_Schrodinger) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [3, 2]
        nessential = [3, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        collapse_type = "none"
        initial_condition = {type = "basis", oscIDs = [1]}
      )",
      logger);
  // For Schrodinger solver, BASIS is converted to DIAGONAL, so n_initial_conditions = nessential[1] = 2
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::DIAGONAL);
  EXPECT_EQ(initcond.osc_IDs.value(), std::vector<size_t>({1}));
  EXPECT_EQ(config.getNInitialConditions(), 2);
}

TEST_F(TomlParserTest, InitialCondition_Basis_Lindblad) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [3, 2]
        nessential = [3, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        collapse_type = "decay"
        initial_condition = {type = "basis", oscIDs = [1]}
      )",
      logger);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::BASIS);
  EXPECT_EQ(initcond.osc_IDs.value(), std::vector<size_t>({1}));
  // For Lindblad solver, n_initial_conditions = nessential[1]^2 = 2^2 = 4
  EXPECT_EQ(config.getNInitialConditions(), 4);
}

TEST_F(TomlParserTest, ParsePiPulseSettings_Structure) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}

        [[apply_pipulse]]
        oscID = 0
        tstart = 0.5
        tstop = 1.0
        amp = 0.8
      )",
      logger);

  const auto& pulses = config.getApplyPiPulses();
  EXPECT_EQ(pulses.size(), 2);

  EXPECT_EQ(pulses[0].size(), 1);
  EXPECT_DOUBLE_EQ(pulses[0][0].tstart, 0.5);
  EXPECT_DOUBLE_EQ(pulses[0][0].tstop, 1.0);
  EXPECT_DOUBLE_EQ(pulses[0][0].amp, 0.8);

  // zero pulse for other oscillator
  EXPECT_EQ(pulses[1].size(), 1);
  EXPECT_DOUBLE_EQ(pulses[1][0].tstart, 0.5);
  EXPECT_DOUBLE_EQ(pulses[1][0].tstop, 1.0);
  EXPECT_DOUBLE_EQ(pulses[1][0].amp, 0.0);
}

TEST_F(TomlParserTest, ParsePiPulseSettings_Multiple) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}

        [[apply_pipulse]]
        oscID = 0
        tstart = 0.5
        tstop = 1.0
        amp = 0.8

        [[apply_pipulse]]
        oscID = 1
        tstart = 0
        tstop = 0.5
        amp = 0.2
      )",
      logger);

  const auto& pulses = config.getApplyPiPulses();
  EXPECT_EQ(pulses.size(), 2);

  EXPECT_EQ(pulses[0].size(), 2);
  EXPECT_DOUBLE_EQ(pulses[0][0].tstart, 0.5);
  EXPECT_DOUBLE_EQ(pulses[0][0].tstop, 1.0);
  EXPECT_DOUBLE_EQ(pulses[0][0].amp, 0.8);
  EXPECT_DOUBLE_EQ(pulses[0][1].tstart, 0.);
  EXPECT_DOUBLE_EQ(pulses[0][1].tstop, 0.5);
  EXPECT_DOUBLE_EQ(pulses[0][1].amp, 0.0);

  EXPECT_EQ(pulses[1].size(), 2);
  EXPECT_DOUBLE_EQ(pulses[1][0].tstart, 0.5);
  EXPECT_DOUBLE_EQ(pulses[1][0].tstop, 1.0);
  EXPECT_DOUBLE_EQ(pulses[1][0].amp, 0.0);
  EXPECT_DOUBLE_EQ(pulses[1][1].tstart, 0.);
  EXPECT_DOUBLE_EQ(pulses[1][1].tstop, 0.5);
  EXPECT_DOUBLE_EQ(pulses[1][1].amp, 0.2);
}

TEST_F(TomlParserTest, ControlSegments_Spline0) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}

        [[control_segments]]
        oscID = 0
        type = "spline0"
        num = 150
        tstart = 0.0
        tstop = 1.0
      )",
      logger);

  const auto& control_segments = config.getControlSegments(0);
  EXPECT_EQ(control_segments.size(), 1);
  EXPECT_EQ(control_segments[0].type, ControlType::BSPLINE0);
  EXPECT_EQ(control_segments[0].nspline.value(), 150);
  EXPECT_DOUBLE_EQ(control_segments[0].tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_segments[0].tstop.value(), 1.0);
}

TEST_F(TomlParserTest, ControlSegments_Spline) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 4.1]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [[control_segments]]
        oscID = 0
        type = "spline"
        num = 10
        [[control_segments]]
        oscID = 1
        type = "spline"
        num = 20
        tstart = 0.0
        tstop = 1.0
        [[control_segments]]
        oscID = 1
        type = "spline"
        num = 30
        tstart = 1.0
        tstop = 2.0
      )",
      logger);

  // Check first oscillator with one segment
  const auto& control_seg0 = config.getControlSegments(0);
  EXPECT_EQ(control_seg0.size(), 1);
  EXPECT_EQ(control_seg0[0].type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg0[0].nspline.value(), 10);
  EXPECT_DOUBLE_EQ(control_seg0[0].tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg0[0].tstop.value(), config.getTotalTime());

  // Check second oscillator with two segments
  const auto& control_seg1 = config.getControlSegments(1);
  EXPECT_EQ(control_seg1.size(), 2);

  EXPECT_EQ(control_seg1[0].type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg1[0].nspline.value(), 20);
  EXPECT_DOUBLE_EQ(control_seg1[0].tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg1[0].tstop.value(), 1.0);

  EXPECT_EQ(control_seg1[1].type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg1[1].nspline.value(), 30);
  EXPECT_DOUBLE_EQ(control_seg1[1].tstart.value(), 1.0);
  EXPECT_DOUBLE_EQ(control_seg1[1].tstop.value(), 2.0);
}

TEST_F(TomlParserTest, ControlSegments_Step) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2,2]
        transfreq = [4.1,4.1]
        rotfreq = [0.0,0.0]
        initial_condition = {type = "basis"}

        [[control_segments]]
        oscID = 0
        type = "step"
        step_amp1 = 0.1
        step_amp2 = 0.2
        tramp = 0.3
        tstart = 0.4
        tstop = 0.5
        [[control_segments]]
        oscID = 1
        type = "step"
        step_amp1 = 0.1
        step_amp2 = 0.2
        tramp = 0.3
      )",
      logger);

  // Check first oscillator
  const auto& control_seg0 = config.getControlSegments(0);
  EXPECT_EQ(control_seg0.size(), 1);
  EXPECT_EQ(control_seg0[0].type, ControlType::STEP);
  EXPECT_EQ(control_seg0[0].step_amp1.value(), 0.1);
  EXPECT_DOUBLE_EQ(control_seg0[0].step_amp2.value(), 0.2);
  EXPECT_DOUBLE_EQ(control_seg0[0].tramp.value(), 0.3);
  EXPECT_DOUBLE_EQ(control_seg0[0].tstart.value(), 0.4);
  EXPECT_DOUBLE_EQ(control_seg0[0].tstop.value(), 0.5);

  // Check second oscillator
  const auto& control_seg1 = config.getControlSegments(1);
  EXPECT_EQ(control_seg1.size(), 1);
  EXPECT_EQ(control_seg1[0].type, ControlType::STEP);
  EXPECT_EQ(control_seg1[0].step_amp1.value(), 0.1);
  EXPECT_DOUBLE_EQ(control_seg1[0].step_amp2.value(), 0.2);
  EXPECT_DOUBLE_EQ(control_seg1[0].tramp.value(), 0.3);
  EXPECT_DOUBLE_EQ(control_seg1[0].tstart.value(), 0.0); // default start time
  EXPECT_DOUBLE_EQ(control_seg1[0].tstop.value(), config.getTotalTime()); // default stop time
}

TEST_F(TomlParserTest, ControlSegments_Defaults) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [[control_segments]]
        oscID = 1
        type = "spline0"
        num = 150
        tstart = 0.0
        tstop = 1.0
        [[control_bounds]]
        oscID = 1
        values = [2.0]
      )",
      logger);

  // Check first oscillator has default settings
  const auto& control_seg0 = config.getControlSegments(0);
  EXPECT_EQ(control_seg0.size(), 1);
  EXPECT_EQ(control_seg0[0].type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg0[0].nspline.value(), 10);
  EXPECT_DOUBLE_EQ(control_seg0[0].tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg0[0].tstop.value(), config.getTotalTime());

  // Check second oscillator has given settings
  const auto& control_seg1 = config.getControlSegments(1);
  const auto& control_bounds1 = config.getControlBounds(1);
  EXPECT_EQ(control_seg1.size(), 1);
  EXPECT_EQ(control_seg1[0].type, ControlType::BSPLINE0);
  EXPECT_EQ(control_bounds1.size(), 1);
  EXPECT_DOUBLE_EQ(control_bounds1[0], 2.0);
  EXPECT_EQ(control_seg1[0].nspline.value(), 150);
  EXPECT_DOUBLE_EQ(control_seg1[0].tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg1[0].tstop.value(), 1.0);

  // Check third oscillator has default settings
  const auto& control_seg2 = config.getControlSegments(2);
  EXPECT_EQ(control_seg2.size(), 1);
  EXPECT_EQ(control_seg2[0].type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg2[0].nspline.value(), 10);
  EXPECT_DOUBLE_EQ(control_seg2[0].tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg2[0].tstop.value(), config.getTotalTime());
}

TEST_F(TomlParserTest, ControlSegments_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [[control_segments]]
        type = "spline0"
        num = 150
        tstart = 0.0
        tstop = 1.0
      )",
      logger);

  // Check first oscillator has given settings
  const auto& control_seg0 = config.getControlSegments(0);
  EXPECT_EQ(control_seg0.size(), 1);
  EXPECT_EQ(control_seg0[0].type, ControlType::BSPLINE0);
  EXPECT_EQ(control_seg0[0].nspline.value(), 150);
  EXPECT_DOUBLE_EQ(control_seg0[0].tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg0[0].tstop.value(), 1.0);

  // Check second oscillator has given settings
  const auto& control_seg1 = config.getControlSegments(1);
  EXPECT_EQ(control_seg1.size(), 1);
  EXPECT_EQ(control_seg1[0].type, ControlType::BSPLINE0);
  EXPECT_EQ(control_seg1[0].nspline.value(), 150);
  EXPECT_DOUBLE_EQ(control_seg1[0].tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg1[0].tstop.value(), 1.0);
}

TEST_F(TomlParserTest, ControlInitialization_Defaults) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2, 2]
        transfreq = [4.1, 4.1, 4.1]
        rotfreq = [0.0, 0.0, 0.0]
        initial_condition = {type = "basis"}

        [[control_initialization]]
        oscID = 1
        type = "random"
        amplitude = 2.0
      )",
      logger);

  // Check first oscillator has default settings
  const auto& control_init0 = config.getControlInitializations(0);
  EXPECT_EQ(control_init0.size(), 1);
  EXPECT_EQ(control_init0[0].type, ControlSegmentInitType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init0[0].amplitude, 0.0);
  EXPECT_DOUBLE_EQ(control_init0[0].phase, 0.0);

  // Check second oscillator has given settings
  const auto& control_init1 = config.getControlInitializations(1);
  EXPECT_EQ(control_init1.size(), 1);
  EXPECT_EQ(control_init1[0].type, ControlSegmentInitType::RANDOM);
  EXPECT_DOUBLE_EQ(control_init1[0].amplitude, 2.0);
  EXPECT_DOUBLE_EQ(control_init1[0].phase, 0.0);

  // Check third oscillator has default settings
  const auto& control_init2 = config.getControlInitializations(2);
  EXPECT_EQ(control_init2.size(), 1);
  EXPECT_EQ(control_init2[0].type, ControlSegmentInitType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init2[0].amplitude, 0.0);
  EXPECT_DOUBLE_EQ(control_init2[0].phase, 0.0);
}

TEST_F(TomlParserTest, ControlInitialization) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2, 2, 2, 2]
        transfreq = [4.1, 4.1, 4.1, 4.1, 4.1]
        rotfreq = [0.0, 0.0, 0.0, 0.0, 0.0]
        initial_condition = {type = "basis"}

        [[control_initialization]]
        oscID = 0
        type = "constant"
        amplitude = 1.0
        phase = 1.1
        [[control_initialization]]
        oscID = 1
        type = "constant"
        amplitude = 2.0
        [[control_initialization]]
        oscID = 2
        type = "random"
        amplitude = 3.0
        phase = 3.1
        [[control_initialization]]
        oscID = 3
        type = "random"
        amplitude = 4.0
        [[control_initialization]]
        oscID = 4
        type = "random"
        amplitude = 5.0
        phase = 5.1
        [[control_initialization]]
        oscID = 4
        type = "constant"
        amplitude = 6.0
        phase = 6.1
      )",
      logger);

  // Check first oscillator
  const auto& control_init0 = config.getControlInitializations(0);
  EXPECT_EQ(control_init0.size(), 1);
  EXPECT_EQ(control_init0[0].type, ControlSegmentInitType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init0[0].amplitude, 1.0);
  EXPECT_DOUBLE_EQ(control_init0[0].phase, 1.1);

  // Check second oscillator
  const auto& control_init1 = config.getControlInitializations(1);
  EXPECT_EQ(control_init1.size(), 1);
  EXPECT_EQ(control_init1[0].type, ControlSegmentInitType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init1[0].amplitude, 2.0);
  EXPECT_DOUBLE_EQ(control_init1[0].phase, 0.0);

  // Check third oscillator
  const auto& control_init2 = config.getControlInitializations(2);
  EXPECT_EQ(control_init2.size(), 1);
  EXPECT_EQ(control_init2[0].type, ControlSegmentInitType::RANDOM);
  EXPECT_DOUBLE_EQ(control_init2[0].amplitude, 3.0);
  EXPECT_DOUBLE_EQ(control_init2[0].phase, 3.1);

  // Check fourth oscillator
  const auto& control_init3 = config.getControlInitializations(3);
  EXPECT_EQ(control_init3.size(), 1);
  EXPECT_EQ(control_init3[0].type, ControlSegmentInitType::RANDOM);
  EXPECT_DOUBLE_EQ(control_init3[0].amplitude, 4.0);
  EXPECT_DOUBLE_EQ(control_init3[0].phase, 0.0);

  // Check fifth oscillator with two segments
  const auto& control_init4 = config.getControlInitializations(4);
  EXPECT_EQ(control_init4.size(), 2);
  EXPECT_EQ(control_init4[0].type, ControlSegmentInitType::RANDOM);
  EXPECT_DOUBLE_EQ(control_init4[0].amplitude, 5.0);
  EXPECT_DOUBLE_EQ(control_init4[0].phase, 5.1);
  EXPECT_EQ(control_init4[1].type, ControlSegmentInitType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init4[1].amplitude, 6.0);
  EXPECT_DOUBLE_EQ(control_init4[1].phase, 6.1);
}

TEST_F(TomlParserTest, ControlInitialization_File) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}

        [[control_initialization]]
        type = "file"
        filename = "params.dat"
      )",
      logger);

  EXPECT_TRUE(config.getControlInitializationFile().has_value());
  EXPECT_EQ(config.getControlInitializationFile().value(), "params.dat");
}

TEST_F(TomlParserTest, ControlInitialization_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 4.1]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [[control_initialization]]
        type = "constant"
        amplitude = 1.0
        phase = 1.1
      )",
      logger);

  // Check first oscillator
  const auto& control_init0 = config.getControlInitializations(0);
  EXPECT_EQ(control_init0.size(), 1);
  EXPECT_EQ(control_init0[0].type, ControlSegmentInitType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init0[0].amplitude, 1.0);
  EXPECT_DOUBLE_EQ(control_init0[0].phase, 1.1);

  // Check second oscillator
  const auto& control_init1 = config.getControlInitializations(1);
  EXPECT_EQ(control_init1.size(), 1);
  EXPECT_EQ(control_init1[0].type, ControlSegmentInitType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init1[0].amplitude, 1.0);
  EXPECT_DOUBLE_EQ(control_init1[0].phase, 1.1);
}

TEST_F(TomlParserTest, ControlBounds) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}

        [[control_segments]]
        oscID = 0
        type = "spline"
        num = 10
        tstart = 0.0
        tstop = 1.0
        [[control_segments]]
        oscID = 0
        type = "step"
        step_amp1 = 0.1
        step_amp2 = 0.2
        tramp = 0.3
        tstart = 0.4
        tstop = 0.5
        [[control_segments]]
        oscID = 0
        num = 20
        type = "spline0"
        tstart = 1.0
        tstop = 2.0
        [[control_bounds]]
        oscID = 0
        values = [1.0, 2.0]
      )",
      logger);

  // Check control bounds for the three segments
  const auto& control_bounds0 = config.getControlBounds(0);
  EXPECT_EQ(control_bounds0.size(), 3);
  EXPECT_EQ(control_bounds0[0], 1.0);
  EXPECT_EQ(control_bounds0[1], 2.0);
  EXPECT_EQ(control_bounds0[2], 2.0); // Use last bound for extra segments
}

TEST_F(TomlParserTest, ControlBounds_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 4.1]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [[control_segments]]
        type = "spline"
        num = 10
        tstart = 0.0
        tstop = 1.0
        [[control_segments]]
        type = "step"
        step_amp1 = 0.1
        step_amp2 = 0.2
        tramp = 0.3
        tstart = 0.4
        tstop = 0.5
        [[control_segments]]
        num = 20
        type = "spline0"
        tstart = 1.0
        tstop = 2.0
        [[control_bounds]]
        values = [1.0, 2.0]
      )",
      logger);

  // Check control bounds for the three segments
  const auto& control_bounds0 = config.getControlBounds(0);
  EXPECT_EQ(control_bounds0.size(), 3);
  EXPECT_EQ(control_bounds0[0], 1.0);
  EXPECT_EQ(control_bounds0[1], 2.0);
  EXPECT_EQ(control_bounds0[2], 2.0); // Use last bound for extra segments

  const auto& control_bounds1 = config.getControlBounds(1);
  EXPECT_EQ(control_bounds1.size(), 3);
  EXPECT_EQ(control_bounds1[0], 1.0);
  EXPECT_EQ(control_bounds1[1], 2.0);
  EXPECT_EQ(control_bounds1[2], 2.0); // Use last bound for extra segments
}

TEST_F(TomlParserTest, CarrierFrequencies) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}

        [[carrier_frequency]]
        oscID = 0
        values = [1.0, 2.0]
      )",
      logger);

  const auto& carrier_freq0 = config.getCarrierFrequencies(0);
  EXPECT_EQ(carrier_freq0.size(), 2);
  EXPECT_EQ(carrier_freq0[0], 1.0);
  EXPECT_EQ(carrier_freq0[1], 2.0);
}

TEST_F(TomlParserTest, CarrierFrequencies_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 4.1]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}

        [[carrier_frequency]]
        values = [1.0, 2.0]
      )",
      logger);

  const auto& carrier_freq0 = config.getCarrierFrequencies(0);
  EXPECT_EQ(carrier_freq0.size(), 2);
  EXPECT_EQ(carrier_freq0[0], 1.0);
  EXPECT_EQ(carrier_freq0[1], 2.0);

  const auto& carrier_freq1 = config.getCarrierFrequencies(1);
  EXPECT_EQ(carrier_freq1.size(), 2);
  EXPECT_EQ(carrier_freq1[0], 1.0);
  EXPECT_EQ(carrier_freq1[1], 2.0);
}

TEST_F(TomlParserTest, OptimTarget_GateType) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}
        optim_target = {target_type = "gate", gate_type = "cnot"}
      )",
      logger);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::GATE);
  EXPECT_EQ(target.gate_type.value(), GateType::CNOT);
}

TEST_F(TomlParserTest, OptimTarget_GateFromFile) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}
        optim_target = {target_type = "gate", gate_type = "file", gate_file = "/path/to/gate.dat"}
      )",
      logger);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::GATE);
  EXPECT_EQ(target.gate_type.value(), GateType::FILE);
  EXPECT_EQ(target.gate_file.value(), "/path/to/gate.dat");
}

TEST_F(TomlParserTest, OptimTarget_PureState) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [3, 3, 3]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}
        optim_target = {target_type = "pure", levels = [0,1,2]}
      )",
      logger);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::PURE);
  const auto& levels = target.levels.value();
  EXPECT_EQ(levels.size(), 3);
  EXPECT_EQ(levels[0], 0);
  EXPECT_EQ(levels[1], 1);
  EXPECT_EQ(levels[2], 2);
}

TEST_F(TomlParserTest, OptimTarget_FromFile) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}
        optim_target = {target_type = "file", filename = "/path/to/target.dat"}
      )",
      logger);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::FROMFILE);
  EXPECT_EQ(target.file.value(), "/path/to/target.dat");
}

TEST_F(TomlParserTest, OptimTarget_DefaultPure) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}
      )",
      logger);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::PURE);
  // For default pure state, levels should be set to ground state
  EXPECT_TRUE(target.levels.has_value());
  EXPECT_FALSE(target.levels.value().empty());
}

TEST_F(TomlParserTest, OptimWeights) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 4.1]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}
        optim_weights = [2.0, 1.0]
      )",
      logger);

  const auto& weights = config.getOptimWeights();
  EXPECT_EQ(weights.size(), 4);
  EXPECT_DOUBLE_EQ(weights[0], 0.4);
  EXPECT_DOUBLE_EQ(weights[1], 0.2);
  EXPECT_DOUBLE_EQ(weights[2], 0.2);
  EXPECT_DOUBLE_EQ(weights[3], 0.2);
}

TEST_F(TomlParserTest, ComprehensiveNonDefaultSettings_PrintConfigValidation) {
  // Create a comprehensive configuration with all non-default settings (single oscillator)
  std::string input_toml = R"(# Configuration settings
# =============================================

nlevels = [3]
nessential = [2]
ntime = 2500
dt = 0.02
transfreq = [5.500000]
selfkerr = [-0.200000]
crosskerr = [0.050000]
Jkl = [0.010000]
rotfreq = [2.100000]
collapse_type = "decay"
decay_time = [25.000000]
dephase_time = [45.000000]
initial_condition = {type = "pure", levels = [1]}
control_enforceBC = false
optim_target = {target_type = "gate", gate_type = "hadamard"}
gate_rot_freq = [1.500000]
optim_objective = "jtrace"
optim_weights = [1.000000]
optim_atol = 1e-06
optim_rtol = 0.001
optim_ftol = 1e-07
optim_inftol = 0.0001
optim_maxiter = 150
optim_regul = 0.002
optim_penalty = 0.1
optim_penalty_param = 0.8
optim_penalty_dpdm = 0.05
optim_penalty_energy = 0.02
optim_penalty_variation = 0.03
optim_regul_tik0 = true
datadir = "/custom/output/path"
output_frequency = 5
optim_monitor_frequency = 25
runtype = "optimization"
usematfree = true
linearsolver_type = "gmres"
linearsolver_maxiter = 20
timestepper = "imr"
rand_seed = 12345

[[apply_pipulse]]
oscID = 0
tstart = 1
tstop = 2
amp = 0.75

[[control_segments]]
oscID = 0
type = "spline"
num = 25
tstart = 0.5
tstop = 1.5

[[control_initialization]]
oscID = 0
type = "random"
amplitude = 1.500000
phase = 0.500000

[[control_bounds]]
oscID = 0
values = [5.000000, 8.000000]

[[carrier_frequency]]
oscID = 0
values = [2.500000, 3.000000]

[[write]]
oscID = 0
type = ["population", "expectedenergy"]

)";

  Config config = Config::fromTomlString(input_toml, logger);

  // Call printConfig and capture the output
  std::stringstream printed_output;
  config.printConfig(printed_output);

  std::string output = printed_output.str();

  // Split both strings into lines
  auto splitIntoLines = [](const std::string& text) {
    std::vector<std::string> lines;
    std::istringstream stream(text);
    std::string line;
    while (std::getline(stream, line)) {
      lines.push_back(line);
    }
    return lines;
  };

  std::vector<std::string> input_lines = splitIntoLines(input_toml);
  std::vector<std::string> output_lines = splitIntoLines(output);

  // Compare line counts
  EXPECT_EQ(input_lines.size(), output_lines.size()) 
    << "Line count mismatch: input has " << input_lines.size() 
    << " lines, output has " << output_lines.size() << " lines";

  // Compare each line
  size_t max_lines = std::max(input_lines.size(), output_lines.size());
  for (size_t i = 0; i < max_lines; ++i) {
    if (i >= input_lines.size()) {
      FAIL() << "Output has extra line " << (i + 1) << ": '" << output_lines[i] << "'";
    } else if (i >= output_lines.size()) {
      FAIL() << "Output missing line " << (i + 1) << ": '" << input_lines[i] << "'";
    } else {
      EXPECT_EQ(input_lines[i], output_lines[i]) << "Line " << (i + 1) << " differs";
    }
  }
}

TEST_F(TomlParserTest, ApplyPipulse_UnknownKey) {
  ASSERT_DEATH({
    Config config = Config::fromTomlString(
        R"(
          nlevels = [2]
          transfreq = [4.1]
          rotfreq = [0.0]
          initial_condition = {type = "basis"}

          [[apply_pipulse]]
          oscID = 0
          tstart = 0.0
          tstop = 1.0
          amp = 0.5
          invalid_param = 123
        )",
        logger);
  }, "ERROR: Unknown key 'invalid_param' in apply_pipulse\\.");
}

TEST_F(TomlParserTest, ControlSegments_UnknownKey) {
  ASSERT_DEATH({
    Config config = Config::fromTomlString(
        R"(
          nlevels = [2]
          transfreq = [4.1]
          rotfreq = [0.0]
          initial_condition = {type = "basis"}

          [[control_segments]]
          type = "spline"
          num = 5
          wrong_field = 99
        )",
        logger);
  }, "ERROR: Unknown key 'wrong_field' in control_segments\\.");
}

TEST_F(TomlParserTest, ControlInitialization_UnknownKey) {
  ASSERT_DEATH({
    Config config = Config::fromTomlString(
        R"(
          nlevels = [2]
          transfreq = [4.1]
          rotfreq = [0.0]
          initial_condition = {type = "basis"}

          [[control_initialization]]
          type = "file"
          filename = "params.dat"
          foo = 42
        )",
        logger);
  }, "ERROR: Unknown key 'foo' in control_initialization\\.");
}

TEST_F(TomlParserTest, ControlBounds_UnknownKey) {
  ASSERT_DEATH({
    Config config = Config::fromTomlString(
        R"(
          nlevels = [2]
          transfreq = [4.1]
          rotfreq = [0.0]
          initial_condition = {type = "basis"}

          [[control_bounds]]
          oscID = 0
          values = [1.0, 2.0]
          extra_key = "not_allowed"
        )",
        logger);
  }, "ERROR: Unknown key 'extra_key' in control_bounds\\.");
}

TEST_F(TomlParserTest, CarrierFrequency_UnknownKey) {
  ASSERT_DEATH({
    Config config = Config::fromTomlString(
        R"(
          nlevels = [2]
          transfreq = [4.1]
          rotfreq = [0.0]
          initial_condition = {type = "basis"}

          [[carrier_frequency]]
          oscID = 0
          values = [4.0]
          forbidden_key = 42
        )",
        logger);
  }, "ERROR: Unknown key 'forbidden_key' in carrier_frequency\\.");
}

TEST_F(TomlParserTest, Write_UnknownKey) {
  ASSERT_DEATH({
    Config config = Config::fromTomlString(
        R"(
          nlevels = [2]
          transfreq = [4.1]
          rotfreq = [0.0]
          initial_condition = {type = "basis"}

          [[write]]
          oscID = 0
          type = ["expectedenergy"]
          unsupported_option = true
        )",
        logger);
  }, "ERROR: Unknown key 'unsupported_option' in write\\.");
}
