#include <gtest/gtest.h>

#include <cstdio>

#include "config.hpp"

class CfgParserTest : public ::testing::Test {};

TEST_F(CfgParserTest, ParseBasicSettings) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        transfreq = 4.1
        rotfreq = 0.0
        ntime = 500
        dt = 0.05
        collapse_type = none
        initialcondition = basis
      )",
      false);

  EXPECT_EQ(config.getNTime(), 500);
  EXPECT_DOUBLE_EQ(config.getDt(), 0.05);
  EXPECT_EQ(config.getDecoherenceType(), DecoherenceType::NONE);
}

TEST_F(CfgParserTest, ParseVectorSettings) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2, 3
        ntime = 1000
        dt = 0.1
        transfreq = 4.1, 4.8, 5.2
        rotfreq = 0.0, 0.0
        initialcondition = basis
      )",
      false);

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

TEST_F(CfgParserTest, ParseOutputSettings) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2, 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1, 4.8
        rotfreq = 0.0, 0.0
        initialcondition = basis
        output0 = population
        output1 = population, expectedEnergy
      )",
      false);

  // Verify output settings
  auto output = config.getOutputObservables();
  EXPECT_EQ(output.size(), 2); // 2: population, expectedEnergy 
  EXPECT_EQ(output[0], OutputType::EXPECTED_ENERGY);
  EXPECT_EQ(output[1], OutputType::POPULATION);
}

TEST_F(CfgParserTest, ParseStructSettings) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        transfreq = 4.1
        rotfreq = 0.0
        ntime = 1000
        dt = 0.1
        optim_target = gate, cnot
        initialcondition = diagonal, 0
      )",
      false);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::GATE);
  EXPECT_EQ(target.gate_type.value(), GateType::CNOT);

  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::DIAGONAL);
  EXPECT_EQ(initcond.subsystem.value(), std::vector<size_t>{0});
}

TEST_F(CfgParserTest, ApplyDefaults) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1
        initialcondition = basis
      )",
      false);

  EXPECT_EQ(config.getDecoherenceType(), DecoherenceType::NONE); // Default
  EXPECT_EQ(config.getRotFreq()[0], 0.0); // Default
}

TEST_F(CfgParserTest, InitialCondition_FromFile) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        transfreq = 4.1
        ntime = 1000
        dt = 0.1
        rotfreq = 0.0
        initialcondition = file, test.dat
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::FROMFILE);
  EXPECT_EQ(initcond.filename.value(), "test.dat");
  EXPECT_EQ(config.getNInitialConditions(), 1);
}

TEST_F(CfgParserTest, InitialCondition_Pure) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 3, 2
        transfreq = 4.1, 4.8
        ntime = 1000
        dt = 0.1
        rotfreq = 0.0, 0.0
        initialcondition = pure, 1, 0
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  // Test backward compatibility: "pure" in CFG should map to PRODUCT_STATE
  EXPECT_EQ(initcond.type, InitialConditionType::PRODUCT_STATE);
  EXPECT_EQ(initcond.levels.value(), std::vector<size_t>({1, 0}));
  EXPECT_EQ(config.getNInitialConditions(), 1);
}

TEST_F(CfgParserTest, InitialCondition_Performance) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        transfreq = 4.1
        ntime = 1000
        dt = 0.1
        rotfreq = 0.0
        initialcondition = performance
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::PERFORMANCE);
  EXPECT_EQ(config.getNInitialConditions(), 1);
}

TEST_F(CfgParserTest, InitialCondition_Ensemble) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 3, 2
        transfreq = 4.1, 4.8
        ntime = 1000
        dt = 0.1
        rotfreq = 0.0, 0.0
        collapse_type = decay
        initialcondition = ensemble, 0, 1
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::ENSEMBLE);
  EXPECT_EQ(initcond.subsystem.value(), std::vector<size_t>({0, 1}));
  EXPECT_EQ(config.getNInitialConditions(), 1);
}

TEST_F(CfgParserTest, InitialCondition_ThreeStates) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 3
        transfreq = 4.1
        ntime = 1000
        dt = 0.1
        rotfreq = 0.0
        collapse_type = decay
        initialcondition = 3states
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::THREESTATES);
  EXPECT_EQ(config.getNInitialConditions(), 3);
}

TEST_F(CfgParserTest, InitialCondition_NPlusOne_SingleOscillator) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 3
        transfreq = 4.1
        rotfreq = 0.0
        ntime = 1000
        dt = 0.1
        collapse_type = decay
        initialcondition = nplus1
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::NPLUSONE);
  // For nlevels = [3], system dimension N = 3, so n_initial_conditions = N + 1 = 4
  EXPECT_EQ(config.getNInitialConditions(), 4);
}

TEST_F(CfgParserTest, InitialCondition_NPlusOne_MultipleOscillators) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2, 3
        transfreq = 4.1, 4.8
        ntime = 1000
        dt = 0.1
        rotfreq = 0.0, 0.0
        collapse_type = decay
        initialcondition = nplus1
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::NPLUSONE);
  // For nlevels = [2, 3], system dimension N = 2 * 3 = 6, so n_initial_conditions = N + 1 = 7
  EXPECT_EQ(config.getNInitialConditions(), 7);
}

TEST_F(CfgParserTest, InitialCondition_Diagonal_Schrodinger) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 3, 2
        ntime = 1000
        dt = 0.1
        nessential = 3, 2
        transfreq = 4.1, 4.8
        rotfreq = 0.0, 0.0
        collapse_type = none
        initialcondition = diagonal, 1
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::DIAGONAL);
  EXPECT_EQ(initcond.subsystem.value(), std::vector<size_t>({1}));
  // For Schrodinger solver (decoherence_type = none), n_initial_conditions = nessential[1] = 2
  EXPECT_EQ(config.getNInitialConditions(), 2);
}

TEST_F(CfgParserTest, InitialCondition_Basis_Schrodinger) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 3, 2
        ntime = 1000
        dt = 0.1
        nessential = 3, 2
        transfreq = 4.1, 4.8
        rotfreq = 0.0, 0.0
        collapse_type = none
        initialcondition = basis, 1
      )",
      false);
  // For Schrodinger solver, BASIS is converted to DIAGONAL, so n_initial_conditions = nessential[1] = 2
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::DIAGONAL);
  EXPECT_EQ(initcond.subsystem.value(), std::vector<size_t>({1}));
  EXPECT_EQ(config.getNInitialConditions(), 2);
}

TEST_F(CfgParserTest, InitialCondition_Basis_Lindblad) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 3, 2
        ntime = 1000
        dt = 0.1
        nessential = 3, 2
        transfreq = 4.1, 4.8
        rotfreq = 0.0, 0.0
        collapse_type = decay
        initialcondition = basis, 1
      )",
      false);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::BASIS);
  EXPECT_EQ(initcond.subsystem.value(), std::vector<size_t>({1}));
  // For Lindblad solver, n_initial_conditions = nessential[1]^2 = 2^2 = 4
  EXPECT_EQ(config.getNInitialConditions(), 4);
}

TEST_F(CfgParserTest, ControlParameterizations_Spline0) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1
        rotfreq = 0.0
        initialcondition = basis
        control_segments0 = spline0, 150, 0.0, 1.0
      )",
      false);

  const auto& control_seg0 = config.getControlParameterizations(0);
  EXPECT_EQ(control_seg0.type, ControlType::BSPLINE0);
  EXPECT_EQ(control_seg0.nspline.value(), 150);
  EXPECT_DOUBLE_EQ(control_seg0.tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg0.tstop.value(), 1.0);
}

TEST_F(CfgParserTest, ControlParameterizations_Spline) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2, 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1, 4.1
        rotfreq = 0.0, 0.0
        initialcondition = basis
        control_segments0 = spline, 10
        control_segments1 = spline, 20, 0.0, 1.0
      )",
      false);

  // Check first oscillator 
  const auto& control_seg0 = config.getControlParameterizations(0);
  EXPECT_EQ(control_seg0.type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg0.nspline.value(), 10);
  EXPECT_FALSE(control_seg0.tstart.has_value());
  EXPECT_FALSE(control_seg0.tstop.has_value());

  // Check second oscillator 
  const auto& control_seg1 = config.getControlParameterizations(1);
  EXPECT_EQ(control_seg1.type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg1.nspline.value(), 20);
  EXPECT_DOUBLE_EQ(control_seg1.tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg1.tstop.value(), 1.0);
}

TEST_F(CfgParserTest, ControlParameterizations_Defaults) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2, 2, 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1, 4.8
        rotfreq = 0.0, 0.0
        initialcondition = basis
        control_segments1 = spline0, 150, 0.0, 1.0
        control_bounds1 = 2.0
      )",
      false);

  // Check first oscillator has default settings
  const auto& control_seg0 = config.getControlParameterizations(0);
  EXPECT_EQ(control_seg0.type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg0.nspline.value(), 10);
  EXPECT_FALSE(control_seg0.tstart.has_value());
  EXPECT_FALSE(control_seg0.tstop.has_value());

  // Check second oscillator has given settings
  const auto& control_seg1 = config.getControlParameterizations(1);
  const double control_bound1 = config.getControlAmplitudeBound(1);
  EXPECT_EQ(control_seg1.type, ControlType::BSPLINE0);
  EXPECT_DOUBLE_EQ(control_bound1, 2.0);
  EXPECT_EQ(control_seg1.nspline.value(), 150);
  EXPECT_DOUBLE_EQ(control_seg1.tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg1.tstop.value(), 1.0);

  // Check third oscillator had default settings
  const auto& control_seg2 = config.getControlParameterizations(2);
  EXPECT_EQ(control_seg2.type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg2.nspline.value(), 10);
  EXPECT_FALSE(control_seg2.tstart.has_value());
  EXPECT_FALSE(control_seg2.tstop.has_value());
}

TEST_F(CfgParserTest, ControlInitialization_Defaults) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2, 2, 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1, 4.1, 4.1
        rotfreq = 0.0, 0.0, 0.0
        initialcondition = basis
        control_initialization1 = random, 2.0
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

  // Check third oscillator defaults to the second's settings
  const auto& control_init2 = config.getControlInitializations(2);
  EXPECT_EQ(control_init2.type, ControlInitializationType::CONSTANT);
  EXPECT_DOUBLE_EQ(control_init2.amplitude.value(), 0.0);
}

TEST_F(CfgParserTest, ControlInitializationSettings) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2, 2, 2, 2, 2
        transfreq = 4.1, 4.1, 4.1, 4.1, 4.1
        ntime = 1000
        dt = 0.1
        rotfreq = 0.0, 0.0, 0.0, 0.0, 0.0
        initialcondition = basis
        control_segments0 = spline_amplitude, 10, 1.0
        control_initialization0 = constant, 1.0, 1.1
        control_segments1 = spline, 10
        control_initialization1 = constant, 2.0
        control_segments2 = spline_amplitude, 10, 1.0
        control_initialization2 = random, 3.0, 3.1
        control_segments3 = spline, 10
        control_initialization3 = random, 4.0
        control_segments4 = spline_amplitude, 10, 1.0
        control_initialization4 = random, 5.0, 5.1
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

  // Check fifth oscillator with two parameterizations
  const auto& control_init4 = config.getControlInitializations(4);
  EXPECT_EQ(control_init4.type, ControlInitializationType::RANDOM);
  EXPECT_DOUBLE_EQ(control_init4.amplitude.value(), 5.0);
  EXPECT_DOUBLE_EQ(control_init4.phase.value(), 5.1);
}

TEST_F(CfgParserTest, ControlInitialization_File) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1
        rotfreq = 0.0
        initialcondition = basis
        control_initialization0 = file, myparams.dat
      )",
      false);

  const auto& control_init0 = config.getControlInitializations(0);
  EXPECT_EQ(control_init0.type, ControlInitializationType::FILE);
  EXPECT_EQ(control_init0.filename.value(), "myparams.dat");
}

TEST_F(CfgParserTest, ControlBounds) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1
        rotfreq = 0.0
        initialcondition = basis
        control_segments0 = spline, 10, 0.0, 1.0, step, 0.1, 0.2, 0.3, 0.4, 0.5, spline0, 20, 1.0, 2.0
        control_bounds0 = 1.5
      )",
      false);

  // Check control bound
  const double control_bound0 = config.getControlAmplitudeBound(0);
  EXPECT_DOUBLE_EQ(control_bound0, 1.5);
  // Check that only the first parameterization is active
  const auto& control_seg0 = config.getControlParameterizations(0);
  EXPECT_EQ(control_seg0.type, ControlType::BSPLINE);
  EXPECT_EQ(control_seg0.nspline.value(), 10);
  EXPECT_DOUBLE_EQ(control_seg0.tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_seg0.tstop.value(), 1.0);
}

TEST_F(CfgParserTest, CarrierFrequencies) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1
        rotfreq = 0.0
        initialcondition = basis
        carrier_frequency0 = 1.0, 2.0
      )",
      false);

  const auto& carrier_freq0 = config.getCarrierFrequencies(0);
  EXPECT_EQ(carrier_freq0.size(), 2);
  EXPECT_EQ(carrier_freq0[0], 1.0);
  EXPECT_EQ(carrier_freq0[1], 2.0);
}

TEST_F(CfgParserTest, OptimTarget_GateType) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1
        rotfreq = 0.0
        initialcondition = basis
        optim_target = gate, cnot
      )",
      false);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::GATE);
  EXPECT_EQ(target.gate_type.value(), GateType::CNOT);
}

TEST_F(CfgParserTest, OptimTarget_GateFromFile) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1
        rotfreq = 0.0
        initialcondition = basis
        optim_target = gate, file, /path/to/gate.dat
      )",
      false);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::GATE);
  EXPECT_EQ(target.gate_type.value(), GateType::FILE);
  EXPECT_EQ(target.filename.value(), "/path/to/gate.dat");
}

TEST_F(CfgParserTest, OptimTarget_PureState) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 3, 3, 3
        ntime = 1000
        dt = 0.1
        transfreq = 4.1
        rotfreq = 0.0
        initialcondition = basis
        optim_target = pure, 0, 1, 2
      )",
      false);

  const auto& target = config.getOptimTarget();
  // Test backward compatibility: "pure" in CFG should map to STATE
  EXPECT_EQ(target.type, TargetType::STATE);
  EXPECT_TRUE(target.levels.has_value());
  const auto& levels = target.levels.value();
  EXPECT_EQ(levels.size(), 3);
  EXPECT_EQ(levels[0], 0);
  EXPECT_EQ(levels[1], 1);
  EXPECT_EQ(levels[2], 2);
}

TEST_F(CfgParserTest, OptimTarget_FromFile) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1
        rotfreq = 0.0
        initialcondition = basis
        optim_target = file, /path/to/targetstate.dat
      )",
      false);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::STATE);
  EXPECT_TRUE(target.filename.has_value());
  EXPECT_EQ(target.filename.value(), "/path/to/targetstate.dat");
}

TEST_F(CfgParserTest, OptimTarget_DefaultNone) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1
        rotfreq = 0.0
        initialcondition = basis
      )",
      false);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::NONE);
}

TEST_F(CfgParserTest, OptimWeights) {
  Config config = Config::fromCfgString(
      R"(
        nlevels = 2, 2
        ntime = 1000
        dt = 0.1
        transfreq = 4.1, 4.1
        rotfreq = 0.0, 0.0
        initialcondition = basis
        optim_weights = 2.0
      )",
      false);

  const auto& weights = config.getOptimWeights();
  EXPECT_EQ(weights.size(), 4);
  EXPECT_DOUBLE_EQ(weights[0], 0.25);
  EXPECT_DOUBLE_EQ(weights[1], 0.25);
  EXPECT_DOUBLE_EQ(weights[2], 0.25);
  EXPECT_DOUBLE_EQ(weights[3], 0.25);
}
