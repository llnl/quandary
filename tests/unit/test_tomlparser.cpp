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
        decoherence = {type = "none"}
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

TEST_F(TomlParserTest, ParseJklAndCrossKerrSettings) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2, 2, 2]
        transfreq = [4.1, 4.8, 5.2]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}
        Jkl = { "0-1" = 0.15, "1-2" = 0.25, "2-3" = 0.35 }
        crosskerr = { "0-1" = 0.1, "2-3" = 0.05 }
      )",
      logger);

  auto Jkl = config.getJkl();
  size_t num_osc = config.getNumOsc();
  size_t num_pairs_osc = (num_osc - 1) * num_osc / 2;
  EXPECT_EQ(Jkl.size(), num_pairs_osc);
  EXPECT_DOUBLE_EQ(Jkl[0], 0.15); // 0-1
  EXPECT_DOUBLE_EQ(Jkl[1], 0.0);  // 0-2
  EXPECT_DOUBLE_EQ(Jkl[2], 0.0);  // 0-3
  EXPECT_DOUBLE_EQ(Jkl[3], 0.25); // 1-2
  EXPECT_DOUBLE_EQ(Jkl[4], 0.0);  // 1-3
  EXPECT_DOUBLE_EQ(Jkl[5], 0.35); // 2-3
  auto crosskerr = config.getCrossKerr();
  EXPECT_EQ(crosskerr.size(), num_pairs_osc);
  EXPECT_DOUBLE_EQ(crosskerr[0], 0.1); // 0-1
  EXPECT_DOUBLE_EQ(crosskerr[1], 0.0); // 0-2
  EXPECT_DOUBLE_EQ(crosskerr[2], 0.0); // 0-3
  EXPECT_DOUBLE_EQ(crosskerr[3], 0.0); // 1-2
  EXPECT_DOUBLE_EQ(crosskerr[4], 0.0); // 1-3
  EXPECT_DOUBLE_EQ(crosskerr[5], 0.05); // 2-3
}

TEST_F(TomlParserTest, ParseOutputSettings_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}
        output_type = ["population", "fullstate"]
      )",
      logger);

  // Verify output settings
  auto output = config.getOutputType();
  EXPECT_EQ(output.size(), 2); // Population
  EXPECT_EQ(output[0], OutputType::POPULATION);
  EXPECT_EQ(output[1], OutputType::FULLSTATE);
}

TEST_F(TomlParserTest, ParseStructSettings) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "diagonal", oscIDs = [0]}
        optim_target = {type = "gate", gate_type = "cnot"}
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

TEST_F(TomlParserTest, InitialCondition_ProductState) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [3, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "product_state", levels = [1, 0]}
      )",
      logger);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::PRODUCT_STATE);
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
        decoherence = {type = "decay"}
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
        decoherence = {type = "decay"}
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
        decoherence = {type = "decay"}
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
        decoherence = {type = "decay"}
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
        decoherence = {type = "none"}
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
        decoherence = {type = "none"}
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
        decoherence = {type = "decay"}
        initial_condition = {type = "basis", oscIDs = [1]}
      )",
      logger);
  const auto& initcond = config.getInitialCondition();
  EXPECT_EQ(initcond.type, InitialConditionType::BASIS);
  EXPECT_EQ(initcond.osc_IDs.value(), std::vector<size_t>({1}));
  // For Lindblad solver, n_initial_conditions = nessential[1]^2 = 2^2 = 4
  EXPECT_EQ(config.getNInitialConditions(), 4);
}

TEST_F(TomlParserTest, ControlParameterizations_Spline0) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}

        control_parameterization = { type = "spline0", num = 150, tstart = 0.0, tstop = 1.0 }
      )",
      logger);

  const auto& control_parameterizations = config.getControlParameterizations(0);
  EXPECT_EQ(control_parameterizations.type, ControlType::BSPLINE0);
  EXPECT_EQ(control_parameterizations.nspline.value(), 150);
  EXPECT_DOUBLE_EQ(control_parameterizations.tstart.value(), 0.0);
  EXPECT_DOUBLE_EQ(control_parameterizations.tstop.value(), 1.0);
}

TEST_F(TomlParserTest, ControlParameterizations_Spline) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 4.1]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}

        control_parameterization = {
          "0" = { type = "spline", num = 10 },
          "1" = { type = "spline", num = 20, tstart = 0.0, tstop = 1.0 }
        }
      )",
      logger);

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
        nlevels = [2, 2, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}

        control_parameterization = {
          "1" = { type = "spline0", num = 150, tstart = 0.0, tstop = 1.0 }
        }
      )",
      logger);

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
        nlevels = [2, 2]
        transfreq = [4.1, 4.8]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}

        control_parameterization = { type = "spline0", num = 150, tstart = 0.0, tstop = 1.0 }
      )",
      logger);

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
        nlevels = [2, 2, 2]
        transfreq = [4.1, 4.1, 4.1]
        rotfreq = [0.0, 0.0, 0.0]
        initial_condition = {type = "basis"}

        control_initialization = { 
          "1" = { type = "random", amplitude = 2.0 }
        }
      )",
      logger);

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

TEST_F(TomlParserTest, ControlInitialization) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2, 2, 2, 2]
        transfreq = [4.1, 4.1, 4.1, 4.1, 4.1]
        rotfreq = [0.0, 0.0, 0.0, 0.0, 0.0]
        initial_condition = {type = "basis"}

        control_initialization = {
          "0" = { type = "constant", amplitude = 1.0, phase = 1.1 },
          "1" = { type = "constant", amplitude = 2.0 },
          "2" = { type = "random", amplitude = 3.0, phase = 3.1 },
          "3" = { type = "random", amplitude = 4.0 },
          "4" = { type = "random", amplitude = 5.0, phase = 5.1 }
        }

        control_parameterization = {
          "0" = { type = "spline_amplitude", num = 10, scaling = 1.0 },
          "2" = { type = "spline_amplitude", num = 10, scaling = 1.0 },
          "4" = { type = "spline_amplitude", num = 10, scaling = 1.0 }
        }
      )",
      logger);

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
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}

        control_initialization = { type = "file", filename = "myparams.dat" }
      )",
      logger);

  const auto& control_init0 = config.getControlInitializations(0);
  EXPECT_EQ(control_init0.type, ControlInitializationType::FILE);
  EXPECT_EQ(control_init0.filename.value(), "myparams.dat");
}

TEST_F(TomlParserTest, ControlInitialization_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 4.1]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}

        control_initialization = { type = "constant", amplitude = 1.0 }
      )",
      logger);

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
        nlevels = [3, 3, 3]
        transfreq = [4.1, 4.1, 4.1]
        rotfreq = [0.0, 0.0, 0.0]
        initial_condition = {type = "basis"}

        control_initialization = { 
          "1" = { type = "random", amplitude = 0.05 }
        }
      )",
      logger);

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
        nlevels = [2, 2]
        transfreq = [4.1, 3.3]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}

        control_parameterization = { type = "spline", num = 10, tstart = 0.0, tstop = 1.0 }
        
        control_bounds = 1.5 
      )",
      logger);

  // Check control bound
  const double control_bound0 = config.getControlBound(0);
  const double control_bound1 = config.getControlBound(1);
  EXPECT_DOUBLE_EQ(control_bound0, 1.5);
  EXPECT_DOUBLE_EQ(control_bound1, 1.5);
}

TEST_F(TomlParserTest, ControlBounds_Defaults) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 3.0]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}
      )",
      logger);

  // Check control bound
  const double control_bound0 = config.getControlBound(0);
  const double control_bound1 = config.getControlBound(1);
  EXPECT_DOUBLE_EQ(control_bound0, 1e12);
  EXPECT_DOUBLE_EQ(control_bound1, 1e12);
}

TEST_F(TomlParserTest, ControlBounds_AllOscillators) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2, 2]
        transfreq = [4.1, 3.3]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}

        control_bounds = {
         "0" = 1.0,
         "1" = 2.0 
      }
      )",
      logger);

  // Check control bound
  const double control_bound0 = config.getControlBound(0);
  const double control_bound1 = config.getControlBound(1);
  EXPECT_DOUBLE_EQ(control_bound0, 1.0);
  EXPECT_DOUBLE_EQ(control_bound1, 2.0);
}

TEST_F(TomlParserTest, CarrierFrequencies) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}

        carrier_frequency = {"0" = [1.0, 2.0]}
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

        carrier_frequency = [1.0, 2.0]
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
        nlevels = [2,2,2]
        transfreq = [4.1, 3.4, 5.6]
        rotfreq = [0.0, 0.0, 0.0]
        initial_condition = {type = "basis"}
        optim_target = {type = "gate", gate_type = "cnot"}
      )",
      logger);

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
        nlevels = [2,2]
        transfreq = [4.1, 3.4]
        rotfreq = [0.0, 0.0]
        initial_condition = {type = "basis"}
        optim_target = {type = "gate", gate_type = "cnot", gate_rot_freq = [1.0, 2.0]}
      )",
      logger);

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
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}
        optim_target = {type = "gate", filename = "/path/to/gate.dat"}
      )",
      logger);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::GATE);
  EXPECT_EQ(target.gate_type.value(), GateType::FILE);
  EXPECT_EQ(target.filename.value(), "/path/to/gate.dat");
}

TEST_F(TomlParserTest, OptimTarget_ProductState) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [3, 3, 3]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}
        optim_target = {type = "state", levels = [0,1,2]}
      )",
      logger);

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
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}
        optim_target = {type = "state", filename = "/path/to/target.dat"}
      )",
      logger);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::STATE);
  EXPECT_TRUE(target.filename.has_value());
  EXPECT_EQ(target.filename.value(), "/path/to/target.dat");
}

TEST_F(TomlParserTest, OptimTarget_DefaultNone) {
  Config config = Config::fromTomlString(
      R"(
        nlevels = [2]
        transfreq = [4.1]
        rotfreq = [0.0]
        initial_condition = {type = "basis"}
      )",
      logger);

  const auto& target = config.getOptimTarget();
  EXPECT_EQ(target.type, TargetType::NONE);
  EXPECT_FALSE(target.gate_type.has_value());
  EXPECT_FALSE(target.levels.has_value());
  EXPECT_FALSE(target.filename.has_value());
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

// TEST_F(TomlParserTest, ComprehensiveNonDefaultSettings_PrintConfigValidation) {
//   // Create a comprehensive configuration with all non-default settings (single oscillator)
//   std::string input_toml = R"(# Configuration settings
// # =============================================

// nlevels = [3]
// nessential = [2]
// ntime = 2500
// dt = 0.02
// transfreq = [5.500000]
// selfkerr = [-0.200000]
// crosskerr = {}
// Jkl = {}
// rotfreq = [2.100000]
// decoherence = {type = "decay"}
// decay_time = [25.000000]
// dephase_time = [45.000000]
// initial_condition = {type = "product_state", levels = [1]}
// control_enforceBC = false
// optim_target = {type = "gate", gate_type = "hadamard"}
// gate_rot_freq = [1.500000]
// optim_objective = "jtrace"
// optim_weights = [1.000000]
// optim_atol = 1e-06
// optim_rtol = 0.001
// optim_ftol = 1e-07
// optim_inftol = 0.0001
// optim_maxiter = 150
// optim_regul = 0.002
// optim_penalty = 0.1
// optim_penalty_param = 0.8
// optim_penalty_dpdm = 0.05
// optim_penalty_energy = 0.02
// optim_penalty_variation = 0.03
// optim_regul_tik0 = true
// datadir = "/custom/output/path"
// output_frequency = 5
// optim_monitor_frequency = 25
// runtype = "optimization"
// usematfree = true
// linearsolver_type = "gmres"
// linearsolver_maxiter = 20
// timestepper = "imr"
// rand_seed = 12345
// output_type = ["population", "expectedenergy"]


// control_parameterization = { type = "spline", num = 25, tstart = 0.5, tstop = 1.5 }

// control_initialization = { type = "random", amplitude = 1.500000, phase = 0.500000 }

// [[control_bounds]]
// oscID = 0
// values = [5.000000, 8.000000]

// [[carrier_frequency]]
// oscID = 0
// values = [2.500000, 3.000000]


// )";

//   Config config = Config::fromTomlString(input_toml, logger);

//   // Test that printed config is valid TOML that can be parsed
//   std::stringstream printed_output;
//   config.printConfig(printed_output);
//   std::string output = printed_output.str();

//   // Verify output is valid TOML by parsing it
//   ASSERT_NO_THROW({
//     Config::fromTomlString(output, logger);
//   });
// }

TEST_F(TomlParserTest, ControlParameterizations_UnknownKey) {
  ASSERT_DEATH({
    Config config = Config::fromTomlString(
        R"(
          nlevels = [2]
          transfreq = [4.1]
          rotfreq = [0.0]
          initial_condition = {type = "basis"}

          control_parameterization = { type = "spline", num = 5, wrong_field = 99 }
        )",
        logger);
  }, "ERROR: Unknown key 'wrong_field' in control_parameterization\\.");
}

TEST_F(TomlParserTest, ControlInitialization_UnknownKey) {
  ASSERT_DEATH({
    Config config = Config::fromTomlString(
        R"(
          nlevels = [2]
          transfreq = [4.1]
          rotfreq = [0.0]
          initial_condition = {type = "basis"}

          control_initialization = { type = "file", filename = "params.dat", foo = 42 }
        )",
        logger);
  }, "ERROR: Unknown key 'foo' in control_initialization\\.");
}

TEST_F(TomlParserTest, ControlBounds_InvalidOscillatorID) {
  ASSERT_DEATH({
    Config config = Config::fromTomlString(
        R"(
          nlevels = [2]
          transfreq = [4.1]
          rotfreq = [0.0]
          initial_condition = {type = "basis"}

          control_bounds = { "5" = 1.0 }
        )",
        logger);
  }, "control_bounds oscillator ID 5 exceeds number of oscillators");
}

TEST_F(TomlParserTest, CarrierFrequency_InvalidOscillatorID) {
  ASSERT_DEATH({
    Config config = Config::fromTomlString(
        R"(
          nlevels = [2]
          transfreq = [4.1]
          rotfreq = [0.0]
          initial_condition = {type = "basis"}

          carrier_frequency = {"5" = [4.0]}
        )",
        logger);
  }, "ERROR: carrier_frequency oscillator ID 5 exceeds number of oscillators");
}

