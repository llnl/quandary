#pragma once

#include <algorithm>
#include <petsc.h>
#include <cstddef>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <toml++/toml.hpp>
#include <vector>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <limits>
#include <random>
#include <type_traits>
#include "cfgparser.hpp"
#include "defs.hpp"
#include "mpi_logger.hpp"
#include "config_defaults.hpp"
#include "config_validators.hpp"

/**
 * @brief Identity type wrapper - passes through the type unchanged.
 *
 * Used with ConfigFieldsT to create ValidatedConfig where all fields are concrete types.
 */
template <typename T>
using Identity = T;

/**
 * @brief Templated configuration data struct for field definitions.
 *
 * This template allows the same field definitions to be instantiated as:
 * - ValidatedConfig (Wrapper = Identity): All fields are concrete types (validated config)
 * - RawConfig (Wrapper = std::optional): All fields are optional (for input from TOML/Python)
 *
 * Fields that are inherently optional (like hamiltonian files) don't use the Wrapper
 * and remain std::optional in both instantiations.
 *
 * @tparam Wrapper Type wrapper template - either Identity or std::optional
 */
template <template <typename> class Wrapper>
struct ConfigFieldsT {
  // System parameters
  Wrapper<std::vector<size_t>> nlevels;
  Wrapper<std::vector<size_t>> nessential;
  Wrapper<size_t> ntime;
  Wrapper<double> dt;
  Wrapper<std::vector<double>> transfreq;
  Wrapper<std::vector<double>> selfkerr;
  Wrapper<std::vector<double>> crosskerr;
  Wrapper<std::vector<double>> Jkl;
  Wrapper<std::vector<double>> rotfreq;
  Wrapper<DecoherenceType> decoherence_type;
  Wrapper<std::vector<double>> decay_time;
  Wrapper<std::vector<double>> dephase_time;
  Wrapper<InitialConditionSettings> initial_condition;

  // Inherently optional - no Wrapper
  std::optional<std::string> hamiltonian_file_Hsys;
  std::optional<std::string> hamiltonian_file_Hc;

  // Control parameters
  Wrapper<bool> control_zero_boundary_condition;
  Wrapper<std::vector<ControlParameterizationSettings>> control_parameterizations;
  Wrapper<std::vector<ControlInitializationSettings>> control_initializations;
  Wrapper<std::vector<double>> control_amplitude_bounds;
  Wrapper<std::vector<std::vector<double>>> carrier_frequencies;

  // Optimization parameters
  Wrapper<OptimTargetSettings> optim_target;
  Wrapper<ObjectiveType> optim_objective;
  Wrapper<std::vector<double>> optim_weights;
  Wrapper<double> optim_tol_grad_abs;
  Wrapper<double> optim_tol_grad_rel;
  Wrapper<double> optim_tol_finalcost;
  Wrapper<double> optim_tol_infidelity;
  Wrapper<size_t> optim_maxiter;
  Wrapper<double> optim_tikhonov_coeff;
  Wrapper<bool> optim_tikhonov_use_x0;
  Wrapper<double> optim_penalty_leakage;
  Wrapper<double> optim_penalty_weightedcost;
  Wrapper<double> optim_penalty_weightedcost_width;
  Wrapper<double> optim_penalty_dpdm;
  Wrapper<double> optim_penalty_energy;
  Wrapper<double> optim_penalty_variation;

  // Output parameters
  Wrapper<std::string> output_directory;
  Wrapper<std::vector<OutputType>> output_observables;
  Wrapper<size_t> output_timestep_stride;
  Wrapper<size_t> output_optimization_stride;

  // Solver parameters
  Wrapper<RunType> runtype;
  Wrapper<bool> usematfree;
  Wrapper<LinearSolverType> linearsolver_type;
  Wrapper<size_t> linearsolver_maxiter;
  Wrapper<TimeStepperType> timestepper_type;
  Wrapper<int> rand_seed;
};

/**
 * @brief Validated configuration data with concrete types.
 *
 * All fields have been validated and contain final values.
 */
using ValidatedConfig = ConfigFieldsT<Identity>;

/**
 * @brief Raw input configuration data with optional fields.
 *
 * Used for receiving configuration from TOML parsing or Python bindings.
 * All fields are optional to allow partial specification with defaults.
 */
using RawConfig = ConfigFieldsT<std::optional>;

/**
 * @brief Configuration class containing all validated settings.
 *
 * Contains validated, typed configuration parameters. All fields have been
 * validated with defaults set. Handles parsing from TOML configuration files
 * and deprecated CFG format, as well as printing log of used configuration.
 * This class is immutable after construction.
 *
 * @note Adding a new toml configuration option:
 *
 * 1) Add new member variable to the Config class below
 *
 * 2) Add parsing logic in constructor Config::Config(toml::table) in src/config.cpp. Either access directly from toml table, or use chains of validators for type-safe extraction and validation. For example:
 *
 *    Parsing simple scalar fields (T foo):
 *      - `foo = toml["foo"].value()`: required scalar field, no default, not validated)
 *      - `foo = toml["foo"].value_or(default_value)`: optional scalar with default, not validated
 *      - `foo = validators::field<T>(toml, "foo").<validation_chain>.value()`: validated required scalar
 *      - `foo = validators::field<T>(toml, "foo").<validation_chain>.valueOr(default_value)`: validated optional scalar with default
 *
 *    Parsing vectors (std::vector<T> foo):
 *      - `foo = toml["foo"].as_array()->as<T>().value_or(default_vector)`: vector with default, not validated
 *      - `foo = validators::vectorField<T>(toml, "foo").<validation_chain>.value()`: validated required vector
 *      - `foo = validators::vectorField<T>(toml, "foo").<validation_chain>.value(default_vector)`: validated optional vector with default
 *
 *    Available validation chain methods can found in include/config_validators.hpp, such as:
 *      `.positive()`, `.greaterThan(val)`, `.lessThan(val)`, `.minLength(len)`, `.hasLength(len)`, etc.
 *
 * 3) Add a getter method to the Config class below (`T getFoo() const {return foo};`)
 *
 * 4) Add printing logic in Config::printConfig() in src/config.cpp
 */
class Config {
 private:
  MPILogger logger; ///< MPI-aware logger for output messages.
  ValidatedConfig data; ///< All validated configuration fields.
  size_t n_initial_conditions; ///< Number of initial conditions (computed, not in ValidatedConfig)

 public:
  /**
   * @brief Constructs a Config from a RawConfig struct
   *
   * This constructor performs the main validation and default value application for all
   * configuration parameters. It uses the validators framework to check constraints
   * and apply defaults in a consistent way.
   *
   * @param input The pre-parsed configuration input data
   * @param quiet_mode Whether to suppress logging output
   */
  Config(const RawConfig& input, bool quiet_mode = false);
  Config(const toml::table& table, bool quiet_mode = false);

  // TODO cfg: delete this when .cfg format is removed.
  Config(const ParsedConfigData& settings, bool quiet_mode = false);

  ~Config() = default;

  static Config fromFile(const std::string& filename, bool quiet_mode = false);
  static Config fromToml(const std::string& toml_filename, bool quiet_mode = false);
  static Config fromTomlString(const std::string& toml_content, bool quiet_mode = false);

  // TODO cfg: delete these when .cfg format is removed.
  static Config fromCfg(const std::string& cfg_filename, bool quiet_mode = false);
  static Config fromCfgString(const std::string& cfg_content, bool quiet_mode = false);

  void printConfig(std::stringstream& log) const;

  /**
   * @brief Configure this Config for control evaluation at a specific sample rate.
   *
   * Modifies time grid based on points_per_ns while keeping total time constant,
   * sets run type to EVALCONTROLS, and configures control initialization to read
   * from the specified file. Caution: this modifies the Config in-place.
   *
   * @param points_per_ns Sample rate [points per nanosecond]
   * @param pcof_file Path to file containing control parameters
   * @param output_dir Output directory for results
   */
  void setupForEvalControls(double points_per_ns, const std::string& pcof_file, const std::string& output_dir);

  // getters
  const std::vector<size_t>& getNLevels() const { return data.nlevels; }
  size_t getNLevels(size_t i_osc) const { return data.nlevels[i_osc]; }
  size_t getNumOsc() const { return data.nlevels.size(); }
  const std::vector<size_t>& getNEssential() const { return data.nessential; }
  size_t getNEssential(size_t i_osc) const { return data.nessential[i_osc]; }
  size_t getNTime() const { return data.ntime; }
  double getDt() const { return data.dt; }
  double getTotalTime() const { return data.ntime * data.dt; }

  const std::vector<double>& getTransFreq() const { return data.transfreq; }
  const std::vector<double>& getSelfKerr() const { return data.selfkerr; }
  const std::vector<double>& getCrossKerr() const { return data.crosskerr; }
  const std::vector<double>& getJkl() const { return data.Jkl; }
  const std::vector<double>& getRotFreq() const { return data.rotfreq; }
  DecoherenceType getDecoherenceType() const { return data.decoherence_type; }
  const std::vector<double>& getDecayTime() const { return data.decay_time; }
  const std::vector<double>& getDephaseTime() const { return data.dephase_time; }
  size_t getNInitialConditions() const { return n_initial_conditions; }
  const InitialConditionSettings& getInitialCondition() const { return data.initial_condition; }
  const std::optional<std::string>& getHamiltonianFileHsys() const { return data.hamiltonian_file_Hsys; }
  const std::optional<std::string>& getHamiltonianFileHc() const { return data.hamiltonian_file_Hc; }

  const ControlParameterizationSettings& getControlParameterizations(size_t i_osc) const { return data.control_parameterizations[i_osc]; }
  bool getControlZeroBoundaryCondition() const { return data.control_zero_boundary_condition; }
  const ControlInitializationSettings& getControlInitializations(size_t i_osc) const {
    return data.control_initializations[i_osc];
  }
  double getControlAmplitudeBound(size_t i_osc) const { return data.control_amplitude_bounds[i_osc]; }
  const std::vector<double>& getCarrierFrequencies(size_t i_osc) const { return data.carrier_frequencies[i_osc]; }
  const OptimTargetSettings& getOptimTarget() const { return data.optim_target; }
  ObjectiveType getOptimObjective() const { return data.optim_objective; }
  const std::vector<double>& getOptimWeights() const { return data.optim_weights; }
  double getOptimTolGradAbs() const { return data.optim_tol_grad_abs; }
  double getOptimTolGradRel() const { return data.optim_tol_grad_rel; }
  double getOptimTolFinalCost() const { return data.optim_tol_finalcost; }
  double getOptimTolInfidelity() const { return data.optim_tol_infidelity; }
  size_t getOptimMaxiter() const { return data.optim_maxiter; }
  double getOptimTikhonovCoeff() const { return data.optim_tikhonov_coeff; }
  bool getOptimTikhonovUseX0() const { return data.optim_tikhonov_use_x0; }
  double getOptimPenaltyLeakage() const { return data.optim_penalty_leakage; }
  double getOptimPenaltyWeightedCost() const { return data.optim_penalty_weightedcost; }
  double getOptimPenaltyWeightedCostWidth() const { return data.optim_penalty_weightedcost_width; }
  double getOptimPenaltyDpdm() const { return data.optim_penalty_dpdm; }
  double getOptimPenaltyEnergy() const { return data.optim_penalty_energy; }
  double getOptimPenaltyVariation() const { return data.optim_penalty_variation; }

  const std::string& getOutputDirectory() const { return data.output_directory; }
  const std::vector<OutputType>& getOutputObservables() const { return data.output_observables; }
  size_t getOutputTimestepStride() const { return data.output_timestep_stride; }
  size_t getOutputOptimizationStride() const { return data.output_optimization_stride; }
  RunType getRuntype() const { return data.runtype; }
  bool getUseMatFree() const { return data.usematfree; }
  LinearSolverType getLinearSolverType() const { return data.linearsolver_type; }
  size_t getLinearSolverMaxiter() const { return data.linearsolver_maxiter; }
  TimeStepperType getTimestepperType() const { return data.timestepper_type; }
  int getRandSeed() const { return data.rand_seed; }

 private:
  /**
   * @brief Finalizes configuration after initial validation
   *
   * This method is called after the constructor finishes basic validation
   * and applies secondary transformations, dependencies between parameters,
   * and final default value application that depends on other parameters.
   */
  void finalize();

  /**
   * @brief Validates final configuration for consistency
   *
   * This method performs additional validation on the finalized configuration
   * to ensure all parameters are consistent with each other and satisfy
   * requirements for simulation. Called after finalize().
   */
  void validate() const;

  size_t computeNumInitialConditions(InitialConditionSettings init_cond_settings, std::vector<size_t> nlevels, std::vector<size_t> nessential, DecoherenceType decoherence_type) const;

  void setRandSeed(int rand_seed_);

  // TODO cfg: delete these when .cfg format is removed.
  template <typename T>
  std::vector<std::vector<T>> parseOscillatorSettingsCfg(const std::optional<std::map<int, std::vector<T>>>& indexed, size_t num_entries, const std::vector<T>& default_values = {}) const;

  std::vector<ControlParameterizationSettings> parseControlParameterizationsCfg(const std::optional<std::map<int, ControlParameterizationData>>& parameterizations_map) const;

  std::vector<ControlInitializationSettings> parseControlInitializationsCfg(const std::optional<std::map<int, ControlInitializationSettings>>& init_configs) const;
};
