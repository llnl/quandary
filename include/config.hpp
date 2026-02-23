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
 * - Setup (Wrapper = std::optional): All fields are optional (for input from TOML/Python)
 *
 * Example expansion for a field like `Wrapper<size_t> ntime`:
 *   - Setup:           std::optional<size_t> ntime;  (user may or may not set it)
 *   - ValidatedConfig: size_t ntime;                 (always has a validated value)
 *
 * Fields that are inherently optional (like hamiltonian files) don't use the Wrapper
 * and remain std::optional in both instantiations (optional even after validation).
 *
 * @tparam Wrapper Type wrapper template - either Identity or std::optional
 */
template <template <typename> class Wrapper>
struct ConfigFieldsT {
  // System parameters
  Wrapper<std::vector<size_t>> nlevels; ///< Number of levels per subsystem
  Wrapper<std::vector<size_t>> nessential; ///< Number of essential levels per subsystem (Default: same as nlevels)
  Wrapper<size_t> ntime; ///< Number of time steps used for time-integration
  Wrapper<double> dt; ///< Time step size (ns). Determines final time: T=ntime*dt
  Wrapper<std::vector<double>> transition_frequency; ///< Fundamental transition frequencies for each oscillator (GHz)
  Wrapper<std::vector<double>> selfkerr; ///< Self-kerr frequencies for each oscillator (GHz)
  Wrapper<std::vector<double>> crosskerr_coupling; ///< Cross-kerr coupling frequencies for each oscillator coupling (GHz)
  Wrapper<std::vector<double>> dipole_coupling; ///< Dipole-dipole coupling frequencies for each oscillator coupling (GHz)
  Wrapper<std::vector<double>> rotation_frequency; ///< Rotational wave approximation frequencies for each subsystem (GHz)
  Wrapper<DecoherenceType> decoherence_type; ///< Switch between Schroedinger and Lindblad solver
  Wrapper<std::vector<double>> decay_time; ///< Time of decay operation (T1) per oscillator (for Lindblad solver)
  Wrapper<std::vector<double>> dephase_time; ///< Time of dephase operation (T2) per oscillator (for Lindblad solver)
  Wrapper<InitialConditionSettings> initial_condition; ///< Initial condition configuration

  // Inherently optional - no Wrapper
  std::optional<std::string> hamiltonian_file_Hsys; ///< File to read the system Hamiltonian from
  std::optional<std::string> hamiltonian_file_Hc; ///< File to read the control Hamiltonian from

  // Control parameters
  Wrapper<bool> control_zero_boundary_condition; ///< Decide whether control pulses should start and end at zero
  Wrapper<std::vector<ControlParameterizationSettings>> control_parameterizations; ///< Control parameterizations for each oscillator
  Wrapper<std::vector<ControlInitializationSettings>> control_initializations; ///< Control initializations for each oscillator
  Wrapper<std::vector<double>> control_amplitude_bounds; ///< Control amplitude bounds for each oscillator
  Wrapper<std::vector<std::vector<double>>> carrier_frequencies; ///< Carrier frequencies for each oscillator

  // Optimization parameters
  Wrapper<OptimTargetSettings> optim_target; ///< Grouped optimization target configuration
  Wrapper<ObjectiveType> optim_objective; ///< Objective function measure
  Wrapper<std::vector<double>> optim_weights; ///< Weights for summing up the objective function
  Wrapper<double> optim_tol_grad_abs; ///< Absolute gradient tolerance
  Wrapper<double> optim_tol_grad_rel; ///< Relative gradient tolerance
  Wrapper<double> optim_tol_final_cost; ///< Final time cost tolerance
  Wrapper<double> optim_tol_infidelity; ///< Infidelity tolerance
  Wrapper<size_t> optim_maxiter; ///< Maximum iterations
  Wrapper<double> optim_tikhonov_coeff; ///< Coefficient of Tikhonov regularization for the design variables
  Wrapper<bool> optim_tikhonov_use_x0; ///< Switch to use Tikhonov regularization with ||x - x_0||^2 instead of ||x||^2
  Wrapper<double> optim_penalty_leakage; ///< Leakage penalty coefficient
  Wrapper<double> optim_penalty_weightedcost; ///< Weighted cost penalty coefficient
  Wrapper<double> optim_penalty_weightedcost_width; ///< Width parameter for weighted cost penalty
  Wrapper<double> optim_penalty_dpdm; ///< Second derivative penalty coefficient
  Wrapper<double> optim_penalty_energy; ///< Energy penalty coefficient
  Wrapper<double> optim_penalty_variation; ///< Amplitude variation penalty coefficient

  // Output parameters
  Wrapper<std::string> output_directory; ///< Directory for output files
  Wrapper<std::vector<OutputType>> output_observables; ///< Specify the desired observables.
  Wrapper<size_t> output_timestep_stride; ///< Output frequency in the time domain: write output every N time-steps
  Wrapper<size_t> output_optimization_stride; ///< Frequency of writing output during optimization iterations

  // Solver parameters
  Wrapper<RunType> runtype; ///< Runtype options: simulation, gradient, or optimization
  Wrapper<bool> usematfree; ///< Use matrix free solver, instead of sparse matrix implementation
  Wrapper<LinearSolverType> linearsolver_type; ///< Solver type for solving the linear system at each time step
  Wrapper<size_t> linearsolver_maxiter; ///< Set maximum number of iterations for the linear solver
  Wrapper<TimeStepperType> timestepper_type; ///< The time-stepping algorithm
  Wrapper<int> rand_seed; ///< Fixed seed for the random number generator for reproducibility
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
using Setup = ConfigFieldsT<std::optional>;

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
   * @brief Constructs a Config from a Setup struct
   *
   * This constructor performs the main validation and default value application for all
   * configuration parameters. It uses the validators framework to check constraints
   * and apply defaults in a consistent way.
   *
   * @param input The pre-parsed configuration input data
   * @param quiet_mode Whether to suppress logging output
   */
  Config(const Setup& input, bool quiet_mode = false);
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

  // toString functions for settings structs (used for printing and Python __repr__)
  static std::string toString(const InitialConditionSettings& initial_condition);
  static std::string toString(const OptimTargetSettings& optim_target);
  static std::string toString(const ControlParameterizationSettings& control_param);
  static std::string toString(const ControlInitializationSettings& control_init);

  // getters
  const std::vector<size_t>& getNLevels() const { return data.nlevels; }
  size_t getNLevels(size_t i_osc) const { return data.nlevels[i_osc]; }
  size_t getNumOsc() const { return data.nlevels.size(); }
  const std::vector<size_t>& getNEssential() const { return data.nessential; }
  size_t getNEssential(size_t i_osc) const { return data.nessential[i_osc]; }
  size_t getNTime() const { return data.ntime; }
  double getDt() const { return data.dt; }
  double getTotalTime() const { return data.ntime * data.dt; }

  const std::vector<double>& getTransitionFrequency() const { return data.transition_frequency; }
  const std::vector<double>& getSelfKerr() const { return data.selfkerr; }
  const std::vector<double>& getCrossKerrCoupling() const { return data.crosskerr_coupling; }
  const std::vector<double>& getDipoleCoupling() const { return data.dipole_coupling; }
  const std::vector<double>& getRotationFrequency() const { return data.rotation_frequency; }
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
  double getOptimTolFinalCost() const { return data.optim_tol_final_cost; }
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
