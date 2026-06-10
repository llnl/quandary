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
#include "defs.hpp"
#include "mpi_logger.hpp"
#include "config_defaults.hpp"
#include "config_validators.hpp"

/**
 * @brief Templated struct containing all configuration data definitions.
 *
 * Each data field is wrapped in a template Wrapper, which is either 
 * - Wrapper = Identity: for ValidatedConfig, where all fields have concrete types 
 *   after validating the user input and applying defaults.
 * - Wrapper = std::optional: for ConfigInput, where all fields are std::optional to allow 
 *   flexible user input from TOML or Python. Validation and defaulting happens in the
 *   Config constructor.
 * 
 * Fields that are inherently optional (like hamiltonian files) don't use the Wrapper
 * and remain std::optional in both instantiations (optional even after validation).
 *
 * @tparam Wrapper Type wrapper template - either Identity or std::optional
 */
template <template <typename> class Wrapper>
struct ConfigDataT {
  // System parameters
  Wrapper<std::vector<size_t>> nlevels; ///< Number of levels per subsystem
  Wrapper<std::vector<size_t>> nessential; ///< Number of essential levels per subsystem (Default: same as nlevels)
  Wrapper<size_t> ntime; ///< Number of time steps used for time-integration
  Wrapper<double> dt; ///< Time step size (ns). Determines final time: T=ntime*dt
  Wrapper<double> total_time; ///< Total evolution time (ns). Alternative to specifying ntime and dt.
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
  Wrapper<std::vector<double>> control_amplitude_bound; ///< Control amplitude bounds for each oscillator
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
template <typename T>
using Identity = T;
using ValidatedConfig = ConfigDataT<Identity>;

/**
 * @brief Raw input configuration data with optional fields.
 *
 * Used for receiving configuration from TOML parsing or Python bindings.
 * All fields are optional to allow partial specification with defaults.
 */
using ConfigInput = ConfigDataT<std::optional>;

/**
 * @brief Configuration class containing all validated settings.
 *
 * Contains the validated, typed configuration. All data fields have been
 * validated with defaults set. Handles parsing from TOML configuration files
 * and printing log of used configuration.
 * This class is immutable after construction.
 *
 * @note Error handling: Use exceptions (std::runtime_error, std::invalid_argument)
 * instead of logger.exitWithError(). The Config class is called from Python via
 * nanobind, and exitWithError terminates the process — crashing Jupyter notebooks
 * and other interactive sessions. Exceptions are caught by nanobind and converted
 * to Python exceptions, allowing graceful error handling.
 *
 * @note Adding a new toml configuration option:
 *
 * 1) Add new member variable to the ConfigDataT struct
 * 2) Add TOML extraction logic in `extractConfigInput()` (src/config.cpp),
 * 3) Add validation/defaulting logic in `Config::Config(const ConfigInput&)`.
 *    Extraction and validation are intentionally split. 
 * 4) Add a getter method to the Config class below 
 * 5) Add printing logic in Config::printConfig() in src/config.cpp
 * 6) Add Python bindings in python/bindings.cpp if the new config option should be accessible from Python. For example, add `.def_prop_ro("foo", &Config::getFoo, "Description of foo")` to the nb::class_<Config> definition.
 */
class Config {
 private:
  MPILogger logger; ///< MPI-aware logger for output messages.
  ValidatedConfig validated_config; ///< All validated configuration fields.
  size_t n_initial_conditions; ///< Number of initial conditions (computed, not in ValidatedConfig)

 public:
  /**
  * @brief Constructs a Config from a ConfigInput struct (primary constructor).
   *
   * This is the main constructor that all other constructors delegate to.
   * It takes a pre-parsed ConfigInput struct where all fields are std::optional, 
   * validates all fields using the validator chains in config_validators.hpp, applies 
   * defaults for missing optional fields, and calls finalize()/validate() for cross-field 
   * consistency checks.
   *
   * Used directly from Python (nanobind) and programmatic C++ use, where the
   * caller populates a ConfigInput struct with std::optional fields.
   *
   * @param input The pre-parsed configuration input data
   * @param quiet_mode Whether to suppress logging output
   */
  Config(const ConfigInput& input, bool quiet_mode = false);

  /**
   * @brief Constructs a Config from a TOML table.
   *
   * Extracts TOML values into a ConfigInput struct (via internal `extractConfigInput()`
   * in src/config.cpp), where all config data fields are std::optional, then delegates 
   * to `Config(const ConfigInput&)` to set up and store the ValidatedConfig instance where 
   * all fields have concrete types after validation and defaulting.
   *
   * @param table Parsed TOML table
   * @param quiet_mode Whether to suppress logging output
   */
  Config(const toml::table& table, bool quiet_mode = false);

  ~Config() = default;

  /**
   * @brief Parses a TOML file and constructs a Config object from the TOML table.
   * @param filename Path to the configuration file
   * @param quiet_mode Whether to suppress logging output
   */
  static Config fromFile(const std::string& filename, bool quiet_mode = false);

  /**
   * @brief Parses a TOML-formatted string and constructs a Config object from the TOML table.
   * @param toml_content String containing TOML configuration
   * @param quiet_mode Whether to suppress logging output
   */
  static Config fromString(const std::string& toml_content, bool quiet_mode = false);

  /**
   * @brief Prints the validated configuration to TOML format.
   * @param log Output stream to write the TOML representation to
   */
  void printConfig(std::stringstream& log) const;

  /**
   * @brief Convert settings structs to human-readable strings.
   *
   * Used for TOML config logging and Python __repr__ methods.
   */
  static std::string toString(const InitialConditionSettings& initial_condition);
  static std::string toString(const OptimTargetSettings& optim_target);
  static std::string toString(const ControlParameterizationSettings& control_param);
  static std::string toString(const ControlInitializationSettings& control_init);

  // getters
  const std::vector<size_t>& getNLevels() const { return validated_config.nlevels; }
  size_t getNLevels(size_t i_osc) const { return validated_config.nlevels[i_osc]; }
  size_t getNumOsc() const { return validated_config.nlevels.size(); }
  const std::vector<size_t>& getNEssential() const { return validated_config.nessential; }
  size_t getNEssential(size_t i_osc) const { return validated_config.nessential[i_osc]; }
  size_t getNTime() const { return validated_config.ntime; }
  double getDt() const { return validated_config.dt; }
  double getTotalTime() const { return validated_config.total_time; }

  const std::vector<double>& getTransitionFrequency() const { return validated_config.transition_frequency; }
  const std::vector<double>& getSelfKerr() const { return validated_config.selfkerr; }
  const std::vector<double>& getCrossKerrCoupling() const { return validated_config.crosskerr_coupling; }
  const std::vector<double>& getDipoleCoupling() const { return validated_config.dipole_coupling; }
  const std::vector<double>& getRotationFrequency() const { return validated_config.rotation_frequency; }
  DecoherenceType getDecoherenceType() const { return validated_config.decoherence_type; }
  const std::vector<double>& getDecayTime() const { return validated_config.decay_time; }
  const std::vector<double>& getDephaseTime() const { return validated_config.dephase_time; }
  size_t getNInitialConditions() const { return n_initial_conditions; }
  const InitialConditionSettings& getInitialCondition() const { return validated_config.initial_condition; }
  const std::optional<std::string>& getHamiltonianFileHsys() const { return validated_config.hamiltonian_file_Hsys; }
  const std::optional<std::string>& getHamiltonianFileHc() const { return validated_config.hamiltonian_file_Hc; }

  const ControlParameterizationSettings& getControlParameterizations(size_t i_osc) const { return validated_config.control_parameterizations[i_osc]; }
  bool getControlZeroBoundaryCondition() const { return validated_config.control_zero_boundary_condition; }
  const ControlInitializationSettings& getControlInitializations(size_t i_osc) const {
    return validated_config.control_initializations[i_osc];
  }
  double getControlAmplitudeBound(size_t i_osc) const { return validated_config.control_amplitude_bound[i_osc]; }
  const std::vector<double>& getCarrierFrequencies(size_t i_osc) const { return validated_config.carrier_frequencies[i_osc]; }
  const OptimTargetSettings& getOptimTarget() const { return validated_config.optim_target; }
  ObjectiveType getOptimObjective() const { return validated_config.optim_objective; }
  const std::vector<double>& getOptimWeights() const { return validated_config.optim_weights; }
  double getOptimTolGradAbs() const { return validated_config.optim_tol_grad_abs; }
  double getOptimTolGradRel() const { return validated_config.optim_tol_grad_rel; }
  double getOptimTolFinalCost() const { return validated_config.optim_tol_final_cost; }
  double getOptimTolInfidelity() const { return validated_config.optim_tol_infidelity; }
  size_t getOptimMaxiter() const { return validated_config.optim_maxiter; }
  double getOptimTikhonovCoeff() const { return validated_config.optim_tikhonov_coeff; }
  bool getOptimTikhonovUseX0() const { return validated_config.optim_tikhonov_use_x0; }
  double getOptimPenaltyLeakage() const { return validated_config.optim_penalty_leakage; }
  double getOptimPenaltyWeightedCost() const { return validated_config.optim_penalty_weightedcost; }
  double getOptimPenaltyWeightedCostWidth() const { return validated_config.optim_penalty_weightedcost_width; }
  double getOptimPenaltyDpdm() const { return validated_config.optim_penalty_dpdm; }
  double getOptimPenaltyEnergy() const { return validated_config.optim_penalty_energy; }
  double getOptimPenaltyVariation() const { return validated_config.optim_penalty_variation; }

  const std::string& getOutputDirectory() const { return validated_config.output_directory; }
  const std::vector<OutputType>& getOutputObservables() const { return validated_config.output_observables; }
  size_t getOutputTimestepStride() const { return validated_config.output_timestep_stride; }
  size_t getOutputOptimizationStride() const { return validated_config.output_optimization_stride; }
  RunType getRuntype() const { return validated_config.runtype; }
  bool getUseMatFree() const { return validated_config.usematfree; }
  LinearSolverType getLinearSolverType() const { return validated_config.linearsolver_type; }
  size_t getLinearSolverMaxiter() const { return validated_config.linearsolver_maxiter; }
  TimeStepperType getTimestepperType() const { return validated_config.timestepper_type; }
  int getRandSeed() const { return validated_config.rand_seed; }

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
};
