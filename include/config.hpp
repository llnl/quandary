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
  // Logging
  MPILogger logger; ///< MPI-aware logger for output messages.

  // General options
  std::vector<size_t> nlevels; ///< Number of levels per subsystem
  std::vector<size_t> nessential; ///< Number of essential levels per subsystem (Default: same as nlevels)
  size_t ntime; ///< Number of time steps used for time-integration
  double dt; ///< Time step size (ns). Determines final time: T=ntime*dt
  std::vector<double> transfreq; ///< Fundamental transition frequencies for each oscillator (GHz)
  std::vector<double> selfkerr; ///< Self-kerr frequencies for each oscillator (GHz)
  std::vector<double> crosskerr; ///< Cross-kerr coupling frequencies for each oscillator coupling (GHz)
  std::vector<double> Jkl; ///< Dipole-dipole coupling frequencies for each oscillator coupling (GHz)
  std::vector<double> rotfreq; ///< Rotational wave approximation frequencies for each subsystem (GHz)
  DecoherenceType decoherence_type; ///< Switch between Schroedinger and Lindblad solver
  std::vector<double> decay_time; ///< Time of decay operation (T1) per oscillator (for Lindblad solver)
  std::vector<double> dephase_time; ///< Time of dephase operation (T2) per oscillator (for Lindblad solver)
  size_t n_initial_conditions; ///< Number of initial conditions
  InitialConditionSettings initial_condition; ///< Initial condition configuration
  std::optional<std::string> hamiltonian_file_Hsys; ///< File to read the system Hamiltonian from
  std::optional<std::string> hamiltonian_file_Hc; ///< File to read the control Hamiltonian from

  // Optimization options
  bool control_zero_boundary_condition; ///< Decide whether control pulses should start and end at zero
  std::vector<ControlParameterizationSettings> control_parameterizations; ///< Control parameterizations for each oscillator
  std::vector<ControlInitializationSettings> control_initializations; ///< Control initializations for each oscillator
  std::vector<double> control_amplitude_bounds; ///< Control amplitude bounds for each oscillator
  std::vector<std::vector<double>> carrier_frequencies; ///< Carrier frequencies for each oscillator
  OptimTargetSettings optim_target; ///< Grouped optimization target configuration
  ObjectiveType optim_objective; ///< Objective function measure
  std::vector<double> optim_weights; ///< Weights for summing up the objective function
  double optim_tol_grad_abs; ///< Absolute gradient tolerance
  double optim_tol_grad_rel; ///< Relative gradient tolerance
  double optim_tol_finalcost; ///< Final time cost tolerance
  double optim_tol_infidelity; ///< Infidelity tolerance
  size_t optim_maxiter; ///< Maximum iterations
  double optim_tikhonov_coeff; ///< Coefficient of Tikhonov regularization for the design variables
  bool optim_tikhonov_use_x0; ///< Switch to use Tikhonov regularization with ||x - x_0||^2 instead of ||x||^2
  double optim_penalty_leakage; ///< Leakage penalty coefficient
  double optim_penalty_weightedcost; ///< Weighted cost penalty coefficient
  double optim_penalty_weightedcost_width; ///< Width parameter for weighted cost penalty
  double optim_penalty_dpdm; ///< Second derivative penalty coefficient
  double optim_penalty_energy; ///< Energy penalty coefficient
  double optim_penalty_variation; ///< Amplitude variation penalty coefficient

  // Output and runtypes
  std::string output_directory; ///< Directory for output files
  std::vector<OutputType> output_observables; ///< Specify the desired observables.
  size_t output_timestep_stride; ///< Output frequency in the time domain: write output every <num> time-step
  size_t output_optimization_stride; ///< Frequency of writing output during optimization iterations
  RunType runtype; ///< Runtype options: simulation, gradient, or optimization
  bool usematfree; ///< Use matrix free solver, instead of sparse matrix implementation
  LinearSolverType linearsolver_type; ///< Solver type for solving the linear system at each time step
  size_t linearsolver_maxiter; ///< Set maximum number of iterations for the linear solver
  TimeStepperType timestepper_type; ///< The time-stepping algorithm
  int rand_seed; ///< Fixed seed for the random number generator for reproducibility

 public:
  Config(const MPILogger& logger, const toml::table& table);

  // TODO cfg: delete this when .cfg format is removed.
  Config(const MPILogger& logger, const ParsedConfigData& settings);

  ~Config() = default;

  static Config fromFile(const std::string& filename, const MPILogger& logger);
  static Config fromToml(const std::string& toml_filename, const MPILogger& logger);
  static Config fromTomlString(const std::string& toml_content, const MPILogger& logger);

  // TODO cfg: delete these when .cfg format is removed.
  static Config fromCfg(const std::string& cfg_filename, const MPILogger& logger);
  static Config fromCfgString(const std::string& cfg_content, const MPILogger& logger);

  void printConfig(std::stringstream& log) const;

  // getters
  const std::vector<size_t>& getNLevels() const { return nlevels; }
  size_t getNLevels(size_t i_osc) const { return nlevels[i_osc]; }
  size_t getNumOsc() const { return nlevels.size(); }
  const std::vector<size_t>& getNEssential() const { return nessential; }
  size_t getNEssential(size_t i_osc) const { return nessential[i_osc]; }
  size_t getNTime() const { return ntime; }
  double getDt() const { return dt; }
  double getTotalTime() const { return ntime * dt; }

  const std::vector<double>& getTransFreq() const { return transfreq; }
  const std::vector<double>& getSelfKerr() const { return selfkerr; }
  const std::vector<double>& getCrossKerr() const { return crosskerr; }
  const std::vector<double>& getJkl() const { return Jkl; }
  const std::vector<double>& getRotFreq() const { return rotfreq; }
  DecoherenceType getDecoherenceType() const { return decoherence_type; }
  const std::vector<double>& getDecayTime() const { return decay_time; }
  const std::vector<double>& getDephaseTime() const { return dephase_time; }
  size_t getNInitialConditions() const { return n_initial_conditions; }
  const InitialConditionSettings& getInitialCondition() const { return initial_condition; }
  const std::optional<std::string>& getHamiltonianFileHsys() const { return hamiltonian_file_Hsys; }
  const std::optional<std::string>& getHamiltonianFileHc() const { return hamiltonian_file_Hc; }

  const ControlParameterizationSettings& getControlParameterizations(size_t i_osc) const { return control_parameterizations[i_osc]; }
  bool getControlZeroBoundaryCondition() const { return control_zero_boundary_condition; }
  const ControlInitializationSettings& getControlInitializations(size_t i_osc) const {
    return control_initializations[i_osc];
  }
  double getControlAmplitudeBound(size_t i_osc) const { return control_amplitude_bounds[i_osc]; }
  const std::vector<double>& getCarrierFrequencies(size_t i_osc) const { return carrier_frequencies[i_osc]; }
  const OptimTargetSettings& getOptimTarget() const { return optim_target; }
  ObjectiveType getOptimObjective() const { return optim_objective; }
  const std::vector<double>& getOptimWeights() const { return optim_weights; }
  double getOptimTolGradAbs() const { return optim_tol_grad_abs; }
  double getOptimTolGradRel() const { return optim_tol_grad_rel; }
  double getOptimTolFinalCost() const { return optim_tol_finalcost; }
  double getOptimTolInfidelity() const { return optim_tol_infidelity; }
  size_t getOptimMaxiter() const { return optim_maxiter; }
  double getOptimTikhonovCoeff() const { return optim_tikhonov_coeff; }
  bool getOptimTikhonovUseX0() const { return optim_tikhonov_use_x0; }
  double getOptimPenaltyLeakage() const { return optim_penalty_leakage; }
  double getOptimPenaltyWeightedCost() const { return optim_penalty_weightedcost; }
  double getOptimPenaltyWeightedCostWidth() const { return optim_penalty_weightedcost_width; }
  double getOptimPenaltyDpdm() const { return optim_penalty_dpdm; }
  double getOptimPenaltyEnergy() const { return optim_penalty_energy; }
  double getOptimPenaltyVariation() const { return optim_penalty_variation; }

  const std::string& getOutputDirectory() const { return output_directory; }
  const std::vector<OutputType>& getOutputObservables() const { return output_observables; }
  size_t getOutputTimestepStride() const { return output_timestep_stride; }
  size_t getOutputOptimizationStride() const { return output_optimization_stride; }
  RunType getRuntype() const { return runtype; }
  bool getUseMatFree() const { return usematfree; }
  LinearSolverType getLinearSolverType() const { return linearsolver_type; }
  size_t getLinearSolverMaxiter() const { return linearsolver_maxiter; }
  TimeStepperType getTimestepperType() const { return timestepper_type; }
  int getRandSeed() const { return rand_seed; }

 private:
  void finalize();
  void validate() const;

  size_t computeNumInitialConditions(InitialConditionSettings init_cond_settings, std::vector<size_t> nlevels, std::vector<size_t> nessential, DecoherenceType decoherence_type) const;

  void setRandSeed(int rand_seed_);

  // Helper function to parse a single control parameterization table
  ControlParameterizationSettings parseControlParameterizationSpecs(const toml::table& param_table) const;

  // Helper function to parse a single control initialization table
  ControlInitializationSettings parseControlInitializationSpecs(const toml::table& table) const;

  /**
   * @brief Parses optimization target settings from TOML table
   *
   * @param table TOML table containing the configuration
   * @param num_osc Number of oscillators
   * @return Parsed optimization target settings
   */
  OptimTargetSettings parseOptimTarget(const toml::table& table, size_t num_osc) const;


  // TODO cfg: delete these when .cfg format is removed.
  template <typename T>
  std::vector<std::vector<T>> parseOscillatorSettingsCfg(const std::optional<std::map<int, std::vector<T>>>& indexed, size_t num_entries, const std::vector<T>& default_values = {}) const;

  std::vector<ControlParameterizationSettings> parseControlParameterizationsCfg(const std::optional<std::map<int, ControlParameterizationData>>& parameterizations_map) const;

  std::vector<ControlInitializationSettings> parseControlInitializationsCfg(const std::optional<std::map<int, ControlInitializationSettings>>& init_configs) const;
};
