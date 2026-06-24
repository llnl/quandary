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

  // General options
  std::optional<std::vector<size_t>> nlevels; ///< Number of levels per subsystem
  std::optional<std::vector<size_t>> nessential; ///< Number of essential levels per subsystem (Default: same as nlevels)
  std::optional<size_t> ntime; ///< Number of time steps used for time-integration
  std::optional<double> dt; ///< Time step size (ns). Determines final time: T=ntime*dt
  std::optional<double> total_time; ///< Total evolution time (ns). Alternative to specifying ntime and dt.
  std::optional<std::vector<double>> transition_frequency; ///< Fundamental transition frequencies for each oscillator (GHz)
  std::optional<std::vector<double>> selfkerr; ///< Self-kerr frequencies for each oscillator (GHz)
  std::optional<std::vector<double>> crosskerr_coupling; ///< Cross-kerr coupling frequencies for each oscillator coupling (GHz)
  std::optional<std::vector<double>> dipole_coupling; ///< Dipole-dipole coupling frequencies for each oscillator coupling (GHz)
  std::optional<std::vector<double>> rotation_frequency; ///< Rotational wave approximation frequencies for each subsystem (GHz)
  std::optional<DecoherenceType> decoherence_type; ///< Switch between Schroedinger and Lindblad solver
  std::optional<std::vector<double>> decay_time; ///< Time of decay operation (T1) per oscillator (for Lindblad solver)
  std::optional<std::vector<double>> dephase_time; ///< Time of dephase operation (T2) per oscillator (for Lindblad solver)
  std::optional<InitialConditionSettings> initial_condition; ///< Initial condition configuration
  std::optional<std::string> hamiltonian_file_Hsys; ///< File to read the system Hamiltonian from
  std::optional<std::string> hamiltonian_file_Hc; ///< File to read the control Hamiltonian from

  // Optimization options
  std::optional<bool> control_zero_boundary_condition; ///< Decide whether control pulses should start and end at zero
  std::optional<std::vector<ControlParameterizationSettings>> control_parameterizations; ///< Control parameterizations for each oscillator
  std::optional<std::vector<ControlInitializationSettings>> control_initializations; ///< Control initializations for each oscillator
  std::optional<std::vector<double>> control_amplitude_bound; ///< Control amplitude bounds for each oscillator
  std::optional<std::vector<std::vector<double>>> carrier_frequencies; ///< Carrier frequencies for each oscillator
  std::optional<OptimTargetSettings> optim_target; ///< Grouped optimization target configuration
  std::optional<ObjectiveType> optim_objective; ///< Objective function measure
  std::optional<std::vector<double>> optim_weights; ///< Weights for summing up the objective function
  std::optional<double> optim_tol_grad_abs; ///< Absolute gradient tolerance
  std::optional<double> optim_tol_grad_rel; ///< Relative gradient tolerance
  std::optional<double> optim_tol_final_cost; ///< Final time cost tolerance
  std::optional<double> optim_tol_infidelity; ///< Infidelity tolerance
  std::optional<size_t> optim_maxiter; ///< Maximum iterations
  std::optional<double> optim_tikhonov_coeff; ///< Coefficient of Tikhonov regularization for the design variables
  std::optional<bool> optim_tikhonov_use_x0; ///< Switch to use Tikhonov regularization with ||x - x_0||^2 instead of ||x||^2
  std::optional<double> optim_penalty_leakage; ///< Leakage penalty coefficient
  std::optional<double> optim_penalty_weightedcost; ///< Weighted cost penalty coefficient
  std::optional<double> optim_penalty_weightedcost_width; ///< Width parameter for weighted cost penalty
  std::optional<double> optim_penalty_dpdm; ///< Second derivative penalty coefficient
  std::optional<double> optim_penalty_energy; ///< Energy penalty coefficient
  std::optional<double> optim_penalty_variation; ///< Amplitude variation penalty coefficient

  // Output and runtypes
  std::optional<std::string> output_directory; ///< Directory for output files
  std::optional<std::vector<OutputType>> output_observables; ///< Specify the desired observables.
  std::optional<size_t> output_timestep_stride; ///< Output frequency in the time domain: write output every <num> time-step
  std::optional<size_t> output_optimization_stride; ///< Frequency of writing output during optimization iterations
  std::optional<RunType> runtype; ///< Runtype options: simulation, gradient, or optimization
  std::optional<bool> usematfree; ///< Use matrix free solver, instead of sparse matrix implementation
  std::optional<LinearSolverType> linearsolver_type; ///< Solver type for solving the linear system at each time step
  std::optional<size_t> linearsolver_maxiter; ///< Set maximum number of iterations for the linear solver
  std::optional<TimeStepperType> timestepper_type; ///< The time-stepping algorithm
  std::optional<int> rand_seed; ///< Fixed seed for the random number generator for reproducibility

  size_t n_initial_conditions; ///< Number of initial conditions (computed, nor parsed)

 public:
  Config(bool quiet_mode); // for python
  Config(const toml::table& toml, bool quiet_mode = false); // Parses the content, sets defaults, and validates the configuration
  static Config fromFile(const std::string& filename, bool quiet_mode = false);
  static Config fromString(const std::string& toml_content, bool quiet_mode = false);

  ~Config() = default;

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
  const std::vector<size_t>& getNLevels() const { return nlevels.value(); }
  size_t getNLevels(size_t i_osc) const { return nlevels.value()[i_osc]; }
  size_t getNumOsc() const { return nlevels.value().size(); }
  const std::vector<size_t>& getNEssential() const { return nessential.value(); }
  size_t getNEssential(size_t i_osc) const { return nessential.value()[i_osc]; }
  size_t getNTime() const { return ntime.value(); }
  double getDt() const { return dt.value(); }
  double getTotalTime() const { return total_time.value(); }

  const std::vector<double>& getTransitionFrequency() const { return transition_frequency.value(); }
  const std::vector<double>& getSelfKerr() const { return selfkerr.value(); }
  const std::vector<double>& getCrossKerrCoupling() const { return crosskerr_coupling.value(); }
  const std::vector<double>& getDipoleCoupling() const { return dipole_coupling.value(); }
  const std::vector<double>& getRotationFrequency() const { return rotation_frequency.value(); }
  DecoherenceType getDecoherenceType() const { return decoherence_type.value(); }
  const std::vector<double>& getDecayTime() const { return decay_time.value(); }
  const std::vector<double>& getDephaseTime() const { return dephase_time.value(); }
  size_t getNInitialConditions() const { return n_initial_conditions; }
  const InitialConditionSettings& getInitialCondition() const { return initial_condition.value(); }
  const std::optional<std::string>& getHamiltonianFileHsys() const { return hamiltonian_file_Hsys; }
  const std::optional<std::string>& getHamiltonianFileHc() const { return hamiltonian_file_Hc; }

  const ControlParameterizationSettings& getControlParameterizations(size_t i_osc) const { return control_parameterizations.value()[i_osc]; }
  bool getControlZeroBoundaryCondition() const { return control_zero_boundary_condition.value(); }
  const ControlInitializationSettings& getControlInitializations(size_t i_osc) const {
    return control_initializations.value()[i_osc];
  }
  double getControlAmplitudeBound(size_t i_osc) const { return control_amplitude_bound.value()[i_osc]; }
  const std::vector<double>& getCarrierFrequencies(size_t i_osc) const { return carrier_frequencies.value()[i_osc]; }
  const OptimTargetSettings& getOptimTarget() const { return optim_target.value(); }
  ObjectiveType getOptimObjective() const { return optim_objective.value(); }
  const std::vector<double>& getOptimWeights() const { return optim_weights.value(); }
  double getOptimTolGradAbs() const { return optim_tol_grad_abs.value(); }
  double getOptimTolGradRel() const { return optim_tol_grad_rel.value(); }
  double getOptimTolFinalCost() const { return optim_tol_final_cost.value(); }
  double getOptimTolInfidelity() const { return optim_tol_infidelity.value(); }
  size_t getOptimMaxiter() const { return optim_maxiter.value(); }
  double getOptimTikhonovCoeff() const { return optim_tikhonov_coeff.value(); }
  bool getOptimTikhonovUseX0() const { return optim_tikhonov_use_x0.value(); }
  double getOptimPenaltyLeakage() const { return optim_penalty_leakage.value(); }
  double getOptimPenaltyWeightedCost() const { return optim_penalty_weightedcost.value(); }
  double getOptimPenaltyWeightedCostWidth() const { return optim_penalty_weightedcost_width.value(); }
  double getOptimPenaltyDpdm() const { return optim_penalty_dpdm.value(); }
  double getOptimPenaltyEnergy() const { return optim_penalty_energy.value(); }
  double getOptimPenaltyVariation() const { return optim_penalty_variation.value(); }

  const std::string& getOutputDirectory() const { return output_directory.value(); }
  const std::vector<OutputType>& getOutputObservables() const { return output_observables.value(); }
  size_t getOutputTimestepStride() const { return output_timestep_stride.value(); }
  size_t getOutputOptimizationStride() const { return output_optimization_stride.value(); }
  RunType getRuntype() const { return runtype.value(); }
  bool getUseMatFree() const { return usematfree.value(); }
  LinearSolverType getLinearSolverType() const { return linearsolver_type.value(); }
  size_t getLinearSolverMaxiter() const { return linearsolver_maxiter.value(); }
  TimeStepperType getTimestepperType() const { return timestepper_type.value(); }
  int getRandSeed() const { return rand_seed.value(); }

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
