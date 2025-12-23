#pragma once

#include <petsc.h>

#include <cstddef>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <toml++/toml.hpp>
#include <vector>

#include "cfgparser.hpp"
#include "defs.hpp"
#include "mpi_logger.hpp"

/**
 * @brief Final validated configuration class.
 *
 * Contains only validated, typed configuration parameters. All fields are required
 * and have been validated with defaults set. This class is immutable after construction.
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
  LindbladType collapse_type; ///< Switch between Schroedinger and Lindblad solver
  std::vector<double> decay_time; ///< Time of decay collapse operation (T1) per oscillator (for Lindblad solver)
  std::vector<double> dephase_time; ///< Time of dephase collapse operation (T2) per oscillator (for Lindblad solver)
  size_t n_initial_conditions; ///< Number of initial conditions
  InitialCondition initial_condition; ///< Initial condition configuration
  std::vector<std::vector<PiPulseSegment>> apply_pipulse; ///< Apply a pi-pulse to oscillator with specified parameters
  std::optional<std::string> hamiltonian_file_Hsys; ///< File to read the system Hamiltonian from
  std::optional<std::string> hamiltonian_file_Hc; ///< File to read the control Hamiltonian from

  // Optimization options
  bool control_enforceBC; ///< Decide whether control pulses should start and end at zero
  std::optional<std::string> control_initialization_file; ///< Global control initialization file for all oscillators
  std::vector<std::vector<ControlSegment>> control_segments; ///< Control segments for each oscillator
  std::vector<std::vector<ControlSegmentInitialization>>
      control_initializations; ///< Control initializations for each oscillator
  std::vector<std::vector<double>> control_bounds; ///< Control bounds for each oscillator
  std::vector<std::vector<double>> carrier_frequencies; ///< Carrier frequencies for each oscillator
  OptimTargetSettings optim_target; ///< Grouped optimization target configuration
  std::vector<double> gate_rot_freq; ///< Frequency of rotation of the target gate, for each oscillator (GHz)
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
  std::string datadir; ///< Directory for output files
  std::vector<OutputType> output_type; ///< Specify the desired output types. 
  size_t output_frequency; ///< Output frequency in the time domain: write output every <num> time-step
  size_t optim_monitor_frequency; ///< Frequency of writing output during optimization iterations
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

  ~Config();

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
  LindbladType getCollapseType() const { return collapse_type; }
  const std::vector<double>& getDecayTime() const { return decay_time; }
  const std::vector<double>& getDephaseTime() const { return dephase_time; }
  size_t getNInitialConditions() const { return n_initial_conditions; }
  const InitialCondition& getInitialCondition() const { return initial_condition; }
  const std::vector<std::vector<PiPulseSegment>>& getApplyPiPulses() const { return apply_pipulse; }
  const std::vector<PiPulseSegment>& getApplyPiPulse(size_t i_osc) const { return apply_pipulse[i_osc]; }
  const std::optional<std::string>& getHamiltonianFileHsys() const { return hamiltonian_file_Hsys; }
  const std::optional<std::string>& getHamiltonianFileHc() const { return hamiltonian_file_Hc; }

  const std::vector<ControlSegment>& getControlSegments(size_t i_osc) const { return control_segments[i_osc]; }
  bool getControlEnforceBC() const { return control_enforceBC; }
  const std::vector<ControlSegmentInitialization>& getControlInitializations(size_t i_osc) const {
    return control_initializations[i_osc];
  }
  const std::optional<std::string> getControlInitializationFile() const { return control_initialization_file; }
  const std::vector<double>& getControlBounds(size_t i_osc) const { return control_bounds[i_osc]; }
  double getControlBound(size_t i_osc, size_t i_seg) const { return control_bounds[i_osc][i_seg]; }
  const std::vector<double>& getCarrierFrequencies(size_t i_osc) const { return carrier_frequencies[i_osc]; }
  double getCarrierFrequency(size_t i_osc, size_t i_seg) const { return carrier_frequencies[i_osc][i_seg]; }
  const OptimTargetSettings& getOptimTarget() const { return optim_target; }
  const std::vector<double>& getGateRotFreq() const { return gate_rot_freq; }
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

  const std::string& getDataDir() const { return datadir; }
  const std::vector<OutputType>& getOutputType() const { return output_type; }
  size_t getOutputFrequency() const { return output_frequency; }
  size_t getOptimMonitorFrequency() const { return optim_monitor_frequency; }
  RunType getRuntype() const { return runtype; }
  bool getUseMatFree() const { return usematfree; }
  LinearSolverType getLinearSolverType() const { return linearsolver_type; }
  size_t getLinearSolverMaxiter() const { return linearsolver_maxiter; }
  TimeStepperType getTimestepperType() const { return timestepper_type; }
  int getRandSeed() const { return rand_seed; }

 private:
  void finalize();
  void validate() const;

  size_t computeNumInitialConditions() const;
  void setRandSeed(int rand_seed_);

  // Table validation helper
  void validateTableKeys(const toml::table& table, const std::set<std::string>& allowed_keys,
                         const std::string& table_name) const;

  // Conversion helper methods
  template <typename T>
  std::vector<std::vector<T>> parseOscillatorSettings(const toml::array& array_of_tables, size_t num_entries,
                                                      std::vector<T> default_values, const std::string& field_name,
                                                      const std::set<std::string>& allowed_keys,
                                                      const std::string& table_name) const;

  InitialCondition parseInitialCondition(std::optional<InitialConditionType> opt_type,
                                         const std::optional<std::string>& filename,
                                         const std::optional<std::vector<size_t>>& levels,
                                         const std::optional<std::vector<size_t>>& osc_IDs) const;

  void addPiPulseSegment(std::vector<std::vector<PiPulseSegment>>& apply_pipulse, size_t oscilID, double tstart,
                         double tstop, double amp) const;

  ControlSegment parseControlSegment(const toml::table& table) const;
  std::vector<std::vector<ControlSegment>> parseControlSegments(const toml::array& array_of_tables,
                                                                size_t num_entries) const;

  std::vector<std::vector<ControlSegmentInitialization>> parseControlInitializations(
      const toml::array& array_of_tables, size_t num_entries, std::optional<std::string>& control_init_file) const;

  OptimTargetSettings parseOptimTarget(TargetType type, const std::optional<GateType>& gate_type,
                                       const std::optional<std::string>& gate_file,
                                       const std::optional<std::vector<size_t>>& levels,
                                       const std::optional<std::string>& file) const;


  /**
   * @brief Parses coupling parameters from TOML table format
   *
   * Converts a TOML table of coupling parameters specified with string keys like "0-1", "1-4"
   * into a flat vector indexed by pair index. Pairs are ordered as (0,1), (0,2), ..., (0,n-1),
   * (1,2), (1,3), ..., (1,n-1), ..., (n-2,n-1).
   *
   * @param toml TOML table containing the configuration
   * @param key Name of the parameter (e.g., "Jkl", "crosskerr")
   * @param num_osc Number of oscillators
   * @param default_value Default value for unspecified couplings
   * @return Vector of coupling values indexed by pair
   */
  std::vector<double> parseCouplingParameters(const toml::table& toml, const std::string& key, size_t num_osc, double default_value) const;


  /**
   * @brief Prints coupling parameters in table format
   *
   * Converts a flat coupling vector back to the table format with "i-j" keys.
   * Only prints non-zero couplings.
   *
   * @param couplings Vector of coupling values indexed by pair
   * @param num_osc Number of oscillators
   * @return String in TOML table format
   */
  std::string printCouplingParameters(const std::vector<double>& couplings, size_t num_osc) const;

  // TODO cfg: delete these when .cfg format is removed.
  template <typename T>
  std::vector<std::vector<T>> parseOscillatorSettingsCfg(const std::optional<std::map<int, std::vector<T>>>& indexed,
                                                         size_t num_entries,
                                                         const std::vector<T>& default_values = {}) const;

  std::vector<std::vector<ControlSegment>> parseControlSegmentsCfg(
      const std::optional<std::map<int, std::vector<ControlSegmentData>>>& segments_opt) const;
  ControlSegment parseControlSegmentCfg(const ControlSegmentData& seg_config) const;
  std::vector<std::vector<ControlSegmentInitialization>> parseControlInitializationsCfg(
      const std::optional<std::map<int, std::vector<ControlInitializationData>>>& init_configs) const;
};
