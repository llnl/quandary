#pragma once

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "defs.hpp"

struct PiPulseData {
  size_t oscil_id; ///< Oscillator ID
  double tstart; ///< Start time
  double tstop; ///< Stop time
  double amp; ///< Amplitude
};

struct ControlSegmentData {
  ControlType control_type; ///< Type of control segment
  std::vector<double> parameters; ///< Parameters for control segment
};

struct ControlInitializationData {
  ControlSegmentInitType init_seg_type; ///< Type of initialization per segment
  std::optional<double> amplitude; ///< Initial amplitude
  std::optional<double> phase; ///< Initial phase (optional)
  std::optional<std::string> filename; ///< Filename (for file init - mutually exclusive with above)
};

/**
 * @brief Configuration settings passed to Config constructor.
 *
 * Contains all optional configuration parameters that can be provided
 * to configure a Config object. Used by CfgParser to pass settings.
 */
struct ParsedConfigData {
  // General parameters
  std::optional<std::vector<size_t>> nlevels;
  std::optional<std::vector<size_t>> nessential;
  std::optional<size_t> ntime;
  std::optional<double> dt;
  std::optional<std::vector<double>> transfreq;
  std::optional<std::vector<double>> selfkerr;
  std::optional<std::vector<double>> crosskerr;
  std::optional<std::vector<double>> Jkl;
  std::optional<std::vector<double>> rotfreq;
  std::optional<LindbladType> collapse_type;
  std::optional<std::vector<double>> decay_time;
  std::optional<std::vector<double>> dephase_time;
  std::optional<InitialCondition> initialcondition;
  std::optional<std::vector<PiPulseData>> apply_pipulse;
  std::optional<std::string> hamiltonian_file_Hsys;
  std::optional<std::string> hamiltonian_file_Hc;

  // Control and optimization parameters
  std::optional<std::map<int, std::vector<ControlSegmentData>>> indexed_control_segments;
  std::optional<bool> control_enforceBC;
  std::optional<std::map<int, std::vector<ControlInitializationData>>> indexed_control_init;
  std::optional<std::map<int, std::vector<double>>> indexed_control_bounds;
  std::optional<std::map<int, std::vector<double>>> indexed_carrier_frequencies;
  std::optional<OptimTargetSettings> optim_target;
  std::optional<std::vector<double>> gate_rot_freq;
  std::optional<ObjectiveType> optim_objective;
  std::optional<std::vector<double>> optim_weights;
  std::optional<double> optim_atol;
  std::optional<double> optim_rtol;
  std::optional<double> optim_ftol;
  std::optional<double> optim_inftol;
  std::optional<size_t> optim_maxiter;
  std::optional<double> optim_regul;
  std::optional<double> optim_penalty;
  std::optional<double> optim_penalty_param;
  std::optional<double> optim_penalty_dpdm;
  std::optional<double> optim_penalty_energy;
  std::optional<double> optim_penalty_variation;
  std::optional<bool> optim_regul_tik0;
  std::optional<bool> optim_regul_interpolate; // deprecated

  // Output parameters
  std::optional<std::string> datadir;
  std::optional<std::map<int, std::vector<OutputType>>> indexed_output;
  std::optional<size_t> output_frequency;
  std::optional<size_t> optim_monitor_frequency;
  std::optional<RunType> runtype;
  std::optional<bool> usematfree;
  std::optional<LinearSolverType> linearsolver_type;
  std::optional<size_t> linearsolver_maxiter;
  std::optional<TimeStepperType> timestepper_type;
  std::optional<int> rand_seed;
};
