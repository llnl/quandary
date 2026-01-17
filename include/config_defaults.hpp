#pragma once

#include <cstddef>
#include <string>

#include "defs.hpp"

/**
 * @brief Centralized configuration defaults for Quandary.
 *
 * This namespace contains all default values used in the configuration system.
 * Simple compile-time constants are defined here, while computed defaults
 * that depend on other settings are provided through functions.
 */
namespace ConfigDefaults {

// General options
const double ROTFREQ = 0.0; ///< Default rotational wave approximation frequency (GHz)
const double SELFKERR = 0.0; ///< Default self-kerr frequency (GHz)
const double CROSSKERR = 0.0; ///< Default cross-kerr frequency (GHz)
const double JKL = 0.0; ///< Default dipole-dipole coupling frequency (GHz)
const DecoherenceType DECOHERENCE_TYPE = DecoherenceType::NONE; ///< Default decoherence type enum
const double DECAY_TIME = 0.0; ///< Default decay time
const double DEPHASE_TIME = 0.0; ///< Default dephase time

// Optimization options
const bool CONTROL_ENFORCE_BC = true; ///< Default control pulse boundary conditions enforcement
const ControlType CONTROL_TYPE = ControlType::BSPLINE; ///< Default control parameterization type
const size_t CONTROL_SPLINE_COUNT = 10; ///< Default number of B-spline basis functions
const ControlInitializationType CONTROL_INIT_TYPE = ControlInitializationType::CONSTANT; ///< Default control initialization amplitude
const double CONTROL_INIT_AMPLITUDE = 0.0; ///< Default control initialization amplitude
const double CONTROL_INIT_PHASE = 0.0; ///< Default control initialization phase

const double CONTROL_BOUND = 1e12; ///< Default amplitude bound for control pulses
const double CARRIER_FREQ = 0.0; ///< Default carrier frequency

const TargetType OPTIM_TARGET = TargetType::NONE; ///< Default optimization target: NONE
const GateType GATE_TYPE = GateType::NONE; ///< Default gate type
const double GATE_ROT_FREQ = 0.0; ///< Default gate rotational frequency
const ObjectiveType OPTIM_OBJECTIVE = ObjectiveType::JTRACE; ///< Default objective function
const double OPTIM_WEIGHT = 1.0; ///< Default optimization weight per initial condition

const double OPTIM_TIKHONOV_COEFF = 1e-4; ///< Default Tikhonov regularization coefficient
const bool OPTIM_TIKHONOV_USE_X0 = false; ///< Default Tikhonov regularization type
const double OPTIM_TOL_GRAD_ABS = 1e-8; ///< Default absolute gradient tolerance
const double OPTIM_TOL_GRAD_REL = 1e-4; ///< Default relative gradient tolerance
const double OPTIM_TOL_FINALCOST = 1e-8; ///< Default final time cost tolerance
const double OPTIM_TOL_INFIDELITY = 1e-5; ///< Default infidelity tolerance
const size_t OPTIM_MAXITER = 200; ///< Default maximum optimization iterations

const double OPTIM_PENALTY_LEAKAGE = 0.0; ///< Default first integral penalty coefficient
const double OPTIM_PENALTY_WEIGHTEDCOST = 0.0; ///< Default weighted cost penalty coefficient
const double OPTIM_PENALTY_WEIGHTEDCOST_WIDTH = 0.5; ///< Default weighted cost penalty width
const double OPTIM_PENALTY_DPDM = 0.0; ///< Default second derivative penalty coefficient
const double OPTIM_PENALTY_ENERGY = 0.0; ///< Default energy penalty coefficient
const double OPTIM_PENALTY_VARIATION = 0.01; ///< Default amplitude variation penalty coefficient

const std::string DATADIR = "./data_out"; ///< Default output directory

const size_t OUTPUT_TIMESTEP_STRIDE = 1; ///< Default output frequency
const size_t OUTPUT_OPTIMIZATION_STRIDE = 10; ///< Default optimization monitoring frequency

const RunType RUNTYPE = RunType::SIMULATION; ///< Default run type
const bool USEMATFREE = false; ///< Default matrix-free solver setting
const LinearSolverType LINEARSOLVER_TYPE = LinearSolverType::GMRES; ///< Default linear solver type
const size_t LINEARSOLVER_MAXITER = 10; ///< Default linear solver max iterations
const TimeStepperType TIMESTEPPER_TYPE = TimeStepperType::IMR; ///< Default time stepper type
const int RAND_SEED = 1; ///< Default random seed

} // namespace ConfigDefaults

/* Structs to group certain configuration settings */

/**
 * @brief Settings for initial conditions. Required, no defaults.
 */
struct InitialConditionSettings {
  InitialConditionType type; ///< Type of initial condition

  // Optional settings - populate based on type
  std::optional<std::string> filename; ///< For FROMFILE: File to read initial condition from
  std::optional<std::vector<size_t>> levels; ///< For PRODUCT_STATE: Quantum level for each oscillator
  std::optional<std::vector<size_t>> subsystem; ///< For ENSEMBLE, DIAGONAL, BASIS: Oscillator IDs
};

/**
 * @brief Settings for optimization targets with defaults
 */
struct OptimTargetSettings {
  TargetType type = ConfigDefaults::OPTIM_TARGET; ///< Type of optimization target
  std::optional<GateType> gate_type = std::nullopt; ///< For GATE: Type of the gate
  std::optional<std::vector<double>> gate_rot_freq = std::nullopt; ///< For GATE: Gate rotation frequencies for each oscillator
  std::optional<std::vector<size_t>> levels = std::nullopt; ///< For STATE: Level occupations for each oscillator
  std::optional<std::string> filename = std::nullopt; ///< For GATE or STATE: File path to target gate or state
};

/**
 * @brief Settings for control parameterizations with defaults.
 */
struct ControlParameterizationSettings {
  ControlType type = ConfigDefaults::CONTROL_TYPE; ///< Type of control parameterization
  std::optional<size_t> nspline = ConfigDefaults::CONTROL_SPLINE_COUNT; ///< Number of basis functions in this parameterization
  std::optional<double> tstart = std::nullopt; ///< Start time of the control parameterization
  std::optional<double> tstop = std::nullopt; ///< Stop time of the control parameterization
  std::optional<double> scaling = std::nullopt; ///< Amplitude scaling factor, only for BSPLINEAMP

};

/**
 * @brief Settings for control initialization with defaults.
 */
struct ControlInitializationSettings {
  ControlInitializationType type = ConfigDefaults::CONTROL_INIT_TYPE; ///< Initialization type
  std::optional<double> amplitude = ConfigDefaults::CONTROL_INIT_AMPLITUDE; ///< Initial control pulse amplitude
  std::optional<double> phase = std::nullopt; ///< Initial control pulse phase
  std::optional<std::string> filename = std::nullopt; ///< Filename for FILE type
};
