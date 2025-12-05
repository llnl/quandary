#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "defs.hpp"

/**
 * @brief Centralized configuration defaults for Quandary.
 *
 * This namespace contains all default values used in the configuration system.
 * Simple compile-time constants are defined here, while computed defaults
 * that depend on other settings are provided through functions.
 */
namespace ConfigDefaults {

// Control system defaults
const double CONTROL_BOUND = 10000.0; ///< Default amplitude bound for control pulses
const double CARRIER_FREQ = 0.0; ///< Default carrier frequency
const size_t SPLINE_COUNT = 10; ///< Default number of B-spline basis functions
const double CONTROL_INIT_AMPLITUDE = 0.0; ///< Default control initialization amplitude
const double CONTROL_INIT_PHASE = 0.0; ///< Default control initialization phase

// Optimization defaults
const double OPTIM_REGUL = 1e-4; ///< Default Tikhonov regularization coefficient
const double OPTIM_ATOL = 1e-8; ///< Default absolute gradient tolerance
const double OPTIM_RTOL = 1e-4; ///< Default relative gradient tolerance
const double OPTIM_FTOL = 1e-8; ///< Default final time cost tolerance
const double OPTIM_INFTOL = 1e-5; ///< Default infidelity tolerance
const size_t OPTIM_MAXITER = 200; ///< Default maximum optimization iterations

// Penalty defaults
const double PENALTY = 0.0; ///< Default first integral penalty coefficient
const double PENALTY_PARAM = 0.5; ///< Default Gaussian variance parameter
const double PENALTY_DPDM = 0.0; ///< Default second derivative penalty coefficient
const double PENALTY_ENERGY = 0.0; ///< Default energy penalty coefficient
const double PENALTY_VARIATION = 0.01; ///< Default amplitude variation penalty coefficient

// Random initialization defaults
const double RANDOM_AMPLITUDE = 0.1; ///< Default amplitude for random control initialization

// Time stepping defaults
const size_t NTIME = 1000; ///< Default number of time steps
const double DT = 0.1; ///< Default time step size (ns)

// Output defaults
const size_t OUTPUT_FREQUENCY = 1; ///< Default output frequency
const size_t OPTIM_MONITOR_FREQUENCY = 10; ///< Default optimization monitoring frequency
const size_t LINEARSOLVER_MAXITER = 10; ///< Default linear solver max iterations

// Boolean defaults
const bool CONTROL_ENFORCE_BC = true; ///< Default boundary conditions enforcement
const bool OPTIM_REGUL_TIK0 = false; ///< Default Tikhonov regularization type
const bool USEMATFREE = false; ///< Default matrix-free solver setting

// Enum defaults
const GateType GATE_TYPE = GateType::NONE; ///< Default gate type
const LindbladType COLLAPSE_TYPE_ENUM = LindbladType::NONE; ///< Default collapse type enum
const ObjectiveType OPTIM_OBJECTIVE = ObjectiveType::JFROBENIUS; ///< Default objective function
const RunType RUNTYPE = RunType::SIMULATION; ///< Default run type
const LinearSolverType LINEARSOLVER_TYPE = LinearSolverType::GMRES; ///< Default linear solver type
const TimeStepperType TIMESTEPPER_TYPE = TimeStepperType::IMR; ///< Default time stepper type

// String defaults
const std::string DATADIR = "./data_out"; ///< Default output directory

// Functions for computed defaults that depend on other settings
// Note: These are declared here but defined in config.cpp to avoid circular dependencies


} // namespace ConfigDefaults
