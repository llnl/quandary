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
const double DT = 0.1; ///< Default time step size (ns)
const size_t NTIME = 1000; ///< Default number of time steps
const double SELFKERR = 0.0; ///< Default self-kerr frequency (GHz)
const double CROSSKERR = 0.0; ///< Default cross-kerr frequency (GHz)
const double JKL = 0.0; ///< Default dipole-dipole coupling frequency (GHz)
const LindbladType COLLAPSE_TYPE_ENUM = LindbladType::NONE; ///< Default collapse type enum
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

const TargetType OPTIM_TARGET = TargetType::PRODUCT_STATE; ///< Default optimization target: product state transfer
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

const size_t OUTPUT_FREQUENCY = 1; ///< Default output frequency
const size_t OPTIM_MONITOR_FREQUENCY = 10; ///< Default optimization monitoring frequency

const RunType RUNTYPE = RunType::SIMULATION; ///< Default run type
const bool USEMATFREE = false; ///< Default matrix-free solver setting
const LinearSolverType LINEARSOLVER_TYPE = LinearSolverType::GMRES; ///< Default linear solver type
const size_t LINEARSOLVER_MAXITER = 10; ///< Default linear solver max iterations
const TimeStepperType TIMESTEPPER_TYPE = TimeStepperType::IMR; ///< Default time stepper type
const int RAND_SEED = 1; ///< Default random seed

} // namespace ConfigDefaults
