#pragma once

#include <cstddef>
#include <string>

#include "config.hpp"
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
const size_t NTIME = 1000; ///< Default number of time steps
const double DT = 0.1; ///< Default time step size (ns)
const double SELFKERR = 0.0; ///< Default self-kerr frequency (GHz)
const double CROSSKERR = 0.0; ///< Default cross-kerr frequency (GHz)
const double JKL = 0.0; ///< Default dipole-dipole coupling frequency (GHz)
const LindbladType COLLAPSE_TYPE_ENUM = LindbladType::NONE; ///< Default collapse type enum
const double DECAY_TIME = 0.0; ///< Default decay time
const double DEPHASE_TIME = 0.0; ///< Default dephase time

// Optimization options
const bool CONTROL_ENFORCE_BC = true; ///< Default boundary conditions enforcement
const size_t CONTROL_SEG_SPLINE_COUNT = 10; ///< Default number of B-spline basis functions
const double CONTROL_SEG_TSTART = 0.0; ///< Default control segment start time
const double CONTROL_INIT_AMPLITUDE = 0.0; ///< Default control initialization amplitude
const double CONTROL_INIT_PHASE = 0.0; ///< Default control initialization phase
const double CONTROL_INIT_RANDOM_AMPLITUDE = 0.1; ///< Default amplitude for random control initialization

const double CONTROL_BOUND = 10000.0; ///< Default amplitude bound for control pulses
const double CARRIER_FREQ = 0.0; ///< Default carrier frequency

const PureOptimTarget OPTIM_TARGET = PureOptimTarget{}; ///< Default optimization target: pure state transfer
const GateType GATE_TYPE = GateType::NONE; ///< Default gate type
const double GATE_ROT_FREQ = 0.0; ///< Default gate rotational frequency
const ObjectiveType OPTIM_OBJECTIVE = ObjectiveType::JFROBENIUS; ///< Default objective function
const double OPTIM_WEIGHT = 1.0; ///< Default optimization weight per initial condition

const double OPTIM_REGUL = 1e-4; ///< Default Tikhonov regularization coefficient
const double OPTIM_ATOL = 1e-8; ///< Default absolute gradient tolerance
const double OPTIM_RTOL = 1e-4; ///< Default relative gradient tolerance
const double OPTIM_FTOL = 1e-8; ///< Default final time cost tolerance
const double OPTIM_INFTOL = 1e-5; ///< Default infidelity tolerance
const size_t OPTIM_MAXITER = 200; ///< Default maximum optimization iterations
const OptimTolerance OPTIM_TOLERANCE = OptimTolerance{
    OPTIM_ATOL,
    OPTIM_RTOL,
    OPTIM_FTOL,
    OPTIM_INFTOL,
    OPTIM_MAXITER
}; ///< Default optimization tolerance settings

const double PENALTY = 0.0; ///< Default first integral penalty coefficient
const double PENALTY_PARAM = 0.5; ///< Default Gaussian variance parameter
const double PENALTY_DPDM = 0.0; ///< Default second derivative penalty coefficient
const double PENALTY_ENERGY = 0.0; ///< Default energy penalty coefficient
const double PENALTY_VARIATION = 0.01; ///< Default amplitude variation penalty coefficient
const OptimPenalty OPTIM_PENALTY = OptimPenalty{
    PENALTY,
    PENALTY_PARAM,
    PENALTY_DPDM,
    PENALTY_ENERGY,
    PENALTY_VARIATION
}; ///< Default optimization penalty settings

const bool OPTIM_REGUL_TIK0 = false; ///< Default Tikhonov regularization type

const std::string DATADIR = "./data_out"; ///< Default output directory

const size_t OUTPUT_FREQUENCY = 1; ///< Default output frequency
const size_t OPTIM_MONITOR_FREQUENCY = 10; ///< Default optimization monitoring frequency

const RunType RUNTYPE = RunType::SIMULATION; ///< Default run type
const bool USEMATFREE = false; ///< Default matrix-free solver setting
const LinearSolverType LINEARSOLVER_TYPE = LinearSolverType::GMRES; ///< Default linear solver type
const size_t LINEARSOLVER_MAXITER = 10; ///< Default linear solver max iterations
const TimeStepperType TIMESTEPPER_TYPE = TimeStepperType::IMR; ///< Default time stepper type
const int RAND_SEED = -1; ///< Default random seed

} // namespace ConfigDefaults
