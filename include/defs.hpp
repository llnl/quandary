#include <map>
#include <string>

/**
 * @file defs.hpp
 * @brief Core type definitions and enumerations for Quandary quantum optimal control.
 *
 * This file contains fundamental type definitions, enumeration classes, and constants
 * used throughout the Quandary quantum optimal control framework. It defines solver
 * types, target specifications, objective functions, and control parameterizations.
 */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SQR(x) (x)*(x)

#pragma once

/**
 * @brief Available Lindblad operator types for open quantum systems, or NONE for closed quantum systems.
 *
 * For open quantum systems, this defines the types of dissipation and decoherence operators 
 * that can be applied to the quantum system in Lindblad master equation. 
 * 
 * @note If this is NONE, the quantum system is considered closed, solving Schroedinger's 
 * equation rather than Lindblad's master equation.
 */
enum class LindbladType {
  NONE,    ///< No Lindblad operators (closed system)
  DECAY,   ///< Decay operators only
  DEPHASE, ///< Dephasing operators only
  BOTH     ///< Both decay and dephasing operators
};

const std::map<std::string, LindbladType> LINDBLAD_TYPE_MAP = {
    {"none", LindbladType::NONE},
    {"decay", LindbladType::DECAY},
    {"dephase", LindbladType::DEPHASE},
    {"both", LindbladType::BOTH}
};

/**
 * @brief Available types of initial conditions that are propagated through the quantum dynamics.
 *
 * Defines how the initial quantum state is specified and prepared
 * for simulation or optimization.
 */
enum class InitialConditionType {
  FROMFILE,    ///< Read initial condition from file
  PURE,        ///< Pure state initial condition
  ENSEMBLE,    ///< Ensemble of states
  DIAGONAL,    ///< Diagonal density matrix
  BASIS,       ///< Basis state
  THREESTATES, ///< Three-state system
  NPLUSONE,    ///< N+1 state system
  PERFORMANCE  ///< Performance test configuration
};

const std::map<std::string, InitialConditionType> INITCOND_TYPE_MAP = {
    {"file", InitialConditionType::FROMFILE},
    {"pure", InitialConditionType::PURE},
    {"ensemble", InitialConditionType::ENSEMBLE},
    {"diagonal", InitialConditionType::DIAGONAL},
    {"basis", InitialConditionType::BASIS},
    {"3states", InitialConditionType::THREESTATES},
    {"nplus1", InitialConditionType::NPLUSONE},
    {"performance", InitialConditionType::PERFORMANCE}
};

/**
 * @brief Initial condition configuration.
 */
struct InitialCondition {
  InitialConditionType type; ///< Type of initial condition

  // Optional fields - populate based on type
  std::optional<std::string> filename; ///< For FROMFILE: File to read initial condition from
  std::optional<std::vector<size_t>> levels; ///< For PURE: Quantum level for each oscillator
  std::optional<std::vector<size_t>> osc_IDs; ///< For ENSEMBLE, DIAGONAL, BASIS: Oscillator IDs
};

/**
 * @brief Types of optimization targets for quantum control.
 *
 * Defines the target quantum state or operation for optimization.
 */
enum class TargetType {
  GATE,      ///< Gate optimization: \f$\rho_{\text{target}} = V\rho(0) V^\dagger\f$
  PURE,      ///< Pure state preparation: \f$\rho_{\text{target}} = e_m e_m^\dagger\f$ for some integer \f$m\f$
  FROMFILE   ///< Target state read from file, vectorized density matrix format
};

const std::map<std::string, TargetType> TARGET_TYPE_MAP = {
    {"gate", TargetType::GATE},
    {"pure", TargetType::PURE},
    {"file", TargetType::FROMFILE}
};

/**
 * @brief Types of objective functions for quantum control optimization.
 *
 * Defines different metrics for measuring the quality of quantum control.
 */
enum class ObjectiveType {
  JFROBENIUS, ///< Weighted Frobenius norm: \f$\frac{1}{2} \frac{\|\rho_{\text{target}} - \rho(T)\|_F^2}{w}\f$, where \f$w\f$ = purity of \f$\rho_{\text{target}}\f$
  JTRACE,     ///< Weighted Hilbert-Schmidt overlap: \f$1 - \frac{\text{Tr}(\rho_{\text{target}}^\dagger \rho(T))}{w}\f$, where \f$w\f$ = purity of \f$\rho_{\text{target}}\f$
  JMEASURE    ///< Pure state measurement: \f$\text{Tr}(O_m \rho(T))\f$ for observable \f$O_m\f$
};

const std::map<std::string, ObjectiveType> OBJECTIVE_TYPE_MAP = {
    {"jfrobenius", ObjectiveType::JFROBENIUS},
    {"jtrace", ObjectiveType::JTRACE},
    {"jmeasure", ObjectiveType::JMEASURE}
};

/**
 * @brief Available types for solving linear systems at each time step.
 *
 * Defines the numerical methods used to solve linear systems
 * arising in quantum dynamics simulations at each time step.
 */
enum class LinearSolverType{
  GMRES,   ///< Uses Petsc's GMRES solver (default)
  NEUMANN  ///< Uses Neuman power iterations
};

const std::map<std::string, LinearSolverType> LINEAR_SOLVER_TYPE_MAP = {
    {"gmres", LinearSolverType::GMRES},
    {"neumann", LinearSolverType::NEUMANN}
};

/**
 * @brief Types of execution modes.
 *
 * Defines what type of computation should be performed.
 */
enum class RunType {
  SIMULATION,   ///< Runs one simulation to compute the objective function (forward)
  GRADIENT,     ///< Runs a simulation followed by the adjoint for gradient computation (forward & backward)
  OPTIMIZATION, ///< Runs optimization iterations
  EVALCONTROLS, ///< Only evaluates the current control pulses (no simulation)
  NONE          ///< Don't run anything
};

const std::map<std::string, RunType> RUN_TYPE_MAP = {
    {"simulation", RunType::SIMULATION},
    {"gradient", RunType::GRADIENT},
    {"optimization", RunType::OPTIMIZATION},
    {"evalcontrols", RunType::EVALCONTROLS},
    {"none", RunType::NONE}
};

/**
 * @brief Types of control parameterizations for quantum control pulses.
 *
 * Defines how control pulses are parameterized for optimization and simulation.
 */
enum class ControlType {
  NONE,       ///< Non-controllable
  BSPLINE,    ///< Control pulses are parameterized with 2nd order BSpline basis functions with carrier waves
  BSPLINEAMP, ///< Paramerizes only the amplitudes of the control pulse with 2nd order BSpline basis functions 
  STEP,       ///< Control parameter is the width of a step function for a given amplitude
  BSPLINE0    ///< Control pulses are parameterized with Zeroth order Bspline (piece-wise constant)
};

const std::map<std::string, ControlType> CONTROL_TYPE_MAP = {
    {"none", ControlType::NONE},
    {"spline", ControlType::BSPLINE},
    {"spline_amplitude", ControlType::BSPLINEAMP},
    {"step", ControlType::STEP},
    {"spline0", ControlType::BSPLINE0}
};

/**
 * @brief Types of control initializations per segment
 */
enum class ControlSegmentInitType {
  CONSTANT, ///< Constant
  RANDOM,   ///< Random
  FILE,     ///< From file
};

const std::map<std::string, ControlSegmentInitType> CONTROL_SEGMENT_INIT_TYPE_MAP = {
    {"constant", ControlSegmentInitType::CONSTANT},
    {"random", ControlSegmentInitType::RANDOM},
    {"file", ControlSegmentInitType::FILE},
};

/**
 * @brief Types of time-stepping methods for evolving quantum states.
 *
 * Defines the numerical methods used to evolve quantum states in time.
 */
enum class TimeStepperType {
  IMR,   ///< Implicit Midpoint Rule (2nd order)
  IMR4,  ///< Implicit Midpoint Rule with 4th order extrapolation
  IMR8,  ///< Implicit Midpoint Rule with 8th order extrapolation
  EE,    ///< Explicit Euler (1st order)
};

const std::map<std::string, TimeStepperType> TIME_STEPPER_TYPE_MAP = {
    {"imr", TimeStepperType::IMR},
    {"imr4", TimeStepperType::IMR4},
    {"imr8", TimeStepperType::IMR8},
    {"ee", TimeStepperType::EE}
};

/**
 * @brief Types of quantum gates used in quantum control.
 */
enum class GateType {
  NONE,     ///< No gate
  XGATE,    ///< X gate (Pauli-X)
  YGATE,    ///< Y gate (Pauli-Y)
  ZGATE,    ///< Z gate (Pauli-Z)
  HADAMARD, ///< Hadamard gate
  CNOT,     ///< CNOT gate
  SWAP,     ///< SWAP gate
  SWAP_0Q,  ///< Multi-qubit SWAP gate with 0 qubit
  CQNOT,    ///< Multi-qubit CQNOT gate
  QFT,      ///< QFT gate
  FILE,     ///< Gate defined in a file
};

const std::map<std::string, GateType> GATE_TYPE_MAP = {
    {"none", GateType::NONE},
    {"xgate", GateType::XGATE},
    {"ygate", GateType::YGATE},
    {"zgate", GateType::ZGATE},
    {"hadamard", GateType::HADAMARD},
    {"cnot", GateType::CNOT},
    {"swap", GateType::SWAP},
    {"swap0q", GateType::SWAP_0Q},
    {"cqnot", GateType::CQNOT},
    {"qft", GateType::QFT},
    {"file", GateType::FILE}
};

/**
 * @brief Types of output files to be written
 */
enum class OutputType {
  EXPECTED_ENERGY,           ///< Expected energy
  EXPECTED_ENERGY_COMPOSITE, ///< Expected energy composite
  POPULATION,                ///< Population
  POPULATION_COMPOSITE,      ///< Population composite
  FULLSTATE,                 ///< Full state
};

const std::map<std::string, OutputType> OUTPUT_TYPE_MAP = {
  {"expectedenergy", OutputType::EXPECTED_ENERGY},
  {"expectedenergycomposite", OutputType::EXPECTED_ENERGY_COMPOSITE},
  {"population", OutputType::POPULATION},
  {"populationcomposite", OutputType::POPULATION_COMPOSITE},
  {"fullstate", OutputType::FULLSTATE},
};
