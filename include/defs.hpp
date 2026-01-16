#include <map>
#include <optional>
#include <string>
#include <vector>

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
enum class DecoherenceType {
  NONE,    ///< No Lindblad operators (closed system)
  DECAY,   ///< Decay operators only
  DEPHASE, ///< Dephasing operators only
  BOTH     ///< Both decay and dephasing operators
};

const std::map<std::string, DecoherenceType> DECOHERENCE_TYPE_MAP = {
    {"none", DecoherenceType::NONE},
    {"decay", DecoherenceType::DECAY},
    {"dephase", DecoherenceType::DEPHASE},
    {"both", DecoherenceType::BOTH}
};

/**
 * @brief Available types of initial conditions that are propagated through the quantum dynamics.
 *
 * Defines how the initial quantum state is specified and prepared
 * for simulation or optimization.
 */
enum class InitialConditionType {
  FROMFILE,      ///< Read initial condition from file
  PRODUCT_STATE, ///< Product state initial condition
  ENSEMBLE,      ///< Ensemble of states
  DIAGONAL,      ///< Diagonal density matrix
  BASIS,         ///< Basis state
  THREESTATES,   ///< Three-state system
  NPLUSONE,      ///< N+1 state system
  PERFORMANCE    ///< Performance test configuration
};

const std::map<std::string, InitialConditionType> INITCOND_TYPE_MAP = {
    {"file", InitialConditionType::FROMFILE},
    {"state", InitialConditionType::PRODUCT_STATE},
    {"ensemble", InitialConditionType::ENSEMBLE},
    {"diagonal", InitialConditionType::DIAGONAL},
    {"basis", InitialConditionType::BASIS},
    {"3states", InitialConditionType::THREESTATES},
    {"nplus1", InitialConditionType::NPLUSONE},
    {"performance", InitialConditionType::PERFORMANCE}
};

/**
 * @brief Types of optimization targets for quantum control.
 *
 * Defines the target quantum state or operation for optimization.
 */
enum class TargetType {
  NONE,          ///< No target specified (no optimization)
  GATE,          ///< Gate optimization: \f$\rho_{\text{target}} = V\rho(0) V^\dagger\f$ for V either read from file or chosen from default set of gates
  STATE,         ///< State preparation: Either read from file, or \f$\rho_{\text{target}} = e_m e_m^\dagger\f$ for some integer \f$m\f$
};

const std::map<std::string, TargetType> TARGET_TYPE_MAP = {
    {"none", TargetType::NONE},
    {"gate", TargetType::GATE},
    {"state", TargetType::STATE},
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
 * @brief Types of control parameterizations.
 *
 * Defines how control pulses are parameterized for optimization and simulation.
 */
enum class ControlType {
  NONE,       ///< Non-controllable
  BSPLINE,    ///< Control pulses are parameterized with 2nd order BSpline basis functions with carrier waves
  BSPLINEAMP, ///< Paramerizes only the amplitudes of the control pulse with 2nd order BSpline basis functions 
  BSPLINE0    ///< Control pulses are parameterized with Zeroth order Bspline (piece-wise constant)
};

const std::map<std::string, ControlType> CONTROL_TYPE_MAP = {
    {"none", ControlType::NONE},
    {"spline", ControlType::BSPLINE},
    {"spline_amplitude", ControlType::BSPLINEAMP},
    {"spline0", ControlType::BSPLINE0}
};

/**
 * @brief Types of control initializations 
 */
enum class ControlInitializationType {
  CONSTANT, ///< Constant
  RANDOM,   ///< Random
  FILE,     ///< From file
};

const std::map<std::string, ControlInitializationType> CONTROL_INITIALIZATION_TYPE_MAP = {
    {"constant", ControlInitializationType::CONSTANT},
    {"random", ControlInitializationType::RANDOM},
    {"file", ControlInitializationType::FILE},
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


/* Structs to group certain configuration settings */

/**
 * @brief Settings for initial conditions
 */
struct InitialConditionSettings {
  InitialConditionType type; ///< Type of initial condition

  // Optional settings - populate based on type
  std::optional<std::string> filename; ///< For FROMFILE: File to read initial condition from
  std::optional<std::vector<size_t>> levels; ///< For PRODUCT_STATE: Quantum level for each oscillator
  std::optional<std::vector<size_t>> subsystem; ///< For ENSEMBLE, DIAGONAL, BASIS: Oscillator IDs
};

/**
 * @brief Settings for optimization targets
 */
struct OptimTargetSettings {
  TargetType type; ///< Type of optimization target (NONE, GATE, STATE)

  // Optional settings - populate based on type
  std::optional<GateType> gate_type; ///< For GATE: Type of the gate
  std::optional<std::vector<double>> gate_rot_freq; ///< For GATE: Gate rotation frequencies for each oscillator 
  std::optional<std::vector<size_t>> levels; ///< For STATE: Level occupations for each oscillator
  std::optional<std::string> filename; ///< For GATE or STATE: File path to target gate or state
};

/**
 * @brief Settings for control parameterizations.
 *
 * Defines a controllable parameterization for an oscillator with corresponding starting and finish times. 
 * Which fields are used depends on the control type.
 */
struct ControlParameterizationSettings {
  ControlType type; ///< Type of control parameterization

  // Common fields for B-spline types (BSPLINE, BSPLINEAMP, BSPLINE0)
  std::optional<size_t> nspline; ///< Number of basis functions in this parameterization
  std::optional<double> tstart; ///< Start time of the control parameterization
  std::optional<double> tstop; ///< Stop time of the control parameterization
  std::optional<double> scaling; ///< Amplitude scaling factor, only for BSPLINEAMP

};

/**
 * @brief Settings for control initialization.
 */
struct ControlInitializationSettings {
  ControlInitializationType type; ///< Initialization type
  std::optional<double> amplitude; ///< Initial control pulse amplitude
  std::optional<double> phase; ///< Initial control pulse phase
  std::optional<std::string> filename; ///< Filename for FILE type
};

