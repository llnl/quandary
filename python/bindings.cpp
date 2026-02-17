#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <sstream>

#include "config.hpp"
#include "quandary_core.hpp"

namespace nb = nanobind;

// Declare extension module "_quandary_impl" (must match name given in CMake)
NB_MODULE(_quandary_impl, m) {
  m.doc() = "Quandary Python bindings.\n\n"
            "Low-level C++ interface for quantum optimal control simulations.\n"
            "Most users should use the high-level interface from quandary.new instead.";

  // Exceptions
  static nb::exception<validators::ValidationError> validation_exc(
    m, "ValidationError", PyExc_RuntimeError);

  // Enums
  nb::enum_<DecoherenceType>(m, "DecoherenceType", "Decoherence model type for Lindblad master equation")
    .value("NONE", DecoherenceType::NONE, "No decoherence (Schrodinger equation)")
    .value("DECAY", DecoherenceType::DECAY, "T1 decay only")
    .value("DEPHASE", DecoherenceType::DEPHASE, "T2 dephasing only")
    .value("BOTH", DecoherenceType::BOTH, "Both T1 decay and T2 dephasing");

  nb::enum_<InitialConditionType>(m, "InitialConditionType", "Type of initial quantum state")
    .value("FROMFILE", InitialConditionType::FROMFILE, "Load initial state from file")
    .value("PRODUCT_STATE", InitialConditionType::PRODUCT_STATE, "Product state of specified levels")
    .value("ENSEMBLE", InitialConditionType::ENSEMBLE, "Ensemble of basis states")
    .value("DIAGONAL", InitialConditionType::DIAGONAL, "Diagonal density matrix")
    .value("BASIS", InitialConditionType::BASIS, "Computational basis states")
    .value("THREESTATES", InitialConditionType::THREESTATES, "Three-state initial condition")
    .value("NPLUSONE", InitialConditionType::NPLUSONE, "N+1 initial condition")
    .value("PERFORMANCE", InitialConditionType::PERFORMANCE, "Performance test initial condition");

  nb::enum_<TargetType>(m, "TargetType", "Type of optimization target")
    .value("NONE", TargetType::NONE, "No target (simulation only)")
    .value("GATE", TargetType::GATE, "Target unitary gate")
    .value("STATE", TargetType::STATE, "Target quantum state");

  nb::enum_<ObjectiveType>(m, "ObjectiveType", "Objective function for optimization")
    .value("JFROBENIUS", ObjectiveType::JFROBENIUS, "Frobenius norm objective")
    .value("JTRACE", ObjectiveType::JTRACE, "Trace fidelity objective")
    .value("JMEASURE", ObjectiveType::JMEASURE, "Measurement-based objective");

  nb::enum_<LinearSolverType>(m, "LinearSolverType", "Linear solver for time stepping")
    .value("GMRES", LinearSolverType::GMRES, "GMRES iterative solver")
    .value("NEUMANN", LinearSolverType::NEUMANN, "Neumann series approximation");

  nb::enum_<RunType>(m, "RunType", "Type of computation to perform")
    .value("SIMULATION", RunType::SIMULATION, "Forward simulation only")
    .value("GRADIENT", RunType::GRADIENT, "Compute objective and gradient")
    .value("OPTIMIZATION", RunType::OPTIMIZATION, "Run full optimization")
    .value("EVALCONTROLS", RunType::EVALCONTROLS, "Evaluate existing control pulses")
    .value("NONE", RunType::NONE, "No computation");

  nb::enum_<ControlType>(m, "ControlType", "Control pulse parameterization type")
    .value("NONE", ControlType::NONE, "No control")
    .value("BSPLINE", ControlType::BSPLINE, "B-spline parameterization")
    .value("BSPLINEAMP", ControlType::BSPLINEAMP, "B-spline with amplitude scaling")
    .value("BSPLINE0", ControlType::BSPLINE0, "B-spline with zero boundary conditions");

  nb::enum_<ControlInitializationType>(m, "ControlInitializationType", "Control pulse initialization method")
    .value("CONSTANT", ControlInitializationType::CONSTANT, "Constant amplitude initialization")
    .value("RANDOM", ControlInitializationType::RANDOM, "Random initialization")
    .value("FILE", ControlInitializationType::FILE, "Load from file");

  nb::enum_<TimeStepperType>(m, "TimeStepperType", "Time integration method")
    .value("IMR", TimeStepperType::IMR, "Implicit midpoint rule (2nd order)")
    .value("IMR4", TimeStepperType::IMR4, "4th order implicit method")
    .value("IMR8", TimeStepperType::IMR8, "8th order implicit method")
    .value("EE", TimeStepperType::EE, "Explicit Euler (1st order)");

  nb::enum_<GateType>(m, "GateType", "Predefined quantum gate types")
    .value("NONE", GateType::NONE, "No gate")
    .value("XGATE", GateType::XGATE, "Pauli X gate")
    .value("YGATE", GateType::YGATE, "Pauli Y gate")
    .value("ZGATE", GateType::ZGATE, "Pauli Z gate")
    .value("HADAMARD", GateType::HADAMARD, "Hadamard gate")
    .value("CNOT", GateType::CNOT, "Controlled-NOT gate")
    .value("SWAP", GateType::SWAP, "SWAP gate")
    .value("SWAP_0Q", GateType::SWAP_0Q, "SWAP 0-Q gate")
    .value("CQNOT", GateType::CQNOT, "Controlled-qNOT gate")
    .value("QFT", GateType::QFT, "Quantum Fourier Transform")
    .value("FILE", GateType::FILE, "Load gate from file");

  nb::enum_<OutputType>(m, "OutputType", "Observable output types")
    .value("EXPECTED_ENERGY", OutputType::EXPECTED_ENERGY, "Expected energy per oscillator")
    .value("EXPECTED_ENERGY_COMPOSITE", OutputType::EXPECTED_ENERGY_COMPOSITE, "Expected energy of composite system")
    .value("POPULATION", OutputType::POPULATION, "Level populations per oscillator")
    .value("POPULATION_COMPOSITE", OutputType::POPULATION_COMPOSITE, "Level populations of composite system")
    .value("FULLSTATE", OutputType::FULLSTATE, "Full quantum state vector");

  // Structs
  nb::class_<InitialConditionSettings>(m, "InitialConditionSettings",
    "Settings for quantum initial state configuration.\n\n"
    "Specify how the initial quantum state is prepared. Different fields are used\n"
    "depending on the condition_type.")
    .def(nb::init<>())
    .def_rw("condition_type", &InitialConditionSettings::type, "(InitialConditionType) Type of initial condition")
    .def_rw("filename", &InitialConditionSettings::filename, "(str | None) File path (for FROMFILE type)")
    .def_rw("levels", &InitialConditionSettings::levels, "(list[int] | None) Level occupations per oscillator (for PRODUCT_STATE)")
    .def_rw("subsystem", &InitialConditionSettings::subsystem, "(list[int] | None) Oscillator IDs (for ENSEMBLE, DIAGONAL, BASIS)")
    .def("__repr__", [](const InitialConditionSettings& s) { return Config::toString(s); });

  nb::class_<OptimTargetSettings>(m, "OptimTargetSettings",
    "Settings for optimization target (gate or state).\n\n"
    "Configure what the optimization should achieve. For GATE targets, specify\n"
    "gate_type or filename. For STATE targets, specify levels or filename.")
    .def(nb::init<>())
    .def_rw("target_type", &OptimTargetSettings::type, "(TargetType) Type of target (GATE or STATE)")
    .def_rw("gate_type", &OptimTargetSettings::gate_type, "(GateType | None) Predefined gate type (for GATE target)")
    .def_rw("gate_rot_freq", &OptimTargetSettings::gate_rot_freq, "(list[float] | None) Gate rotation frequencies [GHz] (for GATE target)")
    .def_rw("levels", &OptimTargetSettings::levels, "(list[int] | None) Target level occupations (for STATE target)")
    .def_rw("filename", &OptimTargetSettings::filename, "(str | None) File path to custom gate or state")
    .def("__repr__", [](const OptimTargetSettings& s) { return Config::toString(s); });

  nb::class_<ControlParameterizationSettings>(m, "ControlParameterizationSettings",
    "Settings for control pulse parameterization.\n\n"
    "Configure how control pulses are represented and parameterized during optimization.")
    .def(nb::init<>())
    .def_rw("control_type", &ControlParameterizationSettings::type, "(ControlType) Parameterization type (BSPLINE, etc.)")
    .def_rw("nspline", &ControlParameterizationSettings::nspline, "(int >= 1 | None) Number of B-spline basis functions")
    .def_rw("tstart", &ControlParameterizationSettings::tstart, "(float | None) Start time of parameterization [ns]")
    .def_rw("tstop", &ControlParameterizationSettings::tstop, "(float | None) Stop time of parameterization [ns]")
    .def_rw("scaling", &ControlParameterizationSettings::scaling, "(float | None) Amplitude scaling factor (BSPLINEAMP only)")
    .def("__repr__", [](const ControlParameterizationSettings& s) { return Config::toString(s); });

  nb::class_<ControlInitializationSettings>(m, "ControlInitializationSettings",
    "Settings for initial control pulse values.\n\n"
    "Configure how control parameters are initialized before optimization.")
    .def(nb::init<>())
    .def_rw("init_type", &ControlInitializationSettings::type, "(ControlInitializationType) Initialization method")
    .def_rw("amplitude", &ControlInitializationSettings::amplitude, "(float | None) Initial amplitude [GHz] (for CONSTANT)")
    .def_rw("phase", &ControlInitializationSettings::phase, "(float | None) Initial phase [rad]")
    .def_rw("filename", &ControlInitializationSettings::filename, "(str | None) File path (for FILE type)")
    .def("__repr__", [](const ControlInitializationSettings& s) { return Config::toString(s); });

  // Setup - mutable configuration with all fields optional
  nb::class_<Setup>(m, "Setup",
    "Mutable configuration builder for Quandary simulations.\n\n"
    "Use this class to build a configuration by setting individual fields.\n"
    "All fields are optional and will use defaults if not specified.\n"
    "Pass to Config() or run() to validate and execute.\n\n"
    "Example:\n"
    "    config = Setup()\n"
    "    config.nlevels = [2, 2]\n"
    "    config.ntime = 100\n"
    "    config.dt = 0.01\n"
    "    results = run(config)\n\n"
    "Note: Use .copy() to create independent copies of Setup objects.")
    .def(nb::init<>())
    .def(nb::init<const Setup&>(),
      nb::arg("other"),
      "Copy constructor - creates a copy of another Setup")
    .def("copy", [](const Setup& self) {
        return Setup(self);
      }, "Create a copy of this Setup")
    // System parameters
    .def_rw("nlevels", &Setup::nlevels, "(list[int]) Number of levels per subsystem")
    .def_rw("nessential", &Setup::nessential, "(list[int]) Number of essential levels per subsystem (default: same as nlevels)")
    .def_rw("ntime", &Setup::ntime, "(int) Number of time steps")
    .def_rw("dt", &Setup::dt, "(float) Time step size [ns]")
    .def_rw("transition_frequency", &Setup::transition_frequency, "(list[float]) Fundamental transition frequencies [GHz]")
    .def_rw("selfkerr", &Setup::selfkerr, "(list[float]) Self-Kerr frequencies [GHz]")
    .def_rw("crosskerr_coupling", &Setup::crosskerr_coupling, "(list[float]) Cross-Kerr coupling frequencies [GHz]")
    .def_rw("dipole_coupling", &Setup::dipole_coupling, "(list[float]) Dipole-dipole coupling frequencies [GHz]")
    .def_rw("rotation_frequency", &Setup::rotation_frequency, "(list[float]) Rotating frame frequencies [GHz]")
    .def_rw("decoherence_type", &Setup::decoherence_type, "(DecoherenceType) Decoherence model type")
    .def_rw("decay_time", &Setup::decay_time, "(list[float]) T1 decay times [ns]")
    .def_rw("dephase_time", &Setup::dephase_time, "(list[float]) T2 dephasing times [ns]")
    .def_rw("initial_condition", &Setup::initial_condition, "(InitialConditionSettings) Initial quantum state configuration")
    // Inherently optional
    .def_rw("hamiltonian_file_Hsys", &Setup::hamiltonian_file_Hsys, "(str | None) Optional file path for system Hamiltonian")
    .def_rw("hamiltonian_file_Hc", &Setup::hamiltonian_file_Hc, "(str | None) Optional file path for control Hamiltonians")
    // Control parameters
    .def_rw("control_zero_boundary_condition", &Setup::control_zero_boundary_condition, "(bool) Enforce zero control amplitude at boundaries")
    .def_rw("control_parameterizations", &Setup::control_parameterizations, "(list[ControlParameterizationSettings]) Control parameterizations per oscillator")
    .def_rw("control_initializations", &Setup::control_initializations, "(list[ControlInitializationSettings]) Control initializations per oscillator")
    .def_rw("control_amplitude_bounds", &Setup::control_amplitude_bounds, "(list[float]) Maximum control amplitude per oscillator [GHz]")
    .def_rw("carrier_frequencies", &Setup::carrier_frequencies, "(list[list[float]]) Carrier frequencies for each control [GHz]")
    // Optimization parameters
    .def_rw("optim_target", &Setup::optim_target, "(OptimTargetSettings) Optimization target configuration")
    .def_rw("optim_objective", &Setup::optim_objective, "(ObjectiveType) Objective function type")
    .def_rw("optim_weights", &Setup::optim_weights, "(list[float]) Weights for summing objective function")
    .def_rw("optim_tol_grad_abs", &Setup::optim_tol_grad_abs, "(float) Absolute gradient tolerance")
    .def_rw("optim_tol_grad_rel", &Setup::optim_tol_grad_rel, "(float) Relative gradient tolerance")
    .def_rw("optim_tol_final_cost", &Setup::optim_tol_final_cost, "(float) Final cost tolerance")
    .def_rw("optim_tol_infidelity", &Setup::optim_tol_infidelity, "(float) Infidelity tolerance")
    .def_rw("optim_maxiter", &Setup::optim_maxiter, "(int) Maximum optimization iterations")
    .def_rw("optim_tikhonov_coeff", &Setup::optim_tikhonov_coeff, "(float) Tikhonov regularization coefficient")
    .def_rw("optim_tikhonov_use_x0", &Setup::optim_tikhonov_use_x0, "(bool) Use initial guess in Tikhonov regularization")
    .def_rw("optim_penalty_leakage", &Setup::optim_penalty_leakage, "(float) Leakage penalty coefficient")
    .def_rw("optim_penalty_weightedcost", &Setup::optim_penalty_weightedcost, "(float) Weighted cost penalty coefficient")
    .def_rw("optim_penalty_weightedcost_width", &Setup::optim_penalty_weightedcost_width, "(float) Weighted cost penalty width parameter")
    .def_rw("optim_penalty_dpdm", &Setup::optim_penalty_dpdm, "(float) Second derivative penalty coefficient")
    .def_rw("optim_penalty_energy", &Setup::optim_penalty_energy, "(float) Energy penalty coefficient")
    .def_rw("optim_penalty_variation", &Setup::optim_penalty_variation, "(float) Amplitude variation penalty coefficient")
    // Output parameters
    .def_rw("output_directory", &Setup::output_directory, "(str) Directory for output files")
    .def_rw("output_observables", &Setup::output_observables, "(list[OutputType]) Observable quantities to output")
    .def_rw("output_timestep_stride", &Setup::output_timestep_stride, "(int) Output frequency in time steps")
    .def_rw("output_optimization_stride", &Setup::output_optimization_stride, "(int) Output frequency in optimization iterations")
    // Solver parameters
    .def_rw("runtype", &Setup::runtype, "(RunType) Type of computation to perform")
    .def_rw("usematfree", &Setup::usematfree, "(bool) Use matrix-free solver")
    .def_rw("linearsolver_type", &Setup::linearsolver_type, "(LinearSolverType) Linear solver type")
    .def_rw("linearsolver_maxiter", &Setup::linearsolver_maxiter, "(int) Maximum linear solver iterations")
    .def_rw("timestepper_type", &Setup::timestepper_type, "(TimeStepperType) Time integration method")
    .def_rw("rand_seed", &Setup::rand_seed, "(int) Random seed for reproducibility");

  // Config - validated configuration with computed values
  nb::class_<Config>(m, "Config",
    "Validated and immutable Quandary configuration.\n\n"
    "This class contains validated configuration parameters with all defaults applied.\n"
    "Instances are immutable after construction. All properties are read-only.\n\n"
    "Create from Setup, TOML file, or TOML string:\n"
    "    config = Config(setup)\n"
    "    config = Config.from_toml('simulation.toml')\n"
    "    config = Config.from_toml_string(toml_content)")
    .def(nb::init<const Setup&, bool>(),
      nb::arg("input"), nb::arg("quiet") = false,
      "Create a validated Config from Setup")
    .def(nb::init<const Config&>(),
      nb::arg("other"),
      "Copy constructor - creates a copy of another Config")
    .def_static("from_file", &Config::fromFile,
      nb::arg("filename"), nb::arg("quiet") = false,
      "Load and validate a Config from a file (auto-detects TOML or .cfg format)")
    .def_static("from_toml", &Config::fromToml,
      nb::arg("toml_filename"), nb::arg("quiet") = false,
      "Load and validate a Config from a TOML file")
    .def_static("from_toml_string", &Config::fromTomlString,
      nb::arg("toml_content"), nb::arg("quiet") = false,
      "Load and validate a Config from a TOML string")
    // System configuration
    .def_prop_ro("nlevels", [](const Config& c) { return c.getNLevels(); },
      "Number of total levels (essential + guard) per oscillator")
    .def_prop_ro("nessential", [](const Config& c) { return c.getNEssential(); },
      "Number of essential levels per oscillator")
    .def_prop_ro("num_osc", &Config::getNumOsc,
      "Number of oscillators")
    // Time configuration
    .def_prop_ro("ntime", &Config::getNTime,
      "Number of time steps")
    .def_prop_ro("dt", &Config::getDt,
      "Time step size [ns]")
    .def_prop_ro("total_time", &Config::getTotalTime,
      "Total simulation time [ns]")
    // Physical parameters
    .def_prop_ro("transition_frequency", [](const Config& c) { return c.getTransitionFrequency(); },
      "Transition frequencies [GHz]")
    .def_prop_ro("selfkerr", [](const Config& c) { return c.getSelfKerr(); },
      "Self-Kerr coefficients [GHz]")
    .def_prop_ro("crosskerr_coupling", [](const Config& c) { return c.getCrossKerrCoupling(); },
      "Cross-Kerr coefficients [GHz]")
    .def_prop_ro("dipole_coupling", [](const Config& c) { return c.getDipoleCoupling(); },
      "Dipole-dipole coupling strengths [GHz]")
    .def_prop_ro("rotation_frequency", [](const Config& c) { return c.getRotationFrequency(); },
      "Rotating frame frequencies [GHz]")
    // Decoherence
    .def_prop_ro("decoherence_type", &Config::getDecoherenceType,
      "Decoherence type (NONE, DECAY, DEPHASE, or BOTH)")
    .def_prop_ro("decay_time", [](const Config& c) { return c.getDecayTime(); },
      "T1 decay times [ns]")
    .def_prop_ro("dephase_time", [](const Config& c) { return c.getDephaseTime(); },
      "T2 dephasing times [ns]")
    // Initial conditions
    .def_prop_ro("n_initial_conditions", &Config::getNInitialConditions,
      "Number of initial conditions")
    .def_prop_ro("initial_condition", [](const Config& c) { return c.getInitialCondition(); },
      "Initial condition settings")
    // Output
    .def_prop_ro("output_directory", &Config::getOutputDirectory,
      "Output directory path")
    // Methods
    .def("to_toml", [](const Config& c) {
        std::stringstream ss;
        c.printConfig(ss);
        return ss.str();
      }, "Serialize the configuration to TOML format")
    .def("__str__", [](const Config& c) {
        std::stringstream ss;
        c.printConfig(ss);
        return ss.str();
      });

  // Run function - accepts Config
  m.def("run", [](const Config& config, bool quiet) {
      return runQuandary(config, quiet);
    },
    nb::arg("config"), nb::arg("quiet") = false,
    "Run a Quandary simulation or optimization from a Config");

  // Run function - accepts Setup, creates Config internally
  m.def("run", [](const Setup& input, bool quiet) {
      Config config(input, quiet);
      return runQuandary(config, quiet);
    },
    nb::arg("config"), nb::arg("quiet") = false,
    "Run a Quandary simulation or optimization from a Setup");

  // Run from file - loads TOML and runs directly
  m.def("run_from_file", [](const std::string& filename, bool quiet) {
      Config config = Config::fromFile(filename, quiet);
      return runQuandary(config, quiet);
    },
    nb::arg("filename"), nb::arg("quiet") = false,
    "Run a Quandary simulation or optimization from a TOML file");
}
