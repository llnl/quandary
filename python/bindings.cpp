#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <sstream>
// Make Config class members public for nanobind access
#define private public  
#include "config.hpp"
#undef private
#include "main.hpp"
#include "petsc.h"
#ifdef WITH_SLEPC
#include <slepceps.h>
#endif

namespace nb = nanobind;

// Module name "_quandary_impl" must match the CMake target in nanobind_add_module().
// The underscore prefix is a Python convention for internal modules -- users import
// the public API through python/quandary/new/__init__.py instead.
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
    "gate_type or filename. For STATE target, specify filename.")
    .def(nb::init<>())
    .def_rw("target_type", &OptimTargetSettings::type, "(TargetType) Type of target (GATE or STATE)")
    .def_rw("gate_type", &OptimTargetSettings::gate_type, "(GateType | None) Predefined gate type (for GATE target)")
    .def_rw("gate_rot_freq", &OptimTargetSettings::gate_rot_freq, "(list[float] | None) Gate rotation frequencies [GHz] (for GATE target)")
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

  // Config - final validated configuration data 
  nb::class_<Config>(m, "Config",
    "Validated Quandary configuration.\n\n"
    "This class contains configuration parameters and can be edited from Python.\n\n"
    "Create from constructors, TOML file, or TOML string:\n"
    "    config = Config(False)\n"
    "    config = Config.from_file('simulation.toml')\n"
    "    config = Config.from_string(toml_content)")
    .def(nb::init<bool>(),
      nb::arg("quiet") = false,
      "Construct an empty Config object")
    .def(nb::init<const toml::table&, bool>(),
      nb::arg("toml"), nb::arg("quiet") = false,
      "Construct and validate Config from TOML table")
    .def(nb::init<const Config&>(),
      nb::arg("other"),
      "Copy constructor - creates a copy of another Config")
    .def_static("from_file", &Config::fromFile,
      nb::arg("filename"), nb::arg("quiet") = false,
      "Load and validate a Config from a TOML file")
    .def_static("from_string", &Config::fromString,
      nb::arg("toml_content"), nb::arg("quiet") = false,
      "Load and validate a Config from a TOML string")
    // System configuration
    .def_rw("nlevels", &Config::nlevels,
      "Number of levels per subsystem")
    .def_rw("nessential", &Config::nessential,
      "Number of essential levels per subsystem")
    .def_rw("ntime", &Config::ntime,
      "Number of time steps used for time-integration")
    .def_rw("dt", &Config::dt,
      "Time step size [ns]. Determines final time: T=ntime*dt")
    .def_rw("total_time", &Config::total_time,
      "Total evolution time [ns]. Alternative to ntime and dt")
    .def_rw("transition_frequency", &Config::transition_frequency,
      "Fundamental transition frequencies for each oscillator [GHz]")
    .def_rw("selfkerr", &Config::selfkerr,
      "Self-kerr frequencies for each oscillator [GHz]")
    .def_rw("crosskerr_coupling", &Config::crosskerr_coupling,
      "Cross-kerr coupling frequencies [GHz]")
    .def_rw("dipole_coupling", &Config::dipole_coupling,
      "Dipole-dipole coupling frequencies [GHz]")
    .def_rw("rotation_frequency", &Config::rotation_frequency,
      "Rotating wave approximation frequencies [GHz]")
    .def_rw("decoherence_type", &Config::decoherence_type,
      "Switch between Schroedinger and Lindblad solver")
    .def_rw("decay_time", &Config::decay_time,
      "Decay time T1 per oscillator")
    .def_rw("dephase_time", &Config::dephase_time,
      "Dephase time T2 per oscillator")
    .def_rw("initial_condition", &Config::initial_condition,
      "Initial condition configuration")
    .def_rw("hamiltonian_file_Hsys", &Config::hamiltonian_file_Hsys,
      "Optional file to read the system Hamiltonian from")
    .def_rw("hamiltonian_file_Hc", &Config::hamiltonian_file_Hc,
      "Optional file to read the control Hamiltonian from")
    // Optimization options
    .def_rw("control_zero_boundary_condition", &Config::control_zero_boundary_condition,
      "Whether control pulses should start and end at zero")
    .def_rw("control_parameterizations", &Config::control_parameterizations,
      "Control parameterizations for each oscillator")
    .def_rw("control_initializations", &Config::control_initializations,
      "Control initializations for each oscillator")
    .def_rw("control_amplitude_bounds", &Config::control_amplitude_bounds,
      "Control amplitude bounds for each oscillator")
    .def_rw("carrier_frequencies", &Config::carrier_frequencies,
      "Carrier frequencies for each oscillator")
    .def_rw("optim_target", &Config::optim_target,
      "Grouped optimization target configuration")
    .def_rw("optim_objective", &Config::optim_objective,
      "Objective function measure")
    .def_rw("optim_weights", &Config::optim_weights,
      "Weights for summing up the objective function")
    .def_rw("optim_tol_grad_abs", &Config::optim_tol_grad_abs,
      "Absolute gradient tolerance")
    .def_rw("optim_tol_grad_rel", &Config::optim_tol_grad_rel,
      "Relative gradient tolerance")
    .def_rw("optim_tol_final_cost", &Config::optim_tol_final_cost,
      "Final time cost tolerance")
    .def_rw("optim_tol_infidelity", &Config::optim_tol_infidelity,
      "Infidelity tolerance")
    .def_rw("optim_maxiter", &Config::optim_maxiter,
      "Maximum optimization iterations")
    .def_rw("optim_tikhonov_coeff", &Config::optim_tikhonov_coeff,
      "Coefficient of Tikhonov regularization")
    .def_rw("optim_tikhonov_use_x0", &Config::optim_tikhonov_use_x0,
      "Use Tikhonov regularization with ||x - x0||^2")
    .def_rw("optim_penalty_leakage", &Config::optim_penalty_leakage,
      "Leakage penalty coefficient")
    .def_rw("optim_penalty_weightedcost", &Config::optim_penalty_weightedcost,
      "Weighted cost penalty coefficient")
    .def_rw("optim_penalty_weightedcost_width", &Config::optim_penalty_weightedcost_width,
      "Width parameter for weighted cost penalty")
    .def_rw("optim_penalty_dpdm", &Config::optim_penalty_dpdm,
      "Second derivative penalty coefficient")
    .def_rw("optim_penalty_energy", &Config::optim_penalty_energy,
      "Energy penalty coefficient")
    .def_rw("optim_penalty_variation", &Config::optim_penalty_variation,
      "Amplitude variation penalty coefficient")
    // Output and run types
    .def_rw("output_directory", &Config::output_directory,
      "Directory for output files")
    .def_rw("output_observables", &Config::output_observables,
      "Desired observables")
    .def_rw("output_timestep_stride", &Config::output_timestep_stride,
      "Write output every N time-step")
    .def_rw("output_optimization_stride", &Config::output_optimization_stride,
      "Write output every N optimization iterations")
    .def_rw("runtype", &Config::runtype,
      "Run type: simulation, gradient, or optimization")
    .def_rw("usematfree", &Config::usematfree,
      "Use matrix free solver")
    .def_rw("linearsolver_type", &Config::linearsolver_type,
      "Linear solver type")
    .def_rw("linearsolver_maxiter", &Config::linearsolver_maxiter,
      "Maximum linear solver iterations")
    .def_rw("timestepper_type", &Config::timestepper_type,
      "Time-stepping algorithm")
    .def_rw("rand_seed", &Config::rand_seed,
      "Random number generator seed")
    .def_rw("n_initial_conditions", &Config::n_initial_conditions,
      "Number of initial conditions")
    // Methods
    .def("printConfig", [](const Config& c) {
        std::stringstream ss;
        c.printConfig(ss);
        return ss.str();
      }, "Return configuration as TOML text")
    .def("__str__", [](const Config& c) {
        std::stringstream ss;
        c.printConfig(ss);
        return ss.str();
      });

  m.def("inputFromFile", [](const std::string& filename, bool quiet) {
      return Config::fromFile(filename, quiet);
    },
    nb::arg("filename"), nb::arg("quiet") = false,
    "Parse a TOML file and construct a Config without validation");

  // Run function - accepts Config
  m.def("run", [](const Config& config, bool quiet) {
      return runQuandary(config, quiet);
    },
    nb::arg("config"), nb::arg("quiet") = false,
    "Run a Quandary simulation or optimization from a Config");

  // Run from file - loads TOML and runs directly
  m.def("run_from_file", [](const std::string& filename, bool quiet) {
      Config config = Config::fromFile(filename, quiet);
      return runQuandary(config, quiet);
    },
    nb::arg("filename"), nb::arg("quiet") = false,
    "Run a Quandary simulation or optimization from a TOML file");

  // Finalize PETSc if initialized. Registered as a Python atexit handler
  // so PETSc is cleaned up before MPI is finalized.
  m.def("_finalize_petsc", []() {
      PetscBool initialized = PETSC_FALSE;
      PetscInitialized(&initialized);
      if (initialized) {
        PetscBool finalized = PETSC_FALSE;
        PetscFinalized(&finalized);
        if (!finalized) {
          #ifdef WITH_SLEPC
          SlepcFinalize();
          #else
          PetscFinalize();
          #endif
        }
      }
    });
}
