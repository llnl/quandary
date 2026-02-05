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
  m.doc() = "Quandary Python bindings";

  // Enums
  nb::enum_<DecoherenceType>(m, "DecoherenceType")
    .value("NONE", DecoherenceType::NONE)
    .value("DECAY", DecoherenceType::DECAY)
    .value("DEPHASE", DecoherenceType::DEPHASE)
    .value("BOTH", DecoherenceType::BOTH);

  nb::enum_<InitialConditionType>(m, "InitialConditionType")
    .value("FROMFILE", InitialConditionType::FROMFILE)
    .value("PRODUCT_STATE", InitialConditionType::PRODUCT_STATE)
    .value("ENSEMBLE", InitialConditionType::ENSEMBLE)
    .value("DIAGONAL", InitialConditionType::DIAGONAL)
    .value("BASIS", InitialConditionType::BASIS)
    .value("THREESTATES", InitialConditionType::THREESTATES)
    .value("NPLUSONE", InitialConditionType::NPLUSONE)
    .value("PERFORMANCE", InitialConditionType::PERFORMANCE);

  nb::enum_<TargetType>(m, "TargetType")
    .value("NONE", TargetType::NONE)
    .value("GATE", TargetType::GATE)
    .value("STATE", TargetType::STATE);

  nb::enum_<ObjectiveType>(m, "ObjectiveType")
    .value("JFROBENIUS", ObjectiveType::JFROBENIUS)
    .value("JTRACE", ObjectiveType::JTRACE)
    .value("JMEASURE", ObjectiveType::JMEASURE);

  nb::enum_<LinearSolverType>(m, "LinearSolverType")
    .value("GMRES", LinearSolverType::GMRES)
    .value("NEUMANN", LinearSolverType::NEUMANN);

  nb::enum_<RunType>(m, "RunType")
    .value("SIMULATION", RunType::SIMULATION)
    .value("GRADIENT", RunType::GRADIENT)
    .value("OPTIMIZATION", RunType::OPTIMIZATION)
    .value("EVALCONTROLS", RunType::EVALCONTROLS)
    .value("NONE", RunType::NONE);

  nb::enum_<ControlType>(m, "ControlType")
    .value("NONE", ControlType::NONE)
    .value("BSPLINE", ControlType::BSPLINE)
    .value("BSPLINEAMP", ControlType::BSPLINEAMP)
    .value("BSPLINE0", ControlType::BSPLINE0);

  nb::enum_<ControlInitializationType>(m, "ControlInitializationType")
    .value("CONSTANT", ControlInitializationType::CONSTANT)
    .value("RANDOM", ControlInitializationType::RANDOM)
    .value("FILE", ControlInitializationType::FILE);

  nb::enum_<TimeStepperType>(m, "TimeStepperType")
    .value("IMR", TimeStepperType::IMR)
    .value("IMR4", TimeStepperType::IMR4)
    .value("IMR8", TimeStepperType::IMR8)
    .value("EE", TimeStepperType::EE);

  nb::enum_<GateType>(m, "GateType")
    .value("NONE", GateType::NONE)
    .value("XGATE", GateType::XGATE)
    .value("YGATE", GateType::YGATE)
    .value("ZGATE", GateType::ZGATE)
    .value("HADAMARD", GateType::HADAMARD)
    .value("CNOT", GateType::CNOT)
    .value("SWAP", GateType::SWAP)
    .value("SWAP_0Q", GateType::SWAP_0Q)
    .value("CQNOT", GateType::CQNOT)
    .value("QFT", GateType::QFT)
    .value("FILE", GateType::FILE);

  nb::enum_<OutputType>(m, "OutputType")
    .value("EXPECTED_ENERGY", OutputType::EXPECTED_ENERGY)
    .value("EXPECTED_ENERGY_COMPOSITE", OutputType::EXPECTED_ENERGY_COMPOSITE)
    .value("POPULATION", OutputType::POPULATION)
    .value("POPULATION_COMPOSITE", OutputType::POPULATION_COMPOSITE)
    .value("FULLSTATE", OutputType::FULLSTATE);

  // Structs
  nb::class_<InitialConditionSettings>(m, "InitialConditionSettings")
    .def(nb::init<>())
    .def_rw("condition_type", &InitialConditionSettings::type)
    .def_rw("filename", &InitialConditionSettings::filename)
    .def_rw("levels", &InitialConditionSettings::levels)
    .def_rw("subsystem", &InitialConditionSettings::subsystem);

  nb::class_<OptimTargetSettings>(m, "OptimTargetSettings")
    .def(nb::init<>())
    .def_rw("target_type", &OptimTargetSettings::type)
    .def_rw("gate_type", &OptimTargetSettings::gate_type)
    .def_rw("gate_rot_freq", &OptimTargetSettings::gate_rot_freq)
    .def_rw("levels", &OptimTargetSettings::levels)
    .def_rw("filename", &OptimTargetSettings::filename);

  nb::class_<ControlParameterizationSettings>(m, "ControlParameterizationSettings")
    .def(nb::init<>())
    .def_rw("control_type", &ControlParameterizationSettings::type)
    .def_rw("nspline", &ControlParameterizationSettings::nspline)
    .def_rw("tstart", &ControlParameterizationSettings::tstart)
    .def_rw("tstop", &ControlParameterizationSettings::tstop)
    .def_rw("scaling", &ControlParameterizationSettings::scaling);

  nb::class_<ControlInitializationSettings>(m, "ControlInitializationSettings")
    .def(nb::init<>())
    .def_rw("init_type", &ControlInitializationSettings::type)
    .def_rw("amplitude", &ControlInitializationSettings::amplitude)
    .def_rw("phase", &ControlInitializationSettings::phase)
    .def_rw("filename", &ControlInitializationSettings::filename);

  // QuandaryConfig - mutable configuration with all fields optional
  nb::class_<RawConfig>(m, "QuandaryConfig")
    .def(nb::init<>())
    // System parameters
    .def_rw("nlevels", &RawConfig::nlevels)
    .def_rw("nessential", &RawConfig::nessential)
    .def_rw("ntime", &RawConfig::ntime)
    .def_rw("dt", &RawConfig::dt)
    .def_rw("transition_frequency", &RawConfig::transition_frequency)
    .def_rw("selfkerr", &RawConfig::selfkerr)
    .def_rw("crosskerr_coupling", &RawConfig::crosskerr_coupling)
    .def_rw("dipole_coupling", &RawConfig::dipole_coupling)
    .def_rw("rotation_frequency", &RawConfig::rotation_frequency)
    .def_rw("decoherence_type", &RawConfig::decoherence_type)
    .def_rw("decay_time", &RawConfig::decay_time)
    .def_rw("dephase_time", &RawConfig::dephase_time)
    .def_rw("initial_condition", &RawConfig::initial_condition)
    // Inherently optional
    .def_rw("hamiltonian_file_Hsys", &RawConfig::hamiltonian_file_Hsys)
    .def_rw("hamiltonian_file_Hc", &RawConfig::hamiltonian_file_Hc)
    // Control parameters
    .def_rw("control_zero_boundary_condition", &RawConfig::control_zero_boundary_condition)
    .def_rw("control_parameterizations", &RawConfig::control_parameterizations)
    .def_rw("control_initializations", &RawConfig::control_initializations)
    .def_rw("control_amplitude_bounds", &RawConfig::control_amplitude_bounds)
    .def_rw("carrier_frequencies", &RawConfig::carrier_frequencies)
    // Optimization parameters
    .def_rw("optim_target", &RawConfig::optim_target)
    .def_rw("optim_objective", &RawConfig::optim_objective)
    .def_rw("optim_weights", &RawConfig::optim_weights)
    .def_rw("optim_tol_grad_abs", &RawConfig::optim_tol_grad_abs)
    .def_rw("optim_tol_grad_rel", &RawConfig::optim_tol_grad_rel)
    .def_rw("optim_tol_final_cost", &RawConfig::optim_tol_final_cost)
    .def_rw("optim_tol_infidelity", &RawConfig::optim_tol_infidelity)
    .def_rw("optim_maxiter", &RawConfig::optim_maxiter)
    .def_rw("optim_tikhonov_coeff", &RawConfig::optim_tikhonov_coeff)
    .def_rw("optim_tikhonov_use_x0", &RawConfig::optim_tikhonov_use_x0)
    .def_rw("optim_penalty_leakage", &RawConfig::optim_penalty_leakage)
    .def_rw("optim_penalty_weightedcost", &RawConfig::optim_penalty_weightedcost)
    .def_rw("optim_penalty_weightedcost_width", &RawConfig::optim_penalty_weightedcost_width)
    .def_rw("optim_penalty_dpdm", &RawConfig::optim_penalty_dpdm)
    .def_rw("optim_penalty_energy", &RawConfig::optim_penalty_energy)
    .def_rw("optim_penalty_variation", &RawConfig::optim_penalty_variation)
    // Output parameters
    .def_rw("output_directory", &RawConfig::output_directory)
    .def_rw("output_observables", &RawConfig::output_observables)
    .def_rw("output_timestep_stride", &RawConfig::output_timestep_stride)
    .def_rw("output_optimization_stride", &RawConfig::output_optimization_stride)
    // Solver parameters
    .def_rw("runtype", &RawConfig::runtype)
    .def_rw("usematfree", &RawConfig::usematfree)
    .def_rw("linearsolver_type", &RawConfig::linearsolver_type)
    .def_rw("linearsolver_maxiter", &RawConfig::linearsolver_maxiter)
    .def_rw("timestepper_type", &RawConfig::timestepper_type)
    .def_rw("rand_seed", &RawConfig::rand_seed);

  // Config - validated configuration with computed values
  nb::class_<Config>(m, "Config")
    .def(nb::init<const RawConfig&, bool>(),
      nb::arg("input"), nb::arg("quiet") = false,
      "Create a validated Config from RawConfig")
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
    .def("setup_for_eval_controls", &Config::setupForEvalControls,
      nb::arg("points_per_ns"), nb::arg("pcof_file"), nb::arg("output_dir"),
      "Configure for control evaluation at specified sample rate");

  // Run function - accepts Config
  m.def("run", [](const Config& config, bool quiet) {
      return runQuandary(config, quiet);
    },
    nb::arg("config"), nb::arg("quiet") = false,
    "Run a Quandary simulation or optimization from a Config");

  // Run function - accepts QuandaryConfig, creates Config internally
  m.def("run", [](const RawConfig& input, bool quiet) {
      Config config(input, quiet);
      return runQuandary(config, quiet);
    },
    nb::arg("config"), nb::arg("quiet") = false,
    "Run a Quandary simulation or optimization from a QuandaryConfig");

  // Run from file - loads TOML and runs directly
  m.def("run_from_file", [](const std::string& filename, bool quiet) {
      Config config = Config::fromFile(filename, quiet);
      return runQuandary(config, quiet);
    },
    nb::arg("filename"), nb::arg("quiet") = false,
    "Run a Quandary simulation or optimization from a TOML file");
}
