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
    .def_rw("type", &InitialConditionSettings::type)
    .def_rw("filename", &InitialConditionSettings::filename)
    .def_rw("levels", &InitialConditionSettings::levels)
    .def_rw("subsystem", &InitialConditionSettings::subsystem);

  nb::class_<OptimTargetSettings>(m, "OptimTargetSettings")
    .def(nb::init<>())
    .def_rw("type", &OptimTargetSettings::type)
    .def_rw("gate_type", &OptimTargetSettings::gate_type)
    .def_rw("gate_rot_freq", &OptimTargetSettings::gate_rot_freq)
    .def_rw("levels", &OptimTargetSettings::levels)
    .def_rw("filename", &OptimTargetSettings::filename);

  nb::class_<ControlParameterizationSettings>(m, "ControlParameterizationSettings")
    .def(nb::init<>())
    .def_rw("type", &ControlParameterizationSettings::type)
    .def_rw("nspline", &ControlParameterizationSettings::nspline)
    .def_rw("tstart", &ControlParameterizationSettings::tstart)
    .def_rw("tstop", &ControlParameterizationSettings::tstop)
    .def_rw("scaling", &ControlParameterizationSettings::scaling);

  nb::class_<ControlInitializationSettings>(m, "ControlInitializationSettings")
    .def(nb::init<>())
    .def_rw("type", &ControlInitializationSettings::type)
    .def_rw("amplitude", &ControlInitializationSettings::amplitude)
    .def_rw("phase", &ControlInitializationSettings::phase)
    .def_rw("filename", &ControlInitializationSettings::filename);

  // QuandaryConfig - mutable configuration with all fields optional
  nb::class_<ConfigInput>(m, "QuandaryConfig")
    .def(nb::init<>())
    // System parameters
    .def_rw("nlevels", &ConfigInput::nlevels)
    .def_rw("nessential", &ConfigInput::nessential)
    .def_rw("ntime", &ConfigInput::ntime)
    .def_rw("dt", &ConfigInput::dt)
    .def_rw("transfreq", &ConfigInput::transfreq)
    .def_rw("selfkerr", &ConfigInput::selfkerr)
    .def_rw("crosskerr", &ConfigInput::crosskerr)
    .def_rw("Jkl", &ConfigInput::Jkl)
    .def_rw("rotfreq", &ConfigInput::rotfreq)
    .def_rw("decoherence_type", &ConfigInput::decoherence_type)
    .def_rw("decay_time", &ConfigInput::decay_time)
    .def_rw("dephase_time", &ConfigInput::dephase_time)
    .def_rw("initial_condition", &ConfigInput::initial_condition)
    // Inherently optional
    .def_rw("hamiltonian_file_Hsys", &ConfigInput::hamiltonian_file_Hsys)
    .def_rw("hamiltonian_file_Hc", &ConfigInput::hamiltonian_file_Hc)
    // Control parameters
    .def_rw("control_zero_boundary_condition", &ConfigInput::control_zero_boundary_condition)
    .def_rw("control_parameterizations", &ConfigInput::control_parameterizations)
    .def_rw("control_initializations", &ConfigInput::control_initializations)
    .def_rw("control_amplitude_bounds", &ConfigInput::control_amplitude_bounds)
    .def_rw("carrier_frequencies", &ConfigInput::carrier_frequencies)
    // Optimization parameters
    .def_rw("optim_target", &ConfigInput::optim_target)
    .def_rw("optim_objective", &ConfigInput::optim_objective)
    .def_rw("optim_weights", &ConfigInput::optim_weights)
    .def_rw("optim_tol_grad_abs", &ConfigInput::optim_tol_grad_abs)
    .def_rw("optim_tol_grad_rel", &ConfigInput::optim_tol_grad_rel)
    .def_rw("optim_tol_finalcost", &ConfigInput::optim_tol_finalcost)
    .def_rw("optim_tol_infidelity", &ConfigInput::optim_tol_infidelity)
    .def_rw("optim_maxiter", &ConfigInput::optim_maxiter)
    .def_rw("optim_tikhonov_coeff", &ConfigInput::optim_tikhonov_coeff)
    .def_rw("optim_tikhonov_use_x0", &ConfigInput::optim_tikhonov_use_x0)
    .def_rw("optim_penalty_leakage", &ConfigInput::optim_penalty_leakage)
    .def_rw("optim_penalty_weightedcost", &ConfigInput::optim_penalty_weightedcost)
    .def_rw("optim_penalty_weightedcost_width", &ConfigInput::optim_penalty_weightedcost_width)
    .def_rw("optim_penalty_dpdm", &ConfigInput::optim_penalty_dpdm)
    .def_rw("optim_penalty_energy", &ConfigInput::optim_penalty_energy)
    .def_rw("optim_penalty_variation", &ConfigInput::optim_penalty_variation)
    // Output parameters
    .def_rw("output_directory", &ConfigInput::output_directory)
    .def_rw("output_observables", &ConfigInput::output_observables)
    .def_rw("output_timestep_stride", &ConfigInput::output_timestep_stride)
    .def_rw("output_optimization_stride", &ConfigInput::output_optimization_stride)
    // Solver parameters
    .def_rw("runtype", &ConfigInput::runtype)
    .def_rw("usematfree", &ConfigInput::usematfree)
    .def_rw("linearsolver_type", &ConfigInput::linearsolver_type)
    .def_rw("linearsolver_maxiter", &ConfigInput::linearsolver_maxiter)
    .def_rw("timestepper_type", &ConfigInput::timestepper_type)
    .def_rw("rand_seed", &ConfigInput::rand_seed);

  // Config - validated configuration with computed values
  nb::class_<Config>(m, "Config")
    .def(nb::init<const ConfigInput&, bool>(),
      nb::arg("input"), nb::arg("quiet") = false,
      "Create a validated Config from ConfigInput")
    .def_prop_ro("n_initial_conditions", &Config::getNInitialConditions,
      "Number of initial conditions")
    .def_prop_ro("output_directory", &Config::getOutputDirectory,
      "Output directory path")
    .def_prop_ro("decoherence_type", &Config::getDecoherenceType,
      "Decoherence type (NONE, DECAY, DEPHASE, or BOTH)")
    .def("to_toml", [](const Config& c) {
        std::stringstream ss;
        c.printConfig(ss);
        return ss.str();
      }, "Serialize the configuration to TOML format");

  // Run function - accepts Config
  m.def("run", [](const Config& config, bool quiet) {
      return runQuandary(config, quiet);
    },
    nb::arg("config"), nb::arg("quiet") = false,
    "Run a Quandary simulation or optimization from a Config");

  // Run function - accepts QuandaryConfig, creates Config internally
  m.def("run", [](const ConfigInput& input, bool quiet) {
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
