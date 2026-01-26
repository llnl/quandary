#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "config.hpp"
#include "quandary_core.hpp"

namespace nb = nanobind;

// Declare extension module "quandary_ext" (must match name given in CMake)
NB_MODULE(quandary_ext, m) {
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

  // Config class
  nb::class_<Config>(m, "Config")
    // fromFile static method
    .def_static("from_file", &Config::fromFile,
      // Arguments
      nb::arg("filename"), nb::arg("quiet") = true,
      // Docstring
      "Load configuration from a TOML file");

  // Run function
  m.def("run", [](const Config& config, bool quiet) { return runQuandary(config, quiet); },
    // Arguments
    nb::arg("config"), nb::arg("quiet") = false,
    // Docstring
    "Run a Quandary simulation or optimization");
}
