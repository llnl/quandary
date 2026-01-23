#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include "config.hpp"
#include "quandary_core.hpp"

namespace nb = nanobind;

// Declare extension module "quandary_ext" (must match name given in CMake)
NB_MODULE(quandary_ext, m) {
  m.doc() = "Quandary Python bindings";

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
