#include <nanobind/nanobind.h>

namespace nb = nanobind;

NB_MODULE(quandary_ext, m) {
  m.doc() = "Quandary Python bindings";
}
