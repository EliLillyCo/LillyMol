#include "pybind11/pybind11.h"

namespace py = pybind11;

PYBIND11_MODULE(lillymol_atom, m) {
  py::module_ lillymol = py::module_::import("lillymol");
  m.attr("Atom") = lillymol.attr("Atom");
}
