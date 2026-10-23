#include "pybind11/pybind11.h"

namespace py = pybind11;

PYBIND11_MODULE(lillymol_set_of_atoms, m) {
  py::module_ lillymol = py::module_::import("lillymol");
  m.attr("Set_of_Atoms") = lillymol.attr("Set_of_Atoms");
}
