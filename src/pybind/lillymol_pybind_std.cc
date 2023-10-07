#include "pybind11/pybind11.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

namespace py = pybind11;

PYBIND11_MODULE(lillymol_standardise, s)
{
  py::class_<Chemical_Standardisation>(s, "Standardise")
    .def(py::init<>())
    .def("activate_all", &Chemical_Standardisation::activate_all, "Activate all transformations")
    .def("set_verbose", &Chemical_Standardisation::set_verbose, "Set verbosity")
    .def("process", &Chemical_Standardisation::process, "Apply active transformations to molecule")
  ;
}
