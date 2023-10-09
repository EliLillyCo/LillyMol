#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "Molecule_Lib/set_of_atoms.h"

namespace py = pybind11;

PYBIND11_MODULE(lillymol_set_of_atoms, s)
{
  py::class_<Set_of_Atoms>(s, "Set_of_Atoms")
    .def(py::init<>())
    .def("empty",
      [](const Set_of_Atoms& s)->bool{
        return s.empty();
      },
      "True if empty"
    )
    .def("size", &Set_of_Atoms::size, "size")
    .def("set_vector",
      [](const Set_of_Atoms& s, std::vector<int>& v, int value) {
        for (auto a : s) {
          v[a] = value;
        }
      }
    )
    .def("__getitem__",
      [](const Set_of_Atoms& me, int ndx) {
        return me.item(ndx);
      })
    .def("__iter__",
      [](const Set_of_Atoms&s) {
        return py::make_iterator(s.begin(), s.end());
      },
      py::keep_alive<0, 1>()
    )
  ;
}
