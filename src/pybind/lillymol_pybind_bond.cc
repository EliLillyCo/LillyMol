#include <iostream>
#include <sstream>
#include <string>

#include "pybind11/pybind11.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"

#include "Foundational/iwstring/iwstring.h"
#include "Molecule_Lib/bond.h"

namespace py = pybind11;

PYBIND11_MODULE(lillymol_bond, b)
{
  py::class_<Bond, std::shared_ptr<Bond>>(b, "Bond")
    .def(py::init<>())
    .def("a1", &Bond::a1, "First atom")
    .def("a2", &Bond::a2, "Second atom")
    .def("other", static_cast<atom_number_t (Bond::*)(int)const>(&Bond::other), "Other connection")
    .def("involves",
      [](const Bond& b, atom_number_t a)->bool{
        return b.involves(a);
      },
      "True if bond involves atom"
    )
    .def("is_single_bond",
      [](const Bond& b)->bool{
        return b.is_single_bond();
      },
      "True if a single bond"
    )
    .def("is_double_bond",
      [](const Bond& b)->bool{
        return b.is_double_bond();
      },
      "True if a double bond"
    )
    .def("is_triple_bond",
      [](const Bond& b)->bool{
        return b.is_triple_bond();
      },
      "True if a triple bond"
    )
    .def("nrings",
      [](const Bond& b){
        int nr;
        if (b.nrings(nr)) {
          return nr;
        }
        throw py::value_error("Bond.nrings:ring membership not computed");
      },
      "Number of rings involving bond"
    )
    .def("__repr__",
      [](const Bond& b) {
        IWString result;
        result << "<Bond " << b.a1();
        if (b.is_single_bond()) {
          result << '-';
        } else if (b.is_double_bond()) {
          result << '=';
        } else if (b.is_triple_bond()) {
          result << '=';
        }
        result << b.a2() << '>';
        return std::string(result.data(), result.size());
      }
    )
  ;
}
