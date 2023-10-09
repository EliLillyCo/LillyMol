#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "Foundational/iwstring/iwstring.h"
#include "Molecule_Lib/chiral_centre.h"

namespace py = pybind11;

static void
AppendChiralComponent(int a,
                      IWString& result) {
  if (a >= 0) {
    result << a;
  } else if (a == CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN) {
    result << 'H';
  } else if (a == CHIRAL_CONNECTION_IS_LONE_PAIR) {
    result << '^';
  } else {
    result << '?';
  }
}

PYBIND11_MODULE(lillymol_chiral_centre, c)
{
  py::class_<Chiral_Centre>(c, "Chiral_Centre")
    .def(py::init<atom_number_t>())
    .def("atom", &Chiral_Centre::a, "Centre atom")
    .def("top_front", &Chiral_Centre::top_front, "Top front")
    .def("top_back", &Chiral_Centre::top_back, "Top back")
    .def("left_down", &Chiral_Centre::left_down, "Left down")
    .def("right_down", &Chiral_Centre::right_down, "Right down")
    .def("invert", &Chiral_Centre::invert, "Invert")
    .def("involves",
      [](const Chiral_Centre& c, atom_number_t zatom)->bool{
        return c.involves(zatom);
      },
      "True if atom is part of chiral centre"
    ) 
    .def("implicit_hydrogens",
      [](const Chiral_Centre&c) {
        return c.implicit_hydrogen_count();
      },
      "Number of implicit hydrogens - can be only 1"
    )
    .def("lone_pairs",
      [](const Chiral_Centre&c) {
        return c.lone_pair_count();
      },
      "Number of lone pairs - can be only 1"
    )
    .def("implicit_hydrogen_is_now_atom_number", &Chiral_Centre::implicit_hydrogen_is_now_atom_number)
    .def("lone_pair_is_now_atom_number", &Chiral_Centre::lone_pair_is_now_atom_number)
    .def("atom_is_now_implicit_hydrogen", &Chiral_Centre::atom_is_now_implicit_hydrogen)
    .def("atom_is_now_lone_pair", &Chiral_Centre::atom_is_now_lone_pair)
    .def("atom_numbers_are_swapped", &Chiral_Centre::atom_numbers_are_swapped)

    .def("__repr__",
      [](const Chiral_Centre& c) {
        IWString result;
        result << "<Chiral_Centre atom " << c.a();
        result << " tf ";
        AppendChiralComponent(c.top_front(), result);
        result << " tb ";
        AppendChiralComponent(c.top_back(), result);
        result << " ld ";
        AppendChiralComponent(c.left_down(), result);
        result << " rd ";
        AppendChiralComponent(c.right_down(), result);
        result << '>';
        return std::string(result.data(), result.size());
      }
    )
  ;

  py::enum_<CahnIngoldPrelog>(c, "CIP")
    .value("R", CahnIngoldPrelog::R)
    .value("S", CahnIngoldPrelog::S)
    .value("Neither", CahnIngoldPrelog::kNeither)
    .value("Unspecified", CahnIngoldPrelog::kUnspecified)
    .export_values();
  ;
}
