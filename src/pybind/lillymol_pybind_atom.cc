#include <iostream>
#include <sstream>
#include <string>

#include "pybind11/pybind11.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"

#include "Molecule_Lib/atom.h"

namespace py = pybind11;

static void
AppendBond(const Bond* b,
           IWString& result) {
  result << b->a1();
  if (b->is_single_bond()) {
    result << '-';
  } else if (b->is_double_bond()) {
    result << '=';
  } else if (b->is_triple_bond()) {
    result << '#';
  }
  result << b->a2();
}

static void
AtomString(Atom& a,
           IWString& result) {
  result << a.atomic_symbol();
  if (a.isotope()) {
    result << " iso " << a.isotope();
  }
  int ih;
  a.compute_implicit_hydrogens(ih);
  result << " ih " << ih;
  result << " ncon " << a.ncon();
  if (a.ncon() == 0) {
    return;
  }

  if (a.ncon() == 1) {
    result << " bond ";
    AppendBond(a[0], result);
    return;
  }

  extending_resizable_array<int> seen;
  for (const Bond* b : a) {
    ++seen[b->a1()];
    ++seen[b->a2()];
  }

  result << " bonded";
  for (int i = 0; i < seen.number_elements(); ++i) {
    if (seen[i] != 1) {
      continue;
    }
    result << ' ' << i;
  }
}

PYBIND11_MODULE(lillymol_atom, a)
{
  py::class_<Atom, std::shared_ptr<Atom>>(a, "Atom")
    .def(py::init<int>())
    .def("atomic_number", static_cast<int (Atom::*)()const>(&Atom::atomic_number), "Atomic Number")
    .def("isotope", static_cast<isotope_t (Atom::*)()const>(&Atom::isotope), "isotope")
    .def("ncon", static_cast<int (Atom::*)()const>(&Atom::ncon), "Number of connections")
    .def("other", static_cast<atom_number_t (Atom::*)(atom_number_t, int)const>(&Atom::other), "Other connection")
    .def("is_organic", &Atom::is_organic, "True if the element is organic")
    .def("atomic_symbol",
      [](const Atom& a)->std::string {
        return a.element()->symbol().AsString();
      },
      "atomic symbol"
    )
    .def("__repr__",
      [](Atom &a) {
        IWString s;
        s << "<Atom " << a.atomic_symbol() << " ncon " << a.ncon();
        if (a.isotope()) {
          s << " iso " << a.isotope();
        }
        s << '>';
        return std::string(s.data(), s.length());
      }
    )
    .def("__getitem__",
      [](const Atom& a, int ndx) {
        return a[ndx];
      },
      py::return_value_policy::reference
    )
    .def("__str__",
      [](Atom& a) {
        IWString result;
        AtomString(a, result);
        return std::string(result.data(), result.size());
      }
    )
    .def("__iter__",
      [](const Atom&a) {
        return py::make_iterator(a.begin(), a.end());
      },
      py::keep_alive<0, 1>()
    )
  ;
}
