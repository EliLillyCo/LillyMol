#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {

void
BindAtomBond(nb::module_& m) {
  nb::class_<Atom>(m, "Atom")
      .def(nb::init<int>(), nb::arg("atomic_number"))
      .def("atomic_number", nb::overload_cast<>(&Atom::atomic_number, nb::const_))
      .def("atomic_symbol", [](const Atom& atom) { return atom.atomic_symbol().AsString(); })
      .def("isotope", &Atom::isotope)
      .def("ncon", nb::overload_cast<>(&Atom::ncon, nb::const_))
      .def("nbonds", nb::overload_cast<>(&Atom::nbonds, nb::const_))
      .def("formal_charge", &Atom::formal_charge)
      .def("atomic_weight", &Atom::atomic_weight)
      .def("exact_mass", &AtomExactMass)
      .def("implicit_hydrogens", &Atom::implicit_hydrogens)
      .def("is_bonded_to",
           [](const Atom& atom, atom_number_t other) {
             return static_cast<bool>(atom.is_bonded_to(other));
           },
           nb::arg("other"))
      .def("valence_ok", &Atom::valence_ok)
      .def("fully_saturated", &Atom::fully_saturated)
      .def("unsaturation", &Atom::unsaturation)
      .def("other", nb::overload_cast<atom_number_t, int>(&Atom::other, nb::const_), nb::arg("atom"), nb::arg("connection"))
      .def("is_organic", &Atom::is_organic)
      .def("x", [](const Atom& atom) { return atom.x(); }, "Return x coordinate")
      .def("y", [](const Atom& atom) { return atom.y(); }, "Return y coordinate")
      .def("z", [](const Atom& atom) { return atom.z(); }, "Return z coordinate")
      .def("distance",
           [](const Atom& atom, const Atom& other) { return atom.distance(other); },
           nb::arg("other"), "Return spatial distance to atom")
      .def("distance",
           [](const Atom& atom, const Coordinates& coords) { return atom.distance(coords); },
           nb::arg("coords"), "Return spatial distance to point")
      .def("connections",
           [](const Atom& atom, atom_number_t atom_number) {
             std::vector<int> result;
             result.reserve(atom.ncon());
             for (const Bond* bond : atom) {
               result.push_back(bond->other(atom_number));
             }
             return result;
           },
           nb::arg("atom_number"))
      .def("__contains__",
           [](const Atom& atom, atom_number_t other) {
             return static_cast<bool>(atom.is_bonded_to(other));
           })
      .def("__len__", [](const Atom& atom) { return atom.ncon(); })
      .def("__getitem__",
           [](const Atom& atom, int index) { return atom[index]; },
           nb::rv_policy::reference_internal)
      .def("__iter__",
           [](const Atom& atom) {
             return nb::make_iterator<nb::rv_policy::reference_internal>(
                 nb::type<Atom>(), "BondIterator", atom.begin(), atom.end());
           },
           nb::keep_alive<0, 1>())
      .def("__sub__",
           [](const Atom& atom, const Atom& other) { return atom.distance(other); },
           nb::arg("other"))
      .def("__repr__", &AtomRepr)
      .def("__str__", &AtomRepr);

  nb::class_<Bond>(m, "Bond")
      .def(nb::init<>())
      .def("a1", &Bond::a1)
      .def("a2", &Bond::a2)
      .def("other", &Bond::other, nb::arg("atom"))
      .def("involves",
           [](const Bond& bond, atom_number_t atom) {
             return static_cast<bool>(bond.involves(atom));
           },
           nb::arg("atom"))
      .def("is_directional", &Bond::is_directional)
      .def("is_single_bond", [](const Bond& bond) { return static_cast<bool>(bond.is_single_bond()); })
      .def("is_double_bond", [](const Bond& bond) { return static_cast<bool>(bond.is_double_bond()); })
      .def("is_triple_bond", [](const Bond& bond) { return static_cast<bool>(bond.is_triple_bond()); })
      .def("is_aromatic_bond", [](const Bond& bond) { return static_cast<bool>(bond.is_aromatic()); })
      .def("is_aromatic", [](const Bond& bond) { return static_cast<bool>(bond.is_aromatic()); })
      .def("btype", &ToBondType)
      .def("nrings", &BondNrings)
      .def("IsInRing", [](const Bond& bond) { return BondNrings(bond) > 0; })
      .def("bond_number_assigned", &Bond::bond_number_assigned)
      .def("bond_number", &Bond::bond_number)
      .def("GetBeginAtomIdx", &Bond::a1)
      .def("GetEndAtomIdx", &Bond::a2)
      .def("GetBondType", &ToBondType)
      .def("__contains__",
           [](const Bond& bond, atom_number_t atom) { return static_cast<bool>(bond.involves(atom)); })
      .def("__repr__", &BondRepr)
      .def("__str__", &BondRepr);


}

}  // namespace lillymol_nb
