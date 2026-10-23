#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {

void
BindSetOfAtomsAndRing(nb::module_& m) {
  nb::class_<Set_of_Atoms>(m, "Set_of_Atoms")
      .def(nb::init<>())
      .def(nb::init<const std::vector<int>&>(), nb::arg("atoms"))
      .def("empty", [](const Set_of_Atoms& atoms) { return atoms.empty(); })
      .def("any_atoms_in_common",
           [](const Set_of_Atoms& lhs, const Set_of_Atoms& rhs) {
             return lhs.any_members_in_common(rhs);
           },
           nb::arg("rhs"))
      .def("first_atom_in_common", &Set_of_Atoms::first_member_in_common,
           nb::arg("rhs"))
      .def("atoms_in_common", &Set_of_Atoms::members_in_common,
           nb::arg("rhs"))
      .def("size", &Set_of_Atoms::size)
      .def("as_list", &Set_of_Atoms::AsVector)
      .def("set_vector",
           [](const Set_of_Atoms& atoms, nb::list values, int value) {
             const Py_ssize_t n = PyList_Size(values.ptr());
             for (atom_number_t atom : atoms) {
               if (atom < 0 || atom >= n) {
                 throw std::out_of_range("set_vector atom index out of range");
               }
               values[atom] = nb::int_(value);
             }
           },
           nb::arg("values"), nb::arg("value"),
           "Set entries in a Python list for atoms in this embedding")
      .def("scatter",
           [](const Set_of_Atoms& atoms, nb::list values, int value) {
             const Py_ssize_t n = PyList_Size(values.ptr());
             for (atom_number_t atom : atoms) {
               if (atom < 0 || atom >= n) {
                 throw std::out_of_range("scatter atom index out of range");
               }
               values[atom] = nb::int_(value);
             }
           },
           nb::arg("values"), nb::arg("value"),
           "Set entries in a Python list for atoms in this embedding")
      .def("increment_vector",
           [](const Set_of_Atoms& atoms, nb::list values, int value) {
             const Py_ssize_t n = PyList_Size(values.ptr());
             for (atom_number_t atom : atoms) {
               if (atom < 0 || atom >= n) {
                 throw std::out_of_range("increment_vector atom index out of range");
               }
               values[atom] = nb::int_(nb::cast<int>(values[atom]) + value);
             }
           },
           nb::arg("values"), nb::arg("value"),
           "Increment entries in a Python list for atoms in this embedding")
      .def("contains_both",
           [](const Set_of_Atoms& atoms, atom_number_t a1, atom_number_t a2) {
             return atoms.contains_atoms(a1, a2);
           },
           nb::arg("a1"), nb::arg("a2"))
      .def("__len__", [](const Set_of_Atoms& atoms) { return atoms.size(); })
      .def("__contains__",
           [](const Set_of_Atoms& atoms, atom_number_t atom) {
             return atoms.contains(atom);
           })
      .def("__getitem__",
           [](const Set_of_Atoms& atoms, int index) { return atoms.item(index); })
      .def("__iter__",
           [](const Set_of_Atoms& atoms) {
             return nb::make_iterator(nb::type<Set_of_Atoms>(), "SetOfAtomsIterator",
                                      atoms.begin(), atoms.end());
           },
           nb::keep_alive<0, 1>())
      .def("__setitem__",
           [](Set_of_Atoms& atoms, int index, atom_number_t atom) {
             atoms[index] = atom;
             return atom;
           })
      .def("__repr__", &SetOfAtomsRepr)
      .def("__str__", &SetOfAtomsStr)
      .def("__eq__",
           [](const Set_of_Atoms& lhs, const std::vector<int>& rhs) {
             return lhs.AsVector() == rhs;
           })
      .def("append", [](Set_of_Atoms& atoms, atom_number_t atom) { atoms.add(atom); },
           nb::arg("atom"))
      .def("extend",
           [](Set_of_Atoms& atoms, const std::vector<int>& extra) {
             for (int atom : extra) {
               atoms << atom;
             }
           },
           nb::arg("atoms"))
      .def("__iadd__",
           [](Set_of_Atoms& lhs, const Set_of_Atoms& rhs) -> Set_of_Atoms& {
             lhs += rhs;
             return lhs;
           },
           nb::arg("rhs"))
      .def("__iadd__",
           [](Set_of_Atoms& lhs, atom_number_t atom) -> Set_of_Atoms& {
             lhs.add(atom);
             return lhs;
           },
           nb::arg("atom"));

  nb::class_<Ring, Set_of_Atoms>(m, "Ring")
      .def(nb::init<>())
      .def("size", &Ring::size)
      .def("ring_number", &Ring::ring_number)
      .def("fragment_membership", &Ring::fragment_membership)
      .def("fused_system_identifier", &Ring::fused_system_identifier)
      .def("is_fused", [](const Ring& ring) { return static_cast<bool>(ring.is_fused()); })
      .def("is_fused_to",
           [](const Ring& lhs, const Ring& rhs) { return static_cast<bool>(lhs.is_fused_to(&rhs)); },
           nb::arg("rhs"))
      .def("fused_ring_neighbours", &Ring::fused_ring_neighbours)
      .def("largest_number_of_bonds_shared_with_another_ring",
           &Ring::largest_number_of_bonds_shared_with_another_ring)
      .def("strongly_fused_ring_neighbours", &Ring::strongly_fused_ring_neighbours)
      .def("is_aromatic", [](const Ring& ring) { return static_cast<bool>(ring.is_aromatic()); },
           "True if the ring is aromatic")
      .def("contains_bond",
           [](const Ring& ring, atom_number_t a1, atom_number_t a2) {
             return static_cast<bool>(ring.contains_bond(a1, a2));
           },
           nb::arg("a1"), nb::arg("a2"))
      .def("contains_both",
           [](const Ring& ring, atom_number_t a1, atom_number_t a2) {
             return static_cast<bool>(ring.contains_both(a1, a2));
           },
           nb::arg("a1"), nb::arg("a2"));


}

}  // namespace lillymol_nb
