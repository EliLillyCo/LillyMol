#include "nanobind/lillymol_nb_internal.h"

#include "Molecule_Lib/mol2graph.h"

#include <queue>

#include <nanobind/ndarray.h>

namespace lillymol_nb {

namespace {

bool
LooksLikeSdfTagRecord(const std::string& line) {
  const size_t open = line.find('<');
  const size_t close = line.rfind('>');
  return !line.empty() && line[0] == '>' && open != std::string::npos &&
         close != std::string::npos && open < close;
}

std::string
SdfTagName(const std::string& line) {
  const size_t open = line.find('<');
  const size_t close = line.rfind('>');
  if (open == std::string::npos || close == std::string::npos || open >= close) {
    return std::string();
  }

  std::string result = line.substr(open + 1, close - open - 1);
  std::replace(result.begin(), result.end(), ' ', '_');
  return result;
}

nb::dict
SdfTags(const Molecule& mol) {
  nb::dict result;

  std::string current_tag;
  std::string current_value;

  auto flush_current_tag = [&]() {
    if (current_tag.empty()) {
      return;
    }
    result[nb::str(current_tag.c_str())] = nb::str(current_value.c_str());
    current_tag.clear();
    current_value.clear();
  };

  for (int i = 0; i < mol.number_records_text_info(); ++i) {
    const std::string line = mol.text_info(i).AsString();
    if (LooksLikeSdfTagRecord(line)) {
      flush_current_tag();
      current_tag = SdfTagName(line);
      continue;
    }

    if (current_tag.empty()) {
      continue;
    }

    if (line.empty()) {
      flush_current_tag();
      continue;
    }

    if (!current_value.empty()) {
      current_value.append("\n");
    }
    current_value.append(line);
  }

  flush_current_tag();
  return result;
}

std::vector<int>
RingAtoms(const Ring& ring) {
  return std::vector<int>(ring.rawdata(), ring.rawdata() + ring.size());
}

std::vector<std::vector<int>>
RingInfoAtomRings(RingInfo& ring_info) {
  Molecule& mol = *ring_info.mol;
  const int nrings = mol.nrings();
  std::vector<std::vector<int>> result;
  result.reserve(nrings);
  for (int i = 0; i < nrings; ++i) {
    result.push_back(RingAtoms(*mol.ringi(i)));
  }
  return result;
}

using FloatNumpyArray1D = nb::ndarray<nb::numpy, float, nb::shape<-1>>;
using IntNumpyArray1D = nb::ndarray<nb::numpy, int, nb::shape<-1>>;

FloatNumpyArray1D
MakeFloatNumpyArray1D(size_t n) {
  float* data = new float[n];
  nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
  return FloatNumpyArray1D(data, {n}, owner);
}

IntNumpyArray1D
MakeIntNumpyArray1D(size_t n) {
  int* data = new int[n];
  nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<int*>(p); });
  return IntNumpyArray1D(data, {n}, owner);
}

FloatNumpyArray1D
FloatArrayFromCoordinates(const float* coords, size_t n) {
  float* data = new float[n];
  std::copy_n(coords, n, data);
  nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
  return FloatNumpyArray1D(data, {n}, owner);
}

FloatNumpyArray1D
GetCoordinatesNumpy(const Molecule& mol) {
  const int matoms = mol.natoms();
  FloatNumpyArray1D result = MakeFloatNumpyArray1D(3 * matoms);
  float* data = result.data();
  for (int i = 0; i < matoms; ++i) {
    data[3 * i + 0] = mol.x(i);
    data[3 * i + 1] = mol.y(i);
    data[3 * i + 2] = mol.z(i);
  }
  return result;
}

void
SetCoordinatesNumpy(Molecule& mol, FloatNumpyArray1D coords) {
  const int expected = 3 * mol.natoms();
  if (coords.ndim() != 1 || static_cast<int>(coords.shape(0)) != expected) {
    throw std::invalid_argument("set_coordinates_numpy requires 3 values per atom");
  }
  mol.SetCoordinates(coords.data());
}

int
ToScaffoldInPlace(Molecule& mol) {
  const int matoms = mol.natoms();
  if (matoms <= 3) {
    return 0;
  }
  if (mol.nrings() == 0) {
    return 0;
  }

  std::unique_ptr<int[]> spinach = std::make_unique<int[]>(matoms);
  mol.identify_spinach(spinach.get());
  if (std::count(spinach.get(), spinach.get() + matoms, 0) == matoms) {
    return 0;
  }

  for (int i = 0; i < matoms; ++i) {
    const Atom& atom = mol[i];
    if (atom.ncon() != 1) {
      continue;
    }
    const Bond* bond = atom[0];
    if (!bond->is_double_bond()) {
      continue;
    }
    const atom_number_t other = bond->other(i);
    if (spinach[other] == 0) {
      spinach[i] = 0;
    }
  }

  return mol.remove_atoms(spinach.get(), 1);
}

Molecule
ScaffoldMolecule(const Molecule& mol) {
  Molecule result(mol);
  ToScaffoldInPlace(result);
  return result;
}

std::vector<int>
RingBondNumbers(Molecule& mol, const Ring& ring) {
  mol.assign_bond_numbers_to_bonds_if_needed();
  std::vector<int> result;
  result.reserve(ring.size());
  const int ring_size = ring.size();
  for (int i = 0; i < ring_size; ++i) {
    const atom_number_t a1 = ring[i];
    const atom_number_t a2 = ring[(i + 1) % ring_size];
    const Bond* bond = mol.bond_between_atoms_if_present(a1, a2);
    if (bond != nullptr && bond->bond_number_assigned()) {
      result.push_back(bond->bond_number());
    }
  }
  return result;
}

std::vector<std::vector<int>>
RingInfoBondRings(RingInfo& ring_info) {
  Molecule& mol = *ring_info.mol;
  const int nrings = mol.nrings();
  std::vector<std::vector<int>> result;
  result.reserve(nrings);
  for (int i = 0; i < nrings; ++i) {
    result.push_back(RingBondNumbers(mol, *mol.ringi(i)));
  }
  return result;
}

int
RingInfoNumBondRings(RingInfo& ring_info, int bond_number) {
  Molecule& mol = *ring_info.mol;
  if (bond_number < 0 || bond_number >= mol.nedges()) {
    throw std::invalid_argument("bond number outside molecule");
  }
  const Bond* bond = mol.bondi(bond_number);
  int result = 0;
  const int nrings = mol.nrings();
  for (int i = 0; i < nrings; ++i) {
    if (mol.ringi(i)->contains_bond(bond->a1(), bond->a2())) {
      ++result;
    }
  }
  return result;
}

bool
RingInfoAreBondsInSameRing(RingInfo& ring_info, int bond1, int bond2) {
  Molecule& mol = *ring_info.mol;
  if (bond1 < 0 || bond1 >= mol.nedges() || bond2 < 0 || bond2 >= mol.nedges()) {
    throw std::invalid_argument("bond number outside molecule");
  }
  const Bond* b1 = mol.bondi(bond1);
  const Bond* b2 = mol.bondi(bond2);
  const int nrings = mol.nrings();
  for (int i = 0; i < nrings; ++i) {
    const Ring* ring = mol.ringi(i);
    if (ring->contains_bond(b1->a1(), b1->a2()) &&
        ring->contains_bond(b2->a1(), b2->a2())) {
      return true;
    }
  }
  return false;
}

std::vector<std::vector<int>>
AtomsByRadius(const Molecule& mol, const Set_of_Atoms& starting_atoms, int max_radius) {
  if (max_radius < 0) {
    throw std::invalid_argument("max_radius must be non-negative");
  }

  std::vector<std::vector<int>> result(max_radius + 1);
  const int matoms = mol.natoms();
  std::vector<int> distance(matoms, -1);
  std::queue<atom_number_t> to_process;

  for (atom_number_t atom : starting_atoms) {
    if (atom < 0 || atom >= matoms) [[unlikely]] {
      throw std::invalid_argument("starting atom is outside the molecule");
    }
    if (distance[atom] >= 0) {
      continue;
    }

    distance[atom] = 0;
    result[0].push_back(atom);
    to_process.push(atom);
  }

  while (!to_process.empty()) {
    const atom_number_t atom = to_process.front();
    to_process.pop();

    const int next_radius = distance[atom] + 1;
    if (next_radius > max_radius) {
      continue;
    }

    for (const Bond* bond : mol[atom]) {
      const atom_number_t other = bond->other(atom);
      if (distance[other] >= 0) {
        continue;
      }

      distance[other] = next_radius;
      result[next_radius].push_back(other);
      to_process.push(other);
    }
  }

  return result;
}

}  // namespace

void
BindMolecule(nb::module_& m) {
  nb::class_<RingInfo>(m, "RingInfo")
      .def("NumRings", [](RingInfo& ring_info) { return ring_info.mol->nrings(); })
      .def("num_rings", [](RingInfo& ring_info) { return ring_info.mol->nrings(); })
      .def("AtomRings", &RingInfoAtomRings)
      .def("atom_rings", &RingInfoAtomRings)
      .def("BondRings", &RingInfoBondRings)
      .def("bond_rings", &RingInfoBondRings)
      .def("NumAtomRings",
           [](RingInfo& ring_info, atom_number_t atom) { return ring_info.mol->nrings(atom); },
           nb::arg("atom"))
      .def("num_atom_rings",
           [](RingInfo& ring_info, atom_number_t atom) { return ring_info.mol->nrings(atom); },
           nb::arg("atom"))
      .def("NumBondRings", &RingInfoNumBondRings, nb::arg("bond"))
      .def("num_bond_rings", &RingInfoNumBondRings, nb::arg("bond"))
      .def("IsAtomInRingOfSize",
           [](RingInfo& ring_info, atom_number_t atom, int ring_size) {
             return ring_info.mol->in_ring_of_given_size(atom, ring_size);
           },
           nb::arg("atom"), nb::arg("ring_size"))
      .def("is_atom_in_ring_of_size",
           [](RingInfo& ring_info, atom_number_t atom, int ring_size) {
             return ring_info.mol->in_ring_of_given_size(atom, ring_size);
           },
           nb::arg("atom"), nb::arg("ring_size"))
      .def("AreAtomsInSameRing",
           [](RingInfo& ring_info, atom_number_t a1, atom_number_t a2) {
             return static_cast<bool>(ring_info.mol->in_same_ring(a1, a2));
           },
           nb::arg("a1"), nb::arg("a2"))
      .def("are_atoms_in_same_ring",
           [](RingInfo& ring_info, atom_number_t a1, atom_number_t a2) {
             return static_cast<bool>(ring_info.mol->in_same_ring(a1, a2));
           },
           nb::arg("a1"), nb::arg("a2"))
      .def("AreBondsInSameRing", &RingInfoAreBondsInSameRing,
           nb::arg("bond1"), nb::arg("bond2"))
      .def("are_bonds_in_same_ring", &RingInfoAreBondsInSameRing,
           nb::arg("bond1"), nb::arg("bond2"));

  nb::class_<Coordinates>(m, "Coordinates")
      .def(nb::init<>())
      .def(nb::init<coord_t, coord_t, coord_t>(), nb::arg("x"), nb::arg("y"), nb::arg("z"))
      .def("x", nb::overload_cast<>(&Coordinates::x, nb::const_))
      .def("y", nb::overload_cast<>(&Coordinates::y, nb::const_))
      .def("z", nb::overload_cast<>(&Coordinates::z, nb::const_))
      .def("setxyz", static_cast<void (Coordinates::*)(coord_t, coord_t, coord_t)>(&Coordinates::setxyz), nb::arg("x"), nb::arg("y"), nb::arg("z"))
      .def("set_x", &Coordinates::set_x, nb::arg("x"))
      .def("set_y", &Coordinates::set_y, nb::arg("y"))
      .def("set_z", &Coordinates::set_z, nb::arg("z"))
      .def("normalise", &Coordinates::normalise)
      .def("norm", &Coordinates::norm)
      .def("length", &Coordinates::length)
      .def("distance",
           [](const Coordinates& coords, const Coordinates& other) { return coords.distance(other); },
           nb::arg("other"))
      .def("dot_product",
           [](const Coordinates& coords, const Coordinates& other) { return coords.dot_product(other); },
           nb::arg("other"))
      .def("__repr__",
           [](const Coordinates& coords) {
             IWString result;
             result << '(' << coords.x() << ',' << coords.y() << ',' << coords.z() << ')';
             return result.AsString();
           })
      .def("__str__",
           [](const Coordinates& coords) {
             IWString result;
             result << '(' << coords.x() << ',' << coords.y() << ',' << coords.z() << ')';
             return result.AsString();
           });

  nb::class_<Molecule>(m, "Molecule")
      .def(nb::init<>())
      .def(nb::init<const Molecule&>(), nb::arg("other"))
      .def("copy", [](const Molecule& mol) { return Molecule(mol); },
           "Return an independent copy of the molecule")
      .def("build_from_smiles",
           [](Molecule& mol, const std::string& smiles) {
             return mol.build_from_smiles(smiles);
           },
           nb::arg("smiles"))
      .def("ok", [](const Molecule& mol) { return static_cast<bool>(mol.ok()); })
      .def("natoms", nb::overload_cast<>(&Molecule::natoms, nb::const_),
           "Return the number of atoms")
      .def("natoms",
           [](const Molecule& mol, atomic_number_t atomic_number) {
             return mol.natoms(atomic_number);
           },
           nb::arg("atomic_number"),
           "number of atoms with the specified atomic number")
      .def("natoms",
           [](const Molecule& mol, const std::string& symbol) {
             const_IWSubstring tmp(symbol);
             const Element* element = get_element_from_symbol_no_case_conversion(tmp);
             if (element == nullptr) {
               throw std::invalid_argument("Unrecognised element");
             }
             return mol.natoms(element);
           },
           nb::arg("symbol"),
           "Number of atoms of type atomic symbol")
      .def("GetNumAtoms", [](const Molecule& mol) { return mol.natoms(); })
      .def("GetNumHeavyAtoms",
           [](const Molecule& mol) { return mol.natoms() - mol.natoms(1); },
           "Return the number of non-hydrogen atoms")
      .def("GetNumBonds", &Molecule::nedges, "Return the number of bonds")
      .def("empty", [](const Molecule& mol) { return static_cast<bool>(mol.empty()); })
      .def("resize", &Molecule::resize, nb::arg("natoms"))
      .def("add_atom",
           [](Molecule& mol, int atomic_number) {
             const Element* element = get_element_from_atomic_number(atomic_number);
             if (element == nullptr) {
               throw std::invalid_argument("Invalid atomic number");
             }
             mol.add(element);
             return mol.natoms() - 1;
           },
           nb::arg("atomic_number"))
      .def("atom", &Molecule::atomi, nb::arg("atom"), nb::rv_policy::reference_internal)
      .def("GetAtomWithIdx", &Molecule::atomi, nb::arg("atom"),
           nb::rv_policy::reference_internal)
      .def("GetAtoms",
           [](const Molecule& mol) {
             return nb::make_iterator<nb::rv_policy::reference_internal>(
                 nb::type<Molecule>(), "AtomIterator", mol.begin(), mol.end());
           },
           nb::keep_alive<0, 1>(), "Return an iterator over atoms")
      .def("nedges", &Molecule::nedges,
           "Return the number of bonds")
      .def("bond", nb::overload_cast<int>(&Molecule::bondi, nb::const_), nb::arg("bond"),
           nb::rv_policy::reference_internal)
      .def("GetBondWithIdx", nb::overload_cast<int>(&Molecule::bondi, nb::const_),
           nb::arg("bond"), nb::rv_policy::reference_internal)
      .def("bonds", &Bonds, nb::rv_policy::reference_internal)
      .def("GetBonds", &Bonds, nb::rv_policy::reference_internal)
      .def("nrings", [](Molecule& mol) { return mol.nrings(); },
           "Return the number of SSSR rings")
      .def("get_ring_membership", &GetRingMembership,
           "Return SSSR ring membership counts for all atoms")
      .def("ring_membership", [](Molecule& mol) { mol.ring_membership(); },
           "Force ring membership perception")
      .def("is_ring_atom", [](Molecule& mol, atom_number_t atom) {
             return static_cast<bool>(mol.is_ring_atom(atom));
           },
           nb::arg("atom"))
      .def("IsInRing", [](Molecule& mol, atom_number_t atom) {
             return mol.ring_bond_count(atom) > 0;
           },
           nb::arg("atom"))
      .def("in_ring_of_size", &Molecule::in_ring_of_given_size,
           nb::arg("atom"), nb::arg("ring_size"))
      .def("IsAtomInRingOfSize", &Molecule::in_ring_of_given_size,
           nb::arg("atom"), nb::arg("ring_size"))
      .def("NumAtomRings", nb::overload_cast<atom_number_t>(&Molecule::nrings),
           nb::arg("atom"))
      .def("ring_bond_count", nb::overload_cast<atom_number_t>(&Molecule::ring_bond_count),
           nb::arg("atom"))
      .def("ring_bond_counts", &RingBondCount)
      .def("fused_system_identifier", &Molecule::fused_system_identifier,
           nb::arg("atom"))
      .def("fused_system_size", &Molecule::fused_system_size, nb::arg("atom"))
      .def("number_ring_systems", nb::overload_cast<>(&Molecule::number_ring_systems))
      .def("ring", nb::overload_cast<int>(&Molecule::ringi), nb::arg("ring"),
           nb::rv_policy::reference_internal)
      .def("non_sssr_rings", &Molecule::non_sssr_rings,
           "Number non SSSR rings")
      .def("non_sssr_ring", &Molecule::non_sssr_ring, nb::arg("ring"),
           nb::rv_policy::reference_internal, "Fetch the ith non SSSR ring")
      .def("rings", &Rings)
      .def("ring_info", [](Molecule& mol) { return RingInfo(&mol); },
           nb::keep_alive<0, 1>(), "Return a RingInfo view")
      .def("GetRingInfo", [](Molecule& mol) { return RingInfo(&mol); },
           nb::keep_alive<0, 1>(), "Return a RingInfo view")
      .def("in_same_ring",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2) {
             return static_cast<bool>(mol.in_same_ring(a1, a2));
           },
           nb::arg("a1"), nb::arg("a2"))
      .def("in_same_ring_system",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2) {
             return static_cast<bool>(mol.in_same_ring_system(a1, a2));
           },
           nb::arg("a1"), nb::arg("a2"))
      .def("largest_ring_size", &Molecule::LargestRingSize)
      .def("is_spiro_fused", nb::overload_cast<atom_number_t>(&Molecule::is_spiro_fused),
           nb::arg("atom"), "True if atom is spiro fused")
      .def("label_atoms_by_atom_number",
           [](Molecule& mol) {
             const int matoms = mol.natoms();
             for (int i = 1; i < matoms; ++i) {
               mol.set_isotope(i, i);
             }
           },
           "Set isotope labels to atom numbers")
      .def("label_atoms_by_ring_system",
           [](Molecule& mol) {
             std::vector<int> result(mol.natoms());
             mol.label_atoms_by_ring_system(result.data());
             return result;
           },
           "Return ring system identifier for each atom")
      .def("label_atoms_by_ring_system_including_spiro_fused",
           [](Molecule& mol) {
             std::vector<int> result(mol.natoms());
             mol.label_atoms_by_ring_system_including_spiro_fused(result.data());
             return result;
           },
           "Return ring system identifier for each atom, including spiro fusion")
      .def("label_atoms_by_ring_system_including_spiro_fused_np",
           [](Molecule& mol) {
             IntNumpyArray1D result = MakeIntNumpyArray1D(mol.natoms());
             mol.label_atoms_by_ring_system_including_spiro_fused(result.data());
             return result;
           },
           "Return ring system identifier for each atom as a NumPy array")
      .def("number_fragments", [](Molecule& mol) { return mol.number_fragments(); },
           "Return the number of fragments")
      .def("fragment_membership", nb::overload_cast<atom_number_t>(&Molecule::fragment_membership),
           nb::arg("atom"), "Return the fragment containing atom")
      .def("get_fragment_membership",
           [](Molecule& mol) {
             std::vector<int> result(mol.natoms());
             mol.fragment_membership(result.data());
             return result;
           },
           "Return fragment membership for each atom")
      .def("atoms_in_fragment", nb::overload_cast<int>(&Molecule::atoms_in_fragment),
           nb::arg("fragment"), "Return number of atoms in fragment")
      .def("delete_fragment", nb::overload_cast<int>(&Molecule::delete_fragment),
           nb::arg("fragment"), "Remove a fragment")
      .def("remove_fragment", nb::overload_cast<int>(&Molecule::delete_fragment),
           nb::arg("fragment"), "Remove a fragment")
      .def("delete_all_fragments_except",
           [](Molecule& mol, int frag) {
             if (mol.empty()) {
               throw std::invalid_argument("molecule contains no atoms");
             }
             const int nf = mol.number_fragments();
             if (nf == 1) {
               throw std::invalid_argument("molecule contains only one fragment");
             }
             if (frag < 0 || frag >= nf) {
               throw std::invalid_argument("fragment index outside [0, number_fragments)");
             }
             return mol.delete_all_fragments_except(frag);
           },
           nb::arg("frag"),
           "Delete all fragments except the zero-based fragment index frag, in place. "
           "Return the number of atoms removed. "
           "Raise ValueError for empty or single-fragment molecules or an invalid fragment index.")
      .def("remove_fragment_containing_atom", &Molecule::remove_fragment_containing_atom,
           nb::arg("atom"), "Remove the fragment containing atom")
      .def("reduce_to_largest_fragment", nb::overload_cast<>(&Molecule::reduce_to_largest_fragment),
           "Strip to the largest fragment")
      .def("reduce_to_largest_fragment_carefully",
           nb::overload_cast<>(&Molecule::reduce_to_largest_fragment_carefully),
           "Strip to the largest fragment with heuristic fragment selection")
      .def("create_components",
           [](Molecule& mol) -> std::optional<std::vector<Molecule>> {
             if (mol.number_fragments() < 2) {
               return std::nullopt;
             }
             resizable_array_p<Molecule> components;
             if (!mol.create_components(components)) {
               return std::nullopt;
             }
             std::vector<Molecule> result;
             result.reserve(components.number_elements());
             for (const Molecule* component : components) {
               result.push_back(*component);
             }
             return result;
           },
           "Split a multi-fragment molecule into component molecules")
      .def("atoms_in_largest_fragment", [](Molecule& mol) { return mol.atoms_in_largest_fragment(); },
           "Return the number of atoms in the largest fragment")
      .def("remove_non_periodic_table_elements", &Molecule::remove_all_non_natural_elements,
           "Remove non periodic table elements")
      .def("organic_only", nb::overload_cast<>(&Molecule::organic_only, nb::const_),
           "True if only organic elements are present")
      .def("non_organic_atom_count", &Molecule::non_organic_atom_count,
           "Return the number of non-organic atoms")
      .def("is_organic",
           [](const Molecule& mol, atom_number_t atom) { return static_cast<bool>(mol.is_organic(atom)); },
           nb::arg("atom"), "True if atom is organic")
      .def("remove_explicit_hydrogens", nb::overload_cast<>(&Molecule::remove_explicit_hydrogens),
           "Remove explicit hydrogens")
      .def("atomic_number", &Molecule::atomic_number, nb::arg("atom"),
           "Return the atomic number for atom")
      .def("set_atomic_number", &Molecule::set_atomic_number, nb::arg("atom"),
           nb::arg("atomic_number"), "Set the atomic number for atom")
      .def("add_bond",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2, BondType bond_type) {
             return mol.add_bond(a1, a2, ToLillyMolBondType(bond_type));
           },
           nb::arg("a1"), nb::arg("a2"), nb::arg("bond_type"), "Add bond")
      .def("set_bond_type_between_atoms",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2, BondType bond_type) {
             return mol.set_bond_type_between_atoms(a1, a2, ToLillyMolBondType(bond_type));
           },
           nb::arg("a1"), nb::arg("a2"), nb::arg("bond_type"), "Set bond type")
      .def("add", [](Molecule& mol, const Molecule& other) { return mol.add_molecule(&other); },
           nb::arg("other"), "Add molecule as a separate fragment")
      .def("assign_bond_numbers_to_bonds", &Molecule::assign_bond_numbers_to_bonds)
      .def("bond_between_atoms",
           [](const Molecule& mol, atom_number_t a1, atom_number_t a2) -> const Bond* {
             return mol.bond_between_atoms_if_present(a1, a2);
           },
           nb::arg("a1"), nb::arg("a2"), nb::rv_policy::reference_internal)
      .def("bond_type_between_atoms",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2) -> std::optional<BondType> {
             mol.compute_aromaticity_if_needed();
             const Bond* bond = mol.bond_between_atoms_if_present(a1, a2);
             if (bond == nullptr) {
               return std::nullopt;
             }
             return ToBondType(*bond);
           },
           nb::arg("a1"), nb::arg("a2"))
      .def("remove_atom", nb::overload_cast<atom_number_t>(&Molecule::remove_atom),
           nb::arg("atom"))
      .def("remove_atoms",
           [](Molecule& mol, const std::vector<int>& to_remove, int flag) {
             return mol.remove_atoms(to_remove.data(), flag);
           },
           nb::arg("to_remove"), nb::arg("flag"))
      .def("remove_atoms",
           [](Molecule& mol, Set_of_Atoms atoms) {
             return mol.remove_atoms(atoms);
           },
           nb::arg("atoms"))
      .def("remove_bond_between_atoms", &Molecule::remove_bond_between_atoms,
           nb::arg("a1"), nb::arg("a2"))
      .def("remove_all_bonds", &Molecule::remove_all_bonds)
      .def("remove_bonds_to_atom",
           [](Molecule& mol, atom_number_t atom) { return static_cast<bool>(mol.remove_bonds_to_atom(atom, 0)); },
           nb::arg("atom"), "Remove all bonds to an atom")
      .def("remove_edge", nb::overload_cast<int>(&Molecule::remove_bond),
           nb::arg("bond"), "Remove a bond by bond number")
      .def("chop", &Molecule::chop, nb::arg("natoms"), "Remove the last atoms")
      .def("remove_all", nb::overload_cast<atomic_number_t>(&Molecule::remove_all),
           nb::arg("atomic_number"))
      .def("atomic_numbers", &AtomicNumbers,
           "Return atomic numbers for all atoms")
      .def("ncon", nb::overload_cast<atom_number_t>(&Molecule::ncon, nb::const_),
           nb::arg("atom"), "Return the connection count for atom")
      .def("connections",
           [](const Molecule& mol, atom_number_t atom) {
             Set_of_Atoms connections = mol.connections(atom);
             return std::vector<int>(connections.rawdata(), connections.rawdata() + connections.size());
           },
           nb::arg("atom"))
      .def("other_atom",
           [](const Molecule& mol, atom_number_t atom, int connection) {
             return mol.other(atom, connection);
           },
           nb::arg("atom"), nb::arg("connection"))
      .def("nbonds", nb::overload_cast<atom_number_t>(&Molecule::nbonds, nb::const_),
           nb::arg("atom"), "Return the bond order sum for atom")
      .def("atomic_symbol",
           [](const Molecule& mol, atom_number_t atom) { return mol.atomic_symbol(atom).AsString(); },
           nb::arg("atom"))
      .def("attached_heteroatom_count", nb::overload_cast<atom_number_t>(&Molecule::attached_heteroatom_count, nb::const_),
           nb::arg("atom"), "Return number of attached heteroatoms")
      .def("is_halogen", &Molecule::is_halogen, nb::arg("atom"), "True if atom is a halogen")
      .def("smarts_equivalent_for_atom",
           [](Molecule& mol, atom_number_t atom) { return mol.smarts_equivalent_for_atom(atom).AsString(); },
           nb::arg("atom"), "Return SMARTS equivalent for atom")
      .def("smarts", [](Molecule& mol) { return mol.smarts().AsString(); },
           "Return molecule as SMARTS")
      .def("is_aromatic",
           [](Molecule& mol, atom_number_t atom) { return static_cast<bool>(mol.IsAromatic(atom)); },
           nb::arg("atom"))
      .def("hybridization",
           [](Molecule& mol, atom_number_t atom) {
             if (!mol.ok_atom_number(atom)) {
               throw std::invalid_argument("hybridization atom number outside [0, natoms)");
             }
             return HybridizationState(mol, atom);
           },
           nb::arg("atom"), "RDKit-like hybridization of atom, computed on demand")
      .def("find_kekule_form",
           [](Molecule& mol, std::vector<int>& atoms) { return static_cast<bool>(mol.find_kekule_form(atoms.data())); },
           nb::arg("atoms"), "Find a Kekule form for atoms")
      .def("pi_electrons",
           [](Molecule& mol, atom_number_t atom) {
             int pi = 0;
             mol.pi_electrons(atom, pi);
             return pi;
           },
           nb::arg("atom"), "Return pi electron count for atom")
      .def("lone_pair_count",
           [](Molecule& mol, atom_number_t atom) {
             int lone_pairs = 0;
             mol.lone_pair_count(atom, lone_pairs);
             return lone_pairs;
           },
           nb::arg("atom"), "Return lone pair count for atom")
      .def("compute_aromaticity_if_needed", &Molecule::compute_aromaticity_if_needed,
           "Ensure aromaticity has been perceived")
      .def("aromatic_atom_count", nb::overload_cast<>(&Molecule::aromatic_atom_count))
      .def("aromatic_ring_count", nb::overload_cast<>(&Molecule::aromatic_ring_count))
      .def("number_chiral_centres", nb::overload_cast<>(&Molecule::chiral_centres, nb::const_))
      .def("remove_all_chiral_centres", &Molecule::remove_all_chiral_centres)
      .def("revert_all_directional_bonds_to_non_directional",
           &Molecule::revert_all_directional_bonds_to_non_directional,
           "Remove directional bond annotations")
      .def("chiral_centre_at_atom",
           [](const Molecule& mol, atom_number_t atom) -> std::optional<Chiral_Centre> {
             Chiral_Centre* centre = mol.chiral_centre_at_atom(atom);
             if (centre == nullptr) {
               return std::nullopt;
             }
             return *centre;
           },
           nb::arg("atom"))
      .def("invert_chirality_on_atom", &Molecule::invert_chirality_on_atom,
           nb::arg("atom"))
      .def("remove_chiral_centre_at_atom", &Molecule::remove_chiral_centre_at_atom,
           nb::arg("atom"))
      .def("chiral_centres",
           [](const Molecule& mol) {
             std::vector<Chiral_Centre> result;
             result.reserve(mol.chiral_centres());
             for (const Chiral_Centre* centre : mol.ChiralCentres()) {
               result.push_back(*centre);
             }
             return result;
           })
      .def("isotope", &Molecule::isotope, nb::arg("atom"),
           "Return the isotope value for atom")
      .def("isotopes", &Isotopes,
           "Return isotope values for all atoms")
      .def("set_isotope",
           [](Molecule& mol, atom_number_t atom, isotope_t isotope) {
             return mol.set_isotope(atom, isotope);
           },
           nb::arg("atom"), nb::arg("isotope"))
      .def("set_isotopes",
           [](Molecule& mol, const Set_of_Atoms& atoms, isotope_t isotope) {
             return mol.set_isotope(atoms, isotope);
           },
           nb::arg("atoms"), nb::arg("isotope"))
      .def("set_isotopes",
           [](Molecule& mol, const std::vector<int>& atoms, isotope_t isotope) {
             return mol.set_isotope(atoms, isotope);
           },
           nb::arg("atoms"), nb::arg("isotope"))
      .def("set_isotopes",
           [](Molecule& mol, IntNumpyArray1D isotopes) {
             if (isotopes.ndim() != 1 || static_cast<int>(isotopes.shape(0)) != mol.natoms()) {
               throw std::invalid_argument("set_isotopes requires one isotope value per atom");
             }
             const int* data = isotopes.data();
             for (int i = 0; i < mol.natoms(); ++i) {
               mol.set_isotope(i, data[i]);
             }
           },
           nb::arg("isotopes"))
      .def("remove_isotopes",
           [](Molecule& mol, int unset_implicit_h) {
             return mol.unset_isotopes(unset_implicit_h);
           },
           nb::arg("unset_implicit_h") = 1)
      .def("transform_to_non_isotopic_form",
           [](Molecule& mol, int unset_implicit_h) {
             return mol.transform_to_non_isotopic_form(unset_implicit_h);
           },
           nb::arg("unset_implicit_h") = 1)
      .def("number_isotopic_atoms", nb::overload_cast<>(&Molecule::number_isotopic_atoms, nb::const_))
      .def("first_atom_with_isotope",
           [](const Molecule& mol, isotope_t isotope) { return mol.atom_with_isotope(isotope); },
           nb::arg("isotope"))
      .def("formal_charge", nb::overload_cast<atom_number_t>(&Molecule::formal_charge, nb::const_),
           nb::arg("atom"))
      .def("set_formal_charge", nb::overload_cast<atom_number_t, formal_charge_t>(&Molecule::set_formal_charge),
           nb::arg("atom"), nb::arg("formal_charge"))
      .def("has_formal_charges", nb::overload_cast<>(&Molecule::has_formal_charges, nb::const_))
      .def("number_formal_charges", nb::overload_cast<>(&Molecule::number_formally_charged_atoms, nb::const_))
      .def("net_formal_charge", nb::overload_cast<>(&Molecule::net_formal_charge, nb::const_))
      .def("molecular_formula",
           [](Molecule& mol) { return mol.molecular_formula().AsString(); },
           "Return the molecular formula")
      .def("amw",
           [](const Molecule& mol) { return lillymol::MolecularWeightIsotopesAsLabels(mol); },
           "Return average molecular weight, treating isotope labels as labels")
      .def("amw_ignore_isotopes", &Molecule::molecular_weight_ignore_isotopes,
           "Return average molecular weight with isotope labels ignored strictly")
      .def("exact_mass", [](Molecule& mol) { return mol.exact_mass(); },
           "Return exact mass")
      .def("implicit_hydrogens", nb::overload_cast<atom_number_t>(&Molecule::implicit_hydrogens),
           nb::arg("atom"))
      .def("explicit_hydrogens", nb::overload_cast<atom_number_t>(&Molecule::explicit_hydrogens, nb::const_),
           nb::arg("atom"))
      .def("hcount",
           [](Molecule& mol, atom_number_t atom) { return mol.hcount(atom); },
           nb::arg("atom"))
      .def("implicit_hydrogens_known", &Molecule::implicit_hydrogens_known,
           nb::arg("atom"))
      .def("make_implicit_hydrogens_explicit", nb::overload_cast<>(&Molecule::make_implicit_hydrogens_explicit))
      .def("unset_all_implicit_hydrogen_information", &Molecule::unset_all_implicit_hydrogen_information,
           "Discard all implicit hydrogen known flags")
      .def("remove_hydrogens_known_flag_to_fix_valence_errors",
           &Molecule::remove_hydrogens_known_flag_to_fix_valence_errors,
           "Remove problematic implicit hydrogen known flags")
      .def("AddHs", [](Molecule& mol) { return mol.make_implicit_hydrogens_explicit(); })
      .def("RemoveHs", [](Molecule& mol) { return mol.remove_all(1); })
      .def("to_scaffold", &ToScaffoldInPlace, "Convert to scaffold in place")
      .def("scaffold", &ScaffoldMolecule, "Return a new Molecule containing the scaffold")
      .def("change_to_graph_form", nb::overload_cast<>(&Molecule::change_to_graph_form),
           "Convert to default graph form in place")
      .def("to_graph",
           [](Molecule& mol, const Mol2Graph& mol2graph) { return mol.change_to_graph_form(mol2graph); },
           nb::arg("mol2graph"), "Convert to graph form under Mol2Graph control")
      .def("valence_ok", nb::overload_cast<>(&Molecule::valence_ok))
      .def("valence_ok", nb::overload_cast<atom_number_t>(&Molecule::valence_ok), nb::arg("atom"))
      .def("lipinski_num_h_donors", &Molecule::LipinskiNumHDonors,
           "Lipinski hydrogen bond donor count")
      .def("lipinski_num_h_acceptors", &Molecule::LipinskiNumHAcceptors,
           "Lipinski hydrogen bond acceptor count")
      .def("rdkit_num_h_donors", &Molecule::RDKitNumHDonors,
           "RDKit compatible hydrogen bond donor count")
      .def("rdkit_num_h_acceptors", &Molecule::RDKitNumHAcceptors,
           "RDKit compatible hydrogen bond acceptor count")
      .def("saturated",
           [](Molecule& mol, atom_number_t atom) { return static_cast<bool>(mol.saturated(atom)); },
           nb::arg("atom"))
      .def("unsaturation", &Molecule::unsaturation, nb::arg("atom"))
      .def("smiles",
           [](Molecule& mol) { return mol.smiles().AsString(); },
           "Return a non-unique SMILES")
      .def("aromatic_smiles",
           [](Molecule& mol) { return mol.aromatic_smiles().AsString(); },
           "Return an aromatic SMILES")
      .def("unique_smiles",
           [](Molecule& mol) { return mol.unique_smiles().AsString(); },
           "Return the unique SMILES")
      .def("random_smiles",
           [](Molecule& mol) { return mol.random_smiles().AsString(); },
           "Return a random SMILES")
      .def("isotopically_labelled_smiles",
           [](Molecule& mol) { return mol.isotopically_labelled_smiles().AsString(); },
           "Return SMILES with isotopes as atom numbers")
      .def("unique_kekule_smiles",
           [](Molecule& mol) { return mol.UniqueKekuleSmiles().AsString(); },
           "Return unique Kekule SMILES")
      .def("smiles_atom_order",
           [](Molecule& mol) {
             std::vector<int> result(mol.natoms());
             mol.smiles_atom_order(result.data());
             return result;
           },
           "Return atom order from the most recent SMILES generation")
      .def("sort_atoms",
           [](Molecule& mol, const std::vector<int>& order) {
             if (static_cast<int>(order.size()) != mol.natoms()) {
               throw std::invalid_argument("sort_atoms requires one value per atom");
             }
             static constexpr int kAscending = 1;
             return mol.sort(order.data(), kAscending);
           },
           nb::arg("order"), "Sort atoms by ascending values in order")
      .def("compute_distance_matrix", &Molecule::recompute_distance_matrix,
           "Compute or recompute the topological distance matrix")
      .def("renumber_atoms",
           [](Molecule& mol, const std::vector<int>& new_number) {
             const int matoms = mol.natoms();
             if (static_cast<int>(new_number.size()) != matoms) {
               throw std::invalid_argument("renumber_atoms requires one entry for each atom");
             }
             std::vector<int> seen(matoms, 0);
             for (int i = 0; i < matoms; ++i) {
               const int destination = new_number[i];
               if (destination < 0 || destination >= matoms) {
                 throw std::invalid_argument("renumber_atoms mapping contains an atom number outside [0, natoms)");
               }
               if (seen[destination]) {
                 throw std::invalid_argument("renumber_atoms mapping contains duplicate atom numbers");
               }
               seen[destination] = 1;
             }
             return mol.renumber_atoms(new_number.data());
           },
           nb::arg("new_number"), "Renumber atoms")
      .def("smiles_starting_with_atom",
           [](Molecule& mol, atom_number_t atom) { return mol.smiles_starting_with_atom(atom).AsString(); },
           nb::arg("atom"), "Return SMILES starting at atom")
      .def("name",
           [](const Molecule& mol) { return mol.name().AsString(); },
           "Return the molecule name")
      .def("set_name",
           [](Molecule& mol, const std::string& name) { mol.set_name(name); },
           nb::arg("name"))
      .def("number_records_text_info", &Molecule::number_records_text_info,
           "Return number of retained SDF/text records")
      .def("text_info",
           [](const Molecule& mol) {
             std::vector<std::string> result;
             result.reserve(mol.number_records_text_info());
             for (int i = 0; i < mol.number_records_text_info(); ++i) {
               result.push_back(mol.text_info(i).AsString());
             }
             return result;
           },
           "Return retained SDF/text records")
      .def("sdf_tags", &SdfTags, "Return retained SDF tags as a dict")
      .def("are_bonded", nb::overload_cast<atom_number_t, atom_number_t>(&Molecule::are_bonded, nb::const_),
           nb::arg("a1"), nb::arg("a2"), "True if atoms are bonded")
      .def("bonds_between", nb::overload_cast<atom_number_t, atom_number_t>(&Molecule::bonds_between),
           nb::arg("a1"), nb::arg("a2"), "Return topological distance between atoms")
      .def("longest_path", nb::overload_cast<>(&Molecule::longest_path),
           "Return the longest topological path")
      .def("most_distant_pair",
           [](Molecule& mol) {
             atom_number_t a1 = INVALID_ATOM_NUMBER;
             atom_number_t a2 = INVALID_ATOM_NUMBER;
             int longest_distance = 0;
             const int matoms = mol.natoms();
             for (int i = 0; i < matoms; ++i) {
               for (int j = i + 1; j < matoms; ++j) {
                 if (mol.fragment_membership(i) != mol.fragment_membership(j)) {
                   continue;
                 }
                 const int distance = mol.bonds_between(i, j);
                 if (distance > longest_distance) {
                   a1 = i;
                   a2 = j;
                   longest_distance = distance;
                 }
               }
             }
             return std::make_tuple(a1, a2);
           },
           "Return the most separated atom pair")
      .def("atoms_on_shortest_path",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2) -> std::optional<std::vector<int>> {
             Set_of_Atoms atoms;
             if (!mol.atoms_between(a1, a2, atoms) || atoms.empty()) {
               return std::nullopt;
             }
             return std::vector<int>(atoms.rawdata(), atoms.rawdata() + atoms.size());
           },
           nb::arg("a1"), nb::arg("a2"), "Return atoms on the shortest path between a1 and a2")
      .def("all_atoms_between",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2) -> std::optional<std::vector<int>> {
             Set_of_Atoms atoms;
             if (!mol.AllAtomsBetween(a1, a2, atoms) || atoms.empty()) {
               return std::nullopt;
             }
             return std::vector<int>(atoms.rawdata(), atoms.rawdata() + atoms.size());
           },
           nb::arg("a1"), nb::arg("a2"), "Return all atoms on shortest paths between a1 and a2")
      .def("down_the_bond",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2) -> std::optional<std::vector<int>> {
             const int matoms = mol.natoms();
             std::unique_ptr<int[]> dtb = std::make_unique<int[]>(matoms);
             std::optional<int> maybe_n = mol.DownTheBond(a1, a2, dtb.get());
             if (!maybe_n) {
               return std::nullopt;
             }
             std::vector<int> result;
             result.reserve(*maybe_n);
             for (int i = 0; i < matoms; ++i) {
               if (i != a2 && dtb[i]) {
                 result.push_back(i);
               }
             }
             return result;
           },
           nb::arg("a1"), nb::arg("a2"), "Return atoms down the bond from a1 to a2")
      .def("atoms_by_radius", &AtomsByRadius, nb::arg("starting_atoms"), nb::arg("max_radius"),
           "Return atom shells by minimum bond distance from starting_atoms")
      .def("reset_atom_map_numbers", &Molecule::reset_all_atom_map_numbers,
           "Reset atom map numbers")
      .def("set_atom_map_number",
           static_cast<void (Molecule::*)(atom_number_t, int)>(&Molecule::set_atom_map_number),
           nb::arg("atom"), nb::arg("map_number"), "Set atom map number")
      .def("atom_map_number",
           nb::overload_cast<atom_number_t>(&Molecule::atom_map_number, nb::const_),
           nb::arg("atom"), "Return atom map number")
      .def("atom_with_atom_map_number", &Molecule::atom_with_atom_map_number,
           nb::arg("map_number"), "Return atom with atom map number")
      .def("canonical_rank", nb::overload_cast<atom_number_t>(&Molecule::canonical_rank),
           nb::arg("atom"), "Return atom canonical rank")
      .def("canonical_ranks",
           [](Molecule& mol) {
             std::vector<int> result(mol.natoms());
             mol.canonical_ranks(result.data());
             return result;
           },
           "Return canonical ranks for all atoms")
      .def("symmetry_class", nb::overload_cast<atom_number_t>(&Molecule::symmetry_class),
           nb::arg("atom"), "Return atom symmetry class")
      .def("number_symmetry_classes", nb::overload_cast<>(&Molecule::number_symmetry_classes),
           "Return number of symmetry classes")
      .def("symmetry_equivalents",
           [](Molecule& mol, atom_number_t atom) {
             Set_of_Atoms tmp;
             mol.symmetry_equivalents(atom, tmp);
             return std::vector<int>(tmp.rawdata(), tmp.rawdata() + tmp.size());
           },
           nb::arg("atom"), "Return atoms symmetry-equivalent to atom")
      .def("partial_charge_type",
           [](const Molecule& mol) { return mol.partial_charge_type().AsString(); },
           "Return the type of partial charges stored")
      .def("invalidate_partial_charges",
           [](Molecule& mol) { mol.invalidate_charges(); },
           "Discard partial charge information")
      .def("partial_charge",
           [](const Molecule& mol, atom_number_t atom) { return mol.partial_charge(atom); },
           nb::arg("atom"), "Return partial charge on atom")
      .def("compute_Abraham_partial_charges",
           [](Molecule& mol) { return mol.compute_Abraham_partial_charges(); },
           "Compute Abraham partial charges")
      .def("compute_Gasteiger_partial_charges",
           [](Molecule& mol) { return mol.compute_Gasteiger_partial_charges(); },
           "Compute Gasteiger partial charges")
      .def("compute_Huckel_partial_charges",
           [](Molecule& mol) { return mol.compute_Huckel_partial_charges(); },
           "Compute Huckel partial charges")
      .def("compute_Gasteiger_Huckel_partial_charges",
           [](Molecule& mol) { return mol.compute_Gasteiger_Huckel_partial_charges(); },
           "Compute Gasteiger-Huckel partial charges")
      .def("gasteiger_partial_charges",
           [](Molecule& mol) -> std::vector<float> {
             mol.compute_Gasteiger_partial_charges();
             std::vector<float> result;
             result.reserve(mol.natoms());
             for (int i = 0; i < mol.natoms(); ++i) {
               result.push_back(mol.partial_charge(i));
             }
             return result;
           },
           "Return list of Gasteiger partial charges")
      .def("move_to_end_of_connection_table", &Molecule::MoveToEndOfConnectionTable,
           nb::arg("atomic_number"), "Move atoms of atomic_number to end of connection table")
      .def("x", nb::overload_cast<atom_number_t>(&Molecule::x, nb::const_),
           nb::arg("atom"), "Return x coordinate")
      .def("y", nb::overload_cast<atom_number_t>(&Molecule::y, nb::const_),
           nb::arg("atom"), "Return y coordinate")
      .def("z", nb::overload_cast<atom_number_t>(&Molecule::z, nb::const_),
           nb::arg("atom"), "Return z coordinate")
      .def("setx", static_cast<void (Molecule::*)(atom_number_t, coord_t)>(&Molecule::setx),
           nb::arg("atom"), nb::arg("x"), "Set x coordinate")
      .def("sety", static_cast<void (Molecule::*)(atom_number_t, coord_t)>(&Molecule::sety),
           nb::arg("atom"), nb::arg("y"), "Set y coordinate")
      .def("setz", static_cast<void (Molecule::*)(atom_number_t, coord_t)>(&Molecule::setz),
           nb::arg("atom"), nb::arg("z"), "Set z coordinate")
      .def("setxyz",
           static_cast<void (Molecule::*)(atom_number_t, coord_t, coord_t, coord_t)>(&Molecule::setxyz),
           nb::arg("atom"), nb::arg("x"), nb::arg("y"), nb::arg("z"),
           "Set atom coordinates")
      .def("get_coordinates",
           [](const Molecule& mol) {
             const int matoms = mol.natoms();
             std::vector<float> result;
             result.reserve(matoms * 3);
             for (int i = 0; i < matoms; ++i) {
               result.push_back(mol.x(i));
               result.push_back(mol.y(i));
               result.push_back(mol.z(i));
             }
             return result;
           },
           "Return coordinates as [x0, y0, z0, x1, y1, z1, ...]")
      .def("get_coordinates_numpy", &GetCoordinatesNumpy,
           "Return coordinates as a flat float32 NumPy array")
      .def("set_coordinates",
           [](Molecule& mol, const std::vector<float>& coords) {
             const int expected = 3 * mol.natoms();
             if (static_cast<int>(coords.size()) != expected) {
               throw std::invalid_argument("set_coordinates requires 3 values per atom");
             }
             mol.SetCoordinates(coords.data());
           },
           nb::arg("coords"), "Set coordinates from [x0, y0, z0, x1, y1, z1, ...]")
      .def("set_coordinates_numpy", &SetCoordinatesNumpy, nb::arg("coords"),
           "Set coordinates from a flat float32 NumPy array")
      .def("dihedral_scan",
           [](Molecule& mol, atom_number_t a2, atom_number_t a3, angle_t angle, float bump_check) {
             std::vector<std::unique_ptr<float[]>> coords = mol.DihedralScan(a2, a3, angle, bump_check);
             std::vector<FloatNumpyArray1D> result;
             result.reserve(coords.size());
             const size_t ncoords = static_cast<size_t>(mol.natoms()) * 3;
             for (const std::unique_ptr<float[]>& coord : coords) {
               result.push_back(FloatArrayFromCoordinates(coord.get(), ncoords));
             }
             return result;
           },
           nb::arg("a2"), nb::arg("a3"), nb::arg("angle"), nb::arg("bump_check") = 0.0f,
           "Scan around the a2-a3 bond and return flat float32 NumPy coordinate arrays")
      .def("translate",
           [](Molecule& mol, coord_t x, coord_t y, coord_t z) {
             mol.translate_atoms(x, y, z);
           },
           nb::arg("x"), nb::arg("y"), nb::arg("z"), "Translate all atoms")
      .def("translate",
           [](Molecule& mol, const Coordinates& delta) { mol.translate_atoms(delta); },
           nb::arg("delta"), "Translate all atoms")
      .def("translate",
           [](Molecule& mol, const std::vector<int>& to_move, int flag,
              coord_t x, coord_t y, coord_t z) {
             if (static_cast<int>(to_move.size()) != mol.natoms()) {
               throw std::invalid_argument("translate selection requires one value per atom");
             }
             const Coordinates delta(x, y, z);
             mol.translate_atoms(delta, to_move.data(), flag);
           },
           nb::arg("to_move"), nb::arg("flag"), nb::arg("x"), nb::arg("y"), nb::arg("z"),
           "Translate selected atoms where to_move[i] == flag")
      .def("rotate",
           [](Molecule& mol, const Coordinates& axis, angle_t angle) {
             return mol.rotate_atoms(axis, angle);
           },
           nb::arg("axis"), nb::arg("angle"), "Rotate all atoms around axis")
      .def("rotate",
           [](Molecule& mol, const std::vector<int>& to_move, int flag,
              const Coordinates& axis, angle_t angle) {
             if (static_cast<int>(to_move.size()) != mol.natoms()) {
               throw std::invalid_argument("rotate selection requires one value per atom");
             }
             return mol.rotate_atoms(axis, angle, to_move.data(), flag);
           },
           nb::arg("to_move"), nb::arg("flag"), nb::arg("axis"), nb::arg("angle"),
           "Rotate selected atoms where to_move[i] == flag")
      .def("bond_length",
           [](const Molecule& mol, atom_number_t a1, atom_number_t a2) -> std::optional<float> {
             if (!mol.are_bonded(a1, a2)) {
               return std::nullopt;
             }
             return mol.bond_length(a1, a2);
           },
           nb::arg("a1"), nb::arg("a2"), "Return distance between bonded atoms")
      .def("bond_angle",
           [](const Molecule& mol, atom_number_t centre, atom_number_t a1, atom_number_t a2) {
             return mol.bond_angle(centre, a1, a2, BondedStatus::kOkNotBonded);
           },
           nb::arg("centre"), nb::arg("a1"), nb::arg("a2"),
           "Return angle defined by three atoms")
      .def("dihedral_angle",
           [](const Molecule& mol, atom_number_t a1, atom_number_t a2,
              atom_number_t a3, atom_number_t a4) {
             return mol.dihedral_angle(a1, a2, a3, a4, BondedStatus::kOkNotBonded);
           },
           nb::arg("a1"), nb::arg("a2"), nb::arg("a3"), nb::arg("a4"),
           "Return dihedral angle defined by four atoms")
      .def("signed_dihedral_angle",
           nb::overload_cast<atom_number_t, atom_number_t, atom_number_t, atom_number_t>(
               &Molecule::signed_dihedral_angle, nb::const_),
           nb::arg("a1"), nb::arg("a2"), nb::arg("a3"), nb::arg("a4"),
           "Return signed dihedral angle defined by four atoms")
      .def("distance_between_atoms",
           nb::overload_cast<atom_number_t, atom_number_t>(
               &Molecule::distance_between_atoms, nb::const_),
           nb::arg("a1"), nb::arg("a2"), "Return spatial distance between atoms")
      .def("longest_intra_molecular_distance", &Molecule::longest_intra_molecular_distance,
           "Return longest spatial distance within the molecule")
      .def("bump_check",
           [](const Molecule& mol, distance_t dist) { return mol.bump_check(dist); },
           nb::arg("dist"), "Return number of non-bonded atom pairs closer than dist")
      .def("highest_coordinate_dimensionality", &Molecule::highest_coordinate_dimensionality,
           "Return highest coordinate dimensionality present")
      .def("discern_chirality_from_3d_structure", &Molecule::discern_chirality_from_3d_structure,
           "Perceive chiral centres from 3D coordinates")
      .def("debug_string", &Molecule::debug_string,
           "Return a dump of internal data structures")
      .def("__repr__",
           [](Molecule& mol) {
             IWString formula;
             mol.isis_like_molecular_formula(formula);
             IWString result;
             result << '<' << mol.name() << " with " << mol.natoms() << " atoms " << formula << '>';
             return result.AsString();
           })
      .def("__str__",
           [](Molecule& mol) {
             IWString result;
             result << mol.smiles() << ' ' << mol.name();
             return result.AsString();
           })
      .def("__iadd__",
           [](Molecule& mol, const Molecule& rhs) -> Molecule& {
             mol += rhs;
             return mol;
           },
           nb::rv_policy::reference_internal)
      .def("__add__",
           [](const Molecule& lhs, const Molecule& rhs) { return lhs + rhs; })
      .def("__eq__",
           [](Molecule& lhs, Molecule& rhs) { return lhs == rhs; },
           nb::arg("rhs"))
      .def("__len__", [](const Molecule& mol) { return mol.natoms(); })
      .def("__getitem__",
           [](const Molecule& mol, int index) { return mol[index]; },
           nb::rv_policy::reference_internal)
      .def("__iter__",
           [](const Molecule& mol) {
             return nb::make_iterator<nb::rv_policy::reference_internal>(
                 nb::type<Molecule>(), "AtomIterator", mol.begin(), mol.end());
           },
           nb::keep_alive<0, 1>())
      .def("__copy__", [](const Molecule& mol) { return Molecule(mol); })
      .def("__contains__", &MoleculeContainsAtomicNumber)
      .def("__contains__", &MoleculeContainsElementSymbol)
      .def("__contains__", &MoleculeContainsSubstructureQuery);

  m.def("set_copy_name_in_molecule_copy_constructor",
        &set_copy_name_in_molecule_copy_constructor, nb::arg("value"),
        "Control whether Molecule copy construction copies the molecule name");
  m.def("Position3D",
        [](Molecule& mol, atom_number_t atom1, float distance, atom_number_t atom2) {
          return lillymol::Position3D(mol, atom1, distance, atom2);
        },
        nb::arg("mol"), nb::arg("atom1"), nb::arg("distance"), nb::arg("atom2"),
        "Move the fragment containing atom2 so atom1 and atom2 are distance apart");
}

}  // namespace lillymol_nb
