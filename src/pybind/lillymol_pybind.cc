#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

#include "pybind11/pybind11.h"
// to convert C++ STL containers to python list
#include "pybind11/stl.h"

#ifdef LILLYMOL_VECTOR_OPAQUE
#include "pybind11/stl_bind.h"
#endif

#include "pybind11/operators.h"

#ifdef LILLYMOL_VECTOR_OPAQUE
// This is not a great idea, because it destroys the otherwise
// easy interoperability between vector<int> and a python List.
// Return to this sometime and see if I can figure it out...
PYBIND11_MAKE_OPAQUE(std::vector<int>);
#endif

#include "Molecule_Lib/chiral_centre.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/mol2graph.pb.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/substructure.h"
#include "pybind/molecule.h"

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

static std::string
RingToString(const Ring& r) {
  IWString result;
  result << "<Ring size " << r.size();
  if (r.is_aromatic()) {
    result << " arom";
  } else if (r.is_non_aromatic()) {
    result << " aliph";
  }
  result << " frag " << r.fragment_membership();
  result << " fsid " << r.fused_system_identifier();
  result << " fused " << r.is_fused();
  result << " atoms";
  for (int i = 0; i < r.number_elements(); ++i) {
    result << ' ' << r[i];
  }

  result << '>';

  return std::string(result.data(), result.size());
}

static int
ToScaffold(Molecule& m) {
  const int matoms = m.natoms();
  if (matoms <= 3) {
    return 0;
  }
  if (m.nrings() == 0) {
    return 0;
  }
  std::unique_ptr<int[]> spinach = std::make_unique<int[]>(matoms);
  m.identify_spinach(spinach.get());
  if (std::count(spinach.get(), spinach.get() + matoms, 0) == matoms) {
    return 0;
  }

  return m.remove_atoms(spinach.get(), 1);
}

enum BondType {
  kUnknown = 0,
  kSingleBond = 1,
  kDoubleBond = 2,
  kTripleBond = 3,
  kAromaticBond = 4
};

PYBIND11_MODULE(lillymol, m) 
{
#ifdef LILLYMOL_VECTOR_OPAQUE
  py::bind_vector<std::vector<int>>(m, "VectorInt");
#endif

  py::class_<Chiral_Centre>(m, "Chiral_Centre")
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
    .def("move_to_end_of_connection_table",
      [](Molecule& m, atomic_number_t z) {
        return m.MoveToEndOfConnectionTable(z);
      },
      "Move atoms of type 'z' to end of connection table"
    )

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

  py::class_<Mol2Graph>(m, "Mol2Graph")
    .def(py::init<>())
    .def("set_exclude_triple_bonds_from_graph_reduction", &Mol2Graph::set_exclude_triple_bonds_from_graph_reduction, "exclude_triple_bonds_from_graph_reduction")
    .def("set_revert_all_directional_bonds_to_non_directional", &Mol2Graph::set_revert_all_directional_bonds_to_non_directional, "revert_all_directional_bonds_to_non_directional")
    .def("set_preserve_cc_double_bonds_no_heteroatoms", &Mol2Graph::set_preserve_cc_double_bonds_no_heteroatoms, "set_preserve_cc_double_bonds_no_heteroatoms")
    .def("set_preserve_cc_double_bonds_saturated", &Mol2Graph::set_preserve_cc_double_bonds_saturated, "preserve_cc_double_bonds_saturated")
    .def("set_append_molecular_formula", &Mol2Graph::set_append_molecular_formula, "set_append_molecular_formula")
    .def("set_aromatic_distinguishing_formula", &Mol2Graph::set_aromatic_distinguishing_formula, "set_aromatic_distinguishing_formula")
    .def("set_remove_chiral_centres", &Mol2Graph::set_remove_chiral_centres, "set_remove_chiral_centres")
    .def("set_active", &Mol2Graph::set_active, "Set active")
    .def("active", &Mol2Graph::active, "True if active")
  ;

         py::class_<Molecule>(m, "Molecule")
		.def(py::init<>())
                .def(py::init([](const Molecule& rhs) {
                  return std::unique_ptr<Molecule>(new Molecule(rhs));
                }))
                .def("ok",
                  [](const Molecule& m)->bool{
                    return m.ok();
                  },
                  "Returns true if internal datastructures ok"
                )
                .def("natoms", static_cast<int (Molecule::*)()const>(&Molecule::natoms), "Number atoms in molecule")
                .def("natoms",
                  [](const Molecule& m, atomic_number_t z) {
                    return m.natoms(z);
                  },
                 "number of atoms with atomic number"
                )
                .def("natoms",
                  [](const Molecule& m, const std::string& asymbol) {
                    const const_IWSubstring tmp(asymbol);
                    const Element* e = get_element_from_symbol_no_case_conversion(tmp);
                    if (e == nullptr) {
                      throw py::value_error("Unrecognised element");;
                    }
                    return m.natoms(e);
                  },
                 "number of atoms with element"
                )
                .def("GetNumAtoms",
                  [](const Molecule& m) {
                    return m.natoms();
                  },
                  "natoms"
                )
                .def("empty",
                  [](const Molecule& m) -> bool {
                    return m.empty();
                  },
                  "True if molecule is empty"
                )
                .def("nedges", static_cast<int (Molecule::*)()const>(&Molecule::nedges), "Number edges in molecule")
                .def("add_atom",
                  [](Molecule& m, int atnum)->bool {
                    const Element* e = get_element_from_atomic_number(atnum);
                    if (e == nullptr) {
                      std::cerr << "Invalid atomic number " << atnum << '\n';
                      return false;
                    }
                    return m.add(e);
                  },
                  "Add an atom with atomic number"
                )
                .def("nrings", static_cast<int (Molecule::*)()>(&Molecule::nrings), "Number rings in molecule")
                .def("nrings", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::nrings), "Rings containing atom")
                .def("non_sssr_rings", &Molecule::non_sssr_rings, "Number non SSSR rings")
                .def("non_sssr_ring", &Molecule::non_sssr_ring, "Fetch the ith non SSSR ring")
                .def("IsInRing",
                  [](Molecule& m, atom_number_t zatom)->bool{
                    return m.ring_bond_count(zatom) > 0;
                  },
                  "True if atom in ring"
                )
                .def("in_ring_of_given_size", &Molecule::in_ring_of_given_size, "True if atom in ring of give size")
                .def("IsAtomInRingOfSize",
                  [](Molecule& m, atom_number_t zatom, int rsize)->bool{
                    return m.in_ring_of_given_size(zatom, rsize);
                  },
                  "True if atom is in a ring size rsize"
                )
                .def("NumAtomRings",
                  [](Molecule& m, atom_number_t zatom)->int{
                    return m.nrings(zatom);
                  },
                  "nrings"
                )

                .def("is_ring_atom",
                  [](Molecule& m, atom_number_t zatom)->bool {
                    return m.is_ring_atom(zatom);
                  },
                  "true if atom is in a ring"
                )
                .def("ring_bond_count", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::ring_bond_count), "Ring bonds at atom")
                .def("get_ring_membership",
                  [](Molecule& m)->std::vector<int> {
                    std::vector<int> rc(m.natoms());
                    m.ring_membership(rc.data());
                    return rc;
                  },
                  "Ring membership for all atoms"
                )
                .def("fused_system_identifier", &Molecule::fused_system_identifier, "Fused system identifier")
                .def("fused_system_size", &Molecule::fused_system_size, "Fused system size")
                .def("number_ring_systems", static_cast<int (Molecule::*)()>(&Molecule::number_ring_systems), "Number ring systems")
                .def("ring", static_cast<const Ring* (Molecule::*)(int)>(&Molecule::ringi), py::return_value_policy::reference, "I'th ring")
                .def("rings",
                  [](Molecule& m)->std::vector<const Ring*>{
                    std::vector<const Ring*> result;
                    // Was not able to make this work passing a vector of the actual rings.
                    // Works once or twice then crashes. Does this create a memory leak?
                    for (const Ring* r : m.sssr_rings()) {
                      result.push_back(new Ring(*r));
                    }
                    return result;
                  }
                )
                //.def("in_same_ring", static_cast<int (Molecule::*)(atom_number_t, atom_number_t)>(&Molecule::in_same_ring), "True if atoms in the same ring")
                .def("in_same_ring",
                  [](Molecule& m, atom_number_t a1, atom_number_t a2)->bool{
                    return m.in_same_ring(a1, a2);
                  },
                  "true if atoms in same ring"
                )
                //.def("in_same_ring_system", static_cast<int (Molecule::*)(atom_number_t, atom_number_t)>(&Molecule::in_same_ring_system), "True if atoms in same ring system")
                .def("in_same_ring_system",
                  [](Molecule& m, atom_number_t a1, atom_number_t a2)->bool{
                    return m.in_same_ring_system(a1, a2);
                  },
                  "True if atoms in same ring system"
                )
                .def("largest_ring_size", static_cast<int (Molecule::*)()>(&Molecule::LargestRingSize), "Largest ring size")
                .def("number_ring_systems", static_cast<int (Molecule::*)()>(&Molecule::number_ring_systems), "Number ring systems")
                .def("is_spiro_fused", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::is_spiro_fused), "True if atom is spiro fused")
                .def("label_atoms_by_ring_system",
                  [](Molecule& m)->std::vector<int>{
                    std::vector<int> rc(m.natoms());
                    m.label_atoms_by_ring_system(rc.data());
                    return rc;
                  },
                  "Return a list of ring system identifiers"
                )
                .def("label_atoms_by_ring_system_including_spiro_fused",
                  [](Molecule& m)->std::vector<int>{
                    std::vector<int> rc(m.natoms());
                    m.label_atoms_by_ring_system_including_spiro_fused(rc.data());
                    return rc;
                  },
                  "Return a list of ring system identifiers"
                )
                .def("amw", static_cast<float (Molecule::*)()const>(&Molecule::molecular_weight), "AMW")
                .def("exact_mass", static_cast<exact_mass_t (Molecule::*)()const>(&Molecule::exact_mass), "Exact Mass")
                .def("ncon", static_cast<int (Molecule::*)(atom_number_t)const>(&Molecule::ncon), "Connections to Atom")
                .def("connections",
                  [](const Molecule& m, atom_number_t zatom)->std::vector<int>{
                    const Set_of_Atoms s = m.connections(zatom);
                    return std::vector<int>(s.rawdata(), s.rawdata() + s.size());
                  },
                  "Atoms connected to atom"
                )
                .def("attached_heteroatom_count", static_cast<int (Molecule::*)(atom_number_t)const>(&Molecule::attached_heteroatom_count), "Number of heteroatoms attached")
                //.def("is_aromatic", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::is_aromatic), "True if atom is aromatic")
                .def("is_aromatic",
                  [](Molecule& m, atom_number_t zatom)->bool{
                    return m.IsAromatic(zatom);
                  },
                  "True if atom is aromatic"
                )
                .def("find_kekule_form",
                  [](Molecule& m, std::vector<int>& atoms)->bool {
                    return m.find_kekule_form(atoms.data());
                  },
                  "Find a Kekule form for 'atoms'"
                )
                .def("pi_electrons",
                  [](Molecule& m, atom_number_t zatom) {
                    int pi;
                    m.pi_electrons(zatom, pi);
                    return pi;
                  },
                  "Pi electrons"
                )
                .def("lone_pair_count", &Molecule::lone_pair_count, "Lone pair count")
                .def("compute_aromaticity_if_needed", &Molecule::compute_aromaticity_if_needed, "Ensure molecule has aromaticity")
                .def("aromatic_atom_count", static_cast<int (Molecule::*)()>(&Molecule::aromatic_atom_count), "Number aromatic atoms")
                .def("aromatic_ring_count", static_cast<int (Molecule::*)()>(&Molecule::aromatic_ring_count), "Number aromatic rings")
                .def("atomic_number", static_cast<atomic_number_t (Molecule::*)(atom_number_t)const>(&Molecule::atomic_number), "atomic number of atom")
                .def("atomic_symbol",
                  [](const Molecule& m, atom_number_t zatom)->std::string{
                    const IWString& s = m.atomic_symbol(zatom);
                    return std::string(s.data(), s.size());
                  },
                  "Atomic symbol for atom"
                )
                .def("is_halogen", &Molecule::is_halogen, "True if atom is Halogen")
                .def("smarts_equivalent_for_atom",
                  [](Molecule& m, atom_number_t zatom)->std::string{
                    const IWString s = m.smarts_equivalent_for_atom(zatom);
                    return std::string(s.data(), s.size());
                  },
                  "smarts for atom"
                )
                .def("smarts",
                  [](Molecule& m)->std::string{
                    return m.smarts().AsString();
                  },
                  "Molecule as smarts - not robust for substructure searching"
                )
                .def("set_atomic_number", static_cast<int (Molecule::*)(atom_number_t, atomic_number_t)>(&Molecule::set_atomic_number), "set atomic number")
                .def("add_bond", static_cast<int (Molecule::*)(atom_number_t, atom_number_t, bond_type_t, int)>(&Molecule::add_bond), "add bond between atoms")
                .def("add_bond", 
                  [] (Molecule& m, atom_number_t a1, atom_number_t a2, BondType bt) {
                    switch (bt) {
                      case BondType::kSingleBond:
                        return m.add_bond(a1, a2, SINGLE_BOND);
                      case BondType::kDoubleBond:
                        return m.add_bond(a1, a2, DOUBLE_BOND);
                      case BondType::kTripleBond:
                        return m.add_bond(a1, a2, TRIPLE_BOND);
                      default:
                        throw py::value_error("add_bond:unrecognised bond type");;
                    }
                  },
                  "add bond between atoms"
                )
                .def("set_bond_type_between_atoms", 
                  [](Molecule& m, atom_number_t a1, atom_number_t a2, BondType bt) {
                    switch (bt) {
                      case BondType::kSingleBond:
                        return m.set_bond_type_between_atoms(a1, a2, SINGLE_BOND);
                      case BondType::kDoubleBond:
                        return m.set_bond_type_between_atoms(a1, a2, DOUBLE_BOND);
                      case BondType::kTripleBond:
                        return m.set_bond_type_between_atoms(a1, a2, TRIPLE_BOND);
                      default:
                        throw py::value_error("Unrecognised bond type");;
                    }
                  },
                  "Set bond type"
                )
                .def("assign_bond_numbers_to_bonds", &Molecule::assign_bond_numbers_to_bonds, "Assign unique id to each bond")
                .def("remove_atom", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::remove_atom), "Remove an atom")
                .def("remove_atoms",
                  [](Molecule& m, const std::vector<int>& to_remove, int flag) {
                    return m.remove_atoms(to_remove.data(), flag);
                  },
                  "Remove atoms value `flag` in `to_remove`"
                )
                .def("remove_atoms",
                  // Note we pass `s` by value since it gets altered.
                  [](Molecule&m, Set_of_Atoms s) {
                    return m.remove_atoms(s);
                  },
                  "Remove a set of atoms"
                )
                .def("sort_atoms",
                  [](Molecule& m, const std::vector<int>& order) {
                    static constexpr int kAscending = 1;
                    return m.sort(order.data(), kAscending);
                  }
                )

                .def("compute_distance_matrix",
                  [](Molecule& m) {
                    return m.recompute_distance_matrix();
                  }
                )
                .def("number_fragments", static_cast<int (Molecule::*)()>(&Molecule::number_fragments), "number fragments")
                .def("fragment_membership", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::fragment_membership), "fragment number for atom")
                .def("delete_fragment", static_cast<int (Molecule::*)(int)>(&Molecule::delete_fragment), "Remove a fragment")
                .def("remove_fragment", static_cast<int (Molecule::*)(int)>(&Molecule::delete_fragment), "Remove a fragment")
                .def("atoms_in_fragment", static_cast<int (Molecule::*)(int)>(&Molecule::atoms_in_fragment), "atoms in fragment")
                .def("remove_fragment_containing_atom", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::remove_fragment_containing_atom), "Remove fragment containing an atom")
                .def("reduce_to_largest_fragment", static_cast<int (Molecule::*)()>(&Molecule::reduce_to_largest_fragment), "Strip to largest fragment")
                .def("reduce_to_largest_fragment_carefully", static_cast<int (Molecule::*)()>(&Molecule::reduce_to_largest_fragment_carefully), "Strip to largest fragment, hueristic driven")
                .def("get_fragment_membership",
                    [](Molecule& m)->std::vector<int> {
                      std::vector<int> rc(m.natoms());
                      m.fragment_membership(rc.data());
                      return rc;
                    },
                    "Fragment membership"
                )
                .def("create_components",
                  [](Molecule& m)->std::optional<std::vector<Molecule*>>{
                    if (m.number_fragments() < 2) {
                      return std::nullopt;
                    }
                    std::vector<Molecule*> res;
                    res.reserve(m.number_fragments());
                    resizable_array_p<Molecule> components;
                    if (! m.create_components(components)) {
                      return std::nullopt;
                    }
                    for (Molecule* c : components) {
                      res.push_back(c);
                    }
                    components.resize_no_delete(0);
                    return res;
                  },
                  "Split into fragments"
                )

                .def("remove_non_periodic_table_elements", static_cast<int (Molecule::*)()>(&Molecule::remove_all_non_natural_elements), "Remove non periodic table elements")
                .def("organic_only", static_cast<int (Molecule::*)()const>(&Molecule::organic_only), "True if only organic elements")

                .def("remove_explicit_hydrogens", static_cast<int (Molecule::*)()>(&Molecule::remove_explicit_hydrogens), "Remove explicit hydrogens")
                .def("remove_all", static_cast<int (Molecule::*)(atomic_number_t)>(&Molecule::remove_all), "Remove all elements with atomic number")
                //.def("remove_bonds_to_atom", static_cast<int (Molecule::*)(atomic_number_t, int)>(&Molecule::remove_bonds_to_atom), "Remove all bonds involving atom")
                .def("remove_bonds_to_atom", 
                  [](Molecule& m, atom_number_t zatom)->bool {
                    return m.remove_bonds_to_atom(zatom, 0);  // Do NOT preserve chirality
                  },
                  "Remove all bonds to an atom"
                )
                .def("remove_edge", static_cast<int (Molecule::*)(int)>(&Molecule::remove_bond), "Remove an edge by number")
                .def("remove_bond_between_atoms", static_cast<int (Molecule::*)(atomic_number_t, atomic_number_t)>(&Molecule::remove_bond_between_atoms), "Remove all bonds involving atom")
                .def("remove_all_bonds", static_cast<int (Molecule::*)()>(&Molecule::remove_all_bonds), "Remove all bonds")
                .def("chop", &Molecule::chop, "Remove the n last atoms")

                .def("implicit_hydrogens", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::implicit_hydrogens), "Implicit Hydrogens on atom")
                .def("explicit_hydrogens", static_cast<int (Molecule::*)(atom_number_t)const>(&Molecule::explicit_hydrogens), "Explicit Hydrogens on atom")
                .def("hcount", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::hcount), "Explicit and implicit Hydrogens on atom")
                .def("saturated",
                  [](Molecule& m, atom_number_t zatom)->bool{
                    return m.saturated(zatom);
                  },
                  "True if atom is fully saturated"
                )
                .def("implicit_hydrogens_known", &Molecule::implicit_hydrogens_known, "True if atom had [] in smiles")
                .def("unset_all_implicit_hydrogen_information", &Molecule::unset_all_implicit_hydrogen_information, "Discard implicit hydrogen known")
                .def("make_implicit_hydrogens_explicit", static_cast<int (Molecule::*)()>(&Molecule::make_implicit_hydrogens_explicit), "Make implicit hydrogens implicit")
                .def("AddHs",
                  [](Molecule& m) {
                    return m.make_implicit_hydrogens_explicit();
                  },
                  "implicit H become implicit"
                )
                .def("RemoveHs",
                  [](Molecule& m) {
                    return m.remove_all(1);
                  },
                  "Remove explicit H"
                )
                .def("to_scaffold",
                  [](Molecule& m) {
                    return ToScaffold(m);
                  },
                  "Convert to scaffold"
                )
                .def("change_to_graph_form",
                  [](Molecule& m) {
                    return m.change_to_graph_form();
                  },
                  "change_to_graph_form - default conditions"
                )
                .def("to_graph",
                  [](Molecule& m, const Mol2Graph& mol2graph) {
                    return m.change_to_graph_form(mol2graph);
                  },
                  "Convert to graph form"
                )

                //.def("valence_ok", static_cast<int (Molecule::*)()>(&Molecule::valence_ok), "True if all atoms ok valence")
                .def("valence_ok",
                  [](Molecule& m)->bool{
                    return m.valence_ok();
                  },
                  "True if all atoms valence ok"
                )

                .def("canonical_rank", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::canonical_rank), "Atom's canonical rank")
                .def("canonical_ranks",
                  [](Molecule& m) ->std::vector<int>{
                    std::vector<int> result(m.natoms());
                    m.canonical_ranks(result.data());
                    return result;
                  },
                  "canonical ranks"
                )
                .def("symmetry_class", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::symmetry_class), "Atom's symmetry class")
                .def("number_symmetry_classes", static_cast<int (Molecule::*)()>(&Molecule::number_symmetry_classes), "Number symmetry classes")
                .def("symmetry_equivalents",
                  [](Molecule& m, atom_number_t zatom)->std::vector<int>{
                    Set_of_Atoms tmp;
                    m.symmetry_equivalents(zatom, tmp);
                    std::vector<int> result;
                    result.reserve(tmp.size());
                    for (atom_number_t a : tmp) {
                      result.push_back(a);
                    }
                    return result;
                  },
                  "Atoms related to zatom by symmetry"
                )

                //.def("build_from_smiles", static_cast<int (Molecule::*)(const std::string&)>(&Molecule::build_from_smiles), "Build from smiles")
                .def("build_from_smiles", 
                  [](Molecule& m, const std::string& s)->bool{
                    return m.build_from_smiles(s);
                  },
                  "Build from smiles"
                )
                .def("smiles", &Molecule::Smiles, "Smiles")
                .def("unique_smiles", &Molecule::UniqueSmiles, "Unique Smiles")
                .def("random_smiles", &Molecule::RandomSmiles, "Random Smiles")
                .def("isotopically_labelled_smiles",
                  [](Molecule& m)->std::string{
                    const IWString s = m.isotopically_labelled_smiles();
                    return std::string(s.data(), s.size());
                  },
                  "Smiles with isotopes as atom number"
                )
                .def("unique_kekule_smiles",
                  [](Molecule& m)->std::string{
                    return m.UniqueKekuleSmiles().AsString();
                  },
                  "Unique Kekule form"
                )
                .def("aromatic_smiles",
                  [](Molecule& m)->std::string{
                    return m.aromatic_smiles().AsString();
                  },
                  "Non unique, aromatic smiles"
                )
                .def("smiles_atom_order",
                  [](Molecule& m)->std::vector<int>{
                    std::vector<int> result(m.natoms());
                    m.smiles_atom_order(result.data());
                    return result;
                  },
                  "atom order in most recent smiles"
                )
                .def("smiles_starting_with_atom",
                  [](Molecule& m, atom_number_t zatom) {
                    const IWString& s = m.smiles_starting_with_atom(zatom);
                    return s.AsString();
                  },
                  "smiles starting at atom"
                )

                .def("remove_hydrogens_known_flag_to_fix_valence_errors",
                     &Molecule::remove_hydrogens_known_flag_to_fix_valence_errors, "Remove problematic square brackets")

                .def("add", static_cast<int (Molecule::*)(const Molecule*)>(&Molecule::add_molecule), "add molecule")
                .def("are_bonded", static_cast<int (Molecule::*)(atom_number_t, atom_number_t)const>(&Molecule::are_bonded), "True if atoms are bonded")

                .def("formal_charge", static_cast<formal_charge_t (Molecule::*)(atom_number_t)const>(&Molecule::formal_charge), "formal charge on atom")
                .def("set_formal_charge", static_cast<void (Molecule::*)(atom_number_t, formal_charge_t)>(&Molecule::set_formal_charge), "set formal charge on atom")
                .def("has_formal_charges", static_cast<int (Molecule::*)()const>(&Molecule::has_formal_charges), "Does the molecule have atoms with formal charges")
                .def("number_formal_charges", static_cast<int (Molecule::*)()const>(&Molecule::number_formally_charged_atoms), "Number of atoms with a formal charge")
                .def("net_formal_charge", static_cast<int (Molecule::*)()const>(&Molecule::net_formal_charge), "Net formal charge")

                .def("number_chiral_centres", static_cast<int (Molecule::*)()const>(&Molecule::chiral_centres), "Number chiral centres")
                .def("remove_all_chiral_centres", static_cast<int (Molecule::*)()>(&Molecule::remove_all_chiral_centres), "Remove all chiral centres")
                .def("chiral_centre_at_atom",
                  [](const Molecule& m, atom_number_t zatom)->std::optional<Chiral_Centre*> {
                    Chiral_Centre* c = m.chiral_centre_at_atom(zatom);
                    if (c == nullptr) {
                      return std::nullopt;
                    }
                    return c;
                  },
                  "Chiral centre at atom",
                  py::return_value_policy::reference
                )
                .def("invert_chirality_on_atom", &Molecule::invert_chirality_on_atom, "Invert chirality")
                .def("remove_chiral_centre_at_atom", &Molecule::remove_chiral_centre_at_atom, "Remove specific chiral centre")
                .def("chiral_centres",
                  [](const Molecule& m)->std::vector<Chiral_Centre*>{
                    std::vector<Chiral_Centre*> result;
                    for (const Chiral_Centre* c : m.ChiralCentres()) {
                      Chiral_Centre* s = new Chiral_Centre(*c);
                      result.push_back(s);
                    }
                    return result;
                  },
                  py::return_value_policy::copy,
                  "List of chiral centres"
                )
                .def("bonds",
                  [](const Molecule& m)->std::vector<const Bond*>{
                    std::vector<const Bond*> result;
                    result.reserve(m.nedges());
                    for (const Bond * b : m.bond_list()) {
                      result.push_back(b);
                    }

                    return result;
                  },
                  py::return_value_policy::move
                )

                .def("remove_isotopes",
                  [](Molecule& m) {
                    return m.transform_to_non_isotopic_form();
                  },
                  "Remove isotopes"
                )
                .def("isotope", &Molecule::isotope, "Isotope on atom")
                .def("set_isotope", 
                  [](Molecule& m, atom_number_t zatom, isotope_t iso) {
                    return m.set_isotope(zatom,iso);
                  },
                  "Set isotope on atom"
                )
                .def("set_isotopes",
                  [](Molecule& m, const Set_of_Atoms& s, isotope_t iso) {
                    return m.set_isotope(s, iso);
                  },
                  "Set isotope for atoms in 's'"
                )
                .def("set_isotopes",
                  [](Molecule& m, const std::vector<atom_number_t>& s, isotope_t iso) {
                    return m.set_isotope(s, iso);
                  },
                  "Set isotope for atoms in 's'"
                )
                .def("number_isotopic_atoms", static_cast<int (Molecule::*)()const>(&Molecule::number_isotopic_atoms), "Number atoms with isotopes")

                .def("bonds_between", static_cast<int (Molecule::*)(atom_number_t, atom_number_t)>(&Molecule::bonds_between), "bonds between atoms")
                .def("longest_path", static_cast<int (Molecule::*)()>(&Molecule::longest_path), "longest path in molecule")
                .def("most_distant_pair",
                  [](Molecule& m)->std::pair<int, int> {
                    atom_number_t a1 = INVALID_ATOM_NUMBER;
                    atom_number_t a2 = INVALID_ATOM_NUMBER;
                    int longest_distance = 0;
                    const int matoms = m.natoms();
                    for (int i = 0; i < matoms; ++i) {
                      for (int j = i + 1; j < matoms; ++j) {
                        if (m.fragment_membership(i) != m.fragment_membership(j)) {
                          continue;
                        }
                        const int d = m.bonds_between(i, j);
                        if (d > longest_distance) {
                          a1 = i;
                          a2 = j;
                          longest_distance = d;
                        }
                      }
                    }
                    return std::make_pair(a1, a2);
                  },
                  "Most separated atoms"
                )

                .def("reset_atom_map_numbers", static_cast<void (Molecule::*)()>(&Molecule::reset_all_atom_map_numbers), "Reset atom map numbers")
                .def("set_atom_map_number", static_cast<void (Molecule::*)(atom_number_t, int)>(&Molecule::set_atom_map_number), "Set atom map number")
                .def("atom_map_number", static_cast<int (Molecule::*)(atom_number_t)const>(&Molecule::atom_map_number), "Set atom map number")
                .def("atom_with_atom_map_number", &Molecule::atom_with_atom_map_number, "Atom with atom map number")

                .def("bond_length", static_cast<distance_t (Molecule::*)(atom_number_t, atom_number_t, int)const>(&Molecule::bond_length), "Bond length between atoms")
                .def("bond_angle", static_cast<distance_t (Molecule::*)(atom_number_t, atom_number_t, atom_number_t, int)const>(&Molecule::bond_angle), "Bond length between atoms")
                .def("bond_angle",
                  [](const Molecule& m, atom_number_t a1, atom_number_t a2, atom_number_t a3)->float {
                    return m.bond_angle(a1, a2, a3, 1);
                  }
                )
                .def("dihedral_angle", static_cast<distance_t (Molecule::*)(atom_number_t, atom_number_t, atom_number_t, atom_number_t, int)const>(&Molecule::dihedral_angle), "Dihedral angle involving atoms")
                .def("signed_dihedral_angle", static_cast<distance_t (Molecule::*)(atom_number_t, atom_number_t, atom_number_t, atom_number_t)const>(&Molecule::signed_dihedral_angle), "Signed dihedral angle involving atoms")
                .def("distance_between_atoms", static_cast<distance_t (Molecule::*)(atom_number_t, atom_number_t)const>(&Molecule::distance_between_atoms), "Spatial distance between atoms")
                .def("longest_intra_molecular_distance", &Molecule::longest_intra_molecular_distance, "longest_intra_molecular_distance")
                .def("translate",
                  [](Molecule& m, distance_t x, distance_t y, distance_t z) {
                    return m.translate_atoms(x, y, z);
                  },
                  "translate coordinates"
                )

                // This did not work, the name always showed up as empty.
                // Besides, I don't think I want to enable mol.name = xxx, when everything else is via functions.
                //.def_property("name", &Molecule::Name, static_cast<void (Molecule::*)(const std::string&)>(&Molecule::set_name), "name")
                .def("name", &Molecule::Name, "Name")
                .def("set_name",
                  [](Molecule& m, const std::string& s) {
                    m.set_name(s);
                  },
                  "name"
                )

                .def("x", static_cast<coord_t (Molecule::*)(atom_number_t)const>(&Molecule::x), "x coordinate")
                .def("y", static_cast<coord_t (Molecule::*)(atom_number_t)const>(&Molecule::y), "y coordinate")
                .def("z", static_cast<coord_t (Molecule::*)(atom_number_t)const>(&Molecule::z), "z coordinate")
                .def("setx", static_cast<void (Molecule::*)(atom_number_t, coord_t)>(&Molecule::setx), "set x coordinate")
                .def("sety", static_cast<void (Molecule::*)(atom_number_t, coord_t)>(&Molecule::sety), "set y coordinate")
                .def("setz", static_cast<void (Molecule::*)(atom_number_t, coord_t)>(&Molecule::setz), "set z coordinate")
                .def("setxyz", static_cast<void (Molecule::*)(atom_number_t, coord_t, coord_t, coord_t)>(&Molecule::setxyz), "Set coordinates")
                .def("discern_chirality_from_3d_structure", &Molecule::discern_chirality_from_3d_structure, "Find chiral centres")

                .def("molecular_formula", static_cast<std::string (Molecule::*)()>(&Molecule::MolecularFormula), "Molecular formula")

                .def("partial_charge_type",
                  [](const Molecule& m)->std::string{
                    return m.partial_charge_type().AsString();
                  },
                  "Type of partial charges stored"
                )
                .def("invalidate_partial_charges",
                  [](Molecule& m) {
                    return m.invalidate_charges();
                  },
                  "Discard any partial charge information"
                )
                .def("partial_charge", &Molecule::partial_charge, "Partial charge on atom")
                .def("compute_Abraham_partial_charges", 
                  [](Molecule& m) {
                    return m.compute_Abraham_partial_charges();
                  },
                  "Abraham partial charges"
                )

                .def("compute_Gasteiger_partial_charges", 
                  [](Molecule& m) {
                    return m.compute_Gasteiger_partial_charges();
                  },
                  "Gasteiger partial charges"
                )
                .def("compute_Huckel_partial_charges",
                  [](Molecule& m) {
                    return m.compute_Huckel_partial_charges();
                  },
                  "Huckel partial charges"
                )
                .def("compute_Gasteiger_Huckel_partial_charges",
                  [](Molecule& m) {
                    return m.compute_Gasteiger_Huckel_partial_charges();
                  },
                  "Gasteiger Huckel partial charges"
                )
                //.def("compute_Del_Re_partial_charges", &Molecule::compute_Del_Re_partial_charges, "Del Re partial charges")
                //.def("compute_Pullman_partial_charges", &Molecule::compute_Pullman_partial_charges, "Pullman partial charges")

                .def("highest_coordinate_dimensionality", static_cast<int (Molecule::*)()const>(&Molecule::highest_coordinate_dimensionality), "highest coordinate dimensionality")
                .def("debug_string", static_cast<std::string (Molecule::*)()const>(&Molecule::debug_string), "Dump of internal data structures")
                .def(py::self += py::self)
                .def(py::self + py::self)
                .def("__repr__",
                  [](Molecule &m) {
                    IWString mf;
                    m.isis_like_molecular_formula(mf);
                    IWString s;
                    s << '<' << m.name() << " with " << m.natoms() << " atoms " << mf << '>';
                    return std::string(s.data(), s.length());
                   }
                 )
                .def("__str__",
                  [](Molecule& m) {
                    IWString s;
                    s << m.smiles() << ' ' << m.name();
                    return std::string(s.data(), s.length());
                  }
                )
                .def("__len__",
                  [](const Molecule& m) {
                    return m.natoms();
                  }
                )
                .def("__getitem__",
                  [](const Molecule& mol, int ndx) {
                    return mol[ndx];
                  },
                  py::return_value_policy::reference
                )
                .def("__iter__",
                  [](const Molecule&m) {
                    return py::make_iterator(m.begin(), m.end());
                  },
                  py::keep_alive<0, 1>()
                )
                .def("__copy__",
                  [](const Molecule& rhs) {
                    return Molecule(rhs);
                  }
                )
                .def("__eq__",
                  [](Molecule& m1, Molecule& m2)->bool{
                    return m1 == m2;
                  },
                  "True if molecules are identical"
                )
                .def("__contains__",
                  [](const Molecule& m, atomic_number_t z)->bool{
                    return m.natoms(z);
                  },
                  "True if molecule contains z"
                )
                .def("__contains__",
                  [](const Molecule& m, const std::string& s)->bool{
                    const_IWSubstring tmp(s);
                    const Element* e = get_element_from_symbol_no_case_conversion(tmp);
                    if (e == nullptr) {
                      throw py::value_error("Unrecognised element type");;
                    }
                    return m.natoms(e);
                  },
                  "True if molecule contains z"
                )
                .def("__contains__",
                  [](Molecule& m, Substructure_Query& q)->bool{
                    return q.substructure_search(&m);
                  },
                  "Substructure search"
                )

  ;

  py::class_<Atom, std::shared_ptr<Atom>>(m, "Atom")
    .def(py::init<int>())
    .def("atomic_number", static_cast<int (Atom::*)()const>(&Atom::atomic_number), "Atomic Number")
    .def("isotope", static_cast<isotope_t (Atom::*)()const>(&Atom::isotope), "isotope")
    .def("ncon", static_cast<int (Atom::*)()const>(&Atom::ncon), "Number of connections")
    .def("nbonds", static_cast<int (Atom::*)()const>(&Atom::nbonds), "Number of bonds - single=1, double=2...")
    .def("formal_charge", &Atom::formal_charge, "formal charge")
    .def("atomic_weight", &Atom::atomic_weight, "atomic weight")
    .def("exact_mass", &Atom::exact_mass, "exact mass")
    .def("is_bonded_to",
      [](const Atom& a, atom_number_t atom)->bool{
          return a.is_bonded_to(atom);
      },
      "True if atom is bonded to other atom"
    )
    .def("valence_ok", &Atom::valence_ok, "True if valence is ok")
    .def("fully_saturated", &Atom::fully_saturated, "True if fully saturated")
    .def("other", static_cast<atom_number_t (Atom::*)(atom_number_t, int)const>(&Atom::other), "Other connection")
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
    .def("connections",
      [](const Atom& a, atom_number_t atom_number)->std::vector<int>{
        std::vector<int> result;
        for (const Bond * b : a) {
          result.push_back(b->other(atom_number));
        }
        return result;
      },
      "List of connected atoms"
    )
    .def("__contains__",
      [](const Atom& a, atom_number_t atom_number)->bool{
        return a.is_bonded_to(atom_number);
      },
      "True if atom_number is bonded"
    )
    .def("__LEN__",
      [](const Atom& a) {
        return a.ncon();
      },
      "Number of connections"
    )


    // This cannot be implemented because the atom does not know its atom number.
    //.def("GetNeighbors",
    //  [](const Atom& s) {
    //  }
    //)
  ;

  py::class_<Bond, std::shared_ptr<Bond>>(m, "Bond")
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
    .def("is_directional", &Bond::is_directional, "True of bond is directional")
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
    .def("is_aromatic_bond",
      [](const Bond& b)->bool{
        return b.is_aromatic();
      },
      "True if aromatic"
    )
    .def("is_aromatic",
      [](const Bond& b)->bool{
        return b.is_aromatic();
      },
      "True if aromatic"
    )
    .def("btype",
      [](const Bond& b)->BondType{
        if (b.is_aromatic()) {
          return BondType::kAromaticBond;
        }
        if (b.is_single_bond()) {
          return BondType::kSingleBond;
        } 
        if (b.is_double_bond()) {
          return BondType::kDoubleBond;
        }
        if (b.is_triple_bond()) {
          return BondType::kTripleBond;
        }
        throw py::value_error("Unrecognised bond type");;
      },
      "btype"
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
    .def("IsInRing",
      [](const Bond& b)->bool{
        int nr;
        if (b.nrings(nr)) {
          return nr;
        }
        throw py::value_error("Bond.nrings:ring membership not computed");
      },
      "True if bond in ring"
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
    .def("GetBeginAtomIdx",
      [](const Bond& b)->int{
        return b.a1();
      },
      "First atom"
    )
    .def("GetEndAtomIdx",
      [](const Bond& b)->int{
        return b.a2();
      },
      "Second atom"
    )
    .def("GetBondType",
      [](const Bond& b)->BondType{
        if (b.is_single_bond()) {
          return kSingleBond;
        }
        if (b.is_double_bond()) {
          return kDoubleBond;
        }
        if (b.is_triple_bond()) {
          return kTripleBond;
        }
        return kUnknown;
      },
      "Bond type"
    )
    .def("__contains__",
      [](const Bond& b, atom_number_t atom)->bool{
        return b.involves(atom);
      },
      "True if atom part of bond"
    )
  ;

  py::class_<std::vector<int>>(m, "IntVector")
    .def(py::init<>())
    .def(py::init([](const std::vector<int>& rhs) {
      return std::make_unique<std::vector<int>>(rhs);
    }))
    .def("clear", &std::vector<int>::clear)
    .def("pop_back", &std::vector<int>::pop_back)
    .def("__len__", [](const std::vector<int> &v) { return v.size(); })
    .def("__iter__", [](std::vector<int> &v) {
       return py::make_iterator(v.begin(), v.end());
    }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
  ;

  py::class_<Set_of_Atoms>(m, "Set_of_Atoms")
    .def(py::init<>())
    .def(py::init([](const std::vector<int>& s) {
        return std::unique_ptr<Set_of_Atoms>(new Set_of_Atoms(s));
    }))
    .def("empty",
      [](const Set_of_Atoms& s)->bool{
        return s.empty();
      },
      "True if empty"
    )
    .def("size", &Set_of_Atoms::size, "size")
#ifdef LILLYMOL_VECTOR_OPAQUE
    .def("scatter",
      [](const Set_of_Atoms& s, std::vector<int>& v, int value) {
        for (auto a : s) {
          v[a] = value;
        }
      }
    )
    .def("increment_vector",
      [](const Set_of_Atoms& s, std::vector<int>& v, int value) {
        for (auto a: s) {
          v[a] += value;
        }
      }
    )
#endif
    .def("contains_both",
      [](const Set_of_Atoms& s, atom_number_t a1, atom_number_t a2)->bool{
        return s.contains_atoms(a1, a2);
      },
      "True if contains both a1 and a2"
    )
    .def("__len__",
      [](const Set_of_Atoms& s) {
        return s.size();
      },
      "Size"
    )
    .def("__contains__",
      [](const Set_of_Atoms& s, atom_number_t atom)->bool{
        return s.contains(atom);
      },
      "Is item included"
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
    .def("__repr__",
      [](const Set_of_Atoms& s) {
        IWString result;
        result << "<Set_of_Atoms";
        for (atom_number_t a : s) {
          result << ' ' << a;
        }
        result << '>';
        return result.AsString();
      }
    )
    .def("__str__",
      [](const Set_of_Atoms& s) {
        IWString result;
        for (atom_number_t a : s) {
          result << ' ' << a;
        }
        return result.AsString();
      }
    )
    // This does not appear to be having any effect.
    .def("__eq__",
      [](const Set_of_Atoms& lhs, const std::vector<int>& rhs)->bool {
        const uint32_t n = lhs.size();
        //std::cerr << "Checking size " << n << " and " << rhs.size() << '\n';
        if (n != rhs.size()) {
          return false;
        }
        for (uint32_t i = 0; i < n; ++i) {
          if (lhs[i] != rhs[i]) {
            return false;
          }
        }
        return true;
      }
    )
    .def("append",
      [](Set_of_Atoms& s, atom_number_t extra) {
        s.add(extra);
      }
    )
    .def("extend",
      [](Set_of_Atoms& s, std::vector<int>& extra) {
        for (int e : extra) {
          s << e;
        }
      }
    )
  ;

  py::class_<Ring, Set_of_Atoms>(m, "Ring")
    .def(py::init<>())
    .def("size", &Ring::size, "size")
    .def("ring_number", &Ring::ring_number, "Unique ring number")
    .def("fragment_membership", &Ring::fragment_membership, "Fragment membership")
    .def("fused_system_identifier", &Ring::fused_system_identifier, "Fused system identifier")
    .def("is_fused",
      [](const Ring& r)->bool{
        return r.is_fused();
      },
      "True if the ring is fused"
    )
    .def("is_fused_to",
      [](const Ring& r1, const Ring* r2)->bool{
        return r1.is_fused_to(r2);
      },
      "True if rings are fused"
    )
    .def("fused_ring_neighbours", &Ring::fused_ring_neighbours, "Number of fused ring neighbours")
    .def("largest_number_of_bonds_shared_with_another_ring", &Ring::largest_number_of_bonds_shared_with_another_ring, "largest_number_of_bonds_shared_with_another_ring")
    .def("strongly_fused_ring_neighbours", &Ring::strongly_fused_ring_neighbours, "strongly_fused_ring_neighbours")
    .def("contains_bond",
      [](const Ring& r, atom_number_t a1, atom_number_t a2)->bool{
        return r.contains_bond(a1, a2);
      },
      "True if a1 and a2 are adjacent"
    )
    .def("contains_both",
      [](const Ring& r, atom_number_t a1, atom_number_t a2)->bool{
        return r.contains_both(a1, a2);
      },
      "True if ring contains both a1 and a2"
    )
#ifdef DELIBERATELY_NOT_IMPLEMENTED
    // THis does not work. The problem is that the Ring object has pointers to
    // other Ring objects, and the copy constructor does not really work.
    // The fused_neighbour array is empty.
    // If anyone ever wants this, we would need to convert back to ring numbers,
    // or do something at the Molecule level.
    .def("fused_neighbours",
      [](const Ring& r)->std::vector<const Ring*>{
        std::vector<const Ring*> result;
        std::cerr << "Ring has " << r.fused_ring_neighbours() << " fused ring nbrs...\n";
        for (int i = 0; i < r.fused_ring_neighbours(); ++i) {
          result.push_back(r.fused_neighbour(i));
        }
        return result;
      },
      py::return_value_policy::move,
      "List of fused rings"
    )
#endif
    .def("is_aromatic",
      [](const Ring& r)->bool{
        return r.is_aromatic();
      },
      "True if the ring is aromatic"
    )
    .def("__getitem__",
      [](const Ring& me, int ndx)->atom_number_t{
        return me.item(ndx);
      })
    .def("__repr__",
      [](const Ring &r) {
        return RingToString(r);
      }
    )
    .def("__str__",
      [](const Ring& r) {
        return RingToString(r);
      }
    )
    .def("__iter__",
      [](const Ring&r) {
        return py::make_iterator(r.begin(), r.end());
      },
      py::keep_alive<0, 1>()
    )
    .def("__len__",
      [](const Ring& r) {
        return r.size();
      },
      "Size"
    )
    .def("__contains__",
      [](const Ring& r, atom_number_t atom)->bool{
        return r.contains(atom);
      },
      "Is item included"
    )
  ;

  py::class_<Element_Transformations>(m, "ElementTransformations")
    .def(py::init<>())
    .def("active",
      [](const Element_Transformations& etrans)->bool{
        return etrans.active();
      },
      "True if active"
    )
    .def("add",
      [](Element_Transformations& etrans, const std::string& directive)->bool{
        const IWString s(directive);
        return etrans.Add(s);
      },
      "Add a transformation 'I=Cl'"
    )
    .def("process",
      [](Element_Transformations& etrans, Molecule& m) {
        return etrans.process(m);
      },
      "Apply transformations"
    )
  ;

  m.def("set_copy_name_in_molecule_copy_constructor", &set_copy_name_in_molecule_copy_constructor, "Copy name in constructor");
  m.def("MolFromSmiles", &MolFromSmiles, "Molecule from smiles");
  m.def("LillyMolFromSmiles", &MolFromSmiles, "Molecule from smiles");
  m.def("set_auto_create_new_elements", &set_auto_create_new_elements, "auto create new elements");
  m.def("set_atomic_symbols_can_have_arbitrary_length", &set_atomic_symbols_can_have_arbitrary_length, "any string is an element");
  m.def("interpret_D_as_deuterium", &element::interpret_d_as_deuterium, "D means '[2H]'");
  m.def("interpret_T_as_deuterium", &element::interpret_t_as_tritium, "T means '[3H]'");
  m.def("set_display_strange_chemistry_messages", &set_display_strange_chemistry_messages, "turn off messages about bad valences");
  m.def("set_auto_create_new_elements", &set_auto_create_new_elements, "Allow arbitrary two letter elements");
  m.def("set_atomic_symbols_can_have_arbitrary_length", &set_atomic_symbols_can_have_arbitrary_length, "Enable elements like 'Ala', 'Gly'");
  m.def("set_display_smiles_interpretation_error_messages", &set_display_smiles_interpretation_error_messages, "Set smiles error messages");
  m.def("count_atoms_in_smiles",
    [](const std::string& smiles) {
      const const_IWSubstring tmp(smiles);
      return count_atoms_in_smiles(tmp);
    }
  );

  py::enum_<BondType>(m, "BondType")
    .value("SINGLE_BOND", kSingleBond)
    .value("DOUBLE_BOND", kDoubleBond)
    .value("TRIPLE_BOND", kTripleBond)
    .value("AROMATIC_BOND", kAromaticBond)
    .export_values();
  ;

  m.def("is_chiral_implicit_hydrogen",
    [](int c)->bool {
      return IsChiralImplicitHydrogen(c);
    },
    "True if chiral connection is an implicit hydrogen"
  );

  m.def("is_chiral_lone_pair",
    [](int c)->bool {
      return IsChiralLonePair(c);
    },
    "True if chiral connection is a lone pair"
  );

}
