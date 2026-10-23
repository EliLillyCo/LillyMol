#include <algorithm>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
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

#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/chiral_centre.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/hybridization.h"
#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/mol2graph.pb.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rotbond_common.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/substructure.h"

#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/xlogp.h"

#include "pybind/molecule.h"

namespace py = pybind11;

// Lifted from test_numpy_dtypes.cpp
template <typename T>
py::array mkarray_via_buffer(size_t n) {
    return py::array(py::buffer_info(
        nullptr, sizeof(T), py::format_descriptor<T>::format(), 1, {n}, {sizeof(T)}));
}


static bool
LooksLikeSdfTagRecord(const std::string& line) {
  const size_t open = line.find('<');
  const size_t close = line.rfind('>');
  return line.size() >= 4 && line[0] == '>' && open != std::string::npos &&
         close != std::string::npos && open < close;
}

static std::string
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

static py::dict
SdfTags(const Molecule& m) {
  py::dict result;

  std::string current_tag;
  std::string current_value;

  auto flush_current_tag = [&]() {
    if (current_tag.empty()) {
      return;
    }
    result[py::str(current_tag)] = py::str(current_value);
    current_tag.clear();
    current_value.clear();
  };

  for (int i = 0; i < m.number_records_text_info(); ++i) {
    const std::string line = m.text_info(i).AsString();
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

    if (! current_value.empty()) {
      current_value.append("\n");
    }
    current_value.append(line);
  }

  flush_current_tag();
  return result;
}

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
  } else if (a == kChiralConnectionIsImplicitHydrogen) {
    result << 'H';
  } else if (a == kChiralConnectionIsLonePair) {
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

  // Add back in singly connected =O and =N
  // Note that we do this for all scaffold atoms,
  // aromatic rings, aliphatic rings and in linker groups.
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m[i];
    if (a.ncon() != 1) {
      continue;
    }
    const Bond* b = a[0];
    if (! b->is_double_bond()) {
      continue;
    }
    atom_number_t o = b->other(i);
    if (spinach[o] == 0) {
      spinach[i] = 0;
    }
  }

  return m.remove_atoms(spinach.get(), 1);
}

static Molecule
Scaffold(const Molecule& m) {
  Molecule result(m);
  ToScaffold(result);
  return result;
}

static std::vector<Set_of_Atoms>
AtomsByRadius(const Molecule& m, const Set_of_Atoms& starting_atoms,
              int max_radius) {
  if (max_radius < 0) {
    throw py::value_error("max_radius must be non-negative");
  }

  std::vector<Set_of_Atoms> result(max_radius + 1);
  const int matoms = m.natoms();
  std::vector<int> distance(matoms, -1);
  std::queue<atom_number_t> to_process;

  for (atom_number_t atom : starting_atoms) {
    if (atom < 0 || atom >= matoms) {
      throw py::value_error("starting atom is outside the molecule");
    }
    if (distance[atom] >= 0) {
      continue;
    }

    distance[atom] = 0;
    result[0].add(atom);
    to_process.push(atom);
  }

  while (!to_process.empty()) {
    const atom_number_t atom = to_process.front();
    to_process.pop();

    const int next_radius = distance[atom] + 1;
    if (next_radius > max_radius) {
      continue;
    }

    for (const Bond* bond : m[atom]) {
      const atom_number_t other = bond->other(atom);
      if (distance[other] >= 0) {
        continue;
      }

      distance[other] = next_radius;
      result[next_radius].add(other);
      to_process.push(other);
    }
  }

  return result;
}

// Ultimately we want kSingleBond to be used as the LillyMol name for
// a single bond. But having these here might make that problematic.
// TODO:ianwatson investigate.
enum BondType {
  kUnknown = 0,
  kSingleBond = 1,
  kDoubleBond = 2,
  kTripleBond = 3,
  kAromaticBond = 4
};

enum ChiralType {
  kChiUnspecified = 0,
  kChiTetrahedralCw = 1,
  kChiTetrahedralCcw = 2,
  kChiOther = 3
};

// BondType is not bond_type_t, so a cast will not do. Both Bond.btype() and
// Molecule.bond_between_atoms report a bond type to python and they must not be
// allowed to disagree, so the mapping lives here.
// Aromatic is tested first because aromatic rings are stored in a Kekule form -
// an aromatic bond is also a single or a double bond, and reporting that would
// be wrong.
static BondType
ToBondType(const Bond& b) {
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
  // kUnknown is deliberately not exported to python, so returning it would hand
  // back an unnamed enum value. An exception is more useful.
  throw py::value_error("Unrecognised bond type");
}

static std::optional<int>
IndexOf(const int* values, int needle) {
  for (int i = 0; i < 4; ++i) {
    if (values[i] == needle) {
      return i;
    }
  }

  return std::nullopt;
}

static int
InversionCount(const std::vector<int>& values) {
  int result = 0;
  for (int i = 0; i < static_cast<int>(values.size()); ++i) {
    for (int j = i + 1; j < static_cast<int>(values.size()); ++j) {
      if (values[i] > values[j]) {
        ++result;
      }
    }
  }

  return result;
}

static int
ChiralConnectionSortKey(int member, int matoms) {
  if (member >= 0) {
    return member;
  }
  if (member == kChiralConnectionIsImplicitHydrogen) {
    return matoms;
  }
  if (member == kChiralConnectionIsLonePair) {
    return matoms + 1;
  }

  return matoms + 2;
}

// Return a LillyMol-defined tetrahedral chiral tag for `zatom`. The tag is
// relative to atom-number order, with implicit Hydrogen after explicit atoms.
static std::optional<ChiralType>
TetrahedralChirality(Molecule& m, atom_number_t zatom, bool check_is_chiral) {
  if (zatom < 0 || zatom >= m.natoms()) {
    throw py::value_error("atom number outside molecule");
  }

  if (check_is_chiral && ! ::is_actually_chiral(m, zatom)) {
    return std::nullopt;
  }

  const Chiral_Centre* c = m.chiral_centre_at_atom(zatom);
  if (c == nullptr) {
    return std::nullopt;
  }

  if (! c->chirality_known() || ! c->complete()) {
    return ChiralType::kChiUnspecified;
  }

  if (c->lone_pair_count() > 0) {
    return ChiralType::kChiOther;
  }

  const int stored[4] = {
    c->top_front(),
    c->top_back(),
    c->left_down(),
    c->right_down()
  };

  std::vector<int> ordered_members(std::begin(stored), std::end(stored));
  const int matoms = m.natoms();
  std::sort(ordered_members.begin(), ordered_members.end(),
    [matoms](int c1, int c2) {
      return ChiralConnectionSortKey(c1, matoms) < ChiralConnectionSortKey(c2, matoms);
    });

  std::vector<int> permutation;
  permutation.reserve(4);
  for (int member : ordered_members) {
    std::optional<int> index = IndexOf(stored, member);
    if (! index) {
      return ChiralType::kChiUnspecified;
    }
    permutation.push_back(*index);
  }

  if (InversionCount(permutation) % 2 == 0) {
    return ChiralType::kChiTetrahedralCcw;
  }

  return ChiralType::kChiTetrahedralCw;
}

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
    .def("turn_on_most_useful_options", &Mol2Graph::TurnOnMostUsefulOptions, "turn on the options you probably want")
    .def("set_active", &Mol2Graph::set_active, "Set active")
    .def("active", &Mol2Graph::active, "True if active")
  ;

  // By adding shared_ptr we avoid copies between C++ and python.
  // It also enables lists of molecules to be passed as std::vector<Molecule*> which is mutable.
  // std::vector<Molecule> will create an array of copies.
         py::class_<Molecule, std::shared_ptr<Molecule>>(m, "Molecule")
		.def(py::init<>())
                .def(py::init([](const Molecule& rhs) {
                  return std::shared_ptr<Molecule>(new Molecule(rhs));
                }))
                .def("ok",
                  [](const Molecule& m)->bool{
                    return m.ok();
                  },
                  "Returns true if internal datastructures ok"
                )
                .def("natoms", static_cast<int (Molecule::*)()const>(&Molecule::natoms), "Number explicit atoms in molecule")
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
                .def("resize", &Molecule::resize, "Change number of atoms - dangerous")
                .def("nedges", static_cast<int (Molecule::*)()const>(&Molecule::nedges), "Number edges in molecule")
                .def("add_atom",
                  [](Molecule& m, int atnum)->atom_number_t {
                    const Element* e = get_element_from_atomic_number(atnum);
                    if (e == nullptr) {
                      std::cerr << "Invalid atomic number " << atnum << '\n';
                      return false;
                    }
                    m.add(e);
                    return m.natoms() - 1;
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
                .def("in_ring_of_size", &Molecule::in_ring_of_given_size, "True if atom in ring of give size")
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
                  "number of rings containing atom"
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
                // A Molecule computes derived properties only when asked, and a
                // Bond holds no back pointer to its Molecule, so until this has
                // run Bond.nrings() reports 0 and Bond.IsInRing() False for a ring
                // bond. This is the force call, for callers that want the bonds to
                // know and have no use for the per atom counts that
                // get_ring_membership allocates and returns.
                .def("ring_membership",
                  [](Molecule& m)->void {
                    m.ring_membership();
                  },
                  "Force ring membership perception, so the bonds know. Returns nothing - get_ring_membership returns the per atom counts"
                )
                .def("fused_system_identifier", &Molecule::fused_system_identifier, "Fused system identifier")
                .def("fused_system_size", &Molecule::fused_system_size, "Fused system size")
                .def("number_ring_systems", static_cast<int (Molecule::*)()>(&Molecule::number_ring_systems), "Number ring systems")
                .def("ring", static_cast<const Ring* (Molecule::*)(int)>(&Molecule::ringi), py::return_value_policy::reference, "I'th ring")
                .def("rings",
                  [](Molecule& m)->std::vector<std::unique_ptr<Ring>>{
                    std::vector<std::unique_ptr<Ring>> result;
                    // Was not able to make this work passing a vector of the actual rings.
                    // Works once or twice then crashes. Does this create a memory leak?
                    for (const Ring* r : m.sssr_rings()) {
                      result.push_back(std::make_unique<Ring>(*r));
                    }
                    return result;
                  },
                  "iterable list of rings in the molecule"
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
                .def("label_atoms_by_atom_number", [](Molecule& m) {
                  const int matoms = m.natoms();
                  for (int i = 1; i < matoms; ++i) {
                    m.set_isotope(i, i);
                  }
                },
                "Isotope becomes atom number"
                )
                .def("label_atoms_by_ring_system",
                  [](Molecule& m)->std::vector<int>{
                    std::vector<int> rc(m.natoms());
                    m.label_atoms_by_ring_system(rc.data());
                    return rc;
                  },
                  "For each atom the ring system identifier"
                )
                .def("label_atoms_by_ring_system_including_spiro_fused",
                  [](Molecule& m)->std::vector<int>{
                    std::vector<int> rc(m.natoms());
                    m.label_atoms_by_ring_system_including_spiro_fused(rc.data());
                    return rc;
                  },
                  "For each atom the ring system identifier"
                )
                .def("label_atoms_by_ring_system_including_spiro_fused_np",
                  [](Molecule& m)-> py::array_t<int> {
                    py::array_t<int> result = mkarray_via_buffer<int>(m.natoms());
                    auto req = result.request();
                    int* ptr = static_cast<int*>(req.ptr);
                    m.label_atoms_by_ring_system_including_spiro_fused(ptr);
                    return result;
                  },
                  "For each atom the ring system identifier"
                )
                .def("amw",
                  [](Molecule& m)->float {
                    if (! m.ContainsIsotopicAtoms()) [[likely]] {
                      return m.molecular_weight();
                    }

                    if (PyErr_WarnEx(PyExc_UserWarning,
                          "amw() erases isotopic labels before weighing, the mass and the "
                          "Hydrogen count both, so [37C]OC weighs the same as COC. LillyMol "
                          "has no table of isotopic masses, so a genuine isotopic molecular "
                          "weight is not available - deuterium is weighed as Hydrogen.", 1) < 0) {
                      throw py::error_already_set();
                    }

                    return lillymol::MolecularWeightIsotopesAsLabels(m);
                  },
R"(Average molecular weight.

Isotopic labels are erased before weighing, both the mass and the Hydrogen
count, so [37C]OC weighs the same as COC. In LillyMol an isotope is usually an
arbitrary atom marker rather than a statement about the nucleus, and a marker
should not change the weight. A strict smiles reading says [37C] has no
Hydrogens; that is almost never what was meant. This matches the amw feature of
molecule_filter.

LillyMol has no table of isotopic masses, so a genuine isotopic molecular weight
is not available from here at all - deuterium is weighed as Hydrogen. If you
need one, use a toolkit that has the table.

A warning is issued the first time an isotope is encountered. Silence it with
warnings.filterwarnings, or turn it into an error with -W error::UserWarning.

See also amw_ignore_isotopes, which differs on one point.)")
                .def("amw_ignore_isotopes",
                     static_cast<float (Molecule::*)()const>(&Molecule::molecular_weight_ignore_isotopes),
R"(Average molecular weight, with an isotopic atom counted at the normal weight of
its element.

Differs from amw in that the declared Hydrogen count is honoured. A bracket atom
states its own Hydrogen count, and [37C] states zero, so [37C]OC weighs 43.045
here where amw gives 46.068. Use this only if you want the strict smiles reading;
amw is what molecule_filter agrees with.

Never refuses, and issues no warning.)")
                .def("exact_mass", static_cast<exact_mass_t (Molecule::*)()const>(&Molecule::exact_mass), "Exact Mass")
                .def("ncon", static_cast<int (Molecule::*)(atom_number_t)const>(&Molecule::ncon), "Connections to Atom")
                .def("connections",
                  [](const Molecule& m, atom_number_t zatom)->std::vector<int>{
                    const Set_of_Atoms s = m.connections(zatom);
                    return std::vector<int>(s.rawdata(), s.rawdata() + s.size());
                  },
                  "Atoms connected to atom"
                )
                .def("other_atom",
                  [](const Molecule& m, atom_number_t zatom, int n)->atom_number_t {
                    return m.other(zatom, n);
                  },
                  "Fetch the atom number of the n'th connection to atom"
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
                  "Pi electrons on atom"
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
                  "smarts for atom. If aromaticity has been computed, will include aromaticity"
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
                .def("remove_atoms",
                  [](Molecule& m, py::array_t<int> to_remove, int flag)->int {
                    //auto req = to_remove.request();
                    int* ptr = static_cast<int*>(to_remove.request().ptr);
                    int rc = 0;
                    for (int i = m.natoms() - 1; i >= 0; --i) {
                      if (ptr[i] == flag) {
                        m.remove_atom(i);
                        ++rc;
                      }
                    }
                    return rc;
                  },
                  "Remove atoms where to_remove[i] == flag"
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
                .def("fragment_membership", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::fragment_membership), "For each atom the fragment number containing that atom")
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
                    "For each atom the fragment membership"
                )
                .def("create_components",
                      [](Molecule& m)->std::optional<std::vector<std::shared_ptr<Molecule>>>{
                    if (m.number_fragments() < 2) {
                      return std::nullopt;
                    }

                    resizable_array_p<Molecule> components;
                    if (! m.create_components(components)) {
                      return std::nullopt;
                    }

                    std::vector<std::shared_ptr<Molecule>> res;
                    res.reserve(components.number_elements());

                    for (Molecule* c : components) {
                      res.emplace_back(c);  // shared_ptr takes ownership
                    }
                    components.resize_no_delete(0);    // array no longer owns them
                    return res;
                  },
                  "Split a multi fragment molecule into fragment molecules"
                )


                .def("remove_non_periodic_table_elements", static_cast<int (Molecule::*)()>(&Molecule::remove_all_non_natural_elements), "Remove non periodic table elements")
                .def("organic_only", static_cast<int (Molecule::*)()const>(&Molecule::organic_only), "True if only organic elements")
                .def("non_organic_atom_count", &Molecule::non_organic_atom_count, "Number of non organic atoms")
                .def("is_organic",
                  [](const Molecule& m, atom_number_t zatom)->bool {
                    return m.is_organic(zatom);
                  },
                  "True if the atom is organic"
                )
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
                .def("lipinski_num_h_donors", &Molecule::LipinskiNumHDonors, "Lipinski hydrogen bond donor count")
                .def("lipinski_num_h_acceptors", &Molecule::LipinskiNumHAcceptors, "Lipinski hydrogen bond acceptor count")
                .def("rdkit_num_h_donors", &Molecule::RDKitNumHDonors, "RDKit compatible hydrogen bond donor count, NOT a Lipinski count")
                .def("rdkit_num_h_acceptors", &Molecule::RDKitNumHAcceptors, "RDKit compatible hydrogen bond acceptor count, NOT a Lipinski count")
                .def("saturated",
                  [](Molecule& m, atom_number_t zatom)->bool{
                    return m.saturated(zatom);
                  },
                  "True if atom is fully saturated"
                )
                .def("unsaturation", &Molecule::unsaturation)
                .def("hybridization",
                  [](Molecule& m, atom_number_t zatom) -> Hybridization {
                    if (! m.ok_atom_number(zatom)) {
                      throw py::value_error("hybridization atom number outside [0, natoms)");
                    }
                    return HybridizationState(m, zatom);
                  },
                  "RDKit-like hybridization of atom, computed on demand"
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
                  "Convert to scaffold in place"
                )
                .def("scaffold",
                  [](const Molecule& m) {
                    return Scaffold(m);
                  },
                  "Return a new Molecule containing the scaffold"
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
                  "Convert to graph form under control of mol2graph"
                )

                //.def("valence_ok", static_cast<int (Molecule::*)()>(&Molecule::valence_ok), "True if all atoms ok valence")
                .def("valence_ok",
                  [](Molecule& m)->bool{
                    return m.valence_ok();
                  },
                  "True if all atoms valence ok"
                )
                .def("valence_ok",
                  [](Molecule& m, atom_number_t a)->bool {
                    return m.valence_ok(a);
                  },
                  "Returns true of atom a has an ok valence"
                )

                .def("canonical_rank", static_cast<int (Molecule::*)(atom_number_t)>(&Molecule::canonical_rank), "Atom's canonical rank")
                .def("canonical_ranks",
                  [](Molecule& m) ->std::vector<int>{
                    std::vector<int> result(m.natoms());
                    m.canonical_ranks(result.data());
                    return result;
                  },
                  "For each atom the canonical rank"
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
                .def("renumber_atoms",
                  [](Molecule& m, const std::vector<int>& new_number)->int {
                    const int matoms = m.natoms();
                    if (static_cast<int>(new_number.size()) != matoms) {
                      throw py::value_error("renumber_atoms requires one entry for each atom");
                    }

                    std::vector<int> seen(matoms, 0);
                    for (int i = 0; i < matoms; ++i) {
                      const int destination = new_number[i];
                      if (destination < 0 || destination >= matoms) {
                        throw py::value_error("renumber_atoms mapping contains an atom number outside [0, natoms)");
                      }
                      if (seen[destination]) {
                        throw py::value_error("renumber_atoms mapping contains duplicate atom numbers");
                      }
                      seen[destination] = 1;
                    }

                    return m.renumber_atoms(new_number.data());
                  },
                  "Renumber atoms. new_number[i] is the new atom number for current atom i"
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

                .def("add", static_cast<int (Molecule::*)(const Molecule*)>(&Molecule::add_molecule), "add molecule, no bonds formed, generates multi-fragment molecule")
                .def("are_bonded", static_cast<int (Molecule::*)(atom_number_t, atom_number_t)const>(&Molecule::are_bonded), "True if atoms are bonded")

                // Both of these use bond_between_atoms_if_present. The plain
                // bond_between_atoms asserts when the atoms are not bonded, and
                // aborting the interpreter is not a reasonable answer to a question.
                //
                // Both perceive aromaticity first. Without it a benzene ring bond
                // reports is_aromatic() False and nrings() 0 - plausible wrong
                // answers rather than errors. Molecule::nrings() does not help, it
                // does not push ring membership down onto the bonds. Once perceived
                // the call is a pointer test.
                .def("bond_between_atoms",
                  [](Molecule& m, atom_number_t a1, atom_number_t a2)->const Bond* {
                    m.compute_aromaticity_if_needed();
                    return m.bond_between_atoms_if_present(a1, a2);
                  },
                  // reference is essential, not decoration. Bond is held by
                  // shared_ptr, so the default policy for a returned pointer,
                  // take_ownership, would have python delete a bond the Molecule
                  // owns. keep_alive stops the molecule being collected while python
                  // still holds one of its bonds.
                  py::return_value_policy::reference,
                  py::keep_alive<0, 1>(),
                  "The Bond between two atoms, None if they are not bonded"
                )
                // Preferred when the type is all that is wanted, which is the usual
                // case. Constructing the python wrapper for a Bond costs far more
                // than the lookup itself, and this avoids it.
                .def("bond_type_between_atoms",
                  [](Molecule& m, atom_number_t a1, atom_number_t a2)->std::optional<BondType> {
                    m.compute_aromaticity_if_needed();
                    const Bond* b = m.bond_between_atoms_if_present(a1, a2);
                    if (b == nullptr) {
                      return std::nullopt;
                    }
                    return ToBondType(*b);
                  },
                  "BondType between two atoms, None if not bonded. Aromatic bonds report AROMATIC_BOND, not their Kekule type"
                )

                .def("formal_charge", static_cast<formal_charge_t (Molecule::*)(atom_number_t)const>(&Molecule::formal_charge), "formal charge on atom")
                .def("set_formal_charge", static_cast<void (Molecule::*)(atom_number_t, formal_charge_t)>(&Molecule::set_formal_charge), "set formal charge on atom")
                .def("has_formal_charges", static_cast<int (Molecule::*)()const>(&Molecule::has_formal_charges), "Does the molecule have atoms with formal charges")
                .def("number_formal_charges", static_cast<int (Molecule::*)()const>(&Molecule::number_formally_charged_atoms), "Number of atoms with a formal charge")
                .def("net_formal_charge", static_cast<int (Molecule::*)()const>(&Molecule::net_formal_charge), "Net formal charge")

                .def("number_chiral_centres", static_cast<int (Molecule::*)()const>(&Molecule::chiral_centres), "Number chiral centres")
                .def("remove_all_chiral_centres", static_cast<int (Molecule::*)()>(&Molecule::remove_all_chiral_centres), "Remove all chiral centres")
                .def("chiral_centre_at_atom",
                  [](const Molecule& m, atom_number_t zatom)->std::optional<Chiral_Centre> {
                    Chiral_Centre* c = m.chiral_centre_at_atom(zatom);
                    if (c == nullptr) {
                      return std::nullopt;
                    }
                    return *c;
                  },
                  "Chiral centre at atom"
//                py::return_value_policy::reference
                )
                .def("invert_chirality_on_atom", &Molecule::invert_chirality_on_atom, "Invert chirality")
                .def("remove_chiral_centre_at_atom", &Molecule::remove_chiral_centre_at_atom, "Remove specific chiral centre")
                .def("chiral_centres",
                  [](const Molecule& m)->std::vector<Chiral_Centre>{
                    std::vector<Chiral_Centre> result;
                    for (const Chiral_Centre* c : m.ChiralCentres()) {
                      result.push_back(*c);
                    }
                    return result;
                  },
//                py::return_value_policy::copy,
                  "List of chiral centres"
                )
                .def("revert_all_directional_bonds_to_non_directional", &Molecule::revert_all_directional_bonds_to_non_directional, "remove bond directionality (cis/trans)")
                .def("bonds",
                  [](const Molecule& m)->std::vector<const Bond*>{
                    std::vector<const Bond*> result;
                    result.reserve(m.nedges());
                    for (const Bond * b : m.bond_list()) {
                      result.push_back(b);
                    }

                    return result;
                  },
                  py::return_value_policy::reference
                )
                .def("bond",
                  [](const Molecule& m, int ndx)->const Bond* {
                    return m.bondi(ndx);
                  },
                  py::return_value_policy::reference,
                  "Return the i'th Bond"
                )

                .def("remove_isotopes",
                  [](Molecule& m) {
                    return m.transform_to_non_isotopic_form();
                  },
                  "Convert isotopes to non isotopic"
                )
                .def("transform_to_non_isotopic_form", &Molecule::transform_to_non_isotopic_form, "Convert isotopes to non isotopic")
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
                .def("set_isotopes",
                  [](Molecule& m, py::array_t<int> iso) {
                    int* ptr = static_cast<int*>(iso.request().ptr);
                    const int matoms = m.natoms();
                    for (int i = 0; i < matoms; ++i) {
                      m.set_isotope(i, ptr[i]);
                    }
                  },
                  "Set isotope for each atom"
                )
                .def("number_isotopic_atoms", static_cast<int (Molecule::*)()const>(&Molecule::number_isotopic_atoms), "Number atoms with isotopes")
                .def("first_atom_with_isotope",
                  [](const Molecule& m, isotope_t iso) -> atom_number_t {
                    return m.atom_with_isotope(iso);
                  },
                  "First atom with isotope:|"
                )

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
                .def("atoms_on_shortest_path",
                  [](Molecule& m, atom_number_t a1, atom_number_t a2) ->std::optional<Set_of_Atoms> {
                    Set_of_Atoms result;
                    if (! m.atoms_between(a1, a2, result)) {
                      return std::nullopt;
                    }

                    if (result.empty()) {
                      return std::nullopt;
                    }

                    return result;
                  },
                  "Return list of atoms on the shortest path between a1 and a2"
                )
                .def("all_atoms_between",
                  [](Molecule& m, atom_number_t a1, atom_number_t a2) ->std::optional<Set_of_Atoms> {
                    Set_of_Atoms result;
                    if (! m.AllAtomsBetween(a1, a2, result)) {
                      return std::nullopt;
                    }

                    if (result.empty()) {
                      return std::nullopt;
                    }

                    return result;
                  },
                  "Return all atoms on shortest paths between a1 and a2"
                )
                .def("down_the_bond",
                  [](Molecule& m, atom_number_t a1, atom_number_t a2)->std::optional<Set_of_Atoms> {
                    const int matoms = m.natoms();
                    std::unique_ptr<int[]> dtb = std::make_unique<int[]>(matoms);
                    std::optional<int> maybe_n = m.DownTheBond(a1, a2, dtb.get());
                    if (! maybe_n) {
                      return std::nullopt;
                    }
                    // std::cerr << "Found " << *maybe_n << " atoms down the " << a1 << " " << a2 << " bond\n";
                    Set_of_Atoms result;
                    result.reserve(*maybe_n);
                    for (int i = 0; i < matoms; ++i) {
                      if (i == a2) {
                        continue;
                      }
                      if (dtb[i]) {
                        result << i;
                      }
                    }
                    // std::cerr << "Returning " << result << '\n';
                    return result;
                  },
                  "Return all the atoms found looking down the bond from a1 to a2"
                )

                .def("atoms_by_radius", &AtomsByRadius,
                  py::arg("starting_atoms"), py::arg("max_radius"),
                  "Atoms grouped by minimum bond distance from starting_atoms"
                )

                .def("reset_atom_map_numbers", static_cast<void (Molecule::*)()>(&Molecule::reset_all_atom_map_numbers), "Reset atom map numbers")
                .def("set_atom_map_number", static_cast<void (Molecule::*)(atom_number_t, int)>(&Molecule::set_atom_map_number), "Set atom map number")
                .def("atom_map_number", static_cast<int (Molecule::*)(atom_number_t)const>(&Molecule::atom_map_number), "Set atom map number")
                .def("atom_with_atom_map_number", &Molecule::atom_with_atom_map_number, "Atom with atom map number")

                .def("bond_length", 
                  [](const Molecule& m, atom_number_t a1, atom_number_t a2)->std::optional<float> {
                    if (! m.are_bonded(a1, a2)) {
                      return std::nullopt;
                    }
                    return m.bond_length(a1, a2);  // no need to check again.
                  },
                  "Distance between two bonded atoms"
                )
                .def("bond_angle",
                  [](const Molecule& m, atom_number_t a1, atom_number_t a2, atom_number_t a3)->float {
                    return m.bond_angle(a1, a2, a3, BondedStatus::kOkNotBonded);
                  },
                  "Angle defined by any three atoms - not necessarily bonded"
                )
                .def("dihedral_angle",
                  [](const Molecule& m, atom_number_t a1, atom_number_t a2,
                     atom_number_t a3, atom_number_t a4) {
                    return m.dihedral_angle(a1, a2, a3, a4, BondedStatus::kOkNotBonded);
                  },
                  "Dihedral angle defined by any four atoms - not necessarily bonded"
                )
                .def("signed_dihedral_angle", static_cast<distance_t (Molecule::*)(atom_number_t, atom_number_t, atom_number_t, atom_number_t)const>(&Molecule::signed_dihedral_angle), "Signed dihedral angle involving atoms")
                .def("distance_between_atoms", static_cast<distance_t (Molecule::*)(atom_number_t, atom_number_t)const>(&Molecule::distance_between_atoms), "Spatial distance between atoms")
                .def("longest_intra_molecular_distance", &Molecule::longest_intra_molecular_distance, "longest_intra_molecular_distance")
                .def("translate",
                  [](Molecule& m, distance_t x, distance_t y, distance_t z) {
                    return m.translate_atoms(x, y, z);
                  },
                  "translate coordinates"
                )
                .def ("translate",
                  [](Molecule& m, std::vector<int>& to_move, int flag, float x, float y, float z)->bool {
                    Coordinates coords(x, y, z);
                    m.translate_atoms(coords, to_move.data(), flag);
                    return 1;
                  },
                  "Move all atoms where 'to_move[i]== flag' by (x,y,z)"
                )
                .def("rotate",
                  [](Molecule& m, const std::vector<int>& to_move, int flag,
                               const Coordinates& axis, float angle) {
                    return m.rotate_atoms(axis, angle, to_move.data(), flag);
                  },
                  "Rotate the atoms for which 'to_move[i]==flag' around `axis` by"
                )

                .def("bump_check",
                  [](const Molecule& m, distance_t dist) {
                    return m.bump_check(dist);
                  }
                  ,
                  "the number of non-bonded atoms within dist of each other"
                )

                // This did not work, the name always showed up as empty.
                // Besides, I don't think I want to enable mol.name = xxx, when everything else is via functions.
                //.def_property("name", &Molecule::Name, static_cast<void (Molecule::*)(const std::string&)>(&Molecule::set_name), "name")
                .def("name", &Molecule::Name, "Name")
                .def("sdf_tags", &SdfTags, "Return retained SDF tags as a dict")
                .def("number_records_text_info", &Molecule::number_records_text_info,
                  "Number of retained SDF/text records")
                .def("text_info",
                  [](const Molecule& m)->std::vector<std::string> {
                    std::vector<std::string> result;
                    result.reserve(m.number_records_text_info());
                    for (int i = 0; i < m.number_records_text_info(); ++i) {
                      result.push_back(m.text_info(i).AsString());
                    }
                    return result;
                  },
                  "Retained SDF/text records")
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

                .def("get_coordinates",
                  [](const Molecule& m)->py::array_t<float> {
                    const int matoms = m.natoms();
                    py::array_t<float> result = mkarray_via_buffer<float>(matoms * 3);
                    auto req = result.request();
                    float* ptr = static_cast<float*>(req.ptr);
                    for (int i = 0; i < matoms; ++i) {
                      ptr[i * 3 + 0] = m.x(i);
                      ptr[i * 3 + 1] = m.y(i);
                      ptr[i * 3 + 2] = m.z(i);
                    }
                    return result;
                  },
                  "Return 3*natoms array of coords. First 3 values hold x,y,z of the first atom..."
                )
                .def("set_coordinates",
                  [](Molecule& m, const pybind11::array_t<float> coords)->void {
                    auto req = coords.request();
                    m.SetXyz(static_cast<float*>(req.ptr));
                  },
                  "Set coordinates. First 3 entries are x,y,z for first atom..."
                )
                .def("dihedral_scan",
                  [](Molecule& m, atom_number_t a2, atom_number_t a3, angle_t angle, float bump_check)->std::vector<py::array_t<float>>{
                    const std::vector<std::unique_ptr<float[]>> coords = m.DihedralScan(a2, a3, angle, bump_check);
                    const int nconf = coords.size();
                    const int matoms = m.natoms();

                    std::vector<py::array_t<float>> result;
                    result.reserve(nconf);

                    for (const std::unique_ptr<float[]>& c : coords) {
                      result.emplace_back(py::array_t<float>(matoms * 3));
                      auto req = result.back().request();
                      float* ptr = static_cast<float*>(req.ptr);
                      std::copy_n(c.get(), matoms * 3, ptr);
                    }
                    return result;
                  },
                  "Scan around `a2`-`a3` bond every `angle` degrees and return list of coordinates at each conformer"
                )

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

                .def("gasteiger_partial_charges",
                  [](Molecule& m) -> std::vector<float> {
                    m.compute_Gasteiger_partial_charges();
                    const int matoms = m.natoms();
                    std::vector<float> result;
                    result.reserve(matoms);
                    for (int i = 0; i < matoms; ++i) {
                      result.push_back(m.partial_charge(i));
                    }
                    return result;
                  },
                  "Return list of Gasteiger partial charges"
                )
                .def("move_to_end_of_connection_table",
                  [](Molecule& m, atomic_number_t z) {
                    return m.MoveToEndOfConnectionTable(z);
                  },
                  "Move atoms of type 'z' to end of connection table"
                )

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
                    // It is an open question whether we should temporarily call set_max_matches_to_find.
                    return q.substructure_search(&m);
                  },
                  "Substructure search"
                )

  ;

  py::class_<Atom, std::shared_ptr<Atom>>(m, "Atom")
    .def(py::init<int>())
    .def("atomic_number", static_cast<int (Atom::*)()const>(&Atom::atomic_number), "Atomic Number")
    .def("atomic_symbol",
      [](const Atom& a)->std::string{
        return a.atomic_symbol().AsString();
      },
      "Atomic symbol"
    )
    .def("isotope", static_cast<isotope_t (Atom::*)()const>(&Atom::isotope), "isotope")
    .def("ncon", static_cast<int (Atom::*)()const>(&Atom::ncon), "Number of connections")
    .def("nbonds", static_cast<int (Atom::*)()const>(&Atom::nbonds), "Number of bonds - single=1, double=2...")
    .def("formal_charge", &Atom::formal_charge, "formal charge")
    .def("atomic_weight", &Atom::atomic_weight, "atomic weight")
    .def("exact_mass", &Atom::exact_mass, "exact mass")
    .def("implicit_hydrogens", &Atom::implicit_hydrogens, "number of implicit Hydrogens attached")
    .def("is_bonded_to",
      [](const Atom& a, atom_number_t atom)->bool{
          return a.is_bonded_to(atom);
      },
      "True if atom is bonded to other atom"
    )
    .def("valence_ok", &Atom::valence_ok, "True if valence is ok")
    .def("fully_saturated", &Atom::fully_saturated, "True if fully saturated")
    .def("unsaturation", &Atom::unsaturation, "nbonds() - ncon()")
    .def("other", static_cast<atom_number_t (Atom::*)(atom_number_t, int)const>(&Atom::other), "Other connection")
    .def("is_organic", &Atom::is_organic, "True if the element is organic")
    .def("x", [](const Atom* a)->float {
        return a->x();
      },
      "x coordinate"
    )
    .def("y", [](const Atom* a)->float {
        return a->y();
      },
      "y coordinate"
    )
    .def("z", [](const Atom* a)->float {
        return a->z();
      },
      "z coordinate"
    )
    .def("distance", [](const Atom* a1, const Atom* a2)->float {
        return a1->distance(*a2);
      },
      "spatial distance between atoms"
    )
    .def("distance", [](const Atom* a1, const Coordinates& c)->float {
        return a1->distance(c);
      },
      "spatial distance between atom and point"
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
    .def("__len__",
      [](const Atom& a) {
        return a.ncon();
      },
      "Number of connections"
    )
    .def("__sub__",
      [](const Atom& a1, const Atom& a2)->float {
        return a1.distance(a2);
      },
      "Spatial distance between two atoms"
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
        return ToBondType(b);
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
          result << '#';
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
    .def("any_atoms_in_common", 
      [](const Set_of_Atoms& s1, const Set_of_Atoms& s2)->bool {
        return s1.any_members_in_common(s2);
      },
      "True if any atom in common"
    )
    .def("first_atom_in_common", &Set_of_Atoms::members_in_common, "Atom number of first atom in common")
    .def("atoms_in_common", &Set_of_Atoms::members_in_common, "Number of atoms shared")
    .def("size", &Set_of_Atoms::size, "size")
    .def("set_vector",
      [](const Set_of_Atoms& s, std::vector<int>& v, int value) {
        for (auto a : s) {
          v[a] = value;
        }
      }
    )
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
    .def("__setitem__",
      [](const Set_of_Atoms& me, int ndx, atom_number_t atom) {
        me[ndx] = atom;
        return atom;
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
    .def("__iadd__",
      [](Set_of_Atoms& lhs, const Set_of_Atoms& rhs)->Set_of_Atoms {
        lhs += rhs;
        return lhs;
      },
      "add contents of RHS to LHS returning lhs"
    )
    .def("__iadd__",
      [](Set_of_Atoms& lhs, atom_number_t a)->Set_of_Atoms {
        lhs.add(a);
        return lhs;
      },
      "Add atom `a` to lhs"
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
      [](const Ring& r1, const Ring& r2)->bool{
        return r1.is_fused_to(&r2);
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
  py::class_<Coordinates>(m, "Coordinates")
    .def(py::init<>())
    .def(py::init<float, float, float>())
    .def("normalise", &Coordinates::normalise)
    .def("__repr__", [](const Coordinates& coords) {
      IWString tmp;
      tmp << '(' << coords.x() << ',' << coords.y() << ',' << coords.z() << ')';
      return tmp.AsString();
    })
    .def("__str__", [](const Coordinates& coords) {
      IWString tmp;
      tmp << '(' << coords.x() << ',' << coords.y() << ',' << coords.z() << ')';
      return tmp.AsString();
    })
  ;

  m.def("tetrahedral_chirality", &TetrahedralChirality,
    py::arg("mol"), py::arg("atom"), py::arg("check_is_chiral") = false,
    "LillyMol-defined tetrahedral chiral tag for an atom, or None. Tags are relative to atom-number order with implicit Hydrogen after explicit atoms. If check_is_chiral is true, first verify the atom is actually chiral");
  m.def("is_actually_chiral",
    [](Molecule& mol, atom_number_t zatom)->bool {
      if (zatom < 0 || zatom >= mol.natoms()) {
        throw py::value_error("atom number outside molecule");
      }
      return ::is_actually_chiral(mol, zatom);
    },
    "True if an atom is actually chiral");
  m.def("set_copy_name_in_molecule_copy_constructor", &set_copy_name_in_molecule_copy_constructor, "Copy name in constructor");
  m.def("NumHAcceptors", [](const Molecule& mol) { return mol.LipinskiNumHAcceptors(); }, "Lipinski hydrogen bond acceptor count");
  m.def("NumHDonors", [](Molecule& mol) { return mol.LipinskiNumHDonors(); }, "Lipinski hydrogen bond donor count");
  m.def("fraction_csp3",
      [](Molecule& mol)->double {
        int carbon = 0;
        int csp3 = 0;

        const int matoms = mol.natoms();
        for (int i = 0; i < matoms; ++i) {
          if (mol.atomic_number(i) != 6) {
            continue;
          }

          ++carbon;
          if (mol.saturated(i)) {
            ++csp3;
          }
        }

        if (carbon == 0) {
          return 0.0;
        }

        return static_cast<double>(csp3) / static_cast<double>(carbon);
      },
      "Fraction of carbon atoms that are fully saturated. The denominator is the "
      "carbon count, and a molecule with no carbon gives 0.0 - the same definition "
      "molecule_filter uses for its sp3_carbon_fraction feature.");

  m.def("RDKitNumHAcceptors", [](Molecule& mol) { return mol.RDKitNumHAcceptors(); }, "RDKit compatible acceptor count, NOT a Lipinski count");
  m.def("RDKitNumHDonors", [](Molecule& mol) { return mol.RDKitNumHDonors(); }, "RDKit compatible donor count, NOT a Lipinski count");
  m.def("LillyMolFromSmiles", &MolFromSmiles, "Molecule from smiles");
  m.def("MolFromSmiles",
    [](const std::string& smiles)->std::optional<Molecule> {
      return MolFromSmiles(smiles);
    },
    "Molecule from smiles"
  );
  m.def("MolFromSmiles",
    [](const std::vector<std::string>& smiles) {
      uint32_t nmols = smiles.size();
      std::vector<Molecule> result(nmols);
      for (uint32_t i = 0; i < nmols; ++i) {
        if (! result[i].build_from_smiles(smiles[i])) {
          std::cerr << "Invalid smiles '" << smiles[i] << "' ignored\n";
          result[i].resize(0);
        }
      }
      return result;
    },
    "Return a list of molecules"
  );
  m.def("set_auto_create_new_elements", &set_auto_create_new_elements, "Allow arbitrary two letter elements");
  m.def("set_atomic_symbols_can_have_arbitrary_length", &set_atomic_symbols_can_have_arbitrary_length, "any string is an element");
  m.def("interpret_D_as_deuterium", &element::interpret_d_as_deuterium, "D means '[2H]'");
  m.def("interpret_T_as_deuterium", &element::interpret_t_as_tritium, "T means '[3H]'");
  m.def("set_display_strange_chemistry_messages", &set_display_strange_chemistry_messages, "turn off messages about bad valences");
  m.def("set_display_smiles_interpretation_error_messages", &set_display_smiles_interpretation_error_messages, "Set smiles error messages");
  m.def("count_atoms_in_smiles",
    [](const std::string& smiles) {
      const const_IWSubstring tmp(smiles);
      return lillymol::count_atoms_in_smiles(tmp);
    }
  );

  py::enum_<BondType>(m, "BondType")
    .value("SINGLE_BOND", kSingleBond)
    .value("DOUBLE_BOND", kDoubleBond)
    .value("TRIPLE_BOND", kTripleBond)
    .value("AROMATIC_BOND", kAromaticBond)
    .export_values();
  ;

  py::enum_<ChiralType>(m, "ChiralType")
    .value("CHI_UNSPECIFIED", kChiUnspecified)
    .value("CHI_TETRAHEDRAL_CW", kChiTetrahedralCw)
    .value("CHI_TETRAHEDRAL_CCW", kChiTetrahedralCcw)
    .value("CHI_OTHER", kChiOther)
    .export_values();
  ;

  py::enum_<Hybridization> hybridization_enum(m, "Hybridization");
  hybridization_enum
    .value("UNSPECIFIED", Hybridization::kUnspecified)
    .value("S", Hybridization::kS)
    .value("SP", Hybridization::kSp)
    .value("SP2", Hybridization::kSp2)
    .value("SP3", Hybridization::kSp3)
    .value("SP2D", Hybridization::kSp2d)
    .value("SP3D", Hybridization::kSp3d)
    .value("SP3D2", Hybridization::kSp3d2)
    .value("OTHER", Hybridization::kOther)
    .export_values();
  m.def("hybridization_name",
    [](Hybridization hybridization) {
      return ToString(hybridization);
    },
    "String name for a Hybridization enum value"
  );

  m.def("hybridization",
    [](Molecule& mol, atom_number_t zatom) -> Hybridization {
      if (! mol.ok_atom_number(zatom)) {
        throw py::value_error("hybridization atom number outside [0, natoms)");
      }
      return HybridizationState(mol, zatom);
    },
    "RDKit-like hybridization of atom, computed on demand"
  );

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

  m.def("Position3D",
    [](Molecule& m, atom_number_t atom1, float distance, atom_number_t atom2) {
      return lillymol::Position3D(m, atom1, distance, atom2);
    },
    "Move atoms in fragment defined by atom2 so that fragment can join atom1"
  );

  m.def("xlogp",
    [](Molecule& m)->std::optional<float> {
      xlogp::XLogPCalc calc;
      std::optional<double> result = calc.LogP(m);
      if (! result) {
        return std::nullopt;
      }
      return *result;
    },
    "xlogp"
  );

  py::class_<alogp::ALogP>(m, "ALogP")
    .def(py::init<>())
    .def("set_rdkit_phoshoric_acid_hydrogen", &alogp::ALogP::set_rdkit_phoshoric_acid_hydrogen,
        "mimic RDKit in how Hydrogens on phosphoric acids are handled")
    .def("set_use_alcohol_for_acid", &alogp::ALogP::set_use_alcohol_for_acid,
        "mimic RDKit in how oxygen atoms in acids are handled")
    .def("logp",
         [](alogp::ALogP& mylogp, Molecule& m) {
          return mylogp.LogP(m);
         },
         "Compute AlogP - or None")
  ;

  // Rotatable bonds.
  py::enum_<quick_rotbond::QuickRotatableBonds::RotBond> (m, "RotBond")
    .value("UNDEFINED", quick_rotbond::QuickRotatableBonds::RotBond::kUndefined)
    .value("QUICK", quick_rotbond::QuickRotatableBonds::RotBond::kQuick)
    .value("EXPENSIVE", quick_rotbond::QuickRotatableBonds::RotBond::kExpensive)
    .export_values();
  ;

  py::class_<quick_rotbond::QuickRotatableBonds>(m, "RotatableBonds")
    .def(py::init<>())
    .def("rotatable_bonds", 
      [](quick_rotbond::QuickRotatableBonds& rotb, Molecule& m)->int {
        return rotb.Process(m, nullptr);
      },
      "Number of rotatable bonds in `m`"
    )
    .def("rotatable_bond_atoms",
      [](quick_rotbond::QuickRotatableBonds& rotb, Molecule& m) {
        std::vector<int> bond_rotatable(m.nedges(), 0);
        rotb.Process(m, bond_rotatable.data());

        std::vector<std::tuple<atom_number_t, atom_number_t>> result;
        result.reserve(m.nedges());
        const int nedges = m.nedges();
        for (int bond_number = 0; bond_number < nedges; ++bond_number) {
          if (! bond_rotatable[bond_number]) {
            continue;
          }
          const Bond* b = m.bondi(bond_number);
          result.emplace_back(b->a1(), b->a2());
        }

        return result;
      },
      "Atom pairs defining the rotatable bonds in `m`"
    )
    .def("set_calculation_type", &quick_rotbond::QuickRotatableBonds::set_calculation_type)
  ;

  py::class_<Charge_Assigner>(m, "ChargeAssigner")
    .def(py::init([]() {
      const char* lillymol_home = getenv("LILLYMOL_HOME");
      if (lillymol_home == NULL) {
        std::cerr << "ChargeAssigner:no LILLYMOL_HOME\n";
        throw std::runtime_error("ChargeAssigner::LILLYMOL_HOME not set");
      }

      auto instance = std::make_unique<Charge_Assigner>();
      IWString data(lillymol_home);
      data << "/data/queries/charges";
      if (!instance->BuildFromDir(data)) {
        std::cerr << "Failed to initialise charge assigner\n";
        return instance;
      }

      return instance;
    }))
    .def("active", [](const Charge_Assigner& chg)->bool {
      return chg.active() > 0;
    })
    .def("set_min_distance_between_charges", &Charge_Assigner::set_min_distance_between_charges,
        "specify minimum bond separation between formal charges assigned")
    .def("process", [](Charge_Assigner& chg, Molecule& m) {
      const int matoms = m.natoms();
      if (matoms < 1) {
        return 0;
      }

      std::unique_ptr<formal_charge_t[]> charges_assigned = std::make_unique<formal_charge_t[]>(matoms);
      if (chg.process(m, charges_assigned.get()) == 0) {
        std::cerr << "No charges assigned\n";
        return 0;
      }

      int rc = 0;
      for (int i = 0; i < matoms; ++i) {
        if (charges_assigned[i] != 0) {
          m.set_formal_charge(i, charges_assigned[i]);
          ++rc;
        }
      }

      return rc;
    }
  )
  ;

  py::class_<Donor_Acceptor_Assigner>(m, "DonorAcceptor")
    .def(py::init([]() {
      const char* lillymol_home = getenv("LILLYMOL_HOME");
      if (lillymol_home == NULL) {
        std::cerr << "DonorAcceptor:no LILLYMOL_HOME\n";
        throw std::runtime_error("DonorAcceptor::LILLYMOL_HOME not set");
      }

      auto instance = std::make_unique<Donor_Acceptor_Assigner>();
      IWString data(lillymol_home);
      data << "/data/queries/hbonds";
      static constexpr int kVerbose = 0;
      if (!instance->BuildFromDir(data, kVerbose)) {
        std::cerr << "Failed to initialise donor acceptor\n";
        return instance;
      }
      instance->set_apply_isotopic_labels(1);

      return instance;
    }))
    .def("active", [](const Donor_Acceptor_Assigner& donor_acceptor)->bool {
      return donor_acceptor.active() > 0;
    })
    .def("process", [](Donor_Acceptor_Assigner& donor_acceptor, Molecule& m) {
      return donor_acceptor.process(m);
    })
    ;
}
