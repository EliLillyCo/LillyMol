#include "dicer_api.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/iwbits/fixed_bit_vector.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/target.h"

#include "dicer_lib.h"

namespace dicer_api {

using std::cerr;

using dicer_lib::USPVPTR;

using fixed_bit_vector::FixedBitVector;

PerMoleculeData::PerMoleculeData(Molecule& m) {
  _matoms = m.natoms();
  _can_break = new_int(_matoms * _matoms);
  _usvptr = new USPVPTR[_matoms];
  _xref = new int[_matoms];
  _atype = nullptr;

  _vector_size = 64;
  // Not really, but we should not be processing large moleules.
  if (_matoms > 64) {
    _vector_size = 128;
  }
}

PerMoleculeData::~PerMoleculeData() {
  delete[] _can_break;
  delete[] _usvptr;
  delete[] _xref;
  delete[] _atype;
}

int
PerMoleculeData::AssignAtomTypes(Molecule& m, Atom_Typing_Specification& atom_typing) {
  _atype = new uint32_t[_matoms];
  std::fill_n(_atype, _matoms, 0);
  return atom_typing.assign_atom_types(m, _atype);
}

int
PerMoleculeData::InitialiseAtomPointers(Molecule& m) {
  return dicer_lib::initialise_atom_pointers(m, _atype, _usvptr);
}

// Note that this uses the current atom numbers - it does NOT
// go back to the parent atom numbers. Used with Recap...
int
PerMoleculeData::CanBreak(const Bond& b) const {
  return _can_break[b.a1() * _matoms + b.a2()];
}

atom_number_t
AtomNumberInParent(const Atom& a) {
  const USPVPTR* v = reinterpret_cast<const USPVPTR*>(a.user_specified_void_ptr());
  return v->atom_number_in_parent();
}

int
PerMoleculeData::CanBreak(const Atom& at1, const Atom& at2) const {
  atom_number_t a1 = AtomNumberInParent(at1);
  atom_number_t a2 = AtomNumberInParent(at2);

  return _can_break[a1 * _matoms + a2];
}

Dicer::Dicer() {
  _max_bonds_to_break = 1;
  _min_fragment_size = 1;
  _max_fragment_size = std::numeric_limits<int>::max();

  _break_amide_bonds = 0;
  _break_cc_bonds = 0;
  _break_ring_chain_bonds = 1;

  _label_join_points = 0;

  _determine_fragment_counts = 1;

  _accumulate_global_fragment_count = 0;

  _work_like_recap = 0;
}

int
Dicer::AddBreakBondSmarts(const std::string& smarts) {
  std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
  if (!q->CreateFromSmarts(smarts)) {
    cerr << "Dicer::AddBreakBondSmarts:invalid smarts '" << smarts << '\n';
    return 0;
  }

  q->set_find_unique_embeddings_only(1);
  _query << q.release();

  return 1;
}

int
Dicer::AddBreakBondQuery(const std::string& fname) {
  const_IWSubstring directive(fname);
  static constexpr int kVerbose = 1;

  if (!process_cmdline_token('*', directive, _query, kVerbose)) {
    cerr << "directive::AddBreakBondQuery:invalid directive '" << directive << "'\n";
    return 0;
  }

  return 1;
}

int
Dicer::AddFragmentRequirementSmarts(const std::string& smarts) {
  std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
  if (!q->CreateFromSmarts(smarts)) {
    cerr << "Dicer::AddBreakBondSmarts:invalid smarts '" << smarts << '\n';
    return 0;
  }

  q->set_find_unique_embeddings_only(1);
  _fragments_must_contain << q.release();

  return 1;
}

int
Dicer::AddFragmentDisqualifierSmarts(const std::string& smarts) {
  std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
  if (!q->CreateFromSmarts(smarts)) {
    cerr << "Dicer::AddBreakBondSmarts:invalid smarts '" << smarts << '\n';
    return 0;
  }

  q->set_find_unique_embeddings_only(1);
  _fragments_cannot_contain << q.release();

  return 1;
}

int
Dicer::set_perceive_symmetry_equivalent_matches(int s) {
  if (_query.empty()) {
    cerr << "Dicer::set_perceive_symmetry_equivalent_matches:no queries\n";
    return 0;
  }

  for (Substructure_Query* q : _query) {
    q->set_perceive_symmetry_equivalent_matches(s);
  }

  return _query.number_elements();
}

int
Dicer::set_atom_type(const std::string& s) {
  const_IWSubstring tmp(s);
  if (!_atom_typing.build(tmp)) {
    cerr << "Dicer::set_atom_type:invalid atom typing '" << s << "'\n";
    return 0;
  }

  return 1;
}

int
Dicer::IdentifyBondsToBreak(Molecule& m, int* can_break) {
  if (_query.size()) {
    return IdentifyBondsToBreakViaQueries(m, can_break);
  }

  return IdentifyBondsToBreakHardCodedRules(m, can_break);
}

int
Dicer::IdentifyBondsToBreakViaQueries(Molecule& m, int* can_break) {
  Molecule_to_Match target(&m);

  const int matoms = m.natoms();

  int rc = 0;
  for (Substructure_Query* q : _query) {
    Substructure_Results sresults;
    if (!q->substructure_search(target, sresults)) {
      continue;
    }
    for (const Set_of_Atoms* e : sresults.embeddings()) {
      atom_number_t a1 = e->item(0);
      atom_number_t a2 = e->item(1);
      can_break[a1 * matoms + a2] = 1;
      can_break[a2 * matoms + a1] = 1;
      ++rc;
    }
  }

  return rc;
}

/*
   The hard coded rules are
   aromatic - *
   heteroatom - *
   unsaturated - *
   heteroatom -C - *
 */

int
Dicer::IdentifyBondsToBreakHardCodedRules(Molecule& m, int* can_break) {
  m.compute_aromaticity_if_needed();

  int rc = 0;

  const int matoms = m.natoms();

  for (const Bond* b : m.bond_list()) {
    if (b->nrings()) {
      continue;
    }

    if (!b->is_single_bond()) {
      continue;
    }

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();

    const int rbc1 = m.ring_bond_count(a1);
    const int rbc2 = m.ring_bond_count(a2);

    if (_break_amide_bonds) {
      ;
    } else if (rbc1 || rbc2) {
      ;
    } else if (dicer_lib::is_amide(m, a1, a2)) {
      continue;
    }

    // Allow biphenyl type bonds to break.
    if (m.is_aromatic(a1) || m.is_aromatic(a2)) {
      can_break[a1 * matoms + a2] = 1;
      can_break[a2 * matoms + a1] = 1;
      rc++;
      continue;
    }

    if (!m.saturated(a1) || !m.saturated(a2)) {
      can_break[a1 * matoms + a2] = 1;
      can_break[a2 * matoms + a1] = 1;
      rc++;
      continue;
    }

    if (!_break_ring_chain_bonds) {
    } else if ((rbc1 == 0 && rbc2 > 0) || (rbc1 > 0 && rbc2 == 0)) {
      can_break[a1 * matoms + a2] = 1;
      can_break[a2 * matoms + a1] = 1;
      rc++;
      continue;
    }

    const atomic_number_t z1 = m.atomic_number(a1);
    const atomic_number_t z2 = m.atomic_number(a2);

    // heteroatoms at both ends of bond
    if (6 != z1 && 6 != z2) {
      can_break[a1 * matoms + a2] = 1;
      can_break[a2 * matoms + a1] = 1;
      rc++;
      continue;
    }

    //  At least one heteroatom. Don't split up CF3 groups
    // This seems wrong, investigate...

    if (m.is_halogen(a1) || m.is_halogen(a2)) {
      continue;
    }

    //  Do allow the bond to a CF3 or t-Butyl to break

    if (dicer_lib::is_cf3_like(m, *b)) {
      can_break[a1 * matoms + a2] = 1;
      can_break[a2 * matoms + a1] = 1;
      rc++;
      continue;
    }

    if (6 != z1 || 6 != z2) {
      can_break[a1 * matoms + a2] = 1;
      can_break[a2 * matoms + a1] = 1;
      rc++;
      continue;
    }

    if (_break_cc_bonds && 6 == z1 && 6 == z2) {
      can_break[a1 * matoms + a2] = 1;
      can_break[a2 * matoms + a1] = 1;
      rc++;
      continue;
    }
  }

// #define ECHO_BONDS_BROKEN
#ifdef ECHO_BONDS_BROKEN
  cerr << m.isotopically_labelled_smiles() << ' ' << m.name() << endl;
  debug_print(cerr);
#endif

  return rc;
}

int
Dicer::Dice(Molecule& m, std::unordered_map<std::string, uint32_t>& fragments) {
  PerMoleculeData pmd(m);

  if (!IdentifyBondsToBreak(m, pmd.can_break())) {
    cerr << "Dicer::Dice:no bonds to break in " << m.smiles() << '\n';
    return 0;
  }

  if (_work_like_recap) {
    return Recap(m, pmd, fragments);
  }

  if (_atom_typing.active()) {
    pmd.AssignAtomTypes(m, _atom_typing);
  }

  pmd.InitialiseAtomPointers(m);

  set_copy_user_specified_atom_void_ptrs_during_create_subset(1);
  set_copy_atom_based_user_specified_void_pointers_during_add_molecule(1);

  const int rc = Dice(m, pmd, 1, fragments);

  if (_accumulate_global_fragment_count) {
    UpdateGlobalFragmentCount(fragments);
  }

  return rc;
}

int
Dicer::Dice(Molecule& m, PerMoleculeData& pmd, int bonds_broken,
            std::unordered_map<std::string, uint32_t>& fragments) {
  const int matoms = m.natoms();

  if (matoms < _min_fragment_size) {
    return 0;
  }

  // cerr << "Begin dicing " << m.smiles() << '\n';

  // Save the bonds that can break. We need to collect them because we alter the
  // bond list during processing.
  Set_of_Atoms b1, b2;
  b1.reserve(matoms);
  b2.reserve(matoms);

  for (const Bond* b : m.bond_list()) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (!pmd.CanBreak(m[a1], m[a2])) {
      continue;
    }

    b1 << b->a1();
    b2 << b->a2();
  }

  std::unique_ptr<int[]> fragment_membership = std::make_unique<int[]>(matoms);

  const int nbonds = b1.number_elements();
  for (int i = 0; i < nbonds; ++i) {
    const atom_number_t a1 = b1[i];
    const atom_number_t a2 = b2[i];

    auto starting_isotopes = BreakTheBond(m, a1, a2);

    m.fragment_membership(fragment_membership.get());

    MakeTwoFragments(m, fragment_membership.get(), pmd, bonds_broken, fragments);

    Restore(m, a1, a2, starting_isotopes);  // add the a1-a2 bond back.
  }

  return fragments.size();
}

int
Dicer::MakeTwoFragments(Molecule& m, const int* fragment_membership, PerMoleculeData& pmd,
                        int bonds_broken,
                        std::unordered_map<std::string, uint32_t>& fragments) {
  assert(m.number_fragments() == 2);
#ifdef DEBUG_MAKE_TWO_FRAGMENTS
  cerr << "Making two fragments from " << m.smiles() << " bb " << bonds_broken << '\n';
#endif
  if (m.number_fragments() != 2) {
    cerr << "Fragment wrong " << m.smiles() << '\n';
  }

  const int matoms = m.natoms();

  // Process each of the fragments
  // TODO:ianwatson  we don't need to form the second bitvector, it is just the complement
  // of the first. Implement that...
  for (int frag_num = 0; frag_num < 2; ++frag_num) {
    if (!_determine_fragment_counts) {  // Suppress duplicates
      std::unique_ptr<FixedBitVector> fp =
          std::make_unique<FixedBitVector>(pmd.vector_size());

      for (int i = 0; i < matoms; ++i) {
        if (fragment_membership[i] == frag_num) {
          fp->set_bit(AtomNumberInParent(m[i]));
        }
      }

      if (pmd.AlreadyFound(fp)) {
#ifdef DEBUG_MAKE_TWO_FRAGMENTS
        Molecule frag;
        m.create_subset(frag, fragment_membership, frag_num, pmd.xref());
        cerr << "From " << m.unique_smiles() << " duplicate " << frag.unique_smiles()
             << '\n';
#endif
        continue;
      }
    }

    Molecule frag;
    m.create_subset(frag, fragment_membership, frag_num, pmd.xref());
    // cerr << "From " << m.smiles() << " created fragment " << frag.smiles() << " bb " <<
    // bonds_broken << '\n';
    MaybeStoreFragment(frag, fragments);

    if (bonds_broken < _max_bonds_to_break) {
      Dice(frag, pmd, bonds_broken + 1, fragments);
    }
  }

  return 1;
}

int
PerMoleculeData::AlreadyFound(std::unique_ptr<FixedBitVector>& fp) {
  for (const FixedBitVector* seen : _found) {
    if (*seen == *fp) {
      // cerr << "Dup " << *seen << " (" << seen->nset() << ") and " << *fp << " (" <<
      // fp->nset() << ")\n";
      return 1;
    }
  }

  _found << fp.release();

  return 0;
}

std::tuple<isotope_t, isotope_t>
Dicer::BreakTheBond(Molecule& m, atom_number_t a1, atom_number_t a2) const {
  m.remove_bond_between_atoms(a1, a2);

  std::tuple<isotope_t, isotope_t> result = {m.isotope(a1), m.isotope(a2)};

  if (_label_join_points) {
    m.set_isotope(a1, _label_join_points);
    m.set_isotope(a2, _label_join_points);
  }

  return result;
}

void
Dicer::Restore(Molecule& m, atom_number_t a1, atom_number_t a2,
               const std::tuple<isotope_t, isotope_t>& iso) const {
  m.add_bond(a1, a2, SINGLE_BOND);
  if (_label_join_points) {
    m.set_isotope(a1, std::get<0>(iso));
    m.set_isotope(a2, std::get<1>(iso));
  }
}

int
Dicer::MaybeStoreFragment(Molecule& m,
                          std::unordered_map<std::string, uint32_t>& fragments) {
  // cerr << "MaybeStoreFragment " << m.unique_smiles() << '\n';
  const int matoms = m.natoms();
  if (matoms < _min_fragment_size) {
    return 0;
  }

  if (matoms > _max_fragment_size) {
    return 0;
  }

  const IWString& usmi = m.unique_smiles();
  const std::string key(usmi.data(), usmi.size());
  // cerr << "key '" << key << "'\n";
  auto iter = fragments.find(key);
  if (iter != fragments.end()) {
    ++iter->second;
    return 0;
  }

  fragments[key] = 1;

  return 1;
}

void
Dicer::UpdateGlobalFragmentCount(const std::unordered_map<std::string, uint32_t>& frags) {
  for (const auto& [usmi, count] : frags) {
    auto iter = _global.find(usmi);
    if (iter == _global.end()) {
      _global[usmi] = 1;
    } else {
      ++iter->second;
    }
  }
}

int
Dicer::Recap(Molecule& m, PerMoleculeData& pmd, 
            std::unordered_map<std::string, uint32_t>& fragments) {
  const int nedges = m.nedges();
  for (int i = nedges - 1; i >= 0; --i) {
    const Bond* b = m.bondi(i);
    if (! pmd.CanBreak(*b)) {
      continue;
    }

    if (_label_join_points) {
      m.set_isotope(b->a1(), _label_join_points);
      m.set_isotope(b->a2(), _label_join_points);
    }

    m.remove_bond_between_atoms(b->a1(), b->a2());
  }

  resizable_array_p<Molecule> components;
  m.create_components(components);

  for (Molecule* frag : components) {
    const IWString& usmi = frag->unique_smiles();

    const std::string as_string(usmi.data(), usmi.size());

    auto [iter, success] = fragments.try_emplace(std::move(as_string), 1);
    if (success) {
    } else {
      ++iter->second;
    }
  }

  return fragments.size();
}

}  // namespace dicer_api
