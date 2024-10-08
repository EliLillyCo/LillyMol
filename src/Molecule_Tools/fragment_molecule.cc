#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/target.h"

#include "fragment_molecule.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace fragment_molecule {

using std::cerr;

MoleculeFragmenter::MoleculeFragmenter() {
  _add_queries_to_default_rules = 0;
  _break_carbon_carbon_bonds = 0;
}

int
MoleculeFragmenter::Initialise(Command_Line& cl) {
  int verbose = cl.option_present('v');
  if (! process_queries(cl, _to_break, verbose, 'q')) {
    cerr << "MoleculeFragmenter::Initialise:canot process to break queries (-q)\n";
    return 0;
  }

  if (! process_queries(cl, _never_break, verbose, 'Q')) {
    cerr << "MoleculeFragmenter::Initialise:canot process never break queries (-Q)\n";
    return 0;
  }

  return 1;
}

constexpr int kUnknown = 0;
constexpr int kCanBreak = 1;
constexpr int kDoNotBreak = 2;

// All bonds implied by the matched atoms are marked with `flag` in `status`.
int
IdentifyMatchingBonds(Molecule_to_Match& target,
                       resizable_array_p<Substructure_Query>& query,
                       int* status,
                       int flag) {
  const int matoms = target.natoms();

  int rc = 0;
  for (Substructure_Query* q : query) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      continue;
    }

    for (const Set_of_Atoms* e : sresults.embeddings()) {
      if (e->size() < 2) {
        continue;
      }
      const int eatoms = e->number_elements();
      for (int i = 0; i < eatoms; ++i) {
        atom_number_t a1 = e->item(i);
        for (int j = i + 1; j < eatoms; ++j) {
          atom_number_t a2 = e->item(j);
          if (! target.molecule()->are_bonded(a1, a2)) {
            continue;
          }
          status[a1 * matoms + a2] = flag;
          status[a2 * matoms + a1] = flag;
          ++rc;
        }
      }
    }
  }

  return rc;
}

int
MoleculeFragmenter::IdentifyBreakableBonds(Molecule& m, resizable_array<int>& bonds_to_break) {
  const int matoms = m.natoms();
  if (matoms < 2) {
    return 0;
  }

  std::unique_ptr<int[]> status(new_int(matoms * matoms, kUnknown));

  if (IdentifyBreakableBonds(m, status.get())== 0) {
    //cerr << "No bonds to be broken\n";
    return 0;
  }

  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      const int s = status[i * matoms + j];

      if (s != kCanBreak) {
        continue;
      }

      int bnum = m.which_bond(i, j);
      if (bnum < 0) {
        continue;
      }

      bonds_to_break << bnum;
      // cerr << "Add " << bnum << " to the results, atoms " << i << " " << j << '\n';
    }
  }

  // now sort
  if (bonds_to_break.size() < 2) {
    return bonds_to_break.size();
  }

  if (bonds_to_break.size() == 2) {
    if (bonds_to_break[0] > bonds_to_break[1]) {
      return 2;
    }
    bonds_to_break.swap_elements(0, 1);
    return 2;
  }

  bonds_to_break.iwqsort_lambda([](int b1, int b2) {
    if (b1 > b2) {
      return -1;
    }
    return 1;
  });

  return bonds_to_break.size();
}

int
SuppressGuainidine(Molecule& m,
                   atom_number_t zatom,
                   int* status) {
  const Atom& a = m[zatom];
  atom_number_t doubly_bonded_nitrogen = INVALID_ATOM_NUMBER;
  atom_number_t singly_bonded_nitrogen1 = INVALID_ATOM_NUMBER;
  atom_number_t singly_bonded_nitrogen2 = INVALID_ATOM_NUMBER;
  for (const Bond* b : a) {
    atom_number_t n = b->other(zatom);
    if (m.atomic_number(n) != 7) {
      continue;
    }
    if (b->is_double_bond()) {
      doubly_bonded_nitrogen = n;
    } else if (singly_bonded_nitrogen1 == INVALID_ATOM_NUMBER) {
      singly_bonded_nitrogen1 = n;
    } else if (singly_bonded_nitrogen2 == INVALID_ATOM_NUMBER) {
      singly_bonded_nitrogen2 = n;
    } else {
      return 0;
    }
  }

  if (doubly_bonded_nitrogen == INVALID_ATOM_NUMBER) {
    return 0;
  }
  if (singly_bonded_nitrogen2 == INVALID_ATOM_NUMBER) {
    return 0;
  }

  const int matoms = m.natoms();

  status[zatom * matoms + doubly_bonded_nitrogen] = status[doubly_bonded_nitrogen * matoms + zatom] = kDoNotBreak;
  status[zatom * matoms + singly_bonded_nitrogen1] = status[singly_bonded_nitrogen1 * matoms + zatom] = kDoNotBreak;
  status[zatom * matoms + singly_bonded_nitrogen2] = status[singly_bonded_nitrogen2 * matoms + zatom] = kDoNotBreak;

  return 1;
}

int
SuppressGuainidine(Molecule& m,
                   int* status) {
  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m[i];
    if (a.atomic_number() != 6) {
      continue;
    }
    if (a.ncon() != 3) {
      continue;
    }

    if (SuppressGuainidine(m, i, status)) {
      ++rc;
    }
  }

  return rc;
}

// `nh` is the Nitrogen in a N=C bond with `carbon` the other end.
int
SuppressAmidine(Molecule& m,
                atom_number_t nh,
                atom_number_t carbon,
                int* status) {
  const Atom& a = m[carbon];
  atom_number_t singly_bonded_nitrogen = INVALID_ATOM_NUMBER;
  for (const Bond * b : a) {
    atom_number_t o = b->other(carbon);
    if (o == nh) {
      continue;
    }
    const Atom& n = m[o];
    if (n.atomic_number() != 7) {
      continue;
    }
    if (n.ncon() != 1) {
      continue;
    }
    singly_bonded_nitrogen = o;
  }

  if (singly_bonded_nitrogen == INVALID_ATOM_NUMBER) {
    return 0;
  }

  const int matoms = m.natoms();

  status[carbon * matoms + singly_bonded_nitrogen] = status[singly_bonded_nitrogen * matoms + carbon] = kDoNotBreak;
  status[carbon * matoms + nh] = status[nh * matoms + carbon] = kDoNotBreak;

  return 1;
}

// Suppress terminal amidine
// Start by looking for the N= atom, then explore around the carbon.
int
SuppressAmidine(Molecule& m,
                int* status) {
  int rc = 0;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m[i];
    if (a.atomic_number() != 7) {
      continue;
    }
    if (a.ncon() != 1) {
      continue;
    }

    const Bond* b = a[0];
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t carbon = b->other(i);
    rc += SuppressAmidine(m, i, carbon, status);
  }

  return rc;
}

// If `carbon` is CF3, CCl3 or t-butyl, update the entries in `status` to `flag`.
int
IsCf3(Molecule& m,
      atom_number_t carbon,
      int* status) {
  assert(m.ncon(carbon) == 4);

  Set_of_Atoms singly_connected;
  atom_number_t not_singly_connected = INVALID_ATOM_NUMBER;
  atomic_number_t other_type = -1;

  for (const Bond* b : m[carbon]) {
    if (! b->is_single_bond()) {
      return 0;
    }

    atom_number_t o = b->other(carbon);
    if (m.ncon(o) == 1) {
      if (singly_connected.empty()) {
        other_type = m.atomic_number(o);
      } else if (m.atomic_number(o) != other_type) {
        return 0;
      }
      singly_connected << o;
    } else {
      if (not_singly_connected >= 0) {
        return 0;
      }
      not_singly_connected = o;
    }
  }

  if (singly_connected.size() != 3) {
    return 0;
  }
  // Probably unnecessary.
  if (not_singly_connected == INVALID_ATOM_NUMBER) {
    return 0;
  }

  const int matoms = m.natoms();
  for (atom_number_t o : singly_connected) {
    status[carbon * matoms + o] = status[o * matoms + carbon] = kDoNotBreak;
  }
  status[carbon * matoms + not_singly_connected] = status[not_singly_connected * matoms + carbon] = kCanBreak;

  // cerr << "IsCf3 set " << carbon << " and " << not_singly_connected << '\n';
  return 1;
}


int
MoleculeFragmenter::IdentifyBreakableBonds(Molecule& m, int* status) {

  int rc = 0;
  if (_to_break.empty() && _never_break.empty()) {
  } else {
    Molecule_to_Match target(&m);
    rc += IdentifyMatchingBonds(target, _to_break, status, kCanBreak);
    (void) IdentifyMatchingBonds(target, _never_break, status, kDoNotBreak);
    if (! _add_queries_to_default_rules) {
      return rc;
    }
  }

  const int matoms = m.natoms();

  // First flag any CF3 type motifs
  int nitrogen_count = 0;

  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m[i];

    if (a.atomic_number() != 6) {
      if (a.atomic_number() == 7) {
        ++nitrogen_count;
      }
      continue;
    }
    if (a.ncon() != 4) {
      continue;
    }

    rc += IsCf3(m, i, status);
  }

  if (nitrogen_count >= 3) {
    int n = SuppressGuainidine(m, status);
    nitrogen_count -= (3 * n); 
  }
  if (nitrogen_count >= 2) {
    SuppressAmidine(m, status);
  }

  m.compute_aromaticity_if_needed();

  for (const Bond* b : m.bond_list()) {
    if (b->is_aromatic()) {
      continue;
    }
    if (! b->is_single_bond()) {
      continue;
    }
    if (b->nrings()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    rc += IdentifyBreakableBond(m, a1, a2, status);
    // cerr << *b << " rc now " << rc << '\n';
  }

  return rc;
}

// Return true if this looks like an amide or an acid
int
IsAmide(Molecule& m,
        atom_number_t nitrogen,
        atom_number_t cs) {
  assert(m.atomic_number(nitrogen) == 7 || m.atomic_number(nitrogen) == 8);

  // cerr << "IsAmide nitrogen " << nitrogen << " " << m.atomic_number(nitrogen) << '\n';
  const Atom& acs = m[cs];
  if (acs.atomic_number() == 6) {
    if (acs.ncon() != 3) {
      return 0;
    }
  } else if (acs.atomic_number() == 16) {
  } else {
    return 0;
  }
  // cerr << "cs " << cs << ' ' << m.smarts_equivalent_for_atom(cs) << '\n';

  // Look for =O or =O groups.
  for (const Bond* b : acs) {
    // cerr << "From " << cs << " to " << *b << " type " << m.smarts_equivalent_for_atom(b->other(cs)) << '\n';
    if (! b->is_double_bond()) {
      continue;
    }
    atom_number_t o = b->other(cs);
    if (o == nitrogen) {
      continue;
    }

    const Atom& ato = m[o];
    if (ato.ncon() != 1) {
      return 0;
    }
    if (ato.atomic_number() == 8) {
      return 1;
    }
    if (ato.atomic_number() == 16) {
      return 1;
    }

    return 0;
  }

  return 0;
}

// We have two non-ring atoms joined by a single bond. Is this bond
// breakable?
int
MoleculeFragmenter::IdentifyBreakableBond(Molecule& m,
                                           atom_number_t a1, atom_number_t a2,
                                           int* status) {
  const int matoms = m.natoms();

  if (status[a1 * matoms + a2] != kUnknown) {
    return 0;
  }

  // cerr << "AToms " << a1 << ',' << a2 << '\n';
  if (m.ring_bond_count(a1) || m.ring_bond_count(a2)) {
    status[a1 * matoms + a2] = status[a2 * matoms + a1] = kCanBreak;
    return 1;
  }

  const Atom& atom1 = m[a1];
  const Atom& atom2 = m[a2];

  const atomic_number_t z1 = atom1.atomic_number();
  const atomic_number_t z2 = atom2.atomic_number();

  const int unsat1 = atom1.unsaturated();
  const int unsat2 = atom2.unsaturated();

  if (z1 == 6 && z2 == 6) {
    if (unsat1 || unsat2) {
    } else if (! _break_carbon_carbon_bonds) {
      return 0;
    }
    status[a1 * matoms + a2] = status[a2 * matoms + a1] = kCanBreak;
    return 1;
  }

  if (z1 == 7 && (z2 == 6 || z2 == 16) && unsat2) {
    if (IsAmide(m, a1, a2)) {
      // cerr << "IsAmide returns true, not breakable\n";
      return 0;
    }
  }

  if (z2 == 7 && (z1 == 6 || z1 == 16) && unsat1) {
    if (IsAmide(m, a2, a1)) {
      // cerr << "IsAmide returns true, not breakable\n";
      return 0;
    }
  }

  // do not break acids
  if (z1 == 8 && (z2 == 6 || z2 == 16) && unsat2) {
    if (IsAmide(m, a1, a2)) {
      // cerr << "IsAmide returns true, not breakable\n";
      return 0;
    }
  }
  if (z2 == 8 && (z1 == 6 || z1 == 16) && unsat1) {
    if (IsAmide(m, a2, a1)) {
      // cerr << "IsAmide returns true, not breakable\n";
      return 0;
    }
  }

  status[a1 * matoms + a2] = status[a2 * matoms + a1] = kCanBreak;

  // std::cerr << "Brekable bond btw " << a1 << " and " << a2 << " " << m.smarts_equivalent_for_atom(a1) << " and " << m.smarts_equivalent_for_atom(a2) << '\n';
  return 1;
}

}  // namespace fragment_molecule
