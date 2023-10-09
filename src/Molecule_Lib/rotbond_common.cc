/*
  Rotatable bonds are computed in iwdescr and in the rotatable_bonds programme.
  In order to keep these in sync as much as possible, we include as many common
  functions here as possible
*/

#include <stdlib.h>

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "molecule.h"
#include "rotbond_common.h"

using std::cerr;

/*
  A bond might be considered rotatable. Is there a triple bond at either end?
*/

int
triple_bond_at_either_end(const Molecule & m,
                          const Bond * b)
{
  Atom * a1 = const_cast<Atom *>(m.atomi(b->a1()));     // nbonds() is non-const

  if (2 == a1->ncon() && 4 == a1->nbonds())
    return 1;

  Atom * a2 = const_cast<Atom *>(m.atomi(b->a2()));     // nbonds() is non-const

  if (2 == a2->ncon() && 4 == a2->nbonds())
    return 1;

  return 0;        // nope, no triple bonds here
}

int
is_non_rotatable_amide(Molecule & m,
                       atom_number_t n,
                       atom_number_t c)
{
  const Atom * nitrogen = m.atomi(n);
  const Atom * carbon = m.atomi(c);

// Swap things around to get the pointers pointing to the right atoms

  if (7 == nitrogen->atomic_number() && 6 == carbon->atomic_number())    // great, got it first time
    ;
  else if (6 == nitrogen->atomic_number() && 7 == carbon->atomic_number())
  {
    std::swap(n, c);
    std::swap(nitrogen, carbon);
  }
  else
    return 0;

  if (2 != nitrogen->ncon())    // only O=C-[NH]- is non-rotatable
    return 0;

  if (1 != m.hcount(n))    // probably not necessary to test
    return 0;

  int ccon = carbon->ncon();
  int cbonds = carbon->nbonds();

  if (3 == ccon && 4 == cbonds)      // [NH]-C(=O)-*
    ;
  else if (2 == ccon && 3 == cbonds)    // [NH]-[CH]=O
    ;
  else
    return 0;

  int doubly_bonded_oxygens = 0;

  for (int i = 0; i < ccon; i++)
  {
    const Bond * b = carbon->item(i);
    if (! b->is_double_bond())
      continue;

    atom_number_t o = b->other(c);

    atomic_number_t z = m.atomic_number(o);
    if (8 == z)
      ;
//  else if (16 == z)    // allow sulphur
//    ;
    else
      return 0;

    doubly_bonded_oxygens++;
  }

  return 1 == doubly_bonded_oxygens;
}

static int
is_cf3_or_t_butyl(Molecule & m,
                  const atom_number_t zatom)
{
  const Atom * a = m.atomi(zatom);

  int n = 0;    // the number of F or C's attached

  atomic_number_t z = kInvalidAtomicNumber;

  for (int i = 0; i < 4; i++)
  {
    atom_number_t j = a->other(zatom, i);

    const Atom * aj = m.atomi(j);

    if (1 != aj->ncon())
      continue;

    if (kInvalidAtomicNumber == z)
    {
      z = aj->atomic_number();
      n = 1;
    }
    else if (z == aj->atomic_number())
      n++;
  }

  return 3 == n;    // we need 3 singly connected C or F neighbours
}

static int
is_cf3_or_t_butyl (Molecule & m,
                   atom_number_t zatom1,
                   atom_number_t zatom2)
{
  const Atom * a1 = m.atomi(zatom1);

  if (4 == a1->ncon() && 4 == a1->nbonds())
  {
    if (is_cf3_or_t_butyl(m, zatom1))
      return 1;
  }

  const Atom * a2 = m.atomi(zatom2);

  if (4 == a2->ncon() && 4 == a2->nbonds())
  {
    if (is_cf3_or_t_butyl(m, zatom2))
      return 1;
  }


  return 0;
}

int
is_non_rotatable_sulphonamide(Molecule & m,
                              atom_number_t zatom1,
                              atom_number_t zatom2)
{
  const Atom * a1 = m.atomi(zatom1);
  const Atom * a2 = m.atomi(zatom2);

  if (16 == a1->atomic_number() && 7 == a2->atomic_number())
    ;
  else if (7 == a1->atomic_number() && 16 == a2->atomic_number())
  {
    std::swap(zatom1, zatom2);
    std::swap(a1, a2);
  }
  else
    return 0;

#ifdef NOT_SURE_IF_THESE_ARE_NEEDED_OR_NOT
  if (2 != a2->ncon())    // only O=C-[NH]- is non-rotatable
    return 0;

  if (1 != m.hcount(zatom2))    // probably not necessary to test
    return 0;
#endif

  if (4 != a1->ncon())
    return 0;

  if (6 != a1->nbonds())
    return 0;

  int doubly_bonded_oxygen = 0;

  for (int i = 0; i < 4; i++)
  {
    const Bond * b = a1->item(i);

    if (! b->is_double_bond())
      continue;

    atom_number_t o = b->other(zatom1);

    if (1 != m.ncon(o))
      continue;

    if (8 != m.atomic_number(o))
      continue;

    doubly_bonded_oxygen++;
  }

//cerr << "doubly_bonded_oxygen " << doubly_bonded_oxygen  << endl;
  return 2 == doubly_bonded_oxygen;
}

/*
  The bond between A1 and A2 might be rotatable. Is it part of a CF3, t-Butyl or O=C-[NH]-
*/

int
part_of_otherwise_non_rotabable_entity(Molecule & m,
                                       atom_number_t a1,
                                       atom_number_t a2)
{
//cerr << " atom " << m.smarts_equivalent_for_atom(a1) << ' ' << m.smarts_equivalent_for_atom(a2) << " is cf3 " << is_cf3_or_t_butyl(m, a1, a2) << endl;

  if (is_cf3_or_t_butyl(m, a1, a2))
    return 1;

  if (is_non_rotatable_amide(m, a1, a2))    // O=C-[NH]- is non-rotatable
    return 1;

  if (is_non_rotatable_sulphonamide(m, a1, a2))    // O=S(=O)-[NH]- is non-rotatable
    return 1;

  return 0;
}

namespace quick_rotbond {

QuickRotatableBonds::QuickRotatableBonds() {
  _calculation = RotBond::kUndefined;
  _isotope = 0;
}

void
DisplayQrbHelp(char flag, std::ostream& output) {
  output << "To control rotatable bond calculation the following directives are recognised\n";
  output << " -" << flag << " fast              all non-ring, [D>1]-[D>1] single bonds are rotatable\n";
  output << " -" << flag << " better            CF3, t-Butyl and amides are excluded\n";
}

int
QuickRotatableBonds::Initialise(Command_Line& cl,
                                char flag) {
  IWString s;
  for (int i = 0; cl.value(flag, s, i); ++i) {
    if (s == "fast") {
      _calculation = RotBond::kQuick;
    } else if (s == "better" ) {
      _calculation = RotBond::kExpensive;
    } else if (s.starts_with("iso=")) {
      s.remove_leading_chars(4);
      if (! s.numeric_value(_isotope)) {
        cerr << "QuickRotatableBonds::Initialise:invalid iso= directive '" << s << "'\n";
        return 0;
      }
    } else if (s == "help") {
      DisplayQrbHelp(flag, cerr);
      return 0;
    } else {
      cerr << "QuickRotatableBonds::Initialise:unrecognised -" << flag << " qualifier '" << s << '\n';
      DisplayQrbHelp(flag, cerr);
      return 0;
    }
  }

  return 1;
}

// Does `zatom` look like the C in a CF3 or t-butyl group.
// Atom `other` is attached to `zatom` but is known to not
// be a terminal group.
// Keep track of the attached atomic numbers, and if there
// are at least 3 all the same, return true.
int
IsCf3(const Molecule& m,
          atom_number_t zatom,
          atom_number_t other) {
  const Atom& a = m.atom(zatom);

  // All atomic numbers attached must be the same.
  atomic_number_t attached = -1;
  int count = 0;

  for (const Bond* b : a) {
    if (! b->is_single_bond()) {
      continue;
    }
    atom_number_t j = b->other(zatom);
    if (j == other) {
      continue;
    }

    atomic_number_t jz = m.atomic_number(j);
    if (attached < 0) {
      attached = jz;
      count = 1;
    } else if (jz != attached) {
      return 0;
    } else {
      ++count;
    }
  }

  // Include CF2 and C(CH3)2
  return count >= 2;
}

// Return true if `carbon` is doubly bonded to one or more =O or =S
// A previous check has verified that `carbon` is attached to a Nitrogen.
int
IsAmide(const Molecule& m,
        atom_number_t carbon) {
  int doubly_bonded_oxygen = 0;

  const Atom& ac = m.atom(carbon);
  for (const Bond * b : ac) {
    if (! b->is_double_bond()) {
      continue;
    }

    const atomic_number_t o = m.atomic_number(b->other(carbon));
    if (o == 8 || o == 16) {
      ++doubly_bonded_oxygen;
    }
  }

  return doubly_bonded_oxygen;
}

int
IsAmide(Molecule& m,
        const Bond* b) {
  const atomic_number_t z1 = m.atomic_number(b->a1());
  const atomic_number_t z2 = m.atomic_number(b->a2());
  if (z2 == 7 && (z1 == 6 || z1 == 16)) {
    return IsAmide(m, b->a1());
  }
  if (z1 == 7 && (z2 == 6 || z2 == 16)) {
    return IsAmide(m, b->a2());
  }

  return 0;
}

int
QuickRotatableBonds::Process(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  const int matoms = m.natoms();
  if (matoms < 2) {
    return 0;
  }

  // Force sssr.
  m.ring_membership();

  switch (_calculation) {
    case RotBond::kUndefined:
      cerr << "QuickRotatableBonds::Process:no calculation defined\n";
      return 0;
    case RotBond::kQuick:
      return Quickest(m);
    case RotBond::kExpensive:
      return Expensive(m);
    default:
      cerr << "QuickRotatableBonds::Process:should not come here\n";
      return 0;
  }
}

std::unique_ptr<int[]>
AtomsWithTripleBonds(const Molecule& m) {
  std::unique_ptr<int[]> result(new_int(m.natoms()));

  const int nedges = m.nedges();
  for (int i = 0; i < nedges; ++i) {
    const Bond* b = m.bondi(i);
    if (b->is_triple_bond()) {
      result[b->a1()] = 1;
      result[b->a2()] = 1;
    }
  }

  return result;
}

int
QuickRotatableBonds::Expensive(Molecule& m) {
  const int matoms = m.natoms();

  resizable_array<const Bond*> candidate_bonds;
  candidate_bonds.resize(12);

  // For each atom, keep track of the number of terminal
  // atoms to which it is attached. This is later used for
  // CF3 and t-butyl
  std::unique_ptr<int[]> terminal_atom_count(new_int(matoms));
  std::unique_ptr<int[]> in_triple_bond = AtomsWithTripleBonds(m);

  const int nedges = m.nedges();
  for (int i = 0; i < nedges; ++i) {
    const Bond* b = m.bondi(i);
    if (! b->is_single_bond()) {
      continue;
    }
    if (b->nrings()) {
      continue;
    }
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (in_triple_bond[a1] || in_triple_bond[a2]) {
      continue;
    }

    if (IsAmide(m, b)) {
      continue;
    }

    const int a1con = m.ncon(b->a1());
    const int a2con = m.ncon(b->a2());
    if (a1con == 1) {
      ++terminal_atom_count[b->a2()];
      continue;
    } else if (a2con == 1) {
      ++terminal_atom_count[b->a1()];
      continue;
    }

    candidate_bonds << b;
  }

  if (candidate_bonds.empty()) {
    return 0;
  }

  // Look for any atoms attached to >=3 terinal atoms for CF3...
  int rc = 0;
  for (const Bond* b : candidate_bonds) {
    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();
    if (in_triple_bond[a1] || in_triple_bond[a2]) {
      continue;
    }

    // If either atom is attached to a terminal atom
    if (terminal_atom_count[a1] >= 3) {
      if (IsCf3(m, a1, a2)) {
        continue;
      }
    } else if (terminal_atom_count[a2] >= 3) {
      if (IsCf3(m, a2, a1)) {
        continue;
      }
    }

    ++rc;
    if (_isotope) {
      m.set_isotope(a1, _isotope);
      m.set_isotope(a2, _isotope);
    }
  }

  return rc;
}

int
QuickRotatableBonds::Quickest(Molecule& m) {
  int rc = 0;

  std::unique_ptr<int[]> in_triple_bond = AtomsWithTripleBonds(m);

  const int nedges = m.nedges();

  for (int i = 0; i < nedges; ++i) {
    const Bond* b = m.bondi(i);
    if (! b->is_single_bond()) {
      continue;
    }
    if (b->nrings()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (in_triple_bond[a1] || in_triple_bond[a2]) {
      continue;
    }

    if (m.ncon(a1) == 1 || m.ncon(a2) == 1) {
      continue;
    }

    ++rc;

    if (_isotope) {
      m.set_isotope(b->a1(), _isotope);
      m.set_isotope(b->a2(), _isotope);
    }
  }

  return rc;
}

}  // namespace quick_rotbond
