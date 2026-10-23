/*
  Rotatable bonds are computed in iwdescr and in the rotatable_bonds programme.
  In order to keep these in sync as much as possible, we include as many common
  functions here as possible
*/

#include <stdlib.h>

#include <algorithm>
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
  output << " -" << flag << " isotope=<iso>     isotope to be placed on all atoms at the end of a rotatable bond\n";
  output << " -" << flag << " uiso              each bond gets a unique isotope\n";
}

int
QuickRotatableBonds::Initialise(Command_Line& cl,
                                char flag) {
  const int verbose = cl.option_present('v');

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
    } else if (s.starts_with("isotope=")) {
      s.remove_leading_chars(8);
      if (! s.numeric_value(_isotope)) {
        cerr << "QuickRotatableBonds::Initialise:invalid isotope\n";
        return 0;
      }
      if (verbose) {
        cerr << "Isotope " << _isotope << " placed on atoms in rotatable bonds\n";
      }
    } else if (s == "uiso") {
      _unique_isotope_each_bond = 1;
      if (verbose) {
        cerr << "Will place a unique isotope on each rotatable bond.\n";
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

    // Explicit Hydrogens are invisible here. Including them would make the
    // attached atoms look inhomogeneous and hide a real CF3 or t-butyl.
    if (jz == 1) {
      continue;
    }

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
/*
  How many double bonds are there from `centre` to a Nitrogen, Oxygen or
  Sulphur?

  Ring double bonds do not count. Aromatic rings are held in a Kekule form, so
  an aromatic carbon really does have a double bond to an aromatic ring
  Nitrogen, and counting those would make the exocyclic bond of every
  aminopyridine look conjugated. The linkages this is used to recognise -
  carbonyl, thiocarbonyl, imine - all have the double bond outside any ring,
  including in a lactam, where the C=O is exocyclic even though the carbon is
  not.
*/

static int
DoublyBondedToHeteroatom(const Molecule& m,
        atom_number_t centre) {
  int rc = 0;

  const Atom& ac = m.atom(centre);
  for (const Bond * b : ac) {
    if (! b->is_double_bond()) {
      continue;
    }

    if (b->nrings()) {
      continue;
    }

    const atomic_number_t z = m.atomic_number(b->other(centre));
    if (z == 7 || z == 8 || z == 16) {
      ++rc;
    }
  }

  return rc;
}

/*
  Is `b` a conjugated linkage - a single bond that is not freely rotatable
  because one end is a carbonyl like centre and the other a heteroatom?

  Covers amides and thioamides, esters and thioesters, amidines, carbamates and
  sulfonamides. This is what RDKit's strict rotatable bond definition excludes
  with

    [CD3](=[N,O,S])-!@[#7,O,S!D1]

  with one deliberate difference. RDKit requires the doubly bonded centre to be
  carbon, so it counts a sulfonamide S-N bond as rotatable. That bond does have
  restricted rotation, so a Sulphur centre is accepted here.

  This was called IsAmide and required the partner to be Nitrogen, so ester,
  thioester and amidine linkages were treated as rotatable. Over
  data/Pubchem1000.smi that put the count above RDKit's on 17% of molecules.
  Agreement went from 754/1000 to 903/1000, and of the 97 that still differ, 84
  are the sulfonamides above.
*/

int
IsConjugatedLinkage(Molecule& m,
        const Bond* b) {
  const atom_number_t a1 = b->a1();
  const atom_number_t a2 = b->a2();
  const atomic_number_t z1 = m.atomic_number(a1);
  const atomic_number_t z2 = m.atomic_number(a2);

  // One end must be the heteroatom, the other the centre.
  if ((z2 == 7 || z2 == 8 || z2 == 16) && (z1 == 6 || z1 == 16)) {
    if (DoublyBondedToHeteroatom(m, a1)) {
      return 1;
    }
  }

  if ((z1 == 7 || z1 == 8 || z1 == 16) && (z2 == 6 || z2 == 16)) {
    if (DoublyBondedToHeteroatom(m, a2)) {
      return 1;
    }
  }

  return 0;
}

int
QuickRotatableBonds::Process(Molecule& m, int* bond_rotatable) {
  if (m.empty()) {
    return 0;
  }

  const int matoms = m.natoms();
  if (matoms < 2) {
    return 0;
  }

  if (bond_rotatable != nullptr) {
    std::fill_n(bond_rotatable, m.nedges(), 0);
    m.assign_bond_numbers_to_bonds_if_needed();
  }

  // Force sssr.
  m.ring_membership();

  switch (_calculation) {
    case RotBond::kUndefined:
      cerr << "QuickRotatableBonds::Process:no calculation defined\n";
      return 0;
    case RotBond::kQuick:
      return Quickest(m, bond_rotatable);
    case RotBond::kExpensive:
      return Expensive(m, bond_rotatable);
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

/*
  ncon() counts an explicit Hydrogen, and the terminal atom test in both
  calculations is a connectivity of one. So writing the Hydrogen of a hydroxyl or
  an amine explicitly would make that heteroatom non terminal and the bond to it
  rotatable - CCO has no rotatable bonds, CCO[H] appeared to have one. Both
  calculations therefore work with the heavy atom degree, and skip bonds to
  Hydrogen outright.

  Nothing here copies or modifies the molecule, which matters because the
  optional bond_rotatable array is indexed by the caller's bond numbering.
*/

static std::unique_ptr<int[]>
HeavyAtomDegree(const Molecule& m) {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> result(new_int(matoms));

  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m.atom(i);
    int heavy = 0;
    for (const Bond* b : a) {
      if (m.atomic_number(b->other(i)) != 1) {
        ++heavy;
      }
    }
    result[i] = heavy;
  }

  return result;
}

int
QuickRotatableBonds::Expensive(Molecule& m, int* bond_rotatable) {
  const int matoms = m.natoms();

  resizable_array<const Bond*> candidate_bonds;
  candidate_bonds.resize(12);

  isotope_t next_isotope = _isotope;

  // For each atom, keep track of the number of terminal
  // atoms to which it is attached. This is later used for
  // CF3 and t-butyl
  std::unique_ptr<int[]> terminal_atom_count(new_int(matoms));
  std::unique_ptr<int[]> in_triple_bond = AtomsWithTripleBonds(m);

  // See HeavyAtomDegree above. Bonds to Hydrogen are skipped entirely below so
  // they neither become candidates nor inflate terminal_atom_count.
  std::unique_ptr<int[]> heavy_ncon = HeavyAtomDegree(m);

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

    if (m.atomic_number(a1) == 1 || m.atomic_number(a2) == 1) {
      continue;
    }

    if (in_triple_bond[a1] || in_triple_bond[a2]) {
      continue;
    }

    if (IsConjugatedLinkage(m, b)) {
      continue;
    }

    const int a1con = heavy_ncon[a1];
    const int a2con = heavy_ncon[a2];
    if (a1con == 1) {
      ++terminal_atom_count[a2];
      continue;
    } else if (a2con == 1) {
      ++terminal_atom_count[a1];
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

    if (bond_rotatable != nullptr) {
      bond_rotatable[b->bond_number()] = 1;
    }

    if (_isotope == 0) {
      continue;
    }

    if (_unique_isotope_each_bond) {
      m.set_isotope(a1, next_isotope);
      m.set_isotope(a2, next_isotope);
      ++next_isotope;
    } else {
      m.set_isotope(a1, _isotope);
      m.set_isotope(a2, _isotope);
    }
  }

  return rc;
}

int
QuickRotatableBonds::Quickest(Molecule& m, int* bond_rotatable) {
  int rc = 0;

  std::unique_ptr<int[]> in_triple_bond = AtomsWithTripleBonds(m);

  // See HeavyAtomDegree above.
  std::unique_ptr<int[]> heavy_ncon = HeavyAtomDegree(m);

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

    if (m.atomic_number(a1) == 1 || m.atomic_number(a2) == 1) {
      continue;
    }

    if (in_triple_bond[a1] || in_triple_bond[a2]) {
      continue;
    }

    if (heavy_ncon[a1] == 1 || heavy_ncon[a2] == 1) {
      continue;
    }

    ++rc;

    if (bond_rotatable != nullptr) {
      bond_rotatable[b->bond_number()] = 1;
    }

    if (_isotope) {
      m.set_isotope(b->a1(), _isotope);
      m.set_isotope(b->a2(), _isotope);
    }
  }

  return rc;
}

int
RotatableBondsBetween(Molecule& m,
                      atom_number_t a1,
                      atom_number_t a2,
                      const int* bond_rotatable) {
  const int d = m.bonds_between(a1, a2);

  //cerr << "RotatableBondsBetween from " << a1 << " to " << a2 << " d " << d << '\n';

  for (const Bond* b : m[a1]) {
    const atom_number_t o = b->other(a1);
    if (m.bonds_between(o, a2) != d - 1) {
      continue;
    }

    int rc = 0;
    if (bond_rotatable[b->bond_number()]) {
      rc = 1;
    }

    if (o == a2) {
      return rc;
    }

    return rc + RotatableBondsBetween(m, o, a2, bond_rotatable);
  }

  return 0;
}

int
QuickRotatableBonds::RotatableBondsBetween(Molecule& m,
                atom_number_t a1,
                atom_number_t a2) {
  const int matoms = m.natoms();
  if (matoms < 2) {
    return 0;
  }

  if (m.fragment_membership(a1) != m.fragment_membership(a2)) {
    cerr << "QuickRotatableBonds::RotatableBondsBetween:inconsistent fragment membership\n";
    return 0;
  }

  std::unique_ptr<int[]> bond_rotatable = std::make_unique<int[]>(m.nedges());

  Process(m, bond_rotatable.get());

  return quick_rotbond::RotatableBondsBetween(m, a1, a2, bond_rotatable.get());
}

std::unique_ptr<int[]>
QuickRotatableBonds::RotatableBondsBetween(Molecule& m) {
  const int matoms = m.natoms();
  std::unique_ptr<int[]> result = std::make_unique<int[]>(matoms * matoms);

  if (matoms <= 1) {
    return result;
  }

  if (matoms == 2) {
    std::fill_n(result.get(), 0, 2*2);
    return result;
  }

  std::unique_ptr<int[]> bond_rotatable = std::make_unique<int[]>(m.nedges() * m.nedges());

  Process(m, bond_rotatable.get());

  for (int i = 0; i < matoms; ++i) {
    result[i * matoms + i] = 0;
    for (int j = i + 1; j < matoms; ++j) {
      int r = quick_rotbond::RotatableBondsBetween(m, i, j, bond_rotatable.get());
      result[i * matoms + j] = r;
      result[j * matoms + i] = r;
    }
  }

  return result;
}

}  // namespace quick_rotbond
