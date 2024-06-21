#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Tools/alogp.h"

namespace alogp {

using std::cerr;

struct PerMoleculeData {
  Molecule& mol;

  std::unique_ptr<atomic_number_t[]> z;
  int* unsaturation;
  int* aromatic;
  int* aryl_count;
  int* attached_heteroatom_count;
  int* single_bond_count;
  int* double_bond_count;
  int* triple_bond_count;
  formal_charge_t* formal_charge;
  int* assigned;
  // Useful for debugging.
  float* atom_value;

  PerMoleculeData(Molecule& m);
  ~PerMoleculeData();
};

PerMoleculeData::PerMoleculeData(Molecule& m) : mol(m) {
  const int matoms = m.natoms();

  z = m.AtomicNumbers();

  unsaturation = new_int(matoms);
  aromatic = new_int(matoms);
  aryl_count = new_int(matoms);
  attached_heteroatom_count = new_int(matoms);
  single_bond_count = new_int(matoms);
  double_bond_count = new_int(matoms);
  triple_bond_count = new_int(matoms);
  formal_charge = new formal_charge_t[matoms];
  std::fill_n(formal_charge, matoms, 0);

  assigned = new_int(matoms);
  atom_value = new float[matoms];
  std::fill_n(atom_value, matoms, 0.0f);

  for (const Bond* b : m.bond_list()) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (z[a1] != 6) {
      ++attached_heteroatom_count[a2];
    }
    if (z[a2] != 6) {
      ++attached_heteroatom_count[a1];
    }

    if (b->is_single_bond()) {
      ++single_bond_count[a1];
      ++single_bond_count[a2];
    } else if (b->is_double_bond()) {
      ++double_bond_count[a1];
      ++double_bond_count[a2];
      ++unsaturation[a1];
      ++unsaturation[a2];
    } else if (b->is_triple_bond()) {
      ++triple_bond_count[a1];
      ++triple_bond_count[a2];
      unsaturation[a1] += 2;
      unsaturation[a2] += 2;
    }
    
    if (b->is_aromatic()) {
      ++aromatic[a1];
      ++aromatic[a2];
    }
  }

  for (int i = 0; i < matoms; ++i) {
    formal_charge[i] = m.formal_charge(i);

    if (aromatic[i] == 0) {
      continue;
    }

    for (const Bond* b : m[i]) {
      atom_number_t o = b->other(i);
      ++aryl_count[o];
    }
  }
}

PerMoleculeData::~PerMoleculeData() {
  delete [] unsaturation;
  delete [] aromatic;
  delete [] aryl_count;
  delete [] attached_heteroatom_count;
  delete [] single_bond_count;
  delete [] double_bond_count;
  delete [] triple_bond_count;
  delete [] formal_charge;
  delete [] assigned;
  delete [] atom_value; 
}

// All error messages include the smarts of the problematic atom and
// the molecule name.
IWString
Diagnostic(PerMoleculeData& pmd, atom_number_t zatom) {
  IWString result;
  result << pmd.mol.smarts_equivalent_for_atom(zatom) << ' ' << pmd.mol.name();

  return result;
}

// Return the Bond that is not aromatic
const Bond*
NonAromaticConnection(Molecule& m, atom_number_t zatom) {
  assert(m.ncon(zatom) == 3);
  assert(m.is_aromatic(zatom));

  for (const Bond* b : m[zatom]) {
    if (b->is_aromatic()) {
      continue;
    }
    return b;
  }

  return nullptr;
}

int
ALogP::AromaticCarbon(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  const Atom& a = pmd.mol[zatom];
  const int acon = a.ncon();

  // C18 [cH]
  if (acon == 2) {
    result += 0.1581;
    pmd.atom_value[zatom] = 0.1581;
    pmd.assigned[zatom] = kC18;
    return 1;
  }

  // C19 aromatic bridghead [c](:a)(:a):a
  if (acon == 3 && pmd.aromatic[zatom] == 3) {
    result += 0.2955;
    pmd.atom_value[zatom] = 0.2955;
    pmd.assigned[zatom] = kC19;
    return 1;
  }

  assert(acon == 3);

  const Bond* b = NonAromaticConnection(pmd.mol, zatom);
  if (b == nullptr) {
    cerr << "AromaticCarbon:no outside ring bond? " << Diagnostic(pmd, zatom) << '\n';
    return 0;
  }

  const atom_number_t o = b->other(zatom);

  // C25
  if (b->is_double_bond() &&
      (pmd.z[o] == 6 || pmd.z[o] == 7 || pmd.z[o] == 8)) {
    result += -0.8186;
    pmd.atom_value[zatom] = -0.8186;
    pmd.assigned[zatom] = kC25;
    return 1;
  }

  // C14
  if (pmd.z[o] == 9) {
    // result += 0.0000;  zero result not added, no-op
    pmd.assigned[zatom] = kC14;
    // pmd.atom_value[zatom] = 0.0;
    return 1;
  }
  // C15
  if (pmd.z[o] == 17) {
    result += 0.2450;
    pmd.assigned[zatom] = kC15;
    pmd.atom_value[zatom] = 0.2450;
    return 1;
  }
  // C16
  if (pmd.z[o] == 35) {
    result += 0.1980;
    pmd.assigned[zatom] = kC16;
    pmd.atom_value[zatom] = 0.1980;
    return 1;
  }
  // C17
  if (pmd.z[o] == 53) {
    result += 0.000;
    pmd.assigned[zatom] = kC17;
    pmd.atom_value[zatom] = 0.00;
    return 1;
  }

  // C20 4-aromatic [c](:a)(:a)-a
  if (pmd.aromatic[o]) {
    result += 0.2713;
    pmd.assigned[zatom] = kC20;
    pmd.atom_value[zatom] = 0.2713;
    return 1;
  }
  // C21 [c](:a)(:a)-C
  if (pmd.z[o] == 6 && b->is_single_bond()) {
    result += 0.1360;
    pmd.assigned[zatom] = kC21;
    pmd.atom_value[zatom] = 0.1360;
    return 1;
  }
  // C22 [c](:a)(:a)-N
  if (pmd.z[o] == 7 && b->is_single_bond()) {
    result += 0.4619;
    pmd.assigned[zatom] = kC22;
    pmd.atom_value[zatom] = 0.4619;
    return 1;
  }
  // C23 [c](:a)(:a)-O
  if (pmd.z[o] == 8 && b->is_single_bond()) {
    result += 0.5437;
    pmd.assigned[zatom] = kC23;
    pmd.atom_value[zatom] = 0.5437;
    return 1;
  }
  // C24 [c](:a)(:a)-S
  if (pmd.z[o] == 16 && b->is_single_bond()) {
    result += 0.1893;
    pmd.assigned[zatom] = kC24;
    pmd.atom_value[zatom] = 0.1893;
    return 1;
  }
  // C13 aromatic heteratoms [cH0]-[!(C,N,O,S,F,Cl,Br,I)]’
  if (pmd.z[o] != 6 && b->is_single_bond()) {
    result +=-0.5443;
    pmd.assigned[zatom] = kC13;
    pmd.atom_value[zatom] = -0.5443;
    return 1;
  }

  result += 0.08129;
  pmd.assigned[zatom] = kCS;
  pmd.atom_value[zatom] = 0.08129;

  return 1;
}

int
ALogP::SaturatedPrimaryCarbon(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  assert(pmd.z[zatom] == 6);
  assert(pmd.mol.ncon(zatom) == 1);

  const atom_number_t o = pmd.mol.other(zatom, 0);

  if (pmd.aromatic[o]) {
    if (pmd.aromatic[o] && pmd.z[o] == 6) {
      result += 0.08452;
      pmd.assigned[zatom] = kC8;
      pmd.atom_value[zatom] = 0.08452;
      return kC8;
    } else {
      result += -0.1444;
      pmd.assigned[zatom] = kC9;
      pmd.atom_value[zatom] = -0.1444;
      return kC9;
    }
  } else {  // aliphatic
    // C1
    if (pmd.z[o] == 6) {
      result += 0.1441;
      pmd.assigned[zatom] = kC1;
      pmd.atom_value[zatom] = 0.1441;
      return kC1;
    } else {  // heteroatom
      // C3 primary heteroatom
      result += -0.2035;
      pmd.assigned[zatom] = kC3;
      pmd.atom_value[zatom] = -0.2035;
      return kC3;
    }
  }
}

int
ALogP::SaturatedSecondaryCarbom(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  assert(pmd.z[zatom] == 6);
  assert(pmd.mol.ncon(zatom) == 2);

  // C10 [CH2X4]a  secondary aromatic
  if (pmd.aryl_count[zatom]) {
    result += -0.0516;
    pmd.assigned[zatom] = kC10;
    pmd.atom_value[zatom] = -0.0516;
    return kC10;
  }

  // C1 secondary aliphatic
  if (pmd.attached_heteroatom_count[zatom] == 0) {
    result += 0.1441;
    pmd.assigned[zatom] = kC1;
    pmd.atom_value[zatom] = 0.1441;
    return 1;
  }


  // C3
  result += -0.2035;
  pmd.assigned[zatom] = kC3;
  pmd.atom_value[zatom] = -0.2035;
  return 1;
}

int
ALogP::SaturatedCarbon(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  const Atom& a = pmd.mol[zatom];
  const int acon = a.ncon();

#ifdef NOW_BEING_HANDLED_EXERNALLY
  // C1
  if (a.ncon() == 0) {
    result += 0.1441;
    pmd.assigned[zatom] = kC1;
    return 1;
  }
#endif

  const int ahc = pmd.attached_heteroatom_count[zatom];

  if (acon == 1) {
    return SaturatedPrimaryCarbon(pmd, zatom, result);
  }

  if (acon == 2) {
    return SaturatedSecondaryCarbom(pmd, zatom, result);
  }

  // C11
  if (acon == 3 && pmd.aryl_count[zatom]) {
    result += 0.1193;
    pmd.assigned[zatom] = kC11;
    pmd.atom_value[zatom] = 0.1193;
    return kC11;
  }

  // C12
  if (acon == 4 && pmd.aryl_count[zatom]) {
    result += -0.0967;
    pmd.assigned[zatom] = kC12;
    pmd.atom_value[zatom] = -0.0967;
    return 1;
  }

  // C2
  if (ahc == 0) {
    // result += 0.0; zero contribution, no-op
    pmd.assigned[zatom] = kC2;
    // pmd.atom_value[zatom] = 0.000;
    return 1;
  }
  // C4
  result += -0.2051;
  pmd.assigned[zatom] = kC4;
  pmd.atom_value[zatom] = -0.2051;
  return 1;
}

// Return the first atom number that is doubly bonded to `zatom`.
atom_number_t
DoublyBondedTo(PerMoleculeData& pmd, atom_number_t zatom) {
  const Atom& a = pmd.mol[zatom];
  for (const Bond* b : a) {
    if (b->is_double_bond()) {
      return b->other(zatom);
    }
  }

  return kInvalidAtomNumber;
}

int
ALogP::UnSaturatedCarbon(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  // C7 [CX2]#A acetylene, nitrile
  if (pmd.triple_bond_count[zatom] && pmd.formal_charge[zatom] == 0) {
    result += 0.00170;
    pmd.assigned[zatom] = kC7;
    pmd.atom_value[zatom] = 0.00170;
    return 1;
  }

  // CS [C-]#N
  if (pmd.formal_charge[zatom] == -1 && pmd.triple_bond_count[zatom]) {
    result += 0.08129;
    pmd.assigned[zatom] = kCS;
    pmd.atom_value[zatom] = 0.08129;
    return 1;
  }

  assert(pmd.double_bond_count[zatom] > 0);

  atom_number_t dbl = DoublyBondedTo(pmd, zatom);
  if (dbl == kInvalidAtomNumber) {
    cerr << "HUH, no doubly bond!! " << Diagnostic(pmd, zatom) << '\n';
    return 0;
  }

  // C5 [C]=[A#X] C = Heteratom
  if (pmd.z[dbl] != 6) {
    result += -0.2783;
    pmd.assigned[zatom] = kC5;
    pmd.atom_value[zatom] = -0.2783;
    return 1;
  }

  // At this stage the doubly bonded atom is definitely a carbon.

  // C6 C=C aliphatic
  if (pmd.aryl_count[zatom] == 0) {
    result += 0.1551;
    pmd.assigned[zatom] = kC6;
    pmd.atom_value[zatom] = 0.1551;
    return 1;
  }

  // C26 C=C aromatic [C]()C)(a)A’, ‘[C]()C)(c)a’, ‘[CH]() C)a’, ‘[C] ) c
  result += 0.2640;
  pmd.assigned[zatom] = kC26;
  pmd.atom_value[zatom] = 0.2640;
  return 1;
}


int
ALogP::Carbon(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  if (pmd.aromatic[zatom]) {
    return AromaticCarbon(pmd, zatom, result);
  }

  if (pmd.unsaturation[zatom] == 0) {
    return SaturatedCarbon(pmd, zatom, result);
  }

  return UnSaturatedCarbon(pmd, zatom, result);
}

// The number of aromatic connections to `zatom`.
int
ArylCount(PerMoleculeData& pmd, atom_number_t zatom) {
  int rc = 0;
  for (const Bond* b : pmd.mol[zatom]) {
    atom_number_t o = b->other(zatom);
    if (pmd.aromatic[o]) {
      ++rc;
    }
  }

  return rc;
}

int
ALogP::AromaticNitrogen(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  // N11 Unprotonated aromatic
  const Atom& a = pmd.mol[zatom];

  if (a.ncon() == 2 && pmd.formal_charge[zatom] == 0) {
    result +=  -0.3239;
    pmd.assigned[zatom] = kN11;
    pmd.atom_value[zatom] = -0.3239;
    return 1;
  }

  // N1(=C2C(=N(=O)C=C1)C=CC=C2)=O CHEMBL2104626
  // which can be written in charge separate form. RDKit
  // assigns this as if it is the charge separated form.
  // Probably that is what is intended in the paper.
  // N12
  if (a.ncon() == 3 && pmd.formal_charge[zatom] == 0 &&
      pmd.double_bond_count[zatom] == 2) {
    result += -1.119;
    pmd.assigned[zatom] = kN12;
    pmd.atom_value[zatom] = -1.119;
    return 1;
  }

  if (pmd.formal_charge[zatom] == 0) {
    result += -0.3239;
    pmd.assigned[zatom] = kN11;
    pmd.atom_value[zatom] = -0.3239;
    return 1;
  }

  // N12 protonated aromatic
  result += -1.119;
  pmd.assigned[zatom] = kN12;
  pmd.atom_value[zatom] = -1.119;
  return 1;
}

// The NH2 atom in an amidine
int
SaturatedNitrogenIsAmidine(PerMoleculeData& pmd, atom_number_t zatom) {
  atom_number_t carbon = pmd.mol.other(zatom, 0);

  if (pmd.attached_heteroatom_count[carbon] != 2) {
    return 0;
  }
  if (pmd.aromatic[carbon]) {
    return 0;
  }
  if (pmd.double_bond_count[carbon] != 1) {
    return 0;
  }
  if (pmd.single_bond_count[carbon] != 2) {
    return 0;
  }

  for (const Bond* b : pmd.mol[carbon]) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t o = b->other(carbon);
    // Should we also check that `o` is singly bonded.
    if (pmd.z[o] == 7) {
      return 1;
    }
  }

  return 0;
}

// The NH2 atom in a guanidine?
int
SaturatedNitrogenIsGuanidine(PerMoleculeData& pmd, atom_number_t zatom) {
  atom_number_t carbon = pmd.mol.other(zatom, 0);

  if (pmd.attached_heteroatom_count[carbon] != 3) {
    return 0;
  }
  if (pmd.aromatic[carbon]) {
    return 0;
  }
  if (pmd.double_bond_count[carbon] != 1) {
    return 0;
  }
  if (pmd.single_bond_count[carbon] != 2) {
    return 0;
  }

  for (const Bond* b : pmd.mol[zatom]) {
    atom_number_t o = b->other(carbon);
    if (o == zatom) {
      continue;
    }
    if (pmd.z[o] != 7) {
      return 0;
    }
  }

  return 1;
}

// Return the number of adjacent unsaturations to `zatom` -
// excludes aromatics.
int
Vinyl(const PerMoleculeData& pmd, atom_number_t zatom) {
  int rc = 0;
  for (const Bond* b : pmd.mol[zatom]) {
    const atom_number_t o = b->other(zatom);
    if (pmd.aromatic[o]) {
      continue;
    }
    if (pmd.unsaturation[o]) {
      ++rc;
    }
  }

  return rc;
}

// Is a terminal nitrogen part of an amide. Allow sulfonamides
int
TerminalNitrogenIsAmide(const PerMoleculeData& pmd, atom_number_t zatom) {
  assert(pmd.mol.ncon(zatom) == 1);

  // Allow either carbon or Sulphur
  const atom_number_t carbon = pmd.mol.other(zatom, 0);
  if (pmd.z[carbon] == 6) {
  } else if (pmd.z[carbon] == 16) {
  } else {
    return 0;
  }

  int doubly_bonded_oxygen_count = 0;
  for (const Bond* b : pmd.mol[carbon]) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t o = b->other(carbon);
    if (pmd.z[o] == 8) {
    } else if (pmd.z[o] == 16) {
    } else {
      continue;
    }

    ++doubly_bonded_oxygen_count;
  }

  return doubly_bonded_oxygen_count;
}

// Is zatom the N of a sulfonamide
int
IsSulfonamide(const PerMoleculeData& pmd, atom_number_t zatom) {
  for (const Bond* b : pmd.mol[zatom]) {
    const atom_number_t sulphur = b->other(zatom);
    if (pmd.z[sulphur] != 16) {
      continue;
    }

    if (pmd.unsaturation[sulphur] == 0) {
      return 0;
    }

    if (pmd.attached_heteroatom_count[sulphur] < 2) {
      return 0;
    }

    for (const Bond* b2 : pmd.mol[sulphur]) {
      if (! b2->is_double_bond()) {
        continue;
      }
      atom_number_t o = b2->other(sulphur);
      if (pmd.z[o] == 8) {
        return 1;
      }
    }
  }

  return 0;
}

int
ALogP::SinglyConnectedSaturatedNitrogen(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  assert(pmd.mol.ncon(zatom) == 1);

  if (pmd.attached_heteroatom_count[zatom] == 0 &&
      SaturatedNitrogenIsGuanidine(pmd, zatom)) {
    result += -1.0190;
    pmd.assigned[zatom] = kN1;
    pmd.atom_value[zatom] = -1.0190;
    return 1;
  }

  if (pmd.attached_heteroatom_count[zatom] == 0 &&
      SaturatedNitrogenIsAmidine(pmd, zatom)) {
    result += -1.0190;
    pmd.assigned[zatom] = kN1;
    pmd.atom_value[zatom] = -1.0190;
    return 1;
  }

  const int vinyl = Vinyl(pmd, zatom);

  // We are not using a charge assigner, so we approximate a primary amine this way.
  if (pmd.attached_heteroatom_count[zatom] == 0 &&
      vinyl == 0) {
    atom_number_t o = pmd.mol.other(zatom, 0);
    if (pmd.aromatic[o]) {  // N3
      result += -1.0270;
      pmd.assigned[zatom] = kN3;
      pmd.atom_value[zatom] = -1.0270;
      return 1;
    } else {  // N1
      result += -1.0190;
      pmd.assigned[zatom] = kN1;
      pmd.atom_value[zatom] = -1.0190;
      return 1;
    }
  }

  // but notice that RDKit assigns N1 to a terminal amide
  if (pmd.attached_heteroatom_count[zatom] == 0 &&
      vinyl && TerminalNitrogenIsAmide(pmd, zatom)) {
    result += -1.0190;
    pmd.assigned[zatom] = kN1;
    pmd.atom_value[zatom] = -1.0190;
    return 1;
  }

  if (pmd.attached_heteroatom_count[zatom] == 1 &&
      vinyl && IsSulfonamide(pmd, zatom)) {
    result += -1.0190;
    pmd.assigned[zatom] = kN1;
    pmd.atom_value[zatom] = -1.0190;
    return 1;
  }

  // N3 maybe not charged, but still adjacent to arom. Aniline
  if (pmd.aryl_count[zatom] == 1) {
    result += -1.0270;
    pmd.assigned[zatom] = kN3;
    pmd.atom_value[zatom] = -1.0270;
    return 1;
  }

  // NS nitrogen supplemental.
  result += -0.4806;
  pmd.assigned[zatom] = kNS;
  pmd.atom_value[zatom] = -0.4806;
  return 1;
}

int
ALogP::SaturatedNitrogen(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  const Atom& a = pmd.mol[zatom];
  // cerr << "SaturatedNitrogen " << pmd.mol.aromatic_smiles() << '\n';

  if (a.ncon() == 1) {
    return SinglyConnectedSaturatedNitrogen(pmd, zatom, result);
  }

  const int vinyl = Vinyl(pmd, zatom);
  const int aryl = ArylCount(pmd, zatom);

  if (a.ncon() == 2 && pmd.attached_heteroatom_count[zatom] == 0 &&
     Vinyl(pmd, zatom) == 0) {
    // N4
    if (aryl) {
      result += -0.5188;
      pmd.assigned[zatom] = kN4;
      pmd.atom_value[zatom] = -0.5188;
      return 1;
    } else {  // N2
      result += -0.7096;
      pmd.assigned[zatom] = kN2;
      pmd.atom_value[zatom] = -0.7096;
      return 1;
    }
  }

  // N4 - not charged
  if (a.ncon() == 2 && pmd.aryl_count[zatom]) {
    result += -0.5188;
    pmd.assigned[zatom] = kN4;
    pmd.atom_value[zatom] = -0.5188;
    return 1;
  }

  // N2 - not charged
  if (a.ncon() == 2 && aryl == 0) {
    result += -0.7096;
    pmd.assigned[zatom] = kN2;
    pmd.atom_value[zatom] = -0.7096;
    return 1;
  }

  if (a.ncon() == 3 && pmd.attached_heteroatom_count[zatom] == 0 &&
      Vinyl(pmd, zatom) == 0) {
    // N8
    if (aryl) {
      result += -0.4458;
      pmd.assigned[zatom] = kN8;
      pmd.atom_value[zatom] = -0.4458;
      return 1;
    } else {  // N7
      result += -0.3187;
      pmd.assigned[zatom] = kN7;
      pmd.atom_value[zatom] = -0.3187;
      return 1;
    }
  }

  // N8 not charged
  if (a.ncon() == 3 && aryl) {
    result += -0.4458;
    pmd.assigned[zatom] = kN8;
    pmd.atom_value[zatom] = -0.4458;
    return 1;
  }

  // N7 but with vinyl connections.
  if (a.ncon() == 3 && pmd.attached_heteroatom_count[zatom] == 0 &&
      pmd.aryl_count[zatom] == 0) {
    result += -0.3187;
    pmd.assigned[zatom] = kN7;
    pmd.atom_value[zatom] = -0.3187;
    return 1;
  }

  // N7 that might be part of a sulfonamide
  if (a.ncon() == 3 && pmd.attached_heteroatom_count[zatom] == 1 &&
      vinyl == 1 && IsSulfonamide(pmd, zatom)) {
    result += -0.3187;
    pmd.assigned[zatom] = kN7;
    pmd.atom_value[zatom] = -0.3187;
    return 1;
  }

  // N7 with very few restrictions
  if (a.ncon() == 3 && pmd.attached_heteroatom_count[zatom] == 1 &&
      vinyl == 0) {
    result += -0.3187;
    pmd.assigned[zatom] = kN7;
    pmd.atom_value[zatom] = -0.3187;
    return 1;
  }

  // N7 with hardly any restriction. Note that these would probably
  // not be charged.
  if (a.ncon() == 3 && pmd.attached_heteroatom_count[zatom] == 1) {
    result += -0.3187;
    pmd.assigned[zatom] = kN7;
    pmd.atom_value[zatom] = -0.3187;
    return 1;
  }

  // N7 with hardly any restriction. Note that these would probably
  // not be charged.
  if (a.ncon() == 3 && pmd.attached_heteroatom_count[zatom] == 2) {
    result += -0.3187;
    pmd.assigned[zatom] = kN7;
    pmd.atom_value[zatom] = -0.3187;
    return 1;
  }

  if (a.ncon() == 4 && pmd.formal_charge[zatom] == 1) {
    // N13
    result += -0.3396;
    pmd.assigned[zatom] = kN13;
    pmd.atom_value[zatom] = -0.3396;
    return 1;
  }

  // N14
  if (pmd.formal_charge[zatom]) {
    result += 0.2887;
    pmd.assigned[zatom] = kN14;
    pmd.atom_value[zatom] = 0.2887;
    return 1;
  }

  // N14
  if (a.ncon() == 2 && pmd.unsaturation[zatom] == 3) {
    result += 0.2887;
    pmd.assigned[zatom] = kN14;
    pmd.atom_value[zatom] = 0.2887;
    return 1;
  }

  // N14
  if (pmd.unsaturation[zatom] == 2 && a.ncon() == 1 &&
      pmd.attached_heteroatom_count[zatom] == 1) {
    result += 0.2887;
    pmd.assigned[zatom] = kN14;
    pmd.atom_value[zatom] = 0.2887;
    return 1;
  }

  // NS nitrogen supplemental.
  result += -0.4806;
  pmd.assigned[zatom] = kNS;
  pmd.atom_value[zatom] = -0.4806;
  return 1;
}

int
IsNitro(const PerMoleculeData& pmd, atom_number_t zatom) {
  int doubly_bonded_oxygen_count = 0;
  for (const Bond * b : pmd.mol[zatom]) {
    if (! b->is_double_bond()) {
      continue;
    }
    atom_number_t o = b->other(zatom);
    if (pmd.z[o] == 8) {
      ++doubly_bonded_oxygen_count;
    }
  }

  return doubly_bonded_oxygen_count == 2;
}

int
ALogP::UnSaturatedNitrogen(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  const Atom& a = pmd.mol[zatom];

  // N9 Nitrile
  if (a.ncon() == 1 && pmd.triple_bond_count[zatom] &&
      pmd.attached_heteroatom_count[zatom] == 0) {
    result += 0.01508;
    pmd.assigned[zatom] = kN9;
    pmd.atom_value[zatom] = 0.01508;
    return 1;
  }

  // N14 in [C-]#[N+]
  if (a.ncon() == 2 && pmd.formal_charge[zatom] == 1 &&
      pmd.triple_bond_count[zatom] == 1) {
    result += 0.2887;
    pmd.assigned[zatom] = kN14;
    pmd.atom_value[zatom] = 0.2887;
    return 1;
  }


  // Nitro is assigned N13
  if (a.ncon() == 3 && pmd.double_bond_count[zatom] == 2 &&
      pmd.attached_heteroatom_count[zatom] >= 2 && IsNitro(pmd, zatom)) {
    result += -0.3396;
    pmd.assigned[zatom] = kN13;
    pmd.atom_value[zatom] = -0.3396;
    return 1;
  }

  // N7 imine
  if (a.ncon() == 1 && pmd.double_bond_count[zatom] == 1) {
    result += 0.0837;
    pmd.assigned[zatom] = kN7;
    pmd.atom_value[zatom] = 0.0837;
    return 1;
  }

  // N6 substituted imine
  if (a.ncon() == 2 && pmd.double_bond_count[zatom] == 1 && pmd.triple_bond_count[zatom] == 0) {
    result += 0.1836;
    pmd.assigned[zatom] = kN6;
    pmd.atom_value[zatom] = 0.1836;
    return 1;
  }

  // N14 N#N=
  if (a.ncon() == 1 && pmd.triple_bond_count[zatom] == 1 &&
      pmd.attached_heteroatom_count[zatom]) {
    result += 0.2887;
    pmd.assigned[zatom] = kN14;
    pmd.atom_value[zatom] = 0.2887;
    return 1;
  }

  // N14 N#N=
  if (a.ncon() == 2 && pmd.triple_bond_count[zatom] == 1 &&
      pmd.double_bond_count[zatom] == 1 &&
      pmd.attached_heteroatom_count[zatom]) {
    result += 0.2887;
    pmd.assigned[zatom] = kN14;
    pmd.atom_value[zatom] = 0.2887;
    return 1;
  }

  // NS nitrogen supplemental.
  result += -0.4806;
  pmd.assigned[zatom] = kNS;
  pmd.atom_value[zatom] = -0.4806;
  return 1;
}

int
ALogP::Nitrogen(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  if (pmd.aromatic[zatom]) {
    return AromaticNitrogen(pmd, zatom, result);
  }

  if (pmd.unsaturation[zatom] == 0) {
    return SaturatedNitrogen(pmd, zatom, result);
  }

  return UnSaturatedNitrogen(pmd, zatom, result);
}

int
ALogP::UnSaturatedOxygen(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  assert(pmd.double_bond_count[zatom] == 1);

  const atom_number_t o = pmd.mol.other(zatom, 0);
  // O5
  if (pmd.z[o] == 7) {
    result += 0.0335;
    pmd.assigned[zatom] = kO5;
    pmd.atom_value[zatom] = 0.0335;
    return 1;
  }

  // O8
  if (pmd.z[o] == 6 && pmd.aromatic[o]) {
    result += 0.1788;
    pmd.assigned[zatom] = kO8;
    pmd.atom_value[zatom] = 0.1788;
    return 1;
  }

  // IAW. Copy behaviour from RDKit which assigns O6 to the oxygens in O=S=O
  // which seems quite wrong.
  if (pmd.z[o] == 16) {
    result += -0.3339;
    pmd.assigned[zatom] = kO6;
    pmd.atom_value[zatom] = -0.3339;
    return 1;
  }

  // O11 carbonyl heteroatom
  if (pmd.z[o] == 6 && pmd.attached_heteroatom_count[o] == 3) {
    result += 0.4833;
    pmd.assigned[zatom] = kO11;
    pmd.atom_value[zatom] = 0.4833;
    return 1;
  }

  // O10 carbonyl aromatic
  if (pmd.attached_heteroatom_count[zatom] == 0 && ArylCount(pmd, o)) {
    result += 0.1129;
    pmd.assigned[zatom] = kO10;
    pmd.atom_value[zatom] = 0.1129;
    return 1;
  }

  // O9 carbonyl aliphatic
  if (pmd.z[o] == 6) {
    result += -0.1526;
    pmd.assigned[zatom] = kO9;
    pmd.atom_value[zatom] = -0.1526;
    return 1;
  }

  // OS
  result += -0.1188;
  pmd.assigned[zatom] = kOS;
  pmd.atom_value[zatom] = -0.1188;
  return 1;
}

int
IsAcid(PerMoleculeData& pmd, atom_number_t zatom) {
  assert (pmd.z[zatom] == 8);

  const atom_number_t carbon = pmd.mol.other(zatom, 0);
  for (const Bond* b : pmd.mol[carbon]) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t o = b->other(carbon);
    if (pmd.z[o] == 8) {
      return 1;
    }
  }

  return 0;
}

int
AttachedAromaticCount(PerMoleculeData& pmd, atom_number_t zatom) {
  int rc = 0;
  for (const Bond* b : pmd.mol[zatom]) {
    const atom_number_t o = b->other(zatom);
    if (pmd.aromatic[o]) {
      ++rc;
    }
  }

  return rc;
}

// Return true if zatom is attached to a fully saturated Nitrogen
int
IsNHydroxy(const PerMoleculeData& pmd, atom_number_t zatom) {
  assert(pmd.mol.ncon(zatom) == 1);
  atom_number_t o = pmd.mol.other(zatom, 0);
  if (pmd.z[o] != 7) {
    return 0;
  }
  if (pmd.unsaturation[o]) {
    return 0;
  }

  return 1;
}

int
ALogP::SaturatedOxygen(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  const Atom& a = pmd.mol[zatom];
  const int acon = a.ncon();

  // cerr << "ALogP::SaturatedOxygen acid " << IsAcid(pmd, zatom) << '\n';
  // O12
  if (acon == 1 && IsAcid(pmd, zatom)) {
    if (_use_alcohol_for_acid) {
      result += -0.2893;
      pmd.assigned[zatom] = kO2;
      pmd.atom_value[zatom] = -0.2893;
    } else {
      result += -1.326;
      pmd.assigned[zatom] = kO12;
      pmd.atom_value[zatom] = -1.326;
    }
    return 1;
  }

  // O2 alcohol
  if (acon == 1 && pmd.attached_heteroatom_count[zatom] == 0) {
    result += -0.2893;
    pmd.assigned[zatom] = kO2;
    pmd.atom_value[zatom] = -0.2893;
    return 1;
  }

  // O2 N-Hydroxy
  if (acon == 1 && pmd.attached_heteroatom_count[zatom] == 1 &&
     pmd.aryl_count[zatom] == 0 && IsNHydroxy(pmd, zatom)) {
    result += -0.2893;
    pmd.assigned[zatom] = kO2;
    pmd.atom_value[zatom] = -0.2893;
    return 1;
  }

  int aac = AttachedAromaticCount(pmd, zatom);
  // O3 Aliphatic ether.
  if (acon == 2 && aac == 0) {
    result += -0.0684;
    pmd.assigned[zatom] = kO3;
    pmd.atom_value[zatom] = -0.0684;
    return 1;
  }

  // O4 Aromatic ether
  if (acon == 2 && aac) {
    result += -0.4195;
    pmd.assigned[zatom] = kO4;
    pmd.atom_value[zatom] = -0.4195;
    return 1;
  }

  if (acon == 1) {
    atom_number_t o = pmd.mol.other(zatom, 0);
    // O6 oxide
    if (pmd.z[o] == 16) {
      result += -0.3339;
      pmd.assigned[zatom] = kO6;
    pmd.atom_value[zatom] = -0.3339;
      return 1;
    }

    // N-oxide Note that this is grouped with other queries in the paper.
    // O5
    if (pmd.z[o] == 7) {
      result += 0.0335;
      pmd.assigned[zatom] = kO5;
      pmd.atom_value[zatom] = 0.0335;
      return 1;
    }

    // O6
    if (pmd.z[o] == 16) {
      result += -0.3339;
      pmd.assigned[zatom] = kO6;
      pmd.atom_value[zatom] = 0.3339;
      return 1;
    }
  }

  // OS oxygen supplemental.
  result += -0.1188;
  pmd.assigned[zatom] = kOS;
  pmd.atom_value[zatom] = -0.1188;
  return 1;
}


int
ALogP::Oxygen(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  // O1 aromatic
  if (pmd.aromatic[zatom]) {
    result += 0.1552;
    pmd.assigned[zatom] = kO1;
    pmd.atom_value[zatom] = 0.1552;
    return 1;
  }

  if (pmd.double_bond_count[zatom]) {
    return UnSaturatedOxygen(pmd, zatom, result);
  }

  return SaturatedOxygen(pmd, zatom, result);
}

int
ALogP::Fluorine(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  result += 0.4202;
  pmd.assigned[zatom] = kF;
  pmd.atom_value[zatom] = 0.4202;
  return 1;
}

int
ALogP::Chlorine(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  result += 0.6895;
  pmd.assigned[zatom] = kCl;
  pmd.atom_value[zatom] = 0.6895;
  return 1;
}

int
ALogP::Bromine(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  result += 0.8456;
  pmd.assigned[zatom] = kBr;
  pmd.atom_value[zatom] = 0.8456;
  return 1;
}

int
ALogP::Iodine(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  result += 0.8857;
  pmd.assigned[zatom] = kI;
  pmd.atom_value[zatom] = 0.8857;
  return 1;
}

int
ALogP::Phosphorus(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  result += 0.8612;
  pmd.assigned[zatom] = kP;
  pmd.atom_value[zatom] = 0.8612;
  return 1;
}

// Or nitrogen to deal with RDKit's handling of
// S(=O)(=O)(N=S(C)CCCC)C1=CC=C(C)C=C1 CHEMBL1533581
int
DoublyBondedToOxygen(const PerMoleculeData& pmd, atom_number_t zatom) {
  for (const Bond* b : pmd.mol[zatom]) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t o = b->other(zatom);
    if (pmd.z[o] == 8) {
      return 1;
    }
    if (pmd.z[o] == 7) {
      return 1;
    }
  }

  return 0;
}

int
ALogP::Sulphur(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  // S3
  if (pmd.aromatic[zatom]) {
    result += 0.6237;
    pmd.assigned[zatom] = kS3;
    pmd.atom_value[zatom] = 0.6237;
    return 1;
  }

  // RDKit seems to apply S2 to S=O groups, which seems wrong
  // Maybe this should be made a settable behaviour.
  if (pmd.double_bond_count[zatom] && DoublyBondedToOxygen(pmd.mol, zatom)) {
    result += -0.0024;
    pmd.assigned[zatom] = kS2;
    pmd.atom_value[zatom] = -0.0024;
    return 1;
  }

  // Not sure this can happen...
  // S2
  if (pmd.formal_charge[zatom]) {
    result += -0.0024;
    pmd.assigned[zatom] = kS2;
    pmd.atom_value[zatom] = -0.0024;
    return 1;
  }

  // S1
  result += 0.6482;
  pmd.assigned[zatom] = kS1;
  pmd.atom_value[zatom] = 0.6482;
  return 1;
}


// [#1]OC=[#6]’, ‘[#1]OC=[#7]’, ‘[#1]OC=O’, ‘[#1]OC=S 
// Also picks up phosphoric acids.
int
ALogP::IsHydrogenAcid(PerMoleculeData& pmd, atom_number_t zatom) {
  if (pmd.z[zatom] != 8) {
    return 0;
  }

  if (pmd.aromatic[zatom]) {
    return 0;
  }

  const atom_number_t carbon = pmd.mol.other(zatom, 0);
  if (pmd.aromatic[carbon]) {
    return 0;
  }

  // RDKit does not consider phosphoric acid H's to be acidic.
  // An argument could be made to only consider 1 of the two
  // to be acidic, but that would slow things down by breaking
  // the paradigm of each atom being processed independently.
  if (pmd.z[carbon] == 15 && _rdkit_phoshoric_acid_hydrogen) {
    return 0;
  }

  for (const Bond* b : pmd.mol[carbon]) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t o = b->other(carbon);
    if (pmd.z[o] == 6 || pmd.z[o] == 7 || pmd.z[o] == 8 || pmd.z[o] == 16) {
      return 1;
    }
  }

  return 0;
}

// ‘[#1]OO’, ‘[#1]OS’
int
IsPeroxide(PerMoleculeData& pmd, atom_number_t zatom) {
  if(pmd.z[zatom] != 8) {
    return 0;
  }

  const atom_number_t os = pmd.mol.other(zatom, 0);
  if (pmd.mol.ncon(os) != 2) {
    return 0;
  }
  if (pmd.z[os] == 8 || pmd.z[os] == 16) {
    return 1;
  }

  return 0;
}

// [#1][#7]’, ‘[#1]O[#7]
int
IsHydrogenAmine(PerMoleculeData& pmd, atom_number_t zatom) {
  if (pmd.z[zatom] == 7) {
    return 1;
  }

  if (pmd.z[zatom] != 8) {
    return 0;
  }

  if (pmd.attached_heteroatom_count[zatom] != 1) {
    return 0;
  }

  atom_number_t o = pmd.mol.other(zatom, 0);
  return pmd.z[o] == 7;
}

// [#1]O[CX4]’, ‘[#1]Oc’, ‘[#1]O[#1]’, ‘[#1]O[#5]’, ‘[#1]O[#14]’, ‘[#1]O[#15]’, ‘[#1]O[#33]’,
// ‘[#1]O[#50]’, ‘[#1][#5]’, ‘[#1][#14]’, ‘[#1][#15]’, ‘[#1][#16]’, ‘[#1][#50]
int
IsHydrogenAlcohol1(PerMoleculeData& pmd, atom_number_t zatom) {
  assert(pmd.z[zatom] == 8);

  const atom_number_t o = pmd.mol.other(zatom, 0);

  // [#1]O[CX4]
  if (pmd.z[o] == 6 && pmd.unsaturation[o] == 0) {
    return 1;
  }

  // [#1][#16]
  if (pmd.z[o] == 16) {
    return 1;
  }

  // [#1][#15]
  if (pmd.z[o] == 15) {
    return 1;
  }

  // [#1]Oc
  if (pmd.z[o] == 6 && pmd.aromatic[o]) {
    return 1;
  }

  // [#1]O[!(C,N,O,S)]
  if (pmd.z[o] == 6) {
  } else if (pmd.z[o] ==7) {
  } else if (pmd.z[o] == 8) {
  } else if (pmd.z[o] == 16) {
  } else {
    return 1;
  }

  return 0;
}

// [#1]O[CX4]’, ‘[#1]Oc’, ‘[#1]O[!(C,N,O,S)]’, ‘[#1][!C,N,O)]
int
IsHydrogenAlcohol(PerMoleculeData& pmd, atom_number_t zatom) {
  if (pmd.z[zatom] == 16) {
    return 1;
  }
  if (pmd.z[zatom] == 15) {
    return 1;
  }

  if (pmd.z[zatom] == 8) {
    return IsHydrogenAlcohol1(pmd, zatom);
  }

  // The smarts !C,N,O really does not make sense, N,O does the same thing.
  if (pmd.z[zatom] == 6) {
    return 0;
  }

  if (pmd.z[zatom] == 7 || pmd.z[zatom] == 8) {
    return 1;
  }

  return 0;
}

int
IsHydroCarbon(PerMoleculeData& pmd, atom_number_t zatom) {
  return pmd.z[zatom] == 6;
}

// #define DEBUG_ADD_HYDROGEN
int
ALogP::AddHydrogenContributions(PerMoleculeData& pmd, float& result) {
  Molecule& m = pmd.mol;
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    int h = m.hcount(i);
    if (_use_alcohol_for_acid && pmd.formal_charge[i] == -1 &&
        pmd.assigned[i] == kO2) {
      h = 1;
    }
    if (h > 0 && _rdkit_charged_nitrogen && pmd.formal_charge[i] == 1 &&
        (pmd.assigned[i] == kN1 || pmd.assigned[i] == kN2 ||
         pmd.assigned[i] == kN7)) {
      --h;
    }

    if (h == 0) {
      continue;
    }
#ifdef DEBUG_ADD_HYDROGEN
    cerr << "atom " << i << " result so far " << result << '\n';
#endif

    // H4
    if (IsHydrogenAcid(pmd, i)) {
      result += h * 0.2980;
#ifdef DEBUG_ADD_HYDROGEN
      cerr << "IsHydrogenAcid\n";
#endif
    } else if (IsPeroxide(pmd, i)) {   // H4
      result += h * 0.2980;
#ifdef DEBUG_ADD_HYDROGEN
      cerr << "IsPeroxide\n";
#endif
    } else if (IsHydrogenAmine(pmd, i)) {  // H3
      result += h * 0.2142;
#ifdef DEBUG_ADD_HYDROGEN
      cerr << "IsHydrogenAmine\n";
#endif
    } else if (IsHydrogenAlcohol(pmd, i)) {  // H2
      result += h * -0.2677;
#ifdef DEBUG_ADD_HYDROGEN
      cerr << "IsHydrogenAlcohol\n";
#endif
    } else if (IsHydroCarbon(pmd, i)) {  // H1
      result += h * 0.1230;
#ifdef DEBUG_ADD_HYDROGEN
      cerr << "IsHydroCarbon\n";
#endif
    } else {  // HS
      result += h * 0.1125;
#ifdef DEBUG_ADD_HYDROGEN
      cerr << "Default\n";
#endif
    }
#ifdef DEBUG_ADD_HYDROGEN
    cerr << " atom " << i << " h " << h << " result " << result << '\n';
#endif
  }

  return 1;
}

int
IsPrimaryAmine(const PerMoleculeData& pmd,  atom_number_t zatom) {
  const atom_number_t carbon = pmd.mol.other(zatom, 0);
  return pmd.unsaturation[carbon] == 0;
}

int
IsAcid(const PerMoleculeData& pmd,  atom_number_t zatom) {
  const atom_number_t carbon = pmd.mol.other(zatom, 0);
  if (pmd.z[carbon] == 6) {
  } else if (pmd.z[carbon] == 16) {
  } else {
    return 0;
  }

  if (pmd.unsaturation[carbon] == 0) {
    return 0;
  }

  for (const Bond* b : pmd.mol[carbon]) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t o = b->other(carbon);
    if (pmd.z[o] == 8) {
      return 1;
    }
    if (pmd.z[o] == 16) {
      return 1;
    }
    return 0;
  }

  return 0;
}

int
ALogP::AddZwitterionCorrection(PerMoleculeData& pmd, float& result) {
  Molecule& m = pmd.mol;
  const int matoms = m.natoms();

  int got_acid = 0;
  int got_nh2 = 0;
  for (int i = 0; i < matoms; ++i) {
    if (pmd.z[i] == 6) {
      continue;
    }

    if (pmd.z[i] == 7 && pmd.mol.ncon(i) == 1 && pmd.unsaturation[i] == 0 &&
        pmd.attached_heteroatom_count[i] == 0 && pmd.aryl_count[i] == 0 &&
        IsPrimaryAmine(pmd, i)) {
      got_nh2 = 1;
      if (got_acid) {
        break;
      }
    } else if (pmd.z[i] == 8 && pmd.mol.ncon(i) == 1 && pmd.aryl_count[i] == 0 &&
               IsAcid(pmd, i)) {
      got_acid = 1;
      if (got_nh2) {
        break;
      }
    }
  }

  if (! got_acid || ! got_nh2) {
    return 0;
  }

  result += -0.5;

  return 1;
}

std::optional<float>
ALogP::LogP(Molecule& m) {

  // Silently remove any explicit Hydrogen atoms.
  m.remove_all(1);

  const int matoms = m.natoms();
  if (matoms == 0) {
    return 0.0f;
  }

  // If we can exclude the case of single atom molecules, we avoid
  // a bunch of logic in the other functions, so in the interests of speed,
  // handle these separately.
  if (matoms == 1) {
    return SingleAtomSpecialCase(m);
  }

  m.compute_aromaticity_if_needed();

  PerMoleculeData pmd(m);

  float result = 0.0;
  for (int i = 0; i < matoms; ++i) {
    const atomic_number_t z = m.atomic_number(i);
    if (z == 6) {
      if (! Carbon(pmd, i, result)) {
        return std::nullopt;
      }
    } else if (z == 7) {
      if (! Nitrogen(pmd, i, result)) {
        return std::nullopt;
      }
    } else if (z == 8) {
      if (! Oxygen(pmd, i, result)) {
        return std::nullopt;
      }
    } else if (z == 9) {
      if (! Fluorine(pmd, i, result)) {
        return std::nullopt;
      }
    } else if (z == 15) {
      if (! Phosphorus(pmd, i, result)) {
        return std::nullopt;
      }
    } else if (z == 16) {
      if (! Sulphur(pmd, i, result)) {
        return std::nullopt;
      }
    } else if (z == 17) {
      if (! Chlorine(pmd, i, result)) {
        return std::nullopt;
      }
    } else if (z == 35) {
      if (! Bromine(pmd, i, result)) {
        return std::nullopt;
      }
    } else if (z == 53) {
      if (! Iodine(pmd, i, result)) {
        return std::nullopt;
      }
    } else {
      if (_display_error_messages) {
        cerr << "ALogP:unclassified atom " << Diagnostic(pmd, i) << '\n';
        return std::nullopt;
      }
    }
  }

#ifdef DEBUG_ADD_HYDROGEN
  std::cerr << "Sum before adding hudrogens " << result << '\n';
#endif

// #define ECHO_ATOMIC_CONTRIBUTIONS
#ifdef ECHO_ATOMIC_CONTRIBUTIONS
  float sum = 0.0f;
  for (int i = 0; i < matoms; ++i) {
    sum += pmd.atom_value[i];
    cerr << i << ' ' << pmd.atom_value[i] << ' ' << m.smarts_equivalent_for_atom(i) << '\t' << sum << '\n';
  }
  cerr << "Sum " << sum << '\n';
#endif

  if (! AddHydrogenContributions(pmd, result)) {
    return std::nullopt;
  }

  if (_label_with_atom_type) {
    m.set_isotopes(pmd.assigned);
  }

  if (_add_zwitterion_correction) {
    AddZwitterionCorrection(pmd, result);
  }

  return result;
}

ALogP::ALogP () {
}

std::optional<double>
ALogP::SingleAtomSpecialCase(Molecule& m) {
  assert(m.natoms() == 1);

  const atomic_number_t z = m.atomic_number(0);
  switch (z) {
    case 6:
      return 0.1441 + 4 * 0.1230;
    default:
      if (_display_error_messages) {
        PerMoleculeData pmd(m);
        cerr << "ALogP::SingleAtomSpecialCase:element not recognised " <<
                Diagnostic(pmd, 0) << '\n';
        return std::nullopt;
      }
      return std::nullopt;
  }
}

}  // namespace alogp


