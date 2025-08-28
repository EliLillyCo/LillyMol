// At some stage this should be updated to be able to read the proto based data files.

#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

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
  m.compute_aromaticity_if_needed();

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

ForFastScoring::ForFastScoring() {
  atype = nullptr;
  hcount = nullptr;
  htype = nullptr;
  zwit = 0;
}

ForFastScoring::~ForFastScoring() {
  if (atype != nullptr) {
    delete [] atype;
  }
  if (hcount != nullptr) {
    delete [] hcount;
  }
  if (htype != nullptr) {
    delete [] htype;
  }
}

int
ForFastScoring::resize(int s) {
  atype = new int[s];
  hcount = new int[s];
  htype = new int [s];

  return 1;
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

void
ALogP::DefaultParameters() {
  _params.value[kC1] = 0.1441;
  _params.value[kC2] = 0.0000;
  _params.value[kC3] = -0.2035;
  _params.value[kC4] = -0.2051;
  _params.value[kC5] = -0.2783;
  _params.value[kC6] = 0.1551;
  _params.value[kC7] = 0.00170;
  _params.value[kC8] = 0.08452;
  _params.value[kC9] = -0.1444;
  _params.value[kC10] = -0.0516;
  _params.value[kC11] = 0.1193;
  _params.value[kC12] = -0.0967;
  _params.value[kC13] = -0.5443;
  _params.value[kC14] = 0.0000;
  _params.value[kC15] = 0.2450;
  _params.value[kC16] = 0.1980;
  _params.value[kC17] = 0.000;
  _params.value[kC18] = 0.1581;
  _params.value[kC19] = 0.2955;
  _params.value[kC20] = 0.2713;
  _params.value[kC21] = 0.1360;
  _params.value[kC22] = 0.4619;
  _params.value[kC23] = 0.5437;
  _params.value[kC24] = 0.1893;
  _params.value[kC25] = -0.8186;
  _params.value[kC26] = 0.2640;
  _params.value[kC27] = 0.2148;   //  ‘[CX4][!(C,N,O,P,S,F,Cl,Br,I)]'. We only process organic elements.
  _params.value[kCS] = 0.08129;
  _params.value[kH1] = 0.1230;
  _params.value[kH2] = -0.2677;
  _params.value[kH3] = 0.2142;
  _params.value[kH4] = 0.2980;
  _params.value[kHS] = 0.1125;
  _params.value[kN1] = -1.0190;
  _params.value[kN2] = -0.7096;
  _params.value[kN3] = -1.0270;
  _params.value[kN4] = -0.5188;
  _params.value[kN5] = 0.08387;   // [NH+0]=A, [NH+0]=a imine. TODO:ianwatson fix
  _params.value[kN6] = 0.1836;
  _params.value[kN7] = -0.3187;
  _params.value[kN8] = -0.4458;
  _params.value[kN9] = 0.01508;
  _params.value[kN10] = -1.950;   // ‘[NH3+*]', ‘[NH2+*]', ‘[NH+*] protonated amine
  _params.value[kN11] = -0.3239;
  _params.value[kN12] = -1.119;
  _params.value[kN13] = -0.3396;
  _params.value[kN14] = 0.2887;
  _params.value[kNS] = -0.4806;
  _params.value[kO1] = 0.1552;
  _params.value[kO2] = -0.2893;
  _params.value[kO3] = -0.0684;
  _params.value[kO4] = -0.4195;
  _params.value[kO5] = 0.0335;
  _params.value[kO6] = -0.3339;
  _params.value[kO7] = -1.189;   // ‘[OX1−*][!(N,S)]' oxide
  _params.value[kO8] = 0.1788;
  _params.value[kO9] = -0.1526;
  _params.value[kO10] = 0.1129;
  _params.value[kO11] = 0.4833;
  _params.value[kO12] = -1.326;
  _params.value[kOS] = -0.1188;
  _params.value[kF] = 0.4202;
  _params.value[kCl] = 0.6895;
  _params.value[kBr] = 0.8456;
  _params.value[kI] = 0.8857;
  _params.value[kHal] = -2.996;   // [#9−*]', ‘[#17−*]', ‘[#35−*]', [#53−*]' ‘[#53+*]' ionic halogens
  _params.value[kP] = 0.8612;
  _params.value[kS1] = 0.6482;
  _params.value[kS2] = -0.0024;
  _params.value[kS3] = 0.6237;
  _params.value[kMe1] = -0.3808;   // All remaining p-block elements
  _params.value[kMe2] = -0.0025;   // All remaining d-block elements
  _params.value[kZwit] = -0.50;
  _params.value[kBias] = 0.0;
}

int
ALogP::AromaticCarbon(PerMoleculeData& pmd, atom_number_t zatom) {
  const Atom& a = pmd.mol[zatom];
  const int acon = a.ncon();

  // C18 [cH]
  if (acon == 2) {
    pmd.assigned[zatom] = kC18;
    return 1;
  }

  // C19 aromatic bridghead [c](:a)(:a):a
  if (acon == 3 && pmd.aromatic[zatom] == 3) {
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
    pmd.assigned[zatom] = kC25;
    return 1;
  }

  // C14
  if (pmd.z[o] == 9) {
    pmd.assigned[zatom] = kC14;
    return 1;
  }
  // C15
  if (pmd.z[o] == 17) {
    pmd.assigned[zatom] = kC15;
    return 1;
  }
  // C16
  if (pmd.z[o] == 35) {
    pmd.assigned[zatom] = kC16;
    return 1;
  }
  // C17
  if (pmd.z[o] == 53) {
    pmd.assigned[zatom] = kC17;
    return 1;
  }

  // C20 4-aromatic [c](:a)(:a)-a
  if (pmd.aromatic[o]) {
    pmd.assigned[zatom] = kC20;
    return 1;
  }
  // C21 [c](:a)(:a)-C
  if (pmd.z[o] == 6 && b->is_single_bond()) {
    pmd.assigned[zatom] = kC21;
    return 1;
  }
  // C22 [c](:a)(:a)-N
  if (pmd.z[o] == 7 && b->is_single_bond()) {
    pmd.assigned[zatom] = kC22;
    return 1;
  }
  // C23 [c](:a)(:a)-O
  if (pmd.z[o] == 8 && b->is_single_bond()) {
    pmd.assigned[zatom] = kC23;
    return 1;
  }
  // C24 [c](:a)(:a)-S
  if (pmd.z[o] == 16 && b->is_single_bond()) {
    pmd.assigned[zatom] = kC24;
    return 1;
  }
  // C13 aromatic heteratoms [cH0]-[!(C,N,O,S,F,Cl,Br,I)]’
  if (pmd.z[o] != 6 && b->is_single_bond()) {
    pmd.assigned[zatom] = kC13;
    return 1;
  }

  pmd.assigned[zatom] = kCS;

  return 1;
}

int
ALogP::SaturatedPrimaryCarbon(PerMoleculeData& pmd, atom_number_t zatom) {
  assert(pmd.z[zatom] == 6);
  assert(pmd.mol.ncon(zatom) == 1);

  const atom_number_t o = pmd.mol.other(zatom, 0);

  if (pmd.aromatic[o]) {
    if (pmd.aromatic[o] && pmd.z[o] == 6) {
      pmd.assigned[zatom] = kC8;
      return kC8;
    } else {
      pmd.assigned[zatom] = kC9;
      return kC9;
    }
  } else {  // aliphatic
    // C1
    if (pmd.z[o] == 6) {
      pmd.assigned[zatom] = kC1;
      return 1;   // kC1 is zero, so cannot return it. All these functions shuold be boolean.
    } else {  // heteroatom
      // C3 primary heteroatom
      pmd.assigned[zatom] = kC3;
      return kC3;
    }
  }
}

int
ALogP::SaturatedSecondaryCarbom(PerMoleculeData& pmd, atom_number_t zatom) {
  assert(pmd.z[zatom] == 6);
  assert(pmd.mol.ncon(zatom) == 2);

  // C10 [CH2X4]a  secondary aromatic
  if (pmd.aryl_count[zatom]) {
    pmd.assigned[zatom] = kC10;
    return kC10;
  }

  // C1 secondary aliphatic
  if (pmd.attached_heteroatom_count[zatom] == 0) {
    pmd.assigned[zatom] = kC1;
    return 1;
  }


  // C3
  pmd.assigned[zatom] = kC3;
  return 1;
}

int
ALogP::SaturatedCarbon(PerMoleculeData& pmd, atom_number_t zatom) {
  const Atom& a = pmd.mol[zatom];
  const int acon = a.ncon();

#ifdef NOW_BEING_HANDLED_EXERNALLY
  // C1
  if (a.ncon() == 0) {
    result += _params.value[q];  // 0.1441;
    pmd.assigned[zatom] = kC1;
    return 1;
  }
#endif

  const int ahc = pmd.attached_heteroatom_count[zatom];

  if (acon == 1) {
    return SaturatedPrimaryCarbon(pmd, zatom);
  }

  if (acon == 2) {
    return SaturatedSecondaryCarbom(pmd, zatom);
  }

  // C11
  if (acon == 3 && pmd.aryl_count[zatom]) {
    pmd.assigned[zatom] = kC11;
    return kC11;
  }

  // C12
  if (acon == 4 && pmd.aryl_count[zatom]) {
    pmd.assigned[zatom] = kC12;
    return 1;
  }

  // C2
  if (ahc == 0) {
    pmd.assigned[zatom] = kC2;
    return 1;
  }
  // C4
  pmd.assigned[zatom] = kC4;
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
ALogP::UnSaturatedCarbon(PerMoleculeData& pmd, atom_number_t zatom) {
  // C7 [CX2]#A acetylene, nitrile
  if (pmd.triple_bond_count[zatom] && pmd.formal_charge[zatom] == 0) {
    pmd.assigned[zatom] = kC7;
    return 1;
  }

  // CS [C-]#N
  if (pmd.formal_charge[zatom] == -1 && pmd.triple_bond_count[zatom]) {
    pmd.assigned[zatom] = kCS;
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
    pmd.assigned[zatom] = kC5;
    return 1;
  }

  // At this stage the doubly bonded atom is definitely a carbon.

  // C6 C=C aliphatic
  if (pmd.aryl_count[zatom] == 0) {
    pmd.assigned[zatom] = kC6;
    return 1;
  }

  // C26 C=C aromatic [C]()C)(a)A’, ‘[C]()C)(c)a’, ‘[CH]() C)a’, ‘[C] ) c
  pmd.assigned[zatom] = kC26;
  return 1;
}


int
ALogP::Carbon(PerMoleculeData& pmd, atom_number_t zatom) {
  if (pmd.aromatic[zatom]) {
    return AromaticCarbon(pmd, zatom);
  }

  if (pmd.unsaturation[zatom] == 0) {
    return SaturatedCarbon(pmd, zatom);
  }

  return UnSaturatedCarbon(pmd, zatom);
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
ALogP::AromaticNitrogen(PerMoleculeData& pmd, atom_number_t zatom) {
  // N11 Unprotonated aromatic
  const Atom& a = pmd.mol[zatom];

  if (a.ncon() == 2 && pmd.formal_charge[zatom] == 0) {
    pmd.assigned[zatom] = kN11;
    return 1;
  }

  // N1(=C2C(=N(=O)C=C1)C=CC=C2)=O CHEMBL2104626
  // which can be written in charge separate form. RDKit
  // assigns this as if it is the charge separated form.
  // Probably that is what is intended in the paper.
  // N12
  if (a.ncon() == 3 && pmd.formal_charge[zatom] == 0 &&
      pmd.double_bond_count[zatom] == 2) {
    pmd.assigned[zatom] = kN12;
    return 1;
  }

  if (pmd.formal_charge[zatom] == 0) {
    pmd.assigned[zatom] = kN11;
    return 1;
  }

  // N12 protonated aromatic
  pmd.assigned[zatom] = kN12;
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
ALogP::SinglyConnectedSaturatedNitrogen(PerMoleculeData& pmd, atom_number_t zatom) {
  assert(pmd.mol.ncon(zatom) == 1);

  if (pmd.attached_heteroatom_count[zatom] == 0 &&
      SaturatedNitrogenIsGuanidine(pmd, zatom)) {
    pmd.assigned[zatom] = kN1;
    return 1;
  }

  if (pmd.attached_heteroatom_count[zatom] == 0 &&
      SaturatedNitrogenIsAmidine(pmd, zatom)) {
    pmd.assigned[zatom] = kN1;
    return 1;
  }

  const int vinyl = Vinyl(pmd, zatom);

  // We are not using a charge assigner, so we approximate a primary amine this way.
  if (pmd.attached_heteroatom_count[zatom] == 0 &&
      vinyl == 0) {
    atom_number_t o = pmd.mol.other(zatom, 0);
    if (pmd.aromatic[o]) {  // N3
      pmd.assigned[zatom] = kN3;
      return 1;
    } else {  // N1
      pmd.assigned[zatom] = kN1;
      return 1;
    }
  }

  // but notice that RDKit assigns N1 to a terminal amide
  if (pmd.attached_heteroatom_count[zatom] == 0 &&
      vinyl && TerminalNitrogenIsAmide(pmd, zatom)) {
    pmd.assigned[zatom] = kN1;
    return 1;
  }

  if (pmd.attached_heteroatom_count[zatom] == 1 &&
      vinyl && IsSulfonamide(pmd, zatom)) {
    pmd.assigned[zatom] = kN1;
    return 1;
  }

  // N3 maybe not charged, but still adjacent to arom. Aniline
  if (pmd.aryl_count[zatom] == 1) {
    pmd.assigned[zatom] = kN3;
    return 1;
  }

  // NS nitrogen supplemental.
  pmd.assigned[zatom] = kNS;
  return 1;
}

int
ALogP::SaturatedNitrogen(PerMoleculeData& pmd, atom_number_t zatom) {
  const Atom& a = pmd.mol[zatom];
  // cerr << "SaturatedNitrogen " << pmd.mol.aromatic_smiles() << '\n';

  if (a.ncon() == 1) {
    return SinglyConnectedSaturatedNitrogen(pmd, zatom);
  }

  const int vinyl = Vinyl(pmd, zatom);
  const int aryl = ArylCount(pmd, zatom);

  if (a.ncon() == 2 && pmd.attached_heteroatom_count[zatom] == 0 &&
     Vinyl(pmd, zatom) == 0) {
    // N4
    if (aryl) {
      pmd.assigned[zatom] = kN4;
      return 1;
    } else {  // N2
      pmd.assigned[zatom] = kN2;
      return 1;
    }
  }

  // N4 - not charged
  if (a.ncon() == 2 && pmd.aryl_count[zatom]) {
    pmd.assigned[zatom] = kN4;
    return 1;
  }

  // N2 - not charged
  if (a.ncon() == 2 && aryl == 0) {
    pmd.assigned[zatom] = kN2;
    return 1;
  }

  if (a.ncon() == 3 && pmd.attached_heteroatom_count[zatom] == 0 &&
      Vinyl(pmd, zatom) == 0) {
    // N8
    if (aryl) {
      pmd.assigned[zatom] = kN8;
      return 1;
    } else {  // N7
      pmd.assigned[zatom] = kN7;
      return 1;
    }
  }

  // N8 not charged
  if (a.ncon() == 3 && aryl) {
    pmd.assigned[zatom] = kN8;
    return 1;
  }

  // N7 but with vinyl connections.
  if (a.ncon() == 3 && pmd.attached_heteroatom_count[zatom] == 0 &&
      pmd.aryl_count[zatom] == 0) {
    pmd.assigned[zatom] = kN7;
    return 1;
  }

  // N7 that might be part of a sulfonamide
  if (a.ncon() == 3 && pmd.attached_heteroatom_count[zatom] == 1 &&
      vinyl == 1 && IsSulfonamide(pmd, zatom)) {
    pmd.assigned[zatom] = kN7;
    return 1;
  }

  // N7 with very few restrictions
  if (a.ncon() == 3 && pmd.attached_heteroatom_count[zatom] == 1 &&
      vinyl == 0) {
    pmd.assigned[zatom] = kN7;
    return 1;
  }

  // N7 with hardly any restriction. Note that these would probably
  // not be charged.
  if (a.ncon() == 3 && pmd.attached_heteroatom_count[zatom] == 1) {
    pmd.assigned[zatom] = kN7;
    return 1;
  }

  // N7 with hardly any restriction. Note that these would probably
  // not be charged.
  if (a.ncon() == 3 && pmd.attached_heteroatom_count[zatom] == 2) {
    pmd.assigned[zatom] = kN7;
    return 1;
  }

  if (a.ncon() == 4 && pmd.formal_charge[zatom] == 1) {
    // N13
    pmd.assigned[zatom] = kN13;
    return 1;
  }

  // N14
  if (pmd.formal_charge[zatom]) {
    pmd.assigned[zatom] = kN14;
    return 1;
  }

  // N14
  if (a.ncon() == 2 && pmd.unsaturation[zatom] == 3) {
    pmd.assigned[zatom] = kN14;
    return 1;
  }

  // N14
  if (pmd.unsaturation[zatom] == 2 && a.ncon() == 1 &&
      pmd.attached_heteroatom_count[zatom] == 1) {
    pmd.assigned[zatom] = kN14;
    return 1;
  }

  // NS nitrogen supplemental.
  pmd.assigned[zatom] = kNS;
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
ALogP::UnSaturatedNitrogen(PerMoleculeData& pmd, atom_number_t zatom) {
  const Atom& a = pmd.mol[zatom];

  // N9 Nitrile
  if (a.ncon() == 1 && pmd.triple_bond_count[zatom] &&
      pmd.attached_heteroatom_count[zatom] == 0) {
    pmd.assigned[zatom] = kN9;
    return 1;
  }

  // N14 in [C-]#[N+]
  if (a.ncon() == 2 && pmd.formal_charge[zatom] == 1 &&
      pmd.triple_bond_count[zatom] == 1) {
    pmd.assigned[zatom] = kN14;
    return 1;
  }


  // Nitro is assigned N13
  if (a.ncon() == 3 && pmd.double_bond_count[zatom] == 2 &&
      pmd.attached_heteroatom_count[zatom] >= 2 && IsNitro(pmd, zatom)) {
    pmd.assigned[zatom] = kN13;
    return 1;
  }

  // N5 imine
  if (a.ncon() == 1 && pmd.double_bond_count[zatom] == 1) {
    pmd.assigned[zatom] = kN5;
    return 1;
  }

  // N6 substituted imine
  if (a.ncon() == 2 && pmd.double_bond_count[zatom] == 1 && pmd.triple_bond_count[zatom] == 0) {
    pmd.assigned[zatom] = kN6;
    return 1;
  }

  // N14 N#N=
  if (a.ncon() == 1 && pmd.triple_bond_count[zatom] == 1 &&
      pmd.attached_heteroatom_count[zatom]) {
    pmd.assigned[zatom] = kN14;
    return 1;
  }

  // N14 N#N=
  if (a.ncon() == 2 && pmd.triple_bond_count[zatom] == 1 &&
      pmd.double_bond_count[zatom] == 1 &&
      pmd.attached_heteroatom_count[zatom]) {
    pmd.assigned[zatom] = kN14;
    return 1;
  }

  // NS nitrogen supplemental.
  pmd.assigned[zatom] = kNS;
  return 1;
}

int
ALogP::Nitrogen(PerMoleculeData& pmd, atom_number_t zatom) {
  if (pmd.aromatic[zatom]) {
    return AromaticNitrogen(pmd, zatom);
  }

  if (pmd.unsaturation[zatom] == 0) {
    return SaturatedNitrogen(pmd, zatom);
  }

  return UnSaturatedNitrogen(pmd, zatom);
}

int
ALogP::UnSaturatedOxygen(PerMoleculeData& pmd, atom_number_t zatom) {
  assert(pmd.double_bond_count[zatom] == 1);

  const atom_number_t o = pmd.mol.other(zatom, 0);
  // O5
  if (pmd.z[o] == 7) {
    pmd.assigned[zatom] = kO5;
    return 1;
  }

  // O8
  if (pmd.z[o] == 6 && pmd.aromatic[o]) {
    pmd.assigned[zatom] = kO8;
    return 1;
  }

  // IAW. Copy behaviour from RDKit which assigns O6 to the oxygens in O=S=O
  // which seems quite wrong.
  if (pmd.z[o] == 16) {
    pmd.assigned[zatom] = kO6;
    return 1;
  }

  // O11 carbonyl heteroatom
  if (pmd.z[o] == 6 && pmd.attached_heteroatom_count[o] == 3) {
    pmd.assigned[zatom] = kO11;
    return 1;
  }

  // O10 carbonyl aromatic
  if (pmd.attached_heteroatom_count[zatom] == 0 && ArylCount(pmd, o)) {
    pmd.assigned[zatom] = kO10;
    return 1;
  }

  // O9 carbonyl aliphatic
  if (pmd.z[o] == 6) {
    pmd.assigned[zatom] = kO9;
    return 1;
  }

  // OS
  pmd.assigned[zatom] = kOS;

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
ALogP::SaturatedOxygen(PerMoleculeData& pmd, atom_number_t zatom) {
  const Atom& a = pmd.mol[zatom];
  const int acon = a.ncon();

  // cerr << "ALogP::SaturatedOxygen acid " << IsAcid(pmd, zatom) << '\n';
  // O12
  if (acon == 1 && IsAcid(pmd, zatom)) {
    if (_use_alcohol_for_acid) {
      pmd.assigned[zatom] = kO2;
    } else {
      pmd.assigned[zatom] = kO12;
    }
    return 1;
  }

  // O2 alcohol
  if (acon == 1 && pmd.attached_heteroatom_count[zatom] == 0) {
    pmd.assigned[zatom] = kO2;
    return 1;
  }

  // O2 N-Hydroxy
  if (acon == 1 && pmd.attached_heteroatom_count[zatom] == 1 &&
     pmd.aryl_count[zatom] == 0 && IsNHydroxy(pmd, zatom)) {
    pmd.assigned[zatom] = kO2;
    return 1;
  }

  int aac = AttachedAromaticCount(pmd, zatom);
  // O3 Aliphatic ether.
  if (acon == 2 && aac == 0) {
    pmd.assigned[zatom] = kO3;
    return 1;
  }

  // O4 Aromatic ether
  if (acon == 2 && aac) {
    pmd.assigned[zatom] = kO4;
    return 1;
  }

  if (acon == 1) {
    atom_number_t o = pmd.mol.other(zatom, 0);
    // O6 oxide
    if (pmd.z[o] == 16) {
      pmd.assigned[zatom] = kO6;
      return 1;
    }

    // N-oxide Note that this is grouped with other queries in the paper.
    // O5
    if (pmd.z[o] == 7) {
      pmd.assigned[zatom] = kO5;
      return 1;
    }

    // O6
    if (pmd.z[o] == 16) {
      pmd.assigned[zatom] = kO6;
      return 1;
    }
  }

  // OS oxygen supplemental.
  pmd.assigned[zatom] = kOS;
  return 1;
}


int
ALogP::Oxygen(PerMoleculeData& pmd, atom_number_t zatom) {
  // O1 aromatic
  if (pmd.aromatic[zatom]) {
    pmd.assigned[zatom] = kO1;
    return 1;
  }

  if (pmd.double_bond_count[zatom]) {
    return UnSaturatedOxygen(pmd, zatom);
  }

  return SaturatedOxygen(pmd, zatom);
}

int
ALogP::Fluorine(PerMoleculeData& pmd, atom_number_t zatom) {
  pmd.assigned[zatom] = kF;
  return 1;
}

int
ALogP::Chlorine(PerMoleculeData& pmd, atom_number_t zatom) {
  pmd.assigned[zatom] = kCl;
  return 1;
}

int
ALogP::Bromine(PerMoleculeData& pmd, atom_number_t zatom) {
  pmd.assigned[zatom] = kBr;
  return 1;
}

int
ALogP::Iodine(PerMoleculeData& pmd, atom_number_t zatom) {
  pmd.assigned[zatom] = kI;
  return 1;
}

int
ALogP::Phosphorus(PerMoleculeData& pmd, atom_number_t zatom) {
  pmd.assigned[zatom] = kP;
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
ALogP::Sulphur(PerMoleculeData& pmd, atom_number_t zatom) {
  // S3
  if (pmd.aromatic[zatom]) {
    pmd.assigned[zatom] = kS3;
    return 1;
  }

  // RDKit seems to apply S2 to S=O groups, which seems wrong
  // Maybe this should be made a settable behaviour.
  if (pmd.double_bond_count[zatom] && DoublyBondedToOxygen(pmd.mol, zatom)) {
    pmd.assigned[zatom] = kS2;
    return 1;
  }

  // Not sure this can happen...
  // S2
  if (pmd.formal_charge[zatom]) {
    pmd.assigned[zatom] = kS2;
    return 1;
  }

  // S1
  pmd.assigned[zatom] = kS1;
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
      result += h* _params.value[kH4];  // h * 0.2980;
#ifdef DEBUG_ADD_HYDROGEN
      cerr << "IsHydrogenAcid\n";
#endif
    } else if (IsPeroxide(pmd, i)) {   // H4
      result += h * _params.value[kH4];  // h * 0.2980;
#ifdef DEBUG_ADD_HYDROGEN
      cerr << "IsPeroxide\n";
#endif
    } else if (IsHydrogenAmine(pmd, i)) {  // H3
      result += h * _params.value[kH3];  // h * 0.2142;
#ifdef DEBUG_ADD_HYDROGEN
      cerr << "IsHydrogenAmine\n";
#endif
    } else if (IsHydrogenAlcohol(pmd, i)) {  // H2
      result += h * _params.value[kH2];  // h * -0.2677;
#ifdef DEBUG_ADD_HYDROGEN
      cerr << "IsHydrogenAlcohol\n";
#endif
    } else if (IsHydroCarbon(pmd, i)) {  // H1
      result += (h * _params.value[kH1]);  // h * 0.1230;
#ifdef DEBUG_ADD_HYDROGEN
      cerr << "IsHydroCarbon " << kH1 << " hcount " << h << " value " << _params.value[kH1] << '\n';
#endif
    } else {  // HS
      result += h * _params.value[kHS];  // h * 0.1125;
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
ALogP::GetHydrogenContributions(Molecule& m, PerMoleculeData& pmd, ForFastScoring& for_fast_scoring) {

  const int matoms = m.natoms();

  int* hcount = for_fast_scoring.hcount;
  int* htype = for_fast_scoring.htype;

  for (int i = 0; i < matoms; ++i) {
    hcount[i] = m.hcount(i);
    htype[i] = -1;

    if (_use_alcohol_for_acid && pmd.formal_charge[i] == -1 &&
        pmd.assigned[i] == kO2) {
      hcount[i] = 1;
    }
    if (hcount[i] > 0 && _rdkit_charged_nitrogen && pmd.formal_charge[i] == 1 &&
        (pmd.assigned[i] == kN1 || pmd.assigned[i] == kN2 ||
         pmd.assigned[i] == kN7)) {
      --hcount[i];
    }

    if (hcount[i] == 0) {
      continue;
    }

    if (IsHydrogenAcid(pmd, i)) {  // H4
      htype[i] = kH4;
    } else if (IsPeroxide(pmd, i)) { [[unlikely]]  // H4
      htype[i] = kH4;
    } else if (IsHydrogenAmine(pmd, i)) {  // H3
      htype[i] = kH3;
    } else if (IsHydrogenAlcohol(pmd, i)) {  // H2
      htype[i] = kH2;
    } else if (IsHydroCarbon(pmd, i)) {  // H1
      htype[i] = kH1;
    } else {  // HS
      htype[i] = kHS;
    }
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

  result += _params.value[kZwit];  // -0.5;   not sure what this is

  return 1;
}

std::optional<float>
ALogP::LogP(Molecule& m) {

  // Silently remove any explicit Hydrogen atoms.
  m.remove_all(1);

  const int matoms = m.natoms();

  if (matoms == 0) [[unlikely]] {
    return 0.0f;
  }

  if (matoms == 1) [[unlikely]] {
    return SingleAtomSpecialCase(m);
  }

  // If the current molecule contains isotopes, we save and restore them.
  std::unique_ptr<isotope_t[]> isosave;

  if (m.ContainsIsotopicAtoms())  [[unlikely]] {
    isosave = m.GetIsotopes();
    m.transform_to_non_isotopic_form();
  }

  PerMoleculeData pmd(m);

  std::optional<float> result = LogPInner(m, pmd);

  if (isosave) [[unlikely]] {
    m.set_isotopes(isosave.get());
  }

  return result;
}

std::optional<float>
ALogP::LogPInner(Molecule& m, PerMoleculeData& pmd) {
  // If we can exclude the case of single atom molecules, we avoid
  // a bunch of logic in the other functions, so in the interests of speed,
  // handle these separately.

  const int matoms = m.natoms();

  if (matoms == 1) [[unlikely]] {
    return SingleAtomSpecialCase(m);
  }

  for (int i = 0; i < matoms; ++i) {
    const atomic_number_t z = m.atomic_number(i);
    if (z == 6) {
      if (! Carbon(pmd, i)) {
        return std::nullopt;
      }
    } else if (z == 7) {
      if (! Nitrogen(pmd, i)) {
        return std::nullopt;
      }
    } else if (z == 8) {
      if (! Oxygen(pmd, i)) {
        return std::nullopt;
      }
    } else if (z == 9) {
      if (! Fluorine(pmd, i)) {
        return std::nullopt;
      }
    } else if (z == 15) {
      if (! Phosphorus(pmd, i)) {
        return std::nullopt;
      }
    } else if (z == 16) {
      if (! Sulphur(pmd, i)) {
        return std::nullopt;
      }
    } else if (z == 17) {
      if (! Chlorine(pmd, i)) {
        return std::nullopt;
      }
    } else if (z == 35) {
      if (! Bromine(pmd, i)) {
        return std::nullopt;
      }
    } else if (z == 53) {
      if (! Iodine(pmd, i)) {
        return std::nullopt;
      }
    } else if (z < 17) {
      pmd.assigned[i] = kMe1;
    } else if (z < 35) {
      pmd.assigned[i] = kMe2;
    } else if (_fail_if_unclassified_atom) {
      if (_display_error_messages) {
        cerr << "ALogP:unclassified atom " << Diagnostic(pmd, i) << '\n';
      }
      return std::nullopt;
    } else {
      // This is wrong, but failures are just not going to happen.
      pmd.assigned[i] = kC1;
    }
  }

  float result = _params.value[kBias];

  for (int i = 0; i < matoms; ++i) {
    int atype = pmd.assigned[i];
//  cerr << " atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " assuigned " << atype << " value " << _params.value[atype] << '\n';
    result += _params.value[atype];
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

  if (_add_zwitterion_correction) {
    AddZwitterionCorrection(pmd, result);
  }

  if (_label_with_atom_type) {
    m.set_isotopes(pmd.assigned);
  }

  return result;
}

ALogP::ALogP () {
  DefaultParameters();
}

std::optional<double>
ALogP::SingleAtomSpecialCase(Molecule& m) {
  assert(m.natoms() == 1);

  const atomic_number_t z = m.atomic_number(0);
  switch (z) {
    case 6:
      return _params.value[kC1] + 4 * _params.value[kH1];
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

int
ALogP::ReadConfiguration(IWString& fname) {
  std::optional<alogp::AlogpConfiguration> maybe_proto =
        iwmisc::ReadTextProto<alogp::AlogpConfiguration>(fname);

  if (! maybe_proto) {
    cerr << "ALogP::ReadTextProto:cannot read textproto '" << fname << "'\n";
    return 0;
  }

  return ConfigFromProto(*maybe_proto);
}

int
ALogP::ConfigFromProto(const alogp::AlogpConfiguration& proto) {
  // cerr << "Alogp::ConfigFromProto: with " << proto.ShortDebugString() << '\n';
  if (! proto.has_parameters()) {
    cerr << "ALogP::ConfigFromProto:no parameters\n";
    return 0;
  }

  if (! ParametersFromProto(proto.parameters())) {
    cerr << "ALogP::ConfigFromProto:cannot transfer parameters from proto\n";
    return 0;
  }

  if (proto.has_add_zwitterion_correction()) {
    _add_zwitterion_correction = proto.add_zwitterion_correction();
  }

  if (proto.has_rdkit_phoshoric_acid_hydrogen()) {
    _rdkit_phoshoric_acid_hydrogen = proto.rdkit_phoshoric_acid_hydrogen();
  }

  if (proto.has_rdkit_charged_nitrogen()) {
    _rdkit_charged_nitrogen = proto.rdkit_charged_nitrogen();
  }

  if (proto.has_use_alcohol_for_acid()) {
    _use_alcohol_for_acid = proto.use_alcohol_for_acid();
  }

  return 1;
}

int
ALogP::ParametersFromProto(const alogp::AlogpParameters& proto) {
  if (! proto.has_c1()) {
    cerr << "No value for c1\n";
    return 0;
  }
  _params.value[kC1] = proto.c1();

  if (! proto.has_c2()) {
    cerr << "No value for c2\n";
    return 0;
  }
  _params.value[kC2] = proto.c2();

  if (! proto.has_c3()) {
    cerr << "No value for c3\n";
    return 0;
  }
  _params.value[kC3] = proto.c3();

  if (! proto.has_c4()) {
    cerr << "No value for c4\n";
    return 0;
  }
  _params.value[kC4] = proto.c4();

  if (! proto.has_c5()) {
    cerr << "No value for c5\n";
    return 0;
  }
  _params.value[kC5] = proto.c5();

  if (! proto.has_c6()) {
    cerr << "No value for c6\n";
    return 0;
  }
  _params.value[kC6] = proto.c6();

  if (! proto.has_c7()) {
    cerr << "No value for c7\n";
    return 0;
  }
  _params.value[kC7] = proto.c7();

  if (! proto.has_c8()) {
    cerr << "No value for c8\n";
    return 0;
  }
  _params.value[kC8] = proto.c8();

  if (! proto.has_c9()) {
    cerr << "No value for c9\n";
    return 0;
  }
  _params.value[kC9] = proto.c9();

  if (! proto.has_c10()) {
    cerr << "No value for c10\n";
    return 0;
  }
  _params.value[kC10] = proto.c10();

  if (! proto.has_c11()) {
    cerr << "No value for c11\n";
    return 0;
  }
  _params.value[kC11] = proto.c11();

  if (! proto.has_c12()) {
    cerr << "No value for c12\n";
    return 0;
  }
  _params.value[kC12] = proto.c12();

  if (! proto.has_c13()) {
    cerr << "No value for c13\n";
    return 0;
  }
  _params.value[kC13] = proto.c13();

  if (! proto.has_c14()) {
    cerr << "No value for c14\n";
    return 0;
  }
  _params.value[kC14] = proto.c14();

  if (! proto.has_c15()) {
    cerr << "No value for c15\n";
    return 0;
  }
  _params.value[kC15] = proto.c15();

  if (! proto.has_c16()) {
    cerr << "No value for c16\n";
    return 0;
  }
  _params.value[kC16] = proto.c16();

  if (! proto.has_c17()) {
    cerr << "No value for c17\n";
    return 0;
  }
  _params.value[kC17] = proto.c17();

  if (! proto.has_c18()) {
    cerr << "No value for c18\n";
    return 0;
  }
  _params.value[kC18] = proto.c18();

  if (! proto.has_c19()) {
    cerr << "No value for c19\n";
    return 0;
  }
  _params.value[kC19] = proto.c19();

  if (! proto.has_c20()) {
    cerr << "No value for c20\n";
    return 0;
  }
  _params.value[kC20] = proto.c20();

  if (! proto.has_c21()) {
    cerr << "No value for c21\n";
    return 0;
  }
  _params.value[kC21] = proto.c21();

  if (! proto.has_c22()) {
    cerr << "No value for c22\n";
    return 0;
  }
  _params.value[kC22] = proto.c22();

  if (! proto.has_c23()) {
    cerr << "No value for c23\n";
    return 0;
  }
  _params.value[kC23] = proto.c23();

  if (! proto.has_c24()) {
    cerr << "No value for c24\n";
    return 0;
  }
  _params.value[kC24] = proto.c24();

  if (! proto.has_c25()) {
    cerr << "No value for c25\n";
    return 0;
  }
  _params.value[kC25] = proto.c25();

  if (! proto.has_c26()) {
    cerr << "No value for c26\n";
    return 0;
  }
  _params.value[kC26] = proto.c26();

  if (! proto.has_c27()) {
    cerr << "No value for c27\n";
    return 0;
  }
  _params.value[kC27] = proto.c27();

  if (! proto.has_cs()) {
    cerr << "No value for cs\n";
    return 0;
  }
  _params.value[kCS] = proto.cs();

  if (! proto.has_h1()) {
    cerr << "No value for h1\n";
    return 0;
  }
  _params.value[kH1] = proto.h1();

  if (! proto.has_h2()) {
    cerr << "No value for h2\n";
    return 0;
  }
  _params.value[kH2] = proto.h2();

  if (! proto.has_h3()) {
    cerr << "No value for h3\n";
    return 0;
  }
  _params.value[kH3] = proto.h3();

  if (! proto.has_h4()) {
    cerr << "No value for h4\n";
    return 0;
  }
  _params.value[kH4] = proto.h4();

  if (! proto.has_hs()) {
    cerr << "No value for hs\n";
    return 0;
  }
  _params.value[kHS] = proto.hs();

  if (! proto.has_n1()) {
    cerr << "No value for n1\n";
    return 0;
  }
  _params.value[kN1] = proto.n1();

  if (! proto.has_n2()) {
    cerr << "No value for n2\n";
    return 0;
  }
  _params.value[kN2] = proto.n2();

  if (! proto.has_n3()) {
    cerr << "No value for n3\n";
    return 0;
  }
  _params.value[kN3] = proto.n3();

  if (! proto.has_n4()) {
    cerr << "No value for n4\n";
    return 0;
  }
  _params.value[kN4] = proto.n4();

  if (! proto.has_n5()) {
    cerr << "No value for n5\n";
    return 0;
  }
  _params.value[kN5] = proto.n5();

  if (! proto.has_n6()) {
    cerr << "No value for n6\n";
    return 0;
  }
  _params.value[kN6] = proto.n6();

  if (! proto.has_n7()) {
    cerr << "No value for n7\n";
    return 0;
  }
  _params.value[kN7] = proto.n7();

  if (! proto.has_n8()) {
    cerr << "No value for n8\n";
    return 0;
  }
  _params.value[kN8] = proto.n8();

  if (! proto.has_n9()) {
    cerr << "No value for n9\n";
    return 0;
  }
  _params.value[kN9] = proto.n9();

  if (! proto.has_n10()) {
    cerr << "No value for n10\n";
    return 0;
  }
  _params.value[kN10] = proto.n10();

  if (! proto.has_n11()) {
    cerr << "No value for n11\n";
    return 0;
  }
  _params.value[kN11] = proto.n11();

  if (! proto.has_n12()) {
    cerr << "No value for n12\n";
    return 0;
  }
  _params.value[kN12] = proto.n12();

  if (! proto.has_n13()) {
    cerr << "No value for n13\n";
    return 0;
  }
  _params.value[kN13] = proto.n13();

  if (! proto.has_n14()) {
    cerr << "No value for n14\n";
    return 0;
  }
  _params.value[kN14] = proto.n14();

  if (! proto.has_ns()) {
    cerr << "No value for ns\n";
    return 0;
  }
  _params.value[kNS] = proto.ns();

  if (! proto.has_o1()) {
    cerr << "No value for o1\n";
    return 0;
  }
  _params.value[kO1] = proto.o1();

  if (! proto.has_o2()) {
    cerr << "No value for o2\n";
    return 0;
  }
  _params.value[kO2] = proto.o2();

  if (! proto.has_o3()) {
    cerr << "No value for o3\n";
    return 0;
  }
  _params.value[kO3] = proto.o3();

  if (! proto.has_o4()) {
    cerr << "No value for o4\n";
    return 0;
  }
  _params.value[kO4] = proto.o4();

  if (! proto.has_o5()) {
    cerr << "No value for o5\n";
    return 0;
  }
  _params.value[kO5] = proto.o5();

  if (! proto.has_o6()) {
    cerr << "No value for o6\n";
    return 0;
  }
  _params.value[kO6] = proto.o6();

  if (! proto.has_o7()) {
    cerr << "No value for o7\n";
    return 0;
  }
  _params.value[kO7] = proto.o7();

  if (! proto.has_o8()) {
    cerr << "No value for o8\n";
    return 0;
  }
  _params.value[kO8] = proto.o8();

  if (! proto.has_o9()) {
    cerr << "No value for o9\n";
    return 0;
  }
  _params.value[kO9] = proto.o9();

  if (! proto.has_o10()) {
    cerr << "No value for o10\n";
    return 0;
  }
  _params.value[kO10] = proto.o10();

  if (! proto.has_o11()) {
    cerr << "No value for o11\n";
    return 0;
  }
  _params.value[kO11] = proto.o11();

  if (! proto.has_o12()) {
    cerr << "No value for o12\n";
    return 0;
  }
  _params.value[kO12] = proto.o12();

  if (! proto.has_os()) {
    cerr << "No value for os\n";
    return 0;
  }
  _params.value[kOS] = proto.os();

  if (! proto.has_f()) {
    cerr << "No value for f\n";
    return 0;
  }
  _params.value[kF] = proto.f();

  if (! proto.has_cl()) {
    cerr << "No value for cl\n";
    return 0;
  }
  _params.value[kCl] = proto.cl();

  if (! proto.has_br()) {
    cerr << "No value for br\n";
    return 0;
  }
  _params.value[kBr] = proto.br();

  if (! proto.has_i()) {
    cerr << "No value for i\n";
    return 0;
  }
  _params.value[kI] = proto.i();

  if (! proto.has_hal()) {
    cerr << "No value for hal\n";
    return 0;
  }
  _params.value[kHal] = proto.hal();

  if (! proto.has_p()) {
    cerr << "No value for p\n";
    return 0;
  }
  _params.value[kP] = proto.p();

  if (! proto.has_s1()) {
    cerr << "No value for s1\n";
    return 0;
  }
  _params.value[kS1] = proto.s1();

  if (! proto.has_s2()) {
    cerr << "No value for s2\n";
    return 0;
  }
  _params.value[kS2] = proto.s2();

  if (! proto.has_s3()) {
    cerr << "No value for s3\n";
    return 0;
  }
  _params.value[kS3] = proto.s3();

  if (! proto.has_me1()) {
    cerr << "No value for me1\n";
    return 0;
  }
  _params.value[kMe1] = proto.me1();

  if (! proto.has_me2()) {
    cerr << "No value for me2\n";
    return 0;
  }
  _params.value[kMe2] = proto.me2();

  if (! proto.has_zwit()) {
    cerr << "No value for zwit\n";
    return 0;
  }
  _params.value[kZwit] = proto.zwit();

  if (! proto.has_bias()) {
    cerr << "No value for bias\n";
    return 0;
  }
  _params.value[kBias] = proto.bias();

  return 1;
}

template <typename T>
int
ALogP::SetWeights(uint32_t n, const T* values) {
  if (n != kLast) {
    cerr << "AlogP::SetWeights:size mismatch, got " << n <<
            " should be " << kLast << '\n';
    return 0;   
  }

  for (uint32_t i = 0; i < n; ++i) {
    _params.value[i] = values[i];
  }

  return 1;   
}

template int ALogP::SetWeights<double>(uint32_t, const double*);

int
ALogP::FillForFastScoring(Molecule& m, ForFastScoring& for_fast_scoring) {
  const int matoms = m.natoms();

  if (m.empty()) [[unlikely]] {
    return 0;
  }

  PerMoleculeData pmd(m);

  std::optional<float> pred = LogPInner(m, pmd);
  if (! pred) {
    return 0;
  }

  for_fast_scoring.resize(matoms);

  for (int i = 0; i < matoms; ++i) {
    for_fast_scoring.atype[i] = pmd.assigned[i];
  }


  return GetHydrogenContributions(m, pmd, for_fast_scoring);
}

float
ALogP::LogP(Molecule& m, const ForFastScoring& for_fast_scoring) {
  float result = _params.value[kBias];

  const int matoms = m.natoms();
  if (matoms == 1) [[ unlikely ]] {
    return *SingleAtomSpecialCase(m);
  }

  const int* atype = for_fast_scoring.atype;
  const int* hcount = for_fast_scoring.hcount;
  const int* htype = for_fast_scoring.htype;

#ifdef DEBUG_SUM_CONTRIBUTIONS
  for (int i = 0; i < matoms; ++i) {
    result += _params.value[atype[i]];
  }
  cerr << "Fast version, sum before H " << result << '\n';
#endif

  result = _params.value[kBias];

  for (int i = 0; i < matoms; ++i) {
    result += _params.value[atype[i]];

    if (hcount[i] == 0) {
      continue;
    }

    const int hti = htype[i];
    // cerr << "H atom " << i << " type " << hti << " count " << hcount[i] << " value " << _params.value[hti] << '\n';

    // cerr << (hcount[i] * _params.value[hti]) << " being added\n";

    if (hcount[i] == 1) {
      result += _params.value[hti];
    } else if (hcount[i] == 2) {
      result += (_params.value[hti] + _params.value[hti]);
    } else {
      result += (hcount[i] * _params.value[hti]);
    }
  }

  return result;
}


}  // namespace alogp


