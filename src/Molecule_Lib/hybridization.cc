#include "Molecule_Lib/hybridization.h"

#include "Molecule_Lib/atom.h"
#include "Molecule_Lib/bond.h"
#include "Molecule_Lib/molecule.h"

namespace {

bool
IsCarbonylLike(Molecule& m, atom_number_t carbon, atom_number_t exclude) {
  if (m.atomic_number(carbon) != 6) {
    return false;
  }

  const Atom* atom = m.atomi(carbon);
  for (const Bond* bond : *atom) {
    if (! bond->is_double_bond()) {
      continue;
    }

    const atom_number_t other = bond->other(carbon);
    if (other == exclude) {
      continue;
    }

    switch (m.atomic_number(other)) {
      case 7:
      case 8:
      case 16:
        return true;
      default:
        break;
    }
  }

  return false;
}

bool
IsAttachedToPiSystem(Molecule& m, atom_number_t zatom) {
  const Atom* atom = m.atomi(zatom);
  for (const Bond* bond : *atom) {
    if (! bond->is_single_bond()) {
      continue;
    }

    const atom_number_t other = bond->other(zatom);
    if (m.is_aromatic(other)) {
      return true;
    }

    const atomic_number_t other_z = m.atomic_number(other);
    if (other_z != 6 && ! (m.atomic_number(zatom) == 8 && other_z == 7)) {
      continue;
    }

    const Atom* other_atom = m.atomi(other);
    for (const Bond* bond2 : *other_atom) {
      if (bond2->other(other) == zatom) {
        continue;
      }
      if (bond2->is_double_bond() || bond2->is_triple_bond() || bond2->is_aromatic()) {
        return true;
      }
    }
  }

  return false;
}

bool
IsAttachedToCarbonyl(Molecule& m, atom_number_t zatom) {
  const atomic_number_t z = m.atomic_number(zatom);
  if (z != 7 && z != 8) {
    return false;
  }

  const Atom* atom = m.atomi(zatom);
  for (const Bond* bond : *atom) {
    if (! bond->is_single_bond()) {
      continue;
    }

    const atom_number_t other = bond->other(zatom);
    if (IsCarbonylLike(m, other, zatom)) {
      return true;
    }
  }

  return false;
}

bool
IsSulfonamideNitrogen(Molecule& m, atom_number_t nitrogen) {
  if (m.atomic_number(nitrogen) != 7) {
    return false;
  }

  const Atom* atom = m.atomi(nitrogen);
  for (const Bond* bond : *atom) {
    if (! bond->is_single_bond()) {
      continue;
    }

    const atom_number_t sulphur = bond->other(nitrogen);
    if (m.atomic_number(sulphur) != 16) {
      continue;
    }

    int doubly_bonded_oxygen = 0;
    const Atom* sulphur_atom = m.atomi(sulphur);
    for (const Bond* sulfur_bond : *sulphur_atom) {
      if (! sulfur_bond->is_double_bond()) {
        continue;
      }
      const atom_number_t other = sulfur_bond->other(sulphur);
      if (m.atomic_number(other) == 8) {
        ++doubly_bonded_oxygen;
      }
    }

    if (doubly_bonded_oxygen >= 2) {
      return true;
    }
  }

  return false;
}

bool
IsTwoConnectedSulfonamideNitrogenAttachedToAromaticCarbon(Molecule& m, atom_number_t nitrogen) {
  if (! IsSulfonamideNitrogen(m, nitrogen)) {
    return false;
  }

  const Atom* atom = m.atomi(nitrogen);
  if (atom->ncon() != 2) {
    return false;
  }

  for (const Bond* bond : *atom) {
    if (! bond->is_single_bond()) {
      continue;
    }
    const atom_number_t other = bond->other(nitrogen);
    if (m.atomic_number(other) == 6 && m.is_aromatic(other)) {
      return true;
    }
  }

  return false;
}

bool
HasTripleBond(const Atom& atom) {
  for (const Bond* bond : atom) {
    if (bond->is_triple_bond()) {
      return true;
    }
  }

  return false;
}

bool
HasDoubleBond(const Atom& atom) {
  for (const Bond* bond : atom) {
    if (bond->is_double_bond()) {
      return true;
    }
  }

  return false;
}

}  // namespace

Hybridization
HybridizationState(Molecule& m, atom_number_t zatom) {
  if (! m.ok_atom_number(zatom)) {
    return Hybridization::kUnspecified;
  }

  const atomic_number_t z = m.atomic_number(zatom);
  if (z == 0) {
    return Hybridization::kUnspecified;
  }

  if (z == 1) {
    return Hybridization::kS;
  }

  const Atom* atom = m.atomi(zatom);
  const int ncon = atom->ncon();
  if (ncon == 0) {
    return Hybridization::kUnspecified;
  }

  if (m.is_aromatic(zatom)) {
    return Hybridization::kSp2;
  }

  if (z == 7 && ncon == 4 && atom->formal_charge() > 0) {
    return Hybridization::kSp3;
  }

  // A two-connected sulfonamide nitrogen attached to an aromatic carbon is
  // reported as SP2 by RDKit. Check this before the general sulfonamide rule.
  if (IsTwoConnectedSulfonamideNitrogenAttachedToAromaticCarbon(m, zatom)) {
    return Hybridization::kSp2;
  }

  // Sulfonamide nitrogen should win over other amide-like environments. Some
  // nitrogens are attached to both SO2 and C=O groups.
  if (IsSulfonamideNitrogen(m, zatom)) {
    return Hybridization::kSp3;
  }

  // RDKit reports common phosphates, sulfones, sulfonamides, and other
  // hypervalent sulfur forms as SP3 even though the connection table can
  // contain P=O or S=O double bonds.
  if (((z == 15 || z == 16) && ncon == 4) ||
      (z == 16 && ncon == 3 && (HasDoubleBond(*atom) || HasTripleBond(*atom)))) {
    return Hybridization::kSp3;
  }

  if (HasTripleBond(*atom)) {
    return Hybridization::kSp;
  }

  if (HasDoubleBond(*atom)) {
    return Hybridization::kSp2;
  }

  if (IsAttachedToCarbonyl(m, zatom)) {
    return Hybridization::kSp2;
  }

  if ((z == 7 || z == 8) && IsAttachedToPiSystem(m, zatom)) {
    return Hybridization::kSp2;
  }

  if (ncon <= 1) {
    return Hybridization::kSp3;
  }

  if (ncon == 2) {
    return Hybridization::kSp3;
  }

  if (ncon == 3) {
    return Hybridization::kSp3;
  }

  if (ncon == 4) {
    return Hybridization::kSp3;
  }

  if (ncon == 5) {
    return Hybridization::kSp3d;
  }

  if (ncon == 6) {
    return Hybridization::kSp3d2;
  }

  return Hybridization::kOther;
}

const char*
ToString(Hybridization hybridization) {
  switch (hybridization) {
    case Hybridization::kUnspecified:
      return "UNSPECIFIED";
    case Hybridization::kS:
      return "S";
    case Hybridization::kSp:
      return "SP";
    case Hybridization::kSp2:
      return "SP2";
    case Hybridization::kSp3:
      return "SP3";
    case Hybridization::kSp2d:
      return "SP2D";
    case Hybridization::kSp3d:
      return "SP3D";
    case Hybridization::kSp3d2:
      return "SP3D2";
    case Hybridization::kOther:
      return "OTHER";
  }

  return "OTHER";
}
