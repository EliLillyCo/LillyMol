#include "smarts_for_atom_subset.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/smiles.h"

namespace lillymol {
namespace {

class ResetBondAromaticity {
 public:
  ResetBondAromaticity() {
    set_include_bond_aromaticity_in_smiles(1);
  }

  ~ResetBondAromaticity() {
    set_include_bond_aromaticity_in_smiles(0);
  }
};

}  // namespace

IWString
SmartsForAtomSubset(Molecule& m, atom_number_t start_atom,
                    uint32_t atom_type, const int* include_atom) {
  if (m.empty() || include_atom == nullptr) {
    return IWString();
  }

  // Should we issue a warning, or perhaps return a '*'?
  if (atom_type == 0) {
    return IWString();
  }

  m.compute_aromaticity_if_needed();
  (void) m.ring_membership();

  Smiles_Information smiles_information(m.natoms());
  smiles_information.set_smiles_is_smarts(1);
  smiles_information.prepare_to_build_ordering(m.natoms());
  smiles_information.set_create_smarts_embedding(1);

  int included_atoms = 0;
  for (atom_number_t i = 0; i < m.natoms(); ++i) {
    if (! include_atom[i]) {
      continue;
    }

    ++included_atoms;
    smiles_information.set_user_specified_atomic_smarts(
        i, SmartsForAtomType(m, i, atom_type));
  }

  if (included_atoms == 0) {
    return IWString();
  }

  ResetBondAromaticity reset_bond_aromaticity;
  if (start_atom == kInvalidAtomNumber) {
    return IWString(m.smiles(smiles_information, include_atom));
  } else {
    return IWString(m.smiles_starting_with_atom(start_atom, smiles_information, include_atom));
  }
}

IWString
SmartsForAtomSubset(Molecule& m, uint32_t atom_type,
                    const int* include_atom) {
  return SmartsForAtomSubset(m, kInvalidAtomNumber, atom_type, include_atom);
}

}  // namespace lillymol
