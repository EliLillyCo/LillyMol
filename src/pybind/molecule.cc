// Some functions to mimic RDKit
#include <optional>

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"

std::optional<Molecule>
MolFromSmiles(const std::string& smi) {
  Molecule m;
  if (! m.build_from_smiles(smi)) {
    return std::nullopt;
  }

  // Think of this as being like sanitize=True
  m.compute_aromaticity_if_needed();

  return m;
}

// This does not work because Substructure_Query lacks a move constructor.
#ifdef NEED_MOVE_CONSTRUCTOR___
std::optional<Substructure_Query>
QueryFromSmarts(const std::string& smarts) {
  Substructure_Query result;
  if (! result.create_from_smarts(smarts)) {
    return std::nullopt;
  }

  return result;
}
#endif
