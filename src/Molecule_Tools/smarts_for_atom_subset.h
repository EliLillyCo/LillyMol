#ifndef MOLECULE_TOOLS_SMARTS_FOR_ATOM_SUBSET_H_
#define MOLECULE_TOOLS_SMARTS_FOR_ATOM_SUBSET_H_

#include <cstdint>

#include "Foundational/iwstring/iwstring.h"
#include "Molecule_Lib/molecule.h"

namespace lillymol {

// Returns SMARTS for the subgraph induced by `include_atom`.
// `atom_type` is applied to every included atom.
// `include_atom` must contain m.natoms() values.
// Temporarily changes process-global SMILES aromatic-bond output state and is
// therefore not safe to call concurrently.
IWString SmartsForAtomSubset(Molecule& m, uint32_t atom_type,
                             const int* include_atom);
// The smarts generated will start with `start_atom` as the first atom
// in the smarts.
IWString SmartsForAtomSubset(Molecule& m, atom_number_t start_atom, uint32_t atom_type,
                             const int* include_atom);

}  // namespace lillymol

#endif  // MOLECULE_TOOLS_SMARTS_FOR_ATOM_SUBSET_H_
