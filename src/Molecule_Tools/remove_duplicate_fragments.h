#ifndef MOLECULE_TOOLS_REMOVE_DUPLICATE_FRAGMENTS_H
#define MOLECULE_TOOLS_REMOVE_DUPLICATE_FRAGMENTS_H

#include "Molecule_Lib/molecule.h"

namespace remove_duplicate_fragments {

// Remove duplicate fragments from `m`. Returns the
// number of atoms removed, and `fragments_removed` will be
// set to the number of fragments removed.
int RemoveDuplicateFragments(Molecule& m, int &fragments_removed);
}  // namespace remove_duplicate_fragments

#endif  // MOLECULE_TOOLS_REMOVE_DUPLICATE_FRAGMENTS_H
