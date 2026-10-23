#ifndef MOLECULE_LIB_HYBRIDIZATION_H_
#define MOLECULE_LIB_HYBRIDIZATION_H_

#include "Molecule_Lib/iwmtypes.h"

class Molecule;

enum class Hybridization {
  kUnspecified = 0,
  kS,
  kSp,
  kSp2,
  kSp3,
  kSp2d,
  kSp3d,
  kSp3d2,
  kOther,
};

// Return an RDKit-like hybridization assignment for `zatom` in `m`.
// This is computed on demand and not stored in the Molecule or Atom.
Hybridization HybridizationState(Molecule& m, atom_number_t zatom);

const char* ToString(Hybridization hybridization);

#endif  // MOLECULE_LIB_HYBRIDIZATION_H_
