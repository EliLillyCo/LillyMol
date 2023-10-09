#ifndef MOLECULE_TOOLS_TOPOLOGICAL_TORSION_H
#define MOLECULE_TOOLS_TOPOLOGICAL_TORSION_H
// Generator for topological torsions.

#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Molecule_Lib/molecule.h"

namespace topological_torsion {

struct TorsionOptions {
  // If the 4th atom in a torsion is the starting atom, we
  // can optionally fingerprint that path.
  bool fingerprint_3_membered_rings = false;

  // If there is a bond between the first and last member
  // of a path, we can optionally fingerprint that path.
  bool fingerprint_4_membered_rings = false;

  bool use_bond_orders = true;
};

Sparse_Fingerprint_Creator
TopologicalTorsion(Molecule& m, 
                   const atom_type_t* atom_type,
                   const int * include_atom,
                   TorsionOptions& options);

}  // namespace topological_torsion

#endif // MOLECULE_TOOLS_TOPOLOGICAL_TORSION_H
