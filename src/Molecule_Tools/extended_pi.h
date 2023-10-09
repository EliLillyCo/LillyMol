#include "Molecule_Lib/molecule.h"

namespace radius_pi_extension {
struct RadiusPiExtensionParams {
  // Include all atoms out to this radius regardless.
  int radius = 3;

  // Extend through aromatic bonds this far beyond `radius`.
  int aromatic_bond_extend = 0;

  // Extend conjugated bonds this far beyond `radius`.
  int conjugated_bond_extend = 0;

  // Atoms that are included because of being within `radius`
  // get assigned this number.
  int in_radius_marker = 1;

  // By default, each layer gets a different number. Breaks if
  // radius > 100.
  int in_radius_step = 1;

  int aromatic_extension_marker = 100;
  int conjugated_extension_marker = 200;

};

// Mark entries in `in_system` with their relationship to the starting
// entries in that array.
// First all atoms within params.radius bonds of initially marked atoms
// are included.
// Then, if params.aromatic_bond_extend is set, aromatic atoms are
// included.
// Then, if params.conjugated_bond_extend is set, conjugated atoms are
// included.
// Note, there is an arbitrary order dependence.
//   Aromatic extensions are perceived first, and then conjugated.
//   Once an atom is marked as visited, it will not be visited again,
//   so if an atom is not visited via an aromatic, it will not be found
//   in the conjugated search, even though it would have been found via
//   a separated conjugated search.
int
RadiusPiExtension(Molecule& m, int * in_system, const RadiusPiExtensionParams& params);

}  // namespace radius_pi_extension
