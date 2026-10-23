#ifndef MOLECULE_LIB_PLANARITY_H_
#define MOLECULE_LIB_PLANARITY_H_

#include <utility>
#include <vector>

#include "Molecule_Lib/molecule.h"

namespace iwplanarity {

enum class PlanarityStatus {
  kPlanar,
  kNonPlanar,
  kError
};

struct PlanarityResult {
  // Status defaults to error.
  PlanarityStatus status = PlanarityStatus::kError;

  // pairs of atom numbers that define the obstructions.
  std::vector<std::pair<atom_number_t, atom_number_t>> obstruction_bonds;
};

// If status is kNonPlanar, then `obstruction_bonds` will be filled.
PlanarityResult Planarity(const Molecule& m);

}  // namespace iwplanarity

#endif // MOLECULE_LIB_PLANARITY_H_
