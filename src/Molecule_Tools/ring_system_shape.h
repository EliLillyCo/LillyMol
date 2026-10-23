#ifndef MOLECULE_TOOLS_RING_SYSTEM_SHAPE_H_
#define MOLECULE_TOOLS_RING_SYSTEM_SHAPE_H_

#include <vector>

#include "Molecule_Lib/molecule.h"

namespace ring_system_shape {

// A bond from a ring-system atom to an atom outside that ring system.
struct RingSystemAttachment {
  atom_number_t ring_atom = INVALID_ATOM_NUMBER;
  atom_number_t outside_atom = INVALID_ATOM_NUMBER;
};

// For one exit point, the farthest ring-system atom(s), using shortest-path
// distances constrained to atoms in the ring system.
struct RingSystemSpan {
  atom_number_t from = INVALID_ATOM_NUMBER;
  int max_separation = -1;
  std::vector<atom_number_t> farthest_atoms;
};

enum class RingSystemShapeClass {
  kNotApplicable,
  kRodLike,
  kNotRodLike,
  kMultiSubstituted,
  kInvalid,
};

struct RingSystemShape {
  std::vector<RingSystemAttachment> attachments;
  std::vector<RingSystemSpan> ring_system_spans;

  // Meaningful only for ring systems with exactly two exit points. The observed
  // separation is the shortest path between the two exit points within the
  // ring system. The rod deficit is the difference between the best possible
  // exit separation and the observed separation.
  int observed_separation = -1;
  int rod_deficit = -1;

  RingSystemShapeClass shape_class = RingSystemShapeClass::kInvalid;
};

// Analyse atoms for which `ring_system_membership[i] == ring_system_id`.
// Non-ring atoms are expected to have a different membership value, often zero.
// Returns false if the requested ring system is absent or internally
// disconnected according to shortest-path traversal within that system.
bool AnalyseRingSystemShape(const Molecule& m, const int* ring_system_membership,
                            int ring_system_id, RingSystemShape& result);

// Isotope-aware variant. When `only_process_isotope` is positive, only ring
// atoms with that isotope can define exit vectors.
bool AnalyseRingSystemShape(const Molecule& m, const int* ring_system_membership,
                            int ring_system_id, isotope_t only_process_isotope,
                            RingSystemShape& result);

// Counts non-ring atoms that appear to be substantial branch points. Atoms in
// any ring system are skipped first. Terminal single-atom decorations and
// terminal multiply bonded atoms, like carbonyl oxygen, are not considered
// evidence of branching.
int NonRingBranchPointCount(const Molecule& m, const int* ring_system_membership);

// Counts ring atoms from which two or more substantial non-ring substituents
// emerge. This is not part of the ring-system geometric classification; it is
// a higher-level shape descriptor useful for deciding whether a molecule is
// globally rod-like.
int RingAtomBranchPointCount(const Molecule& m, const int* ring_system_membership);

// Isotope-aware variant. When `only_process_isotope` is positive, only ring
// atoms with that isotope can contribute ring-atom branch points.
int RingAtomBranchPointCount(const Molecule& m, const int* ring_system_membership,
                             isotope_t only_process_isotope);

}  // namespace ring_system_shape

#endif  // MOLECULE_TOOLS_RING_SYSTEM_SHAPE_H_
