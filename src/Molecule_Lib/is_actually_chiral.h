#ifndef MOLECULE_LIB_IS_ACTUALLY_CHIRAL_H_
#define MOLECULE_LIB_IS_ACTUALLY_CHIRAL_H_

#include "iwmtypes.h"
#include "path_scoring.h"

class Molecule;

// This is not entirely robust. It might be safer, but more expensive
// to do a full symmetry determination.
extern int is_actually_chiral (Molecule & m, atom_number_t);
extern int is_actually_chiral (Molecule & m, atom_number_t, resizable_array_p<Path_Scoring> &);

namespace lillymol {

// Returns true of `zatom` is asymmetric.
int IsActuallyChiralBySymmetry(Molecule& m, atom_number_t zatom);

// Remove invalid chiral centres.

// Returns the number of invalid chiral centres removed
// This version works by using Path_Scoring objects, which might be
// less expensive than using symmetry, but which might not be as accurate.

int do_remove_invalid_chiral_centres (Molecule & m);

// Using symmetry to identify bad chiral centres is more expensive, but
// more reliable.

int RemoveInvalidChiralCentresUsingSymmetry(Molecule & m);

// The states returned by GetChiralityStatus
enum class ChiralStatus {
  // Not marked and not chiral.
  kNotChiral,
  // Marked and invalid.
  kChiral,
  // Marked but not valid.
  kInvalid,
  // Not marked but actually chiral
  kUnmarked
};

std::unique_ptr<ChiralStatus[]> ChiralityStatus(Molecule& m);

}  // namespace lillymol

extern void set_allow_unsaturated_atoms_to_be_chiral (int s);

// This has never been used and is a bad idea.
//extern void set_max_iterations (int m);

#endif  // MOLECULE_LIB_IS_ACTUALLY_CHIRAL_H_
