#ifndef MOLECULE_TOOLS_LADDERANE_H_
#define MOLECULE_TOOLS_LADDERANE_H_

#include "Molecule_Lib/molecule.h"

namespace ladderane {

// Return the largest number number of fused cyclobutanes in a ring system.
// Does NOT hit cubane, or patterns where the fused rings are not in a line.
int
CountLadderane(Molecule& m);

// Returns the largest number of spiro-fused cyclopropanes in a molecule.
int
CountPolySpiroCycloPropane(Molecule& m);

}  // namespace ladderane

#endif  // MOLECULE_TOOLS_LADDERANE_H_
