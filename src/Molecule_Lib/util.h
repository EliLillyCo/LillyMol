
#include <cstdint>
#include <optional>

#include "Foundational/iwstring/iwstring.h"

#include "iwmtypes.h"
#include "istream_and_type.h"
#include "set_of_atoms.h"

class Molecule;

namespace lillymol {
// Given a file name, return the number of molecules in that file.
std::optional<uint64_t> CountMoleculesInFile(const IWString& fname, FileType fype);

// Set entries in `destination` for all atoms in `seeds`, and all atoms within
// `radius` bonds of any seed atom, to `flag`. `destination` must have at least
// m.natoms() entries. Returns the number of entries changed to `flag`.
// A negative radius is a no-op.
int SetAtomsWithinRadius(const Molecule& m, const Set_of_Atoms& seeds, int radius,
                         int flag, int* destination);
}  // namespace lillymol
