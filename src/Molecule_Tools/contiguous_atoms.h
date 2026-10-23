#ifndef MOLECULE_TOOLS_CONTIGUOUS_ATOMS_H_
#define MOLECULE_TOOLS_CONTIGUOUS_ATOMS_H_

namespace contiguous_atoms {

// In iwdescr and MedchemRules we need a measure of the
// maximum size of the largest contiguous set of bonded
// carbon atoms.
int LargestContiguousCarbonGroup(Molecule& m);

// If the caller has a temporary array of at least m.natoms().
int LargestContiguousCarbonGroup(Molecule& m, int* tmp);

};

#endif // MOLECULE_TOOLS_CONTIGUOUS_ATOMS_H_
