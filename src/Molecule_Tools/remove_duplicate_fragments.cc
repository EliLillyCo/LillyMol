#include <iostream>
#include <memory>

#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iwmisc/misc.h"

#include "remove_duplicate_fragments.h"

namespace remove_duplicate_fragments {

using std::cerr;
using std::endl;

int
RemoveDuplicateFragments(Molecule& m,
                         const int* fragment_membership,
                         int& fragments_removed)
{
  const int nf = m.number_fragments();

  const int matoms = m.natoms();

  IWString_STL_Hash_Set usmi;

  Set_of_Atoms atoms_to_be_deleted;
  atoms_to_be_deleted.resize(matoms);

  for (int i = 0; i < nf; i++) {
    Molecule tmp;
    m.create_subset(tmp, fragment_membership, i);

    const IWString & smiles = tmp.unique_smiles();

    if (usmi.contains(smiles)) {
      for (int j = 0; j < matoms; j++) {
        if (i == fragment_membership[j])
          atoms_to_be_deleted.add(j);
      }
      fragments_removed++;
    } else {
      usmi.insert(smiles);
    }
  }

  if (atoms_to_be_deleted.empty())
    return 0;

  return m.remove_atoms(atoms_to_be_deleted);
}

int
RemoveDuplicateFragments(Molecule & m, int & fragments_removed)
{
  fragments_removed = 0;

  if (m.number_fragments() < 2) {
    return 0;
  }

  std::unique_ptr<int[]> fragment_membership = m.fragment_membership();

  return RemoveDuplicateFragments(m, fragment_membership.get(), fragments_removed);
}

}  // namespace remove_duplicate_fragments
