#include <stdlib.h>

#include "iw_stl_hash_set.h"

#include "misc.h"

#include "molecule.h"

#include "do_remove_duplicate_fragments.h"

/*
*/

static int
do_remove_duplicate_fragments (Molecule & m,
                               const int * fragment_membership,
                               int * number_times_fragment_this_size_encountered,
                               int & fragments_removed)
{
  int nf = m.number_fragments ();

  int matoms = m.natoms ();

  IWString_STL_Hash_Set usmi;

  Set_of_Atoms atoms_to_be_deleted;
  atoms_to_be_deleted.resize (matoms);

  for (int i = 0; i < nf; i++)
  {
    int aif = m.atoms_in_fragment (i);

    number_times_fragment_this_size_encountered[aif]++;

    if (1 == number_times_fragment_this_size_encountered[aif])     // first time we have seen a fragment with this many atoms. Therefore it is not a duplicate
      continue;

    Molecule tmp;

    m.create_subset (tmp, fragment_membership, i);

    const IWString & smiles = tmp.unique_smiles ();

    if (usmi.contains (smiles))
    {
      for (int j = 0; j < matoms; j++)
      {
        if (i == fragment_membership[j])
          atoms_to_be_deleted.add (j);
      }
      fragments_removed++;
    }
    else
      usmi.insert (smiles);
  }

  if (0 == atoms_to_be_deleted.number_elements ())
    return 1;

  return m.remove_atoms (atoms_to_be_deleted);
}

int
do_remove_duplicate_fragments (Molecule & m, int & fragments_removed)
{
  fragments_removed = 0;

  int nf = m.number_fragments ();

  if (nf < 2)
  {
    cerr << "do_remove_duplicate_fragments: only " << nf << " fragments in '" << m.name () << "'\n";
    return 1;
  }

  int matoms = m.natoms ();

  int * fragment_membership = new int[matoms];
  int * fragment_size = new_int (matoms);

  m.fragment_membership (fragment_membership);

  int rc = do_remove_duplicate_fragments (m, fragment_membership, fragment_size, fragments_removed);

  delete fragment_membership;
  delete fragment_size;

  return rc;
}
