#ifndef REACTION_DUPLICATE_H
#define REACTION_DUPLICATE_H

/*
  In a couple of the programmes I need to suppress duplicate products from reactions.
  Since these are likely to be generated in a systematic form, I first lookup the smiles
  of the molecule and then the unique smiles. There are counters here to see if that
  is saving time or not.
*/

#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/molecule.h"

class Reaction_Duplicate
{
  private:
    IW_STL_Hash_Set _smiles;
    IW_STL_Hash_Set _unique_smiles;

    int _molecules_submitted;
    int _unique_molecules_encountered;
    int _detected_via_standard_smiles;

  public:
    Reaction_Duplicate ();

    int is_duplicate(Molecule &);

    int report (std::ostream &) const;
};

#endif
