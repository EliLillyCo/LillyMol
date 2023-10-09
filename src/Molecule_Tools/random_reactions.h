#ifndef RANDOM_REACTIONS_H

/*
  Used by smiles_mutation and random_molecular_permutations
*/

#include <random>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/iwreaction.h"

class Random_Reactions
{
  private:
    int _nrxn;
    IWReaction * _rxn;
    std::mt19937_64 _rng;

    Reaction_Iterator * _iter;

    int _attempts_made;
    int _molecules_formed;

//  private functions

  public:
    Random_Reactions ();
    ~Random_Reactions ();

    int build (const Command_Line & cl, char flag, int verbose);

    int active() const { return _nrxn;}

    template <typename T> int report (T &) const;

    int perform_random_reaction (Molecule & m);
    int perform_random_reaction (const IWString & smiles);
};

#endif
