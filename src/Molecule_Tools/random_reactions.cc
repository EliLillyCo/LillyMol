#include "random_reactions.h"

Random_Reactions:: Random_Reactions()
{
  _rxn = nullptr;
  _nrxn = 0;

  _iter = nullptr;

   _attempts_made = 0;
   _molecules_formed = 0;

  return;
}

Random_Reactions::~Random_Reactions()
{
  if (nullptr != _rxn)
    delete [] _rxn;

  if (nullptr != _iter)
    delete [] _iter;

  return;
}

template <typename T>
int
Random_Reactions::report (T & output) const
{
  output << "Random_Reactions::report:attempted " << _attempts_made << " reactions, formed " << _molecules_formed << '\n';
  return 1;
}

int
Random_Reactions::build (const Command_Line & cl,
                         char flag,
                         int verbose)
{
  std::random_device rd;
  _rng.seed(rd());

  assert (0 == _nrxn);
  _nrxn = cl.option_count(flag);

  if (0 == _nrxn)
  {
    cerr << "Random_Reactions::build:no files specified\n";
    return 0;
  }

  Sidechain_Match_Conditions smc;
  smc.set_find_unique_embeddings_only(1);
  smc.set_one_embedding_per_start_atom(1);
  smc.set_ignore_symmetry_related_matches(1);
  smc.set_process_hit_number(0);

  _rxn = new IWReaction[_nrxn];
  _iter = new Reaction_Iterator[_nrxn];

  for (auto i = 0; i < _nrxn; ++i)
  {
    const char * fname = cl.option_value(flag, i);
    if (! _rxn[i].do_read(fname, smc))
    {
      cerr << "Random_Reactions::build:cannot read reaction file '" << fname << "'\n";
      return 0;
    }

    _iter[i].initialise(_rxn[i]);
  }

  if (verbose)
    cerr << "Defined " << _nrxn << " reactions\n";

  return _nrxn;
}

int
Random_Reactions::perform_random_reaction (const IWString & smiles)
{
  Molecule m;

  if (! m.build_from_smiles(smiles))
  {
//  cerr << "Random_Reactions::perform_random_reaction:cannot interpret smiles '" << smiles << "'\n";
    return 0;
  }

  return perform_random_reaction (m);
}

int
Random_Reactions::perform_random_reaction (Molecule & m)
{
  _attempts_made++;

  std::uniform_int_distribution<int> ur(0, _nrxn - 1);

  int r = ur(_rng);

  IWReaction & ri = _rxn[r];

  Substructure_Results sresults;

  int nhits = ri.substructure_search(&m, sresults);

  if (0 == nhits)
    return 0;

  const Set_of_Atoms * e;
  if (1 == nhits)
    e = sresults.embedding(0);
  else    // choose one at random
  {
    std::uniform_int_distribution<int> u(0, nhits-1);
    int j = u(_rng);
    e = sresults.embedding(j);
  }

  int ns = ri.number_sidechains();

  if (0 == ns)                  // no sidechains, just a molecular change
    _iter[r]++;
  else
    _iter[r].randomise(_rng);

  if (! _rxn[r].in_place_transformation(m, e, _iter[r]))
    return 0;

  _molecules_formed++;
  return 1;
}
