#include <random>

#include "misc.h"

#include "iwreaction.h"

Reaction_Iterator::Reaction_Iterator ()
{
  _number_sidechains = 0;

  _reagents_in_sidechain = NULL;

  _reagent = NULL;

  return;
}

Reaction_Iterator::~Reaction_Iterator ()
{
  if (NULL != _reagents_in_sidechain)
    delete [] _reagents_in_sidechain;

  if (NULL != _reagent)
    delete [] _reagent;

  return;
}

int
Reaction_Iterator::ok () const
{
  return 1;
}

int
Reaction_Iterator::debug_print (std::ostream & os) const
{
  os << "Iterator for " << _number_sidechains << " sidechains\n";
  for (int i = 0; i < _number_sidechains; i++)
  {
    os << _reagents_in_sidechain[i] << " reagents, current " << _reagent[i] << endl;
  }

  return os.good ();
}

int
Reaction_Iterator::initialise (const IWReaction & r)
{
  _number_sidechains = r.number_sidechains ();

  assert (_number_sidechains > 0);

//cerr << "Initialising reaction iterator with " << _number_sidechains << " sidechains\n";

  _reagents_in_sidechain = new int[_number_sidechains];

  _reagent = new_int (_number_sidechains);

  for (int i = 0; i < _number_sidechains; i++)
  {
    const Sidechain_Reaction_Site * s = r.sidechain (i);

    _reagents_in_sidechain[i] = s->number_reagents ();
  }

  _active = 1;

  return 1;
}

std::ostream &
operator << (std::ostream & os, const Reaction_Iterator & r)
{
  r.debug_print (os);

  return os;
}

void
Reaction_Iterator::operator++ (int notused)
{
  assert (active ());

  for (int i = _number_sidechains - 1; i >= 0; i--)
  {
    if (_reagent[i] < _reagents_in_sidechain[i] - 1)
    {
      _reagent[i]++;
      return;
    }

    _reagent[i] = 0;
  }

  _active = 0;

  return;
}

#if (GCC_VERSION >= 40900)
void
Reaction_Iterator::randomise (std::mt19937_64 & rng)
{
  for (int i = 0; i < _number_sidechains; ++i)
  {
    std::uniform_int_distribution<int> u(0, _reagents_in_sidechain[i] = 1);

    const int s = u(rng);

    _reagent[i] = s;
  }

  return;
}
#endif
