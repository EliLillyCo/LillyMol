#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "misc.h"

#include "molecule.h"
#include "path_scoring.h"
#include "chiral_centre.h"
#include "is_actually_chiral.h"

static int max_iterations = 0;

void
set_max_iterations (int m)
{
  assert (m > 0);

  max_iterations = m;
}

static int allow_unsaturated_atoms_to_be_chiral = 0;

void
set_allow_unsaturated_atoms_to_be_chiral (int s)
{
  allow_unsaturated_atoms_to_be_chiral = s;
}

/*
  To determine if an atom is chiral or not, we need to perform path tracing
  from that atom.
*/

static int
is_actually_chiral (Molecule & m,
                    atom_number_t zatom,
                    resizable_array_p<Path_Scoring> & ps,
                    int * claimed,
                    Atom * const * atom)
{
  const Atom * a = atom[zatom];

  int acon = a->ncon();

  if (ps.number_elements())
    ps.resize_keep_storage(0);

  ps.resize(acon);

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);

    Path_Scoring * p = new Path_Scoring;

    const Atom * aj = atom[j];

    p->initialise(j, aj);

    if (b->is_single_bond())
      p->set_first_bond(1);
    else if (! allow_unsaturated_atoms_to_be_chiral || a->atomic_number() < 14)
    {
      delete p;
      return 0;
    }
    else if (b->is_double_bond())    // tetrahedral Sulphur types
      p->set_first_bond(2);

    claimed[j] = 1;

    ps.add(p);
  }

  int stopped;
  if (resolved(ps, stopped))
    return 1;

  int iterations = 0;

  while (1)
  {
    for (int i = 0; i < acon; i++)
    {
      if (ps[i]->active())
        ps[i]->advance(atom, claimed);
    }

    int stopped;
    if (resolved(ps, stopped))
      return 1;

    if (stopped)      // not resolved, but cannot go any further
      return 0;

    int number_active = 0;
    for (int i = 0; i < acon; i++)
    {
      if (! ps[i]->active())
        continue;

      ps[i]->update_claimed(claimed);
      number_active++;
    }

    if (number_active < 2)
      return 0;

    iterations++;

    if (max_iterations > 0 && iterations >= max_iterations)
    {
      cerr << "Not resolved by " << max_iterations << " iterations\n";
      return 0;
    }
  }

  return 1;
}

/*
  The query for an asymmetric carbon atom will hit things like the
  carbon in t-butyl. We need to examine the neighbours to make sure
  that this atom actually is an asymmetric centre

  Thought about putting in a more aggressive check on the number of
  connections, but too dangerous. Even 2 == ncon is problematic because
  you could have an atom with a lone-pair and an implicit Hydrogen
*/

int
is_actually_chiral (Molecule & m,
                    atom_number_t zatom)
{
  resizable_array_p<Path_Scoring> ps;

  return is_actually_chiral(m, zatom, ps);
}

int
is_actually_chiral (Molecule & m,
                    atom_number_t zatom,
                    resizable_array_p<Path_Scoring> & ps)
{
  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  if (acon < 2 || acon > 4)
    return 0;

  const int hcount = m.hcount(zatom);

  if (hcount > 1)    // what if isotopic Hydrogen???
    return 0;

  int lp;
  if (m.lone_pair_count(zatom, lp) && lp > 1)
    return 0;

  if (1 == hcount && 1 == lp && 7 == a->atomic_number())    // never
    return 0;

  if (acon < a->nbonds() && ! allow_unsaturated_atoms_to_be_chiral)
    return 0;

  m.compute_aromaticity_if_needed();    // so bonds get aromatic character

  if (m.is_aromatic(zatom))
    return 0;

  const int matoms = m.natoms();

  int * claimed = new_int(matoms); std::unique_ptr<int[]> free_claimed(claimed);

  claimed[zatom] = 1;

  Atom * const * atoms = new Atom *[matoms]; std::unique_ptr<Atom * const []> free_atoms(atoms);

  m.atoms( (const Atom **) atoms);

//cerr << "Detailed calculation on " << m.smarts_equivalent_for_atom(zatom) << endl;

  return is_actually_chiral(m, zatom, ps, claimed, atoms);
}

int
do_remove_invalid_chiral_centres (Molecule & m)
{
  int nc = m.chiral_centres();
  if (0 == nc)
    return 0;

// Removing a chiral centre while we are scanning the set would mess things up,
// so we make a list of the atoms with invalid chiral centres and remove them later

  Set_of_Atoms centres_to_be_removed;

  for (int i = 0; i < nc; i++)
  {
    Chiral_Centre * c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);

    atom_number_t a = c->a();

    if (! is_actually_chiral(m, a))
      centres_to_be_removed.add(a);
  }

  if (centres_to_be_removed.number_elements())
  {
    for (int i = 0; i < centres_to_be_removed.number_elements(); i++)
    {
      m.remove_chiral_centre_at_atom(centres_to_be_removed[i]);
    }
  }

  return nc;
}
