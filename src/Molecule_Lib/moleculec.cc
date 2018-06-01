/*
  Molecule functions which are dependent on the absence of 
  MOLECULE_BOND_LIST
  The whole file is encased in an ifdef for that.
*/

#ifndef MOLECULE_BOND_LIST

#include <stdlib.h>
#include <assert.h>

#include "molecule.h"

#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

template class resizable_array_p <Connection>;        // instantiate the template
template class resizable_array_base <Connection *>;   // instantiate the template

/*
  Return the atom number of the J'th connection to atom I
*/

atom_number_t
Molecule::other (atom_number_t i, int j) const
{
  assert (OK_ATOM_NUMBER (this, i));
  assert (j >= 0);

  return _things[i]->other (j);
}

bond_type_t
Molecule::btype (atom_number_t i, int j) const
{
  assert (OK_ATOM_NUMBER (this, i));
  assert (j >= 0);

  return _things[i]->btype (j);
}

/*
  This function returns true if atoms a1 and a2 (in molecule m)
  are bonded to each other.
*/

int
Molecule::are_bonded (atom_number_t a1, atom_number_t a2) const
{
  assert (ok_2_atoms(a1, a2));

  const Atom * aa1 = _things[a1];
  return aa1->involves (a2);

  return 0;
}

/*
  Function to return the number of bonds in a molecule.
  In the absence of a bond list, this number must be computed.
*/

int
Molecule::nbonds () const
{
  assert (ok ());

  int bonds = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    const Atom *a = _things[i];
    bonds += a->ncon ();
  }

  assert (bonds == bonds / 2 * 2);

  return bonds / 2;
}

/*
  Function to return the I'th bond in THIS molecule.
  Somewhat of a dilemma here. When compiled for a bond list, the
  bondi can return a pointer to a bond, because the bonds exist within
  the molecule. In this mode, there are no bonds stored.
  We fudge by keeping a statically defined bond which is owned by
  this function. After all, we do return a const pointer (famous
  last words).
*/

const Bond *
Molecule::bondi (int bi) const
{
  assert (ok ());
  assert (bi >= 0);

  if (0 == _number_elements)
    return NULL;

  static Bond *b = NULL;

  DELETE_IF_NOT_NULL (b);

  int bcount = 0;      // internal bond count
  for (int i = 0; i < _number_elements; i++)
  {
    Atom *a = _things[i];
    int ancon = a->ncon ();
    for (int j = 0; j < ancon; j++)
    {
      atom_number_t k = a->other (j);
      if (k > i)
      {
        if (bi == bcount++)
	{
	  b = new Bond (i, k, a->btype (j));
	  return b;
	}
      }
    }
  }

  return NULL;
}

/*
  This function returns the bond type of the J'th connection to atom I in THIS
  molecule.
  When compiled with a bond list, the molecule contains bonds, so it is easy to
  return a const pointer to one. When compiled without a bond list, there are
  no bonds in the molecule. We keep one static bond.
*/

static Bond *b = NULL;

const Bond *
Molecule::bondi (atom_number_t i, int j) const
{
  assert (OK_ATOM_NUMBER (this, i));
  assert (j >= 0);

  const Atom *a = _things[i];
  int ancon = a->ncon ();

  assert (j < ancon);

  atom_number_t a1 = a->other (j);
  bond_type_t   bt = a->btype (j);

  DELETE_IF_NOT_NULL (b);

  b = new Bond (i, a1, bt);

  return b;
}

/*
  Function to add all the atoms of molecule m2 to molecule m1
  There is no attempt to actually join the fragments, we
  just copy over connectivity, geometry and perhaps charges.
*/

int
Molecule::add_molecule (const Molecule *m2)
{
  assert (ok ());
  assert (OK_MOLECULE (m2));

  int m2atoms = m2->_number_elements;

  if (m2atoms <= 0)       // wierd
    return 1;

  int m1atoms = _number_elements;

  resize (m1atoms + m2atoms);

  for (int i = 0; i < m2atoms; i++)
  {
    const Atom *a2 = m2->_things[i];
    Atom *tmp = new Atom (a2);
    Molecule::add (tmp);
    set_charge (i, m2->charge_on_atom (i));

    int ncon = a2->ncon ();
    tmp->resize (ncon);
    for (int j = 0; j < ncon; j++)
    {
      atom_number_t other = a2->other (j);
      bond_type_t   btype = a2->btype (j);
      tmp->add_con (m1atoms + other, btype);
    }
  }

  _set_modified ();

  return 1;
}

/*
  When building a molecule from several other molecules, we frequently
  need to exclude one or more atoms from that process. For now, we
  have two functions to do this, one for one exception, the other for
  two. Perhaps this should be generalised.
*/


int
Molecule::add_molecule_without (const Molecule *m2, atom_number_t except)
{
  assert (ok ());
  assert (OK_ATOM_NUMBER (m2, except));

  int m1atoms = _number_elements;
  int m2atoms = m2->_number_elements;
  resize (m1atoms + m2atoms - 1);

// first add the atoms, then the bonds.

  for (int i = 0; i < m2atoms; i++)
    if (i != except)
    {
      Atom *a = m2->_things[i];
      Atom *tmp = new Atom (a);
      Molecule::add (tmp);
      atom_number_t j = ADJUSTED_ATOM_NUMBER_1 (i, except);
      set_charge (j, m2->charge_on_atom (i));
    }

  for (i = 0; i < m2atoms; i++)
    if (i != except)
    {
      int ncon = m2->ncon (i);

      for (int k = 0; k < ncon; k++)
      {
        atom_number_t m = m2->other (i, k);
        if (m != except)
        {
          atom_number_t ai = ADJUSTED_ATOM_NUMBER_1 (i, except);
          atom_number_t am = ADJUSTED_ATOM_NUMBER_1 (m, except);
          if (am > ai)
            add_bond (ai, am, m2->btype (i, k));
	}
      }
    }

  _set_modified ();
  assert (ok ());

  return 1;
}

int
Molecule::add_molecule_without (const Molecule *m2, atom_number_t except1,
                                atom_number_t except2)
{
  assert (ok ());
  assert (OK_2_ATOMS (m2, except1, except2));

  int m1atoms = _number_elements;
  int m2atoms = m2->_number_elements;
  resize (m1atoms + m2atoms - 2);

// first add the atoms, then the bonds.

  for (int i = 0; i < m2atoms; i++)
    if (i != except1 && i != except2)
    {
      Atom *a = m2->_things[i];
      Atom *tmp = new Atom (a);
      Molecule::add (tmp);
      atom_number_t j = ADJUSTED_ATOM_NUMBER_2 (i, except1, except2);
      set_charge (j, m2->charge_on_atom (i));
    }

  for (i = 0; i < m2atoms; i++)
    if (i != except1 && i != except2)
    {
      int ncon = m2->ncon (i);

      for (int k = 0; k < ncon; k++)
      {
        atom_number_t m = m2->other (i, k);
        if ( (m != except1) && (m != except2))
        {
          atom_number_t ai = _ADJUSTED_ATOM_NUMBER_2 (i, except1, except2);
          atom_number_t am = _ADJUSTED_ATOM_NUMBER_2 (m, except1, except2);
          if (am > ai)
            add_bond (ai, am, m2->btype (i, k));
	}
      }
    }

  _set_modified ();
  assert (ok ());

  return 1;
}

/*
  add_bond can be called with a flag to indicate a partially built molecule.
  In that case, bonds to atoms not yet defined can be added.
  At the end of such addition, the bonding information may not be symmetric.
  This function makes all bonding symmetric.
*/

int
Molecule::finished_bond_addition ()
{
  assert (ok ());
  
  if (! _partially_built)
  {
    cerr << "finished bond addition: complete molecule\n";
    return 1;
  }

  _partially_built = 0;

  if (0 == _number_elements)
    return 1;

  for (int i = 0; i < _number_elements; i++)
  {
    Atom *ai = _things[i];
    int ancon = ai->ncon ();
    for (int j = 0; j < ancon; j++)
    {
      atom_number_t k = ai->other (j);
      assert (OK_ATOM_NUMBER (this, k));
      if (! are_bonded (k, i))
      {
        Atom *ak = _things[k];
        ak->add_con (i, ai->btype (j));
      }
    }
  }

  return 1;
}

/*
  This function returns the number of bonds removed
*/

int
Molecule::remove_bonds_to_atom (atom_number_t zatom, int adjust_atom_numbers)
{
  assert (OK_ATOM_NUMBER (this, zatom));

  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    Atom * a = _things[i];
    rc += a->remove_bonds_to_atom (zatom, adjust_atom_numbers);
  }

  return rc;
}

#endif    // MOLECULE_BOND_LIST
