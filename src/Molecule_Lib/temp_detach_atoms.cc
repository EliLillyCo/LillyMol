#include <stdlib.h>

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "temp_detach_atoms.h"
#include "chiral_centre.h"

#include "molecule.h"

Temp_Detach_Atoms::Temp_Detach_Atoms()
{
  _active = 1;

  _remove_hydrogens_no_longer_needed = 1;

  _matoms = 0;

  _connection = NULL;

  _bt = SINGLE_BOND;

  return;
}

Temp_Detach_Atoms::~Temp_Detach_Atoms()
{
  if (NULL != _connection)
    delete [] _connection;

  return;
}

int
Temp_Detach_Atoms::recognised_directive(const const_IWSubstring & token)
{
  if ("noremove" == token)
  {
    _remove_hydrogens_no_longer_needed = 0;
    return 1;
  }
  
  if ("nodetach" == token)
  {
    _active = 0;
    return 1;
  }

  return 0;
}

void
Temp_Detach_Atoms::do_not_reattach_to_atom(atom_number_t a)
{
  assert (a >= 0 && a < _matoms);
  assert (NULL != _connection );

  _connection[a] = -1;

  return;
}

/*
  Break any bond between a singly bonded Z and its attached atom. Record
  the identify of the neighbouring atom in the _CONNECTION array
*/

int
Temp_Detach_Atoms::detach_atoms(Molecule & m, atomic_number_t z)
{
  int matoms = m.natoms();

  if (0 == matoms)
    return 0;

  _chiral_centre.resize_keep_storage(0);

  if (matoms > _matoms)
  {
    if (NULL != _connection)
      delete [] _connection;

    _connection = new int[matoms + matoms];     // just a little easier if we have two arrays
  }

  _matoms = matoms;

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (z == m.atomic_number(i) && 1 == m.ncon(i))
    {
      atom_number_t o = m.other(i, 0);

      _connection[i] = o;

      rc++;

//    m.remove_bond_between_atoms(i, o);

      const Chiral_Centre * c = m.chiral_centre_at_atom(o);
      if (NULL == c)
        continue;

      _chiral_centre.add(m.remove_no_delete_chiral_centre_at_atom(o));
    }
    else
      _connection[i] = INVALID_ATOM_NUMBER;
  }

  _need_to_reattach = rc;

  if (0 == rc)
    return rc;

  int * rmb = _connection + matoms;     // remember, we over-allocated it

  for (int i = 0; i < matoms; ++i)
  {
    if (INVALID_ATOM_NUMBER != _connection[i])
      rmb[i] = 1;
    else
      rmb[i] = 0;
  }

  m.remove_bonds_involving_these_atoms(rmb);

  return rc;
}

int
Temp_Detach_Atoms::reattach_atoms(Molecule & m)
{
  assert (_matoms == m.natoms());

  if (0 == _need_to_reattach)
    return 1;      // don't need to do anything

  Set_of_Atoms atoms_to_be_removed;
  atoms_to_be_removed.resize(_matoms);

  for (int i = 0; i < _matoms; i++)
  {
    atom_number_t o = _connection[i];

    if (o < 0)
      continue;

    int known;

    if (m.implicit_hydrogens(o, known))
    {
      m.add_bond(i, o, _bt);
      if (known)
      {
        m.set_implicit_hydrogens_known(o, 0);
        m.unset_all_implicit_hydrogen_information(o);
      }
    }
    else
      atoms_to_be_removed.add(i);
  }

//cerr << "Temp_Detach_Atoms::reattach_atoms:replacing " << _chiral_centre.number_elements() << " chiral centres\n";

  for (int i = 0; i < _chiral_centre.number_elements(); ++i)
  {
    Chiral_Centre * c = _chiral_centre[i];
    if (! m.valid_chiral_centre(c))
      delete c;
    else
      m.add_chiral_centre(c, 1);
  }

  _chiral_centre.resize_keep_storage(0);

  if (_remove_hydrogens_no_longer_needed)
    m.remove_atoms(atoms_to_be_removed);

  return 1;
}

