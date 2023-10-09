#include <stdlib.h>

#include "Foundational/iwmisc/misc.h"
#include "iwrcb.h"

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

using std::cerr;
using std::endl;

/*
  The _present array is used to keep track of which pairs are present.
*/

Ring_Closure_Bonds::Ring_Closure_Bonds()
{
  _atoms_in_molecule = 0;

  _present = nullptr;

  return;
}

Ring_Closure_Bonds::Ring_Closure_Bonds(const Ring_Closure_Bonds & rhs) : _atoms_in_molecule(rhs._atoms_in_molecule)
{
  _atoms_in_molecule = rhs._atoms_in_molecule;

  if (0 == rhs._atoms_in_molecule)
  {
    _present = nullptr;
    return;
  }

  _present = new int[_atoms_in_molecule];

  copy_vector(_present, rhs._present, _atoms_in_molecule);

  if (rhs._number_elements > 0)
    resizable_array<int>::operator = (rhs);
  
  return;
}

Ring_Closure_Bonds::~Ring_Closure_Bonds()
{
  if (nullptr != _present)
    delete [] _present;

  return;
}

int
Ring_Closure_Bonds::ok() const
{
  if (_atoms_in_molecule < 0)
    return 0;

  if (nullptr == _present && _number_elements > 0 && 0 == _atoms_in_molecule)    // empty
    return 1;

  if (nullptr != _present && _atoms_in_molecule > 0 && 0 == _number_elements)    // the ring openings never touch the same atom twice
    return 1;

  if (nullptr != _present && _number_elements > 0 && _atoms_in_molecule > 0)    // in use
    return 1;

  return 0;    // should not happen
}

Ring_Closure_Bonds &
Ring_Closure_Bonds::operator= (const Ring_Closure_Bonds & rhs)
{
  assert (rhs.ok());

  if (nullptr != _present)
  {
    delete [] _present;
    _present = nullptr;
  }

  if (nullptr == rhs._present)
  {
    _atoms_in_molecule = 0;
    resizable_array<int>::resize(0);
    return *this;
  }

  activate(rhs._atoms_in_molecule);

  for (int i = 0; i < rhs._number_elements; i++)
  {
    resizable_array<int>::add(rhs._things[i]);
  }

  return *this;
}

int 
Ring_Closure_Bonds::write_bonds(std::ostream & output) const
{
  output << "Ring_Closure_Bonds::write_bonds:molecule has " << _atoms_in_molecule << " atoms\n";

  for (int i = 0; i < _atoms_in_molecule; i++)
  {
    if (_present[i] < 0)
      continue;

    atom_number_t a1 = _present[i];
    atom_number_t a2 = _present[a1];

    if (a1 < a2)
      output << " ring closure bond between " << a1 << " and " << a2 << endl;
  }

  for (int i = 0; i < _number_elements; i++)
  {
    int j = _things[i];

    atom_number_t a1 = j / _atoms_in_molecule;
    atom_number_t a2 = j % _atoms_in_molecule;

    output << "  ring closure bond between " << a1 << " and " << a2 << endl;
  }

  return output.good();
}

int
Ring_Closure_Bonds::activate(int s)
{
  _atoms_in_molecule = s;

  if (nullptr != _present)
    delete [] _present;

  _present = new_int(_atoms_in_molecule, -1);

  resizable_array<int>::resize_keep_storage(0);

  return 1;
}

int
Ring_Closure_Bonds::reset()
{
  set_vector(_present, _atoms_in_molecule, -1);

  resizable_array<int>::resize_keep_storage(0);

  return 1;
}

void 
Ring_Closure_Bonds::invalidate()
{
  DELETE_IF_NOT_NULL(_present);

  _atoms_in_molecule = -9;

  resizable_array<int>::resize_keep_storage(0);

  return;
}

int
Ring_Closure_Bonds::_form_corresponding_integer(atom_number_t a1,
                                                atom_number_t a2) const
{
  if (a1 < a2)
    return a1 * _atoms_in_molecule + a2;
  else
    return a2 * _atoms_in_molecule + a1;
}

int
Ring_Closure_Bonds::add(atom_number_t a1, atom_number_t a2)
{
  assert (nullptr != _present);

  if (-1 == _present[a1] && -1 == _present[a2])    // neigher one been seen before
  {
    _present[a1] = a2;
    _present[a2] = a1;

    return 1;
  }

// one or both of these atoms has appeared before. We are going to have
// to add a bond. Use -2 to indicate that the atom is present in the bonds

  if (-1 == _present[a1])
    _present[a1] = -2;
  else if (-1 == _present[a2])
    _present[a2] = -2;

  int b = _form_corresponding_integer(a1, a2);

  resizable_array<int>::add(b);

  return _number_elements;
}

int
Ring_Closure_Bonds::contains(atom_number_t a1, atom_number_t a2) const
{
  assert (nullptr != _present);

  if (-1 == _present[a1] || -1 == _present[a2])
    return 0;

  if (a1 == _present[a2] && a2 == _present[a1])
    return 1;

// Need to go look in the bond list

  int b = _form_corresponding_integer(a1, a2);

  return resizable_array<int>::contains(b);
}

int
Ring_Closure_Bonds::is_the_same(const Ring_Closure_Bonds & rhs) const
{
  if (nullptr == _present && nullptr == rhs._present)
    ;
  else if (nullptr == _present && nullptr != rhs._present)
    return 0;
  else if (nullptr != _present && nullptr == rhs._present)
    return 0;
  else
  {
    for (int i = 0; i < _atoms_in_molecule; i++)
    {
      if (_present[i] != rhs._present[i])
        return 0;
    }
  }

  if (_number_elements != rhs._number_elements)
    return 0;

  if (1 == _number_elements)
    return _things[0] == rhs._things[1];

  for (int i = 0; i < _number_elements; i++)
  {
    if (! rhs.resizable_array<int>::contains(_things[i]))
      return 0;
  }

  return 1;
}

int
Ring_Closure_Bonds::report_differences(const Ring_Closure_Bonds & rhs,
                                       std::ostream & output) const
{
  if (nullptr == _present || nullptr == rhs._present)
    return output.good();

  output << "Ring_Closure_Bonds::report_differences:comparing items for " << _atoms_in_molecule << " atoms\n";

  if (_atoms_in_molecule != rhs._atoms_in_molecule)
  {
    output << "Ring_Closure_Bonds:report_differences:atoms in molecule differs, " << _atoms_in_molecule << " and " << rhs._atoms_in_molecule << endl;
    return 0;
  }

  int lhs_bonds = 0;
  int rhs_bonds = 0;
  int differences_found = 0;

  for (int i = 0; i < _atoms_in_molecule; i++)
  {
    if (_present[i] >= 0)
      lhs_bonds++;

    if (rhs._present[i] > 0)
      rhs_bonds++;

    if (_present[i] == rhs._present[i])
      continue;

    output << "Atom " << i << " present value: lhs " << _present[i] << " rhs " << rhs._present[i] << endl;
    differences_found++;
  }

  if (differences_found)
    output << "LHS has " << (lhs_bonds / 2) << " bonds, rhs has " << (rhs_bonds / 2) << " bonds\n";

  for (int i = 0; i < _number_elements; i++)
  {
    int j = _things[i];

    if (rhs.resizable_array<int>::contains(j))
      continue;

    atom_number_t a1 = j / _atoms_in_molecule;
    atom_number_t a2 = j % _atoms_in_molecule;

    output << "Atoms " << a1 << " to " << a2 << " not in RHS\n";
  }

  for (int i = 0; i < rhs._number_elements; i++)
  {
    int j = rhs._things[i];

    if (resizable_array<int>::contains(j))
      continue;

    atom_number_t a1 = j / _atoms_in_molecule;
    atom_number_t a2 = j % _atoms_in_molecule;

    output << "Atoms " << a1 << " to " << a2 << " not in LHS\n";
  }

  return output.good();
}

/*
  When we re-derive the ring closure bonds in a subset, we want to check to
  ensure that all the bonds we have are also in the parent set
*/

int
Ring_Closure_Bonds::is_subset_of(const Ring_Closure_Bonds & parent) const
{
  assert (_atoms_in_molecule == parent._atoms_in_molecule);

  if (nullptr == _present)    // should not happen
    return 0;

  int rc = 1;

  for (int i = 0; i < _atoms_in_molecule; i++)
  {
    if (_present[i] < 0)
      continue;

    atom_number_t a1 = _present[i];
    atom_number_t a2 = _present[a1];

    if (a2 < a1)    // don't check things twice
      continue;

    if (parent.contains(a1, a2))
      continue;

    cerr << "Ring_Closure_Bonds::is_subset_of:pair " << a1 << " to " << a2 << " not part of parent\n";
    rc = 0;
  }

  for (int i = 0; i < _number_elements; i++)
  {
    int j = _things[i];

    if (parent.resizable_array<int>::contains(j))
      continue;

    atom_number_t a1 = j / _atoms_in_molecule;
    atom_number_t a2 = j % _atoms_in_molecule;
    cerr << "Ring_Closure_Bonds::is_subset_of:pair " << a1 << " to " << a2 << " not part of parent\n";
    rc = 0;
  }

  return rc;
}
