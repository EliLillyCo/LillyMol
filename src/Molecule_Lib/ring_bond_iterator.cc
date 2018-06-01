#include <stdlib.h>

#include "path.h"

Ring_Bond_Iterator::Ring_Bond_Iterator (const Ring & r)
{
  _atoms_in_ring = r.number_elements ();

  _atom_number = r.rawdata ();

  _ndx = 0;

  return;
}

atom_number_t 
Ring_Bond_Iterator::a2 () const
{
  if (_ndx == _atoms_in_ring - 1)
    return _atom_number[0];

  return _atom_number[_ndx + 1];
}

void
Ring_Bond_Iterator::operator-- (int notused)
{
  _ndx--;

  return;
}

atom_number_t
Ring_Bond_Iterator::previous_atom () const
{
  if (_ndx > 0 && _ndx <= _atoms_in_ring)   // even if exhausted
    return _atom_number[_ndx - 1];

  if (0 == _ndx)
    return _atom_number[_atoms_in_ring - 1];

  return INVALID_ATOM_NUMBER;    // not sure this can happen
}

atom_number_t
Ring_Bond_Iterator::next_atom () const
{
  if (_ndx >= 0 && _ndx < _atoms_in_ring - 1)   // even if exhausted
    return _atom_number[_ndx + 1];

  if (_atoms_in_ring - 1 == _ndx)
    return _atom_number[0];

  return INVALID_ATOM_NUMBER;    // not sure this can happen
}

Ring_Atom_Iterator::Ring_Atom_Iterator (const Ring & r)
{
  _atoms_in_ring = r.number_elements ();

  _atom_number = r.rawdata ();

  _ndx = 0;

  return;
}

atom_number_t 
Ring_Atom_Iterator::next () const
{
  if (_ndx == _atoms_in_ring - 1)
    return _atom_number[0];

  return _atom_number[_ndx + 1];
}

atom_number_t 
Ring_Atom_Iterator::prev () const
{
  if (0 == _ndx)
    return _atom_number[_atoms_in_ring - 1];

  return _atom_number[_ndx - 1];
}

int
Ring_Atom_Iterator::is_next_or_previous (atom_number_t zatom) const
{
// First check previous

  atom_number_t c;
  if (0 == _ndx)
    c = _atom_number[_atoms_in_ring - 1];
  else
    c = _atom_number[_ndx - 1];

  if (zatom == c)
    return 1;

// Then check next
    
  if (_atoms_in_ring - 1 == _ndx)
    c = _atom_number[0];
  else
    c = _atom_number[_ndx + 1];

  return c == zatom;
}

atom_number_t
Ring_Atom_Iterator::move_forward()
{
  _ndx++;
  if (_ndx > _atoms_in_ring - 1)
    _ndx = 0;

  return _atom_number[_ndx];
}

atom_number_t
Ring_Atom_Iterator::move_backward()
{
  _ndx--;
  if (_ndx < 0)
    _ndx = _atoms_in_ring - 1;

  return _atom_number[_ndx];
}
