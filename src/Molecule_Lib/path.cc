#include <stdlib.h>
#include <assert.h>

// define this symbol to get the private molecule functions below

#define COMPILING_MOLECULER_CC

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "iwbits.h"

#include "molecule.h"
#include "path.h"

#define RING_FRAGMENT_MEMBERSHIP_UNKNOWN -8

Ring::Ring ()
{
  _is_fused = 0;
  _fused_system_identifier = -1;
  _aromaticity = AROMATICITY_NOT_DETERMINED;
  _ring_number = -9;
  _largest_number_of_bonds_shared_with_another_ring = 0;
  _fragment_membership = RING_FRAGMENT_MEMBERSHIP_UNKNOWN;
  _number_strongly_fused_neighbours = 0;
}

Ring::Ring (const Ring & rhs)
{
  _is_fused = rhs._is_fused;
  _fused_system_identifier = rhs._fused_system_identifier;
  _aromaticity = rhs._aromaticity;
  _ring_number = rhs._ring_number;
  _largest_number_of_bonds_shared_with_another_ring = rhs._largest_number_of_bonds_shared_with_another_ring;
  _fragment_membership = rhs._fragment_membership;
  _number_strongly_fused_neighbours = rhs._number_strongly_fused_neighbours;

  int n = rhs.number_elements();

  resize(n);

  for (int i = 0; i < n; i++)
  {
    add(rhs._things[i]);
  }

  return;
}

/*
  Tried to write the OK function more carefully, but ran into problems
  with the way we build the fusion information, so it is left simple.
*/

int
Ring::ok () const
{
  return resizable_array<atom_number_t>::ok ();
}

int
Ring::set_fused_to (Ring * r, int bic)
{
  assert (ok ());
  assert (r->ok ());

  _fused_neighbours.add_if_not_already_present (r);

  if (bic > _largest_number_of_bonds_shared_with_another_ring)
    _largest_number_of_bonds_shared_with_another_ring = bic;
  if (bic > 1)
    _number_strongly_fused_neighbours++;

  return 1;
}

const Ring *
Ring::fused_neighbour (int n) const
{
  assert (ok ());
  
  return _fused_neighbours[n];
}

int
Ring::propagate_fused_system_identifier (int fsid)
{
  _fused_system_identifier = fsid;
  _is_fused = 1;

  int rc = 1;

  int nf = _fused_neighbours.number_elements ();
  assert (nf > 0);

  for (int i = 0; i < nf; i++)
  {
    Ring * ri = _fused_neighbours[i];
    assert (ri->ok ());

    if (! ri->is_fused ())
      rc += ri->propagate_fused_system_identifier (fsid);
  }

  return rc;
}

int
Ring::set_aromaticity (aromaticity_type_t arom)
{
  assert (AROMATIC == arom || NOT_AROMATIC == arom);

  _aromaticity = arom;

  return 1;
}

/*
  When computing ring membership, we frequently need to increment values in
  an array for each atom in a ring. This provides an efficient way of doing
  that.

  Note that the values in ring_membership are probably initialised with
  a negative value, so the argument FLOOR (which defaults to 0) is used
  to handle that case.
*/

int 
Ring::update_ring_membership (int * ring_membership, int increment,
                              int floor) const
{
  assert (NULL != ring_membership);

  for (int i = 0; i < _number_elements; i++)
  {
    atom_number_t j = _things[i];
    if (ring_membership[j] < floor)
      ring_membership[j] = floor + increment;
    else
      ring_membership[j] += increment;
  }

  return _number_elements;
}

/*
  We need to know if an isolated ring is spiro fused. 
  We look at ring_membership and see if any of the atoms hit are in more
  than 1 rings
*/

int 
Ring::spiro_fused (const int * ring_membership) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    atom_number_t j = _things[i];

    if (ring_membership[j] > 1)
      return 1;
  }

  return 0;   // nope, all ring_membership items must have been 0 or 1
}

/*
  The isolated rings have all called isolated_ring_update_ring_membership,
  and have incremented ring_membership. Now we look at the fused rings
  and see if they touch any atoms that have been set in ring_membership
*/

int
Ring::fused_ring_check_for_spiro_fusion (const int * ring_membership) const
{
  return any_members_set_in_array(ring_membership);

/*for (int i = 0; i < _number_elements; i++)
  {
    atom_number_t j = _things[i];

    if (ring_membership[j] > 0)
      return 1;
  }*/

  return 0;
}

std::ostream &
operator << (std::ostream & os, const Ring & r)
{
  os << "Ring " << r._ring_number;
  if (r.is_fused() || r.fused_system_identifier () >= 0 || r.fused_ring_neighbours ())
  {
    os << " (";
    if (r.is_fused())
      os << "FSysId " << r.fused_system_identifier ();
    if (r.fused_ring_neighbours ())
      os << " fused " << r.fused_ring_neighbours ();
    if (r.strongly_fused_ring_neighbours () > 0)
      os << ", " << r.strongly_fused_ring_neighbours () << " strongly";
    os << ')';
  }

  if (r.is_aromatic ())
    os << " arom";
  else if (r.is_non_aromatic ())
    os << " non-arom";

  if (RING_FRAGMENT_MEMBERSHIP_UNKNOWN != r.fragment_membership ())
    os << " in fragment " << r.fragment_membership ();

  os << " has " << r.number_elements () << " atoms :";

  for (int i = 0; i < r.number_elements (); i++)
  {
    os << " " << r[i];
  }

  return os;
}

int
path_length_comparitor_longer (Path * const * pp1, Path * const * pp2)
{
  assert (NULL != pp1);
  assert (NULL != pp2);

  const Path * p1 = *pp1;
  const Path * p2 = *pp2;

  const int n1 = p1->number_elements ();
  const int n2 = p2->number_elements ();

  if (n1 > n2)
    return 1;
  else if (n1 == n2)
    return 0;
  else
    return -1;
}

int
path_length_comparitor_shorter (const void * p1, const void * p2)
{
  assert (NULL != p1);
  assert (NULL != p2);

  const Path ** pp1 = (const Path **) p1;
  const Path ** pp2 = (const Path **) p2;

  const int n1 = (*pp1)->number_elements ();
  const int n2 = (*pp2)->number_elements ();

  if (n1 < n2)
    return 1;
  else if (n1 == n2)
    return 0;
  else
    return -1;
}

//#define DEBUG_SET_BONDS_IN_RING

int
Molecule::_set_bonds_in_ring (int * zbonds,
                              const Ring * r)
{
#ifdef DEBUG_SET_BONDS_IN_RING
  cerr << "Finding bonds for ring " << (*r) << endl;
  cerr << "Hits bonds";
#endif

  int n = r->number_elements ();
  atom_number_t prev = r->last_item ();
  for (int i = 0; i < n; i++)
  {
    atom_number_t j = r->item (i);
    int k = _bond_list.which_bond (prev, j);
    assert (k >= 0);
    zbonds[k] = 1;

    prev = j;

#ifdef DEBUG_SET_BONDS_IN_RING
    const Bond * b = _bond_list[k];
    cerr << " (" << b->a1 () << "," << b->a2 () << ") " << k;
#endif
  }

#ifdef DEBUG_SET_BONDS_IN_RING
  cerr << endl;
#endif

  return 1;
}

int
Molecule::_set_bonds_in_ring (IW_Bits_Base * bits,
                              const Ring * r)
{
  int n = r->number_elements ();
  atom_number_t prev = r->last_item ();
  for (int i = 0; i < n; i++)
  {
    atom_number_t j = r->item (i);
    int k = _bond_list.which_bond (prev, j);
    bits->set (k, 1);

    prev = j;

#ifdef DEBUG_SET_BONDS_IN_RING
    const Bond * b = _bond_list[k];
    cerr << " (" << b->a1 () << "," << b->a2 () << ") " << k;
#endif
  }

#ifdef DEBUG_SET_BONDS_IN_RING
  cerr << endl;
#endif

  return 1;
}

int
Ring::contains_bond (atom_number_t a1, atom_number_t a2) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (a1 != _things[i])
      continue;

    if (i < _number_elements - 1 && a2 == _things[i + 1])
      return 1;

    if (i > 0 && a2 == _things[i - 1])
      return 1;

    if (i == _number_elements - 1)
      return a2 == _things[0];

    if (0 == i)
      return a2 == _things[_number_elements - 1];

    return 0;
  }

  return 0;
}

int
Ring::contains_both (atom_number_t a1, atom_number_t a2) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i] == a1)
    {
      for (int j = i + 1; j < _number_elements; j++)
      {
        if (_things[j] == a2)
          return 1;
      }
      return 0;
    }
    if (_things[i] == a2)
    {
      for (int j = i + 1; j < _number_elements; j++)
      {
        if (_things[j] == a1)
          return 1;
      }
      return 0;
    }
  }

  return 0;
}

int
Ring::_compute_bonds_shared_with (const Ring & rhs,
                                  int my_direction,
                                  int lhs_ndx,
                                  int rhs_ndx) const
{
  assert (_things[lhs_ndx] == rhs._things[rhs_ndx]);

  atom_number_t a1 = next_after_wrap(lhs_ndx, my_direction);

  int initial_rhs_ndx = rhs_ndx;

  int rhs_direction = 0;

  if (a1 == rhs.next_after_wrap(rhs_ndx, 1))
    rhs_direction = 1;
  else
  {
    rhs_ndx = initial_rhs_ndx;
    if (a1 == rhs.next_after_wrap(rhs_ndx, -1))
      rhs_direction = -1;
  }

  if (0 == rhs_direction)
    return 0;

  int rc = 1;

  while (rc < _number_elements)   // guard against the same ring as both lhs and rhs
  {
    a1 = next_after_wrap(lhs_ndx, my_direction);
    if (a1 != rhs.next_after_wrap(rhs_ndx, rhs_direction))
      return rc;

    rc++;
  }

  return rc;   // probably should not come to here
}

int
Ring::compute_bonds_shared_with (const Ring & rhs) const
{
  atom_number_t lhs_ndx;
  atom_number_t rhs_ndx = -1;

  for (lhs_ndx = 0; lhs_ndx < _number_elements; lhs_ndx++)
  {
    rhs_ndx = rhs.index (_things[lhs_ndx]);
    if (rhs_ndx < 0)
      continue;

    if (lhs_ndx == _number_elements - 1)   // only one atom in common, so obviously zero bonds in common
      return 0;

    break;
  }

  if (rhs_ndx < 0)  // didn't find anything
    return 0;

// At this stage, we have a match, but we don't know the directions.
// We need to match backwards and forwards since we might have a case like
// 10 9 8 7 6            first matching atom is the central one in a sequence
//  3 4 0 9 10 6

  int rc = _compute_bonds_shared_with(rhs, 1, lhs_ndx, rhs_ndx);
  rc += _compute_bonds_shared_with(rhs, -1, lhs_ndx, rhs_ndx);

  return rc;
}

/*
  Someone has made a determination to transfer a non_sssr_ring to
  the sssr ring set
*/

int
Molecule::_transfer_from_non_sssr_to_sssr_ring_set (int nssr_ndx, int sssr_ndx)
{
  Ring * rfrom = _sssr_rings[sssr_ndx];
  Ring * rto = _non_sssr_rings[nssr_ndx];

  for (int i = 0; i < _sssr_rings.number_elements (); i++)
  {
    if (i == sssr_ndx)
      continue;

    _sssr_rings[i]->ring_moving_from_non_sssr_to_sssr (rfrom, rto);
  }

// transfer any extra information

  rto->set_ring_number (sssr_ndx);
  rto->set_fused_system_identifier (rfrom->fused_system_identifier ());
  rto->set_fragment_membership (rfrom->fragment_membership ());

// Do the transfer

  _sssr_rings[sssr_ndx] = rto;
  _non_sssr_rings[nssr_ndx] = rfrom;

  return 1;
}

int
Ring::ring_moving_from_non_sssr_to_sssr (Ring * rfrom,
                                         Ring * rto)
{
  if (! _fused_neighbours.remove_first (rfrom))
    return 0;

  _fused_neighbours.add (rto);

  return 1;
}

Ring_Atom_Iterator
Ring::find (atom_number_t a) const
{
  Ring_Atom_Iterator rc(*this);

  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i] != a)
      continue;

    rc.set_index(i);
    return rc;
  }

  rc.set_index(_number_elements);

  return rc;
}
