#include <iostream>
#include <string.h>
#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

// Define this symbol before including molecule.h

#define COMPILING_MOLECULER_CC
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "iwqsort.h"
#include "iwbits.h"
#include "primes.h"
#include "misc.h"

class Rings_Found;

#include "path.h"
#include "molecule.h"
#include "misc2.h"

void
Molecule::_initialise_ring_membership()
{
  assert (NULL == _ring_membership);

  if (0 == _number_elements)
    return;

  _ring_membership = new_int(_number_elements, IW_RING_MEMBERSHIP_NOT_COMPUTED);
}

void
Molecule::_determine_ring_or_non_ring(atom_number_t a)
{
  int frag_id = _fragment_information.fragment_membership(a);

  assert (frag_id >= 0);

  if (NULL == _ring_membership)
    _initialise_ring_membership();

  assert (IW_RING_MEMBERSHIP_NOT_COMPUTED == _ring_membership[a]);

  int * tmp = new_int(_number_elements); std::unique_ptr<int[]> free_tmp(tmp);

  _find_raw_rings_for_fragment(frag_id, tmp);

  return;
}

int
Molecule::is_non_ring_atom (atom_number_t a)
{
  assert (ok_atom_number(a));

  if (0 == nrings())     // molecule has no rings, all are non ring atoms
    return 1;

  if (NULL == _ring_membership)
    _initialise_ring_membership();

  if (IW_RING_MEMBERSHIP_NOT_COMPUTED == _ring_membership[a])
    _determine_ring_or_non_ring(a);

  if (0 == _ring_membership[a])
    return 1;
  else
    return 0;
}

int
Molecule::is_ring_atom (atom_number_t a)
{
  assert (ok_atom_number(a));

  if (0 == nrings())    // molecule has no rings, no atom is a ring atom
    return 0;

  if (NULL == _ring_membership)
    _initialise_ring_membership();

  if (IW_RING_MEMBERSHIP_NOT_COMPUTED == _ring_membership[a])
  {
    if (_things[a]->number_elements() <= 1)   // single and unconnected atoms are not in rings
    {
      _ring_membership[a] = 0;
      return 0;
    }

    _determine_ring_or_non_ring(a);
  }

  if (0 == _ring_membership[a])
    return 0;
  else
    return 1;
}

int
Molecule::nrings (atom_number_t a)
{
  assert (ok_atom_number(a));

  if (NULL != _ring_membership && _ring_membership[a] >= 0)
    return _ring_membership[a];

  if (0 == nrings())     // no atoms in molecule, atom A not in any
    return 0;

  if (NULL == _ring_membership)
    _initialise_ring_membership();

  if (1 == _things[a]->ncon())
  {
    _ring_membership[a] = 0;
    return 0;
  }

  int frag_id = _fragment_information.fragment_membership(a);

  if (0 == _fragment_information.rings_in_fragment(frag_id))
  {
    _ring_membership[a] = 0;
    return 0;
  }

  if (IW_RING_MEMBERSHIP_NOT_COMPUTED == _ring_membership[a])
    _determine_ring_or_non_ring(a);

  if (_ring_membership[a] >= 0)
    return _ring_membership[a];

//_determine_sssr_ring_membership(a);

  _determine_sssr_for_fragment(frag_id);

  assert (_ring_membership[a] > 0);

  return _ring_membership[a];
}

int
Molecule::ring_bond_count (atom_number_t zatom)
{
  int nr = nrings(zatom);

  if (0 == nr)
    return 0;

  if (1 == nr)
    return 2;

  const Atom * a = _things[zatom];

  int acon = a->ncon();

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->nrings())
      rc++;
  }

  return rc;
}

int
Molecule::is_part_of_fused_ring_system(atom_number_t a)
{
  int nr = nrings();
  if (0 == nr)
    return 0;

  if (0 == nrings(a))
    return 0;

  (void) ring_membership();   // force sssr

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = _sssr_rings[i];
    if (! r->contains(a))
      continue;

    return r->is_fused();
  }

  return 0;
}

/*
  The number of rings in the system that includes atom A
*/

int
Molecule::fused_system_size_no_compute (atom_number_t a) const
{
  for (int i = 0; i < _sssr_rings.number_elements(); i++)
  {
    const Ring * r = _sssr_rings[i];
    if (! r->contains(a))
      continue;

    if (! r->is_fused())      // an isolated ring
      return 1;

//  See how many rings in this system

    int rc = 0;
    int fsid = r->fused_system_identifier();

    for (int j = 0; j < _sssr_rings.number_elements(); j++)
    {
      const Ring * rj = _sssr_rings[j];
      if (fsid == rj->fused_system_identifier())
        rc++;
    }
    return rc;
  }

  return 0;   // not in any rings
}

int
Molecule::fused_system_size (atom_number_t a)
{
  if (0 == nrings(a))
    return 0;

  (void) ring_membership();     // force sssr    
  return fused_system_size_no_compute( a );
}

int
Molecule::rings_with_fused_system_identifier (int f)
{
  if (0 == nrings())
    return 0;

  (void) ring_membership();     // force sssr

  int rc = 0;

  for (int i = 0; i < _sssr_rings.number_elements(); i++)
  {
    const Ring * r = _sssr_rings[i];
    if (r->fused_system_identifier() == f)
      rc++;
  }

  return rc;
}

int
Molecule::_unused_fused_ring_system_identifier()
{
  int nr = _experimental_sssr_rings.number_elements();
  if (0 == nr)
    return 0;

  int rc = -100;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = _experimental_sssr_rings[i];
    if (! r->is_fused())
      continue;

    if (r->fused_system_identifier() > rc)
      rc = r->fused_system_identifier();
  }

  if (rc < 0)
    rc = 0;

  return rc;
}

int
Molecule::print_ring_info (std::ostream & os) const
{
  assert (ok());

  if (IW_NRINGS_NOT_COMPUTED == _number_sssr_rings && IW_NRINGS_NOT_COMPUTED == _nrings)
  {
    os << "Nrings not computed\n";
    return 1;
  }

  if (IW_NRINGS_NOT_COMPUTED != _number_sssr_rings)
    os << "Molecule contains " << _number_sssr_rings << " SSSR rings\n";

  if (IW_NRINGS_NOT_COMPUTED != _nrings)
    os << "Molecule contains " << _nrings << " computed rings\n";

  int nr = _sssr_rings.number_elements();
  if (0 == nr)
    os << "SSSR rings not perceived\n";
  else
  {
    for (int i = 0; i < nr; i++)
    {
      const Ring * r = _sssr_rings[i];
      os << " Ring " << i << " " << (*r) << endl;
    }

    int nq = _non_sssr_rings.number_elements();
    if (nq)
    {
      os << nq << " non sssr rings\n";
      for (int i = 0; i < nq; i++)
      {
        const Ring * r = _non_sssr_rings[i];
        os << " Non SSSR Ring " << i << ' ' << (*r) << endl;
      }
    }
  }

  if (NULL == _ring_membership)
    os << "Ring membership is NULL\n";
  else
  {
    for (int i = 0; i < _number_elements; i++)
    {
//    const Atom * a = _things[i];

      os << "Atom " << i << ' ' << smarts_equivalent_for_atom(i);

      if (_ring_membership[i] >= 0)
        os << " in " << _ring_membership[i] << " rings\n";
      else if (IW_RING_MEMBERSHIP_NOT_COMPUTED == _ring_membership[i])
        os << " ring membership not computed\n";
      else if (IW_RING_MEMBERSHIP_IS_A_RING_ATOM == _ring_membership[i])
        os << " is a ring atom\n";
    }
  }

  return 1;
}

int
Molecule::experimental_print_ring_info (std::ostream & os) const
{
  if (IW_NRINGS_NOT_COMPUTED == _nrings && IW_NRINGS_NOT_COMPUTED == _number_sssr_rings)
  {
    os << "Nrings not computed\n";
    return 1;
  }

  if (IW_NRINGS_NOT_COMPUTED != _number_sssr_rings)
    os << "Molecule contains " << _number_sssr_rings << " SSSR rings\n";

  if (IW_NRINGS_NOT_COMPUTED != _nrings)
    os << "Molecule contains " << _nrings << " rings\n";

  int nr = _experimental_sssr_rings.number_elements();
  if (0 == nr)
    os << "SSSR rings not perceived\n";
  else
  {
    for (int i = 0; i < nr; i++)
    {
      const Ring * r = _experimental_sssr_rings[i];
      os << " Ring " << i << " " << (*r) << endl;
    }
  }

  if (NULL == _ring_membership)
    os << "Ring membership is NULL\n";
  else
  {
    for (int i = 0; i < _number_elements; i++)
    {
      os << "Atom " << i;
      if (_ring_membership[i] >= 0)
        os << " in " << _ring_membership[i] << " rings\n";
      else if (IW_RING_MEMBERSHIP_NOT_COMPUTED == _ring_membership[i])
        os << " ring membership not computed\n";
      else if (IW_RING_MEMBERSHIP_IS_A_RING_ATOM == _ring_membership[i])
        os << " is a ring atom\n";
    }
  }

  return 1;
}

/*
  Audit those rings already perceived
*/

int
Molecule::check_ring_info() const
{
  int nr = _sssr_rings.number_elements();

  int rc = 1;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = _sssr_rings[i];
    if (! ok_ring(r))
      return 0;

//  ok_path does not check that first and last are bonded.

    if (! are_bonded(r->item(0), r->last_item()))
      return 0;

    if (i != r->ring_number())
    {
      cerr << "Molecule::check_ring_info: ring " << i << " bad ring number " << r->ring_number() << endl;
      rc = 0;
    }
  }

  return rc;
}

int
Molecule::_compute_number_sssr_rings_by_eulers_formula()
{
  int nf = number_fragments();

  _number_sssr_rings = _bond_list.number_elements() - _number_elements + nf;

  return _number_sssr_rings;
}

int
Molecule::nrings_no_compute() const
{
  if (_nrings >= 0)               // already computed
    return _nrings;
  else
    return NOT_COMPUTED;          // not computed
}

int
Molecule::nrings()
{
  if (nrings_no_compute() >= 0)     // already computed
    return _nrings;

  if (_number_elements <= 2)
  {
    _nrings = 0;
    _number_sssr_rings = 0;
    return 0;
  }

  if (IW_NRINGS_NOT_COMPUTED != _number_sssr_rings)
    return _number_sssr_rings;

  _nrings = _compute_number_sssr_rings_by_eulers_formula();

#ifdef DEBUG_NRINGS
  cerr << "Molecule::nrings: molecule contains " << _nrings << " rings\n";
#endif

  if (0 == _nrings)    // chain molecule
    return 0;

  assert (number_fragments() > 0);

  assert (_nrings >= 0 && _nrings <= _bond_list.number_elements());

  return _nrings;
}

int
Molecule::number_sssr_rings()
{
  if (IW_NRINGS_NOT_COMPUTED != _number_sssr_rings)
    return _number_sssr_rings;

  if (_number_elements < 3)
  {
    _number_sssr_rings = 0;
    _nrings = 0;
    return 0;
  }

  _compute_number_sssr_rings_by_eulers_formula();

  _nrings = _number_sssr_rings;

  return _number_sssr_rings;
}

/*
  The number of rings of a given size for atom A
*/

int
Molecule::nrings (atom_number_t a, int ring_size)
{
  assert (ok_atom_number(a));
  assert (REASONABLE_RING_SIZE(ring_size));

  (void) ring_membership();    // force SSSR determination

  int nr = _sssr_rings.number_elements();

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = _sssr_rings[i];
    if (ring_size == r->number_elements() && r->contains(a))
      rc++;
  }

  return rc;
}

int
Molecule::nrings_size (int ring_size)
{
  assert (ok());

  (void) ring_membership();     // force sssr

  int nr = _sssr_rings.number_elements();

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = _sssr_rings[i];
    if (ring_size == r->number_elements())
      rc++;
  }

  return rc;
}

/*
  Return the array of ring memberships.
*/

const int *
Molecule::ring_membership()
{
  if (_sssr_rings.number_elements() >= nrings())   // > to allow for non sssr determinations
    return _ring_membership;

  if (NULL == _ring_membership)
    _initialise_ring_membership();

  if (0 == nrings())
    set_vector(_ring_membership, _number_elements, 0);
  else
    _force_complete_sssr_determination();

  return _ring_membership;
}

int 
Molecule::ring_or_non_ring (int * result)
{
  if (0 == nrings())
  {
    set_vector(result, _number_elements, 0);
    return 1;
  }

  for (int i = 0; i < _number_elements; i++)
  {
    result[i] = is_ring_atom(i);
  }

  return 1;
}
int
Molecule::ring_membership (int * rm)
{
  if (0 == nrings())
  {
    set_vector(rm, _number_elements, 0);
    return 1;
  }

  const int * tmp = ring_membership();     // this forces sssr determination

  copy_vector(rm, tmp, _number_elements);

  return _nrings;
}

int
Molecule::in_same_ring_no_compute( const atom_number_t & a1, const atom_number_t & a2 ) const
{
  if ( _nrings <= 0 )
    return 0;

  for (int i = 0; i < _nrings; i++)
  {
    const Ring * r = _sssr_rings[i];
    if (r->contains(a1) && r->contains(a2))
      return 1;
  }
  return 0;
}

int
Molecule::in_same_ring (atom_number_t a1, atom_number_t a2)
{
  assert (ok_2_atoms(a1, a2));

  if (0 == nrings())
    return 0;

  (void) ring_membership();

  int nr = _sssr_rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
//  const Ring * r = _sssr_rings[i];

    if (_sssr_rings[i]->contains_both(a1, a2))
      return 1;
  }

  return 0;
}

int
Molecule::in_same_rings (atom_number_t a1, atom_number_t a2)
{
  assert (ok_2_atoms(a1, a2));

  if (0 == nrings())
    return 0;

  if (! is_ring_atom(a1) || ! is_ring_atom(a2))
    return 0;

  (void) ring_membership();

  int nr = _sssr_rings.number_elements();

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = _sssr_rings[i];
    if (r->contains(a1) && r->contains(a2))
      rc++;
  }

  return rc;
}

int
Molecule::in_same_ring_system (atom_number_t a1, atom_number_t a2)
{
  assert (ok_2_atoms(a1, a2));

  if (0 == nrings())
    return 0;

  if (! is_ring_atom(a1) || ! is_ring_atom(a2))
    return 0;

  (void) ring_membership();

  int nr = _sssr_rings.number_elements();

  if (1 == nr)
    return 1;    // both atoms are ring atoms, and there is only one ring

  if (2 == nr)
  {
    if (_sssr_rings[0]->is_fused())   // both atoms are ring atoms, all rings are fused
      return 1;
    else                // rings in separate ring systems
      return _sssr_rings[0]->contains_both(a1, a2) || _sssr_rings[1]->contains_both(a1, a2);
  }

//#define DEBUG_IN_SAME_RING_SYSTEM
#ifdef DEBUG_IN_SAME_RING_SYSTEM
  cerr << "Checking " << nr << " SSSR rings for atoms " << a1 << " and " << a2 << endl;
#endif

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = _sssr_rings[i];

    if (ri->contains(a1))
      ;
    else if (ri->contains(a2))
      std::swap(a1, a2);
    else 
      continue;

#ifdef DEBUG_IN_SAME_RING_SYSTEM
    if (ri->contains(a2))
      cerr << "Both atoms in same ring\n";
#endif

    if (ri->contains(a2))   // in same ring, definitely in same ring system
      return 1;

    const int fsid = ri->fused_system_identifier();

#ifdef DEBUG_IN_SAME_RING_SYSTEM
    cerr << "First atom in fsid " << fsid << ' ' << *ri << endl;
#endif

    if (! ri->is_fused())     // ring is not fused and A2 is not in it
      continue;

    for (int j = i + 1; j < nr; j++)
    {
      const Ring * rj = _sssr_rings[j];

#ifdef DEBUG_IN_SAME_RING_SYSTEM
      cerr << "Check ring " << j << " fsid " << rj->fused_system_identifier() << ", atom2? " << rj->contains(a2) << ' ' << *rj << endl;
#endif

      if (fsid != rj->fused_system_identifier())
        continue;

      if (rj->contains(a2))
        return 1;
    }
  }

  return 0;
}

/*
  The I'th ring. Determine it if needed. Be careful
*/

const Ring *
Molecule::ringi_no_compute (int i) const
{
  assert (ok());

  if ( nrings_no_compute() <= 0 )
    return NULL;

  assert (i >= 0);

  return _sssr_rings[i];
}

const Ring *
Molecule::ringi (int i)
{
  assert (ok());

  if (0 == nrings())
    return NULL;

  assert (i >= 0);

  (void) ring_membership();   // ensure SSSR rings determined

  assert (_sssr_rings.ok_index(i));

  return _sssr_rings[i];
}

/*
  Return the I'th ring of size RING_SIZE
*/

const Ring *
Molecule::ringi (int which_ring, int ring_size)
{
  assert (ok());

  if (0 == nrings())
    return NULL;

  assert (which_ring >= 0 && which_ring < nrings());

  (void) ring_membership();

  int nr = _sssr_rings.number_elements();

  int found = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * r = _sssr_rings[i];
    if (r->number_elements() == ring_size)
    {
      if (which_ring == found)
        return r;

      found++;
    }
  }

  cerr << "Only " << found << " rings of size " << ring_size << endl;
  return NULL;
}

const Ring *
Molecule::ring_containing_atom (atom_number_t a)
{
  if (is_non_ring_atom(a))
    return NULL;

  (void) ring_membership();   // ensure ring membership available

  int nr = _sssr_rings.number_elements();
  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = ringi(i);
    if (ri->contains(a))
      return ri;
  }

  assert (NULL == "Ring atom in no rings!!!!!");
  return NULL;
}

int
Molecule::ring_sizes_for_atom (atom_number_t a, List_of_Ring_Sizes & ring_sizes)
{
  assert (ok_atom_number(a));

  ring_sizes.resize_keep_storage(0);

  if (0 == nrings())
    return 0;

  (void) ring_membership();

  int nr = _sssr_rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = _sssr_rings[i];
    if (r->contains(a))
      ring_sizes.add_if_not_already_present(r->number_elements());
  }

  return ring_sizes.number_elements();
}

int
Molecule::ring_sizes_for_all_atoms (resizable_array_p<List_of_Ring_Sizes> & ring_sizes)
{
  assert (ok());

  if (0 == nrings())
    return 0;

  assert (ring_sizes.number_elements() == _number_elements);

  (void) ring_membership();

  int nr = _sssr_rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = _sssr_rings[i];
    int n = r->number_elements();
    for (int j = 0; j < n; j++)
    {
      atom_number_t k = r->item(j);
      List_of_Ring_Sizes * rsk = ring_sizes[k];
      rsk->add_if_not_already_present(n);
    }
  }

  return 1;
}

/*
*/

int
Molecule::ring_sizes_for_non_sssr_rings (atom_number_t a,
                                         List_of_Ring_Sizes & ring_sizes,
                                         int include_duplicates)
{
  assert (ok_atom_number(a));

  if (ring_sizes.number_elements())
    ring_sizes.resize(0);

  (void) ring_membership();    // force ring perception

  const int nq = _non_sssr_rings.number_elements();
  for (int i = 0; i < nq; i++)
  {
    const Ring * r = _non_sssr_rings[i];
    if (r->contains(a))
    {
      if (include_duplicates)
        ring_sizes.add(r->number_elements());
      else
        ring_sizes.add_if_not_already_present(r->number_elements());
    }
  }

  return ring_sizes.number_elements();
}

int
Molecule::non_sssr_rings_no_compute() const
{
  return _non_sssr_rings.number_elements();
}

int
Molecule::non_sssr_rings()
{
  assert (ok());

  (void) ring_membership();  // force ring perception

  return _non_sssr_rings.number_elements();
}

const Ring *
Molecule::non_sssr_ring_no_compute (int i) const
{
  return _non_sssr_rings[i];
}

const Ring *
Molecule::non_sssr_ring (int i)
{
  assert (ok());

  (void) ring_membership();  // force ring perception

  return _non_sssr_rings[i];
}

int
Molecule::ring_membership_including_non_sssr_rings (int * result)
{
  int rc = ring_membership(result);
  if (0 == rc)
    return 0;

  int nr = _non_sssr_rings.number_elements();
  if (0 == nr)
    return rc;

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = _non_sssr_rings[i];
    r->increment_vector(result, 1);
  }

  return rc;     // just what should this function return?
}

int
Molecule::nrings_including_non_sssr_rings (atom_number_t a)
{
  int rc = nrings(a);

  if (0 == rc)
    return rc;

  for (int i = 0; i < _non_sssr_rings.number_elements(); i++)
  {
    if (_non_sssr_rings[i]->contains(a))
      rc++;
  }

  return rc;
}

/*int
Molecule::saturated (atom_number_t a) 
{
  assert (ok_atom_number (a));

  formal_charge_t fc = formal_charge (a);

  const Element * e = elementi (a);
  int nv = e->normal_valence() + fc;    // the number of bonds available

  int acon = ncon (a) + implicit_hydrogens (a);

  if (acon < nv)
    return 0;

  if (acon == e->normal_valence())
    return 1;

  if (acon > e->normal_valence())
  {
    int alt = number_alternate_valences;
    for (int i = 0; i < number_alternate_valences; i++)
    {
      int j = e->alternate_valence (i) + fc;
      if (acon < j)
        return 0;
      if (acon == j)
        return 1;
    }
  }

// If we come out here, connectivity is higher than all known valences.

  return 1;    // assume saturated
}*/

int
Molecule::experimental_sssr()
{
  int nr = nrings();

  if (0 == nr)
    return 1;

  if (_sssr_rings.number_elements() == nr)
    return 1;

  return 1;
}

//#define DEBUG_FIND_SSSR_FOR_THESE_FUSED_RAW_RINGS

int
Molecule::_find_sssr_for_these_fused_raw_rings (int fused_sys_id, 
                                                int * tmp)
{
#ifdef DEBUG_FIND_SSSR_FOR_THESE_FUSED_RAW_RINGS
  cerr << "Finding SSSR for fsid " << fused_sys_id << endl;
#endif

  set_vector(tmp, _number_elements, 0);

  int nraw = _raw_rings.number_elements();

  int rings_being_processed = 0;
  for (int i = 0; i < nraw; i++)
  {
    const Ring * r = _raw_rings[i];
    if (fused_sys_id == r->fused_system_identifier())
      rings_being_processed++;
  }

#ifdef DEBUG_FIND_SSSR_FOR_THESE_FUSED_RAW_RINGS
  cerr << "Processing " << rings_being_processed << " rings with id " << fused_sys_id << endl;
#endif

  if (2 == rings_being_processed)
    return _easy_case_two_rings(fused_sys_id, tmp);

// The more difficult and more general case

  for (int i = nraw - 1; i >= 0; i--)
  {
    const Ring * r = _raw_rings[i];
    if (fused_sys_id != r->fused_system_identifier())
      continue;

    r->set_vector(tmp, 1);
    _raw_rings.remove_item(i);
  }

  return _pearlman_sssr(tmp, 1);
}

int
Molecule::_determine_sssr_for_fragment (int f)
{
  int nraw = _raw_rings.number_elements();

  resizable_array<int> fsid;

  for (int i = 0; i < nraw; i++)
  {
    const Ring * r = _raw_rings[i];

    atom_number_t a = r->item(0);

    if (f == _fragment_information.fragment_membership(a))
      fsid.add_if_not_already_present(r->fused_system_identifier());
  }

  int * tmp = new int[_number_elements]; std::unique_ptr<int[]> free_tmp(tmp);

  for (int i = 0; i < fsid.number_elements(); i++)
  {
    _find_sssr_for_these_fused_raw_rings(fsid[i], tmp);
  }

  return 1;
}

/*
  Determine the sssr ring membership for a given atom

  WRONG. Because of spiro fusions this doesn't work. If someone
  asks for ring membership of an atom that isn't at the spiro
  fusion point, we'll never do the spiro fusion properly.
*/

#ifdef FIX_WHEN_FSID_CHANGES

int
Molecule::_determine_sssr_ring_membership (atom_number_t a)
{
  assert (IW_RING_MEMBERSHIP_IS_A_RING_ATOM == _ring_membership[a]);

// First identify the fsid's of all raw rings that contain the atom

  int nraw = _raw_rings.number_elements();
  assert(nraw);

// Because of spiro fusions, there can be as many as two fused
// system identifiers

  int fsid1 = -100;
  int fsid2 = -100;

  for (int i = 0; i < nraw; i++)
  {
    const Ring * r = _raw_rings[i];
    if (! r->contains(a))
      continue;

    int f = r->fused_system_identifier();

    if (fsid1 < 0)
      fsid1 = f;
    else if (f == fsid1)
      continue;
    else if (fsid2 < 0)   // found spiro fusion. We are done
    {
      fsid2 = f;
      break;
    }
  }

  assert (fsid1 >= 0);

// Process all raw rings with the same fused system identifier

  int * tmp = new int[_number_elements]; std::unique_ptr<int[]> free_tmp(tmp);

  _find_sssr_for_these_fused_raw_rings(fsid1, tmp);

  if (fsid2 >= 0)
    _find_sssr_for_these_fused_raw_rings(fsid2, tmp);

  return 1;
}

#endif

/*int
Molecule::qnrings (atom_number_t a)
{
  assert (ok_atom_number (a));

  if (NULL != _ring_membership && _ring_membership[a] >= 0)
    return _ring_membership[a];

  int nr = nrings();
  if (0 == nr)
    return 0;

  (void) number_fragments();

  int frag_id = _things[a]->fragment_membership();
  if (0 == _bonds_in_fragment[frag_id] - _atoms_in_fragment[frag_id] + 1)
    return 0;

  if (IW_RING_MEMBERSHIP_NOT_COMPUTED == _ring_membership[a])
  {
    _determine_ring_or_non_ring (a);
    if (_ring_membership[a] >= 0)
      return _ring_membership[a];
  }

  assert (IW_RING_MEMBERSHIP_IS_A_RING_ATOM == _ring_membership[a]);
    
  _determine_sssr_ring_membership (a);

  assert (_ring_membership[a] >= 0);

  return _ring_membership[a];
}*/

int
Molecule::_sssr_for_all_raw_rings (int * tmp)
{
  while (_raw_rings.number_elements())
  {
    const Ring * r = _raw_rings[0];
    int fused_sys_id = r->fused_system_identifier();

    _find_sssr_for_these_fused_raw_rings (fused_sys_id, tmp);
  }

  _nrings = _sssr_rings.number_elements();   // in case a non-sssr determination

  return 1;
}

/*
  Many times a procedure will need to ensure that a full sssr has been done
  First determine all the ring/non-ring atoms. This may also
  find the sssr, if the ring system is not too complex.
*/

int
Molecule::_force_complete_sssr_determination()
{
  int nr = nrings();

  if (0 == nr)
    return 1;

  if (NULL == _ring_membership)
  {
    _initialise_ring_membership();
    _find_raw_rings();     // gets all of them
  }
  else
  {
    for (int i = 0; i < _number_elements; i++)
    {
      if (IW_RING_MEMBERSHIP_NOT_COMPUTED == _ring_membership[i])
      {
        _determine_ring_or_non_ring(i);

        assert (IW_RING_MEMBERSHIP_NOT_COMPUTED != _ring_membership[i]);
      }
    }
  }

// If all the rings were non-fused, we are done

  int nraw = _raw_rings.number_elements();

  if (0 == nraw)
  {
    _sort_by_ring_size();

    _assign_ring_numbers(nr);

    return 1;
  }

  int * tmp = new int[_number_elements]; std::unique_ptr<int[]> free_tmp(tmp);

  _sssr_for_all_raw_rings(tmp);

  _sort_by_ring_size();

  _assign_ring_numbers(nr);

  return 1;
}

/*
  Fetch all the atoms in a given ring system into RESULT
  We return the number of rings in the system
*/

int
Molecule::get_fused_system (int fsid, Set_of_Atoms & result)
{
  (void) ring_membership();     // force sssr determination

  int nr = _sssr_rings.number_elements();
  if (nr <= 1)
  {
    cerr << "Molecule::get_fused_system: only " << nr << " rings in the molecule\n";
    iwabort();
    return 0;
  }

  if (0 == result.elements_allocated())
    result.resize(_number_elements);

  int rc = 0;
  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = _sssr_rings[i];
    if (fsid == ri->fused_system_identifier())
    {
      rc++;
      result.add_non_duplicated_elements(*ri);
    }
  }

  if (0 == rc)
    cerr << "Molecule::get_fused_system: no rings with fsid = " << fsid << endl;

  return rc;
}

int
Molecule::label_atoms_by_ring_system_no_compute (int * r) const
{
  set_vector(r, _number_elements, 0);

  int nr = _sssr_rings.number_elements();

  if (0 == nr)
    return 0;

  _sssr_rings[0]->set_vector(r, 1);   // always

  if (1 == nr)
    return 1;

  if (2 == nr)
  {
    if (_sssr_rings[0]->is_fused())
    {
      _sssr_rings[1]->set_vector(r, 1);
      return 1;
    }
    else
    {
      _sssr_rings[1]->set_vector(r, 2);
      return 2;
    }
  }

// We have more than two rings. We don't know if ring[0] is fused or not, so we restart

  int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  int f = 0;       // the number we will assign

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])
      continue;

    const Ring * ri = _sssr_rings[i];

    f++;

    if (i > 0)
      ri->set_vector(r, f);

    if (! ri->is_fused())
      continue;

    for (int j = i + 1; j < nr; j++)
    {
      if (ring_already_done[j])
        continue;

      const Ring * rj = _sssr_rings[j];

      if (ri->fused_system_identifier() != rj->fused_system_identifier())
        continue;

      rj->set_vector(r, f);

      ring_already_done[j] = 1;
    }
  }

  return f;
}

int
Molecule::label_atoms_by_ring_system (int * r)
{
  (void) ring_membership();

  return label_atoms_by_ring_system_no_compute (r);
}

int
Molecule::label_atoms_by_ring_system_including_spiro_fused(int * r)
{
  (void) ring_membership();

  set_vector(r, _number_elements, 0);

  const int nr = _sssr_rings.number_elements();

  if (0 == nr)
    return 0;

  if (1 == nr)
  {
    const Ring * ri = _sssr_rings[0];
    ri->set_vector(r, 1);
    return 1;
  }

  int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  int f = 0;

  int rings_processed = 0;

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])
      continue;

    const Ring * ri = _sssr_rings[i];

    f++;

    ri->set_vector(r, f);

    rings_processed++;

    if (nr == rings_processed)
      break;

    if (ri->is_fused())
    {
      int fused_ring_neighbours = ri->fused_ring_neighbours();
      for (int j = 0; j < fused_ring_neighbours; j++)
      {
        const Ring * rj = ri->fused_neighbour(j);
        rj->set_vector(r, f);
        ring_already_done[rj->ring_number()] = 1;
        rings_processed++;
      }

      if (nr == rings_processed)
        return f;
    }

//  Add any spiro fusions, and maybe other ring systems

    while (1)
    {
      int added_another_ring = 0;

      for (int j = i + 1; j < nr; j++)
      {
        if (ring_already_done[j])
          continue;

        const Ring * rj = _sssr_rings[j];

        if (! rj->any_members_set_in_array(r))
          continue;

        rj->set_vector(r, f);

        ring_already_done[j] = 1;

        rings_processed++;

        if (nr == rings_processed)
          return f;

        added_another_ring = 1;
      }

      if (0 == added_another_ring)
        break;
    }
  }

  return f;   // the number of different ring systems we have identified
}

int
Molecule::in_ring_of_given_size (atom_number_t a, int ring_size)
{
  if (! is_ring_atom(a))
    return 0;

  (void) ring_membership();     // force sssr determination

  int nr = _sssr_rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = _sssr_rings[i];

    if (ri->number_elements() != ring_size)
      continue;

    if (ri->contains(a))
      return 1;
  }

  return 0;      // didn't find a match
}  

class
Ring_Size_Comparator
{
  private:
  public:
    int operator () (const Ring * r1, const Ring * r2) const
    {
      int n1 = r1->number_elements();
      int n2 = r2->number_elements();
      if (n1 < n2)
        return -1;
      if (n1 > n2)
        return 1;

      return 0;
    }
};

static Ring_Size_Comparator rsc;

void
Molecule::_sort_by_ring_size ()
{
  int n = _sssr_rings.number_elements();

  if (n < 2)
    return;

  _sssr_rings.iwqsort(rsc);
}

void
Molecule::_assign_ring_numbers (int nr)
{
  nr = _sssr_rings.number_elements();   // may be different from NR with -K esssr

  for (int i = 0; i < nr; i++)
  {
    _sssr_rings[i]->set_ring_number(i);
  }

  return;
}

/*
  Very detailed examination of R to see whether or not it is consistent with being part of our molecule
*/

int
Molecule::ok_ring (const Ring * r) const
{
  assert (r->ok());

  int ring_size = r->number_elements();

  if (0 == ring_size)    // I guess the empty set is OK
    return 1;

  if (r->item(0) < 0 || r->item(0) >= _number_elements)
    return 0;

  if (1 == _number_elements)      // not very interesting
    return 1;

  for (int i = 1; i < ring_size; i++)
  {
    atom_number_t j = r->item(i - 1);
    atom_number_t k = r->item(i);

    if (k < 0 || k >= _number_elements)
      return 0;
    
    if (! _things[j]->is_bonded_to(k))
      return 0;
  }

  return 1;
}

/*
  We are processing a fused system with just two rings
  The TMP array has already been set with the atoms involved.

  We have just two cases:
    The rings are already SSSR rings.
    One of the rings includes the other

  But we only process the case of one shared bond, so unless all atoms
  are shared, we go to _pearlman_sssr
*/

//#define DEBUG_EASY_CASE_TWO_RINGS

int
Molecule::_easy_case_two_rings(int fused_sys_id,
                               int * tmp)
{
//cerr << "Via _easy_case_two_rings\n";

  int fid = _unused_fused_system_identifier();    // we need to assign a new fused system identifier to these rings

  int nraw = _raw_rings.number_elements();

#ifdef DEBUG_EASY_CASE_TWO_RINGS
  cerr << nraw << " raw rings, looking for fsed " << fused_sys_id << endl;
#endif

  Ring * r1 = NULL;
  Ring * r2 = NULL;
  for (int i = nraw - 1; i >= 0; i--)
  {
    Ring * r = _raw_rings[i];

    if (fused_sys_id != r->fused_system_identifier())
      continue;

    r->increment_vector(tmp, 1);
    r->set_fused_system_identifier(fid);
    r->set_is_fused(1);
    if (NULL == r1)
    {
      r1 = r;
      _raw_rings.remove_no_delete(i);
    }
    else if (NULL == r2)
    {
      r2 = r;
      _raw_rings.remove_no_delete(i);
      break;
    }
  }

  assert (NULL != r2);

// If there are only two atoms in two rings, we have found the SSSR

  if (2 == count_occurrences_of_item_in_array(2, _number_elements, tmp))
  {
    r1->update_ring_membership(_ring_membership, 1);
    r2->update_ring_membership(_ring_membership, 1);

    r1->set_fused_to(r2, 1);
    r2->set_fused_to(r1, 1);

    if (r1->number_elements() < r2->number_elements())
    {
      _add_ring_to_sssr(r1);
      _add_ring_to_sssr(r2);
    }
    else
    {
      _add_ring_to_sssr(r2);
      _add_ring_to_sssr(r1);
    }

#ifdef NEVER_WANT_THESE
    if (accumulate_non_sssr_rings())
    {
      Ring * r = new Ring;
      r->add_non_duplicated_elements(*r1);
      r->add_non_duplicated_elements(*r2);
      r->set_fragment_membership(r1->fragment_membership());
      _non_sssr_rings.add(r);
    }
#endif

    return 2;
  }

// We have more than two atoms shared. Is there only one bond shared? 
// Find an atom in 1 ring that is bonded to an atom in 2 rings

  Ring other_ring;

  other_ring.resize(_number_elements);
  other_ring.set_fused_system_identifier(fid);
  other_ring.set_is_fused(1);
  other_ring.set_fragment_membership(r1->fragment_membership());

  atom_number_t first_2 = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _number_elements; i++)
  {
    if (1 != tmp[i])
      continue;

    const Atom * a = _things[i];
    int acon = a->ncon();
    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);
      if (2 != tmp[k])
        continue;

      other_ring.add(i);
      first_2 = k;
      break;
    }

    if (INVALID_ATOM_NUMBER != first_2)
      break;
  }

#ifdef DEBUG_EASY_CASE_TWO_RINGS
  cerr << "Initial start atom " << other_ring[0] << " first_2 is " << first_2 << endl;
#endif

// Build up the rest of the tmp[] == 1 ring in bonded order

  atom_number_t start_atom = other_ring[0];
  atom_number_t second_2 = INVALID_ATOM_NUMBER;
  while (1)
  {
#ifdef DEBUG_EASY_CASE_TWO_RINGS
    cerr << "start atom is " << start_atom << endl;
#endif

    tmp[start_atom] = 0;                   // don't want to turn back on ourselves
    const Atom * a = _things[start_atom];
    int acon = a->ncon();
    int found_connection_with_tmp1 = 0;    // is the ring continuing
    atom_number_t connected_atom_with_tmp2 = INVALID_ATOM_NUMBER;
    for (int i = 0; i < acon; i++)
    {
      atom_number_t j = a->other(start_atom, i);
      if (0 == tmp[j])
        continue;

      if (1 == tmp[j])
      {
        other_ring.add(j);
        found_connection_with_tmp1 = 1;
        start_atom = j;
        break;
      }

      if (2 == tmp[j])
        connected_atom_with_tmp2 = j;
    }

    if (found_connection_with_tmp1)
      continue;

    second_2 = connected_atom_with_tmp2;
    break;
  }

// We have identified two atoms at the ends of a sequence of atoms in just 1 ring. 
// Are these bonded?

#ifdef DEBUG_EASY_CASE_TWO_RINGS
  cerr << "first_2 is " << first_2 << " second_2 is " << second_2 << endl;
#endif
  if (! _things[first_2]->is_bonded_to(second_2))
  {
    r1->set_vector(tmp, 1);
    r2->set_vector(tmp, 1);
    delete r1;
    delete r2;
//  cerr << "Must be done via SSSR\n";
    return _pearlman_sssr(tmp, 1);
  }

// The atoms in the Ring must be in bonded order...

  other_ring.add(second_2);
  other_ring.add(first_2);

#ifdef DEBUG_EASY_CASE_TWO_RINGS
  cerr << "R1    " << (*r1) << endl;
  cerr << "R2    " << (*r2) << endl;
  cerr << "Other " << other_ring << endl;
#endif

// The two smallest rings are the SSSR rings.

  int nr1 = r1->number_elements();
  int nr2 = r2->number_elements();
  int notr = other_ring.number_elements();

  if (nr1 <= notr && nr2 < notr)    // R1 and R2 are the smallest
    ;
  else if (notr <= nr1 && nr2 < nr1)
    *r1 = other_ring;
  else if (notr <= nr2 && nr1 < nr2)
    *r2 = other_ring;

#ifdef DEBUG_EASY_CASE_TWO_RINGS
    cerr << "Smallest R1 " << (*r1) << endl;
    cerr << "Smallest R2 " << (*r2) << endl;
#endif

  r1->update_ring_membership(_ring_membership, 1);
  r2->update_ring_membership(_ring_membership, 1);

  r1->set_fused_to(r2, 1);
  r2->set_fused_to(r1, 1);

  if (r1->number_elements() < r2->number_elements())
  {
    _add_ring_to_sssr(r1);
    _add_ring_to_sssr(r2);
  }
  else
  {
    _add_ring_to_sssr(r2);
    _add_ring_to_sssr(r1);
  }

#ifdef NEVER_WANT_THESE
  not sure what I was thinking of here, we never want that outside ring

  if (accumulate_non_sssr_rings())
  {
    Ring * r = new Ring;
    r->add_non_duplicated_elements(*r1);
    r->add_non_duplicated_elements(*r2);
    r->set_fragment_membership(r1->fragment_membership());
    _non_sssr_rings.add(r);
  }
#endif

  return 2;
}

int
Molecule::_just_one_unclassified_spinach_connection (atom_number_t zatom,
                                                     const int * spinach) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  int unclassified_connections_found = 0;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (1 == spinach[j])
      continue;

    unclassified_connections_found++;
  }

//cerr << "From atom " << zatom << " there are " << unclassified_connections_found << " unclassified connections\n";

  return 1 == unclassified_connections_found;
}

/*
  The preset values in SPINACH must be something other than 1
*/

int
Molecule::identify_spinach_preset (int * spinach) const
{
  Set_of_Atoms stack;
  stack.resize(_number_elements);

  for (int i = 0; i < _number_elements; i++)
  {
    if (1 == spinach[i])
    {
      cerr << "Molecule::identify_spinach_preset:cannot use preset value of 1\n";
      return 0;
    }

    if (spinach[i])
      continue;

    if (1 == _things[i]->ncon())
      stack.add(i);
  }

  if (0 == stack.number_elements())
    return 0;    // maybe we should return something different

  return _identify_spinach(spinach, stack);
}

int
Molecule::identify_spinach (int * spinach)
{
  if (0 == nrings())
  {
    set_vector(spinach, _number_elements, 1);
    return _number_elements;
  }

  if (0 == _number_elements)
    return 0;

  set_vector(spinach, _number_elements, 0);

  Set_of_Atoms stack;
  stack.resize(_number_elements);

  for (int i = 0; i < _number_elements; i++)
  {
    if (1 == _things[i]->ncon())
    {
      stack.add(i);
      spinach[i] = 1;
    }
  }

  return _identify_spinach(spinach, stack);
}

int
Molecule::_identify_spinach (int * spinach,
                             Set_of_Atoms & stack) const
{
  int rc = 0;

  while (stack.number_elements() > 0)
  {
    atom_number_t zatom = stack.pop();

//  cerr << "Processing atom " << zatom << " spinach? " << spinach[zatom] << endl;

    const Atom * a = _things[zatom];

    int acon = a->ncon();

    if (1 == acon)
      ;
    else if (spinach[zatom])
      continue;

    spinach[zatom] = 1;

    rc++;

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(zatom, j);

      if (spinach[k] > 0)
        continue;

      if (2 == _things[k]->ncon())   // will definitely be part of spinach
        ;
      else if (! _just_one_unclassified_spinach_connection(k, spinach))
        continue;

      stack.add(k);
    }
//  cerr << "After atom " << zatom << " stack contains " << stack.number_elements() << endl;
  }

  return rc;
}

int
Molecule::rings_with_strongly_fused_ring_neighbours()
{
  int nr = nrings();

  if (nr < 3)
    return 0;

  (void) ring_membership();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = _sssr_rings[i];

    if (ri->strongly_fused_ring_neighbours())
      rc++;
  }

  return rc;
}

static Int_Comparator_Larger icl;

#if defined(__GNUG__)
template void iwqsort<int, Int_Comparator_Larger>(int*, int, Int_Comparator_Larger&);
template void iwqsort<int, Int_Comparator_Larger>(int*, int, Int_Comparator_Larger&, void*);
template void move_in_from_left<int, Int_Comparator_Larger>(int*, int&, int&, int, Int_Comparator_Larger&, void*);
template void move_in_from_right<int, Int_Comparator_Larger>(int*, int&, int&, Int_Comparator_Larger&);
template void compare_two_items<int, Int_Comparator_Larger>(int*, Int_Comparator_Larger&, void*);
template void swap_elements<int>(int&, int&, void*);
#endif

int
Molecule::is_spiro_fused (atom_number_t zatom)
{
  if (nrings() < 2)
    return 0;

  if (is_non_ring_atom(zatom))
    return 0;

  const Atom * a = _things[zatom];

  if (4 != a->ncon())
    return 0;

  (void) ring_membership();   // force sssr

  Set_of_Atoms connections;
  connections.resize(4);

  for (int i = 0; i < 4; i++)
  {
    const Bond * b = a->item(i);
    if (0 == b->nrings())
      return 0;

    connections.add(b->other(zatom));
  }

// Ugly, we get things like this
// C12(N3CC4(C)C(=O)C(C)(C3)CN1C4)C1=C(C=CC(=C1)Br)NC2=O 417829 PBCHM1097623:PBCHM6966522
// where there are 3 or even 4 rings incident on the spiro fusion.
// For now, just gather all rings incident on this atom

  resizable_array<const Ring *> rings_here;

  for (int i = 0; i < _sssr_rings.number_elements(); i++)
  {
    const Ring * ri = _sssr_rings[i];

    if (ri->contains(zatom))
      rings_here.add(ri);
  }

  int nr = rings_here.number_elements();

#ifdef DEBUG_IS_SPIRO_FUSED
  cerr << "Got " << nr << " rings\n";
#endif

  if (nr < 2)
  {
    cerr << "Molecule::is_non_ring_atom:strange, atom " << zatom << " not in two rings\n";
    return 0;
  }

//atom_number_t a0 = connections[0];
//atom_number_t a1 = connections[1];
//atom_number_t a2 = connections[2];
//atom_number_t a3 = connections[3];

#ifdef DEBUG_IS_SPIRO_FUSED
  cerr << "Centre atom " << zatom << endl;
  cerr << "atom " << a0 << " ring_membership " << _ring_membership[a0] << endl;
  cerr << "atom " << a1 << " ring_membership " << _ring_membership[a1] << endl;
  cerr << "atom " << a2 << " ring_membership " << _ring_membership[a2] << endl;
  cerr << "atom " << a3 << " ring_membership " << _ring_membership[a3] << endl;
#endif

// Form a number for each atom that indicates which of the rings it is in

  int rc[4] = {1, 1, 1, 1};

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = rings_here[i];

    for (int j = 0; j < 4; j++)
    {
      atom_number_t aj = connections[j];

      if (ri->contains(aj))
        rc[j] *= primes[i];
    }
  }

// In a spiro system, the rc array should contain just two different values
// Do a sort. Maybe should do this in place...

  ::iwqsort(&(rc[0]), 4, icl);

  if (rc[0] == rc[1] && rc[1] != rc[2] && rc[2] == rc[3])
    return 1;

  return 0;

// DO the easy case first
// count the number of times the atoms occur in each of the rings. 

#ifdef QOWERU
  if (2 == nrz)
  {
    int sro0 = spiro_ring_occupancy (r1, r2, a0);
    if (0 == sro0 || 11 == sro0)
      return 0;

    int sro1 = spiro_ring_occupancy (r1, r2, a1);
    if (0 == sro1 || 11 == sro1)
      return 0;

    int sro2 = spiro_ring_occupancy (r1, r2, a2);
    if (0 == sro2 || 11 == sro2)
      return 0;

    int sro3 = spiro_ring_occupancy (r1, r2, a3);
    if (0 == sro3 || 11 == sro3)
      return 0;

// Within the sro numbers, we must have two 1's and two 10's

    cerr << sro0 << sro1 << sro2 << sro3 << endl;

    return (22 == sro0 + sro1 + sro2 + sro3);
  }
#endif

// Now the more difficult case of 3 or 4 rings incident on the fusion point

  return 0;
}
int
Molecule::fused_system_identifier (atom_number_t a)
{
  const int nr = nrings();
  if (0 == nr)
    return -1;

  if (! is_ring_atom(a))
    return -1;

  ring_membership();

  const Ring * r = ring_containing_atom(a);

  if (NULL == r)    // should not happen
    return -1;

  return r->fused_system_identifier();
}

const Ring * const *
Molecule::cbeginRing()
{
  (void) ring_membership();

  return _sssr_rings.cbegin();
}

const Ring * const * 
Molecule::cendRing()
{
  (void) ring_membership();

  return _sssr_rings.cend();
}
