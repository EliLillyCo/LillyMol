#include <stdlib.h>
#include <assert.h>
#include <iostream>

/*
  Class for substructure matching. There will be one of these atoms for each
  one in the target molecule. Attributes are computed as they are required.
  All values are assigned TARGET_ATOM_NOT_COMPUTED initially

  Each atom has a pointer 
*/

#define TARGET_ATOM_NOT_COMPUTED -1923

/*
  The value of _nrings can be
   TARGET_ATOM_NOT_COMPUTED
   TARGET_IS_RING_ATOM
   0,1,2,3.....   
*/

#define TARGET_IS_RING_ATOM     -1902

#include "misc.h"

#include "misc2.h"
#include "path.h"
#include "target.h"
#include "substructure.h"

#ifdef FINGERPRINT_SUBSTRUCTURE_SEARCHES
#include "iwmfingerprint.h"
#endif

//#define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

/*
  Apr 99. When a bond is aromatic we can make it no longer match it's
  underlying Kekule form
  Jan 2000 Change this to being default behaviour for better match with daylight
  May 2000, lots of problems, change back

  Jun 2008. Introduce a new value for this parameter, 2. In that case, we
  allow aromatic rings that don't have alternating to retain their
  Kekule forms, but not rings that have alternating forms.
*/

static int _aromatic_bonds_lose_kekule_identity = 0;

void
set_aromatic_bonds_lose_kekule_identity (int s)
{
  _aromatic_bonds_lose_kekule_identity = s;
}

int
aromatic_bonds_lose_kekule_identity ()
{
  return _aromatic_bonds_lose_kekule_identity;
}

/*
  It is optional as to whether we initialise the element counts functionality during a search
*/

static int initialise_element_counts = 0;

void
set_initialise_element_counts (int s)
{
  initialise_element_counts = s;
}

static int global_setting_nrings_includes_non_sssr_rings = 0;

void
set_global_setting_nrings_includes_non_sssr_rings(int s)
{
  global_setting_nrings_includes_non_sssr_rings = s;
}

static int global_setting_nbonds_includes_implicit_hydrogens = 0;

void set_global_setting_nbonds_includes_implicit_hydrogens(int s)
{
  global_setting_nbonds_includes_implicit_hydrogens = s;
}



Bond_and_Target_Atom & 
Target_Atom::other (int i)
{
  if (_ncon && NULL == _other)
    _allocate_other();

  assert (i >= 0 && i < _ncon);

  return _other[i];
}

/*
  When time comes, my guess is that efficiency will improve by making
  _ring_sizes a resizable_array<int> rather than allocating it as is
  done now.
*/

void
Target_Atom::_default_values ()
{
  _other = NULL;

  _nbonds                      = TARGET_ATOM_NOT_COMPUTED;
  _formal_charge               = TARGET_ATOM_NOT_COMPUTED;
  _nrings                      = TARGET_ATOM_NOT_COMPUTED;
  _ring_bond_count             = TARGET_ATOM_NOT_COMPUTED;
  _ncon2                       = TARGET_ATOM_NOT_COMPUTED;
  _hcount                      = TARGET_ATOM_NOT_COMPUTED;
  _aromaticity                 = TARGET_ATOM_NOT_COMPUTED;
  _multiple_bond_to_heteroatom = TARGET_ATOM_NOT_COMPUTED;

  _vinyl                       = TARGET_ATOM_NOT_COMPUTED;
  _aryl                        = TARGET_ATOM_NOT_COMPUTED;
  _lone_pair_count = TARGET_ATOM_NOT_COMPUTED;

  _chirality_fetched = 0;
  _chiral_centre = NULL;

  _isotope                     = TARGET_ATOM_NOT_COMPUTED;
  _userAtomType                = TARGET_ATOM_NOT_COMPUTED;
  _attached_heteroatom_count   = TARGET_ATOM_NOT_COMPUTED;
  _fused_system_size           = TARGET_ATOM_NOT_COMPUTED;

  _all_rings_kekule            = TARGET_ATOM_NOT_COMPUTED;

  _heteroatoms_in_ring = TARGET_ATOM_NOT_COMPUTED;

  return;
};

Target_Atom::Target_Atom ()
{
  _default_values();

  return;
}

/*
  Make an atom never match anything
*/

#define TARGET_ATOM_NOTHING_CAN_MATCH -55

void
Target_Atom::invalidate ()
{
  _ncon2                       = TARGET_ATOM_NOTHING_CAN_MATCH;
  _nbonds                      = TARGET_ATOM_NOTHING_CAN_MATCH;
  _nrings                      = TARGET_ATOM_NOTHING_CAN_MATCH;
  _ring_bond_count             = TARGET_ATOM_NOTHING_CAN_MATCH;
  _formal_charge               = TARGET_ATOM_NOTHING_CAN_MATCH;
  _hcount                      = TARGET_ATOM_NOTHING_CAN_MATCH;
  _aromaticity                 = TARGET_ATOM_NOTHING_CAN_MATCH;
  _attached_heteroatom_count   = TARGET_ATOM_NOTHING_CAN_MATCH;
  _fused_system_size           = TARGET_ATOM_NOTHING_CAN_MATCH;
  _multiple_bond_to_heteroatom = TARGET_ATOM_NOTHING_CAN_MATCH;

  _vinyl                       = TARGET_ATOM_NOTHING_CAN_MATCH;
  _aryl                        = TARGET_ATOM_NOTHING_CAN_MATCH;
  _isotope                     = TARGET_ATOM_NOTHING_CAN_MATCH;
  _userAtomType                = TARGET_ATOM_NOTHING_CAN_MATCH;

  _lone_pair_count = TARGET_ATOM_NOTHING_CAN_MATCH;

  _heteroatoms_in_ring = TARGET_ATOM_NOTHING_CAN_MATCH;

  _chirality_fetched = 0;

  _chiral_centre = NULL;

  return;
}

void
Target_Atom::initialise (Molecule * m, atom_number_t man, Atom * a,
                         Molecule_to_Match * target)
{
  _m = m;
  _my_atom_number = man;
  _my_atom = a;
  _target = target;

  _element = a->element();

  _ncon = a->ncon();

  return;
}

Target_Atom::~Target_Atom ()
{
  if (-8 == _ncon2)
    cerr << "Freeing already deleted Target_Atom\n";
  _ncon2 = -8;

#ifdef USE_IWMALLOC
  iwmalloc_check_all_malloced(stderr);
#endif

  if (NULL != _other)
    delete [] _other;

  return;
}

int
Target_Atom::ok () const
{
  return 1;
//return (OK_ATOM_NUMBER (_m, _my_atom_number));     takes too much time
}

int
Target_Atom::debug_print (std::ostream & os) const
{
  os << "Target atom for atom " << _my_atom_number << " atomic number " << _element->atomic_number() << endl;

  if (TARGET_ATOM_NOT_COMPUTED != _ncon)
  {
    os << "ncon = " << _ncon << endl;

    if (NULL != _other)
    {
      for (int i = 0; i < _ncon; i++)
      {
        Bond_and_Target_Atom & a = _other[i];
        os << "Bond " << i << " of type " << a.bond()->btype() << " to " << a.other()->atom_number() << endl;
      }
    }
  }

  return 1;
}

void
Target_Atom::discard_chirality()
{
  if (NULL != _chiral_centre)
    _chiral_centre = NULL;

  _chirality_fetched = 1;   // don't fetch it again

  return;
}

int
Target_Atom::nbonds ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _nbonds)
  {
    _nbonds = _my_atom->nbonds();
    if (global_setting_nbonds_includes_implicit_hydrogens)
      _nbonds += _my_atom->implicit_hydrogens();
  }

  return _nbonds;
}

int
Target_Atom::hcount ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _hcount)
  {
    if (_ncon && NULL == _other)
      _allocate_other();

    _hcount = _my_atom->implicit_hydrogens();

    for (int i = 0; i < _ncon; i++)
    {
      const Target_Atom * a = _other[i].other();
      if (1 == a->atomic_number())
        _hcount++;
    }
  }

  return _hcount;
}

int
Target_Atom::daylight_x ()
{
  return _my_atom->implicit_hydrogens() + _ncon;
}

int
Target_Atom::userAtomType ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _userAtomType)
    _userAtomType = _my_atom->userAtomType();

  return _userAtomType;
}

int
Target_Atom::isotope ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _isotope)
    _isotope = _my_atom->isotope();

  return _isotope;
}

int
Target_Atom::nrings ()
{
  if (_nrings >= 0)
    return _nrings;

  if (global_setting_nrings_includes_non_sssr_rings)
    _nrings = _m->nrings_including_non_sssr_rings(_my_atom_number);
  else
    _nrings = _m->nrings(_my_atom_number);

//cerr << "Target_Atom::nrings:atom " << _my_atom_number << " nrings " << _nrings << endl;

  return _nrings;
}

int
Target_Atom::ring_bond_count ()
{
//cerr << "On entry _ring_bond_count " << _ring_bond_count << endl;
  if (_ring_bond_count >= 0)
    return _ring_bond_count;

  nrings();

  if (0 == _nrings)
  {
    _ring_bond_count = 0;
    return 0;
  }

  if (_nrings < 2)
  {
    _ring_bond_count = _nrings + 1;
    return _ring_bond_count;
  }

  _ring_bond_count = _m->ring_bond_count(_my_atom_number);

  return _ring_bond_count;
}

/*
  Because of similarities in the code, attached_heteroatom_count also
  computes _multiple_bond_to_heteroatom
*/

int
Target_Atom::attached_heteroatom_count ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _attached_heteroatom_count)
  {
    if (0 == ncon())
      return 0;

    if (NULL == _other)
      _allocate_other();

    _attached_heteroatom_count = 0;
    _multiple_bond_to_heteroatom = 0;
    for (int i = 0; i < _ncon; i++)
    {
      Bond_and_Target_Atom & q = _other[i];

      const Target_Atom * a = q.other();
      atomic_number_t za = a->atomic_number();
      if (6 == za || 1 == za)
        continue;

      _attached_heteroatom_count++;

      const Bond * b = q.bond();
      if (! b->is_single_bond())
        _multiple_bond_to_heteroatom++;
    }
  }

  return _attached_heteroatom_count;
}

/*
  This relies on attached_heteroatom_count () to do all the work
*/

int
Target_Atom::multiple_bond_to_heteroatom ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _multiple_bond_to_heteroatom)
    (void) attached_heteroatom_count();

  return _multiple_bond_to_heteroatom;
}

int
Target_Atom::is_ring_atom ()
{
  if (_nrings >= 0)
    return _nrings;

  if (TARGET_IS_RING_ATOM == _nrings)
    return 1;

  if (_ncon <= 1)
  {
    _nrings = 0;
    return 0;
  }

  if (TARGET_ATOM_NOT_COMPUTED == _nrings)
  {
    if (_m->is_ring_atom(_my_atom_number))
    {
      _nrings = TARGET_IS_RING_ATOM;
      return 1;
    }

    _nrings = 0;
    return 0;
  }

  cerr << "Target_Atom::is_ring_atom: should not come to here " << _nrings << endl;
  iwabort();
  return _nrings;
}

int
Target_Atom::is_non_ring_atom ()
{
  return ! is_ring_atom();
}

int
Target_Atom::all_rings_kekule()
{
//cerr << "Target_Atom::all_rings_containing_atom_are_kekule" << _all_rings_kekule << endl;

  if (TARGET_ATOM_NOT_COMPUTED == _all_rings_kekule)
    _all_rings_kekule = _m->all_rings_containing_atom_are_kekule(_my_atom_number);

  return _all_rings_kekule;
}

int
Target_Atom::ncon2 ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _ncon2)
    _compute_ncon2();

  return _ncon2;
}

int
Target_Atom::ncon2_value_set() const
{
  return TARGET_ATOM_NOT_COMPUTED != _ncon2;
}

/*int
Target_Atom::match_nbonds (int nbonds)
{
  assert (nbonds >= 0);

  if (TARGET_ATOM_NOT_COMPUTED == _nbonds)
    _nbonds = _my_atom->nbonds();

  return _nbonds == nbonds;
}*/

int
Target_Atom::formal_charge ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _formal_charge)
    _formal_charge = _my_atom->formal_charge();

  return _formal_charge;
}

int
Target_Atom::_compute_ncon2 ()
{
  _ncon2 = ncon();   // also forces computation of _ncon

  for (int i = 0; i < _ncon; i++)
  {
    Target_Atom * a = other(i).other();   // Bond_and_Target_Atom->other()
    _ncon2 += a->ncon() - 1;
  }

  return _ncon2;
}

const List_of_Ring_Sizes *
Target_Atom::sssr_ring_sizes ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _nrings)
    _nrings = _m->nrings(_my_atom_number);

  if (0 == _nrings)
  {
    assert (0 == _sssr_ring_sizes.number_elements());
    return &_sssr_ring_sizes;
  }

  if (TARGET_IS_RING_ATOM == _nrings)
    _nrings = _m->nrings(_my_atom_number);

// We are dealing with an atom in a ring

  if (0 == _sssr_ring_sizes.number_elements())
  {
    (void) _m->ring_sizes_for_atom(_my_atom_number, _sssr_ring_sizes);
    if (0 == _sssr_ring_sizes.number_elements())
    {
      cerr << "Atom " << _my_atom_number << " in " << _nrings << " rings, " << _sssr_ring_sizes.number_elements() << endl;
      cerr << "Molecule says " << _m->nrings(_my_atom_number) << endl;
      _m->debug_print(cerr);
      cerr << _m->smiles() << endl;
      assert (_sssr_ring_sizes.number_elements());
    }
  }

  return &_sssr_ring_sizes;
}

void
Target_Atom::_get_ring_sizes_and_aromaticity ()
{
  (void) _m->compute_aromaticity_if_needed();
  int rings_in_target = _m->nrings();

  for (int i = 0; i < rings_in_target; i++)
  {
    const Ring * r = _m->ringi(i);
    if (! r->contains(_my_atom_number))
      continue;

    if (r->is_aromatic())
      _aromatic_ring_sizes.add(r->number_elements());
    else if (r->is_non_aromatic())
      _aliphatic_ring_sizes.add(r->number_elements());
    else
      cerr << "Target_Atom::_get_ring_sizes_and_aromaticity: undetermined ring aromaticity\n";
  }

  return;
}

const List_of_Ring_Sizes *
Target_Atom::aromatic_ring_sizes ()
{
  if (0 == _nrings)
  {
    assert (0 == _sssr_ring_sizes.number_elements());
    return &_sssr_ring_sizes;
  }

// At this stage, we know it has rings.

  if (0 == _aromatic_ring_sizes.number_elements() &&
      0 == _aliphatic_ring_sizes.number_elements())
    _get_ring_sizes_and_aromaticity();

  return &_aromatic_ring_sizes;
}

const List_of_Ring_Sizes *
Target_Atom::aliphatic_ring_sizes ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _nrings)
    _nrings = _m->nrings(_my_atom_number);

  if (0 == _nrings)
  {
    assert (0 == _sssr_ring_sizes.number_elements());
    return &_sssr_ring_sizes;
  }

// At this stage, we know it has rings.

  if (0 == _aromatic_ring_sizes.number_elements() &&
      0 == _aliphatic_ring_sizes.number_elements())
    _get_ring_sizes_and_aromaticity();

  return &_aliphatic_ring_sizes;
}

//#define DEBUG_RING_SIZES_INCLUDING_NON_SSSR

const List_of_Ring_Sizes *
Target_Atom::ring_sizes_including_non_sssr ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _nrings)
    _nrings = _m->nrings(_my_atom_number);

  if (0 == _nrings)
  {
    assert (0 == _sssr_ring_sizes.number_elements());
    return &_sssr_ring_sizes;
  }

// At this stage, we know it has rings. Are there non sssr rings
// associated with this atom too. 

  if (0 == _ring_sizes_including_non_sssr.number_elements())
  {
    (void) sssr_ring_sizes();    // fill the _sssr_ring_sizes array
    
    _ring_sizes_including_non_sssr += _sssr_ring_sizes;

    List_of_Ring_Sizes tmp;
    _m->ring_sizes_for_non_sssr_rings(_my_atom_number, tmp, 0);   // the 0 means no duplicates

    _ring_sizes_including_non_sssr += tmp;
    _ring_sizes_including_non_sssr.sort(int_comparitor_larger);

#ifdef DEBUG_RING_SIZES_INCLUDING_NON_SSSR
    cerr << "Atom " << _my_atom_number << " has ring sizes";
    for (int i = 0; i < _ring_sizes_including_non_sssr.number_elements(); i++)
    {
      cerr << ' ' << _ring_sizes_including_non_sssr[i];
    }
    cerr << endl;
#endif

  }

  return & _ring_sizes_including_non_sssr;
}

/*
  Mar 2000.
  Whenever aromaticity could not be determined, I used to put the special value
  AROMATICITY_NOT_DETERMINED. This seems too difficult, as this special case
  then needs to be checked in substructure_spec.cc
  Instead, change the logic so anything that isn't aromatic is aliphatic;
*/

aromaticity_type_t
Target_Atom::aromaticity ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _aromaticity)
  {
    if (! _m->aromaticity(_my_atom_number, _aromaticity))
    {
//    cerr << "Cannot determine aromaticity of atom " << _my_atom_number << 
//            " (" << _my_atom->atomic_symbol() << ' ' << ncon() << " connections)\n";
      SET_ALIPHATIC_ATOM(_aromaticity);
    }
  }

  return _aromaticity;
}

int
Target_Atom::is_aromatic ()
{
  return IS_AROMATIC_ATOM (aromaticity());
}

Chiral_Centre *
Target_Atom::chiral_centre ()
{
  if (! _chirality_fetched)
  {
    _chiral_centre = _m->chiral_centre_at_atom(_my_atom_number);
    _chirality_fetched = 1;
  }

  return _chiral_centre;
}

/*
  Ran into problems trying to detect vinyl-nitro groups.
  We define a vinyl group as being a SINGLE bond to an unsaturated atom
*/

int
Target_Atom::vinyl ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _vinyl)
  {
    if (0 == ncon())
      return 0;

    if (NULL == _other)
      _allocate_other();

    _vinyl = 0;

    for (int i = 0; i < _ncon; i++)
    {
      const Bond * b = _other[i].bond();
      if (! b->is_single_bond())
        continue;

      Target_Atom * o = _other[i].other();

      if (o->nbonds() == o->ncon())     // fully saturated
        continue;

      if (IS_AROMATIC_ATOM(o->aromaticity()))
        continue;

      _vinyl++;
    }
  }

  return _vinyl;
}

int
Target_Atom::aryl ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _aryl)
  {
    if (0 == ncon())
      return 0;

    if (NULL == _other)
      _allocate_other();

    _aryl = 0;

    for (int i = 0; i < _ncon; i++)
    {
      Target_Atom * o = _other[i].other();

      if (IS_AROMATIC_ATOM(o->aromaticity()))
        _aryl++;
    }
  }

  return _aryl;
}

int
Target_Atom::_count_heteroatoms_in_ring (const Ring * r) const
{
  int rs = r->number_elements();

  int rc = 0;

  for (int i = 0; i < rs; i++)
  {
    atom_number_t j = r->item(i);

#ifdef DEBUG_CARBOCYCLE
    cerr << "Atom number for atom " << j << " is " << _all_atoms[j]->atomic_number() << endl;
#endif

    if (6 != _target->operator[](j).atomic_number())
      rc++;
  }

#ifdef DEBUG_CARBOCYCLE
  cerr << "Found " << rc << " heteroatoms in ring\n";
#endif

  return rc;
}

/*
  Design decision about what does "heteroatoms_in_ring" mean?
  For now, I am defining it as the largest number of heteroatoms
  in a ring that contains our atom.

  Basically this is a bad idea, because of SSSR considerations
  and the common case of atoms being in more than one ring.
  Should have thought of this before...

  Jul 2015. Change behaviour.
  First accumulate all rings that contain the atom.

  If only 1 ring, we are done, the answer is clear.

  If multiple rings, we return the largest number of heteroatoms
  containing in a ring containing the atom, but we stop looking
  once the rings get too large. Massive kludge, should probably
  deprecate this feature, but it is useful...
*/

int
Target_Atom::heteroatoms_in_ring ()
{
#ifdef DEBUG_CARBOCYCLE
  cerr << "Testing carbocycle attribute for atom " << _my_atom_number << endl;
#endif

  if (TARGET_ATOM_NOT_COMPUTED != _heteroatoms_in_ring)     // already determined
    return _heteroatoms_in_ring;

  _heteroatoms_in_ring = 0;

  if (! is_ring_atom())      // a chain atom does not match a ring heteroatom specification
    return 0;

  int nr = _m->nrings();

  resizable_array<const Ring *> rings_this_atom;

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = _m->ringi(i);
    if (! r->contains(_my_atom_number))
      continue;

    rings_this_atom.add(r);
  }

  nr = rings_this_atom.number_elements();

  _heteroatoms_in_ring = _count_heteroatoms_in_ring(rings_this_atom[0]);

  if (1 == nr)    // the easy case
    return _heteroatoms_in_ring;

  for (int i = 1; i < nr; i++)
  {
    const Ring * r = rings_this_atom[i];

#ifdef DEBUG_CARBOCYCLE
    cerr << "Ring " << i << " contains atom " << _my_atom_number << endl;
#endif
    if (r->number_elements() >= 2 * rings_this_atom[0]->number_elements())
      break;

    int h = _count_heteroatoms_in_ring(r);

    if (h > _heteroatoms_in_ring)
      _heteroatoms_in_ring = h;
  }

  return _heteroatoms_in_ring;
}

int
Target_Atom::lone_pair_count ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _lone_pair_count)
  {
    if (! _my_atom->lone_pair_count(_lone_pair_count))
      _lone_pair_count = 0;
  }

  return _lone_pair_count;
}

int
Target_Atom::_allocate_other ()
{
  if (0 == _ncon)
    return 0;

  _other = new Bond_and_Target_Atom[_ncon];

  for (int i = 0; i < _ncon; i++)
  {
    const Bond * b = _my_atom->item(i);

    atom_number_t o = b->other(_my_atom_number);

    Target_Atom & other = _target->operator[](o);

    _other[i].initialise(b, other);
  }

  return _ncon;
}

/*
  A Substructure_Atom is checking a ring closure.
  We check to see whether or not we have a bond to atom A2.
*/

Bond_and_Target_Atom *
Target_Atom::bond_to_atom (atom_number_t a2)
{
  if (_ncon && NULL == _other)
    _allocate_other();

  for (int i = 0; i < _ncon; i++)
  {
    Bond_and_Target_Atom & b = _other[i];

    const Target_Atom * o = b.other();
    if (a2 == o->atom_number())
    {
      return &b;
    }
  }

  return NULL;
}


int
Target_Atom::in_same_rings (Target_Atom * a)
{
  if (0 == a->nrings())
    return 0;

  if (0 == nrings())
    return 0;

  return _m->in_same_rings(_my_atom_number, a->atom_number());
}

int
Target_Atom::is_bonded_to (Target_Atom * other) const
{
  for (int i = 0; i < _ncon; i++)
  {
    if (other == _other[i].other())
      return 1;
  }

  return 0;
}

/*int
Target_Atom::isolated_ring ()
{
  if (! is_ring_atom())
    return 0;

  if (nrings() > 1)
    return 0;

  if (TARGET_ATOM_NOT_COMPUTED == _isolated_ring)
  {
    if (_m->is_part_of_fused_ring_system(_my_atom_number))
      _isolated_ring = 0;
    else
      _isolated_ring = 1;
  }

  return _isolated_ring;
}*/

int
Target_Atom::fused_system_size ()
{
  if (! is_ring_atom())
    return 0;

  if (TARGET_ATOM_NOT_COMPUTED == _fused_system_size)
    _fused_system_size = _m->fused_system_size(_my_atom_number);

  return _fused_system_size;
}

int
Target_Atom::fragment_membership ()
{
  return _m->fragment_membership(_my_atom_number);
}

void
Target_Atom::establish_aromatic_bonds ()
{
  if (_ncon > 0 && NULL != _other)
  {
    for (int i = 0; i < _ncon; i++)
    {
      _other[i].establish_aromatic_bonds();
    }
  }

  return;
}

/*void
Target_Atom::establish_aromatic_bonds (int matoms,
                                       const int * in_fixed_kekule_ring)
{
  if (_ncon && NULL == _other)
    _allocate_other();

  cerr << "Atom " << _my_atom_number << " establish_aromatic_bonds, ncon " << _ncon << endl;
  for (int i = 0; i < _ncon; i++)
  {
    _other[i].establish_aromatic_bonds(matoms, _my_atom_number, in_fixed_kekule_ring);
  }

  return;
}*/

Bond_and_Target_Atom::Bond_and_Target_Atom ()
{
  _nbonds = _nrings = _aromatic = TARGET_ATOM_NOT_COMPUTED;

  _bond = NULL;
  _other = NULL;
}

void
Bond_and_Target_Atom::initialise (const Bond * b,
                                  Target_Atom & other)
{
  _bond = b;
  _other = &other;

  establish_aromatic_bonds();

  return;
}

/*void
Bond_and_Target_Atom::initialise (const Bond * b,
                                  Target_Atom & other,
                                  int nr)
{
//cerr << "Bond_and_Target_Atom::initialise: atoms " << b->a1() << " to " << b->a2() << " type " << b->btype() << endl;
  _bond = b;
  _other = &other;
  _nrings = nr;

  establish_aromatic_bonds();

  return;
}*/

/*
  After the Molecule_to_Match object has been created, and after initialise ()
  has been run, the caller may decide that we need to get aromaticity information
  from the molecule again - for example the query uses aromatic bond
  specifications. That's why this is public
*/

void
Bond_and_Target_Atom::establish_aromatic_bonds ()
{
  assert (NULL != _bond);

  if (_bond->is_aromatic())
  {
    _aromatic = 1;

    if (0 == _aromatic_bonds_lose_kekule_identity)
      _bond_types = (BOND_TYPE_ONLY_MASK & (_bond->btype()));
    else if (2 == _aromatic_bonds_lose_kekule_identity)
      _bond_types = (BOND_TYPE_ONLY_MASK & (_bond->btype()));
    else if (1 == _aromatic_bonds_lose_kekule_identity)
      _bond_types = AROMATIC_BOND;
    else
      _bond_types = (BOND_TYPE_ONLY_MASK & (_bond->btype()));
  }
  else
  {
    _aromatic = 0;
    _bond_types = (BOND_TYPE_ONLY_MASK & (_bond->btype()));
  }

//cerr << " bond " << _bond->a1() << " to " << _bond->a2() << " types " << _bond_types << endl;

  return;
}

/*
  More complex form where we are dealing with fixed kekule forms.
  The fixed_kekule_form array is an nbonds*nbonds array. Since bonds
  don't know their atom number, we get passed the atom number at
  one end of the bond, and ask the _other object to get the atom
  at the other end
*/

/*
void
Bond_and_Target_Atom::establish_aromatic_bonds (int matoms, 
                                                atom_number_t a1,
                                                const int * fixed_kekule_form)
{
  assert (NULL != _bond);
  assert (NULL != fixed_kekule_form);
  cerr << " Initialising bond, " << _aromatic_bonds_lose_kekule_identity << endl;

  if (! _bond->is_aromatic())
  {
    _bond_types = (BOND_TYPE_ONLY_MASK & (_bond->btype()));
    _aromatic = 0;
    return;
  }

  _aromatic = 1;

  if (0 == _aromatic_bonds_lose_kekule_identity)
    _bond_types = (BOND_TYPE_ONLY_MASK & (_bond->btype()));
  else if (1 == _aromatic_bonds_lose_kekule_identity)
    _bond_types = AROMATIC_BOND;
  else if (2 == _aromatic_bonds_lose_kekule_identity)
  {
    int ndx = a1 * matoms + _other->atom_number();
    cerr << "Atoms " << a1 << " and " << _other->atom_number() << " ndx " << ndx << " fixed? " << fixed_kekule_form[ndx] << endl;
    if (fixed_kekule_form[ndx])
      _bond_types = (BOND_TYPE_ONLY_MASK & (_bond->btype()));
    else
      _bond_types = AROMATIC_BOND;
  }

//cerr << " bond " << _bond->a1() << " to " << _bond->a2() << " types " << _bond_types << endl;

  return;
}
*/

int
Bond_and_Target_Atom::ok () const
{
  if (! _bond->ok())
    return 0;

  return _other->ok();
}

int
Bond_and_Target_Atom::nrings ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _nrings)
    _nrings = _bond->nrings();

  return _nrings;
}

/*boolean
Bond_and_Target_Atom::aromatic ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _aromatic)
    _aromatic = (0 != _bond->is_aromatic());    // ensure 0 or 1

  return _aromatic;
}*/

/*static void
initialise_bond_types (Molecule & m, 
                       bond_type_t * bt)
{
  m.compute_aromaticity_if_needed();

  int nb = m.nedges();

  int matoms = m.natoms();

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    bt[a1 * matoms + a2] = (b->btype() & BOND_TYPE_ONLY_MASK);
    bt[a2 * matoms + a1] = (b->btype() & BOND_TYPE_ONLY_MASK);
  }

  int nr = m.nrings();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    if (! is_fixed_kekule_form(m, *ri))
      continue;

    for (Ring_Bond_Iterator j(*ri); j != ri->end(); j++)
    {
      atom_number_t a1 = j.a1();
      atom_number_t a2 = j.a2();

      const Bond * b = m.bond_between_atoms(a1, a2);

      bond_type_t to_use;
      if (b->is_single_bond())
        to_use = SINGLE_BOND | AROMATIC_BOND;
      else
        to_use = DOUBLE_BOND | AROMATIC_BOND;

      bt[a1 * matoms + a2] = to_use;
      bt[a2 * matoms + a1] = to_use;
    }
  }
}*/

/*
  Sept 2000. Ran into problems with the code for matching aromatic
  bonds. I wanted lazy evaluation of aromaticity, but the bonds get
  their aromaticity in the constructor, so I give up and just ask
  the molecule to compute its aromaticity
*/

void
Molecule_to_Match::_initialise_molecule (Molecule * m)
{
  assert (m->ok());

  _m = m;

  _establish_aromatic_bonds_called = 0;

  if (_aromatic_bonds_lose_kekule_identity)
    m->compute_aromaticity_if_needed();

  _nrings = TARGET_ATOM_NOT_COMPUTED;
  _aromatic_ring_count     = TARGET_ATOM_NOT_COMPUTED;
  _non_aromatic_ring_count = TARGET_ATOM_NOT_COMPUTED;
  _fused_ring_count        = TARGET_ATOM_NOT_COMPUTED;
  _strongly_fused_ring_count = TARGET_ATOM_NOT_COMPUTED;
  _isolated_ring_count     = TARGET_ATOM_NOT_COMPUTED;
  _number_isotopic_atoms   = TARGET_ATOM_NOT_COMPUTED;
  _ring_object_count       = TARGET_ATOM_NOT_COMPUTED;

  _natoms = m->natoms();

  _target_atom = new Target_Atom[_natoms];

  _atom = new Atom *[_natoms];

  m->atoms((const Atom **) _atom);   // fetch the atoms

//_first = new_int(HIGHEST_ATOMIC_NUMBER + 1, INVALID_ATOM_NUMBER);
//_last  = new_int(HIGHEST_ATOMIC_NUMBER + 1, INVALID_ATOM_NUMBER);

  set_vector(_first, HIGHEST_ATOMIC_NUMBER + 1, INVALID_ATOM_NUMBER);   // no need to initialis last

  if (initialise_element_counts)
    set_vector(_count, HIGHEST_ATOMIC_NUMBER + 1, 0);

// We have two separate loops rather than putting the test for
// initialise_element_counts inside the loop. Warning, potential
// maintenance problem - code maintainability sacrificed for speed

  if (initialise_element_counts)
  {
    for (int i = 0; i < _natoms; i++)
    {
      Atom * a = _atom[i];
      _target_atom[i].initialise(m, i, a, this);

      if (! a->element()->is_in_periodic_table())
        continue;

      atomic_number_t z = a->atomic_number();

      if (INVALID_ATOM_NUMBER == _first[z])
        _first[z] = i;
      _last[z] = i;

      _count[z]++;
    }
  }
  else
  {
    for (int i = 0; i < _natoms; i++)
    {
      Atom * a = _atom[i];
      _target_atom[i].initialise(m, i, a, this);

      if (! a->element()->is_in_periodic_table())
        continue;

      atomic_number_t z = a->atomic_number();

      if (INVALID_ATOM_NUMBER == _first[z])
        _first[z] = i;
      _last[z] = i;
    }
  }

  _spinach_or_between_rings = NULL;

  _start_matching_at = INVALID_ATOM_NUMBER;

  if (use_fingerprints_for_screening_substructure_searches())
  {
#ifdef FINGERPRINT_SUBSTRUCTURE_SEARCHES
    IWMFingerprint iwmfingerprint;
    iwmfingerprint.construct_fingerprint(*m);
    _fingerprint = new IW_Bits_Base;
    (*_fingerprint) = iwmfingerprint;
#endif
  }
  else
    _fingerprint = NULL;

  return;
}

Molecule_to_Match::Molecule_to_Match (Molecule * m) : _m(m)
{
  _initialise_molecule(m);

  return;
}

Molecule_to_Match::Molecule_to_Match()
{
  _m = NULL;
  _target_atom = NULL;
  _atom = NULL;
  _spinach_or_between_rings = NULL;
  _fingerprint = NULL;

  _natoms = -1;

  _start_matching_at = INVALID_ATOM_NUMBER;

  return;
}

Molecule_to_Match::~Molecule_to_Match ()
{
  if (NULL != _target_atom)
  {
    delete [] _target_atom;
  }
    
  if (NULL != _atom)
    delete [] _atom;

  if (NULL != _spinach_or_between_rings)
    delete [] _spinach_or_between_rings;

  if (NULL != _fingerprint)
    delete _fingerprint;

  return;
}

int
Molecule_to_Match::ok () const
{
  return 1;    // for now
}

int
Molecule_to_Match::debug_print (std::ostream & os) const
{
  os << "Molecule_to_Match:natoms " << _natoms << endl;

  for (int i = 0; i < _natoms; i++)
  {
    _target_atom[i].debug_print(os);
//  os << " atom " << i << " type " << _target_atom[i].atomic_number() << " ncon " << _target_atom[i].ncon() << endl;
  }

  return os.good();
}
int
Molecule_to_Match::initialise_molecule (Molecule * m)
{
  assert (NULL == _m);

  _initialise_molecule(m);

  return 1;
}

/*
  In the fragmenter application, there is a need to pre-compute all
  aromaticity values, as the molecule from which the target was being
  derived is being changed, which would destroy aromaticity
*/

int
Molecule_to_Match::compute_aromaticity ()
{
  int rc = 1;
  for (int i = 0; i < _natoms; i++)
  {
    if (AROMATICITY_NOT_DETERMINED ==  _target_atom[i].aromaticity())
      rc = 0;
  }

  return rc;
}

/*static void
determine_in_fixed_kekule_ring(Molecule & m,
                               int * in_fixed_kekule_ring)
{
  int nr = m.nrings();

  m.compute_aromaticity_if_needed();

  int matoms = m.natoms();

  cerr << "Examining '" << m.name() << "' for fixed kekule ring forms\n";

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    if (! is_fixed_kekule_form(m, *ri))
      continue;

    for (Ring_Bond_Iterator j(*ri); j != ri->end(); j++)
    {
      atom_number_t a1 = j.a1();
      atom_number_t a2 = j.a2();

      in_fixed_kekule_ring[a1* matoms + a2] = 1;
      in_fixed_kekule_ring[a2* matoms + a1] = 1;
      cerr << "Atoms " << a1 << " and " << a2 << " in fixed kekule ring, size " << ri->number_elements() << endl;

    }
  }

  return;
}*/


/*
  We are being matched by a query that requires aromatic bond specifications.
*/

void
Molecule_to_Match::establish_aromatic_bonds ()
{
  if (_establish_aromatic_bonds_called)
    return;

  _m->compute_aromaticity_if_needed();

  for (int i = 0; i < _natoms; i++)
  {
    _target_atom[i].establish_aromatic_bonds();
  }

  _establish_aromatic_bonds_called = 1;

  return;
}

int
Molecule_to_Match::first (atomic_number_t z) const
{
  assert (_first);

  if (z < 0 || z > HIGHEST_ATOMIC_NUMBER)
    return -1;
  
  return _first[z];
}

int
Molecule_to_Match::last (atomic_number_t z) const
{
  assert (_last);
  
  return _last[z];
}

int
Molecule_to_Match::nrings ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _nrings)
    _nrings = _m->nrings();

  return _nrings;
}

int
Molecule_to_Match::_determine_ring_counts ()
{
  (void) nrings();

  _aromatic_ring_count  = 0;
  _non_aromatic_ring_count = 0;
  _isolated_ring_count = 0;
  _fused_ring_count = 0;
  _strongly_fused_ring_count = 0;

  if (0 == _nrings)
    return 1;

  (void) _target_atom[0].aromaticity();    // force aromaticity determination

  _rings.resize(_nrings);
  for (int i = 0; i < _nrings; i++)
  {
    const Ring * r = _m->ringi(i);
    if (r->is_aromatic())
      _aromatic_ring_count++;
    else
      _non_aromatic_ring_count++;

    if (r->is_fused())
    {
      _fused_ring_count++;
      if (r->largest_number_of_bonds_shared_with_another_ring() > 1)
        _strongly_fused_ring_count++;
    }
    else
      _isolated_ring_count++;

    _rings.add(r);
  }

  return 1;
}

/*
  A ring object is either an isolated ring, or a ring system
*/

int
Molecule_to_Match::ring_object_count ()
{
  if (TARGET_ATOM_NOT_COMPUTED != _ring_object_count)
    return _ring_object_count;

  _ring_object_count = 0;

  if (nrings() < 2)
    return _nrings;

  _determine_ring_counts();

  if (_isolated_ring_count == _nrings)   // all isolated
    return _nrings;

// Deal with the cases where there is at least some fustion

  assert (_nrings - _isolated_ring_count >= 2);

  if (2 == _nrings)     // must be a single two member ring system
    return 1;

  if (3 == _nrings)
  {
    if (0 == _isolated_ring_count)    // must be 3
      return 1;
    return 2;                         // must be 1 + 2
  }

  if (4 == _nrings)
  {
    if (1 == _isolated_ring_count)    // must be 1 and 3
      return 2;
    else if (2 == _isolated_ring_count)    // must be 2 and 1 and 1
      return 3;
  }

// Now we need to count them

  int * ring_already_done = new_int(_nrings);

//cerr << "Checking " << _nrings << " rings\n";
  for (int i = 0; i < _nrings; i++)
  {
    if (ring_already_done[i])
      continue;

    const Ring * r = _rings[i];

    _ring_object_count++;

//  cerr << "Ring " << i << " fused is " << r->is_fused() << endl;
    if (! r->is_fused())
      continue;

    for (int j = i + 1; j < _nrings; j++)
    {
      if (_rings[j]->fused_system_identifier() == r->fused_system_identifier())
      {
        ring_already_done[j] = 1;
      }
    }
  }

  delete ring_already_done;

  return _ring_object_count;
}

int
Molecule_to_Match::aromatic_ring_count ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _aromatic_ring_count)
    _determine_ring_counts();

  return _aromatic_ring_count;
}

int
Molecule_to_Match::non_aromatic_ring_count ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _non_aromatic_ring_count)
    _determine_ring_counts();

  return _non_aromatic_ring_count;
}

const Ring *
Molecule_to_Match::ringi (int i)
{
  if (TARGET_ATOM_NOT_COMPUTED == _aromatic_ring_count)
    _determine_ring_counts();

  assert (_rings.ok_index(i));

  return _rings[i];
}

int
Molecule_to_Match::fused_ring_count ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _fused_ring_count)
    _determine_ring_counts();

  return _fused_ring_count;
}

int
Molecule_to_Match::strongly_fused_ring_count ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _strongly_fused_ring_count)
    _determine_ring_counts();

  return _strongly_fused_ring_count;
}

int
Molecule_to_Match::isolated_ring_count ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _isolated_ring_count)
    _determine_ring_counts();

  return _isolated_ring_count;
}

int
Molecule_to_Match::number_isotopic_atoms ()
{
  if (TARGET_ATOM_NOT_COMPUTED == _number_isotopic_atoms)
    _number_isotopic_atoms = _m->number_isotopic_atoms();

  return _number_isotopic_atoms;
}

int
Molecule_to_Match::number_isotopic_atoms (int iso)
{
  return _m->number_isotopic_atoms(iso);
}

int
Molecule_to_Match::heteroatom_count () const
{
  int rc = 0;

  for (int i = 0; i < _natoms; i++)
  {
    atomic_number_t z = _target_atom[i].atomic_number();

    if (6 == z || 1 == z)
      ;
    else
      rc++;
  }

  return rc;
}

int
Molecule_to_Match::atoms_with_atomic_number (atomic_number_t z) const
{
  if (INVALID_ATOM_NUMBER == _first[z])
    return 0;

  if (initialise_element_counts)
    return _count[z];

  int rc = 0;


/*cerr << "First " << z << " at " << _first[z] << " to " << _last[z] << endl;

  for (int i = 0; i < _natoms; i++)
  {
    cerr << " i = " << i << " atomic number " << _target_atom[i].atomic_number() << endl;
  }*/

  for (int i = _first[z]; i <= _last[z]; i++)
  {
    if (z == _target_atom[i].atomic_number())
      rc++;
  }

  return rc;
}

void
Molecule_to_Match::invalidate (const Set_of_Atoms & e)
{
  for (int i = 0; i < e.number_elements(); i++)
  {
    atom_number_t j = e[i];

    _target_atom[j].invalidate();
  }

  return;
}

int
Molecule_to_Match::is_spinach (atom_number_t a)
{
  if (NULL == _spinach_or_between_rings)
    _initialise_spinach_or_between_rings();

  if (_spinach_or_between_rings[a] <= 0)
    return 0;

  return _spinach_or_between_rings[a];
}

int
Molecule_to_Match::is_between_rings (atom_number_t a)
{
  if (NULL == _spinach_or_between_rings)
    _initialise_spinach_or_between_rings();

  return TARGET_BETWEEN_RING == _spinach_or_between_rings[a];
}

int
Molecule_to_Match::_initialise_spinach_or_between_rings ()
{
  assert (NULL == _spinach_or_between_rings);

  _spinach_or_between_rings = new_int(_natoms, TARGET_UNSET);

  int nr = _m->nrings();

  if (0 == nr)    // no spinach, no between rings
    return 1;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = _m->ringi(i);
    ri->set_vector(_spinach_or_between_rings, TARGET_IS_RING);
  }

// First identify the between rings atoms.  This is inefficient
// because we may end up looking down each bond twice, but the
// searching is fundamentally different - for spinach we gather all
// atoms down the bond, for between ring features, we take a path to a
// ring atom

  for (int i = 0; i < _natoms; i++)
  {
    if (TARGET_IS_RING != _spinach_or_between_rings[i])   // bonds between rings start with ring atoms
      continue;

    const Atom * a = _atom[i];

    int acon = a->ncon();

    if (2 == acon)    // 2 connections, no branches off this atom
      continue;

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

      if (TARGET_IS_RING == _spinach_or_between_rings[k])
        continue;

      _identify_between_rings(i, k);
    }
  }

  for (int i = 0; i < _natoms; i++)
  {
    if (TARGET_IS_RING == _spinach_or_between_rings[i])    // spinach starts with ring atoms
      ;
    else if (TARGET_BETWEEN_RING == _spinach_or_between_rings[i])    // between ring atoms start on rings
      ;
    else               // no spinach starting here
      continue;

//  We have an atom that is between rings

    const Atom * a = _m->atomi(i);

    int acon = a->ncon();

    if (2 == acon)    // only > 2 connections can grow out spinach
      continue;

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

      if (TARGET_IS_RING == _spinach_or_between_rings[k])
        continue;

      if (TARGET_BETWEEN_RING == _spinach_or_between_rings[k])
        continue;

      int atoms_in_spinach = _grow_spinach(k, TARGET_SPTMP);

//    cerr << "From atom " << i << " find " << atoms_in_spinach << " atoms_in_spinach\n";

      _spinach_or_between_rings[k] = atoms_in_spinach;
      for (int l = 0; l < _natoms; l++)
      {
        if (TARGET_SPTMP == _spinach_or_between_rings[l])
          _spinach_or_between_rings[l] = atoms_in_spinach;
      }
    }
  }

//#define DEBUG_INITIALISE_SPINACH_OR_BETWEEN_RINGS
#ifdef DEBUG_INITIALISE_SPINACH_OR_BETWEEN_RINGS
  for (int i = 0; i < _natoms; i++)
  {
    cerr << "Atom " << i << " value " << _spinach_or_between_rings[i] << endl;
  }
#endif

  return 1;
}

int
Molecule_to_Match::_grow_spinach (atom_number_t zatom,
                                  int flag)
{
  _spinach_or_between_rings[zatom] = flag;

  const Atom * a = _atom[zatom];

  int acon = a->ncon();

  int rc = 1;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (TARGET_UNSET == _spinach_or_between_rings[j])
      rc += _grow_spinach(j, flag);
  }

  return rc;
}

/*
  Note that we don't find the shortest distance between the rings, just the first one we find
*/

int
Molecule_to_Match::_identify_between_rings (atom_number_t aprev,
                                            atom_number_t zatom)
{
  const Atom * a = _atom[zatom];

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (j == aprev)    // forward only
      continue;

    if (TARGET_IS_RING == _spinach_or_between_rings[j] || _identify_between_rings(zatom, j))
    {
      _spinach_or_between_rings[zatom] = TARGET_BETWEEN_RING;
      return 1;
    }
  }

  return 0;
}

int
Molecule_to_Match::is_superset (const IW_Bits_Base & fp) const
{
  assert (NULL != _fingerprint);

  return fp.is_subset(*_fingerprint);
}

void
Molecule_to_Match::discard_chirality_information()
{
  for (int i = 0; i < _natoms; i++)
  {
    _target_atom[i].discard_chirality();
  }

  return;
}
int
Target_Atom::is_spinach () 
{ 
  return _target->is_spinach(_my_atom_number);
}

Target_Atom &
Target_Atom::atom (int i) const
{
  return _target->operator[](i);
}
