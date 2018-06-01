#include <memory>
#include <iostream>
using std::cerr;
using std::endl;

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "misc.h"

/*
  We need some private molecule functions
*/

#define COMPILING_SMILES_CC

#include "molecule.h"
#include "path.h"
#include "chiral_centre.h"

static int copy_user_specified_atom_void_ptrs_during_create_subset = 0;

void
set_copy_user_specified_atom_void_ptrs_during_create_subset (int s)
{
  copy_user_specified_atom_void_ptrs_during_create_subset = s;
}

/*
  Helper function for transferring directional bonding.
  A will be an unsaturated atom at one end of a directional double bond.
  Identify the atoms atttached to the other end of the double bond
*/

#ifdef OLD_TRANSFER_DIRECTIONAL_BOND_INFO
int
Molecule::_transfer_directional_bond_info (const Atom * a, 
                                           atom_number_t zatom,
                                           int * transfer_bond_directionality) const
{
  int a1con = a->ncon();

  for (int i = 0; i < a1con; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    atom_number_t a2 = b->other(zatom);   // atom at other end of double bond

    const Atom * aa2 = _things[a2];

    int a2con = aa2->ncon();

    int rc = 0;

    for (int j = 0; j < a2con; j++)
    {
      const Bond * b2 = aa2->item(j);

      if (! b2->is_directional())
        continue;

//    If the bond number is available, use it

      int k = b2->bond_number();

      if (k < 0)
        k = which_bond(a2, b2->other(a2));

      transfer_bond_directionality[k] = 2;
      rc++;
    }
    
    return rc;
  }

  return 0;     // will never come here
}
#endif

/*
  zatom
*/

/*
  one or both ends of a cis-trans bond will not be in the subset.
  Unset any directional bonds attached. But beware of bonds that
  might be part of a valid adjacent cis-trans bond
*/

void
Molecule::_invalidate_attatched_directional_bonds (Molecule & subset,
                                                   const int * xref,
                                                   const Bond * b) const
{
  cerr << "Molecule::_invalidate_attatched_directional_bonds:not implemented, see Ian\n";

  atom_number_t a1 = b->a1();
  atom_number_t a2 = b->a2();

  if (xref[a1] >= 0)
  {
  }

  if (xref[a2] >= 0)
  {
  }

  return;
}

void
Molecule::_set_directional_bonds (Molecule & subset,
                                  atom_number_t zatom,
                                  const int * xref) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_directional())
      continue;

    atom_number_t j = b->other(zatom);
    if (xref[j] < 0)
      continue;

    if (b->is_directional_up())
      subset._set_bond_directionality(xref[zatom], xref[j], IW_BOND_DIRECTIONAL_UP);
    else
      subset._set_bond_directionality(xref[zatom], xref[j], IW_BOND_DIRECTIONAL_DOWN);
  }

  return;
}

int
Molecule::_count_directional_attachments (atom_number_t zatom) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_directional())
      continue;

    rc++;
  }

  return rc;
}

/*
  We are creating a subset and there are some bonds which are directional
  which have been included in the subset. Which directional bonds need to
  be transferred
*/


int
Molecule::_transfer_directional_bond_info (Molecule & subset,
                                           const int * xref) const
{
// We need to first place all the bonds that are known to be OK.
// Then we check any that might need to be invalidated

  int nb = _bond_list.number_elements();

  resizable_array<int> directional_bonds_to_invalidate;

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];

    if (! b->is_double_bond())
      continue;

    if (! b->part_of_cis_trans_grouping())
      continue;

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

//  cerr << "Examining double directional bond between " << a1 << " (" << xref[a1] << ") and " << a2 << " (" << xref[a2] << ")\n";

    if (xref[a1] < 0 || xref[a2] < 0)
    {
      directional_bonds_to_invalidate.add(i);
      continue;
    }

    if (0 == _count_directional_attachments(a1) || 0 == _count_directional_attachments(a2))
    {
      directional_bonds_to_invalidate.add(i);
      continue;
    }

    _set_directional_bonds(subset, a1, xref);
    _set_directional_bonds(subset, a2, xref);

    Bond * bs = const_cast<Bond *>(subset.bond_between_atoms(xref[a1], xref[a2]));
    bs->set_part_of_cis_trans_grouping(1);
  }

  for (int i = 0; i < directional_bonds_to_invalidate.number_elements(); i++)
  {
    int j = directional_bonds_to_invalidate[i];

    const Bond * b = _bond_list[j];

    _invalidate_attatched_directional_bonds(subset, xref, b);
  }

  return 1;
}

int
Molecule::_transfer_wedge_bond_info (Molecule & subset,
                                     const int * xref) const
{
// We need to first place all the bonds that are known to be OK.
// Then we check any that might need to be invalidated

  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];

    if (! b->is_single_bond())
      continue;

    atom_number_t a1 = b->a1();

    if (xref[a1] < 0)
      continue;

    atom_number_t a2 = b->a2();

    if (xref[a2] < 0)
      continue;

    if (! b->is_wedge_any())
      continue;

    if (b->is_wedge_up())
      subset.set_wedge_bond_between_atoms(xref[a1], xref[a2], 1);
    else if (b->is_wedge_down())
      subset.set_wedge_bond_between_atoms(xref[a1], xref[a2], -1);
  }

  return 1;
}

static void
create_ring_with_cross_referenced_atom_number (const Ring & rfrom,
                                               Ring & rto,
                                               const int * xref)
{
  int n = rfrom.number_elements();

  const atom_number_t * a = rfrom.rawdata();

  for (int i = 0; i < n; i++)
  {
    rto.add(xref[a[i]]);
  }

  return;
}

int
Molecule::create_components (int bond_number,
                             Molecule & m1,
                             Molecule & m2)
{
  assert(1 == number_fragments());
  assert(0 == m1.natoms());
  assert(0 == m2.natoms());

  compute_aromaticity_if_needed();

  const Bond * b = _bond_list[bond_number];
  assert(0 == b->nrings());

  atom_number_t break1 = b->a1();
  atom_number_t break2 = b->a2();

//cerr << "Bond " << bond_number << " atoms " << break1 << " and " << break2 << endl;

// Lazy, we need 3 arrays, just allocate one

  int * tmp = new_int(_number_elements * 2); std::unique_ptr<int[]> free_tmp(tmp);

  identify_side_of_bond(tmp, break1 , 1, break2);

  int * xref = tmp + _number_elements;

  int ndx1 = 0;
  int ndx2 = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == tmp[i])
    {
      xref[i] = ndx1;
      m1.add(new Atom(_things[i]));
      ndx1++;
    }
    else
    {
      xref[i] = ndx2;
      m2.add(new Atom(_things[i]));
      ndx2++;
    }
  }

  for (int i = 0; i < m1._number_elements; i++)
  {
    m1._things[i]->set_modified();
  }

  for (int i = 0; i < m2._number_elements; i++)
  {
    m2._things[i]->set_modified();
  }

  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; i++)
  {
    if (i == bond_number)
      continue;

    const Bond * b = _bond_list[i];

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (0 == tmp[a1])
      m1.add_bond(xref[a1], xref[a2], b->btype(), 1);
    else
      m2.add_bond(xref[a1], xref[a2], b->btype(), 1);
  }

  m1._fragment_information.all_atoms_in_one_fragment(ndx1, m1._bond_list.number_elements());
  m2._fragment_information.all_atoms_in_one_fragment(ndx2, m2._bond_list.number_elements());

// Do ring data

  int nr = _sssr_rings.number_elements();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = _sssr_rings[i];

    atom_number_t r0 = ri->item(0);

    Ring * r = new Ring;
    create_ring_with_cross_referenced_atom_number(*ri, *r, xref);

    if (0 == tmp[r0])
      m1._sssr_rings.add(r);
    else
      m2._sssr_rings.add(r);
  }

  m1._nrings = m1._sssr_rings.number_elements();
  m2._nrings = m2._sssr_rings.number_elements();
  m1._number_sssr_rings = m1._sssr_rings.number_elements();
  m2._number_sssr_rings = m2._sssr_rings.number_elements();

  m1._ring_membership = new int[ndx1];
  m2._ring_membership = new int[ndx2];
  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == tmp[i])
      m1._ring_membership[xref[i]] = _ring_membership[i];
    else
      m2._ring_membership[xref[i]] = _ring_membership[i];
  }

//cerr << "m1 has " << m1._nrings << " rings, m2 has " << m2._nrings << endl;

  if (m1._nrings)
  {
    m1._aromaticity = new aromaticity_type_t[ndx1];
    for (int i = 0; i < _number_elements; i++)
    {
      if (0 == tmp[i])
        m1._aromaticity[xref[i]] = _aromaticity[i];
    }
  }

  if (m2._nrings)
  {
    m2._aromaticity = new aromaticity_type_t[ndx2];
    for (int i = 0; i < _number_elements; i++)
    {
      if (1 == tmp[i])
        m2._aromaticity[xref[i]] = _aromaticity[i];
    }
  }

  assert(m1.ok());
  assert(m2.ok());

  return 1;
}

/*static void
translate_atom_numbers (Chiral_Centre & c,
                        const int initial_natoms,
                         const int * xref)
{
  c.set_centre(xref[c.a()] - initial_natoms);

  atom_number_t a = c.top_front();
  if (a >= 0)
    c.set_top_front(xref[a] - initial_natoms);

  a = c.top_back();
  if (a >= 0)
    c.set_top_back(xref[a] - initial_natoms);

  a = c.left_down();
  if (a >= 0)
    c.set_left_down(xref[a] - initial_natoms);

  a = c.right_down();
  if (a >= 0)
    c.set_right_down(xref[a] - initial_natoms);

  return;
}*/

/*
  There are lots of use cases of create_components where all we really want
  is to separate the fragments. Leaving one fragment in the initial molecule
  and just creating a separate, single, molecule is fine.

  ZATOM will be an atom in the fragment that will remain in THIS
*/

//#define DEBUG_SPLIT_OFF_FRAGMENTS

int
Molecule::split_off_fragments (const atom_number_t zatom, Molecule & frags)
{
#ifdef DEBUG_SPLIT_OFF_FRAGMENTS
  cerr << "Molecule::split_off_fragments start, parent contains " << nrings() << " rings\n";
#endif
  if (_charges || _atom_type)
  {
    cerr << "Molecule::split_off_fragments:does not work with molecules having partial charges and/or atom types. See Ian\n";    // laziness
    return 0;
  }

  assert(number_fragments() > 1);
  assert(0 == frags.natoms());

  int * tmp = new_int(_number_elements); std::unique_ptr<int[]> free_tmp(tmp);

  identify_side_of_bond(tmp, zatom, 1, INVALID_ATOM_NUMBER);

  const int initial_natoms = _number_elements;

#ifdef DEBUG_SPLIT_OFF_FRAGMENTS
  for (auto i = 0; i < initial_natoms; ++i)
  {
    cerr << " atom " << i << " atomic number " << atomic_number(i) << " side " << tmp[i] << endl;
  }
#endif

  int atoms_in_residual = 0;
  int atoms_in_fragments = 0;

  for (int i = 0; i < initial_natoms; ++i)
  {
    if (1 == tmp[i])
    {
      tmp[i] = atoms_in_residual;
      atoms_in_residual++;
    }
    else
    {
      tmp[i] = initial_natoms + atoms_in_fragments;
      atoms_in_fragments++;
    }

#ifdef DEBUG_SPLIT_OFF_FRAGMENTS
    cerr << " xref " << i << " " << tmp[i] << endl;
#endif
  }

  return _common_fragment_extraction(tmp, initial_natoms, atoms_in_residual, atoms_in_fragments, frags);
}

int
Molecule::excise_fragment (const atom_number_t zatom, Molecule & frag)
{
  assert(number_fragments() > 1);
  assert(0 == frag.natoms());

  int * tmp = new_int(_number_elements); std::unique_ptr<int[]> free_tmp(tmp);

  identify_side_of_bond(tmp, zatom, 1, INVALID_ATOM_NUMBER);    // just doing a cheap fragment membership. tmp = 1 means being removed

  const int initial_natoms = _number_elements;

#ifdef DEBUG_EXCISE_FRAGMENT
  for (auto i = 0; i < initial_natoms; ++i)
  {
    cerr << " atom " << i << " atomic number " << atomic_number(i) << " side " << tmp[i] << endl;
  }
#endif

  int atoms_in_residual = 0;
  int atoms_in_fragment = 0;

  for (int i = 0; i < initial_natoms; ++i)
  {
    if (0 == tmp[i])     // in the residual
    {
      tmp[i] = atoms_in_residual;
      atoms_in_residual++;
    }
    else
    {
      tmp[i] = initial_natoms + atoms_in_fragment;
      atoms_in_fragment++;
    }

#ifdef DEBUG_EXCISE_FRAGMENT
    cerr << " xref " << i << " " << tmp[i] << endl;
#endif
  }

  return _common_fragment_extraction(tmp, initial_natoms, atoms_in_residual, atoms_in_fragment, frag);
}

int
Molecule::_common_fragment_extraction (int * tmp,
                                       const int initial_natoms,
                                       const int atoms_in_residual,
                                       const int atoms_in_fragment,
                                       Molecule & f)
{
  assert(atoms_in_residual + atoms_in_fragment == initial_natoms);

  if (0 == atoms_in_fragment || 0 == atoms_in_residual)
  {
    cerr << "Molecule::excise_fragment:all atoms in one frgament\n";
    return 0;
  }

  f.resize(atoms_in_fragment);

  for (int i = initial_natoms - 1; i >= 0; --i)
  {
    if (tmp[i] < initial_natoms)    // being retained
      continue;

    Atom * a = this->remove_no_delete(i);
    f.add(a);
  }

  if (f._number_elements > 1)    // reverse the array to preserve atom numbering
  {
    for (int i = 0, j = f._number_elements - 1; i < j; ++i, --j)
    {
      std::swap(f._things[i], f._things[j]);
    }
  }

  for (int i = _bond_list.number_elements() - 1; i >= 0; --i)
  {
    Bond * b = _bond_list[i];

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();
//  cerr << " bond in starting molecule btw " << a1 << " (" << tmp[a1] << ") and " << a2 << " (" << tmp[a2] << ")\n";

    if (tmp[a1] < initial_natoms)    // being retained, just renumber
    {
      b->set_a1a2(tmp[a1], tmp[a2]);
//    cerr << "  bond in residual btw " << b->a1() << " and " << b->a2() << endl;
    }
    else              // needs to go to the new molecule
    {
      b = _bond_list.remove_no_delete(i);
      b->set_a1a2(tmp[a1] - initial_natoms, tmp[a2] - initial_natoms);
      f._bond_list.add(b);
//    cerr << "  bond in fragment btw " << b->a1() << " and " << b->a2() << endl;
    }
  }

// for chiral centres, we need numbering in each molecule

  const int nc = _chiral_centres.number_elements();
  if (nc)
  {
    for (int i = nc - 1; i >= 0; --i)
    {
      const atom_number_t c = _chiral_centres[i]->a();

      if (tmp[c] >= initial_natoms)    // must go to f
      {
        Chiral_Centre * c = _chiral_centres.remove_no_delete(i);
        f._chiral_centres.add(c);
      }
      else
        _chiral_centres[i]->adjust_atom_numbers(tmp);
    }

    const auto ncf = f._chiral_centres.number_elements();
    if (ncf)
    {
      for (int i = 0; i < initial_natoms; ++i)
      {
        if (tmp[i] >= initial_natoms)
          tmp[i] -= initial_natoms;
      }

      for (int i = 0; i < ncf; ++i)
      {
        Chiral_Centre * c = f._chiral_centres[i];
        c->adjust_atom_numbers(tmp);
      }
    }
  }

  _set_modified();     // rings and aromaticity too hard, just recompute

#ifdef DEBUG_EXCISE_FRAGMENT
  cerr << "excise_fragment:residual contains " << _number_elements << " atoms and " << _bond_list.size() << " bonds. Frags " << f._number_elements << " and " << f._bond_list.size() << " bonds\n";
  this->debug_print(cerr);
  cerr << "FRAG\n";
  f.debug_print(cerr);
#endif

  assert(f.ok());
  assert(ok());

  return atoms_in_fragment;
}

/*
  We don't play nicely with the molecular consistency here, in that during
  the intermediate steps, the Molecule is NOT in a consistent state, for
  example, bonds are removed, but we have not called set_modified()
*/

int
Molecule::create_components_across_bonds (const int * rmbond,
                                          resizable_array_p<Molecule> & components)
{
  int ne = nedges();

  resizable_array<Bond *> add_back_later;

  _bond_list.assign_bond_numbers_to_bonds_if_needed();

  ring_membership();

  for (int i = ne - 1; i >=0; i--)
  {
    if (0 == rmbond[i])
      continue;

    Bond * b = _bond_list[i];

    _bond_list.remove_no_delete(i);

    add_back_later.add(b);

    _things[b->a1()]->remove_first(b);
    _things[b->a2()]->remove_first(b);
  }

  int n = add_back_later.number_elements();

  if (0 == n)
    return 0;

//_fragment_information.debug_print(cerr);
  Fragment_Information frag_info;

  _compute_fragment_information(frag_info, 0);   // 0 means do NOT update ring info

  int nf = frag_info.number_fragments();

  const int * in_frag = frag_info.fragment_membership();
//_fragment_information.debug_print(cerr);

  for (int i = 1; i <= nf; i++)
  {
    Molecule * m = new Molecule;
    create_subset(*m, in_frag, i);
    components.add(m);
  }

//_fragment_information.debug_print(cerr);
  for (int i = n - 1; i >= 0; i--)
  {
    Bond * abl = add_back_later[i];

    int bn = abl->bond_number();

    if (bn == _bond_list.number_elements())
      _bond_list.add(abl);
    else
      _bond_list.insert_before(bn, abl);

    _things[abl->a1()]->add(abl);
    _things[abl->a2()]->add(abl);
  }

  return nf;
}

int
Molecule::create_subset (Molecule & subset,
                         const int * process_these,
                         int id,
                         int * xref) const
{
  int ndx = subset.natoms();    // we can start with a non-empty molecule
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (id != process_these[i])
    {
      xref[i] = -1;
      continue;
    }

    const Atom * a = _things[i];
    Atom * as = new Atom(a);
    subset.Molecule::add(as);

    if (copy_user_specified_atom_void_ptrs_during_create_subset)    // actually it always happens in the Atom copy constructor...
      as->set_user_specified_void_ptr(const_cast<void *>(a->user_specified_void_ptr()));

    xref[i] = ndx;

    ndx++;
    rc++;

//  cerr << "Atom " << i << " becomes subset atom " << xref[i] << endl;

    if (_charges)
    {
      charge_t q = charge_on_atom(i);
      if (q)
        subset.set_charge(xref[i], q);
    }
  }

  int directional_bonds_may_be_present = 0;
  int wedge_bonds_may_be_present = 0;

  int nb = _bond_list.number_elements();
  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];

    atom_number_t a1 = b->a1();
    if (xref[a1] < 0)
      continue;

    atom_number_t a2 = b->a2();
    if (xref[a2] < 0)
      continue;

    subset.add_bond(xref[a1], xref[a2], b->btype());

    if (b->is_directional())
      directional_bonds_may_be_present++;
    if (b->is_wedge_any())
      wedge_bonds_may_be_present++;
  }

//if (0 == subset.natoms())
//  cerr << "Molecule::create_subset: warning empty molecule created\n";

// What if there are chiral centres

  int nc = _chiral_centres.number_elements();
  if (0 == nc && ! directional_bonds_may_be_present && ! wedge_bonds_may_be_present)
  {
    subset.recompute_implicit_hydrogens();
    return rc;
  }

// Need to process any chiral centres for which the centre atom is in the subset

#ifdef DEBUG_CREATE_SUBSET
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << " i = " << i << " process_these " << process_these[i] << " xref " << xref[i] << endl;
  }
#endif

  for (int i = 0; i < nc; i++)
  {
    const Chiral_Centre * c = _chiral_centres[i];
//  cerr << "Chiral centre at atom " << c->a() << " xref " << xref[c->a()] << endl;
    if (xref[c->a()] < 0)       // centre atom not a chiral centre in original molecule
      continue;

    atom_number_t j = xref[c->a()];     // atom number of the chiral centre in the subset molecule

//  cerr << "in subset " << c->all_atoms_in_subset(process_these, id) << endl;
//  c->debug_print(cerr);

    if (c->all_atoms_in_subset(process_these, id))
    {
      Chiral_Centre * nc = new Chiral_Centre(j);
      if (! nc->make_copy(*c, xref))
      {
        cerr << "Molecule::create_subset: yipes, make copy chiral centre failed\n";
        nc->debug_print(cerr);
        return 0;
      }

      subset._chiral_centres.add(nc);
    }
    else     // only some of the atoms are in the subset. No chiral centre in the subset molecule
    {
      Atom * a = subset._things[j];
      a->set_implicit_hydrogens_known(0);
      a->recompute_implicit_hydrogens(j);
    }
  }

// cis trans bonds also need to be added

  if (directional_bonds_may_be_present)
    _transfer_directional_bond_info(subset, xref);

  if (wedge_bonds_may_be_present)
    _transfer_wedge_bond_info(subset, xref);

  subset.recompute_implicit_hydrogens();
//subset.debug_print(cerr);
  return rc;
}

int
Molecule::create_subset (Molecule & subset,
                         const int * process_these,
                         int id) const
{
  assert(ok());
  assert(subset.ok());

  int * tmp = new_int(_number_elements, -1); std::unique_ptr<int[]> free_tmp(tmp);

  return create_subset(subset, process_these, id, tmp);
}

int
Molecule::create_components (const int * fragID,
                             resizable_array_p<Molecule> & components ) const
{
  if ( _number_elements <= 1 )
    return 0;

  int     ntmp = _number_elements;
  int     cntr = ntmp;
  int     fcntr = 0;
  while( cntr > 0 )
  {
    if (locate_item_in_array(fcntr, ntmp, fragID) < 0)
    {
      fcntr++;
      continue;
    }

    Molecule * m = new Molecule;
    create_subset( *m, fragID, fcntr );

    if ( 0 == m->natoms( ) )
    {
      delete m;
      continue;
    }                 // backup null check

    cntr -= m->natoms();
    fcntr ++;
    components.add( m );
  }

  return fcntr;
}

template <typename T>
int
Molecule::create_components (resizable_array_p<T> & components)
{
  assert(ok());

  if (0 == _number_elements)
    return 0;

  if (_number_elements <= 1)
  {
    cerr << "Molecule::create_components: molecule one or fewer atoms\n";
    return 0;
  }

  int nf = number_fragments();
  if (1 == nf)
  {
    cerr << "Molecule::create_components: molecule contains only one component\n";
    return 0;
  }

  int * tmp = new int[_number_elements]; std::unique_ptr<int[]> free_tmp(tmp);

  fragment_membership(tmp);

  assign_bond_numbers_to_bonds_if_needed();

  for (int i = 0; i < nf; i++)
  {
//  cerr << "Creating component from fragment " << i << endl;

    T * m = new T;
    create_subset(*m, tmp, i);
    components.add(m);
  }

  assert(nf == components.number_elements());

  return nf;
}

template int Molecule::create_components(resizable_array_p<Molecule>&);

int
Molecule::fragment_membership (int * f)
{
  assert(ok());
  assert(NULL != f);

  if (! _fragment_information.contains_valid_data())
    (void) number_fragments();

  copy_vector(f, _fragment_information.fragment_membership(), _number_elements);

  return number_fragments();
}

int
Molecule::fragment_membership (atom_number_t a)
{
  assert (ok_atom_number(a));

  if (! _fragment_information.contains_valid_data())
    (void) number_fragments();

  return _fragment_information.fragment_membership(a);
}

int
Molecule::number_fragments()
{
  assert(ok());

  if (_fragment_information.contains_valid_data())
    return _fragment_information.number_fragments();

  if (0 == _number_elements)
    return 0;

  return _compute_fragment_information(_fragment_information);
}

int
Molecule::atoms_in_fragment (int zfrag) 
{
  assert(zfrag >= 0);

  assert(zfrag < number_fragments());

  return _fragment_information.atoms_in_fragment(zfrag);
}

int
Molecule::atoms_in_fragment (Set_of_Atoms & aif,
                             int f)
{
  assert(ok());
  assert(0 == aif.number_elements());

  if (0 == _number_elements)
    return 0;

  int nf;

  if (! _fragment_information.contains_valid_data())
    nf = number_fragments();
  else
    nf = _fragment_information.number_fragments();

  if (1 == nf)
  {
    aif.resize(_number_elements);
    for (int i = 0; i < _number_elements; i++)
    {
      aif.add(i);
    }

    return _number_elements;
  }

  assert(f >= 0 && f < nf);

  _fragment_information.atoms_in_fragment(_number_elements, f, aif);

  return _fragment_information.atoms_in_fragment(f);
}

int
Molecule::atoms_in_largest_fragment()
{
  if (0 == _number_elements)
    return 0;

  const int f = largest_fragment();

  return _fragment_information.atoms_in_fragment(f);
}

int
Molecule::largest_fragment ()
{
  if (0 == _number_elements)
    return 0;

  const int nf = number_fragments();

  if (1 == nf)
    return _number_elements;

  int lf = _fragment_information.atoms_in_fragment(0);
  int rc = 0;

  if (2 == nf)    // a common case
  {
    if (_fragment_information.atoms_in_fragment(1) > lf)
      return 1;
    else
      return 0;
  }


  for (int i = 1; i < nf; i++)
  {
    if (_fragment_information.atoms_in_fragment(i) > lf)
    {
      lf = _fragment_information.atoms_in_fragment(i);
      rc = i;
    }
  }

  return rc;
}

/*
  Note that it might be a good idea to allow for duplicate detection in this
  function, or at least offer it as an option. Is it OK or not for there to
  be duplicate atom numbers in AIF. Perhaps the underlying Set_of_Atoms object
  could do that some day...
*/

int
Molecule::add_atoms_in_fragment (Set_of_Atoms & aif,
                                 int f)
{
  Set_of_Atoms tmp;

  atoms_in_fragment(tmp, f);

  aif += tmp;

  return atoms_in_fragment(f);
}

int
Molecule::reduce_to_largest_fragment()
{
  int nf = number_fragments();

  if (nf <= 1)
    return 1;

  int maxat = _fragment_information.atoms_in_fragment(0);    // largest number of atoms in a fragment

  int imax  = 0;                        // which fragment

// If the number of atoms in fragment 0 is > half the atoms in the molecule, we
// know this must be the largest fragment

  if (maxat <= _number_elements / 2)
  {
    for (int i = 1; i < nf; i++)
    {
      int j = _fragment_information.atoms_in_fragment(i);

      if (j > maxat)
      {
        maxat = j;
        imax = i;
      }
    }
  }

  Set_of_Atoms atoms_to_be_removed;
  atoms_to_be_removed.resize(maxat);

  const int * fragment_membership = _fragment_information.fragment_membership();

  for (int i = 0; i < _number_elements; i++)
  {
    if (fragment_membership[i] != imax)
      atoms_to_be_removed.add(i);
  }

  return remove_atoms(atoms_to_be_removed);
}

atom_number_t
Molecule::first_atom_in_fragment (int f)
{
  assert(f >= 0);

  int nf = number_fragments();

  if (1 == nf)
    return 0;

  assert(f < nf);

  const int * fragment_membership = _fragment_information.fragment_membership();

  for (int i = 0; i < _number_elements; i++)
  {
    if (fragment_membership[i] == f)
      return i;
  }

  cerr << "Molecule::first_atom_in_fragment:no atom found in fragment " << f << "\n";    // how could this happen
  
  return INVALID_ATOM_NUMBER;
}

int
Molecule::rings_in_fragment (int f)
{
  assert(f >= 0);

  assert(f < number_fragments());

  return _fragment_information.rings_in_fragment(f);
}

/*
  We are building a bond subset and have placed a bond. Turn off that bond
  in BOND_LOOKUP_TABLE
*/

static void
turn_off_bond (atom_number_t a1, 
               atom_number_t a2, 
               int * bond_lookup_table,
               int natoms)
{
  bond_lookup_table[a1 * natoms + a2] = 0;
  bond_lookup_table[a2 * natoms + a1] = 0;

  return;
}

/*
  Dangerous, we allocate an natoms*natoms array to put a 
  quick-access bond lookup table
*/

int
Molecule::create_subset_by_bond (Molecule & subset,
                                 const int * these_bonds_only,
                                 int flag) const
{
  subset.resize(0);

  int bond_lookup_table_size = _number_elements * (_number_elements - 1) + _number_elements + 1;

  int * bond_lookup_table = new_int(bond_lookup_table_size);

  if (NULL == bond_lookup_table)
  {
    cerr << "Molecule::create_subset_by_bond:cannot allocate bond lookup table " << _number_elements << endl;
    return 0;
  }

  std::unique_ptr<int[]> free_bond_lookup_table(bond_lookup_table);

  int * xref = new_int(_number_elements, -1);
  if (NULL == xref)
  {
    cerr << "Molecule::create_subset_by_bond:cannot allocate xref array " << _number_elements << endl;
    return 0;
  }

  std::unique_ptr<int[]> free_xref(xref);

  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; i++)
  {
    if (flag != these_bonds_only[i])   // bond I is being omitted
      continue;
  
    const Bond * b = _bond_list[i];
  
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    bond_lookup_table[a1 * _number_elements + a2] = 1;
    bond_lookup_table[a2 * _number_elements + a1] = 1;
  }

  for (int i = 0; i < bond_lookup_table_size; i++)
  {
    if (0 == bond_lookup_table[i])
      continue;

    atom_number_t a1 = i % _number_elements;

    _create_bond_subset_starting_with(subset, a1, bond_lookup_table, xref);
  }

  return subset._number_elements;
}

int
Molecule::_create_bond_subset_starting_with (Molecule & subset,
                                             atom_number_t zatom,
                                             int * bond_lookup_table,
                                             int * xref) const
{
  Set_of_Atoms atoms_to_be_processed;
  atoms_to_be_processed.resize(_number_elements);

  int x = subset._number_elements;    // the atom number it will get

  Atom * a = new Atom(_things[zatom]);
  a->set_implicit_hydrogens_known(0);
  subset.add(a);

//cerr << "Starting subset with atom " << zatom << endl;

  xref[zatom] = x;

  atoms_to_be_processed.add(zatom);

  while (atoms_to_be_processed.number_elements() > 0)
  {
    atom_number_t current_atom = atoms_to_be_processed.pop();

    assert(xref[current_atom] >= 0);

    const Atom * a = _things[current_atom];

    int acon = a->ncon();

    for (int i = 0; i < acon; i++)
    {
      const Bond * b = a->item(i);

      atom_number_t a1 = b->a1();
      atom_number_t a2 = b->a2();

      if (0 == bond_lookup_table[a1 * _number_elements + a2])
        continue;

      atom_number_t o = b->other(current_atom);

      if (xref[o] >= 0)   // already done that, need to put bond between them
      {
        if (! subset.are_bonded(xref[current_atom], xref[o]))
          subset.add_bond(xref[current_atom], xref[o], b->btype());
        turn_off_bond(current_atom, o, bond_lookup_table, _number_elements);
        continue;
      }

      int x = subset._number_elements;    // the atom number it will get in the new molecule

      Atom * t = new Atom(_things[o]);
      t->set_implicit_hydrogens_known(0);
      subset.add(t);

      subset.add_bond(xref[current_atom], x, b->btype());

      turn_off_bond(current_atom, o, bond_lookup_table, _number_elements);

      xref[o] = x;

      atoms_to_be_processed.add(o);
    }
  }

  for (int i = 0; i < subset._number_elements; i++)
  {
    int notused;
    subset._things[i]->recompute_implicit_hydrogens(notused);
  }

  return subset._number_elements;
}

Fragment_Information::Fragment_Information()
{
  _number_fragments = -1;

  _fragment_membership = NULL;

  return;
}

Fragment_Information::~Fragment_Information()
{
  if (NULL != _fragment_membership)
    delete [] _fragment_membership;

  _number_fragments = -2;

  return;
}

int
Fragment_Information::contains_valid_data() const
{
  return (-1 != _number_fragments);
}

int
Fragment_Information::debug_print (std::ostream & os) const
{
  if (_number_fragments >= 0)
    os << "Molecule contains " << _number_fragments << " components\n";

  if (_number_fragments > 0 && _atoms_in_fragment.number_elements())
  {
    for (int i = 0; i < _number_fragments; i++)
    {
      os << "Fragment " << i << ", " << _atoms_in_fragment[i] << " atoms, " << _bonds_in_fragment[i] << " bonds";

      os << " " << rings_in_fragment(i) << " rings\n";
    }
  }

  return os.good();
}

/*
  Build fragment membership without reference to smiles ordering
*/

int
Fragment_Information::initialise (int matoms)
{
  if (NULL == _fragment_membership)
  {
    _fragment_membership = new_int(matoms, FRAGMENT_MEMBERSHIP_NOT_SET);
    if (NULL == _fragment_membership)
    {
      cerr << "Fragment_Information::initialise:cannot allocate " << matoms << " atoms\n";
      return 0;
    }
  }
  else
    set_vector(_fragment_membership, matoms, FRAGMENT_MEMBERSHIP_NOT_SET);

  _number_fragments = -1;

  _atoms_in_fragment.resize_keep_storage(0);
  _bonds_in_fragment.resize_keep_storage(0);

  return 1;
}

int
Molecule::invalidate_fragment_membership()
{
// do NOT assert OK here, _set_modified_no_ok() will die.

  _fragment_information.invalidate();

  return 1;
}

void
Fragment_Information::invalidate()
{
  _number_fragments = -1;

  if (NULL != _fragment_membership)
  {
    delete [] _fragment_membership;
    _fragment_membership = NULL;
  }

  _atoms_in_fragment.resize_keep_storage(0);
  _bonds_in_fragment.resize_keep_storage(0);

  return;
}

int
Fragment_Information::all_atoms_in_one_fragment (int natoms, int nbonds)
{
  _number_fragments = 1;

  if (NULL != _fragment_membership)
    delete [] _fragment_membership;

  _fragment_membership = new_int(natoms, 0);

  _atoms_in_fragment.resize_keep_storage(0);
  _bonds_in_fragment.resize_keep_storage(0);

  _atoms_in_fragment.add(natoms);
  _bonds_in_fragment.add(nbonds);

  return 1;
}

int
Fragment_Information::set_number_fragments (int nf)
{
  assert(nf > 0);

  _number_fragments = nf;

  _atoms_in_fragment.extend(nf);
  _bonds_in_fragment.extend(nf);

  return 1;
}

int
Fragment_Information::atoms_in_fragment (int matoms,
                                         int f,
                                         Set_of_Atoms & s) const
{
  s.resize_keep_storage(0);

  assert(f >= 0 && f < _number_fragments);

  for (int i = 0; i < matoms; i++)
  {
    if (f == _fragment_membership[i])
      s.add(i);
  }

  return s.number_elements();
}

int
Molecule::_compute_fragment_information (Fragment_Information & fragment_information,
                                         int update_ring_info)
{
  if (fragment_information.contains_valid_data())
    return fragment_information.number_fragments();

  if (0 == _number_elements)
    return 1;

  fragment_information.initialise(_number_elements);

  int * fragment_membership = fragment_information.fragment_membership();

  Set_of_Atoms mystack;
  mystack.resize(_number_elements);

  int nf = 0;

  for (int i = 0; i < _number_elements; ++i)
  {
    if (FRAGMENT_MEMBERSHIP_NOT_SET != fragment_membership[i])
      continue;

    mystack.add(i);

    while (mystack.number_elements() > 0)
    {
      atom_number_t zatom = mystack.pop();

      if (fragment_membership[zatom] >= 0)
        continue;

      fragment_membership[zatom] = nf;

      const Atom * a = _things[zatom];

      int acon = a->ncon();
  
      for (int i = 0; i < acon; i++)
      {
        atom_number_t j = a->other(zatom, i);

        if (FRAGMENT_MEMBERSHIP_NOT_SET == fragment_membership[j])
          mystack.add(j);
      }
    }

    nf++;
  }

  assert(nf > 0);

  if (! fragment_information.set_number_fragments(nf))
    return 0;

  resizable_array<int> & atoms_in_fragment = fragment_information.atoms_in_fragment();
  resizable_array<int> & bonds_in_fragment = fragment_information.bonds_in_fragment();

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom * a = _things[i];

    int f = fragment_membership[i];

    atoms_in_fragment[f]++;

    bonds_in_fragment[f] += a->ncon();    // counts the bonds twice
  }

  int fragments_with_no_rings = 0;

  for (int i = 0; i < nf; i++)
  {
    bonds_in_fragment[i] = bonds_in_fragment[i] / 2;

    if (0 == fragment_information.rings_in_fragment(i))
      fragments_with_no_rings++;
  }

//#define DEBUG_COMPUTE_FRAGMENT_MEMBERSHIP
#ifdef DEBUG_COMPUTE_FRAGMENT_MEMBERSHIP
  cerr << "Fragment membership computed '" << _molecule_name << "' with " << _number_elements << " atoms has " << _number_fragments << " fragments\n";
  for (int i = 0; i < _number_elements; i++)
  {
    cerr << "Atom " << i << " in fragment " << fragment_membership[i] << endl;
  }
  for (int i = 0; i < nf; i++)
  {
    cerr << "Fragment " << i << " " << atoms_in_fragment[i] << " atoms and " << bonds_in_fragment[i] << " bonds, nr = " << _fragment_information.rings_in_fragment(i) << endl;
  }
#endif

// We can do some very quick ring membership stuff

  if (0 == fragments_with_no_rings)
    return nf;

  if (! update_ring_info)
    return nf;

  if (NULL == _ring_membership)
    _initialise_ring_membership();

  for (int i = 0; i < _number_elements; i++)
  {
    int f = fragment_membership[i];

    if (0 == fragment_information.rings_in_fragment(f))
      _ring_membership[i] = 0;
  }

  return nf;
}

static int 
next_available_atom (int needle,
                     int n,
                     const int * haystack,
                     const int * include_atom,
                     int & istart)
{
  assert(NULL != include_atom);

  for ( ; istart < n; istart++)
  {
    if (0 == include_atom[istart])
      continue;

    if (needle == haystack[istart])
    {
      istart++;
      return istart - 1;
    }
  }

  return -1;
}

static int
next_available_atom (int needle,
                     int n,
                     const int * haystack,
                     int & istart)
{
  while (istart < n)
  {
    if (needle == haystack[istart])
    {
      istart++;
      return istart - 1;
    }

    istart++;
  }

  return -1;
}

/*
  Split off to separate function if doing a subset
*/

int
Molecule::compute_fragment_information(Fragment_Information & fragment_information,
                                       const int * include_atom) const
{
  if (0 == _number_elements)
    return 1;

  fragment_information.initialise(_number_elements);

  if (NULL != include_atom)
    return _compute_fragment_information_subset(fragment_information, include_atom);

  int * fragment_membership = fragment_information.fragment_membership();

  Set_of_Atoms mystack;
  mystack.resize(_number_elements);

  int nf = 0;

  int istart = 0;

  while (1)
  {
    const atom_number_t i = next_available_atom(FRAGMENT_MEMBERSHIP_NOT_SET, _number_elements, fragment_membership, istart);

    if (i < 0)
      break;

    mystack.add(i);

    while (mystack.number_elements() > 0)
    {
      const atom_number_t zatom = mystack.pop();

      if (fragment_membership[zatom] >= 0)
        continue;

      fragment_membership[zatom] = nf;

      const Atom * a = _things[zatom];

      const int acon = a->ncon();
  
      for (int i = 0; i < acon; i++)
      {
        atom_number_t j = a->other(zatom, i);

        if (FRAGMENT_MEMBERSHIP_NOT_SET == fragment_membership[j])
          mystack.add(j);
      }
    }

    nf++;
  }

  if (0 == nf)
  {
    cerr << "Molecule::compute_fragment_information:subset contains no atoms\n";
    return 0;
  }
  assert(nf > 0);

  if (! fragment_information.set_number_fragments(nf))
    return 0;

  resizable_array<int> & atoms_in_fragment = fragment_information.atoms_in_fragment();
  resizable_array<int> & bonds_in_fragment = fragment_information.bonds_in_fragment();

  int * aif = atoms_in_fragment.rawdata();   // access raw data for speed
  for (int i = 0; i < _number_elements; i++)
  {
    const int f = fragment_membership[i];
    aif[f]++;
  }

  const int nb = nedges();

  int * bif = bonds_in_fragment.rawdata();    // for speed
  for (int i = 0; i < nb; ++i)
  {
    const Bond * b = _bond_list[i];

    const atom_number_t a = b->a1();    // no subsetting, no need to check other end of bond

    const int f = fragment_membership[a];

    bif[f]++;
  }

  return nf;
}

int
Molecule::_compute_fragment_information_subset(Fragment_Information & fragment_information,
                                               const int * include_atom) const
{
  int * fragment_membership = fragment_information.fragment_membership();

  Set_of_Atoms mystack;
  mystack.resize(_number_elements);

  int nf = 0;

  int istart = 0;

  while (1)
  {
    const atom_number_t i = next_available_atom(FRAGMENT_MEMBERSHIP_NOT_SET, _number_elements, fragment_membership, include_atom, istart);

    if (i < 0)
      break;

    mystack.add(i);

    while (mystack.number_elements() > 0)
    {
      const atom_number_t zatom = mystack.pop();

      if (fragment_membership[zatom] >= 0)
        continue;

      fragment_membership[zatom] = nf;

      const Atom * a = _things[zatom];

      const int acon = a->ncon();
  
      for (int i = 0; i < acon; i++)
      {
        atom_number_t j = a->other(zatom, i);

        if (0 == include_atom[j])
          continue;

        if (FRAGMENT_MEMBERSHIP_NOT_SET == fragment_membership[j])
          mystack.add(j);
      }
    }

    nf++;
  }

  if (0 == nf)
  {
    cerr << "Molecule::compute_fragment_information:subset contains no atoms\n";
    return 0;
  }
  assert(nf > 0);

  if (! fragment_information.set_number_fragments(nf))
    return 0;

  const int nb = nedges();

  resizable_array<int> & atoms_in_fragment = fragment_information.atoms_in_fragment();
  resizable_array<int> & bonds_in_fragment = fragment_information.bonds_in_fragment();

  for (int i = 0; i < _number_elements; i++)
  {
    if (0 == include_atom[i])
      continue;

    const int f = fragment_membership[i];
    atoms_in_fragment[f]++;
  }

  for (int i = 0; i < nb; ++i)
  {
    const Bond * b = _bond_list[i];

    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();

    const int f1 = fragment_membership[a1];
    const int f2 = fragment_membership[a2];

    if (f1 == f2)
      bonds_in_fragment[f1]++;
  }

  return nf;
}
