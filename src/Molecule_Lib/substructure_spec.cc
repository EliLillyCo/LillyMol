#include <cctype>

#include "misc.h"
#include "misc2.h"
#include "smiles.h"
#include "path.h"

#include "substructure.h"
#include "target.h"
#include "tokenise_atomic_smarts.h"

/*
  Feb 2005. Want to be able to modify the behaviour of the H
  directive in a smarts. The Daylight rule is that H means
  exactly one implicit hydrogen. I want it to optionally be
  minimum of 1 hydrogen
*/

static int _h_means_exactly_one_hydrogen = 1;

int 
h_means_exactly_one_hydrogen()
{
  return _h_means_exactly_one_hydrogen;
}

void
set_h_means_exactly_one_hydrogen (int s)
{
  _h_means_exactly_one_hydrogen = s;
}

int
Substructure_Atom_Specifier::atomic_number (atomic_number_t & z) const
{
  if (_element.number_elements())
  {
    z = _element[0]->atomic_number();
    return 1;
  }

  return 0;
}

void
Substructure_Atom_Specifier::_default_values()
{
  _aromaticity = SUBSTRUCTURE_NOT_SPECIFIED;
  _chirality   = SUBSTRUCTURE_NOT_SPECIFIED;

  _preference_value = 0;

  _attributes_specified = 0;

  _match_spinach_only = -1;    // -1 means not specified

  _all_rings_kekule = SUBSTRUCTURE_NOT_SPECIFIED;

  _symmetry_group = 0;

  return;
}

Substructure_Atom_Specifier::Substructure_Atom_Specifier()
{
  _default_values();

  return;
}

Substructure_Atom_Specifier::Substructure_Atom_Specifier (const Element * e)
{
  _default_values();

  _element.add(e);

  return;
}

Substructure_Atom_Specifier::~Substructure_Atom_Specifier()
{
  if (-8 == _attributes_specified)
    cerr << "Deleting already deleted Substructure_Atom\n";

  assert (ok());

  _attributes_specified = -8;

  return;
}

int
Substructure_Atom_Specifier::ok() const
{
  if (-8 == _attributes_specified)
    return 0;

  return 1;
}

int
Substructure_Atom_Specifier::debug_print(std::ostream & os,
                               const IWString & indentation) const
{
  assert (os.good());

  os << indentation << "Substructure Atom Specifier, " << _attributes_specified << " attributes specified\n";

  if (_preference_value)
    os << indentation << "  Preference " << _preference_value << endl;

  if (! ok())
    os << "Warning, OK fails\n";

  os << indentation << "  atomic number ";
  if (_element.number_elements())
  {
    for (int i = 0; i < _element.number_elements(); i++)
    {
      os << " " << _element[i]->atomic_number();
    }
    os << endl;
  }
  else
    os << " not specified\n";

  if (_ncon.is_set())
    os << indentation << "  ncon " << _ncon << endl;
  if (_nbonds.is_set())
    os << indentation << "  nbonds " << _nbonds << endl;
  if (_ncon2.is_set())
    os << indentation << "  ncon2 " << _ncon2 << endl;
  if (_formal_charge.is_set())
    os << indentation << "  Formal charge " << _formal_charge << endl;
  if (_nrings.is_set())
    os << indentation << "  Nrings " << _nrings << endl;
  if (_ring_bond_count.is_set())
    os << indentation << "  RingBondCount " << _ring_bond_count << endl;
  if (_hcount.is_set())
    os << indentation << "  Hcount " << _hcount << endl;
  if (_isotope.is_set())
    os << indentation << "  Isotope " << _isotope << endl;
  if (_lone_pair_count.is_set())
    os << indentation << "  Lone Pair " << _lone_pair_count << endl;
  if (_unsaturation.is_set())
    os << indentation << "  unsaturation " << _unsaturation << endl;
  if (_attached_heteroatom_count.is_set())
  os << indentation << "  Attached heteroatom count " << _attached_heteroatom_count << endl;
  if (SUBSTRUCTURE_NOT_SPECIFIED != _aromaticity)
    os << indentation << "  Aromaticity = " << _aromaticity << endl;

  if (SUBSTRUCTURE_NOT_SPECIFIED != _chirality)
    os << indentation << "  chiral " << _chirality << endl;
  if (_symmetry_degree.is_set())
    os << indentation << "  symmd " << _symmetry_degree << endl;
  if (_symmetry_group > 0)
    os << indentation << "  symmg " <<  _symmetry_group << endl;

  os << endl;

  return 1;
}

int
Substructure_Atom_Specifier::terse_details (std::ostream & os,
                               const IWString & indentation) const
{
  assert (os.good());

  os << indentation << "Details on substructure Atom Specifier\n";

  if (! ok())
    os << "Warning, OK fails\n";

  os << indentation << "Substructure atom specifications\n";
  if (_preference_value)
    os << indentation << " preference " << _preference_value << endl;
  if (_element.number_elements())
  {
    os << indentation << " atomic number";
    for (int i = 0; i < _element.number_elements(); i++)
    {
      os << ' ' << _element[i]->atomic_number();
    }
    os << endl;
  }
  if (_ncon.is_set())
    os << indentation << " ncon " << _ncon << endl;
  if (_nbonds.is_set())
    os << indentation << " nbonds " << _nbonds << endl;
  if (_ncon2.is_set())
    os << indentation << " ncon2 " << _ncon2 << endl;
  if (_formal_charge.is_set())
    os << indentation << " Formal charge " << _formal_charge << endl;
  if (_nrings.is_set())
    os << indentation << " Nrings " << _nrings << endl;
  if (_ring_bond_count.is_set())
    os << indentation << " RingBondCount " << _ring_bond_count << endl;
  if (_hcount.is_set())
    os << indentation << " Hcount " << _hcount << endl;
  if (_lone_pair_count.is_set())
    os << indentation << " Lone Pair " << _lone_pair_count << endl;
  if (SUBSTRUCTURE_NOT_SPECIFIED != _aromaticity)
    os << indentation << "Aromaticity = " << _aromaticity << endl;

  return 1;
}

int
Substructure_Atom_Specifier::involves_rings() const
{
  assert (NULL == "This is not working");

  if (_nrings.is_set())
  {
  }

  return 0;
}

/*
  A molecule has been made from the query. Each component of the query
  needs to ensure that it does not violate any aspect of the molecule.

  Note that the Molecule may have 0 == nrings for any atom, but the
  query itself may specify a positive ring value, so all we really
  can do is to check that any max_nrings specified is at least as
  large as nr
*/

int
Substructure_Atom_Specifier::_adjust_nrings (int nr)
{
//  There is nothing we can do in this case. Even if the atom in the query
//  does not appear in a ring, the query can ultimately be embedded in a ring.

  if (0 == nr && _nrings.is_set())
    return 1;
  
  int mxnr;
  if (! _nrings.max(mxnr))
    return 1;

  if (mxnr < nr)
  {
    cerr << "Substructure_Atom_Specifier::_adjust_nrings: nrings violation\n";
    cerr << "Perceived " << nr << " rings in query, max is " << mxnr << endl;
    if (! _nrings.adjust_to_accommodate(nr))
    {
      cerr << "Could not adjust nr " << _nrings << endl;
    }
  }


  return 1;
}

/*
  The query has been found to be cyclic. Each Substructure_Atom must ensure
  that its own _ring_sizes specifications are consistent with the query.
*/

int
Substructure_Atom_Specifier::_adjust_ring_sizes (const List_of_Ring_Sizes & ring_sizes_perceived)
{
  if (0 == ring_sizes_perceived.number_elements())
    return 1;

  if (! _ring_size.is_set())
  {
    for (int i = 0; i < ring_sizes_perceived.number_elements(); i++)
    {
      _ring_size.add(ring_sizes_perceived[i]);
    }

    return 1;
  }

  _ring_size.add_non_duplicated_elements(ring_sizes_perceived);

  return 1;
}

/*
  A ring has been detected in the query. Each substructure atom must check to 
  ensure its own consistency with the ring system in the query.

  There are two parts of ring specification
    _nrings and _ring_sizes
*/

int
Substructure_Atom_Specifier::adjust_ring_specifications (int nr, 
                                  const List_of_Ring_Sizes & ring_sizes)
{
  assert (nr >= 0);

  if (! _adjust_nrings(nr))
    return 0;

  if (! _adjust_ring_sizes(ring_sizes))
    return 0;

  return 1;
}

/*
  Substructure_Atom_Specifier::min_ncon is only used by ::create_molecule. It returns an
  estimate of the minimum connectivity associated with an atom.

  A floor on the return code is the actual number of connections present in the
  query.

  If a minimum has been specified, return that value.
  Otherwise, check to see if there is an explicit connectivity available.
*/

int
Substructure_Atom_Specifier::min_ncon() const
{
  int tmp;
  if (_ncon.min(tmp))       // an explicit min value is known
    return tmp;

  iwmin<int> rc(0);

  for (int i = 0; i < _ncon.number_elements(); i++)
    rc.extra(_ncon[i]);

  return rc.minval();
}


int
Substructure_Atom_Specifier::min_nbonds() const
{
  int minimum_possible = min_ncon();

  int nbonds;
  if (_nbonds.min(nbonds))
  {
    assert (nbonds >= minimum_possible);
    return nbonds;
  }

  if (_nbonds.is_set())
  {
    int nbonds = _nbonds[0];
    assert (nbonds >= minimum_possible);
    return nbonds;
  }

  return minimum_possible;
}

int
Substructure_Atom_Specifier::set_element (const Element * e)
{
  assert (ok());
  assert (e->ok());

  if (!  _element.add_if_not_already_present(e))
    return 0;

  _element_unique_id.add(e->unique_id());

  return 1;
}

int
Substructure_Atom_Specifier::set_ncon (int ncon)
{
  assert (ok());
  assert (ncon >= 0);

  return _ncon.add(ncon);
}

int
Substructure_Atom_Specifier::set_min_ncon (int ncon)
{
  assert (ok());
  assert (ncon >= 0);

  return _ncon.set_min(ncon);
}

int
Substructure_Atom_Specifier::set_max_ncon (int ncon)
{
  assert (ok());
  assert (ncon >= 0);

  return _ncon.set_max(ncon);
}

int
Substructure_Atom_Specifier::set_min_nbonds (int nbonds)
{
  assert (ok());
  assert (nbonds >= 0);

  return _nbonds.set_min(nbonds);
}

int
Substructure_Atom_Specifier::set_max_nbonds (int nbonds)
{
  assert (ok());
  assert (nbonds >= 0);

  return _nbonds.set_max(nbonds);
}

int
Substructure_Atom_Specifier::set_nbonds (int nbonds)
{
  assert (ok());
  assert (nbonds >= 0);

  return _nbonds.add(nbonds);
}

int
Substructure_Atom_Specifier::set_ncon2 (int ncon2)
{
  assert (ok());
  assert (ncon2 >= 0);

  return _ncon2.add(ncon2);
}

int
Substructure_Atom_Specifier::set_nrings (int nrings)
{
  assert (ok());
  assert (nrings >= 0);

  return _nrings.add(nrings);
}

int
Substructure_Atom_Specifier::set_min_nrings (int nr)
{
  assert (ok());
  assert (nr > 0);

  return _nrings.set_min(nr);
}

int
Substructure_Atom_Specifier::set_formal_charge (formal_charge_t fc)
{
  assert (ok());

  return _formal_charge.add(fc);
}

int
Substructure_Atom_Specifier::set_attached_heteroatom_count (int ahc)
{
  assert (ok());

  assert (ahc >= 0);

  return _attached_heteroatom_count.add(ahc);
}

int
Substructure_Atom_Specifier::set_chirality (int c)
{
  assert (ok());

  assert (c >= 0 && c <= 2);

  _chirality = c;

  return 1;
}

int
Substructure_Atom_Specifier::formal_charge (formal_charge_t & fc) const
{
  assert (ok());

  if (! _formal_charge.is_set())
    return 0;

  if (_formal_charge.number_elements())
  {
    fc = _formal_charge[0];
    return 1;
  }

//  Should probably do a check for min and max. Have never used them, so...

  return 0;
}

/*
  The difference between set_aromaticity and update_aromaticity is
  that set always does an overwrite.
*/

int
Substructure_Atom_Specifier::set_aromaticity (aromaticity_type_t arom)
{
  assert (OK_ATOM_AROMATICITY(arom));

  _aromaticity = arom;

  return 1;
}

int
Substructure_Atom_Specifier::update_aromaticity (aromaticity_type_t arom)
{
  assert (OK_ATOM_AROMATICITY(arom));

  if (IS_AROMATIC_ATOM(arom))
  {
    if (SUBSTRUCTURE_NOT_SPECIFIED == _aromaticity)
      _aromaticity = AROMATIC;
    else
      SET_AROMATIC_ATOM(_aromaticity);
  }

  if (IS_ALIPHATIC_ATOM(arom))
  {
    if (SUBSTRUCTURE_NOT_SPECIFIED == _aromaticity)
      _aromaticity = NOT_AROMATIC;
    else
      SET_ALIPHATIC_ATOM(_aromaticity);
  }

  return 1;
}

//#define DEBUG_RING_SIZE_MATCHES

/*
  All the ring size matching functions look alike, so they are done in
  this function.
*/

static int
match_ring_sizes(const Min_Max_Specifier<int> & ring_sizes_in_query,
                 const List_of_Ring_Sizes * ring_sizes_in_molecule,
                 const char * comment)
{
  (void) comment;     // suppress compiler warnings of unused arg

#ifdef DEBUG_RING_SIZE_MATCHES
  cerr << "Checking '" << comment << "' ring sizes. Atom is in rings of size:\n";
  for (int i = 0; i < ring_sizes_in_molecule->number_elements(); i++)
    cerr << "   " << ring_sizes_in_molecule->item(i) <<
         ", match = " << ring_sizes_in_query.matches(ring_sizes_in_molecule->item(i)) << endl;
#endif

  for (int i = 0; i < ring_sizes_in_molecule->number_elements(); i++)
  {
    int nr = ring_sizes_in_molecule->item(i);
    if (ring_sizes_in_query.matches(nr))
      return 1;
  }

  return 0;    // no match found
}

int
Substructure_Atom_Specifier::_match_scaffold_bonds_attached_to_ring(Target_Atom & target)
{
  Molecule & m = *(target.m());

  m.ring_membership();    // force just in case

  int nr = m.nrings();

  for (auto i = 0; i < nr; ++i)
  {
    const Ring * ri = m.ringi(i);

    if (! ri->contains(target.atom_number()))
      continue;

    if (_match_scaffold_bonds_attached_to_ring(target, *ri))
      return 1;
  }

  return 0;
}

int
Substructure_Atom_Specifier::_match_symmetry_degree(Target_Atom & target_atom) const
{
  Molecule & m = *(target_atom.m());

  const int matoms = m.natoms();

  const auto a = target_atom.atom_number();

  const int * s = m.symmetry_classes();

  const int sa = s[a];

  int symmetric_atoms = 0;
  for (int i = 0; i < matoms; ++i)
  {
    if (sa == s[i])
      symmetric_atoms++;
  }

  return _symmetry_degree.matches(symmetric_atoms);
}

int
Substructure_Atom_Specifier::_match_scaffold_bonds_attached_to_ring(Target_Atom & target_atom,
                                                                    const Ring & r) const
{
  Molecule & m = *(target_atom.m());

  int ring_size = r.number_elements();

  int exocyclic_scaffold_bonds = 0;

  for (auto i = 0; i < ring_size; ++i)
  {
    atom_number_t j = r[i];

    const Atom * aj = m.atomi(j);

    int jcon = aj->ncon();

    if (2 == jcon)
      continue;

    for (int k = 0; k < jcon; ++k)
    {
      const Bond * b = aj->item(k);

      atom_number_t l = b->other(j);

      if (r.contains(l))
        continue;

      if (b->nrings())    // fused to another ring. By convention this is an exocyclic scaffold bond
      {
        exocyclic_scaffold_bonds++;
        continue;
      }

      Target_Atom & al = target_atom.atom(l);

      if (al.is_spinach())
        continue;

      exocyclic_scaffold_bonds++;
    }
  }

//cerr << "Ring " << r << " has " << exocyclic_scaffold_bonds <<  " exocyclic_scaffold_bonds, match? " << _scaffold_bonds_attached_to_ring.matches(exocyclic_scaffold_bonds) << endl;

  return _scaffold_bonds_attached_to_ring.matches(exocyclic_scaffold_bonds);
}

/*
  In order to match a query atom, we first check any properties of the Substructure_Atom
  itself, followed by the individual Substructure_Atom_Specifiers, and then check
  any rejections
*/

//#define DEBUG_ATOM_MATCHES

int
Substructure_Atom_Specifier::_matches (Target_Atom & target)
{
// Check all attributes

  int attributes_checked = 0;

#ifdef DEBUG_ATOM_MATCHES
  cerr << "Need to check " << _attributes_specified << " attributes\n";
  if (_element.number_elements())
  {
    cerr << "Check element, target is " << target.atomic_number() <<
            " check " << _element.number_elements() << " values, match = " << (_element.contains(target.element())) << endl;
  }
#endif

#ifdef USING_ELEMENT
   if (_element.number_elements())
   {
     if ( ! _element.contains(target.element()))
       return 0;

     attributes_checked++;
     if (attributes_checked == _attributes_specified)
       return 1;
   }
#endif

#ifdef DEBUG_ATOM_MATCHES
   if (_element_unique_id.number_elements())
     cerr << "check ASHV, target " << target.element_unique_id() << " check " << _element_unique_id.number_elements() << " match " << _element_unique_id.contains(target.element_unique_id()) << endl;
   for (int i = 0; i < _element_unique_id.number_elements(); ++i)
   {
     cerr << ' ' << _element_unique_id[i] << endl;
   }
#endif

   if (_element_unique_id.number_elements())
   {
     if ( ! _element_unique_id.contains(target.element_unique_id()))
       return 0;

     attributes_checked++;
     if (attributes_checked == _attributes_specified)
       return 1;
   }

#ifdef DEBUG_ATOM_MATCHES
   if (_ncon.is_set())
     cerr << "Check ncon, target is " << target.ncon() << " match = " << _ncon.matches(target.ncon()) << endl;
#endif

   if (_ncon.is_set())
   {
     if (! _ncon.matches(target.ncon()))
       return 0;

     attributes_checked++;
     if (attributes_checked == _attributes_specified)
       return 1;
   }


#ifdef DEBUG_ATOM_MATCHES
   if (_nbonds.is_set())
     cerr << "Check nbonds, target is " << target.nbonds() << " match = " << _nbonds.matches(target.nbonds()) << endl;
#endif

   if (_nbonds.is_set())
   {
     if (! _nbonds.matches(target.nbonds()))
       return 0;

     attributes_checked++;
     if (attributes_checked == _attributes_specified)
       return 1;
   }


#ifdef DEBUG_ATOM_MATCHES
   if (_formal_charge.is_set())
     cerr << "Check formal charges, target = " << target.formal_charge() << endl;
#endif

  if (_formal_charge.is_set())
  {
    if (! _formal_charge.matches(target.formal_charge()))
       return 0;

     attributes_checked++;
     if (attributes_checked == _attributes_specified)
       return 1;
  }


#ifdef DEBUG_ATOM_MATCHES
  if (_nrings.is_set())
  {
    cerr << "Check nrings. ";
    if (1 == _nrings.number_elements() && 0 == _nrings[0])
      cerr << "Specified non ring, target is " << target.is_ring_atom() << endl;
    else 
      cerr << "Nrings for target is " << target.nrings() << ", match? " << _nrings.matches(target.nrings()) << endl;
  }
#endif

  if (_nrings.is_set())
  {
    if (1 == _nrings.number_elements() && 0 == _nrings[0])
    {
      if (target.is_ring_atom())
        return 0;
    }
    else if (! _nrings.matches(target.nrings()))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }

#ifdef DEBUG_ATOM_MATCHES
  if (_ring_bond_count.is_set())
    cerr << "Check _ring_bond_count, target " << target.ring_bond_count() << " matches " << _ring_bond_count.matches(target.ring_bond_count()) << endl;
#endif
  if (_ring_bond_count.is_set())
  {
    if (1 == _ring_bond_count.number_elements() && 0 == _ring_bond_count[0])
    {
      if (target.is_ring_atom())
        return 0;
    }
    else if (! _ring_bond_count.matches(target.ring_bond_count()))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }


#ifdef DEBUG_ATOM_MATCHES
  if (_hcount.is_set())
    cerr << "Check hcount: target " << target.hcount() << " matches " << _hcount.matches(target.hcount()) << endl;
#endif

  if (_hcount.is_set())
  {
    if (! _hcount.matches(target.hcount()))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }

#ifdef DEBUG_ATOM_MATCHES
  if (SUBSTRUCTURE_NOT_SPECIFIED != _aromaticity)
     cerr << "Let's try aromaticity " << _aromaticity << " vs " << target.aromaticity() << endl;
#endif

  if (SUBSTRUCTURE_NOT_SPECIFIED != _aromaticity)
  {
    aromaticity_type_t tmp = target.aromaticity();

    if (tmp == _aromaticity)
      ;
    else if (IS_AROMATIC_ATOM(_aromaticity) && IS_AROMATIC_ATOM(tmp))
      ;
    else if (IS_ALIPHATIC_ATOM(_aromaticity) && IS_ALIPHATIC_ATOM(tmp))
      ;
    else
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }


#ifdef DEBUG_ATOM_MATCHES
  if (_attached_heteroatom_count.is_set())
    cerr << "Try attached heteroatom count: target =" << target.attached_heteroatom_count() << 
            " match is " << _attached_heteroatom_count.matches(target.attached_heteroatom_count()) << endl;
#endif

  if (_attached_heteroatom_count.is_set())
  {
    int tahc = target.attached_heteroatom_count();
    if (! _attached_heteroatom_count.matches(tahc))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }


#ifdef DEBUG_ATOM_MATCHES
  if (_unsaturation.is_set())
  {
    int tmp = target.nbonds() - target.ncon();
    cerr << "Try unsaturation: target =" << tmp << 
            " match is " << _unsaturation.matches(tmp) << endl;
  }
#endif

  if (_unsaturation.is_set())
  {
    int tmp = target.nbonds() - target.ncon();
    if (! _unsaturation.matches(tmp))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }


  if (_vinyl.is_set())
  {
    if (! _vinyl.matches(target.vinyl()))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }


  if (_aryl.is_set())
  {
    if (! _aryl.matches(target.aryl()))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }


  if (_lone_pair_count.is_set())
  {
    int lpc = target.lone_pair_count();
    if (! _lone_pair_count.matches(lpc))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }


  if (_heteroatoms_in_ring.is_set())
  {
    int hir = target.heteroatoms_in_ring();
    if (! _heteroatoms_in_ring.matches(hir))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }

#ifdef DEBUG_ATOM_MATCHES
  if (_ring_size.is_set())
    cerr << "Test ring sizes SET rs = " << _ring_size.is_set() << endl;
#endif       

  if (_ring_size.is_set())
  {
    const List_of_Ring_Sizes * ring_sizes_in_molecule = target.sssr_ring_sizes();
    if (! match_ring_sizes(_ring_size, ring_sizes_in_molecule, "ring_size"))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }

#ifdef DEBUG_ATOM_MATCHES
  if (_match_spinach_only >= 0)
  {
    cerr << "Checking _match_spinach_only " << _match_spinach_only << " target.is_spinach() " << target.is_spinach() << endl;
  }
#endif

  if (_match_spinach_only < 0)    // not specified
    ;
  else if (_match_spinach_only && ! target.is_spinach())
    return 0;
  else if (0 == _match_spinach_only && target.is_spinach())
    return 0;
  else
  {
    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }

  if (_fused_system_size.is_set())
  {
    int tmp = target.fused_system_size();
    if (! _fused_system_size.matches(tmp))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }


  if (_daylight_x.is_set())
  {
    int tmp = target.daylight_x();
    if (! _daylight_x.matches(tmp))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }

#ifdef DEBUG_ATOM_MATCHES
  if (_aromatic_ring_sizes.is_set())
    cerr << " aromatic ring sizes = " << _aromatic_ring_sizes.is_set();
#endif

  if (_aromatic_ring_sizes.is_set())
  {
    const List_of_Ring_Sizes * ring_sizes_in_molecule = target.aromatic_ring_sizes();
    if (! match_ring_sizes(_aromatic_ring_sizes, ring_sizes_in_molecule, "aromatic ring size"))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }

#ifdef DEBUG_ATOM_MATCHES
  if (_aliphatic_ring_sizes.is_set())
    cerr << "Aliphatic ring sizes " << _aliphatic_ring_sizes.is_set() << endl;
#endif

  if (_aliphatic_ring_sizes.is_set())
  {
    const List_of_Ring_Sizes * ring_sizes_in_molecule = target.aliphatic_ring_sizes();
    if (! match_ring_sizes(_aliphatic_ring_sizes, ring_sizes_in_molecule, "aliphatic ring size"))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }

#ifdef DEBUG_ATOM_MATCHES
  if (_ncon2.is_set())
    cerr << "Check hcount, target is " << target.nrings() << endl;
#endif
  if (_ncon2.is_set())
  {
    if (! _ncon2.matches(target.ncon2()))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }

  if (_isotope.is_set())
  {
    if (! _isotope.matches(target.isotope()))
      return 0;

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }

// We can't really check chirality now because the adjoining atoms are
// not matched yet. If our query says chiral of some kind and the
// matched atom has no chirality, that's a non-match

#ifdef DEBUG_ATOM_MATCHES
  if (SUBSTRUCTURE_NOT_SPECIFIED != _chirality)
    cerr << "Checking chirality " << _chirality << " vs " << target.chiral_centre() << endl;
#endif

  if (SUBSTRUCTURE_NOT_SPECIFIED != _chirality)
  {
    if (NULL == target.chiral_centre())    // no chiral centre on matched atom
    {
      if (_chirality)
        return 0;
    }
    else     // matched atom has chiral centre
    {
      if (0 == _chirality)
        return 0;
    }

    attributes_checked++;
    if (attributes_checked == _attributes_specified)
      return 1;
  }

  if (_scaffold_bonds_attached_to_ring.is_set())
  {
    if (! target.is_ring_atom())    // cannot match
      return 0;

    if (! _match_scaffold_bonds_attached_to_ring(target))
      return 0;

    attributes_checked++;

    if (attributes_checked == _attributes_specified)
      return 1;
  }

  if (_symmetry_degree.is_set())
  {
    if (! _match_symmetry_degree(target))
      return 0;

    attributes_checked++;

    if (attributes_checked == _attributes_specified)
      return 1;
  }

  if (SUBSTRUCTURE_NOT_SPECIFIED != _all_rings_kekule)
  {
    attributes_checked++;

    const auto k = target.all_rings_kekule();
    if (0 == _all_rings_kekule && 0 == k)
      ;
    else if (_all_rings_kekule > 0 && k > 0)
      ;
    else
      return 0;

    if (attributes_checked == _attributes_specified)
      return 1;
  }

  if (_attributes_specified != attributes_checked)
  {
    cerr << "Oops, attributes specified = " << _attributes_specified << " processed = " << attributes_checked << endl;
    debug_print(cerr);
    assert (NULL == "This should not happen");
  }

#ifdef DEBUG_ATOM_MATCHES
  cerr << "All " << _attributes_specified << " attributes checked, returning match\n";
#endif

  return 1;
}


/*
  Public interface for _matches
*/

int
Substructure_Atom_Specifier::matches (Target_Atom & target)
{
  int rc = _matches(target);

#ifdef DEBUG_ATOM_MATCHES
  cerr << "Substructure_Atom_Specifier::matches: to atom " << target.atom_number() << " match " << rc << endl;
#endif

  return rc;
}

/*
  This function was written for the fingerprint routines.
  They need to know whether or not this query atom is specified 
  as being in a ring or not
*/

int
Substructure_Atom_Specifier::determine_ring_or_non_ring (int & result) const
{
  if (! _nrings.is_set() && ! _ring_bond_count.is_set())
    return 0;

  if (1 == _nrings.number_elements())
  {
    result = _nrings[0];
    return 1;
  }

  if (1 == _ring_bond_count.number_elements())
  {
    result = _ring_bond_count[0];
    return 0;
  }

// Come back to this sometime and finish the logic for _ring_bond_count

  if (! _nrings.is_set())
    return 0;

  int tmp;
  if (_nrings.min(tmp) && tmp > 0)
  {
    result = tmp;
    return 1;
  }

  if (0 == _nrings.number_elements())    // must just be a max value specified
    return 0;

// There are multiple values for nrings. If they are all > 0, then ok.
// Fingerprint just needs to know ring or non ring.

  for (int i = 0; i < _nrings.number_elements(); i++)
  {
    if (0 == _nrings[i])
      return 0;
  }

// All _nrings values were > 0

  result = _nrings[0];

  return 1;
}

/*
  This function was written for the fingerprint routines.
  Does an explicit ncon value exist for this atom.
*/

int
Substructure_Atom_Specifier::ncon_specification (int & result) const
{
  if (! _ncon.is_set())
    return 0;

  if (_ncon.number_elements() > 1)
    return 0;

  if (1 == _ncon.number_elements())
  {
    result = _ncon[0];
    return 1;
  }

  return 0;
}

/*
  When making a molecule from a query, we need to know whether or not there
  is a unique hcount specification for an atom
*/

int
Substructure_Atom_Specifier::hcount_specification (int & result) const
{
  if (! _hcount.is_set())
    return 0;

  if (_hcount.number_elements() > 1)
    return 0;

  if (1 == _hcount.number_elements())
  {
    result = _hcount[0];
    return 1;
  }

  return 0;
}

int
Substructure_Atom_Specifier::formal_charge_specification (formal_charge_t & result) const
{
  if (! _formal_charge.is_set())
    return 0;

  if (_formal_charge.number_elements() > 1)
    return 0;

  if (1 == _formal_charge.number_elements())
  {
    result = _formal_charge[0];
    return 1;
  }

  return 0;
}

/*
  A Substructure_Query object can do a consistency check upon itself.
  All components must support a check_internal_consistency function.
*/

int
Substructure_Atom_Specifier::check_internal_consistency (int connections) const
{
  (void) connections;
  return 1;
}

int
Substructure_Atom_Specifier::ring_sizes_specified (resizable_array<int> & ring_sizes) const
{
  ring_sizes.add_non_duplicated_elements(_aromatic_ring_sizes);

  ring_sizes.add_non_duplicated_elements(_aliphatic_ring_sizes);

  return ring_sizes.number_elements();
}

int
Substructure_Atom_Specifier::_set_implicit_hydrogens (Molecule & m,
                                            atom_number_t a) const
{
  int ih;
  if (! hcount_specification(ih))
  {
//  cerr << "Setting implicit hcount for atom " << i << " to zero\n";
    m.set_implicit_hydrogens(a, 0);
  }
  else
    m.set_implicit_hydrogens(a, ih);

  return 1;
}

int
Substructure_Atom_Specifier::_fill_min_ncon (Molecule & m,
                                   atom_number_t a) const
{
  int connection_shortage = min_ncon() - m.ncon(a);
  if (0 == connection_shortage)
    return 1;

  assert (connection_shortage > 0);

  int bond_shortage = min_nbonds() - m.nbonds(a);
  assert (bond_shortage >= 0);

// cerr << "Atom " << a << " cs = " << connection_shortage << " bs = " << bond_shortage << endl;
  for (int i = 0; i < connection_shortage; i++)
  {
    Atom * h = new Atom(0);
    m.add(h);
    bond_type_t bt = SINGLE_BOND;
    if (bond_shortage > connection_shortage)
    {
      bt = DOUBLE_BOND;
      bond_shortage--;
    }
    m.add_bond(a, m.natoms() - 1, bt);
  }

  return 1;
}

Atom *
Substructure_Atom_Specifier::create_atom() const
{
  if (_element.number_elements())
    return new Atom(_element[0]);

  return new Atom(0);
}

int
Substructure_Atom::_parse_smarts_environment (const Atomic_Smarts_Component & env)
{
  if (! env.starts_with("$(") || ! env.ends_with(")"))
  {
    cerr << "Substructure_Atom::_parse_smarts_environment: environments must start with '$(' and end with ')'\n";
    return 0;
  }

  return _environment.create_from_smarts(env);
}

/*
  Does the smarts end with : followed by digits. The caller must have checked to make
  sure there is a colon in here.
  We trim MYSMARTS to remove what we consume
*/

int
Substructure_Atom::_extract_initial_atom_number (const_IWSubstring & mysmarts)
{
  if (mysmarts.ends_with(':'))
    return 0;

  int lastcolon = mysmarts.rindex(':');
  assert (lastcolon > 0);

  int znumber = 0;
  for (int i = lastcolon + 1; i < mysmarts.length(); i++)
  {
    char c = mysmarts[i];

    int d = c - '0';
    if (d < 0 || d > 9)
      return 0;

    znumber = 10 * znumber + d;
  }

// All characters between the last ':' and the end were digits.

  mysmarts.iwtruncate(lastcolon);

//_initial_atom_number = _unique_id = _atom_map_number = znumber;  NO! wrong

  _unique_id = _atom_map_number = znumber;

//cerr << "Substructure_Atom::_extract_initial_atom_number: set " << _initial_atom_number << endl;

  return 1;
}

/*
  Someone has found a colon. We are to extract the number following the colon.
  We must find a closing square bracket.
*/

/*static int
fetch_atom_number (const_IWSubstring & smarts,
                   int istart, 
                   int & unique_id)
{
  assert (':' == smarts[istart]);

  istart++;     // skip over the ':'

  unique_id = 0;

  int rc = 0;      // the number of characters beyond the colon we process

  while (istart < smarts.length())
  {
    if (']' == smarts[istart])
    {
      if (0 == rc)
        cerr << "fetch_atom_number: must specify atom number after colon\n";

      return rc;
    }

    int tmp = smarts[istart] - '0';

    if (tmp >= 0 && tmp <= 9)
      unique_id = 10 * unique_id + tmp;
    else
    {
      cerr << "Invalid atom number specifier '" << smarts << "'\n";
      return 0;
    }

    istart++;
    rc++;
  }

  return 1;
}*/

/*
  We have encountered a '$(' in an atomic smarts. Examine the string to
  find the closing ).
*/

static int
fetch_environment (const_IWSubstring & env)
{
  assert (env.starts_with ("$("));

  int paren_level = 1;

  for (int i = 2; i < env.nchars(); i++)
  {
    if ('(' == env[i])
      paren_level++;
    else if (')' == env[i])
    {
      paren_level--;
      if (0 == paren_level)
      {
        env.iwtruncate(i + 1);
        return i + 1;
      }
    }
  }

  return 0;    // yipes, no closing paren found!
}

/*static int
fetch_environment (const const_IWSubstring & env, int & characters_processed,
                   IWString & result)
{
  assert (0 == result.nchars());
  assert ('$' == env[characters_processed]);

  int chars_remaining = env.nchars() - characters_processed;

  if (chars_remaining < 4)
  {
    cerr << "fetch_environment: the environment must have at least $(.)\n";
    return 0;
  }

  result.resize (chars_remaining);

  if ('(' != env[characters_processed + 1])
  {
    cerr << "fetch_environment: environment must start with '$('\n";
    return 0;
  }

  int paren_level = 1;

  result = "$(";

  for (int i = characters_processed + 2; i < env.nchars(); i++)
  {
    result += env[i];
//  cerr << "Env is now '" << result << "'\n";

    if ('(' == env[i])
      paren_level++;
    else if (')' == env[i])
    {
      paren_level--;
      if (0 == paren_level)
      {
        if (3 == result.nchars())
        {
          cerr << "fetch_environment: empty environment\n";
          return 0;
        }

        characters_processed += result.nchars();
        return 1;
      }
    }
  }

  cerr << "fetch_environment: mismatched parentheses\n";
  return 0;
}*/

/*
  Jul 99. Change the default value to 1. Why did I ever do it differently?
*/

static int respect_aliphatic_smarts = 1;

void
set_respect_aliphatic_smarts (int s)
{
  respect_aliphatic_smarts = s;
}

static void
truncate_after_digits (const const_IWSubstring & ignore_these,
                       const_IWSubstring & s)
{
  int nchars = s.length();

  for (int i = 0; i < nchars; i++)
  {
    char c = s[i];
    if (isdigit(c))
      continue;

    if (ignore_these.contains(c))
      continue;

    s.iwtruncate(i);
    return;
  }

  return;     // looked good all the way out
}

//#define DEBUG_ATOM_CONSTRUCT_FROM_SMARTS_TOKEN

/*
  this turned out to be surprisingly difficult, and revealed some interesting
  aspects of how I'm doing substructures. I quickly realised that my data
  structure cannot represent something like 
    c,n;R,H
  which (I think) is grouped as
    (c,n);(R,H)
  which means "a carbon or nitrogen" AND "either a ring atom or an attached H"

  Oct 97, well wait a minute, of course I can store this. Store "carbon or nitrogen"
  in the object's Substructure_Atom_Specifier, and then two components - remember
  components are matched with the OR operator.

  Therefore I must place a limitation, only one "and group" may have a
  comma operator, and to make things even easier for myself, this must
  be the first grouping

  Dec 2005. Major problems. Cannot correctly parse [O,S;H] or the
  other way round [H;O,S]. These must be interpreted as "Oxygen or Sulphur,
  and having one implicit Hydrogen"

  But, we tokenise the smarts, and we present 'H' by itself to ->construct_from_smarts_token()
  and it will interpret the H as elemental Hydrogen. We need to pre-screen
  the tokens and see if any of them match elemental hydrogen

  Jan 2006: This is still broken. Not sure how to fix it. Get back to this...
*/

static IW_Regular_Expression elemental_hydrogen_rx("^(H|H\\+|H-|[0-9]+H)$");

#define SMARTS_PREVIOUS_TOKEN_UNSPECIFIED            0

#define SMARTS_PREVIOUS_TOKEN_OPERATOR_AND_SEMICOLON 1
#define SMARTS_PREVIOUS_TOKEN_OPERATOR_AND_AMPER     2
#define SMARTS_PREVIOUS_TOKEN_OPERATOR_NOT           4    
#define SMARTS_PREVIOUS_TOKEN_OPERATOR_OR            8

#define SMARTS_PREVIOUS_TOKEN_MASS                  16
#define SMARTS_PREVIOUS_TOKEN_ENVIRONMENT           32
#define SMARTS_PREVIOUS_TOKEN_ELEMENT               64
#define SMARTS_PREVIOUS_TOKEN_V                     128
#define SMARTS_PREVIOUS_TOKEN_CHARGE                256
#define SMARTS_PREVIOUS_TOKEN_CHIRALITY             512
#define SMARTS_PREVIOUS_TOKEN_RING                 1024
#define SMARTS_PREVIOUS_TOKEN_X                    2048
#define SMARTS_PREVIOUS_TOKEN_T                    4096
#define SMARTS_PREVIOUS_TOKEN_RBC                  8192

int
Substructure_Atom::construct_from_smarts_token(const const_IWSubstring & smarts)
{
#ifdef DEBUG_ATOM_CONSTRUCT_FROM_SMARTS_TOKEN
  cerr << "Atom parsing smarts '" << smarts << "' nchars = " << smarts.length() << "\n";
#endif

  if ('[' != smarts[0])
    return construct_from_smiles_token(smarts);

  int characters_to_process = smarts.length();

// We have a more complex specification.

  const char * initial_smarts_ptr = smarts.rawchars();

// The most important decision is whether or not the specification can
// all fit in the parent, or do we need to create components and logical
// operators.

  int right_square_bracket = -1;
  int nsemicolons = 0;
  int ncommas = 0;
  int ncarets = 0;
  int nnot = 0;
  int environment_present = 0;
  int environment_level = 0;
  int paren_level = 0;      // parentheses must balance
  int square_bracket_level = 1;     // the square brackets must balance
  int ncolon = 0;
  int curly_brace_level = 0;

  for (int i = 1; i < characters_to_process; i++)     // smarts[0] is the opening square bracket
  {
    if ('[' == smarts[i])
    {
      square_bracket_level++;
      continue;
    }

    if ('{' == smarts[i])
    {
      curly_brace_level++;
      continue;
    }

    if ('}' == smarts[i])
    {
      curly_brace_level--;
      continue;
    }

    if (']' == smarts[i])
    {
      square_bracket_level--;
      if (0 == paren_level && 0 == square_bracket_level)
      {
        right_square_bracket = i;
        break;
      }
    }
    else if (';' == smarts[i])
      nsemicolons++;
    else if (',' == smarts[i])
      ncommas++;
    else if ('^' == smarts[i])
      ncarets++;
    else if ('!' == smarts[i])
      nnot++;
    else if ('$' == smarts[i] && i < characters_to_process - 1 && '(' == smarts[i + 1])
    {
      if (0 == environment_level)
        environment_present++;
      environment_level++;
    }
    else if ('(' == smarts[i])
      paren_level++;
    else if (')' == smarts[i])
      paren_level--;
    else if (':' == smarts[i])
      ncolon++;

#ifdef DEBUG_ATOM_CONSTRUCT_FROM_SMARTS_TOKEN
    cerr << "After '" << smarts[i] << "' paren level " << paren_level << endl;
#endif
  }

  assert (0 == paren_level);

  if (right_square_bracket < 0)
  {
    smiles_error_message(initial_smarts_ptr, characters_to_process, 0,
                  "Unterminated bracket specifier");
    return 0;
  }

// I guess we could interpret '[]' as match anything....

  if (1 == right_square_bracket)
  {
    smiles_error_message(initial_smarts_ptr, characters_to_process, 0, "Empty bracket specifier");
    return 0;
  }

// Make our own local copy

  const_IWSubstring mysmarts = smarts;

  mysmarts.remove_leading_chars(1);    // get rid of the leading '['

  mysmarts.iwtruncate(right_square_bracket - 1);

#ifdef DEBUG_ATOM_CONSTRUCT_FROM_SMARTS_TOKEN
  cerr << "After getting rid of rsb '" << mysmarts << "'\n";
#endif

  characters_to_process = right_square_bracket + 1;

// If the smarts ends with ':nn', that is the atom number. We need to be careful because
// colon's can be part of an environment [xxx$(a:a)]

  if (ncolon)
  {
    _extract_initial_atom_number(mysmarts);
    ncolon--;
  }

// Any further IW customisations to atomic smarts

  int nchars = mysmarts.length() - 4;    // we want to look for /IW.

  for (int i = 0; i < nchars; ++i)
  {
    if (! mysmarts.matches_at_position(i, "/IW"))
      continue;

    const_IWSubstring c(mysmarts.rawchars() + i + 3, mysmarts.length() - i - 3);

    if ('x' == c[0])
    {
      _include_in_embedding = 0;
    }
    else if (c.starts_with("fsid"))
    {
      c.remove_leading_chars(4);
      if (! isdigit(c[0]))    // only single digit ring ids are allowed in smarts
      {
        cerr << "Substructure_Atom::construct_from_smarts_token:invalid fsid qualifier '" << c << "'\n";
        return 0;
      }
      _fused_system_id = c[0] - '0';
    }
    else if (c.starts_with ("rid"))
    {
      c.remove_leading_chars(3);
      if (! isdigit(c[0]))    // only single digit ring ids are allowed in smarts
      {
        cerr << "Substructure_Atom::construct_from_smarts_token:invalid rid qualifier '" << c << "'\n";
        return 0;
      }
      _ring_id = c[0] - '0';
    }
    else if (c.starts_with("fss"))
      ;
    else if (c.starts_with("Vy"))
      ;
    else if (c.starts_with("Ar"))
      ;
    else if (c.starts_with("spch"))
      ;
    else if (c.starts_with("hr") || c.starts_with("rh"))
      ;
    else if (c.starts_with("rscb"))
      ;
    else if (c.starts_with("symd"))
      ;
    else if (c.starts_with("symg"))
      ;
    else if (c.starts_with("Kl"))
      ;
    else if (c.starts_with("Nv"))
    {
      c.remove_leading_chars(2);
      int nv;
      if (! isdigit(c[0]))
      {
        cerr << "Substructure_Atom::construct_from_smarts_token:invalid Numeric Value specifier '" << c << "'\n";
        return 0;
      }
      fetch_numeric(c, nv, c.length());
      _numeric_value = static_cast<double>(nv);
    }
    else
    {
      cerr << "Substructure_Atom::construct_from_smiles_token: unrecognised /IW qualifier '" << c << "'\n";
      return 0;
    }
  }

  if (0 == mysmarts.length())
    return characters_to_process;

// If only the ';' operator is present, we can process in place.
// Apr 2000. WRONG! The smarts '[r4;r5]' fails - both ring sizes get put
// into the _ring_size archive and it becomes an OR condition.

// We have one or more operators.

// Scan the smarts looking for tokens followed by operators

#ifdef DEBUG_ATOM_CONSTRUCT_FROM_SMARTS_TOKEN
  cerr << "Processing " << characters_to_process << " characters\n";
#endif

  Atomic_Smarts_Component tokens;

  if (! tokens.parse(mysmarts))
  {
    cerr << "Cannot parse smarts '";
    for (int i = 0; i < characters_to_process; i++)
    {
      cerr << smarts[i];
    }
    cerr << "'\n";
    return 0;
  }

#ifdef DEBUG_ATOM_CONSTRUCT_FROM_SMARTS_TOKEN
  cerr << "After parsing tokens\n";
  tokens.debug_print(cerr);
#endif

// We need an index for placing unary operators

  int uopindex = 0;

  Atomic_Smarts_Component * asc = & tokens;

  do
  {
#ifdef DEBUG_ATOM_CONSTRUCT_FROM_SMARTS_TOKEN
    cerr << "Building component from '" << *asc << "'\n";
#endif

    if (asc->starts_with("$("))
    {
      if (! _parse_smarts_environment(*asc))
      {
        cerr << "Cannot parse environment token '" << (*asc) << "'\n";
        return 0;
      }
    }
    else
    {
      Substructure_Atom_Specifier * a = new Substructure_Atom_Specifier;
      if (! a->construct_from_smarts_token(*asc))
      {
        cerr << "Substructure_Atom::construct_from_smarts_token: cannot parse component '" << *(asc) << "'\n";
        delete a;
        return 0;
      }

      _components.add(a);

      if (IW_LOGEXP_UNDEFINED != asc->op())
        _operator.add_operator(asc->op());

      _operator.set_unary_operator(uopindex, asc->unary_operator());
      uopindex++;

#ifdef DEBUG_ATOM_CONSTRUCT_FROM_SMARTS_TOKEN
      cerr << "Unary operator is " << asc->unary_operator() << ", number results = " << _operator.number_results() << endl;
      _operator.debug_print(cerr);
#endif
    }

  } while (NULL != (asc = asc->next()));

#ifdef DEBUG_ATOM_CONSTRUCT_FROM_SMARTS_TOKEN
  cerr << "After building, operator is\n";
  _operator.debug_print(cerr);
#endif

  assert (ok());

  return right_square_bracket + 1;
}

/*
  The '@' character is capital 2 on the keyboard
*/

static int
snarf_capital_two (const char * smarts, int nchars)
{
  int rc = 1;
  nchars--;
  smarts++;
  while (nchars > 0)
  {
    if ('@' != *smarts)
      return rc;
    
    rc++;
    smarts++;
    nchars--;
  }

  return rc;
}

int
Substructure_Atom::construct_from_smarts_token (const char * smarts, int nchars)
{
  const_IWSubstring tmp;
  tmp.set(smarts, nchars);

  return construct_from_smarts_token(tmp);
}

/*
  This is a variant of the fetch_numeric in misc2.cc
  The only difference is that we allow a leading '>' or '<' before the digit

  We return the number of characters processed

  When gcc fixes namespaces, change this to get the leading '>' or '<' and
  then call the function from misc2

  Need to be careful here too, parsing just the '>' or '<' is an error
*/

int
smarts_fetch_numeric (const char * string, int & value, 
                      int & qualifier)
{
  int rc;

  if ('<' == *string)
  {
    qualifier = -1;
    string++;
    rc = 1;
  }
  else if ('>' == *string)
  {
    qualifier = 1;
    string++;
    rc = 1;
  }
  else
  {
    qualifier = 0;
    rc = 0;
  }

  value = 0;
  while (isdigit(*string))
  {
    int tmp = *string - '0';
    if (tmp < 0 || tmp > 9)    // how could this happen?
      return 0;

    value = value * 10 + tmp;

    rc++;

    string++;
  }

  if (1 == rc && 0 != qualifier)    // just a '>' or '<' by itself is an error
    return 0;

  if (qualifier < 0 && 0 == value)
  {
    cerr << "smarts_fetch_numeric: < 0 not allowed\n";
    return 0;
  }

  return rc;
}

//#define DEBUG_GET_ATOMIC_NUMBER_OR_SYMBOL

int
Substructure_Atom_Specifier::_get_atomic_number_or_symbol (const char * smarts,
                                                           const int characters_to_process) 
{
  if (0 == characters_to_process)
    return 0;

#ifdef DEBUG_GET_ATOMIC_NUMBER_OR_SYMBOL
  cerr << "Atomic number specification '";
  for (int i = 0; i < characters_to_process; ++i)
  {
    cerr << smarts[i];
  }
  cerr << "'\n";
#endif

  int nchars = 0;

  const Element * e = NULL;

  if ('{' == smarts[0] && characters_to_process > 2)
  {
    const_IWSubstring s(smarts+1, characters_to_process-1);   // skip over '{'

    int close_brace = s.index('}');
    if (close_brace <= 0)
      return 0;
    s.iwtruncate(close_brace);
#ifdef DEBUG_GET_ATOMIC_NUMBER_OR_SYMBOL
    cerr << "getting element for '" << s << "'\n";
#endif
    const_IWSubstring token;
    for (int i = 0; s.nextword(token, i, ','); )
    {
      e = get_element_from_symbol_no_case_conversion(token);

      if (NULL != e)
        _add_element(e);
      else if (auto_create_new_elements())
      {
        e = create_element_with_symbol(token);
        if (NULL == e)
          return 0;

        _add_element(e);
      }
      else
      {
        cerr << "Substructure_Atom_Specifier::_get_atomic_number_or_symbol:no element for '" << token << "'\n";
        return 0;
      }
    }
    nchars = 1 + close_brace + 1;
  }
  else
  {
    int z;
//  cerr << "FETCHING NUMBER FROM '" << smarts[0] << "', characters_to_process " << characters_to_process << endl;
    nchars = fetch_numeric_char(smarts, z, characters_to_process);
    if (0 == nchars)
      return 0;

    e = get_element_from_atomic_number(z);
    _add_element(e);
  }

#ifdef DEBUG_GET_ATOMIC_NUMBER_OR_SYMBOL
  cerr << "After _get_atomic_number_or_symbol element unique id\n";
  for (int i = 0; i < _element_unique_id.number_elements(); ++i)
  {
    cerr << " i = " << i << ' ' << _element_unique_id[i] << endl;
  }
#endif

  return nchars;
}

int 
Substructure_Atom_Specifier::_add_element (const atomic_number_t z)
{
  const Element * e = get_element_from_atomic_number(z);

  if (NULL == e)
    return 0;

  return _add_element(e);
}

int
Substructure_Atom_Specifier::_add_element (const Element * e)
{
  _element.add(e);
  _element_unique_id.add(e->unique_id());

//cerr << "Substructure_Atom_Specifier::_add_element:added " << e->symbol() << " unique_id " << e->unique_id() << endl;
   
  return 1;
}

//#define DEBUG_CONSTRUCT_FROM_SMARTS_TOKEN

/*
  An atomic smarts has been tokenised for us. Parse it.

  zop will contain any operator that has come before the token.
  This is needed for dealing with H. For example [O,S;H].
  Note that this is still a kludge
*/

int
Substructure_Atom_Specifier::construct_from_smarts_token (const const_IWSubstring & zsmarts)
{
  int not_operator = 0;
  int characters_to_process = zsmarts.length();

#ifdef DEBUG_CONSTRUCT_FROM_SMARTS_TOKEN
  cerr << "Specifier parsing smarts '" << zsmarts << "' " << characters_to_process << " chars\n";
#endif

  if (0 == characters_to_process)
  {
    cerr << "Substructure_Atom_Specifier::construct_from_smarts_token: empty smarts\n";
    return 0;
  }

// Special case of H by itself. Treating it out here keeps the loop below a little
// simpler

  if (1 == characters_to_process && 'H' == zsmarts[0])
    return _add_element(1);

  const char * initial_smarts_ptr = zsmarts.rawchars();
  const char * smarts = zsmarts.rawchars();

  int characters_processed = 0;

  int previous_token_was = SMARTS_PREVIOUS_TOKEN_UNSPECIFIED;

// We keep track of the formal charge. 
// Feb 99. In order to recognise '+0' and '-0' we need to keep track of 
// whether or not a formal charge specifier has been encountered;

  int fc = 0;
  int fc_encountered = 0;

// Dec 97, communication from Dave Weininger.

// To make SMILES and SMARTS rules as similar as possible, there is an
// important exception to the above rule: H in brackets is taken to be a
// hydrogen atomic symbol if it is the *first* elemental primitive in the

// May 98. Look at how v4.52 software manual, I'm not sure this is correct.
// Change to make hydrogen never recognised as an element (use #1 if you want it)

  int first_elemental_primitive_encountered = 0;

  while (characters_processed < characters_to_process)
  {
    const char s = *smarts;

#ifdef DEBUG_CONSTRUCT_FROM_SMARTS_TOKEN
    cerr << "Examining smarts character '" << s << "', " << characters_to_process << " characters to process\n";
#endif

    int nchars = 0;        // how many 'extra' characters consumed

    const Element * e;

//  When parsing this, it is convenient to know whether or not the next
//  character is a lowercase letter

    int next_char_is_lowercase_letter = 0;
    int next_char_is_digit = 0;
    int next_char_is_relational = 0;
    int next_char_is_charge = 0;

    if (characters_processed < characters_to_process)
    {
      char cnext = smarts[1];

      if (islower(cnext))
        next_char_is_lowercase_letter = 1;
      else if (isdigit(cnext))
        next_char_is_digit = 1;
      else if ('>' == cnext || '<' == cnext)
        next_char_is_relational = 1;
      else if ('+' == cnext || '-' == cnext)
        next_char_is_charge = 1;
    }

//  Oct 97. Change parsing rules for better consistency with Daylight.
//  Try to consume leading characters as an element specifier
//  Quickly ran into 'H' and 'D', which would otherwise be considered Hydrogen and Deuterium
//  Make the change that D and T are no longer recognised as elements
//  But then, what about He and Hf (and Ha if that's what Hahnium ends up as)

    if ('H' == s && (next_char_is_digit || next_char_is_relational))
    {
      int hh;
      int gtlt;
      nchars = smarts_fetch_numeric(smarts + 1, hh, gtlt);
      if (0 == nchars)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "Bad H qualifier");
        return 0;
      }
      if (1 == gtlt) //// '>' found
        _hcount.set_min(hh + 1);
      else if (-1 == gtlt)
        _hcount.set_max(hh - 1);
      else
        _hcount.add(hh);
    }
    else if ('H' == s && next_char_is_lowercase_letter && (nchars = element_from_smarts_string(smarts, characters_to_process - characters_processed, e)))
    {
      nchars--;
      _element.add(e);
      _element_unique_id.add(e->unique_id());
      first_elemental_primitive_encountered = 1;
    }

//  May 98, remove this conditional

//  else if ('H' == s && 0 == first_elemental_primitive_encountered)
//  {
//    _atomic_number.add(1);
//    first_elemental_primitive_encountered = 1;
//  }

//  The Daylight site says that H+ and H- mean Explicit Hydrogen

    else if ('H' == s && next_char_is_charge && ! first_elemental_primitive_encountered)
    {
      _add_element(1);
      first_elemental_primitive_encountered = 1;
    }

//  Isotopic Hydrogen must be detected

    else if ('H' == s && SMARTS_PREVIOUS_TOKEN_MASS == previous_token_was)
    {
      _add_element(1);
      first_elemental_primitive_encountered = 1;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
    }

    else if ('H' == s)
    {
      if (_h_means_exactly_one_hydrogen)
        _hcount.add(1);
      else
        _hcount.set_min(1);
    }
    else if (isupper(s) && characters_to_process >= 2 && (nchars = element_from_smarts_string(smarts, characters_to_process - characters_processed, e)))   // beware autocreate
    {
      nchars--;
      _add_element(e);
      first_elemental_primitive_encountered = 1;
      if (respect_aliphatic_smarts && e->organic())
        _aromaticity = NOT_AROMATIC;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
    }
    else if (isupper(s) && 'A' != s && 'D' != s && 'H' != s && 'R' != s && 'G' != s && 'T' != s && (nchars = element_from_smarts_string(smarts, characters_to_process - characters_processed, e)))   // beware autocreate
    {

#ifdef DEBUG_CONSTRUCT_FROM_SMARTS_TOKEN
      cerr << "Character is element " << e->symbol() << " nchars = " << nchars << "\n";
#endif

      nchars--;    // remember, we are looking at the base character
      _add_element(e);
      first_elemental_primitive_encountered = 1;
      if (respect_aliphatic_smarts && e->organic())
        _aromaticity = NOT_AROMATIC;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
    }

//  Beware the 'x' directive. Fortunately no elements end in x - unless someone
//  does autocreate, in which case Cx will be ambiguous...

//  Beware of the 'r' directive. Consider Or in a smarts.
//  If the element Or has been created, it will have been picked up
//  in the previous test. If not, we need to look for it here.

    else if (isupper(s) && next_char_is_lowercase_letter && 
             ('r' == smarts[1] || 'v' == smarts[1] || 'x' == smarts[1]) &&
             (NULL != (e = get_element_from_symbol_no_case_conversion(s))))
    {
      _add_element(e);
      first_elemental_primitive_encountered = 1;
      _aromaticity = NOT_AROMATIC;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
    }
    else if ('D' == s)       // degree specifier (note that Deuterium is not valid in a smarts)
    {
      int ncon;
      int gtlt;
      nchars = smarts_fetch_numeric(smarts + 1, ncon, gtlt);
      if (0 == nchars || ncon < 0)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "D specifier has no default");
        return 0;
      }

      if (1 == gtlt)    // '>' found
      {
        _ncon.set_min(ncon + 1);
      }
      else if (-1 == gtlt)    // '<' found
      {
        _ncon.set_max(ncon - 1);
      }
      else
        _ncon.add(ncon);
    }
    else if (isdigit(s))       // atomic mass specifier
    {
      int msdif;
      nchars = fetch_numeric_char(smarts, msdif, zsmarts.length() - characters_processed);
      previous_token_was = SMARTS_PREVIOUS_TOKEN_MASS;
      _isotope.add(msdif);
      nchars--;     // remember, nchars is really the number of extra characters consumed
    }
    else if ('<' == s)   // upper bound on isotope
    {
      int iso;
      nchars = fetch_numeric_char(smarts + 1, iso, zsmarts.length() - characters_processed - 1);
      if (0 == nchars)
      {
        cerr << "Substructure_Atom_Specifier::construct_from_smarts_token:invalid isotope <\n";
        return 0;
      }
      previous_token_was = SMARTS_PREVIOUS_TOKEN_MASS;
      _isotope.set_max(iso - 1);
      assert(_isotope.ok());
    }
    else if ('>' == s)   // lower bound on isotope
    {
      int iso;
      nchars = fetch_numeric_char(smarts + 1, iso, zsmarts.length() - characters_processed - 1);
      if (0 == nchars)
      {
        cerr << "Substructure_Atom_Specifier::construct_from_smarts_token:invalid isotope >\n";
        return 0;
      }
      previous_token_was = SMARTS_PREVIOUS_TOKEN_MASS;
      _isotope.set_min(iso + 1);
      assert (_isotope.ok());
    }
    else if ('*' == s)     // any atom specifier, don't do anything
    {
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
      first_elemental_primitive_encountered = 1;
    }
    else if ('a' == s)     // an aromatic atom
    {
      if (SUBSTRUCTURE_NOT_SPECIFIED == _aromaticity)
        _aromaticity = AROMATIC;
      else
        SET_AROMATIC_ATOM(_aromaticity);

      first_elemental_primitive_encountered = 1;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
    }
    else if ('A' == s)     // an aliphatic atom
    {
      if (SUBSTRUCTURE_NOT_SPECIFIED == _aromaticity)
        _aromaticity = NOT_AROMATIC;
      else
        SET_ALIPHATIC_ATOM(_aromaticity);
      first_elemental_primitive_encountered = 1;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
    }
    else if ('c' == s)
    {
      _add_element(6);
      _aromaticity = AROMATIC;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
    }
    else if ('n' == s)
    {
      _add_element(7);
      _aromaticity = AROMATIC;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
    }
    else if ('o' == s)
    {
      _add_element(8);
      _aromaticity = AROMATIC;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
    }
    else if ('h' == s)     // implicit h count specifier, ignored
    {
    }
    else if ('s' == s)
    {
      _add_element(16);
      _aromaticity = AROMATIC;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
    }
    else if ('p' == s)
    {
      _add_element(15);
      _aromaticity = AROMATIC;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
    }
    else if ('R' == s)     // number of rings specifier
    {
      int rr;
      int ltgt;
      nchars = smarts_fetch_numeric(smarts + 1, rr, ltgt);
      if (0 == nchars)
        _nrings.set_min(1);
      else if (ltgt > 0)
        _nrings.set_min(rr + 1);
      else if (ltgt < 0)
        _nrings.set_max(rr - 1);
      else
        _nrings.add(rr);

      previous_token_was = SMARTS_PREVIOUS_TOKEN_RING;
    }

//  March 2007. Beware of things like [rR1] and [R1r]

    else if ('r' == s)     // ring size specifier
    {
      int rr;
      int ltgt;
      nchars = smarts_fetch_numeric(smarts + 1, rr, ltgt);
      if (0 == nchars)      // just 'r' by itself means in a ring
        _nrings.set_min(1);
      else if (0 == rr)     // 'r0' means not in a ring
        _nrings.add(0);
      else if (rr + ltgt < 3)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "ring sizes must be >= 3");
        return 0;
      }
      else if (ltgt > 0)
        _ring_size.set_min(rr + 1);
      else if (ltgt < 0)
      {
        assert (rr > 3);      // no such thing as a 2 membered ring, so 'r<3' is nonsensical
        _ring_size.set_max(rr - 1);
      }
      else
        _ring_size.add(rr);
      previous_token_was = SMARTS_PREVIOUS_TOKEN_RING;
    }
    else if ('X' == s)     // connectivity - total connections
    {
      int xx;
      int ltgt;
      nchars = smarts_fetch_numeric(smarts + 1, xx, ltgt);
      if (0 == xx)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "the 'X' specifier must be followed by a whole number");
        return 0;
      }
      if (ltgt > 0)
        _daylight_x.set_min(xx + 1);
      else if (ltgt < 0)
        _daylight_x.set_max(xx - 1);
      else
        _daylight_x.add(xx);

      previous_token_was = SMARTS_PREVIOUS_TOKEN_X;
    }
    else if ('x' == s)     // ring bond count
    {
      int xx;
      int ltgt;
      nchars = smarts_fetch_numeric(smarts + 1, xx, ltgt);
      if (0 == nchars)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "the 'x' specifier must be followed by a whole number");
        return 0;
      }

      if (ltgt > 0)
        _ring_bond_count.set_min(xx + 1);
      else if (ltgt < 0)
        _ring_bond_count.set_max(xx - 1);
      else
        _ring_bond_count.add(xx);

      previous_token_was = SMARTS_PREVIOUS_TOKEN_RBC;
    }
    else if ('v' == s)     // total valence (nbonds)
    {
      int vv;
      int ltgt;
      nchars = smarts_fetch_numeric(smarts + 1, vv, ltgt);
      if (0 == vv)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "the 'v' specifier must be followed by a whole number");
        return 0;
      }
      else if (ltgt > 0)
        _nbonds.set_min(vv + 1);
      else if (ltgt < 0)
      {
        _nbonds.set_max(vv - 1);
      }
      else
        _nbonds.add(vv);

      previous_token_was = SMARTS_PREVIOUS_TOKEN_V;
    }
    else if ('-' == s)     // negative charge specifier
    {
      if (fc > 0)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "+ and - cannot be combined");
        return 0;
      }

      fc_encountered++;

      int ff;
      nchars = fetch_numeric_char (smarts + 1, ff, zsmarts.length() - characters_processed - 1);
      if (nchars)             // should do more error checking
        fc = -ff;
      else
        fc--;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_CHARGE;
    }
    else if ('+' == s)     // positive charge specifier
    {
      if (fc < 0)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "+ and - cannot be combined");
        return 0;
      }

      fc_encountered++;

      int ff;
      nchars = fetch_numeric_char(smarts + 1, ff, zsmarts.length() - characters_processed - 1);
      if (nchars)
        fc = ff;
      else
        fc++;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_CHARGE;
    }
    else if ('#' == s)     // atomic number specifier
    {
      nchars = _get_atomic_number_or_symbol(smarts+1, zsmarts.length() - characters_processed - 1);
      if (0 == nchars)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "# specifier has no default");
        return 0;
      }

      first_elemental_primitive_encountered = 1;
      previous_token_was = SMARTS_PREVIOUS_TOKEN_ELEMENT;
    }
    else if ('@' == s)     // chirality specifier
    {
      _chirality = snarf_capital_two(smarts, characters_to_process - characters_processed);
      if (_chirality > 3)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "too many @'s");
        return 0;
      }
      nchars = _chirality - 1;     // there was no leading letter before the @
//    cerr << "Chirality count " << (_chirality + 1) << ", nchars = " << nchars << endl;

      previous_token_was = SMARTS_PREVIOUS_TOKEN_CHIRALITY;
    }
    else if ('G' == s)    // unsaturation, iaw extension to smarts. Note we cannot use 'U' because that would be potentially ambiguous with Uranium
    {
      int u;
      int ltgt;
      nchars = smarts_fetch_numeric(smarts + 1, u, ltgt);
      if (0 == nchars)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "The unsaturation specifier 'G' has no default");
        return 0;
      }
      if (0 == ltgt)
        _unsaturation.add(u);
      else if (ltgt > 0)
        _unsaturation.set_min(u + 1);
      else if (ltgt < 0)
        _unsaturation.set_max(u - 1);
      previous_token_was = SMARTS_PREVIOUS_TOKEN_CHIRALITY;
    }
    else if ('T' == s)   // attached heteroatom count, iaw extension to smarts
    {
      int a;
      int ltgt;
      nchars = smarts_fetch_numeric(smarts + 1, a, ltgt);
      if (0 == nchars)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "The attached heteroatom count specifier 'T' has no default");
        return 0;
      }
      if (0 == ltgt)
        _attached_heteroatom_count.add(a);
      else if (ltgt > 0)
        _attached_heteroatom_count.set_min(a + 1);
      else if (ltgt < 0)
        _attached_heteroatom_count.set_max(a - 1);
      previous_token_was = SMARTS_PREVIOUS_TOKEN_T;
    }

//  Now the operators. Note the comma operator is not allowed

    else if ('!' == s)     // negation (tight bonding)
    {
      previous_token_was = SMARTS_PREVIOUS_TOKEN_OPERATOR_NOT;
      not_operator = 1;
      smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "! operator not supported");
      return 0;
    }
    else if ('&' == s)     // tight binding and
    {
      previous_token_was = SMARTS_PREVIOUS_TOKEN_OPERATOR_AND_AMPER;
    }
    else if (';' == s)     // low priority and
    {
      previous_token_was = SMARTS_PREVIOUS_TOKEN_OPERATOR_AND_AMPER;
    }
    else if (',' == s)
    {
      smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "Comma operator not allowed in specifier");
      return 0;
    }
    else if ('$' == s && characters_processed < characters_to_process && '(' == smarts[1])
    {
      const_IWSubstring environment;
      environment.set(smarts, characters_to_process - characters_processed);

      nchars = fetch_environment(environment);   // tuncates at end of environment
      if (0 == nchars)
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "invalid environment specifier");
        return 0;
      }

      previous_token_was = SMARTS_PREVIOUS_TOKEN_ENVIRONMENT;
    }
    else if ('/' == s && (characters_processed + 3) < characters_to_process && 'I' == smarts[1] && 'W' == smarts[2])
    {
      const_IWSubstring c(smarts + 3, characters_to_process - characters_processed - 3);    // handy here

      if ('x' == c[0])
        nchars = 3 + 1 - 1;
      else if (c.length() > 3 && c.starts_with("fss"))
      {
        c.remove_leading_chars(3);
        truncate_after_digits("=><,", c);
        if (! _fused_system_size.initialise(c))
        {
          smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "Invalid fss qualifier");
          return 0;
        }
        nchars = 3 + 3 + c.length() - 1;
      }
      else if (c.length() > 2 && c.starts_with("Vy"))
      {
        c.remove_leading_chars(2);
        truncate_after_digits("=<>", c);
        if (! _vinyl.initialise(c))
        {
          smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "Invalid Vy qualifier");
          return 0;
        }
        nchars = 3 + 2 + c.length() - 1;
      }
      else if (c.length() > 2 && c.starts_with("Ar"))
      {
        c.remove_leading_chars(2);
        truncate_after_digits("=<>", c);
        if (! _aryl.initialise(c))
        {
          smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "Invalid Ar qualifier");
          return 0;
        }

        nchars = 3 + 2 + c.length() - 1;
      }
      else if (c.length() > 3 && c.starts_with("rid"))
        nchars = 3 + 3 + 1 - 1;                     // will fail if more than two digits for rid
      else if (c.length() > 4 && c.starts_with("fsid"))
        nchars = 3 + 4 + 1 - 1;                     // will fail if more than two digits for rid
      else if (c.length() > 4 && c.starts_with("spch"))
      {
        c.remove_leading_chars(4);
        if ('1' == c[0])
          _match_spinach_only = 1;
        else if ('0' == c[0])
          _match_spinach_only = 0;
        else
        {
          smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "Invalid spch qualifier");
          return 0;
        }
        nchars = 3 + 4 + 1 - 1;
      }
      else if (c.length() > 2 && (c.starts_with("rh") || c.starts_with("hr")))
      {
        c.remove_leading_chars(2);
        truncate_after_digits ("=<>", c);
        if (! _heteroatoms_in_ring.initialise(c))
        {
          smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "Invalid rh/hr qualifier");
          return 0;
        }
        nchars = 3 + 2 + c.length() - 1;
      }
      else if (c.length() > 4 && c.starts_with("rscb"))
      {
        c.remove_leading_chars(4);
        truncate_after_digits("=<>", c);
        if (! _scaffold_bonds_attached_to_ring.initialise(c))
        {
          smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "Invalid rscb/hr qualifier");
          return 0;
        }
        nchars = 3 + 4 + c.length() - 1;
      }
      else if (c.starts_with("Nv"))
      {
        c.remove_leading_chars(2);
        nchars = 3 + 2;
        while (c.length() && isdigit(c[0]))
        {
          c += 1;
          nchars++;
        }
        nchars--;
      }
      else if (c.starts_with("symd"))
      {
        c.remove_leading_chars(4);
        truncate_after_digits("=<>", c);
        if (! _symmetry_degree.initialise(c))
        {
          smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "Invalid symd qualifier");
          return 0;
        }
        nchars = 3 + 4 + c.length() - 1;
      }
      else if (c.starts_with("symg"))
      {
        c.remove_leading_chars(4);
        truncate_after_digits("=<>", c);
        if (1 != c.length() || ! isdigit(c[0]) || '0' == c[0])
        {
          smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "Invalid symg qualifier");
          return 0;
        }
        _symmetry_group = c[0] - '0';
        nchars = 3 + 4 + c.length() - 1;
      }
      else if (c.starts_with("Kl"))
      {
        c.remove_leading_chars(2);
        if ('0' == c[0])
          _all_rings_kekule = 0;
        else if ('1' == c[0])
          _all_rings_kekule = 1;
        else
        {
          smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "Invalid KL qualifier");
          return 0;
        }
        nchars = 3 + 2 + 1 - 1;
//      _attributes_specified++;
      }
      else
      {
        smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "unrecognised /IW qualifier");
        return 0;
      }
    }
    else
    {
      smiles_error_message(initial_smarts_ptr, characters_to_process, characters_processed, "unrecognised SMARTS character");
      return 0;
    }

    if (nchars)
    {
      characters_processed += 1 + nchars;
      smarts += 1 + nchars;
    }
    else
    {
      characters_processed++;
      smarts++;
    }

//  cerr << "_match_spinach_only " << _match_spinach_only << endl;
  }

  if (fc_encountered)
    _formal_charge.add(fc);

  if (ignore_chirality_in_smarts_input())
    _chirality = 0;

#ifdef DEBUG_CONSTRUCT_FROM_SMARTS_TOKEN
  cerr << "Consumed " << characters_processed << " characters\n";
  cerr << "ASKLDJHASD _element_unique_id array contains " << _element_unique_id.number_elements() << " items\n";
#endif

  return characters_processed;
}

/*
  A non-bracketed element has been encountered in a smarts
*/

int
Substructure_Atom_Specifier::construct_from_smiles_token (const char * smi, int nchars)
{
  const_IWSubstring tmp;
  tmp.set(smi, nchars);

  return construct_from_smiles_token(tmp);
}

int
Substructure_Atom_Specifier::construct_from_smiles_token (const const_IWSubstring & smiles)
{
  _attributes_specified++;

  if ('a' == smiles[0])
  {
    _aromaticity = AROMATIC;
    return 1;
  }

  if ('A' == smiles[0])
  {
    _aromaticity = NOT_AROMATIC;
    return 1;
  }

  const Element * e;

  int nchars = parse_smiles_token(smiles.rawchars(), smiles.length(), e, _aromaticity);
  if (0 == nchars || NULL == e)
  {
    cerr << "Substructure_Atom_Specifier::construct_from_smiles_token: cannot parse '" << smiles << "'\n";
//  iwabort();
    return 0;
  }

  if (! respect_aliphatic_smarts)
    _aromaticity = SUBSTRUCTURE_NOT_SPECIFIED;
  else if (AROMATIC == _aromaticity)
    ;
  else
    _aromaticity = NOT_AROMATIC;

  _add_element(e);

  _attributes_specified++;
    
  return nchars;
}

/*
  during parsing a smarts, we split the smarts into and blocks (via the ';' char).
  These must be reconciled.

  Note that this is broken.  For each attribute, we really need to
  check our components
*/

int
Substructure_Atom_Specifier::reconcile_and_conditions (const Substructure_Atom_Specifier * s)
{
  if (s->_element.number_elements())
  {
    if (_element.number_elements())
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: element conflict\n";
      return 0;
    }

    _element += s->_element;
  }

  if (s->_ncon.is_set())
  {
    if (_ncon.is_set())
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: ncon conflict\n";
      return 0;
    }

    _ncon = s->_ncon;
  }

  if (s->_ncon2.is_set())
  {
    if (_ncon2.is_set())
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: ncon2 conflict\n";
      return 0;
    }

    _ncon2 = s->_ncon2;
  }

  if (s->_nbonds.is_set())
  {
    if (_nbonds.is_set())
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: nbonds conflict\n";
      return 0;
    }

    _nbonds = s->_nbonds;
  }

  if (s->_formal_charge.is_set())
  {
    if (_formal_charge.is_set())
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: formal_charge conflict\n";
      return 0;
    }

    _formal_charge = s->_formal_charge;
  }

  if (s->_nrings.is_set())
  {
    if (_nrings.is_set())
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: nrings conflict\n";
      return 0;
    }

    _nrings = s->_nrings;
  }

  if (s->_ring_size.is_set())
  {
    if (_ring_size.is_set())
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: ring_size conflict\n";
      return 0;
    }

    _ring_size = s->_ring_size;
  }

  if (s->_hcount.is_set())
  {
    if (_hcount.is_set())
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: hcount conflict\n";
      return 0;
    }

    _hcount = s->_hcount;
  }

  if (SUBSTRUCTURE_NOT_SPECIFIED != s->_aromaticity)
  {
    if (SUBSTRUCTURE_NOT_SPECIFIED != _aromaticity)
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: aromaticity conflict\n";
      return 0;
    }

    _aromaticity = s->_aromaticity;
  }

// do chirality sometime

  if (s->_attached_heteroatom_count.is_set())
  {
    if (_attached_heteroatom_count.is_set())
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: attached_heteroatom_count conflict\n";
      return 0;
    }

    _attached_heteroatom_count = s->_attached_heteroatom_count;
  }

  if (s->_aromatic_ring_sizes.is_set())
  {
    if (_aromatic_ring_sizes.is_set())
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: aromatic_ring_sizes conflict\n";
      return 0;
    }

    _aromatic_ring_sizes = s->_aromatic_ring_sizes;
  }

  if (s->_aliphatic_ring_sizes.is_set())
  {
    if (_aliphatic_ring_sizes.is_set())
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: aliphatic_ring_sizes conflict\n";
      return 0;
    }

    _aliphatic_ring_sizes = s->_aliphatic_ring_sizes;
  }

  if (s->_unsaturation.is_set())
  {
    if (_unsaturation.is_set())
    {
      cerr << "Substructure_Atom_Specifier::reconcile_and_conditions: unsaturation conflict\n";
      return 0;
    }

    _unsaturation = s->_unsaturation;
  }

#ifdef DEBUG_ATOM_CONSTRUCT_FROM_SMARTS_TOKEN
  cerr << "And reconciliaton done\n";
#endif

  return 1;
}

/*int
Substructure_Atom_Specifier::smarts (IWstring & s) const
{
  aromaticity_type_t arom = _aromaticity;
  if (AROMATICITY_NOT_DETERMINED == arom)
    arom = NOT_AROMATIC;

  int na = _atomic_number.number_elements();
  for (int i = 0; i < na; i++)
  {
    atomic_number_t z = _atomic_number[i];
    if (i > 0)
      s += ',';


    const Element * e = get_element_from_atomic_number (z);
    assert (e);

    int needs_sqb = e->needs_square_brackets();

    if (needs_sqb)
      s += '[';

    IWString tmp = e->symbol();
    if (AROMATIC == arom)
    {
    }

    s += tmp;

    if (needs_sqb)
      s += ']';
  }

  if (! _ncon.is_set())
    ;
  else if (0 == _ncon.number_elements())
    s += "{ncon min or max"};
  else
  {
    for (int i = 0; i < _ncon.number_elements(); i++)
    {
      if (i > 0)
        s += ',';

      append_int (s, _ncon[i]);
    }
  }

  return 1;
}*/

/*
  Before initiating a search, we need to determin the number of attributes
  which have been specified.
*/

int
Substructure_Atom_Specifier::attributes_specified()
{
  int rc = 0;
//if (_element.number_elements())
//  rc++;
  if (_element_unique_id.number_elements())
    rc++;
  if (_ncon.is_set())
    rc++;
  if (_ncon2.is_set())
    rc++;
  if (_nbonds.is_set())
    rc++;
  if (_formal_charge.is_set())
    rc++;
  if (_nrings.is_set())
    rc++;
  if (_ring_bond_count.is_set())
    rc++;
  if (_ring_size.is_set())
    rc++;
  if (_hcount.is_set())
    rc++;
  if (SUBSTRUCTURE_NOT_SPECIFIED != _aromaticity)
    rc++;
  if (SUBSTRUCTURE_NOT_SPECIFIED != _chirality)
    rc++;
  if (_aromatic_ring_sizes.is_set())
    rc++;
  if (_aliphatic_ring_sizes.is_set())
    rc++;
  if (_attached_heteroatom_count.is_set())
    rc++;
  if (_lone_pair_count.is_set())
    rc++;
  if (_unsaturation.is_set())
    rc++;
  if (_daylight_x.is_set())
    rc++;
  if (_isotope.is_set())
    rc++;
  if (_aryl.is_set())
    rc++;
  if (_match_spinach_only >= 0)
    rc++;
  if (_fused_system_size.is_set())
    rc++;
  if (_vinyl.is_set())
    rc++;
  if (_heteroatoms_in_ring.is_set())
    rc++;
  if (_scaffold_bonds_attached_to_ring.is_set())
    rc++;  
  if (_symmetry_degree.is_set())
    rc++;
  if (SUBSTRUCTURE_NOT_SPECIFIED != _all_rings_kekule)
    rc++;

//#define DEBUG_ATTRIBUTES_SPECIFIED
#ifdef DEBUG_ATTRIBUTES_SPECIFIED
  if (_element.number_elements())
    cerr << "ele is specified \n";
  if (_ncon.is_set())
    cerr << "nc is specified \n";
  if (_ncon2.is_set())
    cerr << "nc2 is specified \n";
  if (_nbonds.is_set())
    cerr << "nb is specified \n";
  if (_formal_charge.is_set())
    cerr << "fc is specified \n";
  if (_nrings.is_set())
    cerr << "nr is specified \n";
  if (_ring_bond_count.is_set())
    cerr << "rbc is specified \n";
  if (_ring_size.is_set())
    cerr << "rs is specified \n";
  if (_hcount.is_set())
    cerr << "hc is specified \n";
  if (SUBSTRUCTURE_NOT_SPECIFIED != _aromaticity)
    cerr << "ar is specified \n";
  if (SUBSTRUCTURE_NOT_SPECIFIED != _chirality)
    cerr << "chir is specified \n";
  if (_aromatic_ring_sizes.is_set())
    cerr << "ars is specified \n";
  if (_aliphatic_ring_sizes.is_set())
    cerr << "alrs is specified \n";
  if (_attached_heteroatom_count.is_set())
    cerr << "ahc is specified \n";
  if (_lone_pair_count.is_set())
    cerr << "lp is specified \n";
  if (_unsaturation.is_set())
    cerr << "us is specified \n";
  if (_daylight_x.is_set())
    cerr << "Daylight X is specified \n";
  if (_isotope.is_set())
    cerr << "is is specified\n";
  if (_aryl.is_set())
    cerr << "aryl is specified\n";
  if (_vinyl.is_set())
    cerr << "vinyl is specified\n";
  if (_heteroatoms_in_ring.is_set())
    cerr << "heteroatoms in ring specified\n";
  if (_fused_system_size.is_set())
    cerr << "fused system size is specified\n";

  cerr << "Atom has " << rc << " attributes specified\n";
#endif

  _attributes_specified = rc;


  return rc;
}
