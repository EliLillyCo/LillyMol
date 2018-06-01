#include <stdlib.h>
#include <memory>

#include "misc.h"

//#define TESTING_STUFF
#ifdef TESTING_STUFF
#include "cmdline.h"
#include "aromatic.h"
#endif

#include "iwreaction.h"
#include "molecule_to_query.h"

static int smirks_lost_atom_means_remove_frgment = 0;

void
set_smirks_lost_atom_means_remove_frgment(int s)
{
  smirks_lost_atom_means_remove_frgment = s;
}

int
IWReaction::construct_from_smirks(const const_IWSubstring & smirks)
{
  if (! smirks.starts_with("F:"))
    return _construct_from_smirks(smirks);

  IWString fname(smirks);
  fname.remove_leading_chars(2);

  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good())
  {
    cerr << "IWReaction::construct_from_smirks:cannot open '" << fname << "'\n";
    return 0;
  }

  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "IWReaction::const_IWSubstring:cannot read smirks from '" << fname << "'\n";
    return 0;
  }

  return _construct_from_smirks(buffer);
}

int
IWReaction::_construct_from_smirks(const const_IWSubstring & smirks)
{
  int square_bracket_level = 0;
  const int n = smirks.length();

  int first_open_angle = -1;
  int second_open_angle = -1;

  for (int i = 0; i < n; ++i)
  {
    const char c = smirks[i];
    if (']' == c)
      square_bracket_level--;
    else if ('[' == c)
      square_bracket_level++;
    else if (0 == square_bracket_level && '>' == c)
    {
      if (first_open_angle < 0)
        first_open_angle = i;
      else
      {
        second_open_angle = i;
        break;
      }
    }
  }

  if (second_open_angle < 0)
  {
    cerr << "IWReaction::build_from_smirks:yipes, did not find >> construct '" << smirks << "'\n";
    return 0;
  }

  const_IWSubstring reagents, products;

  smirks.from_to(0, first_open_angle - 1, reagents);

//cerr << "Reagents " << reagents << endl;

  smirks.from_to(second_open_angle + 1, smirks.length() - 1, products);

//cerr << "products " << products << endl;

  return construct_from_smirks(reagents, products);
}

void
fetch_closing_paren(const_IWSubstring & buffer,
                     IWString & token)
{
  int paren_level = 1;
  for (auto i = 1; i < buffer.length(); ++i)
  {
    const auto c = buffer[i];
    if ('(' == c)
      paren_level++;
    else if (')' == c)
    {
      paren_level--;

      if (0 == paren_level)
      {
        buffer.from_to(1, i-1, token);
        buffer += i;
        return;
      }
    }
  }

  cerr << "fetch_closing_paren:did not find closing paren '" << buffer << "'\n";   // should not happen

  return;
}

/*
  We try to match the Daylight component level grouping
  thing, WITH THE EXCEPTION OF THE ZERO LEVEL CONCEPT

  So, (smarts1.smarts2).smarts3

  means smarts1 and smarts2 must be in the scaffold, and smarts3 is
  in the sidechain.
*/

static int
tokenise_smarts_into_components(const_IWSubstring smarts,    // pass by value
                                resizable_array_p<IWString> & components)
{
//cerr << "Smarts '" << smarts << "'\n";

  while (smarts.length() > 0)
  {
//  cerr << "  smarts '" << smarts << "'\n";

    if (smarts.starts_with('('))
    {
      IWString * c = new IWString;
      fetch_closing_paren(smarts, *c);
      components.add(c);
      smarts++;
      if (smarts.starts_with('.'))
        smarts++;
    }
    else
    {
      IWString * s = new IWString;

//    cerr << "Check " << smarts.length() << " characters\n";
      for (int i = 0; i < smarts.length(); ++i)
      {
        const char c = smarts[i];

        if ('.' != c)
        {
          s->add(c);
          continue;
        }

        if (smarts.matches_at_position(i, "...", 3))
        {
          *s << "...";
          i += 2;
          continue;
        }

//      cerr << "Found '" << *s << "'\n";
        components.add(s);

        smarts += i;
//      cerr << "Smarts updated to " << smarts << endl;
        break;
      }
      if (s->length())
      {
        components.add(s);
        smarts += s->length();
      }
    }
  }

  return components.number_elements();
}

#ifdef NO_LONGER_USED_JJJJJ
static void
file_scope_identify_atom_map_numbers(const Substructure_Query & q,
                          extending_resizable_array<int> & numbers_in_use)
{
  for (int i = 0; i < q.number_elements(); ++i)
  {
    q[i]->identify_atom_map_numbers(numbers_in_use);
  }

  return;
}

static void
file_scope_identify_atom_map_numbers(const Molecule & m,
                                     extending_resizable_array<int> & numbers_in_use)
{
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    const auto a = m.atom_map_number(i);

    if (0 == a)
      continue;

    numbers_in_use[a]++;
  }

  return;
}
#endif

int
IWReaction::construct_from_smirks(const const_IWSubstring & reagents,
                                  const const_IWSubstring & products)
{
  if (0 != reagents.balance('(', ')'))
  {
    cerr << "IWReaction::construct_from_smirks:unbalanced parentheses in reagent specification '" << reagents << "'\n";
    return 0;
  }

  resizable_array_p<IWString> smarts;

  if (! tokenise_smarts_into_components(reagents, smarts))
  {
    cerr << "IWReaction::construct_from_smirks:cannot separate reagent smarts\n";
    return 0;
  }

//#define DEBUG_SMIRKS_TOKENISING
#ifdef DEBUG_SMIRKS_TOKENISING
  cerr << "From '" << reagents << "' generate\n";
  for (auto i = 0; i < smarts.number_elements(); ++i)
  {
    const IWString & s = *(smarts[i]);

    cerr << ' ' << s << '\n';
  }
#endif

  const IWString & s = *(smarts[0]);
  if (! this->create_from_smarts(s))
  {
    cerr << "IWReaction::construct_from_smirks:cannot interpret scaffold smarts '" << s << "'\n";
    return 0;
  }

  for (auto i = 1; i < smarts.number_elements(); ++i)
  {
    Sidechain_Reaction_Site * r = new Sidechain_Reaction_Site;

    if (! r->create_from_smarts(*smarts[i]))
    {
      cerr << "IWReaction::construct_from_smirks:invalid sidechain smarts '" << (*smarts[i]) << "'\n";
      delete r;
      return 0;
    }

    _sidechains.add(r);
  }

  Substructure_Query product_molecule;

  if (! product_molecule.create_from_smarts(products))
  {
    cerr << "IWReaction::construct_from_smirks:invalid product " << products << endl;
    return 0;
  }

// the RHS cannot contain ambiguous elements

  if (product_molecule.any_query_atom([] (const Substructure_Atom * q) { if (q->element().number_elements() > 1)
                                                                   return 1;
                                                                 else
                                                                   return 0; }))
  {
    cerr << "IWReaction::_identify_changes_from_smirks:RHS has multiple elements, impossible\n";
    return 0;
  }

// GAther up all the atom map numbers that are in use on both sides

  extending_resizable_array<int> ramn;     // reagent atom map numbers
  this->identify_atom_map_numbers(ramn);

  const auto ns = _sidechains.number_elements();

  for (int i = 0; i < ns; ++i)
  {
    _sidechains[i]->identify_atom_map_numbers(ramn);
  }

#ifdef DEBUG_CONSTRUCT_FROM_SMIRKS
  cerr << "Reagent atom numbers";
  for (int i = 0; i < ramn.number_elements(); ++i)
  {
    if (ramn[i] > 0)
      cerr << ' ' << i;
  }
  cerr << endl;
#endif

// On the LHS, each atom map number must be used once

  int h = 0;    // highest atom man number in use

  for (int i = 1; i < ramn.number_elements(); ++i)
  {
    if (0 == ramn[i])
      continue;

    h = i;

    if (1 == ramn[i])
      continue;

    cerr << "IWReaction::construct_from_smirks:atom map number " << i << " found " << ramn[i] << " occurrences\n";
    return 0;
  }

  const int h2 = product_molecule.highest_atom_map_number();
  if (h2 > h)
    h = h2;

#ifdef DEBUG_CONSTRUCT_FROM_SMIRKS
  cerr << "Highest atom map number in use " << h << endl;
#endif

// We need to assign atom map numbers to unmapped atoms in the product. These are, by definition
// orphan atoms.

  product_molecule.assign_atom_map_numbers(h);

#ifdef DEBUG_CONSTRUCT_FROM_SMIRKS
  cerr << "After assigning atom map number to product atoms, h = " << h << endl;
#endif

// Ensure that all atom numbers are present. Identify any orphans - atoms in the products that
// do not have a corresponding atom in the reagents

  extending_resizable_array<int> pamn(0);
  product_molecule.identify_atom_map_numbers(pamn);

  if (0 == pamn.number_elements())
  {
    cerr << "IWReaction::construct_from_smirks:No mapped atoms on RHS, cannot continue\n";
    return 0;
  }

  int istop = h + 1;

  resizable_array<int> orphan_atoms;

  for (int i = 1; i < istop; i++)
  {
#ifdef DEBUG_CONSTRUCT_FROM_SMIRKS
    cerr << "i = " << i << " ran " << ramn[i] << " pan " << pamn[i] << endl;
#endif

    if (1 == ramn[i] && 1 == pamn[i])    // one occurrence on LHS and RHS
      continue;

    if (0 == ramn[i] && 0 == pamn[i])    // atom map number not used either side, OK
      continue;

    if (ramn[i] == pamn[i])    // equal, but not 1, bad
    {
      cerr << "IWReaction::construct_from_smirks:atom map " << i << " used " << ramn[i] << " in reagents and " << pamn[i] << " in products\n";
      return 0;
    }

    if (1 == ramn[i] && 0 == pamn[i])    // atom disappears, this is OK
      continue;

    if (0 == ramn[i] && 1 == pamn[i])
    {
      orphan_atoms.add(i);
      continue;
    }

    cerr << "IWReaction::construct_from_smirks:atom number mismatch, atom " << i << " in reagents " << ramn[i] << " in products " << pamn[i] << endl;
    return 0;
  }

  Molecule orphan_molecule;

  if (orphan_atoms.number_elements())
  {
    Sidechain_Reaction_Site * x = new Sidechain_Reaction_Site();

    if (! _create_orphan_molecule(orphan_atoms, orphan_molecule, product_molecule, *x))
    {
      cerr << "IWReaction::_identify_changes_from_smirks:cannot construct orphan molecule, has " << orphan_atoms.number_elements() << " orphan atoms\n";
      delete x;
      return 0;
    }
    
    _sidechains.add(x);
  }

  for (int i = 0; i < _sidechains.number_elements(); ++i)
  {
    _sidechains[i]->set_sidechain_number(i);
  }

  IWString foo ("FOO.qry");
  product_molecule.write_msi(foo);

  return _identify_changes_from_smirks(ramn, h+1, orphan_molecule, product_molecule);
}

static int
add_bond_to_orphan(Molecule & orphan_molecule,
                   const atom_number_t a1,
                   const atom_number_t a2,
                   const Substructure_Bond * b)
{
  if (orphan_molecule.are_bonded(a1, a2))    // this process will find most bonds twice
    return 0;

  const bond_type_t bt = b->types_matched();
  if (IS_SINGLE_BOND(bt))
    orphan_molecule.add_bond(a1, a2, SINGLE_BOND);
  else if (IS_DOUBLE_BOND(bt))
    orphan_molecule.add_bond(a1, a2, DOUBLE_BOND);
  else if (IS_TRIPLE_BOND(bt))
    orphan_molecule.add_bond(a1, a2, TRIPLE_BOND);
  else
  {
    cerr << "IWReaction::_create_orphan_molecule:bond type not unitary form " << bt << ", assuming single bond!!\n";
    orphan_molecule.add_bond(a1, a2, SINGLE_BOND);
  }

  return 1;
}

int
IWReaction::_create_orphan_molecule(const resizable_array<int> & orphan_atoms,    // atom map numbers
                                    Molecule & orphan_molecule,
                                    const Substructure_Query & product_molecule,
                                    Sidechain_Reaction_Site & sc)
{
  const int n = orphan_atoms.number_elements();
  const int mx = orphan_atoms.max_val();

// convenient to have two cross reference arrays

  int * mol2qry = new int[n];std::unique_ptr<int[]> free_mol2qry(mol2qry);
  int * amap2mol = new_int(mx + 1, -1); std::unique_ptr<int[]> free_amap2mol(amap2mol);

  const Substructure_Atom ** query_atoms = new const Substructure_Atom *[n]; std::unique_ptr<const Substructure_Atom *[]> free_query_atoms(query_atoms);
//resizable_array<const Substructure_Atom *> query_atoms;

  for (int i = 0; i < orphan_atoms.number_elements(); ++i)
  {
    const int oi = orphan_atoms[i];

    const Substructure_Atom * qi = product_molecule.query_atom_with_atom_map_number(oi);

    const Substructure_Atom_Specifier * si;
    if (qi->ncomponents() > 0)              // was complex smarts
      si = qi->component(0);
    else                                    // was simple
      si = qi;

    assert (NULL != qi);

    const auto & e = si->element();
    if (0 == e.number_elements())
    {
      cerr << "IWReaction::_create_orphan_molecule:RHS orphan atom has no element specified, impossible. Atom map " << oi << endl;
      return 0;
    }

    const int x = orphan_molecule.natoms();

    query_atoms[x] = qi;     // the X'th query atom

    amap2mol[oi] = x;     // atom map OI is atom X in molecule

    orphan_molecule.add(e[0]);
    orphan_molecule.set_atom_map_number(x, oi);

    const auto & iso = si->isotope();
    if (iso.number_elements() > 0)
      orphan_molecule.set_isotope(x, iso[0]);
    const auto & fc = si->formal_charge();
    if (fc.number_elements() > 0)
      orphan_molecule.set_formal_charge(x, fc[0]);
  }

// Now the atoms have been added do the bonds. We need to handle both bonds within the orphan
// atoms, and bonds to other parts of the products

  for (int i = 0; i < orphan_atoms.number_elements(); ++i)
  {
    const Substructure_Atom * qi = query_atoms[i];

    const int nc = qi->number_children();

    for (int j = 0; j < nc; ++j)   // scan all atoms attached to the corresponding query atom
    {
      const Substructure_Atom * c = qi->child(j);
      const Substructure_Bond * b = c->bond_to_parent();

      const int atom_number_in_orphan = amap2mol[c->atom_map_number()];

      if (atom_number_in_orphan >= 0)
        add_bond_to_orphan(orphan_molecule, i, atom_number_in_orphan, b);
    }

    const auto & bonds = qi->bonds();

    for (int j = 0; j < bonds.number_elements(); ++j)
    {
      const Substructure_Bond * b = bonds[j];
      const Substructure_Atom * a = b->a();
      const int atom_number_in_orphan = amap2mol[a->atom_map_number()];
      if (atom_number_in_orphan >= 0)
        add_bond_to_orphan(orphan_molecule, i, atom_number_in_orphan, b);
    }
  }

  sc.set_single_reagent(orphan_molecule);

#ifdef DEBUG_CONSTRUCT_FROM_SMIRKS
  cerr << "Build orphan molecule " << orphan_molecule.smiles() << endl;
#endif

  return 1;
}

/*
  There appear to be two fundamentally different type of products.

  if all product atoms are mapped to a query atom on the LHS, this reaction can be used
  for enumerations of numerous reagents.

  If there are orphan atoms on the RHS, not present on the LHS, then we must assume this is a
  fixed reagent sidechain. 

  As we look at the molecules on the RHS, there will be two kinds of atoms.
  Those with a corresponding atom number on the reagent side.
    For those we just set the element, charge and isotope.
  Those which do not appear on the reagent side.
    for these, we must figure out how they are attached to those that appear on the LHS
*/

int
IWReaction::_identify_changes_from_smirks(const extending_resizable_array<int> & atom_map_numbers_in_reagents,
                                          const int istop,
                                          const Molecule & orphan_molecule,
                                          const Substructure_Query & product_molecule)
{
#ifdef DEBUG_CONSTRUCT_FROM_SMIRKS
  cerr << "IWReaction::_identify_changes_from_smirks:must examine " << istop << " mapped atoms\n";

  product_molecule.print_connectivity_graph(cerr);
#endif

  resizable_array<int> atoms_lost;

  for (int i = 0; i < istop; i++)
  {
//  cerr << "Atom map i " << i << " mapped in reagents " << atom_map_numbers_in_reagents[i] << endl;
    if (atom_map_numbers_in_reagents[i] <= 0)   // atom map number not present on LHS, orphan
      continue;

    Reaction_Site * r1 = _reaction_site_with_atom_map_number(i);
    assert (NULL != r1);

    const Substructure_Atom * q1 = r1->query_atom_with_atom_map_number(i);
    assert (NULL != q1);

    const Substructure_Atom * q2 = product_molecule.query_atom_with_atom_map_number(i);

    if (NULL == q2)     // atom eliminated in products
      atoms_lost.add(i);
    else
      _discern_atomic_changes(*r1, *q1, *q2);
  }

  if (atoms_lost.number_elements())
  {
    for (int i = 0; i < atoms_lost.number_elements(); ++i)
    {
      const int j = atoms_lost[i];

      Reaction_Site * r = _reaction_site_with_atom_map_number(j);
      assert (NULL != r);

      const Substructure_Atom * s = r->query_atom_with_atom_map_number(j);
      assert (NULL != s);

      const int x = s->initial_atom_number();

      if (smirks_lost_atom_means_remove_frgment)
        r->add_fragment_to_be_removed(x);
      else
        r->add_atom_to_be_removed(x);
    }
  }

#ifdef DEBUG_CONSTRUCT_FROM_SMIRKS
  for (int i = 1; i < istop; ++i)
  {
    cerr << i << " atom_map_numbers_in_reagents " << atom_map_numbers_in_reagents[i] << endl;
  }
#endif

// Any bond changes. this must be done on a per sidechain basis because there may be bonding
// changes that occur just within a given sidechain
// e.g. ([CH3:1].[NH2D1:2]).[OH:3]-[C:4]=[O:5]>>[CH2:1]-[NH:2]-[C:4]=[O:5]

  for (int i = 0; i < istop; i++)
  {
//  cerr << i << " atom_map_numbers_in_reagents " << atom_map_numbers_in_reagents[i] << endl;

    Reaction_Site * i_site = _reaction_site_with_atom_map_number(i);    // may be NULL

    for (int j = i + 1; j < istop; j++)
    {
//    cerr << i << ',' << j << " atom_map_numbers_in_reagents " << atom_map_numbers_in_reagents[i] << " j " << ' ' << atom_map_numbers_in_reagents[j] << endl;
      if (atom_map_numbers_in_reagents[i] <= 0 && atom_map_numbers_in_reagents[j] <= 0)   // numbers not used
        continue;

      Reaction_Site * j_site = _reaction_site_with_atom_map_number(j);    // may be NULL

      if (NULL == i_site && NULL == j_site)          // should not happen
        continue;

      const Substructure_Bond * bonded_in_product = product_molecule.bond_between_atom_map_numbers(i, j);  
//    cerr << "Between mapped atoms " << i << " and " << j << " bnd " << bonded_in_product << " Product\n";

      const Substructure_Bond * bonded_in_reagent;
      if (i_site == j_site)           // two mapped atoms in same reagent
        bonded_in_reagent = i_site->bond_between_atom_map_numbers(i, j);
      else           // not in same reagent, therefore not bonded
        bonded_in_reagent = NULL; 

//    cerr << "bonded? " << bonded_in_reagent << ' ' << bonded_in_product << endl;

      if (NULL == bonded_in_reagent && NULL == bonded_in_product)    // no bond btw I and J either side
        continue;

      bond_type_t btr = NULL == bonded_in_reagent ? INVALID_BOND_TYPE : BOND_TYPE_ONLY(bonded_in_reagent->types_matched());
      bond_type_t btp = NULL == bonded_in_product ? INVALID_BOND_TYPE : BOND_TYPE_ONLY(bonded_in_product->types_matched());

      if (btr == btp)     // same bond type, no need to change anything
        continue;

//    the bonding changes. Is it a removal, an addition or a change? What components are involved?

//    cerr << "Between atoms " << i << " and " << j << " reagent? " << bonded_in_reagent << " (" << btr << ") products " << bonded_in_product << " (" << btp << ")\n";
//    cerr << "Sites " << i_site << " and " << j_site << endl;

      if (i_site == j_site)     // new/changed or removed bond in the same component
      {
        const Substructure_Atom * ai = i_site->query_atom_with_atom_map_number(i);
        const Substructure_Atom * aj = j_site->query_atom_with_atom_map_number(j);
//      cerr << "Atom numbers " << ai->initial_atom_number() << " and " << aj->initial_atom_number() << " xBTP " << btp << endl;
        if (INVALID_BOND_TYPE == btp)       // bond is removed
          i_site->add_bond_to_be_broken(ai->initial_atom_number(), aj->initial_atom_number());
        else                    // new or changed
          i_site->add_bond_to_be_made(ai->initial_atom_number(), aj->initial_atom_number(), btp);

        continue;
      }

      if (NULL != i_site && NULL != j_site)     // new inter partile bond btw existing reagents
      {
//      cerr << "Adding inter particle bond btw " << i << " and " << j << endl;
        if (! _from_smirks_add_inter_particle_bond(i, i_site, j, j_site, btp))
          return 0;
        continue;
      }

//    One of these atom map numbers is not present on the LHS

//    cerr << "Bond involves different sites. I " << i_site << " J " << j_site << " type " << btp << endl;
      assert (INVALID_BOND_TYPE != btp);

      if (NULL == i_site)
      {
        const Substructure_Atom * q2 = j_site->query_atom_with_atom_map_number(j);
        const atom_number_t x2 = q2->initial_atom_number();

        const atom_number_t x1 = orphan_molecule.atom_with_atom_map_number(i);
        const int r = _sidechain_with_mapped_atom(j);

        _sidechains.last_item()->add_inter_particle_bond(r, x2, x1, btp);
      }
      else if (NULL == j_site)
      {
        const Substructure_Atom * q1 = i_site->query_atom_with_atom_map_number(i);
        const atom_number_t x1 = q1->initial_atom_number();

        const atom_number_t x2 = orphan_molecule.atom_with_atom_map_number(j);
        const int r = _sidechain_with_mapped_atom(i);
        _sidechains.last_item()->add_inter_particle_bond(r, x1, x2, btp);
      }
      else    // HUH!
        return 0;
    }
  }

// deal with chirality sometime - hard...

  return 1;
}

/*
  This is complicated, because we do not know which reaction component
  contains the two atoms.
  It might be from the scaffold to one of the sidechains, or it might
  involve two sidechains
*/

int
IWReaction::_from_smirks_add_inter_particle_bond(int amap1, 
                                                 Reaction_Site * a1site,
                                                 int amap2,
                                                 Reaction_Site * a2site,
                                                 bond_type_t bt)
{
  assert (NULL != a1site);
  assert (NULL != a2site);

  const Substructure_Atom * a1 = a1site->query_atom_with_atom_map_number(amap1);
  const Substructure_Atom * a2 = a2site->query_atom_with_atom_map_number(amap2);
  assert (NULL != a1);
  assert (NULL != a2);

// If either one of the sites is the scaffold, handle specially

  if (a1site == this)
  {
    const int r2 = _sidechain_with_mapped_atom(amap2);
    _sidechains[r2]->add_inter_particle_bond(-1, a1->initial_atom_number(), a2->initial_atom_number(), bt);
  }
  else if (a2site == this)
  {
    const int r1 = _sidechain_with_mapped_atom(amap1);
    _sidechains[r1]->add_inter_particle_bond(-1, a2->initial_atom_number(), a1->initial_atom_number(), bt);
  }
  else         // different sidechains
  {
    const int r1 = _sidechain_with_mapped_atom(amap1);
    const int r2 = _sidechain_with_mapped_atom(amap2);
    _sidechains[r1]->add_inter_particle_bond(r2, a2->initial_atom_number(), a1->initial_atom_number(), bt);
  }
  
  return 1;
}

/*
  Atom A1 is in the scaffold. We need to find which sidechain has a2
*/

int
IWReaction::_from_smirks_add_inter_particle_bond_involves_scaffold(int a1, int a2,
                                                                    bond_type_t bt)
{
  int r2 = -1;

  for (auto i = 0; i < _sidechains.number_elements(); ++i)
  {
    Substructure_Atom * x = _sidechains[i]->query_atom_with_atom_map_number(a2);
    if (NULL == x)
      continue;

    r2 = i;
    break;
  }

  cerr << "Atom " << a1 << " in scaffold, atom " << a2 << " in sidechain " << r2 << endl;

  if (r2 < 0)
  {
    cerr << "IWReaction::_from_smirks_add_inter_particle_bond_involves_scaffold:no sidechain for " << a2 << endl;
    return 0;
  }

  _sidechains[r2]->add_inter_particle_bond(0, a1, a2, bt);

  return 1;
}

/*
  Due to how the Substructure_Atom object is built, the attributes will be in the
  first Substructure_Atom_Specifier component. Fetch that if it is present
*/

int
IWReaction::_discern_atomic_changes(Reaction_Site & r,
                                    const Substructure_Atom & q1,
                                    const Substructure_Atom & q2)
{
  const int nc1 = q1.ncomponents();
  const int nc2 = q2.ncomponents();

  const auto a = q1.initial_atom_number();

  if (nc1 > 0 && nc2 > 0)
    return _discern_atomic_changes_specifier(r, a, *q1.component(0), *q2.component(0));
  if (nc1 > 0)
    return _discern_atomic_changes_specifier(r, a, *q1.component(0), q2);
  if (nc2 > 0)
    return _discern_atomic_changes_specifier(r, a, q1, *q2.component(0));

  return _discern_atomic_changes_specifier(r, a, q1, q2);
}

int
IWReaction::_discern_atomic_changes_specifier(Reaction_Site & r,
                                    const atom_number_t a,
                                    const Substructure_Atom_Specifier & q1,
                                    const Substructure_Atom_Specifier & q2)
{

  const auto & e2 = q2.element();
  if (e2.number_elements() > 1)
  {
    cerr << "IWReaction::_discern_atomic_changes:RHS has multiple elements\n";
    for (int i = 0; i < e2.number_elements(); ++i)
    {
      cerr << ' ' << e2[i];
    }
    cerr << endl;
    return 0;
  }

  if (e2.number_elements() > 0)
  {
    Reaction_Change_Element * rce = new Reaction_Change_Element;
    rce->set_atom(a);
    rce->set_element(e2[0]);
    r.add_element_to_be_changed(rce);
  }

  const auto & fc1 = q1.Substructure_Atom_Specifier::formal_charge();
  const auto & fc2 = q2.Substructure_Atom_Specifier::formal_charge();

  int need_to_set_formal_charge = 99;

  if (0 == fc1.number_elements() && 0 == fc2.number_elements())   // not specified either side
    ;
  else if (fc1.number_elements() > 0 && 0 == fc2.number_elements())    // set on LHS, but not on right
    need_to_set_formal_charge = 0;
  else if (fc2.number_elements() > 0)        // set on RHS, must set
    need_to_set_formal_charge = fc2[0];

  if (99 != need_to_set_formal_charge)
    r.add_formal_charge_to_assign(a, need_to_set_formal_charge);

  const auto & iso1 = q1.isotope();
  const auto & iso2 = q2.isotope();

  if (iso2.number_elements() > 1)
  {
    cerr << "IWReaction::_discern_atomic_changes:RHS has multiple isotopes\n";
    return 0;
  }

  int need_to_set_isotope = 2776704;

  if (0 == iso1.number_elements() && 0 == iso2.number_elements())
    ;
  else if (iso1.number_elements() > 0 && 0 == iso2.number_elements())
    need_to_set_isotope = 0;
  else if (iso2.number_elements())
    need_to_set_isotope = iso2[0];

  if (2776704 != need_to_set_isotope)
    r.add_isotope_to_be_placed(a, need_to_set_isotope);

//cerr << "need_to_set_isotope " << need_to_set_isotope << endl;

  const aromaticity_type_t arom1 = q1.Substructure_Atom_Specifier::aromaticity();
  const aromaticity_type_t arom2 = q2.Substructure_Atom_Specifier::aromaticity();

  if (arom1 != arom2)    // finish this sometime...
  {
  }

  return 1;
}


const Element *
Substructure_Atom::first_specified_element() const
{
  if (_element.number_elements())
    return _element[0];

  for (int i = 0; i < _components.number_elements(); i++)
  {
    const Substructure_Atom_Specifier * c = _components[i];

    if (c->element().number_elements())
      return c->element()[0];
  }

  return NULL;
}

int
Substructure_Atom::first_specified_formal_charge() const
{
  if (_formal_charge.number_elements())
    return _formal_charge[0];

  for (int i = 0; i < _components.number_elements(); i++)
  {
    const Substructure_Atom_Specifier * c = _components[i];

    const Min_Max_Specifier<int> & fc = c->formal_charge();

    if (fc.number_elements())
      return fc[0];
  }

  return 0;
}

int
Substructure_Atom::first_specified_isotope() const
{
  if (_isotope.number_elements())
    return _isotope[0];

  for (int i = 0; i < _components.number_elements(); i++)
  {
    const Substructure_Atom_Specifier * c = _components[i];

    const Min_Max_Specifier<int> & iso = c->isotope();

    if (iso.number_elements())
      return iso[0];
  }

  return 0;
}

#ifdef TESTING_STUFF

static int verbose = 0;

static void
usage(int rc)
{
  exit(rc);
}

template <typename T>
int
test_reaction_from_smirks_record(const const_IWSubstring & buffer,
                                  T & output)
{
  IWReaction rxn;

  if (! rxn.construct_from_smirks(buffer))
  {
    cerr << "Invalid smirks\n";
    return 0;
  }

  return 1;
}

template <typename T>
int
test_reaction_from_smirks(iwstring_data_source & input,
                           T & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! test_reaction_from_smirks_record(buffer, output))
    {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

template <typename T>
int
test_reaction_from_smirks(const char * fname,
                           T & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_reaction_from_smirks(input, output);
}

static int
test_reaction_from_smirks(const int argc,
                           char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:l");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }
  else 
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  return 0;
}

int
main(int argc, char **argv)
{
  int rc = test_reaction_from_smirks(argc, argv);

  return rc;
}
#endif
