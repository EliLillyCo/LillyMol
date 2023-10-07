#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <memory>

#include "assert.h"

#define COMPILING_CTB

#include "Foundational/iwmisc/misc.h"

#include "molecule.h"
#include "smiles.h"
#include "misc2.h"

using std::cerr;
using std::endl;

static int File_Scope_discard_directional_bonds_on_input = 0;

void
set_discard_directional_bonds_on_input (int s)
{
  File_Scope_discard_directional_bonds_on_input = s;
}

int
Molecule::cis_trans_bonds_present() const
{
  return _bond_list.cis_trans_bonds_present();
}

int
Molecule::revert_all_directional_bonds_to_non_directional()
{
  int rc = 0;

  for (int i = 0; i < _bond_list.number_elements(); i++)
  {
    Bond * b = _bond_list[i];

    if (b->is_directional())
    {
      b->set_not_directional();
      rc++;
    }
    else if (b->part_of_cis_trans_grouping())
      b->set_part_of_cis_trans_grouping(0);
  }

  return rc;
}

/*
  When a molecule is being read in without its invalid directional info, we append something
  to the name
*/

int
Molecule::_append_bad_cis_trans_input_text_to_name()
{
  _molecule_name << " (no CTB)";

  return _molecule_name.length();
}

/*
  An atom is going to be unbonded - probably as a prelude to being
  removed competely.  We need to invalidate any cis-trans
  stereochemistry that depends on that atom Note that if there is
  still one directional bonds at our end of the double bond, we can
  leave things in place
*/

int
Molecule::_atom_being_unbonded_check_directional_bonds (atom_number_t zatom,
                                                       int preserve_chirality)
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  resizable_array<Bond *> directional_bonds_found;

  Bond * part_of_cis_trans_grouping = nullptr;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_single_bond())
    {
      if (b->is_directional())
        directional_bonds_found.add(const_cast<Bond *>(b));
    }
    else if (b->is_double_bond() && b->part_of_cis_trans_grouping())
      part_of_cis_trans_grouping = const_cast<Bond *>(b);
  }

  if (nullptr != part_of_cis_trans_grouping)      // breaking up the central part of the cis-trans grouping
    return _invalidate_directional_double_bond(*part_of_cis_trans_grouping);

  if (directional_bonds_found.empty())
    return 1;

// The atom is not the centre of a cis-trans bond, but is bonded to such a bond

  for (int i = 0; i < directional_bonds_found.number_elements(); i++)
  {
    Bond * b = directional_bonds_found[i];

    atom_number_t j = b->other(zatom);          // Atom J may or may not be part of a cis-trans bond

    const Atom * aj = _things[j];

    if (! aj->is_centre_of_cis_trans_bond())   // directional single bond bonded to double bond that isn't part of a cis-trans system
      continue;

    int hj = hcount(j);

#ifdef DEBUG_ATOM_BEING_UNBONDED_CHECK_DIRECTIONAL_BONDS
    cerr << "Directional bond from " << zatom << " to " << j << " preserve? " << preserve_chirality << ", hj " << hj << endl;
#endif

//  If we are removing a non-H from something that already has an H, we
//  must invalidate

    if (1 == hj && 1 != a->atomic_number())
    {
      b->set_not_directional();
      continue;
    }

//  are we removing an explicit H from an atom that has another non-H connection.
//  That can stay intact

    if (1 == hj && 1 == a->atomic_number() && 3 == aj->ncon())
      continue;

#ifdef DEBUG_ATOM_BEING_UNBONDED_CHECK_DIRECTIONAL_BONDS
    cerr << "Invalidating\n";
#endif

//  The whole cis-trans system needs to be invalidated

    for (int k = 0; k < aj->ncon(); k++)
    {
      Bond * b = const_cast<Bond *>(aj->item(k));

      if (! b->is_double_bond())
        continue;

      _invalidate_directional_double_bond(*b);
      break;
    }
  }

  return 1;
}

/*
  The Bond B is going to disappear, invalidate any related cis-trans info
*/

int
Molecule::_invalidate_directional_double_bond (Bond & b)
{
  assert (b.is_double_bond());

  b.set_part_of_cis_trans_grouping(0);

  _invalidate_directional_bonds_at_end_of_double_bond(b.a1());
  _invalidate_directional_bonds_at_end_of_double_bond(b.a2());

  return 1;
}

int
Molecule::_invalidate_directional_bonds_at_end_of_double_bond (atom_number_t zatom)
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    Bond * b = const_cast<Bond *>(a->item(i));

    if (b->is_double_bond())
      continue;

    if (! b->is_directional())
    {
      cerr << "Molecule::_invalidate_directional_bonds_at_end_of_double_bond:strange, non directional bond " << (*b) << endl;
      continue;
    }

    atom_number_t j = b->other(zatom);

    if (_things[j]->is_centre_of_cis_trans_bond())
      continue;

    b->set_not_directional();
  }

  return 1;
}

/*
  Some possibly large changes have happened to the molecule.  Do we
  need to invalidate any directional bonds
*/

int
Molecule::_remove_directionality_from_bonds_not_actually_directional()
{
  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; i++)
  {
    Bond * b = _bond_list[i];

    if (! b->part_of_cis_trans_grouping())
      continue;

    if (_adjacent_directional_bonds_ok(*b))
      continue;

    _invalidate_directional_double_bond(*b);
  }

// Now look for any directional bonds that may be un-attached to anything

  for (int i = 0; i < nb; i++)
  {
    Bond * b = _bond_list[i];

    if (! b->is_directional())
      continue;

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (_things[a1]->is_centre_of_cis_trans_bond() ||
        _things[a2]->is_centre_of_cis_trans_bond())
      continue;

    b->set_not_directional();
  }

  return 1;
}

/*
  Accord smiles contain things like

c1(c2c([nH]c1)cc(cc2)F)\C=1\CCN(C/C1)CC\C=1\C(N/2\C(=N/C1/C)\C=C/C=C2)=O.Cl PBCHM11689904

  where directional bonds are used around rings, and may not touch
  any double bond.
*/

int
Molecule::_unset_directional_bonds_not_adjacent_to_a_double_bond()
{
  int nb = _bond_list.number_elements();

  int rc = 0;

  for (int i = 0; i < nb; i++)
  {
    Bond * b = _bond_list[i];

    if (! b->is_directional())
      continue;

    const Atom * a = _things[b->a1()];
    if (a->ncon() < a->nbonds())
      continue;

    a = _things[b->a2()];
    if (a->ncon() < a->nbonds())
      continue;

    b->set_not_directional();
    rc++;
  }

  return rc;
}

int
Molecule::_adjacent_directional_bonds_ok (const Bond & b) const
{
  assert (b.is_double_bond());

  if (! _adjacent_directional_bonds_ok(b.a1()))
    return 0;

  if (! _adjacent_directional_bonds_ok(b.a2()))
    return 0;

  return 1;
}

int
Molecule::_adjacent_directional_bonds_ok (atom_number_t zatom) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  if (acon < 2)
    return 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_double_bond())
      continue;

    if (! b->is_directional())
      return 0;
  }

  return 1;
}

/*
  quick check on whether or not a given atom can be one of the
  doubly bonded atoms in a cis-trans bond
*/

int
Molecule::_can_be_end_of_directional_bond (atom_number_t zatom) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  if (acon < 2)
    return 0;

  if (acon > 3)
    return 0;

  if (a->nbonds() != acon + 1)
    return 0;

  return 1;
}

/*
  We have just read a smiles.

  We need to make sure all the directional bonds are OK and assign any double bonds
  to part_of_cis_trans_grouping()
  
  If there are bonds like

  C\C(N)=C(/F)Cl 

  fill in the unspecified, but unambiguous directional bonds

*/

int
Molecule::__finished_reading_smiles_assign_and_check_directional_bonds()
{
//_unset_directional_bonds_not_adjacent_to_a_double_bond();  seems to be not needed, watch to see if problems

  int nb = _bond_list.number_elements();

// First identify those bonds that are initially specified as directional

  resizable_array<const Bond *> directional_bonds;

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];

    if (b->is_directional())
      directional_bonds.add(b);
  }

  int nd = directional_bonds.number_elements();

  if (0 == nd)
    return 1;    // nothing to check

// Things like F\C=C(/C)F

  _set_any_unambiguous_unset_directional_bonds(directional_bonds);

  nd = directional_bonds.number_elements();

  if (1 == nd)
  {
#ifdef DEBUG_CIS_TRANS_BOND_X
    if (display_smiles_interpretation_error_messages())
      cerr << "Molecule::_finished_reading_smiles_assign_and_check_directional_bonds:only one directional bond, impossible\n";
#endif

    return 0;
  }

// Now make sure that there are no unmatched directional bonds

  resizable_array<const Bond *> coupled;

  for (int i = 0; i < nd; i++)
  {
    const Bond * b = directional_bonds[i];

    if (! _identify_directional_bonds_across_double_bonds(b, coupled))
    {
#ifdef DEBUG_CIS_TRANS_BOND_X
      if (display_smiles_interpretation_error_messages())
        cerr << "Molecule::_finished_reading_smiles_assign_and_check_directional_bonds:mismatched directional bond " << b->a1() << " " << b->a2() << endl;
#endif
      return 0;
    }
  }

  if (coupled.number_elements() != nd)
  {
#ifdef DEBUG_CIS_TRANS_BOND_X
    if (display_smiles_interpretation_error_messages())
      cerr << "Molecule::_finished_reading_smiles_assign_and_check_directional_bonds:incomplete cis trans bonds\n";
#endif
    return 0;
  }

// Scan those bonds initially specified as directional

#ifdef DEBUG_CIS_TRANS_BOND_X
  cerr << "Scanning " << nd << " directional bonds\n";
#endif

  for (int i = 0; i < nd; i++)
  {
    const Bond * b = directional_bonds[i];

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (in_same_ring(a1, a2))
      continue;

    int valid_cis_trans_form_found = 0;

    if (! _process_directional_system(a1, a2, valid_cis_trans_form_found))
      return 0;

    if (valid_cis_trans_form_found)
      ;
    else if (! _process_directional_system(a2, a1, valid_cis_trans_form_found))
      return 0;

    if (valid_cis_trans_form_found)
      continue;

//#define DEBUG_CHECK_INCONCISTENT_CIS_TRANS
#ifdef DEBUG_CHECK_INCONCISTENT_CIS_TRANS
    if (display_smiles_interpretation_error_messages())
      cerr << "Molecule::_finished_reading_smiles_assign_and_check_directional_bonds:invalid or incomplete cis trans bond, atoms " << a1 << " and " << a2 << endl;
#endif

    return 0;
  }

  return 1;
}

/*
  I have had so many problems with directional bonds that there are
  plenty of invalid molecules out there, so be permissive
*/

int
Molecule::_finished_reading_smiles_assign_and_check_directional_bonds()
{
  if (File_Scope_discard_directional_bonds_on_input)
  {
    revert_all_directional_bonds_to_non_directional();
    return 1;
  }

  if (__finished_reading_smiles_assign_and_check_directional_bonds())
    return 1;

  revert_all_directional_bonds_to_non_directional();
  return 1;
}

static int
discern_bond_directionality (const Bond * b,
                             atom_number_t a)
{
  if (! b->is_directional())
    return 0;

  if (a == b->a1())
  {
    if (b->is_directional_up())
      return  1;
    else
      return -1;
  }
  else if (a == b->a2())
  {
    if (b->is_directional_up())
      return -1;
    else
      return  1;
  }
  else
    return 0;
}

/*
  We have a potential cis-trans bonding system

*/

int
Molecule::_process_directional_system (atom_number_t lhs1,
                                       atom_number_t db1,
                                       int & valid_cis_trans_form_found)
{
//cerr << "Examining possible cis-trans system involving atoms " << lhs1 << " and " << db1 << endl;

  const Atom * a = _things[db1];

  int acon = a->ncon();

  if (2 == acon)     // only 2 or 3 connections possible
    ;
  else if (3 == acon)
    ;
  else
    return 1;      // return code zero is reserved for failures

  if (acon + 1 != a->nbonds())   // the doubly bonded item must be unsaturated
    return 1;

  const Bond * lhsb1 = nullptr;
  const Bond * lhsb2 = nullptr;
  const Bond * double_bond = nullptr;
  atom_number_t db2 = INVALID_ATOM_NUMBER;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t a2 = b->other(db1);

    if (b->is_double_bond())
    {
      db2 = a2;
      double_bond = b;
    }
    else if (a2 == lhs1)
      lhsb1 = b;
    else
      lhsb2 = b;
  }

  assert (nullptr != lhsb1);
  assert (lhsb1->is_directional());
  assert (INVALID_ATOM_NUMBER != db2);

  a = _things[db2];

  acon = a->ncon();

  if (2 == acon)
    ;
  else if (3 == acon)
    ;
  else
    return 1;   // not a valid cis-trans system here

  const Bond * rhsb1 = nullptr;
  const Bond * rhsb2 = nullptr;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_double_bond())
      continue;

    if (nullptr == rhsb1)
      rhsb1 = b;
    else
    {
      rhsb2 = b;
      break;
    }
  }

  if (nullptr == rhsb1)    // no valid directional bond here
    return 1;

  if (rhsb1->is_directional())   // good
    ;
  else if (nullptr == rhsb2)        // not a valid directional bond here
    return 1;
  else if (rhsb2->is_directional())   // ensure rhsb1 is directional
    std::swap(rhsb1, rhsb2);
  else                           // not a valid directional bond here
    return 1;

// Now we need to see if the bonding is complete

  int lhs1_direction = discern_bond_directionality(lhsb1, db1);

  if (nullptr == lhsb2)   // no second bond on LHS
    ;
  else if (lhsb2->is_directional())   // already marked, make sure consistent
  {
    int lhsb2_directionality = discern_bond_directionality(lhsb2, db1);
    if (lhsb2_directionality == lhs1_direction)
    {
      if (display_smiles_interpretation_error_messages())
        cerr << "Molecule::_process_directional_system:inconsistent directional bonding on bond " << (*lhsb2) << endl;
      return 0;
    }
  }
  else if (lhs1_direction < 0)
    const_cast<Bond *>(lhsb2)->set_directional_up (db1, lhsb2->other(db1));
  else
    const_cast<Bond *>(lhsb2)->set_directional_down(db1, lhsb2->other(db1));

// Now that we have straightened out lhsb2, look at the right hand side

  int rhsb1_direction = discern_bond_directionality(rhsb1, db2);

  if (nullptr == rhsb2)
    ;
  else if (rhsb2->is_directional())
  {
    int rhsb2_directionality = discern_bond_directionality(rhsb2, db2);
    if (rhsb2_directionality == rhsb1_direction)
    {
      if (display_smiles_interpretation_error_messages())
      {
        cerr << "Molecule::_process_directional_system:inconsistent directional bonding on bonds\n";
        cerr << (*rhsb2) << endl;
        cerr << (*rhsb1) << endl;
      }
      return 0;
    }
  }
  else if (rhsb1_direction < 0)
    const_cast<Bond *>(rhsb2)->set_directional_up (db2, rhsb2->other(db2));
  else if (rhsb1_direction > 0)
    const_cast<Bond *>(rhsb2)->set_directional_down(db2, rhsb2->other(db2));

  valid_cis_trans_form_found = 1;

  const_cast<Bond *>(double_bond)->set_part_of_cis_trans_grouping(1);

  return 1;
}

/*
  We have found a directional bond B.
  We need to identify other directional bonds that are at the end of
  a double bond from it.
*/

int
Molecule::_identify_directional_bonds_across_double_bonds (const Bond * b,
                                resizable_array<const Bond *> & coupled) const
{
  return _identify_directional_bonds_across_double_bonds(b->a1(), coupled) +
         _identify_directional_bonds_across_double_bonds(b->a2(), coupled);
}

int
Molecule::_identify_directional_bonds_across_double_bonds (atom_number_t zatom,
                                                           resizable_array<const Bond *> & coupled) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  if (acon == a->nbonds())   // no double bonds here
    return 0;

  const Bond * double_bond = nullptr;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_double_bond())
    {
      double_bond = b;
      break;
    }
  }

  if (nullptr == double_bond)   // maybe someone drew a directional bond to a triple bond
    return 0;

  atom_number_t db2 = double_bond->other(zatom);

  a = _things[db2];
  acon = a->ncon();

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);
    if (b->is_double_bond())
      continue;

    if (! b->is_directional())
      continue;

    coupled.add_if_not_already_present(b);
    rc++;
  }

  return rc;
}

static int
directional_bond_attached (const Atom * a)
{
  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_directional())
      return 1;
  }

  return 0;
}

/*
  This is complicated because we could get inconsistent specifications.
  We ignore that...
*/

int
Molecule::_set_any_unambiguous_unset_directional_bonds(resizable_array<const Bond *> & directional_bonds)
{
  int ne = _bond_list.number_elements();

  for (int i = 0; i < ne; i++)
  {
    const Bond * b = _bond_list[i];

    if (! b->is_double_bond())
      continue;

    atom_number_t a1 = b->a1();

    if (! directional_bond_attached(_things[a1]))
      continue;

    atom_number_t a2 = b->a2();

    if (! directional_bond_attached(_things[a2]))
      continue;

    _fill_in_missing_directional_bond_specification(a1, directional_bonds);
    _fill_in_missing_directional_bond_specification(a2, directional_bonds);
  }

  return 1;
}

/*
  Fix any unspecified directional bonds
*/

int
Molecule::_fill_in_missing_directional_bond_specification (atom_number_t zatom,
                                        resizable_array<const Bond *> & directional_bonds)
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  if (3 != acon)   // the only case we can process
    return 0;

  const Bond * directional = nullptr;
  Bond * not_directional = nullptr;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_double_bond())
      continue;

    if (b->is_directional())
    {
      if (nullptr != directional)   // already got a directional bond, must be OK
        return 1;

      directional = b;
    }
    else
    {
      if (nullptr != not_directional)   // nothing directional here
        return 1;

      not_directional = const_cast<Bond *>(b);
    }
  }

  if (nullptr == directional)
    return 0;

  if (nullptr == not_directional)
    return 0;

// The not_directional bond must now be set

  if (directional->is_directional_up())
  {
    if (zatom == directional->a1())
    {
      if (zatom == not_directional->a1())
        not_directional->set_directional_down();
      else
        not_directional->set_directional_up();
    }
    else
    {
      if (zatom == not_directional->a1())
        not_directional->set_directional_up();
      else
        not_directional->set_directional_down();
    }
  }
  else    // is a down bond
  {
    if (zatom == directional->a1())
    {
      if (zatom == not_directional->a1())
        not_directional->set_directional_up();
      else
        not_directional->set_directional_down();
    }
    else
    {
      if (zatom == not_directional->a1())
        not_directional->set_directional_down();
      else
        not_directional->set_directional_up();
    }
  }

#ifdef DEBUG_CHECK_INCONCISTENT_CIS_TRANS
  if (not_directional->is_directional_up())
    cerr << "Non directional bond " << not_directional->a1() << ' ' << not_directional->a2() << " now up\n";
  else if (not_directional->is_directional_down())
    cerr << "Non directional bond " << not_directional->a1() << ' ' << not_directional->a2() << " now down\n";
  else
    cerr << "HUH\n";
#endif

  directional_bonds.add(not_directional);

  return 1;
}

/*
  An atom has one or more directional bonds attached, and is unsaturated.
  Identify the doubly bonded other atom and mark the bond as being part of
  a cis-trans system
*/

/*int
Molecule::_mark_adjacent_double_bond_with_directional_atoms_at_other_end (atom_number_t a3,
                                        int * process_these_atoms)
{
  const Atom * a = _things[a3];

  int acon = a->ncon();

  Bond * doubly_bonded = nullptr;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_double_bond())
    {
      doubly_bonded = const_cast<Bond *>(b);
      break;
    }
  }

  if (nullptr == doubly_bonded)     // should be impossible
  {
    cerr << "Molecule::_mark_adjacent_double_bond_with_directional_atoms_at_other_end:no doubly bonded connection\n";
    return 0;
  }

#ifdef DEBUG_MARK_ADJACENT_DOUBLE_BOND_WITH
  cerr << "Doubly bonded bond from " <<a3 << " is " << doubly_bonded->a1() << " to " << doubly_bonded->a2() << endl;
#endif

  if (doubly_bonded->part_of_cis_trans_grouping())    // already processed
    return 1;

  atom_number_t a4 = doubly_bonded->other(a3);

  if (0 == _things[a4]->number_directional_bonds_attached())
    return 0;

  if (in_same_ring(a3, a4))
    return 0;

  doubly_bonded->set_part_of_cis_trans_grouping(1);

  process_these_atoms[a4] = 1;

  return 1;
}*/

/*int
Molecule::_adjacent_directional_bonds_mutually_consistent (atom_number_t zatom)
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  if (acon < 3)    // must be consistent
    return 1;

  int up_bond = 0;
  int down_bond = 0;

  Bond * non_directional_bond_encountered = nullptr;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_double_bond())
      continue;

    if (b->is_directional_up())
    {
      if (zatom == b->a1())
        up_bond++;
      else
        down_bond++;
    }
    else if (b->is_directional_down())
    {
      if (zatom == b->a1())
        down_bond++;
      else
        up_bond++;
    }
    else
      non_directional_bond_encountered = const_cast<Bond *> (b);
  }

//cerr << "At atom " << zatom << " found " << up_bond << " up bonds and " << down_bond << " down bonds\n";

  if (1 == up_bond && 1 == down_bond)
    return 1;

  if (0 == up_bond && 0 == down_bond)   // how could that happen?
    return 0;

  if (nullptr == non_directional_bond_encountered)    // seems unlikely
    return 0;

//cerr << "Non directional bond from " << zatom << " is " << non_directional_bond_encountered->a1() << " to " << non_directional_bond_encountered->a2() << endl;

  if (1 == up_bond && 0 == down_bond)     // need a down bond
  {
    if (zatom == non_directional_bond_encountered->a1())
      non_directional_bond_encountered->set_directional_down();
    else
      non_directional_bond_encountered->set_directional_up();
  }
  else if (0 == up_bond && 1 == down_bond)    // need an up bond
  {
    if (zatom == non_directional_bond_encountered->a1())
      non_directional_bond_encountered->set_directional_up();
    else
      non_directional_bond_encountered->set_directional_down();
  }
  else
    return 0;

  return 1;
}*/

/*
  ZATOM is at the end of a double bond.
  Identify all the directional bonds attached, and if possible,
  the unprocessed atom at the other end of the double bond
*/

//#define DEBUG_IDENTIFY_LINKED_CIS_TRANS_BONDS 1

int
Molecule::_identify_linked_cis_trans_bonds(resizable_array<const Bond *> & bonds_to_be_flipped,
                                           atom_number_t previous_atom,
                                           atom_number_t zatom,
                                           int * bond_already_done) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

#ifdef DEBUG_IDENTIFY_LINKED_CIS_TRANS_BONDS
  cerr << "Atom " << zatom << " is attached to directional bond, ncon " << acon << endl;
#endif

  if (acon < 2)
    return 0;

  atom_number_t atom_at_end_of_double_bond = INVALID_ATOM_NUMBER;

  int double_bonds_encountered = 0;

  Set_of_Atoms coupled_atoms;   // unsaturated atoms linked

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);

    if (j == previous_atom)
      continue;

    if (b->is_double_bond())
    {
      atom_at_end_of_double_bond = j;
      double_bonds_encountered++;
      continue;
    }

    if (! b->is_directional())
      continue;

    int bn = b->bond_number();

    if (bond_already_done[bn])
      continue;

    bond_already_done[bn] = 1;
    bonds_to_be_flipped.add(b);

    const Atom * aj = _things[j];

    if (aj->ncon() + 1 == aj->nbonds())
      coupled_atoms.add(j);
  }

#ifdef DEBUG_IDENTIFY_LINKED_CIS_TRANS_BONDS
  cerr << "Found " << coupled_atoms.number_elements() << " possibly coupled atoms, double bond " << atom_at_end_of_double_bond << endl;
#endif

  for (int i = 0; i < coupled_atoms.number_elements(); i++)
  {
    atom_number_t j = coupled_atoms[i];

    _identify_linked_cis_trans_bonds(bonds_to_be_flipped, zatom, j, bond_already_done);
  }

  if (INVALID_ATOM_NUMBER == atom_at_end_of_double_bond || double_bonds_encountered > 1)
    return 1;

  return _identify_linked_cis_trans_bonds(bonds_to_be_flipped, zatom, atom_at_end_of_double_bond, bond_already_done);
}

int
Molecule::_identify_linked_cis_trans_bonds(resizable_array<const Bond *> & bonds_to_be_flipped,
                                           const Bond * current_bond,
                                           int * bond_already_done) const
{
  bonds_to_be_flipped.add(current_bond);
  bond_already_done[current_bond->bond_number()] = 1;

  atom_number_t a1 = current_bond->a1();
  atom_number_t a2 = current_bond->a2();

#ifdef DEBUG_IDENTIFY_LINKED_CIS_TRANS_BONDS
  cerr << "Starting with directional bond " << a1 << " to " << a2 << endl;
#endif

  const Atom * aa1 = _things[a1];
  const Atom * aa2 = _things[a2];

// One end must be a double bond

  if (aa1->ncon() + 1 == aa1->nbonds())
    _identify_linked_cis_trans_bonds(bonds_to_be_flipped, a2, a1, bond_already_done);

  if (aa2->ncon() + 1 == aa2->nbonds())
    _identify_linked_cis_trans_bonds(bonds_to_be_flipped, a1, a2, bond_already_done);

  return bonds_to_be_flipped.number_elements();
}

int
Molecule::flip_cis_trans_bond_and_all_related_directional_bonds(const Bond * b)
{
  if (! b->is_directional())
  {
    cerr << "Molecule::flip_cis_trans_bond_and_all_related_directional_bonds:bond is not directional\n";
    cerr << "Atoms " << b->a1() << " and " << b->a2() << endl;
    return 0;
  }

  _bond_list.assign_bond_numbers_to_bonds_if_needed();

  int * bond_already_done = new_int(_bond_list.number_elements()); std::unique_ptr<int[]> free_bond_already_done(bond_already_done);

  resizable_array<const Bond *> bonds_to_be_flipped;

  _identify_linked_cis_trans_bonds(bonds_to_be_flipped, b, bond_already_done);

  int n = bonds_to_be_flipped.number_elements();

  for (int i = 0; i < n; i++)
  {
    Bond * b = const_cast<Bond *>(bonds_to_be_flipped[i]);

    if (b->is_directional_up())
      b->set_directional_down();
    else if (b->is_directional_down())
      b->set_directional_up();
    else
      cerr << "Molecule::flip_cis_trans_bond_and_all_related_directional_bonds:non directional bond\n";
  }

  return 1;
}

int
Molecule::_adjust_cis_trans_bonds_to_canonical_form(const int * canonical_rank)
{
  assign_bond_numbers_to_bonds_if_needed();

  int nb = _bond_list.number_elements();

  if (0 == nb)
    return 1;

  int * bond_already_done = new_int(nb); std::unique_ptr<int[]> free_bond_already_done(bond_already_done);

  for (int i = 0; i < nb; i++)
  {
    if (bond_already_done[i])
      continue;

    const Bond * bi = _bond_list[i];

    if (! bi->is_directional())
      continue;

    resizable_array<const Bond *> bonds_in_grouping;

    _identify_linked_cis_trans_bonds(bonds_in_grouping, bi, bond_already_done);

    _canonicalise_linked_group_of_cis_trans_bonds(bonds_in_grouping, canonical_rank);
  }

  return 1;
}

/*
  Compute a canonical number associated with a bond.
  Arbitrary computation, just so long as it will be unique
*/

static int
compute_bond_canonical_rank(const Bond * b,
                            int matoms,
                            const int * canonical_rank)
{
  atom_number_t a1 = b->a1();
  atom_number_t a2 = b->a2();

  int c1 = canonical_rank[a1];
  int c2 = canonical_rank[a2];

  if (c1 > c2)
    std::swap(c1, c2);

  return c1 * (matoms + matoms) + c2;
}

int
Molecule::_canonicalise_linked_group_of_cis_trans_bonds(const resizable_array<const Bond *> &bonds_in_grouping,
                                        const int * canonical_rank)
{
  int n = bonds_in_grouping.number_elements();

  int highest_score = compute_bond_canonical_rank(bonds_in_grouping[0], _number_elements, canonical_rank);
  int ndx_of_highest_score = 0;

  for (int i = 1; i < n; i++)
  {
    int s = compute_bond_canonical_rank(bonds_in_grouping[i], _number_elements, canonical_rank);
    if (s > highest_score)
    {
      highest_score = s;
      ndx_of_highest_score = i;
    }
  }

  const Bond * b = bonds_in_grouping[ndx_of_highest_score];

//cerr << "Thereare " << n << " bonds in the cis trans grouping\n";

  int c1 = canonical_rank[b->a1()];
  int c2 = canonical_rank[b->a2()];

  if (c1 < c2 && b->is_directional_down())
    ;
  else if (c1 > c2 && b->is_directional_up())
    ;
  else
  {
//  cerr << "FLIPPING\n";
    for (int i = 0; i < n; i++)
    {
      Bond * b = const_cast<Bond *>(bonds_in_grouping[i]);

      if (b->is_directional_up())
        b->set_directional_down();
      else if (b->is_directional_down())
        b->set_directional_up();
    }
  }

  return 1;
}

/*
  If a single bond has lost its directionality, we must remove any cis-trans
  bonding information that depends on it.
*/

int
Molecule::_bond_is_no_longer_directional (const Bond * b)
{
  const Bond * db = _identify_double_bond(b->a1());

  if (nullptr != db)
    _remove_directional_bonding_associated_with_bond(db);

  db = _identify_double_bond(b->a2());

  if (nullptr != db)
    _remove_directional_bonding_associated_with_bond(db);

  return 1;
}

/*
  B describes a bond that may have had directional bonds associated with it.
  That bond will no longer have any directionality, so remove any
  directionality data associated with it
*/

int
Molecule::_remove_directional_bonding_associated_with_bond (const Bond * b)
{
  return 1;
}

const Bond *
Molecule::_identify_double_bond (atom_number_t zatom) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  if (acon == a->nbonds())
    return nullptr;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_double_bond())
      return b;
  }

  return nullptr;
}

int
Molecule::_set_bond_directionality(atom_number_t a1,
                                   atom_number_t a2,
                                   int dir)
{
  _set_modified();

  Bond * b = const_cast<Bond *>(_things[a1]->bond_to_atom(a2));

  assert (nullptr != b);
  if (b == nullptr)
  {
    cerr << "Molecule::could not find bond between atoms " << a1 << " and " << a1   << endl;
    return 0;
  }	
  if (IW_BOND_DIRECTIONAL_UP == dir)
    b->set_directional_up();
  else if (IW_BOND_DIRECTIONAL_DOWN == dir)
    b->set_directional_down();
  else
  {
    cerr << "Molecule::_set_bond_directionality:not sure what to do with " << dir << endl;
    return 0;
  }

  return 1;
}
