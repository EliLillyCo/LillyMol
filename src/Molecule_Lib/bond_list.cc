#include <iostream>
#include <assert.h>
#include <memory>

using std::cerr;
using std::endl;

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "bond_list.h"

Bond_list::Bond_list()
{
  _magic = BOND_LIST_MAGIC;
}

Bond_list::~Bond_list()
{
  _magic = 0;
}

int
Bond_list::ok() const
{
  if (BOND_LIST_MAGIC != _magic)
    return 0;

  if (! resizable_array_p<Bond>::ok())
    return 0;

  return 1;
}

int
Bond_list::debug_print(std::ostream & os) const
{
  assert(os.good());

  os << "Bond list contains " << _number_elements  << " members\n";

  for (int i = 0; i < _number_elements; i++)
  {
    os << " Bond " << i << " ";
    const Bond *b = _things[i];
    b->debug_print(os);
  }

  return 1;
}

int
Bond_list::which_bond (atom_number_t a1, atom_number_t a2) const
{
  assert(ok());
  assert(a1 >= 0 && a2 >= 0 && a1 != a2);

  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->involves(a1, a2))
      return i;
  }

  return -1;
}

Bond_list &
Bond_list::operator = (Bond_list && rhs)
{
  resizable_array_p<Bond> & prhs = rhs;
  resizable_array_p<Bond> & plhs = *this;

  plhs = std::move(prhs);

  return *this;
}

int
Bond_list::remove_bond_between_atoms (atom_number_t a1, atom_number_t a2)
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->involves(a1, a2))
    {
      remove_item(i);
      return 1;
    }
  }

  cerr << "Bond_list::remove_bond_between_atoms: no bond between atoms " << a1 << " and " << a2 << endl;
  assert(NULL == "this should not happen");
  return 0;
}

int
Bond_list::remove_bonds_involving_these_atoms (const int * rm)
{
  int rc = 0;

  for (int i = _number_elements - 1; i >= 0; --i)
  {
    const Bond * b = _things[i];

    if (rm[b->a1()] || rm[b->a2()])
    {
      remove_item(i);
      rc++;
    }
  }

  return rc;
}

/*
  Change the bond list to reflect the fact that atoms A1 and A2
  have been swapped. 
*/

int
Bond_list::swap_atoms (atom_number_t a1, atom_number_t a2)
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    Bond * b = _things[i];
    rc += b->swap_atoms(a1, a2);
  }

  return rc;
}

int
Bond_list:: move_atom_to_end_of_atom_list (atom_number_t zatom, int atoms_in_molecule)
{
  for (int i = 0; i < _number_elements; i++)
  {
    Bond * b = _things[i];

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (a1 == zatom)
      a1 = atoms_in_molecule - 1;
    else if (a1 > zatom)
      a1--;

    if (a2 == zatom)
      a2 = atoms_in_molecule - 1;
    else if (a2 > zatom)
      a2--;

    b->set_a1a2(a1, a2);
  }

  return 1;
}

#ifdef BONDS_KNOW_RING_MEMBERSHIP


/*
  Note that the strategy of just checking the first bond does not really work
  because we are not properly invalidating the bond ring info - we are setting
  it to zero rather than unknown. Must fix this sometime...
*/

int
Bond_list::invalidate_ring_info()
{
  if (_number_elements > 0 && ! _things[0]->nrings_known())   // /work on the assumption that if the first bond does not have ring info. then none of them do.
    return 0;

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->invalidate_nrings();
  }

  return 1;
}

#endif

int
Bond_list::copy_bond_types (bond_type_t * b) const
{
  assert(nullptr != b);

  for (int i = 0; i < _number_elements; i++)
  {
    b[i] = BOND_TYPE_ONLY(_things[i]->btype());
  }

  return _number_elements;
}

/*
  Optimisation!! If the first bond doesn't have a number assigned, we assume
  that none do. Dangerous, but for efficiency
*/

void
Bond_list::invalidate_bond_numbers()
{
  if (0 == _number_elements)
    return;

  if (! _things[0]->bond_number_assigned())
    return;

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->invalidate_bond_number();
  }

  return;
}

int
Bond_list::set_all_bond_types (bond_type_t bt)
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_not_directional();   // if all bonds are being set to one type, directionality will be lost

    if (bt == _things[i]->btype())
      continue;

    _things[i]->set_bond_type(bt);
    rc++;
  }

  return rc;
}

int
Bond_list::cis_trans_bonds_present() const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->is_directional())
      return 1;
  }

  return 0;
}

int
Bond_list::remove_bonds_involving_these_atoms (const int * r,
                                               int adjust_atom_numbers)
{
  int rc = 0;

  for (auto i = _number_elements - 1; i >= 0; --i)
  {
    if (_things[i]->either_atom_set_in_array(r))    // bond must be removed
    {
      remove_item(i);
      rc++;
    }
    else if (adjust_atom_numbers)
      _things[i]->adjust_for_loss_of_atom(i);
  }

  return rc;
}

int
Bond_list::adjust_atom_numbers_for_loss_of_atom (const atom_number_t zatom)
{
  for (int i = 0; i < _number_elements; ++i)
  {
    _things[i]->adjust_for_loss_of_atom(zatom);
  }

  return 1;
}

int
Bond_list::new_atom_numbers (const int * xref)
{
  for (int i = 0; i < _number_elements; ++i)
  {
    _things[i]->new_atom_numbers(xref);
  }

  return 1;
}
