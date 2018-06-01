#ifndef BOND_LIST_H
#define BOND_LIST_H

#include "iwaray.h"

#include "bond.h"

class Beep;

#define BOND_LIST_MAGIC -134211

class Bond_list : public resizable_array_p<Bond>
{
  private:
    magic_number_t _magic;

// private functions

    int _maximum_connectivity (int *, int) const;

  public:
    Bond_list ();
    ~Bond_list ();

    int ok () const;

    int debug_print (std::ostream &) const;

    Bond_list & operator = (Bond_list && rhs);

    int nbonds () const { return _number_elements;}

    int     which_bond (atom_number_t, atom_number_t) const;

    int     remove_bonds_to_atom (atom_number_t, int = 0);
    int     remove_bond_between_atoms (atom_number_t, atom_number_t);
    int     remove_bonds_involving_these_atoms (const int * rm);

    Bond *  bond_between_atoms (atom_number_t, atom_number_t) const;

    int     swap_atoms   (int, int);
    int     move_atom_to_end_of_atom_list (atom_number_t, int);

//  In case someone wants rapid access to the bond types

    int     copy_bond_types (bond_type_t *) const;

    int     set_modified ();

#ifdef BONDS_KNOW_RING_MEMBERSHIP
    int     invalidate_ring_info ();
    int     assign_ring_membership_to_bonds (const resizable_array_p<Beep> & beeps);
#endif

    int     assign_bond_numbers (int istart);
    int     assign_bond_numbers_to_bonds_if_needed();
    void    invalidate_bond_numbers ();

    int     set_all_bond_types (bond_type_t);

    int     unset_all_permanent_aromatic_bonds ();

    int     cis_trans_bonds_present() const;

    int     remove_bonds_involving_these_atoms (const int * r, int adjust_atom_numbers);

    int     adjust_atom_numbers_for_loss_of_atom (const atom_number_t zatom);

    int     new_atom_numbers (const int * xref);   // makeing a subset, adjust the bonds to reflect the new numbering
};

#define OK_BOND_LIST(b) ( NULL != (b) && b->ok () )

#endif
