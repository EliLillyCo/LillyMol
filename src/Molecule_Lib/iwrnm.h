#ifndef IW_RING_NUMBER_MANAGER_H
#define IW_RING_NUMBER_MANAGER_H

/*
  This class is used when constructing a smiles.
  When placing a ring opening digit, we must keep track of
    (a) The atom which completes the ring number
    (b) The type of the bond
    (c) The atom which created the ring number
*/

#include "chiral_centre.h"

class Ring_Number_Manager
{
  private:
    int _nr;
    int * _ring_id;
//  bond_type_t * _bt;
    const Bond ** _bond;
    atom_number_t * _from_atom;

    int _include_aromaticity_in_smiles;    // initialised during constructor and never changed
    int _include_cis_trans_in_smiles;      // initialised during constructor and never changed

//  private functions

    void _append_ring_closure_digits (IWString & smiles,
                            int ring_closure_number,
                            const Bond * b,
                            atom_number_t ato) const;

    int _process_ring (IWString & smiles, int ring, atom_number_t afrom);

    int _place_ring_closure (IWString & smiles,
                   atom_number_t a,
                   atom_number_t afrom);
    int _append_ring_closures_for_chiral_atom (IWString & smiles,
                    atom_number_t a,
                    const resizable_array<atom_number_t> & ring_closures_found);
    int _generate_ring_opening_chars (IWString & ring_opening_chars,
                                                  const resizable_array<const Bond *> & ring_opening_bonds,
                                                  atom_number_t ato);

    void _default_values ();

  public:
    Ring_Number_Manager ();
    Ring_Number_Manager (int);
    ~Ring_Number_Manager ();

    int debug_print (std::ostream &) const;

    int ok () const;

    int activate (int);

//  int store_ring (IWString &, const Bond *, atom_number_t);
    int append_ring_closing_and_opening_digits (IWString & smiles, atom_number_t zatom, const resizable_array<atom_number_t> & ring_closures_found, const resizable_array<const Bond *> & ring_opening_bonds, const Chiral_Centre *);

    int append_ring_closures_for_atom (IWString &,
                atom_number_t, 
                const resizable_array<atom_number_t> &,
                const Chiral_Centre *);
};

#endif
