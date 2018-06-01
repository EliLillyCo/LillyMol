#ifndef IW_RING_CLOSURE_BONDS_H
#define IW_RING_CLOSURE_BONDS_H

#include "bond.h"
#include "iwaray.h"

/*
  We need to keep track of whether or not given pairs of atoms are
  stored. Some trickery with the _present array does that if both
  atoms are hit just once. Otherwise, we store integers
  (a1 * _atoms_in_molecule + a2) in the resizable_array<int>
*/

class Ring_Closure_Bonds : public resizable_array<int>
{
  private:
    int _atoms_in_molecule;

    int * _present;

//  private functions

    int _form_corresponding_integer (atom_number_t a1, atom_number_t a2) const;

  public:
    Ring_Closure_Bonds ();
    Ring_Closure_Bonds (const Ring_Closure_Bonds &);
    ~Ring_Closure_Bonds ();

    Ring_Closure_Bonds & operator= (const Ring_Closure_Bonds &);

    int ok () const;

    int write_bonds (std::ostream & output) const;

    int reset ();

    void invalidate ();

    int activate (int);

    int add (atom_number_t, atom_number_t);
    int contains (atom_number_t, atom_number_t) const;

//  I didn't call this operator== because it the order of the integers in the resizable_array<int> may be different

    int is_the_same (const Ring_Closure_Bonds &) const;

    int report_differences (const Ring_Closure_Bonds &, std::ostream &) const;

    int is_subset_of (const Ring_Closure_Bonds &) const;
};

#endif
