#ifndef TMP_DETACH_ATOMS
#define TMP_DETACH_ATOMS

/*
  When dealing with explicit hydrogens, both the donor_acceptor and
  charge_assigner objects may need to temporarily remove explicit
  Hydrogens in order for the queries to work
*/

#include "molecule.h"

class Temp_Detach_Atoms
{
  private:

    int _active;

    int _verbose;

// When detaching-reattaching explicit Hydrogens, by default we remove
// any Hydrogen no longer needed (on a Carboxyllic acid for example). 
// We ran into cases where we needed to preserve the number of atoms
// in the molecule. We break the bond, but just leave the atom.

    int _remove_hydrogens_no_longer_needed;

//  As a small consistency check, we record the number of atoms in
//  the molecule we most recently broke apart. Woe to anyone who
//  passes different molecules to different calles to detach_atoms and reattach_atoms!!

    int _matoms;

//  if we didn't detach any atoms, we don't need to reattach any

    int _need_to_reattach;

    int * _connection;

//  If the atom helps define a chiral centre, then we must make a copy of that chiral centre
//  so we can re-create it when re-attaching. No need to worry with non isotopic Hydrogen atoms

    resizable_array<Chiral_Centre *> _chiral_centre;

    bond_type_t _bt;

  public:
    Temp_Detach_Atoms ();
    ~Temp_Detach_Atoms ();

    void set_verbose (int v) { _verbose = v;}

//  The object recognised directives like "nodetach" and "noremove"

    int recognised_directive (const const_IWSubstring &);

//  We may no longer want a specific atom re-attached

    void do_not_reattach_to_atom (atom_number_t);

    int active () const { return _active;}
    int natoms () const { return _matoms;}

    int detach_atoms (Molecule &, atomic_number_t z = 1);
    int reattach_atoms (Molecule &);
};

#endif
