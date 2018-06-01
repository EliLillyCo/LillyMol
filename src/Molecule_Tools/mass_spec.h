#ifndef IW_FRAGMENTATION_H
#define IW_FRAGMENTATION_H 1

#include <iostream>

class Molecule;

class Fragmentation
{
  private:
    int _atoms_lost;
    int _heavy_atoms_lost;
    int _atoms_in_molecule;
    int * _lost;
    int * _rings_lost;
    int _nr;               // rings in the molecule - never changes.
    int _nr_remaining;     // rings remaining in the molecule
    int * _nrings;         // nrings for each atom - changes as molecule is destroyed
    int * _ncon;           // ncon for each atom  - never changes
    atomic_mass_t _mass_lost;
    int _max_atoms_to_lose;

    resizable_array<atomic_mass_t> _masses_lost;

    int * _is_exposed_ring_atom;

//  Private functions

    int _remaining_ncon (Molecule *, atom_number_t) const;

    int _compute_mass_and_lose (Molecule *, atom_number_t, atomic_mass_t &);
    int _lose_atom             (Molecule *, atom_number_t);
    int _lose_atoms            (Molecule *, const Set_of_Atoms &);

    int _break_longest_ring_with_exposed_atoms (Molecule *);
    int _break_longest_ring                    (Molecule *);
    int _break_rings                           (Molecule *);

    int _invalidate_ring      (Molecule *, int);
    int _strip_exposed_ring   (Molecule *, int);
    int _strip_exposed_rings  (Molecule *);

    int _identify_anchor             (Molecule *, atom_number_t);
    int _strip_exposed_hetero_group  (Molecule *, atom_number_t);
    int _strip_exposed_hetero_groups (Molecule *);

    int _strip_halogens (Molecule *);

  public:
    Fragmentation (Molecule *);
    ~Fragmentation ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int perform_fragmentation (Molecule *, std::ostream &);

    int atoms_lost () const { return _atoms_lost;}
    int heavy_atoms_lost () const { return _heavy_atoms_lost;}
    atomic_mass_t mass_lost () const { return _mass_lost;}
};

#endif
