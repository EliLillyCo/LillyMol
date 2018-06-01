#ifndef IW_PATH_H
#define IW_PATH_H 1

#include <iostream>

#include "molecule.h"

class Path : public resizable_array<atom_number_t>
{
  private:
  public:
};

class Ring_Atom_Iterator;

class Ring : public Set_of_Atoms
{
  private:
    int                _is_fused;

    int                _fused_system_identifier;

//  Often had a need to keep track of which ring was which, so in Aug 97
//  I added a unique identifier for each ring

    int                _ring_number;

    aromaticity_type_t _aromaticity;

    resizable_array<Ring *> _fused_neighbours;

//  Aug 2000. When computing complexity descriptors it is handy to
//  know the highest degree of fusion between this ring and another
    
    int _largest_number_of_bonds_shared_with_another_ring;

//  The number of neighbours that are strongly fused

    int _number_strongly_fused_neighbours;

//  Nov 2000. Handy to know which fragment each ring is in

    int _fragment_membership;

//  private functions

    int _compute_bonds_shared_with (const Ring & rhs, int my_direction,
                                  int lhs_ndx,
                                  int rhs_ndx) const;
  friend
    std::ostream & operator << (std::ostream &, const Ring &);
  friend
    std::ostream & operator << (std::ostream & os, const Ring * r) { return os << (*r);}

  public:
    Ring ();
    Ring (const Ring &);

    int ok () const;

    int ring_number () const { return _ring_number;};
    void set_ring_number (int n) { _ring_number = n;};

    int fragment_membership () const { return _fragment_membership;};
    void set_fragment_membership (int f) { _fragment_membership = f;};

    int  fused_system_identifier () const { return _fused_system_identifier;}
    void set_fused_system_identifier (int i) { _fused_system_identifier = i;}    // not intended for public use
    int  propagate_fused_system_identifier (int);                                // not intended for public use

    int  is_fused () const { return _is_fused;}
    void set_is_fused (int s) { _is_fused = s;}                                  // not intended for public use

    int set_aromaticity (aromaticity_type_t);
    int set_aromaticity_to_not_determined () {return _aromaticity = AROMATICITY_NOT_DETERMINED;};

    int is_aromatic () const { return AROMATIC == _aromaticity;}
    int is_non_aromatic () const { return NOT_AROMATIC == _aromaticity;}
    int undetermined_aromaticity () const
          { return AROMATICITY_NOT_DETERMINED == _aromaticity;}

    void invalidate_aromaticity() { _aromaticity = AROMATICITY_NOT_DETERMINED;}

    int  fused_ring_neighbours () const { return _fused_neighbours.number_elements ();}
    const Ring * fused_neighbour (int) const;
    int  set_fused_to (Ring * r, int);
    int  is_fused_to (const Ring * r) const { return _fused_neighbours.contains ((Ring *) r);}   // loss of const OK

    int largest_number_of_bonds_shared_with_another_ring () const { return _largest_number_of_bonds_shared_with_another_ring;}
    int strongly_fused_ring_neighbours () const { return _number_strongly_fused_neighbours;}

    int contains_bond (atom_number_t, atom_number_t) const;
    int contains_both (atom_number_t, atom_number_t) const;

//  update_ring_membership is only used internally (pearlman.cc)

    int update_ring_membership (int *, int = 1, int = 0) const;

//  Used when finding raw rings, probably not useful in general use

    int  spiro_fused (const int * ring_membership) const;
    int  fused_ring_check_for_spiro_fusion (const int *) const;

    int  end () const { return _number_elements;}    // used by the iterator

// In aromatic.cc I need a means of computing whether an SSSR ring is
// fused to a non-sssr ring.  The non-sssr ring will not be in the
// fused_neighbour list

    int compute_bonds_shared_with (const Ring &) const;

//  We are moving a ring from the non_sssr set. Rings know their fused neighbours
//  so each ring that might have a reference to the old ring needs to be updated

    int ring_moving_from_non_sssr_to_sssr (Ring * rfrom, Ring * rto);

    Ring_Atom_Iterator find (atom_number_t) const;
};

class Ring_Bond_Iterator
{
  private:
    int _atoms_in_ring;
    int _ndx;
    const atom_number_t * _atom_number;    // extracted from the Ring during constructor

  public:
    Ring_Bond_Iterator (const Ring &);

    atom_number_t a1 () const { return _atom_number[_ndx];}
    atom_number_t a2 () const;

    atom_number_t previous_atom () const;
    atom_number_t next_atom () const;

    void set_index (int s) { _ndx = s;}   // no checking, no nothing!!

    void operator++ (int notused) { _ndx++;}
    void operator++ () { _ndx++;}
    void operator-- (int notused);
    bool operator != (int n) const { return _ndx != n;}
    bool operator == (int n) const { return _ndx == n;}
};

class Ring_Atom_Iterator
{
  private:
    int _atoms_in_ring;
    int _ndx;
    const atom_number_t * _atom_number;   // extracted from the Ring during constructor

  public:
    Ring_Atom_Iterator (const Ring &);

    void set_index (int s) { _ndx = s;}   // no checking, no nothing!!
    int  ndx () const { return _ndx;}    // just for debugging

    atom_number_t prev () const;
    atom_number_t current () const { return _atom_number[_ndx];}
    atom_number_t next () const;     // will point to Ring[0] even if iterator is at the end

    void operator++  (int notused) { _ndx++;}
    void operator--  (int notused) { _ndx--;}    // no checking
    bool operator != (int n) const { return _ndx != n;}
    bool operator == (int n) const { return _ndx == n;}

//  whereas the operators ++ and -- just increment and decrement, these two
//  can be used to deliberately move in directions. In each case, they return
//  the new atom

    atom_number_t move_forward ();
    atom_number_t move_backward ();

//  A frequent operation is to look for an atom that is bonded to a ring atom, but outside the ring
//  All it does is check a0 and a2

    int is_next_or_previous (atom_number_t) const;
};

extern int path_length_comparitor_longer  (Path * const *, Path * const *);
extern int path_length_comparitor_shorter (Path * const *, Path * const *);
extern int fused_system (const resizable_array_p<Path> & , int, int *, int);

/*
  This function is intended for the case of looking at two fused
  rings and identifying the fused atoms at the join
*/

extern int find_adjacent_atoms_in_common_between_two_rings (const Ring & ri, 
                      const Ring & rj,
                      atom_number_t & a1,
                      atom_number_t & a2);

#endif
