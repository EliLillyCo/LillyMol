#ifndef IW_CHIRAL_CENTRE_H
#define IW_CHIRAL_CENTRE_H

#include "iwmtypes.h"

#include "iwaray.h"

class IWString;
class Molecule;
class Bond;

#define CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN (INVALID_ATOM_NUMBER - 4)
#define CHIRAL_CONNECTION_IS_LONE_PAIR         (INVALID_ATOM_NUMBER - 5)

/*
  There is an unusual twist to the Chiral_Centre class.
  It is used to not only describe chiral centres who's chirality
  is known, but also to identify chiral centres who's chirality
  is undetermined.

  The flag _chirality_known indicates this. When 0 == _chirality_known,
  the object will not add any chirality specifiers to a smiles

*/

class Chiral_Centre
{
  private:
    atom_number_t _a;

    atom_number_t _top_front, _top_back, _left_down, _right_down;

    int _chirality_known;

// Private functions

    int _discern_clockwise (IWString & smiles,
                            const int * zorder,
                            atom_number_t south_west,
                            atom_number_t north,
                            atom_number_t south_east,
                            const resizable_array<atom_number_t> & ring_openings,
                            const resizable_array<atom_number_t> & ring_closures) const;

    void _default_values ();

  public:
    Chiral_Centre (atom_number_t);
    Chiral_Centre (const Chiral_Centre &);
    ~Chiral_Centre ();

    int ok () const;
    int debug_print (std::ostream &) const;

    atom_number_t a () const { return _a;}
    void set_centre (atom_number_t s) { _a = s;}     // no checking

    int make_copy (const Chiral_Centre &, const int *);

    int invert ();

//  Note that mdl.cc uses the value of _chirality_known as a
//  temporary storage for the type of the chiral centre, so
//  make sure we don't ever change this feature. Could get
//  around it by making mdl.cc allocate an array of chirality
//  types....

    int  chirality_known () const { return _chirality_known;}
    void set_chirality_known (int i) { _chirality_known = i;}

    int complete () const;

//  Here connections include things like implicit hydrogens

    int number_connections_specified() const;

//  Here atoms mean actual atom numbers

    int number_atoms_specified () const;

    int involves (atom_number_t at) const;

    int involves (atom_number_t, atom_number_t) const;

    int implicit_hydrogen_count () const;
    int lone_pair_count () const;

    int centre_atom_has_a_lone_pair ();

    void adjust_for_loss_of_atom (atom_number_t);
    int  adjust_atom_numbers (const int * xref);

    int  convert_to_implicit_hydrogen (atom_number_t);

    int top_front () const { return _top_front;}
    int top_back  () const { return _top_back;}
    int left_down () const { return _left_down;}
    int right_down () const { return _right_down;}

    int set_top_front (atom_number_t ntf);
    int set_top_back  (atom_number_t ntb);
    int set_left_down (atom_number_t nld);
    int set_right_down (atom_number_t nrd);

    int got_ring_opening_bond (int ring_number, int chiral_count);
    int got_ring_closure_bond (int ring_number, atom_number_t a);

    int mdl_stereo_centre_value () const;

    int mdl_stereo_centre_value (atom_number_t sw,
                                 atom_number_t north,
                                 atom_number_t se) const;


    int append_smiles_chirality_symbol (IWString & smiles, const int * zorder,
                         atom_number_t previous_atom,
                         const resizable_array<const Bond *> & ring_opening_bonds,
                         const resizable_array<atom_number_t> & ring_closures) const;

//  add_molecule uses this function to copy chiral centres from one molecule to another

    int make_copy (const Chiral_Centre *);

//  When we make implicit hydrogens explicit, we need a means of telling
//  a chiral centre what the new atom number is

    int implicit_hydrogen_is_now_atom_number (atom_number_t);

    int lone_pair_is_now_atom_number (atom_number_t);

//  Similarly, when we break a bond, we may have an atom which now becomes an
//  implicit hydrogen or a lone pair

    int atom_is_now_implicit_hydrogen (atom_number_t);
    int atom_is_now_lone_pair         (atom_number_t);

    int change_atom_number (atom_number_t, atom_number_t);
    int move_atom_to_end_of_atom_list (atom_number_t, int);

    int atom_numbers_are_swapped (atom_number_t, atom_number_t);

    int make_top_front (atom_number_t);

//  The Molecule::create_subset often needs to know whether or not the
//  atoms in this chiral centre are all in the subset

    int all_atoms_in_subset (const int *, int id) const;

//  Unique smiles with chirality need a means of providing a -1,0,1 ranking
//  for the atoms involved in a chiral centre

    int orientation (const unsigned int *) const;

    int influence (const unsigned int *, atom_number_t, atom_number_t) const;
    int influence (const unsigned int *, atom_number_t) const;

//  A convenient way of cycling through all the atoms in a Chiral_Centre object. Note that
//  it does not return any implicit hydrogens or lone pairs

    atom_number_t next_atom (int &) const;

    int set_vector (int *, int) const;

    void new_atom_numbers (const int *);     // someone is changing the atoms in the molecule
};

#endif

