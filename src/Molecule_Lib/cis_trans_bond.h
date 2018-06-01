#ifndef CIS_TRANS_BOND_H
#define CIS_TRANS_BOND_H

#include <iostream>
#include "iwmtypes.h"

class IWString;

class Cis_Trans_Bond;

class Cis_Trans_Bond_Leaf
{
  private:
    atom_number_t    _a;
    Cis_Trans_Bond * _coupled;

//  private functions

    void _default_values ();
    int _is_coupled (atom_number_t my_root, 
                     atom_number_t other_root,
                     atom_number_t other);

  public:
    Cis_Trans_Bond_Leaf ();
    Cis_Trans_Bond_Leaf (const Cis_Trans_Bond_Leaf &, int = 0);

    atom_number_t atom () const { return _a;}
    void set_atom (atom_number_t aa) { _a = aa;}

    void swap (Cis_Trans_Bond_Leaf &);

    Cis_Trans_Bond * coupled () { return _coupled;}
    void  set_coupled (Cis_Trans_Bond * c) { _coupled = c;}

    int is_coupled (atom_number_t root, Cis_Trans_Bond * ctb, int & orientation);

    int renumber (const int *);

    Cis_Trans_Bond_Leaf & operator = (atom_number_t);

    void adjust_for_loss_of_atom (atom_number_t);
};

inline int operator == (const Cis_Trans_Bond_Leaf & l, const atom_number_t a)
{
  return l.atom () == a;
}

inline int operator == (atom_number_t a, const Cis_Trans_Bond_Leaf & l)
{
  return a == l.atom ();
}

inline int operator != (const atom_number_t a, const Cis_Trans_Bond_Leaf & l)
{
  return l.atom () != a;
}

inline std::ostream & operator << (std::ostream & os, const Cis_Trans_Bond_Leaf & l)
{
  return os << l.atom ();
}

class Bond;

class Cis_Trans_Bond
{
  private:
    atom_number_t _left_root;

    Cis_Trans_Bond_Leaf _left_up;
    Cis_Trans_Bond_Leaf _left_down;

    atom_number_t _right_root;

    Cis_Trans_Bond_Leaf _right_up;
    Cis_Trans_Bond_Leaf _right_down;

    int _couple_id;

//  private functions

    void _default_values ();

  public:
    Cis_Trans_Bond ();
    Cis_Trans_Bond (const Cis_Trans_Bond &, int = 0);

    int ok () const;
    int debug_print (std::ostream &) const;

    void set_left_root (atom_number_t a) { _left_root = a;}
    atom_number_t left_root () const { return _left_root;}
    atom_number_t left_up   () const { return _left_up.atom ();}
    atom_number_t left_down () const { return _left_down.atom ();}

    void set_left_up   (atom_number_t a) { _left_up = a;}
    void set_left_down (atom_number_t a) { _left_down = a;}

    int set_bond_directionality (Bond &) const;

    int  extra_lhs (atom_number_t, int);

    atom_number_t right_root () const { return _right_root;}
    atom_number_t right_up   () const { return _right_up.atom ();}
    atom_number_t right_down () const { return _right_down.atom ();}

    void set_right_root (atom_number_t a) { _right_root = a;}
    void set_right_up   (atom_number_t a) { _right_up = a;}
    void set_right_down (atom_number_t a) { _right_down = a;}

    int  extra_rhs (atom_number_t, int);

    int involves (atom_number_t) const;
    int involves (atom_number_t, atom_number_t) const;

    int  couple_id () const { return _couple_id;}
    void set_couple_id (int s) { _couple_id = s;}

    int invert_left_right ();
    int invert_up_down ();

    int reorient (int orientation1, int orientation2);

    int make_consistent (atom_number_t, atom_number_t, int);

    int renumber (const int *);

    int adjust_for_loss_of_atom (atom_number_t);

    int change_atom_number (atom_number_t, atom_number_t);

    int process_directional_bond_for_smiles (IWString & smiles,
                                              atom_number_t a, atom_number_t anchor) const;

    int discern_coupling (Cis_Trans_Bond *, int &);

    int on_same_side (atom_number_t, atom_number_t) const;
//  int append_bond_symbol (IWString &, atom_number_t, atom_number_t) const;

    int atom_is_now_implicit_hydrogen (atom_number_t);
};

#endif
