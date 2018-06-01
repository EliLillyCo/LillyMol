#ifndef BOND_H
#define BOND_H

#include <iostream>
#include <functional>

#include "iwmtypes.h"

#ifndef MOLECULE_H
class Molecule;
#endif

class IWString;

class Connection
{
  protected:
    atom_number_t _a2;
    bond_type_t   _btype;

  public:
    Connection ();
    Connection (atom_number_t, bond_type_t);
    ~Connection ();

    int           debug_print (std::ostream &);

    atom_number_t a2 () const { return _a2; }

    bond_type_t   btype () const { return _btype; }

    void        set_bond_type   (bond_type_t);

    void        set_aromatic ();
    void        set_non_aromatic ();

    void        set_permanent_aromatic (int);    // 0 or 1 values only
};

#define BOND_PROPERTY_UNKNOWN -3

#define BONDS_KNOW_RING_MEMBERSHIP
#ifdef BONDS_KNOW_RING_MEMBERSHIP
#define UNKNOWN_BOND_NRINGS -3
#endif

/*
  The possible things to go in _directional

  First, participation in a cis-trans bonding scheme
*/

#define IW_BOND_DIRECTIONAL_UP 0x1
#define IW_BOND_DIRECTIONAL_DOWN 0x2

/*
  It is handy if we know if a double bond is part of a cis-trans bonding arrangement
*/

#define IW_BOND_DIRECTIONAL_DOUBLE_BOND 0x4

/*
  These will be read from and written to MDL connection tables
  WEDGE_UP means an up wedge from _a1 to _a2
*/

#define IW_BOND_DIRECTION_WEDGE_UP 0x8
#define IW_BOND_DIRECTION_WEDGE_DOWN 0x10

#define IW_BOND_DIRECTION_WEDGE_EITHER 0x20

/*
  This next one corresponds to directionality 3 in the MDL file format - cis-trans double bond either
  Not used for smiles
*/

#define IW_BOND_CIS_TRANS_DIRECTIONAL_EITHER 0x40

class Bond: public Connection
{
  friend
    std::ostream & operator << (std::ostream &, const Bond &);

  private:
    atom_number_t    _a1;

//  This single variable holds info on wedge bonding and directional
//  bonding for cis-trans bonds.

    int _directional;

#ifdef BONDS_KNOW_RING_MEMBERSHIP
    int             _nrings;
#endif

//  Sept 2003, decided that bonds may know their bond number

    int _bond_number;

// private functions

    void _default_values ();

  public:
    Bond (atom_number_t, atom_number_t, bond_type_t);
    Bond (const Molecule *, atom_number_t, atom_number_t, bond_type_t);
    Bond (const Bond &);
    Bond ();
    ~Bond ();

    int ok () const;
#ifdef BONDS_SHOULD_NOT_KNOW_ABOUT_MOLECULES
    int ok (const Molecule *) const;
#endif

    void copy_directionality_specifications (const Bond *);

    int debug_print (std::ostream &) const;

    atom_number_t a1 ()    const { return _a1; }
    atom_number_t a2 ()    const { return _a2; }
    bond_type_t   btype () const { return _btype; }
    atom_number_t other (atom_number_t i) const 
    {
      if (i == _a1)
        return _a2;

      return _a1;
    }

//  If both bonds aromatic they are the same, else compare by bond type

    int same_bond_type (const Bond & b) const;

    int           involves (atom_number_t q) const { return q == _a1 || q == _a2;}

//  Involves these two atoms

    int           involves (atom_number_t, atom_number_t) const;

//  Do two bonds share a common end

    int           joins (const Bond * b) const { return _a1 == b->_a1 || _a1 == b->_a2 ||
                                                        _a2 == b->_a1 || _a2 == b->_a2;}

//  Do two bonds share a common end. If so, put the common atom in A

    int           joins (const Bond * b, atom_number_t & a) const;

//  If the bond involves a given atom, return 1 and the identity of the other atom

    int           involves_and_what_is_other (atom_number_t, atom_number_t &) const;

    void          adjust_for_loss_of_atom (atom_number_t);

    void          set_bond_type (bond_type_t);

    int           number_of_bonds () const;

    void          set_directional_up   ();
    void          set_directional_down ();

//  Set the directionality according to the atoms. The direction is up/down from s1 to s2

    void          set_directional_up   (atom_number_t s1, atom_number_t s2);
    void          set_directional_down (atom_number_t s1, atom_number_t s2);

    void          set_not_directional () { _directional = 0;}
    int           is_directional () const { return (_directional & (IW_BOND_DIRECTIONAL_UP|IW_BOND_DIRECTIONAL_DOWN)); }
    int           is_directional_up () const { return (_directional & IW_BOND_DIRECTIONAL_UP);}
    int           is_directional_down () const { return (_directional & IW_BOND_DIRECTIONAL_DOWN);}

//  Double bonds can be part of a cis-trans grouping

    int           part_of_cis_trans_grouping () const { return (_directional & IW_BOND_DIRECTIONAL_DOUBLE_BOND);}
    int           set_part_of_cis_trans_grouping (int s);

    void          set_wedge_up ();
    int           is_wedge_up () const { return (_directional & IW_BOND_DIRECTION_WEDGE_UP);}
    void          set_wedge_down ();
    int           is_wedge_down () const { return (_directional & IW_BOND_DIRECTION_WEDGE_DOWN);}
    void          set_wedge_either ();
    int           is_wedge_either () const { return (_directional & IW_BOND_DIRECTION_WEDGE_EITHER);}

    void          set_cis_trans_either_double_bond () { _directional = (_directional | IW_BOND_CIS_TRANS_DIRECTIONAL_EITHER);}
    int           is_cis_trans_either_double_bond () const { return _directional & IW_BOND_CIS_TRANS_DIRECTIONAL_EITHER;}

//  Definitive means either UP or DOWN but not EITHER

    int           is_wedge_definitive () const { return 0 != (_directional & (IW_BOND_DIRECTION_WEDGE_UP | IW_BOND_DIRECTION_WEDGE_DOWN)) ;}
    int           is_wedge_any () const { return 0 != (_directional & (IW_BOND_DIRECTION_WEDGE_UP | IW_BOND_DIRECTION_WEDGE_DOWN | IW_BOND_DIRECTION_WEDGE_EITHER)) ;}

    int           swap_atoms (atom_number_t, atom_number_t);

    int           set_a1 (atom_number_t);
    int           set_a2 (atom_number_t);
    int           set_a1a2 (atom_number_t, atom_number_t);

//  Used with VDOM - max efficiency

    void          set_a1a2_btype (atom_number_t newa1, atom_number_t newa2, bond_type_t newbt) {
                                  _a1 = newa1; _a2 = newa2; _btype = newbt;}

    int           is_single_bond () const { return IS_SINGLE_BOND (_btype);}
    int           is_double_bond () const { return IS_DOUBLE_BOND (_btype);}
    int           is_triple_bond () const { return IS_TRIPLE_BOND (_btype);}
    int           is_aromatic () const { return IS_AROMATIC_BOND (_btype);}
    int           is_permanent_aromatic () const { return IS_PERMANENT_AROMATIC_BOND (_btype);}

    void          append_bond_type (IWString &, atom_number_t afrom, int inc_arom) const;   // used in smiles construction
    void          append_bond_type_space_for_nothing(IWString & smiles, atom_number_t ato, int include_aromaticity_in_smiles) const;

#ifdef BONDS_KNOW_RING_MEMBERSHIP
    int set_nrings (int);
    int nrings (int &) const;
    int nrings () const;
    int in_another_ring ();
    int nrings_known () const;
    void invalidate_nrings ();
#endif

    int            bond_number_assigned () const { return _bond_number >= 0;}
    int            bond_number () const { return _bond_number;}
    void           set_bond_number (int s) { _bond_number = s;}
    void           invalidate_bond_number () { _bond_number = -1;}

    int            either_atom_set_in_array (const int *) const;

    int            either_atom_true (const int *, std::function<int(int, int)>) const;            // an array of integers,

    int            new_atom_numbers (const int * xref);   // makeing a subset, adjust the bonds to reflect the new numbering
};

#endif
