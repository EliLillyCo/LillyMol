#ifndef IW_TARGET_H
#define IW_TARGET_H

#include <iostream>

#include "iwbits.h"

#include "iwmtypes.h"
#include "molecule.h"

class Target_Atom;

class Bond_and_Target_Atom
{
  protected:

    int _nbonds;      // the  number of bonds 1-4

    bond_type_t _bond_types;     // bits for each bond type

    int _nrings;
    int _aromatic;

    const Bond * _bond;
    Target_Atom * _other;

  public:
    Bond_and_Target_Atom ();

    int ok () const;

    void initialise (const Bond *, Target_Atom &);
    void initialise (const Bond *, Target_Atom &, int nr);

    void establish_aromatic_bonds ();
    void establish_aromatic_bonds (int matoms, atom_number_t a1, const int * fixed_kekule_form);

    const Bond  * bond  () const { return _bond;}
    Target_Atom * other () const { return _other;}


    bond_type_t btype () const { return _bond_types;}
    int nrings ();
    boolean aromatic () const { return _aromatic;}
};

class Molecule_to_Match;

class Target_Atom
{
  private:
    Molecule * _m;
    atom_number_t _my_atom_number;   // to which atom number in _m do I correspond.
    Atom *        _my_atom;          // the actual atom in _m

//  Target_Atom * _all_atoms;
    Molecule_to_Match * _target;    // object that owns this atom

    int _ncon;
    Bond_and_Target_Atom * _other;

    const Element * _element;
    int _nbonds;
    formal_charge_t _formal_charge;
    int _nrings;
    int _ring_bond_count;
    int _ncon2;
    int _hcount;
    aromaticity_type_t _aromaticity;
    int _multiple_bond_to_heteroatom;

    int _vinyl;
    int _aryl;

    int _lone_pair_count;

//  We store a pointer to a chiral centre if present. We need to be
//  able to distinguish between not yet fetched and NULL (which means
//  no chiral centre on this atom)

    int                _chirality_fetched;
    Chiral_Centre *    _chiral_centre;

    int                _all_rings_kekule;

    List_of_Ring_Sizes _sssr_ring_sizes;
    List_of_Ring_Sizes _aromatic_ring_sizes;
    List_of_Ring_Sizes _aliphatic_ring_sizes;

    List_of_Ring_Sizes _ring_sizes_including_non_sssr;
    int _attached_heteroatom_count;
    int _fused_system_size;
    int _isotope;
    int _userAtomType;
    int _heteroatoms_in_ring;

//  private functions

    void _default_values ();
    int  _compute_ncon2  ();
    int  _allocate_other ();

    void _get_ring_sizes_and_aromaticity ();
    int  _count_heteroatoms_in_ring (const Ring * r) const;

  public:
    Target_Atom ();
    ~Target_Atom ();

    int ok () const;
    int debug_print (std::ostream &) const;

    void initialise (Molecule *, atom_number_t, Atom *, Molecule_to_Match *);
    void establish_aromatic_bonds ();
//  void establish_aromatic_bonds (int, const int * in_fixed_kekule_form);

    const Bond_and_Target_Atom & operator [] (int s) const { return _other[s];}

    atomic_number_t atomic_number () const { return _element->atomic_number();}
    int atomic_symbol_hash_value () const { return _element->atomic_symbol_hash_value();}
    int element_unique_id () const { return _element->unique_id();}
    const Element * element () const { return _element;}

    int ncon () const { return _ncon;}
    void set_ncon (int s) { _ncon = s;}   // could get out of sync with _number_elements
    int nbonds ();
    void set_nbonds(int s) { _nbonds = s;}
    int nrings ();
    void set_nrings (int s) { _nrings = s;}
    int ring_bond_count ();
    void set_ring_bond_count(int s) { _ring_bond_count = s;}
    int is_ring_atom ();
    int is_non_ring_atom ();
    int ncon2 ();
    void set_ncon2(int s) { _ncon2 = s;}
    int  ncon2_value_set() const;   // ask whether or not the value has been computed/set
    formal_charge_t formal_charge ();
    void set_formal_charge (formal_charge_t s) { _formal_charge = s;}
    atom_number_t atom_number () const { return _my_atom_number;}
    int hcount ();
    void set_hcount (int s) { _hcount = s;}
    int daylight_x ();
    aromaticity_type_t aromaticity ();
    int is_aromatic ();
    void set_is_aromatic (int s);
    int vinyl ();
    int aryl ();
    void set_aryl (int s) { _aryl = s;}

    int all_rings_kekule();
    void set_all_rings_kekule(const int s);

//  If we want to exclude chirality in substructure searches, or
//  during unique determinations, we can suppress chirality by
//  marking it as already fetched

    void discard_chirality();

    Chiral_Centre * chiral_centre ();

    int attached_heteroatom_count ();
//  int isolated_ring ();
    int fused_system_size ();
    int multiple_bond_to_heteroatom ();
    int lone_pair_count ();
    int isotope ();
    int userAtomType ();
    int heteroatoms_in_ring ();

    int fragment_membership ();

//  Because ring_sizes is a multi-valued thing, that is handled differently

    const List_of_Ring_Sizes * sssr_ring_sizes ();
    const List_of_Ring_Sizes * aromatic_ring_sizes ();
    const List_of_Ring_Sizes * aliphatic_ring_sizes ();
    const List_of_Ring_Sizes * ring_sizes_including_non_sssr ();

    int in_same_rings (Target_Atom *);
    int is_bonded_to (Target_Atom *) const;

    Bond_and_Target_Atom & other (int);
    Bond_and_Target_Atom * bond_to_atom (atom_number_t);

    int atoms_in_target_molecule () const { return _m->natoms ();}

    Molecule * m () const { return _m;}

    const Atom * atom () const { return _my_atom;}

//  We can sometimes speed things up by ensuring that a given atom doesn't match anything

    void invalidate ();

    void set_isotope (int i) { _isotope = i;}
    void set_userAtomType (int i) { _userAtomType = i;}
    void set_element (const Element * e) { _element = e;}

    void set_atom_number(int s) { _my_atom_number = s;}

    int  is_spinach ();

    Target_Atom & atom (int i) const;

};
#define TARGET_IS_RING 0
#define TARGET_BETWEEN_RING -1
#define TARGET_SPTMP -3
#define TARGET_UNSET -4

class Molecule_to_Match
{
  private:
    Molecule * _m;
    int _natoms;
    Atom ** _atom;

    Target_Atom * _target_atom;

    int _nrings;
    int _aromatic_ring_count;
    int _non_aromatic_ring_count;
    int _fused_ring_count;
    int _strongly_fused_ring_count;
    int _isolated_ring_count;
    int _ring_object_count;
    int _number_isotopic_atoms;

//  We need to keep track of whether or not _establish_aromatic_bonds has been called

    int _establish_aromatic_bonds_called;

    resizable_array<const Ring *> _rings;

//  sometimes a substructure search can gain some speed from
//  searching only a range of the atoms. For each atomic
//  number, we record the first and last occurrence

    int _first[HIGHEST_ATOMIC_NUMBER + 1];
    int _last[HIGHEST_ATOMIC_NUMBER + 1];
    int _count[HIGHEST_ATOMIC_NUMBER + 1];

//  Apr 03. Need queries dealing with spinach and atoms between rings. If the atom is in a ring
//  value will be TARGET_IS_RING. If it is between rings, value will be TARGET_BETWEEN_RING. Otherwise
//  it is the number of atoms in this particular piece of spinach

    int * _spinach_or_between_rings;

//  Oct 2003. Want to be able to start a query at a particular atom

    atom_number_t _start_matching_at;

    IW_Bits_Base * _fingerprint;


//  private functions


    int _determine_ring_counts ();

    int _grow_spinach (atom_number_t zatom, int flag);
    int _identify_between_rings (atom_number_t aprev, atom_number_t zatom);

    int _initialise_spinach_or_between_rings ();

    void _initialise_molecule (Molecule * m);

  public:
    Molecule_to_Match ();
    Molecule_to_Match (Molecule *);
    ~Molecule_to_Match ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int initialise_molecule (Molecule * m);

    int natoms () const { return _natoms;}

    int number_fragments () const { return _m->number_fragments ();}

    Target_Atom & operator [] (int i) const { return _target_atom[i];}

    void discard_chirality_information();

    int first (atomic_number_t) const;
    int last  (atomic_number_t) const;
    int first_atomic_symbol_hash_value (int) const;
    int last_atomic_symbol_hash_value  (int) const;

    int atoms_with_atomic_number (atomic_number_t) const;

    int is_superset (const IW_Bits_Base &) const;

    atom_number_t start_matching_at () const { return _start_matching_at;}
    void set_start_matching_at (atom_number_t s) { _start_matching_at = s;}

    Molecule * molecule () const { return _m;}

//  this is dangerous. It was added for use with Valhalla, where something
//  that isn't a MOlecule at all is passed here. Beware...

    void set_molecule (Molecule * s) { _m = s;}

    int compute_aromaticity ();

    void establish_aromatic_bonds ();

    int nrings ();
    int aromatic_ring_count ();
    int non_aromatic_ring_count ();

    int fused_ring_count ();
    int strongly_fused_ring_count ();
    int isolated_ring_count ();

    int ring_object_count ();

    const Ring * ringi (int i);

    int number_isotopic_atoms ();
    int number_isotopic_atoms (int);

    int heteroatom_count () const;

    int net_formal_charge() const { return _m->net_formal_charge();}

    void invalidate (const Set_of_Atoms &);   // after getting a match you can invalidate all the matched atoms

    int is_spinach (atom_number_t);
    int is_between_rings (atom_number_t);

};

extern void set_aromatic_bonds_lose_kekule_identity (int);
extern int  aromatic_bonds_lose_kekule_identity ();

/*
  Computing the number of instances of each element may speed up
  some searches
*/

extern void set_initialise_element_counts (int);
extern void set_global_setting_nrings_includes_non_sssr_rings(int s);
extern void set_global_setting_nbonds_includes_implicit_hydrogens(int s);

#endif
