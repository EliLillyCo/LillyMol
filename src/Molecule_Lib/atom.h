#ifndef MOLECULE_LIB_ATOM_H_
#define MOLECULE_LIB_ATOM_H_

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Foundational/iwaray/iwaray.h"

#include "coordinates.h"
#include "element.h"
#include "iwmtypes.h"

#define ATOM_PROPERTY_UNKNOWN -476

#include "bond.h"
#include "set_of_atoms.h"

struct BondAndAtom {
  const Bond* bond;
  atom_number_t atom;
};

class __attribute__((visibility("default"))) Atom : public resizable_array<Bond*>,
                                                    public Coordinates,
                                                    public std::enable_shared_from_this<Atom>{
 private:
  const Element* _element;

  //  Normally the implicit hydrogen count is a computed value. Sometimes
  //  however a known hcount value is read from a file.

  short _implicit_hydrogens;
  short _implicit_hydrogens_known;  // known value from file

  short _formal_charge;
  short _radical;
  short _nbonds;

  //  Some atoms may be classified as permanently aromatic

  short _permanent_aromatic;

  isotope_t _isotope;
  int _userAtomType;

  void* _user_specified_void_ptr;

  int _atom_map;

  //  private functions

  void _default_values(const Element*);
  void _constructor_copy_atom_attributes(const Atom& other_atom);
  int _compute_implicit_hydrogens(int& result);

  template <typename T>
  int _common_saturation(const T& comparitor) const;

 public:
  Atom(const char*);
  Atom(const Element*);
  Atom(atomic_number_t);
  Atom(const Atom*);  // make a copy of an atom - connections are not copied
  Atom(const Atom&);
  ~Atom();

  int debug_print(std::ostream&) const;
  std::string debug_string() const;

  int ok() const {
    return nullptr != _element;
  }

  int audit() const;

  isotope_t isotope() const {
    return _isotope;
  }

  void set_isotope(isotope_t i) {
    _isotope = i;
  }

  int is_isotope() const;

  int is_organic() const {
    return _element->organic();
  }

  int userAtomType() const {
    return _userAtomType;
  }

  void set_userAtomType(int i) {
    _userAtomType = i;
  }

  int add(Bond*);

  // Number of atoms connected.
  int ncon() const {
    return _number_elements;
  }

  // how many of the atoms set in the array are connected
  int ncon(atom_number_t, const int* include_atom) const;

  int nbonds();
  int nbonds() const;
  int recompute_nbonds();

  int permanent_aromatic() const {
    return _permanent_aromatic;
  }

  void set_permanent_aromatic(int s) {
    _permanent_aromatic = s;
  }

  int molecule_being_resized(int);

  const Bond* bond_to_atom(atom_number_t) const;
  const Bond* bond_to_atom(atom_number_t myAtomId, atom_number_t otherAtomId) const;

  int is_halogen() const;

  const IWString& atomic_symbol() const;

  const Element* element() const {
    return _element;
  }

  const Element& elementq() const {
    return *_element;
  }

  atomic_number_t atomic_number() const;

  int atomic_symbol_hash_value() const {
    return _element->atomic_symbol_hash_value();
  }

  void set_element(const Element*);  // dangerous to have public

  int implicit_hydrogens_computed() const {
    return _implicit_hydrogens >= 0;
  }

  int implicit_hydrogens();  // not const as it stores the result
  int recompute_implicit_hydrogens(int&);
  int compute_implicit_hydrogens(int&);  // note this does not store the result, it is non
                                         // const because it may compute nbonds
  int set_implicit_hydrogens(int, int = 0);

  void unset_all_implicit_hydrogen_information();

  //  during Kekule determinations I may change the bonds to an atom and need a quick
  //  means of having that atom reset anything that may be dependent on its bonding

  void invalidate_computed_values_after_bond_change();

  int implicit_hydrogens_known() const {
    return _implicit_hydrogens_known;
  }

  void set_implicit_hydrogens_known(int i);

  formal_charge_t formal_charge() const {
    return formal_charge_t(_formal_charge);
  }

  int set_formal_charge(formal_charge_t);

  int is_radical() const {
    return _radical;
  }

  void set_radical(const int s) {
    _radical = s;
  }

  void set_modified();

  // In order to facilitate reproducible things, it can be useful
  // to place the bonds attached into a canonical order.
  int CanonicaliseBonds();

  //  void set_cb(int, atom_number_t, bond_type_t);
  //  int  add_con(atom_number_t, bond_type_t);
  int remove_bonds_to_atom(atom_number_t, int = 0);

  atom_number_t other(atom_number_t, int) const;  // atom number of i'th connection
  bond_type_t btype_to_connection(int) const;     // bond order  to i'th connection
  bond_type_t btype_to_atom(atom_number_t myAtomId, atom_number_t otherAtomId)
      const;  // bond order of bond to otherAtomId

  int other_and_type(atom_number_t, int, atom_number_t&, bond_type_t&) const;

  //  void set_bond_type_to_connection(int, bond_type_t);    // set bond order of i'th
  //  connection
  int set_bond_type_to_atom(atom_number_t, bond_type_t);

  //  void set_atom_number(int, int);        // set atom number of i'th connection
  //  void set_directional_bond(int, int);    // set the I'th bond to be directional

  int is_centre_of_cis_trans_bond() const;

  int number_directional_bonds_attached() const;

  int connections(atom_number_t, atom_number_t*, bond_type_t* = nullptr) const;
  int connections(atom_number_t, Set_of_Atoms&) const;
  Set_of_Atoms connections(atom_number_t zatom) const;

  // Enable structured bindings for the Bonds and atoms connected.
  //   const Atom& atom = ...
  //   for (const auto [bond, atom] : atom.BondsAndConnections(zatom)) {
  //      if (bond->is_single_bond())...
  //      if (atom == previous_atom)....
  //   }
  // Convenient to use, but slower than accessing Bonds sequentially.
  std::vector<BondAndAtom> BondsAndConnections(atom_number_t zatom) const;

  // Not implemented yet, easy to do, not sure needed...
  // A more efficient, but less convenient version
  //   BondAndAtom bond_and_atom;
  //   for (int i = 0; atom.BondAndConnection(zatom, i, bond_and_atom); ) {
  //   }
  int BondAndConnection(atom_number_t zatom, int& iter, BondAndAtom& result) const;

  int connections_and_types(atom_number_t, Set_of_Atoms&,
                            resizable_array<bond_type_t>&) const;
  int connections_and_types(atom_number_t, atom_number_t*, bond_type_t*) const;

  int bond_types(resizable_array<bond_type_t>&) const;
  int bond_types(bond_type_t*) const;

  int is_bonded_to(atom_number_t) const;  // are we bonded to atom number J
  int is_bonded_to(atom_number_t, bond_type_t&) const;

  atomic_mass_t atomic_weight() const;

  int exact_mass(exact_mass_t&) const;

  int lone_pair_count(int&);
  int pi_electrons(int&);

  int valence_ok();

  //  When making a smiles, it is often convenient to have multiple bonds first

  void multiple_bonds_first();

  //  various things for smiles

  int next_atom_for_smiles(atom_number_t my_atom_number, int* already_done,
                           int* canonical_order, atom_number_t& next_atom) const;
  int next_atom_for_unique_smiles(atom_number_t my_atom_number, int* already_done,
                                  int* canonical_order, atom_number_t& next_atom) const;
  int next_atom_for_random_smiles(atom_number_t my_atom_number, int* already_done,
                                  int* canonical_order, atom_number_t& next_atom) const;

  //  Several molecule formats need coordinates written in a common format

  int write_coordinates(std::ostream&, int = 0) const;

  int fully_saturated() const;  // nbonds() == ncon)
  int unsaturated() const;      // nbonds() < ncon()
  int unsaturation() const;

  const void* user_specified_void_ptr() const {
    return _user_specified_void_ptr;
  }

  void* user_specified_void_ptr() {
    return _user_specified_void_ptr;
  }

  void set_user_specified_void_ptr(void* v) {
    _user_specified_void_ptr = v;
  }

  int remove_connections_to_any_of_these_atoms(const int*);

  int atom_map() const {
    return _atom_map;
  }

  void set_atom_map(const int s) {
    _atom_map = s;
  }

  // A common operation is to ask if this atom `zatom` is
  //   a. singly connected
  //   b. with a double bond.
  // In that case return the other atom,
  // otherwise return kInvalidAtomNumber.
  atom_number_t SingleConnectionDoubleBondTo(atom_number_t zatom) const;
};

extern void reset_atom_file_scope_variables();

// Returns the dihedral angle implied by a1-a2-a3-a4
extern angle_t angle_between_atoms(const Atom& a1, const Atom& a2, const Atom& a3,
                                   const Atom& a4);

#define OK_ATOM(a) (nullptr != (a) && (a)->ok())

// Deprecated, not useful.
// extern int how_many_atoms();

extern Coordinates form_unit_vector(const Atom&, const Atom&);

extern void set_display_abnormal_valence_messages(int);
extern int display_abnormal_valence_messages();

extern void set_copy_implicit_hydrogen_count_in_atom_copy_constructor(int s);

extern int set_reasonable_formal_charge_range(formal_charge_t, formal_charge_t);
extern int reasonable_formal_charge_value(formal_charge_t c);
extern void set_reset_implicit_hydrogens_known_on_bond_removal(const int s);
extern void set_four_connected_neutral_nitrogen_has_h(const int s);

extern void set_alternate_valences_give_hcount(const int s);
#endif  // MOLECULE_LIB_ATOM_H_
