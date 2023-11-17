#ifndef MOLECULE_LIB_TOGGLE_KEKULE_FORM_H_
#define MOLECULE_LIB_TOGGLE_KEKULE_FORM_H_

#include "Molecule_Lib/toggle_kekule_form.pb.h"

#include "bond.h"
#include "molecule.h"

class Command_Line;
class msi_attribute;

/*
  We want to fix a number of molecules into a single Kekule form
*/

class Toggle_Kekule_Form_Temporary_Arrays {
 private:
  //  At the end of the computation we recompute the aromaticity to see if
  //  all our atoms are still aromatic. Store the starting aromaticity.

  aromaticity_type_t* _arom;

  //  For efficiency, we get a copy of the atoms in the molecule

  Atom** _atom;

  // We keep track of which atoms are mentioned in a directive to
  // set a given bond.

  int* _specified_by_bonds;

  //  We keep track of whether or not each atom already has a double bond

  int* _has_double_bond;

  //  We also keep track of which bonds can change. This is matoms*matoms
  int* _can_change_bond;

  // A temporary array.
  int* _process_these;

  //  We keep track of which bonds are OK already
  int* _correct;

  //  When we get an atom that cannot be changed, all rings that contain
  //  that atom cannot change
  int* _ring_can_vary;

  int* _atom_can_change;

  // If multiple bonds in a ring are covered by the atoms specified
  // in _bonds, we mark that ring as being definitely a candidate for change.
  int* _ring_specified_as_changeable;

 public:
  Toggle_Kekule_Form_Temporary_Arrays(Molecule&);
  ~Toggle_Kekule_Form_Temporary_Arrays();

  int*
  can_change_bond() {
    return _can_change_bond;
  }

  aromaticity_type_t*
  arom() {
    return _arom;
  }

  Atom**
  atom() {
    return _atom;
  }

  int*
  has_double_bond() {
    return _has_double_bond;
  }

  int*
  process_these() {
    return _process_these;
  }

  int*
  correct() {
    return _correct;
  }

  void
  set_ring_can_toggle(int r, int s) {
    _ring_can_vary[r] = s;
  }

  int
  ring_can_toggle(int r) const {
    return _ring_can_vary[r];
  }

  void
  set_atom_can_change(int a, int s) {
    _atom_can_change[a] = s;
  }

  int
  atom_can_change(int a) const {
    return _atom_can_change[a];
  }

  int*
  ring_specified_as_changeable() {
    return _ring_specified_as_changeable;
  }

  int
  ring_specified_as_changeable(int rnum) const {
    return _ring_specified_as_changeable[rnum];
  }

  int*
  specified_by_bonds() {
    return _specified_by_bonds;
  }

  int
  specified_by_bonds(int s) {
    return _specified_by_bonds[s];
  }
};

class Toggle_Kekule_Form {
 private:
  // Each bond specifies a pair of (bonded) matched atoms and the bond between then.
  // Even though a Bond holds atom_number_t type, they are interpreted as indices into
  // an embedding.
  resizable_array_p<Bond> _bond;

  //  By default, we do NOT allow a pyrrole Nitrogen to toggle

  int _allow_pyrrole_to_change;

  int _display_error_messages;

  //  If we get in molecules that have come from reactions, there may be issues with the
  //  Hydrogen count

  int _unset_unnecessary_implicit_hydrogens_known_values;

  // If we have made a change to the ring system already, it may no longer be aromatic,
  // so make it possible to NOT check whether all the bonds requested start in an aromatic
  // ring.
  int _check_all_bonds_aromatic;

  //  private functions

  void _no_changes_to_atom(Molecule& m, atom_number_t zatom,
                           Toggle_Kekule_Form_Temporary_Arrays&) const;

  void _set_all_bonds_to_single(Molecule& m, int id,
                                Set_of_Atoms& double_bonds_to_be_restored,
                                Toggle_Kekule_Form_Temporary_Arrays&) const;

  int _restore_previous_bonding(Molecule& m,
                                const Set_of_Atoms& bonds_to_be_restored) const;

  int _bond_is_correct(const Molecule& m, const Set_of_Atoms& embedding,
                       const Bond* b) const;
  int _all_bonds_correct(const Molecule& m, const Set_of_Atoms& embedding,
                         Toggle_Kekule_Form_Temporary_Arrays&) const;
  int _all_bonds_aromatic(Molecule& m, const Set_of_Atoms& embedding) const;
  int _all_atoms_aromatic(Molecule& m, int, int,
                          Toggle_Kekule_Form_Temporary_Arrays&) const;
  int MultipleBondsChanged(const Ring& r, const Set_of_Atoms& embedding) const;

  int grow_ring_system(Molecule& m, const Ring* r, int id,
                       Toggle_Kekule_Form_Temporary_Arrays& tkfta,
                       int* ring_already_done);

  int _get_ring_system_atoms(resizable_array<int>& atoms_to_process, int rid,
                             atom_number_t zatom,
                             Toggle_Kekule_Form_Temporary_Arrays&) const;
  int _set_our_bonds(Molecule& m, const Set_of_Atoms& embedding, int id,
                     Toggle_Kekule_Form_Temporary_Arrays&) const;

  int _do_not_process_rings_containing(Molecule& m, atom_number_t zatom,
                                       Toggle_Kekule_Form_Temporary_Arrays& tkfta) const;

  void IdentifyChangeableRings(Molecule& m, const Set_of_Atoms& embedding,
                     Toggle_Kekule_Form_Temporary_Arrays&) const;
  int TurnOffBondsAlreadyCorrect(Molecule& m,
                const Set_of_Atoms& embedding,
                Toggle_Kekule_Form_Temporary_Arrays& tkfta) const;
  void _do_chemistry_aromatic_ring(Molecule& m, const Ring& r,
                                   Toggle_Kekule_Form_Temporary_Arrays& tkfta) const;
  int _process_ring_system(Molecule& m, resizable_array<int>& atoms_to_process, int rid,
                           int zitem, const atom_number_t prev,
                           Toggle_Kekule_Form_Temporary_Arrays&) const;
  int _process_ring_system(Molecule& m, const Set_of_Atoms& embedding, int rid,
                           Toggle_Kekule_Form_Temporary_Arrays&);
  int _process(Molecule& m, const int* process_these, int* already_done);
  int _process(Molecule& m, const Set_of_Atoms& embedding, int* process_these);
  int _process(Molecule& m, const Set_of_Atoms& embedding,
               Toggle_Kekule_Form_Temporary_Arrays&);
  int _process_single_ring(Molecule& m, const Set_of_Atoms& embedding, int id,
                           Toggle_Kekule_Form_Temporary_Arrays&);
  int _process_single_ring2(Molecule& m, const Set_of_Atoms& embedding, int id,
                            Toggle_Kekule_Form_Temporary_Arrays&);
  int IdentifyRingsSpecifiedByDirectives(Molecule& m, const Set_of_Atoms& embedding,
                                         Toggle_Kekule_Form_Temporary_Arrays& tkfta);
  int IdentifyAtomsInvolvedInBonds(const Set_of_Atoms& embedding,
                                   Toggle_Kekule_Form_Temporary_Arrays& tkfta) const;

 public:
  Toggle_Kekule_Form();
  ~Toggle_Kekule_Form();

  int ok() const;
  int debug_print(std::ostream&) const;

  int
  active() const {
    return _bond.number_elements();
  }

  void
  set_display_error_messages(int s) {
    _display_error_messages = s;
  }

  void
  set_unset_unnecessary_implicit_hydrogens_known_values(const int s) {
    _unset_unnecessary_implicit_hydrogens_known_values = s;
  }

  void
  set_check_all_bonds_aromatic(int s) {
    _check_all_bonds_aromatic = s;
  }

  const Bond* contains_bond(atom_number_t a1, atom_number_t a2) const;
  int will_change_ring(const Ring* r, const Set_of_Atoms&) const;

  int construct_from_command_line(Command_Line&, char, int = 0);
  int add_bond_from_msi_attribute(const msi_attribute&);
  int add_bond(int, int, bond_type_t);
  int add_bond(Bond* b);

  int ConstructFromProto(const ToggleKekuleForm::ToggleKekuleForm& proto);

  void
  set_allow_pyrrole_to_change(int s) {
    _allow_pyrrole_to_change = s;
  }

  int write_msi(std::ostream&, const IWString&, const char*) const;

  int ok_embedding(const Set_of_Atoms& embedding) const;

  int process(Molecule&, const Set_of_Atoms&, int&);

  int process(Molecule&, atom_number_t a1, atom_number_t a2, bond_type_t bt,
              int& changed);

  void adjust_for_subset(const int* include_these_atoms);

  void
  reset() {
    _bond.resize(0);
    return;
  }
};

#endif  // MOLECULE_LIB_TOGGLE_KEKULE_FORM_H_
