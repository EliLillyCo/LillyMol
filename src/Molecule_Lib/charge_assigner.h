#ifndef MOLECULE_LIB_CHARGE_ASSIGNER_H_
#define MOLECULE_LIB_CHARGE_ASSIGNER_H_

#include <iostream>
#include <optional>
#include <vector>

#include "Foundational/iwaray/iwaray.h"

#include "iwmtypes.h"
#include "qry_wstats.h"
#include "temp_detach_atoms.h"

class Substructure_Hit_Statistics;
class Single_Substructure_Query;
class Substructure_Atom;
class Molecule_to_Match;
class Command_Line;

struct ChargeAndQuery {
  atom_number_t atom;
  formal_charge_t formal_charge;
  int query_number;
  bool operator==(const ChargeAndQuery& rhs) const;
};

extern std::ostream& operator<<(std::ostream& output, const ChargeAndQuery& afq);

class Charge_Assigner : public resizable_array_p<Substructure_Hit_Statistics>
{
  private:
    int _verbose;

    int _overwrite_existing_formal_charges;

    int _molecules_examined;
    int _molecules_changed;
    int _negative_charges_assigned;
    int _positive_charges_assigned;

//  We need to make sure that we don't assign charges too close to each other

    int _min_distance_between_charges;
    
//  Another possibility is to replace the atoms with a different element

    const Element * _positive_element;
    isotope_t _positive_isotope;
    const Element * _negative_element;
    isotope_t _negative_isotope;

//  We can isotopically label the positive and negative atoms with
//  the atomic number of the atom from which they were derived

    int _apply_isotopic_labels;

//  When changing the element type, we can run into problems with implicit
//  hydrogens, chiral centres and such. Every time we assign a different element,
//  we can optionally remove any chiral centre info from that atom

    int _remove_chiral_centres_from_changed_atoms;

//  Another possibility is to preserve the implicit hydrogen count

    int _preserve_implicit_hydrogen_count;

//  Sometimes we will want to add an explicit Hydrogen atom to an atom that 
//  has just been assigned a formal positive charge.

    int _attach_explicit_hydrogen_to_positive_atom;

    int _molecules_receiving_negative_charges;
    int _molecules_receiving_positive_charges;

//  If we have a simple query - just one set of query atoms, we can
//  compute one time which atom in the embedding bears the charge, and
//  what that charge is.

    int * _which_atom;
    formal_charge_t * _charge_to_assign;

    Temp_Detach_Atoms _temp_detach_hydrogens;

//  Mar 2004. I want to be able to apply charges, but not change the molecule, just
//  get the results in the array passed to process ()

    int _apply_charges_to_molecule;

//  private functions

    int BuildFromDefaultEnvs();
    int BuildFromEnv(const IWString& env);
    int BuildFromEnvValue(const IWString& dir);

    void InitialiseWhichAtom();

    void _increment_global_counters (const Set_of_Atoms & positive_charges_assigned,
                                             const Set_of_Atoms & negative_charges_assigned);
    void _increment_global_counters(std::vector<ChargeAndQuery>& results);

    int _enumerate_possibilities0 (Molecule & m,
                                          const Set_of_Atoms & s,
                                          resizable_array_p<Set_of_Atoms> & possibilities) const;
    int _enumerate_possibilities1 (Molecule & m,
                                          const Set_of_Atoms & s,
                                          int istart,
                                          resizable_array_p<Set_of_Atoms> & possibilities) const;
    void _remove_positive_charge_hits_on_chiral_atoms (Molecule & m,
                        Set_of_Atoms & positive_charges_assigned,
                        formal_charge_t * charges_assigned) const;

    int _add_charge_to_atom (Molecule_to_Match & target,
                           Set_of_Atoms & positive_charges_assigned,
                           Set_of_Atoms & negative_charges_assigned,
                           formal_charge_t * charges_assigned,
                           atom_number_t zatom,
                           formal_charge_t fc) const;
    int _add_charge_to_result(Molecule_to_Match& target,
                        int query_number,
                        atom_number_t zatom,
                        formal_charge_t fc,
                        formal_charge_t* charges_assigned,
                        std::vector<ChargeAndQuery>& result);

    int _identify_charged_atoms_too_close (Molecule & m,
                                           const Set_of_Atoms & s,
                                           int * times_too_close) const;

    int AtomsTooClose(Molecule& m, atom_number_t a1, atom_number_t a2) const;

    void _remove_hits_too_close (Molecule & m,
                                         Set_of_Atoms & s,
                                         formal_charge_t * charges_assigned) const;
    void _remove_hits_too_close_isolation_score (Molecule & m,
                                         Set_of_Atoms & s,
                                         formal_charge_t * charges_assigned) const;
    void RemoveSitesTooClose(Molecule & m,
                formal_charge_t* charges_assigned,
                std::vector<ChargeAndQuery>& result);
    int _remove_lower_preference_hits(Molecule & m,
                                Substructure_Results & sresults) const;
    int _do_make_implicit_hydrogens_explicit (Molecule & m, const int * charges_assigned) const;
    int _numeric_value_present (const Substructure_Atom * a) const;
    int _numeric_value_present (const Single_Substructure_Query * q) const;
    int _all_queries_have_numeric_value () const;

    int MaybeUnscaled(formal_charge_t fc) const;

    int _process (Molecule &, resizable_array_p<Molecule> & charged_forms);
    int _process (Molecule_to_Match &, Set_of_Atoms &, Set_of_Atoms &, formal_charge_t *);
    int _process (Molecule &, formal_charge_t *);
    int _process(Molecule& m, formal_charge_t* charges_assigned, std::vector<ChargeAndQuery>& result);

  public:
    Charge_Assigner ();
    ~Charge_Assigner ();

    int construct_from_command_line (Command_Line& cl, int verbose = 0, char flag = 'G');
    int build (const const_IWSubstring &);

    void set_verbose (int v) { _verbose = v;}

    int report (std::ostream &) const;

    void set_min_distance_between_charges (int s) { _min_distance_between_charges = s;}

//  array is actual charges assigned

    int process (Molecule &, formal_charge_t * = nullptr);
    int process (Molecule &, resizable_array_p<Molecule> & charged_forms);

    int Process(Molecule& m, std::vector<ChargeAndQuery>& afq);

    int number_queries () const { return _number_elements;}

    int active () const { return _number_elements;}
    void deactivate () { resize(0);}

    void set_apply_charges_to_molecule (int s) { _apply_charges_to_molecule = s;}
};

extern void display_standard_charge_assigner_options (std::ostream &, char);

#endif  // MOLECULE_LIB_CHARGE_ASSIGNER_H_
