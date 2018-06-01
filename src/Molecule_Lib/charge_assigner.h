#ifndef CHARGE_ASSIGNER_H
#define CHARGE_ASSIGNER_H

/*
*/

#include "iwaray.h"

#include "iwmtypes.h"
#include "temp_detach_atoms.h"
#include "qry_wstats.h"

class Substructure_Hit_Statistics;
class Single_Substructure_Query;
class Substructure_Atom;
class Molecule_to_Match;
class Command_Line;

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
    int _positive_isotope;
    const Element * _negative_element;
    int _negative_isotope;

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

//  Feb 2006. I want to know which queries were responsible for placing the
//  charges. To do this, we can assign each formal charge to a multiplier
//  times the query number. So, if query 0 assigned a +1 charge, that would
//  still be +1, but if query 2 assigned +1, that would be 201. If query 2
//  assigned a negative charge, it would be -201

    int _assigned_charge_multiplier;

//  private functions

    void _increment_global_counters (const Set_of_Atoms & positive_charges_assigned,
                                             const Set_of_Atoms & negative_charges_assigned);

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

    int _identify_charged_atoms_too_close (Molecule & m,
                                           const Set_of_Atoms & s,
                                           int * times_too_close) const;
    void _remove_hits_too_close (Molecule & m,
                                         Set_of_Atoms & s,
                                         formal_charge_t * charges_assigned) const;
    void _remove_hits_too_close_isolation_score (Molecule & m,
                                         Set_of_Atoms & s,
                                         formal_charge_t * charges_assigned) const;
    int _remove_lower_preference_hits(Molecule & m,
                                Substructure_Results & sresults) const;
    int _do_make_implicit_hydrogens_explicit (Molecule & m, const int * charges_assigned) const;
    int _numeric_value_present (const Substructure_Atom * a) const;
    int _numeric_value_present (const Single_Substructure_Query * q) const;
    int _all_queries_have_numeric_value () const;

    int _process (Molecule &, resizable_array_p<Molecule> & charged_forms);
    int _process (Molecule_to_Match &, Set_of_Atoms &, Set_of_Atoms &, formal_charge_t *);
    int _process (Molecule &, formal_charge_t *);

  public:
    Charge_Assigner ();
    ~Charge_Assigner ();

    int construct_from_command_line (Command_Line &, int = 0, char = 'G');
    int build (const const_IWSubstring &);

    void set_verbose (int v) { _verbose = v;}

    int report (std::ostream &) const;

    void set_assigned_charge_multiplier (int s) { _assigned_charge_multiplier = s;}

    void set_min_distance_between_charges (int s) { _min_distance_between_charges = s;}

//  array is actual charges assigned

    int process (Molecule &, formal_charge_t * = NULL);
    int process (Molecule &, resizable_array_p<Molecule> & charged_forms);

    int number_queries () const { return _number_elements;}

    int active () const { return _number_elements;}
    void deactivate () { resize(0);}

    void set_apply_charges_to_molecule (int s) { _apply_charges_to_molecule = s;}
};

extern void display_standard_charge_assigner_options (std::ostream &, char);

#endif
