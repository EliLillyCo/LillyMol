#ifndef IWFIX_STRUCTURES_H
#define IWFIX_STRUCTURES_H

#include <iostream>

#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/substructure.h"

/*
  Correct various structural defects that are most likely just 
  transcription errors.
*/

class Structure_Fixing
{
  private:
    int _active;

    int _directive_fix_almost_nitro;     // We convert -N(=O)-O to a Nitro group

//  We put a formal charge on -N(-*)= when the Nitrogen is in a ring that looks like a Pyridine

    int _directive_fix_substituted_pyridine;

//  If we see a neutral Nitrogen atom with 4 bonds, put a charge on it

    int _directive_fix_four_valent_nitrogens_charged;

//  The atom might be OK, but it has a strange hydrogen count

    int _directive_remove_hydrogens_known_flag_to_fix_valence_errors;

//  maybe if a charge is removed, the atom will now be OK.

    int _directive_try_removing_charges_on_atoms_with_bad_valences;

//  SD3G0 is a valence error, but with a +ve charge it is OK

    int _directive_add_positive_charge_to_three_valent_sulphur;

//  ID2G0 is a valence error, but with a -ve charge it is OK

    int _directive_add_positive_charge_to_two_valent_halogen;

//  Almost any formal charge on a carbon should be driven out.

    int _directive_remove_charges_from_carbon_atoms;

    int _directive_remove_charges_from_halogen_atoms;

    // de-novo molecules may have bad implicit Hydrogens =[OH]
    int _directive_fix_obviously_wrong_implicit_hydrogens;

    Substructure_Query _pyridine_query;

    int _molecules_examined;
    int _molecules_changed;

    int _partial_nitros_changed;
    int _substituted_pyridines_changed;
    int _four_valent_nitrogens_charged;
    int _three_valent_sulphurs_changed;
    int _two_valent_halogens_changed;
    int _charges_removed;
    int _charged_carbons_neutralised;
    int _charged_halogens_neutralised;

    IWString _append_to_changed_structures;

//  private functions

    int _do_try_removing_charges_on_atoms_with_bad_valences (Molecule & m, const Atom ** atom);
    int _do_fix_four_valent_nitrogens_charged (Molecule & m, const Atom ** atom);
    int _do_fix_almost_nitro (Molecule & m, const Atom ** atom,
                     int & has_3_connected_4_bonds_uncharged_nitrogen);
    int _do_fix_substituted_pyridine (Molecule & m, const Atom ** atom);
    int _do_add_positive_charge_to_three_valent_sulphur (Molecule & m, const Atom ** atom);
    int _do_add_positive_charge_to_two_valent_halogen (Molecule & m, const Atom ** atom);
    int _do_remove_charges_from_halogen_atoms(Molecule & m, const Atom ** atom);
    int _do_remove_charges_from_carbon_atoms(Molecule & m, const Atom ** atom);
    int _process (Molecule & m, const Atom ** atom);

  public:
    Structure_Fixing ();

    int initialise (Command_Line & cl, char c, int vb);

    int usage (char flag, std::ostream &) const;

    int active () const { return _active;}

    int report (std::ostream &) const;

    void set_append_to_changed_structures (const IWString & s) { _append_to_changed_structures = s;}

    int process (Molecule &);
};

#endif
