#ifndef IW_AROMATIC_H
#define IW_AROMATIC_H

#include <iostream>

/*
  Aromaticity can be either by Daylight rules, or Pearlman's
*/

#define Simple_4n_plus_2 1
#define Daylight 2
#define Pearlman 3

/*
  Dec 98. The Wang Fu Lai clogp paper. They use rules which look mostly like
  Pearlman rules, but they aromatise furan and the sulphur analogue
*/

#define WangFuLai 4

#define Vijay_Gombar 5

#define EVERYTHING_HAS_A_PI_ELECTRON 6

#define Pipeline_Pilot 7

#define PUBCHEM_AROMATICITY 8

#define ANY_EVEN_NUMBER_OF_PI_ELECTRONS 9

extern int set_global_aromaticity_type (int);
extern int global_aromaticity_type ();

extern int display_standard_aromaticity_options (std::ostream &);

class Command_Line;

extern int process_standard_aromaticity_options (Command_Line &, int = 0, char = 'A');

extern int input_aromatic_structures ();
extern void set_input_aromatic_structures (int);

extern int allow_input_without_valid_kekule_form  ();
extern void set_allow_input_without_valid_kekule_form (int);

extern int  allow_delocalised_carbonyl_bonds ();
extern void set_allow_delocalised_carbonyl_bonds (int);

extern int  discard_non_aromatic_kekule_input ();
extern void set_discard_non_aromatic_kekule_input (int);

extern int convert_chain_aromatic_bonds ();
extern void set_convert_chain_aromatic_bonds (int);

extern void set_aromatic_chain_bonds_are_ok (int s);

/*
  When outputting forms other than SMILES, are aromatic bonds written as
  aromatic, or as kekule forms
*/

extern void set_write_aromatic_bonds (int);
extern int write_aromatic_bonds ();

extern void set_warn_aromatic_chain_atoms (int);

extern void set_kekule_try_positive_nitrogen (int s);

extern void set_all_bonds_in_aromatic_ring_must_be_aromatic (int s);

extern void set_display_no_kekule_form_message(int s);
extern int  display_no_kekule_form_message();

extern void set_allow_pipeline_pilot_aromaticity_on_input(int s);

/*
  The largest ring that can be aromatic
*/

extern void set_max_aromatic_ring_size (int s);

extern void set_perform_kekule_perception (int s);

extern void reset_aromatic_file_scope_variables ();
extern void  reset_mdl_file_scope_variables ();

extern void set_strongly_fused_rings_can_be_aromatic(int);
extern void set_aromatic_rings_must_contain_unsaturation(const int s);

#endif

/* arch-tag: cd7a479d-ca70-4047-b1d1-1d7beb0bf194 */
