#ifndef MOLECULE_LIB_MOLECULEIO_H_
#define MOLECULE_LIB_MOLECULEIO_H_

// Functions controlling aspects of reading writing and building molecules.
// This is WIP, there are more functions still in molecule.h that need to
// be migrated here, and some of the functions listed here have not been
// implemented.

#include <sys/types.h>

#include "Foundational/iwstring/iwstring.h"

#include "iwmtypes.h"

namespace moleculeio {
int set_store_pdb_atom_information(int s);
int set_store_pdb_atom_information(int);
int set_use_stored_atom_information_when_writing_pdb_files(int);

const IWString & newline_string();

char input_file_delimiter();
void set_record_delimiter(char s);

int input_is_dos_mode();
void set_dos_mode(int);
void set_write_DOS_records(int);
int write_DOS_records();

// Ultimately these should move to the istream_and_type class.
int number_connection_table_errors_to_skip();
void set_number_connection_table_errors_to_skip(int s);

int skip_first_molecules();
void set_skip_first_molecules(int s);

int do_only_n_molecules();
void set_do_only_n_molecules(int s);

//int valid_file_type(FileType);
////int suffix_for_file_type(FileType);

// Too disruptive to move this here.
//FileType discern_file_type_from_name(IWString const&);

void set_seek_to(off_t);
off_t seek_to_from_command_line();

off_t max_offset_from_command_line();
void set_max_offset_from_command_line(off_t);

void set_discern_cis_trans_bonds(int);
int discern_cis_trans_bonds();
void set_discern_chirality_from_3d_coordinates(int);
int  discern_chirality_from_3d_coordinates();


void set_ignore_all_chiral_information_on_input(int s);
void set_ignore_incorrect_chiral_input(int s);
int ignore_incorrect_chiral_input();

void set_ignore_bad_cis_trans_input(int s);
int ignore_bad_cis_trans_input();

void set_discard_directional_bonds_on_input(int s);
int discard_directional_bonds_on_input();

int convert_from_mdl_charge(int);
int int3d(const_IWSubstring const&, int&, int&, int*);
int ignore_all_chiral_information_on_input();
int read_extra_text_info();
void set_read_extra_text_info(int s);
//int create_file_with_appropriate_name(const_IWSubstring const&, IWString&, FileType, int);
void set_flush_files_after_writing_each_molecule(int);
int flush_files_after_writing_each_molecule();
void set_write_extra_text_info(int);
int write_extra_text_info();

FileType string_to_file_type(const_IWSubstring const&);

int set_mol2_write_formal_charge_as_partial_charge(int);
int set_pdb_number_within_sequence(int);
int set_pdb_number_by_element_count(int);
int set_write_pdb_files_in_fragment_order(int);

int unconnect_covalently_bonded_non_organics_on_read();
void set_unconnect_covalently_bonded_non_organics_on_read(int);

//void set_tdt_identifier_dataitem(const const_IWSubstring &);
//void set_tdt_append_dataitem(const const_IWSubstring &);

void ResetFileScopeVariables();
}  // namespace moleculeio
#endif
