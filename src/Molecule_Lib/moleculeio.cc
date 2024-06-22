// Holds the global settings defined in moleculeio.h

#include <iostream>
#include <limits>

#include "moleculeio.h"

namespace moleculeio {

using std::cerr;

int _ignore_all_chiral_information_on_input = 0;

void
set_ignore_all_chiral_information_on_input(int s) {
  _ignore_all_chiral_information_on_input = s;

  return;
}

int
ignore_all_chiral_information_on_input() {
  return _ignore_all_chiral_information_on_input;
}

// Another way for fileconv to function is to use chiral information
// only if it is correct

int _ignore_incorrect_chiral_input = 0;

void
set_ignore_incorrect_chiral_input(int i) {
  _ignore_incorrect_chiral_input = i;

  return;
}

int
ignore_incorrect_chiral_input() {
  return _ignore_incorrect_chiral_input;
}

int _flush_files_after_writing_each_molecule = 0;

int
flush_files_after_writing_each_molecule() {
  return _flush_files_after_writing_each_molecule;
}

void
set_flush_files_after_writing_each_molecule(int i) {
  _flush_files_after_writing_each_molecule = i;
}

// When inputting MOLFILE's or TDT's we can optionally save all
// the non-connection-table records in the molecule's _extra_info
// array.

static int _read_extra_text_info = 0;

void
set_read_extra_text_info(int r)
{
  _read_extra_text_info = r;
}

int
read_extra_text_info() {
  return _read_extra_text_info;
}

// Does the extra text info get written or not

static int _write_extra_text_info = 0;

void set_write_extra_text_info(int w) {
  _write_extra_text_info = w;
}

int
write_extra_text_info() {
  return _write_extra_text_info;
}

// Especially with 3rd party molecules we can get files without a newline

char record_delimiter = '\n';

char 
input_file_delimiter() {
  return record_delimiter;
}

void
set_record_delimiter(char s) {
  record_delimiter = s;
}

int dos_mode = 1;    // Mar 2005. Change to default

int 
input_is_dos_mode() {
  return dos_mode;
}

void
set_dos_mode(int s) {
  dos_mode = s;
}

IWString file_scope_newline_string('\n');

void
generate_newline_string(IWString & newline_string)
{
  if (write_DOS_records()) {
    newline_string.resize_keep_storage(0);
    newline_string << static_cast<char>(13) << '\n';   // cannot put newline in src
  }
  else 
    newline_string = '\n';

  return;
}

const IWString & 
newline_string() {
  return file_scope_newline_string;
}

static int _write_DOS_records = 0;

void 
set_write_DOS_records(int s) {
  _write_DOS_records = s;

  generate_newline_string(file_scope_newline_string);
}

int
write_DOS_records() {
  return _write_DOS_records;
}

/*
  When reading MDL files we can look at the coordinates to perceive
  cis-trans bonds
*/

int _discern_cis_trans_bonds = 0;

int
discern_cis_trans_bonds()
{
  return _discern_cis_trans_bonds;
}

void
set_discern_cis_trans_bonds(int s)
{
  _discern_cis_trans_bonds = s;
}

int _discern_chirality_from_3d_coordinates = 0;

void
set_discern_chirality_from_3d_coordinates(int s)
{
  _discern_chirality_from_3d_coordinates = s;
}

int
discern_chirality_from_3d_coordinates()
{
  return _discern_chirality_from_3d_coordinates;
}

int _ignore_bad_cis_trans_input = 0;
int 
ignore_bad_cis_trans_input() {
  return _ignore_bad_cis_trans_input;
}

void
set_ignore_bad_cis_trans_input(int s) {
  _ignore_bad_cis_trans_input = s;
}

int _discard_directional_bonds_on_input = 0;

void
set_discard_directional_bonds_on_input(int s) {
  _discard_directional_bonds_on_input = s;
}

int discard_directional_bonds_on_input() {
  return _discard_directional_bonds_on_input;
}

int _number_connection_table_errors_to_skip = 0;
void
set_number_connection_table_errors_to_skip(int s) {
  _number_connection_table_errors_to_skip = s;
}

int number_connection_table_errors_to_skip() {
  return _number_connection_table_errors_to_skip;
}

int _skip_first_molecules = 0;
void
set_skip_first_molecules(int s) {
  _skip_first_molecules = s;
}
int
skip_first_molecules() {
  return _skip_first_molecules;
}

int _do_only_n_molecules = 0;

void
set_do_only_n_molecules(int s) {
  _do_only_n_molecules = s;
}

int
do_only_n_molecules() {
  return _do_only_n_molecules;
}

static off_t _seek_to_from_command_line = 0;

void set_seek_to(off_t o)
{
  _seek_to_from_command_line = o;
}

off_t
seek_to_from_command_line()
{
  return _seek_to_from_command_line;
}

static off_t _max_offset_from_command_line = std::numeric_limits<off_t>::max();

off_t
max_offset_from_command_line()
{
  return _max_offset_from_command_line;
}

void
set_max_offset_from_command_line(off_t s)
{
  _max_offset_from_command_line = s;
}

/*
  Optionally remove all bonds from non-organic elements.
  Possible values:
    0  do nothing
    1  remove all bonds to metals
    2  remove all bonds to a non-organic if all atoms attached are heteroatoms.

  TODO: make this an enumeration
*/

int _unconnect_covalently_bonded_non_organics_on_read = 0;

int
unconnect_covalently_bonded_non_organics_on_read()
{
  return _unconnect_covalently_bonded_non_organics_on_read;
}

void
set_unconnect_covalently_bonded_non_organics_on_read(int s)
{
  _unconnect_covalently_bonded_non_organics_on_read = s;
}

//  Converts a string to one of our known file types. 
//  Does its work silently.
//  Returns FILE_TYPE_INVALID if a match is not found.

FileType
_string_to_file_type(const const_IWSubstring & file_type)
{
  assert (file_type.nchars());

  if ("mdl" == file_type)
    return FILE_TYPE_MDL;
  if ("pdb" == file_type)
    return FILE_TYPE_PDB;
  if ("mmod" == file_type)
    return FILE_TYPE_MMOD;
  if ("smi" == file_type)
    return FILE_TYPE_SMI;
  if ("usmi" == file_type)
    return FILE_TYPE_USMI;
  if ("msi" == file_type)
    return FILE_TYPE_MSI;
  if ("tdt" == file_type)
    return FILE_TYPE_TDT;
  if ("gfp" == file_type)
    return FILE_TYPE_TDT;
  if ("utdt" == file_type)
    return FILE_TYPE_UTDT;
  if ("tdtnausmi" == file_type)
    return FILE_TYPE_TDT_NAUSMI;
  if ("rdf" == file_type)
    return FILE_TYPE_RDF;
  if ("qry" == file_type)
    return FILE_TYPE_QRY;
  if ("rsmi" == file_type)
    return FILE_TYPE_RSMI;
  if ("sdf" == file_type)
    return FILE_TYPE_SDF;
  if ("mol" == file_type)
    return FILE_TYPE_SDF;
  if ("mol2" == file_type)
    return FILE_TYPE_MOL2;
  if ("chm" == file_type)
    return FILE_TYPE_CHM;
  if ("moe" == file_type)
    return FILE_TYPE_MOE;
  if ("mrk" == file_type)
    return FILE_TYPE_MRK;
  if ("wchm" == file_type)
    return FILE_TYPE_WCHM;
  if ("nausmi" == file_type)
    return FILE_TYPE_NAUSMI;
  if ("cif" == file_type)
    return FILE_TYPE_CIF;
  if ("smt" == file_type)
    return FILE_TYPE_SMT;
  if ("mrv" == file_type)
    return FILE_TYPE_MRV;
  if ("inchi" == file_type)
    return FILE_TYPE_INCHI;
  if ("csv" == file_type)
    return FILE_TYPE_CSV;
  if ("textproto" == file_type ||
      "txtproto" == file_type ||
      "txt" == file_type)
    return FILE_TYPE_TXTPROTO;
  if ("xyz" == file_type)
    return FILE_TYPE_XYZ;
  
  return FILE_TYPE_INVALID;
}

// Complains if unrecognised
FileType
string_to_file_type(const const_IWSubstring & file_type)
{
  assert (file_type.nchars());

  FileType rc = _string_to_file_type(file_type);
  if (rc != FILE_TYPE_INVALID) {
    return rc;
  }

  cerr << "string_to_file_type: unrecognised type '" << file_type << "'\n";
  
  return rc;
}

void
ResetFileScopeVariables() {
  _number_connection_table_errors_to_skip = 0;

  _ignore_all_chiral_information_on_input = 0;
  _ignore_incorrect_chiral_input = 0;
  _flush_files_after_writing_each_molecule = 0;
  _read_extra_text_info = 0;
  _write_extra_text_info = 0;
  record_delimiter = '\n';
  dos_mode = 1;
  _write_DOS_records = 0;
  _discern_cis_trans_bonds = 0;
  _discern_chirality_from_3d_coordinates = 0;
  file_scope_newline_string = '\n';

  _discern_cis_trans_bonds = 0;
  _ignore_bad_cis_trans_input = 0;

  _discard_directional_bonds_on_input = 0;
  _skip_first_molecules = 0;
  _do_only_n_molecules = std::numeric_limits<int>::max();

  _seek_to_from_command_line = 0;
  _max_offset_from_command_line = std::numeric_limits<off_t>::max();

  _unconnect_covalently_bonded_non_organics_on_read = 0;
}

}  // namespace moleculeio
