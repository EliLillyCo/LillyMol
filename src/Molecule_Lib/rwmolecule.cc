/*
  Various functions associated with reading and writing structure files
*/

#include <iostream>
#include <limits>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/string_data_source.h"

#include "element.h"
#include "mdl.h"
#include "misc2.h"
#include "molecule.h"
#include "moleculeio.h"
#include "smiles.h"

using std::cerr;
using std::endl;

/*
  Fileconv can be run so as to ignore all chiral information on input
*/

#ifdef NOW_IN_MOLECULEIO
static int _ignore_all_chiral_information_on_input = 0;

void
set_ignore_all_chiral_information_on_input(int i)
{
  _ignore_all_chiral_information_on_input = i;

  return;
}

int
ignore_all_chiral_information_on_input()
{
  return _ignore_all_chiral_information_on_input;
}

/*
  Another way for fileconv to function is to use chiral information
  only if it is correct
*/

static int _ignore_incorrect_chiral_input = 0;

void
set_ignore_incorrect_chiral_input(int i)
{
  _ignore_incorrect_chiral_input = i;

  return;
}

int
ignore_incorrect_chiral_input()
{
  return _ignore_incorrect_chiral_input;
}

static int _flush_files_after_writing_each_molecule = 0;

int
flush_files_after_writing_each_molecule()
{
  return _flush_files_after_writing_each_molecule;
}

void
set_flush_files_after_writing_each_molecule(int i)
{
  _flush_files_after_writing_each_molecule = i;
}

/*
  When inputting MOLFILE's or TDT's we can optionally save all
  the non-connection-table records in the molecule's _extra_info
  array.
*/

static int _read_extra_text_info = 0;

void
set_read_extra_text_info(int r)
{
  _read_extra_text_info = r;
}

/*
  When reading MDL files we can look at the coordinates to perceive
  cis-trans bonds
*/

static int _discern_cis_trans_bonds = 0;

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

static int _discern_chirality_from_3d_coordinates = 0;

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

static int _ignore_bad_cis_trans_input = 0;

void
set_ignore_bad_cis_trans_input(int s)
{
  _ignore_bad_cis_trans_input = s;
}

int
ignore_bad_cis_trans_input()
{
  return _ignore_bad_cis_trans_input;
}

int
read_extra_text_info()
{
  return _read_extra_text_info;
}

/*
  Does the extra text info get written or not
*/

static int _write_extra_text_info = 0;

void set_write_extra_text_info(int w)
{
  _write_extra_text_info = w;
}

int
write_extra_text_info()
{
  return _write_extra_text_info;
}

/*
  Especially with 3rd party molecules we can get files without a newline
*/

static char record_delimiter = '\n';

char 
input_file_delimiter()
{
  return record_delimiter;
}

static int dos_mode = 1;    // Mar 2005. Change to default

int 
input_is_dos_mode()
{
  return dos_mode;
}

static int _skip_first_molecules = 0;

void
set_skip_first_molecules(int s)
{
  _skip_first_molecules = s;
}

int 
skip_first_molecules()
{
  return _skip_first_molecules;
}

static int _do_only_n_molecules = 0;

void
set_do_only_n_molecules(int d)
{
  _do_only_n_molecules = d;
}

int
do_only_n_molecules()
{
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
#endif

#ifdef NOW_IN_MOLECULEIO
static int _number_connection_table_errors_to_skip = 0;

int
number_connection_table_errors_to_skip()
{
  return _number_connection_table_errors_to_skip;
}

void
set_number_connection_table_errors_to_skip(int s)
{
  _number_connection_table_errors_to_skip = s;
}

/*
  Optionally remove all bonds from non-organic elements.
  Possible values:
    0  do nothing
    1  remove all bonds to metals
    2  remove all bonds to a non-organic if all atoms attached are heteroatoms.

  TODO: make this an enumeration
*/

static int _unconnect_covalently_bonded_non_organics_on_read = 0;

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
#endif

static int _put_formal_charges_on_neutral_ND3v4 = 0;

void
set_put_formal_charges_on_neutral_ND3v4(int s)
{
  _put_formal_charges_on_neutral_ND3v4 = s;
}

int
put_formal_charges_on_neutral_ND3v4()
{
  return _put_formal_charges_on_neutral_ND3v4;
}

// If coordinates are read in in some other units (Bohr?) we can 
// apply a scaling upon reading.
static float coordinate_scaling = 0.0f;

void set_coordinate_scaling(float scale) {
  coordinate_scaling = scale;
}

static int _ignore_self_bonds = 0;
void
set_ignore_self_bonds(int s) {
  _ignore_self_bonds = s;
}

int
ignore_self_bonds() {
  return _ignore_self_bonds;
}

namespace smiles {
// When reading a smiles, if the 2nd token looks like |...|, process
// that as a Chemaxon smiles extension
int _discern_chemaxon_smiles_extensions = 0;

void
SetDiscernChemaxonSmilesExtensions(int s) {
  _discern_chemaxon_smiles_extensions = s;
}

int
DiscernChemaxonSmilesExtensions() {
  return _discern_chemaxon_smiles_extensions;
}

int process_quoted_smiles = 0;
void
SetProcessQuotedSmiles(int s) {
  process_quoted_smiles = s;
}
int ProcessQuotedSmiles() {
  return process_quoted_smiles;
}

}  // namespace smiles

#ifdef NOW_IN_MOLECULEIO
static IWString file_scope_newline_string('\n');

void
generate_newline_string(IWString & newline_string)
{
  if (write_DOS_records())
  {
    newline_string.resize_keep_storage(0);
    newline_string << static_cast<char>(13) << '\n';   // cannot put newline in src
  }
  else 
    newline_string = '\n';

  return;
}

const IWString & 
newline_string()
{
  return file_scope_newline_string;
}

static int _write_DOS_records = 0;

void 
set_write_DOS_records(int s)
{
  _write_DOS_records = s;

  generate_newline_string(file_scope_newline_string);
}

int
write_DOS_records()
{
  return _write_DOS_records;
}
#endif

/*
  Return the appropriate suffix for a given file type.
*/

const char *
suffix_for_file_type(FileType file_type)
{
  
  switch (file_type)
  {
    case FILE_TYPE_MDL:
      return "mdl";
      break;

    case FILE_TYPE_PDB:
      return "pdb";
      break;

    case FILE_TYPE_MMOD:
      return "mmod";
      break;

    case FILE_TYPE_SMI:
      return "smi";
      break;

    case FILE_TYPE_USMI:
      return "smi";
      break;

    case FILE_TYPE_MSI:
      return "msi";
      break;

    case FILE_TYPE_TDT:
      return "tdt";
      break;

    case FILE_TYPE_UTDT:
      return "tdt";
      break;

    case FILE_TYPE_RDF:
      return "rdf";
      break;

    case FILE_TYPE_QRY:
      return "qry";
      break;

    case FILE_TYPE_RSMI:
      return "smi";
      break;

    case FILE_TYPE_SDF:
      return "sdf";
      break;

    case FILE_TYPE_MOL2:
      return "mol2";
      break;

    case FILE_TYPE_PSF:
      return "psf";
      break;

    case FILE_TYPE_CRD:
      return "crd";
      break;

    case FILE_TYPE_MRK:
      return "mrk";
      break;

    case FILE_TYPE_NAUSMI:
      return "smi";
      break;

    case FILE_TYPE_WCHM:
      return "chm";
      break;

    case FILE_TYPE_CIF:
      return "cif";
      break;

    case FILE_TYPE_SMT:
      return "smt";
      break;

    case FILE_TYPE_MOE:
      return "moe";
      break;

    case FILE_TYPE_MRV:
      return "mrv";
      break;

    case FILE_TYPE_INCHI:
      return "inchi";
      break;

    case FILE_TYPE_TDT_NAUSMI:
      return "tdt";
      break;

    case FILE_TYPE_CSV:
        return "csv";
        break;

    case FILE_TYPE_TXTPROTO:
        return "textproto";
        break;

    case FILE_TYPE_XYZ:
        return "xyz";
        break;

    default:
      return nullptr;
      break;
  }
  
  iwabort();
  return nullptr;    // should never come here
}

#ifdef NOLONGER_NEEDED_ASDASD
const char *
suffix_for_file_type(int file_type)
{
  switch (file_type)
  {
    case MDL:
      return "mdl";
      break;

    case PDB:
      return "pdb";
      break;

    case MMOD:
      return "mmod";
      break;

    case SMI:
      return "smi";
      break;

    case USMI:
      return "smi";
      break;

    case MSI:
      return "msi";
      break;

    case TDT:
      return "tdt";
      break;

    case UTDT:
      return "tdt";
      break;

    case RDF:
      return "rdf";
      break;

    case QRY:
      return "qry";
      break;

    case RSMI:
      return "smi";
      break;

    case SDF:
      return "sdf";
      break;

    case MOL2:
      return "mol2";
      break;

    case IWMTYPE_PSF:
      return "psf";
      break;

    case IWMTYPE_CRD:
      return "crd";
      break;

    case IWMTYPE_MOE:
      return "moe";
      break;

    case IWMTYPE_MRK:
      return "mrk";
      break;

    case IWMTYPE_NAUSMI:
      return "smi";
      break;

    case IWMTYPE_WCHM:
      return "chm";
      break;

    case IWMTYPE_CIF:
      return "cif";
      break;

    case IWMTYPE_SMT:
      return "smt";
      break;

    case IWMTYPE_MRV:
      return "mrv";
      break;

    case IWMTYPE_INCHI:
      return "inchi";
      break;

    case IWMTYPE_TDT_NAUSMI:
      return "tdt";
      break;

    default:
      return nullptr;
      break;
  }
  
  iwabort();
  return nullptr;    // should never come here
}
#endif

/*
  We must be careful stripping off old prefixes in case the directory name
  contains a period.

  Look for the last period and check to see if there are any directory
  separators out there
*/

static void
remove_suffix(IWString & fname)
{
  if (fname.ends_with(".gz"))
    fname.chop(3);

  int period = fname.rindex('.');
  if (period < 0)
    return;

  if (fname.ends_with('.'))
  {
    fname.chop();
    return;
  }

// only remove the suffix if it is a known suffix

  const_IWSubstring suffix = fname.substr(period + 1);

  if (FILE_TYPE_INVALID == moleculeio::string_to_file_type(suffix)) {   // not a known suffix
    return;
  }

  for (int i = period + 1; i < fname.nchars(); i++) {
    if ('/' == fname[i]) {
      return;
    }
  }

// No directory separators found

  fname.iwtruncate(period);

  return;
}

/*
  Given an old file name, what would be the corresponding name for a
  file with a different type.

  By default, we strip off any existing suffix, but any suffix can
  be retained if keep_existing_suffix is set.
*/

int
create_file_with_appropriate_name(const const_IWSubstring & old_name,
                                   IWString & new_name,
                                   FileType file_type,
                                   int keep_existing_suffix)
{
  const char * new_suffix = suffix_for_file_type(file_type);
  if (nullptr == new_suffix) {
    cerr << "create_file_with_appropriate_name: unrecognised type " << 
            file_type << endl;

    new_name = "UNK_TYPE";
    return 0;
  }

  IWString dot_suffix;
  dot_suffix << '.' << new_suffix;

  new_name = old_name;

  if (new_name.ends_with(dot_suffix)) {
    return 1;
  }

  if (! keep_existing_suffix) {     // get rid of old suffix
    remove_suffix(new_name);
    new_name += '.';
  } else {
    if (! new_name.ends_with('.'))
      new_name += '.';
  }

  new_name += new_suffix;

  return 1;
}

int
append_appropriate_suffix(IWString & fname, FileType file_type)
{
  const char * new_suffix = suffix_for_file_type(file_type);
  if (nullptr == new_suffix)
  {
    cerr << "append_appropriate_suffix: unrecognised type " << 
            file_type << endl;

    return 0;
  }

  if (! fname.ends_with('.'))
    fname += '.';

  fname += new_suffix;

  return 1;
}

/*
  Tries to identify file type from a file name
  Returns 0 if unsuccessful
*/

FileType
discern_file_type_from_name(const IWString & file_name)
{
  if (0 == file_name.length())
    return FILE_TYPE_INVALID;

  if (file_name.ends_with(".gz"))
  {
    const_IWSubstring tmp(file_name);
    tmp.chop(3);
    return discern_file_type_from_name(tmp);
  }

  if (file_name.ends_with(".smi"))
    return FILE_TYPE_SMI;
  if (file_name.ends_with(".sdf"))
    return FILE_TYPE_SDF;
  if (file_name.ends_with(".mdl"))
    return FILE_TYPE_SDF;
  if (file_name.ends_with(".mol"))
    return FILE_TYPE_SDF;
  if (file_name.ends_with(".pdb"))
    return FILE_TYPE_PDB;
  if (file_name.ends_with(".mmod"))
    return FILE_TYPE_MMOD;
  if (file_name.ends_with(".msi"))
    return FILE_TYPE_MSI;
  if (file_name.ends_with(".tdt"))
    return FILE_TYPE_TDT;
  if (file_name.ends_with(".rdf"))
    return FILE_TYPE_RDF;
  if (file_name.ends_with(".qry"))
    return FILE_TYPE_QRY;
  if (file_name.ends_with(".mol2"))
    return FILE_TYPE_MOL2;
  if (file_name.ends_with(".chm"))
    return FILE_TYPE_CHM;
  if (file_name.ends_with(".moe"))
    return FILE_TYPE_MOE;
  if (file_name.ends_with(".mrk"))
    return FILE_TYPE_MRK;
  if (file_name.ends_with(".wchm"))
    return FILE_TYPE_WCHM;
  if (file_name.ends_with(".gfp"))
    return FILE_TYPE_TDT;
  if (file_name.ends_with(".mrv"))
    return FILE_TYPE_MRV;
  if (file_name.ends_with(".inchi"))
    return FILE_TYPE_INCHI;
  if (file_name.ends_with(".cif"))
    return FILE_TYPE_CIF;
  if (file_name.ends_with(".csv"))
    return FILE_TYPE_CSV;
  if (file_name.ends_with(".textproto") ||
      file_name.ends_with(".txtproto") ||
      file_name.ends_with(".txt"))
    return FILE_TYPE_TXTPROTO;
  if (file_name.ends_with(".xyz"))
    return FILE_TYPE_XYZ;

  return FILE_TYPE_INVALID;
}

static int
process_option(const Command_Line & cl, FileType & type, const char c, const char * type_name)
{

  if (! cl.option_present(c))
  {
    cerr << "Must specify " << type_name << " type via -" << c << " option\n";
    return 0;
  }

  const_IWSubstring type_string;
  int i = 0;
  cl.value(c, type_string, i);

  FileType tmp = moleculeio::string_to_file_type(type_string);
  if (0 == tmp)
  {
    cerr << "Unrecognised " << type_name << " type '" << type_string << "'\n";
    return 0;
  }

  type = tmp;

  return 1;
}

int
process_output_type(const Command_Line & cl, FileType & output_type)
{
  return process_option(cl, output_type, 'o', "output");
}

int
all_files_recognised_by_suffix(const Command_Line & cl)
{
  int rc = 1;

  int nf = cl.number_elements();

  for (int i = 0; i < nf; i++)
  {
    const char * fname = cl[i];

    if (! discern_file_type_from_name(fname))
    {
      cerr << "all_files_recognised_by_suffix: cannot discern file type '" << fname << "'\n";
      rc = 0;
    }
  }

  return rc;
}

int valid_file_type(FileType ftype) {
  return ftype != FILE_TYPE_INVALID;
}

#ifdef NOLONGERUSED_QWE
int
valid_file_type(int ftype)
{
  switch (ftype)
  {
    case MDL:
      return 1;
      break;

    case PDB:
      return 1;
      break;

    case MMOD:
      return 1;
      break;

    case SMI:
      return 1;
      break;

    case USMI:
      return 1;
      break;

    case MSI:
      return 1;
      break;

    case TDT:
      return 1;
      break;

    case UTDT:
      return 1;
      break;

    case RDF:
      return 1;
      break;

    case RSMI:
      return 1;
      break;

    case SDF:
      return 1;
      break;

    case MOL2:
      return 1;
      break;

    case FILE_TYPE_PSF:
      return 1;
      break;

    case FILE_TYPE_CRD:
      return 1;
      break;

    case FILE_TYPE_CHM:
      return 1;
      break;

    case FILE_TYPE_MOE:
      return 1;
      break;

    case FILE_TYPE_MRK:
      return 1;
      break;

    case FILE_TYPE_WCHM:
      return 1;
      break;

    case FILE_TYPE_NAUSMI:
      return 1;
      break;

    case FILE_TYPE_CIF:
      return 1;
      break;

    case FILE_TYPE_SMT:
      return 1;
      break;

    case FILE_TYPE_MRV:
      return 1;
      break;

    case FILE_TYPE_INCHI:
      return 1;
      break;

    case FILE_TYPE_TDT_NAUSMI:
      return 1;
      break;

    case FILE_TYPE_XYZ:
      return 1;
      break;

    default:
      return 0;
      break;
  }
}
#endif


/*
  We have read a molecule with aromatic atoms. Before doing Kekule perception
  we need to remove any covalently bonded non-organics
  Nov 2007. Change this to only un-connect them when they are bonded
  only to heteroatoms
  Oct 2023. More useful to make the default always remove all bonds.
  The old behaviour can be obtained by setting _unconnect_covalently_bonded_non_organics_on_read=2
*/

int
Molecule::_do_unconnect_covalently_bonded_non_organics()
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom * a = _things[i];

    if (a->element()->organic() || 0 == a->ncon())
      continue;

    if (moleculeio::unconnect_covalently_bonded_non_organics_on_read() == 1) {
      remove_bonds_to_atom(1);
      ++rc;
    } else if (attached_heteroatom_count(i) == a->ncon()) {
      remove_bonds_to_atom(i);
      rc++;
    }
  }

//cerr << "Broke bonds to " << rc << " covalent non-organics\n";

  return rc;
}

int
Molecule::read_molecule_ds(iwstring_data_source & input, FileType input_type)
{
  assert (ok());
  assert (input.good());

  int rc = 0;
  if (FILE_TYPE_MDL == input_type || FILE_TYPE_SDF == input_type)
    rc = read_molecule_mdl_ds(input);

  else if (FILE_TYPE_PDB == input_type)
    rc = read_molecule_pdb_ds(input);

  else if (FILE_TYPE_MMOD == input_type)
    rc = read_molecule_mmod_ds(input);

  else if (FILE_TYPE_SMI == input_type || FILE_TYPE_USMI == input_type)  // why would we input USMI?
    rc = read_molecule_smi_ds(input);

  else if (FILE_TYPE_MSI == input_type)
    rc = read_molecule_msi_ds(input);

  else if (FILE_TYPE_TDT == input_type)
    rc = read_molecule_tdt_ds(input);

  else if (FILE_TYPE_RDF == input_type)
    rc = read_molecule_rdf_ds(input);

  else if (FILE_TYPE_MOL2 == input_type)
    rc = read_molecule_mol2_ds(input);

  else if (FILE_TYPE_MOE == input_type)
    rc = read_molecule_moe_ds(input);

  else if (FILE_TYPE_MRK == input_type)
    rc = read_molecule_mrk_ds(input);

  else if (FILE_TYPE_MRV == input_type)
    rc = read_molecule_mrv_ds(input);

  else if (FILE_TYPE_INCHI == input_type)
    rc = read_molecule_inchi_ds(input);

  else if (FILE_TYPE_CIF == input_type)
    rc = read_molecule_cif_ds(input);

  else if (FILE_TYPE_CSV == input_type)
    rc = read_molecule_csv_ds(input);

  else if (FILE_TYPE_TXTPROTO == input_type)
    rc = read_molecule_textproto_ds(input);

  else if (FILE_TYPE_XYZ == input_type)
    rc = read_molecule_xyz_ds(input);

  else
  {
    cerr << "read_molecule_ds: Unknown type " << input_type << "\n";
    iwabort();
  }
  //cerr << "Having read molecule name is " << _molecule_name << "'\n";
  //cerr << Name() << '\n';

  if (rc == 0) {
    return 0;
  }

  if (coordinate_scaling > 0.0f) {
    ScaleCoordinates(coordinate_scaling);
  }

  if (moleculeio::unconnect_covalently_bonded_non_organics_on_read()) {
    _do_unconnect_covalently_bonded_non_organics();
  }

  return rc;
}

static FileType type_to_write_with_operator = FILE_TYPE_SMI;

void
set_type_to_write_with_operator(FileType s)
{
  type_to_write_with_operator = s;
}

std::ostream &
operator << (std::ostream & os, Molecule & m)
{
  m.write_molecule(os, type_to_write_with_operator);

  return os;
}

int
Molecule::write_molecule(std::ostream & os, FileType output_type, 
                         const IWString & comments)
{
  assert (ok());
  assert (os.good());

  int rc = 0;
  if (FILE_TYPE_MDL == output_type || FILE_TYPE_SDF == output_type)
    rc = write_molecule_mdl(os, comments);
  else if (FILE_TYPE_PDB == output_type)
    rc = write_molecule_pdb(os, comments);
  else if (FILE_TYPE_MMOD == output_type)
    rc = write_molecule_mmod(os);    // seems that no function defined with a comments field?
  else if (FILE_TYPE_SMI == output_type)
    rc = write_molecule_smi(os, name());
  else if (FILE_TYPE_USMI == output_type)
    rc = write_molecule_usmi(os, name());
  else if (FILE_TYPE_MSI == output_type)
    rc = write_molecule_msi(os, name());
  else if (FILE_TYPE_TDT == output_type)
    rc = write_molecule_tdt(os, name());
  else if (FILE_TYPE_UTDT == output_type)
    rc = write_molecule_tdt_unique(os, name());
  else if (FILE_TYPE_TDT_NAUSMI == output_type)
    rc = write_molecule_tdt_nausmi(os, name());
  else if (FILE_TYPE_RSMI == output_type)
    rc = write_molecule_rsmi(os, name());
  else if (FILE_TYPE_MOL2 == output_type)
    rc = write_molecule_mol2(os);
  else if (FILE_TYPE_PSF == output_type)
    rc = write_molecule_psf(os);
  else if (FILE_TYPE_CRD == output_type)
    rc = write_molecule_crd(os);
  else if (FILE_TYPE_MOE == output_type) {
    cerr << "Molecule::write_molecule:cannot write MOE type\n";
    rc = 0;
  }
  else if (FILE_TYPE_MRK == output_type)
    rc = write_molecule_mrk(os);
  else if (FILE_TYPE_WCHM == output_type)
    rc = write_molecule_wchm(os);
  else if (FILE_TYPE_NAUSMI == output_type)
    rc = write_molecule_nausmi(os, name());
  else if (FILE_TYPE_CIF == output_type)
    rc = write_molecule_cif(os);
  else if (FILE_TYPE_SMT == output_type)
    rc = write_molecule_smarts(os);
  else if (FILE_TYPE_MRV == output_type)
    rc = write_molecule_mrv(os);
  else if (FILE_TYPE_INCHI == output_type)
    rc = write_molecule_inchi(os);
  else if (FILE_TYPE_CSV == output_type)
    rc = write_molecule_csv(os);
  else if (FILE_TYPE_TXTPROTO == output_type)
    rc = write_molecule_textproto(os);
  else if (FILE_TYPE_XYZ == output_type)
    rc = write_molecule_xyz(os);
  else
    cerr << "Molecule::write_molecule: unrecognised type " << output_type << "\n";
    
  return rc;
}

template <typename T>
int
rwmolecule_error(const char * message, T & input)
{
  assert (nullptr != message);

  cerr << message << ", line " << input.lines_read() << "\n";
  if (input.at_eof())
    return 0;

  IWString buffer;
  input.most_recent_record(buffer);
  cerr << "Buffer '" << buffer << "'\n";

  return 0;
}

template int rwmolecule_error<String_Data_Source>(char const*, String_Data_Source&);
template int rwmolecule_error<iwstring_data_source>(char const*, iwstring_data_source&);

static void
display_input_help(std::ostream & os)
{
  os << R"( -i sdf                  SDF input\n";
 -i smi                  smiles input
 -i tdt                  TDT input
 -i mdl                  MDL format (generally, use 'sdf' instead)
 -i textproto            certain text proto file formats, see 'prototag' below
 -i ICTE=<nn>            ignore as many as <nn> connection table errors
 -i info                 collect extra text records in input (SDF and TDT only)
 -i skip=nn              skip over the first MM molecules in the input
 -i do=nn                only process NN molecules
 -i seek=offset          seek to byte offset OFFSET before starting reading
 -i stop=offset          stop reading once file is at byte offset OFFSET
 -i maxq=<charge>        set maximum plausible atomic partial charge
 -i minq=<charge>        set minimum plausible atomic partial charge
 -i mq=<charge>          set min (-charge) and max (+charge) plausible atomic partial charge
 -i mfc=<charge>         set min and max plausible formal charges
 -i d@3d                 Discern chirality from 3d coordinates
 -i ignore_bad_chiral    ignore obviously incorrect chirality specifications
 -i discard_chiral       discard all chiral information input
 -i mdlRincch            read incorrect chirality flags - explicit Hydrogen problem
 -i addmih               add implicit hydrogens to H deficient chiral centres
 -i dctb                 Discern Cis-Trans Bonds from geometry
 -i dwedge               Discern chirality from Wedge bonds only
 -i xctb                 Discard all cis-trans information on input
 -i ibctb                Ignore erroneous cis-trans bond input
 -i ignore_self_bonds    ignore bonds with the same atom at each end
 -i uccvno               immediately unconnect covalently bonded non organics on input
 -i Hiso                 allow implicit Hydrogens on aromatic isotopically labeled atoms [9c]
 -i rmhknown             remove implicit Hydrogens known property if possible
 -i RMHKNOWN             remove ALL implicit Hydrogens known attributes
 -i fixnd3v4             put formal charge on 3 connected, 4 bonded neutral Nitrogen
 -i mscale=<scale>       multiply all coordinates by <scale> upon reading
 -i mdlD                 allow D to be recognised as [2H]
 -i mdlT                 allow T to be recognised as [3H]
 -i mdl3chop             chop long elements in MDL files to first 2 chars
 -i mdl3=XX              change all long element symbols to XX
 -i ignore_bad_m         ignore unrecognised 'M' records in SDF files
 -i ignore_fatal_m       ignore otherwise fatal errors in M records
 -i mdlquiet             don't report unrecognised records to stdout
 -i SDFID:XX             name follows '> <XX>' record in SDF file
 -i SDFNONAME            discard name in first record of SD file
 -i SDFNOPREPEND         do not prepend the SDF identifier to the identifier value
 -i SDFTAG2JSON          append all sdf tags as JSON form\
 -i NAME2JSON            when writing in JSON form, convert the name to 'name: "...."'
 -i firstsdftag          take first tag in an sd file as the name
 -i allsdfid             concatenate all sdf identifiers
 -i ISISEXTREG           try to discern the ISIS external registry number in SDF files
 -i mdlatomalias         replace elements by their alias symbols, A records in MDL files
 -i Galias               replace elements by their alias symbols, G records in MDL files
 -i gsubsdf=<c>          gsub all spaced values in sdf data records with <c>
 -i mdlRisonum           read isotopes as numbers rather than diffs from normal mass
 -i MDLIBT:num1=num2     translate input bond type 'num1' to 'num2'
 -i mdlustere            accumulate unclassified MDL chiral centres (type 3)
 -i mdlsep=<..>          separator between tags when reading mdl files 
 -i mdlnwce              do NOT write messages about chirality problems in sdf files
 -i sasge                MDL V30: convert single atom SGROUP labels to elements
 -i RDFID:XX             detect RDFILE info in tag <XX>
 -i RDFSTART:XX          records in RDFILES start with tag <XX>
 -i mol2fc               try to assign formal changes when reading MOL2 files
 -i mol2rfc              interpret mol2 charges as formal charges
 -i TDTSMI:XX            smiles is in tag XX rather than $SMI<>
 -i TDTID:XX             identifier is XX rather than PCN
 -i TDTAPPEND:XX         append dataitem XX to name
 -i TDTINCTAG            when appending TDT items, keep tags
 -i ignoretdtnosmi       ignore TDT's with no '$SMI' dataitem
 -i DOS                  remove '^M' characters from the end of each record
 -i delim=<char>         record delimiter (default '\\n', use ^M for DOS)
 -i chemaxon             parse 'smiles |chemaxon| name' forms
 -i squoted              smiles records may have double quotes (typically \"smiles |chemaxon|\", name)
 -i csvsmicol=<col>      when reading csv files, smiles is in column number <col>
 -i csvsminam=<col>      when reading csv files, smiles is in column name <col>
 -i csvidcol=<col>       when reading csv files, identifier is in column number <col>
 -i csvidnam=<col>       when reading csv files, identifier is in column name <col>
 -i csvallcol            when reading csv files, all extra columns become the name
 -i prototag=<tag>       when reading textptoto files, the smiles tag
 -i MOEFRAGS=<fname>     when reading Moe files, write named fragments to <fname>
)";

  // Unimplemented or deprecated features. TODO:ianwatson cleanup/implement...
  //os << " -i mdlRisoMnum          read isotopes as numbers rather than diffs from normal mass
  // Feature seems unimplemented.
  //os << " -i IDM:TAG              the name is in the M TAG record\n";
  // Feature seems unimplemented
  // os << " -i RPSDFTAG=XX          replace the first sdf tag with <XX>\n";
  // Very strange option, remove from help message.
  // os << " -i d@3d=nn              Discern chirality from 3d coordinates, only if < nn atoms\n";

  exit(0);
}

int
process_input_type(const Command_Line & cl, FileType & input_type)
{
  MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();

  input_type = FILE_TYPE_INVALID;

  int i = 0;
  const_IWSubstring optval;
  while (cl.value('i', optval, i++))
  {
//  cerr << "Parsing token '" << optval << "'\n";

    if ("info" == optval)
    {
      moleculeio::set_read_extra_text_info(1);
    }
    else if ("ignore_bad_m" == optval)
    {
      mdlfos->set_ignore_unrecognised_mdl_m_records(1);
    }
    else if ("ignore_fatal_m" == optval)
    {
       mdlfos->set_die_on_erroneous_m_input(0);
    }
    else if ("mdlquiet" == optval)
    {
      mdlfos->set_report_unrecognised_records(0);
    }
//  else if (optval.starts_with("mdlquiet="))     seems like a good idea, implement it sometime
//  {
//    optval.remove_leading_chars(9);
//    mdlfor->set_ignore_unrecognised_mdl_tag(optval);
//  }
    else if ("ignore_self_bonds" == optval)
    {
      mdlfos->set_ignore_self_bonds(1);
      set_add_same_bond_twice_fatal(0);
    }
    else if ("uccvno" == optval)
    {
      moleculeio::set_unconnect_covalently_bonded_non_organics_on_read(1);
    }
    else if ("Hiso" == optval)
    {
      set_add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens(1);
    }
    else if (optval.starts_with("MDLIBT:"))
    {
      optval.remove_leading_chars(7);
      if (! mdlfos->process_mdl_bond_translation(optval))
      {
        cerr << "Invalid MDL input bond type translation '" << optval << "'\n";
        return 0;
      }
    }
    else if (optval.starts_with("SDFID:"))
    {
      IWString sdfid = optval.substr(6);

      if (sdfid.empty()) {
        cerr << "The 'SDFID:' specifier must be followed by the SDF identifier\n";
        return 0;
      }

      sdfid.unhtml();

      if (! mdlfos->set_sdf_identifier(sdfid))    // bad pattern
        return 0;
    }
    else if (optval == "SDFTAG2JSON") {
      mdlfos->set_sdf_tags_to_json(1);
    }
    else if (optval == "NAME2JSON") {
      mdlfos->set_name_to_json(1);
    } else if ("ICTE" == optval)
    {
      moleculeio::set_number_connection_table_errors_to_skip(std::numeric_limits<int>::max());
    }
    else if (optval.starts_with("ICTE="))
    {
      optval.remove_leading_chars(5);
      int t;
      if (! optval.numeric_value(t) || t < 0)
      {
        cerr << "Invalid ICTE qualifier '" << optval << "'\n";
        return 0;
      }
      moleculeio::set_number_connection_table_errors_to_skip(t);
    }
    else if ("allsdfid" == optval || "ALLSDFID" == optval)
    {
      mdlfos->set_fetch_all_sdf_identifiers(1);
    }
    else if (optval.starts_with("gsubsdf="))
    {
      optval.remove_leading_chars(8);
      mdlfos->set_gsub_mdl_file_data(optval[0]);
    }
    else if ("SDFNONAME" == optval)
    {
      mdlfos->set_discard_sdf_molecule_name(1);
    }
    else if ("SDFNOPREPEND" == optval)
    {
      mdlfos->set_prepend_sdfid(0);
    }
    else if ("firstsdftag" == optval)
    {
      mdlfos->set_take_first_tag_as_name(1);
    }
    else if (optval.starts_with("RPSDFTAG="))
    {
      optval.remove_leading_chars(9);
      mdlfos->set_replace_first_sdf_tag(optval);
    }
    else if (optval.starts_with("RDFID:"))
    {
      optval.remove_leading_chars(6);
      mdlfos->add_rdfile_identifier(optval);
    }
    else if (optval.starts_with("RDFSTART:"))
    {
      optval.remove_leading_chars(9);
      mdlfos->set_rdfile_start_of_record(optval);
    }
    else if ("mdlatomalias" == optval)
    {
      mdlfos->set_set_elements_based_on_atom_aliases(1);
      set_auto_create_new_elements(1);
    }
    else if ("Galias" == optval)
    {
      mdlfos->set_mdl_g_records_hold_atom_symbols(1);
      set_auto_create_new_elements(1);
    }
    else if ("ignoretdtnosmi" == optval || "IGNORETDTNOSMI" == optval)
    {
      set_ignore_tdts_with_no_smiles(1);
    }
    else if ("rmhknown" == optval)
    {
      set_unset_implicit_hydrogens_known_if_possible(1);
    }
    else if ("RMHKNOWN" == optval)
    {
      set_unset_all_implicit_hydrogens_known_attributes(1);
    }
    else if ("ISISEXTREG" == optval)
    {
      mdlfos->set_extract_isis_extregno(1);
    }
    else if (optval.starts_with("IDM:"))
    {
      optval.remove_leading_chars(4);
      mdlfos->set_mdl_name_in_m_tag(optval);
    }
    else if ("mdl3chop" == optval)
    {
      mdlfos->set_truncate_long_symbols(1);
    }
    else if (optval.starts_with("mdl3="))
    {
      optval.remove_leading_chars(5);
      mdlfos->set_mdl_change_long_symbols_to(optval);
    }
    else if ("mdlD" == optval)
    {
      mdlfos->set_allow_deuterium(1);
      element::set_interpret_d_as_deuterium(1);
    }
    else if ("mdlT" == optval)
    {
      mdlfos->set_allow_tritium(1);
      element::set_interpret_t_as_tritium(1);
    }
    else if ("mdlRisonum" == optval)
    {
      mdlfos->set_read_isotopes_as_numbers_rather_than_differences_from_normal(1);
    }
    else if ("mdlRisoMnum" == optval)
    {
      mdlfos->set_read_M_isotopes_as_numbers_rather_than_differences_from_normal(1);
    }
    else if ("mdlRincch" == optval)
    {
      mdlfos->set_mdl_read_h_correct_chiral_centres(0);
    }
    else if ("ignore_bad_chiral" == optval)
    {
      moleculeio::set_ignore_incorrect_chiral_input(1);
    }
    else if ("discard_chiral" == optval)
    {
      moleculeio::set_ignore_all_chiral_information_on_input(1);
    }
    else if ("addmih" == optval)
    {
      set_automatically_add_implicit_hydrogen_to_incomplete_chiral_centre(1);
    }
    else if (optval.starts_with("TDTID:"))
    {
      const const_IWSubstring tdtid = optval.substr(6);

      set_tdt_identifier_dataitem(tdtid);
    }
    else if (optval.starts_with("TDTSMI:"))
    {
      const const_IWSubstring tdtsmi = optval.substr(7);

      set_smiles_tag(tdtsmi);
    }
    else if (optval.starts_with("TDTAPPEND:"))
    {
      const const_IWSubstring a = optval.substr(10);

      set_tdt_append_dataitem(a);
    }
    else if ("TDTINCTAG" == optval)
    {
      set_tdt_append_dataitem_content(0);
    }
    else if ("mdlustere" == optval)
    {
      mdlfos->set_accumulate_mdl_chirality_features(1);
    }
    else if ("dctb" == optval)
    {
      moleculeio::set_discern_cis_trans_bonds(1);
    }
    else if ("d@3d" == optval || "d3d" == optval)
    {
      moleculeio::set_discern_chirality_from_3d_coordinates(1);
    }
    else if (optval.starts_with("d@3d=") || optval.starts_with("d3d="))
    {
      optval.remove_up_to_first('=');
      int d;
      if (! optval.numeric_value(d) || d < 3)
      {
        cerr << "The 'd@3d=' qualifier must be a valid atom count\n";
        return 0;
      }

      moleculeio::set_discern_chirality_from_3d_coordinates(d);
    }
    else if ("dwedge" == optval)
    {
      mdlfos->set_discern_chirality_from_wedge_bonds(1);
    }
    else if ("ibctb" == optval)
    {
      moleculeio::set_ignore_bad_cis_trans_input(1);
    }
    else if ("xctb" == optval || "rmctb" == optval)
    {
      moleculeio::set_discard_directional_bonds_on_input(1);
    }
    else if (optval.starts_with("delim="))
    {
      if (7 != optval.length())
      {
        cerr << "The delim= qualifier must be followed by a single character\n";
        return 0;
      }

      moleculeio::set_record_delimiter(optval.last_item());
    }
    else if ("DOS" == optval || "dos" == optval)
    {
      moleculeio::set_dos_mode(1);
    }
    else if (optval.starts_with("skip="))
    {
      optval.remove_leading_chars(5);
      int n;
      if (! optval.numeric_value(n) || n < 0)
      {
        cerr << "The skip first molecules directive 'skip=nn' must be followed by a whole positive number\n";
        return 0;
      }
      moleculeio::set_skip_first_molecules(n);
    }
    else if (optval.starts_with("do="))
    {
      int n;
      optval.remove_leading_chars(3);
      if (! optval.numeric_value(n) || n < 0)
      {
        cerr << "The do only molecules directive 'do=nn' must be followed by a whole positive number\n";
        return 0;
      }
      moleculeio::set_do_only_n_molecules(n);
    }
    else if (optval.starts_with("seek="))
    {
      optval.remove_leading_chars(5);

      off_t tmp;
      if (! optval.numeric_value(tmp)) {
        cerr << "Invalid seek to specifier 'seek=" << optval << "'\n";
        return 0;
      }

      moleculeio::set_seek_to(tmp);
    }
    else if (optval.starts_with("stop="))
    {
      optval.remove_leading_chars(5);

      off_t tmp;
      if (! optval.numeric_value(tmp)) {
        cerr << "Invalid stop specifier 'stop=" << optval << "'\n";
        return 0;
      }

      moleculeio::set_max_offset_from_command_line(tmp);
    }
    else if (optval.starts_with("maxq="))
    {
      optval.remove_leading_chars(5);
      charge_t q;
      if (! optval.numeric_value(q))
      {
        cerr << "INvalid maximum charge value '" << optval << "'\n";
        return 0;
      }

      if (! set_min_reasonble_atomic_partial_charge_value(q))
        return 0;
    }
    else if (optval.starts_with("minq="))
    {
      optval.remove_leading_chars(5);
      charge_t q;
      if (! optval.numeric_value(q))
      {
        cerr << "Invalid minimum charge value '" << optval << "'\n";
        return 0;
      }

      if (! set_min_reasonble_atomic_partial_charge_value(q))
        return 0;
    }
    else if (optval.starts_with("mq="))
    {
      optval.remove_leading_chars(3);
      charge_t q;
      if (! optval.numeric_value(q))
      {
        cerr << "Invalid minimum charge value '" << optval << "'\n";
        return 0;
      }

      if (! set_reasonable_atomic_partial_charge_range(-q, q))
        return 0;
    }
    else if (optval.starts_with("mfc="))
    {
      optval.remove_leading_chars(4);
      formal_charge_t q;
      if (! optval.numeric_value(q) || q < 0)
      {
        cerr << "Invalid max formal charge specification '" << optval << "'\n";
        return 0;
      }

      set_reasonable_formal_charge_range(-q, q);
    }
    else if ("mol2fc" == optval)
    {
      tripos::set_mol2_assign_default_formal_charges(1);
    }
    else if ("mol2rfc" == optval)
    {
      tripos::set_mol2_read_charge_column_contains_formal_charges(1);
    }
    else if ("fixnd3v4" == optval)
    {
       set_put_formal_charges_on_neutral_ND3v4(1);
    }
    else if (optval.starts_with("csvsmicol=")) {
      optval.remove_leading_chars(10);
      int c;
      if (! optval.numeric_value(c) || c < 1) {
        cerr << "process_input_type:invalid csv smiles column " << c << '\n';
        return 0;
      }
      if (! lillymol_csv::set_smiles_column(c - 1)) {
        return 0;
      }
    }
    else if (optval.starts_with("csvidcol=")) {
      optval.remove_leading_chars(9);
      int c;
      if (! optval.numeric_value(c) || c < 1) {
        cerr << "process_input_type:invalid csv id column " << c << '\n';
        return 0;
      }
      if (! lillymol_csv::set_id_column(c - 1)) {
        return 0;
      }
    }
    else if (optval.starts_with("csvsminam=")) {
      optval.remove_leading_chars(10);
      lillymol_csv::set_smiles_column_name(optval);
    }
    else if (optval.starts_with("csvidnam=")) {
      optval.remove_leading_chars(9);
      lillymol_csv::set_id_column_name(optval);
    }
    else if (optval == "csvallcol") {
      lillymol_csv::set_all_tokens_to_name(1);
    }
    else if (optval.starts_with("prototag=")) {
      optval.remove_leading_chars(9);
      lillymol_textproto::set_smiles_tag(optval);
    }
    else if (optval.starts_with("MOEFRAGS=")) {
      optval.remove_leading_chars(9);
      if (! moe::SetupMoeFragmentStream(optval)) {
        cerr << "Cannot initialise fragment stream when reading Moe '" << optval << "'\n";
        return 0;
      }
    }
    else if (optval.starts_with("mdlsep=")) {
      optval.remove_leading_chars(7);
      mdlfos->set_mdl_insert_between_sdf_name_tokens(optval);
    } else if (optval == "mdlnwce") {
      mdlfos->set_mdl_display_invalid_chiral_connectivity(0);
    } else if ("sasge" == optval) {
      mdlfos->set_convert_single_atom_sgroup_to_element(1); 
    } else if (optval == "chemaxon") {
      smiles::SetDiscernChemaxonSmilesExtensions(1);
    } else if (optval == "squoted") {
      smiles::SetProcessQuotedSmiles(1);
    } else if (optval.starts_with("mscale=")) {
      optval.remove_leading_chars(7);
      if (! optval.numeric_value(coordinate_scaling) || coordinate_scaling < 0.0f) {
        cerr << "process_input_type:invalid mscale directive " << optval << '\n';
        return 0;
      }
    } else if ("help" == optval) {
      display_input_help(cerr);
    } else if (input_type) {
      cerr << "Only one -i option is allowed, '" << optval << "' unrecognised\n";
      return 0;
    } else if (0 == (input_type = moleculeio::string_to_file_type(optval))) {
      cerr << "Unrecognised input type or directive '" << optval << "'\n";
      return 0;
    }
  }

  if (0 != input_type)
    ;
  else if (! all_files_recognised_by_suffix(cl))
  {
    cerr << "Unrecognised input type(s)\n";
    return 0;
  }

  return 1;
}

void
reset_rwmolecule_file_scope_variables()
{
//  _ignore_all_chiral_information_on_input = 0;
//  _ignore_incorrect_chiral_input = 0;
//  _flush_files_after_writing_each_molecule = 0;
//  _read_extra_text_info = 0;
//  _discern_cis_trans_bonds = 0;
//  _discern_chirality_from_3d_coordinates = 0;
//  _ignore_bad_cis_trans_input = 0;
//  _write_extra_text_info = 0;
//    record_delimiter = '\n';
//    dos_mode = 1;    // Mar 2005. Change to default
//  _skip_first_molecules = 0;
//  _do_only_n_molecules = 0;
//  _seek_to_from_command_line = 0;
//  _max_offset_from_command_line = std::numeric_limits<off_t>::max();
//  _number_connection_table_errors_to_skip = 0;
//  _unconnect_covalently_bonded_non_organics_on_read = 0;
    _put_formal_charges_on_neutral_ND3v4 = 0;
//    file_scope_newline_string = '\n';
//    _write_DOS_records = 0;
    type_to_write_with_operator = FILE_TYPE_SMI;

    return;
}
