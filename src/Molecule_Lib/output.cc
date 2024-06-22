#include <stdlib.h>

#include <iostream>

#if (__GNUC_MINOR__ == 95)
#include <strstream>
#else
#include <sstream>
#endif

#include "Foundational/cmdline/cmdline.h"

#include "mdl.h"
#include "molecule.h"
#include "moleculeio.h"
#include "output.h"
#include "smiles.h"

using std::cerr;
using std::endl;

void
Molecule_Output_Object::_default_values() {
  _verbose = 0;
  _molecules_per_file = 0;
  _name_token_for_fname = -1;
  _molecules_written = 0;
  _files_created = 0;

  _use_stdout = 0;

  return;
}

Molecule_Output_Object::Molecule_Output_Object() {
  _default_values();

  return;
}

Molecule_Output_Object::~Molecule_Output_Object() {
  assert(ok());

  if (-1028 == _molecules_per_file) {
    cerr << "Deleting already freed Molecule_Output_Object\n";
  }
  _molecules_per_file = -1028;

  return;
}

int
Molecule_Output_Object::ok() const {
  if (0 == _number_elements) {
    return 1;
  }

  if (_molecules_per_file) {
    ;
  } else if (_output_types.number_elements() != _number_elements) {
    return 0;
  }

  // For now, if we have streams, we cannot also be writing to individual files

  if (_molecules_per_file) {
    return 1;
  }

  return 1;
}

int
Molecule_Output_Object::debug_print(std::ostream& os) const {
  assert(os.good());

  os << "Molecule output object with " << _number_elements << " streams\n";
  if (_output_types.number_elements()) {
    os << "Output types are";
    for (int i = 0; i < _output_types.number_elements(); i++) {
      os << " " << suffix_for_file_type(_output_types[i]);
    }
    os << endl;
  }

  if (_name_token_for_fname >= 0) {
    os << "Name derived from token " << (_name_token_for_fname + 1) << " of the name\n";
  }

  if (_molecules_per_file > 0) {
    os << _molecules_per_file << " molecules per file, stem is '" << _stem << "'\n";
  }

  os << "Have written " << molecules_written() << " molecules\n";

  return 1;
}

int
Molecule_Output_Object::active() const {
  if (_number_elements) {
    return _number_elements;
  }

  if (_use_stdout) {
    return 1;
  }

  if (_name_token_for_fname >= 0) {
    return 1;
  }

  return _molecules_per_file;
}

/*
  Cannot be good() if there are no streams defined
*/

int
Molecule_Output_Object::good() const {
  if (_use_stdout) {
    return std::cout.good();
  }

  if (0 == _number_elements) {
    if (1 == _molecules_per_file) {  // who knows if things are OK or not, no files open
      return 1;
    }

    return 0;
  }

  for (int i = 0; i < _number_elements; i++) {
    if (!_things[i]->good()) {
      return 0;
    }
  }

  return 1;
}

int
Molecule_Output_Object::do_close() {
  resizable_array_p<ofstream_and_type>::resize(0);
  _output_types.resize(0);

  return 1;
}

int
Molecule_Output_Object::do_flush() {
  for (int i = 0; i < _number_elements; i++) {
    _things[i]->flush();
  }

  return 1;
}

/*
  For now, single_file_per_molecule and output to streams are mutually
  exclusive, so if we turn on single_file_per_molecule, we destroy (close)
  any of the streams we had.
*/

int
Molecule_Output_Object::set_molecules_per_file(int newv) {
  _molecules_per_file = newv;

  // if (newv)
  //   return _turn_on_single_file_per_molecule();
  // else
  //   return _turn_off_single_file_per_molecule();

  return 1;
}

// #define DEBUG_NEW_STEM

/*
  A new stem has been specified. Construct all the files with the
  appropriate suffix.
*/

int
Molecule_Output_Object::new_stem(const const_IWSubstring& nstem, int keep_suffix) {
  assert(ok());
  assert(nstem.length());

  resize(0);

#ifdef DEBUG_NEW_STEM
  cerr << "Molecule_Output_Object::new_stem: setting stem to '" << nstem
       << "' molecules per file = " << _molecules_per_file
       << " ntypes = " << _output_types.number_elements() << endl;
#endif

  if (_molecules_per_file) {
    if (_name_token_for_fname >= 0) {
      cerr << "Molecule_Output_Object::new_stem: names derived from molecule names\n";
      return 1;  // we ignore this
    }

    _stem = nstem;
    _molecules_written = 0;
    return 1;
  }

  const int nt = _output_types.number_elements();
  if (nt == 0) {
    cerr << "Molecule_Output_Object::new_stem: no output types defined\n";
    return 0;
  }

  if ('-' == nstem) {
    if (1 != nt) {
      cerr << "Molecule_Output_Object::new_stem: request for stdout, but " << nt
           << " output types\n";
      return 0;
    }

    _use_stdout = 1;

    return _write_any_header_records();
  }

  resize(nt);

  return _open_new_stems(nstem, keep_suffix);
}

int
Molecule_Output_Object::write(Molecule& m) {
  assert(ok());
  assert(m.ok());

  if (_use_stdout) {
    assert(1 == _output_types.number_elements());

    _molecules_written++;

    if (!m.write_molecule(std::cout, _output_types[0])) {
      return 0;
    }

    if (0 == _number_elements) {  // just stdout
      return 1;
    }
  }

  if (0 == _number_elements && 0 == _molecules_per_file && _name_token_for_fname < 0) {
    cerr << "Molecule_Output_Object::write: object has no outputs defined\n";
    return 0;
  }

  int need_new_stem = 0;
  if (_name_token_for_fname >= 0) {
    need_new_stem = 1;
  } else if (0 == _molecules_per_file) {  // no limit
    ;
  } else if (0 == _molecules_written % _molecules_per_file) {
    need_new_stem = 1;
  }

  if (need_new_stem) {
    IWString nstem;      // a new file name stem if needed
    IWString save_stem;  // restore if needed

    if (_name_token_for_fname >= 0) {
      if (!m.name().word(_name_token_for_fname, nstem)) {
        cerr << "Molecule_Output_Object::_write: cannot access " << _name_token_for_fname
             << " token in '" << m.name() << "'\n";
        return 0;
      }
    } else {
      save_stem = _stem;

      nstem = _stem;
      nstem << _files_created;
    }

    if (!_open_new_stems(nstem)) {
      return 0;
    }

    if (save_stem.length()) {
      _stem = save_stem;
    }
  }

  // cerr << "Setup to write " << _number_elements << " instances of '" << m.name() << "
  // to '" << _stem << "'\n";

  int rc = 1;

  for (int i = 0; i < _number_elements; i++) {
    if (!_things[i]->write_molecule(&m)) {
      rc = 0;
    }
  }

  _molecules_written++;

  return rc;
}

int
Molecule_Output_Object::write(Molecule* m) {
  return write(*m);
}

int
Molecule_Output_Object::_open_new_stems(const const_IWSubstring& nstem, int keep_suffix) {
  if (_number_elements) {
    resize_keep_storage(0);
  }

  // cerr << "Opening stem '" << nstem << "'\n";

  int nt = _output_types.number_elements();

  for (int i = 0; i < nt; i++) {
    IWString tmp;
    FileType output_type = _output_types[i];

    if (!create_file_with_appropriate_name(nstem, tmp, output_type, keep_suffix)) {
      cerr << "Cannot process file type " << output_type << endl;
      return 0;
    }

    ofstream_and_type* n = new ofstream_and_type(output_type, tmp);
    if (!n->good()) {
      cerr << "Molecule_Output_Object::new_stem: cannot open '" << tmp << "'\n";
      delete n;
    } else {
      add(n);
    }

#ifdef DEBUG_NEW_STEM
    cerr << "Opened file '" << tmp << endl;
#endif
  }

  _files_created++;

  _stem = nstem;

  _write_any_header_records();

  return _number_elements;
}

// Currently does not respect the output separator setting
// defined in csv.cc
int
Molecule_Output_Object::_write_any_header_records() const {
  for (int i = 0; i < _output_types.number_elements(); ++i) {
    if (_output_types[i] == FILE_TYPE_CSV) {
      if (_use_stdout) {
        std::cout << "SMILES,ID\n";
      } else {
        *_things[i] << "SMILES,ID\n";
      }
    }
  }

  return 1;
}

/*int
Molecule_Output_Object::_write_single_file_per_molecule (Molecule & m)
{
  IWString name_stem;

  if (_name_token_for_fname >= 0)
  {
    if (! m.name().word(_name_token_for_fname, name_stem))
    {
      cerr << "Molecule_Output_Object::_write_single_file_per_molecule: cannot access " <<
_name_token_for_fname << " token in '" << m.name() << "'\n"; return 0;
    }
  }
  else
  {
    name_stem = _stem;
    name_stem << _molecules_written;
  }

  int rc = 0;

  int nt = _output_types.number_elements();

  for (int i = 0; i < nt; i++)
  {
    IWString fname = name_stem;

    int output_type = _output_types[i];

    fname << '.' << suffix_for_file_type(output_type);

    ofstream output(fname.null_terminated_chars(), ios::out);

    if (! output.good())
    {
      cerr << "write_structure_to_new_file: Cannot open '" << fname << "'\n";
      return 0;
    }
  }

  for (int i = 0; i < _output_types.number_elements(); i++)
  {
    if (! m.write_molecule(output, output_type))
      return 0;

    rc++;
  }

  return rc;
}*/

/*
  We must make sure that we don't allow something like '-o smi -o usmi'
  where the two types use the same suffix
*/

int
Molecule_Output_Object::_suffix_already_used(const char* suffix) {
  for (int i = 0; i < _output_types.number_elements(); i++) {
    const char* si = suffix_for_file_type(_output_types[i]);
    if (0 == ::strcmp(suffix, si)) {
      return 1;
    }
  }

  return 0;  // not already used
}

int
Molecule_Output_Object::_add_output_type(FileType otype, const const_IWSubstring& c,
                                         MDL_File_Supporting_Material& mdlfos) {
  if (FILE_TYPE_CHM == otype) {  // need to add two different types
    if (!_add_output_type(FILE_TYPE_PSF, c, mdlfos)) {
      return 0;
    }
    if (!_add_output_type(FILE_TYPE_CRD, c, mdlfos)) {
      return 0;
    }

    return 1;
  }

  if (_output_types.contains(otype)) {
    cerr << "Molecule_Output_Object::determine_output_types: duplicate output type '" << c
         << "' ignored\n";
    return 1;  // not fatal, just ignore these
  }

  const char* suffix = suffix_for_file_type(otype);
  if (_suffix_already_used(suffix)) {
    cerr << "Molecule_Output_Object::determine_output_types: suffix '" << suffix
         << "' already specified\n";
    return 0;
  }

  _output_types.add(otype);
  if (FILE_TYPE_SDF == otype) {
    mdlfos.set_write_mdl_charges_as_m_chg(1);
  }

  return 1;
}

static int
display_output_help_screen(std::ostream& os, char opt) {
  // clang-format off
  os << "The following qualifiers are recognised by the -" << opt << " flag\n";

  os << " -" << opt << " smi3d       append coordinates after smiles\n";
  os << " -" << opt << " flush       flush files after writing each molecule\n";
  os << " -" << opt << " info        write any associated text info (if possible)\n";
  os << " -" << opt << " MULT        create multiple output files, one per molecule\n";
  os << " -" << opt << " MULT=nn     create multiple output files, <nn> molecules per file\n";
  os << " -" << opt << " TOKEN=n     with MULT, use token N from the name as the file name\n";
  os << " -" << opt << " ISIS        create as close as possible to an ISIS file\n";
  os << " -" << opt << " AROM        when writing MDL files, write aromatic bonds as '4'\n";
  os << " -" << opt << " V30         when writing MDL files, write V30 file type\n";
  os << " -" << opt << " DOS         include ^M characters in output\n";
  os << " -" << opt << " mdlWisonum  write isotopes as numbers rather than diffs from normal mass\n";
  os << " -" << opt << " mdlWisoMnum write isotopes as numbers rather than diffs from normal mass\n";
  os << " -" << opt << " mdlWincch   write stereo flags (1,2) without special consideration for\n";
  os << "                explicit Hydrogens - may give incorrect values\n";
  os << " -" << opt << " mdlMEND     always add the M  END record to MDL files\n";
  os << " -" << opt << " RxSymbol    when writing R1 elements to MDL files, do NOT translate to R# form\n";
  os << " -" << opt << " pdbso       write the atoms in PDB files in sequence (fragment) order\n";
  os << " -" << opt << " pdbnws      label atoms in pdb files by number within sequence\n";
  os << " -" << opt << " pdbec       label atoms in pdb files by element C1 C2 C3 N1 O1 O2 C4\n";
  os << " -" << opt << " pdbsname    use stored atom names when writing PDB files\n";
  os << " -" << opt << " mol2wfc     write the formal charge in the charge column of MOL2 files\n";
  os << " -" << opt << " nochiral    exclude chirality info from smiles and mdl outputs\n";
  os << " -" << opt << " nochiralflag don't write the chiral flag info to mdl files\n";
  os << " -" << opt << " NOCT        exclude any CIS/TRANS information from smiles output\n";
  os << " -" << opt << " smisep=<c>  separator between smiles and name: 'smisep=tab' for example\n";
  os << " -" << opt << " <type>      specify one or more output types (smi,usmi,nausmi,rsmi,sdf,tdt,mol2,marvin)\n";
  // clang-format off

  return os.good();
}

int
Molecule_Output_Object::determine_output_types(const Command_Line& cl, char opt)
{
  assert(_output_types.empty());

  MDL_File_Supporting_Material* mdlfos = global_default_MDL_File_Supporting_Material();

  _output_types.resize(cl.option_count(opt));

  int i = 0;
  const_IWSubstring c;
  while (cl.value(opt, c, i++))
  {
    if ("flush" == c)
    {
      moleculeio::set_flush_files_after_writing_each_molecule(1);
      continue;
    }

    if ("info" == c)
    {
      moleculeio::set_write_extra_text_info(1);
      continue;
    }

    if (c.starts_with("MULT=") || c.starts_with("mult="))
    {
      c.remove_leading_chars(5);

      if (! c.numeric_value(_molecules_per_file) || _molecules_per_file < 1)
      {
        cerr << "Invalid molecules per file specification '" << c << "'\n";
        return 0;
      }

      continue;
    }

    if (c.starts_with("MULT") || c.starts_with("mult"))
    {
      _molecules_per_file = 1;
      continue;
    }

    if (c.starts_with("TOKEN=") || c.starts_with("token="))
    {
      c.remove_leading_chars(6);
      if (! c.numeric_value(_name_token_for_fname) || _name_token_for_fname < 1)
      {
        cerr << "The TOKEN= qualifier must be followed by a valid word number\n";
        return 0;
      }

      _name_token_for_fname--;    // convert to C numbering
      continue;
    }

    if (c.starts_with("STEM=") || c.starts_with("stem="))
    {
      c.remove_leading_chars(5);
      _stem = c;
      continue;
    }

    if ("NOCT" == c)
    {
      set_include_cis_trans_in_smiles(0);
      continue;
    }

    if ("ISIS" == c)
    {
      _output_types.add_if_not_already_present(FILE_TYPE_SDF);
      mdlfos->set_write_isis_standard(1);
      mdlfos->set_write_mdl_charges_as_m_chg(1);
      continue;
    }

    if ("AROM" == c)
    {
      mdlfos->set_mdl_write_aromatic_bonds(1);
      continue;
    }

    if ("V30" == c)
    {
      mdlfos->set_write_v30_mdl_files(1);
      _output_types.add_if_not_already_present(FILE_TYPE_SDF);
      continue;
    }

    if ("mdlWisonum" == c)
    {
      mdlfos->set_write_isotopes_as_numbers_rather_than_differences_from_normal(1);
      continue;
    }

    if ("mdlWisoMnum" == c)
    {
      mdlfos->set_write_M_isotopes_as_numbers_rather_than_differences_from_normal(1);
      continue;
    }

    if ("mdlWincch" == c)
    {
      mdlfos->set_mdl_write_h_correct_chiral_centres(0);
      continue;
    }

    if ("mdlMEND" == c)
    {
      mdlfos->set_write_mdl_m_end_record(
          2);    // the 2 indicates high priority, over-rides defaults
      continue;
    }

    if ("RxSymbol" == c)
    {
      mdlfos->set_write_Rn_groups_as_element(1);
      continue;
    }

    if ("nochiral" == c)
    {
      set_include_chiral_info_in_smiles(0);
      mdlfos->set_include_chiral_info_in_mdl_outputs(0);
      continue;
    }

    if ("nochiralflag" == c)
    {
      mdlfos->set_write_mdl_chiral_flags(0);
      continue;
    }

    if ("DOS" == c)
    {
      moleculeio::set_write_DOS_records(1);
      continue;
    }

    if ("pdbso" == c || "pdbfo" == c)
    {
      set_write_pdb_files_in_fragment_order(1);
      continue;
    }

    if ("pdbec" == c)
    {
      set_pdb_number_by_element_count(1);
      continue;
    }

    if ("pdbnws" == c)
    {
      set_pdb_number_within_sequence(1);
      continue;
    }

    if ("pdbsname" == c)
    {
      set_store_pdb_atom_information(1);
      set_use_stored_atom_information_when_writing_pdb_files(1);
      continue;
    }

    if ("smi3d" == c)
    {
      set_append_coordinates_after_each_atom(1);
      continue;
    }

    if ("mol2wfc" == c)
    {
      tripos::set_mol2_write_formal_charge_as_partial_charge(1);
      continue;
    }

    if (c.starts_with("smisep=")) {
      c.remove_leading_chars(7);
      IWString tmp(c);
      if (! char_name_to_char(tmp)) {
        cerr << "Molecule_Output_Object::determine_output_types:invalid smisep '" << c << "'\n";
        return 0;
      }
      smiles::set_smiles_output_separator(tmp[0]);
      continue;
    }

    if ("help" == c)
    {
      display_output_help_screen(cerr, opt);
      exit(1);
    }

    FileType tmp = moleculeio::string_to_file_type(c);
    if (FILE_TYPE_INVALID == tmp)
    {
      cerr << "Molecule_Output_Object::determine_output_types:Unrecognised output type '" << c
           << "'\n";
      return 0;
    }

    if (! _add_output_type(tmp, c, *mdlfos))
      return 0;
  }

  if (_output_types.empty())
    _output_types.add(FILE_TYPE_SMI);

  if (_name_token_for_fname >= 0 && _molecules_per_file <= 0)
    _molecules_per_file = 1;

  return _output_types.number_elements();
}

int
Molecule_Output_Object::add_output_type(FileType ot)
{
  if (_number_elements)
  {
    cerr
        << "Molecule_Output_Object::add_output_type: cannot add types after initialised - fix this some day\n";
    return 0;
  }

  if (_output_types.contains(ot))
  {
    cerr << "Molecule_Output_Object::add_output_type:type " << ot
         << " already specified. Ignored\n";
    return 0;
  }

  const char* suffix = suffix_for_file_type(ot);
  if (nullptr == suffix)
  {
    cerr << "Molecule_Output_Object::add_output_type: unrecognised type " << ot << endl;
    return 0;
  }

  if (_suffix_already_used(suffix))
  {
    cerr << "Molecule_Output_Object::add_output_type: suffix '" << suffix
         << "' already specified\n";
    return 0;
  }

  _output_types.add(ot);

  return _output_types.number_elements();
}

/*
  Is this object writing to a particular file?
*/

int
Molecule_Output_Object::would_use_name(const char* name) const
{
  int no = _output_types.number_elements();

  for (int i = 0; i < no; i++)
  {
    IWString tmp;

    FileType output_type = _output_types[i];
    if (! create_file_with_appropriate_name(name, tmp, output_type))
    {
      cerr << "Cannot process file type " << output_type << endl;
      return 0;
    }

    if (tmp == name)
      return 1;
  }

  return 0;
}

/*
  When using one of these objects, we frequently need to know whether
  or not, given a proposed stem, will it overwrite a given file name.
*/

int
Molecule_Output_Object::would_use_name(const const_IWSubstring& proposed_stem,
                                       const char* name_to_match) const
{
  int no = _output_types.number_elements();

  for (int i = 0; i < no; i++)
  {
    IWString tmp;

    FileType output_type = _output_types[i];
    if (! create_file_with_appropriate_name(proposed_stem, tmp, output_type))
    {
      cerr << "Cannot process file type " << output_type << endl;
      return 0;
    }

    if (tmp == name_to_match)
      return 1;
  }

  return 0;
}

int
Molecule_Output_Object::would_overwrite_input_files(const Command_Line& cl,
                                                    const const_IWSubstring& proposed_stem) const
{
  if (_output_types.empty())
  {
    cerr << "Molecule_Output_Object::would_overwrite_input_files: no output types\n";
    return 0;
  }

  for (int i = 0; i < cl.number_elements(); i++)
  {
    const char* fname = cl[i];
    if (would_use_name(proposed_stem, fname))
      return 1;
  }

  return 0;
}

std::ostream&
Molecule_Output_Object::stream_for_type(int ftype) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (ftype == _output_types[i])
      return *(_things[i]);
  }

  cerr << "Molecule_Output_Object::stream_for_type:no stream for type " << ftype << endl;
  return std::cout;
}
