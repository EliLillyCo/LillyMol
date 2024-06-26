// Do lookups in a structure database.

#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/standardise.h"

#include "structure_database.h"

namespace structure_database {

namespace fs = std::filesystem;

using std::cerr;

// By convention the Usage function tells how to use the tool.
void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << R"(Looks up molecules in a structure database built by structure-database_load.
 -d <dbname>    database name - generated by structure_database_load. 
                Both <dbname>.smi.bdb and <dbname>.graph.bdb database files must be present.
 -l             strip input molecule to the largest fragment before doing the lookup.
 -c             remove chirality from the input molecules before doing database lookup.
 -t             transform the input molecule to graph form before doing the database lookup.
 -s             remove invalid chiral centres before lookup
 -F <file>      stream for molecules found in the database
 -p             append database identifier to name
 -e <sep>       search all databases, insert <sep> string between database finds, def , def |
 -U <file>      stream for molecules NOT found in the database
 -M <file>      stream for molecules with database lookup status
 -W <fname>     write per atom count lookup rates to <fname>
 -v             verbose output
)";
// clang-format on

  ::exit(rc);
}

// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    // If doing a tautomer lookup, needs special handling.
    int _tautomer_lookup = 0;

    int _remove_invalid_chiral_centres = 0;

    // We can have search any number of databases.
    resizable_array_p<StructureDatabase> _database;

    // One or more or'd Lookup enum values.
    uint32_t _mask;

    // When we find a molecule in a database, we can append the database
    // value to the name of the molecule.
    int _append_database_id;

    // If there are multiple databases, by default we stop looking once
    // we have a match, but we can continue searching.
    int _search_all_databases;

    Molecule_Output_Object _stream_for_found;
    Molecule_Output_Object _stream_for_not_in_database;

    Chemical_Standardisation _chemical_standardisation;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    int _molecules_read = 0;
    int _molecules_in_database = 0;

    extending_resizable_array<int> _atoms_in_input;
    extending_resizable_array<int> _atoms_in_found;

    // The -W option.
    IWString _write_per_atom_lookup_rates;

    // The -M option
    IWString_and_File_Descriptor _stream_for_lookup_status;

    // If we are concatenating results from databases, the separator
    // between them. We use the first character only.
    IWString _between_db;

  // private functions
    int CommonOpenFile(const Command_Line& cl, const IWString& stem, 
                Molecule_Output_Object& output);
    int OpenDatabase(IWString& dbname);
    int HandleInDatabase(Molecule& m, const IWString& fromdb);
    int HandleNotInDatabase(Molecule& m);
    int WritePerAtomLookupRates(const char* fname);
    int WritePerAtomLookupRates(IWString_and_File_Descriptor& output);
    int MaybeWriteLookupStatus(const IWString& id, const resizable_array<int>& status);

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    // After each molecule is read, but before any processing
    // is attempted, do any preprocessing transformations.
    int Preprocess(Molecule& m);

    int Process(Molecule& mol);

    // If any summary reports have been requested.
    int WriteSummaryData();

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _remove_invalid_chiral_centres = 0;
  _mask = Lookup::kExact;
  _tautomer_lookup = 0;
  _append_database_id = 0;
  _between_db = '|';
  _molecules_read = 0;
  _molecules_in_database = 0;
  _search_all_databases = 0;
}

std::string
Basename(const IWString& fname) {
  std::string tmp(fname.data(), fname.size());
  return fs::path(tmp).filename();
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
      Usage(8);
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove all chirality\n";
    }
  }

  if (cl.option_present('s')) {
    if (_remove_chirality) {
      cerr << "Using the -s option with the -c option does not make sense\n";
      return 0;
    }
    _remove_invalid_chiral_centres = 1;
    if (_verbose) {
      cerr << "Will remove invalid chiral centres before lookup\n";
    }
  }

  if (cl.option_present('t')) {
    _tautomer_lookup = 1;
    _mask = Lookup::kGraph;
  }

  if (! cl.option_present('d')) {
    cerr << "Must specify one o more databases via the -d option\n";
    Usage(1);
  }

  if (cl.option_present('W')) {
    cl.value('W', _write_per_atom_lookup_rates);
    if (_verbose) {
      cerr << "Per atom lookup rates written to '" << _write_per_atom_lookup_rates << '\n';
    }
  }

  resizable_array<IWString> dbnames;

  if (cl.option_present('d')) {
    IWString fname;
    for (int i = 0; cl.value('d', fname, i); ++i) {
      dbnames << fname;
      if (! OpenDatabase(fname)) {
        cerr << "Cannot open database '" << fname << "'\n";
        return 0;
      }
    }

    if (_verbose) {
      cerr << "Opened " << _database.size() << " structure databases\n";
    }
  } 

  // If we are doing chemical standardisation here, turn off in the databases.
  if (cl.option_present('g')) {
    for (StructureDatabase* d : _database) {
      d->TurnOffChemicalStandardisation();
    }
  }

  if (cl.option_present('e')) {
    _search_all_databases = 1;
    if (_database.size() == 1) {
      cerr << "Searching all databases when there is just one database does not make much sense\n";
    }
    _between_db = cl.option_value('e');
    if (! char_name_to_char(_between_db)) {
      cerr << "Invalid character name '" << _between_db << '\n';
      return 0;
    }
    if (_verbose) {
      cerr << "Will search all databases '" << _between_db << "' inserted between db's\n";
    }
  }

  if (cl.option_present('p')) {
    if (!cl.option_present('F') && !cl.option_present('M')) {
      cerr << "The -p option only makes sense with the -F or -M option\n";
      return 0;
    }

    _append_database_id = 1;
    if (_verbose) {
      cerr << "Will append database identifiers\n";
    }
  }

  if (cl.option_present('M')) {
    IWString fname = cl.option_value('M');
    if (! _stream_for_lookup_status.open(fname.null_terminated_chars())) {
      cerr << "Cannot open stream for lookup status '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Per database lookup status written to '" << fname << "'\n";
    }
    static constexpr char kSep = ' ';

    _stream_for_lookup_status << "Id";
    for (int i = 0; i < _database.number_elements(); ++i) {
      _stream_for_lookup_status << kSep << Basename(dbnames[i]);
    }
    _stream_for_lookup_status << '\n';
  }

  if (cl.option_present('F')) {
    IWString fname = cl.string_value('F');
    if (! CommonOpenFile(cl, fname, _stream_for_found)) {
      cerr << "Cannot open stream for in database (-F) '" << fname << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Molecules found in the database written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('U')) {
    IWString fname = cl.string_value('U');
    if (! CommonOpenFile(cl, fname, _stream_for_not_in_database)) {
      cerr << "Cannot open stream for not in database (-U) '" << fname << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Molecules not found in the database written to '" << fname << "'\n";
    }
  }

  return 1;
}

int
Options::CommonOpenFile(const Command_Line& cl, const IWString& stem, 
                Molecule_Output_Object& output) {
  if (! cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SMI);
  } else if (! output.determine_output_types(cl, 'o')) {
    cerr << "Options::CommonOpenFile:cannot determing output type information\n";
    return 0;
  }

  if (output.would_overwrite_input_files(cl, stem)) {
    cerr << "Options::CommonOpenFile:cannot overwrite input file '" << stem << "'\n";
    return 0;
  }

  return output.new_stem(stem);
}

int
Options::OpenDatabase(IWString& dbname) {
  std::unique_ptr<StructureDatabase> db = std::make_unique<StructureDatabase>();
  if (! db->OpenForReading(dbname)) {
    cerr << "Options::OpenDatabase:cannot open '" << dbname << "'\n";
    return 0;
  }

  _database << db.release();

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << "Found " << _molecules_in_database << " molecules ";
  output << iwmisc::Fraction<float>(_molecules_in_database, _molecules_read);
  output << '\n';

  return 1;
}

int
Options::WriteSummaryData() {
  if (_write_per_atom_lookup_rates.empty()) {
    return 1;
  }

  return WritePerAtomLookupRates(_write_per_atom_lookup_rates);
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_remove_invalid_chiral_centres) {
    lillymol::RemoveInvalidChiralCentresUsingSymmetry(m);
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

int
Options::Process(Molecule& m) {
  ++_molecules_read;

  ++_atoms_in_input[m.natoms()];

  int found_match = 0;

  // A per-database lookup result, for the -M file.
  resizable_array<int> status;

  IWString fromdb;
  for (StructureDatabase* db : _database) {
    IWString thisdb;
    if (! db->Lookup(m, _mask, thisdb)) {
      status << 0;
      continue;
    }

    status << 1;
    if (_append_database_id) {
      fromdb.append_with_spacer(thisdb, _between_db[0]);
    }

    ++found_match;
    if (! _search_all_databases) {
      break;
    }
  }

  MaybeWriteLookupStatus(m.name(), status);

  if (found_match) {
    return HandleInDatabase(m, fromdb);
  }

  return HandleNotInDatabase(m);
}

int
Options::HandleNotInDatabase(Molecule& m) {
  if (_stream_for_not_in_database.active()) {
    return _stream_for_not_in_database.write(m);
  }

  return 1;
}

int
Options::HandleInDatabase(Molecule& m, const IWString& fromdb) {
  ++_molecules_in_database;

  ++_atoms_in_found[m.natoms()];

  if (! _stream_for_found.active()) {
    return 1;
  }

  if (! fromdb.empty()) {
    m << ' ' << fromdb;
  }

  return _stream_for_found.write(m);
}

int
Options::MaybeWriteLookupStatus(const IWString& id, const resizable_array<int>& status) {
  if (! _stream_for_lookup_status.active()) {
    return 1;
  }

  static constexpr char kSep = ' ';

  _stream_for_lookup_status << id;
  for (int i = 0; i < _database.number_elements(); ++i) {
    _stream_for_lookup_status << kSep << '0';
  }
  _stream_for_lookup_status << '\n';

  _stream_for_lookup_status.write_if_buffer_holds_more_than(8192);

  return 1;
}

int
Options::WritePerAtomLookupRates(IWString_and_File_Descriptor& output) {
  static constexpr char kSep = ' ';

  output << "Natoms" << kSep << "NMols" << kSep << "Ratio\n";

  for (int i = 1; i < _atoms_in_input.number_elements(); i++) {
    if (0 == _atoms_in_input[i]) {
      continue;
    }

    int matches;
    if (i >= _atoms_in_found.number_elements()) {
      matches = 0;
    } else {
      matches = _atoms_in_found[i];
    }

    float ratio = iwmisc::Fraction<float>(matches, _atoms_in_input[i]);

    output << i << kSep <<_atoms_in_input[i] << kSep << matches << kSep << ratio << '\n';
  }

  return 1;
}

int
Options::WritePerAtomLookupRates(const char* fname) {
  IWString_and_File_Descriptor output;

  if (!output.open(fname)) {
    cerr << "Cannot open per atom count file name '" << fname << "'\n";
    return 0;
  }

  return WritePerAtomLookupRates(output);
}


int
StructureDatabaseLookup(Options& options,
                Molecule& m) {
  // cerr << "Looking up " << m.unique_smiles() << ' ' << m.name() << '\n';
  return options.Process(m);
}

int
StructureDatabaseLookup(Options& options,
                data_source_and_type<Molecule>& input) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    // cerr << "Just read " << m->unique_smiles() << ' ' << m->name() << '\n';
    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! StructureDatabaseLookup(options, *m)) {
      return 0;
    }
  }

  return 1;
}

int
StructureDatabaseLookup(Options& options,
             const char * fname,
             FileType input_type) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "StructureDatabaseLookup:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return StructureDatabaseLookup(options, input);
}

int
StructureDatabaseLookup(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:d:F:U:pe:M:W:cltT:s");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }
  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }


  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  for (const char * fname : cl) {
    if (! StructureDatabaseLookup(options, fname, input_type)) {
      cerr << "StructureDatabaseLookup::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  options.WriteSummaryData();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace structure_database

int
main(int argc, char ** argv) {

  int rc = structure_database::StructureDatabaseLookup(argc, argv);

  return rc;
}
