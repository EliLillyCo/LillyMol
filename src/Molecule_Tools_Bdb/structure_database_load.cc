// Load a BerkeleyDB database of unique smiles

#include <iostream>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "structure_database.h"

namespace structure_database {

using std::cerr;

int molecules_read = 0;

int remove_invalid_chiral_centres = 0;

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
  cerr << "Loads molecules into a set of structure databases\n";
  cerr << " -d <stem>   name step for the two databases\n";
  cerr << "             <stem>.bdb and <stem>.graph.bdb will be formed\n";
  cerr << " -s          remove invalid chiral centres from input molecules\n";
  cerr << " -v          verbose output\n";
// clang-format on

  ::exit(rc);
}

int
StructureDatabaseLoad(data_source_and_type<Molecule>& input,
                StructureDatabase& database) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    ++molecules_read;

    if (remove_invalid_chiral_centres) {
      //lillymol::RemoveInvalidChiralCentresUsingSymmetry(*m);
      lillymol::do_remove_invalid_chiral_centres(*m);
    }

    if (! database.Store(*m)) {
      return 0;
    }
  }

  return 1;
}

int
StructureDatabaseLoad(const char * fname,
             FileType input_type,
             StructureDatabase& database) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "StructureDatabaseLoad:cannot open '" << fname << "'\n";
    return 0;
  }

  return StructureDatabaseLoad(input, database);
}

int
StructureDatabaseLoad(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:d:s");

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

  if (cl.option_present('s')) {
    remove_invalid_chiral_centres = 1;
    if (verbose) {
      cerr << "Will remove invalid chiral centres\n";
    }
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

  StructureDatabase database;
  if (! cl.option_present('d')) {
    cerr << "Must specify database stem via the -d option\n";
    Usage(1);
  }

  if (cl.option_present('d')) {
    IWString stem = cl.string_value('d');
    IWString dbname;
    dbname << stem << ".smi.bdb";
    if (! database.OpenStructureDatabaseForWriting(dbname)) {
      cerr << "Cannot open '" << dbname << "'\n";
      return 1;
    }

    dbname = stem;
    dbname << ".graph.bdb";
    if (! database.OpenGraphDatabaseForWriting(dbname)) {
      cerr << "Cannot open graph db '" << dbname << "'\n";
      return 1;
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  for (const char * fname : cl) {
    if (! StructureDatabaseLoad(fname, input_type, database)) {
      cerr << "StructureDatabaseLoad::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << "Processed " << molecules_read << " molecules\n";
  }

  return 0;
}

}  // namespace structure_database

int
main(int argc, char ** argv) {

  int rc = structure_database::StructureDatabaseLoad(argc, argv);

  return rc;
}
