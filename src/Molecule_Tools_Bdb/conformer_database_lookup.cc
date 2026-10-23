// Looks up molecules in BerkeleyDB databases built by conformer_database_load.

#include "sys/stat.h"
#include "sys/types.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include "db_cxx.h"

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mdl.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools_Bdb/conformer_database.pb.h"
#else
#include "conformer_database.pb.h"
#endif


namespace conformer_database {

using std::cerr;

void 
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << R"(Looks up molecules in a BerkeleyDB database built by conformer_database_store.
conformer_database_lookup -d /path/to/database.bdb -o ISIS -S - file.smi > conformer.sdf

 -d <dbname>            name of database built by conformer_database_store.
 -U <fname>             write molecules not found to <fname>.smi.
 -g ...                 chemical standardisation applied to the lookup key.
                          conformers retain any explicit Hydrogen atoms stored.
 -v                     verbose output
)";
  // clang-format on
  // clang-format off

  ::exit(rc);
}

class ConformerDB{
  private:
    std::unique_ptr<Db> _db;

    uint64_t _molecules_read;
    uint64_t _molecules_found;

    int _max_conformers;

    extending_resizable_array<uint32_t> _number_conformers;

    Chemical_Standardisation _chemical_standardisation;
    int _remove_chirality;
    int _reduce_to_largest_fragment;

    IWString_and_File_Descriptor _stream_for_not_found;

    int _verbose;

    // Private functions

    int HandleMoleculeNotFound(Molecule& m);

  public:
    ConformerDB();
    ~ConformerDB();

    int WriteConformers(Molecule& m, const conformer_database::Conformers& conformers);
    int AddConformer(Molecule& m,
        const std::unique_ptr<int[]>& atom_order_in_smiles,
        conformer_database::Conformers& proto);

    int Initialise(Command_Line& cl);

    int Process(Molecule& m, Molecule_Output_Object& output);

    int Report(std::ostream& output) const;
};

ConformerDB::ConformerDB() {
  _verbose = 0;

  _molecules_read = 0;
  _molecules_found = 0;

  _max_conformers = std::numeric_limits<int>::max();

  _remove_chirality = 0;
  _reduce_to_largest_fragment = 0;
}

ConformerDB::~ConformerDB() {
}

int
ConformerDB::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (! cl.option_present('d')) {
    cerr << "Must specify database via the -d option\n";
    Usage(1);
  }

  if (cl.option_present('d')) {
    const char *dbname = cl.option_value('d');

    int flags;
    DBTYPE dbtype;
    int mode;

    dbtype = DB_UNKNOWN;
    flags = 0;
    mode = 0;

    _db.reset(new Db(NULL, DB_CXX_NO_EXCEPTIONS));

    if (int rc = _db->open(NULL, dbname, NULL, dbtype, flags, mode); rc != 0) {
      cerr << "ConformerDB::Initialise:cannot open database '" << dbname << "'\n";
      _db->err(rc, "");
      _db.reset();
      return 0;
    }

    if (_verbose) {
      cerr << "Opened " << dbname << "'\n";
    }
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "ConformerDB::Initialise:cannot initialise chemical standardisation (-g)\n";
      return 0;
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from all molecules\n";
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to the largest fragment\n";
    }
  }

  if (cl.option_present('U')) {
    IWString fname = cl.string_value('U');
    fname.EnsureEndsWith(".smi");
    if (! _stream_for_not_found.open(fname)) {
      cerr << "ConformerDB::initialise:cannot open stream for molecules not found '" << fname << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Molecules not found in db written to '" << fname << "'\n";
    }
  }

  return 1;
}

int
ConformerDB::Process(Molecule& m,
                Molecule_Output_Object& output) {
  ++_molecules_read;

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }
  
  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  lillymol::set_include_coordinates_with_smiles(0);
  const IWString& usmi = m.unique_smiles();
  lillymol::set_include_coordinates_with_smiles(1);

  Dbt zkey((void*)usmi.data(), usmi.length());
  Dbt zdata;

  if (int rc = _db->get(NULL, &zkey, &zdata, 0); rc == 0) {
    ++_molecules_found;
  } else if (rc == DB_NOTFOUND) {
    return HandleMoleculeNotFound(m);
  } else {
    _db->err(rc, "Unspecified error on lookup");
    return 0;
  }

  conformer_database::Conformers proto;
  if (! proto.ParseFromArray(zdata.get_data(), zdata.get_size())) {
    cerr << "Invalid proto " << m.name() << '\n';
    cerr << usmi << '\n';
    return 0;
  }

  int nwrite = 0;
  if (proto.conformers_size() > _max_conformers) {
    nwrite = _max_conformers;
  } else {
    nwrite = proto.conformers_size();
  }

  ++_number_conformers[nwrite];

  Molecule fromdb;
  if (! fromdb.build_from_smiles(proto.smiles())) {
    cerr << "ConformerDB::Process:cannot parse database smiles\n";
    cerr << proto.smiles() << '\n';
    cerr << m.smiles() << " lookup key " << m.name() << '\n';
    return 0;
  }

  const int matoms = fromdb.natoms();

  int ndx = 0;
  for (int i = 0; i < nwrite; ++i) {
    const conformer_database::Conformer& c = proto.conformers(i);
    for (int j = 0; j < matoms; ++j) {
      fromdb.setxyz(j, c.xyz(ndx), c.xyz(ndx + 1), c.xyz(ndx + 2));
      ndx += 3;
    }
    fromdb.invalidate_smiles();
    output.write(fromdb);
  }

  return 1;
}

int
ConformerDB::HandleMoleculeNotFound(Molecule& m) {
  if (_verbose > 1) {
    cerr << m.name() << " not found\n";
  }

  if (_stream_for_not_found.is_open()) {
    _stream_for_not_found << m.smiles() << ' ' << m.name() << '\n';
    _stream_for_not_found.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
ConformerDB::Report(std::ostream& output) const {
  output << "ConformerDB::Report:read " << _molecules_read << " molecules\n";
  output << "Found " << _molecules_found << '\n';
  for (int i = 0; i < _number_conformers.number_elements(); ++i) {
    if (_number_conformers[i] > 0) {
      output << _number_conformers[i] << " molecules had " << i << " conformers\n";
    }
  }

  return 1;
}

int
ConformerDatabaseLookup(ConformerDB& cdb,
        data_source_and_type<Molecule>& input,
        Molecule_Output_Object& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (!cdb.Process(*m, output)) {
      cerr << "ConformerDatabaseLookup:error processing " << m->name() << '\n';
      return 0;
    }
  }

  return 1;
}

int
ConformerDatabaseLookup(ConformerDB& cdb,
        const char* fname,
        FileType input_type,
        Molecule_Output_Object& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "ConformerDB:Cannot open '" << fname << "'\n";
    return 0;
  }

  return ConformerDatabaseLookup(cdb, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:g:d:i:o:S:clU:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_present('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(0);
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (!cl.option_present('i') && 1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {
    input_type = FILE_TYPE_SMI;
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(6);
    }
  }

  if (0 == input_type && !all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  ConformerDB cdb;
  if (! cdb.Initialise(cl)) {
    cerr << "Cannot initialise database\n";
    return 1;
  }

  Molecule_Output_Object output;
  if (cl.option_present('o')) {
    MDL_File_Supporting_Material* mdlfos = global_default_MDL_File_Supporting_Material();
    output.add_output_type(FILE_TYPE_SDF);
    mdlfos->set_write_isis_standard(1);
    mdlfos->set_write_mdl_charges_as_m_chg(1);
  } else {
    output.add_output_type(FILE_TYPE_SMI);
    lillymol::set_include_coordinates_with_smiles(1);
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (output.would_overwrite_input_files(cl, fname)) {
      cerr << "ConformerDatabaseLookup:cannot overwrite input files (-S) " << fname << "'\n";
      return 1;
    }
    if (! output.new_stem(fname)) {
      cerr << "ConformerDatabaseLookup:cannot open output '" << fname << "'\n";
      return 1;
    }
  } else {
    output.new_stem("-");
  }

  for (const char* fname : cl) {
    if (! ConformerDatabaseLookup(cdb, fname, input_type, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    output.do_flush();
    cdb.Report(cerr);
  }

  return 0;
}

}  // namespace conformer_database

int
main(int argc, char **argv) {
  int rc = conformer_database::Main(argc, argv);

  return rc;
}
