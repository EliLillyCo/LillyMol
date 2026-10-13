// loading a BerkeleyDB of conformers. Use conform_database_lookup for retrievals.

#include "sys/stat.h"
#include "sys/types.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>

#include "db_cxx.h"

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iwstring/string_with_commas.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/moleculeio.h"
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

  // clang-format off
  cerr << R"(Builds a berkeleyDB database of conformers.
Generate 3D structures, with the conformers for each molecule adjacent to each other.

conformer_database_store -d /path/to/database.bdb file.sdf

The key to the database is the (Hydrogen suppressed) unique smiles of the molecule, and
the data is a serialized conformer_database.Conformers proto that contains the smiles
with Hydrogens and the coordinates.

 -d <dbname>            name of database to be created.
 -r <n>                 report progress every <n> molecules processed.
 -g ...                 chemical standardisation. Applied to the database key, stored conformers keep H's.
 -v                     verbose output.
)";
  // clang-format on

  ::exit(rc);
}

// If we have specified a tag for conformer energies, each molecule will have
// an associated energy.
class MoleculeWithEnergy : public Molecule {
  private:
    float _energy;
  public:
    MoleculeWithEnergy() {
      _energy = 0.0f;
    }

    float energy() const {
      return _energy;
    }
    void set_energy(float s) {
      _energy = s;
    }
};

class ConformerDB{
  private:
    std::unique_ptr<Db> _db;

    uint64_t _molecules_read;
    uint64_t _duplicates;

    extending_resizable_array<uint32_t> _number_conformers;

    Chemical_Standardisation _chemical_standardisation;

    Report_Progress _report_progress;

    // If the energy of the conformer is in the .sdf file,
    // this is the tag line.
    IWString _energy_tag;

    int _verbose;

    // Private functions

    int NextSetOfConformers(data_source_and_type<MoleculeWithEnergy>& input,
                    resizable_array_p<MoleculeWithEnergy>& molecules);
    int GetEnergy(MoleculeWithEnergy& m);

  public:
    ConformerDB();
    ~ConformerDB();

    int WriteConformers(MoleculeWithEnergy& m, const conformer_database::Conformers& conformers);
    int AddConformer(MoleculeWithEnergy& m,
        const std::unique_ptr<int[]>& atom_order_in_smiles,
        conformer_database::Conformers& proto);

    int Initialise(Command_Line& cl);

    int Process(data_source_and_type<MoleculeWithEnergy>& input);

    int Report(std::ostream& output) const;
};

ConformerDB::ConformerDB() {
  _molecules_read = 0;
  _duplicates = 0;
  _verbose = 0;
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

    dbtype = DB_BTREE;
    flags = DB_CREATE;
    mode = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;

    _db.reset(new Db(NULL, DB_CXX_NO_EXCEPTIONS));

    if (int rc = _db->open(NULL, dbname, NULL, dbtype, flags, mode); rc != 0) {
      cerr << "ConformerDB::Initialise:cannot open database '" << dbname << "'\n";
      _db->err(rc, "");
      _db.reset();
      return 0;
    }

    if (_verbose) {
      cerr << "Smiles will be written to database '" << dbname << "'\n";
    }
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "ConformerDB::Initialise:cannot initialise chemical standardisation (-g)\n";
      return 0;
    }
  }

  if (cl.option_present('r')) {
    if (! _report_progress.initialise(cl, 'r', _verbose)) {
      cerr << "ConformerDB::Initialise:cannot initialise progress reporting (-r)\n";
      return 0;
    }
  }

  if (cl.option_present('G')) {
    IWString g = cl.string_value('G');
    _energy_tag << ">  <" << g << '>';
    if (_verbose) {
      cerr << "Energy in sdf tag '" << _energy_tag << "'\n";
    }
    moleculeio::set_read_extra_text_info(1);
  }

  return 1;
}

int
ConformerDB::AddConformer(MoleculeWithEnergy& m,
        const std::unique_ptr<int[]>& atom_order_in_smiles,
        conformer_database::Conformers& proto) {
  conformer_database::Conformer* c  = proto.mutable_conformers()->Add();
  c->set_energy(m.energy());
  const int matoms = m.natoms();
  c->mutable_xyz()->resize(matoms * 3, 0.0f);
  for (int i = 0; i < matoms; ++i) {
    int j = atom_order_in_smiles[i];
    const Atom& a = m[j];
    c->set_xyz(3 * i, a.x());
    c->set_xyz(3 * i + 1, a.y());
    c->set_xyz(3 * i + 2, a.z());
  }

  return 1;
}

int
Reset(MoleculeWithEnergy& m, std::unique_ptr<int[]>& atom_order_in_smiles,
      conformer_database::Conformers& proto) {
  proto.Clear();
  const IWString& n = m.name();
  proto.set_id(n.rawdata(), n.length());
  const IWString& s = m.smiles();
  proto.set_smiles(s.rawdata(), s.length());

  const resizable_array<atom_number_t>& aois = m.atom_order_in_smiles();
  const int matoms = m.natoms();
  //for (int i = 0; i < matoms; ++i) {
  //  cerr << i << ' ' << m.smarts_equivalent_for_atom(aois[i]) << " aois " << aois[i] << '\n';
  //}
  atom_order_in_smiles.reset(new int[matoms]);

  std::copy_n(aois.rawdata(), matoms, atom_order_in_smiles.get());

  return 1;
}

int
ConformerDB::WriteConformers(MoleculeWithEnergy& m, const conformer_database::Conformers& conformers) {
  // This will happen on the first call.
  if (conformers.conformers_size() == 0) [[ unlikely]] {
    return 1;
  }

  if (_report_progress()) {
    cerr << "Stored conformers for " << iwstring::with_commas(_molecules_read)
         << " molecules\n";
  }

  _number_conformers[conformers.conformers_size()]++;

  // Even though conformers may have explicit Hydrogens, the database key does not.
  m.remove_all(1);

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  const IWString& usmi = m.unique_smiles();

  // For now, duplicates are discarded. Too hard...
  Dbt dkey((void *) usmi.rawdata(), usmi.length());
  Dbt zdata;
  if (int rc = _db->get(NULL, &dkey, &zdata, 0); rc == 0) {
    ++_duplicates;
    if (_verbose > 1) {
      cerr << "ConformerDB::WriteConformers:existing data found\n";
      cerr << m.smiles() << ' ' << m.name() << '\n';
    }
    return 0;
  }

  std::string serialized;
  if (!conformers.SerializeToString(&serialized)) {
    cerr << "SerializeToString failed\n";
    return 0;
  }

  zdata.set_data(serialized.data());
  zdata.set_size(serialized.size());

  if (int rc = _db->put(NULL, &dkey, &zdata, 0); rc == 0) {
    return 1;
  } else {
    cerr << "Error storing " << m.name() << '\n';
    _db->err(rc, "");
    return 0;
  }

  return 1;
}

int
ConformerDB::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << "Skipped " << _duplicates << " duplicate data\n";
  for (int i = 0; i < _number_conformers; ++i) {
    if (_number_conformers[i]) {
      output << _number_conformers[i] << " molecules had " << i << " conformers\n";
    }
  }

  return 1;
}

// Look through the text_info in `m` and find the _energy_tag record.
// The next record should hold the energy.
// Set m.energy with the resulting value.
int
ConformerDB::GetEnergy(MoleculeWithEnergy& m) {
  // We should not have been called if the tag was empty.
  if (_energy_tag.empty()) [[unlikely]] {
    return 0;
  }

  const int n = m.number_records_text_info();
  if (n == 0) {
    cerr << "ConformerDB::GetEnergy:no text records in " << m.name() << '\n';
    return 0;
  }

  for (int i = 0; i < n; ++i) {
    const IWString& s = m.text_info(i);
    if (! s.starts_with(_energy_tag)) {
      continue;
    }

    if (i == (n - 1)) [[unlikely]] {
      cerr << "ConformerDB::GetEnergy:truncated text info " << m.name() << '\n';
      return 0;
    }

    const IWString& e = m.text_info(i + 1);
    float tmp;
    if (! e.numeric_value(tmp)) {
      cerr << "ConformerDB::GetEnergy:invalid energy value '" << e << "' for " << m.name() << '\n';
      return 0;
    }

    m.set_energy(tmp);
    return 1;
  }

  cerr << "ConformerDB::GetEnergy:did not find " << _energy_tag << " in " << m.name() << '\n';
  return 0;
}

// Return a set of conformers - share same name and number of atoms.
int
ConformerDB::NextSetOfConformers(data_source_and_type<MoleculeWithEnergy>& input,
                    resizable_array_p<MoleculeWithEnergy>& molecules) {
  molecules.resize_keep_storage(0);

  // We detect end of a group when we read a molecule with a different name.
  // We retain that molecule here and it becomes the first member of the next group.
  static MoleculeWithEnergy* next_molecule = nullptr;

  IWString mname;
  // all molecules in a set of conformers must have the same number of atoms.
  int natoms = 0;

  if (next_molecule != nullptr) [[ likely ]] {
    molecules << next_molecule;
    mname = next_molecule->name();
    natoms = next_molecule->natoms();
    next_molecule = nullptr;
  }

  MoleculeWithEnergy * m;
  while ((m = input.next_molecule()) != nullptr) {
    if (_energy_tag.size() > 0 && ! GetEnergy(*m)) {
      cerr << "ConformerDB::NextSetOfConformers:cannot retrieve energy\n";
      return 0;
    }

    if (mname.empty()) [[ unlikely ]] {  // first time called only.
      mname = m->name();
      natoms = m->natoms();
      molecules << m;
    } else if (mname == m->name()) {
      if (m->natoms() != natoms) {
        cerr << "NextSetOfConformers:atom count mismatch: expected " << natoms << 
                " got " << m->natoms() << ' ' << m->name() << '\n';
        return 0;
      }

      molecules << m;
    } else {
      next_molecule = m;
      return molecules.number_elements();
    }
  }

  // Normal EOF
  return molecules.number_elements();
}

int
ConformerDB::Process(data_source_and_type<MoleculeWithEnergy>& input) {
  resizable_array_p<MoleculeWithEnergy> conformers;

  while (NextSetOfConformers(input, conformers)) {
    _molecules_read += conformers.size();

    std::unique_ptr<int[]> atom_order_in_smiles;
    conformer_database::Conformers proto;
    Reset(*conformers[0], atom_order_in_smiles, proto);

    for (MoleculeWithEnergy* m : conformers) {
      AddConformer(*m, atom_order_in_smiles, proto);
    }
    WriteConformers(*conformers[0], proto);
  }

  return 1;
}

int
ConformerDatabaseStore(ConformerDB& cdb, const char* fname, FileType input_type) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<MoleculeWithEnergy> input(input_type, fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return cdb.Process(input);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:g:d:i:r:G:");
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

  for (const char* fname : cl) {
    if (! ConformerDatabaseStore(cdb, fname, input_type)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
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
