// Read serialized protos from dicer and load a BerkeleyDB database'
// consisting of
//  key: fragment unique smiles
//  value: dicer textproto
// The key must have an isotopic label that matches the isotope in the value

#include "sys/stat.h"
#include "sys/types.h"

#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include "google/protobuf/text_format.h"

#include "absl/container/flat_hash_map.h"

#include "db_cxx.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/string_with_commas.h"

#include "Molecule_Lib/molecule.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/dicer_fragments.pb.h"
#else
#include "dicer_fragments.pb.h"
#endif

namespace dicer_to_complement_db {

using std::cerr;

int
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(Reads dicer output and constructs a BerkeleyDB database for the complementary fragment tool.
A typical usage might be 

dicer -B serialized_proto -S /tmp/rand.dc -C iso=9 -I 9 -m 4 -M 17 -k 2 -X 500 file.smi
dicer_to_complement_db -d /tmp/rand.bdb /tmp/rand.dc

 -d <dbname>            The database to generate.
 -f <natoms>            Minimum fragment size.
 -F <natoms>            maximum fragment size.
 -c <natoms>            Minimum complement size.
 -C <natoms>            maximum complement size.
 -x <natoms>            When comparing the fragment with the complement, the complement can lose at most <natoms>
 -X <natoms>            When comparing the fragment with the complement, the complement can gain at most <natoms>
 -w                     Write the stored data before writing to db. Beware may be very large.
 -r <n>                 Report progress every <n> molecules processed.
 -v                     Verbose output.
                
)";
  // clang-format on
  ::exit(rc);
}

// As fragments are gathered, we collect a set of complementary fragments
// associated with that fragment. 
class Complements {
  private:
    // Even though we are storing the complementary fragments, this is the
    // number of atoms in the LHS. Just because there is no other place
    // to easily and efficiently store it.
    int _natoms;

    absl::flat_hash_map<std::string, dicer_data::DicerFragment> _complement;

  public:
    Complements();

    int Initialise(const std::string& parent_name, const dicer_data::DicerFragment& proto);

    int natoms() const {
      return _natoms;
    }

    int Extra(const std::string& parent_name, const dicer_data::DicerFragment& proto);

    int ToProto(dicer_data::ComplementaryFragments& proto) const;

    int WriteAsText(const std::string& lhs, IWString_and_File_Descriptor& output) const;
};

Complements::Complements() {
  _natoms = 0;
}

int
Complements::Initialise(const std::string& parent_name, const dicer_data::DicerFragment& proto) {
  _natoms = proto.nat();

  // Make a copy since we need to set the parent attribute.
  dicer_data::DicerFragment tmp(proto);
  tmp.set_par(parent_name);
  tmp.set_n(1);
  _complement[proto.comp()] = std::move(tmp);

  return 1;
}

int
Complements::Extra(const std::string& parent_name, const dicer_data::DicerFragment& proto) {
  auto iter = _complement.find(proto.comp());
  if (iter != _complement.end()) {
    uint32_t n = iter->second.n();
    iter->second.set_n(n + 1);
    return 1;
  }

  dicer_data::DicerFragment frag;
  frag.set_smi(proto.comp());
  frag.set_n(1);
  frag.set_par(parent_name);
  cerr << "New fragment " << frag.ShortDebugString() << '\n';

  _complement[proto.comp()] = std::move(frag);

  return 1;
}

int
Complements::WriteAsText(const std::string& lhs, IWString_and_File_Descriptor& output) const {
  google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(1);

  std::string buffer;

  for (const auto& [k, v] : _complement) {
    output << lhs << ' ';
    if (! printer.PrintToString(v, &buffer)) {
      cerr << "PrintToString failed\n";
      return 0;
    }

    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
Complements::ToProto(dicer_data::ComplementaryFragments& proto) const {
 
  for (const auto& [smi, v] : _complement) {
    dicer_data::ComplementaryFragment* f = proto.mutable_frag()->Add();
    f->set_smi(smi);
    f->set_id(v.par());
    f->set_n(v.n());
  }

  return 1;
}

class Options {
  private:
    int _verbose;

    std::unique_ptr<Db> _database;

    uint64_t _molecules_read = 0;

    // The size of fragments can be controlled by dicer, but that is expensive.
    // We can cheaply filter things here.
    uint64_t _atom_count_fail = 0;

    // In order to keep database sizes down, impose size limits on what gets stored.
    int _min_fragment_size = 0;
    int _max_fragment_size = std::numeric_limits<int>::max();
    int _min_complement_size = 0;
    int _max_complement_size = std::numeric_limits<int>::max();

    // When storing complementary fragments, we control how different
    // the atom count can be. For example if we examine C as a fragment,
    // there will be a huge number of complementary fragments available
    // for it.
    int _max_extra_atoms;
    int _max_fewer_atoms;

    Report_Progress _report_progress;

//  absl::flat_hash_map<std::string, dicer_data::DicerFragment> _fragment;
    absl::flat_hash_map<std::string, Complements> _fragment;

    // Useful for debugging and development.
    uint32_t _stop_after_processing = 0;

  // Private functions.
    int Process(const dicer_data::DicedMolecule& proto, const dicer_data::DicerFragment& frag);

    int OkAtomCounts(int n1, int n2) const;

    int StoreInDb(const IWString& smiles, const Complements& complements);

  public:
    Options();
    ~Options();

    int Initialise(Command_Line& cl);

    int Process(const dicer_data::DicedMolecule& proto);

    int WriteAsText(IWString_and_File_Descriptor& output) const;

    // Write the contents of _fragment t- _database.
    int StoreInDb();
};

Options::Options() {
  _verbose = 0;

  _molecules_read = 0;

  _max_extra_atoms = 0;
  _max_fewer_atoms = 0;
}

Options::~Options() {
  if (_database) {
    _database->close(0);
  }
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (! cl.option_present('d')) {
    cerr << "Must specify name of database to build via the -d option\n";
    return 0;
  }

  if (cl.option_present('f')) {
    if (! cl.value('f', _min_fragment_size) || _min_fragment_size < 1) {
      cerr << "The minimum fragment atom count (-c) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Min fragment size " << _min_fragment_size << '\n';
    }
  }

  if (cl.option_present('F')) {
    if (! cl.value('F', _max_fragment_size) || _max_fragment_size < _min_fragment_size) {
      cerr << "The maximum fragment atom count (-C) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Max fragment size " << _max_fragment_size << '\n';
    }
  }

  if (cl.option_present('c')) {
    if (! cl.value('c', _min_complement_size) || _min_complement_size < 1) {
      cerr << "The minimum complement atom count (-c) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Min complement size " << _min_complement_size << '\n';
    }
  }

  if (cl.option_present('C')) {
    if (! cl.value('C', _max_complement_size) || _max_complement_size < _min_complement_size) {
      cerr << "The maximum complement atom count (-C) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Max complement size " << _max_complement_size << '\n';
    }
  }

  if (cl.option_present('x')) {
    if (! cl.value('x', _max_fewer_atoms) || _max_fewer_atoms < 0) {
      cerr << "The maximum fewer atoms in the complement (-c) option must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will not store a complementary fragment if it has " << _max_fewer_atoms << " fewer atoms\n";
    }
  }

  if (cl.option_present('X')) {
    if (! cl.value('X', _max_extra_atoms) || _max_extra_atoms < 0) {
      cerr << "The maximum examine atoms in the complement (-C) option must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will not store a complementary fragment if it has " << _max_extra_atoms << " extra atoms\n";
    }
  }

  if (cl.option_present('r')) {
    if (! _report_progress.initialise(cl, 'r', _verbose)) {
      cerr << "Cannot initialise progress reporting (-r)\n";
      return 0;
    }
  }

  if (cl.option_present('d')) {
    _database = std::make_unique<Db>(nullptr, DB_CXX_NO_EXCEPTIONS);

    const char *dbname = cl.option_value('d');

    int flags;
    DBTYPE dbtype;
    int mode;

    if (dash_s(dbname)) {
      dbtype = DB_UNKNOWN;
      flags = 0;
      mode = 0;
    } else {
      dbtype = DB_BTREE;
      flags = DB_CREATE;
      mode = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;
    }

    int rc = _database->open(NULL, dbname, NULL, dbtype, flags, mode);

    if (0 != rc) {
      cerr << "Cannot open database '" << dbname << "'\n";
      _database.get()->err(rc, "");
      return 2;
    }

    if (_verbose) {
      cerr << "Smiles will be written to database '" << dbname << "'\n";
    }
  }

  return 1;
}

int
Options::Process(const dicer_data::DicedMolecule& proto) {

  ++_molecules_read;

  for (const dicer_data::DicerFragment& frag : proto.fragment()) {
    Process(proto, frag);
  }

  return 1;
}

int
Options::Process(const dicer_data::DicedMolecule& proto, const dicer_data::DicerFragment& frag) {
#ifdef THIS_NEVER_HAPPENS
  if (! frag.has_comp()) {
    cerr << "Options::Process:no complementary fragment\n";
    return 0;
  }
#endif

  if (_report_progress()) {
    cerr << "Read " << iwstring::with_commas(_molecules_read)
         << " molecules, hash contains " << iwstring::with_commas(_fragment.size())
         << " complements\n";

  }

  const int atoms_in_parent = proto.natoms();
  const int atoms_in_frag = frag.nat();
  if (! OkAtomCounts(atoms_in_parent, atoms_in_frag)) {
    ++_atom_count_fail;
    return 1;
  }

  auto iter = _fragment.find(frag.smi());
  if (iter != _fragment.end()) {
    iter->second.Extra(proto.name(), frag);
  } else {
    Complements comp;
    comp.Initialise(proto.name(), frag);
    _fragment[frag.smi()] = std::move(comp);
  }

  return 1;
}

// We have a fragmentation with `atoms_in_parent` atoms in the parent
// and `atoms_in_fragment` atoms in the fragment.
// Return true if this is an OK scenario
int
Options::OkAtomCounts(int atoms_in_parent, int atoms_in_fragment) const {
  assert (atoms_in_parent > atoms_in_fragment);

  int atoms_in_complement = atoms_in_parent - atoms_in_fragment;

  if (atoms_in_fragment < _min_fragment_size) {
    return 0;
  }
  if (atoms_in_fragment > _max_fragment_size) {
    return 0;
  }

  if (atoms_in_complement < _min_complement_size) {
    return 0;
  }
  if (atoms_in_complement > _max_complement_size) {
    return 0;
  }

#ifdef NOT_IMPLEMENTED
  if (_max_extra_atoms > 0) {
    if (n1 < n2) {
      if ((n2 - n1) > _max_extra_atoms) {
        return 0;
      }
    }
  }

  if (_max_fewer_atoms > 0) {
    if (n1 > n2) {
      if ((n1 - n2) > _max_fewer_atoms) {
        return 0;
      }
    }
  }
#endif

  return 1;
}

int
Options::WriteAsText(IWString_and_File_Descriptor& output) const {
  cerr << "Options::WriteAsText called\n";
  for (const auto& [k, v] : _fragment) {
    v.WriteAsText(k, output);
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
Options::StoreInDb() {
  if (_report_progress.active()) {
    _report_progress.reset();
  }

  for (const auto& [k, v] : _fragment) {
    if (! StoreInDb(k, v)) {
      cerr << "Error processing '" << k << "'\n";
      return 0;
    }

    if (_report_progress()) {
      cerr << "Stored " << iwstring::with_commas(_report_progress.times_called())
           << " items\n";
    }
  }

  if (_verbose) {
    cerr << "Stored " << iwstring::with_commas(_fragment.size()) << " fragments\n";
  }

  return 1;
}

int
Options::StoreInDb(const IWString& smiles, const Complements& complements) {
  Dbt zkey((char*)smiles.data(), smiles.length());

  dicer_data::ComplementaryFragments proto;
  complements.ToProto(proto);

  std::string buffer;
  if (! proto.SerializeToString(&buffer)) {
    cerr << "SerializeToString failed\n";
    return 0;
  }

  Dbt zdata((char*)buffer.data(), buffer.size());

  if (int rc = _database->put(NULL, &zkey, &zdata, 0); rc != 0) {
    cerr << "Options::StoreInDb:store failed ";
    _database->err(rc, "");
    return 0;
  }

  return 1;
}

int
DicerToComplmentDb(iw_tf_data_record::TFDataReader& input,
                   uint32_t stop_after_processing,
                   Options& options) {
  uint32_t items_read = 0;
  while (1) {
    std::optional<dicer_data::DicedMolecule> maybe_proto = input.ReadProto<dicer_data::DicedMolecule>();
    if (! maybe_proto) {
      return 1;
    }

    // cerr << maybe_proto->ShortDebugString() << '\n';
    options.Process(*maybe_proto);
    ++items_read;
    if (items_read > stop_after_processing) {
      cerr << "Processing stopped after " << items_read << " items\n";
      return 1;
    }
  }

  return 1;
}

int
DicerToComplmentDb(const char* fname,
                   uint32_t stop_after_processing,
                   Options& options) {
  iw_tf_data_record::TFDataReader input(fname);

  if (! input.good()) {
    cerr << "DicerToComplmentDb:cannot open '" << fname << "'\n";
    return 0;
  }

  return DicerToComplmentDb(input, stop_after_processing, options);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vd:f:F:c:Cx:X::r:ws:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments, must specify file generated by -S option from dicer\n";
    Usage(1);
  }

  uint32_t stop_after_processing = std::numeric_limits<uint32_t>::max();

  if (cl.option_present('s')) {
    if (! cl.value('s', stop_after_processing)) {
      cerr << "The stop processing option (-s) must be a whole +ve number\n";
      return 0;
    }
  }

  for (const char* fname: cl) {
    if (! DicerToComplmentDb(fname, stop_after_processing, options)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  IWString_and_File_Descriptor output(1);
  if (cl.option_present('w')) {
    options.WriteAsText(output);
  }

  if (! options.StoreInDb()) {
    cerr << "DB store failed\n";
    return 1;
  }

  return 0;
}

}  // namespace dicer_to_complement_db

int
main(int argc, char **argv) {
  int rc = dicer_to_complement_db::Main(argc, argv);

  return rc;
}
