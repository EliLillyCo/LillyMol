// Read serialized protos from dicer and load a BerkeleyDB database'
// consisting of
//  key: fragment unique smiles
//  value: dicer textproto
// The key must have an isotopic label that matches the isotope in the value

#include "sys/stat.h"
#include "sys/types.h"

#include <iostream>
#include <memory>
#include <string>

#include "absl/container/flat_hash_map.h"

#include "db_cxx.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"

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

  _complement[proto.comp()] = std::move(frag);

  return 1;
}

int
Complements::WriteAsText(const std::string& lhs, IWString_and_File_Descriptor& output) const {
  for (const auto& [k, v] : _complement) {
    output << lhs << ' ';
    output << v.ShortDebugString() << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

class Options {
  private:
    int _verbose;

    std::unique_ptr<Db> _database;

    uint64_t _molecules_read;

    // When storing complementary fragments, we control how different
    // the atom count can be. For example if we examine C as a fragment,
    // there will be a huge number of complementary fragments available
    // for it.
    int _max_extra_atoms;
    int _max_fewer_atoms;

//  absl::flat_hash_map<std::string, dicer_data::DicerFragment> _fragment;
    absl::flat_hash_map<std::string, Complements> _fragment;

  // Private functions.
    int Process(const std::string& parent_name, const dicer_data::DicerFragment& frag);
    int OkAtomCountDifference(int n1, int n2) const;

  public:
    Options();
    ~Options();

    int Initialise(Command_Line& cl);

    int Process(const dicer_data::DicedMolecule& proto);

    int WriteAsText(IWString_and_File_Descriptor& output) const;
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

  if (cl.option_present('c')) {
    if (! cl.value('c', _max_fewer_atoms) || _max_fewer_atoms < 0) {
      cerr << "The maximum fewer atoms in the complement (-c) option must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will not store a complementary fragment if it has " << _max_fewer_atoms << " fewer atoms\n";
    }
  }

  if (cl.option_present('C')) {
    if (! cl.value('C', _max_extra_atoms) || _max_extra_atoms < 0) {
      cerr << "The maximum examine atoms in the complement (-C) option must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will not store a complementary fragment if it has " << _max_extra_atoms << " extra atoms\n";
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

  for (const dicer_data::DicerFragment& frag : proto.fragment()) {
    Process(proto.name(), frag);
  }

  return 1;
}

int
Options::Process(const std::string& parent_name, const dicer_data::DicerFragment& frag) {
#ifdef THIS_NEVER_HAPPENS
  if (! frag.has_comp()) {
    cerr << "Options::Process:no complementary fragment\n";
    return 0;
  }
#endif

  auto iter = _fragment.find(frag.smi());
  if (iter != _fragment.end()) {
    return iter->second.Extra(parent_name, frag);
  }

  Complements comp;
  comp.Initialise(parent_name, frag);

  _fragment[frag.smi()] = std::move(comp);

  return 1;
}

// We are processing a fragment with `n1` atoms and are thinking about
// adding a complementary fragment with `n2` atoms.
// Return 1 if this is consistent with the atom count difference constraints.
int
Options::OkAtomCountDifference(int n1, int n2) const {
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

  return 1;
}

int
Options::WriteAsText(IWString_and_File_Descriptor& output) const {
  for (const auto& [k, v] : _fragment) {
    v.WriteAsText(k, output);
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
DicerToComplmentDb(iw_tf_data_record::TFDataReader& input,
                   Options& options) {
  while (1) {
    std::optional<dicer_data::DicedMolecule> maybe_proto = input.ReadProto<dicer_data::DicedMolecule>();
    if (! maybe_proto) {
      return 1;
    }

    // cerr << maybe_proto->ShortDebugString() << '\n';
    options.Process(*maybe_proto);
  }

  return 1;
}

int
DicerToComplmentDb(const char* fname, Options& options) {
  iw_tf_data_record::TFDataReader input(fname);

  if (! input.good()) {
    cerr << "DicerToComplmentDb:cannot open '" << fname << "'\n";
    return 0;
  }

  return DicerToComplmentDb(input, options);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vd:c:C:");

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

  for (const char* fname: cl) {
    if (! DicerToComplmentDb(fname, options)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  IWString_and_File_Descriptor output(1);
  options.WriteAsText(output);

  return 0;
}

}  // namespace dicer_to_complement_db

int
main(int argc, char **argv) {
  int rc = dicer_to_complement_db::Main(argc, argv);

  return rc;
}
