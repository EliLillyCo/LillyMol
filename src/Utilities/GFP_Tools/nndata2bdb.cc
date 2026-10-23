// Convert a file of nnbr::NearNeighbours protos so a BerkeleyDB database

#include "sys/stat.h"
#include "sys/types.h"

#include <cstdint>
#include <iostream>
#include <optional>
#include <string_view>

#include "db_cxx.h"

#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/data_source/tfdatarecord.h"

#ifdef BUILD_BAZEL
#include "Utilities/GFP_Tools/nearneighbours.pb.h"
#else
#include "nearneighbours.pb.h"
#endif

namespace nndata2bdb {

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
  cerr << R"(Creates a BerkeleyDB database of nnbr::NearNeighbours protos from a TFDataRecord
file of serialized protos.

gfp_nearneighbours_single_file_tbb -n 16 -S haystack.dat -n 10 haystack.gfp
nndata2bdb -d haystack.bdb haystack.dat

 -d <dbname>            name of database to create.
 -rpt <n>               report progress every <n> items stored.
 -v                     verbose output.
)";
  // clang-format on

  ::exit(rc);
}

class Worker {
  private:
    int _verbose;

    Db* _db;

    uint64_t _report_progress;
    uint64_t _next_report;

    uint64_t _items_read;

  public:
    Worker();
    ~Worker();

    int Initialise(Command_Line_v2& cl);

    int Store(const nnbr::NearNeighbours& proto, const std::string_view& serialized);

    int Report(std::ostream& output) const;
};

Worker::Worker() {
  _verbose = 0;
  _db = nullptr;
  _report_progress = 0;
  _next_report = 0;
  _items_read = 0;
}

Worker::~Worker() {
  if (_db != nullptr) {
    _db->close(0);
    delete _db;
  }
}

int
Worker::Initialise(Command_Line_v2& cl) {
  _verbose = cl.option_present('v');

  if (! cl.option_present('d')) {
    cerr << "Worker::Initialise:must specify name of database via the -d option\n";
    Usage(1);
  }

  if (cl.option_present('d')) {
    IWString dbname = cl.string_value('d');
    dbname.EnsureEndsWith(".bdb");

    int flags;
    DBTYPE dbtype;
    int mode;

    if (dash_s(dbname.null_terminated_chars())) {
      dbtype = DB_UNKNOWN;
      flags = 0;
      mode = 0;
    } else {
      dbtype = DB_BTREE;
      flags = DB_CREATE;
      mode = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;
    }

    _db = new Db(NULL, DB_CXX_NO_EXCEPTIONS);

    if (int rc = _db->open(NULL, dbname, NULL, dbtype, flags, mode); rc != 0) {
      cerr << "Cannot open database '" << dbname << "'\n";
      _db->err(rc, "");
      return 2;
    }

    if (_verbose) {
      cerr << "Opened '" << dbname << "' for writing\n";
    }
  }

  if (cl.option_present("rpt")) {
    if (! cl.value("rpt", _report_progress) || _report_progress == 0) {
      cerr << "The report progress option -rpt must be a whole +ve integer\n";
      Usage(1);
    }

    if (_verbose) {
      cerr << "Will report progress every " << _report_progress << " values processed\n";
    }
  }

  return 1;
}

// If `name` contains multiple words, truncate to the first word.
// Kind of fails if `name` starts with a space...
void
MaybeTruncateName(IWString& name) {
  name.remove_leading_chars(' ');
  name.truncate_at_first(' ');
}

int
Worker::Store(const nnbr::NearNeighbours& proto,
              const std::string_view& serialized) {
  IWString id = proto.name();
  if (id.empty()) {
    cerr << "NNData2Bdb:empty id\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  MaybeTruncateName(id);

  ++_items_read;

  Dbt dbkey((void*)id.rawdata(), id.size());

  Dbt zdata((void*)serialized.data(), serialized.length());

  // By default, we replace any existing content.
  if (int rc = _db->put(NULL, &dbkey, &zdata, 0); rc != 0) {
    cerr << "NNData2Bdb:cannot store '" << id << " data size " << serialized.length() << '\n';
    _db->err(rc, "");
    return 0;
  }

  return 1;
}

int
Worker::Report(std::ostream& output) const {
  output << "Read " << _items_read << " items\n";

  return 1;
}

int
NNData2Bdb(const nnbr::NearNeighbours& proto,
           const std::string_view& serialized, Worker& worker) {
  return worker.Store(proto, serialized);
}

int
NNData2Bdb(iw_tf_data_record::TFDataReader& reader, Worker& worker) {
  while (1) {
    std::optional<const_IWSubstring> serialized = reader.Next(); 
    if (! serialized) {
      if (reader.eof()) {
        return 1;
      }

      cerr << "NNData2Bdb:error return\n";
      return 0;
    }

    const std::string_view as_string(serialized->data(), serialized->length());

    nnbr::NearNeighbours proto;
    if (! proto.ParseFromString(as_string)) {
      cerr << "NNData2Bdb:invalid serialized proto\n";
      return 0;
    }

    if (! NNData2Bdb(proto, as_string, worker)) {
      return 0;
    }
  }

  return 1;
}

int
NNData2Bdb(const char* fname, Worker& worker) {
  iw_tf_data_record::TFDataReader reader(fname);
  if (! reader.good()) {
    cerr << "NNData2Bdb:cannot open '" << fname << "'\n";
    return 0;
  }

  return NNData2Bdb(reader, worker);
}

int
Main(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-d=s-rpt=ipos");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  Worker worker;
  if (! worker.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  for (const char* fname : cl) {
    if (! NNData2Bdb(fname, worker)) {
      cerr << "Cannot process '" << fname << "'\n";
      return 1;
    }
  }

  if (cl.option_present('v')) {
    worker.Report(cerr);
  }

  return 0;
}

}  // namespace nndata2bdb

int
main(int argc, char **argv) {
  int rc = nndata2bdb::Main(argc, argv);

  return rc;
}
