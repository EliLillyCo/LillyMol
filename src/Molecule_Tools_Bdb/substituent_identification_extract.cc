// List data from a database built by substituent_identification.

#include <string.h>

#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>

#include "db_cxx.h"

#include "google/protobuf/text_format.h"

#include "absl/strings/string_view.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iwstring.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools_Bdb/substituent_identification.pb.h"
#else
#include "substituent_identification.pb.h"
#endif

namespace substituent_identification {

using std::cerr;

struct DBKey {
  int radius;
  uint32_t b;
};

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << R"(Lists the contents of a BerkeleyDB built by substitutent_identification.
 -v             verbose output.
)";
  // clang-format off
  // clang-format on

  ::exit(rc);
}

class Options {
  private:
    int _verbose;

    Db _db;

    int _min_radius;
    int _max_radius;

    uint32_t _fail_radius;

    uint32_t _min_count;
    uint32_t _max_count;

    uint32_t _fail_count;

    uint64_t _records_read;
    uint64_t _records_written;

    google::protobuf::TextFormat::Printer _printer;

    Accumulator_Int<uint64_t> _bytes;
    Accumulator_Int<uint64_t> _number_replacements;

    char _output_separator;

  // Private functions
    int List(const Dbt& zkey, const Dbt& zdata, IWString_and_File_Descriptor& output);
    int List(const DBKey& key, const substituent_identification::Replacements& proto,
              IWString_and_File_Descriptor& output);

  public:
    Options();
    ~Options();

    int Initialise(Command_Line& cl);

    int List(IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() : _db(NULL, DB_CXX_NO_EXCEPTIONS) {
  _verbose = 0;

  _min_radius = 0;
  _max_radius = std::numeric_limits<int>::max();

  _fail_radius = 0;

  _min_count = 0;
  _max_count = std::numeric_limits<uint32_t>::max();

  _records_read = 0;
  _records_written = 0;

  _fail_count = 0;

  _output_separator = ' ';

  _printer.SetSingleLineMode(true);
}

Options::~Options() {
  // What happens if the db has not been opened??
  _db.close(0);
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Must specify name of database to list as a command line option\n";
    Usage(1);
  }

  const char* dbname = cl[0];

  if (int rc = _db.open(NULL, dbname, NULL, DB_UNKNOWN, DB_RDONLY, 0); rc != 0) {
    cerr << "Cannot open database '" << dbname << "' '";
    _db.err(rc, "");
    return 0;
  }

  if (cl.option_present('r')) {
    if (! cl.value('r', _min_radius) || _min_radius < 0) {
      cerr << "Invalid min radius (-r)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only list fragments with radius > " << _min_radius << " bonds\n";
    }
  }

  if (cl.option_present('R')) {
    if (! cl.value('R', _max_radius) || _max_radius < _min_radius) {
      cerr << "Invalid max radius (-R)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only list fragments with radius < " << _max_radius << " bonds\n";
    }
  }

  if (cl.option_present('o')) {
    IWString o;
    cl.value('o', o);

    static constexpr int kNoMessage = 0;

    if (! char_name_to_char(o, kNoMessage)) {
      cerr << "Invalid character name '" << o << "'\n";
      return 0;
    }
    _output_separator = o[0];
  }

  return 1;
}

int
Options::List(IWString_and_File_Descriptor& output) {
  Dbc* cursor = NULL;

  if (int rc = _db.cursor(NULL, &cursor, 0); rc != 0) {
    _db.err(rc, "cannot acquire cursor");
    return 0;
  }

  Dbt zkey, zdata;

  int rc;
  while (0 == (rc = cursor->get(&zkey, &zdata, DB_NEXT))) {
    _records_read++;
    if (! List(zkey, zdata, output)) {
      cerr << "Error\n";
      return 0;
    }
  }

  cursor->close();

  return 1;
}

// Return true if `zkey` is '_RADIUS'
bool
IsRadius(const Dbt& zkey) {
  if (zkey.get_size() != 7) {
    return false;
  }

  return ::strncmp((const char*)(zkey.get_data()), "_RADIUS", 7) == 0;
}

int
Options::List(const Dbt& zkey, const Dbt& zdata, IWString_and_File_Descriptor& output) {
  if (IsRadius(zkey)) {
    return 1;
  }

  if (zkey.get_size() != 8) {
    return 1;
  }

//std::string_view data((const char*) zdata.get_data(), zdata.get_size());
  const absl::string_view data(reinterpret_cast<const char*>(zdata.get_data()), zdata.get_size());

  const DBKey* key = reinterpret_cast<const DBKey*>(zkey.get_data());

  substituent_identification::Replacements replacements;
  if (! replacements.ParseFromString(data)) {
    cerr << "Options::List:cannot parse proto, length " << data.size() << '\n';
    cerr.write((const char*)zkey.get_data(), zkey.get_size());
//  return 0;
  }

  if (_verbose) {
    _bytes.extra(zdata.get_size());
    _number_replacements.extra(replacements.replacement_size());
  }

  return List(*key, replacements, output);
}

int
Options::List(const DBKey& key, const substituent_identification::Replacements& proto,
              IWString_and_File_Descriptor& output) {
  if (key.radius < _min_radius || key.radius > _max_radius) {
    ++_fail_radius;
//  cerr << "Fail radius\n";
    return 1;
  }

  if (uint32_t c = proto.replacement_size(); c < _min_count || c > _max_count) {
    ++_fail_count;
//  cerr << "Fail count\n";
    return 1;
  }

  ++_records_written;

  std::string buffer;
  if (! _printer.PrintToString(proto, &buffer)) {
    cerr << "Options::List cannot convert to text'" << proto.ShortDebugString() << "'\n";
    return 0;
  }

  output << key.radius << _output_separator << key.b << ' ' << buffer << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _records_read << " records, wrote " << _records_written << '\n';
  if (_min_radius > 0) {
    output << _fail_radius << " records failed radius " << _min_radius << " to " << _max_radius << '\n';
  }
  if (_min_count > 0) {
    output << _fail_count << " records with fewer than " << _min_count << " items\n";
  }
  output << "Data had btw " << _bytes.minval() << " and " << _bytes.maxval() <<
          " bytes, ave " << static_cast<float>(_bytes.average()) << '\n';

  output << "Replacements had btw " << _number_replacements.minval() << " and " <<
            _number_replacements.maxval() << " items ave " <<
            static_cast<float>(_number_replacements.average()) << '\n';
  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vr:R:o:c:C:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_present('v');

  if (cl.empty()) {
    cerr << "Must specify database to list as a command line argument\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);
  if (! options.List(output)) {
    cerr << "Error processing '" << cl[0] << "'\n";
    return 1;
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace substituent_identification

int
main(int argc, char** argv) {

  return substituent_identification::Main(argc, argv);
}
