// Collate multiple independently generated ring replacement sets.

#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

#include "absl/container/flat_hash_map.h"

#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/text_format.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"

#include "Molecule_Tools/replacement_ring.pb.h"

namespace ring_replacement_collate {

using std::cerr;
namespace fs = std::filesystem;

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
  cerr << R"(
)";
  // clang-format on
  ::exit(rc);
}

// We keep track of rings by type. And within each type there are protos
// describing the rings found.

using RingData = absl::flat_hash_map<std::string, absl::flat_hash_map<std::string, RplRing::ReplacementRing>>;

class Options {
  private:
    // Regular expression describing the file names to be processed.
    std::unique_ptr<RE2> _rx;

    int _ignore_errors;
    int _errors_encountered;

    std::string _prefix;

  // private functions
    int AggregateRings(const char* dirname,
                        const std::filesystem::directory_entry& fname,
                        RingData& rings);
    int AggregateRings(iwstring_data_source& input,
                        absl::flat_hash_map<std::string, RplRing::ReplacementRing>& rings);
    int AggregateRings(const RplRing::ReplacementRing& ring,
                        absl::flat_hash_map<std::string, RplRing::ReplacementRing>& rings);
    int WriteRings(const absl::flat_hash_map<std::string, RplRing::ReplacementRing>& rings,
                    const std::string name_stem,
                    const std::string ring_type);
    int WriteRings(const absl::flat_hash_map<std::string, RplRing::ReplacementRing>& rings,
                    IWString_and_File_Descriptor& output);

  public:
    Options();

    int Initialise(Command_Line_v2& cl);

    int AggregateRings(const char* dirname, 
                RingData& rings);

    // The name stem for where the output files are written.
    // If this includes a directory path, that must have already been
    // created.
    // The file name will be   prefix/_prefix_ringtype.smi
    int WriteRings(const RingData& rings, const std::string& name_stem);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _ignore_errors = 0;
  _errors_encountered = 0;

  const_IWSubstring tmp(R"(_([0-9,a,A]+)\.smi$)");
  iwre2::RE2Reset(_rx, tmp);

  _prefix = "ring_";
}

int
Options::Initialise(Command_Line_v2& cl) {
  int verbose = cl.option_present('v');

  if (cl.option_present("ignore_errors")) {
    _ignore_errors = 1;
    if (verbose) {
      cerr << "Will ignore errors\n";
    }
  }

  if (cl.option_present("prefix")) {
    cl.value("prefix", _prefix);
    if (verbose) {
      cerr << "Prefix for output files '" << _prefix << "'\n";
    }
  }

  return 1;
}

int
Options::AggregateRings(const char* dirname,
               RingData& rings) {
  //    for (auto it{fs::directory_iterator("sandbox")}; it != fs::directory_iterator(); ++it)

  for (auto& fname : fs::directory_iterator(dirname)) {
    if (! AggregateRings(dirname, fname, rings)) {
      cerr << "Error processing '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::AggregateRings(const char* dirname,
                        const std::filesystem::directory_entry& fname,
                        RingData& rings) {
  std::string ring_type;
  const std::string x = fname.path();
  if (! RE2::PartialMatch(x, *_rx, &ring_type)) {
    cerr << "Options::AggregateRings:cannot determine ring type '" << fname << 
            "' pattern '" << _rx->pattern() << "'\n";
    ++_errors_encountered;
    return _ignore_errors;
  }

  fs::path full_path_name(dirname);
  full_path_name /= fname;

  iwstring_data_source input(full_path_name.c_str());
  if (! input.good()) {
    cerr << "Options::AggregateRings:cannot open '" << full_path_name << "'\n";
    return 0;
  }

  auto iter = rings.find(ring_type);
  if (iter != rings.end()) {
     return AggregateRings(input, iter->second);
  }

  // New ring type, create;
  absl::flat_hash_map<std::string, RplRing::ReplacementRing> proto;
  auto iter2 = rings.emplace(ring_type, std::move(proto));

  return AggregateRings(input, std::get<1>(*iter2.first));
}

int
Options::AggregateRings(iwstring_data_source& input,
                        absl::flat_hash_map<std::string, RplRing::ReplacementRing>& rings) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    RplRing::ReplacementRing proto;
    google::protobuf::io::ArrayInputStream zero_copy_array(buffer.data(), buffer.nchars());
    if (!google::protobuf::TextFormat::Parse(&zero_copy_array, &proto)) {
      cerr << "Options::AggregateRings:cannot parse " << buffer << "'\n";
      return 0;
    }

    if (! AggregateRings(proto, rings)) {
      cerr << "Options::AggregateRings:cannot parse proto " << buffer << '\n';
      return 0;
    }
  }

  // cerr << "At end of " << input.fname() << " hrve " << rings.size() << " lines_read " << input.lines_read() << '\n';

  return 1;
}

int
Options::AggregateRings(const RplRing::ReplacementRing& ring,
                        absl::flat_hash_map<std::string, RplRing::ReplacementRing>& rings) {
  if (ring.usmi().size() == 0) {
    cerr << "EMpty usmi " << ring.ShortDebugString() << '\n';
    return 0;
  }

  auto iter = rings.find(ring.usmi());
  if (iter != rings.end()) {
    const auto n = ring.n();
    iter->second.set_n(n + ring.n());
    return 1;
  }

  rings.emplace(ring.usmi(), ring);

  return 1;
}

int
Options::Report(std::ostream& output) const {
  return 1;
}

int
Options::WriteRings(const RingData& rings, const std::string& name_stem) {
  for (const auto& [ring_type, myrings] : rings) {
    if (! WriteRings(myrings, name_stem, ring_type)) {
      cerr << "Options::WriteRings:cannot write " << name_stem << ' ' << ring_type << '\n';
      return 0;
    }
  }

  return 1;
}

int
Options::WriteRings(const absl::flat_hash_map<std::string, RplRing::ReplacementRing>& rings,
                    const std::string destdir,
                    const std::string ring_type) {
  IWString full_path_name;
  full_path_name << destdir << '/' << _prefix << ring_type << ".smi";

  IWString_and_File_Descriptor output;
  if (! output.open(full_path_name)) {
    cerr << "Options::WriteRings:cannot open '" << full_path_name << "'\n";
    return 0;
  }

  return WriteRings(rings, output);
}

int
Options::WriteRings(const absl::flat_hash_map<std::string, RplRing::ReplacementRing>& rings,
                    IWString_and_File_Descriptor& output) {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;

  for (const auto& [_, proto] : rings) {
    if (! printer.PrintToString(proto, &buffer)) {
      cerr << "Options::WriteRings write '" << proto.ShortDebugString() << "'\n";
      return 0;
    }

    output << buffer;
    output << '\n';

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}
                        

int
Main(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v--destdir=s-ignore_errors-prefix=s");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }
  
  const int verbose = cl.option_present('v');

  Options options;
  if ( !options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (! cl.option_present("destdir")) {
    cerr << "Must specify destination directory via the -destdir option\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  RingData rings;

  for (const char* dir : cl) {
    if (! options.AggregateRings(dir, rings)) {
      cerr << "Error processing directory '" << dir << "'\n";
      return 1;
    }
  }

  if (rings.empty()) {
    cerr << "NO rings found\n";
    return 1;
  }

  if (verbose) {
    cerr << "Read data on " << rings.size() << " ring types\n";
  }

  IWString tmp;
  cl.value("destdir", tmp);
  std::string destdir(tmp.AsString());
  fs::path path_name(destdir);
  if (! fs::is_directory(path_name)) {
    if (! fs::create_directories(path_name)) {
      cerr << "Cannot create '" << path_name << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << "Writing " << rings.size() << " ring types\n";
    Accumulator_Int<int32_t> acc;
    for (const auto& [_, ring_type] : rings) {
      for (const auto& [_, proto] : ring_type) {
        acc.extra(proto.n());
      }
    }
    cerr << "Counts btw " << acc.minval() << " and " << acc.maxval() <<
            " ave " << acc.average() << '\n';
  }

  if (! options.WriteRings(rings, destdir)) {
    cerr << "Cannot write rings to " << destdir << "'\n";
    return 1;
  }

  return 0;
}

}  // namespace ring_replacement_collate

int
main(int argc, char** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  return ring_replacement_collate::Main(argc, argv);
}
