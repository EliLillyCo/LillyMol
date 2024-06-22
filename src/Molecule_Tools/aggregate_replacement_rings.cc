// Aggregate a bunch of RplRing::ReplacementRing protos

#include <iostream>
#include <string>
#include <unordered_map>

#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
// #include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Tools/replacement_ring.pb.h"

namespace aggregate_replacement_rings {
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
  cerr << "Aggregates RplRing::sReplacementRing textproto files\n";
  cerr << " -S <fname>  resulting textproto file containing the aggregated result\n";
  cerr << " -v          verbose output\n";
  // clang-format on

  ::exit(rc);
}

int
Write(const std::unordered_map<std::string, RplRing::ReplacementRing>& smi2proto,
                IWString_and_File_Descriptor& output) {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;

  for (const auto& [_, proto] : smi2proto) {
    printer.PrintToString(proto, &buffer);
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
Write(const std::unordered_map<std::string, RplRing::ReplacementRing>& smi2proto,
      IWString& fname) {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "Write:cannot open '" << fname << "'\n";
    return 0;
  }

  return Write(smi2proto, output);
}

int
AggregateReplacementRingsLine(const RplRing::ReplacementRing& proto,
                std::unordered_map<std::string, RplRing::ReplacementRing>& smi2proto) {
  const std::string& usmi = proto.usmi();
  auto iter = smi2proto.find(usmi);

  if (iter != smi2proto.end()) {
    RplRing::ReplacementRing& existing = iter->second;
    auto n = existing.n();
    existing.set_n(n + proto.n());
    return 1;
  }

  // We could think about moving this instead of copying if efficiency matters.
  smi2proto[usmi] = proto;
  return 1;
}

int
AggregateReplacementRingsLine(const const_IWSubstring& buffer,
                std::unordered_map<std::string, RplRing::ReplacementRing>& smi2proto) {
  absl::string_view tmp(buffer.data(), buffer.length());

  RplRing::ReplacementRing proto;
  if (!google::protobuf::TextFormat::ParseFromString(tmp, &proto)) {
    cerr << "AggregateReplacementRingsLine:invalid proto\n";
    cerr << buffer << '\n';
    return 0;
  }

  return AggregateReplacementRingsLine(proto, smi2proto);

}

int
AggregateReplacementRings(iwstring_data_source& input, std::unordered_map<std::string,
                RplRing::ReplacementRing>& smi2proto) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! AggregateReplacementRingsLine(buffer, smi2proto)) {
      cerr << "AggregateReplacementRingsLine:error processing " << buffer << '\n';
      return 0;
    }
  }

  return smi2proto.size();
}

int
AggregateReplacementRings(const char* fname, std::unordered_map<std::string,
                RplRing::ReplacementRing>& smi2proto) {
  iwstring_data_source input(fname);

  if (! input.good()) {
    cerr << "AggregateReplacementRings:cannot open '" << fname << "'\n";
    return 0;
  }

  return AggregateReplacementRings(input, smi2proto);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vS:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered \n";
    Usage(1);
  }

  const int verbose = cl.option_present('v');

  if (! cl.option_present('S')) {
    cerr << "Must specify output file via the -S option\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  std::unordered_map<std::string, RplRing::ReplacementRing> smi2proto;

  for (const char* fname : cl) {
    if (! AggregateReplacementRings(fname, smi2proto)) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
    if (verbose) {
      cerr << "After reading '" << fname << " have " << smi2proto.size() <<
              " replacement rings\n";
    }
  }

  if (verbose) {
    cerr << "Read " << smi2proto.size() << " replacement rings from " << cl.size() << " files\n";
  }

  IWString fname = cl.string_value('S');

  if (! Write(smi2proto, fname)) {
    cerr << "Cannot write '" << fname << "'\n";
    return 1;
  }

  return 0;
}

}  // namespace aggregate_replacement_rings

int
main(int argc, char** argv) {
  int rc = aggregate_replacement_rings::Main(argc, argv);

  return rc;
}
