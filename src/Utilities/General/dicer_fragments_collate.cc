// Aggregate multiple dicer_data::DicerFragment protos.
// There is a working python version of this but it is too slow.

#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>

#include "absl/container/flat_hash_map.h"
#include "absl/log/initialize.h"

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace dicer_fragments_collate {

using std::cerr;

// Found that when storing a hash from unique smiles
// to proto, the memory requirements were very large,
// 9.5GB for what were 4.4GB on disk and which ended
// up 3.1GB when combined. I expected protos to be
// memory efficient. Let's see if this is better.
// Note that we do not store the unique smiles, it is the key.
// This does not make much of an improvement. The job now
// takes 7.5GB of RAM

struct DicerData {
  IWString smiles;
  std::string par;
  uint32_t n;
  uint32_t natoms;
};

class DicerFragmentsCollate {
  private:
    // When writing the result, a minumum n: value in order to be written.
    uint32_t _support;

    // Our objective is to generate a combined hash that counts
    // the instances of each unique smiles.
    absl::flat_hash_map<std::string, dicer_data::DicerFragment> _hash;

    // hopefully more memory efficient.
    absl::flat_hash_map<std::string, DicerData> _hash_data;

    // The default is to use the full proto.
    int _hash_value_is_proto;

    // Is the input a valid textproto, or does it have a leading
    // non-unique smiles.
    int _has_leading_non_unique_smiles;

    // It can be convenient to limit the number of atoms
    uint32_t _max_atoms;
    int _discarded_for_too_many_atoms;

    // If it has a leading non unique smiles, keep track of them and
    // write them at the end.
    absl::flat_hash_map<std::string, std::string> _usmi_to_non_unique;

    // If the input is serialized protos.
    int _input_is_tfdata;

    int _verbose;

    Report_Progress _report_progress;

    int _items_read;
    int _suppressed_by_support;
    uint32_t _max_count;
    uint32_t _singleton_count;

  // private functions

    int AccumulateRecord(const_IWSubstring buffer);
    int AccumulateRecord(dicer_data::DicerFragment& proto,
                const IWString& non_unique_smiles);

    int AccumulateTFDataRecord(iw_tf_data_record::TFDataReader& input);

    int InsertIntoDataHash(dicer_data::DicerFragment& proto,
                const IWString& non_unique_smiles);
    int InsertIntoProtoHash(dicer_data::DicerFragment& proto,
                const IWString& non_unique_smiles);

    int WriteFromDataHash(IWString_and_File_Descriptor& output);

    int WriteFromProtoHash(IWString_and_File_Descriptor& output);
    int WriteNonUniqueSmiles(const std::string& usmi,
                             IWString_and_File_Descriptor& output) const;

  public:
    DicerFragmentsCollate();

    int Initialise(Command_Line_v2& cl);

    int Accumulate(const char* fname);
    int Accumulate(iwstring_data_source& input);
    int AccumulateTFDataRecord(const char* fname);

    int Write(IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

DicerFragmentsCollate::DicerFragmentsCollate() {
  _has_leading_non_unique_smiles = 1;

  _verbose = 0;

  _input_is_tfdata = 0;

  _hash_value_is_proto = 1;

  _max_atoms = std::numeric_limits<uint32_t>::max();
  _discarded_for_too_many_atoms = 0;

  _items_read = 0;
  _support = 0;
  _suppressed_by_support = 0;
  _max_count = 0;
  _singleton_count = 0;
}

int
DicerFragmentsCollate::Initialise(Command_Line_v2& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('p')) {
    if (! cl.value('p', _support) || _support < 1) {
      cerr << "The support for writing option (-p) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only write examples with at least " << _support << " exemplars\n";
    }
  }
  
  if (cl.option_present("nosmi")) {
    _has_leading_non_unique_smiles = 0;
  }

  if (cl.option_present("maxat")) {
    cl.value("maxat", _max_atoms);
    if (_max_atoms < 1) {
      cerr << "Invalid -maxat value\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will discard input with more than " << _max_atoms << " atoms\n";
    }
  }

  if (cl.option_present('r')) {
    int rpt;
    if (! cl.value('r', rpt) || rpt < 1) {
      cerr << "The report every option (-r) must be a whole +ve number\n";
      return 0;
    }
    _report_progress.set_report_every(rpt);
    if (_verbose) {
      cerr << "Will report progress every " << rpt << " items read\n";
    }
  }

  if (cl.option_present("minimal")) {
    _hash_value_is_proto = 0;
    if (_verbose) {
      cerr << "Will extract only essential items from the proto\n";
    }
  }

  if (cl.option_present("tfdata")) {
    _input_is_tfdata = 1;
    if (_verbose) {
      cerr << "Input assumed to be TFDataRecord serialized protos\n";
    }
  }

  return 1;
}

int
DicerFragmentsCollate::Accumulate(const char* fname) {
  if (_input_is_tfdata) {
    return AccumulateTFDataRecord(fname);
  }

  iwstring_data_source input(fname);

  if (! input.good()) {
    cerr << "DicerFragmentsCollate::Accumulate:cannot open '" << fname << "'\n";
    return 0;
  }

  if (_verbose) {
    cerr << "DicerFragmentsCollate::Accumulate:start '" << fname << "'\n";
  }

  return Accumulate(input);
}

int
DicerFragmentsCollate::Accumulate(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! AccumulateRecord(buffer)) {
      cerr << "DicerFragmentsCollate::Accumulate:cannot process line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }
  }

  if (_hash_value_is_proto) {
    return _hash.size();
  } else {
    return _hash_data.size();
  }
}

// Note local copy of argument.
int
DicerFragmentsCollate::AccumulateRecord(const_IWSubstring buffer) {
  int i = 0;
  const_IWSubstring token;

  IWString non_unique_smiles;
  if (_has_leading_non_unique_smiles) {
    buffer.nextword(non_unique_smiles, i);
    buffer += i;
    i = 0;
  }

  using google::protobuf::io::ArrayInputStream;
  google::protobuf::io::ArrayInputStream zero_copy_array(buffer.data(), buffer.nchars());

  dicer_data::DicerFragment proto;
  if (!google::protobuf::TextFormat::Parse(&zero_copy_array, &proto)) {
    cerr << "DicerFragmentsCollate:AccumulateRecord:cannot parse proto " << buffer << '\n';
    return 0;
  }

  return AccumulateRecord(proto, non_unique_smiles);
}

int
DicerFragmentsCollate::AccumulateTFDataRecord(const char* fname) {
  iw_tf_data_record::TFDataReader input(fname);
  if (! input.good()) {
    cerr << "DicerFragmentsCollate::AccumulateTFDataRecord:cannot open '" << fname << "'\n";
    return 0;
  }

  return AccumulateTFDataRecord(input);
}

int
DicerFragmentsCollate::AccumulateTFDataRecord(iw_tf_data_record::TFDataReader& input) {

  // Never used here.
  static IWString non_unique_smiles;

  while (true) {
    std::optional<dicer_data::DicerFragment> maybe_proto =
       input.ReadProto<dicer_data::DicerFragment>();

    if (maybe_proto) {
      AccumulateRecord(*maybe_proto, non_unique_smiles);
    } else if (input.eof()) {
      return 1;
    } else {
      cerr << "DicerFragmentsCollate::AccumulateTFDataRecord: error reading\n";
      return 0;
    }
  }

  return 1;
}

// Note `proto` is NOT const, we move it.
int
DicerFragmentsCollate::AccumulateRecord(dicer_data::DicerFragment& proto,
                const IWString& non_unique_smiles) {
  ++_items_read;
  if (_report_progress()) {
    cerr << _items_read << " items read, hash size ";
    if (_hash_value_is_proto) {
      cerr << _hash.size();
    } else {
      cerr << _hash_data.size();
    }
    cerr << '\n';
  }

  if (proto.nat() > _max_atoms) {
    ++_discarded_for_too_many_atoms;
    return 1;
  }

  if (_hash_value_is_proto) {
    return InsertIntoProtoHash(proto, non_unique_smiles);
  } else {
    return InsertIntoDataHash(proto, non_unique_smiles);
  }
}

int
DicerFragmentsCollate::InsertIntoProtoHash(dicer_data::DicerFragment& proto,
                const IWString& non_unique_smiles) {
  auto iter = _hash.find(proto.smi());
  if (iter != _hash.end()) {
    const int n = iter->second.n();
    iter->second.set_n(n + 1);
    return 1;
  }

  if (! non_unique_smiles.empty()) {
    std::string tmp(non_unique_smiles.data(), non_unique_smiles.size());

    _usmi_to_non_unique.emplace(std::make_pair(proto.smi(), std::move(tmp)));
  }

  // The key is the unique smiles, so clear from `proto`.
  std::string key = proto.smi();

  proto.clear_smi();

  _hash.emplace(std::make_pair(std::move(key), std::move(proto)));

  return 1;
}

int
DicerFragmentsCollate::InsertIntoDataHash(dicer_data::DicerFragment& proto,
                const IWString& non_unique_smiles) {
  auto iter = _hash_data.find(proto.smi());
  if (iter != _hash_data.end()) {
    const int n = iter->second.n;
    iter->second.n = n + 1;
    iter->second.natoms = proto.nat();
    return 1;
  }

  DicerData data;
  if (! non_unique_smiles.empty()) {
    data.smiles = non_unique_smiles;
  }
  data.par = proto.par();
  data.n = proto.n();

  _hash_data.insert(std::make_pair(proto.smi(), std::move(data)));

  return 1;
}

int
DicerFragmentsCollate::Write(IWString_and_File_Descriptor& output) {
  if (_hash_value_is_proto) {
    return WriteFromProtoHash(output);
  } else {
    return WriteFromDataHash(output);
  }
}

int
DicerFragmentsCollate::WriteFromProtoHash(IWString_and_File_Descriptor& output) {

  static google::protobuf::TextFormat::Printer printer;  
  printer.SetSingleLineMode(true);

  std::string buffer;
  _max_count = 0;
  _singleton_count = 0;

  for (auto& [usmi, proto] : _hash) {
    if (proto.n() == 1) {
      ++_singleton_count;
    }

    if (_support > 0 && proto.n() < _support) {
      ++_suppressed_by_support;
      continue;
    }

    if (proto.n() > _max_count) {
      _max_count = proto.n();
    }

    proto.set_smi(usmi);

    if (! printer.PrintToString(proto, &buffer)) {
      cerr << "DicerFragmentsCollate::cannot print '" << proto.ShortDebugString() << "'\n";
      continue;
    }

    if (buffer.ends_with(' ')) {
      buffer.pop_back();
    }

    if (_has_leading_non_unique_smiles) {
      if (! WriteNonUniqueSmiles(proto.smi(), output)) {
        return 0;
      }
    }

    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(8192);
  }

  output.flush();

  return 1;
}

int
DicerFragmentsCollate::WriteFromDataHash(IWString_and_File_Descriptor& output) {

  static google::protobuf::TextFormat::Printer printer;  
  printer.SetSingleLineMode(true);

  std::string buffer;
  _max_count = 0;
  _singleton_count = 0;

  for (const auto& [usmi, data] : _hash_data) {
    if (data.n == 1) {
      ++_singleton_count;
    }

    if (_support > 0 && data.n < _support) {
      ++_suppressed_by_support;
      continue;
    }

    if (data.n > _max_count) {
      _max_count = data.n;
    }

    if (! data.smiles.empty()) {
      output << data.smiles;
    }
    output << " smi: " << usmi;
    output << " par: " << data.par;
    output << " n: " << data.n;
    output << " nat: " << data.natoms;
    output << '\n';

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

int
DicerFragmentsCollate::WriteNonUniqueSmiles(const std::string& usmi,
                                            IWString_and_File_Descriptor& output) const {
  const auto iter = _usmi_to_non_unique.find(usmi);
  if (iter == _usmi_to_non_unique.end()) {
    cerr << "DicerFragmentsCollate::WriteNonUniqueSmiles:no smiles for '" << usmi << "'\n";
    return 0;
  }

  output << iter->second << ' ';

  return 1;
}

int
DicerFragmentsCollate::Report(std::ostream& output) const {
  output << "DicerFragmentsCollate::Report: read " << _items_read <<
            " items, ";
  if (_max_atoms != std::numeric_limits<uint32_t>::max()) {
    output << "discarded " << _discarded_for_too_many_atoms <<
              " inputs for having more than " << _max_atoms << " atoms\n";
  }
  if (_hash_value_is_proto) {
    output << _hash.size();
  } else {
    output << _hash_data.size();
  }
  output << " unique smiles\n";

  if (_support > 0) {
    output << _suppressed_by_support << " suppressed by support requirement " << _support << '\n';
  }
  output << _singleton_count << " singletons, max count " << _max_count << '\n';

  return 1;
}

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
  cerr << "Aggregates multiple dicer_data::DicerFragment text proto files\n";
  cerr << " -p <support>        minimum support level (n: value) for inclusion\n";
  cerr << " -nosmi              each record does not contain a leading non-unique smiles\n";
  cerr << " -maxat <n>          discard input that contains more than <maxat> atoms\n";
  cerr << " -r <n>              report progress every <n> items read\n";
  cerr << " -minimal            extract only the essential information from the protos\n";
  cerr << " -tfdata             data is TFDataRecord serialized protos\n";
  cerr << " -v                  verbose output\n";
  // clang-format off

  ::exit(rc);
}

int
Main(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-p=ipos-nosmi-r=ipos-minimal-tfdata-maxat=ipos");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options present\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  DicerFragmentsCollate doit;
  if(! doit.Initialise(cl)) {
    cerr << "Cannot initialise DicerFragmentsCollate\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(0);
  }

  for (const char* fname : cl) {
    if (! doit.Accumulate(fname)) {
      cerr << "Cannot process '" << fname << "'\n";
      return 1;
    }
  }

  IWString_and_File_Descriptor output(1);
  if (! doit.Write(output)) {
    cerr << "Cannot write\n";
  }

  output.flush();

  if (verbose) {
    doit.Report(cerr);
  }

  return 0;
}

}  // namespace dicer_fragments_collate

int
main(int argc, char** argv) {
  absl::InitializeLog();

  int rc = dicer_fragments_collate::Main(argc, argv);

  return rc;
}
