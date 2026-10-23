// Converts from nnbr::NearNeighbours to nnbr::NearNeighboursIndices.

#include <iostream>
#include <optional>
#include <string>

#include "absl/container/flat_hash_map.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/absl_hash.h"

#ifdef BUILD_BAZEL
#include "Utilities/GFP_Tools/nearneighbours.pb.h"
#else
#include "nearneighbours.pb.h"
#endif

namespace nn_id_to_ndx {

using std::cerr;
using iw_tf_data_record::TFDataReader;
using iw_tf_data_record::TFDataWriter;

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
  cerr << R"(Converts a file or serialised nnbr::NearNeighbours protos to nnbr::NearNeighboursIndices
  -smiles <fname>       smiles and ids are in separate file. Must be same order.
  -T <distance>         truncate neighbour lists to distance <distance>.
  -S <fname>            output file.
  -v                    verbose output.
)";

  exit(rc);
}

class Data {
  private:
    absl::flat_hash_map<IWString, uint32_t> _id_to_ndx;

    resizable_array_p<IWString> _smiles;
    resizable_array_p<IWString> _id;

    // During this transformaion, we can chop the neighbour lists.
    float _max_distance;
    uint64_t _discarded_for_max_distance;

    Report_Progress _report_progress;

    // we only collect statistics if verbose is set.
    int _verbose;

    uint64_t _items_read;

    // The neighbour count upon input.
    Accumulator_Int<uint64_t> _acc_nbrs_initial;
    // If we have a maximum distance, the number of neighbours written.
    Accumulator_Int<uint64_t> _acc_nbrs_final;

  // private functions.
    uint32_t ReadSmiles(iwstring_data_source& input);
    bool ReadSmilesRecord(const const_IWSubstring& buffer);
    bool GetSmilesIds(const nnbr::NearNeighbours& proto);
    uint32_t GetSmilesIds(TFDataReader& input);
    bool EstablishCrossReferences(const IWString& smiles,
                                 IWString& id);

    int Process(const nnbr::NearNeighbours& proto, TFDataWriter& writer);
    uint32_t Process(TFDataReader& reader, TFDataWriter& writer);

  public:
    Data();

    int Initialise(Command_Line_v2& cl);

    uint32_t ReadSmilesIds(const IWString& fname);
    uint32_t GetSmilesIds(const IWString& fname);

    uint32_t size() const {
      return _id_to_ndx.size();
    }

    uint32_t Process(const IWString& input_fname, const IWString& output_fname);

    int Report(std::ostream& output) const;
};

Data::Data() {
  _items_read = 0;
  _max_distance = 1.1f;  // Something larger than 1.0
  _discarded_for_max_distance = 0;
  _verbose = 0;
}

int
Data::Initialise(Command_Line_v2& cl) {
  _verbose = cl.option_present('v');

  if (cl.option_present('T')) {
    if (! cl.value('T', _max_distance) || _max_distance <= 0.0f ||
        _max_distance > 1.0f) {
      cerr << "Data::Initialise:the max distance option -T must be a valid distance\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will truncate distances to " << _max_distance << '\n';
    }
  }

  if (cl.option_present('r')) {
    uint64_t rpt;
    if (! cl.value('r', rpt)) {
      cerr << "Data::Initialise:invalid report progress option (-r)\n";
      return 0;
    }
    _report_progress.set_report_every(rpt);

    if (_verbose) {
      cerr << "Will report progress every " << rpt << " items processed\n";
    }
  }

  return 1;
}

uint32_t
Data::ReadSmilesIds(const IWString& fname) {
  iwstring_data_source input(fname);

  if (! input.good()) {
    cerr << "Data::ReadSmiles:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadSmiles(input);
}

uint32_t
Data::ReadSmiles(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! ReadSmilesRecord(buffer)) {
      cerr << "Data::ReadSmiles:cannot parse '" << buffer << "'\n";
      return 0;
    }
  }

  return _id_to_ndx.size();
}

bool
Data::ReadSmilesRecord(const const_IWSubstring& buffer) {
  IWString smiles, id;
  int i = 0;
  if (! buffer.nextword(smiles, i) || ! buffer.nextword(id, i)) {
    cerr << "Data::ReadSmilesRecord:invalid input\n";
    return 0;
  }

  return EstablishCrossReferences(smiles, id);
}

bool
Data::EstablishCrossReferences(const IWString& smiles,
                               IWString& id) {
  id.truncate_at_first(' ');

  if (auto iter = _id_to_ndx.find(id); iter != _id_to_ndx.end()) {
    cerr << "Spread::EstablishCrossReferences:duplicate identifier '" << id << "'\n";
    return false;
  }

  auto s = _id_to_ndx.size();

  assert(s == _smiles.size());

  _id << new IWString(id);
  _smiles << new IWString(smiles);
  
  _id_to_ndx[id] = s;

  return true;
}

bool
Data::GetSmilesIds(const nnbr::NearNeighbours& proto) {
  IWString smiles = proto.smiles();
  IWString id = proto.name();

  return EstablishCrossReferences(smiles, id);
}

uint32_t
Data::GetSmilesIds(TFDataReader& input) {
  uint32_t items_read = 0;

  while (1) {
    std::optional<nnbr::NearNeighbours> maybe_proto =
      input.ReadProto<nnbr::NearNeighbours>();
    if (! maybe_proto) {
      return items_read;
    }

//  if (maybe_proto->nbr_size() == 0) {
//    continue;
//  }

    if (! GetSmilesIds(*maybe_proto)) {
      return 0;
    }

    ++items_read;
  }

  return items_read;
}

uint32_t
Data::GetSmilesIds(const IWString& fname) {
  TFDataReader input(fname);
  if (! input.good()) {
    cerr << "Data::GetSmilesIds:cannot open '" << fname << "'\n";
    return 0;
  }

  return GetSmilesIds(input);
}

uint32_t
Data::Process(const IWString& input_fname, const IWString& output_fname) {

  TFDataReader reader(input_fname);
  if (! reader.good()) {
    cerr << "Data::Process:cannot open input '" << input_fname << "'\n";
    return 0;
  }

  TFDataWriter writer;
  if (! writer.Open(output_fname)) {
    cerr << "Data::Process:cannot open output '" << output_fname << "'\n";
    return 0;
  }

  return Process(reader, writer);
}

uint32_t
Data::Process(TFDataReader& reader, TFDataWriter& writer) {
  uint32_t items_read = 0;

  while (1) {
    std::optional<nnbr::NearNeighbours> maybe_proto =
      reader.ReadProto<nnbr::NearNeighbours>();
    if (! maybe_proto) {
      return items_read;
    }

    ++items_read;
    if (_report_progress()) {
      cerr << "Processed " << items_read << " items\n";
    }

//  if (maybe_proto->nbr_size() == 0) {
//    continue;
//  }

    if (! Process(*maybe_proto, writer)) {
      return 0;
    }
  }

  return items_read;
}

int
Data::Process(const nnbr::NearNeighbours& proto, TFDataWriter& writer) {
  ++_items_read;

  nnbr::NearNeighboursIndices to_write;
  to_write.set_smiles(proto.smiles());
  to_write.set_name(proto.name());
  to_write.mutable_nbr()->Reserve(proto.nbr().size());

  if (_verbose) {
    _acc_nbrs_initial.extra(proto.nbr_size());
  }

  for (const auto& nbr : proto.nbr()) {
    if (nbr.dist() > _max_distance) {
      ++_discarded_for_max_distance;
      continue;
    }

    auto* s = to_write.mutable_nbr()->Add();
    const auto iter = _id_to_ndx.find(nbr.id());
    if (iter == _id_to_ndx.end()) {
      cerr << "Data::Process:cannot find '" << nbr.id() << "' in hash\n";
      return 0;
    }
    s->set_id(iter->second);
    s->set_dist(nbr.dist());
  }

  if (_verbose) {
    _acc_nbrs_final.extra(to_write.nbr_size());
  }

  if (! writer.WriteSerializedProto<nnbr::NearNeighboursIndices>(to_write)) {
    cerr << "Data::Process:error writing " << to_write.ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}

int
Data::Report(std::ostream& output) const {
  output << "Contains " << _id_to_ndx.size() << " items\n";
  output << "Read " << _items_read << " items\n";
  output << "Initial nbr counts btw " << _acc_nbrs_initial.minval() << " and " <<
            _acc_nbrs_initial.maxval() << " ave " <<
            static_cast<float>(_acc_nbrs_initial.average()) << '\n';
  if (_max_distance < 1.0f) {
    cerr << "Omitted " << _discarded_for_max_distance << 
            " neighbours longer than " << _max_distance << '\n';
    output << "Final   nbr counts btw " << _acc_nbrs_final.minval() << " and " <<
            _acc_nbrs_final.maxval() << " ave " <<
            static_cast<float>(_acc_nbrs_final.average()) << '\n';
  }

  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-smiles=sfile-S=s-T=float-r=ipos");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (! cl.option_present('S')) {
    cerr << "Must specify name of output file via the -S option\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Must specify input TFDataRecord file as command line argument\n";
    Usage(1);
  }

  Data data;

  if (! data.Initialise(cl)) {
    cerr << "Cannot initialise data\n";
    return 1;
  }

  const IWString input_fname(cl[0]);

  const auto tzero = std::chrono::system_clock::now();

  if (cl.option_present("smiles")) {
    IWString fname = cl.string_value("smiles");
    if (! data.ReadSmilesIds(fname)) {
      cerr << "nn_id_to_ndx:cannot read smiles from '" << fname << "\n";
      return 1;
    }
  } else {
    if (!data.GetSmilesIds(input_fname)) {
      cerr << "nn_id_to_ndx:cannot read ids and smiles from '" << input_fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << "Read " << data.size() << " ids\n";
    auto now = std::chrono::system_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(now - tzero);
    cerr << "Reading smiles/ids took " << elapsed_seconds.count() << " seconds\n";
  }

  IWString output_fname = cl.string_value('S');

  if (! data.Process(input_fname, output_fname)) {
    cerr << "nn_id_to_ndx::error processing '" << input_fname << "'\n";
    return 1;
  }

  if (verbose) {
    data.Report(cerr);
  }

  return 0;
}

}  // namespace nn_id_to_ndx

int
main(int argc, char** argv) {
  const int rc = nn_id_to_ndx::Main(argc, argv);

  return rc;
}
