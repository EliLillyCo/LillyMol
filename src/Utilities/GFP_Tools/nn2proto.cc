// Convert .nn files to proto form

#include <iostream>
#include <limits>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION 1
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Utilities/GFP_Tools/nearneighbours.pb.h"
#include "Utilities/GFP_Tools/nndata.h"

namespace nn2proto {
using std::cerr;

IWString smiles_tag = "$SMI<";
IWString identifier_tag = "PCN<";
IWString distance_tag = "DIST<";

// Some tags are ignored here.
IWString cluster_tag = "CLUSTER<";
IWString csize_tag = "CSIZE<";

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
  cerr << "Converts a .nn file to proto form\n";
  cerr << " -o <sep>      set token separator (default ,)\n";
  cerr << " -n <nbrs>     only include <nbrs> neighbours\n";
  cerr << " -s ...        sort the targets, enter '-s help' for info\n";
  cerr << " -z            do not write molecules with no neighbours\n";
  cerr << " -l            input has come from leader\n";
  cerr << " -S <fname>    write serialized ListOfNearNeighbours proto to <fname>\n";
  cerr << " -T <fname>    write TFDataRecord file of serialized Nbr protos to <fname>\n";
  cerr << " -v            verbose output\n";
// clang-format on

  ::exit(rc);
}

class Neighbour {
 private:
  IWString _smiles;
  IWString _id;
  float _distance;

 public:
  Neighbour();

  int Build(iwstring_data_source& input);

  float distance() const {
    return _distance;
  }

  int Write(char sep, IWString_and_File_Descriptor& output) const;

  const IWString& smiles() const {
    return _smiles;
  }

  const IWString& id() const {
    return _id;
  }
};

Neighbour::Neighbour() {
  _distance = -1.0f;
}

int
Neighbour::Write(char sep, IWString_and_File_Descriptor& output) const {
  output << _smiles;
  output << sep << _id;
  output << sep << _distance;
  return 1;
}

// We assume that `input` is pointing at the start of a neighbour
// and we expect to read smiles,id,distance, but we do not assume
// that order.
int
Neighbour::Build(iwstring_data_source& input) {
  const_IWSubstring buffer;

  // We need 3 components initialised, and then we are done.
  int components_initialised = 0;
  while (input.next_record(buffer)) {
    // cerr << "Neighbour::Build:reading '" << buffer << "'\n";
    if (_smiles.empty() && buffer.starts_with(smiles_tag)) {
      buffer.remove_leading_chars(smiles_tag.length());
      buffer.chop();
      _smiles = buffer;
      ++components_initialised;
    } else if (_id.empty() && buffer.starts_with(identifier_tag)) {
      buffer.remove_leading_chars(identifier_tag.length());
      buffer.chop();
      _id = buffer;
      ++components_initialised;
    } else if (_distance < 0.0f && buffer.starts_with(distance_tag)) {
      buffer.remove_leading_chars(distance_tag.length());
      buffer.chop();
      if (!buffer.numeric_value(_distance) || _distance < 0.0f || _distance > 1.0f) {
        cerr << "Neighbour::Build:invalid distance '" << buffer << "'\n";
        return 0;
      }
      ++components_initialised;
    } else {
      cerr << "Neighbour::Build:invalid input, line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }
    // cerr << "components_initialised " <<components_initialised << '\n';
    if (components_initialised == 3) {
      return 1;
    }
  }

  cerr << "Neighbour::Build:read beyond nndata\n";
  return 0;
}

// An individual entry in the NN file.
// smiles, id and a list of neighbours.
class NNData {
 private:
  IWString _smiles;
  IWString _id;
  resizable_array_p<Neighbour> _neighbour;

 public:
  int Build(int from_leader, iwstring_data_source& input);

  int number_neighbours() const {
    return _neighbour.number_elements();
  }

  const Neighbour* operator[](int ndx) const {
    return _neighbour[ndx];
  }

  nnbr::NearNeighbours ToProto() const;

  int Write(IWString_and_File_Descriptor& output) const;
};

int
FetchRecord(iwstring_data_source& input, const IWString& tag, IWString& destination) {
  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    cerr << "FetchRecord:cannot fetch\n";
    return 0;
  }
  if (!buffer.starts_with(tag)) {
    cerr << "FetchRecord:tag mismatch, expected '" << tag << "' got " << buffer << '\n';
    return 0;
  }

  buffer.remove_leading_chars(tag.length());
  buffer.chop();
  destination = buffer;

  return 1;
}

// Read from `input` until we get a CSIZE tag. Then parse
// that to set `cluster_size`.
int
DetermineClusterSize(iwstring_data_source& input, int& cluster_size) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with('|')) {
      cerr << "DetermineClusterSize:csize data not found\n";
      return 0;
    }
    if (!buffer.starts_with(csize_tag)) {
      continue;
    }

    buffer.remove_leading_chars(csize_tag.length());
    buffer.chop();
    if (!buffer.numeric_value(cluster_size) || cluster_size < 0) {
      cerr << "DetermineClusterSize:invalid cluster size '" << buffer << "'\n";
      return 0;
    }

    return 1;
  }

  // Should never come here.
  return 0;
}

int
NNData::Build(int from_leader, iwstring_data_source& input) {
  if (!FetchRecord(input, smiles_tag, _smiles) ||
      !FetchRecord(input, identifier_tag, _id)) {
    cerr << "NNData::Build:cannot read first two records\n";
    return 0;
  }

  if (from_leader) {
    int cluster_size = 0;
    if (!DetermineClusterSize(input, cluster_size)) {
      cerr << "NNData::Build:cannot determine cluster size\n";
      return 0;
    }

    if (cluster_size == 0) {
      return 1;
    }
  }

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer == '|') {
      return 1;
    }
    input.push_record();
    std::unique_ptr<Neighbour> nbr = std::make_unique<Neighbour>();
    if (!nbr->Build(input)) {
      cerr << "NNData::Build:fatal error at line " << input.lines_read() << '\n';
      return 0;
    }
    _neighbour << nbr.release();
  }

  return 1;
}

// Write the neighbour data. Tokens separated by `sep`
// and we must produce `max_nbrs` columns.
int
NNData::Write(IWString_and_File_Descriptor& output) const {
  nnbr::NearNeighbours proto = ToProto();

  return gfp::WriteNNData(proto, output);
}

nnbr::NearNeighbours
NNData::ToProto() const {
  nnbr::NearNeighbours proto;  // to be returned.
  proto.set_name(_id.data(), _id.length());
  proto.set_smiles(_smiles.data(), _smiles.length());

  for (const Neighbour* neighbour : _neighbour) {
    nnbr::Nbr* nbr = proto.add_nbr();
    const IWString& smi = neighbour->smiles();
    const IWString& id = neighbour->id();
    nbr->set_smi(smi.data(), smi.length());
    nbr->set_id(id.data(), id.length());
    nbr->set_dist(neighbour->distance());
  }

  return proto;
}

enum class SortType {
  kNone = 0,
  kShortestDistance = 1,
  kLongestDistance = 2,
  kNumberNeighbours = 3,
};

class NN2Proto {
 private:
  int _verbose;
  // The output separator.
  char _sep;

  // We can limit the number of neighbours output
  int _neighbours_to_write;

  // All the nearest neighbour data read from the input.
  resizable_array_p<NNData> _nn_data;

  SortType _sort_type;

  // If the -z option is present, do not write a record
  // if there are zero neighbours.
  int _write_molecules_with_zero_neighbours;

  // If the data has come from leader, we may have zero
  // sized clusters, which are hard to parse.
  int _from_leader;

  // private functions.

  int Sort(SortType& s);
  int SortByNumberNeighbours();
  int SortByClosestDist();

 public:
  NN2Proto();

  int Initialise(Command_Line& cl);

  int Accumulate(int from_leader, const char* fname);
  int Accumulate(int from_leader, iwstring_data_source& input);

  int number_targets() const {
    return _nn_data.number_elements();
  }

  int SortIfRequested() {
    return Sort(_sort_type);
  }

  int from_leader() const {
    return _from_leader;
  }

  // Writes text_format to `output`.
  int Write(IWString_and_File_Descriptor& output) const;

  int WriteBinary(IWString& fname) const;

  int WriteTFDataRecord(iw_tf_data_record::TFDataWriter& output) const;
};

NN2Proto::NN2Proto() {
  _sep = ',';
  _sort_type = SortType::kNone;
  _neighbours_to_write = std::numeric_limits<int>::max();
  _write_molecules_with_zero_neighbours = 1;
  _from_leader = 0;
}

void
DisplaySortOptions(std::ostream& output) {
  output << " -s nbrs      sort by the number of neighbours\n";
  output << " -s dist      sort by the closest distance\n";
  exit(0);
}

int
NN2Proto::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('o')) {
    IWString o = cl.string_value('o');
    char_name_to_char(o);
    _sep = o[0];
    if (_verbose) {
      cerr << "Token separator '" << _sep << "'\n";
    }
  }

  if (cl.option_present('s')) {
    const IWString s = cl.string_value('s');
    if (s == "none") {
      _sort_type = SortType::kNone;
    } else if (s == "nbrs") {
      _sort_type = SortType::kNumberNeighbours;
    } else if (s == "dist") {
      _sort_type = SortType::kShortestDistance;
    } else if (s == "maxdist") {
      _sort_type = SortType::kLongestDistance;
    } else if (s == "help") {
      DisplaySortOptions(cerr);
    } else {
      cerr << "NN2Proto:Initialise:unrecognised sort type '" << s << "'\n";
      DisplaySortOptions(cerr);
    }
  }

  if (cl.option_present('l')) {
    _from_leader = 1;
    if (_verbose) {
      cerr << "Assuming data is from leader\n";
    }
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', _neighbours_to_write) || _neighbours_to_write < 0) {
      cerr << "NN2Proto::Initialise:the number of neighbours to write (-n) must be a "
              "whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will write a max of " << _neighbours_to_write << " neighbours\n";
    }
  }

  if (cl.option_present('z')) {
    _write_molecules_with_zero_neighbours = 0;
    if (_verbose) {
      cerr << "MOlecules with zero neighbours will not be written\n";
    }
  }

  return 1;
}

int
NN2Proto::Accumulate(int from_leader, const char* fname) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "NN2Proto::Accumulate:cannot open '" << fname << "'\n";
    return 0;
  }

  return Accumulate(from_leader, input);
}

int
NN2Proto::Accumulate(int from_leader, iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    std::unique_ptr<NNData> nn = std::make_unique<NNData>();
    input.push_record();

    if (!nn->Build(from_leader, input)) {
      cerr << "Cannot fetch start of nbr list\n";
      return 0;
    }

    _nn_data << nn.release();
  }

  if (_nn_data.empty()) {
    cerr << "NN2Proto:Accumulate:no data\n";
    return 0;
  }

  return 1;
}

int
NN2Proto::Sort(SortType& sort_type) {
  switch (sort_type) {
    case SortType::kNone:
      return 1;
    case SortType::kNumberNeighbours:
      return SortByNumberNeighbours();
    case SortType::kShortestDistance:
      return SortByClosestDist();
    default:
      cerr << "Unimplemented sort\n";
      return 0;
  }
}

int
NN2Proto::SortByNumberNeighbours() {
  _nn_data.iwqsort_lambda([](const NNData* n1, const NNData* n2) {
    if (n1->number_neighbours() < n2->number_neighbours()) {
      return -1;
    } else if (n1->number_neighbours() > n2->number_neighbours()) {
      return 1;
    } else {
      return 0;
    }
  });

  return 1;
}

// Sorting by closest neighbour is complicate by the fact
// that there may be no neighbours.
int
NN2Proto::SortByClosestDist() {
  _nn_data.iwqsort_lambda([](const NNData* n1, const NNData* n2) {
    float d1 = std::numeric_limits<float>::max();
    if (n1->number_neighbours() > 0) {
      d1 = (*n1)[0]->distance();
    }
    float d2 = std::numeric_limits<float>::max();
    if (n2->number_neighbours() > 0) {
      d2 = (*n2)[0]->distance();
    }
    if (d1 < d2) {
      return -1;
    } else if (d1 > d2) {
      return 1;
    } else {
      return 0;
    }
  });

  return 1;
}

int
NN2Proto::Write(IWString_and_File_Descriptor& output) const {
  int targets_with_no_neighbours = 0;

  for (const NNData* nn_data : _nn_data) {
    if (!_write_molecules_with_zero_neighbours && nn_data->number_neighbours() == 0) {
      ++targets_with_no_neighbours;
      continue;
    }
    nn_data->Write(output);
  }

  if (_verbose && !_write_molecules_with_zero_neighbours) {
    cerr << "Skipped " << targets_with_no_neighbours << " targets with no neighbours\n";
  }

  return 1;
}

int
NN2Proto::WriteTFDataRecord(iw_tf_data_record::TFDataWriter& output) const {
  for (const NNData* nn_data : _nn_data) {
    const nnbr::NearNeighbours proto = nn_data->ToProto();
    output.WriteSerializedProto<nnbr::NearNeighbours>(proto);
  }

  return 1;
}

int
NN2Proto::WriteBinary(IWString& fname) const {
  int targets_with_no_neighbours = 0;

  nnbr::ListOfNearNeighbours proto;

  for (const NNData* nn_data : _nn_data) {
    if (!_write_molecules_with_zero_neighbours && nn_data->number_neighbours() == 0) {
      ++targets_with_no_neighbours;
      continue;
    }

    nnbr::NearNeighbours* extra = proto.add_nearneighbours();
    // Potentially expensive operation.
    *extra = std::move(nn_data->ToProto());
  }

  if (_verbose && !_write_molecules_with_zero_neighbours) {
    cerr << "Skipped " << targets_with_no_neighbours << " targets with no neighbours\n";
  }

  return iwmisc::WriteBinaryProto(proto, fname);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vo:s:n:zlS:T:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_present('v');

  if (cl.empty()) {
    cerr << "INsufficient arguments\n";
    Usage(1);
  }

  NN2Proto nn_to_csv;
  if (!nn_to_csv.Initialise(cl)) {
    cerr << "Cannot initialise NN2Proto conditions\n";
    return 1;
  }

  for (const char* fname : cl) {
    if (!nn_to_csv.Accumulate(nn_to_csv.from_leader(), fname)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << "Read data on " << nn_to_csv.number_targets() << " targets\n";
  }

  nn_to_csv.SortIfRequested();

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (! nn_to_csv.WriteBinary(fname)) {
      cerr << "Binary write to '" << fname << "' failed\n";
      return 1;
    }

    return 0;
  }

  if (cl.option_present('T')) {
    IWString fname = cl.string_value('T');
    iw_tf_data_record::TFDataWriter writer;
    if (! writer.Open(fname)) {
      cerr << "Cannot open '" << fname << "' for TFDataRecord proto output\n";
      return 1;
    }

    if (! nn_to_csv.WriteTFDataRecord(writer)) {
      cerr << "TFDataRecord write to '" << fname << "' failed\n";
      return 1;
    }

    return 0;
  }

  IWString_and_File_Descriptor output(1);

  nn_to_csv.Write(output);

  return 0;
}

}  // namespace nn2proto

int
main(int argc, char** argv) {
  int rc = nn2proto::Main(argc, argv);
  return rc;
}
