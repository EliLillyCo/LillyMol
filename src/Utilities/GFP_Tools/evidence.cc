// Examines values associated with structures and looks for
// consistency among the neighbours.
// Needs:
//   An activity file via the -A option.
//   A file of nnbr::NearNeighbours serialized protos

#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <ranges>
#include <algorithm>

#include "absl/container/flat_hash_map.h"

#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Utilities/GFP_Tools/nearneighbours.pb.h"
#include "Utilities/GFP_Tools/evidence.pb.h"

namespace evidence {

using std::cerr;
using iw_tf_data_record::TFDataReader;

// when dealing with values that need to be sorted, use this to
// keep track of the id of the computed value.
template <typename T>
struct NdxValue {
  uint32_t ndx;
  T value;
  // After a set of values have been assembled and sorted, an
  // external entity can assign a percentily to each item.
  int percentile = -1;
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
  // clang-format off
  cerr << "Computes measures of internal consistency for a series of measured values based on how\n";
  cerr << "consistent the values are across the neighbours\n";
  cerr << "Takes a single argument, a TFDataRecord serialized nnbr::NearNeighbours protos, such as\n";
  cerr << "might be produced by the -T option to nn2proto\n";
  cerr << "The following arguments are recognised\n";
  cerr << " -config <fname>     an EvidenceData::Options textproto file with options\n";
  cerr << " -A <fname>          descriptor file containing activity values for each molecule\n";
  cerr << " -smiles <fname>     smiles for each molecule, will be included in the output\n";
  cerr << " -C                  data is classlfication type (not implemented)\n";
  cerr << " -diff               for each column generated, insert an extra column with difference from actual\n";
  cerr << " -v                  verbose output\n";
  // clang-format on

  ::exit(rc);
}

template <typename T>
class Nbr {
  private:
    // As read from the proto;
    float _distance;

    // This will be read from _id_to_activity
    T _activity;

    // Used as an array index into the _item array in the Evidence class.
    int _ndx;

  public:
    Nbr();

    void Set(float d, T a, int n) {
      _distance = d;
      _activity = a;
      _ndx = n;
    }

    float distance() const  {
      return _distance;
    }
    T activity() const  {
      return _activity;
    }
    int ndx() const  {
      return _ndx;
    }

};

template <typename T>
Nbr<T>::Nbr() {
  _distance = 0.0;
  _activity = {};
  _ndx = -1;
}

// Values get stored along with their percentile in the set.
template <typename T>
struct ValuePercentile {
  T value;
  int percentile;
};

template <typename T>
class Item {
  private:
    IWString _id;
    T _activity;

    int _number_nbrs;
    Nbr<T>* _nbrs;

    // As we do calculations, we append those to this array and at the 
    // end of calculations, we will write these results.
    resizable_array<float> _results;

    // Most of the results will have an associated percentile.
    // TODO:ianwatson should make a struct to hold result and percentile
    resizable_array<int> _percentile;

    // Some results may be undefined - zero nbr list for example.
    // It will be externally specified to some value that will not
    // collide with any actual value.
    T _undefined_value;
    static constexpr char kMissing = '.';

  // private functions

    void MaybeAppendDiff(bool append_diff);
    void AppendUndefined(bool append_diff);

  public:
    Item();
    ~Item();

    void set_undefined_value(T s) {
      _undefined_value = s;
    }

    void set(const IWString& s, T v) {
      _id = s;
      _activity = v;
    }

    const IWString& id() const {
      return _id;
    }
    const T activity() const {
      return _activity;
    }
    int number_nbrs() const {
      return _number_nbrs;
    }

    int PrintNbrs(std::ostream& output) const;

    int AllocateNbrs(int nbrs);

    void AddPercentile(int p) {
      _percentile << p;
    }

    // Note that we make no attempt to shorten _nbrs. Could save some
    // memory that way...
    int TruncateNbrList(int nbrs) {
      _number_nbrs = nbrs;
      return 1;
    }

    int AddNbr(int ndx, float distance, T activity, int uid);

    const Nbr<T> nbr(int ndx) const {
      return _nbrs[ndx];
    }

    // The value for `k` nearest neighbours.
    void AddKnn(bool append_diff, int k);
    // Among the k neares neighbours, the closest value, regardless of
    // how far apart the items are wrt distance.
    void AddClosest(bool append_diff, int k);
    // Add the shortest distance
    void AddShortestDistance(bool append_diff);

    void AddAverageWithinDistance(bool append_diff, float distance);

    // Add a piecewise linear weighted model.
    int AddPiecewiseLinear(bool append_diff, const EvidenceData::PiecewiseLinear& piecewise_linear);

    // Weighted average value, for the closest `nbrs` neighbours.
    int AddWeightedAverage(bool append_diff, int nbrs);

    // Return the `ndx` item from _results.
    float result(int ndx) const {
      return _results[ndx];
    }

    // Once all calculations have been accumulated, write the _results
    // array to `output`.
    int WriteResults(char sep, IWString_and_File_Descriptor& output) const;
};

template <typename T>
Item<T>::Item() {
  _activity = {};
  _number_nbrs = 0;
  _nbrs = nullptr;
  _undefined_value = 0;
}

template <typename T>
Item<T>::~Item() {
  if (_nbrs != nullptr) {
    delete [] _nbrs;
  }
}

template <typename T>
int
Item<T>::AllocateNbrs(int nbrs) {
  if (nbrs == 0) {
    return 1;
  }

  _nbrs = new Nbr<T>[nbrs];
  _number_nbrs = nbrs;

  return 1;
}

template <typename T>
int
Item<T>::AddNbr(int ndx, float distance, T activity, int uid) {
  _nbrs[ndx].Set(distance, activity, uid);

  return 1;
}

template <typename T>
void
Item<T>::AppendUndefined(bool append_diff) {
  _results << _undefined_value;

  if (append_diff) {
    _results << _undefined_value;
  }
}

// If `append_diff` is true, append the diff between the last item in _results
// to _results.
template <typename T>
void
Item<T>::MaybeAppendDiff(bool append_diff) {
  if (! append_diff) {
    return;
  }

  const float last_result = _results.back();
  if (last_result == _undefined_value) {
    _results << _undefined_value;
    return;
  }

  _results << (last_result - _activity);
}

template <typename T>
void
Item<T>::AddShortestDistance(bool append_diff) {
  if (_number_nbrs == 0) {
    AppendUndefined(append_diff);
    return;
  }

  _results << _nbrs[0].distance();

  if (! append_diff) {
    return;
  }

  _results << (_nbrs[0].activity() - _activity);
}

template <typename T>
void
Item<T>::AddKnn(bool append_diff, int k) {

  if (_number_nbrs == 0) {
    AppendUndefined(append_diff);
    return;
  }

  if (k > _number_nbrs) {
    k = _number_nbrs;
  }

  float result = 0.0;
  for (int i = 0; i < k; ++i) {
    result += _nbrs[i].activity();
  }

  result = result / static_cast<float>(k);

  _results << result;

  MaybeAppendDiff(append_diff);
}

template <typename T>
void
Item<T>::AddClosest(bool append_diff, int k) {
  if (_number_nbrs == 0) {
    AppendUndefined(append_diff);
    return;
  }

  if (k > _number_nbrs) {
    k = _number_nbrs;
  }

  float min_diff = abs(_activity - _nbrs[0].activity());
  float closest_value = _nbrs[0].activity();

  for (int i = 1; i < k; ++i) {
    float d = abs(_activity - _nbrs[i].activity());
    if (d < min_diff) {
      min_diff = d;
      closest_value = _nbrs[i].activity();
    }
  }

  _results << closest_value;

  MaybeAppendDiff(append_diff);
}

template <typename T>
void
Item<T>::AddAverageWithinDistance(bool append_diff, float distance) {
  if (_number_nbrs == 0) {
    AppendUndefined(append_diff);
    return;
  }

  float tot = 0.0;
  for (int i = 0; i < _number_nbrs; ++i) {
    tot += _nbrs[i].activity();
  }

  _results << tot / static_cast<float>(_number_nbrs);

  MaybeAppendDiff(append_diff);
}

template <typename T>
int
Item<T>::AddPiecewiseLinear(bool append_diff,
                            const EvidenceData::PiecewiseLinear& piecewise_linear) {
  if (_number_nbrs == 0) {
    AppendUndefined(append_diff);
    return 1;
  }

  float tot = 0.0;
  float sum_weights = 0.0;

  for (int i = 0; i < _number_nbrs; ++i) {
    if (_nbrs[i].distance() <= piecewise_linear.min()) {
      tot += _nbrs[i].activity();
      sum_weights += 1.0f;
      continue;
    }

    if (_nbrs[i].distance() >= piecewise_linear.max()) {
      break;
    }

    // Need to interpolate.
    float w = 1.0f - (_nbrs[i].distance() - piecewise_linear.min()) /
              (piecewise_linear.max() - piecewise_linear.min());
    tot += _nbrs[i].activity() * w;
    sum_weights += w;
  }

  if (sum_weights == 0.0f) {
    _results << _undefined_value;
  } else {
    _results << (tot / sum_weights);
  }

  MaybeAppendDiff(append_diff);

  return 1;
}

template <typename T>
int
Item<T>::AddWeightedAverage(bool append_diff, int nbrs) {
  if (_number_nbrs == 0) {
    AppendUndefined(append_diff);
    return 1;
  }

  if (nbrs > _number_nbrs) {
    nbrs = _number_nbrs;
  }

  float sum_weights = 0.0f;
  float tot = 0.0f;
  for (int i = 0; i < nbrs; ++i) {
    float w = 1.0f - _nbrs[i].distance();
    tot += _nbrs[i].activity() * w;
    sum_weights += w;
  }

  _results << (tot / sum_weights);

  MaybeAppendDiff(append_diff);

  return 1;
}

template <typename T>
int
Item<T>::WriteResults(char sep, IWString_and_File_Descriptor& output) const {
  for (float r : _results) {
    output << sep;
    if (r == _undefined_value) {
      output << kMissing;
    } else {
      output << r;
    }
  }

  return 1;
}

class Evidence {
  private:
    int _verbose;

    EvidenceData::Options _config;

    // We need some hashes while ingesting the data.
    //IW_STL_Hash_Map_float _id_to_activity;
    //IW_STL_Hash_Map_int _id_to_ndx;
    absl::flat_hash_map<std::string, float> _id_to_activity;
    absl::flat_hash_map<std::string, int> _id_to_ndx;

    // A vector of all the Items in the input. Constructed from the -A file
    Item<float>* _item;
    int _number_items;

    int _is_classification;

    // Our output is a descriptor file.
    resizable_array_p<IWString> _header;

    // We can append after each calculated result its difference from the true
    // value.
    bool _append_differences;

    // If a smiles file has been specified.
    absl::flat_hash_map<IWString, IWString> _id_to_smiles;

    char _input_separator;
    char _output_separator;

  // private functions
   int ReadActivityDataRecord(const const_IWSubstring& buffer,
                IW_STL_Hash_Map_String& id_to_activity_string);
    int ReadActivityData(IWString& fname,
                IW_STL_Hash_Map_String& id_to_activity_string);
    int ReadActivityData(iwstring_data_source& input,
                IW_STL_Hash_Map_String& id_to_activity_string);

    int ReadSmiles(IWString& fname);
    int ReadSmiles(iwstring_data_source& input);
    int ReadSmilesRecord(const const_IWSubstring& buffer);

    int ToFloat(const IW_STL_Hash_Map_String& id_to_activity_string);

    int SetupItems();

    int ReadNbrList(TFDataReader& input);

    int Initialise(const nnbr::NearNeighbours& proto);

    int AddToHeader(const IWString& s);
    int AddToHeader(const char* stem, int ndx);

    int PrintNbrs(std::ostream& output) const;

    void MaybeWriteSmiles(const IWString& id, IWString_and_File_Descriptor& output) const;

    int Knn();
    int Knn(int k);
    int ClosestValue();
    int ClosestValue(int k);
    int ShortestDistance();
    int Closest();
    int AverageWithinDistance();
    int AverageWithinDistance(float distance);
    int PiecewiseWeightedAverage();
    int PiecewiseWeightedAverage(const EvidenceData::PiecewiseLinear& piecewise_linear);
    int WeightedAverage();
    int WeightedAverage(int nbrs);

    int Process(float* values);

  public:
    Evidence();
    ~Evidence();

    int Initialise(Command_Line_v2& cl);

    int ReadNbrList(const char* fname);

    int Process();

    int WriteResults(IWString_and_File_Descriptor& output) const;
};

Evidence::Evidence() {
  _verbose = 0;
  _is_classification = 0;

  _append_differences = false;

  _item = nullptr;
  _number_items = 0;

  _input_separator = ' ';
  _output_separator = ' ';
}

Evidence::~Evidence() {
  if (_item != nullptr) {
    delete [] _item;
  }
}

int
Evidence::Initialise(Command_Line_v2& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present("config")) {
    IWString fname = cl.string_value("config");
    std::optional<EvidenceData::Options> maybe_proto = 
                iwmisc::ReadTextProtoCommentsOK<EvidenceData::Options>(fname);
    if (! maybe_proto) {
      cerr << "Evidence::Initialise:cannot read options from '" << fname << "'\n";
      return 0;
    }
    _config = std::move(*maybe_proto);
    if (_verbose) {
      cerr << "Options read from '" << fname << "'\n";
    }
  }

  if (cl.option_present('C')) {
    _is_classification = 1;
    if (_verbose) {
      cerr << "Data treated as classification\n";
    }
  }

  if (cl.option_present("diff")) {
    _append_differences = true;
    if (_verbose) {
      cerr << "Will append a column of differences after each computed value\n";
    }
  }

  if (! cl.option_present('A')) {
    cerr << "Evidence::Initialise:must specify the activity file via the -A option\n";
    Usage(1);
  }

  IW_STL_Hash_Map_String id_to_activity_string;

  if (cl.option_present('A')) {
    IWString fname = cl.string_value('A');
    if (! ReadActivityData(fname, id_to_activity_string)) {
      cerr << "Cannot read activity data '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << id_to_activity_string.size() << " activity/data values from '" << fname << "'\n";
    }
  }

  if (_is_classification) {
  } else if (! ToFloat(id_to_activity_string)) {
    cerr << "Cannot convert activity values to float form\n";
    return 0;
  }

  if (! SetupItems()) {
    cerr << "Evidence::Initialise:cannot initialise internal arrays\n";
    return 0;
  }

  if (cl.option_present("smiles")) {
    IWString fname = cl.string_value("smiles");
    if (! ReadSmiles(fname)) {
      cerr << "Evidence::Initialise:cannot read smiles file '" << fname << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Read " << _id_to_smiles.size() << " id->smiles relationships from '" << fname << "'\n";
    }
  }

  return 1;
}

int
Evidence::ReadActivityData(IWString& fname,
                IW_STL_Hash_Map_String& id_to_activity_string) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Evidence::ReadActivityData:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadActivityData(input, id_to_activity_string);
}

int
Evidence::ReadActivityData(iwstring_data_source& input,
                IW_STL_Hash_Map_String& id_to_activity_string) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Evidence::ReadActivityData:empty file\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! ReadActivityDataRecord(buffer, id_to_activity_string)) {
      cerr << "Evidence::ReadActivityData:cannot process '" << buffer << "'\n";
      return 0;
    }
  }

  return id_to_activity_string.size();
}

int
Evidence::ReadActivityDataRecord(const const_IWSubstring& buffer,
                IW_STL_Hash_Map_String& id_to_activity_string) {
  IWString id, activity;
  int i = 0;
  if (! buffer.nextword(id, i, _input_separator) ||
      ! buffer.nextword(activity, i, _input_separator) ||
      id.empty() || activity.empty()) {
    cerr << "Evidence::ReadActivityDataRecord:empty or invalid record\n";
    return 0;
  }

  // No checking, so if there are duplicate values the last one matters.
  id_to_activity_string[id] = activity;

  return 1;
}

int
Evidence::ReadSmiles(IWString& fname) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Evidence::ReadSmiles:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadSmiles(input);
}

int
Evidence::ReadSmiles(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! ReadSmilesRecord(buffer)) {
      cerr << "Evidence::ReadSmiles:cannot process '" << buffer << "'\n";
      return 0;
    }
  }

  return _id_to_smiles.size();
}

int
Evidence::ReadSmilesRecord(const const_IWSubstring& buffer) {
  IWString smiles;
  int i = 0;
  buffer.nextword(smiles, i);

  IWString id;
  if (! buffer.nextword(id, i) || id.empty()) {
    return 0;
  }

  _id_to_smiles.emplace(std::make_pair(std::move(id), std::move(smiles)));

  return 1;
}

// Convert the values in id_to_activity_string to _id_to_activity
int
Evidence::ToFloat(const IW_STL_Hash_Map_String& id_to_activity_string) {
  for (const auto& [id, as_string] : id_to_activity_string) {
    float v;
    if (! as_string.numeric_value(v)) {
      cerr << "Evidence::ToFloat:invalid float '" << id << "' value '" << as_string << "'\n";
      return 0;
    }

    std::string s(id.data(), id.size());
    _id_to_activity.emplace(std::make_pair(std::move(s), v));
  }

  return _id_to_activity.size();
}

// allocate the _item array and populate it with the values from
// _id_to_activity.
// While doing that, populate _id_to_ndx.
int
Evidence::SetupItems() {
  _number_items = _id_to_activity.size();
  _item = new Item<float>[_number_items];

  for (int ndx = 0; const auto& [id, value] : _id_to_activity) {
    _item[ndx].set(id, value);
    _id_to_ndx[id] = ndx;
    ++ndx;
  }

  float min_activity = std::numeric_limits<float>::max();
  float max_activity = -std::numeric_limits<float>::max();

  for (const auto& [id, value] : _id_to_activity) {
    if (value < min_activity) {
      min_activity = value;
    }
    if (value > max_activity) {
      max_activity = value;
    }
  }

  float range = max_activity - min_activity;
  if (range == 0.0f) {
    cerr << "Evidence::SetupItems:HUH range is zero " << min_activity << " to " << max_activity << '\n';
  }

  const float undefined_value = min_activity - range;
  if (_verbose) {
    cerr << "Range " << min_activity << ',' << max_activity << " undefined value " <<
            undefined_value << '\n';
  }
  for (int i = 0; i < _number_items; ++i) {
    _item[i].set_undefined_value(undefined_value);
  }

  return _id_to_ndx.size();
}

int
Evidence::AddToHeader(const IWString& s) {
  _header << new IWString(s);

  if (_append_differences) {
    IWString*  d = new IWString(s);
    *d << ".diff";
    _header << d;
  }

  return 1;
}

int
Evidence::AddToHeader(const char* stem, int ndx) {
  std::unique_ptr<IWString> tmp = std::make_unique<IWString>(stem);
  *tmp << ndx;
  _header << tmp.release();

  if (_append_differences) {
    std::unique_ptr<IWString> tmp = std::make_unique<IWString>(stem);
    *tmp << ndx << ".diff";
    _header << tmp.release();
  }

  return 1;
}

int
Evidence::ReadNbrList(const char* fname) {
  TFDataReader input(fname);
  if (! input.good()) {
    cerr << "Evidence::Process:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadNbrList(input);
}

int
Evidence::ReadNbrList(TFDataReader& input) {
  while (1) {
    std::optional<nnbr::NearNeighbours> proto = input.ReadProto<nnbr::NearNeighbours>();
    if (! proto) {
      break;
    }
    if (! Initialise(*proto)) {
      cerr << "Evidence::Process:cannot build " << proto->ShortDebugString() << '\n';
      return 0;
    }
  }

#ifdef ECHO_NBRS
  if (_verbose) {
    PrintNbrs(cerr);
  }
#endif

  return 1;
}

int
Evidence::Initialise(const nnbr::NearNeighbours& proto) {
  const auto iter_ndx = _id_to_ndx.find(proto.name());
  if (iter_ndx == _id_to_ndx.end()) {
    cerr << "Evidence::Initialise:no index for '" << proto.name() << "'\n";
    return 0;
  }

  Item<float>& item = _item[iter_ndx->second];

  const uint32_t nbrs = proto.nbr().size();

  item.AllocateNbrs(nbrs);

  for (uint32_t i = 0; i < nbrs; ++i) {
    const nnbr::Nbr& nbr = proto.nbr(i);

    if (_config.max_distance() > 0.0f && nbr.dist() > _config.max_distance()) {
      item.TruncateNbrList(i);
      return 1;
    }

    const auto iter_ndx = _id_to_ndx.find(nbr.id());
    if (iter_ndx == _id_to_ndx.end()) {
      cerr << "Evidence::Initialise:cannot find '" << nbr.id() << "'\n";
      return 0;
    }
    const auto iter_activity = _id_to_activity.find(nbr.id());
    if (iter_activity == _id_to_activity.end()) {
      cerr << "Evidence::Initialise:no activity for '" << nbr.id() << "'\n";
      return 0;
    }

    item.AddNbr(i, nbr.dist(), iter_activity->second, iter_ndx->second);
  }

  return 1;
}

int
Evidence::PrintNbrs(std::ostream& output) const {
  static const char kSep = ' ';

  output << "Evidence:data on " << _number_items << " items\n";
  for (int i = 0; i < _number_items; ++i) {
    const int n = _item[i].number_nbrs();
    output << _item[i].id() << kSep << _item[i].number_nbrs() << '\n';
    for (int j = 0; j < n; ++j) {
      const Nbr<float>& nbr = _item[i].nbr(j);
      output << j << kSep << _item[nbr.ndx()].id() << kSep << nbr.distance() << '\n';
    }
  }

  return 1;
}

int
Evidence::Process() {
  // The first result is always the shortest distance.
  ShortestDistance();

  if (_config.knn().size()) {
    Knn();
  }

  if (_config.closest_value().size()) {
    ClosestValue();
  }

  if (_config.average_within_distance().size()) {
    AverageWithinDistance();
  }

  if (_config.piecewise_linear_size() > 0) {
    PiecewiseWeightedAverage();
  }

  if (_config.weighted_average_size() > 0) {
    WeightedAverage();
  }

  // This is not needed...
  std::unique_ptr<float[]> values = std::make_unique<float[]>(_header.size());

  return Process(values.get());
}

// `ndx_and_value` must be sorted by value.
// assign the percentile value to each item based on the position
// in the sorted list.
// In the interests of simplicity, we do not handle identical values
// just do the simplest/fastest thing possible.
void
AssignPercentileRanks(NdxValue<float>* ndx_and_value,
                      int nvalues) {
  ndx_and_value[0].percentile = 0;


  for (int i = 1; i < nvalues; ++i) {
    int p = static_cast<int>(iwmisc::Fraction<float>(i, nvalues) * 100.0f);
    ndx_and_value[i].percentile = p;
  }
  ndx_and_value[nvalues - 1].percentile = 100;

  return;
}

int
Evidence::Process(float* values) {
  std::unique_ptr<NdxValue<float>[]> ndx_and_value = std::make_unique<NdxValue<float>[]>(_number_items);

  int number_columns = _header.size();
  for (int i = 0; i < number_columns; ++i) {
    for (int j = 0; j < _number_items; ++j) {
      ndx_and_value[j].ndx = j;
      ndx_and_value[j].value = _item[j].result(i);
    }

    std::sort(ndx_and_value.get(), ndx_and_value.get() + _number_items,
                        [](const NdxValue<float>& nv1,
                                const NdxValue<float>& nv2) {
                return nv1.value < nv2.value;
              });
    AssignPercentileRanks(ndx_and_value.get(), _number_items);
    for (int j = 0; j < _number_items; ++j) {
      Item<float>& item = _item[ndx_and_value[j].ndx];
      item.AddPercentile(ndx_and_value[j].percentile);
    }
  }

  for (int i = 0; i < _number_items; ++i) {
  }

  return 1;
}

int
Evidence::Knn() {
  for (int k : _config.knn()) {
    AddToHeader("KNN", k);

    Knn(k);
  }

  return 1;
}

int
Evidence::Knn(int k) {
  for (int i = 0; i < _number_items; ++i) {
    _item[i].AddKnn(_append_differences, k);
  }

  return 1;
}

int
Evidence::AverageWithinDistance() {
  // c++ 23 
  // for (auto const [ndx, distance] : std::views::enumerate(_config.average_within_distance())) {

  for (int ndx = 0; float distance : _config.average_within_distance()) {
    AddToHeader("AVDist", ndx);
    ++ndx;

    AverageWithinDistance(distance);
  }

  return 1;
}

int
Evidence::AverageWithinDistance(float distance) {
  for (int i = 0; i < _number_items; ++i) {
    _item[i].AddAverageWithinDistance(_append_differences, distance);
  }

  return 1;
}

int
Evidence::ClosestValue() {
  for (int k : _config.knn()) {
    AddToHeader("Nearest", k);

    ClosestValue(k);
  }

  return 1;
}

int
Evidence::ClosestValue(int k) {
  for (int i = 0; i < _number_items; ++i) {
    _item[i].AddClosest(_append_differences, k);
  }

  return 1;
}

int
Evidence::ShortestDistance() {
  AddToHeader("Dmin");

  for (int i = 0; i < _number_items; ++i) {
    _item[i].AddShortestDistance(_append_differences);
  }

  return 1;
}

int
Evidence::PiecewiseWeightedAverage() {

  for (int ndx = 0; const EvidenceData::PiecewiseLinear& pwl : _config.piecewise_linear()) {
    AddToHeader("Pwlinear", ndx);
    ++ndx;
    Evidence::PiecewiseWeightedAverage(pwl);
  }

  return 1;
}

int
Evidence::PiecewiseWeightedAverage(const EvidenceData::PiecewiseLinear& piecewise_linear) {
  for (int i = 0; i < _number_items; ++i) {
    _item[i].AddPiecewiseLinear(_append_differences, piecewise_linear);
  }

  return 1;
}

int
Evidence::WeightedAverage() {
  for (int ndx = 0; uint32_t nbrs : _config.weighted_average()) {
    AddToHeader("WtAve", ndx);
    ++ndx;
    Evidence::WeightedAverage(nbrs);
  }

  return 1;
}
int
Evidence::WeightedAverage(int nbrs) {
  for (int i = 0; i < _number_items; ++i) {
    _item[i].AddWeightedAverage(_append_differences, nbrs);
  }

  return 1;
}

int
Evidence::WriteResults(IWString_and_File_Descriptor& output) const {
  if (_id_to_activity.size() > 0) {
    output << "Smiles" << _output_separator;
  }

  output << "ID" << _output_separator << "Activity";
  for (const IWString* h : _header) {
    output << _output_separator << *h;
  }
  output << '\n';

  for (int i = 0; i < _number_items; ++i) {
    MaybeWriteSmiles(_item[i].id(), output);

    output << _item[i].id() << _output_separator << _item[i].activity();

    _item[i].WriteResults(_output_separator, output);
    output << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

void
Evidence::MaybeWriteSmiles(const IWString& id, IWString_and_File_Descriptor& output) const {
  if (_id_to_smiles.empty()) {
    return;
  }

  const auto iter = _id_to_smiles.find(id);
  if (iter == _id_to_smiles.end()) {
    cerr << "Evidence::MaybeWriteSmiles:no smiles for '" << id << "'\n";
    output << '*' << _output_separator;
    return;
  }

  output << iter->second << _output_separator;

  return;
}

int
Main(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-A=s-C-config=sfile-diff-smiles=sfile");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  Evidence analysis;
  if (! analysis.Initialise(cl)) {
    cerr << "Cannot initialise Evidence\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(0);
  }

  if (cl.size() > 1) {
    cerr << "Extra arguments ignored\n";
  }

  if (! analysis.ReadNbrList(cl[0])) {
    cerr << "Cannot read neighbour list '" << cl[0] << "'\n";
    return 1;
  }

  if (! analysis.Process()) {
    cerr << "Analysis failed '" << cl[0] << "'\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);
  analysis.WriteResults(output);

  return 0;
}

}  // namespace evidence

int
main(int argc, char **argv) {
  int rc = evidence::Main(argc, argv);

  return rc;
}
