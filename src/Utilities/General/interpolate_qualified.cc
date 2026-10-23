// Use a model to get values for qualified activity values

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <optional>
#include <random>
#include <ranges>
#include <string>
#include <vector>

#define USE_OMP
#ifdef USE_OMP
#include "omp.h"
#endif

#include "google/protobuf/util/json_util.h"

#include "absl/container/flat_hash_map.h"

#include "xgboost/c_api.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"

#ifdef BUILD_BAZEL
#include "Utilities/General/interpolate_qualified.pb.h"
#else
#include "interpolate_qualified.pb.h"
#endif

namespace model_qualified_values {

using std::cerr;

#define safe_xgboost(call) {  \
  int err = (call); \
  if (err != 0) { \
    fprintf(stderr, "%s:%d: error in %s: %s\n", __FILE__, __LINE__, #call, XGBGetLastError());  \
    exit(1); \
  } \
}

class QualifiedData {
  private:
    // The initial value read from the activity file.
    float _initial;

    // The qualifier, -1 means <, 0 means = and +1 means >
    int _qualifier;

    // The number of etimates we have that satisfy the qualifier.
    uint32_t _successes;

    // We keep track of all imputed scores.
    std::vector<float> _score;

  public:
    QualifiedData(float a, int q);

    int DebugPrint(std::ostream& output) const;

    float initial_activity() const {
      return _initial;
    }

    void AnotherEstimate(float s) {
      _score.push_back(s);
      if (_qualifier < 0) {
        if (s < _initial) {
          ++_successes;
        }
      } else {
        if (s > _initial) {
          ++_successes;
        }
      }
    }

    std::optional<float> BestEstimate() const;

    int successes() const {
      return _successes;
    }

    // Return true of an estimate of `s` is a successful extrapolation
    // beyond our qualified value.
    bool SuccessfulEstimate(float s)  const;

};

QualifiedData::QualifiedData(float a, int q) {
  _initial = a;
  _qualifier = q;
  _successes = 0;
}

int
QualifiedData::DebugPrint(std::ostream& output) const {
  output << "QualifiedData, initial " << _initial << " with " << _score.size() << " predicted values";
  for (float v : _score) {
    output << ' ' << v;
  }
  output << '\n';

  return output.good();
}

std::optional<float>
QualifiedData::BestEstimate() const {
  if (_score.empty()) {
    return std::nullopt;
  }

  if (_qualifier < 0) {
    auto iter = std::min_element(_score.cbegin(), _score.cend());
    return *iter;
  } else if (_qualifier > 0) {
    auto iter = std::max_element(_score.cbegin(), _score.cend());
    return *iter;
  } else {
    return std::nullopt;  // not sure this can happen
  }
}

bool
QualifiedData::SuccessfulEstimate(float s) const {
  if (_qualifier < 0) {
    return s < _initial;
  }

  if (_qualifier > 0) {
    return s > _initial;
  }

  // Should not come here.
  return false;
}

class Measurement {
  private:
    IWString _id;

    // If this point is either measured, or has been imputed, this will be set.
    std::optional<float> _y;

    // If this measurement was initially qualified, this will be set.
    std::unique_ptr<QualifiedData> _qualified;

    // We do not own this, it will be transferred from the part of the code
    // that reads the descriptor file.
    const float* _x;

  public:
    Measurement(const IWString& s, float a, int q, const float* x);

    const IWString& id() const {
      return _id;
    }

    bool HasY() const {
      return _y.has_value();
    }

    std::optional<float> y() const {
      return _y;
    }

    bool HasQualified() const {
      return _qualified != nullptr;
    }

    const std::unique_ptr<QualifiedData>& qualified() const {
      return _qualified;
    }

    int ConvertToUnqualified();

    void CopyX(float* destination, uint32_t ncols) const {
      std::copy_n(_x, ncols, destination);
    }

    float InitialActivity() const {
      return _qualified->initial_activity();
    }

    void AnotherEstimate(float s) {
      _qualified->AnotherEstimate(s);
    }

    // If we are unqualified, return the activity.
    // If we are qualified, get the best estimate.
    float BestEstimate() const;
};

Measurement::Measurement(const IWString& s, float a, int q, const float* x) {
  _id = s;
  _x = x;

  if (q == 0) {
    _y = a;
    return;
  }

  _qualified = std::make_unique<QualifiedData>(a, q);

  return;
}

float
Measurement::BestEstimate() const {
  if (_qualified == nullptr) {
    return *_y;
  }

  if (std::optional<float> q = _qualified->BestEstimate(); q) {
    return *q;
  }

  return _qualified->initial_activity();
}

int
Measurement::ConvertToUnqualified() {
  _y = _qualified->initial_activity();;
  _qualified.release();

  return 1;
}

class Data {
  private:
    int _verbose;

    // As we read the descriptor file, fill this hash.
    absl::flat_hash_map<IWString, float*> _id_to_descriptors;

    uint32_t _ncols;
    uint32_t _nrows;

    // As we read the activity data, we place it into one of these
    // arrays. One contains unqualified data and the other those
    // measurements that are qualified.
    resizable_array_p<Measurement> _activity;
    resizable_array_p<Measurement> _qualified_activity;

    uint32_t _initial_qualified_count;

    // This will be passed to the trainer so it must be dimensioned to
    // _ncols*_nrows. 
    float* _x;
    // This will be passed to the learner as the response.
    float* _y;

    // The features for the test set.
    float* _xtest;

    Accumulator<double> _acc_activity;

    // When predictions are made, we can optionally truncate those values to
    // the observed range of values. Probably a good idea.
    int _truncate_predicted_values_to_experimental_range;

    // As we load the _x and _y arrays, this is the number of values for
    // in the training set.
    int _training_set_size;
    // The number of items in the test set. These will be loaded into the
    // _x and _y arrays just after the first _training_set_size rows.
    int _test_set_size;

    // The name of the second column in the activity file.
    IWString _response;

    // The number of items used in each test set.
    uint32_t _batch_size;

    char** _feature_names;

    interpolate_qualified::InterpolateQualifiedConfig _config;

    std::string _config_json;

    // By default, we just cycle through the qualified values in batches.
    // But instead, we can choose the batches randomly.
    int _choose_batches_randomly;
    // When selecting random subsets, keep track of which items alread selected.
    int* _selected;

    int _num_iterations;

    int _ignore_no_descriptors;

    IWString _descriptor_file_header;

    IWString _stem;

    // When we write the new activity file, we can write an extra column with the
    // initial data.
    int _add_initial_expt;

    char _separator;

  // Private functions
    int ReadDescriptors(iwstring_data_source& input);
    int ReadDescriptorRecord(const const_IWSubstring& buffer);
    int ReadActivity(iwstring_data_source& input);
    int ReadActivityRecord(const const_IWSubstring& buffer);
    int StoreFeatureNames(const const_IWSubstring& buffer);

    DMatrixHandle CreateDmatrixHandle(float* xstart, uint32_t nrows, float* y);
    BoosterHandle CreateBooster(DMatrixHandle dmats[]);

    int UnconsiderAtExtrema();
    float MaybeTruncateToRange(float v) const;

    int GetNextBatch();
    int GetNextRandomBatch();
    int GeneratePredictions(DMatrixHandle& training_data, DMatrixHandle& test_data);

    int UpdatePredictions(uint64_t const* out_shape, uint64_t out_dim, const float* out_result);

    int CountSuccessfulExtrapolations() const;

    int BuildAndScore();

    int WriteHeader(IWString_and_File_Descriptor& output) const;
    int WriteTestSet();
    int WriteTrainingSet();
    int WriteExtrapolatedResults() const;

  public:
    Data();
    ~Data();

    int Initialise(Command_Line_v2& cl);

    int ReadDescriptors(IWString& fname);

    int ReadActivity(IWString& fname);

    int Optimise();
};

Data::Data() {
  _separator = ' ';
  _nrows = 0;
  _ncols = 0;
  _training_set_size = 0;
  _batch_size = 10;
  _feature_names = nullptr;
  _ignore_no_descriptors = 0;
  _x = nullptr;
  _y = nullptr;
  _xtest = nullptr;
  _truncate_predicted_values_to_experimental_range = 0;
  _num_iterations = 500;
  _stem = "/tmp/intplq";
  _add_initial_expt = 0;
  _choose_batches_randomly = 0;
  _selected = nullptr;
}

Data::~Data() {
  if (_feature_names != nullptr) {
    for (uint32_t i = 0; i < _ncols; ++i) {
      // cerr << "Deleting column " << i << " '" << _feature_names[i] << "'\n";
      delete [] _feature_names[i];
    }
    delete [] _feature_names;
  }

  for (auto [k, v] : _id_to_descriptors) {
    delete [] v;
  }

  delete [] _x;
  delete [] _y;
  delete [] _xtest;
  if (_selected != nullptr) {
    delete [] _selected;
  }
}

int
Data::Initialise(Command_Line_v2& cl) {
  _verbose = cl.option_present('v');

  if (cl.option_present('h')) {
    int nthreads;
    if (! cl.value('h', nthreads)) {
      cerr << "Invalid -h (number threads) value\n";
      return 0;
    }
    cerr << "set  LGBM_SetMaxThreads \n";
//  if (LGBM_SetMaxThreads(nthreads) < 0) {
//    cerr << "Data::Initialise:cannot set LGBM_SetMaxThreads " << nthreads << '\n';
//    return 0;
//  }
    if (_verbose) {
      cerr << "Set number of threads to " << nthreads << '\n';
    }
  }

  if (cl.option_present("config")) {
    IWString fname = cl.string_value("config");
    std::optional<interpolate_qualified::InterpolateQualifiedConfig> tmp =
                iwmisc::ReadTextProto<interpolate_qualified::InterpolateQualifiedConfig>(fname) ;
    if (! tmp) {
      cerr << "Data::Initialise:cannot read config '" << fname << "'\n";
      return 0;
    }
    _config = std::move(*tmp);

    google::protobuf::util::JsonPrintOptions options;
    auto status = google::protobuf::util::MessageToJsonString(_config, &_config_json, options);

    if (!status.ok()) {
      cerr << "Data::Initialise:cannot proto convert to JSON\n";
      return 0;
    }
  }

  if (cl.option_present('X')) {
    IWString fname = cl.string_value('X');
    if (! ReadDescriptors(fname)) {
      cerr << "Cannot read descriptors from '" << fname << "'\n";
      return 0;
    }
  }

  if (_verbose) {
    cerr << "Read " << _nrows << " rows of X values each with " << _ncols << " columns\n";
  }

  _activity.reserve(_nrows);
  _qualified_activity.reserve(_nrows);

  if (cl.option_present('b')) {
    _ignore_no_descriptors = 1;
    if (_verbose) {
      cerr << "Will ignore activity values with no descriptors\n";
    }
  }

  if (cl.option_present('Y')) {
    IWString fname = cl.string_value('Y');
    if (! ReadActivity(fname)) {
      cerr << "Data::Initialise:cannot read activity from '" << fname << "'\n";
      return 0;
    }
  }

  if (_verbose) {
    cerr << "Read " << _qualified_activity.size() << " items with qualified data\n";
  }

  _initial_qualified_count = _qualified_activity.size();

  if (_activity.empty() || _qualified_activity.empty()) {
    cerr << "Either no qualified data or all data qualified\n";
    return 0;
  }

  _x = new float[_nrows * _ncols];
  cerr << "Y array sized to " << _nrows << " rows\n";
  _y = new float[_nrows];
  _xtest = new float[_initial_qualified_count * _ncols];

  if (_nrows != _id_to_descriptors.size()) {
    cerr << "Data::Initialise:size mismatch, _nrows " << _nrows << " but " <<
             _id_to_descriptors.size() << " descriptors\n";
    return 0;
  }

  if (_verbose) {
    cerr << "Read data on " << _nrows << " molecules\n";
  }

  for (int i = 0; i < _activity.number_elements(); ++i) {
    std::optional<float> s = _activity[i]->y();
    _y[i] = *s;
  }

  for (uint32_t i = 0; i < _activity.size(); ++i) {
    _activity[i]->CopyX(_x + i * _ncols, _ncols);
  }

  if (cl.option_present("trunc")) {
    _truncate_predicted_values_to_experimental_range = 1;
    if (_verbose) {
      cerr << "Interpolated values will be truncated to the experimental range\n";
    }
  }

  if (cl.option_present("batch")) {
    if (! cl.value("batch", _batch_size) || _batch_size == 0) {
      cerr << "Invalid batch size\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Batch size " << _batch_size << '\n';
    }
  }

  if (cl.option_present("rand")) {
    _choose_batches_randomly = 1;
    if (_verbose) {
      cerr << "Subsets will be chosen randomly\n";
    }
    _selected = new int[_qualified_activity.size()];
  }

  if (cl.option_present('S')) {
    cl.value('S', _stem);
    if (_verbose) {
      cerr << "Interpolated values written to '" << _stem << "'\n";
    }
  }

  if (cl.option_present("initial")) {
    _add_initial_expt = 1;
    if (_verbose) {
      cerr << "Output file will contain the initial activity\n";
    }
  }

  return 1;
}

int
Data::ReadDescriptors(IWString& fname) {
  iwstring_data_source input;
  if (! input.open(fname)) {
    cerr << "Data::ReadDescriptors:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadDescriptors(input);
}

int
Data::ReadDescriptors(iwstring_data_source& input) {
  const_IWSubstring buffer;

  if (! input.next_record(buffer)) {
    cerr << "Data::ReadDescriptors:cannot read header record\n";
    return 0;
  }

  _nrows = input.records_remaining();
  if (_nrows < 10) {
    cerr << "Data:;ReadDescriptors:not enough data in file\n";
    return 0;
  }
  cerr << "Find " << _nrows << " rows of data\n";

  _ncols = buffer.nwords(_separator) - 1;
  if (_ncols < 2) {
    cerr << "Data::ReadDescriptors:invalid header record\n";
    cerr << buffer << '\n';
    return 0;
  }

  cerr << "_ncols " << _ncols << '\n';

  StoreFeatureNames(buffer);
  cerr << "Begin reading\n";

  while (input.next_record(buffer)) {
    if (! ReadDescriptorRecord(buffer)) {
      cerr << "Data::ReadDescriptors:error processing line " << input.lines_read() << '\n';
      return 0;
    }
  }

  return _id_to_descriptors.size();
}

int
Data::StoreFeatureNames(const const_IWSubstring& buffer) {
  IWString token;
  int i = 0;
  buffer.nextword(token, i, _separator);

  _feature_names = new char*[_ncols];
  cerr << "Allocated " << _ncols << " columns\n";

  for (int col = 0; buffer.nextword(token, i, _separator); ++col) {
    _feature_names[col] = new char[token.length() + 1];
    ::strcpy(_feature_names[col], token.null_terminated_chars());
    cerr << "Filling column " << col << " value " << _feature_names[col] << "'\n";
  }

  return 1;
}

int
Data::ReadDescriptorRecord(const const_IWSubstring& buffer) {
  IWString id;
  int i = 0;

  std::unique_ptr<float[]> descriptors = std::make_unique<float[]>(_ncols);

  const_IWSubstring token;
  if (! buffer.nextword(id, i, _separator)) {
    cerr << "Data::ReadDescriptorRecord:cannot read id\n";
    return 0;
  }
//cerr << "id is " << id << '\n';

  for (int col = 0; buffer.nextword(token, i, _separator); ++col) {
    if (! token.numeric_value(descriptors[col])) {
      cerr << "Data::ReadDescriptorRecord:invalid numeric '" << token << "'\n";
      return 0;
    }
  }

  _id_to_descriptors[id] = descriptors.release();

  return _id_to_descriptors.size();
}

int
Data::ReadActivity(IWString& fname) {
  iwstring_data_source input;
  if (! input.open(fname)) {
    cerr << "Data:;ReadActivity:cannot open '" <<fname << "'\n";
    return 0;
  }

  return ReadActivity(input);
}

int
Data::ReadActivity(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if(! input.next_record(buffer)) {
    cerr << "Data::ReadActivity:cannot read header\n";
    return 0;
  }


  int i = 0;
  buffer.nextword(_response, i, _separator);
  buffer.nextword(_response, i, _separator);

  while (input.next_record(buffer)) {
    if (! ReadActivityRecord(buffer)) {
      cerr << "Data::ReadActivity:invalid input\n";
      cerr << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
Data::ReadActivityRecord(const const_IWSubstring& buffer) {
  IWString id;
  int i = 0;
  buffer.nextword(id, i, _separator);

  auto iter = _id_to_descriptors.find(id);
  if (iter == _id_to_descriptors.end()) {
    cerr << "Data::ReadActivityRecord:no descriptors for '" << id << "'\n";
    return _ignore_no_descriptors;
  }

  const_IWSubstring token;
  if (! buffer.nextword(token, i, _separator)) {
    cerr << "Data::ReadActivityRecord:too few tokens in input\n";
    return 0;
  }

  int qualifier = 0;
  if (token.starts_with('<')) {
    qualifier = -1;
    ++token;
  } else if (token.starts_with('=')) {
    ++token;
  } else if (token.starts_with('>')) {
    qualifier = 1;
    ++token;
  } 

  float activity;
  if (! token.numeric_value(activity)) {
    cerr << "Data:;ReadActivityRecord:invalid activity '" << token << "'\n";
    return 0;
  }

  _acc_activity.extra(activity);

  std::unique_ptr<Measurement> m = std::make_unique<Measurement>(id, activity, qualifier, iter->second);
  if (qualifier == 0) {
    _activity << m.release();
  } else {
    _qualified_activity << m.release();
  }
  
  return 1;
}

int
Data::WriteHeader(IWString_and_File_Descriptor& output) const {
  output << "id";
  for (uint32_t i = 0; i < _ncols; ++i) {
    output << ' ' << _feature_names[i];
  }
  output << '\n';

  return 1;
}

int
Data::WriteTrainingSet() {
  static int ndx = 0;

  IWString fname;
  fname << "/tmp/train" << ndx << ".x";

  IWString_and_File_Descriptor output;
  if (! output.open(fname.null_terminated_chars())) {
    cerr << "Data::WriteTrainingSet:cannot open '" << fname << "'\n";
    return 0;
  }

  WriteHeader(output);

  for (uint32_t i = 0; i < _activity.size(); ++i) {
    output << _activity[i]->id();
    for (uint32_t j = 0; j < _ncols; ++j) {
      output << ' ' << _x[i * _ncols + j];
    }
    output << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  fname.resize_keep_storage(0);
  fname << "/tmp/train" << ndx << ".y";
  ++ndx;

  IWString_and_File_Descriptor yfile;
  if (! yfile.open(fname.null_terminated_chars())) {
    cerr << "Data:;WriteTrainingSet:cannot open yfioe " << fname << "'\n";
    return 0;
  }

  yfile << "ID " << _response << '\n';
  for (uint32_t i = 0; i < _activity.size(); ++i) {
    std::optional<float> tmp = _activity[i]->y();
    yfile << _activity[i]->id() << ' ' << *tmp << '\n';
  }

  return 1;
}

int
Data::WriteTestSet() {
  static int ndx = 0;

  IWString fname;
  fname << "/tmp/test" << ndx << ".x";

  IWString_and_File_Descriptor xfile;
  if (! xfile.open(fname.null_terminated_chars())) {
    cerr << "Data::WriteTestSet:cannot open '" << fname << "'\n";
    return 0;
  }

  WriteHeader(xfile);

  for (uint32_t i = 0; i < _qualified_activity.size(); ++i) {
    xfile << _qualified_activity[i]->id();
    for (uint32_t j = 0; j < _ncols; ++j) {
      xfile << ' ' << _x[_activity.size() * _ncols + i * _ncols + j];
    }
    xfile << '\n';
    xfile.write_if_buffer_holds_more_than(4096);
  }

  fname.resize_keep_storage(0);
  fname << "/tmp/test" << ndx << ".y";
  ++ndx;
  IWString_and_File_Descriptor yfile;
  if (! yfile.open(fname.null_terminated_chars())) {
    cerr << "Data:;WriteTestSet:cannot open '" << fname << "'\n";
    return 0;
  }

  yfile << "ID " << _response << '\n';
  for (uint32_t i = 0; i < _qualified_activity.size(); ++i) {
    yfile << _qualified_activity[i]->id() << ' ' << _qualified_activity[i]->InitialActivity() << '\n';
  }

  return 1;
}

int
WriteDMatrixHandle(DMatrixHandle& dmatrix, IWString& fname) {
  IWString_and_File_Descriptor output;
  if (! output.open(fname.null_terminated_chars())) {
    cerr << "WriteDMatrixHandle:cannot open '" << fname << "'\n";
    return 0;
  }

  bst_ulong nrows;
  safe_xgboost(XGDMatrixNumRow(dmatrix, &nrows));
  bst_ulong ncols;
  safe_xgboost(XGDMatrixNumCol(dmatrix, &ncols));

  return 1;
}

DMatrixHandle
Data::CreateDmatrixHandle(float* xstart, uint32_t nrows, float* y) {
  static constexpr float kMissingValue = std::numeric_limits<float>::max();

  DMatrixHandle rc;
  safe_xgboost(XGDMatrixCreateFromMat(xstart, nrows, _ncols, kMissingValue, &rc));
  if (y != nullptr) {
    safe_xgboost(XGDMatrixSetFloatInfo(rc, "label", y, nrows));
  }

  return rc;

  // Useful for debugging.
  if (y != nullptr) {
    WriteTrainingSet();
  } else {
    WriteTestSet();
  }

  return rc;
}

BoosterHandle
Data::CreateBooster(DMatrixHandle dmats[]) {
  BoosterHandle booster;
  safe_xgboost(XGBoosterCreate(dmats, 1, &booster));
  // Never could get this to work
  // safe_xgboost(XGDMatrixSetStrFeatureInfo(booster, "feature_name", (const char**)_feature_names, _ncols));
  safe_xgboost(XGBoosterSetParam(booster, "objective", "reg:squarederror"));
  safe_xgboost(XGBoosterSetParam(booster, "eta", "0.4"));
  safe_xgboost(XGBoosterSetParam(booster, "max_depth", "5"));
  safe_xgboost(XGBoosterSetParam(booster, "n_estimators", "500"));

  return booster;
}

// qualified values that are at the min or max values do not need to be processed.
int
Data::UnconsiderAtExtrema() {
  float minval = *_activity[0]->y();
  float maxval = *_activity[0]->y();
  for (const Measurement* m : _activity) {
    const float v = *m->y();
    if (v < minval) {
      minval = v;
    } else if (v > maxval) {
      maxval = v;
    }
  }

  int rc = 0;
  for (int i = _qualified_activity.number_elements() - 1; i>= 0; --i) {
    const std::unique_ptr<QualifiedData>& q = _qualified_activity[i]->qualified();
    if (q->initial_activity() <= minval ||
        q->initial_activity() >= maxval) {
      Measurement* m = _qualified_activity[i];
      _qualified_activity.remove_no_delete(i);
      m->ConvertToUnqualified();
      _activity << m;
      ++rc;
    }
  }

  if (_verbose) {
    cerr << "Turned off " << rc << " qualified items at extrema\n";
  }

  return rc;
}


int
Data::Optimise() {
  UnconsiderAtExtrema();

  // First step is to build a model on unqualified and predict qualified.

  for (int i = 0; i < _qualified_activity.number_elements(); ++i) {
    _qualified_activity[i]->CopyX(_xtest + i * _ncols, _ncols);
  }

  BuildAndScore();

  if (_verbose) {
    int successful = CountSuccessfulExtrapolations();
    cerr << "After initial model, " << successful << " successful extrapolations\n";
  }

  uint32_t nsteps = _qualified_activity.size() / _batch_size;
  if (nsteps * _batch_size < _qualified_activity.size()) {
    ++nsteps;
  }

  for (uint32_t i = 0; i < nsteps; ++i) {
    if (_verbose) {
      cerr << "Begin step " << i << " _activity contains " << _activity.size() << '\n';
    }

    int bsize = GetNextBatch();
    if (bsize == 0) {
      break;
    }

    BuildAndScore();
  }

  if (_verbose) {
    int successful = CountSuccessfulExtrapolations();
    cerr << "After batched models, " << successful << " successful extrapolations\n";
  }

  int improved = 0;
  int qualified = 0;
  for (uint32_t i = 0; i < _activity.size(); ++i) {
    const Measurement* m = _activity[i];
    if (! m->HasQualified()) {
      continue;
    }
    ++qualified;
    const std::unique_ptr<QualifiedData>& q = m->qualified();
    cerr << " i " << i << ' ';
    q->DebugPrint(cerr);
    float b = *q->BestEstimate();
    if (q->SuccessfulEstimate(b)) {
      ++improved;
    }
    cerr << "BestEstimate " << *q->BestEstimate() << '\n';
  }

  if (_verbose) {
    cerr << "Improved " << improved << " of " << qualified << " qualified values\n";
  }

  // The last items in the _activity array are the qualified.
  uint32_t istart = _activity.size() - _initial_qualified_count;

  for (uint32_t i = _activity.size() - 1; i >= istart; --i) {
    const std::unique_ptr<QualifiedData>& q = _activity[i]->qualified();
    if (! q) [[unlikely]] {
      cerr << "Data::Optimise:no qualified value for " << i << '\n';
      return 0;
    }

    // If this has successfully extrapolated, keep it in train.
    if (q->successes()) {
      continue;
    }
    // No successes, move back to qualified.
    _qualified_activity << _activity[i];
    _activity.remove_no_delete(i);
  }

  if (_qualified_activity.empty()) {
    cerr << "Data::Optimise:all values successfully extrapolated\n";
    return WriteExtrapolatedResults();
  }

  // to be safe, just copy everything

  for (uint32_t i = 0; i < _activity.size(); ++i) {
    _activity[i]->CopyX(_x + i * _ncols, _ncols);
    if (_activity[i]->HasQualified()) {
      const std::unique_ptr<QualifiedData>& q = _activity[i]->qualified();
      _y[i] = *q->BestEstimate();
    } else {
      std::optional<float> a = _activity[i]->y();
      _y[i] = *a;
    }
  }

  for (uint32_t i = 0; i < _qualified_activity.size(); ++i) {
    _qualified_activity[i]->CopyX(_xtest + i * _ncols, _ncols);
  }

  BuildAndScore();

  improved = 0;
  for (const Measurement* m : _qualified_activity) {
    const std::unique_ptr<QualifiedData>& q = m->qualified();
    float b = *q->BestEstimate();
    if (q->SuccessfulEstimate(b)) {
      ++improved;
    }
  }

 if (_verbose) {
  cerr << "Final round finds " << improved << " of " << _qualified_activity.size() << " now extrapolated\n";
 }

 if (_verbose) {
   uint32_t successes = CountSuccessfulExtrapolations();
   cerr << "Of " << _initial_qualified_count << " starting qualified measurements, extrapolated " <<
           successes << " fractional " << iwmisc::Fraction<float>(successes, _initial_qualified_count)
           << '\n';
 }

  return WriteExtrapolatedResults();
}

int
Data::BuildAndScore() {
  DMatrixHandle training_data = CreateDmatrixHandle(_x, _activity.size(), _y);

  BoosterHandle booster = CreateBooster(&training_data);
  for (int i = 0; i < _num_iterations; ++i) {
    safe_xgboost(XGBoosterUpdateOneIter(booster, i, training_data));
  }

  DMatrixHandle dtest = CreateDmatrixHandle(_xtest, _qualified_activity.size(), nullptr);

  /* Shape of output prediction */
  uint64_t const* out_shape;
  /* Dimension of output prediction */
  uint64_t out_dim;
  /* Pointer to a thread local contigious array, assigned in prediction function. */
  float const* out_result = NULL;

  static const char* config = R"({
  "type": 0,
  "training": false,
  "iteration_begin": 0,
  "iteration_end": 100,
  "strict_shape": false
})";

  safe_xgboost(XGBoosterPredictFromDMatrix(booster, dtest, config, &out_shape, &out_dim, &out_result));

  UpdatePredictions(out_shape, out_dim, out_result);

  XGDMatrixFree(training_data);
  XGDMatrixFree(dtest);
  XGBoosterFree(booster);

  return 1;
}

int
Data::CountSuccessfulExtrapolations() const {
  int rc = 0;

  for (const Measurement* m : _qualified_activity) {
    const std::unique_ptr<QualifiedData>& q = m->qualified();
    if (q->successes()) {
      ++rc;
    }
  }

  for (const Measurement* m : _activity | std::views::reverse) {
    if (! m->HasQualified()) {
      break;
    }
    const std::unique_ptr<QualifiedData>& q = m->qualified();
    if (q->successes()) {
      ++rc;
    }
  }
    
  return rc;
}

int
Data::GetNextRandomBatch() {
  uint32_t nsel = std::min(_batch_size, _qualified_activity.size());

  std::fill_n(_selected, _qualified_activity.size(), 0);

  static std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint32_t> u(0, _qualified_activity.size() - 1);

  uint32_t count = 0;
  for (uint32_t i = 0; i < 2 * nsel; ++i) {
    int j = u(gen);
    if (_selected[j]) {
      continue;
    }
    _selected[j] = 1;
    ++count;
    if (count == nsel) {
      break;
    }
  }

  uint32_t x_index = _activity.size() * _ncols;
  uint32_t y_index = _activity.size();

  for (int i = _qualified_activity.number_elements() - 1; i >= 0; --i) {
    if (_selected[i] == 0) {
      continue;
    }

    Measurement* m = _qualified_activity[i];
    _qualified_activity.remove_no_delete(i);

    m->CopyX(_x + x_index, _ncols);
    x_index += _ncols;
    _y[y_index] = m->BestEstimate();
    ++y_index;
    _activity << m;
  }

  x_index = 0;
  for (uint32_t i = 0; i < _qualified_activity.size(); ++i) {
    _qualified_activity[i]->CopyX(_xtest + x_index, _ncols);
    x_index += _ncols;
  }

  return 1;
}

int
Data::WriteExtrapolatedResults() const {
  IWString fname(_stem);
  fname << ".dat";

  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "Data::WriteExtrapolatedResults:cannot open '" << fname << "'\n";
    return 0;
  }
  cerr << "Just opened '" << fname << "'\n";

  static constexpr char kSep = ' ';

  output << "ID" << kSep << _response;
  if (_add_initial_expt) {
    output << kSep << _response << "_init";
  }
  output << '\n';

  for (const Measurement* m : _activity) {
    output << m->id() << kSep;
    if (m->HasQualified()) {
      const std::unique_ptr<QualifiedData>& q = m->qualified();
      output << *q->BestEstimate();
    } else {
      std::optional<float> y = m->y();
      output << *y;
    }

    if (_add_initial_expt) {
      output << kSep;
      if (m->HasQualified()) {
        const std::unique_ptr<QualifiedData>& q = m->qualified();
        output << q->initial_activity();
      } else {
        std::optional<float> a = m->y();
        output << *a;
      }
    }
    output << '\n';
  }

  // All of these will have qualified data.
  for (const Measurement* m : _qualified_activity) {
    output << m->id() << kSep;
    const std::unique_ptr<QualifiedData>& q = m->qualified();
    std::optional<float> b = *q->BestEstimate();
    if (q->SuccessfulEstimate(*b)) {
      output << *b;
    } else {
      output << q->initial_activity();
    }

    if (_add_initial_expt) {
      output << kSep << q->initial_activity();
    }
    output << '\n';
  }

  return 1;
}

// Remove _batch_size values from _qualified_activity to _activity.
// For each item moved into _activity, copy their X values into the _x array.
// For each item moved into _activity, copy their Y values into the _y array.
int
Data::GetNextBatch() {
  uint32_t n = std::min(_batch_size, _qualified_activity.size());
  if (n == 0) {
    return 0;
  }

  if (_choose_batches_randomly) {
    return GetNextRandomBatch();
  }

  // Start copying into the _x array here.
  uint32_t initial_x_index = _activity.size() * _ncols;
  uint32_t initial_y_index = _activity.size();

  for (uint32_t i = 0; i < n; ++i) {
    _activity << _qualified_activity.pop();
    _activity.back()->CopyX(_x + initial_x_index + i * _ncols, _ncols);
    _y[initial_y_index + i] = _activity.back()->BestEstimate();
  }

  // Copy the remaining qualified values to their proper location.
  for (uint32_t i = 0; i < _qualified_activity.size(); ++i) {
    _qualified_activity[i]->CopyX(_xtest + i * _ncols, _ncols);
  }

  return n;
}

// build a model on the data in _x
int
Data::GeneratePredictions(DMatrixHandle& training_data, DMatrixHandle& test_data) {
  return 1;
}

float
Data::MaybeTruncateToRange(float v) const {
  if (! _truncate_predicted_values_to_experimental_range) {
    return v;
  }

  if (v < _acc_activity.minval()) {
    return _acc_activity.minval();
  }

  if (v > _acc_activity.maxval()) {
    return _acc_activity.maxval();
  }

  return v;
}

int
Data::UpdatePredictions(uint64_t const* out_shape, uint64_t out_dim, const float* out_result) {
  // cerr << "UpdatePredictions out_dim " << out_dim << " shape " << *out_shape << '\n';
  Accumulator<double> acc_expt;
  Accumulator<double> acc_pred;
  Accumulator<double> diffs;

  for (uint64_t i = 0; i < *out_shape; ++i) {
    float expt = _qualified_activity[i]->InitialActivity();
    float pred = MaybeTruncateToRange(out_result[i]);
    acc_pred.extra(pred);
    // cerr << "  i " << i << ' ' << _qualified_activity[i]->id() << " expt " << expt << " pred " << pred << '\n';
    diffs.extra(std::abs(expt - pred));
    acc_expt.extra(expt);
    _qualified_activity[i]->AnotherEstimate(pred);
  }

  if (acc_pred.empty()) {
    return 1;
  }

  // This does not work for some reason, R2 is always wrong.
#ifdef DOES_NOT_WORK
  if (_verbose) {
    cerr << " ave expt " << static_cast<float>(acc_expt.average()) <<
            " ave mean " << static_cast<float>(acc_pred.average()) << '\n';
    float rms = std::sqrt(diffs.sum_of_squares() / diffs.n());
    float mean = acc_expt.average();
    float sst = 0.0f;
    for (uint64_t i = 0; i < out_shape[0]; ++i) {
      float y = _qualified_activity[i]->InitialActivity();
      float d = mean - y;
      sst += d * d;
    }
    float r2 = 1.0f - diffs.sum_of_squares() / sst;
    cerr << "RMS " << rms << " R2 " << r2 << '\n';
  }
#endif // DOES_NOT_WORK

  // Cut and paste from iwstats.
  if (_verbose) {
    const int len = out_shape[0];

    Accumulator<double> acc_expt, acc_pred, acc_diffs;
    for (int i = 0; i < len; ++i) {
      acc_expt.extra(_qualified_activity[i]->InitialActivity());
      acc_pred.extra(out_result[i]);
      acc_diffs.extra(_qualified_activity[i]->InitialActivity() - out_result[i]);
    }

    float obar = acc_expt.average();
    float pbar = acc_pred.average();
    // cerr << " obar " << obar << " pbar " << pbar << '\n';

    float n = 0.0;    // numerator
    float dno = 0.0;  // denominator for observed
    float dnp = 0.0;  // denominator for predicted
    for (int i = 0; i < len; ++i) {
      float expt = _qualified_activity[i]->InitialActivity();
      float pred = out_result[i];
      n += (expt - obar) * (pred - pbar);
      dno += (expt - obar) * (expt - obar);
      dnp += (pred - pbar) * (pred - pbar);
    }

    float rms = std::sqrt(acc_diffs.sum_of_squares() / acc_diffs.n());
    float r2 = n * n / (dno * dnp);
    // cerr << "dno " << dno << " dnp " << dnp << '\n';
    cerr << "RMS " << rms << " R2 " << r2 << '\n';
  }

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
  cerr << R"(Imputes qualified data.
Given measured data like <3.14 or >10, build models on the rest of the data and
see if predicted values for the qualified data can give better values. If a model
prediction is within the qualified range, less than 3.14 or greater than 10, that
imputed value is used.
 -X <fname>             File of descriptor values. Used for model building.
 -Y <fname>             Activity file.
 -i <sep>               Input file separator, default ' '.
 -b                     Ignore activity values (Y) with no descriptors (X).
 -batch <size>          Number of qualified values per batch.
 -rand                  Form batches randomly - rather than sequentially.
 -trunc                 Truncate any predicted value to the range of experimental values.
 -S <stem>              Write updated activity file to <stem>.activity.
 -initial               Also write the initial (qualified) values to the output file.
 -v                     Verbose output.
)";

  ::exit(rc);
}

int
Main(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-X=sfile-Y=sfile-i=s-batch=ipos-h=int-b-xgboost=sfile-S=s-initial-rand-trunc");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unknown_options_encountered\n";
    Usage(1);
  }

  if (! cl.option_present("X")) {
    cerr << "Must specify the descriptor file via the -X option\n";
    Usage(1);
  }

  if (! cl.option_present("Y")) {
    cerr << "Must specify the activity file via the -Y option\n";
    Usage(1);
  }

  if (! cl.empty()) {
    cerr << "Other command line arguments are ignored\n";
  }

  Data data;
  if (! data.Initialise(cl)) {
    cerr << "Cannot initialise\n";
    return 1;
  }
  
  data.Optimise();
 // Dummy();

  return 0;
}

}  // namespace model_qualified_values

int
main(int argc, char ** argv)
{
  int rc = model_qualified_values::Main(argc, argv);

  return rc;
}
