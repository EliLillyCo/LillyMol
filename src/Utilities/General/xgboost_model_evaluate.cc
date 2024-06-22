// Evaluate an xgboost model

#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <optional>
#include <string>

#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "xgboost/json_io.h"
#include "xgboost/c_api.h"

#include "Utilities/General/xgboost_model.pb.h"

namespace lillymol_xgboost {

using std::cerr;

// We are using the C api, since I was not able to find
// good examples of how to use the C++ API.
// TODO:ianwatson figure this out sometime

#define safe_xgboost(call) {  \
  int err = (call); \
  if (err != 0) { \
    fprintf(stderr, "%s:%d: error in %s: %s\n", __FILE__, __LINE__, #call, XGBGetLastError());  \
    exit(1); \
  } \
}

class Options {
  private:
    uint64_t _molecules_read;

    // What to do if the model does not return a valid prediction. By default, we fail.
    int _ignore_failed_predictions;

    // If we are doing batch scoring - which is the default.
    // It will never make sense to do one-at-a-time scoring,
    // since we observe that it is very slow - presumably because
    // of the overhead of decoding the JSON input.
    // zero is interpreted as no batching
    uint32_t _batch_size;

    char _input_separator;
    char _output_separator;

  // private functions

  public:
    Options();

    int Initialise(Command_Line_v2& cl);

    int ignore_failed_predictions() const {
      return _ignore_failed_predictions;
    }

    uint32_t batch_size() const {
      return _batch_size;
    }

    void set_input_separator(char s) {
      _input_separator = s;
    }
    char input_separator() const {
      return _input_separator;
    }
    void set_output_separator(char s) {
      _output_separator = s;
    }
    char output_separator() const {
      return _output_separator;
    }

    int WriteId(const const_IWSubstring& buffer, IWString_and_File_Descriptor& output) const;

    int WriteResult(float value, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _molecules_read = 0;
  _ignore_failed_predictions = 0;

  _batch_size = 1000;

  _input_separator = ' ';
  _output_separator = ' ';
}

int
Options::Initialise(Command_Line_v2& cl) {
  const int verbose = cl.option_count('v');

  if (cl.option_present("batch")) {
    if (! cl.value("batch", _batch_size) || _batch_size < 1) {
      cerr << "Options::Initialise:invalid batch size (-b)\n";
      return 0;
    }

    if (verbose) {
      cerr << "Batch size " << _batch_size << '\n';
    }
  }

  return 1;
}

int
Options::WriteId(const const_IWSubstring& buffer,
        IWString_and_File_Descriptor& output) const {
  const_IWSubstring token;
  int i = 0;
  buffer.nextword_single_delimiter(token, i, _input_separator);
  output << token;

  return 1;
}

int
Options::WriteResult(float value,
                     IWString_and_File_Descriptor& output) {
  output << _output_separator;
  output << value;
  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

void
Usage() {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(
Evaluate an xgboost descriptor model that has been built with xgbd_make.
 -mdir <directory>          model directory created with xgbd_make (do this)
                            undocumented -XJSON and -META options are also available - see the source code.
 -batch <n>                 do batched evaluation - default is 1000 items per batch.
 -precision <n>             number of significant digits in the output
 -report_unused             report unused features in the input
 -i <char>                  input delimiter - default space. Character names recognised.
 -o <char>                  output default - default space. Character names recognised.
 -v                         verbose output
)";

  // clang-format on

  exit(0);
}

class XGBDModel {
  private:
    BoosterHandle _booster;

    // The model has a column order and we may need to rearrange the
    // column order before passing to the model.

    // The number of features in the model.
    bst_ulong _ncol;

    int _columns_in_input;
    int* _column_xref;

    // We translate our input buffer to an array of floats which is then
    // sent to...
    float* _value;

    // The first call must pass a header record, and we use that to setup
    // _column_xref.
    int _first_call;

    // It can be noisy to report unused features, and quite harmless
    // TODO:ianwatson enable setting this (why?)
    int _report_unused_features;

    // Read with the -P option.
    xgboost_model::XGBoostModel _proto;

    IWString _missing_value;

    char _input_separator;

  // Private functions

  public:
    XGBDModel();
    ~XGBDModel();

    void set_input_separator(char s) {
      _input_separator = s;
    }

    const std::string& response() const {
      return _proto.response();
    }

    // This function assumes the names of the xgboost file and the
    // model metadata file.
    int Build(const IWString& mdir);

    // Alternatively an external caller can specifically initialise the 
    // two parts of the model.
    int ReadModelProto(IWString& fname);
    int ReadModel(IWString& fname);

    uint32_t ncol() const {
      return _ncol;
    }
  
    // Given a descriptor file record, extract the id and floating point
    // values, in the correct column order according to _column_xref.
    template <typename T> int GetValues(const const_IWSubstring& buffer,
                     T& id,
                     float* values) const;

    // Must be called before Score.
    int EstablishColumnXref(const const_IWSubstring& header);

#ifdef NO_LONGER_UASED____
    int StringToVector(const const_IWSubstring& buffer,
               int j, 
               float* values);
#endif

    std::optional<double> Score(const const_IWSubstring& buffer);
    int ScoreBatch(const float* values, uint32_t items_in_chunk, float* pred);
};

XGBDModel::XGBDModel() {
  _columns_in_input = 0;
  _column_xref = nullptr;

  _value = nullptr;

  safe_xgboost(XGBoosterCreate(NULL, 0, &_booster));

  _first_call = true;

  _report_unused_features = 0;

  _missing_value = ".";

  _input_separator = ' ';
}

XGBDModel::~XGBDModel() {
  XGBoosterFree(_booster);

  if (_column_xref) {
    delete [] _column_xref;
  }
  if (_value) {
    delete [] _value;
  }
}

int
XGBDModel::ReadModelProto(IWString& fname) {
  std::optional<xgboost_model::XGBoostModel> proto =
    iwmisc::ReadTextProto<xgboost_model::XGBoostModel>(fname);
  if (! proto) {
    cerr << "XGBDModel::ReadModelProto:cannot read '" << fname << "'\n";
    return 0;
  }

  _proto = std::move(*proto);

  return 1;
}

int
XGBDModel::ReadModel(IWString& fname) {
  safe_xgboost(XGBoosterLoadModel(_booster, fname.null_terminated_chars()));

  safe_xgboost(XGBoosterGetNumFeature(_booster, &_ncol));

  _value = new float[_ncol];

  return 1;
}

int
XGBDModel::Build(const IWString& mdir) {
  IWString fname;
  fname << mdir << '/' << "model_metadata.txt";
  if (! ReadModelProto(fname)) {
    cerr << "Cannot read model proto '" << fname << "'\n";
    return 0;
  }

  fname = mdir;
  fname << '/' << "xgboost.json";
  if (! ReadModel(fname)) {
    cerr << "Cannot read model file '" << fname << "'\n";
    return 0;
  }

  return 1;
}

#ifdef XGBOOST_CPP_NOT_WORKING_YET
  iwstring_data_source input(fname);
  if (! input.good())  {
    cerr << "XGBDModel::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  return Build(input);
}

int
XGBDModel::Build(iwstring_data_source& input) {
  IWString file_contents;
  file_contents.resize(input.file_size());
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    file_contents << buffer;
  }

  xgboost::StringView tmp(file_contents.data(), file_contents.size());

  xgboost::JsonReader reader(tmp);

  BoosterHandle x;
  BoosterCreate(0,0,&x);

  const std::string model_filename(input.file_name().data(), input.file_name().size());
  int y = XGBoosterLoadModel(x,model_filename);


  return 1;
}
#endif

int
XGBDModel::EstablishColumnXref(const const_IWSubstring& header) {
  _columns_in_input = header.nwords_single_delimiter(_input_separator);
  if (_columns_in_input < 2) {
    cerr << "XGBDModel::EstablishColumnXref:invalid header '" << header << "'\n";
    return 0;
  }

  _column_xref = new int[_columns_in_input];
  std::fill_n(_column_xref, _columns_in_input, -1);

  IW_STL_Hash_Map_int name_to_col;

  for (const auto&[name, col] : _proto.name_to_col()) {
    IWString tmp(name);
    name_to_col[tmp] = col;
  }

  if (name_to_col.empty()) {
    cerr << "XGBDModel::EstablishColumnXref:no column data in proto\n";
    return 0;
  }

  int i = 0;
  IWString token;

  // keep track of how many features are found. Must be the same as _ncol.
  // Note that we do not handle duplicate column names - too unusual to worry about.
  uint32_t features_found = 0;

  for (int col = 0; header.nextword_single_delimiter(token, i, _input_separator); ++col) {
    if (col == 0) {
      continue;
    }

    auto iter = name_to_col.find(token);
    // Skip features not required for the model.
    if (iter == name_to_col.end()) {
      if (_report_unused_features) {
        cerr << "Feature '" << token << "' not part of model\n";
      }
      continue;
    }

    _column_xref[col] = iter->second;
    ++features_found;
  }

  if (features_found != _ncol) {
    cerr << "XGBDModel::EstablishColumnXref:only found " << features_found << " of " << 
            _ncol << " columns\n";
    cerr << header << '\n';
    return 0;
  }

  return 1;
}

// Convert `buffer` to `id` and an array of floats.
template <typename T>
int
XGBDModel::GetValues(const const_IWSubstring& buffer,
                     T& id,
                     float* values) const {
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; buffer.nextword_single_delimiter(token, i, _input_separator); ++col) {
    if (col == 0) {
      id = token;
      continue;
    }

    if (_column_xref[col] < 0) {
      continue;
    }

    float v;
    if (token.numeric_value(v)) {
      values[_column_xref[col]] = v;
    } else if (token == _missing_value) {
      values[_column_xref[col]] = std::numeric_limits<float>::quiet_NaN();
    } else {
      cerr << "XGBDModel::Score:invalid value '" << token << "'\n";
      return 0;
    }
  }
  
  return 1;
}

std::optional<double>
XGBDModel::Score(const const_IWSubstring& buffer) {
  const_IWSubstring id;

  if (! GetValues(buffer, id, _value)) {
    return std::nullopt;
  }

  static constexpr bst_ulong kNrow = 1;

  DMatrixHandle dmatrix;
  XGDMatrixCreateFromMat(_value, kNrow, _ncol, std::numeric_limits<float>::quiet_NaN(), &dmatrix);

  static char const config[] =
      "{\"training\": false, \"type\": 0, "
        "\"iteration_begin\": 0, \"iteration_end\": 0, \"strict_shape\": false}";

  /* Shape of output prediction */
  uint64_t const* out_shape;
  /* Dimension of output prediction */
  uint64_t out_dim;
  /* Pointer to a thread local contiguous array, assigned in prediction function. */
  float const* out_result = NULL;
  safe_xgboost(
      XGBoosterPredictFromDMatrix(_booster, dmatrix, config, &out_shape, &out_dim, &out_result));

#ifdef IF_DOING_MATRIX
  for (unsigned int i = 0; i < out_dim; i++){
    printf("prediction[%i] = %f \n", i, out_result[i]);
  }
#endif

  XGDMatrixFree(dmatrix);

//free(out_result);
//delete [] out_result;

  return out_result[0];
}

int
XGBDModel::ScoreBatch(const float* values, uint32_t items_in_chunk, float* pred) {
  DMatrixHandle dmatrix;
  // cerr << "Creating dmatrix from " << items_in_chunk << " rows and " << _ncol << " columns\n";
  safe_xgboost(
    XGDMatrixCreateFromMat(values, items_in_chunk, _ncol, std::numeric_limits<float>::quiet_NaN(), &dmatrix)
  );

  static char const config[] =
      "{\"training\": false, \"type\": 0, "
        "\"iteration_begin\": 0, \"iteration_end\": 0, \"strict_shape\": false}";

  /* Shape of output prediction */
  uint64_t const* out_shape;
  /* Dimension of output prediction */
  uint64_t out_dim;
  /* Pointer to a thread local contiguous array, assigned in prediction function. */
  float const* out_result = NULL;
  safe_xgboost(
      XGBoosterPredictFromDMatrix(_booster, dmatrix, config, &out_shape, &out_dim, &out_result));

  std::copy_n(out_result, items_in_chunk, pred);

  XGDMatrixFree(dmatrix);

  return 1;
}

#ifdef NO_LONGER_UASED____
int
XGBDModel::StringToVector(const const_IWSubstring& buffer,
               int j, 
               float* values) {
  const_IWSubstring token;
  for (int col = 0; buffer.nextword_single_delimiter(token, j, _input_separator); ++col) {
    float v;
    if (! token.numeric_value(v)) {
      cerr << "XGBDModel::Score:invalid value '" << token << "'\n";
      return 0;
    }
    // cerr << "colum " << col << " value " << v << " xref " << _column_xref[col] << '\n';
    values[_column_xref[col]] = v;
  }

  return 1;
}
#endif

struct TempArrays {
  public:
    float* values;
    float* pred;
    IWString* id;

  public:
    TempArrays();
    ~TempArrays();

    int Initialise(uint32_t ncols, uint32_t batch_size);
};

TempArrays::TempArrays() {
  values = nullptr;
  pred = nullptr;
  id = nullptr;
}

TempArrays::~TempArrays() {
  delete [] values;
  delete [] pred;
  delete [] id;
}

int
TempArrays::Initialise(uint32_t ncols, uint32_t batch_size) {
  values = new float[ncols * batch_size];
  pred = new float[batch_size];
  id = new IWString[batch_size];

  return 1;
}

int
ScoreAndWrite(XGBDModel& model,
              Options& options,
              TempArrays& arrays,
              int items_in_chunk,
              IWString_and_File_Descriptor& output) {
  if (! model.ScoreBatch(arrays.values, items_in_chunk, arrays.pred)) {
    cerr << "ScoreAndWrite:cannot score\n";
    return 0;
  }

  for (int i = 0; i < items_in_chunk; ++i) {
    output << arrays.id[i];
    options.WriteResult(arrays.pred[i], output);

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
XGBDModelEvaluateBatched(iwstring_data_source& input,
                  Options& options,
                  XGBDModel& model,
                  IWString_and_File_Descriptor& output) {
  assert(options.batch_size() > 0);

  static TempArrays arrays;

  static int first_call = 1;
  if (first_call) {
    arrays.Initialise(model.ncol(), options.batch_size());
    first_call = 0;
  }

  const_IWSubstring buffer;
  uint32_t items_in_chunk = 0;
  while (input.next_record(buffer)) {
    if (! model.GetValues(buffer, arrays.id[items_in_chunk], arrays.values + model.ncol() * items_in_chunk)) {
      cerr << "XGBDModelEvaluateBatched:bad data\n";
      cerr << buffer << '\n';
      return 0;
    }

    ++items_in_chunk;
    if (items_in_chunk == options.batch_size()) {
      if (! ScoreAndWrite(model, options, arrays, items_in_chunk, output)) {
        return 0;
      }
      items_in_chunk = 0;
    }
  }

  if (items_in_chunk) {
    return ScoreAndWrite(model, options, arrays, items_in_chunk, output);
  }

  return 1;
}

int
XGBDModelEvaluate(iwstring_data_source& input,
                  Options& options,
                  XGBDModel& model,
                  IWString_and_File_Descriptor& output) {
  static int first_call = 1;

  // On the first call, we store the header record in this variable.
  // All subsequent files must have the same header.
  static IWString header;

  const_IWSubstring buffer;

  if (first_call) {
    if (! input.next_record(buffer)) {
      cerr << "XGBDModelEvaluate:cannot read header\n";
      return 0;
    }

    if (! model.EstablishColumnXref(buffer)) {
      cerr << "XGBDModelEvaluate::cannot establish column cross reference\n";
      return 0;
    }

    header = buffer;
    output << "Name" << options.output_separator() << model.response() << '\n';
  } else if (buffer != header) {
    cerr << "XGBDModelEvaluate:header mismatch\n";
    cerr << header << '\n';
    cerr << buffer << '\n';
    return 0;
  }

  if (options.batch_size()) {
    return XGBDModelEvaluateBatched(input, options, model, output);
  }

  while (input.next_record(buffer)) {
    std::optional<float> v = model.Score(buffer);
    if (v) {
    } else if (options.ignore_failed_predictions()) {
      continue;
    } else {
      cerr << "XGBDModelEvaluate::no score '" << buffer << "'\n";
      return 0;
    }

    options.WriteId(buffer, output);
    options.WriteResult(*v, output);
  }

  return 1;
}

int
XGBDModelEvaluate(const char* fname,
                  Options& options,
                  XGBDModel& model,
                  IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "XGBDModelEvaluate:cannot open '" << fname << "'\n";
    return 0;
  }

  return XGBDModelEvaluate(input, options, model, output);
}

int
Main(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-mdir=dir-XJSON=sfile-META=sfile-i=s-o=s-p=ipos-batch=ipos-report_unused");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage();
  }

  const int verbose = cl.option_count('v');

  if (cl.option_present("mdir")) {
  } else if (cl.option_present("XJSON") &&
             cl.option_present("META")) {
  } else {
    cerr << "Must specity the model via either -mdir or both -XJSON -META\n";
    Usage();
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage();
  }

  XGBDModel model;

  if (cl.option_present("mdir")) {
    IWString mdir = cl.option_value("mdir");
    if (! model.Build(mdir)) {
      cerr << "Cannot initialise model -mdir '" << mdir << "'\n";
      return 1;
    }
  } else {
    IWString fname = cl.string_value("META");
    if (! model.ReadModelProto(fname)) {
      cerr << "Cannot read model proto '" << fname << "'\n";
      return 1;
    }

    fname = cl.string_value("XJSON");
    if (! model.ReadModel(fname)) {
      cerr << "Cannot read model file '" << fname << "'\n";
      return 1;
    }
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise local options\n";
    return 1;
  }

  if (cl.option_present('i')) {
    IWString s = cl.string_value('i');
    if (! char_name_to_char(s)) {
      cerr << "Invalid input delimiter specification '" << s << "'\n";
      return 1;
    }

    options.set_input_separator(s[0]);
    model.set_input_separator(s[0]);
  }

  if (cl.option_present('o')) {
    IWString s = cl.string_value('o');
    if (! char_name_to_char(s)) {
      cerr << "Invalid output delimiter specification '" << s << "'\n";
      return 1;
    }

    options.set_output_separator(s[0]);
  }

  if (cl.option_present("precision")) {
    int p;
    if (! cl.value("precision", p) || p < 2) {
      cerr << "The output float precision mus be a whole +ve number > 2\n";
      return 1;
    }

    set_default_iwstring_float_concatenation_precision(p);
    if (verbose) {
      cerr << "Default float precision " << p << '\n';
    }
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (! XGBDModelEvaluate(fname, options, model, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  return 0;
}

}  // namespace lillymol_xgboost

int
main(int argc, char **argv) {
  int rc = lillymol_xgboost::Main(argc, argv);

  return rc;
}
