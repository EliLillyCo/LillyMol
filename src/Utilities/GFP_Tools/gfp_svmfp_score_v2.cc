// Score an svmfp model.
// A model is described by a set of weighted fingerprints
// and a bias term.
// The support vectors, and weights, are read, and then reduced
// to the subset used by the model.
// Model metadata is stored in protos.

#include <fstream>
#include <tuple>

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#ifdef BUILD_BAZEL
#include "Utilities/General/class_label_translation.pb.h"
#else
#include "class_label_translation.pb.h"
#endif

#define FEATURE_SCALER_IMPLEMENTATION
#include "Utilities/General/scaler.h"

#include "Utilities/GFP_Tools/bit_subset.h"
#include "Utilities/GFP_Tools/gfp.h"

#ifdef BUILD_BAZEL
#include "Utilities/GFP_Tools/gfp_model.pb.h"
#include "Utilities/GFP_Tools/gfp_to_svm_lite.pb.h"
#else
#include "gfp_model.pb.h"
#include "gfp_to_svm_lite.pb.h"
#endif

namespace gfp_svmfp_evaluate {

using std::cerr;

int verbose = 0;

// This tag is hard coded, should be command line settable.
IWString support_vector_tag = "NCSV<";
// The tag holding the support vector weight.
IWString weight_tag = "WEIGHT<";

char output_separator = ' ';

// If we have a classification problem, what do we write?
// By default, we write just the numeric score, but we can
// write other things.
bool write_label = true;
bool write_score = false;

Fraction_as_String fraction_as_string;

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
  cerr << "Evaluates svmfp model(s) built with svmfp_make\n";
  cerr << " -mdir <dir>     one or more model directories\n";
  cerr << " -cwrite ...     what to write for classification models, '-cwrite help' for info\n";
  cerr << " -sas <n>        write the score as string with <n> digits of accuracy\n";
  cerr << " -v              verbose output\n";
  // clang-format on

  exit(rc);
}

void
DisplayCWriteOptions(std::ostream& output) {
  output << "The -cwrite option controls output from classification models\n";
  output << " -cwrite label      write the label\n";
  output << " -cwrite score      write the score\n";

  exit(0);
}

// Read a binary encoded proto from `fname` into `proto`.
// TODO:ianwatson transition to utilities in iwmisc
template <typename Proto>
int
ReadBinaryProto(IWString& fname,
                Proto& proto) {
  using std::ios;
  std::fstream input(fname.null_terminated_chars(), ios::in | ios::binary);
  if (! input) {
    cerr << "ReadBinaryProto:cannot open '" << fname << "'\n";
    return 0;
  }
  if (! proto.ParseFromIstream(&input)) {
    cerr << "ReadBinaryProto::cannot decode '" << fname << "'\n";
    return 0;
  }

  return 1;
}

// Read a single gfp from train.gfp in order that the global
// gfp conditions get set up.
// This does not really work for multiple models since they
// may have different fingerprints. Need to wait till we can 
// handle multiple gfp configs at once.
int
ReadSingleGfp(const IWString& mdir, const std::string& fname) {
  IWString path;
  path << mdir << '/' << fname;
  iwstring_data_source input(path);
  if (! input.good()) {
    cerr << "ReadSingleGfp:cannot open '" << path << "'\n";
    return 0;
  }

  IW_TDT tdt;
  if (! tdt.next(input)) {
    cerr << "ReadSingleGfp:cannot read '" << path << "'\n";
    return 0;
  }

  IW_General_Fingerprint notused;
  int fatal = 0;
  if (! notused.construct_from_tdt(tdt, fatal)) {
    cerr << "ReadSingleGfp:bad gfp data in '" << path << "'\n";
    return 0;
  }

  return 1;
}

IWString
PathName(const IWString& dirname, const IWString& fname) {
  IWString result;
  result << dirname << '/' << fname;
  return result;
}

// Read a textproto of type `Proto` from `dir/fname`.
template <typename Proto>
std::optional<Proto>
ReadTextProto(const IWString& dirname, const std::string& fname) {
  IWString path_name = PathName(dirname, fname);
  return iwmisc::ReadTextProto<Proto>(path_name);
}

// Read a binary encoded proto of type `Proto` from `dirname/fname`.
template <typename Proto>
std::optional<Proto>
ReadBinaryProto(const IWString& dirname, const std::string& fname) {
  IWString path_name = PathName(dirname, fname);
  return iwmisc::ReadBinaryProto<Proto>(path_name);
}

// The support vectors are stored in gfp form with a weight.
class WeightedFingerprint {
  private:
    IW_General_Fingerprint _fp;
    double _weight;

  public:
    WeightedFingerprint();

    double weight() const { return _weight;}

    int Build(IW_TDT& tdt);

    IW_General_Fingerprint& fp() { return _fp;}

    double Tanimoto(IW_General_Fingerprint& fp);
};

WeightedFingerprint::WeightedFingerprint() {
  _weight = 0.0;
}

int
WeightedFingerprint::Build(IW_TDT& tdt) {
  int fatal = 0;
  if (! _fp.construct_from_tdt(tdt, fatal)) {
    cerr << "WeightedFingerprint::Build:cannot parse " << tdt << "'\n";
    return 0;
  }

  if (! tdt.dataitem_value(weight_tag, _weight)) {
    cerr << "WeightedFingerprint::Built:invalid weight\n";
    return 0;
  }

  return 1;
}

double
WeightedFingerprint::Tanimoto(IW_General_Fingerprint& fp) {
  return _fp.equal_weight_tanimoto(fp) * _weight;
}

// Model class populated from a svmfp_model proto.
class SvmModel {
  private:
    int _nfingerprints;
    WeightedFingerprint* _support_vector;

    // If directed by the proto, discard count information in sparse fingerprints.
    bool _flatten_counts;

    double _threshold_b;

    // In a weighted average, the weight assigned to this model.
    float _weight;

    // Read from the proto.
    IWString _response_name;

    // A quick way to check if this is a classification problem.
    bool _is_regression;

    // If we have a classification problem, the class labels to apply.
    IWString _class_labels[2];

    // If a regression model, scale to [0,1] and back.
    feature_scaler::FeatureScaler<float> _response_scaler;

    // Enable subsetting a gfp to just those features needed by the model.
    bit_subset::BitSubset _bit_subset;

    // It can be interesting to keep track of the scores computed.
    Accumulator<double> _acc_scores;

  // private functions

    int ReadSupportVectors(const IWString& dir, const IWString& fname);
    int ReadSupportVectors(iwstring_data_source& input);
    int ReadClassLabelTranslation(const IWString& dir, const std::string& fname);
    int ReadRegressionScaling(const IWString& dir, const std::string& fname);

    // Make any changes needed to `gfp` before it is used for scoring.
    // This includes subsetting the bits, and possibly flattening sparse fingerprints.
    int PreProcess(IW_General_Fingerprint& gfp);

    // Once the support vectors are read, preprocess them the same as
    // the fingerprints to be scored.
    int PreprocessSupportVectors();

  public:
    SvmModel();
    ~SvmModel();

    // Initialise an SvmModel based on a SvmfpModel proto.
    // Argument `directory` is needed so that the files
    // specified in `model_proto` can be found.
    int Initialise(const GfpModel::SvmfpModel& model_proto,
                   const IWString& directory);

    bool is_regression() const { return _is_regression;}

    const IWString& response_name() const { return _response_name;}

    // The numeric score for this model.
    // Should be const, but for updating _acc_scores.
    double Score(IW_General_Fingerprint& fp);

    // Given a numeric score, translate to the appropriate label.
    const IWString& ClassLabelForScore(double score) const;

    // Report on the scores assigned.
    int Report(std::ostream& output) const;
};

SvmModel::SvmModel() {
  _nfingerprints = 0;
  _support_vector = nullptr;
  _threshold_b = 0.0;
  _flatten_counts = false;
  _weight = 1.0;
  _is_regression = true;
}

SvmModel::~SvmModel() {
  if (_support_vector != nullptr) {
    delete [] _support_vector;
    _support_vector = nullptr;
  }

  _nfingerprints = 0;
  _weight = 0.0;
}

// Read information from `model_proto` to build state.
int
SvmModel::Initialise(const GfpModel::SvmfpModel& model_proto,
                     const IWString& dir) {
  if (! model_proto.has_bit_subset()) {
    cerr << "SvmModel::Initialise:missing bit_subset\n";
    return 0;
  }

  if (! model_proto.has_train_gfp()) {
    cerr << "SvmModel::Initialise:missing train gfp file\n";
    return 0;
  }

  if (! model_proto.has_support_vectors()) {
    cerr << "SvmModel::Initialise:missing support vectors";
    return 0;
  }

  // Establish gfp state.
  if (! ReadSingleGfp(dir, model_proto.train_gfp())) {
    cerr << "SvmModel::Initialise:cannot read " << model_proto.train_gfp() << '\n';
    return 0;
  }

  if (! model_proto.has_threshold_b()) {
    cerr << "SvmModel::Initialise:no threshold_b\n";
    return 0;
  }

  _threshold_b = model_proto.threshold_b();

  _response_name = model_proto.metadata().response_name();

  std::string fname = model_proto.bit_subset();
  std::optional<GfpBitSubset::GfpBitSubset> subset_proto = 
                ReadBinaryProto<GfpBitSubset::GfpBitSubset>(dir, fname);
  if (! subset_proto) {
    cerr << "Cannot read bit subset proto '" << fname << "'\n";
    return 0;
  }

  if (! _bit_subset.Build(*subset_proto)) {
    cerr << "SvmModel::Initialise:cannot process subset proto\n";
    return 0;
  }

  fname = model_proto.support_vectors();
  if (! ReadSupportVectors(dir, fname)) {
    cerr << "SvmModel::Initialise:cannot read support vectors '" << fname << "'\n";
    return 0;
  }

  PreprocessSupportVectors();

  _flatten_counts = model_proto.metadata().flatten_sparse_fingerprints();

  if (model_proto.metadata().has_class_label_translation()) {
    const std::string fname = model_proto.metadata().class_label_translation();
    if (!ReadClassLabelTranslation(dir, fname)) {
      cerr << "SvmModel:cannot read class_label_translation " << fname << '\n';
      return 0;
    }
    _is_regression = false;
  }

  if (model_proto.metadata().has_response_scaling()) {
    if (! _is_regression) {
      cerr << "SvmModel::Initialise:classification model cannot also have response scaling\n";
      return 0;
    }
    const std::string fname = model_proto.metadata().response_scaling();
    if (! ReadRegressionScaling(dir, fname)) {
      cerr << "SvmModel::Initialise:cannot read response scaling '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

// `fname` contains the name of the file containing the weighted support vectors.
// Read that file and populate _support_vector.
int
SvmModel::ReadSupportVectors(const IWString& dir, const IWString& fname) {
  IWString path_name = PathName(dir, fname);
  iwstring_data_source input(path_name.null_terminated_chars());
  if (! input.good()) {
    cerr << "SvmModel::ReadSupportVectors:cannot open '" << path_name << "'\n";
    return 0;
  }

  std::unique_ptr<RE2> vbar = std::make_unique<RE2>("^\\|$");
  _nfingerprints = input.grep(*vbar);
  if (_nfingerprints == 0) {
    cerr << "SvmModel::ReadSupportVectors:no fingerprints in '" << path_name << "'\n";
    return 0;
  }

  _support_vector = new WeightedFingerprint[_nfingerprints];

  if (verbose) {
    cerr << "Reading " << _nfingerprints << " support vectors\n";
  }

  return ReadSupportVectors(input);
}

// Populate the already allocated _support_vector array with the contents of `input`.
int
SvmModel::ReadSupportVectors(iwstring_data_source& input) {
  IW_TDT tdt;
  for (int i = 0; tdt.next(input); ++i) {
    if (! _support_vector[i].Build(tdt)) {
      cerr << "SvmModel::ReadSupportVectors:cannot read " << tdt << '\n';
      return 0;
    }
  }

  if (verbose > 1) {
    cerr << "Support vectors initialized\n";
  }

  return 1;
}

void FlattenSparseFingerprint(IW_General_Fingerprint& gfp) {
  const int nsparse = number_sparse_fingerprints();
  for (int i = 0; i < nsparse; ++i) {
    Sparse_Fingerprint& sfp = gfp.sparse_fingerprint(i);
    sfp.truncate_counts_at(1);
  }
}

int
SvmModel::PreProcess(IW_General_Fingerprint& gfp) {
  int rc = _bit_subset.MakeSubset(gfp);
  if (_flatten_counts) {
    FlattenSparseFingerprint(gfp);
  }

  return rc;
}

// Returns the number of bits removed.
int
SvmModel::PreprocessSupportVectors() {
  int rc = 0;

  for (int i = 0; i < _nfingerprints; ++i) {
    rc += PreProcess(_support_vector[i].fp());
  }
  return rc;
}

// The file `dir/fname` is a ClassLabelTranslation proto. Read that file
// and use it to establish the _class_labels array.
int
SvmModel::ReadClassLabelTranslation(const IWString& dir, const std::string& fname) {
  std::optional<ClassLabelTranslation::ClassLabelTranslation> mapping =
    ReadBinaryProto<ClassLabelTranslation::ClassLabelTranslation>(dir, fname);
  if (! mapping) {
    cerr << "SvmModel::ReadClassLabelTranslation:cannot read '" << fname << "'\n";
    return 0;
  }

  for (auto [key, value] : mapping->to_numeric()) {
    if (value == -1) {
      _class_labels[0] = key;
    } else if (value == 1) {
      _class_labels[1] = key;
    } else {
      cerr << "SvmModel::ReadClassLabelTranslation:unrecognised pair '" << key << "' value " << value << '\n';
      return 0;
    }
  }

  if (_class_labels[0].empty() || _class_labels[1].empty()) {
    cerr << "SvmModel::ReadClassLabelTranslation:incomplete specification " << mapping->ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}

int
SvmModel::ReadRegressionScaling(const IWString& dir, const std::string& fname) {
  std::optional<FeatureScaling::FeatureScaling> proto =
    ReadBinaryProto<FeatureScaling::FeatureScaling>(dir, fname);
  if (! proto) {
    cerr << "SvmModel::ReadRegressionScaling:cannot read '" << fname << "'\n";
    return 0;
  }

  if (! _response_scaler.Initialise(*proto)) {
    cerr << "SvmModel::ReadRegressionScaling:cannot parse " << proto->ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}

// A common operation is to fetch the number of each kind of fingerprint.
// A convenience function to enable this with one line of code rather than 2.
std::tuple<int, int>
GetNumberFingerprints() {
  return {number_fingerprints(), number_sparse_fingerprints()};
}

int
SvmModel::Report(std::ostream& output) const {
  if (_acc_scores.n() == 0) {
    output << "SvmModel::Report: " << _response_name << " no scores\n";
    return 1;
  }

  output << "SvmModel::Report: " << _response_name << " scored " << _acc_scores.n() << " fingerprints";
  output << " scores btw " << _acc_scores.minval() << " and " << _acc_scores.maxval() << " mean " << _acc_scores.average() << '\n';

  return 1;
}

// Score a single fingerprint. The weighted similarity to our support vectors.
double
SvmModel::Score(IW_General_Fingerprint& gfp) {
  PreProcess(gfp);

  double result = - _threshold_b;
  for (int i = 0; i < _nfingerprints; ++i) {
    result += _support_vector[i].Tanimoto(gfp) ;
  }

  // If needed, convert back to the original range of the model data.
  //cerr << "Before rescaling " << result << '\n';
  if (_response_scaler.Active()) {
    result = _response_scaler.ScaleBackToOrignalRange(result);
  }
  //cerr << "After rescaling " << result << '\n';

  _acc_scores.extra(result);

  return result;
}

const IWString&
SvmModel::ClassLabelForScore(double score) const {
  // Classification problem, translate to class label.
  if (score < 0.0) {
    return _class_labels[0];
  } else {
    return _class_labels[1];
  }
}

// Append a string representation of `score` to `output`.
// If fraction_as_string is active, use it.
void
AppendScore(float score,
            IWString_and_File_Descriptor& output) {

  if (fraction_as_string.active()) {
    fraction_as_string.append_number(output, score);
  } else {
    output << output_separator << score;
  }
}

// Copy from a protobuf::Map to an std::unordered_map.
// There is a typing mismatch, ignored.
void
CopyXref(const google::protobuf::Map<long unsigned int, long unsigned int>& from,
         std::unordered_map<uint32_t, uint32_t>& to) {
  for (auto [b, f] : from) {
    to.emplace(b, f);
  }
}
         
int
GfpSvmfpEvaluate(IW_General_Fingerprint& fp,
                 SvmModel* models,
                 int nmodels,
                 IWString_and_File_Descriptor& output) { 
  append_first_token_of_name(fp.id(), output);

  for (int i = 0; i < nmodels; ++i) {
    const auto score = models[i].Score(fp);
    if (models[i].is_regression()) {
      AppendScore(score, output);
    } else {
      if (write_label) {
        output << output_separator << models[i].ClassLabelForScore(score);
      }
      if (write_score) {
        AppendScore(score, output);
      }
    }
  }

  output << '\n';
  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

int
GfpSvmfpEvaluate(IW_TDT& tdt,
                 SvmModel* models,
                 int nmodels,
                 IWString_and_File_Descriptor& output) { 
  IW_General_Fingerprint fp;
  int fatal;
  if (! fp.construct_from_tdt(tdt, fatal)) {
    cerr << "Cannot decode TDT\n";
    return 0;
  }

  return GfpSvmfpEvaluate(fp, models, nmodels, output);
}

int
GfpSvmfpEvaluate(iwstring_data_source& input,
                 SvmModel* models,
                 int nmodels,
                 IWString_and_File_Descriptor& output) { 
  IW_TDT tdt;
  while (tdt.next(input)) {
    if (! GfpSvmfpEvaluate(tdt, models, nmodels, output)) {
      cerr << "Error processing " << tdt << '\n';
      return 0;
    }
  }

  return 1;
}

int
GfpSvmfpEvaluate(const char * fname,
                 SvmModel* models,
                 int nmodels,
                 IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return GfpSvmfpEvaluate(input, models, nmodels, output);
}

int
GfpSvmfpEvaluate(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-M=sfile-cwrite=s-sas=ipos");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('M')) {
    cerr << "Must specify model proto file (-M)\n";
    Usage(1);
  }

  const int nmodels = cl.option_count("M");
  std::unique_ptr<SvmModel[]> models(new SvmModel[nmodels]);
  for (int i = 0; i < nmodels; ++i) {
    IWString fname = cl.string_value('M', i);
    const IWString dir_name = iwmisc::IWDirname(fname);
    GfpModel::SvmfpModel model;
    // cerr << "fname " << fname << " dir_name '" << dir_name << "'\n";
    if (! ReadBinaryProto(fname, model)) {
      cerr << "Cannot read model proto file '" << fname << "'\n"; 
      return 1;
    }
    if (! models[i].Initialise(model, dir_name)) {
      cerr << "Cannot initialise model " << model.ShortDebugString() << '\n';
      return 1;
    }
  }

  if (verbose && nmodels > 1) {
    cerr << "Read " << nmodels << " model protos\n";
  }

  // Either csv or separate options.
  if (cl.option_present("cwrite")) {
    write_label = false;
    write_score = false;
    IWString cwrite;
    for (int i = 0; cl.value("cwrite", cwrite, i); ++i) {
      int j = 0;
      const_IWSubstring token;
      while (cwrite.nextword(token, j, ',')) {
        if (token == "label") {
          write_label = true;
        } else if (token == "score") {
          write_score = true;
        } else if (token == "help") {
          DisplayCWriteOptions(cerr);
        } else {
          cerr << "Unrecognised -cwrite qualifier '" << cwrite << "\n";
          Usage(1);
        }
      }
    }
  }

  if (cl.option_present("sas")) {
    int npoints;
    cl.value("sas", npoints);
    fraction_as_string.set_leading_string(IWString(output_separator));
    fraction_as_string.initialise(-1.2, 1.2, npoints);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  output << "ID";
  for (int i = 0; i < nmodels; ++i) {
    if (models[i].is_regression()) {
      output << output_separator << models[i].response_name();
      continue;
    }
    output << output_separator << models[i].response_name() << "_label";
    if (write_label && write_score) {
      output << output_separator << models[i].response_name() << "_score";
    }
  }
  output << '\n';
  for (const char * fname : cl) {
    if (! GfpSvmfpEvaluate(fname, models.get(), nmodels, output)) {
      cerr << "GfpSvmfpEvaluate:cannot process '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    for (int i = 0; i < nmodels; ++i) {
      models[i].Report(cerr);
    }
  }

  return 0;
}

}  // namespace gfp_svmfp_evaluate

int
main(int argc, char ** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  return gfp_svmfp_evaluate::GfpSvmfpEvaluate(argc, argv);
}
