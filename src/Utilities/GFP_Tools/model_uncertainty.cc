// Make an estimate of the uncertainty associated with a model prediction.

// Inputs are:
// The output from mispredicted - which should probably also include
//   scoring the test set.
// The fingerprints of the training set - although this might also be
// the support vectors.

// For each molecule for which we want an uncertainty estimate
// find the nearest neighbour(s) in the training set.
// Use the errors associated with those predictions to infer a likely
// error for the molecule being predicted.

#include <algorithm>
#include <iostream>
#include <memory>

#include "absl/container/flat_hash_map.h"
#include "absl/strings/string_view.h"

#include "google/protobuf/text_format.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"

#ifdef BUILD_BAZEL
#include "Utilities/General/mispredicted.pb.h"
#else
#endif

#include "gfp.h"

namespace model_uncertainty {

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
  cerr << R"(Generates an estimate of model uncertainty.
Reads a fingerprint file of the molecules being scored - command line argument(s).
This file must contain the predicted value as part of the name field - see below.
Reads a fingerprint file of the training set - might be just the support vectors.
Reads the output from mispredicted to get an idea of how well training set members are predicted,
or how 'hard' it is to make predictions in a region of structure space.

gfp_make train.smi > train.gfp
mispredicted -E train.activity -T Tfile PRED*

some_kind_of_prediction test.smi > test.pred

# need a .gfp file with the prediced value incorporated into the name field.
fetch_smiles_quick test.pred test.smi > test_with_pred.smi
gfp_make test_with_pred.smi > test.gfp

model_uncertainty -T grain.gfp -M Tfile test.gfp > test.uncertainty

The following options are recognised.
 -T <fname>     fingerprint file of the training set.
 -M <fname>     output from the -T option of mispredicted.
 -p <col>       predicted values are in column <col> of the id field.

)";
  // clang-format on

  ::exit(rc);
}

// We use this for finding nearest neighbours.
struct IdDist {
  uint32_t id;
  float distance;
};

// These percentiles are stored here as fractional.
// One for each member of the training set.
struct MPStats {
  float percentile_ave_error;
  float percentile_std_error;
};


class Options {
  private:
    IW_General_Fingerprint* _train;
    uint32_t _ntrain;

    uint32_t _ntest;
    mispredicted::MispredictedContinuous* _mispredicted;

    MPStats* _stats;

    int _number_nbrs;

    // Distances to training set.
    Accumulator<double> _acc_dist;
    Accumulator<double> _acc_score;

    int _verbose;

    int _predicted_column;
    char _input_separator;
    char _output_separator;

  // private functions.
    uint32_t ReadFingerprints(IWString& fname);
    uint32_t ReadFingerprints(iwstring_data_source& fname);

    int ReadMispredictedData(IWString& fname);
    int ReadMispredictedData(iwstring_data_source& fname);
    int ReadMispredictedData(const const_IWSubstring& buffer,
                              mispredicted::MispredictedContinuous& destination);

    int Process(iwstring_data_source& input, IWString_and_File_Descriptor& output);
    int Process(IW_TDT& tdt, IWString_and_File_Descriptor& output);
    int Process(IW_General_Fingerprint& gfp, IWString_and_File_Descriptor& output);
    int Process(const IWString& id,
                 const IdDist* id_distance,
                 IWString_and_File_Descriptor& output);
    int Process(const IWString& id,
                 float activity,
                 const IdDist* id_distance,
                 IWString_and_File_Descriptor& output);

    int GatherMispredictedStatistics();

  public:
    Options();
    ~Options();

    int Initialise(const Command_Line& cl);

    int Process(const char* fname, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _train = nullptr;
  _ntrain = 0;
  _stats = nullptr;

  _ntest = 0;
  _mispredicted = nullptr;

  _number_nbrs = 1;

  _verbose = 0;

  _predicted_column = 1;

  _input_separator = ' ';
  _output_separator = ' ';
}

Options::~Options() {
  if (_train != nullptr) {
    delete [] _train;
  }
  if (_mispredicted != nullptr) {
    delete [] _mispredicted;
  }
  if (_stats != nullptr) {
    delete [] _stats;
  }
}

int
Options::Initialise(const Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('n')) {
    if (! cl.value('n', _number_nbrs) || _number_nbrs < 1) {
      cerr << "Options::Initialise:the number of neighbours (-n) must be a whole +ve number\n";
      Usage(1);
    }

    if (_verbose) {
      cerr << "Will examine the " << _number_nbrs << " nearest neighbours\n";
    }
  }

  if (! cl.option_present('T')) {
    cerr << "Must specify training set fingerprint file via the -T option\n";
    Usage(1);
  }

  if (! cl.option_present('M')) {
    cerr << "Must specify output from mispredicted via the -M option\n";
    Usage(1);
  }

  IWString fname = cl.string_value('T');
  if (! ReadFingerprints(fname)) {
    cerr << "Options::Initialise:cannot read training set fingerprints (-T) '" << fname << "'\n";
    return 0;
  }

  if (cl.option_present('p')) {
    if (! cl.value('p', _predicted_column) || _predicted_column < 1) {
      cerr << "Options::Initialise:the predicted column (-p) must be a whole +ve column number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Predicted values read from column " << _predicted_column << '\n';
    }
    --_predicted_column;
  }

  fname = cl.string_value('M');
  if (! ReadMispredictedData(fname)) {
    cerr << "Options::Initialise:cannot initialise mispredicted data '" << fname << "'\n";
    return 0;
  }

  if (! GatherMispredictedStatistics()) {
    cerr << "Options::Initialise:cannot gather mispredicted statistics\n";
    return 0;
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _ntest << " predictions\n";
  output << "Nearest dists between " << _acc_dist.minval() << " and " << _acc_dist.maxval() <<
            " ave " << static_cast<float>(_acc_dist.average()) << '\n';
  output << "Scores between " << _acc_score.minval() << " and " << _acc_score.maxval() <<
            " ave " << static_cast<float>(_acc_score.average()) << '\n';

  return output.good();
}

uint32_t
Options::ReadFingerprints(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadFingerprints:cannot open training set fingerprint file '" << fname << "'\n";
    return 0;
  }

  return ReadFingerprints(input);
}

uint32_t
Options::ReadFingerprints(iwstring_data_source& input) {
  _ntrain = input.grep("^PCN<");

  if (_ntrain == 0) {
    cerr << "Options::ReadFingerprints:no fingerprints in file\n";
    return 0;
  }

  if (_verbose) {
    cerr << "Options::ReadFingerprints:find " << _ntrain << " fingerprints in training set fingerprint file\n";
  }

  _train = new IW_General_Fingerprint[_ntrain];

  for (uint32_t i = 0; i < _ntrain; ++i) {
    IW_TDT tdt;
    if (! tdt.next(input)) [[unlikely]] {
      cerr << "Options::ReadFingerprints:cannot read fingerprint " << i << "'\n";
      return 0;
    }

    int fatal = 0;
    if (_train[i].construct_from_tdt(tdt, fatal)) {
      continue;
    }

    if (! fatal) {
      continue;
    }
    cerr << "Options::ReadFingerprints:error reading " << tdt << '\n';
    return 0;
  }

  cerr << "Returning " << _ntrain << '\n';
  return _ntrain;
}

int
Options::ReadMispredictedData(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadMispredictedData:cannot open mispredicted data\n";
    return 0;
  }

  return ReadMispredictedData(input);
}

int
Options::ReadMispredictedData(iwstring_data_source& input) {
  _ntest = input.records_remaining();
  if (_ntest == 0) {
    cerr << "Options::ReadMispredictedData:no data\n";
    return 0;
  }

  _mispredicted = new mispredicted::MispredictedContinuous[_ntest];

  const_IWSubstring buffer;

  for (uint32_t i = 0; i < _ntest; ++i) {
    if (! input.next_record(buffer)) {
      cerr << "Options::ReadMispredictedData:cannot read line " << i << '\n';
      return 0;
    }
    if (! ReadMispredictedData(buffer, _mispredicted[i])) {
      cerr << "Options::ReadMispredictedData:invalid data\n";
      cerr << buffer << '\n';
      return 0;
    }
  }

  return _ntest;
}

int
Options::ReadMispredictedData(const const_IWSubstring& buffer,
                              mispredicted::MispredictedContinuous& destination) {

  absl::string_view s(buffer.data(), buffer.length());
  if (! google::protobuf::TextFormat::ParseFromString(s, &destination)) {
    cerr << "Options::ReadMispredictedData:cannot parse textproto\n";
    return 0;
  }

  return 1;
}

struct IdValue {
  uint32_t id;
  float value;
};

int
Options::GatherMispredictedStatistics() {
  if (_ntrain < 100) {
    cerr << "Options::GatherMispredictedStatistics:not enough training data " << _ntrain << '\n';
    return 0;
  }

  _stats = new MPStats[_ntrain];

  std::unique_ptr<IdValue[]> id_value = std::make_unique<IdValue[]>(_ntrain);
  for (uint32_t i = 0; i < _ntrain; ++i) {
    id_value[i].id = i;
    id_value[i].value = _mispredicted[i].average_difference();
  }

  // We ultimately want high percentile values to be good values,
  // so sort in descending order  high->low   these are errors.
  std::sort(id_value.get(), id_value.get() + _ntrain,
       [](const IdValue idv1, const IdValue idv2) {
        return idv1.value > idv2.value;
    });


  uint32_t next_cutpoint = _ntrain / 100;
  uint32_t percentile = 0;
  for (uint32_t i = 0; i < _ntrain; ++i) {
    if (i == next_cutpoint) {
      ++percentile;
      next_cutpoint = (percentile + 1) * _ntrain / 100;
    }

    uint32_t id = id_value[i].id;

    _stats[id].percentile_ave_error = static_cast<float>(percentile) / 100.0f;
  }

  for (uint32_t i = 0; i < _ntrain; ++i) {
    id_value[i].id = i;
    id_value[i].value = _mispredicted[i].std_difference();
  }
  cerr << "Replacing std with max diff\n";
  for (uint32_t i = 0; i < _ntrain; ++i) {
    id_value[i].id = i;
    float minpred = _mispredicted[i].minpred();
    float maxpred = _mispredicted[i].maxpred();
    //float d1 = abs(_mispredicted[i].activity() - minpred);
    //float d2 = abs(_mispredicted[i].activity() - maxpred);

    id_value[i].value = (maxpred - minpred);
  }

  // Sort so high percentile values are the good values - low.
  std::sort(id_value.get(), id_value.get() + _ntrain,
       [](const IdValue& idv1, const IdValue idv2) {
        return idv1.value > idv2.value;
    });


  percentile = 0;
  next_cutpoint = _ntrain / 100;
  for (uint32_t i = 0; i < _ntrain; ++i) {
    if (i == next_cutpoint) {
      ++percentile;
      next_cutpoint = (percentile + 1) * _ntrain / 100;
    }

    uint32_t id = id_value[i].id;
    _stats[id].percentile_std_error = static_cast<float>(percentile) / 100.0f;
    cerr << " std " << i << " ndx " << id << " percentile " << percentile << " value " << _mispredicted[id].std_difference() << '\n';
  }

  return 1;
}

int
Options::Process(const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::Process:cannot open '" << fname << "'\n";
    return 0;
  }

  return Process(input, output);
}

int
Options::Process(iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    if (! Process(tdt, output)) {
      cerr << "Options::Process:fatal error processing\n";
      cerr << tdt << '\n';
      return 0;
    }
  }

  return 1;
}

int
Options::Process(IW_TDT& tdt, IWString_and_File_Descriptor& output) {
  IW_General_Fingerprint fp;
  int fatal = 0;
  if (fp.construct_from_tdt(tdt, fatal)) {
    return Process(fp, output);
  }

  if (! fatal) {
    return 1;
  }

  cerr << "Options::Process:cannot construct fingerprint\n";
  cerr << tdt << '\n';
  return 0;
}

// We have compared two fingerprints and they are identical.
// If the ID's are the same, we do not process that molecule.
// Chances are that `id1` will be of the form
// 123455 3.14
// and `id2` will be of the form
// 123455
int
SameIds(const IWString& id1, const IWString& id2) {
  if (id1 == id2) [[unlikely]] {
    return 1;
  }

  if (id1.starts_with(id2)) {
    return 1;
  }

  // This probably will never happen.
  return id2.starts_with(id1);
}

int
Options::Process(IW_General_Fingerprint& fp, IWString_and_File_Descriptor& output) {
  static std::unique_ptr<IdDist[]> id_distance;

  if (! id_distance) {
    id_distance.reset(new IdDist[_ntrain]);
  }

  uint32_t nsort = 0;
  for (uint32_t i = 0; i < _ntrain; ++i) {
    float d = fp.tanimoto(_train[i]);
    if (d == 1.0f && SameIds(fp.id(), _train[i].id())) {
      continue;
    }

    id_distance[i].id = i;
    id_distance[i].distance = 1.0f - d;
    ++nsort;
  }

  std::partial_sort(id_distance.get(), id_distance.get() + _number_nbrs,
                    id_distance.get() + nsort,
    [](const IdDist& idd1, const IdDist& idd2) {
      return idd1.distance < idd2.distance;
    });

  for (int i = 0; i < _number_nbrs; ++i) {
    cerr << i << " distance " << id_distance[i].distance << '\n';
  }

  return Process(fp.id(), id_distance.get(), output);
}

int
Options::Process(const IWString& id,
                 const IdDist* id_distance,
                 IWString_and_File_Descriptor& output) {
  const_IWSubstring token;
  if (! id.word(_predicted_column, token)) {
    cerr << "Options::Process:cannot extract column " << (_predicted_column + 1) << 
            " from '" << id << "'\n";
    return 0;
  }

  float activity;
  if (! token.numeric_value(activity)) {
    cerr << "Options::Process:invalid numeric '" << token << "'\n";
    return 0;
  }

  return Process(id, activity, id_distance, output);
}

int
Options::Process(const IWString& id,
                 float activity,
                 const IdDist* id_distance,
                 IWString_and_File_Descriptor& output) {
  float result = 0.0f;

  cerr << "Working on '" << id << "' activity " << activity << '\n';
  float d0 = 0.0f;
  for (int i = 0; i < _number_nbrs; ++i) {
    if (i == 0) {
      d0 = id_distance[i].distance;
      _acc_dist.extra(d0);
    }

    float s = 1.0f - id_distance[i].distance;

    uint32_t nbr_id = id_distance[i].id;

    cerr << "   nbr " << i << " nbr_id " << nbr_id << " sim " << s << " ave " << 
            _stats[nbr_id].percentile_ave_error << " std " << _stats[nbr_id].percentile_std_error << '\n';

    float p_ave = _stats[nbr_id].percentile_ave_error;
    float p_std = _stats[nbr_id].percentile_std_error;

    float score = s * p_ave * p_std;

    result += score;
    cerr << "   result incrmented to " << result << '\n';
  }

  result /= static_cast<float>(_number_nbrs);

  _acc_score.extra(result);

//if (d0 > 0.2) {
//  return 1;
//}

  output << id << _output_separator << d0 << _output_separator << result << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vT:M:p:n:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  output << "ID Pred Dist Uncertainty\n";

  for (const char* fname : cl) {
    if (! options.Process(fname, output)) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace model_uncertainty

int
main(int argc, char **argv) {
  int rc = model_uncertainty::Main(argc, argv);

  return rc;
}
