// Optimise a gfp linear model with NLOpt.

#include <iostream>
#include <random>
#include <vector>

#include "absl/container/flat_hash_map.h"

#include "nlopt.hpp"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.h"

#ifdef BUILD_BAZEL
#include "Obsolete/gfp_linear_model.pb.h"
#else
#include "gfp_linear_model.pb.h"
#endif

namespace gfp_linear_model {

using std::cerr;

void Usage(int rc) {
  // clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(Evaluates a gfp linear model build with gfp_linear_model_optimise.
 -C <fname>     a gfp_linear_model::GfpLinearModel textproto model build by gfp_linear_model_optimise.
 -F ...         standard gfp options, enter '-F help' for info.
 -v             verbose output.
)";
  // clang-format on
}

class Optimise {
  private:
    int _number_fingerprints;
    IW_General_Fingerprint* _fp;
    
    int _number_bits;
    // For each fingerprint, a _number_bits fixed size array of the feature counts for that fingerprint.
    int** _feature_count;

    double _bias;

    // For each fingerprint the experimental value.
    double* _expt;

    uint64_t _write_every;

    // The -S option.
    IWString _result_fname;

  // Private functions.
    int BuildFingerprints(const char* fname);
    int BuildFingerprints(iwstring_data_source& fname);
    int ExtractFeatures(const IWDYFP& fp, int* feature_count);
    int ExtractFeatures(const Sparse_Fingerprint& fp, int* feature_count);
    int ExtractFeatures(const IW_General_Fingerprint& fp, int* feature_count);

    int ReadExperimentalValues(IWString& fname);
    int ReadExperimentalValues(iwstring_data_source& input);


  public:
    Optimise();
    ~Optimise();

    int Initialise(Command_Line& cl);

    // Will be the objective function for NLOpt;
    double Score(const std::vector<double>& x);

    int number_bits() const {
      return _number_bits;
    }

    uint64_t write_every() const {
      return _write_every;
    }

    // Write results to _result_fname.
    int WriteTextProto(const std::vector<double>& weight);
};

Optimise::Optimise() {
  _number_fingerprints = 0;
  _fp = nullptr;

  _number_bits = 1024;

  _feature_count = nullptr;

  _bias = 0.0;

  _expt = nullptr;

  _write_every = 0;
}

Optimise::~Optimise() {
  if (_fp != nullptr) {
    delete [] _fp;
  }

  for (int i = 0; i < _number_fingerprints; ++i) {
    delete [] _feature_count[i];
  }
  delete [] _feature_count;

  delete [] _expt;
}

int
Optimise::Initialise(Command_Line& cl) {
  const int verbose = cl.option_present('v');

  if (! BuildFingerprints(cl[0])) {
    cerr << "Optimise::Initialise:cannot read fingerrpints '" << cl[0] << "'\n";
    return 0;
  }

  if (cl.option_present('b')) {
    if (! cl.value('b', _number_bits) || _number_bits < 100) {
      cerr << "The umber of bits to use (-b) must be a reasonably large +ve number\n";
      Usage(1);
    }
    if (verbose) {
      cerr << "Will use " << _number_bits << " bits\n";
    }
  }

  if (! cl.option_present('E')) {
    cerr << "Must specify training set activity file via the -E option\n";
    Usage(1);
  }

  if (cl.option_present('E')) {
    IWString fname = cl.string_value('E');
    if (! ReadExperimentalValues(fname)) {
      cerr << "Cannot read experimental data from '" << fname << "'\n";
      return 0;
    }
  }

  if (! cl.option_present('S')) {
    cerr << "Must provide an output file name via the -S option\n";
    Usage(1);
  }

  if (cl.option_present('S')) {
    _result_fname = cl.string_value('S');
    if (verbose) {
      cerr << "Resulting proto written to " << _result_fname << "'\n";
    }
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', _write_every) || _write_every < 100) {
      cerr << "The write a result every (-s) option must be a reasonable large +ve number\n";
      return 0;
    }

    if (verbose) {
      cerr << "Will write the result file (-S) every " << _write_every << " stems\n";
    }
  }

  return 1;
}

int
Optimise::BuildFingerprints(const char* fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return BuildFingerprints(input);
}

int
Optimise::BuildFingerprints(iwstring_data_source& input) {
  _number_fingerprints = input.count_records_starting_with("PCN<");
  if (_number_fingerprints == 0) {
    cerr << "Optimise::BuildFingerprints:no fingerprints in input\n";
    return 0;
  }

  _fp = new IW_General_Fingerprint[_number_fingerprints];
  _feature_count = new int*[_number_fingerprints];

  for (int i = 0; i < _number_fingerprints; ++i) {
    _feature_count[i] = new_int(_number_bits);
  }

  IW_TDT tdt;
  for (int ndx = 0; tdt.next(input); ++ndx) {
    int fatal;

    if (!_fp[ndx].construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        return 0;
      }
      continue;
    }

    if (! ExtractFeatures(_fp[ndx], _feature_count[ndx])) {
      return 0;
    }
  }

  return 1;
}

int
Optimise::ExtractFeatures(const IW_General_Fingerprint& fp, int* feature_count) {
  for (int i = 0; i < number_fingerprints(); ++i) {
    ExtractFeatures(fp[i], feature_count);
  }

  for (int i = 0; i < number_sparse_fingerprints(); ++i) {
    ExtractFeatures(fp.sparse_fingerprint(i), feature_count);
  }

  return 1;
}

int
Optimise::ExtractFeatures(const IWDYFP& fp, int* feature_count) {
  int i = 0;
  int b;
  while ((b = fp.next_on_bit(i)) >= 0) {
    // If we knew that we were the only fingerprint this could be an assignment...
    feature_count[b % _number_bits] += 1;
  }

  return 1;
}

int
Optimise::ExtractFeatures(const Sparse_Fingerprint& fp, int* feature_count) {
  int i = 0;
  uint32_t bit;
  int count;
  while (fp.next_bit_set(i, bit, count)) {
    feature_count[bit % _number_bits] += count;
  }
  return 1;
}

int
Optimise::ReadExperimentalValues(IWString& fname) {
  iwstring_data_source input(fname);

  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadExperimentalValues(input);
}

int
Optimise::ReadExperimentalValues(iwstring_data_source& input) {
  // for each ID, where is it in the _fp array.
  absl::flat_hash_map<IWString, uint32_t> id_to_ndx;

  _expt = new double[_number_fingerprints];

  cerr << "Reading expt for " << _number_fingerprints << " fingerprints\n";
  for (int i = 0; i < _number_fingerprints; ++i) {
    const IWString& id = _fp[i].id();
    cerr << "id '" << id << "' i " << i << '\n';
    if (id.contains(' ')) {
      IWString tmp(id);
      tmp.truncate_at_first(' ');
      id_to_ndx[tmp] = i;
    } else {
      id_to_ndx[id] = i;
    }
  }
  cerr << "id_to_ndx has " << id_to_ndx.size() << " entries\n";

  static constexpr char kSep = ' ';

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (input.lines_read() == 1) {
      continue;
    }

    IWString id, val;
    int i = 0;
    if (! buffer.nextword(id, i, kSep) || id.empty() ||
        ! buffer.nextword(val, i, kSep) || val.empty()) {
      cerr << "ReadExperimentalValues:invalid input '" << buffer << "'\n";
      return 0;
    }

    auto f = id_to_ndx.find(id);
    if (f == id_to_ndx.end()) {
      cerr << "Extra expt value " << id << " ignored\n";
      continue;  // extra items in the activity file is not a problem.
//    cerr << "ReadExperimentalValues:no molecule for '" << id << "'\n";
//    cerr << "Remove training set molecules without data and try again\n";
//    cerr << "fetch_smiles_quick train.dat train.smi > newtrain.smi\n";
//    return 0;
    }

    double v;
    if (! val.numeric_value(v)) {
      cerr << "Invalid numeric '" << buffer << "'\n";
      return 0;
    }
    cerr << "id " << id << " goes to " << f->second << '\n';

    _expt[f->second] = v;
  }

#ifdef ECHO_OBJECTIVE
  for (int i = 0; i< _data.expt.number_elements(); ++i) {
    cerr << _data.molecules[i]->name() <<  " expt " << _data.expt[i] << '\n';
  }
#endif

  return 1;
}

int
Optimise::WriteTextProto(const std::vector<double>& weight) {
  gfp_linear_model::GfpLinearModel proto;
  const int n = _number_bits + 1;

  proto.mutable_weight()->Reserve(n);
  for (int i = 0; i < n; ++i) {
    proto.set_weight(i, weight[i]);
  }
  return iwmisc::WriteTextProto<gfp_linear_model::GfpLinearModel>(proto, _result_fname);
}

double
Optimise::Score(const std::vector<double>& x) {
  const double bias = x[_number_bits + 1];

  Accumulator<double> acc;
  for (int i = 0; i < _number_fingerprints; ++i) {
    double result = bias;
    const int* f = _feature_count[i];
    for (int j = 0; j < _number_bits; ++j) {
      if (f[j] == 0) {
        continue;
      }
      if (f[j] == 1) {
        result += x[j];
      } else {
        result += f[j] * x[j];
      }
    }

    double diff = result - _expt[i];
    // cerr << i << " result " << result << " expt " << _expt[i] << '\n';
    acc.extra(diff * diff);
  }

  return std::sqrt(acc.average());
}

double
Objective(const std::vector<double>& x,
          std::vector<double>& grad, void* my_func_data) {
  static uint64_t function_evaluations = 0;

  Optimise* data = reinterpret_cast<Optimise*>(my_func_data);

  const double result = data->Score(x);

  ++function_evaluations;
  cerr << "Objective " << result << ' ' << function_evaluations << '\n';

  if (data->write_every() > 0 && (function_evaluations % data->write_every()) == 0) {
    data->WriteTextProto(x);
  }

  return result;
}

int
DoOptimisation(Optimise& opt) {
  const int n = opt.number_bits() + 1;
  std::vector<double> x(n);

  std::random_device rd;
  std::mt19937 gen(rd()); 
  std::uniform_real_distribution<float> u(-0.5, 0.5);

  for (int i = 0; i < n; ++i) {
    x[i] = u(gen);
  }

  nlopt::opt myopt(nlopt::LN_COBYLA, n);
  myopt.set_min_objective(Objective, reinterpret_cast<void*>(&opt));

  std::vector<double> bounds(n, 2.0);
  myopt.set_upper_bounds(bounds);

  bounds.assign(n, -2.0);
  myopt.set_lower_bounds(bounds);

  myopt.set_maxeval(100000);
  myopt.set_ftol_abs(1.0e-04);

  double opt_f;
  const int rc = myopt.optimize(x, opt_f);

  cerr << "Return code " << rc << " opt_f " << opt_f << '\n';

  if (rc == nlopt::SUCCESS) {
    cerr << "Success " << opt_f << '\n';
  } else if (rc == nlopt::MAXEVAL_REACHED) {
    cerr << "Max number of evaluations reached " << opt_f << '\n';
  } else if (rc == nlopt::STOPVAL_REACHED) {
     cerr << "Stopping value achieved " << opt_f << '\n';
  } else if (rc == nlopt::FTOL_REACHED) {
     cerr << "FTOL reached " << opt_f << '\n';
  } else if (rc == nlopt::MAXTIME_REACHED) {
     cerr << "Maxtime reached " << opt_f << '\n';
  } else {
    cerr << "Return code " << rc << '\n';
  }

  cerr << "Final parameters\n";
  for (int i = 0; i < n; ++i) {
    cerr << i << ' ' << x[i] << '\n';
  }

  return opt.WriteTextProto(x);
}

int
Main(int argc, char** argv) {
   
  Command_Line cl(argc, argv, "vF:S:b:E:s:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Insufficient arguments, must specify training set fingerprint input file\n";
    Usage(2);
  }

  if (need_to_call_initialise_fingerprints(cl)) {
    if (!initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise GFP options\n";
      Usage(23);
    }
  }

  Optimise optimise;

  if (! optimise.Initialise(cl)) {
    cerr << "Cannot initialise optimisation conditions\n";
    return 1;
  }

  DoOptimisation(optimise);

  return 0;
}

}  // namespace gfp_linear_model

int
main(int argc, char** argv) {
  int rc = gfp_linear_model::Main(argc, argv);

  return rc;
}
