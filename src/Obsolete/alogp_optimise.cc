// Optimise an alogp model with nlopt.

#include <iostream>
#include <optional>
#include <random>
#include <vector>

#include "google/protobuf/descriptor.h"
#include "google/protobuf/message.h"

#include "absl/container/flat_hash_map.h"

#include "nlopt.hpp"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/absl_hash.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/alogp.h"

namespace alogp_optimise {

using std::cerr;

// By convention the Usage function tells how to use the tool.
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
  cerr << R"(Optimises ALogP model parameters for an arbitrary response.
 -O <fname>             the objective - usually measured values for each member of the training set.
 -C <fname>             An alogp::AlogpConfiguration textproto with starting parameter values (optional).
                        Values are initialised randomly if not specified.
 -r <float>             Add random values in (0,<float>) to the starting values - from the -C file.
 -P ...                 Set optimisation parameters. Enter '-P help' for info.
 -S <fname>             Write the result to <fname>. By default 'alogp_optimise.textproto'.
 -W <fname>             Using optimised parameters, write training set predictions to <fname>.
 -w <n>                 Write the -W file every <n> steps.
 -b <bounds>            NLOpt needs an upper and lower bound for each parameter. Default is 2.0
 -v                     verbose output
)";
// clang-format on

  ::exit(rc);
}

// Fill the parameters in `proto` with random values.
void
FillRandomly(alogp::AlogpParameters& proto) {
  std::random_device rd;
  std::mt19937 gen(rd()); 
  std::uniform_real_distribution<float> u(0.0, 1.0);

  const /*Reflection*/ auto* reflection = proto.GetReflection();

  using google::protobuf::Descriptor;
  const Descriptor* descriptor = proto.GetDescriptor();
  cerr << "Randomly filling " << descriptor->field_count() << " fields\n";
  for (int i = 0; i < descriptor->field_count(); ++i) {
    reflection->SetFloat(&proto, descriptor->field(i), u(gen));
  }
}

//
void
AddNoise(float noise, alogp::AlogpParameters& proto) {
  std::random_device rd;
  std::mt19937 gen(rd()); 
  std::uniform_real_distribution<float> u(0.0, noise);

  const /*Reflection*/ auto* reflection = proto.GetReflection();

  using google::protobuf::Descriptor;
  const Descriptor* descriptor = proto.GetDescriptor();
  cerr << "Randomly incrementing " << descriptor->field_count() << " fields\n";
  for (int i = 0; i < descriptor->field_count(); ++i) {
    float v = reflection->GetFloat(proto, descriptor->field(i)) + u(gen);
    reflection->SetFloat(&proto, descriptor->field(i), v);
  }
}

// Transfer parameters from `x` to the fields in `proto`.
void
TransferParameters(const std::vector<double>& x,
                alogp::AlogpParameters& proto) {
  const /*Reflection*/ auto* reflection = proto.GetReflection();

  using google::protobuf::Descriptor;
  const Descriptor* descriptor = proto.GetDescriptor();

  assert(descriptor->field_count() == x.size());

  for (int i = 0; i < descriptor->field_count(); ++i) {
    reflection->SetFloat(&proto, descriptor->field(i), x[i]);
  }
}

// Transfer the parameters in `proto` to `destination`.
void
TransferParameters(const alogp::AlogpParameters& proto,
                   std::vector<double>& destination) {
  const /*Reflection*/ auto* reflection = proto.GetReflection();

  using google::protobuf::Descriptor;
  const Descriptor* descriptor = proto.GetDescriptor();

  destination.resize(descriptor->field_count());

  for (int i = 0; i < descriptor->field_count(); ++i) {
    destination[i] = reflection->GetFloat(proto, descriptor->field(i));
  }
}

struct Data {
  // The training set molecules.
  resizable_array_p<Molecule> molecules;

  // For each molecule, the data needed for fast scoring.
  alogp::ForFastScoring* for_fast_scoring;

  // For each training set molecule the experimental value.
  // This is used for the RMS computation.
  resizable_array<double> expt;

  // The instance that does all the evaluations.
  alogp::ALogP my_alogp;

  uint32_t function_calls = 0;

  // Periodically write a result file.
  uint32_t write_result_every = 0;
  // the name of the result file being periodically written.
  IWString fname;

  Data();
  ~Data();
};

Data::Data() {
  for_fast_scoring = nullptr;
}

Data::~Data() {
  if (for_fast_scoring != nullptr) {
    delete [] for_fast_scoring;
  }
}

int
WriteResultFile(const std::vector<double>& x,
                IWString& fname) {
  alogp::AlogpConfiguration proto;
  TransferParameters(x, *proto.mutable_parameters());

  // alogp::ALogP my_alogp;
  // my_alogp.ConfigFromProto(proto);

  return iwmisc::WriteTextProto<alogp::AlogpConfiguration>(proto, fname);
}

// The objective function called by nlopt;
// We are not using a derivative, so grad will always remain empty.
double
Objective(const std::vector<double>& x,
          std::vector<double>& grad, void* my_func_data) {
  Data* data = reinterpret_cast<Data*>(my_func_data);

  data->function_calls += 1;
  if (data->write_result_every > 0 && data->function_calls % data->write_result_every == 0) {
    WriteResultFile(x, data->fname);
  }

  // const uint32_t n = x.size();
  // cerr << "(OBJECTIVE) DImension " << n << '\n';

  data->my_alogp.SetWeights(x.size(), x.data());

  Accumulator<double> acc;

#ifdef SLOW_VERSION
  // SOme checks to make sure fast and slow versions genererate the same results.
  float max_diff = 0.0;
#endif
  uint32_t nmolecules = data->molecules.size();
  for (unsigned int i = 0; i < nmolecules; ++i) {
#ifdef SLOW_VERSION
    std::optional<double> a = data->my_alogp.LogP(*data->molecules[i]);
    if (! a) {
      continue;
    }
#endif

    float fast = data->my_alogp.LogP(*data->molecules[i], data->for_fast_scoring[i]);
#ifdef SLOW_VERSION
    if (abs(fast - *a) > 0.01) {
      cerr << "eror between methods, comp " << *a << " fast " << fast << '\n';
    }
    if (abs(fast - *a) > max_diff) {
      max_diff = abs(fast - *a);
    }
#endif

    // cerr << ' ' << i << " pred " << *a << " expt " << data->expt[i] << '\n';
    // double diff = *a - data->expt[i];
    double diff = fast - data->expt[i];
    acc.extra(diff * diff);
  }

  if (acc.empty()) {
    cerr << "No data\n";
    return 0.0;
  }

  // cerr << "Objective " << std::sqrt(acc.average()) << ' ' << data->function_calls << " max diff " << max_diff << '\n';
  double ave = acc.average();
  if (ave <= 0.0) {
    cerr << "NEgative average encountered!\n";
    return 9.0;
  }
  cerr << "Objective " << std::sqrt(acc.average()) << ' ' << data->function_calls << '\n';
  return std::sqrt(acc.average());
}

// A subset of the parameters that we process on the command line
class NLOptParameters {
  private:

    // Maximum number of objective function evaluations.
    int _maxeval;

    // A stopping criterion
    double _ftol_rel;

    double _ftol_abs;

    // Maximum elapsed time in seconds.
    double _maxtime;

  public:
    NLOptParameters();

    int Initialise(Command_Line& cl, char flag);

    // Before optimisation starts, transfer our setting to `optimiser`.
    int SetOptimisationParameters(nlopt::opt& optimiser);
};

NLOptParameters::NLOptParameters() {
  _maxeval = 10000;
  _ftol_rel = 0.0f;
  _ftol_abs = 1.0e-05;
  _maxtime = 0.0;
}

int
DisplayDashPOptions(std::ostream& output) {
  output << R"(Sets certain nlopt parameters.
 -P maxeval=nnn         maximum number of funcation evaluations (default 10k).
 -P maxtime=seconds     set maximum execution time (in seconds), 3600 = 1 hour, 43200 = 12 hours...
 -P ftol_abs=<float>    absolute tolerance for objective function evaluations.
)";
  ::exit(0);
}

int
NLOptParameters::Initialise(Command_Line& cl, char flag) {
  const int verbose = cl.option_present('v');

  IWString p;
  for (int i = 0; cl.value(flag, p, i); ++i) {
    if (p.starts_with("maxeval=")) {
      p.remove_leading_chars(8);
      if (! p.numeric_value(_maxeval) || _maxeval < 100) {
        cerr << "The maximum number of objective function evaluations must be a non-small positive number\n";
        return 0;
      }
      if (verbose) {
        cerr << "maxeval " << _maxeval << '\n';
      }
    } else if (p.starts_with("maxtime=")) {
      p.remove_leading_chars(8);
      uint64_t t;
      if (! p.numeric_value(t) || t < 60) {
        cerr << "The maxtime directive must be a valid number of seconds\n";
        return 0;
      }
      _maxtime = static_cast<double>(t);
    } else if (p.starts_with("ftol_abs=")) {
      p.remove_leading_chars(9);
      if (! p.numeric_value(_ftol_abs) || _ftol_abs <= 0.0) {
        cerr << "The function absolute tolerance ftol_abs must > 0.0\n";
        return 0;
      }
      if (verbose) {
        cerr << "Function absolute tolerance " << _ftol_abs << '\n';
      }
    } else if (p == "help") {
      DisplayDashPOptions(cerr);
    } else {
      cerr << "Unrecognised '" << flag << " qualifier '" << p << "'\n";
      DisplayDashPOptions(cerr);
    }
  }

  return 1;
}

int
NLOptParameters::SetOptimisationParameters(nlopt::opt& optimiser) {
  if (_maxeval > 0) {
    optimiser.set_maxeval(_maxeval);
  }
  if (_maxtime > 0.0) {
    optimiser.set_maxtime(_maxtime);
  }
  if (_ftol_rel > 0.0) {
    // implement this...
  }
  if (_ftol_abs > 0) {
    optimiser.set_ftol_abs(_ftol_abs);
  }

  return 1;
}

// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
  private:
    int _verbose = 0;

    Data _data;

    // The file with the objective values might have missing values.
    // If we choose to ignore such cases, we will have to adjust the
    // array of molecules. Decided that would be too complicated here
    // and easy for the user to fix. So while this variable exists
    // it is never set true.
    int _ignore_missing_experimental_values;

    // The -P option.
    NLOptParameters _nlopt_params;

    // We read this with the -C option or fill it randomly.
    // The -C option.
    alogp::AlogpConfiguration _proto;

    // Given a starting configuration, add a measure of random noise to that.
    // The -r option.
    float _add_noise_to_starting_config;

    // We box all paramters between in [-bound,_bound]. The -b option.
    double _bound;

    // Once we have a solution, we can write te training set predictions
    // to this file.
    IWString _training_set_predictions_fname;

    // During long running optimisations it can be convenient to have
    // result files written periodically. The -s option.
    uint32_t _write_result_every;

    // Private functions
    int ReadMolecules(const char* fname, FileType input_type);
    int ReadMolecules(data_source_and_type<Molecule>& input);

    int ReadExperimentalValues(IWString& fname);
    int ReadExperimentalValues(iwstring_data_source& input);

    int ProcessDashCOption(Command_Line& cl, char flag);

    int GatherAssignedTypes();

    uint32_t GetParameterDimension() const;

    int WriteResult(const std::vector<double>& x, IWString & fname);
    int WriteTrainingSetPredictions(std::vector<double>& x);

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    // Optimise the problem and write an alogp::AlogpParameters textproto
    // to `fname`. Writing happens even if the computation has not completed
    // successfully - the failure is reported on stderr.
    int Optimise(IWString& fname);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _bound = 2.0;
  _ignore_missing_experimental_values = 0;
  _add_noise_to_starting_config = 0.0;
  _write_result_every = 0;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('P')) {
    if (! _nlopt_params.Initialise(cl, 'P')) {
      cerr << "Cannot initialise nlopt options (-P)\n";
      return 0;
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments, must specify smiles file\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (! ReadMolecules(cl[0], input_type)) {
    cerr << "Cannot read molecules from '" << cl[0] << "'\n";
    return 0;
  }

  if (! cl.option_present('O')) {
    cerr << "Must specify objective values via the -O option\n";
    Usage(1);
  }

  if (cl.option_present('O')) {
    IWString fname = cl.string_value('O');
    if (! ReadExperimentalValues(fname)) {
      cerr << "Cannot read experimental values from '" << fname << "'\n";
      return 0;
    }
  }

  if (! ProcessDashCOption(cl, 'C')) {
    return 0;
  }

  if (cl.option_present('b')) {
    if (! cl.value('b', _bound) || _bound <= 0.0) {
      cerr << "The bounds value (-b) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Bound on all values " << _bound << '\n';
    }
  }

  if (! GatherAssignedTypes()) {
    cerr << "Cannot gather assigned atom types\n";
  }

  if (cl.option_present('W')) {
    cl.value('W', _training_set_predictions_fname);
    if (_verbose) {
      cerr << "Training set predictions will be written to '" << _training_set_predictions_fname << "'\n";
    }
  }

  if (cl.option_present('s')) {
    uint32_t s;
    if (!cl.value('s', s) || s < 100) {
      cerr << "The write a result every (-s) option must be a reasonable large +ve number\n";
      return 0;
    }

    _data.write_result_every = s;

    if (_verbose) {
      cerr << "Will write the result file (-S) every " << s << " stems\n";
    }
  }

  return 1;
}

int
Options::ProcessDashCOption(Command_Line& cl, char flag) {

  if (cl.option_present(flag)) {
    IWString fname = cl.string_value(flag);
    std::optional<alogp::AlogpConfiguration> maybe_proto = 
        iwmisc::ReadTextProto<alogp::AlogpConfiguration>(fname);
    if (! maybe_proto) {
      cerr << "Options::ProcessDashCOption:cannot read textproto '" << fname << "'\n";
      return 0;
    }

    if (! _data.my_alogp.ConfigFromProto(*maybe_proto)) {
      cerr << "Cannot build alogp from proto config file '" << fname << "'\n";
      return 0;
    }
    _proto = *maybe_proto;
  } else {
    FillRandomly(*_proto.mutable_parameters());
    _data.my_alogp.ConfigFromProto(_proto);
  }

  if (cl.option_present('r')) {
    float r;
    if (! cl.value('r', r) || r <= 0.0f) {
      cerr << "The randomise start config option (-r) must be positive\n";
      return 0;
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Trained on " << _data.molecules.size() << " molecules\n";
  output << "Performed " << _data.function_calls << " model evaluations\n";

  return 1;
}

int
Options::ReadMolecules(data_source_and_type<Molecule>& input) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    _data.molecules << m;
  }

  return _data.molecules.size();
}

int
Options::ReadMolecules(const char* fname,
              FileType input_type) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "ReadMolecules:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadMolecules(input);
}

int
Options::ReadExperimentalValues(IWString& fname) {
  iwstring_data_source input(fname);

  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadExperimentalValues(input);
}

int
Options::ReadExperimentalValues(iwstring_data_source& input) {
  // for each ID, where is it in the _molecules array.
  absl::flat_hash_map<IWString, uint32_t> id_to_ndx;
  const int n = _data.molecules.size();

  _data.expt.extend(n, 0.0);

  for (int i = 0; i < n; ++i) {
    const IWString& id = _data.molecules[i]->name();
    if (id.contains(' ')) {
      IWString tmp(id);
      tmp.truncate_at_first(' ');
      id_to_ndx[tmp] = i;
    } else {
      id_to_ndx[id] = i;
    }
  }

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

    _data.expt[f->second] = v;
  }

#ifdef ECHO_OBJECTIVE
  for (int i = 0; i< _data.expt.number_elements(); ++i) {
    cerr << _data.molecules[i]->name() <<  " expt " << _data.expt[i] << '\n';
  }
#endif

  return _data.expt.size();
}

// For each atom perform a computation, requesting that the assigned atom
// type be returned as the isotopic label.
// For successful calculations store those assigned types.
int
Options::GatherAssignedTypes() {
  assert(_data.atype == nullptr);

  const uint32_t n = _data.molecules.size();

  _data.for_fast_scoring = new alogp::ForFastScoring[n];

  uint32_t failures = 0;

  for (uint32_t i = 0; i < n; ++i) {
    Molecule& m = *_data.molecules[i];
    std::optional<double> a = _data.my_alogp.LogP(m);
    if (! a) {
      ++failures;
      continue;
    }

    _data.my_alogp.FillForFastScoring(m, _data.for_fast_scoring[i]);
  }

  if (failures) {
    cerr << "Options::GatherAssignedTypes:encountered " << failures << " failed computations\n";
  }

  return 1;
}

// Given the contents of _proto, set lower bounds in `bounds`
uint32_t
Options::GetParameterDimension() const {
  using google::protobuf::Descriptor;
  const Descriptor* descriptor = _proto.parameters().GetDescriptor();

  return descriptor->field_count();
}

int
Options::Optimise(IWString& fname) {

  // If we are writing files, transfer the file name.
  if (_data.write_result_every > 0) {
    _data.fname = fname;
  }

  const uint32_t n = GetParameterDimension();
  cerr << "Dimension " << n << '\n';

  nlopt::opt myopt(nlopt::LN_COBYLA, n);
  myopt.set_min_objective(Objective, reinterpret_cast<void*>(&_data));

  std::vector<double> bounds(n, -_bound);
  myopt.set_lower_bounds(bounds);

  bounds.assign(n, _bound);
  myopt.set_upper_bounds(bounds);

  _nlopt_params.SetOptimisationParameters(myopt);

  std::vector<double> x;
  x.reserve(n);
  TransferParameters(_proto.parameters(), x);

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

//cerr << "Final parameters\n";
//for (uint32_t i = 0; i < n; ++i) {
//  cerr << i << ' ' << x[i] << '\n';
//}
  TransferParameters(x, *_proto.mutable_parameters());

  if (_training_set_predictions_fname.length() > 0) {
    WriteTrainingSetPredictions(x);
  }

  alogp::ALogP foo;
  foo.ConfigFromProto(_proto);
  const int nmolecules = _data.molecules.number_elements();
  for (int i = 0; i < nmolecules; ++i) {
    Molecule& m = *_data.molecules[i];
    std::optional<float> a = foo.LogP(m);
    cerr << m.name() << ' ' << *a << " FROMPROTO\n";
  }

  return WriteResult(x, fname);
}

int
Options::WriteResult(const std::vector<double>& x, IWString & fname) {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "Options::WriteResult:cannot open '" << fname << "'\n";
    return 0;
  }

#ifdef DEBUG_WRITE_RESULT
  cerr << "Optimisation parameters\n";
  const int n = GetParameterDimension();
  for (int i = 0; i < n; ++i) {
    cerr << i << " result " << x[i] << '\n';
  }
#endif

  TransferParameters(x, *_proto.mutable_parameters());

  return iwmisc::WriteTextProto<alogp::AlogpConfiguration>(_proto, fname);
}

int
Options::WriteTrainingSetPredictions(std::vector<double>& x) {
  std::vector<double> grad;
  double q = Objective(x, grad, &_data);
  cerr << "Final objective " << q << '\n';

  IWString_and_File_Descriptor output;
  if (! output.open(_training_set_predictions_fname)) {
    cerr << "Options::WriteTrainingSetPredictions:cannot open '" << _training_set_predictions_fname << "'\n";
    return 0;
  }

  constexpr char kSep = ' ';

  Accumulator<double> acc;
  uint32_t nmolecules = _data.molecules.size();
  for (uint32_t i = 0; i < nmolecules; ++i) {
    Molecule& m = *_data.molecules[i];
    std::optional<double> q = _data.my_alogp.LogP(m);
    if (! q) {
      continue;
    }

    append_first_token_of_name(m.name(), output);
    output << kSep << *q << '\n';
    output.write_if_buffer_holds_more_than(4096);

    double diff = (*q - _data.expt[i]);
    acc.extra(diff * diff);
  }

  cerr << "Rescored " << std::sqrt(acc.average()) << '\n';

  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:C:O:S:P:W:b:r:s:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }
  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  IWString fname;
  if (cl.option_present('S')) {
    cl.value('S', fname);
  } else {
    fname = "alogp_opt.textproto";
  }
  options.Optimise(fname);

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace alogp_optimise

int
main(int argc, char ** argv) {

  int rc = alogp_optimise::Main(argc, argv);

  return rc;
}
