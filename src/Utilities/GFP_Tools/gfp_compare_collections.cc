// Compare two collections via fingerprint distances.
// Stop when the distribution of distances converges.

#include <stdlib.h>

#include <algorithm>
#include <iostream>

#define USE_OMP
#ifdef USE_OMP
#include "omp.h"
#endif

#define REPORT_PROGRESS_IMPLEMENTATION
#include "re2/re2.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Utilities/GFP_Tools/gfp_standard.h"

namespace gfp_compare_collections {

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
  cerr << R"(Compares a set of fingerprints against another, sampling the profile of distances until converged.
 -p <fname>     fingerprint file of the comparison set - likely Chembl or some other reference set.
 -s <n>         number of items in the -p file. This tool will count if not specified.
 -T <dist>      upper distance threshold - ignore distances larger than <dist>.
 -n <number>    check convergence and write an output file every <number> fingerprints (default 1000).
 -f <number>    first time to check for convergence. If not specified it will be 1000.
                Convergence is assessed either by absolute (-a) or relative (-r) tolerance.
                A comparison is made between the previous distribution and the current values.
 -r <tol>       Relative tolerance. Converged if every bucketised part of the distribution is within <tol>.
 -a <tol>       Absolute tolerance. Converged if every fractional difference is less than <tol>.
 -S <stem>      output files are written as <stem>0.csv, <stem>1.csv...
 -b <size>      batch size. Query fingerprints are processed in groups - default 10.
 -h <threads>   number of OMP threads to use.
 -v             verbose output
)";

  exit(rc);
}

constexpr int kNBins = 1001;

class Options {
  private:
    int _verbose;
    uint64_t _fingerprints_read;

    float _max_distance;

    uint64_t _next_sample;
    uint64_t _sample_interval;

    double _absolute_tolerance;
    double _relative_tolerance;

    GFP_Standard* _fp;
    uint32_t _number_fingerprints;

    // From ChatGPT we find that efficiency improves if we can amortize the
    // overhead of thread creation across multiple fingerprints.
    // The initial version processed fingerprints one at a time and parallel
    // processing happened with each fingerprint.
    uint64_t _batch_size = 10;

    IWString _stem;
    int _next_fname_index;

    uint64_t _acc[kNBins];
    double _distribution[kNBins];
    double _previous_distribution[kNBins];

  // Private functions
    int ReadComparisonSet(IWString& fname);
    int ReadComparisonSet(iwstring_data_source& input);
    void AddToBinnedDistances(float d);
    int Process(iwstring_data_source& input, int& converged);
    int Process(GFP_Standard& gfp_standard, int& converged);
    int WriteDistribution(uint64_t tot);
    int WriteDistribution();

    int IsConverged();
    int AtolConverged() const;
    int RtolConverged() const;
    void ProcessBatch(const GFP_Standard* queries, uint32_t n,
                      int& converged);

  public:
    Options();
    ~Options();

    int Initialise(Command_Line& cl);

    int Process(const char* fname, int& converged);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _fingerprints_read = 0;
  _next_sample = 1000;
  _sample_interval = 1000;
  _max_distance = 1.0f;
  _absolute_tolerance = 0.0;
  _relative_tolerance = 0.0;
  _fp = nullptr;
  _number_fingerprints = 0;
  _next_fname_index = 0;

  std::fill_n(_acc, kNBins, 0);
  std::fill_n(_distribution, kNBins, 0.0);
  std::fill_n(_previous_distribution, kNBins, 0.0);
}

Options::~Options() {
  if (_fp != nullptr) {
    delete [] _fp;
  }
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('s')) {
    if (! cl.value('s', _number_fingerprints) || _number_fingerprints == 0) {
      cerr << "The number of fingerprints in the -p file (-s option) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Fingerprint array sized to " << _number_fingerprints << " fingerprints\n";
    }
  }

  if (cl.option_present('p')) {
    IWString fname = cl.string_value('p');
    if (! ReadComparisonSet(fname)) {
      cerr << "Cannot read comparison set '" << fname << "'\n";
      return 0;
    }
  }

  if (cl.option_present('S')) {
    cl.value('S', _stem);
    if (_verbose) {
      cerr << "output files written with stem '" << _stem << "'\n";
    }
  }

  if (cl.option_present('f')) {
    if (! cl.value('f', _next_sample) || _next_sample < 10) {
      cerr << "The first check option (-f) must be a reasonably large +ve integer\n";
      return 0;
    }

    if (_verbose) {
      cerr << "First check will be after " << _next_sample << " fingerprints\n";
    }
  }

  if (cl.option_present('n'))  {
    if (! cl.value('n', _sample_interval) || _sample_interval == 0) {
      cerr << "The sample interval (-n) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will check for convergence every " << _sample_interval << " fingerprints\n";
    }
  }

  if (cl.option_present('a') && cl.option_present('r')) {
    cerr << "Only one of absolute (-a) or relative (-r) tolerance can be specified\n";
    return 0;
  }

  if (cl.option_present('a')) {
    if (! cl.value('a', _absolute_tolerance) || _absolute_tolerance <= 0.0) {
      cerr << "The convergence tolerance (-a) option must be a +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will be converged if distribution remains within " << _absolute_tolerance << '\n';
    }
  }

  if (cl.option_present('r')) {
    if (! cl.value('r', _relative_tolerance) || _relative_tolerance <= 0.0) {
      cerr << "The relative tolerance convergence criterion (-r) must be a positive value\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Relative tolerance " << _relative_tolerance << '\n';
    }
  }

  if (_absolute_tolerance == 0.0 && _relative_tolerance == 0.0) {
    _relative_tolerance = 0.01;
  }

  if (cl.option_present('T')) {
    if (! cl.value('T', _max_distance) || _max_distance <= 0.0f || _max_distance > 1.0f) {
      cerr << "The maximum distance value (-T) must be a valid distance\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will only record distances <= " << _max_distance << '\n';
    }
  }

  if (cl.option_present('b')) {
    if (! cl.value('b', _batch_size) || _batch_size < 1) {
      cerr << "The batch size parameter must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Query fingerprings processed in batches of size " << _batch_size << '\n';
    }
  }

  return 1;
}

int
Options::ReadComparisonSet(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadComparisonSet:cannot open '" << fname << "'\n";
    return 0;
  }

  if (_number_fingerprints == 0) {
    RE2 vbar("^|");
    _number_fingerprints = input.grep(vbar);
    if (_number_fingerprints == 0) {
      cerr << "Options::ReadComparisonSet:no fingerprints\n";
      return 0;
    }
  }

  _fp = new GFP_Standard[_number_fingerprints];

  return ReadComparisonSet(input);
}

int
FromGfp(const IW_General_Fingerprint& gfp,
        GFP_Standard& destination) {

  destination.build_molecular_properties(gfp.molecular_properties_integer());
  destination.build_iw(gfp[0]);
  destination.build_mk(gfp[1]);
  destination.build_mk2(gfp[2]);

  return 1;
}

int
FromGfp(IW_TDT& tdt,
        GFP_Standard& destination,
        int& fatal) {
  IW_General_Fingerprint gfp;
  fatal = 0;
  if (!gfp.construct_from_tdt(tdt, fatal)) {
    cerr << "Options::ReadComparisonSet:cannot read fingerprint\n";
    return 0;
  }

  return FromGfp(gfp, destination);
}

int
Options::ReadComparisonSet(iwstring_data_source& input) {
  IW_TDT tdt;
  int fatal;
  uint32_t ndx = 0;
  for (; ndx < _number_fingerprints && tdt.next(input); ++ndx) {
    if (! FromGfp(tdt, _fp[ndx], fatal)) {
      if (fatal) {
        cerr << "Options::ReadComparisonSet:error reading fingerprint, line " << input.lines_read() << '\n';
        return 0;
      }
      break;
    }
  }

  _number_fingerprints = ndx;

  if (_verbose) {
    cerr << "Read " << _number_fingerprints << " fingerprints\n";
  }

  return 1;
}

int
Options::Process(const char* fname, int& converged) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::Process:cannot open '" << fname << "'\n";
    return 0;
  }

  return Process(input, converged);
}

int
Options::Process(iwstring_data_source& input, int& converged) {
  std::unique_ptr<GFP_Standard[]> batch = std::make_unique<GFP_Standard[]>(_batch_size + 1);

  IW_TDT tdt;
  int fatal = 0;
  uint32_t ndx = 0;
  while (tdt.next(input)) {
    if (FromGfp(tdt, batch[ndx], fatal)) {
    } else if (fatal) {
      cerr << "Options::Process:error reading fingerprint, line "
           << input.lines_read() << '\n';
      return 0;
    } else {
      break;
    }

    ++ndx;

    if (ndx == _batch_size) {
      ProcessBatch(batch.get(), ndx, converged);
      if (converged) {
        ndx = 0;   // batch has been fully processed, now must appear empty.
        break;
      }

      ndx = 0;
    }
  }

  if (ndx > 0) {
    ProcessBatch(batch.get(), ndx, converged);
  }

  if (converged) {
    cerr << "Converged\n";
    WriteDistribution();
  }

  return 1;
}

void
Options::ProcessBatch(const GFP_Standard* queries, uint32_t n,
                      int& converged) {
  if (_verbose > 1) {
    cerr << "ProcessBatch, processed " << _fingerprints_read << " fingerprints\n";
  }
  
  // Cache member variables - ChatGPT suggestion.
  const float max_distance = _max_distance;
  const uint32_t nfp = _number_fingerprints;
  GFP_Standard* fp = _fp;

#pragma omp parallel
  {
    uint64_t local_acc[kNBins] = {0};

#pragma omp for schedule(static)
    for (uint32_t q = 0; q < n; ++q) {
      const GFP_Standard& query = queries[q];

      for (uint32_t i = 0; i < nfp; ++i) {
        const float d = query.tanimoto_distance(fp[i]);
        if (d > max_distance) {
          continue;
        }

        const int ndx = static_cast<int>(d * 1000.0f);
        ++local_acc[ndx];
      }
    }

#pragma omp critical
    {
      for (int i = 0; i < kNBins; ++i) {
        _acc[i] += local_acc[i];
      }
    }
  }

  _fingerprints_read += n;

  if (_fingerprints_read < _next_sample) {
    return;
  }

  if (IsConverged()) {
    cerr << "IsConverged passed\n";
    converged = 1;
  }

  _next_sample = _fingerprints_read + _sample_interval;
}

// Returns true if we have achieved convergence.
// Not const because it updates _distribution and other variables.
int
Options::IsConverged() {
  if (_next_fname_index > 0) {
    std::copy_n(_distribution, kNBins, _previous_distribution);
  }

  uint64_t tot = 0;
  for (int i = 0; i < kNBins; ++i) {
    tot += _acc[i];
  }

  if (tot == 0) [[ unlikely ]] { 
    cerr << "Options::IsConverged:no distances\n";
    return 0;
  }

  const double multiplier = 1.0 / static_cast<double>(tot);
  for (int i = 0; i < kNBins; ++i) {
    _distribution[i] = static_cast<double>(_acc[i]) * multiplier;
  }

  if (_next_fname_index == 0) {
    // No previous against which to compare.
  } else if (_absolute_tolerance > 0.0) {
    if (AtolConverged()) {
      return 1;
    }
  } else if (_relative_tolerance > 0.0) {
    if (RtolConverged()) {
      return 1;
    }
  }

  cerr << "Not converged, WriteDistribution\n";
  WriteDistribution(tot);

  return 0;
}

// Return true if the distribution is converted via absolute diff.
// This is slightly inefficient, since we could return false once we
// have found a difference greater than _absolute_tolerance.
int
Options::AtolConverged() const {
  double max_diff = 0.0;

  for (int i = 0; i < kNBins; ++i) {
    const double d = std::abs(_distribution[i] - _previous_distribution[i]);
    if (d > max_diff) {
      max_diff = d;
    }
  }

  if (max_diff <= _absolute_tolerance) {
    return 1;
  }

  if (_verbose) {
    cerr << _fingerprints_read << " fingerprints read, rtol max_diff " << max_diff;
  }

  return 0;
}

// Returns true of the distribution is converged via the relative
// tolerance criteria.
// Slightly inefficient because once a failure in rtol is detected
// we should return 0 immediately, but we keep going to if _verbose
// is set, we can return the max_rtol.
int
Options::RtolConverged() const {
  int failed = 0;
  double max_rtol = 0.0;

  for (int i = 0; i < kNBins; ++i) {
    if (_previous_distribution[i] == 0.0 && _distribution[i] > 0.0) {
      failed = 1;
      continue;
    }

    const double d = std::abs(_distribution[i] - _previous_distribution[i]);
    if (d == 0.0) {
      continue;
    }

    const double ave = (_distribution[i] + _previous_distribution[i]) * 0.5;
    const double rtol = d / ave;

    if (rtol > max_rtol) {
      max_rtol = rtol;
    }

    if (rtol > _relative_tolerance) {
      failed = 1;
    }
  }

  if (! failed) {
    return 1;
  }

  if (_verbose) {
    cerr << _fingerprints_read
         << " fingerprints read, max relative tolerance " << max_rtol << '\n';
  }

  return 0;
}

int
Options::WriteDistribution() {
  uint64_t tot = 0;
  for (int i = 0; i < kNBins; ++i) {
    tot += _acc[i];
  }

  return WriteDistribution(tot);
}

int
Options::WriteDistribution(uint64_t tot) {
  IWString fname;
  fname << _stem << _next_fname_index << ".csv";
  cerr << "WriteDistribution '" << fname << "'\n";
  ++_next_fname_index;
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "Options::WriteDistribution:cannot open '" << fname << "'\n";
    return 0;
  }

  static constexpr char kSep = ',';

  Accumulator<double> stats;

  output << "Dist" << kSep << "Fraction\n";
  uint64_t sum = 0;
  for (int i = 0; i < kNBins; ++i) {
    float dist = static_cast<float>(i / 1000.0f);
    output << dist << kSep << static_cast<float>(_distribution[i]) << '\n';

    sum += _acc[i];

    if (_acc[i] > 0) {
      stats.extra(dist, _acc[i]);
    }

    if (dist > _max_distance) {
      break;
    }
  }

  if (_verbose) {
    cerr << "Average " << static_cast<float>(stats.average()) << '\n';
  }

  return 1;
}

void
Options::AddToBinnedDistances(float d) {
  int ndx = static_cast<int>(d * 1000.0f + 0.49999f);
  ++_acc[ndx];
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _fingerprints_read << " fingerprints\n";
  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vs:S:p:a:r:n:f:h:T:b:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (! cl.option_present('S')) {
    cerr << "Must specify output file name stem via the -S option\n";
    Usage(1);
  }

  if (! cl.option_present('p')) {
    cerr << "Must specify comparison set via the -p option\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Must specify input file as command line argument\n";
    Usage(1);
  }

  if (cl.number_elements() > 1) {
    cerr << "Only one input query fingerprint file can be specified\n";
    Usage(1);
  }

#ifdef USE_OMP
  if (cl.option_present('h')) {
    int h;
    if (!cl.value('h', h) || h <= 0) {
      cerr << "The maximum number of threads to use (-h) must be a valid whole +ve "
              "number\n";
      Usage(2);
    }

    omp_set_num_threads(h);
  }
#endif

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  int converged = 0;
  if (! options.Process(cl[0], converged)) {
    cerr << "Error processing '" << cl[0] << "'\n";
    return 1;
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace gfp_compare_collections

int
main(int argc, char **argv) {
  int rc = gfp_compare_collections::Main(argc, argv);

  return rc;
}
