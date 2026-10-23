// Incrementally retain fingerprints that are sufficiently distant from all
// previously retained fingerprints.

#include <atomic>
#include <cstdint>
#include <iostream>
#include <memory>

#include "omp.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iwstring/string_with_commas.h"

#include "Utilities/GFP_Tools/gfp.h"

namespace gfp_diverse_subset {

using std::cerr;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  cerr << R"(Incrementally selects an input-order-dependent diverse fingerprint subset.
For each fingerprint, discard it if its distance from any retained fingerprint
is less than the threshold. Output is SMILES and identifier.

Note that if reading from standard input you must fully specify the fingerprints
being used. By default, most GFP tools peek at the input file in order to determine
what fingerprints are present. That does not work with standard input, so fully
specify the fingerprints to be used.
... generate fingerprint | gfp_diverse_subset -F FPIW -F FPMK -F FPMK2 -P MPR ... -

 -t <dist>      minimum distance from every retained fingerprint (required)
 -s <n>         maximum number of fingerprints to retain
 -x <n>         use OpenMP when <n> fingerprints have been retained (default 1000)
 -h <n>         number of OpenMP threads to use
 -r <n>         report progress every <n> fingerprints read
 -F ...         standard fingerprint options, enter '-F help'
 -v             verbose output
)";
// clang-format on

  ::exit(rc);
}

class Options {
 private:
  int _verbose = 0;
  similarity_type_t _threshold = 0.0f;
  int _parallel_crossover = 1000;
  int _number_threads = 1;
  uint32_t _maximum_retained = 0;

  IWString _smiles_tag = "$SMI<";
  IWString _identifier_tag = "PCN<";

  resizable_array_p<IW_General_Fingerprint> _retained;

  Report_Progress _report_progress;

  uint64_t _fingerprints_read = 0;
  uint64_t _fingerprints_retained = 0;
  uint64_t _fingerprints_discarded = 0;
  uint64_t _comparisons = 0;

  // Private functions.

  bool TooCloseSerial(IW_General_Fingerprint& fp);
  bool TooCloseParallel(IW_General_Fingerprint& fp);
  bool TooClose(IW_General_Fingerprint& fp);

  int Process(iwstring_data_source& input,
              IWString_and_File_Descriptor& output);

  void MaybeReportProgress();

 public:
  int Initialise(Command_Line& cl);

  int Process(const char* fname, IWString_and_File_Descriptor& output);

  bool Full() const {
    return _maximum_retained > 0 &&
           _retained.size() >= _maximum_retained;
  }

  int Report(std::ostream& output) const;
};

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (! cl.option_present('t') ||
      ! cl.value('t', _threshold) ||
      _threshold < 0.0f || _threshold > 1.0f) {
    cerr << "The distance threshold (-t) must be between 0 and 1\n";
    return 0;
  }

  if (cl.option_present('x')) {
    if (! cl.value('x', _parallel_crossover) || _parallel_crossover < 1) {
      cerr << "The parallel crossover (-x) must be a positive integer\n";
      return 0;
    }
  }

  if (cl.option_present('h')) {
    if (! cl.value('h', _number_threads) || _number_threads < 1) {
      cerr << "The number of OpenMP threads (-h) must be a positive integer\n";
      return 0;
    }
    omp_set_num_threads(_number_threads);
  }

  if (cl.option_present('s')) {
    if (! cl.value('s', _maximum_retained) || _maximum_retained < 1) {
      cerr << "The maximum retained size (-s) must be a positive integer\n";
      return 0;
    }
    _retained.reserve(_maximum_retained);
  } else {
    _retained.reserve(100000);
  }

  if (cl.option_present('r')) {
    if (! _report_progress.initialise(cl, 'r', _verbose)) {
      cerr << "Cannot initialise progress reporting (-r)\n";
      return 0;
    }
  }

  if (_verbose) {
    cerr << "Will retain fingerprints at least " << _threshold
         << " from every previously retained fingerprint\n";
    cerr << "Will use OpenMP when " << _parallel_crossover
         << " fingerprints have been retained\n";
    if (_number_threads > 0) {
      cerr << "Will use " << _number_threads << " OpenMP threads\n";
    }
    if (_maximum_retained > 0) {
      cerr << "Will retain at most " << _maximum_retained << " fingerprints\n";
    }
  }

  return 1;
}

bool
Options::TooCloseSerial(IW_General_Fingerprint& fp) {
  const int n = _retained.number_elements();
  for (int i = 0; i < n; ++i) {
    if (! can_be_compared(fp, *_retained[i])) {
      continue;
    }

    ++_comparisons;
    if (fp.distance(*_retained[i]) < _threshold) {
      return true;
    }
  }

  return false;
}

bool
Options::TooCloseParallel(IW_General_Fingerprint& fp) {
  const int n = _retained.number_elements();
//std::atomic<bool> too_close(false);
  bool too_close = false;

#pragma omp parallel for schedule(static) reduction(|| : too_close)
  for (int i = 0; i < n; ++i) {
//  if (too_close.load(std::memory_order_relaxed)) {
//    continue;
//  }

    if (! can_be_compared(fp, *_retained[i])) {
      continue;
    }

//  ++comparisons;
    if (fp.distance(*_retained[i]) < _threshold) {
      too_close = true;
//    too_close.store(true, std::memory_order_relaxed);
    }
  }

//_comparisons += comparisons;
//return too_close.load(std::memory_order_relaxed);
  return too_close;
}

bool
Options::TooClose(IW_General_Fingerprint& fp) {
  if (_retained.number_elements() < _parallel_crossover ||
      _number_threads == 1) {
    return TooCloseSerial(fp);
  }

  return TooCloseParallel(fp);
}

int
Options::Process(iwstring_data_source& input,
                 IWString_and_File_Descriptor& output) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    const_IWSubstring smiles;
    if (! tdt.dataitem_value(_smiles_tag, smiles)) {
      cerr << "Options::Process:no smiles '" << _smiles_tag << "' in TDT\n";
      cerr << tdt;
      return 0;
    }

    const_IWSubstring identifier;
    if (! tdt.dataitem_value(_identifier_tag, identifier)) {
      cerr << "Options::Process:no identifier '" << _identifier_tag << "' in TDT\n";
      cerr << tdt;
      return 0;
    }

    auto fp = std::make_unique<IW_General_Fingerprint>();
    int fatal = 0;
    if (! fp->construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        cerr << "Options::Process:cannot parse fingerprint\n";
        return 0;
      }
      continue;
    }

    ++_fingerprints_read;
    MaybeReportProgress();
    if (TooClose(*fp)) {
      ++_fingerprints_discarded;
      continue;
    }

    output << smiles << ' ' << identifier << '\n';
    // The calculation slows down over time so give more feedback on progress.
    if (_retained.size() > 500000) {  // an arbitrarily chosen number.
      output.flush();
    } else {
      output.write_if_buffer_holds_more_than(8192);
    }

    _retained << fp.release();
    ++_fingerprints_retained;
    if (Full()) {
      return output.good();
    }
  }

  return output.good();
}

void
Options::MaybeReportProgress() {
  if (! _report_progress()) {
    return;
  }

  cerr << "Read " << iwstring::with_commas(_fingerprints_read) << " fingerprints, retained "
       << iwstring::with_commas(_fingerprints_retained) << ", discarded "
       << _fingerprints_discarded << ", comparisons " << _comparisons << '\n';
}

int
Options::Process(const char* fname,
                 IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::Process:cannot open '" << fname << "'\n";
    return 0;
  }

  return Process(input, output);
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _fingerprints_read << " fingerprints, retained "
         << _fingerprints_retained << ", discarded "
         << _fingerprints_discarded << '\n';
//output << _comparisons << " fingerprint comparisons performed\n";
  cerr << "Retained " << _retained.size() << " fingerprints\n";

  return output.good();
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vF:P:t:s:x:h:r:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');
  if (need_to_call_initialise_fingerprints(cl)) {
    if (! initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise fingerprints\n";
      return 1;
    }
  } else if (! initialise_fingerprints(cl[0], verbose)) {
    cerr << "Cannot initialise fingerprints from '" << cl[0] << "'\n";
    return 1;
  }

  Options options;
  if (! options.Initialise(cl)) {
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (! options.Process(fname, output)) {
      return 1;
    }
    if (options.Full()) {
      break;
    }
  }
  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace gfp_diverse_subset

int
main(int argc, char** argv) {
  return gfp_diverse_subset::Main(argc, argv);
}
