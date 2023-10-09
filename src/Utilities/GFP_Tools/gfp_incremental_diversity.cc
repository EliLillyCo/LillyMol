/*
  Determines the incremental diversity added across a series of molecules
*/

#include <stdlib.h>

#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Utilities/GFP_Tools/gfp.h"

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

static char output_separator = ' ';

static int flush_output = 0;

static Accumulator<double> acc;

static Report_Progress report_progress;

static int window_present = 0;

static int rescan_window_fail = 0;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Computes incremental diversity as molecules are added to a collection\n";
  cerr << " -s <n>         number of fingerprints in input\n";
  cerr << " -b <n>         only consider windows of <n> previous molecules when assessing diversity\n";
  cerr << " -f             flush output after each molecule processed\n";
  cerr << " -r <n>         report progress every <n> molecules processed\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static void
compare_no_window(IW_General_Fingerprint* pool, const int ndx, const int lookback,
                  int& id_maxsim, float& maxsim)
{
  int jstart = 0;
  if (lookback > 0 && lookback < ndx) {
    jstart = ndx - lookback;
  }

  for (int j = jstart; j < ndx; ++j) {
    const float s = pool[ndx].tanimoto(pool[j]);
    if (s > maxsim) {
      maxsim = s;
      id_maxsim = j;
    }
  }

  return;
}

static void
compare_with_window(IW_General_Fingerprint* pool, const int ndx, const int lookback,
                    int& id_maxsim, float& maxsim)
{
  int jstart = 0;
  if (lookback > 0 && lookback < ndx) {
    jstart = ndx - lookback;
  }

  for (int j = jstart; j < ndx; ++j) {
    if (!can_be_compared(pool[ndx], pool[j])) {
      continue;
    }

    const float s = pool[ndx].tanimoto(pool[j]);
    if (s > maxsim) {
      maxsim = s;
      id_maxsim = j;
    }
  }

  return;
}

static int
gfp_incremental_diversity(IW_General_Fingerprint* pool, const int istart,
                          const int pool_size, const int lookback,
                          IWString_and_File_Descriptor& output)
{
  const time_t tzero = time(NULL);
  for (int i = istart; i < pool_size; ++i) {
    float maxsim = 0.0f;
    int id_maxsim = -1;

    if (window_present) {
      compare_with_window(pool, i, lookback, id_maxsim, maxsim);
      if (id_maxsim < 0) {
        rescan_window_fail++;
        compare_no_window(pool, i, lookback, id_maxsim, maxsim);
      }
    } else {
      compare_no_window(pool, i, lookback, id_maxsim, maxsim);
    }

    output << pool[i].id() << output_separator << (1.0f - maxsim) << output_separator
           << pool[id_maxsim].id() << output_separator << (i - id_maxsim) << '\n';

    if (verbose) {
      acc.extra(1.0f - maxsim);
    }

    if (flush_output) {
      output.flush();
    } else {
      output.write_if_buffer_holds_more_than(8192);
    }

    if (report_progress()) {
      const time_t tnow = time(NULL);
      cerr << "Processed " << i << " fingerprints, elapsed " << (tnow - tzero) << '\n';
    }
  }

  return 1;
}

static int
build_pool(iwstring_data_source& input, IW_General_Fingerprint*& pool, int& pool_size)
{
  if (0 == pool_size) {
    pool_size = input.count_records_starting_with("$SMI<");

    if (pool_size < 2) {
      cerr << "Not enough fngerprints in input\n";
      return 0;
    }
  }

  pool = new IW_General_Fingerprint[pool_size];

  IW_TDT tdt;
  int fatal = 0;

  int ndx = 0;
  while (tdt.next(input) && ndx < pool_size) {
    if (!pool[ndx].construct_from_tdt(tdt, fatal)) {
      return 0;
    }
    ndx++;
  }

  pool_size = ndx;

  return pool_size;
}

static int
build_pool(const char* fname, IW_General_Fingerprint*& pool, int& pool_size)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return build_pool(input, pool, pool_size);
}

static int
gfp_incremental_diversity(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vF:b:s:fr:W:a:o:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (need_to_call_initialise_fingerprints(cl)) {
    if (!initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise GFP options\n";
      usage(23);
    }
  }

  int lookback = 0;

  if (cl.option_present('b')) {
    if (!cl.value('b', lookback) || lookback < 1) {
      cerr << "The lookback value (-b) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will only look at molecules " << lookback << " from the current\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise progress reporter (-r)\n";
      usage(1);
    }
  }

  window_present = cl.option_present('W');

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IW_General_Fingerprint* pool = nullptr;

  int pool_size = 0;

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 2) {
      cerr << "The pool size option (-s) must be a valid pool size\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Input set assumed to have " << pool_size << " fingerprints\n";
    }
  }

  if (!build_pool(cl[0], pool, pool_size)) {
    cerr << "Cannot build pool '" << cl[0] << "'\n";
    return 1;
  }

  int istart = 1;

  if (cl.option_present('a')) {
    if (!cl.value('a', istart) || istart < 1 || istart > pool_size) {
      cerr << "The start scanning value (-a) must be a whole +ve number >= 1 and <= "
           << pool_size << '\n';
      return 1;
    }

    if (verbose) {
      cerr << "Will start scanning from " << istart << '\n';
    }
  }

  if (cl.option_present('o')) {
    int o;
    if (!cl.value('o', o) || o < istart || o > pool_size) {
      cerr << "The stop value (-o) just be a +ve whole number >= " << istart
           << " and <= " << pool_size << '\n';
      return 1;
    }

    pool_size = o;

    if (verbose) {
      cerr << "Will perform comparisons to " << pool_size << '\n';
    }
  }

  IWString_and_File_Descriptor output(1);

  output << "ID" << output_separator << "Dist" << output_separator << "Prev"
         << output_separator << "Sep\n";

  gfp_incremental_diversity(pool, istart, pool_size, lookback, output);

  output.flush();

  if (verbose) {
    if (window_present) {
      cerr << rescan_window_fail << " rescans due to window fail\n";
    }
    cerr << acc.n() << " incremental distances btw " << acc.minval() << " and "
         << acc.maxval() << " ave " << static_cast<float>(acc.average()) << '\n';
  }

  return 0;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = gfp_incremental_diversity(argc, argv);

  return rc;
}
