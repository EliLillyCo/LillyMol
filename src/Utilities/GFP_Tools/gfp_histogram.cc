/* We want to know the distribution of distances in a set of molecules
 */

#include <math.h>
#include <stdlib.h>

#include <iomanip>
#include <memory>
#include <random>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/histogram/iwhistogram.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/primes.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/tversky.h"

using std::cerr;
using std::endl;

static int verbose = 0;

static IW_General_Fingerprint* pool = nullptr;

static int pool_size = 0;

static int do_intra_library_histogram = 0;

static int nfingerprints = 0;

static int output_distances = 1;

static int report_progress_intra = 0;
static int report_progress_inter = 0;

static resizable_array<similarity_type_t> report_number_shorter_than;

static Tversky tversky;

static int write_bare_data = 0;

static int write_bare_data_as_fractions = 0;

static int only_write_to_last_non_zero_value = 0;

static float fraction_to_sample = static_cast<float>(0.0);

static float inter_collection_fraction_to_sample = static_cast<float>(0.0);

static double convergence_criterion = 0.0;

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
  cerr << "Computes the distance matrix for a pool of fingerprints\n";
  cerr << " -m               output as similarities rather than distances\n";
  cerr << " -I               do intra-library comparisons\n";
  cerr << " -b               write minimal data - should go directly into a plotting programme\n";
  cerr << " -f               with the -b option, write data as fractions rather than counts\n";
  cerr << " -n               nearest neighbour distances only\n";
  cerr << " -t <dist>        report the number of distances shorter than <dist>\n";
  cerr << " -r <number>      report intra-library progress every <number> molecules\n";
  cerr << " -R <number>      report inter-library progress every <number> molecules\n";
  cerr << " -s <size>        specify the number of molecules\n";
  cerr << " -p <file>        optional pool file for inter-collection distributions\n";
  cerr << " -H <min,max,dx>  histogram initialisation conditions (default 0.0,1.0,0.01)\n";
  cerr << " -F, -P, ...      standard fingerprint options, enter -F help for info\n";
  cerr << " -V ...           Tversky specifications\n";
  cerr << " -M samp=<frac>   fraction of the set to sample\n";
  cerr << " -M convg=<frac>  convergence criterion for random sampling\n";
  cerr << " -v               verbose output\n";
  // clang-format on

  exit(rc);
}

static int
do_write_bare_data(const IWHistogram& histogram, std::ostream& output)
{
  const unsigned int* r = histogram.raw_counts();

  float multiplier;
  if (write_bare_data_as_fractions) {
    multiplier = static_cast<float>(1.0) / static_cast<float>(histogram.nsamples());
  } else {
    multiplier = static_cast<float>(1.0);
  }

  int nb = histogram.nbuckets();
  if (only_write_to_last_non_zero_value) {
    nb = histogram.highest_filled_bucket();
  }

  for (int i = 0; i < nb; i++) {
    unsigned int c = r[i];

    float x = histogram.minval() + static_cast<float>(i) * histogram.delta();

    if (write_bare_data_as_fractions) {
      output << x << ' ' << static_cast<float>(c) * multiplier << endl;
    } else {
      output << x << ' ' << c << endl;
    }
  }

  return output.good();
}

static int
allocate_pool()
{
  assert(pool_size > 0 && NULL == pool);

  pool = new IW_General_Fingerprint[pool_size];

  if (verbose) {
    cerr << "Pool sized for " << pool_size << " molecules\n";
  }

  assert(NULL != pool);

  return 1;
}

static int
compute_the_distance(IW_General_Fingerprint& fp1, IW_General_Fingerprint& fp2,
                     similarity_type_t& result)
{
  if (!can_be_compared(fp1, fp2)) {
    return 0;
  }

  if (tversky.active()) {
    result = fp1.tversky(fp2, tversky);
  } else {
    result = fp1.tanimoto(&fp2);
  }

  if (output_distances) {
    result = static_cast<similarity_type_t>(1.0) - result;
  }

  return 1;
}

static void
fill_inter_collection_histogram_randomly(IW_General_Fingerprint& fp,
                                         IWHistogram& histogram,
                                         Accumulator<similarity_type_t>& stats)
{
  int istep = static_cast<int>(nfingerprints * inter_collection_fraction_to_sample) + 1;

  for (int i = 0; i < nfingerprints; i += istep) {
    similarity_type_t result;

    if (!compute_the_distance(fp, pool[i], result)) {
      continue;
    }

    stats.extra(result);
    histogram.extra(result);
  }

  return;
}

static void
fill_inter_collection_histogram(IW_General_Fingerprint& fp, IWHistogram& histogram,
                                Accumulator<similarity_type_t>& stats)
{
  for (int i = 0; i < pool_size; i++) {
    similarity_type_t result;

    if (!compute_the_distance(fp, pool[i], result)) {
      continue;
    }

    stats.extra(result);
    histogram.extra(result);
  }

  return;
}

static int inter_collection_fingerprints_read = 0;

static int
fill_inter_collection_histogram(iwstring_data_source& input, IWHistogram& histogram,
                                Accumulator<similarity_type_t>& stats)
{
  IW_TDT tdt;
  while (tdt.next(input)) {
    IW_General_Fingerprint fp;
    int fatal;
    if (!fp.construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        cerr << "Cannot read inter-collection TDT, line " << input.lines_read() << endl;
        return 0;
      }

      continue;
    }

    inter_collection_fingerprints_read++;

    if (report_progress_inter &&
        0 == inter_collection_fingerprints_read % report_progress_inter) {
      cerr << "Read " << inter_collection_fingerprints_read
           << " inter collection fingerprints\n";
    }

    if (static_cast<float>(0.0) != inter_collection_fraction_to_sample) {
      fill_inter_collection_histogram_randomly(fp, histogram, stats);
    } else {
      fill_inter_collection_histogram(fp, histogram, stats);
    }
  }

  return 1;
}

static int
fill_inter_collection_histogram(const const_IWSubstring& fname, IWHistogram& histogram,
                                Accumulator<similarity_type_t>& stats)
{
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open inter colleciton file '" << fname << "'\n";
    return 0;
  }

  return fill_inter_collection_histogram(input, histogram, stats);
}

static int
fill_intra_pool_histogram_nearest_neighbour(IWHistogram& histogram,
                                            Accumulator<similarity_type_t>& stats)
{
  if (verbose) {
    cerr << "fill_intra_pool_histogram_nearest_neighbour " << nfingerprints << endl;
  }

  similarity_type_t* minval = new float[nfingerprints];
  std::unique_ptr<similarity_type_t> free_minval(minval);

  set_vector(minval, nfingerprints, static_cast<similarity_type_t>(1.0));

  for (int i = 0; i < nfingerprints; i++) {
    if (report_progress_intra && i > 0 && 0 == i % report_progress_intra) {
      cerr << "Processed " << i << " fingerprints\n";
    }

    IW_General_Fingerprint& fpi = pool[i];

    for (int j = i + 1; j < nfingerprints; j++) {
      similarity_type_t result;

      if (!compute_the_distance(fpi, pool[j], result)) {
        continue;
      }

      //    cerr << "Between " << i << " and " << j << " value " << result << endl;

      if (result < minval[i]) {
        minval[i] = result;
      }
      if (result < minval[j]) {
        minval[j] = result;
      }
    }
  }

  for (int i = 0; i < nfingerprints; i++) {
    stats.extra(minval[i]);
    histogram.extra(minval[i]);
  }

  return 1;
}

static int
is_converged(const IWHistogram& histogram, double* current_value)
{
  int nsamples = histogram.nsamples();

  double multiplier = 1.0 / static_cast<double>(nsamples);

  const unsigned int* raw_count = histogram.raw_counts();

  int unconverged_buckets = 0;

#ifdef DEBUG_IS_CONVERGED
  cerr << "Checking for convergence, nsamples " << histogram.nsamples() << endl;
#endif

  double largest_difference = 0.0;

  for (int i = 0; i < histogram.nbuckets(); i++) {
    double new_value = static_cast<double>(raw_count[i]) * multiplier;

#ifdef DEBUG_IS_CONVERGED
    cerr << "i = " << i << " compare " << new_value << " with previous "
         << current_value[i] << endl;
#endif

    if (fabs(new_value - current_value[i]) > convergence_criterion) {
      unconverged_buckets++;
    }

    if (verbose > 2) {
      if (fabs(new_value - current_value[i]) > largest_difference) {
        largest_difference = fabs(new_value - current_value[i]);
      }
    }

    current_value[i] = new_value;
  }

  if (verbose > 2) {
    double total_computations_to_do =
        static_cast<double>(nfingerprints) * static_cast<double>(nfingerprints - 1) / 2.0;
    double percent_done =
        static_cast<double>(histogram.nsamples()) * 200.0 / total_computations_to_do;

    cerr << histogram.nsamples() << " samples, " << static_cast<float>(percent_done)
         << "%, " << unconverged_buckets << " uncoverged buckets, largest "
         << largest_difference << endl;
  }

  return 0 == unconverged_buckets;
}

static int
find_prime_index(int nfingerprints, int& idelta, int& jdelta)
{
  unsigned int tmp = static_cast<unsigned int>(nfingerprints * fraction_to_sample);

  idelta = 0;

  for (int i = 1; i < IWNPRIMES; i++) {
    if (primes[i] < tmp) {
      continue;
    }

    if (0 == primes[i] % nfingerprints) {
      continue;
    }

    idelta = primes[i];
    break;
  }

  if (0 == idelta) {
    idelta = 1;
  }

  for (int i = 0; i < IWNPRIMES; i++) {
    if (primes[i] <= static_cast<unsigned int>(idelta)) {
      continue;
    }

    if (0 == primes[i] % nfingerprints) {
      continue;
    }

    jdelta = primes[i];
    return 1;
  }

  jdelta = 1;

  return 1;
}

static int
fill_intra_pool_histogram_till_converged(IWHistogram& histogram,
                                         Accumulator<similarity_type_t>& stats)
{
  double* current_value = new double[histogram.nbuckets()];
  std::unique_ptr<double> free_current_value(current_value);

  set_vector(current_value, histogram.nbuckets(), static_cast<double>(0.0));

  int idelta, jdelta;

  find_prime_index(nfingerprints, idelta, jdelta);

  if (verbose) {
    cerr << "For pool of size " << nfingerprints << " idelta " << idelta << " and jdelta "
         << jdelta << endl;
  }

  int* jstart = new int[nfingerprints];
  std::unique_ptr<int> free_jstart(jstart);

  for (int i = 0; i < nfingerprints; i++) {
    jstart[i] = i + 1;
  }
  jstart[nfingerprints - 1] = 0;

  int itimes = nfingerprints / jdelta + 1;

  for (int outer_loop = 0; outer_loop < itimes; outer_loop++) {
    //  cerr << "outer_loop " << outer_loop << endl;
    int i = 0;
    while (1) {
      IW_General_Fingerprint& fpi = pool[i];
      //    cerr << " i = " << i << " starting with " << jstart[i] << endl;
      for (int j = 0; j < jdelta; j++) {
        if (i < jstart[i]) {
          similarity_type_t d;
          if (compute_the_distance(fpi, pool[jstart[i]], d)) {
            histogram.extra(d);
            stats.extra(d);
            //          cerr << "Distance between " << i << " and " << jstart[i] << " is "
            //          << d << endl;
          }
        }

        jstart[i] += jdelta;
        //      cerr << "jstart incremented to " << jstart[i] << endl;
        if (jstart[i] >= nfingerprints) {
          jstart[i] = jstart[i] % nfingerprints;
        }

        if (i + 1 == jstart[i]) {  // back where we started
          break;
        }
      }

      i += idelta;
      if (i == nfingerprints) {
        break;
      } else if (i > nfingerprints) {
        i = i % nfingerprints;
      }
    }

    if (is_converged(histogram, current_value)) {
      if (verbose) {
        cerr << "Converged " << convergence_criterion << " after " << histogram.nsamples()
             << " samples "
             << static_cast<float>(histogram.nsamples() * 100) /
                    static_cast<float>(nfingerprints * (nfingerprints - 1)) * 2.0
             << "%\n";
      }
      return 1;
    }
  }

  cerr << "Full evaluation\n";
  return 1;
}

static int
fill_intra_pool_histogram_randomly(IWHistogram& histogram,
                                   Accumulator<similarity_type_t>& stats)
{
  int* idone = new_int(nfingerprints);
  std::unique_ptr<int> free_idone(idone);

  int fingerprints_sampled = 0;

  int fingerprints_to_sample = static_cast<int>(nfingerprints * fraction_to_sample) + 1;

  int jdelta = nfingerprints / fingerprints_to_sample;

  if (0 == jdelta) {  // very hard to imagine
    jdelta = 1;
  }

  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> unfingerprints(0, nfingerprints - 1);
  while (fingerprints_sampled < fingerprints_to_sample) {
    const int i = unfingerprints(rng);

    if (idone[i]) {
      continue;
    }

    fingerprints_sampled++;

    if (report_progress_intra && i > 0 &&
        0 == fingerprints_sampled % report_progress_intra) {
      cerr << "Processed " << fingerprints_sampled << " fingerprints\n";
    }

    IW_General_Fingerprint& fpi = pool[i];

    std::uniform_int_distribution<int> u(0, jdelta);
    int jstart = u(rng);

    if (0 == (i - jstart) % jdelta) {  // avoid comparing I with J
      jstart++;
    }

    for (int j = jstart; j < nfingerprints; j += jdelta) {
      similarity_type_t result;

      if (!compute_the_distance(fpi, pool[j], result)) {
        continue;
      }

      stats.extra(result);
      histogram.extra(result);
    }
  }

  return 1;
}

static int
fill_intra_pool_histogram(IWHistogram& histogram, Accumulator<similarity_type_t>& stats)
{
  for (int i = 0; i < nfingerprints; i++) {
    if (report_progress_intra && i > 0 && 0 == i % report_progress_intra) {
      cerr << "Processed " << i << " fingerprints\n";
    }

    IW_General_Fingerprint& fpi = pool[i];

    for (int j = i + 1; j < nfingerprints; j++) {
      similarity_type_t result;

      if (!compute_the_distance(fpi, pool[j], result)) {
        continue;
      }

      stats.extra(result);
      histogram.extra(result);
    }
  }

  return 1;
}

static int
build_pool(IW_TDT& tdt)
{
  IW_General_Fingerprint& fp = pool[nfingerprints];
  int fatal;

  if (!fp.construct_from_tdt(tdt, fatal)) {
    if (!fatal) {
      return 1;
    }

    return 0;
  }

  nfingerprints++;

  return 1;
}

static int
build_pool(iwstring_data_source& input)
{
  int tdts_read = 0;

  IW_TDT tdt;
  while (tdt.next(input)) {
    tdts_read++;

    if (!build_pool(tdt)) {
      return 0;
    }

    if (nfingerprints >= pool_size) {
      cerr << "Pool is full, " << nfingerprints << endl;
      break;
    }
  }

  if (verbose) {
    cerr << "Read " << tdts_read << " TDT's, pool contains " << nfingerprints
         << " fingerprints\n";
  }

  return 1;
}

static int
build_pool(const char* fname)
{
  iwstring_data_source input(fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (0 == pool_size) {
    std::unique_ptr<re2::RE2> pcn = std::make_unique<RE2>("^PCN<");

    pool_size = input.grep(*pcn);

    if (0 == pool_size) {
      cerr << "zero occurrences of '" << pcn->pattern() << "' in input\n";
      return 0;
    }

    if (!allocate_pool()) {
      return 0;
    }
  }

  return build_pool(input);
}

static int
distance_histogram(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vs:V::F:P:W:Q:mO:r:R:p:IH:bft:nM:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!initialise_fingerprints(cl, verbose)) {
    cerr << "Cannot initialise fingerprint specifications\n";
    usage(23);
  }

  if (cl.option_present('V')) {
    if (!tversky.parse_command_line(cl, 'V', verbose)) {
      cerr << "Cannot get Tversky A specifications\n";
      usage(18);
    }
  }

  if (cl.option_present('m')) {
    output_distances = 0;
    if (verbose) {
      cerr << "Will output similarity values rather than distances\n";
    }
  }

  int nearest_neighbour_distances_only = 0;  // by default, we do all

  if (cl.option_present('n')) {
    nearest_neighbour_distances_only = 1;

    if (verbose) {
      cerr << "Will report only the nearest neighbour distance of each molecule\n";
    }
  }

  if (cl.option_present('I')) {
    do_intra_library_histogram = 1;
    if (verbose) {
      cerr << "Will compute intra-library histogram\n";
    }
  } else if (!cl.option_present('p')) {
    do_intra_library_histogram = 1;
    cerr << "Will compute intra-library histogram\n";
  }

  if (cl.option_present('b')) {
    write_bare_data = 1;

    if (verbose) {
      cerr << "Just barebones data will be written\n";
    }

    if (cl.option_present('f')) {
      write_bare_data_as_fractions = 1;

      if (verbose) {
        cerr << "Barebones data written as fractions\n";
      }
    }
  }

  if (cl.option_present('t')) {
    int i = 0;
    const_IWSubstring t;
    while (cl.value('t', t, i++)) {
      similarity_type_t tmp;
      if (!t.numeric_value(tmp) || tmp < 0.0 || tmp > 1.0) {
        cerr << "The report number shorter than (-t) option must be followed by a valid "
                "distance\n";
        usage(4);
      }

      report_number_shorter_than.add(tmp);

      if (verbose) {
        cerr << "Will report the number of histogram samples closer than " << tmp << endl;
      }
    }
  }

  if (cl.option_present('r')) {
    if (!cl.option_present('I')) {
      cerr << "The report progress intra (-r) option only makes sense with the -I "
              "option\n";
      usage(32);
    }

    if (!cl.value('r', report_progress_intra) || report_progress_intra < 1) {
      cerr << "The report progress (intra) (-r) option must be followed by a positive "
              "whole number\n";
      usage(11);
    }

    if (verbose) {
      cerr << "Will report intra-collection progress every " << report_progress_intra
           << " fingerprints\n";
    }
  }

  if (cl.option_present('R')) {
    if (!cl.option_present('p')) {
      cerr << "The report progress (inter) option only makes sense with the -p option\n";
      usage(14);
    }

    if (!cl.value('R', report_progress_inter) || report_progress_inter < 1) {
      cerr << "The report progress (inter) (-R) option must be followed by a positive "
              "whole number\n";
      usage(11);
    }

    if (verbose) {
      cerr << "Will report inter-collection progress every " << report_progress_inter
           << " fingerprints\n";
    }
  }

  if (cl.option_present('M')) {
    const_IWSubstring m;
    int i = 0;
    while (cl.value('M', m, i++)) {
      if (m.starts_with("convg=")) {
        m.remove_leading_chars(6);
        if (!m.numeric_value(convergence_criterion) || convergence_criterion <= 0.0 ||
            convergence_criterion >= 1.0) {
          cerr << "Invalid convergence criterion '" << m << "'\n";
          usage(4);
        }

        if (verbose) {
          cerr << "Iterative convergence criterion " << convergence_criterion << endl;
        }
      } else if (m.starts_with("samp=")) {
        m.remove_leading_chars(5);
        if (!m.numeric_value(fraction_to_sample) || fraction_to_sample <= 0.0 ||
            fraction_to_sample >= 1.0) {
          cerr << "INvalid fraction to sample '" << m << "'\n";
          usage(4);
        }

        if (verbose) {
          cerr << "Will scan the pool in " << fraction_to_sample << " increments\n";
        }
      } else {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        usage(4);
      }
    }
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 2) {
      cerr << "The -s option must be followed by a whole number > 2\n";
      usage(5);
    }

    if (!allocate_pool()) {
      return 12;
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (1 != cl.number_elements()) {
    cerr << "Takes only one command line argument\n";
    usage(3);
  }

  for (int i = 0; i < cl.number_elements(); i++) {
    if (!build_pool(cl[i])) {
      return i + 1;
    }
  }

  if (verbose) {
    cerr << "Read " << nfingerprints << " fingerprints\n";
  }

  if (0 == nfingerprints) {
    cerr << "No fingerprints read\n";
    return 32;
  }

  if (do_intra_library_histogram) {
    Accumulator<similarity_type_t> stats;

    IWHistogram histogram;

    if (cl.option_present('H')) {
      const_IWSubstring h = cl.string_value('H');

      if (!histogram.initialise(h)) {
        cerr << "Invalid histogram initialisation '" << h << "'\n";
        usage(5);
      }

      if (verbose) {
        cerr << "Histogram initialised ";
        histogram.debug_print(cerr);
      }
    } else {
      (void)histogram.initialise(0.0, 1.0, 0.01);
    }

    if (nearest_neighbour_distances_only) {
      fill_intra_pool_histogram_nearest_neighbour(histogram, stats);
    } else if (convergence_criterion > 0.0) {
      fill_intra_pool_histogram_till_converged(histogram, stats);
    } else if (fraction_to_sample > static_cast<float>(0.0)) {
      fill_intra_pool_histogram_randomly(histogram, stats);
    } else {
      fill_intra_pool_histogram(histogram, stats);
    }

    if (write_bare_data) {
      do_write_bare_data(histogram, std::cout);
    } else {
      std::cout << "Intra collection distribution\n";
      histogram.write(std::cout);
    }

    if (report_number_shorter_than.number_elements() && histogram.nsamples()) {
      std::cout.precision(2);
      for (int i = 0; i < report_number_shorter_than.number_elements(); i++) {
        int n = histogram.number_samples_less_than(report_number_shorter_than[i]);
        float fraction = static_cast<float>(n) / static_cast<float>(histogram.nsamples());
        std::cout << n << " intra collection samples (fraction " << fraction
                  << " ) less than " << report_number_shorter_than[i] << endl;
      }
    }

    if (verbose) {
      cerr << "Computed " << stats.n() << " values, between " << stats.minval() << " and "
           << stats.maxval() << endl;
      if (stats.n() > 1) {
        cerr << "Average " << stats.average() << " variance " << stats.variance() << endl;
      }
    }
  }

  if (cl.option_present('p')) {
    IWHistogram inter_histogram;

    (void)inter_histogram.initialise(0.0, 1.0, 0.01);

    Accumulator<similarity_type_t> inter_stats;

    int i = 0;
    const_IWSubstring p;
    while (cl.value('p', p, i++)) {
      if (!fill_inter_collection_histogram(p, inter_histogram, inter_stats)) {
        cerr << "Cannot process inter-collection file '" << p << "'\n";
        return i;
      }
    }

    if (verbose) {
      cerr << "Read " << inter_collection_fingerprints_read
           << " inter collection fingerprints\n";
      cerr << "Computed " << inter_stats.n() << " values, between "
           << inter_stats.minval() << " and " << inter_stats.maxval() << endl;
      if (inter_stats.n() > 1) {
        cerr << "Average " << inter_stats.average() << " variance "
             << inter_stats.variance() << endl;
      }
    }

    if (write_bare_data) {
      do_write_bare_data(inter_histogram, std::cout);
    } else {
      std::cout << "Inter collection distribution\n";
      inter_histogram.write(std::cout);
    }
  }

  delete[] pool;

  return 0;
}

int
main(int argc, char** argv)
{
  int rc = distance_histogram(argc, argv);

  return rc;
}
