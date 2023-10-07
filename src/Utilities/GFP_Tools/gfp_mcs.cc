/*
  Looks for bits commonly found across the molecules in the input.
*/

#include <stdlib.h>

#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/histogram/iwhistogram.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Utilities/GFP_Tools/gfp.h"

const char* prog_name = nullptr;

using std::cerr;
using std::endl;

static int verbose = 0;

static unsigned int count_threshold = 0;

static float fractional_threshold = 0.0f;

static IWString identifier_tag("PCN<");

static int write_fingerprint_header = 1;

static Fraction_as_String fraction_as_string;

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
  cerr << "Max common bits\n";
  cerr << " -c <n>         discard bits set in <n> or fewer molecules\n";
  cerr << " -s <n>         specify number of fingerprints (optional)\n";
  cerr << " -p <ratio>     specify number of fingerprints as fraction\n";
  cerr << " -S <stem>      write each fingerprint to own file\n";
  cerr << " -w             omit fingerprint header in ouput if only one fingerprint\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static IW_General_Fingerprint* pool = nullptr;

static unsigned int pool_size = 0;

static float reciprocal_pool_size = 0.0F;

static IWHistogram histogram;

class Bit_and_Count
{
 private:
  unsigned int _bit;
  unsigned int _count;

 public:
  Bit_and_Count();

  unsigned int bit() const {
    return _bit;
  }

  unsigned int count() const {
    return _count;
  }

  void set_bit(unsigned int b) {
    _bit = b;
  }

  void set_count(unsigned int c) {
    _count = c;
  }
};

Bit_and_Count::Bit_and_Count()
{
  _bit = 0;
  _count = 0;

  return;
}

class Bit_and_Count_Comparator
{
 private:
 public:
  int operator()(const Bit_and_Count& bac1, const Bit_and_Count& bzc2) const;
};

int
Bit_and_Count_Comparator::operator()(const Bit_and_Count& bac1,
                                     const Bit_and_Count& bac2) const
{
  if (bac1.count() < bac2.count()) {
    return 1;
  }

  if (bac1.count() > bac2.count()) {
    return -1;
  }

  return 0;
}

static Bit_and_Count_Comparator bacc;

static int
write_bits(const Bit_and_Count* bac, int n, IWString_and_File_Descriptor& output)
{
  int mcs_found = 0;

  int singletons = 0;
  for (int i = 0; i < n; i++) {
    unsigned int c = bac[i].count();
    if (c == 1) {
      ++singletons;
    }

    output << bac[i].bit() << ' ' << c;

    float ratio = static_cast<float>(c) * reciprocal_pool_size;

    fraction_as_string.append_number(output, ratio);

    if (c == pool_size) {
      output << " *";
      mcs_found++;
    }
    output << "\n";

    if (histogram.active()) {
    }
  }

  if (write_fingerprint_header) {
    output << "|\n";
  }

  output.write_if_buffer_holds_more_than(32768);

  if (verbose) {
    output.flush();
    cerr << "Encountered " << n << " bits and " << singletons << " singletons "
         << iwmisc::Fraction<float>(singletons, n) << '\n';
  }

  return mcs_found;
}

static int
gfp_mcs_fp(int f, IWString_and_File_Descriptor& output)
{
  int nb = pool[0][f].nbits();

  int* count = new_int(nb);
  std::unique_ptr<int[]> free_count(count);

  for (unsigned int i = 0; i < pool_size; i++) {
    pool[i].item(f).increment_vector(count);
  }

  Bit_and_Count* bac = new Bit_and_Count[nb];
  std::unique_ptr<Bit_and_Count[]> free_bac(bac);

  int ndx = 0;

  for (int i = 0; i < nb; i++) {
    if (count[i] <= static_cast<int>(count_threshold)) {
      continue;
    }

    bac[ndx].set_bit(i);
    bac[ndx].set_count(count[i]);
    ndx++;
  }

  if (0 == ndx) {
    cerr << "None of " << nb << " bits in fingerprint " << f << " above threshold "
         << count_threshold << endl;
    return 1;
  }

  if (verbose && ndx < nb) {
    cerr << "Fixed fingerprint " << fixed_fingerprint_tag(f) << " " << ndx << " of " << nb
         << " bits above threshold\n";
  }

  if (ndx > 1) {
    iwqsort(bac, ndx, bacc);
  }

  if (write_fingerprint_header) {
    output << fixed_fingerprint_tag(f) << "\n";
  }

  int mcs_found = write_bits(bac, ndx, output);

  if (verbose) {
    cerr << "Fingerprint '" << fixed_fingerprint_tag(f) << " found " << mcs_found
         << " mcs bits\n";
  }

  return 1;
}

static int
gfp_mcs_nc(int f, IWString_and_File_Descriptor& output)
{
  IW_Hash_Map<unsigned int, unsigned int> count;

  for (unsigned int i = 0; i < pool_size; i++) {
    const Sparse_Fingerprint& s = pool[i].sparse_fingerprint(f);

    int j = 0;
    unsigned int b;
    int c;

    while (s.next_bit_set(j, b, c)) {
      count[b]++;
    }
  }

  int nb = count.size();

  if (0 == nb) {
    cerr << "Sparse fingerprint " << f << " '" << sparse_fingerprint_tag(f)
         << "' not enough bits set\n";
    return 1;
  }

  Bit_and_Count* bac = new Bit_and_Count[nb];
  std::unique_ptr<Bit_and_Count[]> free_bac(bac);

  unsigned int ndx = 0;

  for (IW_Hash_Map<unsigned int, unsigned int>::const_iterator i = count.begin();
       i != count.end(); ++i) {
    if ((*i).second <= count_threshold) {
      continue;
    }

    bac[ndx].set_bit((*i).first);
    bac[ndx].set_count((*i).second);
    ndx++;
  }

  if (0 == ndx) {
    cerr << "None of " << nb << " bits in sparse fingerprint " << f << " above threshold "
         << count_threshold << endl;
    return 1;
  }

  if (verbose && ndx < count.size()) {
    cerr << "Sparise fingerprint " << sparse_fingerprint_tag(f) << " " << ndx << " of "
         << count.size() << " bits above threshold\n";
  }

  if (ndx > 1) {
    iwqsort(bac, ndx, bacc);
  }

  if (write_fingerprint_header) {
    output << sparse_fingerprint_tag(f) << "\n";
  }

  int mcs_found = write_bits(bac, ndx, output);

  if (verbose) {
    cerr << "Fingerprint '" << sparse_fingerprint_tag(f) << " found " << mcs_found
         << " mcs bits\n";
  }

  return 1;
}

int
allocate_pool(int n, IW_General_Fingerprint*& fp)
{
  if (NULL != fp) {
    cerr << "Warning, fingerprint array already allocated!!\n";
    delete[] fp;
  }

  fp = new IW_General_Fingerprint[n];

  if (NULL == fp) {
    cerr << "Memory failure, could not allocate " << n << " fingerprints\n";
    return 0;
  }

  return n;
}

static int
read_pool(iwstring_data_source& input)
{
  unsigned int items_in_pool = 0;

  IW_TDT tdt;

  int fatal;  // scope here for efficiency

  while (tdt.next(input)) {
    if (!pool[items_in_pool].construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        cerr << "Cannot build pool item\n";
        cerr << tdt;
        return 0;
      }
      continue;
    }

    items_in_pool++;

    if (items_in_pool == pool_size) {
      break;
    }
  }

  pool_size = items_in_pool;

  reciprocal_pool_size = 1.0f / static_cast<float>(pool_size);

  return items_in_pool;
}

static int
read_pool(const char* fname)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open fingerprint file '" << fname << "'\n";
    return 0;
  }

  if (0 == pool_size) {
    pool_size = count_tdts_in_file(input, identifier_tag);
    if (0 == pool_size) {
      return 0;
    }

    if (!allocate_pool(pool_size, pool)) {
      return 0;
    }
  }

  return read_pool(input);
}

static int
gfp_mcs(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vF:P:s:H:c:wS:p:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('H')) {
    const_IWSubstring h = cl.string_value('H');

    if (!histogram.initialise(h)) {
      cerr << "Cannot initialise histogram from '" << h << "'\n";
      return 3;
    }
  }

  if (cl.option_present('c') && cl.option_present('p')) {
    cerr << "The -c and -p options are mutually incompatible\n";
    usage(2);
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', count_threshold) || count_threshold < 1) {
      cerr << "The count threshold (-c) option must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will drop bits present in " << count_threshold << " or fewer molecuoes\n";
    }
  }

  if (cl.option_present('p')) {
    if (!cl.value('p', fractional_threshold) || fractional_threshold <= 0.0f ||
        fractional_threshold > 1.0f) {
      cerr << "The fractional threshold specification (-p) must be a valid fraction\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Fractional threshold " << fractional_threshold << endl;
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (need_to_call_initialise_fingerprints(cl)) {
    if (!initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise general fingerprint options\n";
      usage(17);
    }
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 1) {
      cerr << "The pool size value (-s) must be a whole +ve number\n";
      usage(4);
    }

    if (!allocate_pool(pool_size, pool)) {
      cerr << "Cannot allocate pool " << pool_size << '\n';
      return 3;
    }
  }

  if (cl.number_elements() > 1) {
    cerr << "Sorry, cannot handle multiple input files\n";
    exit(2);
  }

  if (!read_pool(cl[0])) {
    cerr << "Cannot read fingerprints from '" << cl[0] << "'\n";
    return 3;
  }

  if (verbose) {
    cerr << "Read " << pool_size << " fingerprints from '" << cl[0] << "'\n";
  }

  if (1 == pool_size) {
    cerr << "Need multiple fingerprints in order to do MCS\n";
    return 2;
  }

  if (fractional_threshold > 0.0f) {
    count_threshold =
        static_cast<int>(static_cast<float>(pool_size) * fractional_threshold + 0.4999f);

    if (verbose) {
      cerr << "Fractional threshold " << fractional_threshold << " corresponds to "
           << count_threshold << " items\n";
    }
  } else if (count_threshold > pool_size) {
    cerr << "The count threshold (-c) option must not be greater than the number of "
            "fingerprints\n";
    cerr << count_threshold << " too large for " << pool_size << " fingerprints\n";
    return 3;
  }

  int nfp = number_fingerprints() + number_sparse_fingerprints();

  if (0 == nfp) {
    cerr << "No fingerprints to process, cannot continue\n";
    return 2;
  }

  fraction_as_string.set_include_leading_space(1);
  fraction_as_string.initialise(0.0, 1.0, 3);

  if (cl.option_present('w')) {
    if (1 == nfp) {  // good, no confusion
      ;
    } else if (cl.option_present('S')) {  // each goes to separate file
      ;
    } else {
      cerr << "More than one fingerprint present, cannot suppress fingerprint headers\n";
      usage(1);
    }

    write_fingerprint_header = 0;

    if (verbose) {
      cerr << "Fingerprint headers suppressed\n";
    }
  }

  if (cl.option_present('S')) {
    IWString output_stem = cl.string_value('S');

    for (int i = 0; i < number_fingerprints(); i++) {
      IWString fname(output_stem);
      fname << ".FP" << i;
      IWString_and_File_Descriptor output;
      if (!output.open(fname.null_terminated_chars())) {
        cerr << "Cannot open output file '" << fname << "'\n";
        return 3;
      }
      gfp_mcs_fp(i, output);
    }

    for (int i = 0; i < number_sparse_fingerprints(); i++) {
      IWString fname(output_stem);
      fname << ".NC" << i;
      IWString_and_File_Descriptor output;
      if (!output.open(fname.null_terminated_chars())) {
        cerr << "Cannot open output file '" << fname << "'\n";
        return 3;
      }
      gfp_mcs_nc(i, output);
    }
  } else  // output to stdout
  {
    IWString_and_File_Descriptor output(1);

    for (int i = 0; i < number_fingerprints(); i++) {
      gfp_mcs_fp(i, output);
    }

    for (int i = 0; i < number_sparse_fingerprints(); i++) {
      gfp_mcs_nc(i, output);
    }

    output.flush();
  }

  if (verbose) {
  }

  return 0;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = gfp_mcs(argc, argv);

  return rc;
}
