#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <cmath>
#include <iomanip>
using std::setw;

/*
  Reads a file of tdt's and generates statistics on the fingerprints set.
*/

#include "iwbits.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwminmax.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/smiles.h"

char * prog_name = nullptr;
int verbose = 0;

static int fingerprints_read = 0;

/*
  The counter for bits.
*/

#define DEFAULT_FINGERPRINT_LENGTH 2048
static int max_fingerprint_length = DEFAULT_FINGERPRINT_LENGTH;

static int * population = nullptr;

/*
  We also keep track of the lengths of the fingerprints.
*/

static int * fingerprint_length_count = nullptr;
static int longest_fingerprint_read = 0;
static unsigned int total_bits_read = 0;

static int
iw_log2(int x)
{
  return int(log(x) / log(2.0) + 0.1);
}

static int
fppop(const IW_Bits_Base & fp)
{
  assert(fp.nbits() <= max_fingerprint_length);

  if (fp.nbits() > longest_fingerprint_read)
    longest_fingerprint_read = fp.nbits();

  total_bits_read += fp.nbits();

  int tmp = iw_log2(fp.nbits());
  fingerprint_length_count[tmp]++;

  return fp.increment_vector(population);
}

static int
fppop(iwstring_data_source & input)
{
  fingerprints_read = 0;    // reset on a per file basis.
  IWString buffer;
  while (input.next_record(buffer))
  {
    if (! buffer.starts_with("FP<"))
      continue;

    IW_Bits_Base fp;
    if (! fp.construct_from_tdt_record(buffer))
    {
      cerr << "test_iwfp: cannot construct fingerprint from '" << buffer << "'\n";
      return 0;
    }

    fp.set_id(fingerprints_read++);

    (void) fppop(fp);
  }

  cerr << "Read " << fingerprints_read << " fingerprints, longest = " << longest_fingerprint_read << endl;

  int tmp = iw_log2(max_fingerprint_length) + 1;
  int j = 1;
  int recount_fingerprints_read = 0;
  for (int i = 0; i < tmp; i++)
  {
    cerr << setw(4) << fingerprint_length_count[i] << " fingerprints of length " << setw(4) << j << endl;
    recount_fingerprints_read += fingerprint_length_count[i];
    j *= 2;
  }

  if (recount_fingerprints_read != fingerprints_read)
  {
    cerr << "Major internal error, counted " << fingerprints_read << 
            ", accumulators show " << recount_fingerprints_read << endl;
    cerr << "Use results with caution\n";
  }

  for (int i = 0; i < longest_fingerprint_read; i++)
  {
    float tmp = float(population[i]) / float(fingerprints_read);
    cerr << "Bit " << setw(4) << i << " " << setw(5) << population[i] << " hits, probability " << tmp << endl;
  }

// Report highest and lowest probabilities

  unsigned int bits_with_zero_hits = 0;
  unsigned int total_hits = 0;

  iwminmax<int> hits(fingerprints_read, 0);
  for (int i = 0; i < longest_fingerprint_read; i++)
  {
    if (0 == population[i])
      bits_with_zero_hits++;
    else
      total_hits += population[i];

    hits.try_this(population[i]);
  }

  cerr << "Highest nhits = " << hits.maxval() << ", lowest = " << 
          hits.minval() << endl;
  if (bits_with_zero_hits)
    cerr << bits_with_zero_hits << " bits had zero hits\n";

  return 1;
}

/*
  Even though this is a standard library thing for me, I reproduce it here to
  avoid linkage problems
*/

static int *
new_int(int size, int initial_value = 0)
{
  assert(size > 0);

  int * rc = new int[size];

  for (int i = 0; i < size; i++)
    rc[i] = initial_value;

  return rc;
}

#define DELETE_AND_SET_NULL(p) { assert(nullptr != (p)); delete(p); p = nullptr;}

static int
fppop(const char * fname)
{
  if (nullptr == strstr(fname, ".tdt"))
    cerr << "fppop: warning, name of tdt file does not look like tdt '" << fname << "'\n";

  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 1;
  }

  assert(nullptr == population);

  population = new_int(max_fingerprint_length);
  fingerprint_length_count = new_int(iw_log2(max_fingerprint_length) + 1);

  int rc = fppop(input);

  DELETE_AND_SET_NULL(population);
  DELETE_AND_SET_NULL(fingerprint_length_count);

  return rc;
}

static void
usage(int rc)
{
  cerr << "Usage: " << prog_name << " <options> <pool file> <file1> <file2>...\n";
  cerr << " -c                  check fingerprint bit counts\n";
  cerr << " -v                  verbose output\n";

  exit(rc);
}

static int
fppop(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vftT:n:c");
  verbose = cl.option_count('v');

  if (cl.option_present('c'))
    set_check_fingerprint_bit_counts(1);

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  for (int i = 0; i < cl.number_elements(); i++)
    (void) fppop(cl[i]);

  return 0;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = fppop(argc, argv);

  return rc;
}
