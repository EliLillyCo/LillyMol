/*
  Test the Daylight string representation stuff
*/

#include <stdlib.h>
#include <unistd.h>

#include "iwbits.h"
#include "iwrandom.h"
#include "cmdline.h"

static int verbose = 0;

static void
usage (int rc)
{
  cerr << " -n <ntest>        how many tests to perform\n";
  cerr << " -m <nbits>        max size for fingerprints\n";
  exit (rc);
}

static int
tdaylight (int bsize, const int * v)
{
  IW_Bits_Base fp;

  (void) fp.construct_from_array_of_ints (v, bsize);

  for (int i = 0; i < bsize; i++)
  {
    int iset = fp.is_set (i);
    if (v[i] && iset)
      ;
    else if (0 == v[i] && 0 == iset)
      ;
    else
    {
      cerr << "Bit set mismatch, bit " << i << " input " << v[i] << " set " << iset << endl;
      return 0;
    }
  }

  IWString ascii;
  fp.daylight_ascii_representation (ascii);

//cerr << "Ascii representation '" << ascii << "'\n";

  IW_Bits_Base fp2;

  if (! fp2.construct_from_daylight_ascii_bit_rep (ascii.rawchars (), ascii.length ()))
  {
    cerr << "Yipes, cannot parse ASCII '" << ascii << "'\n";
    return 0;
  }

  if (fp2 == fp)
    return 1;

  cerr << "Yipes, fingerprints differ\n";

  if (fp.nbits() != fp2.nbits())
    cerr << "nb " << fp.nbits() << " vs " << fp2.nbits() << endl;
  if (fp.nset() != fp2.nset())
    cerr << "nset " << fp.nset() << " vs " << fp2.nset() << endl;

  fp.printon(cerr);
  cerr << endl;
  fp2.printon(cerr);
  cerr << endl;
  cerr << ascii << endl;
  fp2.daylight_ascii_representation(ascii);
  cerr << ascii << endl;

  return 0;
}

static int
tdaylight (int max_size)
{
  int bsize = intbtwij (1, max_size);

  if (0 != bsize % 8)
    bsize = bsize + 8 - bsize % 8;

  int * tmp = new int[bsize];

  for (int i = 0; i < bsize; i++)
  {
    tmp[i] = intbtwij (0, 1);
  }

#define CHECK_DISTRIBUTION
#ifdef CHECK_DISTRIBUTION
  int n1 = 0;
  for (int i = 0; i < bsize; i++)
  {
    if (tmp[i])
      n1++;
  }

  cerr << n1 << " of " << bsize << " 1 bits (" << int (float (100 * n1) / float (bsize)) << "%)\n";
#endif

  int rc = tdaylight (bsize, tmp);

  delete tmp;

  return rc;
}

static int
tdaylight (int ntest, int max_size)
{
  for (int i = 0; i < ntest; i++)
  {
    if (! tdaylight (max_size))
    {
      cerr << "Test " << i << " failed\n";
      return 0;
    }
  }

  return ntest;
}

static int
tdaylight (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:m:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  int ntest = 1;

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', ntest) || ntest < 1)
    {
      cerr << "The -n option must be followed by a whole positive number\n";
      usage (3);
    }

    if (verbose)
      cerr << "Will perform " << ntest << " tests\n";
  }

  int max_size = 2048;

  if (cl.option_present ('m'))
  {
    if (! cl.value ('m', max_size) || max_size < 1)
    {
      cerr << "The -m option must be followed by a whole positive number\n";
      usage (4);
    }

    if (verbose)
      cerr << "Fingerprints will be a max length of " << max_size << " bits\n";
  }

//iw_random_seed ();     // initialise random number generator with random number

  if (! tdaylight (ntest, max_size))
  {
    cerr << "Tests failed\n";
    return 1;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = tdaylight (argc, argv);

  return rc;
}
