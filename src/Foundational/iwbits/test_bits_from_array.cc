/*
  Scans a descriptor file for similarity to a given vector
*/

#include <random>
#include <memory>

#include "cmdline.h"
#include "iwbits.h"

const char * prog_name = NULL;

static int verbose = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Tester for construct_from_array_of_ints\n";
  cerr << " -n <niter>     number of iterations to perform\n";
  cerr << " -b <nbits>     smallest size bitvector to test\n";
  cerr << " -B <nbits>     largest  size bitvector to test\n";
  cerr << " -d <dens>      smal;est bit density\n";
  cerr << " -D <dens>      largest  bit density\n";
  cerr << " -r <n>         for each iteration (choice of value above) number of replicates\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
fill_array (int * tmp,
            const int nbits,
            const double density,
            std::mt19937_64 & rng,
            std::uniform_real_distribution<double> & d01)
{
  int rc = 0;

  for (int i = 0; i < nbits; ++i)
  {
    const double r = d01(rng);

    if (r < density)
    {
      tmp[i] = 1;
      rc++;
    }
    else
      tmp[i] = 0;
  }

  return rc;
}

static int
test_bits_from_array(const int replicates,
                     std::uniform_int_distribution<int> & bit_length,
                     std::uniform_real_distribution<double> & bit_density,
                     std::uniform_real_distribution<double> & d01,
                     std::mt19937_64 & rng,
                     int * tmp)
{
  const int nbits = bit_length(rng);
  const double density = static_cast<double>(bit_density(rng));

  for (int i = 0; i < replicates; ++i)
  {
    const int nset = fill_array(tmp, nbits, density, rng, d01);
    IW_Bits_Base b;

    b.construct_from_array_of_ints(tmp, nbits);

    if (b.nset() != nset)
    {
      cerr << "Nset mismatch, array " << nset << " bitvector " << b.nset() << endl;
      return 0;
    }
  }

  return 1;
}


static int
test_bits_from_array (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vn:b:B:d:D:r:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  int niter = 1000;

  if (cl.option_present('n'))
  {
    if (! cl.value('n', niter) || niter < 1)
    {
      cerr << "The number of iterations (-n) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will perform " << niter << " iterations\n";
  }

  int replicates = 1;
  if (cl.option_present('r'))
  {
    if (! cl.value('r', replicates) || replicates < 1)
    {
      cerr << "The number of replicates (-r) option must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "For each choice of parameters, will perform " << replicates << " replicates\n";
  }

  int min_bit_length = 32;
  int max_bit_length = 4096;

  if (cl.option_present('b'))
  {
    if (! cl.value('b', min_bit_length) || min_bit_length < 8)
    {
      cerr << "The min bit length (-b) option must be a valid bit length\n";
      usage(1);
    }

    if (verbose)
      cerr << "min_bit_length " << min_bit_length << endl;
  }

  if (cl.option_present('B'))
  {
    if (! cl.value('B', max_bit_length) || max_bit_length < min_bit_length)
    {
      cerr << "The max bit length (-B) option must be a valid bit length greater than min\n";
      usage(1);
    }

    if (verbose)
      cerr << "max_bit_length " << max_bit_length << endl;
  }

  double min_density = 0.1;
  double max_density = 0.9;

  if (cl.option_present('d'))
  {
    if (! cl.value('d', min_density) || min_density < 0.0001f)
    {
      cerr << "The min bit density (-b) option must be a valid bit density\n";
      usage(1);
    }

    if (verbose)
      cerr << "min_density " << min_density << endl;
  }

  if (cl.option_present('D'))
  {
    if (! cl.value('D', max_density) || max_density < min_density)
    {
      cerr << "The max bit density (-D) option must be a valid bit density greater than min\n";
      usage(1);
    }

    if (verbose)
      cerr << "max_density " << max_density << endl;
  }

  std::random_device rd;
  std::mt19937_64 rng(rd());

  std::uniform_int_distribution<int> bit_length(min_bit_length, max_bit_length);
  std::uniform_real_distribution<double> bit_density(min_density, max_density);

  int * tmp = new int[max_bit_length]; std::unique_ptr<int[]> free_tmp(tmp);

  std::uniform_real_distribution<double> d01;

  int rc = 0;

  for (int i = 0; i < niter; ++i)
  {
    if (! test_bits_from_array(replicates, bit_length, bit_density, d01, rng, tmp))
    {
      cerr << "Failure at iteration " << i << endl;
      rc = 1;
      break;
    }
  }

  if (verbose)
  {
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_bits_from_array(argc, argv);

  return rc;
}
