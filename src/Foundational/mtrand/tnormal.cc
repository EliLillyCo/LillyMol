/*
  Tester for normal distribution generators
*/

#include <stdlib.h>
#include <iostream>
using std::cout;

#include "cmdline.h"
#include "iwrandom.h"
#include "iwhistogram.h"
#include "accumulator.h"

const char * prog_name = NULL;

static int verbose = 0;

static int random_numbers_to_generate = 100;

static unsigned int seed = 0;

static Accumulator<double> acc;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Tests random number generators\n";
  cerr << " -H min,max,dx  histogram specifications\n";
  cerr << " -n <n>         number of random numbers to generate\n";
  cerr << " -s <int>       random number seed\n";
  cerr << " -R uniform     test uniform random number generator\n";
  cerr << " -R BM          test Box Muller polar random number generator\n";
  cerr << " -R inverse     test inverse normal random number generator\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
write_histogram (const IWHistogram & histogram,
                 std::ostream & os)
{
  return histogram.write_terse(os);
}

static int
test_uniform_random_distribution (IWHistogram & histogram,
                                  std::ostream & os)
{
  Random_Number_Working_Storage rnws;

  if (0 != seed)
    rnws.seed(seed);

  for (int i = 0; i < random_numbers_to_generate; i++)
  {
    double r = rnws.random_number();

    acc.extra(r);

    histogram.extra(r);
  }

  return write_histogram(histogram, os);
}

static int
test_box_muller_distribution (IWHistogram & histogram,
                              std::ostream & os)
{
  IW_Box_Muller_Normally_Distributed bmnd;

  if (0 != seed)
    bmnd.seed(seed);

  for (int i = 0; i < random_numbers_to_generate; i++)
  {
    double r = bmnd();

    acc.extra(r);

    histogram.extra(r);
  }

  return write_histogram(histogram, os);
}

static int
test_inverse_normal_distribution (IWHistogram & histogram,
                                  std::ostream & os)
{
  IW_Inverse_Normal_Normally_Distributed innd;

  if (0 != seed)
    innd.seed(seed);

  for (int i = 0; i < random_numbers_to_generate; i++)
  {
    double r = innd();

    acc.extra(r);

    histogram.extra(r);
  }

  return write_histogram(histogram, os);
}

static int
tnormal (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vH:n:s:R:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('n'))
  {
    if (! cl.value('n', random_numbers_to_generate) || random_numbers_to_generate < 2)
    {
      cerr << "Must specify at least two random numbers to generate via the -n option\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will generate " << random_numbers_to_generate << " random numbers\n";
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', seed))
    {
      cerr << "Invalid seed (-s)\n";
      usage(2);
    }

    if (verbose)
      cerr << "Seed set to " << seed << endl;
  }
  else
  {
    seed = iw_random_seed();

    if (verbose)
      cerr << "Default seed " << seed << endl;
  }

  IWHistogram histogram;

  if (cl.option_present('H'))
  {
    const_IWSubstring h = cl.string_value('H');

    if (! histogram.initialise(h))
    {
      cerr << "Invalid histogram initialisation conditions '" << h << "'\n";
      usage(3);
    }
  }
  else
    histogram.initialise(0.0, 1.0, 0.01);

  if (cl.number_elements())
    cerr << "Command line arguments ignored\n";

  if (! cl.option_present('R'))
  {
    cerr << "Must specify random number type via the -R option\n";
    usage(2);
  }

  if (cl.option_present('R'))
  {
    const_IWSubstring r = cl.string_value('R');

    if (r.starts_with("unif"))
      test_uniform_random_distribution(histogram, cout);
    else if ("BM" == r)
      test_box_muller_distribution(histogram, cout);
    else if (r.starts_with("inv"))
      test_inverse_normal_distribution(histogram, cout);
    else
    {
      cerr << "Unrecognised distribution type '" << r << "'\n";
      usage(2);
    }
  }

  if (verbose && acc.n() > 1)
  {
    cerr << "Recorded " << acc.n() << " values btw " << static_cast<float>(acc.minval()) << " and " << static_cast<float>(acc.maxval()) << " ave " << static_cast<float>(acc.average()) << " var " << static_cast<float>(acc.variance()) << endl;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tnormal(argc, argv);

  return rc;
}
