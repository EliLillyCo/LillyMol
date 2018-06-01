#include <stdlib.h>
#include <iostream>
using std::cerr;
using std::endl;

/*
  Tester for random number stuff
*/

#include <iostream>

const char * prog_name = NULL;

static int verbose = 0;

#include "iwrandom.h"
#include "iwaray.h"

static int
test_one_zero (int n,
               random_number_seed_t seed)
{
  MTRand_int32 r1(seed);
  MTRand_int32 r2(seed);

  int nzero = 0;
  int none = 0;

  int rc = 0;
  for (int i = 0; i < n; i++)
  {
    int i1 = (r1.closed_closed() > 0.5);
    int i2 = (r2.closed_closed() > 0.5);
    if (i1 != i2)
    {
      std::cerr << "Error, iteration " << i << " i1 = " << i1 << " i2 = " << i2 << std::endl;
      rc++;
    }
    else if (0 == i1)
      nzero++;
    else if (1 == i1)
      none++;
    else
    {
      std::cerr << "Error, iteration " << i << " value not 0 or 1 " << i1 << std::endl;
//    rnws0.debug_print (std::cerr);
//    rnws1.debug_print (std::cerr);
      rc++;
    }
  }

  if (verbose)
    std::cerr << "After " << n << " iterations, " << nzero << " 0 and " << none << " 1 values\n";

  if (0 == nzero || 0 == none || nzero / none > 1 || none / nzero > 1)
  {
    std::cerr << "There seems to be an imbalance of 0 (" << nzero<< ") and 1 (" << none << ")\n";
    rc++;
  }

  return rc;
}

static int
test_uniform (int n,
              random_number_seed_t seed)
{
  int count[1000];
  for (int i = 0; i < 1000; ++i)
  {
    count[i] = 0;
  }

  MTRand_int32 rng(seed);

  for (int i = 0; i < n; ++i)
  {
    double r = rng.open_open();
    const int j = static_cast<int>(r * 1000.0);
    count[j]++;
  }

  int min_count = count[0];
  int max_count = count[0];
  for (int i = 1; i < 1000; ++i)
  {
    if (count[i] < min_count)
      min_count = count[i];
    else if (count[i] > max_count)
      max_count = count[i];
  }

  cerr << "Counts btw " << min_count << " and " << max_count << endl;

  for (int i = 0; i < 1000; ++i)
  {
    std::cout << i << ' ' << count[i] << endl;
  }

  return 1;
}

static void
usage (int rc)
{
  std::cerr << "Usage: " << prog_name << " <options>\n";
  std::cerr << "  -n <number>   the number of iterations to perform\n";
  std::cerr << "  -s RAND       use a random randon number seed\n";
  std::cerr << "  -s <number>   use <number> as random number seed\n";
  std::cerr << "  -v            verbose output\n";
  exit (rc);
}

#include "cmdline.h"

static int
tiwrandom (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:s:");

  int rc = 0;
  verbose = cl.option_count ('v');

  int n = 300;

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', n) || n <= 0)
    {
      std::cerr << "The -n option requires a positive int\n";
      usage (2);
    }

    if (verbose)
      std::cerr << "Will perform " << n << " iterations\n";
  }

  random_number_seed_t seed;

  if (cl.option_present ('s'))
  {
    IWString tmp;
    if (cl.value ('s', tmp) && tmp.starts_with ("RAND"))
    {
      seed = iw_random_seed();
    }
    else if (! cl.value ('s', seed) || 0 == seed)
    {
      std::cerr << "The -s option must be followed by 'RAND' or a non zero number\n";
      usage (3);
    }
  }
  else
    seed = iw_random_seed();

  if (verbose)
    std::cerr << "Using seed " << seed << std::endl;

  rc += test_one_zero (n, seed);

  rc += test_uniform(n, seed);

  if (0 == verbose)
    return rc;

  if (0 == rc)
    std::cerr << "All tests successful\n";
  else
    std::cerr << rc << " failures\n";
    
  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tiwrandom (argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status (stderr);
#endif

  return rc;
}
