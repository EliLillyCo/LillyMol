#include <stdlib.h>

/*
  Tester for random number stuff
*/

#include <iostream>

#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

const char * prog_name = NULL;

static int verbose = 0;

#include "iwrandom.h"
#include "iwaray.h"

static int
test_one_zero (int n,
               random_number_seed_t seed)
{
  Random_Number_Working_Storage rnws0;
  Random_Number_Working_Storage rnws1;

  rnws0.set_seed (seed);
  rnws1.set_seed (seed);

  int nzero = 0;
  int none = 0;

  int rc = 0;
  for (int i = 0; i < n; i++)
  {
    int i0 = rnws0.random_one_or_zero ();
    int i1 = rnws1.random_one_or_zero ();
    if (i0 != i1)
    {
      std::cerr << "Error, iteration " << i << " i0 = " << i0 << " i1 = " << i1 << std::endl;
//    rnws0.debug_print (std::cerr);
//    rnws1.debug_print (std::cerr);
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
      seed = iw_random_seed ();
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
