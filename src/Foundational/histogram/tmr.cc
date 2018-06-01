/*
  Tester for the most recent object
*/

#include <stdlib.h>

#include "cmdline.h"
#include "most_recent.h"
using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

static int verbose = 0;

static int errors = 0;

static void
usage (int rc)
{
  cerr << "Tester for the most recent object\n";
  cerr << " -n <number>     the number of samples to test\n";
  cerr << " -s <number>     the number of samples to save\n";
  cerr << " -t <number>     how often to test\n";
  cerr << " -v              verbose output\n";

  exit (rc);
}

static int 
test_most_recent (int nsamples,
                  int ntest,
                  int mrsize,
                  ostream & os)
{
  IWMost_Recent<int> mr (mrsize);

  for (int i = 1; i < nsamples; i++)
  {
    mr.extra (i);

    if (i > 0 && (0 == i % ntest || i == nsamples - 1))
    {
      int nstored = mr.items_stored ();

      if (verbose)
        os << "Checking after " << mr.nsamples () << " samples, " << nstored << " items stored\n";

      int should_be = i - nstored + 1;

      for (int i = 0; i < nstored; i++)
      {
        IWMost_Recent<int>::const_reference j = mr[i];
        if (verbose)
          os << j;

        if (j != should_be)
        {
          if (0 == verbose)
            os << j;
          os << "   yipes, should be " << should_be;
          errors++;
        }

        if (verbose)
          os << endl;

        should_be++;
      }
    }
  }

// Test the copy operator

  IWMost_Recent<int> mr2 (mr);

  int nstored = mr2.items_stored ();

  if (nstored != mr.items_stored ())
  {
    cerr << "Copy constructor did not preserve items stored " << mr2.items_stored () << " vs " << mr.items_stored () << endl;
    errors++;

    return 0;
  }

  for (int i = 0; i < nstored; i++)
  {
    IWMost_Recent<int>::const_reference j1 = mr[i];
    IWMost_Recent<int>::const_reference j2 = mr2[i];
    if (j1 != j2)
    {
      cerr << "Value mis-match in copy constructor, i = " << i << " values " << j1 << " and " << j2 << endl;
      errors++;
    }
  }

  if (verbose)
    cerr << "Copy constructor OK\n";

  return os.good ();
}

int
test_most_recent (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:s:t:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised_options_encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  int nsamples;
  if (! cl.option_present ('n'))
  {
    nsamples = 1000;
    cerr << "The -n option (number of samples) absent, default " << nsamples << endl;
  }
  else if (! cl.value ('n', nsamples) || nsamples <= 0)
  {
    cerr << "The -n option must be followed by a positive whole number\n";
    usage (5);
  }
  else if (verbose)
    cerr << "Will test " << nsamples << " samples\n";

  int mrsize;
  if (! cl.option_present ('s'))
  {
    mrsize = 11;
    cerr << "The -s option (number of samples to save) absent, default " << mrsize << endl;
  }
  else if (! cl.value ('s', mrsize) || mrsize <= 0)
  {
    cerr << "The -s option must be followed by a positive whole number\n";
    usage (5);
  }
  else if (verbose)
    cerr << "Will save the " << mrsize << " most recent samples\n";

  int ntest;
  if (! cl.option_present ('t'))
  {
    ntest = 13;
    cerr << "The -t option (frequency of testing) absent, default " << ntest << endl;
  }
  else if (! cl.value ('t', ntest) || ntest <= 0)
  {
    cerr << "The -t option must be followed by a positive whole number\n";
    usage (5);
  }
  else if (verbose)
    cerr << "Will test results every " << ntest << " steps\n";

  if (cl.number_elements ())
    cerr << "Arguments ignored\n";

  test_most_recent (nsamples, ntest, mrsize, cout);

  cerr << "Encountered " << errors << " errors\n";

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = test_most_recent (argc, argv);

  return rc;
}
