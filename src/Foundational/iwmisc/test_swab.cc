#include <stdlib.h>
#include <limits>

#ifdef __WIN32__
#include <winsock2.h>
#else
#include <netinet/in.h>
#endif

#include "cmdline.h"
#include "iwrandom.h"
#include "misc.h"

using std::numeric_limits;

static int verbose = 0;

static int test_htonl = 0;

static void
usage (int rc)
{
  cerr << "Speed benchmarking for byte swapping\n";
  cerr << " -h          use htonl\n";
  cerr << " -c          check for correctness - does inversion\n";
  cerr << " -l <len>    length of vector to test\n";
  cerr << " -n <num>    number of tests to do\n";
  cerr << " -v          verbose output\n";

  exit (rc);
}

static void
do_test_htonl (int nw, unsigned int * b)
{
  for (int i = 0; i < nw; i++)
  {
    b[i] = htonl (b[i]);
  }

  return;
}

static int
test_swab (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vl:n:ch");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised_options_encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  int vector_length;

  if (cl.option_present ('l'))
  {
    if (! cl.value ('l', vector_length) || vector_length < 1)
    {
      cerr << "Invalid vector length (-l option)\n";
      usage (2);
    }
  }
  else
  {
    vector_length = 100;
  }

  if (verbose)
    cerr << "Will test vector length " << vector_length << endl;

  int number_tests;

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', number_tests) || number_tests < 1)
    {
      cerr << "Invalid number_tests (-n option)\n";
      usage (2);
    }
  }
  else
  {
    number_tests = 100;
  }

  if (verbose)
    cerr << "Will test " << number_tests << " times\n";

  if (cl.option_present ('h'))
  {
    test_htonl = 1;

    if (verbose)
      cerr << "Will test htonl\n";
  }

  int check_correctness;

  if (cl.option_present ('c'))
  {
    check_correctness = 1;
  }
  else
  {
    check_correctness = 0;
  }

  unsigned int * b = new unsigned int[vector_length];

  for (int i = 0; i < vector_length; i++)
  {
    b[i] = intbtwij (0, numeric_limits<unsigned int>::max());
  }

  if (check_correctness)
  {
    unsigned int * bsave = new unsigned int[vector_length];

    for (int i = 0; i < vector_length; i++)
    {
      bsave[i] = b[i];
    }

    if (test_htonl)
    {
      do_test_htonl (vector_length, b);
      do_test_htonl (vector_length, b);
    }
    else
    {
      rick_higgs_byte_swap (vector_length, b);
      rick_higgs_byte_swap (vector_length, b);
    }

    for (int i = 0; i < vector_length; i++)
    {
      if (bsave[i] != b[i])
      {
        cerr << "Yipes, mismatch i = " << i << " b[i] " << b[i] << " saved " << bsave[i] << endl;
      }
    }

    delete bsave;

    if (verbose)
      cerr << "Multiple swap OK\n";
  }

  if (test_htonl)
  {
    for (int i = 0; i < number_tests; i++)
    {
      do_test_htonl (vector_length, b);
    }
  }
  else
  {
    for (int i = 0; i < number_tests; i++)
    {
      rick_higgs_byte_swap (vector_length, b);
    }
  }

  delete b;

  return 0;
}

int
main (int argc, char ** argv)
{

  int rc = test_swab (argc, argv);

  return rc;
}
