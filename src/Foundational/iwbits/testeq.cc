/*
  I was interested in whether memcmp was faster than my explicit loops
*/

#include <stdlib.h>

#include "cmdline.h"
#include "iwbits.h"
#include "iwrandom.h"

const char * prog_name = NULL;

static int verbose = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Time tester for bit equality operator\n";
  cerr << " -n <iter>      number of iterations to do\n";
  cerr << " -b <number>    bit size to use\n";
  cerr << " -e             do the equality test\n";
  cerr << " -r <number>    re-initialise the bits every <number> iterations\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static void
initialise_randomly (IW_Bits_Base & b,
                     int * tmp,
                     int bsize)
{
  for (int i = 0; i < bsize; i++)
  {
    tmp[i] = intbtwij (0, 1);
  }

  b.construct_from_array_of_ints (tmp, bsize);

  return;
}


static int
testbeq (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:b:er:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  int iterations;

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', iterations) || iterations < 1)
    {
      cerr << "The number of iterations must be a positive whole number\n";
      usage (1);
    }

    if (verbose)
      cerr << "Will do " << iterations << " iterations\n";
  }
  else
  {
    iterations = 1000;

    if (verbose)
      cerr << "doing " << iterations << " iterations by default\n";
  }

  int bsize;

  if (cl.option_present ('b'))
  {
    if (! cl.value ('b', bsize) || bsize < 1)
    {
      cerr << "The size of the bit vectors (-b option) must be a positive whole number\n";
      usage (1);
    }
  }
  else
  {
    bsize = 1024;
  }

  if (verbose)
    cerr << "Will test bits vectors of length " << bsize << " bits\n";

  int * tmp = new int[bsize];

  int do_equality_check = 0;

  if (cl.option_present ('e'))
  {
    do_equality_check = 1;

    if (verbose)
      cerr << "Will do the equality check\n";
  }

  int nchange = 0;

  if (cl.option_present ('r'))
  {
    if (! cl.value ('r', nchange) || nchange < 1)
    {
      cerr << "The change bit option must be a valid whole number\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will re-initialise the bits every " << nchange << " iterations\n";
  }

  IW_Bits_Base b1, b2;

  initialise_randomly (b1, tmp, bsize);
  initialise_randomly (b2, tmp, bsize);

  for (int i = 0; i < iterations; i++)
  {
    if (nchange && 0 == i % nchange)
    {
      initialise_randomly (b1, tmp, bsize);
      initialise_randomly (b2, tmp, bsize);
    }

    (void) (b1 == b2);

    if (do_equality_check)
    {
      (void) (b1 == b2);
      (void) (b2 == b1);
    }
  }

  if (verbose)
  {
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = testbeq (argc, argv);

  return rc;
}
