#include <stdlib.h>

#include "cmdline.h"
#include "iwstring_data_source.h"

#include "running_average.h"

using std::cout;
using std::ostream;

static int verbose = 0;

static int nsample = 10;

typedef double mytype;
typedef Running_Average<mytype> my_running_average;

static void
usage (int rc)
{
  cerr << "Tester for the running average object\n";
  cerr << " -n <number>     how many samples to retain\n";
  cerr << " -v              verbose output\n";

  exit (rc);
}

static int
tra (iwstring_data_source & input, ostream & output)
{
  my_running_average mra (nsample);

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    const_IWSubstring token;
    buffer.word (0, token);

    mytype t;
    if (! token.numeric_value (t))
    {
      cerr << "Yipes, non numeric value '" << buffer << "', line " << input.lines_read () << endl;
      return 0;
    }

    mra.extra (t);

    if (mra.nsamples () > 1)
      output << "After " << mra.nsamples () << " samples, average " << mra.average () << endl;
  }

  return output.good ();
}

static int
tra (const char * fname, ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return tra (input, output);
}

static int
tra (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised_options_encountered\n";
    usage (3);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', nsample) || nsample < 2)
    {
      cerr << "The -n option must be followed by a valid number of samples\n";
      usage (7);
    }

    if (verbose)
      cerr << "Will compute a " << nsample << " point running average\n";
  }
  else
  {
    cerr << "Default sample size " << nsample << endl;
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (3);
  }

  int rc = tra (cl[0], cout);

  return ! rc;
}

int
main (int argc, char ** argv)
{
  int rc = tra (argc, argv);

  return rc;
}
