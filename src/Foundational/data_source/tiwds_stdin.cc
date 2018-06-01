/*
  Tests the iwstring_data_source object reading from stdin
*/

#include <stdlib.h>

#include "cmdline.h"
#include "iwstring_data_source.h"

using std::cerr;
using std::cout;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id$\n";
  cerr << "Echo stdin input to stdout\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
tiwds_stdin (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "v");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.number_elements ())
    cerr << "Command line arguments ignored\n";

  iwstring_data_source input ("-");

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    cout << buffer << endl;
  }

  if (verbose)
    cerr << "Read " << input.lines_read () << " lines from stdin\n";

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tiwds_stdin (argc, argv);

  return rc;
}
