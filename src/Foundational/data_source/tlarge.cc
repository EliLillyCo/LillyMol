#include <stdlib.h>
#include <unistd.h>

/*
  Tester for large file source
*/

#include "cmdline.h"
#include "Large_File.h"

static int verbose = 0;

static void
usage (int rc)
{
  cerr << "Tester for large file reader\n";
  cerr << " -e              echo input file\n";
  cerr << " -v              verbose output\n";

  exit (rc);
}

static int
do_test_echo (IW_Large_File & input,
              int output)
{
  IWString buffer;

  char newline = '\n';

  while (input.next_record (buffer))
  {
    (void) write (output, buffer.rawchars (), buffer.length ());
    write (output, &newline, sizeof (newline));
  }

  return 1;
}

static int
do_test_echo (const char * fname, int output)
{
  IW_Large_File input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return do_test_echo (input, output);
}

static int
tlarge (IW_Large_File & input)
{
  IWString buffer;

  while (input.next_record (buffer))
  {
  }

  cerr << "Read " << input.records_read () << " records\n";
  return 1;
}

static int
tlarge (const char * fname)
{
  IW_Large_File input;

  if (! input.open (fname))
  {
    cerr << "Cannot open via open () '" << fname << "'\n";
    return 0;
  }

  return tlarge (input);
}

static int
tlarge (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "ve");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised options encountered\n";
    usage (3);
  }

  verbose = cl.option_count ('v');

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (3);
  }

  int test_echo = 0;

  if (cl.option_present ('e'))
  {
    test_echo = 1;
  }

  int rc;

  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (test_echo)
      rc = do_test_echo (cl[i], 1);
    else 
      rc = tlarge (cl[i]);

    if (0 == rc)
    {
      rc = i + 1;
      break;
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = tlarge (argc, argv);

  return rc;
}
