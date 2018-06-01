/*
  Test various scenarios for file writing
*/

#include <stdlib.h>
#include <unistd.h>

#include "cmdline.h"
#include "iwstring_data_source.h"

const char * prog_name = NULL;

static int verbose = 0;

static int just_write_block = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id$\n";
  cerr << "Tests writing speed from an IWString\n";
  cerr << " -b             just write a block of data, rather than whole string\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

#define IW_BLK_SIZE 4096

static int
test_write_speed (iwstring_data_source & input,
                  int output)
{
  IWString output_buffer;

  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    output_buffer << buffer << '\n';

    if (output_buffer.length () < IW_BLK_SIZE)
      continue;

    if (just_write_block)
    {
      write (output, output_buffer.rawchars (), IW_BLK_SIZE);

      output_buffer.remove_leading_chars (IW_BLK_SIZE);
    }
    else
    {
      output_buffer.write (output);
      output_buffer.resize_keep_storage (0);
    }
  }

  if (output_buffer.length ())
    output_buffer.write (output);

  return 1;
}

static int
test_write_speed (const char * fname,
                  int output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_write_speed (input, output);
}


static int
test_write_speed (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vb");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('b'))
  {
    just_write_block = 1;

    if (verbose)
      cerr << "Will just write a block\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! test_write_speed (cl[i], 1))
    {
      rc = i + 1;
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

  int rc = test_write_speed (argc, argv);

  return rc;
}
