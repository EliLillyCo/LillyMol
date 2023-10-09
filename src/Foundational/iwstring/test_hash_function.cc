/*
  I have some strings that seem to be doing poorly with the hash
  function.
*/

#include <stdlib.h>

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "iwhash.h"

const char * prog_name = nullptr;

static int verbose = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id$\n";
  cerr << "What does this programme do?\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
test_hash_function (iwstring_data_source & input,
                    IWString_and_File_Descriptor & output)
{
  IWString buffer;

  IWStringHash hasher;

  while (input.next_record(buffer))
  {
    buffer.truncate_at_first(' ');

    size_t h = hasher(buffer);
    output << buffer << ' ' << h << '\n';

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
test_hash_function (const char * fname,
     IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_hash_function(input, output);
}


static int
test_hash_function (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "v");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! test_hash_function(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_hash_function(argc, argv);

  return rc;
}
