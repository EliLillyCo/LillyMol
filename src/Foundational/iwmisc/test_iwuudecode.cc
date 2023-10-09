/*
*/

#include <stdlib.h>
#include <unistd.h>
#include <memory>

#include "cmdline.h"
#include "iwstring_data_source.h"

#include "misc.h"

const char * prog_name = nullptr;

static int verbose = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id$\n";
  cerr << "What does this programme do?\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
test_iwuudecode (const const_IWSubstring & buffer,
                 ostream & output)
{
  int nbytes = IWuudecode_bytes_needed (buffer.length ());

//cerr << "String '" << buffer << "' of length " << buffer.length () << " needs " << nbytes << " bytes\n";

  unsigned char * b = new unsigned char[nbytes]; std::unique_ptr<unsigned char> free_b (b);

  IWuudecode (buffer, b);

  write (1, b, nbytes);
//output.write (b, nbytes);

  return output.good ();
}

static int
test_iwuudecode (iwstring_data_source & input,
                 ostream & output)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (! test_iwuudecode (buffer, output))
      return 0;
  }

  return output.good ();
}

static int
test_iwuudecode (const char * fname,
                 ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_iwuudecode (input, output);
}

static int
test_iwuudecode (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "v");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! test_iwuudecode (cl[i], cout))
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

  int rc = test_iwuudecode (argc, argv);

  return rc;
}
