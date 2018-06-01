/*
  Scans a descriptor file for similarity to a given vector
*/

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "cmdline.h"
#include "iwstring_data_source.h"

#include "misc.h"

using std::cout;
using std::ostream;

const char * prog_name = NULL;

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

#define BUFFER_SIZE 6

static unsigned char buffer[BUFFER_SIZE];

static int
test_iwuuencode (const unsigned char * buffer,
                 int bytes_read,
                 ostream & output)
{
  IWString encoded_form;

  if (! IWuuencode_append (buffer, bytes_read, encoded_form))
  {
    cerr << "Yipes, cannot encode " << bytes_read << " bytes\n";
    return 0;
  }

  output << encoded_form;

  return output.good ();
}

static int
test_iwuuencode (int fd,
                 ostream & output)
{
  int bytes_read;

  while ((bytes_read = read (fd, buffer, BUFFER_SIZE)) > 0)
  {
    if (! test_iwuuencode (buffer, bytes_read, output))
      return 0;
  }

  return 1;
}

static int
test_iwuuencode (const char * fname,
                 ostream & output)
{
  int oflag = O_RDONLY;

#ifdef IWSUN
  int fd = open64 (fname, oflag);
#else
  int fd = open   (fname, oflag);
#endif

  if (fd < 0)
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_iwuuencode (fd, output);
}

static int
test_iwuuencode (int argc, char ** argv)
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
    if (! test_iwuuencode (cl[i], cout))
    {
      rc = i + 1;
      break;
    }
  }

  cout << endl;

  if (verbose)
  {
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_iwuuencode (argc, argv);

  return rc;
}
