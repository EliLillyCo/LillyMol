/*
  Tester for iwzlib
*/

#include <stdlib.h>

#include "cmdline.h"
#include "iwzlib.h"

using std::cerr;
using std::cout;
using std::ostream;

static void
usage (int rc)
{
  exit (rc);
}

static int
test_iwzlib (IW_ZLib_Wrapper & iwzlwrap,
             ostream & output)
{
  IWString buffer;

  while (iwzlwrap.next_record (buffer))
  {
    output << buffer << endl;
  }

  return output.good ();
}

static int
test_iwzlib (const char * fname,
             ostream & output)
{
  IW_ZLib_Wrapper iwzlwrap;

  if (! iwzlwrap.open_file (fname))
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_iwzlib (iwzlwrap, output);
}

static int
test_iwzlib (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "v");

  if (0 == cl.number_elements ())
  {
    cerr << "Must specify file(s) to check on the command line\n";
    usage (4);
  }

  int rc = 0;

  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! test_iwzlib (cl[i], cout))
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
  int rc = test_iwzlib (argc, argv);

  return rc;
}
