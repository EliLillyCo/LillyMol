#include <stdlib.h>
#include <string>
using namespace std;

#include "cmdline.h"

static void
usage(int rc)
{
  cerr << "Tester for std::string funcion with Command_Line object\n";
  cerr << " -x <...>     will be echo'd\n";

  exit(rc);
}

static int
test_std_string_cmdline(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vx:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "unrecognised_options_encountered\n";
    usage(4);
  }

  if (! cl.option_present('x'))
  {
    cerr << "Tester uses the -x option\n";
    usage(4);
  }

  std::string xvalue;
  if (! cl.value('x', xvalue))
  {
    cerr << "Cannot extract -x value\n";
    return 4;
  }

  cout << "-x option value '" << xvalue << "'\n";

  if (1 == cl.option_count('x'))
  {
    cerr << "Only one -x option, cannot test std_string_value\n";
    return 0;
  }

  xvalue = cl.std_string_value('x', 1);

  cout << "Second -x option '" << xvalue << "'\n";

  return 0;

  return 0;
}

int
main(int argc, char ** argv)
{
  int rc = test_std_string_cmdline(argc, argv);

  return rc;
}
