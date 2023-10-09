/*
  Tester for the command line object
*/

#include <stdlib.h>

#include "cmdline_v2.h"
#include "iwstring_data_source.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

const char * prog_name = nullptr;

static int tests_performed = 0;

static int failures = 0;

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
test_iwcmdline (int argc,
                char ** argv,
                IWString constructor_string,
                int expected_good,
                int expected_non_option_arguments,
                ostream & output)
{
  cerr << "Doing test on '" << constructor_string << "', argc = " << argc << endl;

  tests_performed++;

  Command_Line_v2 cl (argc, argv, constructor_string.null_terminated_chars ());

  if (cl.good () && expected_good)
    ;
  else if (! cl.good () && ! expected_good)
    ;
  else
  {
    cerr << "Good state mismatch, expected " << expected_good << " got " << cl.good () << endl;
    return 0;
  }

  if (cl.number_elements () != expected_non_option_arguments)
  {
    cerr << "Non option argument count mismatch, expected " << expected_non_option_arguments << " got " << cl.number_elements () << endl;
    return 0;
  }

  return output.good ();
}

static int
test_iwcmdline (const const_IWSubstring & buffer,
                int nw,
                char ** argv,
                ostream & output)
{
  assert (nw > 2);

  int argc = nw - 2;

  IWString constructor_string;
  int i = 0;

  (void) buffer.nextword (constructor_string, i);

  const_IWSubstring e;   // the expected result
  (void) buffer.nextword (e, i);

  int expected_good;
  if ('1' == e)
    expected_good = 1;
  else if ('0' == e)
    expected_good = 0;
  else
  {
    cerr << "Unrecognised expected result string '" << e << "'\n";
    return 0;
  }

  (void) buffer.nextword (e, i);

  int expected_non_option_arguments;
  if (! e.numeric_value (expected_non_option_arguments) || expected_non_option_arguments < 0)
  {
    cerr << "Number of expected non option arguments must be a whole non negative number\n";
    return 0;
  }

  IWString token;
  int ndx = 1;
  while (buffer.nextword (token, i))
  {
    char * t = new char[token.length () + 1];
    strcpy (t, token.null_terminated_chars ());
    argv[ndx] = t;
    ndx++;
  }

  return test_iwcmdline (argc, argv, constructor_string, expected_good, expected_non_option_arguments, output);
}

/*
  A test record consists of the pattern, the good state of the object
  once the pattern is parsed and the number of non-option arguments
  expected.
*/

static int
test_iwcmdline (const const_IWSubstring & buffer,
                ostream & output)
{
  if (buffer.starts_with ('#') || 0 == buffer.length ())
    return 1;

  int nw = buffer.nwords ();

  if (nw < 3)
  {
    cerr << "cmdline test cases must be at least three tokens\n";
    return 0;
  }

  char ** argv = new char *[nw - 2];

  argv[0] = "programme_name";

  int rc = test_iwcmdline (buffer, nw, argv, output);

  for (int i = 1; i < nw - 2; i++)
  {
    delete argv[i];
  }

  delete argv;

  return rc;
}

static int
test_iwcmdline (iwstring_data_source & input,
                ostream & output)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (! test_iwcmdline (buffer, output))
    {
      cerr << "Failed test '" << buffer << "', line " << input.lines_read () << endl;
      failures++;
    }
  }

  return output.good ();
}

static int
test_iwcmdline (const char * fname,
                ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_strip_leading_blanks (1);

  return test_iwcmdline (input, output);
}

static int
test_iwcmdline (int argc, char ** argv)
{
  cerr << "On entry, argc " << argc << endl;

  if (1 == argc)
  {
    cerr << "test_iwcmdline:must specify command line options\n";
    usage (4);
  }

  const char * fname = argv[1];

  (void) test_iwcmdline (fname, cout);

  cerr << "Performed " << tests_performed << " tests, " << failures << " failed\n";

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_iwcmdline (argc, argv);

  return rc;
}
