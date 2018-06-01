#include <stdlib.h>

/*
  Tester for compiled regular expressions
*/

#include <iostream>
using namespace std;

// #define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "iwcrex.h"
#include "iwlrx.h"

const char * prog_name = NULL;

static int verbose = 0;

static int failures = 0;

typedef IW_Regular_Expression_Template<IW_lrx> IW_Regular_Expression_Test;

static int
tiwcrex ()
{

  IWString s = "Hello world";

  IW_Regular_Expression_Test foo ("^Hello");

  if (foo.matches (s))
    ;
  else
  {
    cerr << "Did not match 'Hello world'\n";
    foo.debug_print (cerr);
    failures++;
  }

  s = "Hel";
  if (foo.matches (s))
  {
    cerr << "Matched something too short\n";
    foo.debug_print (cerr);
    failures++;
  }

  assert (foo.set_pattern ("world$"));
  s = "Hello world";
  if (foo.matches (s))
    ;
  else
  {
    cerr << "Missed match at end\n";
    foo.debug_print (cerr);
    failures++;
  }

  s = "worldq";
  if (foo.matches (s))
  {
    cerr << "Got erroneous match at end\n";
    foo.debug_print (cerr);
    failures++;
  }

  assert (foo.set_pattern ("lo wo"));
  s = "Hello world";
  if (foo.matches (s))
    ;
  else
  {
    cerr << "Missed match anywhere\n";
    foo.debug_print (cerr);
    failures++;
  }

  s = "Helloworld";
  if (foo.matches (s))
  {
    cerr << "Got erroneous match to '" << s << "'\n";
    foo.debug_print (cerr);
    failures++;
  }

  if (! foo.set_pattern ("o"))
  {
    cerr << "Cannot parse pattern 'o'\n";
    failures++;
  }

  if (! foo.matches (s))
  {
    cerr << "Cannot match single character pattern\n";
    foo.debug_print (cerr);
    failures++;
  }

//if (! foo.set_pattern ("\\([[:alnum:]]\\+\\) \\([0-9]\\{6\\}\\)"))
  if (! foo.set_pattern ("([[:alnum:]]+) ([0-9]{6})"))
  {
    cerr << "Cannot set parse complex pattern\n";
    foo.debug_print (cerr);
    failures++;
  }

  if (verbose)
    cerr << "Testing saved subexpressions\n";

  const_IWSubstring ftarget ("HELLO 123456 xyxy");
  if (! foo.matches_save_subexpressions (ftarget))
  {
    cerr << "Hmmm, the complex pattern did not match, to '" << ftarget << "'\n";
    foo.debug_print (cerr);
    failures++;
  }

  if (verbose)
    cerr << "Continuing, line " << __LINE__ << endl;

  const_IWSubstring dollar;
  foo.dollar (1, dollar);

  if ("HELLO" != dollar)
  {
    cerr << "Yipes, complex pattern $1 invalid '" << dollar << "'\n";
    foo.debug_print (cerr);
    failures++;
  }

  dollar = foo.dollar (1);
  if ("HELLO" != dollar)
  {
    cerr << "Yipes, dollar (int) seems to have failed\n";
    foo.debug_print (cerr);
    failures++;
  }

  foo.dollar (2, dollar);
  if ("123456" != dollar)
  {
    cerr << "Yipes, complex pattern $2 invalid '" << dollar << "'\n";
    foo.debug_print (cerr);
    failures++;
  }

  int i;
  if (! foo.dollar (2, i) || i != 123456)
  {
    cerr << "Extracting dollar 2 as numeric failed\n";
    failures++;
  }

  if (verbose)
    cerr << "Line " << __LINE__ << endl;

  IWString tochange = "GLERFL 987654 barf";
  foo.s (tochange, "wow");
  if ("wow wow barf" != tochange)
  {
    cerr << "Yipes, iwcrex_s failed, result '" << tochange << "'\n";
    failures++;
  }

  if (verbose)
    cerr << "Line " << __LINE__ << endl;

//if (! foo.set_pattern ("\\([0-9]\\+\\)-\\1 .*\\([[:alnum:]]\\{4\\}\\)"))
  if (! foo.set_pattern ("([0-9]+)-\\1 .*([[:alnum:]]{4})"))
  {
    cerr << "Cannot compile 2nd complex regexp\n";
    failures++;
  }

  tochange = "987-987 foo barq";

  if (verbose)
    cerr << "Line " << __LINE__ << endl;

  if (! foo.matches (tochange))
  {
    cerr << "Hmmm, doesn't match tochange\n";
    failures++;
  }

  foo.s (tochange, "wobble");
  if ("wobble-987 foo wobble" != tochange)
  {
    cerr << "Yipes, iwcrex_s failed, result '" << tochange << "'\n";
    foo.debug_print (cerr);
    failures++;
  }

  if (0 == failures)
    cerr << prog_name << " all tests successful\n";

  return failures;
}

static int
tiwcrex (IW_Regular_Expression_Test & rx,
         istream & input,
         ostream & output)
{
  int records_matching = 0;
  int records_read = 0;

  IWString buffer;
  while (buffer.getline (input))
  {
    records_read++;
    if (rx.matches (buffer))
    {
      output << buffer << endl;
      records_matching++;
    }
  }

  if (verbose)
    cerr << records_matching << " of " << records_read << " records matched\n";

  return output.good ();
}

static int
tiwcrex (const const_IWSubstring & regexp,
         istream & input,
         ostream & output)
{
  IW_Regular_Expression_Test foo;

  if (! foo.set_pattern (regexp))
  {
    cerr << "Cannot compile regular expression pattern '" << regexp << "'\n";
    return 0;
  }

  return tiwcrex (foo, input, output);
}

static int
tiwcrex (const const_IWSubstring & regexp,
         const char * fname,
         ostream & output)
{
  ifstream input;
  input.open (fname, ios::in);

  if (! input.good ())
  {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 0;
  }

  return tiwcrex (regexp, input, output);
}

static void
usage (int rc)
{
  cerr << "Usage: " << prog_name << " <options>\n";
  cerr << "  -n <ntimes>     the number of times to run each test\n";
  cerr << "  -g <regexp>     grep-like test using <regexp> on a file\n";
  cerr << "  -v              verbose output\n";

  exit (rc);
}

#include "cmdline.h"

static int
tiwcrex (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:g:");

  verbose = cl.option_count ('v');

  int ntimes = 1;

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', ntimes) || ntimes < 1)
    {
      cerr << "The -n switch requires a positive whole number\n";
      usage (4);
    }

    if (verbose)
      cerr << "Will perform " << ntimes << " iterations\n";
  }

  for (int i = 0; i < ntimes; i++)
  {
    (void) tiwcrex ();
  }

  if (cl.option_present ('g'))
  {
    if (0 == cl.number_elements ())
    {
      cerr << "The -g option (grep like test) requires a file name argument\n";
      usage (5);
    }

    const_IWSubstring r = cl.string_value ('g');

    if (! tiwcrex (r, cl[0], cout))
    {
      cerr << "Grep-like test failed\n";
    }
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

#ifdef USE_IWMALLOC
  (void) iwmalloc_initialise_memory_tracking (256);
#endif

  int rc = tiwcrex (argc, argv);

#ifdef USE_IWMALLOC
  iwmalloc_terse_malloc_status (stderr);
//iwmalloc_malloc_status (stderr);
#endif

  return rc;
}
