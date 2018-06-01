#include <stdlib.h>

/*
  Tester for compiled regular expressions
*/

//#define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "iwcrex.h"
#include "iwgrep-2.5.h"
#include "iwstring_data_source.h"

typedef IW_Regular_Expression_Template<IW_grep_25_regex> IW_Regular_Expression_Test;

const char * prog_name = NULL;

static int verbose = 0;

static int failures = 0;

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

  if (! foo.set_pattern ("^FP[A-Z]*<"))
  {
    cerr << "Cannot set pattern to '^FP[A-Z]*<'\n";
    failures++;
  }

  s = "FPQQ<";
  if (! foo.matches (s))
  {
    cerr << "Did not match '" << s << "' to '" << foo.source () << "'\n";
    failures++;
  }

  if (! foo.set_pattern ("^(FOO|BAR)$"))
  {
    cerr << "Cannot set composite pattern\n";
    failures++;
  }

  s = "FOO";
  if (! foo.matches (s))
  {
    cerr << "Cannot recognised alternate '" << s << "' source '" << foo.source () << "'\n";
    failures++;
  }

  s = "BAR";
  if (! foo.matches (s))
  {
    cerr << "Cannot recognised alternate '" << s << "' source '" << foo.source () << "'\n";
    failures++;
  }

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
    cerr << "Hmmm, complex pattern did not match, to '" << ftarget << "'\n";
    foo.debug_print (cerr);
    failures++;
  }

  const_IWSubstring dollar;
  foo.dollar (1, dollar);
  if ("HELLO" != dollar)
  {
    cerr << "Yipes, complex pattern $1 invalid '" << dollar << "' from '" << ftarget << "'\n";
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

//cerr << "Line " << __LINE__ << endl;

  IWString tochange = "GLERFL 987654 barf";
  foo.s (tochange, "wow");
  if ("wow wow barf" != tochange)
  {
    cerr << "Yipes, s() failed, result '" << tochange << "'\n";
    failures++;
  }

//if (! foo.set_pattern ("\\([0-9]\\+\\)-\\1 .*\\([[:alnum:]]\\{4\\}\\)"))
  if (! foo.set_pattern ("([0-9]+)-\\1 .*([[:alnum:]]{4})"))
  {
    cerr << "Cannot compile 2nd complex regexp line " << __LINE__ << "\n";
    failures++;
  }

  tochange = "987-987 foo barq";

  if (! foo.matches (tochange))
  {
    cerr << "Yipes, regexp does not match change target, line " << __LINE__ << endl;
    cerr << "Rx '" << foo.source () << "', target '" << tochange << "'\n";
    failures++;
  }

//cerr << "Line " << __LINE__ << endl;

  if (! foo.matches (tochange))
  {
    cerr << "Hmmm, doesn't match tochange, line " << __LINE__ << "\n";
    failures++;
  }

  foo.s (tochange, "wobble");
  if ("wobble-987 foo wobble" != tochange)
  {
    cerr << "Yipes, iwcrex_s failed, result '" << tochange << "' line " << __LINE__ << "\n";
    foo.debug_print (cerr);
    failures++;
  }

  if (! foo.set_pattern ("^#"))
  {
    cerr << "Cannot compile '^#'\n";
    failures++;
  }

  if (! foo.matches ("# hello"))
  {
    cerr << "RX '" << foo.source () << "' cannot match '# hello'\n";
    failures++;
  }

  if (! foo.set_pattern ("^[a-m,3]{2}[^b]"))
  {
    cerr << "Cannot compile '^[a-m,3]{2}[^b]'\n";
    failures++;
  }

  const char * ss;

  ss = "f";

  if (foo.matches (ss))
  {
    cerr << "Erroneous match of single character '" << ss << "' with '" << foo.source () << "'\n";
    failures++;
  }

  ss = "3";
  if (foo.matches (ss))
  {
    cerr << "Erroneous match of single character '" << ss << "' with '" << foo.source () << "'\n";
    failures++;
  }

  ss = "33b";

  const_IWSubstring tmp (ss);

  if (foo.matches (tmp))
  {
    cerr << "Erroneous match '" << tmp << "' with '" << foo.source () << "'\n";
    failures++;
  }

  tmp = "ggd";

  if (! foo.matches (tmp))
  {
    cerr << "Cannot match '" << tmp << "' with '" << foo.source () << "'\n";
    failures++;
  }

  if (! foo.set_pattern ("([a-e]).{4,7}\\1"))
  {
    cerr << "Cannot compile '([a-e]).{4,7}\\1'\n";
    failures++;
  }

  if (! foo.matches ("bxxxxb"))
  {
    cerr << "'" << foo.source () << "' failed to match 'bxxxxb'\n";
    failures++;
  }

  if (foo.matches ("b...b"))
  {
    cerr << "'" << foo.source () << "' erroneous match to 'b...b'\n";
    failures++;
  }
  return failures;
}

static int
greptest (iwstring_data_source & input, IW_Regular_Expression_Test & rx)
{
  int matches = 0;

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (rx.matches (buffer))
    {
      matches++;
      if (verbose)
        cout << buffer << endl;
    }
  }

  return matches;
}

static int
greptest (const char * fname, IW_Regular_Expression_Test & rx)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  int rc = greptest (input, rx);

  cerr << rc << " matches to '" << rx.source () << "' in '" << fname << "'\n";

  return rc;
}

static void
usage (int rc)
{
  cerr << "Usage: " << prog_name << " <options>\n";
  cerr << "  -n <ntimes>     the number of times to run each test\n";
  cerr << "  -t <regexp>     search files on command line for <regexp>\n";
  cerr << "  -v              verbose output\n";

  exit (rc);
}

#include "cmdline.h"

static int
tiwcrex (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:t:");

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

  if (0 == failures)
    cerr << "All tests successful\n";

  if (cl.option_present ('t'))
  {
    const_IWSubstring t;
    cl.value ('t', t);

    IW_Regular_Expression_Test rx;
    if (! rx.set_pattern (t))
    {
      cerr << "Could not compile regular expression '" << t << "'\n";
      return 73;
    }

    for (int i = 0; i < cl.number_elements (); i++)
    {
      (void) greptest (cl[0], rx);
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
