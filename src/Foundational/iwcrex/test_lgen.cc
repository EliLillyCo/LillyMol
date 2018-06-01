#include <stdlib.h>
/*
  Tester for compiled regular expressions
*/

// #define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "iwcrex.h"
#include "iwlgen.h"

typedef IW_Regular_Expression_Template<IW_lgen> IW_Regular_Expression_Test;

const char * prog_name = NULL;

static int verbose = 0;

static int failures = 0;

static int skip_known_problems = 0;

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

  if (! foo.set_pattern ("world$"))
  {
    cerr << "Huh, cannot reset pattern\n";
    abort ();
  }

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

  if (! foo.set_pattern ("lo wo"))
  {
    cerr << "Huh, cannot reset pattern\n";
    abort ();
  }

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

#ifdef USE_IWMALLOC
  iwmalloc_check_all_malloced (stderr);
#endif

//return 0;

// Aug 97, this next step gives a memory corruption error!

  if (! skip_known_problems)
  {
    if (! foo.set_pattern ("o"))
    {
      cerr << "Cannot set pattern to 'o'\n";
      foo.debug_print (cerr);
      failures++;
    }

#ifdef USE_IWMALLOC
    iwmalloc_check_all_malloced (stderr);
#endif

    if (foo.matches (s))
    {
    }
    else
    {
      cerr << "Could not match single character pattern\n";
      foo.debug_print (cerr);
      failures++;
    }
  }

  if (0 == failures)
    cerr << prog_name << " all tests successful\n";

  return failures;
}


static void
usage (int rc)
{
  cerr << "Usage: " << prog_name << " <options>\n";
  cerr << "  -n <ntimes>     the number of times to run each test\n";
  cerr << "  -k              skip known problems\n";
  cerr << "  -v              verbose output\n";

  exit (rc);
}

#include "cmdline.h"

static int
tiwcrex (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn:k");

  verbose = cl.option_count ('v');

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

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

  if (cl.option_present ('k'))
  {
    skip_known_problems = 1;
    if (verbose)
      cerr << "Will skip know problematic cases\n";
  }

  for (int i = 0; i < ntimes; i++)
    (void) tiwcrex ();

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
  iwmalloc_malloc_status (stderr);
#endif

  return rc;
}
