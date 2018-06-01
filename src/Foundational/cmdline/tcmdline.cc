#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"

#include "cmdline.h"

const char * prog_name;

static int
tcmdline (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "ab:c:F:g");

  if (cl.option_present ('a'))
  {
    cerr << "Option a present\n";
    if (NULL == cl.option_value ('a'))
      cerr << "value is correct\n";
    else
      cerr << "OOps, option value a is not null!\n";
  }
  else
    cerr << "OOps, option a is missing\n";

  IWString a = cl.string_value ('a');

  cerr << "A option value is '" << a << "'\n";

  if (! cl.option_present ('b'))
    cerr << "OOPS, option b not present\n";
  else
  {
    if (NULL != cl.option_value ('b'))
      cerr << "Option 'b' value is '" << cl.option_value ('b') << "'\n";
    else
      cerr << "Option 'b' is null\n";
  }

  IWString b = cl.string_value ('b');

  cerr << "B option value is '" << b << "'\n";

  if (cl.option_present ('F'))
  {
    int i = 0;
    const_IWSubstring f;
    while (cl.value ('F', f, i++))
    {
      cerr << "-F option is '" << f << "'\n";
    }
  }

  cerr << "Found " << cl.number_elements () << " other arguments\n";
  for (int i = 0; i < cl.number_elements (); i++)
    cerr << "Arg " << i << " is '" << cl[i] << "'\n";

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = tcmdline (argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status (stderr);
#endif

  return rc;
}
