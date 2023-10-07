#include "stdlib.h"
#include "iostream.h"

#include "msi_object.h"
#include "cmdline.h"

const char * prog_name = nullptr;

static int verbose = 0;

static int
echomsi (iwstring_data_source & input)
{
  msi_object msi;

  if (! msi.read (input))
  {
    cerr << "Cannot read msi object\n";
    return 0;
  }

  cout << msi;

  return 1;
}

static int
echomsi (const char * fname)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "cannot open '" << fname << "' for input\n";
    return 0;
  }

  return echomsi (input);
}

static int
echomsi (Command_Line & cl)
{
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! echomsi (cl[i]))
      return 999;
  }

  return 0;
}

static int
echomsi (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "v");

  verbose = cl.option_count ('v');

  if (0 == cl.number_elements ())
  {
    cerr << "No inputs specified\n";
    return 1;
  }

  return echomsi (cl);
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = echomsi (argc, argv);

  return rc;
}
