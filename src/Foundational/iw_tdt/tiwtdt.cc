/*
  Tester for TDT's
*/

#include <stdlib.h>

#include "iw_tdt.h"
#include "cmdline.h"
#include "iwcrex.h"

static int verbose = 0;

static int tdts_read = 0;

static IW_Regular_Expression iwrx;

static int
tiwtdt (IW_TDT & tdt)
{
  if (iwrx.active ())
  {
    int c = tdt.number_of_dataitems (iwrx);

    cerr << "Found " << c << " dataitems matching " << iwrx.source () << endl;
  }

  return 1;
}

static int
tiwtdt (iwstring_data_source & input)
{
  IW_TDT tdt;
  while (tdt.next (input))
  {
    if (! tiwtdt (tdt))
      return 0;
  }

  return 1;
}

static int
tiwtdt (const char * fname)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return tiwtdt (input);
}

static int
tiwtdt (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vr:");

  verbose = cl.option_count ('v');

  if (cl.option_present ('r'))
  {
    const_IWSubstring r;
    cl.value ('r', r);

    if (! iwrx.set_pattern (r))
    {
      cerr << "Cannot parse regular expression '" << r << "'\n";
      return 0;
    }
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
//  usage (1);
  }

  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! tiwtdt (cl[i]))
      return i + 1;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = tiwtdt (argc, argv);

  return rc;
}
