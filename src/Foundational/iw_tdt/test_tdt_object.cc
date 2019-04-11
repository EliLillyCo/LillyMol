/*
  Tester for TDT objects
*/

#include <stdlib.h>

#include "cmdline.h"
#include "iwcrex.h"

#include "iw_tdt.h"

const char * prog_name = NULL;

static int verbose = 0;

static int tdts_read = 0;

static int echo_tdts = 0;

static IWString tag_to_remove;

static IW_Regular_Expression rx_to_remove;

static int dataitems_removed = 0;

static extending_resizable_array<int> items_in_tdt;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "test_tdt_object.cc,v 1.1.1.1 2003/09/25 12:35:55 \n";
  cerr << "Tester for TDT objects\n";
  cerr << " -x <tag>       remove all dataitems with <tag>\n";
  cerr << " -X <tag>       remove all dataitems whose tags match regular expression <tag>\n";
  cerr << " -e             echo the input\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

typedef IW_TDT TDT;

static int
do_remove_tag (TDT & tdt,
               IW_Regular_Expression & rx_to_remove)
{
  if (0 == tdt.remove_all (rx_to_remove))
  {
    cerr << "No match to '" << rx_to_remove.source () << "' in\n";
    cerr << tdt;
    return 0;
  }

  return 1;
}

static int
do_remove_tag (TDT & tdt,
               const IWString & tag_to_remove)
{
  int ndx = tdt.index_of_dataitem (tag_to_remove);

  if (ndx < 0)
  {
    cerr << "Cannot find '" << tag_to_remove << "' in\n";
    cerr << tdt;
    return 0;
  }

  tdt.remove_item (ndx);

  return 1;
}

static int
test_tdt_object (iwstring_data_source & input,
                 ostream & output)
{
  TDT tdt;

  while (tdt.next (input))
  {
    tdts_read++;

    int n = tdt.number_elements ();

    items_in_tdt[n]++;

    if (tag_to_remove.length ())
      do_remove_tag (tdt, tag_to_remove);

    if (rx_to_remove.active ())
      do_remove_tag (tdt, rx_to_remove);

    if (echo_tdts)
      output << tdt;
  }

  return output.good ();
}

static int
test_tdt_object (const char * fname,
                 ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_tdt_object (input, output);
}

static int
test_tdt_object (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vex:X:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('e'))
  {
    echo_tdts = 1;

    if (verbose)
      cerr << "Will echo input TDT's\n";
  }

  if (cl.option_present ('x'))
  {
    tag_to_remove = cl.string_value ('x');

    echo_tdts = 1;

    if (verbose)
      cerr << "Will remove tag '" << tag_to_remove << "'\n";

    if (! tag_to_remove.ends_with ('<'))
      tag_to_remove += '<';
  }

  if (cl.option_present ('X'))
  {
    const_IWSubstring x = cl.string_value ('X');

    if (! rx_to_remove.set_pattern (x))
    {
      cerr << "Possibly invalid regular expression '" << x << "'\n";
      return 9;
    }

    echo_tdts = 1;

    if (verbose)
      cerr << "Will remove all tags that match '" << x << "'\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! test_tdt_object (cl[i], cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << tdts_read << " TDT's\n";

    for (int i = 0; i < items_in_tdt.number_elements (); i++)
    {
      if (items_in_tdt[i])
        cerr << items_in_tdt[i] << " TDT's had " << i << " items\n";
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_tdt_object (argc, argv);

  return rc;
}
