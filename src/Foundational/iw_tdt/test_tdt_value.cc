/*
  tester for tdt object.
  Depends on an existing TDT file. Reads the file and writes
  information to stdout. Compare that with a reference value
*/

#include <stdlib.h>

#include "cmdline.h"
#include "iwstring_data_source.h"

#include "iw_tdt.h"

const char * prog_name = NULL;

static int verbose = 0;

static IWString string_tag ("STRING<");
static IWString int_tag ("INT<");
static IWString float_tag ("FLOAT<");

static int tdts_read = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id$\n";
  cerr << "Echo's values from a TDT file\n";
  cerr << " -n             include newlines in TDT's\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
test_tdt_value (const IW_TDT & tdt,
                ostream & output)
{
  IWString s1;
  if (! tdt.dataitem_value (string_tag, s1))
  {
    cerr << "Yipes, no '" << string_tag << "' tag present\n";
    return 0;
  }

  const_IWSubstring s1a;
  if (! tdt.dataitem_value (string_tag, s1a) || s1a != s1)
  {
    cerr << "Fatal error with const_IWSubstring and tag '" << string_tag << "'\n";
    cerr << "s1 '" << s1 << "', s1a '" << s1a << "'\n";
    return 0;
  }

  const_IWSubstring ctag (string_tag);

  if (! tdt.dataitem_value (ctag, s1a) || s1a != s1)
  {
    cerr << "Fatal error with const_IWSubstring tag '" << string_tag << "'\n";
    cerr << "s1 '" << s1 << "', s1a '" << s1a << "'\n";
    return 0;
  }

  if (! tdt.dataitem_value (ctag, s1) || s1 != s1a)
  {
    cerr << "Fatal error with const_IWSubstring tag '" << string_tag << "' into const_IWSubstring\n";
    cerr << "s1 '" << s1 << "', s1a '" << s1a << "'\n";
    return 0;
  }

  if (s1.contains ('\n'))
  {
    cerr << "Newline in data, impossible\n";
    return 0;
  }
  
  else if (s1.ends_with ('\n'))
  {
    cerr << "Dataitem values should not end with newline, but '" << s1 << "' does\n";
    return 0;
  }

  output << "'" << string_tag << "' '";
  if (include_newlines_in_tdt ())
    output.write (s1.rawchars (), s1.length ());
  else
    output << s1;

  output << "'\n";

  int i;
  if (! tdt.dataitem_value (int_tag, i))
  {
    cerr << "Fatal error, no integer tag present '" << int_tag << "'\n";
    return 0;
  }

  ctag = int_tag;

  int j;
  if (! tdt.dataitem_value (ctag, j) || j != i)
  {
    cerr << "Fatal error, cannot extract int using const_IWSubstring '" << ctag << "'\n";
    cerr << " i = " << i << " j = " << j << endl;
    return 0;
  }

  output << "'" << int_tag << "' " << i << endl;

  float f1;

  if (! tdt.dataitem_value (float_tag, f1))
  {
    cerr << "Fatal error, cannot extract float '" << float_tag << "'\n";
    return 0;
  }

  ctag = float_tag;
  float f2;
  if (! tdt.dataitem_value (ctag, f2) || f1 != f2)
  {
    cerr << "Fatal error, cannot extract float with const_IWSubstring '" << ctag << "'\n";
    cerr << " f1 = " << f1 << " f2 = " << f2 << endl;
    return 0;
  }

  output << "'" << float_tag << "' " << f1 << endl;

  return output.good ();
}

static int
test_tdt_value (iwstring_data_source & input,
                ostream & output)
{
  IW_TDT tdt;

  while (tdt.next (input))
  {
    tdts_read++;

    if (! test_tdt_value (tdt, output))
    {
      cerr << "Fatal error processing tdt\n";
      cerr << tdt;
      return 0;
    }
  }

  return output.good ();
}

static int
test_tdt_value (const char * fname,
                ostream & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_tdt_value (input, output);
}

static int
test_tdt_value (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vn");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('n'))
  {
    set_include_newlines_in_tdt (1);

    if (verbose)
      cerr << "TDT's will include newlines\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! test_tdt_value (cl[i], cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "read " << tdts_read << " TDT's\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_tdt_value (argc, argv);

  return rc;
}
