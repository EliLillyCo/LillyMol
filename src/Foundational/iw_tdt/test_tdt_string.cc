/*
  Tester for TDT fromstring
*/

#include <stdlib.h>

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "iw_auto_array.h"

#include "iw_tdt.h"

const char * prog_name = NULL;

static int verbose = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Tester for tdt_from_string functionality\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
test_tdt_string (iwstring_data_source & input,
                 ostream & output)
{
  IW_Regular_Expression vbar("^\\|");

  int ntdt = input.grep(vbar);

  if (ntdt < 1)
  {
    cerr << "No TDTs in input\n";
    return 0;
  }

  if (verbose)
    cerr << "Input contains " << ntdt << " TDT's\n";

  IW_TDT * tdt = new IW_TDT[ntdt]; iw_auto_array<IW_TDT> free_tdt(tdt);

  for (int i = 0; i < ntdt; i++)
  {
    if (! tdt[i].next(input))
    {
      cerr << "Cannot build tdt " << i << " from file descriptor\n";
      return 0;
    }
  }

  if (! input.seekg(0))
  {
    cerr << "Cannot seek back to start of file\n";
    return 0;
  }

  int ndx = 0;

  IWString concatenated_records;

  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    concatenated_records << buffer << '\n';
    if ('|' != buffer)
      continue;

//  cerr << "Just read '" << buffer << "'\n";
    IW_TDT j;
    if (! j.build(concatenated_records))
    {
      cerr << "Cannot build tdt from string\n";
      cerr << concatenated_records;
      return 0;
    }

    if (tdt[ndx] == j)   // great
      ;
    else
    {
      cerr << "Mismatch on tdt " << ndx << endl;
      cerr << concatenated_records;
      return 0;
    }

    ndx++;
    concatenated_records.resize_keep_storage(0);
  }

  return 1;
}

static int
test_tdt_string (const char * fname,
                 ostream & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_tdt_string(input, output);
}


static int
test_tdt_string (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "v");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! test_tdt_string(cl[i], cout))
    {
      cerr << "Test failure processing '" << cl[i] << "'\n";
      rc = i + 1;
      break;
    }
  }

  if (0 == rc)
    cerr << "All tests successful\n";

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_tdt_string(argc, argv);

  return rc;
}
