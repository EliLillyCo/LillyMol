#include <stdlib.h>

/*
  Tester for histogram objects
*/

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "iwhistogram.h"
#include "accumulator.h"

using std::cout;
using std::ostream;

static void
usage (int rc)
{
  cerr << "Tester for histogram objects\n";
  cerr << " -D <hist>       initialise histogram, <hist> = xmin,xmax,dx'\n";
  cerr << " -x              test the extending histogram\n";
  cerr << " -s <nskip>      number of records at top of file to skip\n";
  cerr << " -c <number>     column(s) from which to extract values\n";
  cerr << " -j              treat as a descriptor file\n";
  cerr << " -m <string>     missing value string (ignore)\n";
  cerr << " -r              ignore values out of range of the histogram\n";
  cerr << " -v              verbose output\n";

  exit (rc);
}

static int verbose = 0;

static int records_to_skip = 0;

static int column_to_process = -1;

static IWString missing_value;

typedef float mytype;

static Accumulator<mytype> accumulator;

static int values_out_of_range = 0;

static int ignore_values_out_of_range = 0;

template <typename H>
int
thistogram_token (const const_IWSubstring & token, H & d)
{
  if (missing_value == token)
    return 1;

  mytype tmp;
  if (! token.numeric_value (tmp))
  {
    cerr << "Very bad news, non-numeric token '" << token << "'\n";
    return 0;
  }

  if (verbose)
    accumulator.extra (tmp);

  if (! d.extra (tmp))
  {
    values_out_of_range++;
    if (verbose > 1)
      cerr << "Yipes, value " << tmp << " failed extra\n";

    return ignore_values_out_of_range;
  }

  return 1;
}

#ifdef __GNUG__
template int thistogram_token (const const_IWSubstring &, IWHistogram &);
template int thistogram_token (const const_IWSubstring &, Resizable_Histogram &);
#endif

template <typename H>
int
thistogram (const const_IWSubstring & buffer, H & d)
{
  const_IWSubstring token;
  if (column_to_process >= 0)
  {
    if (! buffer.word (column_to_process, token))
    {
      cerr << "Cannot extract column " << column_to_process << " from '" << buffer << "'\n";
      return 0;
    }

    return thistogram_token (token, d);
  }

// otherwise do all tokens on the line

  int i = 0;
  while (buffer.nextword (token, i))
  {
    if (! thistogram_token (token, d))
      return 0;
  }

  return 1;
}

#ifdef __GNUG__
template int thistogram (const const_IWSubstring &, IWHistogram &);
template int thistogram (const const_IWSubstring &, Resizable_Histogram &);
#endif

template <typename H>
int
thistogram (iwstring_data_source & input, H & d)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (! thistogram (buffer, d))
    {
      cerr << "Yipes, error processing line " << input.lines_read () << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

#ifdef __GNUG__
template int thistogram (iwstring_data_source &, IWHistogram &);
template int thistogram (iwstring_data_source &, Resizable_Histogram &);
#endif

template <typename H>
int
thistogram (const char * fname, H & d)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_strip_leading_blanks ();
  input.set_skip_blank_lines ();

  for (int i = 0; i < records_to_skip; i++)
  {
    const_IWSubstring notused;
    (void) input.next_record (notused);    // no check for premature EOF
  }

  return thistogram (input, d);
}

#ifdef __GNUG__
template int thistogram (const char *, IWHistogram &);
template int thistogram (const char *, Resizable_Histogram &);
#endif

template <typename H>
int
thistogram (Command_Line & cl, H & d, ostream & output)
{
  if (cl.option_present ('D'))
  {
    const_IWSubstring D = cl.string_value ('D');

    if (! d.initialise (D))
    {
      cerr << "Cannot parse the -D option '" << D << "'\n";
      return 32;
    }
  }

  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! thistogram (cl[i], d))
      return i + 1;
  }

  d.debug_print (output);

  return 1;
}

#ifdef __GNUG__
template int thistogram (Command_Line &, IWHistogram &, ostream &);
template int thistogram (Command_Line &, Resizable_Histogram &, ostream &);
#endif

static int
thistogram (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vD:xc:s:m:r");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options present\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (! cl.option_present ('D'))
  {
    cerr << "Must initialise the distribution object via the -D option\n";
    usage (13);
  }

  if (cl.option_present ('c'))
  {
    if (! cl.value ('c', column_to_process) || column_to_process < 1)
    {
      cerr << "Columns to process (-c option) must be positive whole numbers\n";
      usage (9);
    }

    if (verbose)
      cerr << "Will extract column " << column_to_process << " from the input\n";

    column_to_process--;
  }

  if (cl.option_present ('s'))
  {
    if (! cl.value ('s', records_to_skip) || records_to_skip < 0)
    {
      cerr << "The skip records (-s) option requires a whole non-negative number\n";
      usage (7);
    }

    if (verbose)
      cerr << "Will skip the first " << records_to_skip << " records in each file\n";
  }

  if (cl.option_present ('j'))
  {
    records_to_skip = 1;

    if (verbose)
      cerr << "Will process as a descriptor file\n";
  }

  if (cl.option_present ('m'))
  {
    missing_value = cl.string_value ('m');

    if (verbose)
      cerr << "Using missing value '" << missing_value << "'\n";
  }

  if (cl.option_present ('r'))
  {
    ignore_values_out_of_range = 1;

    if (verbose)
      cerr << "Will ignore values out of range of the histogram\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  int rc;
  if (cl.option_present ('x'))
  {
    Resizable_Histogram h;
    if (verbose)
      cerr << "Testing resizable histogram\n";

    rc = thistogram (cl, h, cout);
  }
  else  
  {
    IWHistogram h;
    if (verbose)
      cerr << "Testing regular histogram\n";
    rc = thistogram (cl, h, cout);
  }

  if (verbose)
  {
    cerr << "Read " << accumulator.n () << " values";
    if (accumulator.n ())
    {
      cerr << " between " << accumulator.minval () << " and " << accumulator.maxval ();
      if (accumulator.n () > 1)
        cerr << " ave " << accumulator.average ();
    }
    cerr << endl;

    if (ignore_values_out_of_range && values_out_of_range)
      cerr << values_out_of_range << " values out of range of the histogram\n";
  }

  if (0 == rc)
    return 19;

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = thistogram (argc, argv);

  return rc;
}
