/*
  Tester for the filter object
*/

#include <stdlib.h>

#include <fstream>

#include "iw_tdt.h"

typedef IW_TDT TDT;

#include "iw_tdt_filter.h"
#include "iwstring_data_source.h"
#include "cmdline.h"
#include "iwcrex.h"
#include "iwrandom.h"

static int verbose = 0;

static int tdts_read = 0;

static int tdts_matching = 0;

static int write_matching_tdts = 0;

static int tdts_to_write = -1;

static IW_TDT_Filter filter;

static IW_Regular_Expression rx_nowrite;

static ofstream stream_for_non_matches;

static float probability_of_being_written = 0.0;

static char truncate_before = '\0';

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << "\n";
  cerr << "Tester for the TDT filter\n";
  cerr << " -O <expr>    filter expression\n";
  cerr << " -w           write matching TDT's to stdout\n";
  cerr << " -N <file>    write non-matching TDT's to <file>\n";
  cerr << " -X <regexp>  if writing, don't write dataitems which match <regexp>\n";
  cerr << " -p <prob>    when a record matches, only write it with probability <prob>\n";
  cerr << " -n <number>  maximum number of TDT's to write\n";
//cerr << " -t <char>    truncate data before first <char> - useful for processing clogp\n";   // not implemented
  cerr << " -v           verbose output\n";
  cerr << endl;

  display_tdt_filter_syntax (cerr);

  exit (rc);
}

static int
tfilter (iwstring_data_source & input, 
         int output_fd)
{
  IWString output_buffer;

  TDT tdt;
  while (tdt.next (input))
  {
    tdts_read++;

    if (! filter.matches (tdt))
    {
      if (stream_for_non_matches.rdbuf ()->is_open ())
        stream_for_non_matches << tdt;

      continue;
    }

    tdts_matching++;

    if (! write_matching_tdts)
      continue;

    if (rx_nowrite.active ())
      tdt.remove_all (rx_nowrite);

    if (probability_of_being_written > 0.0)
    {
      if (iwrandom () < probability_of_being_written)
      {
        if (stream_for_non_matches.rdbuf ()->is_open ())
          stream_for_non_matches << tdt;

        continue;
      }
    }

    output_buffer << tdt;

    if (output_buffer.length () > 32768)
    {
      output_buffer.write_whole_blocks_shift_unwritten (output_fd);
//    output_buffer.resize_keep_storage (0);
    }

    if (tdts_to_write > 0 && tdts_matching >= tdts_to_write)
      break;
  }

  if (output_buffer.length ())
    output_buffer.write (output_fd);

  return 1;
}

static int
tfilter (const char * fname, 
         int output_fd)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 0;
  }

  return tfilter (input, output_fd);
}

static int
tfilter (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vO:wX:N:p:n:t:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (! cl.option_present ('O'))
  {
    cerr << "Must specify filter via the -O option\n";
    usage (3);
  }

  if (cl.option_present ('w'))
  {
    write_matching_tdts = true;
    if (verbose)
      cerr << "Matching TDT's will be written to stdout\n";
  }

  if (cl.option_present('t'))
  {
    const_IWSubstring t = cl.string_value('t');

    if (1 != t.length())
    {
      cerr << "Sorry, the -t option takes only a single character\n";
      usage(4);
    }

    if (verbose)
      cerr << "Input data truncated to first '" << t << "'\n";

    truncate_before = t[0];
  }

  if (cl.option_present ('p'))
  {
    if (! cl.value ('p', probability_of_being_written) || probability_of_being_written <= 0.0 || probability_of_being_written >= 1.0)
    {
      cerr << "Invalid probability specifier (-p option)\n";
      usage (18);
    }

    iw_random_seed ();

    if (verbose)
      cerr << "Records that match will be written with probability " << probability_of_being_written << endl;
  }

  if (cl.option_present ('X'))
  {
    if (! cl.option_present ('w'))
    {
      cerr << "The -X option (don't write dataitems) only makes sense with the -w option\n";
      usage (5);
    }

    const_IWSubstring x = cl.string_value ('X');

    if (! rx_nowrite.set_pattern (x))
    {
      cerr << "Cannot parse -X regular expression '" << x << "'\n";
      usage (3);
    }

    if (verbose)
      cerr << "Will delete TDT records matching '" << rx_nowrite.source () << "'\n";
  }

  if (cl.option_present ('O'))
  {
    const_IWSubstring f = cl.string_value ('O');

    if (! filter.build_from_string (f))
    {
      cerr << "Cannot parse tdt filter specifier '" << f << "'\n";
      return 81;
    }
  }

  if (cl.option_present ('N'))
  {
    IWString fname = cl.string_value ('N');

    stream_for_non_matches.open (fname.null_terminated_chars (), ios::out);

    if (! stream_for_non_matches.good ())
    {
      cerr << "Cannot open stream for non-matches '" << fname << "'\n";
      return 12;
    }

    if (verbose)
      cerr << "Will write non-matches to '" << fname << "'\n";
  }

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', tdts_to_write) || tdts_to_write < 1)
    {
      cerr << "The maximum number of tdt's to write (-n option) must be a +ve whole number\n";
      usage (5);
    }

    if (verbose)
      cerr << "Will write a maximum of " << tdts_to_write << " tdt's\n";

    write_matching_tdts = 1;
  }

  if (verbose)
    filter.debug_print (cerr);

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! tfilter (cl[i], 1))
    {
      return i + 1;
    }
  }

  if (verbose)
  {
    cerr << "Read " << tdts_read << " TDT's, " << tdts_matching << " matched the filter\n";
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = tfilter (argc, argv);

  return rc;
}
