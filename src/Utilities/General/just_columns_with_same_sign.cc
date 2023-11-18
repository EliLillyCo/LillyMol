/*
  We want to extract just those columns from a file where every value in
  each column is of the same sign. This is used when generating fingerprints
  and that can only process positive values.
*/

#include <stdlib.h>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iwstring.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static int header_records_to_skip = 0;

static int is_descriptor_file = 0;

static int columns_in_input = 0;

static IWString missing_value('.');

static int flush_output_after_each_molecule = 0;

/*
  This array holds the fate of each column. 0 means unknown, 1 means all positive values,
  -1 means all negative values, and a special value for discarded
*/

#define JCSS_DISCARDED -5

static int * sign_of_column = 0;

static int just_positive_columns = 0;

static int ok_not_same_column_count = 0;

static int not_same_column_count = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Extracts just columns from a file where every row is the same sign\n";
  cerr << " -s <n>         number header records to skip during determination\n";
  cerr << " -j             is a descriptor file\n";
  cerr << " -p             restrict selections to just positive values\n";
  cerr << " -M <...>       missing value token\n";
  cerr << " -X ...         more options, enter '-X help' for info\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
just_columns_with_same_sign (const_IWSubstring & buffer,
                             IWString_and_File_Descriptor& output) {
  int i = 0;
  const_IWSubstring token;

  if (is_descriptor_file)
  {
    buffer.nextword (token, i);
    output << token;
  }

  for (int col = 0; buffer.nextword (token, i); col++)
  {
    if (JCSS_DISCARDED == sign_of_column[col])
      ;
    else if (0 == sign_of_column[col])
      ;
    else
    {
      output.append_with_spacer(token);
    }
  }

  output << '\n';
  if (flush_output_after_each_molecule) {
    output.flush();
  } else {
    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
just_columns_with_same_sign (iwstring_data_source & input,
                             IWString_and_File_Descriptor& output)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (! just_columns_with_same_sign (buffer, output))
    {
      return 0;
    }
  }

  return 1;
}

static int
just_columns_with_same_sign (const char * fname,
                       IWString_and_File_Descriptor& output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return just_columns_with_same_sign (input, output);
}

static int
determine_sign_of_each_column (const const_IWSubstring & buffer)
{
  const_IWSubstring token;
  int i = 0;

  if (is_descriptor_file)
    (void) buffer.nextword (token, i);

  for (int col = 0; buffer.nextword (token, i); col++)
  {
    if (JCSS_DISCARDED == sign_of_column[col])
      continue;

    if (missing_value == token)
      continue;

//  In what may be a significant design flaw, we try to skip numeric evaluations

    if (! token.starts_with ('-'))    // positive or zero
    {
      if (1 == sign_of_column[col])
        continue;
    }
    else if (-1 == sign_of_column[col])
      continue;
    else if (just_positive_columns || 1 == sign_of_column[col])
    {
      sign_of_column[col] = JCSS_DISCARDED;
      continue;
    }
    else if (0 == sign_of_column[col])
    {
      sign_of_column[col] = -1;
      continue;
    }

    float x;
    if (! token.numeric_value (x))
    {
      cerr << "Invalid numeric token '" << token << "'\n";
      return 0;
    }
    
    if (static_cast<float> (0.0) == x)   // cannot decide anything
      ;
    else if (x > static_cast<float> (0.0))
    {
      if (1 == sign_of_column[col])
        ;
      else if (0 == sign_of_column[col])
        sign_of_column[col] = 1;
      else
        sign_of_column[col] = JCSS_DISCARDED;
    }
    else
    {
      if (-1 == sign_of_column[col])
        ;
      else if (0 == sign_of_column[col])
        sign_of_column[col] = -1;
      else
        sign_of_column[col] = JCSS_DISCARDED;
    }
  }

  return 1;
}

static int
initialise_columns_stuff ()
{
  assert (NULL == sign_of_column);
  assert (columns_in_input > 0);

  sign_of_column = new_int (columns_in_input);
  if (NULL == sign_of_column)
  {
    cerr << "Gack, cannot allocate " << columns_in_input << " data points\n";
    return 0;
  }

  if (verbose)
    cerr << "Input contains " << columns_in_input << " columns\n";

  return 1;
}

static int
determine_sign_of_each_column (iwstring_data_source & input)
{
  const_IWSubstring buffer;

  for (int i = 0; i < header_records_to_skip; i++)
  {
    input.next_record (buffer);
  }

  int records_read_this_file = 0;

  while (input.next_record (buffer))
  {
    records_read_this_file++;

    int ncol = buffer.nwords ();

    if (ncol <= 0)
    {
      cerr << "determine_sign_of_each_column:empty row, cannot process\n";
      return 0;
    }

    if (ncol == columns_in_input)
      ;
    else if (0 == columns_in_input)
    {
      columns_in_input = ncol;
      initialise_columns_stuff ();
    }
    else if (ok_not_same_column_count)
    {
      not_same_column_count++;
      if (ncol < columns_in_input)
        columns_in_input = ncol;
    }
    else
    {
      cerr << "determine_sign_of_each_column:column count mis-match, got " << ncol << " expected " << columns_in_input << endl;
      return 0;
    }

    if (! determine_sign_of_each_column (buffer))
    {
      cerr << "Fatal error processing line " << input.lines_read () << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  if (0 == records_read_this_file)
    cerr << "Warning, empty file encountered\n";

  return 1;
}

static int
determine_sign_of_each_column (const char * fname)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return determine_sign_of_each_column (input);
}

static void
DisplayDashXOptions(std::ostream& output) {
  output << "-X flush       flush output after each molecule\n";

  ::exit(0);
}

static int
just_columns_with_same_sign (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vjs:M:pX:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('j'))
  {
    is_descriptor_file = 1;
    header_records_to_skip = 1;
    if (verbose)
      cerr << "Will treat as a descriptor file\n";
  }

  if (cl.option_present ('s'))
  {
    if (! cl.value ('s', header_records_to_skip) || header_records_to_skip < 0)
    {
      cerr << "The number of header records to skip option (-s) must be a whole +ve number\n";
      usage (4);
    }

    if (verbose)
      cerr << "Will skip the first " << header_records_to_skip << " records in each file\n";
  }

  if (cl.option_present ('M'))
  {
    cl.value ('M', missing_value);

    if (verbose)
      cerr << "Missing value string set to '" << missing_value << "'\n";
  }

  if (cl.option_present ('p'))
  {
    just_positive_columns = 1;

    if (verbose)
      cerr << "Only columns where every item is +ve will be output\n";
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "flush") {
        flush_output_after_each_molecule = 1;
        if (verbose) {
          cerr << "Will flush after each molecule\n";
        }
      } else if (x == "help") {
        DisplayDashXOptions(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(cerr);
      }
    }
  }

  if (cl.empty ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  for (const char* fname : cl) {
    if (! determine_sign_of_each_column (fname))
    {
      cerr << "Fatal error profiling '" << fname << "'\n";
      return 1;
    }
  }

  int positive_columns_to_be_output = 0;
  int negative_columns_to_be_output = 0;

  for (int i = 0; i < columns_in_input; i++)
  {
    if (-1 == sign_of_column[i])
      negative_columns_to_be_output++;
    else if (1 == sign_of_column[i])
      positive_columns_to_be_output++;
  }

  if (verbose)
    cerr << "Found " << negative_columns_to_be_output << " columns of -ve numbers, and " << positive_columns_to_be_output << " columns of positive numbers\n";

  if (0 == negative_columns_to_be_output && 0 == positive_columns_to_be_output)
  {
    cerr << "None of " << columns_in_input << " columns selected for output\n";
    return 3;
  }

  IWString_and_File_Descriptor output(1);
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! just_columns_with_same_sign (cl[i], output))
    {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  output.flush();

  if (NULL != sign_of_column)
    delete sign_of_column;

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = just_columns_with_same_sign (argc, argv);

  return rc;
}
