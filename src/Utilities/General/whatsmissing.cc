/*
  We often want to know what's missing from a descriptor file
*/

#include <stdlib.h>
#include <limits>
using namespace std;

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "misc.h"

static int verbose = 0;

static int records_read = 0;

/*
  Include spaces in the missing value string so we can find it
*/

static IWString missing_value_string (" . ");
static IWString missing_value_end_of_record (" .");
static IWString missing_value ('.');

static int max_ok_missing_values = -1;

static int work_quietly = 0;

static IWString_and_File_Descriptor stream_for_discarded_records;

/*
  By default we write the header to the missing value stream, but maybe
  we don't always want that.
*/

static int write_header_to_missing_record_stream = 1;

/*
  We keep the header of the input file here at file scope to make
  processing the -M file easier
*/

static IWString header;

static int records_written = 0;

static int header_written_this_file = 0;
static int header_written_to_missing_values_stream_this_file = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Identifies which descriptors are missing from a descriptor file\n";
  cerr << " -m <string>    missing value (default '" << missing_value << "')\n";
  cerr << " -x <int>       maximum number of times to report each column\n";
  cerr << " -b             break after finding first difference\n";
  cerr << " -n <number>    work as a filter. Records having <number> or fewer\n";
  cerr << "                missing values written to stdout\n";
  cerr << " -M <file>      option file for records discarded by -n option\n";
  cerr << "                to append to an existing file use '>>file'\n";
  cerr << " -q             quiet mode (only recognised with -n)\n";
  cerr << " -u             DO NOT write header record to -M file\n";
  cerr << " -r             insist upon rectangular files - enforce column count check\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static IWString * descriptor_names;
static int ncols = 0;

static int allow_non_rectangular_files = 1;

static int * items_missing_in_column = NULL;

static int * items_missing_per_id = NULL;

/*
  To cut down the amount of stuff coming to the screen, we can limit the number
  of times any descriptor is reported
*/

static int max_report = numeric_limits<int>::max();

static int break_after_finding_first_missing_value = 0;

static int records_with_missing_values = 0;

static int
write_rejected_record (const const_IWSubstring & buffer,
                       IWString_and_File_Descriptor & os)
{
  if (! write_header_to_missing_record_stream)    // don't worry about the header
    ;
  else if (! header_written_to_missing_values_stream_this_file)
  {
    os << header << '\n';
    header_written_to_missing_values_stream_this_file = 1;
  }

  os << buffer << '\n';

  os.write_if_buffer_holds_more_than(32768);

  return os.good ();
}

static int
write_ok_record (const const_IWSubstring & buffer,
                 IWString_and_File_Descriptor & output)
{
  if (! header_written_this_file)
  {
    output << header << '\n';
    header_written_this_file = 1;
  }

  records_written++;
  output << buffer << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return output.good ();
}

static int
whatsmissing (const const_IWSubstring & buffer,
              int line_number,
              int & missing_values_this_id,
              IWString_and_File_Descriptor & output)
{
  missing_values_this_id = 0;

  if (buffer.contains (missing_value_string))
    ;
  else if (buffer.ends_with (missing_value_end_of_record))
    ;
  else
    return 1;

  int i = 0;
  const_IWSubstring id;
  const_IWSubstring token;
  int col = 0;

  while (buffer.nextword (token, i))
  {
    if (0 == col)
      id = token;
    else if (missing_value != token)
      ;
    else
    {
      missing_values_this_id++;
      items_missing_in_column[col]++;
      if (verbose < 2)
        ;
      else if (items_missing_in_column[col] > max_report)
        ;
      else
        output << id << ", line " << line_number << " missing '" << descriptor_names[col] << "', column " << (col + 1) << '\n';
    }
    
    col++;
  }

  items_missing_per_id[missing_values_this_id]++;

  if (max_ok_missing_values < 0)
    output << missing_values_this_id << " of " << (ncols - 1) << " descriptors missing for '" << id << "' line " << line_number << '\n';

//cerr << "After '" << buffer << "', col = " << col << " expected " << ncols << endl;
  if (col == ncols)
    ;
  else if (allow_non_rectangular_files)
    ;
  else
  {
    cerr << "Column count mismatch, got " << col << " expected " << ncols << ", line " << line_number << '\n';
    cerr << buffer << '\n';
    return 0;
  }

  if (missing_values_this_id)
    records_with_missing_values++;

  return output.good ();
}

static int
whatsmissing (iwstring_data_source & input, 
              IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer) && output.good ())
  {
    records_read++;

    int missing_values_this_id = 0;

    if (! whatsmissing (buffer, input.lines_read (), missing_values_this_id, output))   // some kind of fatal error
    {
      return 0;
    }

    if (max_ok_missing_values < 0)     // not functioning as a filter
      ;
    else if (missing_values_this_id <= max_ok_missing_values)    // good record
      write_ok_record (buffer, output);
    else if (stream_for_discarded_records.is_open ())   // too many missing values
      write_rejected_record (buffer, stream_for_discarded_records);

    if (break_after_finding_first_missing_value && records_with_missing_values)
      return 1;
  }

  return output.good ();
}

static int
establish_column_identities (const const_IWSubstring & header)
{
  ncols = header.nwords ();

  if (ncols <= 0)
  {
    cerr << "Yipes, only " << ncols << " columns in the header\n";
    return 0;
  }

  if (verbose)
    cerr << "File has " << ncols << " columns\n";

  descriptor_names = new IWString[ncols];
  items_missing_in_column = new_int (ncols);
  items_missing_per_id = new_int (ncols);

  const_IWSubstring d;
  int i = 0;
  int col = 0;
  while (header.nextword (d, i))
  {
    descriptor_names[col] = d;
    col++;
  }

  if (verbose)
  {
    if (verbose > 2)
    {
      for (int i = 1; i < ncols; i++)
      {
        cerr << " column " << (i + 1) << " is '" << descriptor_names[i] << "'\n";
      }
    }
  }

  return 1;
}

static int
report_overall_column_statistics (ostream & os,
                        int ncols,
                        const int * items_missing_in_column)
{
  for (int i = 0; i < ncols; i++)
  {
    if (items_missing_in_column[i])
      os << " descriptor '" << descriptor_names[i] << "' (column " << (i + 1) << ") missing on " << items_missing_in_column[i] << " records\n";
  }

  return os.good ();
}

static int
whatsmissing (const char * fname, 
              IWString_and_File_Descriptor & output)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (! input.next_record (header))
  {
    cerr << "Cannot read header record from '" << fname << "'\n";
    return 0;
  }

  header_written_this_file = 0;
  header_written_to_missing_values_stream_this_file = 0;

  if (! establish_column_identities (header))
  {
    cerr << "Cannot determine descriptor names, file '" << fname << "'\n";
    return 0;
  }

  int rc = whatsmissing (input, output);

  output.flush();

  if (work_quietly)
    ;
  else 
    report_overall_column_statistics (cerr, ncols, items_missing_in_column);

  if (verbose)
  {
    cerr << "Read " << records_read << " records from '" << fname << "'\n";
    for (int i = 0; i < ncols; i++)
    {
      if (items_missing_per_id[i])
        cerr << items_missing_per_id[i] << " records had " << i << " missing values\n";
    }
  }

  delete [] descriptor_names;

  delete [] items_missing_in_column;
  delete [] items_missing_per_id;

  ncols = 0;

  return rc;
}

static int
whatsmissing (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vx:m:bn:M:qur");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised_options_encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('m'))
  {
    missing_value = cl.string_value ('m');

    missing_value_string = ' ';
    missing_value_string += missing_value;

    missing_value_end_of_record = missing_value_string;
    missing_value_string += ' ';

    if (verbose)
      cerr << "MIssing value set to '" << missing_value << "'\n";
  }

  if (cl.option_present ('x'))
  {
    if (! cl.value ('x', max_report) || max_report < 1)
    {
      cerr << "The max report option (-x) must be a whole positive number\n";
      usage (8);
    }

    if (verbose)
      cerr << "Will report a maximum of " << max_report << " missing values in each column\n";
  }

  if (cl.option_present ('b'))
  {
    break_after_finding_first_missing_value = 1;

    if (verbose)
      cerr << "Will break after finding first missing value\n";
  }

  if (cl.option_present('r'))
  {
    allow_non_rectangular_files = 0;

    if (verbose)
      cerr << "All files must be rectangular\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', max_ok_missing_values) || max_ok_missing_values < 0)
    {
      cerr << "The maximum number of missing values to tolerate (-n option) must be non-negative\n";
      usage (6);
    }

    if (verbose)
    {
      cerr << "Input file(s) will be echo'd, records having more than " << max_ok_missing_values << " missing values discarded\n";
    }

    if (cl.option_present ('q'))
    {
      work_quietly = 1;

      if (verbose)
        cerr << "Working as filter, will work quietly\n";
    }
  }

  if (cl.option_present ('M'))
  {
    if (max_ok_missing_values < 0)
    {
      max_ok_missing_values = 0;
      cerr << "By default, records with any missing values discarded\n";
    }

    IWString m = cl.string_value ('M');

    if (! stream_for_discarded_records.open (m.null_terminated_chars()))
    {
      cerr << "Cannot open stream for discarded records '" << m << "'\n";
      return 7;
    }

    if (verbose)
      cerr << "Discarded records written to '" << m << "'\n";

    if (cl.option_present ('u'))
    {
      write_header_to_missing_record_stream = 0;

      if (verbose)
        cerr << "Will not write header records to missing record stream (-M)\n";
    }
  }

  if (break_after_finding_first_missing_value && max_ok_missing_values >= 0)
  {
    cerr << "The -b option is mutually exclusive of filtering\n";
    usage (5);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  int nfiles = cl.number_elements ();
  for (int i = 0; i < nfiles; i++)
  {
    if (! whatsmissing (cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose && nfiles > 1)
  {
    cerr << "Read " << records_read << " records from " << nfiles << " files, " << records_with_missing_values << " records had missing values\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = whatsmissing (argc, argv);

  return rc;
}
