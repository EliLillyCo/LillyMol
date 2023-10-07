/*
  Processes a descriptor file and removes records which contain values
  greater than N sigma from the mean of that column
*/

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fstream>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION 

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwqsort/iwqsort.h"

using std::cerr;
using std::endl;

static extending_resizable_array<int> ignore_column;
static extending_resizable_array<int> process_column;
static extending_resizable_array<int> down_column_missing_value_count;
static extending_resizable_array<int> across_row_missing_value_count;
static extending_resizable_array<int> records_discarded_for_being_too_high;
static extending_resizable_array<int> records_discarded_for_being_too_low;

/*
  Once we have profiled the columns we will be processing, we can compute their upper and
  lower cutoffs
*/

static double * lower_cutoff = NULL;
static double * upper_cutoff = NULL;

static resizable_array_p<IWString> descriptors_to_process;

static int is_descriptor_file = 0;

static IWString_and_File_Descriptor stream_for_discards;

static int only_write_selected_columns_to_discard_stream = 0;

static int verbose = 0;
static double sigma = 0.0;
static int header_records_to_skip = 0;
static char missing_value = '.';

static int records_read = 0;

static int write_to_output = 1;

const char * prog_name = NULL;

/*
  We can skip any record containing missing values
*/

static int skip_records_with_missing_values = -1;

static int rows_with_missing_values = 0;

static IWString * column_name = NULL;

static int use_median = 0;

/*
  Sept 2008.
  I want this to be able to operate on the Median. But since the
  programme is already written, and I don't want to disturb much,
  I'll write a class that can compute the median, which has 
  signatures just like an accumulator
*/

template <typename T>
class Median_Determination
{
  private:
    resizable_array<T> _values;
    Accumulator<T> _acc;

    int _median_valid;
    double _median;

  public:
    Median_Determination();

    void extra (T & t) { _values.add(t); _acc.extra(t);}

    double average () { return _acc.average();}
    double variance () { return _acc.variance();}
    T minval () const { return _acc.minval();}
    T maxval () const { return _acc.maxval();}

    double median ();
};

template <typename T>
Median_Determination<T>::Median_Determination()
{
  _median_valid = 0;
  _median = 0.0;

  return;
}

template <typename T>
class IWComparator
{
  private:
  public:
    int operator () (const T &, const T &);
};

template <typename T>
int
IWComparator<T>::operator() (const T & p1, const T & p2)
{
  if (p1 < p2)
    return -1;

  if (p1 > p2)
    return 1;

  return 0;
}

template <typename T>
double
Median_Determination<T>::median ()
{
  if (_median_valid == _acc.n())
    return _median;

  IWComparator<T> cptr;
  _values.iwqsort(cptr);

/*for (int i = 0; i < _values.number_elements(); i++)
  {
    cerr << " i = " << i << " value " << _values[i] << endl;
  }*/

  _median = _values[_values.number_elements() / 2];

  _median_valid = _acc.n();

  return _median;
}

#ifdef __GNUG
template class Median_Determination<double>;
#endif

static int
determine_columns_to_process (int columns_in_input,
                              const IWString * column_name,
                              extending_resizable_array<int> & process_column)
{
  process_column.resize (columns_in_input);

  int n = descriptors_to_process.number_elements ();

  int rc = 1;

  for (int i = 0; i < n; i++)
  {
    int found_match = 0;

    const IWString & d = *(descriptors_to_process[i]);

    for (int j = 0; j < columns_in_input; j++)
    {
      if (d != column_name[j])
        continue;

      process_column[j] = 1;
      found_match = 1;
      break;
    }

    if (! found_match)
    {
      cerr << "Could not find descripror '" << d << "'\n";
      rc = 0;
    }
  }

  for (int i = 0; i < columns_in_input; i++)
  {
    if (0 == process_column[i])
      ignore_column[i] = 1;
  }

  return rc;
}

void
compute_cutoffs (Accumulator<double> & acc,
                 double sigma,
                 double & lower_cutoff,
                 double & upper_cutoff)
{
  lower_cutoff = acc.average () - sigma * sqrt (acc.variance ());

  upper_cutoff = acc.average () + sigma * sqrt (acc.variance ());

  return;
}

void
compute_cutoffs (Median_Determination<double> & md,
                 double sigma,
                 double & lower_cutoff,
                 double & upper_cutoff)
{
  lower_cutoff = md.median () - sigma * sqrt (md.variance ());

  upper_cutoff = md.median () + sigma * sqrt (md.variance ());

  return;
}

void
write_central_value (Accumulator<double> & acc,
                     std::ostream & os)
{
  os << " ave " << acc.average();
}

void
write_central_value (Median_Determination<double> & md,
                     std::ostream & os)
{
  os << " median " << md.median();
}

static int
do_only_write_selected_columns_to_discard_stream (const const_IWSubstring & buffer,
                                                  IWString_and_File_Descriptor & stream_for_discards)
{
  int i = 0;
  int col = 0;
  const_IWSubstring token;

  while (buffer.nextword (token, i))
  {
    if (is_descriptor_file && 0 == col)
      ;
    else if (0 == process_column[col])
    {
      col++;
      continue;
    }

    stream_for_discards.append_with_spacer (token);

    col++;
  }

  stream_for_discards += '\n';

  stream_for_discards.write_if_buffer_holds_more_than(32768);

  return stream_for_discards.good ();
}

static int
initialise_file (iwstring_data_source & input,
                 extending_resizable_array<int> & process_column)
{
  input.seekg (0);

  const_IWSubstring buffer;
  for (int i = 0; i < header_records_to_skip; i++)
  {
    if (! input.next_record (buffer))
    {
      cerr << "Cannot read record " << i << " from input\n";
      return 0;
    }

    if (i > 0)
      continue;

//  cerr << "PC has " << process_column.number_elements() << " items\n";

    if (is_descriptor_file && 0 == i && descriptors_to_process.number_elements () > 0)
    {
      int columns_in_input = buffer.nwords ();

      if (columns_in_input < 2)
      {
        cerr << "Yipes, not enough columns in descriptor file\n";
        return 0;
      }

      column_name = new IWString[columns_in_input];

      if (NULL == column_name)
      {
        cerr << "Yipes, cannot allocate " << columns_in_input << " column names\n";
        return 0;
      }

      buffer.split (column_name);

      if (! determine_columns_to_process (columns_in_input, column_name, process_column))
      {
        cerr << "Cannot determine which columns to process from descriptor names\n";
        return 0;
      }
    }
  }

  return 1;
}

/*
  Decide whether or not a record should be kept or discarded.
*/

static int
rmsigma (const const_IWSubstring & buffer,
         const int columns_in_input)
{
  int missing_values_this_row = 0;

  int j = 0;
  for (int i = 0; i < columns_in_input; i++)
  {
    const_IWSubstring token;

    if (! buffer.nextword (token, j))
    {
      cerr << "Not enough columns\n";
      return 0;
    }

    if (ignore_column[i] || ! process_column[i])
      continue;

    if (missing_value == token)
    {
      missing_values_this_row++;
      if (skip_records_with_missing_values >= 0 && missing_values_this_row >= skip_records_with_missing_values)
      {
        rows_with_missing_values++;
        return 0;
      }

      continue;
    }

    double d;
    if (! token.numeric_value (d))      // should never be a problem, this is 2nd time we read the file
    {
      cerr << "Invalid numeric token '" << token << "'\n";
      return 0;
    }

//  cerr << "Value " << d << " in column " << i << ", lc " << lower_cutoff[i] << " uc " << upper_cutoff[i] << endl;

    if (d < lower_cutoff[i])
    {
      records_discarded_for_being_too_low[i]++;
      return 0;
    }

    if (d > upper_cutoff[i])
    {
      records_discarded_for_being_too_high[i]++;
      return 0;
    }
  }

// If we get to here, the record is OK

  return 1;
}

template <typename T>
int
rmsigma (iwstring_data_source & input,
         const int columns_in_input,
         T * acc,
         IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    int nw = buffer.nwords ();
    if (columns_in_input != nw)
    {
      cerr << "Column count mismatch, " << columns_in_input << " vs " <<
              nw << " line " << input.lines_read () << endl;
      return 2;
    }

    int missing_values_this_row = 0;
    int j = 0;
    for (int i = 0; i < nw; i++)
    {
      const_IWSubstring tmp;
      if (! buffer.nextword (tmp, j))
      {
        cerr << "Cannot extract word " << i << " on line " << input.lines_read () << endl;
        return 3;
      }

      if (ignore_column[i])
        continue;

      if (missing_value == tmp)
      {
        down_column_missing_value_count[i]++;
        missing_values_this_row++;
        continue;
      }

      double d;
      if (! tmp.numeric_value (d))
      {
        cerr << "Non numeric value found '" << tmp << "' line " << input.lines_read () << endl;
        return 4;
      }

      acc[i].extra (d);
    }

    across_row_missing_value_count[missing_values_this_row]++;
  }

// Now that we have profiled the input file, establish the upper and lower cutoffs

  for (int i = 0; i < columns_in_input; i++)
  {
    if (ignore_column[i] || ! process_column[i])
      continue;

    compute_cutoffs (acc[i], sigma, lower_cutoff[i], upper_cutoff[i]);
  }

  if (verbose)
  {
    for (int i = 0; i < columns_in_input; i++)
    {
      if (ignore_column[i] || ! process_column[i])
        continue;

      T & acci = acc[i];

      cerr << "Col " << (i + 1) << " min " << acci.minval () << " max " << acci.maxval ();
      write_central_value(acci, cerr);

      cerr << " [" << lower_cutoff[i] << ',' << upper_cutoff[i] << "]";

      cerr << ' ';

      if (lower_cutoff[i] > acci.minval ())
        cerr << '*';

      cerr << ',';

      if (upper_cutoff[i] < acci.maxval ())    // items definitely removed
        cerr << '*';

      cerr << endl;
    }
  }

  records_read = input.lines_read ();

  if (! input.seekg (0))
  {
    cerr << "Cannot seek back to beginning of file\n";
    return 0;
  }

  for (int i = 0; i < header_records_to_skip; i++)
  {
    input.next_record (buffer);
    if (write_to_output)
      output << buffer << '\n';

    if (! stream_for_discards.is_open())
      ;
    else if (only_write_selected_columns_to_discard_stream)
      do_only_write_selected_columns_to_discard_stream (buffer, stream_for_discards);
    else
      stream_for_discards << buffer << '\n';
  }

  output.write_if_buffer_holds_more_than(32768);

  while (input.next_record (buffer) && output.good ())
  {
    if (rmsigma (buffer, columns_in_input))
    {
      if (write_to_output)
      {
        output << buffer << '\n';
        output.write_if_buffer_holds_more_than(32768);
      }
    }
    else if (! stream_for_discards.is_open ())
      ;
    else if (only_write_selected_columns_to_discard_stream)
      do_only_write_selected_columns_to_discard_stream (buffer, stream_for_discards);
    else
      stream_for_discards << buffer << '\n';

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

template <typename T>
int 
rmsigma (iwstring_data_source & input,
         int columns_in_input,
         IWString_and_File_Descriptor & output)
{
  T * acc = new T [columns_in_input]; std::unique_ptr<T[]> free_acc(acc);

  if (NULL == acc)
  {
    cerr << "Sorry, cannot initialise " << columns_in_input << " accumulators\n";
    return 0;
  }

  return rmsigma (input, columns_in_input, acc, output);
}

static int
rmsigma (iwstring_data_source & input,
         IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  if (! input.next_record (buffer))
  {
    cerr << "Cannot read first record of input\n";
    return 1;
  }

  int columns_in_input = buffer.nwords ();
  if (0 == columns_in_input)
  {
    cerr << "First record in file contains no columns\n";
    return 4;
  }

// Reconcile the two methods of specify which columns to process

  int columns_being_processed = 0;

  for (int i = 0; i < columns_in_input; i++)
  {
    if (ignore_column[i] && process_column[i])
    {
      cerr << "Inconsistent specifications for column " << (i + 1) << " ignore " << ignore_column[i] << " process " << process_column[i] << endl;
      return 0;
    }

    if (ignore_column[i])
      process_column[i] = 0;
    else if (process_column[i])
      ignore_column[i] = 0;

    if (process_column[i])
      columns_being_processed++;
  }

  if (0 == columns_being_processed)
  {
    for (int i = is_descriptor_file; i < columns_in_input; i++)
    {
      process_column[i] = 1;
      ignore_column[i] = 0;
    }
  }

// At some stage put in a check to make sure that all columns specified
// with the -c and -i options are actually in range.

  if (process_column.number_elements () > 0)
  {
    ignore_column.resize (columns_in_input);

    for (int i = 0; i < process_column.number_elements (); i++)
    {
      if (0 == process_column[i])
        ignore_column[i]++;
    }
  }
  else if (0 == ignore_column.number_elements())
  {
    process_column.resize(columns_in_input);

    for (int i = 1; i < columns_in_input; i++)
    {
      process_column[i] = 1;
    }
  }

  if (! initialise_file (input, process_column))
  {
    cerr << "Cannot initialise input\n";
    return 0;
  }

  if (verbose > 2)
  {
    for (int i = 0; i < columns_in_input; i++)
    {
      cerr << " col " << i << " ignore " << ignore_column[i] << " proc " << process_column[i] << endl;
    }
  }

  lower_cutoff = new double[columns_in_input];
  upper_cutoff = new double[columns_in_input];

  int rc;
  if (use_median)
    rc = rmsigma<Median_Determination<double> > (input, columns_in_input, output);
  else
    rc = rmsigma<Accumulator<double> > (input, columns_in_input, output);

  output.flush();

  if (verbose)
  {
    int total_records_discarded = 0;
    for (int i = 0; i < columns_in_input; i++)
    {
      if (! process_column[i])
        continue;

      if (0 == records_discarded_for_being_too_high[i] && 0 == records_discarded_for_being_too_low[i] && 0 == down_column_missing_value_count[i])
        continue;

      cerr << "Column " << i;
      if (NULL != column_name)
        cerr << " (" << column_name[i] << ")";

      if (records_discarded_for_being_too_high[i] > 0)
      {
        cerr << ' ' << records_discarded_for_being_too_high[i] << " records too high";
        total_records_discarded += records_discarded_for_being_too_high[i];
      }

      if (records_discarded_for_being_too_low[i] > 0)
      {
        cerr << ' ' << records_discarded_for_being_too_low[i] << " records too low";
        total_records_discarded += records_discarded_for_being_too_low[i];
      }

      if (down_column_missing_value_count[i])
        cerr << ' ' << down_column_missing_value_count[i] << " missing values";

      cerr << endl;
    }

    cerr << "Read " << records_read << " records";
    if (total_records_discarded)
      cerr << ", discarded " << total_records_discarded;
    cerr << endl;

    for (int i = 0; i < across_row_missing_value_count.number_elements (); i++)
    {
       if (across_row_missing_value_count[i])
         cerr << across_row_missing_value_count[i] << " rows had " << i << " missing values\n";
    }
  }

  return rc;
}

static int
rmsigma (const char * fname,
         IWString_and_File_Descriptor & output)
{
  iwstring_data_source input (fname);
  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 1;
  }

  return rmsigma (input, output);
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << "\n";
  cerr << "Removes rows where a given descriptor falls beyond number of sigma of the mean\n";
  cerr << " -g <sigma>     specify how many sigma to trim\n";
  cerr << " -e             use offset from median rather than average\n";
  cerr << " -i <column>    specify columns to ignore\n";
  cerr << " -s <skip>      specify records at top of file to skip\n";
  cerr << " -j             process as a descriptor file (-s 1 -i 1)\n";
  cerr << " -c <column>    specify columns to process\n";
  cerr << " -d <dname>     specify descriptor(s) to process\n";
  cerr << " -D <fname>     write discarded records to <fname>\n";
  cerr << " -k             only write columns being processed to the -D file\n";
  cerr << " -q             suppress normal output\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
rmsigma (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "g:D:vs:i:jd:qkec:");

  verbose = cl.option_count ('v');

  if (! cl.option_present ('g'))
  {
    cerr << "The -g option must be present\n";
    usage (1);
  }

  if (cl.option_present ('i'))
  {
    if (cl.option_present ('c'))
    {
      cerr << "The -i and -c options are mutually exclusive\n";
      usage (2);
    }

    IWString tmp;
    int i = 0;
    while (cl.value ('i', tmp, i++))
    {
      int j;
      if (! tmp.numeric_value (j) || j < 1)
      {
        cerr << "Column specifiers (-i option) must be whole non-negative numbers\n";
        usage (3);
      }

      ignore_column[j - 1]++;
      if (verbose)
        cerr << "Will ignore column " << j << endl;
    }
  }

  if (cl.option_present ('c'))
  {
    if (cl.option_present ('i'))
    {
      cerr << "The -c and -i options are mutually exclusive\n";
      usage (2);
    }

    IWString tmp;
    int i = 0;
    while (cl.value ('c', tmp, i++))
    {
      int j;
      if (! tmp.numeric_value (j) || j < 1)
      {
        cerr << "Column specifiers (-c option) must be whole non-negative numbers\n";
        usage (3);
      }

      process_column[j - 1]++;
      if (verbose)
        cerr << "Column " << j << " will be processed\n";
    }
  }

  if (cl.option_present ('j'))
  {
    is_descriptor_file = 1;

    ignore_column[0] = 1;
    header_records_to_skip = 1;

    if (verbose)
      cerr << "Will process as a descriptor file\n";
  }

  if (! cl.value ('g', sigma) || sigma <= 0.0)
  {
    cerr << "The multiple of sigma must be a positive whole number\n";
    usage (2);
  }

  if (cl.option_present ('s'))
  {
    if (! cl.value ('s', header_records_to_skip) || header_records_to_skip <= 0)
    {
      cerr << "The -s option requires a positive whole number\n";
      usage (4);
    }
    if (verbose)
      cerr << "The first " << header_records_to_skip << " record(s) will be skipped\n";
  }

  if (cl.option_present ('d'))
  {
    int i = 0;
    const_IWSubstring d;
    while (cl.value ('d', d, i++))
    {
      IWString * tmp = new IWString (d);

      descriptors_to_process.add (tmp);
    }

    is_descriptor_file = 1;

    header_records_to_skip = 1;
    ignore_column[0] = 1;
  }

  if (cl.option_present ('k'))
  {
    if (! cl.option_present ('D'))
    {
      cerr << "The -k option only makes sense when used with the -D option\n";
    }

    only_write_selected_columns_to_discard_stream = 1;

    if (verbose)
      cerr << "Will only write selected columns to the discard stream\n";
  }

  if (cl.option_present ('q'))
  {
    if (! cl.option_present ('D'))
      cerr << "Warning, no output\n";

    write_to_output = 0;

    if (verbose)
      cerr << "Normal output suppressed\n";
  }

  if (cl.option_present('e'))
  {
    use_median = 1;

    if (verbose)
      cerr << "Will use offset from median\n";
  }

  if (cl.option_present ('D'))
  {
    IWString tmp;
    cl.value ('D', tmp);
    stream_for_discards.open (tmp.chars ());
    if (! stream_for_discards.good ())
    {
      cerr << "Cannot open '" << tmp << "' for discards\n";
      return 9;
    }
    if (verbose)
      cerr << "Discarded records written to '" << tmp << "'\n";

  }

  if (1 != cl.number_elements ())
  {
    cerr << prog_name << " processes exactly one file\n";
    usage (4);
  }

  IWString_and_File_Descriptor output(1);

  if (! rmsigma(cl[0], output))
  {
    cerr << "Error processing file '" << cl[0] << "'\n";
    return 4;
  }

  output.flush();

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = rmsigma (argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status (stderr);
#endif

  return rc;
}
