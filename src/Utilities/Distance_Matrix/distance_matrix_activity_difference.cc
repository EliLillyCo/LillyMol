/*
  Compares distances with activity differences
*/

#include <stdlib.h>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"

#include "Utilities/Distance_Matrix/IWDistanceMatrixBase.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

const char * prog_name = NULL;

static int verbose = 0;

static IWString not_tested_string('*');

#define NOT_TESTED_VALUE -1.0e+27

static int strip_leading_zeros_from_identifiers = 0;

/*
  By default identifier in column 1, activity in column 2
*/

static int activity_column = 1;

static int skip_header_record_of_activity_file = 0;

static Accumulator<double> * acc = NULL;

static Accumulator<double> global_difference;

static int number_accumulators = 100;

static double max_distance = 1.0;

static double dx = 0.01;

static IWString_and_File_Descriptor stream_for_all_pairs;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id$\n";
  cerr << "Compares distances with activity differences\n";
  cerr << " -A <fname>     file with activities\n";
  cerr << " -c <col>       activity data in column of the -A file\n";
  cerr << " -j             skip first record of activity file\n";
  cerr << " -m <value>     minimum activity value to consider\n";
  cerr << " -z             remove leading zero's from identifiers\n";
  cerr << " -I <string>    string to designate inactive molecules\n";
  cerr << " -N <string>    string to designate molecules not tested\n";
  cerr << " -W <fname>     write all pairs\n";
  cerr << " -n <number>    number of buckets in the distance range\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

template <typename D>
int
allocate_accumulators(const IWDistanceMatrixBase<D> & dm)
{
  D m = dm.maxval();

  max_distance = static_cast<double> (m);

  dx = max_distance / static_cast<double>(number_accumulators);

  assert (number_accumulators > 0);

  acc = new Accumulator<double>[number_accumulators];

  if (NULL == acc)
  {
    cerr << "Cannot allocate " << number_accumulators << " accumulator objects\n";
    return 0;
  }

  return 1;
}

template int allocate_accumulators(const IWDistanceMatrixBase<float> &);

template <typename D>
int
write_accumulator_data(D dist,
                       const Accumulator<double> & acc,
                       std::ostream & output)
{
  if (0 == acc.n())
    return 1;

  output << dist << ' ';

  if (1 == acc.n())
  {
    output << acc.minval() << ' ' << acc.minval() << ' ' << acc.minval() << " 1 0\n";
    return 0;
  }

  output << static_cast<float> (acc.average()) << ' ' << acc.minval() << ' ' << acc.maxval() << ' ' << acc.n() << ' ' << acc.variance() << '\n';

  return 1;
}

template int write_accumulator_data(double, const Accumulator<double> &, std::ostream &);

static int
write_accumulator_data(std::ostream & output)
{
  for (int i = 0; i < number_accumulators; i++)
  {
    write_accumulator_data(i * dx, acc[i], output);
  }

  return 1;
}

template <typename D, typename A>
int 
place_into_accumulator (A activity,
                        D distance)
{
  int ndx = static_cast<int> (distance / dx + 0.4999);

  if (number_accumulators == ndx)
    ndx--;

//cerr << "distance " << distance << " dx = " << dx << " ndx " << ndx << endl;

  assert (ndx >= 0 && ndx < number_accumulators);

  acc[ndx].extra(activity);

  return 1;
}

template int place_into_accumulator(double, float);

template int place_into_accumulator(double, double);

template <typename D, typename A>
int
distance_matrix_activity_difference(const IWDistanceMatrixBase<D> & dm,
                                    int ndx,
                                    const A * activity,
                                    std::ostream & output)
{
  int n = dm.number_molecules();

  const A activity_ndx = activity[ndx];

  if (NOT_TESTED_VALUE == static_cast<A>(activity_ndx))
    return 1;

  for (int i = 0; i < n; i++)
  {
    if (i == ndx)
      continue;

    if (NOT_TESTED_VALUE == static_cast<A>(activity[i]))
      continue;

    D dist = dm.zvalue(ndx, i);

    A activity_difference = activity[ndx] - activity[i];
    if (activity_difference < static_cast<A>(0.0))
      activity_difference = - activity_difference;

    place_into_accumulator(activity_difference, dist);

    global_difference.extra(activity_difference);

    if (stream_for_all_pairs.is_open())
    {
      stream_for_all_pairs << dist << ' ' << activity_difference << '\n';
      stream_for_all_pairs.write_if_buffer_holds_more_than(32768);
    }
  }

  return 1;
}

template int distance_matrix_activity_difference(const IWDistanceMatrixBase<float> & dm,
                                    int ndx,
                                    const double * activity,
                                    std::ostream & output);

template <typename D, typename A>
int
distance_matrix_activity_difference(const IWDistanceMatrixBase<D> & dm,
                                    const A * activity,
                                    std::ostream & output)
{
  int n = dm.number_molecules();

  for (int i = 0; i < n; i++)
  {
    distance_matrix_activity_difference(dm, i, activity, output);
  }

  return 1;
}

template int distance_matrix_activity_difference(const IWDistanceMatrixBase<float> & dm,
                                    const double * activity,
                                    std::ostream & output);

template <typename T>
int
read_distance_matrix (const char * fname,
                      IWDistanceMatrixBase<T> & dm)
{
  return dm.do_read(fname);
}

template int read_distance_matrix (const char * fname,
                                   IWDistanceMatrixBase<float> & dm);

template <typename D, typename A>
int
read_activity_data_record(const const_IWSubstring & buffer,
                   IWDistanceMatrixBase<D> & dm,
                   A * activity)
{
  IWString id = buffer.word(0);

  IW_STL_Hash_Map_int::const_iterator f = dm.find(id);

  if (f == dm.IW_STL_Hash_Map_int::end())
  {
    if (strip_leading_zeros_from_identifiers)
    {
      id.remove_leading_chars('0');
      f = dm.find(id);
    }

    if (f == dm.IW_STL_Hash_Map_int::end())
    {
      cerr << "No distance data for '" << id << "'\n";
      return 0;
    }
  }

  int ndx = (*f).second;

  if (buffer.nwords() <= activity_column)
  {
    cerr << "No column " << (activity_column + 1) << " in activity data\n";
    return 0;
  }

  const_IWSubstring token = buffer.word(activity_column);

  if (not_tested_string == token)
  {
    activity[ndx] = static_cast<A>(NOT_TESTED_VALUE);
    return 1;
  }

  if (! token.numeric_value(activity[ndx]))
  {
    cerr << "Invalid numeric '" << token << "'\n";
    return 0;
  }

  return 1;
}

template int read_activity_data_record(const const_IWSubstring & buffer,
                   IWDistanceMatrixBase<float> & dm,
                   double * activity);

template <typename T, typename A>
int
read_activity_data(iwstring_data_source & input,
                   IWDistanceMatrixBase<T> & dm,
                   A * activity)
{
  input.set_dos(1);

  const_IWSubstring buffer;

  if (skip_header_record_of_activity_file)
    input.next_record(buffer);

  int records_successfully_processed = 0;

  while (input.next_record(buffer))
  {
    if (! read_activity_data_record(buffer, dm, activity))
    {
      cerr << "Cannot extract activity data, line " << input.lines_read() << endl;
      cerr << "'" << buffer << "'\n";
    }
    else
      records_successfully_processed++;
  }

  if (records_successfully_processed < dm.number_molecules())
  {
    cerr << "Incomplete activity data, distance matrix contains " << dm.number_molecules() << " items, only read " << records_successfully_processed << " activity records\n";
    if (records_successfully_processed < 2)
      return 0;
  }

  if (verbose)
    cerr << "Read " << records_successfully_processed << " valid activity records, N = " << dm.number_molecules() << endl;

  return 1;
}

template int read_activity_data(iwstring_data_source & input,
                   IWDistanceMatrixBase<float> & dm,
                   double * activity);

template <typename T, typename A>
int
read_activity_data(const const_IWSubstring & fname,
                   IWDistanceMatrixBase<T> & dm,
                   A * activity)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_activity_data (input, dm, activity);
}

template int read_activity_data(const const_IWSubstring & fname,
                   IWDistanceMatrixBase<float> & dm,
                   double * activity);

static int
distance_matrix_activity_difference (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vc:zn:N:A:jm:W:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present('z'))
  {
    strip_leading_zeros_from_identifiers = 1;
    if (verbose)
      cerr << "Leading zero's will be stripped from identifiers\n";
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', activity_column) || activity_column < 2)
    {
      cerr << "The activity column value (-c) must be a valid column number > 1\n";
      usage(3);
    }

    if (verbose)
      cerr << "Activity data in column " << activity_column << " of -A file\n";

    activity_column--;
  }

  if (cl.option_present('j'))
  {
    skip_header_record_of_activity_file = 1;
    if (verbose)
      cerr << "Will skip the first record of the activity file\n";
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', number_accumulators) || number_accumulators < 2)
    {
      cerr << "The number of distance bins must be a whole number > 1\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will divide the distance range into " << number_accumulators << " equidistant bins\n";
  }

  if (cl.option_present('N'))
  {
    cl.value('N', not_tested_string);
    if (verbose)
      cerr << "Items not tested marked as '" << not_tested_string << "'\n";
  }

  if (! cl.option_present('A'))
  {
    cerr << "Must specify the activity file via the -A option\n";
    usage(4);
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWDistanceMatrixBase<float> dm;

  if (! read_distance_matrix(cl[0], dm))
  {
    cerr << "Cannot read distance matrix '" << cl[0] << "'\n";
    return 4;
  }

  int n = dm.number_molecules();

  if (verbose)
    cerr << "Distance matrix contains data on " << n << " items\n";

  if (! allocate_accumulators(dm))
  {
    cerr << "Cannot allocate global structures\n";
    return 4;
  }

  double * activity = new double[n]; std::unique_ptr<double[]> free_activity(activity);

  set_vector(activity, n, static_cast<double>(NOT_TESTED_VALUE));

  if (cl.option_present('A'))
  {
    const_IWSubstring a = cl.string_value('A');
    if (! read_activity_data(a, dm, activity))
    {
      cerr << "Cannot read activity data from '" << a << "'\n";
      return 4;
    }
  }

  if (cl.option_present('m'))
  {
    double min_activity_to_consider;
    if (! cl.value('m', min_activity_to_consider) || min_activity_to_consider <= 0.0)
    {
      cerr << "The min activity to consider option (-m) must be a +ve real number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will ignore activities less than " << min_activity_to_consider << endl;

    int values_masked = 0;
    for (int i = 0; i < n; i++)
    {
      if (activity[i] < min_activity_to_consider)
      {
        activity[i] = NOT_TESTED_VALUE;
        values_masked++;
      }
    }

    if (verbose)
      cerr << "Masked " << values_masked << " of " << n << " values\n";

    if (values_masked >= (n - 2))
    {
      cerr << "Too many values masked " << values_masked << " of " << n << endl;
      return 5;
    }
  }

  if (cl.option_present('W'))
  {
    const char * w = cl.option_value('W');

    if (! stream_for_all_pairs.open(w))
    {
      cerr << "Cannot open file for all pairs '" << w << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "All pairs written to '" << w << "'\n";

    stream_for_all_pairs << "Dist Diff\n";
  }

  if (! distance_matrix_activity_difference(dm, activity, cout))
  {
    cerr << "Computation failed\n";
    return 4;
  }

  if (NULL != acc)
  {
    cout << "Dist ave_diff min_diff max_diff var\n";
    write_accumulator_data(cout);
  }

  if (verbose)
  {
    if (global_difference.n() > 1)
      cerr << "Whole set activity difference " << global_difference.average() << endl;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = distance_matrix_activity_difference (argc, argv);

  return rc;
}
