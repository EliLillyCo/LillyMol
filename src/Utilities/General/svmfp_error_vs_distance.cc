/*
  We have run svmfp_evaluate.v2.rb with the -c option.
  What is the relationship between distance to nearest
  support vector and prediction error
*/

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#define ACTIVITY_DATA_IMPLEMENATION_H
#include "Foundational/iwmisc/activity_data_from_file.h"
#include "Foundational/iwmisc/misc.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static int classification = false;

static Activity_Data_From_File<float> activity;

static int identifier_column = 0;

static int predicted_activity_column = 1;

static int distance_column = 2;

static int nacc = 10;

static int fault_tolerant = 0;

static int take_absolute_value_of_activity_differences = 0;

static float min_activity_to_consider = std::numeric_limits<float>::max();

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Parses the output of svmfp_evaluate with the -c option\n";
  cerr << "This includes distance to nearest support vector.\n";
  cerr << "This programme summarises the prediction error as a function of distance\n";
  cerr << " -A <fname>     activity file\n";
  cerr << " -C             classification model - only numeric labels -1 and 1 supported\n";
  cerr << " -a             record absolute value of obs-pred for continuous response\n";
  cerr << " -m <min>       minimum obs activity to consider\n";
  cerr << " -n <buckets>   divide the 0-1 range into <buckets> ranges (def 100)\n";
  cerr << " -f             work in a fault tolerant mode - ignores empty files, bad input\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

/*
  We need two classes that accumulate the errors for a classification
  and a regression problem
*/

class Classification_Errors
{
  private:
    int _n;
    int _correct;

  public:
    Classification_Errors();

    void extra (int obs, int pred);

    int  n () const { return _n;}

    int  report (IWString_and_File_Descriptor &) const;
};

class Continuous_Errors
{
  private:
    Accumulator<float> _acc;

  public:

    void extra (float obs, float pred);

    int  n () const { return _acc.n();}

    int  report (IWString_and_File_Descriptor &) const;
};

Classification_Errors::Classification_Errors()
{
  _n = 0;

  _correct = 0;

  return;
}

void
Continuous_Errors::extra (float obs, float pred)
{
  if (obs < min_activity_to_consider)
    return;

  if (take_absolute_value_of_activity_differences)
    _acc.extra(fabs(obs - pred));
  else
    _acc.extra(obs - pred);
}

void
Classification_Errors::extra (int obs, int pred)
{
  _n++;

  if (obs == pred)
    _correct++;

  return;
}

int
Continuous_Errors::report (IWString_and_File_Descriptor & output) const
{
  output << _acc.n();

  if (0 == _acc.n())
  {
    output << " 0 0 0\n";
    return 1;
  }

  float minval = _acc.minval();

  if (fabs(minval) < 1.0e-03)
    minval = 0.0;

  if (1 == _acc.n())
    output << ' ' << minval << ' ' << _acc.maxval() << ' ' << _acc.minval();
  else 
    output << ' ' << minval << ' ' << _acc.maxval() << ' ' << static_cast<float>(_acc.average());

  output << '\n';

  return 1;
}

int
Classification_Errors::report (IWString_and_File_Descriptor & output) const
{
  output << _n << ' ';
  
  if (0 == _n)
    output << "0 0\n";
  else
    output << _correct << ' ' << static_cast<float>(_correct) / static_cast<float>(_n) << "\n";

  return 1;
}

template <typename T, typename A>
int
gather_data_record (const const_IWSubstring & buffer,
                    const Activity_Data_From_File<T> & activity,
                    A * acc)
{
  int i = 0;
  const_IWSubstring token;

  IWString id, predicted_activity_string, distance_string;

  for (int col = 0; buffer.nextword(token, i); col++)
  {
    if (identifier_column == col)
      id = token;
    else if (col == predicted_activity_column)
      predicted_activity_string = token;
    else if (col == distance_column)
      distance_string = token;
  }

//cerr << "Brom'" << buffer << "' get '" << id << "'\n";

  if (0 == id.length() || 0 == predicted_activity_string.length() || 0 == distance_string.length())
  {
    cerr << "Incomplete data!\n";
    return 0;
  }

  T obs;

  if (! activity.get_activity(id, obs))
  {
    cerr << "No experimental data for '" << id << "'\n";
    return 0;
  }

  T pred;

  if (! predicted_activity_string.numeric_value(pred))
  {
    cerr << "Invalid activity value\n";
    return 0;
  }

  float d;

  if (! distance_string.numeric_value(d) || d < 0.0 || d > 1.0)
  {
    cerr << "INvalid distance '" << distance_string << "'\n";
    return 0;
  }

  int ndx = static_cast<int>(nacc * d);

//cerr << id << " obs " << obs << " pred " << pred << " * " << (obs == pred) << " " << d << endl;

  acc[ndx].extra(obs, pred);

  return 1;
}

template <typename T, typename A>
int
gather_data (iwstring_data_source & input,
             const Activity_Data_From_File<T> & activity,
             A * acc)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read header record\n";
    return 0;
  }

  while (input.next_record(buffer))
  {
    if (! gather_data_record (buffer, activity, acc))
    {
      cerr << "Error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return input.lines_read() > 1;   // header record does not count
}

template <typename T, typename A>
int
gather_data (const char * fname,
             const Activity_Data_From_File<T> & activity,
             A * acc)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" <<fname << "'\n";
    return 0;
  }

  return gather_data (input, activity, acc);
}

template <typename T, typename A>
int
svmfp_error_vs_distance (const Command_Line & cl,
                         Activity_Data_From_File<T> & activity,
                         A * acc,
                         IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (dash_s(cl[i]))
      ;
    else if (fault_tolerant)
    {
      cerr << "Skipping empty file '" << cl[i] << "'\n";
      continue;
    }

    if (gather_data(cl[i], activity, acc))
      ;
    else if (fault_tolerant)
    {
      cerr << "Skipping file with problems '" << cl[i] << "'\n";
      continue;
    }
    else
    {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      return 0;
    }
  }

  int n = 0;

  for (int i = 0; i <= nacc; i++)
  {
    n += acc[i].n();
  }

  if (0 == n)
  {
    cerr << "No data available\n";
    return 0;
  }

  if (verbose)
    cerr << "Read data on " << n << " observations\n";

  if (classification)
    output << "Dist N correct fraction\n";
  else
    output << "Dist N min max ave\n";

  for (int i = 0; i <= nacc; i++)
  {
    float d = static_cast<float>(i) / static_cast<float>(nacc);

    output << d << ' ';
    acc[i].report(output);
  }

  return 1;
}

static int
do_classification_models (Command_Line & cl,
                          IWString_and_File_Descriptor & output)
{
  Activity_Data_From_File<int> activity;

  if (! activity.construct_from_command_line(cl, 'A', verbose))
  {
    cerr << "Cannot read activity data (-A)\n";
    return 0;
  }

  Classification_Errors * acc = new Classification_Errors[nacc + 1];

  return svmfp_error_vs_distance(cl, activity, acc, output);
}

static int
do_regression_models (Command_Line & cl,
                      IWString_and_File_Descriptor & output)
{
  Activity_Data_From_File<float> activity;

  if (! activity.construct_from_command_line(cl, 'A', verbose))
  {
    cerr << "Cannot read activity data (-A)\n";
    return 0;
  }

  Continuous_Errors * acc = new Continuous_Errors[nacc + 1];

  return svmfp_error_vs_distance (cl, activity, acc, output);
}

static int
svmfp_error_vs_distance (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:Cn:fam:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('A'))
  {
    cerr << "Must specify activity data via the -A option\n";
    usage(3);
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', nacc) || nacc < 2)
    {
      cerr << "The number of distance buckets (-n) must be a valid +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will report distances in " << nacc << " buckets\n";
  }

  if (cl.option_present('f'))
  {
    fault_tolerant = 1;

    if (verbose)
      cerr << "Will ignore most errors\n";
  }

  if (cl.option_present('a'))
  {
    take_absolute_value_of_activity_differences = 1;

    if (verbose)
      cerr << "Will record the absolute error for continuous response\n";
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', min_activity_to_consider))
    {
      cerr << "The min obs activity to consider (-m) must be a valid float\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will only consider items where obs activity > " << min_activity_to_consider << endl;
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  if (cl.option_present('C'))
  {
    classification = 1;

    if (verbose)
      cerr << "Data is classification data\n";

    predicted_activity_column = 2;
    distance_column = 3;

    do_classification_models (cl, output);
  }
  else
    do_regression_models (cl, output);

  output.flush();

  if (verbose)
  {
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = svmfp_error_vs_distance(argc, argv);

  return rc;
}
