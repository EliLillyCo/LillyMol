/*
  Evaluate a training set of near neighbours and experimental
  values in order to determine optimum values for nn_predictions.
*/

#include <stdlib.h>
#include <math.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <unordered_map>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/iwmisc/misc.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "Utilities/GFP_Tools/bsquared.h"

#include "nn_results_.h"
#include "prediction_weight_function.h"
#include "nn_stuff.h"
#include "nn_specification.h"
#include "iwpvalue.h"
#include "extract_from_tdt_form.h"

using std::cerr;
using std::endl;
using std::dec;
using std::hex;

static int verbose = 0;

/*
  If we are to invoke gfp_predictions, we need to know where it likely is...
*/

static int experimental_results_identifier_column = 0;
static int experimental_results_experimental_value_column =  1;

static IWString missing_value ('.');

/*
  When doing a paired T test, we need a threshold. Should really be part of
  the Cross_Validation_Conditions object, but one of the qsort functions needs
  the value
*/

static double pvalue_threshold = 0.0;

/*
  Once the P value is known, we can then compute the corresponding cutoff value for T
*/

static double tvalue_threshold = 0.0;

/*
  We can control how many predictions are printed after cross validation
*/

static int prediction_types_to_print = 10;


/*
  We can also decide that two B squared values are indistinguishable if they differ by less
  than a certain amount
*/

static double min_meaningful_bsquared_difference = 0.0;

/*
  We may find that something like 0.03-2N is the highest ranked model,
  but it isn't terribly useful - a very complex model.  Therefore we
  can impose a minimum useful distance
*/

static double shortest_useful_distance = -1.0;

/*
  It is instructive to see how the careful comparisons are done
*/

static int careful_comparisons_done = 0;
static int careful_comparisons_resolved_by_complexity = 0;

/*
  We can help automation by writing best results to a file
*/

static std::ofstream stream_for_best_models;

/*
  if we are going to make predictions, we can change the number of predictions here
*/

static int number_predictions_to_print = 1;

static int outlier_neighbours_to_print = 0;

/*
  When writing the outlier data, if all items in the dataset have the same number
  of tokens per identifier, then we don't need to translate spaces to underscores
  when we write the names
*/

static int tokens_per_id = 0;

typedef float experimental_value_t;

/*
  Sometimes we need to strip leading zero's in order to match up experimental data
*/

static int strip_leading_zeros_from_identifiers = 0;

typedef int category_type_t;

typedef double probability_type_t;

static IWString bsquared_or_error_rate;       // either 'BSquared' or 'ERROR'

static int classification_model = 0;

/*
  Sometimes categories are related, like good, better and best
*/

static int good_better_best_model = 0;

/*
  On classification models, there are some global items that are generally helpful
*/

static int number_categories = 0;
static int highest_category_number = 0;

static int * number_in_category = nullptr;

/*
  When doing a classification model, we just sum the weights of
  each neighbour from each category
*/

static weight_t * category_sumw = nullptr;

/*
  In computing the overall score, each category can have a different weight assigned
*/

static weight_t * category_weight = nullptr;

/*
  Sometimes we just want to maximise the classification accuracy of just
  one particular class - typically when there are a very small number
  of one class present in a larger dataset
*/

static int overall_error_rate_is_error_rate_for_particular_class = -1;

/*
  It is often useful to write out the test set predictions
*/

static IWString file_name_for_test_set_predictions;

static int nn_file_is_distance_matrix = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Builds KNN models. Optionally enter pre-computed training set memberships as arguments\n";
  display_standard_nn_stuff (cerr);
  cerr << " -N <condition>   specify neighbour conditions to test\n";
  cerr << " -N file=<fname>  file of conditions\n";
  cerr << " -M <string>      missing value specifier (default " << missing_value << ")\n";
  cerr << " -E file=<fname>  file with training set experimental data\n";
  cerr << " -E col=<col>     experimental data is always column <col> in the training set names\n";
  cerr << "                  'col=last' for the last column in the identifier\n";
  cerr << " -E YYY=          experimental data is a token in the name starting with 'YYY='\n";
  cerr << " -j               skip first line in experimental data file\n";
  cerr << " -C ...           classification model, 'def' for default class weights\n";
  cerr << " -C ...           class weight specification, '-C 0=0.1,1=3.0,2=0.5,...'\n";
  cerr << " -L ...           error balancing specifications\n";
  cerr << " -V ...           options for cross validation. Enter '-V help' for details\n";
  cerr << " -g <value>       B squared values <= <value> are ignored\n";
  cerr << " -p <number>      number of predictions to print\n";
  cerr << " -P <fname>       write leave-one-out training set predicted values to <fname>\n";
  cerr << " -O <fname>       write molecules that are poorly predicted to <fname>\n";
  cerr << " -R <fname>       write model descriptions to <fname> - can be used by nn_predictions\n";
  cerr << " -u <dist>        smallest useful distance. All distances specified will be computed\n";
  cerr << "                  but the one printed will have a distance of at least <dist>\n";
  cerr << " -z               strip leading zero's from identifiers before matching data\n";
  cerr << " -d               input file is a distance matrix\n";
  cerr << " -D max=<dist>    specify maximum distance in input (default 1.0)\n";
  cerr << " -D nwt=<num>     specify number of points in the weight arrays (default 101)\n";
  cerr << " -J ...           various other options, enter '-J help' for info\n";
  cerr << " -v               verbose output\n";

  exit (rc);
}

/*
  We want to be able to do the test/train splits either randomly, or from test/train splits
  specified on the command line
*/

class Training_Sets_From_Files : public resizable_array_p<const_IWSubstring>
{
  private:

//  private functions

    int _external_split (const const_IWSubstring &, IWString_STL_Hash_Set &) const;
    int _external_split (iwstring_data_source & input, IWString_STL_Hash_Set & zresult) const;

  public:
    Training_Sets_From_Files ();

    int initialise_external_splits (const Command_Line &);

    int external_split (int, IWString_STL_Hash_Set &) const;
};

Training_Sets_From_Files::Training_Sets_From_Files ()
{
  return;
}

/*
  All arguments after the first one are file names of training set memberships
*/

int
Training_Sets_From_Files::initialise_external_splits (const Command_Line & cl)
{
  int n = cl.number_elements ();

  if (n < 2)
  {
    cerr << "Training_Sets_From_Files::initialise_external_splits: must specify one or more split files on the command line\n";
    return 0;
  }

  if (verbose)
    cerr << n << " test/train assignments being read from files\n";

  resize (n - 1);

  for (int i = 1; i < n; i++)
  {
    const char * fname = cl[i];

    if (! dash_s (fname))
    {
      cerr << "Training_Sets_From_Files::missing or empty file '" << fname << "'\n";
      return 0;
    }

    const_IWSubstring * tmp = new const_IWSubstring (fname);   // cl stays in scope!

    add (tmp);
  }

  return n - 1;
}

int
Training_Sets_From_Files::external_split (int ndx,
                                          IWString_STL_Hash_Set & zresult) const
{
  if (! ok_index (ndx))
    return 0;

  zresult.clear();

  int rc = _external_split (*(_things[ndx]), zresult);

  if (0 == rc)
    cerr << "Training_Sets_From_Files::external_split: bad news, cannot read '" << (*(_things[ndx])) << "'\n";

  return rc;
}

int
Training_Sets_From_Files::_external_split (const const_IWSubstring & fname,
                                           IWString_STL_Hash_Set & zresult) const
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Training_Sets_From_Files::_external_split: canot open '" << fname << "'\n";
    return 0;
  }

  return _external_split (input, zresult);
}

int
Training_Sets_From_Files::_external_split (iwstring_data_source & input,
                                           IWString_STL_Hash_Set & zresult) const
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    buffer.strip_leading_blanks ();
    buffer.strip_trailing_blanks ();

    if (buffer.starts_with ('#'))
      continue;

    if (strip_leading_zeros_from_identifiers)
      buffer.remove_leading_chars ('0');

    buffer.truncate_at_first (' ');

    zresult.insert (buffer);
  }

  return zresult.size ();
}

/*
  May 2002. When I introduced the distance matrix, I noticed that I got different answers.
  Turns out that the problem is with the ordering of neighbours. Imagine we have 3 neighbours,
  ID    DIST    ACTIVITY
  1     0.32    9.5
  2     0.32    5.5
  3     0.32    1.0

  In any sorting by distance, the ordering of these neighbours is arbitrary. Therefore, during
  cross validation, we offer the option of randomly shuffling the order of equi-distant
  neighbours
*/


class Cross_Validation_Conditions : public Training_Sets_From_Files
{
  private:
    int _percent_training_set;
    int _number_splits;
    int _stratified_sampling;
    int _do_paired_t_test;
    int _reorder_equidistant_neighbours;

  public:
    Cross_Validation_Conditions ();

    int construct_from_command_line (Command_Line & cl, char c);

    int  percent_training_set () const { return _percent_training_set;}
    void set_percent_training_set (int s) { _percent_training_set = s;}

    int  number_splits () const { return _number_splits;}
    void set_number_splits (int s) { _number_splits = s;}

    int  stratified_sampling () const { return _stratified_sampling;}
    void set_stratified_sampling (int s) { _stratified_sampling = s;}

    int  do_paired_t_test () const { return _do_paired_t_test;}
    void set_do_paired_t_test (int s) { _do_paired_t_test = s;}

    int  reorder_equidistant_neighbours () const { return _reorder_equidistant_neighbours;}
    void set_reorder_equidistant_neighbours (int s) { _reorder_equidistant_neighbours = s;}
};

Cross_Validation_Conditions::Cross_Validation_Conditions ()
{
  _percent_training_set = 25;
  _number_splits = 0;
  _stratified_sampling = 0;
  _do_paired_t_test = 0;
  _reorder_equidistant_neighbours = 0;

  return;
}

/*
  This class holds both all the information needed to define the computation as well as
  the resulting measure of the quality of the estimate - Bsquared or Error rate

  I haven't made a separate class for classification and regression. Probably should do
  that, but haven't just in the interests of keeping the number of templates down
*/

class BSquared
{
  private:
    Accumulator<double> _acc;

    resizable_array<double> _bsquared_values;   // one for each split

    int _number_neighbours;

    double _radius;

//  The conditions of an estimate type include the weight function used
   
    NN_Weight_Function _weight_function;

//  when doing classification problems, we may have a probability threshold
//  for the special class

    category_type_t _special_category;
    probability_type_t _special_category_threshold;

//  I tried to use bsquared_comparitor_careful to sort an array of Bsquared values,
//  but it fails. this is because some comparisons are done on the basis of average
//  bsquared values, some comparisons on the basis of complexity. Therefore do all
//  pair-wise comparisons and keep track of the number of times we are ranked best

    int _number_times_better;

//  A human readable description of what kind of model we are

    IWString _type;

//  It can be informative to assign each condition a number based on the
//  initial Bsquared value. 

    int _sequence;

  public:
    BSquared ();

    int build_from_text_string (const IWString &);

    void remove_all_previous_prediction_data (int s);

    void set_sequence (int s) { _sequence = s;}
    int  sequence () const { return _sequence;}

    category_type_t special_category () const { return _special_category;}
    void set_special_category (category_type_t c) { _special_category = c;}

    probability_type_t special_category_threshold () const { return _special_category_threshold;}
    void set_special_category_threshold (probability_type_t p) { _special_category_threshold = p;}

    int  number_neighbours () const { return _number_neighbours;}
    void set_number_neighbours (int n) { _number_neighbours = n;}

    double radius () const { return _radius;}
    void  set_radius (double n) { _radius = n;}

    const IWString & estimate_type () const { return _type;}
    int   compute_estimate_type ();

    const NN_Weight_Function & weight_function () const { return _weight_function;}

    double bsquared_value (int i) const { return _bsquared_values[i];}

//  The companion programme nn_predictions, needs to receive its estimate type in a standard form

    int write_nn_prediction_options (ostream & os) const;

    void extra (double);

    int  report (ostream & os, int = 0) const;

    double average_bsquared () const { return _acc.average ();}   // no checking to make sure we have at least 2 values

    const resizable_array<double> & bsquared_values () const { return _bsquared_values;}

    int number_times_better () const { return _number_times_better;}
    void you_are_better () { _number_times_better++;}
};

BSquared::BSquared ()
{
  _number_neighbours = -1;
  _radius = -1.0;

  _sequence = -1;

  _number_times_better = 0;

  _special_category_threshold = static_cast<probability_type_t> (0.0);
  _special_category = -1;

  return;
}

void
BSquared::extra (double b)
{
  _acc.extra (b);

  _bsquared_values.add (b);

  return;
}

int
BSquared::report (ostream & os, int insert_comment_char) const
{
  if (insert_comment_char)
    os << "# ";

  if (_sequence >= 0)
    os << _sequence << ' ';

  os << bsquared_or_error_rate << " for ";
  if (_type.length ())
    os << _type;
  else
  {
    if (_radius >= 0.0)
      os << "radius " << _radius << ' ';

    if (_number_neighbours >= 0)
      os << _number_neighbours << " neighbours";
  }

  if (_number_times_better)
    os << ". Better than " <<_number_times_better << " others";

  os << '\n';

  if (_acc.n ())
  {
    if (insert_comment_char)
      os << "# ";

    os << _acc.n () << " " << bsquared_or_error_rate << " values between " << _acc.minval () << " and " << _acc.maxval ();
    if (_acc.n () > 1)
    {
      os << " ave " << _acc.average ();
      os << " std err " <<  static_cast<float> (sqrt (_acc.variance ()) / sqrt (static_cast<double> (_acc.n ())));
    }
    os << '\n';

    assert (_bsquared_values.size() == _acc.n());
  }

  return os.good ();
}

/*
  Our conditions have been specified, we now need to determine our type string.
*/

int
BSquared::compute_estimate_type ()
{
  if (_type.length ())
    _type.resize_keep_storage (0);

  if (_radius > 0.0)
  {
    _type += "radius ";
    _type.append_number (_radius, 4);
  }

  if (_number_neighbours > 0)
  {
    if (_type.length ())
      _type += ", ";
    _type << _number_neighbours << " neighbours";
  }

  _type << ' ' << _weight_function.weight_function_type () << " weight";

  if (_special_category >= 0)
  {
    _type << ", ebal " << _special_category << " prob " << _special_category_threshold;
  }

  return 1;
}

/*
  The description of a model must be parsable by parse_model_specifications. We
  use ';' as the separator between tokens, so make sure that character never
  gets used within a specification
*/

int
BSquared::write_nn_prediction_options (ostream & os) const
{
  if (_radius >= 0.0)
  {
    os << _radius;
    if (_number_neighbours > 0)
      os << '-';
  }

  if (_number_neighbours > 0)
    os << _number_neighbours << 'N';

  os << MODEL_COMPONENT_SEPARATOR << "W:" << _weight_function.weight_function_type ();

  if (_special_category >= 0)
    os << MODEL_COMPONENT_SEPARATOR << "S:" << _special_category << '=' << _special_category_threshold;

  return os.good ();
}

void
BSquared::remove_all_previous_prediction_data (int s)
{
  _bsquared_values.resize_keep_storage (0);
  if (s > 0)
    _bsquared_values.resize (s);

  _acc.reset ();

  return;
}

int
BSquared::build_from_text_string (const IWString & n)
{
  const_IWSubstring weight_function;

//cerr << "Before building radius " << _radius << " neighbours " << _number_neighbours << " string '" << n << "'\n";

  if (! parse_model_specifications (n, weight_function, _special_category, _special_category_threshold, _radius, _number_neighbours))
  {
    cerr << "BSquared::build_from_text_string: cannot parse specification '" << n << "'\n";
    return 0;
  }

//cerr << "After building radius " << _radius << " neighbours " << _number_neighbours << endl;

  if (_special_category >= 0 && ! classification_model)
  {
    cerr << "BSquared::build_from_text_string: cannot have special class in non-classification model '" << n << "'\n";
    return 0;
  }

// Exceptionally lazy bit of programming, we never delete these weight functions

  if (0 == weight_function.length ())
  {
    _weight_function.initialise_default_weight_function ();
  }
  else
  {
    if (! _weight_function.construct_from_command_line_token (weight_function, 'W'))
    {
      cerr << "BSquared::build_from_text_string: cannot parse weight function '" << weight_function << "'\n";
      return 0;
    }
  }

  compute_estimate_type ();

  return 1;
}

static int
average_bsquared_comparitor (const void * b1, const void * b2)
{
  double a1 = ((const BSquared *) b1)->average_bsquared ();
  double a2 = ((const BSquared *) b2)->average_bsquared ();

  if (a1 > a2)
    return -1;
  if (a1 < a2)
    return 1;

  return 0;
}

/*
  Sometimes we want to print out the details of the comparisons
*/

static int print_comparison_details = 0;

/*
  We have found that at our given level of confidence, we cannot resolve two distributions.
  See if they can be resolved based on model complexity
*/

static int
resolve_by_complexity (const BSquared & b1, const BSquared & b2,
                       double variance_between_bsquared)
{
  careful_comparisons_resolved_by_complexity++;

// According to the t-test, the two distributions are not different, so 
// resolve them by their characteristics

  double r1 = b1.radius ();
  double r2 = b2.radius ();

  if (r1 > r2)
    return -1;

  if (r1 < r2)
    return 1;

// Radii are the same. Can we resolve by neighbours

  int n1 = b1.number_neighbours ();
  int n2 = b2.number_neighbours ();

// If the variance is zero, the radii are the same, but the number of
// required neighbours different, then we go with the one with the
// fewest number of required neighbours, because probably the larger
// number of neighbours was never found

  if (0.0 == variance_between_bsquared && r1 > 0.0)
  {
    if (n1 > n2)
      return 1;
    else if (n1 < n2)
      return -1;
  }

// Normally, we favour the one with the most neighbours

  if (n1 > n2)
    return -1;

  if (n1 < n2)
    return 1;

  return average_bsquared_comparitor (&b1, &b2);
}

/*
  The person did just one test/train split!!
  We can't test a distribution of Bsquared values, just two numbers
*/

static int
check_just_two_bsquared_values (const BSquared & bs1,
                                const BSquared & bs2)
{
  double b1 = bs1.bsquared_values ()[0];
  double b2 = bs2.bsquared_values ()[0];

  if (fabs (b1 - b2) < 0.002)
    return resolve_by_complexity (bs1, bs2, 0.0);

  if (min_meaningful_bsquared_difference > 0.0 && fabs (b1 - b2) < min_meaningful_bsquared_difference)
    return resolve_by_complexity (bs1, bs2, 0.0);

  if (b1 > b2)
    return -1;

  if (b1 < b2)
    return 1;

  return resolve_by_complexity (bs1, bs2, 0.0);
}

/*
  We have a more careful way of comparing BSquared objects. Each one has
  a range of B2 values encountered. We need to know whether or not they are
  significantly different. Use a paired t-test
*/

static int
bsquared_comparitor_careful (const void * pb1, const void * pb2)
{
//cerr << "Comparing " << pb1 << " with " << pb2 << endl;

  const BSquared & b1 = *((const BSquared *) pb1);
  const BSquared & b2 = *((const BSquared *) pb2);

//cerr << "Comparing '" << b1.estimate_type () << "' " << b1.average_bsquared () << " with '" << b2.estimate_type () << "' " << b2.average_bsquared () << "\n";

  careful_comparisons_done++;

  const resizable_array<double> & bs1 = b1.bsquared_values ();
  const resizable_array<double> & bs2 = b2.bsquared_values ();

  int n = bs1.number_elements ();

  if (n != bs2.number_elements ())
  {
    cerr << "Yipes, comparing " << n << " and " << bs2.number_elements () << " items\n";
    cerr << b1.radius () << " and " << b2.radius () << endl;
    cerr << b1.estimate_type () << " and " << b2.estimate_type () << endl;
    cerr << hex << pb1 << " and " << hex << pb2 << dec << endl;
  }
  assert (n == bs2.number_elements ());

  if (1 == n)
    return check_just_two_bsquared_values (b1, b2);

  Accumulator<double> acc;

  for (int i = 0; i < n; i++)
  {
    double d = bs1[i] - bs2[i];

    acc.extra (d);
  }

  if (0.0 == acc.variance ())    // distributions are identical
    return resolve_by_complexity (b1, b2, 0.0);

// The distributions are not identical. Are they so different that differences are not meaningful?

  if (min_meaningful_bsquared_difference > 0.0)
  {
    double a1 = b1.average_bsquared ();
    double a2 = b2.average_bsquared ();

    if (fabs (a1 - a2) <= min_meaningful_bsquared_difference)
      return resolve_by_complexity (b1, b2, acc.variance ());
  }

// Use the paired t-test to see if the distributions are different or not

  double t = acc.average () / (sqrt (acc.variance ()) / sqrt (static_cast<double> (n)));

//#define USE_FULL_PVALUE_COMPUTATION
#ifdef USE_FULL_PVALUE_COMPUTATION
  double p = iwpvalue (n, fabs (t));
//cerr << "P value for " << n << " " << fabs (t) << " " << p << endl;
//cerr << "Ave " << acc.average () << " variance " << acc.variance () << " n = " << n << endl;

  if (print_comparison_details)
    cerr << "Compare '" << b1.estimate_type () << "' (" << b1.average_bsquared () << ") and '" << b2.estimate_type () << "' (" << b2.average_bsquared () << ") T = " << t << " P = " << p << endl;

  if (p <= pvalue_threshold)    // they are different at this level
    return average_bsquared_comparitor (pb1, pb2);
#else
  if (fabs (t) <= tvalue_threshold)
    return average_bsquared_comparitor (pb1, pb2);
#endif

// Not different at this significance level

  return resolve_by_complexity (b1, b2, acc.variance ());
}

static void
do_all_pair_wise_comparisons (BSquared * bsquared,
                              int number_estimates)
{
#ifdef TIME_COMPARISONS
  time_t t0 = time (NULL);
  cerr << "Starting comparisons at " << t0 << endl;
#endif

  for (int i = 0; i < number_estimates; i++)
  {
    BSquared & bsi = bsquared[i];

    for (int j = i + 1; j < number_estimates; j++)
    {
      BSquared & bsj = bsquared[j];

      int c = bsquared_comparitor_careful (&bsi, &bsj);
      if (0 == c)
        ;
      else if (c < 0)
        bsi.you_are_better ();
      else
        bsj.you_are_better ();
    }
  }

#ifdef TIME_COMPARISONS
  cerr << "Finished comparisons at " << time (NULL) << ", delta " << (time (NULL) - t0) << endl;
#endif

  return;
}

static int
bsquared_comparitor_number_better (const void * pb1, const void * pb2)
{
  const BSquared & b1 = *((const BSquared *) pb1);
  const BSquared & b2 = *((const BSquared *) pb2);

  int nb1 = b1.number_times_better ();
  int nb2 = b2.number_times_better ();

  if (nb1 > nb2)
    return -1;
  else if (nb1 < nb2)
    return 1;

  double avbs1 = b1.average_bsquared ();
  double avbs2 = b2.average_bsquared ();

  double average_bsquared_difference = fabs (avbs1 - avbs2);

  if (min_meaningful_bsquared_difference > 0.0 && average_bsquared_difference > min_meaningful_bsquared_difference)
    return average_bsquared_comparitor (pb1, pb2);

  return resolve_by_complexity (b1, b2, 0.0);    // just pass 0.0 rather than recomputing stuff
}

/*
  When we read experimental data from a file we need a hash to store that data
*/

template <typename E>
class ID_Experimental_Value_Hash : public IW_STL_Hash_Map<IWString, E>
{
  private:
  public:
};

//typedef ID_Experimental_Value_Hash<experimental_value_t> ID_Activity_Hash;
//typedef ID_Experimental_Value_Hash<category_type_t> ID_Category_Hash;

class ID_Activity_Hash : public ID_Experimental_Value_Hash<experimental_value_t>
{
  private:
  public:
    typedef ID_Experimental_Value_Hash<experimental_value_t>::const_iterator const_iterator;
};

class ID_Category_Hash : public ID_Experimental_Value_Hash<category_type_t>
{
  private:
  public:
    typedef IW_STL_Hash_Map<IWString, category_type_t>::const_iterator const_iterator;
};

/*
  Each NN_Item will have a resizable_array of these objects.
  The near neighbour results will be derived from running the training
  set against itself, so everything that appears as a neighbour, will
  also be a target somewhere. 

  This comes from gfp_nearneighbours with the -o option, so the id will be
  a numeric index into the pool

  As we subdivide the training set, we need to know whether or not a
  given neighbour is in the training or test set. For that reason,
  we have a pointer to our corresponding NN_Item object
*/

template <typename E, typename N>
class ID_Dist_Observed 
{
  protected:

    int _ndx;

    similarity_type_t _distance;

// Weights are computed from the distance.  By default, the computation
// will take a float and convert it to an int in order to index into the
// weight function's array.  Much time saved if we pre-compute the index

    int _weight_function_index;

    E _experimental_result;

//  When doing an estimate, each NN_Item needs to know which of its neighbours
//  are in the training set. Therefore each neighbour has a pointer to its
//  corresponding NN_Item and can get the training set membership from there.

    const N * _nn_item;

  public:
    ID_Dist_Observed ();

    int build (iwstring_data_source &, int &);

    int in_training_set () const { return _nn_item->in_training_set ();}

    int index_in_pool () const { return _ndx;}
    void set_index_in_pool (int i) { _ndx = i;}

    const IWString & id () const { return _nn_item->id ();}

    similarity_type_t distance () const { return _distance;}
    void set_distance (similarity_type_t d) { _distance = d;}

    weight_t weight (const NN_Weight_Function & wfn) const { return wfn.weight (_weight_function_index);}

//  void set_nn_item (const N * n) { _nn_item = n;}
//  Never got the templates to work, so we cast the address to void *
    void set_nn_item (const void * n);

    int set_experimental_value_from_id_token (int token);
    int set_experimental_value_from_token (const const_IWSubstring & zprefix);

    E experimental_value () const { return _experimental_result;}
    void set_experimental_value (E e) { _experimental_result = e;}
};

template <typename E, typename N>
ID_Dist_Observed<E, N>::ID_Dist_Observed ()
{
  _ndx = -1;

  _nn_item = nullptr;

  _weight_function_index = -1;

  _experimental_result = static_cast<E> (0);

  return;
}

/*
  As soon as we know the identity of our parent object in the pool, we fetch our
  observed activity from them. Make sure that the pool objects have their experimental
  values determined first.
*/

template <typename E, typename N>
void
ID_Dist_Observed<E, N>::set_nn_item (const void * n)
{
  assert (NULL != n);

  _nn_item = reinterpret_cast<const N *> (n);

  _experimental_result = _nn_item->experimental_value ();

  return;
}

/*
  Common code for extracting the experimental value from a token in a string

  We can speed things up a lot by maintaining a cache of string representations of
  the experimental data - implement this sometime...
*/

static IW_STL_Hash_Map_float expt_to_float;
static IW_STL_Hash_Map_int expt_to_int;

template <typename B, typename A>
int
get_experimental_value (const B & buffer,
                        int data_column,
                        A & expt)
{
//cerr << "Fetching word " << data_column << " from '" << buffer << "'\n";

  if (data_column < 0)
  {
    int d = buffer.nwords () - data_column;
    if (d < 0)
    {
      cerr << "Cannot extract column " << (-data_column) << " from the end of '" << buffer << "'\n";
      return 0;
    }

    data_column = d;
  }

  IWString sexpt;
  if (! buffer.word (data_column, sexpt))
  {
    cerr << "Cannot extract experimental value, column " << data_column << " in '" << buffer << "'\n";
    return 0;
  }

  if (! sexpt.numeric_value (expt))
  {
    cerr << "Invalid numeric value for experimental activity '" << sexpt << "'\n";
    return 0;
  }

  return 1;
}

template <typename E, typename N>
int
ID_Dist_Observed<E, N>::set_experimental_value_from_id_token (int c)
{
  const IWString & parent_id = _nn_item->id ();

  if (! get_experimental_value (parent_id, c, _experimental_result))
  {
    cerr << "ID_Dist_Observed::set_experimental_value_from_id_token: cannot set activity from id token\n";
    return 0;
  }

  return 1;
}

/*
  Determine the ID and the DISTANCE
*/

template <typename E, typename N>
int
ID_Dist_Observed<E, N>::build (iwstring_data_source & input, int & fatal)
{
  fatal = 0;

  int got_distance = 0;

  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (buffer.starts_with (neighbour_tag))
    {
      if (_ndx >= 0)
      {
        cerr << "ID_Dist::build: duplicate id on line " << input.lines_read () << endl;
        return 0;
      }

      if (! extract_from_tdt_form (buffer, neighbour_tag, _ndx))
      {
        cerr << "ID_Dist::build: invalid identifier record on line " << input.lines_read () << endl;
        return 0;
      }
    }
    else if (buffer.starts_with (distance_tag))
    {
      if (got_distance)
      {
        cerr << "ID_Dist::build: duplicate distance on line " << input.lines_read () << endl;
        return 0;
      }

      if (! extract_from_tdt_form (buffer, distance_tag, _distance) || _distance < 0.0)
      {
        cerr << "ID_Dist_Observed::build: invalid distance form '" << buffer << "', line " << input.lines_read () << endl;
        return 0;
      }

      got_distance = 1;

      if (! compute_weight_function_index (_distance, _weight_function_index))
        return 0;
    }
    else if ('|' == buffer)
    {
      if (_ndx >= 0 || got_distance)
      {
        cerr << "ID_Dist::build: incomplete specification\n";
        if (_ndx >= 0)
          cerr << "ndx " << _ndx << endl;
        if (got_distance)
          cerr << "dist " << _distance << endl;
        fatal = 1;

        return 0;
      }

      return 0;     // got to the end of the TDT, no information
    }
    else
      continue;

    if (_ndx >= 0 && got_distance)
      return 1;
  }

  cerr << "ID_Dist_Observed::build: unexpected EOF\n";
  return 0;
}

template <typename E>
int
extract_experimental_value_from_token (const IWString & buffer,
                                       const const_IWSubstring & zprefix,
                                       E & zresult)
{
  int i = 0;
  const_IWSubstring token;

  while (buffer.nextword (token, i))
  {
    if (! token.starts_with (zprefix))
      continue;

    token.remove_leading_chars (zprefix.length ());

    if (! token.numeric_value (zresult))
    {
      cerr << "Invalid '" << zprefix << "' activity '" << buffer << "'\n";
      return 0;
    }

    return 1;
  }

  cerr << "extract_experimental_value_from_token: no identifier tokens start with '" << zprefix << "'\n";
  return 0;
}

template <typename E, typename N>
int
ID_Dist_Observed<E, N>::set_experimental_value_from_token (const const_IWSubstring & zprefix)
{
  const IWString parent_id = _nn_item->id ();

  return extract_experimental_value_from_token (parent_id, zprefix, _experimental_result);
}

template int get_experimental_value (const const_IWSubstring &, int, experimental_value_t &);
template int get_experimental_value (const IWString &, int, experimental_value_t &);

/*
  It is optional as to whether or not the target molecules have an experimental activity.
*/

template <typename E, typename N, typename H>
class NN_Item : public NN_Item_Base<N>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using NN_Item_Base<N>::_id;
    using NN_Item_Base<N>::_number_elements;
    using NN_Item_Base<N>::_things;
//  using resizable_array_base<N *>::swap_elements(int, int);
#endif

  protected:
    E _experimental_result;     // category or continuous response

//  for each prediction, we need a prediced value

    int _number_predictions;

    Set_or_Unset<E> * _predicted;

//  When doing cross validation, we need to assign items to test and training sets,
//  even though this whole set of items is a training set

    int _in_training_set;

//  If we are doing the shuffling of equidistant neighbours, we can save time by
//  recording whether or not we have equidistant neighbours

    int _has_equidistant_neighbours;

    int _times_in_training_set;

//  private functions

    void _determine_whether_I_have_equidistant_neighbours ();

  public:
    NN_Item ();
    ~NN_Item ();

    int debug_print (ostream & os, int) const;

    int times_in_training_set () const { return _times_in_training_set;}

    E experimental_value () const { return _experimental_result;}

    int predicted_value (int i, E &) const;

    int write_predicted_values (ostream & output) const;

    int size_predicted_array (int p);

    void unset_predicted_values ();

    int shuffle_equidistant_neighbours ();

/*
    I never could get the templates for this to work out, so we cheat.

    int initialise_nn_pointers (const IW_STL_Hash_Map<IWString, NN_Item_Ptr> &);
*/

    int initialise_nn_pointers (void * *, int);

    similarity_type_t max_neighbour_distance () const;

    int in_training_set () const { return _in_training_set;}
    void set_training_set_membership (int s);

    int get_experimental_values_from_identifier_token (int col);
    int get_experimental_value_from_token (const const_IWSubstring & zprefix);

    int assign_experimental_result_from_hash (const H &);
};

template <typename E, typename N, typename H>
NN_Item<E, N, H>::NN_Item ()
{
  _predicted = nullptr;

  _in_training_set = 0;

  _has_equidistant_neighbours = -1;

  _times_in_training_set = 0;

  return;
}

template <typename E, typename N, typename H>
NN_Item<E, N, H>::~NN_Item ()
{
  if (NULL != _predicted)
    delete [] _predicted;

  return;
}

template <typename E, typename N, typename H>
int
NN_Item<E, N, H>::debug_print (ostream & os, int neighbours_to_print) const
{
  os << "NN_Item::experimental value " << _experimental_result;
  if (_in_training_set)
    os << " In training set";
  else
    os << " Not in training set";

  os << ", has " << _number_elements << " near neighbours\n";

  int nprint;
  if (neighbours_to_print > _number_elements)
    nprint = _number_elements;
  else
    nprint = neighbours_to_print;

  for (int i = 0; i < nprint; i++)
  {
    const N & sidda = *(_things[i]);

    os << ' ' << i << " '" << sidda.id () << "' dist " << sidda.distance () << " expt " << sidda.experimental_value () << '\n';
  }

  return os.good ();
}

template <typename E, typename N, typename H>
void
NN_Item<E, N, H>::set_training_set_membership (int s)
{
  _in_training_set = s;

  if (s)
    _times_in_training_set++;

  return;
}

template <typename E, typename N, typename H>
int
NN_Item<E, N, H>::predicted_value (int i, E & e) const
{
  assert (NULL != _predicted);

  return _predicted[i].value (e);
}

template <typename E, typename N, typename H>
int
NN_Item<E, N, H>::size_predicted_array (int p)
{
  assert (p > 0);
  assert (NULL == _predicted);

  _number_predictions = p;
  _predicted = new Set_or_Unset<E>[p];

  return NULL != _predicted;
}

template <typename E, typename N, typename H>
void
NN_Item<E, N, H>::unset_predicted_values ()
{
  assert (NULL != _predicted);

  for (int i = 0; i < _number_predictions; i++)
  {
    _predicted[i].unset ();
  }

  return;
}

template <typename E, typename N, typename H>
similarity_type_t
NN_Item<E, N, H>::max_neighbour_distance () const
{
  similarity_type_t rc = static_cast<similarity_type_t> (0.0);

  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->distance () > rc)
      rc = _things[i]->distance ();
  }

  return rc;
}

template <typename E, typename N, typename H>
int
NN_Item<E, N, H>::shuffle_equidistant_neighbours ()
{
  if (0 == _has_equidistant_neighbours)
    return 0;

  if (-1 == _has_equidistant_neighbours)
  {
    _determine_whether_I_have_equidistant_neighbours ();
    if (0 == _has_equidistant_neighbours)
      return 0;
  }

  similarity_type_t prevd = _things[0]->distance ();

  for (int i = 1; i < _number_elements; i++)
  {
    similarity_type_t d = _things[i]->distance ();

    if (d == prevd)
      this->swap_elements (i, i - 1);
    else
      prevd = d;
  }

  return 1;
}

template <typename E, typename N, typename H>
void
NN_Item<E, N, H>::_determine_whether_I_have_equidistant_neighbours ()
{
  similarity_type_t prevd = _things[0]->distance ();

  for (int i = 1; i < _number_elements; i++)
  {
    similarity_type_t d = _things[i]->distance ();

    if (d == prevd)
    {
      _has_equidistant_neighbours = 1;

      return;
    }

    prevd = d;
  }

  _has_equidistant_neighbours = 0;

  return;
}

/*
  The experimental values are tokens in the name fields.
*/

template <typename E, typename N, typename H>
int
NN_Item<E, N, H>::get_experimental_values_from_identifier_token (int col)
{
  if (! get_experimental_value (_id, col, _experimental_result))
  {
    cerr << "Training set member '" << _id << "' has no experimental value.\n";
    return 0;
  }

  return 1;
}

class ID_Dist_Activity;

class NN_Item_Activity : public NN_Item<experimental_value_t, ID_Dist_Activity, ID_Activity_Hash>
{
  private:

//  We keep track of all the estimates made

    Accumulator<experimental_value_t> _acc_pred;

//  private functions

    int _compute_predicted_values (int neighbours_to_use, int ndx, const NN_Weight_Function &);
    int _compute_predicted_values (similarity_type_t threshold_to_use, int min_neighbours_required, int ndx, const NN_Weight_Function & wfn);

  public:

    int compute_predicted_values (const BSquared & bs, int which_prediction, const NN_Weight_Function & wfn);

    int write_range_of_predicted_values (double, ostream & output) const;

//  We can identify molecules that are poorly predicted

    int average_absolute_prediction_error (experimental_value_t &) const;
};

typedef IW_STL_Hash_Map<IWString, NN_Item_Activity> ID_NN_Item_Activity_Hash;

class ID_Dist_Category;

class NN_Item_Category : public NN_Item<category_type_t, ID_Dist_Category, ID_Category_Hash>
{
  private:

//  For classification problems, we need to keep track of the probability of each item
//  being placed in each category. Kind of awkward, we have an array of arrays

    probability_type_t ** _probability;

//  We can also keep track of how often an item is mis-classified

    int _misclassifications;

//  private functions

    int _determine_category_with_highest_weight (int);
    int _determine_category_with_highest_weight (category_type_t, probability_type_t, int);

    int _compute_predicted_values (int neighbours_to_use, int ndx, const NN_Weight_Function & wfn);
    int _compute_predicted_values (similarity_type_t threshold_to_use, int min_neighbours_required, int ndx, const NN_Weight_Function & wfn);
    int _compute_predicted_values (int neighbours_to_use, category_type_t special_category, probability_type_t special_category_threshold, int ndx, const NN_Weight_Function & wfn);
    int _compute_predicted_values (similarity_type_t threshold_to_use, int min_neighbours_required, category_type_t special_category, probability_type_t special_category_threshold, int ndx, const NN_Weight_Function & wfn);

  public:
    NN_Item_Category ();
    ~NN_Item_Category ();

    int allocate_probability_arrays ();

    category_type_t category () const { return _experimental_result;}
    void set_category (category_type_t s) { _experimental_result = s;}

    int compute_predicted_values (const BSquared & bs, int which_prediction, const NN_Weight_Function & wfn);

    void another_misclassification () { _misclassifications++;}
    int  misclassifications () const { return _misclassifications;}
};

class ID_Dist_Activity : public ID_Dist_Observed<similarity_type_t, NN_Item_Activity>
{
  private:
    
  public:
    ID_Dist_Activity ();

    void increment_prediction_and_weight (experimental_value_t & pred,
                                          weight_t & sumw,
                                          const NN_Weight_Function & wfn) const;
};

class ID_Dist_Category: public ID_Dist_Observed<category_type_t, NN_Item_Category>
{
  private:

  public:
    ID_Dist_Category ();

    category_type_t category () const;
};

ID_Dist_Activity::ID_Dist_Activity ()
{
  _experimental_result = static_cast<experimental_value_t> (-9.9);

  return;
}

/*
  Someone is doing an activity model and needs the weight and weight * result
*/

void
ID_Dist_Activity::increment_prediction_and_weight (experimental_value_t & pred,
                                                   weight_t & sumw,
                                                   const NN_Weight_Function & wfn) const
{
  weight_t w = wfn.weight (_weight_function_index);

  sumw += w;
  pred += w * _experimental_result;

  return;
}

ID_Dist_Category::ID_Dist_Category ()
{
  return;
}

template <typename E, typename N, typename H>
int
NN_Item<E, N, H>::get_experimental_value_from_token (const const_IWSubstring & zprefix)
{
  if (! extract_experimental_value_from_token (_id, zprefix, _experimental_result))
  {
    cerr << "NN_Item::get_experimental_value_from_token: no experimental data for '" << _id << "'\n";
    return 0;
  }

  return 1;
}

/*
  Would have been more elegant with

  H::const_iterator

  but that would not compile - not sure why
*/

template <typename E, typename N, typename H>
int
NN_Item<E, N, H>::assign_experimental_result_from_hash (const H & idah)
{
  if (idah.contains (_id))
  {
    _experimental_result = (*(idah.find (_id))).second;
    return 1;
  }

  if (_id.contains (' '))
  {
    IWString tmp (_id);
    tmp.truncate_at_first (' ');
    if (idah.contains (tmp))
    {
      _experimental_result = (*(idah.find (tmp))).second;
      return 1;
    }
  }

  if (strip_leading_zeros_from_identifiers)
  {
    IWString tmp (_id);
    tmp.remove_leading_chars ('0');
    tmp.truncate_at_first (' ');

    if (idah.contains (tmp))
    {
      _experimental_result = (*(idah.find (tmp))).second;
      return 1;
    }
  }

  cerr << "NN_Item::assign_experimental_result_from_hash: no data for '" << _id << "'\n";
  return 0;
}

/*
  Write the results of the last set of predictions
*/

template <typename E, typename N, typename H>
int
NN_Item<E, N, H>::write_predicted_values (ostream & output) const
{
  assert (NULL != _predicted);

  for (int i = 0; i < _number_predictions; i++)
  {
    E e;

    if (_predicted[i].value (e))
    {
      output << ' ' << e;
    }
    else
    {
      output << ' ' << missing_value;
    }
  }

  return output.good ();
}

int
NN_Item_Activity::write_range_of_predicted_values (double average_error, ostream & os) const
{
  os << _acc_pred.n () << ' ' << _experimental_result << ' ';

  if (0 == _acc_pred.n ())
  {
    os << missing_value << ' ' << missing_value << ' ' << missing_value;
    return os.good ();
  }

  os.precision (4);

  if (_acc_pred.n () > 1)
    os << _acc_pred.average ();
  else
    os << missing_value;

  os << ' ' << _acc_pred.minval () << "  " << _acc_pred.maxval () << "  " << average_error << "  " << _things[0]->distance ();

  return os.good ();
}

NN_Item_Category::NN_Item_Category ()
{
  _misclassifications = 0;

  _probability = nullptr;
}

NN_Item_Category::~NN_Item_Category ()
{
  if (NULL != _probability)
  {
    for (int i = 0; i < _number_predictions; i++)
    {
      if (NULL != _probability[i])
        delete _probability[i];
    }

    delete _probability;

    _probability = nullptr;
  }
  
  return;
}

int
NN_Item_Category::allocate_probability_arrays ()
{
  assert (NULL == _probability);
  assert (highest_category_number > 0);

  _probability = new probability_type_t *[_number_predictions];

  for (int i = 0; i < _number_predictions; i++)
  {
    _probability[i] = new probability_type_t[highest_category_number + 1];
  }

  return 1;
}

/*
  Each of the classification predicters need to figure out which category
  got the highest total
*/

int
NN_Item_Category::_determine_category_with_highest_weight (int ndx)
{
  assert (NULL != _predicted);

  weight_t highest_weight = category_sumw[0];
  int zprediction = 0;

  for (int c = 1; c <= highest_category_number; c++)
  {
    weight_t s = category_sumw[c];

    if (s > highest_weight)
    {
      highest_weight = s;
      zprediction = c;
    }
  }

  _predicted[ndx] = zprediction;

  return 1;
}

/*
  We are doing error balancing and one class has been designed as special, with a probability threshold
*/

int
NN_Item_Category::_determine_category_with_highest_weight (category_type_t special_category,
                                                           probability_type_t special_category_threshold,
                                                           int ndx)
{
  assert (NULL != _predicted);
  assert (special_category >= 0 && special_category <= highest_category_number);

  weight_t sumw = sum_vector (category_sumw, highest_category_number + 1);

  probability_type_t p = category_sumw[special_category] / sumw;

#ifdef ECHO_DECISIONS
  for (int i = 0; i <= highest_category_number; i++)
  {
    cerr << "class " << i << " weight " << category_sumw[i] << '\n';
  }
  cerr << "Class " << special_category << " p = " << p << " threshold " << special_category_threshold << '\n';
#endif

  if (p >= special_category_threshold)
  {
    _predicted[ndx] = special_category;
    return 1;
  }

// Did not pass the threshold, just do things the normal way

  return _determine_category_with_highest_weight (ndx);
}

int
NN_Item_Category::compute_predicted_values (const BSquared & bs,
                                   int which_prediction,
                                   const NN_Weight_Function & wfn)
{
  set_vector (category_sumw, highest_category_number + 1, static_cast<weight_t> (0.0));

  int neighbours = bs.number_neighbours ();
  double radius = bs.radius ();
  category_type_t special_class = bs.special_category ();
  
// do the easy case first

  if (special_class < 0)
  {
    if (neighbours > 0 && radius < 0.0)
      return _compute_predicted_values (neighbours, which_prediction, wfn);
    else
      return _compute_predicted_values (radius, neighbours, which_prediction, wfn);
  }

  if (neighbours > 0 && radius < 0.0)
    return _compute_predicted_values (neighbours, special_class, bs.special_category_threshold (), which_prediction, wfn);
  else
    return _compute_predicted_values (radius, neighbours, special_class, bs.special_category_threshold (), which_prediction, wfn);
}

int
NN_Item_Category::_compute_predicted_values (int neighbours_to_use,
                                             int which_prediction,
                                             const NN_Weight_Function & wfn)
{
//cerr << "neighbours_to_use " << neighbours_to_use << " _number_elements " << _number_elements << endl;
  if (neighbours_to_use > _number_elements)
    neighbours_to_use = _number_elements;

  int items_used_in_prediction = 0;

  for (int i = 0; i < neighbours_to_use; i++)
  {
    const ID_Dist_Category & sida = *(_things[i]);

    if (! sida.in_training_set ())     // we want to use training set molecules to predict test set molecules
      continue;

    int c = sida.category ();

    category_sumw[c] += sida.weight (wfn);

    items_used_in_prediction++;
    if (items_used_in_prediction >= neighbours_to_use)
      break;
  }

  return _determine_category_with_highest_weight (which_prediction);
}

int
NN_Item_Category::_compute_predicted_values (similarity_type_t threshold_to_use,
                                   int min_neighbours_required,
                                   int which_prediction,
                                   const NN_Weight_Function & wfn)
{
  int items_used_in_prediction = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    const ID_Dist_Category & sida = *(_things[i]);

    if (! sida.in_training_set ())    // we want to use training set molecules to predict test set molecules
      continue;

    similarity_type_t d = sida.distance ();

    if (d <= threshold_to_use)
      ;
    else if (min_neighbours_required > 0 && items_used_in_prediction < min_neighbours_required)
      ;
    else
      break;

    int c = sida.category ();

    category_sumw[c] += sida.weight (wfn);

    items_used_in_prediction++;
  }

//#define DEBUG_ITEMS_USED
#ifdef DEBUG_ITEMS_USED
  if (0 == items_used_in_prediction)
    cerr << "No neighbours in training set, ndx = " << which_prediction << "\n";
  else
    cerr << "1 Predicted with " << items_used_in_prediction << '\n';
#endif

  return _determine_category_with_highest_weight (which_prediction);
}

int
NN_Item_Category::_compute_predicted_values (int neighbours_to_use,
                                             category_type_t special_category,
                                             probability_type_t special_category_threshold,
                                             int ndx,
                                             const NN_Weight_Function & wfn)
{
  if (neighbours_to_use > _number_elements)
    neighbours_to_use = _number_elements;

  int items_used_in_prediction = 0;

  for (int i = 0; i < neighbours_to_use; i++)
  {
    const ID_Dist_Category & sida = *(_things[i]);

    if (! sida.in_training_set ())     // we want to use training set molecules to predict test set molecules
      continue;

    int c = sida.category ();

    category_sumw[c] += sida.weight (wfn);

    items_used_in_prediction++;
    if (items_used_in_prediction >= neighbours_to_use)
      break;
  }

  return _determine_category_with_highest_weight (special_category, special_category_threshold, ndx);
}

int
NN_Item_Category::_compute_predicted_values (similarity_type_t threshold_to_use,
                                             int min_neighbours_required,
                                             category_type_t special_category,
                                             probability_type_t special_category_threshold,
                                             int ndx,
                                             const NN_Weight_Function & wfn)
{
  int items_used_in_prediction = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    const ID_Dist_Category & sida = *(_things[i]);

    if (! sida.in_training_set ())    // we want to use training set molecules to predict test set molecules
      continue;

    similarity_type_t d = sida.distance ();

    if (d > threshold_to_use)
    {
      if (min_neighbours_required > 0 && items_used_in_prediction < min_neighbours_required)    // need to keep going
        ;
      else
        break;
    }

    int c = sida.category ();

    category_sumw[c] += sida.weight (wfn);

    items_used_in_prediction++;
  }

#ifdef DEBUG_ITEMS_USED
  if (0 == items_used_in_prediction)
    cerr << "No neighbours in training set, ndx = " << ndx << "\n";
  else
    cerr << "2 Predicted with " << items_used_in_prediction << '\n';
#endif

  return _determine_category_with_highest_weight (special_category, special_category_threshold, ndx);
}

int
NN_Item_Activity::compute_predicted_values (const BSquared & bs,
                                            int which_prediction,
                                            const NN_Weight_Function & wfn)
{
  int neighbours = bs.number_neighbours ();
  double radius = bs.radius ();

  if (neighbours > 0 && radius < 0.0)
    return _compute_predicted_values (neighbours, which_prediction, wfn);
  else
    return _compute_predicted_values (radius, neighbours, which_prediction, wfn);
}

int
NN_Item_Activity::_compute_predicted_values (int neighbours_to_use,
                                             int ndx,
                                             const NN_Weight_Function & wfn)
{
  assert (NULL != _predicted);

  if (neighbours_to_use > _number_elements)
    neighbours_to_use = _number_elements;

  int items_used_in_prediction = 0;

  weight_t sumw = static_cast<weight_t> (0.0);
  experimental_value_t pred = static_cast<experimental_value_t> (0.0);

  for (int i = 0; i < _number_elements; i++)
  {
    const ID_Dist_Activity & sida = *(_things[i]);

    if (! sida.in_training_set ())     // we want to use training set molecules to predict test set molecules
      continue;

    sida.increment_prediction_and_weight (pred, sumw, wfn);

    items_used_in_prediction++;

    if (items_used_in_prediction >= neighbours_to_use)
      break;
  }

#ifdef DEBUG_ITEMS_USED
  if (0 == items_used_in_prediction)
    cerr << "No neighbours in training set, ndx = " << ndx << '\n';
  else
    cerr << "3 Predicted with " << items_used_in_prediction << '\n';
#endif

  if (0 == items_used_in_prediction)
    return 1;

  if (static_cast<weight_t> (0.0) != sumw)
  {
    experimental_value_t zresult = pred / sumw;

    _predicted[ndx] = zresult;
    _acc_pred.extra (zresult);
  }

  return 1;
}

int
NN_Item_Activity::_compute_predicted_values (similarity_type_t threshold_to_use,
                                            int min_neighbours_required,
                                            int ndx,
                                            const NN_Weight_Function & wfn)
{
  assert (NULL != _predicted);

  int items_used_in_prediction = 0;

  experimental_value_t pred = static_cast<experimental_value_t> (0.0);
  weight_t sumw = static_cast<weight_t> (0.0);
  for (int i = 0; i < _number_elements; i++)
  {
    const ID_Dist_Activity & sida = *(_things[i]);

    if (! sida.in_training_set ())    // we want to use training set molecules to predict test set molecules
      continue;

    similarity_type_t d = sida.distance ();

    if (d > threshold_to_use)
    {
      if (min_neighbours_required > 0 && items_used_in_prediction < min_neighbours_required)    // need to keep going
        ;
      else
        break;
    }

    sida.increment_prediction_and_weight (pred, sumw, wfn);

    items_used_in_prediction++;
  }

#ifdef DEBUG_ITEMS_USED
  if (0 == items_used_in_prediction)
    cerr << "No neighbours in training set, ndx = " << ndx << "\n";
  else
    cerr << "4 Predicted with " << items_used_in_prediction << '\n';
#endif

  if (0 == items_used_in_prediction)
    return 1;

  if (sumw > static_cast<weight_t> (0.0))
  {
    experimental_value_t zresult = pred / sumw;

    _predicted[ndx] = zresult;
    _acc_pred.extra (zresult);
  }

  return 1;
}

/*
  Since I could not get the templates to work as I wanted, we pass a hash of
  IWString and pointers to int's. They are actually pointers to NN_Item objects
*/

template <typename E, typename N, typename H>
int
NN_Item<E, N, H>::initialise_nn_pointers (void ** xref,
                                          int items_in_pool)
{
  int rc = 1;
  for (int i = 0; i < _number_elements; i++)
  {
    N * sida = _things[i];

    int ndx = sida->index_in_pool ();

    if (ndx < 0 || ndx >= items_in_pool || NULL == xref[ndx])
    {
      cerr << "initialise_nn_pointers: neighbour " << i << " ndx " << ndx << " invalid (max " << items_in_pool << ") or xref is NULL\n";
      return 0;
    }

    sida->set_nn_item (xref[ndx]);
  }

  return rc;
}

int
NN_Item_Activity::average_absolute_prediction_error (experimental_value_t & rc) const
{
  if (0 == _acc_pred.n ())
    return 0;

  rc = static_cast<experimental_value_t> (_acc_pred.average ()) - _experimental_result;

  if (rc < static_cast<experimental_value_t> (0.0))
    rc = -rc;

  return 1;
}


category_type_t
ID_Dist_Category::category () const
{
  return _nn_item->category ();
}

template <typename N, typename H, typename E>
class NN_Results : public NN_Results_Base<N>
{
  protected:

#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using NN_Results_Base<N>::_number_items;
    using NN_Results_Base<N>::_results;
#endif

//  If we are reading test/train splits from an external file, we want a hash that holds the relationship
//  between text id and index in the _results array

    IW_STL_Hash_Map_int _id_to_index;

//  protected functions

    int _read_experimental_value (const const_IWSubstring & buffer, H & activity);
    int _read_experimental_values (iwstring_data_source & input, H & activity);

    int _assign_training_test_set_membership_external (const Cross_Validation_Conditions & cvc, int i);
  private:

//  private functions

    int _compute_predicted_values (int, int);
    int _compute_predicted_values (float, int, int);

    int _initialise_id_index_structures ();

  public:
    NN_Results ();
    ~NN_Results ();

    int debug_print (int, ostream &) const;

    int determine_experimental_values_file (iwstring_data_source & input);
    int determine_experimental_values_file (const const_IWSubstring & experimental_results_file, int skip_header_line_in_experimental_values_file);
    int determine_experimental_values_hash (H & activity);
    int determine_experimental_values_column (int c);
    int determine_experimental_values_token (const const_IWSubstring &);

    int any_items_have_no_neighbours () const;

    void set_training_set_membership (int);

    int initialise_nn_pointers ();

    int shuffle_equidistant_neighbours ();

//#define IW_DASH_S
#ifdef IW_DASH_S
    int size_the_problem (int) { return 1;}
#endif

    void unset_predicted_values ();

    int compute_leave_one_out_values (const BSquared * bsquared, int number_predictions);
    int compute_leave_one_out_values (const BSquared & bsquared, int which_prediction);

    int compute_predicted_values (const BSquared &, int);
    int compute_predicted_values (const BSquared *, int);

    int write_predicted_values (ostream &, int) const;
};

template <typename N, typename H, typename E>
NN_Results<N, H, E>::NN_Results ()
{
  return;
}

template <typename N, typename H, typename E>
NN_Results<N, H, E>::~NN_Results ()
{
  return;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::debug_print (int depth, ostream & os) const
{
  os << "NN_Results::results for " << _number_items << " items\n";

  for (int i = 0; i < _number_items; i++)
  {
    const N & ni = _results[i];

    ni.debug_print (os, depth);
  }

  return os.good ();
}


template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::any_items_have_no_neighbours () const
{
  int items_with_no_neighbours = 0;

  for (int i = 0; i < _number_items; i++)
  {
    if (0 == _results[i].number_elements ())
    {
      items_with_no_neighbours++;
    }
  }

  return items_with_no_neighbours;
}

template <typename N, typename H, typename E>
void
NN_Results<N, H, E>::set_training_set_membership (int s)
{
  for (int i = 0; i < _number_items; i++)
  {
    _results[i].set_training_set_membership (s);
  }

  return;
}

template <typename N, typename H, typename E>
void
NN_Results<N, H, E>::unset_predicted_values ()
{
  for (int i = 0; i < _number_items; i++)
  {
    _results[i].unset_predicted_values ();
  }

  return;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::_initialise_id_index_structures ()
{
  _id_to_index.reserve (_number_items);

  for (int i = 0; i < _number_items; i++)
  {
    const IWString & id = _results[i].id ();

    if (id.contains (' '))
    {
      IWString tmp (id);
      tmp.truncate_at_first (' ');
      _id_to_index[tmp] = i;
    }
    else

    _id_to_index[id] = i;
  }

  return 1;
}

/*
  We need to let each neighbour know where its parent item is
*/

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::initialise_nn_pointers ()
{
  N ** xref = new N *[_number_items];

  if (verbose > 1)
    cerr << "Initialising pool with " << _number_items << " items\n";

  for (int i = 0; i < _number_items; i++)
  {
    xref[i] = _results + i;
  }

  int rc = 1;
  for (int i = 0; i < _number_items; i++)
  {
    if (! _results[i].initialise_nn_pointers ((void **) xref, _number_items))
    {
      rc = 0;
    }
  }

  delete [] xref;

  return rc;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::shuffle_equidistant_neighbours ()
{
  for (int i = 0; i < _number_items; i++)
  {
    _results[i].shuffle_equidistant_neighbours ();
  }

  return 1;
}

/*
  When identifying outliers, we need to keep track of which NN_Item and its
  associated average error
*/

class NN_Item_and_Error
{
  private:
    const NN_Item_Activity * _nn_item;

    experimental_value_t _error;

  public:
    NN_Item_and_Error ();

    void set_nn_item (const NN_Item_Activity * i) { _nn_item = i;}
    const NN_Item_Activity * nn_item () const { return _nn_item;}

    void set_error (experimental_value_t e) { _error = e;}
    experimental_value_t zerror () const { return _error;}
};

NN_Item_and_Error::NN_Item_and_Error ()
{
  _nn_item = nullptr;

  _error = -1.0;

  return;
}

class NN_Results_Activity : public NN_Results<NN_Item_Activity, ID_Activity_Hash, experimental_value_t>
{
  private:

//  When doing Bsquared determinations, we need a list of items sorted by observed activity

    const NN_Item_Activity ** _esorted;

//  The Bsquared needs a couple of arrays

    experimental_value_t * _pred;
    experimental_value_t * _t1;

// private functions

    int _identify_poorly_predicted_items (ostream & output, NN_Item_and_Error * nie) const;
    int _assign_training_test_set_membership_stratified (const Cross_Validation_Conditions & cvc, int);
    int _assign_training_test_set_membership_random (int items_needed);

  public:
    NN_Results_Activity ();
    ~NN_Results_Activity ();

    int initialise_temporary_arrays ();

    int assign_training_test_set_membership (const Cross_Validation_Conditions &, int = 0);

    int compute_model_quality (int number_estimates, BSquared * bsquared, int check_training_set_membership = 1);

    int identify_poorly_predicted_items (ostream &) const;
    int identify_poorly_predicted_items (const char *) const;
};

class NN_Results_Category : public NN_Results<NN_Item_Category, ID_Category_Hash, category_type_t>
{
  private:

// Various arrays needed during training set selection and predictions

  int * _category_count;
  int * _category_errors;
  int * _number_needed_from_category;
  int * _number_selected_from_category;

// If someone enters negative categories, we offset them

  int _category_offset;

// When doing a good_better_best model, we need to know the worst possible error for each
// class. To do this, we need to know the furthest distance from this category to all other categories

  int * _category_range;

// private functions

   int _identify_poorly_predicted_items (ostream & output, const NN_Item_Category ** nim) const;

  public:
    NN_Results_Category ();
    ~NN_Results_Category ();

    int initialise_temporary_arrays ();

    int establish_categories ();

    int assign_number_from_each_category (int percent_training_set);

    int assign_training_test_set_membership (const Cross_Validation_Conditions &, int);

    int compute_model_quality (int number_estimates, BSquared * bsquared, int check_training_set_membership = 1);

    int identify_poorly_predicted_items (ostream &) const;
    int identify_poorly_predicted_items (const char *) const;

    int report_misclassifications_by_class (ostream &) const;
};

NN_Results_Activity::NN_Results_Activity ()
{
  _esorted = nullptr;
  _pred = nullptr;
  _t1 = nullptr;

  return;
}

NN_Results_Activity::~NN_Results_Activity ()
{
  if (NULL != _esorted)
    delete [] _esorted;

   if (NULL != _pred)
     delete [] _pred;

   if (NULL != _t1)
     delete [] _t1;

  return;
}

NN_Results_Category::NN_Results_Category ()
{
  _category_count = nullptr;
  _category_errors = nullptr;
  _number_needed_from_category = nullptr;
  _number_selected_from_category = nullptr;

  _category_range = nullptr;

  return;
}

NN_Results_Category::~NN_Results_Category ()
{
  if (NULL != _category_count)
    delete [] _category_count;

  if (NULL != _category_errors)
    delete [] _category_errors;

  if (NULL != _number_needed_from_category)
    delete [] _number_needed_from_category;

  if (NULL != _number_selected_from_category)
    delete [] _number_selected_from_category;

  if (NULL != _category_range)
    delete [] _category_range;

  return;
}

/*
  All the work of temporary arrays is done in establish_categories ()
*/

int
NN_Results_Category::initialise_temporary_arrays ()
{
  assert (NULL == _category_count);

  return 1;
}

/*
  We need to sort the training set by experimental value
*/

static int 
experimental_value_comparitor (const void * ppv1, const void * ppv2)
{
  const NN_Item_Activity * nni1 = * ((const NN_Item_Activity **) ppv1);
  const NN_Item_Activity * nni2 = * ((const NN_Item_Activity **) ppv2);

  experimental_value_t pv1 = nni1->experimental_value ();
  experimental_value_t pv2 = nni2->experimental_value ();

//cerr << "Comparing '" << nni1->id () << "'  " << pv1 << " with '" << nni2->id () << "' " << pv2 << endl;

  if (pv1 < pv2)
    return 1;

  if (pv1 > pv2)
    return -1;

  return 0;
}

int
NN_Results_Activity::initialise_temporary_arrays ()
{
  assert (NULL == _esorted);

  _esorted = new const NN_Item_Activity * [_number_items];

  for (int i = 0; i < _number_items; i++)
  {
    _esorted[i] = &(_results[i]);
  }

// Sort by experimental value

  qsort (_esorted, _number_items, sizeof (NN_Item_Activity *), experimental_value_comparitor);

#ifdef SHOW_SORTED_VALUES
  for (int i = 0; i < items_in_training_set; i++)
  {
    const NN_Item * nni = _esorted[i];

    experimental_value_t p;
    nni->predicted_value (0, p);

    cerr << " i = " << i << " expt = " << nni->experimental_value () << " pred " << p << '\n';
  }
#endif

  _pred = new experimental_value_t[_number_items];
  _t1 = new experimental_value_t[_number_items];

  return 1;
}

template class resizable_array_p<ID_Dist_Activity>;
template class resizable_array_base<ID_Dist_Activity *>;

template class resizable_array_p<ID_Dist_Category>;
template class resizable_array_base<ID_Dist_Category *>;

/*
  Sort within category and then by misclassifications
*/

static int
nim_comparitor (const void * v1, const void * v2)
{
  const NN_Item_Category * n1 = *((const NN_Item_Category **) v1);
  const NN_Item_Category * n2 = *((const NN_Item_Category **) v2);

  if (n1->category () < n2->category ())
    return -1;

  if (n1->category () > n2->category ())
    return 1;

  if (n1->misclassifications () > n2->misclassifications ())
    return -1;

  if (n1->misclassifications () < n2->misclassifications ())
    return 1;

  return 0;
}

template <typename T>
int
write_outlier_neighbours (const T & t,
                          int n,
                          ostream & output)
{
  if (n > t.number_elements ())
    n = t.number_elements ();

  for (int i = 0; i < n; i++)
  {
    output << "NBR " << t[i]->id () << ' ' << t[i]->distance () << '\n';
  }

  return output.good ();
}

int
NN_Results_Category::_identify_poorly_predicted_items (ostream & output,
                                                const NN_Item_Category ** nim) const
{
  for (int i = 0; i < _number_items; i++)
  {
    nim[i] = &(_results[i]);
  }

  qsort (nim, _number_items, sizeof (NN_Item_Category *), nim_comparitor);

  output << "id misclassifications\n";

  for (int i = 0; i < _number_items; i++)
  {
    const NN_Item_Category * ni = nim[i];

    if (tokens_per_id > 0)
      output << ni->id ();
    else
      write_space_suppressed_string (ni->id (), output);

    output << ' ' << ni->misclassifications () << '\n';

    write_outlier_neighbours (*ni, outlier_neighbours_to_print, output);
  }

  return output.good ();
}

int
NN_Results_Category::identify_poorly_predicted_items (const char * fname) const
{
  std::ofstream output (fname, std::ios::out);

  if (! output.good ())
  {
    cerr << "NN_Results_Category::identify_poorly_predicted_items: cannot open '" << fname << "'\n";
    return 0;
  }

  return identify_poorly_predicted_items (output);
}

int
NN_Results_Category::identify_poorly_predicted_items (ostream & output) const
{
  const NN_Item_Category * * nim = new const NN_Item_Category * [_number_items];

  int rc = _identify_poorly_predicted_items (output, nim);

  delete [] nim;

  return rc;
}

static int
nie_comparitor (const void * v1, const void * v2)
{
  const NN_Item_and_Error * n1 = (const NN_Item_and_Error *) v1;
  const NN_Item_and_Error * n2 = (const NN_Item_and_Error *) v2;

  if (n1->zerror () > n2->zerror ())
    return -1;

  if (n1->zerror () < n2->zerror ())
    return 1;

  return 0;
}

/*
  Sometimes it is helpful to know which molecules are being consistenly mis-predicted
*/

int
NN_Results_Activity::_identify_poorly_predicted_items (ostream & output,
                                        NN_Item_and_Error * nie) const
                                      
{
  for (int i = 0; i < _number_items; i++)
  {
    const NN_Item_Activity & ni = _results[i];

    NN_Item_and_Error & niei = nie[i];

    niei.set_nn_item (&ni);

    experimental_value_t a;

    if (ni.average_absolute_prediction_error (a))
      niei.set_error (a);
    else
      niei.set_error (0.0);
  }

  qsort (nie, _number_items, sizeof (NN_Item_and_Error), nie_comparitor);

  output << "id times_predicted actual ave min max ave ave_error* nndist\n";

  for (int i = 0; i < _number_items; i++)
  {
    const NN_Item_and_Error & niei = nie[i];

    const NN_Item_Activity * ni = niei.nn_item ();

    write_space_suppressed_string (ni->id (), output);
    output.precision (4);
    output << ' ';
    ni->write_range_of_predicted_values (niei.zerror (), output);
    output << '\n';
    write_outlier_neighbours (*ni, outlier_neighbours_to_print, output);
  }

  return 1;
}

int
NN_Results_Activity::identify_poorly_predicted_items (ostream & output) const
{
  NN_Item_and_Error * nie = new NN_Item_and_Error[_number_items];

  int rc = _identify_poorly_predicted_items (output, nie);

  delete [] nie;

  return rc;
}

int
NN_Results_Activity::identify_poorly_predicted_items (const char * fname) const
{
  std::ofstream output (fname, std::ios::out);

  if (! output.good ())
  {
    cerr << "NN_Results_Activity::identify_poorly_predicted_items: cannot open '" << fname << "'\n";
    return 0;
  }

  return identify_poorly_predicted_items (output);
}

template <typename E>
int
echo_id_experimental_value_hash (const ID_Experimental_Value_Hash<E> & idah)
{
  typename ID_Experimental_Value_Hash<E>::const_iterator f;
  for (f = idah.begin (); f != idah.end (); f++)
  {
    cerr << "'" << (*f).first << "' expt = " << (*f).second << '\n';
  }

  return 1;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::_read_experimental_value (const const_IWSubstring & buffer,
                                H & activity)

{
  IWString id;

  if (! buffer.word (experimental_results_identifier_column, id))
  {
    cerr << "Cannot extract identifier from experimental results file\n";
    return 0;
  }

  E expt;

  if (! get_experimental_value (buffer, experimental_results_experimental_value_column, expt))
  {
    if (0 == activity.size())
      return 1;

    return 0;
  }

  if (strip_leading_zeros_from_identifiers)
    id.remove_leading_chars ('0');

  activity[id] = expt;

  return 1;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::determine_experimental_values_token (const const_IWSubstring & zprefix)
{
  int rc = 1;

  for (int i = 0; i < _number_items; i++)
  {
    if (! _results[i].get_experimental_value_from_token (zprefix))
    {
      rc = 0;
    }
  }

  return rc;
}

/*
  Associate names with experimental values
*/

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::_read_experimental_values (iwstring_data_source & input,
                                                H & activity)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (! _read_experimental_value (buffer, activity))
    {
      if (1 == input.lines_read ())    // probably it has a header
        continue;

      cerr << "Cannot process experimental data on line " << input.lines_read () << '\n';
      cerr << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::determine_experimental_values_hash (H & activity)
{
  int rc = _number_items;
  for (int i = 0; i < _number_items; i++)
  {
    if (! _results[i].assign_experimental_result_from_hash (activity))
    {
      rc = 0;
    }
  }

  return rc;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::determine_experimental_values_file (iwstring_data_source & input)
{
  H activity;

  activity.reserve (_number_items);

  if (! _read_experimental_values (input, activity))
  {
    cerr << "Cannot read activities from experimental results file\n";
    return 0;
  }

  if (verbose > 2)
    echo_id_experimental_value_hash (activity);

  return determine_experimental_values_hash (activity);
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::determine_experimental_values_file (const const_IWSubstring & experimental_results_file,
                                    int skip_header_line_in_experimental_values_file)
{
  iwstring_data_source input (experimental_results_file);

  if (! input.ok ())
  {
    cerr << "Cannot open experimental results file '" << experimental_results_file << "'\n";
    return 0;
  }

  if (skip_header_line_in_experimental_values_file)
  {
    const_IWSubstring buffer;
    (void) input.next_record (buffer);
  }

  return determine_experimental_values_file (input);
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::determine_experimental_values_column (int col)
{
  int rc = 1;

  for (int i = 0; i < _number_items; i++)
  {
    if (! _results[i].get_experimental_values_from_identifier_token (col))
    {
      rc = 0;
    }
  }

  if (0 == rc)
    return 0;

  return 1;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::_assign_training_test_set_membership_external (const Cross_Validation_Conditions & cvc,
                                                                    int i)
{
  if (0 == _id_to_index.size ())
    _initialise_id_index_structures ();

  IWString_STL_Hash_Set its;

  if (! cvc.external_split (i, its))
  {
    cerr << "NN_Results::_assign_training_test_set_membership_external: cannot process split " << i << '\n';
    return 0;
  }

  int nmissing = 0;

  for (IWString_STL_Hash_Set::const_iterator i = its.begin (); i != its.end (); ++i)
  {
    const IWString & id = (*i);

    IW_STL_Hash_Map_int::const_iterator f = _id_to_index.find (id);

    if (f == _id_to_index.end ())
    {
      cerr << "NN_Results::_assign_training_test_set_membership_external: no data on '" << id << "'\n";
      nmissing++;
    }
    else
    {
      int j = (*f).second;
      _results[j].set_training_set_membership (1);
    }
  }

  if (nmissing)
  {
    cerr << nmissing << " of " << its.size () << " items missing. Set contains " << _number_items << " items, id index hash " << _id_to_index.size () << '\n';
    return 0;
  }
  else 
    return 1;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::compute_predicted_values (const BSquared & bs,
                                               int which_prediction)
{
  const NN_Weight_Function & wfn = bs.weight_function ();

  for (int i = 0; i < _number_items; i++)
  {
    N & nni = _results[i];

    if (nni.in_training_set ())
      continue;

    nni.compute_predicted_values (bs, which_prediction, wfn);
  }

  return 1;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::compute_predicted_values (const BSquared * bsquared,
                                               int number_predictions)
{
  int rc = 1;
  for (int i = 0; i < number_predictions; i++)
  {
    if (! compute_predicted_values (bsquared[i], i))
    {
      cerr << "Cannot compute predicted values for '" << bsquared[i].estimate_type () << "'\n";
      rc = 0;
    }
  }

  return rc;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::compute_leave_one_out_values (const BSquared * bs,
                                                   int number_predictions)
{
  set_training_set_membership (1);     // everything in the training set

  for (int i = 0; i < number_predictions; i++)
  {
    compute_leave_one_out_values (bs[i], i);
  }

  return 1;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::compute_leave_one_out_values (const BSquared & bs,
                                                   int which_prediction)
{
  const NN_Weight_Function & wfn = bs.weight_function ();

  for (int i = 0; i < _number_items; i++)
  {
    _results[i].compute_predicted_values (bs, which_prediction, wfn);
  }

  return 1;
}

/*
  When doing categorical studies, we need to establish how many categories there
  are, and assign category markers. 

  We need to determine the number of instances of each category and then
  set up the NUMBER_NEEDED_FROM_CATEGORY array

  We are setting up both member variables and some file scope variables

  Categories must be non negative, so translate them
*/

int
NN_Results_Category::establish_categories ()
{
  int lowest_category = 0;

  for (int i = 0; i < _number_items; i++)
  {
    category_type_t c = _results[i].category ();

    if (c < lowest_category)
      lowest_category = c;
  }

  _category_offset = lowest_category;

  if (verbose)
    cerr << "Catogories offset by " << _category_offset << '\n';

  std::unordered_map<int, int> number_of_instances_of_each_category;

  for (int i = 0; i < _number_items; i++)
  {
    category_type_t c = _results[i].category ();

    _results[i].set_category(c - _category_offset);

    c = _results[i].category ();

    assert (c >= 0);

    number_of_instances_of_each_category[c]++;

    if (c > highest_category_number)
      highest_category_number = c;
  }

  number_categories = number_of_instances_of_each_category.size ();

  if (number_categories < 2)     // turn on verbose so they can see the category info
    verbose = 1;

  if (verbose)
    cerr << "Found " << number_categories << " separate categories, highest category number " << highest_category_number << '\n';

  number_in_category = new_int (highest_category_number + 1);

  for (std::unordered_map<int, int>::const_iterator i = number_of_instances_of_each_category.begin (); i != number_of_instances_of_each_category.end (); ++i)
  {
    int category = (*i).first;

    if (verbose)
      cerr << "Found " << (*i).second << " instances of category " << category << '\n';

    number_in_category[category] = (*i).second;
  }

  if (number_categories < 2)
  {
    cerr << "Yipes, only found " << number_categories << " different categories, cannot continue\n";
    return 0;
  }

  category_weight = new weight_t[highest_category_number + 1];

  _category_count = new int[highest_category_number + 1];
  _category_errors = new int[highest_category_number + 1];

  _number_needed_from_category = new_int (highest_category_number + 1);
  _number_selected_from_category = new int[highest_category_number + 1];

  category_sumw = new weight_t[highest_category_number + 1];

  weight_t w = 1.0 / static_cast<weight_t> (number_categories);

  set_vector (category_weight, highest_category_number + 1, w);

  if (good_better_best_model)
  {
    _category_range = new_int (highest_category_number + 1);

    category_type_t first_category = -1;
    for (int i = 0; i <= highest_category_number; i++)
    {
      if (0 == number_in_category[i])
        continue;

      first_category = i;
      break;
    }

    assert (first_category >= 0);

    for (int i = 0; i <= highest_category_number; i++)
    {
      if (0 == number_in_category[i])
        continue;

      int r1 = i - first_category;
      int r2 = highest_category_number - i;

      if (r1 >= r2)
        _category_range[i] = r1;
      else
        _category_range[i] = r2;
    }
  }

  return 1;
}

int
NN_Results_Category::assign_number_from_each_category (int percent_training_set)
{
  assert (percent_training_set > 0 && percent_training_set < 100);
  assert (NULL != number_in_category);
  assert (highest_category_number > 0);

  for (int i = 0; i <= highest_category_number; i++)
  {
    if (0 == number_in_category[i])
      continue;

    _number_needed_from_category[i] = int ( float (number_in_category[i] * percent_training_set) / 100.0);

    if (0 == _number_needed_from_category[i])
    {
      _number_needed_from_category[i] = 1;
      cerr << "Warning, at " << percent_training_set << " percent training set class " << i << " gets only one item in the training set\n";
    }

    if (verbose > 1)
      cerr << "At " << percent_training_set << " percent, category " << i << " will have " << _number_needed_from_category[i] << " of " << number_in_category[i] << " items\n";
  }

  return 1;
}

/*
  Pull the appropriate number of samples from each category
*/

int
NN_Results_Category::assign_training_test_set_membership (const Cross_Validation_Conditions & cvc,
                                                          int i)
{
  if (cvc.number_elements ())
    return _assign_training_test_set_membership_external (cvc, i);

  set_vector (_number_selected_from_category, highest_category_number + 1, 0);

  int categories_completed = 0;

  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<int> u(0, _number_items - 1);

  while (1)
  {
    int i = u(rng);

    NN_Item_Category & ri = _results[i];

    if (ri.in_training_set ())
      continue;

    category_type_t c = ri.experimental_value ();

    if (_number_selected_from_category[c] == _number_needed_from_category[c])
      continue;

    ri.set_training_set_membership (1);
    _number_selected_from_category[c]++;

    if (_number_selected_from_category[c] == _number_needed_from_category[c])
    {
      categories_completed++;
      if (categories_completed == number_categories)
        break;
    }
  }

  return 1;
}

/*
  We do stratified sampling if we can. But, if the interval (delta) is 1,
  then we revert to random sampling
*/

int
NN_Results_Activity::_assign_training_test_set_membership_stratified (const Cross_Validation_Conditions & cvc,
                                   int items_needed)
{
  int delta = _number_items / items_needed;

  if (1 == delta)
  {
    cerr << "Cannot do stratified sampling, doing random...\n";
    return _assign_training_test_set_membership_random (items_needed);
  }

  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<int> u(0, delta - 1);

  for (int i = 0; i < _number_items; i += delta)
  {
    //int j = i + intbtwij (0, delta - 1);
    int j = i + u(rng);

    if (j >= _number_items)
      break;

//  need to use the _esorted array here - this really doesn't work...

    _results[i].set_training_set_membership (1);
  }

  return 1;
}

/*
  We want to assign molecules to the training or test set. 
*/

int
NN_Results_Activity::assign_training_test_set_membership (const Cross_Validation_Conditions & cvc,
                                                          int i)
{
  if (cvc.number_elements ())
    return _assign_training_test_set_membership_external (cvc, i);

  int items_needed = _number_items * cvc.percent_training_set () / 100;

  if (0 == items_needed)
  {
    cerr << "Very small training set\n";
    items_needed = 1;
  }

  if (cvc.stratified_sampling ())
    return _assign_training_test_set_membership_stratified (cvc, items_needed);

  return _assign_training_test_set_membership_random (items_needed);
}

int
NN_Results_Activity::_assign_training_test_set_membership_random (int items_needed)
{
  int max_iterations = _number_items;
  if (items_needed > _number_items / 2)
    max_iterations = 4 * _number_items;

  int items_selected = 0;

  int iterations = 0;

  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<int> u(0, _number_items - 1);

  while (items_selected < items_needed)
  {
    iterations++;

    if (iterations > max_iterations)
    {
      cerr << "Cannot find " << items_needed << " items from " << _number_items << " set\n";
      return 0;
    }

    // int i = intbtwij (0, _number_items - 1);
    int i = u(rng);

    if (_results[i].in_training_set ())
      continue;

    _results[i].set_training_set_membership (1);
    items_selected++;
  }

//#define CHECK_TS_MEMBERSHIP
#ifdef CHECK_TS_MEMBERSHIP
  cerr << items_selected << " of " << items_in_whole_set << " items in training set\n";
  int count = 0;
  for (int i = 0; i < _number_items; i++)
  {
    if (_result[i].in_training_set ())
    {
      count++;
    }
  }

  if (count != items_selected)
    cerr << "Mismatch in assign_training_test_set_membership, count = " << count << " nsel = " << items_selected << '\n';
#endif

  return items_selected;
}

/*
  Compute the model quality for each of the estimates

  We invoke this two times. Once when we are doing leave-one-out
  estimates, and that time we don't check training set membership.
  Later when doing cross validation, we need to check training
  set membership
*/

int
NN_Results_Category::compute_model_quality (int number_estimates,
                                            BSquared * bsquared,
                                            int check_training_set_membership)
{
  for (int i = 0; i < number_estimates; i++)
  {
    set_vector (_category_count,  highest_category_number + 1, 0);
    set_vector (_category_errors, highest_category_number + 1, 0);

    int number_predicted_values = 0;
    for (int j = 0; j < _number_items; j++)
    {
      NN_Item_Category & nnj = _results[j];

      if (! check_training_set_membership)
        ;
      else if (nnj.in_training_set ())
        continue;

      category_type_t predicted_category;
      if (! nnj.predicted_value (i, predicted_category))    // maybe all neighbours also in training set
        continue;

      number_predicted_values++;

      int experimental_category = nnj.category ();

      assert (experimental_category >= 0 && experimental_category <= highest_category_number);

      _category_count[experimental_category]++;

      if (predicted_category == experimental_category)    // great, correctly predicted it
        continue;

      nnj.another_misclassification ();

      if (good_better_best_model)
      {
        if (predicted_category > experimental_category)
          _category_errors[experimental_category] += predicted_category - experimental_category;
        else
          _category_errors[experimental_category] += experimental_category - predicted_category;
      }
      else
        _category_errors[experimental_category]++;
    }

    if (number_predicted_values <= 2)
    {
      cerr << "For " << bsquared[i].estimate_type () << " only " << number_predicted_values << " predicted values. B2 set to 0.0\n";
      bsquared[i].extra (0.0);
      continue;
    }

    if (overall_error_rate_is_error_rate_for_particular_class >= 0)
    {
      int c = overall_error_rate_is_error_rate_for_particular_class;

      double overall_score = 1.0 - static_cast<double> (_category_errors[c]) / static_cast<double> (_category_count[c]);

      bsquared[i].extra (overall_score);

      continue;
    }

    double overall_score = 1.0;

    for (int j = 0; j <= highest_category_number; j++)
    {
      if (0 == _category_errors[j])
        continue;

      double e;
      if (good_better_best_model)
        e = static_cast<double> (_category_errors[j]) / static_cast<double> (number_predicted_values * _category_range[j]);
      else
        e = static_cast<double> (_category_errors[j]) / static_cast<double> (_category_count[j]);

//#define DEBUG_COMPUTE_CATEGORY_ERROR_RATES
#ifdef DEBUG_COMPUTE_CATEGORY_ERROR_RATES
      cerr << "Category " << j << " " << _category_count[j] << " items in category, " << _category_errors[j] << " incorrect " << e << '\n';
#endif

      overall_score = overall_score - e * category_weight[j];
    }

    bsquared[i].extra (overall_score);
  }

  return 1;
}

int
NN_Results_Category::report_misclassifications_by_class (ostream & os) const
{
  os << "Errors by category\n";

  for (int i = 0; i <= highest_category_number; i++)
  {
    if (0 == _category_count[i])
      continue;

    float tmp = 1.0 - static_cast<float> (_category_errors[i]) / static_cast<float> (_category_count[i]);
    os << " category " << i << " contains " << _category_count[i] << " items " << (_category_count[i] - _category_errors[i]) << " correct " << tmp << '\n';
  }

  return os.good ();
}

/*
  We have some predicted values, compute the Bsquared values
*/

int
NN_Results_Activity::compute_model_quality (int number_estimates,
                      BSquared * bsquared,
                      int check_training_set_membership)
{
  for (int i = 0; i < number_estimates; i++)
  {
    int myn = 0;                          // the number of values we get
    for (int j = 0; j < _number_items; j++)
    {
      const NN_Item_Activity * nnj = _esorted[j];

      if (! check_training_set_membership)
        ;
      else if (nnj->in_training_set ())
        continue;

      experimental_value_t e;
      if (! nnj->predicted_value (i, e))    // maybe all neighbours also in training set
        continue;

      _pred[myn] = e;
      myn++;
    }

    double myb2;
    if (myn <= 2)
    {
      cerr << "For " << bsquared[i].estimate_type () << " only " << myn << " predicted values. B2 set to 0.0\n";
      myb2 = 0.0;
    }
    else
      compute_b_squared (_pred, _t1, myn, myb2);

    bsquared[i].extra (myb2);
  }

  return 1;
}

template <typename N, typename H, typename E>
int
NN_Results<N, H, E>::write_predicted_values (ostream & output,
                                             int ndx) const
{
  output << "Name True Predicted\n";

  for (int i = 0; i < _number_items; i++)
  {
    const N & nni = _results[i];

    write_space_suppressed_string (nni.id (), output);
    output << ' ' << nni.experimental_value () << ' ';

    E e;
    if (! nni.predicted_value (ndx, e))
      output << missing_value;
    else
      output << e;

    output << '\n';
  }

  return output.good ();
}

template <typename R>
int
write_test_set_predictions (R & nnr,
                            BSquared & bs, 
                            ostream & output)
{
  nnr.compute_leave_one_out_values (bs, 0);

  return nnr.write_predicted_values (output, 0);
}

template <typename R>
int
write_test_set_predictions (R & nnr,
                            BSquared & bs, 
                            IWString & file_name_for_test_set_predictions)
{
  std::ofstream output (file_name_for_test_set_predictions.null_terminated_chars (), std::ios::out);

  if (! output.good ())
  {
    cerr << "Cannot open '" << file_name_for_test_set_predictions << "'\n";
    return 0;
  }

  return write_test_set_predictions (nnr, bs, output);
}

/*
  Roughly where is the T value that gives a P value of pvalue_threshold
*/

static int
determine_tvalue_threshold (int n,
                            double pvalue_threshold)
{
  double right = 8.0;
  double rp = iwpvalue (n, right);

  if (rp > pvalue_threshold)
  {
    cerr << "Huh, for N = " << n << " T = 8.0, pvalue " << rp << " which is > " << pvalue_threshold << '\n';
    cerr << "Too strict!\n";
    return 0;
  }

  double left = 0.01;
  // double lp = iwpvalue (n, left);
  //cerr << " lp " << lp << " threshold " << pvalue_threshold << " rp " << rp << '\n';
  //assert (lp >= pvalue_threshold && pvalue_threshold >= rp);

  while (1)
  {
    double xmid = (left + right) * 0.5;

    double xp = iwpvalue (n, xmid);

    if (verbose > 2)
      cerr << " X = " << xp << " P " << xp << '\n';

    if (fabs (xp - pvalue_threshold) < 0.00001)
    {
      tvalue_threshold = xp;
      if (verbose)
        cerr << "For N = " << n << " P = " << pvalue_threshold << " T threshold " << tvalue_threshold << '\n';

      return 1;
    }

    if (xp < pvalue_threshold)
    {
      right = xmid;
    }
    else
    {
      left = xmid;
    }
  }

  return 1;
}

/*
  Special case when shortest_useful_distance has been set
*/

static int
print_best_models (const BSquared * bsquared,
                   int number_predictions,
                   ostream & os,
                   double min_dist)
{
  assert (min_dist == shortest_useful_distance);

  int found_model = 0;
  for (int i = 0; i < number_predictions; i++)
  {
    const BSquared & bs = bsquared[i];
    if (bs.radius () >= shortest_useful_distance)
    {
      os << "Shortest useful distance\n";
      bs.report (os);
      found_model = 1;
      break;
    }
  }

  if (! found_model)
  {
    os << "No models with radius >= " << shortest_useful_distance << '\n';
    cerr << "No models with radius >= " << shortest_useful_distance << '\n';

    return os.good ();
  }

  int models_printed = 0;

  for (int i = 0; i < number_predictions; i++)
  {
    const BSquared & bs = bsquared[i];
    if (bs.radius () >= shortest_useful_distance)
    {
      os << KNN_MODEL_PREFIX << ' ';
      bs.report (os);

      models_printed++;

      if (models_printed >= number_predictions_to_print)
        break;
    }
  }

  return os.good ();
}

static int
print_best_models (const BSquared * bsquared,
                   int number_predictions,
                   ostream & os)
{
  bsquared[0].report (os, 1);   // extra arg means comment prefix

  if (shortest_useful_distance >= 0.0)     // needs special treatment
    return print_best_models (bsquared, number_predictions, os, shortest_useful_distance);

  for (int i = 0; i < number_predictions_to_print; i++)
  {
    os << KNN_MODEL_PREFIX << ' ';
    bsquared[i].write_nn_prediction_options (os);
    os << '\n';
  }

  return os.good ();
}

template <typename R>
int
do_cross_validations (R & nnr,
                      int number_estimates,
                      BSquared * bsquared,
                      const Cross_Validation_Conditions & cvc)
{
  int items_in_training_set = nnr.number_elements ();

  if (items_in_training_set < 2)
  {
    cerr << "Sorry, cannot do cross validation on " << items_in_training_set << " items\n";
    return 0;
  }
  else if (items_in_training_set < 16)
    cerr << "Warning!, only " << items_in_training_set << " items in your training set!!!\n";

  for (int i = 0; i < number_estimates; i++)
  {
    bsquared[i].remove_all_previous_prediction_data (0);
  }

  for (int i = 0; i < cvc.number_splits (); i++)
  {
    if (verbose > 1)
      cerr << "Beginning test/train split " << i << '\n';

    if (0 == i % cvc.reorder_equidistant_neighbours ())
      nnr.shuffle_equidistant_neighbours ();

    nnr.set_training_set_membership (0);     // everything out of the training set

    if (! nnr.assign_training_test_set_membership (cvc, i))
      return 0;

    nnr.unset_predicted_values ();

    (void) nnr.compute_predicted_values (bsquared, number_estimates);   // based on the new training set membership

    nnr.compute_model_quality (number_estimates, bsquared);
  }

  if (cvc.number_splits () > 1)     // so we have average predictions
  {
    qsort (bsquared, number_estimates, sizeof (BSquared), average_bsquared_comparitor);

    for (int i = 0; i < number_estimates; i++)
    {
      bsquared[i].set_sequence (i);
    }
  }

  if (verbose)
  {
    int nprint;
    if (cvc.do_paired_t_test ())
    {
      nprint = prediction_types_to_print;
      if (nprint > number_estimates)
        nprint = number_estimates;

      cerr << "Best " << nprint << " estimates, without paired T test\n";
    }
    else
    {
      nprint = number_estimates;

      cerr << "Best estimates based on average Bsquared value\n";
    }

    for (int i = 0; i < nprint; i++)
    {
      const BSquared & bs = bsquared[i];
      bs.report (cerr);
    }

    cerr << '\n';
  }

  if (cvc.do_paired_t_test ())
  {
#ifdef PRINT_ARRAY
    cerr << "Starting sort\n";
    for (int i = 0; i < number_estimates; i++)
    {
      cerr << "Estimate " << i << " type " << bsquared[i].estimate_type () << " at " << hex << &(bsquared[i]) << dec << '\n';
    }
    cerr << "Addresses between " << hex << bsquared << " and " << hex << &(bsquared[number_estimates - 1]) << dec << '\n';
#endif

    do_all_pair_wise_comparisons (bsquared, number_estimates);

    qsort (bsquared, number_estimates, sizeof (BSquared), bsquared_comparitor_number_better);

    if (verbose)
      cerr << "After paired T test of Bsquared values, threshold " << pvalue_threshold << '\n';

    if (verbose > 2)
      print_comparison_details = 1;

    int nprint;
    if (verbose <= 1)
    {
      nprint = prediction_types_to_print;
      if (nprint > number_estimates)
        nprint = number_estimates;
    }
    else
    {
      nprint = number_estimates;
    }

    for (int i = 0; i < nprint; i++)
    {
      const BSquared & bs = bsquared[i];
      bs.report (cerr);

      if (verbose > 1 && i > 0)
        bsquared_comparitor_careful (&(bsquared[i - 1]), &bs);
    }

    if (verbose && careful_comparisons_done)
      cerr << "Performed " << careful_comparisons_done << " comparisons, " << careful_comparisons_resolved_by_complexity << " resolved by complexity\n";
  }

  if (stream_for_best_models.rdbuf ()->is_open ())
    print_best_models (bsquared, number_estimates, stream_for_best_models);

  std::cout << "After cross validation\n";
  print_best_models (bsquared, number_estimates, std::cout);
  cerr << "After cross validation\n";
  print_best_models (bsquared, number_estimates, cerr);

  if (file_name_for_test_set_predictions.length ())
    write_test_set_predictions (nnr, bsquared[0], file_name_for_test_set_predictions);

#ifdef DEBUG_TRAINING_SET_MEMBERSHIP
  for (int i = 0; i < items_in_training_set; i++)
  {
    cerr << "In training set " << nnr[i].times_in_training_set () << '\n';
  }
#endif

  return 1;
}

static int
write_best_results (ostream & os,
                    double maxb2,
                    const BSquared & best_estimate,
                    const Accumulator<double> & bvals)
{
  os << "Leave One Out " << bsquared_or_error_rate << " " << maxb2 << " " << best_estimate.estimate_type () << '\n';
  os << "values between " << bvals.minval () << " and " << bvals.maxval ();
  if (bvals.n () > 1)
    os << " average " << bvals.average ();
  os << '\n';

  return os.good ();
}

static void
display_cross_validation_options ()
{
  cerr << " -V times=nn      perform NN training/test set splits\n";
  cerr << " -V <n>           use N% of the dataset as the training set\n";
//cerr << " -V strat         do stratified sampling throughout activity range\n";
  cerr << " -V seed=nn       use nn as the seed for the random number generator\n";
  cerr << " -V pvalue=xx     do paired T test with P value threshold xx when comparing Bsquared values\n";
  cerr << " -V eqshuf=nn     every <nn> splits, randomly shuffle equidistant neighbours\n";

  return;
}

int
Cross_Validation_Conditions::construct_from_command_line (Command_Line & cl,
                                                          char c)
{
#ifdef SEEDS_NO_LONGER_SUPPORTED
  int seed_set = 0;
#endif

  int i = 0;
  const_IWSubstring v;
  while (cl.value (c, v, i++))
  {
    if (v.starts_with ("times="))
    {
      v.remove_leading_chars (6);
      if (! v.numeric_value (_number_splits) || _number_splits < 1)
      {
        cerr << "Invalid number of splits '" << v << "'\n";
        usage (5);
      }

      if (verbose)
        cerr << "Will do " << _number_splits << " splits of the training set\n";
    }
    else if (v.starts_with ("strat"))
    {
      _stratified_sampling = 1;

      if (verbose)
        cerr << "Will do stratified sampling if possible\n";
    }
    else if (v.starts_with ("seed="))
    {
      v.remove_leading_chars (5);

      long s;      // unfortunately, we don't have a numeric_value () method for type random_number_seed_t
      if (! v.numeric_value (s))
      {
        cerr << "Invalid seed value '" << v << "'\n";
        return 0;
      }

#ifdef SEEDS_NO_LONGER_SUPPORTED
      iw_set_rnum_seed (s);
      seed_set = 1;
#endif
    }
    else if (v.starts_with ("pvalue="))
    {
      v.remove_leading_chars (7);
      if (! v.numeric_value (pvalue_threshold) || pvalue_threshold <= 0.0 || pvalue_threshold >= 1.0)
      {
        cerr << "Invalid P value specifier '" << v << "'\n";
        usage (4);
      }

      if (verbose)
        cerr << "Will compare estimates via paired T test of Bsquared values, threshold " << pvalue_threshold << '\n';

      _do_paired_t_test = 1;
    }
    else if (v.starts_with ("eqshuf="))
    {
      v.remove_leading_chars (7);

      int tmp;
      if (! v.numeric_value (tmp) || tmp < 1)
      {
        cerr << "The reorder equidistant neigbhbours directive must be a whole positive number\n";
        usage (6);
      }

      if (verbose)
        cerr << "Will reorder equidistant neighbours every " << tmp << " splits\n";

      _reorder_equidistant_neighbours = tmp;
    }
    else if ("help" == v)
    {
      display_cross_validation_options ();
      exit (1);
    }
    else
    {
      if (v.ends_with ('%'))
        v.chop ();

      if (! v.numeric_value (_percent_training_set) || _percent_training_set < 1 || _percent_training_set > 99)
      {
        cerr << "INvalid training set percentage '" << v << "'\n";
        usage (5);
      }

      if (verbose)
        cerr << "Will use " << _percent_training_set << " percent of the input for training\n";
    }
  }

  int splits_on_command_line = cl.number_elements () - 1;

  if (_number_splits && splits_on_command_line)
  {
    cerr << "Cross_Validation_Conditions::construct_from_command_line: cannot specify both 'times=' and training set files\n";
    return 0;
  }

  if (splits_on_command_line)
  {
    if (! Training_Sets_From_Files::initialise_external_splits (cl))
    {
      cerr << "Cross_Validation_Conditions::construct_from_command_line:cannot initialise external splits\n";
      return 0;
    }

    _number_splits = _number_elements;
  }
  else
  {
    if (_do_paired_t_test && 0 == _number_splits)
    {
      _number_splits = 30;

      cerr << "Paired T test requested, but nsplits not specified. Set to " << _number_splits << '\n';
    }

#ifdef SEEDS_NO_LONGER_SUPPORTED
    if (! seed_set)
      iw_random_seed ();
#endif

//  if (_stratified_sampling)
//    use_the_categories_for_sampling ();
  }

  if (_do_paired_t_test)  
    determine_tvalue_threshold (_number_splits, pvalue_threshold);

  if (0 == _reorder_equidistant_neighbours)
    _reorder_equidistant_neighbours = _number_splits + 1; // so the modulus test above never happens

  return 1;
}

template <typename R>
int
do_leave_one_out (R & nnr,
                  int number_estimates,
                  BSquared * bsquared)
{
  assert (number_estimates > 0);
  if (verbose > 1)
    cerr << "Doing " << number_estimates << " estimates\n";

  for (int i = 0; i < nnr.number_results (); i++)
  {
    nnr[i].size_predicted_array (number_estimates);
  }

  if (! nnr.compute_leave_one_out_values (bsquared, number_estimates))
  {
    cerr << "Cannot compute leave one out values\n";
    return 1;
  }

  int items_in_training_set = nnr.number_elements ();

  cerr << "Leave one out B squared values for the whole dataset - best case\n";
  cerr << items_in_training_set << " items in the training set\n";

  nnr.compute_model_quality (number_estimates, bsquared, 0);

  Accumulator<double> accb;

  double maxb2 = -2.0;
  int    imaxb2 = -1;

  for (int i = 0; i < number_estimates; i++)
  {
    double bsq = bsquared[i].bsquared_value (0);

    cerr << i << ' ' << bsquared_or_error_rate << " for " << bsquared[i].estimate_type () << " : " << bsq << '\n';
  
    accb.extra (bsq);
  
    if (bsq > maxb2)
    {
      maxb2 = bsq;
      imaxb2 = i;
    }
  }

  if (imaxb2 < 0)
  {
    cerr << "No predicted values found\n";
    return 3;
  }

  write_best_results (cerr, maxb2, bsquared[imaxb2], accb);
  write_best_results (std::cout, maxb2, bsquared[imaxb2], accb);

  nnr.unset_predicted_values ();

  return 1;
}

static int
read_neighbour_conditions_file (iwstring_data_source & input,
                                resizable_array_p<IWString> & model_conditions)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (0 == buffer.length ())
      continue;

    if (buffer.starts_with ('#'))
      continue;

    buffer.truncate_at_first (' ');

    IWString * tmp = new IWString (buffer);
    model_conditions.add (tmp);
  }

  return model_conditions.number_elements ();
}

static int
read_neighbour_conditions_file (const const_IWSubstring & n,
                           resizable_array_p<IWString> & model_conditions)
{
  const_IWSubstring fname (n);
  fname.remove_up_to_first ('=');

  if (0 == fname.length ())
  {
    cerr << "Must specify file name after file=\n";
    return 0;
  }

  iwstring_data_source input (fname);
  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (! input.skip_records (1))
    cerr << "Cannot read first record\n";

  return read_neighbour_conditions_file (input, model_conditions);
}

static int
read_neighbour_conditions (const const_IWSubstring & n,
                           resizable_array_p<IWString> & model_conditions)
{
  if (n.starts_with ("file="))
    return read_neighbour_conditions_file (n, model_conditions);

  IWString * tmp = new IWString (n);
  model_conditions.add (tmp);

  return 1;
}

/*
  Parse the -N option - a list of model conditions to build
*/

static int
read_models_to_be_built (Command_Line & cl,
                         int & number_estimates,
                         BSquared * & bsquared)
{
  int nN = cl.option_count ('N');

  assert (nN > 0);

  resizable_array_p<IWString> model_conditions;

  for (int i = 0; i < nN; i++)
  {
    const_IWSubstring n;
    cl.value ('N', n, i);

    if (! read_neighbour_conditions (n, model_conditions))
    {
      cerr << "Cannot read model specifications from '" << n << "'\n";
      return 0;
    }
  }

  number_estimates = model_conditions.number_elements ();

  if (0 == number_estimates)
  {
    cerr << "Yipes, no models read from -N option(s)\n";
    return 0;
  }

  bsquared = new BSquared [number_estimates];

  int rc = 1;
  for (int i = 0; i < number_estimates; i++)
  {
    const IWString & nc = *(model_conditions[i]);

    if (! bsquared[i].build_from_text_string (nc))
    {
      cerr << "Cannot parse neighbour specifications from '" << nc << "'\n";
      rc = 0;
    }
  }

  return rc;
}

static int
initialise_dm_to_nn_conditions (DM_to_NN_Conditions<float> & dmc,
                                const BSquared * bsquared,
                                int number_estimates,
                                int percent_training_set)
{
  similarity_type_t max_distance_encountered = 0.0;
  int max_neighbours_needed = 0;

  for (int i = 0; i < number_estimates; i++)
  {
    const BSquared & bs = bsquared[i];

    if (bs.radius () > max_distance_encountered)
      max_distance_encountered = bs.radius ();

    if (bs.number_neighbours () > max_neighbours_needed)
      max_neighbours_needed = bs.number_neighbours ();
  }

  if (verbose)
    cerr << "Neighbours required by models to " << max_neighbours_needed << '\n';

  dmc.set_min_neighbours (1);
  if (max_neighbours_needed > 1)
  {
    max_neighbours_needed = 100 * max_neighbours_needed / percent_training_set;
    if (max_neighbours_needed < 2)   // hard to imagine
      max_neighbours_needed = 2;

    dmc.set_max_neighbours (max_neighbours_needed);
  }

  if (verbose)
    cerr << "Max distance among model conditions " << max_distance_encountered << '\n';

  if (max_distance_encountered > 0.0)
    dmc.set_max_distance (max_distance_encountered);

  return 1;
}

template <typename R>
int
establish_tokens_per_id (const R & nnr)
{
  tokens_per_id = 0;

  int n = nnr.number_elements ();

  for (int i = 0; i < n; i++)
  {
    const IWString & id = nnr[i].id ();

    int nw = id.nwords ();

    if (0 == i)
      tokens_per_id = nw;
    else if (nw == tokens_per_id)
      ;
    else
    {
      tokens_per_id = 0;
      return 0;
    }
  }

  return 1;
}


template <typename R>
int
initialisation_nnr (Command_Line & cl,
                    R & nnr,
                    const BSquared * bsquared,
                    int number_estimates,
                    int percent_training_set)
{
#ifdef IW_DASH_S
  if (cl.option_present ('s'))
  {
    int nprocess;
    if (! cl.value ('s', nprocess) || nprocess < 2)
    {
      cerr << "The number of molecules to process (-s option) must be a whole positive number\n";
      usage (13);
    }

    if (! nnr.size_the_problem (nprocess))
    {
      cerr << "Bad news, cannot size the problem for " << nprocess << " molecules\n";
      return 14;
    }

    if (verbose)
      cerr << "Problem sized for " << nprocess << " molecules\n";
  }
#endif

  int rc;
  if (nn_file_is_distance_matrix)
  {
    DM_to_NN_Conditions<float> dmc;
    initialise_dm_to_nn_conditions (dmc, bsquared, number_estimates, percent_training_set);

    IWString_STL_Hash_Set only_use;    // not used

    rc = nnr.build_from_distance_matrix (cl[0], dmc, only_use);
  }
  else
    rc = nnr.build (cl[0]);

  if (0 == rc)
  {
    cerr << "Cannot build pool from '" << cl[0] << "'\n";
    return 0;
  }

  if (nnr.any_items_have_no_neighbours ())
  {
    cerr << nnr.any_items_have_no_neighbours () << " members of the training set have no neighbours, impossible\n";
    return 0;
  }

  if (verbose)
    cerr << "Read near neighbour results for " << nnr.number_results () << " molecules\n";

  establish_tokens_per_id (nnr);

  if (cl.option_present ('E'))
  {
    const_IWSubstring e = cl.string_value ('E');

    if (stream_for_best_models.rdbuf ()->is_open ())
      stream_for_best_models << EXPT_SOURCE << ' ' << e << '\n' << '\n';

    int rc = 1;

    if (e.starts_with ("col="))
    {
      e.remove_leading_chars (4);
      int experimental_column;

      if ("last" == e)
        experimental_column = -1;
      else if (! e.numeric_value (experimental_column) || experimental_column < 1)
      {
        cerr << "Invalid experimental column 'col=" << e << "'\n";
        usage (7);
      }

      experimental_column--;    // convert to our notation

      rc = nnr.determine_experimental_values_column (experimental_column);
    }
    else if (e.starts_with ("file="))
    {
      e.remove_leading_chars (5);

      rc = nnr.determine_experimental_values_file (e, cl.option_present ('j'));
    }
    else
    {
      rc = nnr.determine_experimental_values_token (e);
    }

    if (0 == rc)
    {
      cerr << "Cannot determine experimental values (-E option)\n";
      return 0;
    }
  }

// Very important that initialise_nn_pointers () be called after the experimental
// values are associated with the pool objects

  if (! nnr.initialise_nn_pointers ())
  {
    cerr << "Bad news, must be bad identifiers in input. Get help\n";
    return 0;
  }

  return 1;
}

/*
  If we have error balancing in effect, the computations get initialised differently
*/

/*static int
initialise_bsquared_array (int number_weight_functions,
                           const NN_Weight_Function * weight_function,
                           const resizable_array<int> & neighbours_to_use,
                           const resizable_array<int> & min_neighbours_required,
                           const resizable_array<float> & distances_to_use,
                           category_type_t special_class,
                           const resizable_array<probability_type_t> & probabilities_to_use,
                           int & number_estimates,
                           BSquared * & bsquared)
{
  assert (probabilities_to_use.number_elements () > 0);

  int number_neighbour_types = neighbours_to_use.number_elements () + min_neighbours_required.number_elements () * distances_to_use.number_elements ();
  number_estimates = number_neighbour_types * number_weight_functions;
  number_estimates = number_estimates * probabilities_to_use.number_elements ();

  if (verbose && number_estimates > 1)
    cerr << "A total of " << number_estimates << " estimates will be made\n";

  assert (number_estimates > 0);

  bsquared = new BSquared[number_estimates];

  int ndx = 0;
  for (int i = 0; i < number_weight_functions; i++)
  {
    const NN_Weight_Function & wi = weight_function[i];

    for (int j = 0; j < neighbours_to_use.number_elements (); j++)
    {
      for (int k = 0; k < probabilities_to_use.number_elements (); k++)
      {
        probability_type_t p = probabilities_to_use[k];

        BSquared & bs = bsquared[ndx];

        bs.set_weight_function (&wi);

        bs.set_number_neighbours (neighbours_to_use[j]);

        bs.set_special_category (special_class);
        bs.set_special_category_threshold (p);

        bs.compute_estimate_type ();

        ndx++;
      }
    }

    for (int j = 0; j < min_neighbours_required.number_elements (); j++)
    {
      for (int k = 0; k < distances_to_use.number_elements (); k++)
      {
        for (int l = 0; l < probabilities_to_use.number_elements (); l++)
        {
          probability_type_t p = probabilities_to_use[l];

          BSquared & bs = bsquared[ndx];

          bs.set_weight_function (&wi);
  
          bs.set_radius (distances_to_use[k]);
          bs.set_number_neighbours (min_neighbours_required[j]);

          bs.set_special_category (special_class);
          bs.set_special_category_threshold (p);

          bs.compute_estimate_type ();
  
          ndx++;
        }
      }
    }
  }

  assert (ndx == number_estimates);

  if (verbose > 1)
  {
    for (int i = 0; i < number_estimates; i++)
    {
      const BSquared & bs = bsquared[i];

      cerr << "number " << i << " is " << bs.estimate_type () << '\n';
    }
  }

  assert (ndx == number_estimates);

  return 1;
}*/

/*
  If there is no error balancing going on
*/

/*static int
initialise_bsquared_array (int number_weight_functions,
                           const NN_Weight_Function * weight_function,
                           const resizable_array<int> & neighbours_to_use,
                           const resizable_array<int> & min_neighbours_required,
                           const resizable_array<float> & distances_to_use,
                           int & number_estimates,
                           BSquared * & bsquared)
{
  int number_neighbour_types = neighbours_to_use.number_elements () + min_neighbours_required.number_elements () * distances_to_use.number_elements ();
  number_estimates = number_neighbour_types * number_weight_functions;

  if (verbose)
    cerr << "A total of " << number_estimates << " estimates will be made\n";

  bsquared = new BSquared[number_estimates];

  int ndx = 0;
  for (int i = 0; i < number_weight_functions; i++)
  {
    const NN_Weight_Function & wi = weight_function[i];

    for (int j = 0; j < neighbours_to_use.number_elements (); j++)
    {
      BSquared & bs = bsquared[ndx];

      bs.set_number_neighbours (neighbours_to_use[j]);
      bs.set_weight_function (&wi);

      bs.compute_estimate_type ();

      ndx++;
    }

    for (int j = 0; j < min_neighbours_required.number_elements (); j++)
    {
      for (int k = 0; k < distances_to_use.number_elements (); k++)
      {
        BSquared & bs = bsquared[ndx];

        bs.set_weight_function (&wi);

        bs.set_radius (distances_to_use[k]);
        bs.set_number_neighbours (min_neighbours_required[j]);

        bs.compute_estimate_type ();

        ndx++;
      }
    }
  }

  assert (ndx == number_estimates);

  if (verbose > 1)
  {
    for (int i = 0; i < number_estimates; i++)
    {
      const BSquared & bs = bsquared[i];

      cerr << "number " << i << " is " << bs.estimate_type () << '\n';
    }
  }

  assert (ndx == number_estimates);

  return 1;
}*/

/*static int
initialise_neighbour_conditions (Command_Line & cl,
                                 resizable_array<int> & neighbours_to_use,
                                 resizable_array<float> & distances_to_use,
                                 resizable_array<int> & min_neighbours_required)
{
  if (cl.option_present ('n'))
  {
    if (! parse_dash_n (cl, 'n', neighbours_to_use, verbose))
    {
      cerr << "Cannot initialise number of neighbours to use (-n option)\n";
      usage (11);
    }
  }

  if (cl.option_present ('t'))
  {
    double max_distance = 0.0;     // no longer using that feature

    if (! parse_dash_t (cl, 't', distances_to_use, min_neighbours_required, max_distance, verbose))
    {
      cerr << "Cannot parse radius options (-t)\n";
      usage (12);
    }

    if (0 == min_neighbours_required.number_elements ())
      min_neighbours_required.add (1);
  }

  if (verbose > 1)
  {
    cerr << neighbours_to_use.number_elements () << " neighbours to use\n";
    cerr << distances_to_use.number_elements () << " distances to use\n";
    cerr << min_neighbours_required.number_elements () << " min neighbours required\n";
  }

  return 1;
}*/

/*
  We can save a lot by having all the conditions with common weight functions
  grouped together
*/

#ifdef NOT_USED_j
static int
bsquared_weight_function_comparitor (const void * bs1, const void * bs2)
{
  const BSquared * b1 = (const BSquared *) (bs1);
  const BSquared * b2 = (const BSquared *) (bs2);

  const NN_Weight_Function & b1wf = b1->weight_function ();
  const NN_Weight_Function & b2wf = b2->weight_function ();

  const IWString & b1type = b1wf.weight_function_type ();
  const IWString & b2type = b2wf.weight_function_type ();

  return b1type.strcmp (b2type);
}
#endif

template <typename T>
int
do_outlier_stuff (T & nnr,
                  Command_Line & cl,
                  char flag)
{
  int i = 0;
  const_IWSubstring o;
  IWString fname;
  while (cl.value (flag, o, i++))
  {
    if (o.starts_with ("nbr="))
    {
      o.remove_leading_chars (4);
      if (! o.numeric_value (outlier_neighbours_to_print) || outlier_neighbours_to_print < 0)
      {
        cerr << "Invalid number of outlier neighbours to write '" << o << "'\n";
        usage (5);
      }

      if (verbose)
        cerr << "Will write " << outlier_neighbours_to_print << " neighbours with each outlier record\n";
    }
    else if (0 == fname.length ())
    {
      fname = o;
    }
    else
    {
      cerr << "Unrecognised -O qualifier '" << o << "'\n";
      usage (5);
    }
  }

  if (0 == fname.length ())
  {
    cerr << "No file specified for -O\n";
    usage (3);
  }

  return nnr.identify_poorly_predicted_items (fname.null_terminated_chars ());
}
static int
nn_training_continuous (Command_Line & cl,
                        int number_weight_functions,
                        NN_Weight_Function * weight_function)
{
  BSquared * bsquared;
  int number_estimates;

  if (! read_models_to_be_built (cl, number_estimates, bsquared))
  {
    cerr << "Cannot initialise neighbour conditions from -N option\n";
    return 0;
  }

// Sort the conditions so weight functions are grouped together

//#define echo_sorted_bsquared
#ifdef echo_sorted_bsquared
  cerr << "Before sorting\n";
  for (int i = 0; i < number_estimates; i++)
  {
    cerr << "Condition " << bsquared[i].estimate_type () << '\n';
  }
#endif

//qsort (bsquared, number_estimates, sizeof (BSquared), bsquared_weight_function_comparitor);

#ifdef echo_sorted_bsquared
  cerr << "After sorting\n";
  for (int i = 0; i < number_estimates; i++)
  {
    cerr << "Condition " << bsquared[i].estimate_type () << '\n';
  }
#endif

  Cross_Validation_Conditions cvc;

  if (! cvc.construct_from_command_line (cl, 'V'))
  {
    cerr << "Cannot initialise cross validation conditions\n";
    return 0;
  }

  NN_Results_Activity nnr;

  if (! initialisation_nnr (cl, nnr, bsquared, number_estimates, cvc.percent_training_set ()))
  {
    cerr << "Cannot initialise near neighbours pool\n";
    return 0;
  }

  nnr.initialise_temporary_arrays ();

//nnr.debug_print (10, cerr);

  do_leave_one_out (nnr, number_estimates, bsquared);

  if (cl.option_present ('V'))
  {
    do_cross_validations (nnr, number_estimates, bsquared, cvc);
  }

  delete [] bsquared;

  if (cl.option_present ('O'))
    do_outlier_stuff (nnr, cl, 'O');

  return 1;
}

/*
  The error balancing specification must look like

  class,min,max,delta
*/

static int
initialise_error_balancing (const_IWSubstring c,             // passed by value
                            category_type_t & special_class,
                            resizable_array<probability_type_t> & probabilities_to_use)
{
  if (4 != c.nwords (','))
    return 0;

  const_IWSubstring token;
  int i = 0;

  c.nextword (token, i, ',');

  if (! token.numeric_value (special_class) || special_class < 0 || special_class > highest_category_number)
  {
    cerr << "INvalid special class '" << token << "'\n";
    return 0;
  }

  c.nextword (token, i, ',');

  probability_type_t pmin;
  if (! token.numeric_value (pmin) || pmin <= 0.0 || pmin >= 1.0)
  {
    cerr << "Invalid minimum probability specification '" << token << "'\n";
    return 0;
  }

  c.nextword (token, i, ',');

  probability_type_t pmax;
  if (! token.numeric_value (pmax) || pmax <= pmin || pmax >= 1.0)
  {
    cerr << "Invalid max probability specification '" << token << "'\n";
    return 0;
  }

  c.nextword (token, i, ',');

  probability_type_t pdelta;
  if (! token.numeric_value (pdelta) || pdelta <= 0.0 || pdelta >= 1.0)
  {
    cerr << "Invalid delta probability specification '" << token << "'\n";
    return 0;
  }

  for (probability_type_t p = pmin; p <= pmax; p += pdelta)
  {
    probabilities_to_use.add (p);
    if (verbose)
      cerr << " Will do error balancing for class " << special_class << " at P = " << p << '\n';
  }

  if (1 == probabilities_to_use.number_elements ())
  {
    cerr << "Only one probability specified, impossible\n";
    return 0;
  }

  return 1;
}

static void
display_standard_classification_options (ostream & os)
{
  os << " -C def             default classification model\n";
  os << " -C ggb             model is a good, better, best classification model\n";
  os << " -C mac=<cls>       error rate is just the rate for class<cls>\n";
  os << " -C <cat>=<w>       assign category <cat> weight <w>\n";

  exit (0);
}

static int
initialise_classification_conditions (Command_Line & cl,
                                      char flag,
                                      int * weight_specified)
{
  assert (highest_category_number > 0);

  int weights_specified = 0;    // how many classes have a weight specified
  int default_specified = 0;
  int gbb_specified = 0;

  const_IWSubstring c;
  int i = 0;
  while (cl.value (flag, c, i++))
  {
    if ("def" == c)
    {
      default_specified = 1;
      continue;
    }

    if ("gbb" == c)
    {
      gbb_specified = 1;
      continue;
    }

    if (c.starts_with ("mac="))
    {
      c.remove_leading_chars (4);
      if (! c.numeric_value (overall_error_rate_is_error_rate_for_particular_class) || overall_error_rate_is_error_rate_for_particular_class < 0)
      {
        cerr << "The maximum accuracy class specification, must be a valid class number\n";
        usage (4);
      }
      if (verbose)
        cerr << "Will try to maximise accuracy for class " << overall_error_rate_is_error_rate_for_particular_class << " only\n";
      continue;
    }

    if ("help" == c)
    {
      display_standard_classification_options (cerr);
    }

//  We want to allow individual cat=w, and cat=w,cat=w,cat=w

    int j = 0;
    const_IWSubstring token;
    while (c.nextword (token, j, ','))
    {
      const_IWSubstring cat, w;
      if (! token.split (cat, '=', w))
      {
        cerr << "Invalid category weight specification '" << token << "'\n";
        return 0;
      }

      int category;
      if (! cat.numeric_value (category) || category < 0 || category > highest_category_number)
      {
        cerr << "Invalid category number '" << cat << "'\n";
        return 0;
      }

      weight_t weight;
      if (! w.numeric_value (weight) || weight < 0.0)
      {
        cerr << "Invalid category weight '" << w << "'\n";
        return 0; 
      }

      category_weight[category] = weight;
      weight_specified[category] = 1;
      weights_specified++;
    }
  }

  if (gbb_specified && default_specified)
  {
    cerr << "Cannot do both 'gbb' and 'def'\n";
    return 0;
  }

  if ((gbb_specified || default_specified) && weights_specified)
  {
    cerr << "Can choose either default category weights or default, not both\n";
    return 0;
  }

  if (weights_specified)
  {
    weight_t sumw = 0.0;
    weights_specified = 0;

    for (int i = 0; i <= highest_category_number; i++)
    {
      if (weight_specified[i])
      {
        weights_specified++;
        sumw += category_weight[i];
      }
      else       // unspecified weights are 1.0
        sumw += 1.0;
    }

    for (int i = 0; i <= highest_category_number; i++)
    {
      category_weight[i] = category_weight[i] / sumw;
    }
  }

  if (overall_error_rate_is_error_rate_for_particular_class >= 0 &&
      overall_error_rate_is_error_rate_for_particular_class > highest_category_number)
  {
    cerr << "The most accurately predicted class number must be between 0 and " << highest_category_number << ", " << overall_error_rate_is_error_rate_for_particular_class << " is invalid\n";
    return 0;
  }

  if (verbose > 1)
  {
    cerr << "Category weights\n";
    for (int i = 0; i <= highest_category_number; i++)
    {
      cerr << "Category " << i << " " << number_in_category[i] << " items, weight " << category_weight[i] << '\n';
    }
  }

  return 1;
}

static int
initialise_classification_conditions (Command_Line & cl,
                                      char flag)
{
  int * tmp = new_int (highest_category_number + 1);

  std::unique_ptr<int> free_tmp (tmp);

  int rc = initialise_classification_conditions (cl, flag, tmp);

  return rc;
}

static int
nn_training_category (Command_Line & cl,
                      int number_weight_functions,
                      NN_Weight_Function * weight_function)
{
  resizable_array<probability_type_t> probabilities_to_use;
  category_type_t special_class = -1;

  if (cl.option_present ('L'))
  {
    const_IWSubstring l = cl.string_value ('L');
    if (! initialise_error_balancing (l, special_class, probabilities_to_use))
    {
      cerr << "Cannot process error balancing conditions '" << l << "' (-L option)\n";
      return 0;
    }
  }

  BSquared * bsquared;
  int number_estimates;

  if (cl.option_present ('N'))
  {
    if (! read_models_to_be_built (cl, number_estimates, bsquared))
    {
      cerr << "Cannot initialise neighbour conditions from -N option\n";
      return 0;
    }
  }

  Cross_Validation_Conditions cvc;

  if (! cvc.construct_from_command_line (cl, 'V'))
  {
    cerr << "Cannot initialise cross validation conditions\n";
    return 0;
  }

  NN_Results_Category nnr;

  if (! initialisation_nnr (cl, nnr, bsquared, number_estimates, cvc.percent_training_set ()))
  {
    cerr << "Cannot initialise near neighbours pool\n";
    return 0;
  }

  if (! nnr.establish_categories ())
  {
    cerr << "Cannot establish categories from input\n";
    return 0;
  }

  if (! initialise_classification_conditions (cl, 'C'))
  {
    cerr << "Cannot parse the -C option\n";
    return 0;
  }

  do_leave_one_out (nnr, number_estimates, bsquared);

  if (cl.option_present ('V'))
  {
    nnr.assign_number_from_each_category (cvc.percent_training_set ());

    do_cross_validations (nnr, number_estimates, bsquared, cvc);
  }

  if (cl.option_present ('O'))
    do_outlier_stuff (nnr, cl, 'O');

  if (file_name_for_test_set_predictions.length () > 0)
  {
    nnr.set_training_set_membership (1);     // everything in the training set
    nnr.unset_predicted_values ();
    nnr.compute_leave_one_out_values (bsquared[0], 0);
    nnr.compute_model_quality (1, bsquared, 0);
    nnr.report_misclassifications_by_class (cerr);
  }

  delete [] category_weight;
  delete [] bsquared;

  return 1;
}

static int
display_dash_j_options (ostream & os)
{
  os << " -J trgfp=<file>   specify training set gfp file (for written to -R file)\n";

  exit (0);

  return 0;          // keep the SGI compiler quiet
}

static int
nn_predictions (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vE:jdW:M:V:T:O:g:u:p:C:s:L:N:zP:R:D:J:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "unrecognised options encountered\n";
    usage (4);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('C'))
  {
    classification_model = 1;

    bsquared_or_error_rate = "ERROR";

    const_IWSubstring c = cl.string_value ('C');

    if (verbose)
      cerr << "Will treat as a classification model\n";
  }
  else
  {
    bsquared_or_error_rate = "BSquared";
  }

  if (cl.option_present ('z'))
  {
    strip_leading_zeros_from_identifiers = 1;

    if (verbose)
      cerr << "Leading zero's will be removed from identifiers when doing experimental data match-ups\n";
  }

  if (cl.option_present ('g'))
  {
    if (! cl.value ('g', min_meaningful_bsquared_difference) || min_meaningful_bsquared_difference <= 0.0 || min_meaningful_bsquared_difference >= 1.0)
    {
      cerr << "The minimum meaningful B squared difference (-g) option must be a positive number\n";
      usage (8);
    }

    if (verbose)
      cerr << "Differences in " << bsquared_or_error_rate << " less than " << min_meaningful_bsquared_difference << " will be ignored\n";
  }

  if (cl.option_present ('u'))
  {
    if (! cl.value ('u', shortest_useful_distance) || shortest_useful_distance <= 0.0 || shortest_useful_distance >= maximum_distance_possible ())
    {
      cerr << "INvalid shortest useful distance (-u) option\n";
      usage (7);
    }

    if (verbose)
      cerr << "The best distance printed will be at least " << shortest_useful_distance << '\n';
  }

  if (cl.option_present ('M'))
  {
    missing_value = cl.string_value ('M');

    if (verbose)
      cerr << "Missing values written as '" << missing_value << "'\n";
  }

  if (cl.option_present ('D'))
  {
    if (! initialise_distances (cl, 'D', verbose))
    {
      cerr << "Cannot initialise distances (-D option)\n";
      usage (5);
    }
  }

  if (cl.option_present ('p'))
  {
    if (! cl.value ('p', prediction_types_to_print) || prediction_types_to_print < 1)
    {
      cerr << "Invalid number of prediction types to print (-p option)\n";
      usage (21);
    }

    if (verbose)
      cerr << "Will print the best " << prediction_types_to_print << " predictions\n";
  }

  if (cl.option_present ('P'))
  {
    file_name_for_test_set_predictions = cl.string_value ('P');

    if (verbose)
      cerr << "Will write best training set prediction to '" << file_name_for_test_set_predictions << "'\n";
  }

  if (! cl.option_present ('N'))
  {
    cerr << "Must specify models to be built with the -N option\n";
    usage (5);
  }

// We don't use a gfp file, but if we have that info, we can write it to the -R file

  IWString trgfp;

  if (cl.option_present ('J'))
  {
    int i = 0;
    const_IWSubstring j;
    while (cl.value ('J', j, i++))
    {
      if (j.starts_with ("trgfp="))
      {
        trgfp = j;
        trgfp.remove_leading_chars (6);
      }
      else if ("help" == j)
      {
        display_dash_j_options (cerr);
      }
      else
      {
        cerr << "Unrecognised -J qualifier '" << j << "'\n";
        display_dash_j_options (cerr);
      }
    }
  }

  if (cl.option_present ('d'))
  {
    nn_file_is_distance_matrix = 1;

    if (verbose)
      cerr << "NN file is a distance matrix\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (8);
  }

  if (cl.option_present ('R'))
  {
    IWString fname = cl.string_value ('R');

    stream_for_best_models.open (fname.null_terminated_chars (), std::ios::out);

    if (! stream_for_best_models.good ())
    {
      cerr << "Cannot open file for best models '" << fname << "'\n";
      return 8;
    }

    if (verbose)
      cerr << "Best results written to '" << fname << "'\n";

    stream_for_best_models << KNN_MODEL_FILE_VERSION_1 << '\n';

//  training set gfp file will always be in the model directory, so this
//  record appears as a comment

    stream_for_best_models << '#' << TR_GFP_DIRECTIVE << ' ';
    if (trgfp.length ())
      stream_for_best_models << trgfp << '\n';
    else
      stream_for_best_models << TR_GFP_UNKNOWN << "\n";

    stream_for_best_models << '\n';

    stream_for_best_models << KNN_MODEL_TYPE << ' ';
    if (classification_model)
      stream_for_best_models << KNN_CLASSIFICATION << '\n';
    else
      stream_for_best_models << KNN_CONTINUOUS << '\n';

    stream_for_best_models << '\n';
  }

  int number_weight_functions = cl.option_count ('W');

  NN_Weight_Function * weight_function;
  if (0 == number_weight_functions)
  {
    weight_function = new NN_Weight_Function[1];
    weight_function[0].initialise_default_weight_function ();
    number_weight_functions = 1;
  }
  else
  {
    weight_function = new NN_Weight_Function[number_weight_functions];

    for (int i = 0; i < number_weight_functions; i++)
    {
      const_IWSubstring w;
      cl.value ('W', w, i);

      if (! weight_function[i].construct_from_command_line_token (w, 'W'))
      {
        cerr << "Cannot initialise weight function specification '" << w << "'\n";
        usage (14);
      }

      if (verbose)
      {
        cerr << "Weight function " << i << " is '" << weight_function[i].weight_function_type () << "'\n";
        if (verbose > 1)
          weight_function[i].echo_weight_function (cerr);
       }
    }
  }

  if (verbose)
    cerr << "Initialised " << number_weight_functions << " weight functions\n";

// We must have a source of experimental data

  if (! cl.option_present ('E'))
  {
    cerr << "Must specify source of experimental data via the -E option\n";
    usage (4);
  }

  int rc;
  if (classification_model)
  {
    rc = nn_training_category (cl, number_weight_functions, weight_function);
  }
  else
  {
    rc = nn_training_continuous (cl, number_weight_functions, weight_function);
  }

  return 0 == rc;
}

int
main (int argc, char ** argv)
{
  int rc = nn_predictions (argc, argv);

  return rc;
}
