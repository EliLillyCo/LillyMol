#include <iostream>
#include <math.h>
#include <memory>
using std::cerr;
using std::endl;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwstring.h"
#include "iwrandom.h"
#include "accumulator.h"
#include "set_or_unset.h"
#include "iwcrex.h"

/*
  As the buffers get larger, it may be advantageous to sort the stored values
  Jun 2000, always sort the buffers. Needed for the quartile determinations
  for scaleA
*/

static int sort_stored_values = 1;

/*
  We treat binary and near binary descriptors specially.
  If a descriptor has >= continuous_threshold different values, then it is continuous
*/

static int continuous_threshold = 5;

static float binary_collapse_threshold = static_cast<float>(0.0);

/*
  Nov 2000. Just ignore columns with invalid data....
*/

static int ignore_columns_with_non_numeric_data = 0;

/*
  When reducing to a given number of columns, we have various ways we can choose
  the most desirable columns;
*/

#define HV_MOST_VALUES 1
#define HV_CV 2
#define HV_SCALEA 3

static int chosen_highest_variance_measure = HV_MOST_VALUES;

static float adjusted_median_threshold = static_cast<float>(0.01);

static int nkeep = 0;

/*
  Oct 2014. Enable buffring the whole file
*/

static int buffer_file = 0;

static resizable_array_p<IWString> all_data;

/*
  Once the buffer is full, we can store new values with their closest
  neighbour
*/

static int store_with_closest_neighbour = 0;

static IWString missing_value('.');

static char input_separator(' ');
static char output_separator(' ');

/*
  We keep track of all the different values encountered in a column
*/

template <typename T>
class value_count : public resizable_array<T>, public Accumulator<T>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using resizable_array_base<T>::_number_elements;
    using resizable_array_base<T>::_elements_allocated;
    using resizable_array_base<T>::_things;
#endif

  private:
    int * _count;

    int _missing_values;

    float _fraction_of_most_common_value;

// private functions

    int _extra_sorted(T);
    int _can_collapse_to_binary(T & v1, int & c1, T & v2, int & c2) const;
    T _scaleA() const;
    T _interpolate(int number_needed, int number_encountered, int i) const;
    T _find_quantile(int number_needed) const;

    T _adjusted_cv() const;
    T _rescaled_scaleA() const;
    T _unscaled_scaleA() const;

    T _quantile_difference(int, int) const;

    int _total_number_values_stored() const;

  public:
    value_count();
    ~value_count();

    int count(int i) const { return _count[i];};

    int missing_values() const { return _missing_values;}

    int initialise(int);
    int looks_like_fingerprint_data() const;

    int extra(const const_IWSubstring &);

    int well_distributed(float) const;

    int only_one_value_significantly_hit(float) const;
    float fraction_of_most_common_value();

    int several_values_significantly_hit(int nvalues, float p) const;
    int median_too_common(float p, int column) const;

    T rescaled_variance() const;
    T best_variance_measure();
    T chosen_variance_measure();

    int all_values_greater_than_zero() const;
};

template <typename T>
value_count<T>::value_count()
{
  _count = NULL;

  _missing_values = 0;

  _fraction_of_most_common_value = -1.0;
}

template <typename T>
value_count<T>::~value_count()
{
  delete [] _count;
  _count = NULL;
}

template <typename T>
int
value_count<T>::initialise(int nkeep)
{
  assert (NULL == _count);

  resizable_array_base<T>::resize(nkeep);

  _count = new int[nkeep];
  for (int i = 0; i < nkeep; i++)
  {
    _count[i] = 0;
  }

  return 1;
}

template <typename T>
int
value_count<T>::looks_like_fingerprint_data() const
{
  if (_number_elements > 2)
    return 0;

  if (0 == _number_elements)
    return 0;

  if (1 == _number_elements)
  {
    if (0.0 == _things[0] || 1.0 == _things[0])
      return 1;
    else
      return 0;
  }

  if (0.0 == _things[0] && 1.0 == _things[1])
    return 1;
  if (1.0 == _things[0] && 0.0 == _things[1])
    return 1;

  return 0;
}

template <typename T>
int
value_count<T>::_total_number_values_stored() const
{
  int rc = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    rc += _count[i];
  }

  return rc;
}

//#define DEBUG_EXTRA_SORTED

template <typename T>
int
value_count<T>::_extra_sorted(T v)
{
#ifdef DEBUG_EXTRA_SORTED
  cerr << "Adding " << v << " to";
  for (int i = 0; i < _number_elements; i++)
    cerr << ' ' << _things[i];
  cerr << endl;
#endif

  int left = 0;

  if (v < _things[left])
  {
    if (_number_elements >= nkeep)    // array is full
      return 1;

    resizable_array<T>::insert_at_beginning(v);
    for (int i = _number_elements - 1; i > 0; i--)
      _count[i] = _count[i - 1];

    _count[0] = 1;

    return 1;
  }
  else if (v == _things[left])
  {
    _count[left]++;
    return 1;
  }

  int right = _number_elements - 1;

  if (v > _things[right])
  {
    if (_number_elements >= nkeep)      // array already full
      return 1; 

    this->add(v);
    _count[_number_elements - 1] = 1;

    return 1;
  }
  else if (v == _things[right])
  {
    _count[right]++;
    return 1;
  }

// At this stage, it is somewhere in the middle

  int middle;

  while ((middle = (right + left) / 2) != left)
  {
    if (v < _things[middle])
      right = middle;
    else if (v > _things[middle])
      left = middle;
    else
    {
      _count[middle]++;
      return 1;
    }
  }

#ifdef CHECK_BINARY_SEARCH
  for (int i = 0; i < _number_elements; i++)
  {
    if (v == _things[i])
    {
      cerr << "Binary search failed\n";
    }
  }
#endif

// Looks like this is a new value

  if (_number_elements >= nkeep)    // no room for new value
  {
    if (0 == store_with_closest_neighbour)
      return 1;
    if (v - _things[middle] < _things[middle + 1] - v)
      _count[middle]++;
    else
      _count[middle + 1]++;

    return 1;
  }

  resizable_array<T>::insert_after(middle, v);
  for (int i = _number_elements - 1; i > middle + 1; i--)
    _count[i] = _count[i - 1];

  _count[middle + 1] = 1;

#ifdef DEBUG_EXTRA_SORTED
  cerr << "After insertion:";
  for (int i = 0; i < _number_elements; i++)
    cerr << ' ' << _things[i];
  cerr << endl;
#endif

  return 1;
}

template <typename T>
int
value_count<T>::extra(const const_IWSubstring & token)
{
  if (missing_value == token)
  {
    _missing_values++;
    return 1;
  }

  T v;
  if (! token.numeric_value(v))
    return 0;

  Accumulator<T>::extra(v);

  if (0 == _number_elements)
  {
    resizable_array<T>::add(v);
    _count[0] = 1;

    return 1;
  }

  if (sort_stored_values)
    return _extra_sorted(v);

  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i] == v)
    {
      _count[i]++;
      return 1;
    }
  }

  if (_number_elements < nkeep)
  {
    resizable_array<T>::add(v);
    _count[_number_elements - 1] = 1;
  }

  return 1;
}

/*
  Well, what is well distributed?

  We compute the most hit and the next most hit values.
  If the ratio is less than P, then this is not well distributed
*/

template <typename T>
int
value_count<T>::well_distributed(float p) const
{
  if (1 == _number_elements)
    return 0;

  int max_pop1 = 0;
  int max_pop2 = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (_count[i] > max_pop1)
    {
      max_pop2 = max_pop1;
      max_pop1 = _count[i];
    }
    else if (_count[i] > max_pop2)
      max_pop2 = _count[i];
  }

  if (static_cast<float>(max_pop2) / static_cast<float>(max_pop1) < p)
    return 0;
  else
    return 1;
}

/*
  If the most common value comprises P or more of the samples,
  that means not well distributed
*/

template <typename T>
int
value_count<T>::only_one_value_significantly_hit(float p) const
{
  int most_common = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (_count[i] > most_common)
      most_common = _count[i];
  }

  if (static_cast<float>(most_common) / static_cast<float>(Accumulator<T>::n()) > p)
    return 1;
  else
    return 0;
}

template <typename T>
float
value_count<T>::fraction_of_most_common_value()
{
  if (_fraction_of_most_common_value > 0.0)
    return _fraction_of_most_common_value;

  int most_common = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (_count[i] > most_common)
      most_common = _count[i];
  }

  _fraction_of_most_common_value = static_cast<float>(most_common) / static_cast<float>(Accumulator<T>::n());

  return _fraction_of_most_common_value;
}

template <typename T>
int
value_count<T>::median_too_common(float p, int column) const
{
  if (! sort_stored_values)
  {
    cerr << "Median only valid with sorted values\n";
    abort();
  }

  int sum = 0;    // the number counted so far
  int target = (Accumulator<T>::n() + 1) / 2;    // the threshold for median

  int median = -1;
  int highest_count = 0;
  int value_with_highest_count = -1;

  for (int i = 0; i < _number_elements; i++)
  {
    if (_count[i] > highest_count)
    {
      highest_count = _count[i];
      value_with_highest_count = i;
    }

    sum += _count[i];
    if (median < 0 && sum >= target)    // we have the median
    {
      median = i;
    }
  }

  if (median != value_with_highest_count)
  {
    cerr << "Warning, column " << column << " median " << median << " different from most hit " << value_with_highest_count << endl;
    cerr << "Median has " << _count[median] << " hits, most common is " << _count[value_with_highest_count] << endl;
  }

  if (_count[median] >= p)
    return 1;
  else
    return 0;
}

/*
  Another test is to require that NVALUES or more values are
  hit more than P of the time.
  For example, we might require that at least two values are
  hit 20% of the time or more.
*/

template <typename T>
int
value_count<T>::several_values_significantly_hit(int nvalues, float p) const
{
  assert (p >= 0.0 && p < 1.0);

  int cutoff = int(float(Accumulator<T>::n()) * p);
//cerr << "Cutoff for " << n() << " values with " << p << " is " << cutoff << endl;

  int number_above_cutoff = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (_count[i] >= cutoff)
    {
      number_above_cutoff++;
      if (number_above_cutoff >= nvalues)
        return 1;
    }
  }

  return 0;
}

/*
  Given C1 instances of V1 and C2 instances of V2, what is the binary variance
*/

static double
binary_descriptor_variance(int c1, int c2)
{
  double p = static_cast<double>(c1) / static_cast<double>(c1 + c2);

  return sqrt(p * (1.0 - p));
}

//#define DEBUG_BEST_VARIANCE_MEASURE 1

template <typename T>
T
value_count<T>::best_variance_measure()
{
#ifdef DEBUG_BEST_VARIANCE_MEASURE
  if (_number_elements >= continuous_threshold)
    cerr << "best_variance_measure:: > continuous threshold, scaleA\n";
#endif

  if (_number_elements >= continuous_threshold)
    return _scaleA();

  if (_number_elements < 2)
    return static_cast<T>(0.0);

#ifdef DEBUG_BEST_VARIANCE_MEASURE
  if (2 == _number_elements)
    cerr << "best_variance_measure:: just two different values\n";
#endif

  if (2 == _number_elements)
    return static_cast<T>(binary_descriptor_variance(_count[0], _count[1]));

  T v1, v2;
  int c1, c2;

#ifdef DEBUG_BEST_VARIANCE_MEASURE
  if (binary_collapse_threshold > 0.0 && _can_collapse_to_binary(v1, c1, v2, c2))
    cerr << "best_variance_measure:: collapsed to binary\n";
#endif

  if (binary_collapse_threshold > 0.0 && _can_collapse_to_binary(v1, c1, v2, c2))
    return static_cast<T> (binary_descriptor_variance(c1, c2));

#ifdef DEBUG_BEST_VARIANCE_MEASURE
  cerr << "best_variance_measure::default case, using scaleA\n";
#endif

  return _scaleA();
}

//#define DEBUG_CHOSEN_VARIANCE_MEASURE

template <typename T>
T
value_count<T>::chosen_variance_measure()
{
  if (HV_MOST_VALUES == chosen_highest_variance_measure)
    return static_cast<T>(_number_elements);

  if (HV_CV == chosen_highest_variance_measure)
    return _adjusted_cv();

  if (HV_SCALEA == chosen_highest_variance_measure)
    return _rescaled_scaleA();

  abort();

  return 0.0;
}

template <typename T>
T
value_count<T>::_adjusted_cv() const
{
  T sa = _scaleA();

// Now find the median.
// Interesting issue with number needed. Should we use Accumulator::n () or count the
// number of items we sampled and stored. For now we use stored values...

  int number_present = _total_number_values_stored();

  int number_needed = number_present / 2;

  T median = _find_quantile(number_needed);

  T rc;

#ifdef DEBUG_CHOSEN_VARIANCE_MEASURE
  cerr << "Adjusted CV: scaleA " << sa  << " median " << median;
  if (median >= - adjusted_median_threshold && median <= adjusted_median_threshold)
    cerr << " median below threshold";
  cerr << endl;
#endif

  if (median <= adjusted_median_threshold && median >= -adjusted_median_threshold)    // denominator would be set to 1.0
    rc = sa;
  else
    rc= sa / median;

#ifdef DEBUG_CHOSEN_VARIANCE_MEASURE
  cerr << "  adjusted cv = " << rc << endl;
#endif

  return rc;
}

template <typename T>
T
value_count<T>::_rescaled_scaleA() const
{
  int number_stored = _total_number_values_stored();

  int number_needed = number_stored * 2 / 100;

  T q2 = _find_quantile(number_needed);

  number_needed = number_stored * 98 / 100;

  T q98 = _find_quantile(number_needed);

  T sa = _scaleA();

  T denominator = q98 - q2;

#ifdef DEBUG_CHOSEN_VARIANCE_MEASURE
  cerr << "Rescaled scaleA: q2 = " << q2 << " q98 = " << q98 << " sa = " << sa << endl;
#endif

  if (0.0 == denominator)
    return sa;

  return sa / denominator;
}

template <typename T>
T
value_count<T>::rescaled_variance() const
{
  if (_number_elements < 2)
    return static_cast<T>(0.0);

  if (2 == _number_elements)
    return Accumulator<T>::variance();

  T xmin = Accumulator<T>::minval();
  T range = _quantile_difference(98, 2);

  cerr << "rescaled_variance: between " << Accumulator<T>::minval() << " and " <<  Accumulator<T>::maxval() << " range " << (Accumulator<T>::maxval() - Accumulator<T>::minval()) << " q98-q2 " << range << endl;

  Accumulator<T> tmp;

  for (int i = 0; i < _number_elements; i++)
  {
    T x = _things[i];

    x = (x - xmin) / range;

    tmp.extra(x, _count[i]);
  }

  return tmp.variance();
}

//#define DEBUG_SCALEA

/*
  Jun 2000. Designed to limit
  the influence of outliers. The formula is

  min (sd(dat), quantile(dat, 0.75) - quantile(dat, 0.25))
*/

template <typename T>
T
value_count<T>::_scaleA() const
{
  if (2 == _number_elements)     // must be binary data
    return Accumulator<T>::variance();

  int number_stored = _total_number_values_stored();

// We need to find the 25, and 75'th percentiles

  int number_needed = number_stored / 4;
  if (0 == number_needed)
    number_needed = 1;

  T q25 = _find_quantile(number_needed);

  number_needed = number_stored * 3 / 4;

  T q75 = _find_quantile(number_needed);

  assert (q25 <= q75);

  T sd = rescaled_variance();

  T denominator = _quantile_difference(98, 2);

#ifdef DEBUG_SCALEA
  cerr << "q25 = " << q25 << " q75 = " << q75 << ", denominator(q98 - q2) = " << denominator << endl;
#endif

  if (0.0 == denominator)
    return sd;

  T sa = (q75 - q25) / denominator / 1.38;

#ifdef DEBUG_SCALEA
  cerr << "sa = " << sa << " sd = " << sd << endl;
#endif

  if (0.0 == sa)
    return sd;

  if (sa < sd)
    return sa;
  else
    return sd;
}

/*
  No scaling of anything
*/

template <typename T>
T
value_count<T>::_unscaled_scaleA() const
{
  assert(_number_elements > 2);

  int number_stored = _total_number_values_stored();

// We need to find the 25, and 75'th percentiles

  int number_needed = number_stored / 4;
  if (0 == number_needed)
    number_needed = 1;

  T q25 = _find_quantile(number_needed);

  number_needed = number_stored * 3 / 4;

  T q75 = _find_quantile(number_needed);

#ifdef DEBUG_SCALEA
  cerr << "q25 = " << q25 << " q75 = " << q75 << endl;
#endif

  assert(q25 <= q75);

  T sd = Accumulator<T>::variance();

  T sa = (q75 - q25) / 1.38;

#ifdef DEBUG_SCALEA
  cerr << "sa = " << sa << " sd = " << sd << endl;
#endif

  if (0.0 == sa)
    return sd;

  if (sa < sd)
    return sa;
  else
    return sd;
}

//#define DEBUG_INTERPOLATE

/*
  In doing the determination of quarties, we are looking for NUMBER_NEEDED, items.
  We have scanned through I and have encountered NUMBER_ENCOUNTERED samples. Return
  our best estimate of the real quartile position
*/

template <typename T>
T
value_count<T>::_interpolate(int number_needed, int number_encountered, int i) const
{
  assert (i > 0);
  assert (number_encountered > number_needed);

#ifdef DEBUG_INTERPOLATE
  cerr << "Interpolating between " << _things[i - 1] << " and " << _things[i] << " encountered " << number_encountered << " needed " << number_needed << endl;
#endif

// what was the count at (i - 1)

  int number_encountered_previously = number_encountered - _count[i];

  assert (number_encountered_previously > 0);

  T ratio = static_cast<T>(number_needed - number_encountered_previously) / static_cast<T>(_count[i]);

#ifdef DEBUG_INTERPOLATE
  cerr << "Number needed " << number_needed << " interpolation ratio " << ratio << endl;
#endif

  T rc = _things[i - 1] + (_things[i] - _things[i - 1]) * ratio;

#ifdef DEBUG_INTERPOLATE
  cerr << "From " << _things[i - 1] << " to " << _things[i] << " value " << rc << endl;
#endif

  return rc;
}

template <typename T>
T
value_count<T>::_find_quantile(int number_needed) const
{
  int number_encountered = 0;

  for (int i = 0; i < _number_elements; i++)
  {
    number_encountered += _count[i];

    if (number_encountered < number_needed)     // haven't found enough yet
      continue;

    if (number_encountered == number_needed)    // got it exactly
      return _things[i];

    if (0 == i)         // cannot interpolate
      return _things[i];

    return _interpolate(number_needed, number_encountered, i);    // found too many, need to interpolate back
  }

  assert (NULL == "Should not come here");

  return 0.0;
}

template <typename T>
T
value_count<T>::_quantile_difference(int i1, int i2) const
{
  assert (i1 >= i2);

  int number_stored = _total_number_values_stored();

  int n1 = i1 * number_stored / 100;
  int n2 = i2 * number_stored / 100;

  T q1 = _find_quantile(n1);
  T q2 = _find_quantile(n2);

//cerr << "_quantile_difference: q1 = " << q1 << " q2 = " << q2 << " returning " << (q1 - q2) << endl;

  return q1 - q2;
}

/*
  Can a distribution be collapsed to binary

  c1, c2, c3 and c4 are the counts of the most frequent values - in order
*/

template <typename T>
int
value_count<T>::_can_collapse_to_binary(T & v1, int & c1, T & v2, int & c2) const
{
  assert (_number_elements > 2);     // should not be called in this case

  if (_number_elements >= continuous_threshold)    // cannot be collapsed
    return 0;

  T v3;
  int c3;
  T v4;
  int c4;

  c1 = _count[0];
  v1 = _things[0];

  c2 = c3 = c4 = 0;

  for (int i = 1; i < _number_elements; i++)
  {
    if (_count[i] > c1)
    {
      c4 = c3;
      v4 = v3;
      c3 = c2;
      v3 = v2;
      c2 = c1;
      v2 = v1;

      c1 = _count[i];
      v1 = _things[i];
    }
    else if (_count[i] > c2)
    {
      c4 = c3;
      v4 = v3;
      c3 = c2;
      v3 = v2;

      c2 = _count[i];
      v2 = _things[i];
    }
    else if (_count[i] > c3)
    {
      c4 = c3;
      v4 = v3;

      c3 = _count[i];
      v3 = _things[i];
    }
    else if (_count[i] > c4)
    {
      c4 = _count[i];
      v4 = _things[i];
    }
  }

  assert (c1 >= c2 && c2 >= c3 && c3 >= c4);

  float binary_density = static_cast<float>(c1 + c2) / static_cast<float>(c1 + c2 + c3 + c4);

#ifdef DEBUG_BEST_VARIANCE_MEASURE
  cerr << "Binary density " << binary_density << endl;
#endif

  if (binary_density < binary_collapse_threshold)    // cannot be collapsed
    return 0;

// We must now combine the one or two minor values with one of the two major values

  return 1;
}

template <typename T>
int
value_count<T>::all_values_greater_than_zero() const
{
  assert (_number_elements > 0);

// the _things array is sorted in ascending order

  if (_things[0] > 0 && Accumulator<T>::minval() > static_cast<T>(0.0))
    return 1;

  return 0;
}

class value_count_float : public value_count<float>
{
  private:
  public:
    value_count_float();
    ~value_count_float();
};

value_count_float::value_count_float()
{
}

value_count_float::~value_count_float()
{
}

#ifdef __GNUG__
template class resizable_array_p<value_count_float>;
template class resizable_array_base<value_count_float *>;

template class value_count<float>;
#endif

const char * prog_name = NULL;

static int verbose = 0;

static int print_statistics = 0;

static value_count_float * counters = NULL;

static int translate_tabs = 0;

/*
  Dealing with the ignore column is a little difficult, because
  the values come from the command line, before we know how many
  records there are in the file. So, we collect them in a
  resizable_array<int> and once the file size is known, we fill
  the array
*/

static resizable_array<int> ignore_columns_from_cl;

static int * ignore_column = NULL;

/*
  Alternatively, we can process only columns which match a regexp (descriptor files only)
*/

static IW_Regular_Expression process_column_regexp;

static extending_resizable_array<int> missing_value_counter;

static int is_descriptor_file = 0;

static int nskip = 0;

static float probability = 2.0;    // probability of choosing a record

static int print_distributions = 0;

static int columns_in_input = 0;

static int min_different_values = 2;

static int records_to_process = 0;

/*
  If all the columns are in the same units (like a Comfa grid) we
  can eliminate columns below a given variance
*/

static float min_variance = -1.0;

/*
  Alternatively, we can keep the N columns with the highest variance
*/

static int highest_variance_nkeep = 0;

/*
  Although probably more useful is to keep the columns with the lowest
  value of the item most commonly hit
*/

static int lowest_x_nkeep = 0;

/*
  By default, we cease processing if we find that the file needs no changes.
  We can optionally write the file
*/

static int write_file_even_if_no_changes = 0;

/*
  most_common_max_pop and max_fraction_for_median are very similar.
  One works on the most common value, the other on the median
*/

static float most_common_max_pop = 0.0;

static float max_fraction_for_median = 0.0;

/*
  Variables for the -f and -b options
*/

static int min_values_above_population_threshold = 0;
static float population_threshold = static_cast<float>(0.0);

/*
  Given that we know the average and variance of each column, we can
  optionally normalise what we write
*/

static int normalise_output = 0;

/*
  Aug 97. The -q option allows rejecting a column of the fraction
  of missing values is > a threshold
*/

static float maximum_missing_value_fraction = static_cast<float>(-1.0);

/*
  Aug 2000.
  Chris Waller was using Coner Keys, and wanted to eliminate columns for which
  every member was greater than 0 - that is, the key was set in every row
*/

static int remove_if_all_greater_than_zero = 0;

/*
  Sept 2002.
  We want to get rid of descriptors where everything is "small"
*/

static float remove_if_all_values_less_than = static_cast<float>(0.0);

#include "iwstring_data_source.h"

static IWString * column_title = NULL;

/*
  Generic routine for writing column info.
  Note that the invoker considers columns to start with column 1, so
  we increment the value written by 1 to accommodate that.
*/

static void
write_column(int c, std::ostream & os)
{
  os << "column " << (c + 1);
  if (is_descriptor_file)
    os << " '" << column_title[c] << "' ";

  return;
}

static int
pair_comparitor(const void * pv1, const void * pv2)
{
  const std::pair<int, float> * p1 = (std::pair<int, float> *)pv1;
  const std::pair<int, float> * p2 = (std::pair<int, float> *)pv2;

  if (p1->second > p2->second)
    return -1;
  else if (p1->second < p2->second)
    return 1;
  else
    return 0;
}

static int
do_highest_variance_nkeep(value_count_float * counters,
                           int * keep,
                           int & columns_suppressed,
                           std::pair<int, float> * cv)
{
  int nprocess = 0;

  for (int i = 0; i < columns_in_input; i++)
  {
    if (0 == keep[i])
      continue;

    if (ignore_column && ignore_column[i])
      continue;

    std::pair<int, float> & c = cv[nprocess];
    nprocess++;

    value_count_float & v = counters[i];

    c.first = i;

#ifdef DEBUG_BEST_VARIANCE_MEASURE
    cerr << "Getting chosen variance for column " << (i + 1) << " '" << column_title[i] << "', " << v.number_elements() << " values between " << v.minval() << " and " << v.maxval() << endl;
#endif

    if (v.n() <= 1)
      c.second = 0.0;
    else if (highest_variance_nkeep)
      c.second = counters[i].chosen_variance_measure();
    else if (lowest_x_nkeep)
      c.second = 1.0 - counters[i].fraction_of_most_common_value();   // we want to keep the columns with the lowest fraction

#ifdef DEBUG_BEST_VARIANCE_MEASURE
    cerr << "Column " << (i + 1) << " '" << column_title[i] << "' chosen variance " << c.second << endl;
#endif

  }

  int nkeep;     // the number of columns requested to be kept

  if (highest_variance_nkeep)
    nkeep = highest_variance_nkeep;
  else if (lowest_x_nkeep)
    nkeep = lowest_x_nkeep;
  else
    abort();

  if (nprocess <= nkeep)    // already fewer than nkeep columns still alive
  {
    if (verbose)
      cerr << "Requested " << nkeep << " most variable columns. Only " << nprocess << " still active\n";
    return 1;
  }

  qsort(cv, nprocess, sizeof(std::pair<int, float>), pair_comparitor);

//#define DEBUG_QSORT
#ifdef DEBUG_QSORT
  cerr << "After sorting " << nprocess << " columns\n";
  for (int i = 0; i < nprocess; i++)
  {
    const std::pair<int, float> & v = cv[i];
    cerr << i << ' ';
    write_column(v.first, cerr);
    cerr << " value " << v.second << endl;
  }
#endif

// Anything beyond nkeep will be suppressed

  for (int i = nkeep; i < nprocess; i++)
  {
    const std::pair<int, float> & v = cv[i];

    keep[v.first] = 0;

    columns_suppressed++;
  }

  return 1;
}

static int
do_highest_variance_nkeep(value_count_float * counters,
                           int * keep,
                           int & columns_suppressed)
{
  std::pair<int, float> * cv = new std::pair<int, float> [columns_in_input]; std::unique_ptr<std::pair<int, float>[]> free_cv(cv);

  return do_highest_variance_nkeep(counters, keep, columns_suppressed, cv);
}

static int
establish_column_titles(const const_IWSubstring & buffer)
{
  assert (columns_in_input > 0);
  assert (NULL != ignore_column);

  column_title = new IWString[columns_in_input + 1];    // need an extra one at the end for the last call

  int i = 0;
  int column = 0;
  while (buffer.nextword(column_title[column], i, input_separator))
  {
    if (ignore_column[column])
      ;
    else if (! process_column_regexp.active())
      ;
    else if (! process_column_regexp.matches(column_title[column]))
      ignore_column[column] = 1;
    else if (verbose)
    {
      write_column(column, cerr);
      cerr << " matches the process column regular expression\n";
    }

    column++;
  }

  return 1;
}

static int
read_all_data(iwstring_data_source & input,
               resizable_array_p<IWString> & all_data)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    all_data.add(new IWString(buffer));
  }

  return all_data.number_elements();
}

static int
echo_file_buffered(std::ostream & output)
{
  const int n = all_data.size();

  for (int i = 0; i < n; ++i)
  {
    output << *(all_data[i]) << '\n';
  }

  return output.good();
}

static int
echo_file(iwstring_data_source & input,
           std::ostream & output)
{
  if (buffer_file)
    return echo_file_buffered(output);

  input.seekg(0);

  const_IWSubstring buffer;
  while (input.next_record(buffer) && output.good())
  {
    output << buffer << '\n';
  }

  return output.good();
}

static int
get_next_record (iwstring_data_source & input,
                 const resizable_array_p<IWString> & all_data,
                 int & ndx,
                 const_IWSubstring & buffer)
{
  if (buffer_file)
  {
    if (ndx >= all_data.number_elements())
      return 0;

    buffer = *(all_data[ndx]);
    ndx++;
    return 1;
  }

  return input.next_record(buffer);
}

/*
  Note that if we have too many columns in a record, we will spill off
  the end of the counters array. they should have used tcount
*/

static int
notenoughvariance(iwstring_data_source & input,
                   std::ostream & output)
{
  const_IWSubstring buffer;

  int ndx = 0;    // index into all_data

  if (buffer_file)
    ndx = nskip;
  else
  {
    for (int i = 0; i < nskip; i++)
    {
      input.next_record(buffer);
    }
  }

  int records_examined = 0;    // different from records read with probability

  while (get_next_record(input, all_data, ndx, buffer))
  {
    if (probability < 1.0  && iwrandom() > probability)
      continue;

    records_examined++;
    if (records_to_process > 0 && records_examined >= records_to_process)
      break;

    int i = 0;
    int nwords_this_record = 0;
    int missing_values_this_record = 0;
    const_IWSubstring token;
    while (buffer.nextword(token, i, input_separator))
    {
      if (columns_in_input > 0 && nwords_this_record >= columns_in_input)
      {
        cerr << "Column count mismatch, should be only " << columns_in_input << " columns, at " << input.lines_read() << endl;
        return 71;
      }

      value_count_float & v = counters[nwords_this_record];

      nwords_this_record++;
      if (ignore_column && ignore_column[nwords_this_record - 1])
        continue;

      if (missing_value == token)
        missing_values_this_record++;

//    cerr << "Token " << nwords_this_record << " on line " << input.lines_read() << " is '" << token << "'\n";
      if (! v.extra(token))
      {
        cerr << "Non numeric value, token " << nwords_this_record << " '" << token << "' on line " << input.lines_read() << endl;
        if (ignore_columns_with_non_numeric_data)
        {
          ignore_column[nwords_this_record - 1] = 1;
          continue;
        }
        return 87;
      }
    }

    if (columns_in_input <= 0)
      columns_in_input = nwords_this_record;
    else if (nwords_this_record != columns_in_input)
    {
      cerr << "Column count mismatch, " << columns_in_input << " vs " << nwords_this_record << " at " << input.lines_read() << endl;
      return 71;
    }

    missing_value_counter[missing_values_this_record]++;
  }

// We have now scanned the file, and are ready to make decisions on what to
// remove or keep

  if (verbose)
  {
    cerr << "Read " << records_examined << " records";
    if (records_examined != input.lines_read())
      cerr << " (" << input.lines_read() << " lines read)";
    cerr << ", each with " << columns_in_input << " columns\n";
    if (verbose > 1 || print_distributions)
    {
      for (int i = 0; i < columns_in_input; i++)
      {
        if (ignore_column && ignore_column[i])
          continue;

        value_count_float & v = counters[i];
        write_column(i, cerr);
        cerr << v.number_elements() << " values between ";
        cerr << v.minval() << " and " << v.maxval();
        if (v.n() > 1)
          cerr << " ave " << v.average() << " variance " << v.variance() << " -x " << v.fraction_of_most_common_value();
        if (v.missing_values())
          cerr << ", " << v.missing_values() << " missing values";
        cerr << endl;
        if (print_distributions)
        {
          for (int j = 0; j < v.number_elements(); j++)
          {
            cerr << ' ' << v[j] << '(' << v.count(j) << ')';
          }
          cerr << endl;
        }
      }
    }
  }

// If we are printing statistics, we also gather statistics on the variances, both raw and scaled

  Accumulator<float> vstats, vstats_scaled;
  
  int columns_suppressed = 0;

  int * keep = new int[columns_in_input]; std::unique_ptr<int[]> free_keep(keep);

  for (int i = 0; i < columns_in_input; i++)
  {
    keep[i] = 1;    // by default

    if (ignore_column && ignore_column[i])
      continue;

    value_count_float & v = counters[i];
    if (print_statistics)
    {
      write_column(i, cerr);
      cerr << ' ' << v.number_elements() << " values";
      if (v.number_elements() > 1)
      {
        vstats.extra (v.variance());
        float r = v.rescaled_variance();
        vstats_scaled.extra(r);
        cerr << " between " << v.minval() << " and " << v.maxval() << " ave " << v.average() << " variance " << v.variance();
        if (! v.looks_like_fingerprint_data())     // rescaled variances are the same as the variance for fingerprints
          cerr << " (rescaled " << r << ")";
      }
      if (v.missing_values())
        cerr << ", " << v.missing_values() << " missing values";
      cerr << endl;
    }

    if (1 == v.number_elements())
    {
      keep[i] = 0;
      if (verbose)
      {
        write_column(i, cerr);
        cerr << " is constant " << v.minval() << " - suppressed\n";
      }
      columns_suppressed++;
      continue;
    }

    if (v.number_elements() < min_different_values)
    {
      keep[i] = 0;
      if (verbose)
      {
        write_column(i, cerr);
        cerr << " only has " << v.number_elements() << " distinct values - suppressed\n";
      }
      columns_suppressed++;
      continue;
    }

    if (most_common_max_pop > 0.0 && v.fraction_of_most_common_value() > most_common_max_pop)
    {
      keep[i] = 0;
      if (verbose)
      {
        write_column(i, cerr);
        cerr << "(" << v.number_elements() << " different values) fails fraction test " << v.fraction_of_most_common_value() << " - suppressed\n";
      }
      columns_suppressed++;
      continue;
    }

    if (max_fraction_for_median > 0.0 && v.median_too_common(max_fraction_for_median, i))
    {
      keep[i] = 0;
      if (verbose)
      {
        write_column(i, cerr);
        cerr << "(" << v.number_elements() << " different values) fails median test " << max_fraction_for_median << " - suppressed\n";
      }
      columns_suppressed++;
      continue;
    }

    if (min_values_above_population_threshold > 0 && 
        ! v.several_values_significantly_hit(min_values_above_population_threshold, population_threshold))
    {
      keep[i] = 0;
      if (verbose)
      {
        write_column(i, cerr);
        cerr << "(" << v.number_elements() << " different values) fails population test - suppressed\n";
      }
      columns_suppressed++;
      continue;
    }

//  Remember, the case of 0.0 for maximum_missing_value_fraction is special

    if (0.0 == maximum_missing_value_fraction && v.missing_values())
    {
      keep[i] = 0;
      if (verbose)
      {
        write_column(i, cerr);
        cerr << "(" << v.number_elements() << " different values) contains " << v.missing_values() << " missing values - suppressed\n";
      }
      columns_suppressed++;
      continue;
    }

    if (maximum_missing_value_fraction > 0.0 &&
        static_cast<float>(v.missing_values()) / static_cast<float>(v.n() + v.missing_values()) >= maximum_missing_value_fraction)
    {
      keep[i] = 0;
      if (verbose)
      {
        write_column(i, cerr);
        cerr << "(" << v.number_elements() << " different values, " << v.missing_values() << " missing values) fails missing value fraction test " <<
              (static_cast<float>(v.missing_values()) / static_cast<float>(v.n() + v.missing_values())) << " - suppressed\n";
      }
      columns_suppressed++;
      continue;
    }

    if (min_variance > 0.0 && v.number_elements() > 1 && v.best_variance_measure() < min_variance)
    {
      keep[i] = 0;
      if (verbose)
      {
        write_column(i, cerr);
        cerr << "(" << v.number_elements() << " different values) variance " << v.best_variance_measure() << " (rescaled " << v.rescaled_variance() << " ) below threshold - suppressed\n";
      }

      columns_suppressed++;
      continue;
    }

    if (remove_if_all_greater_than_zero && v.all_values_greater_than_zero())
    {
      keep[i] = 0;

      if (verbose)
      {
        write_column(i, cerr);
        cerr << " all values > 0 - suppressed\n";
      }

      columns_suppressed++;
      continue;
    }

    if (remove_if_all_values_less_than > static_cast<float>(0.0) && (v.maxval() - v.minval()) < remove_if_all_values_less_than)
    {
      keep[i] = 0;

      if (verbose)
      {
        write_column(i, cerr);
        cerr << " all values within fuzz - suppressed\n";
      }

      columns_suppressed++;
      continue;
    }
  }

  if (print_statistics)
  {
    cerr << "Column variances between " << vstats.minval() << " and " << vstats.maxval();
    if (vstats.n() > 1)
      cerr << ". Average variance " << vstats.average();
    cerr << endl;

    cerr << "Rescaled variances between " << vstats_scaled.minval() << " and " << vstats_scaled.maxval();
    if (vstats_scaled.n() > 1)
      cerr << ". Average rescaled variance " << vstats_scaled.average();
    cerr << endl;
  }

  if (highest_variance_nkeep || lowest_x_nkeep)
  {
    do_highest_variance_nkeep(counters, keep, columns_suppressed);
  }

  if (verbose)
    cerr << "Columns suppressed = " << columns_suppressed << endl;

  if (columns_suppressed)    // we need to change the input
    ;
  else if (normalise_output)   // no columns suppressed, but need to process the input
    ;
  else if (write_file_even_if_no_changes)   // no changes, but they want it echo'd
    return echo_file(input, output);
  else       // no changes, no echo requested
  {
    cerr << "No columns suppressed\n";
    return 0;
  }

  int columns_to_be_output = 0;
  for (int i = 0; i < columns_in_input; i++)
  {
    if (keep[i])
      columns_to_be_output++;
  }

  if (0 == columns_to_be_output)
  {
    cerr << "No columns selected for output\n";
    return 0;
  }

  if (verbose)
    cerr << "Selected " << columns_to_be_output << " columns for output\n";

  if (! buffer_file)
    (void) input.seekg(0);

  IWString output_buffer;
  output_buffer.resize(buffer.nchars());

  int nlines = 0;     // need to keep track of whether we are in the nskip lines at the top

  ndx = 0;    // restart all_data

  while (get_next_record(input, all_data, ndx, buffer))
  {
    nlines++;

    int i = 0;
    const_IWSubstring token;
    int iword = 0;
    while (buffer.nextword(token, i, input_separator))
    {
      if (! keep[iword])
      {
        iword++;
        continue;
      }

      if (iword > 0)
        output_buffer += output_separator;

//    figure out all the cases in which we need to convert to an actual number

      if (normalise_output && (nlines > nskip) &&
          (NULL == ignore_column || 0 == ignore_column[iword])
          && missing_value != token)
      {
        float tmp;
        (void) token.numeric_value(tmp);
        value_count_float & v = counters[iword];
        tmp = (tmp - v.average()) / v.variance();
        output_buffer += tmp;
      }
      else
        output_buffer += token;

      iword++;
    }

    output << output_buffer << '\n';

    output_buffer.resize_keep_storage(0);
  }

  return 0;
}


static int
notenoughvariance(const char * fname,
                   std::ostream & output)
{
  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 99;
  }

  if (translate_tabs)
    input.set_translate_tabs(1);

  if (buffer_file)
  {
    if (! read_all_data(input, all_data))
      return 0;

    if (verbose)
      cerr << "Read " << all_data.size() << " records from '" << fname << "'\n";
  }

// Determine the number of tokens per record

  IWString buffer;

  if (buffer_file)
    buffer = *(all_data[0]);
  else if (! input.next_record(buffer))
  {
    cerr << "Cannot read header from '" << fname << "'\n";
    return 61;
  }

  columns_in_input = buffer.nwords();
  if (verbose)
    cerr << "Input contains " << columns_in_input << " columns\n";

  if (0 == columns_in_input)
  {
    cerr << "Very bad news, header record contains zero columns\n";
    return 51;
  }

  if (buffer_file)
    ;
  else if (! input.seekg(0))
  {
    cerr << "Yipes, cannot seek back to beginning of file\n";
    return 13;
  }

  if (ignore_columns_from_cl.number_elements() || process_column_regexp.active() || 
      ignore_columns_with_non_numeric_data)
  {
    ignore_column = new int[columns_in_input];
    for (int i = 0; i < columns_in_input; i++)
    {
      ignore_column[i] = 0;
    }
  }

  for (int i = 0; i < ignore_columns_from_cl.number_elements(); i++)
  {
    int j = ignore_columns_from_cl[i];
    ignore_column[j] = 1;
  }

  if (is_descriptor_file)
    (void) establish_column_titles(buffer);

  counters = new value_count_float[columns_in_input]; std::unique_ptr<value_count_float[]> free_counters(counters);

  for (int i = 0; i < columns_in_input; i++)
  {
    value_count_float & v = counters[i];

    if (NULL == ignore_column || ! ignore_column[i])
      v.initialise(nkeep);
  }

  if (! buffer_file)
    input.reset_record_counter();

  return notenoughvariance(input, output);
}

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Usage:  <options> <input_file>\n";
  cerr << " -n <buffer>    specify buffer size - number of different values stored (mandatory)\n";
  cerr << " -m <string>    specify missing value string (default '" << missing_value << "')\n";
  cerr << " -p <prob>      sample records with probability <p>\n";
  cerr << " -d <number>    minimum number of different values in a column (default 2)\n";
  cerr << " -x <fraction>  highest fraction of hits for most common value\n";
  cerr << "                -x 0.8 means reject that column if the most common\n";
  cerr << "                value is 80% or more of the values found\n";
  cerr << endl;
  cerr << " -x <keep=n>    keep only N columns with the lowest -x value\n";
  cerr << " -f <fraction>  specify fraction for -b option\n";
  cerr << " -b <number>    used with -f option\n";
  cerr << "                -b 3 -f 0.1 means reject columns for which fewer than three\n";
  cerr << "                values are hit more than 10% of the time\n";
  cerr << "                Note that -b 3 -f 0.5 is impossible - cannot have\n";
  cerr << "                three values each hitting 50% of the time\n";
  cerr << endl;
  cerr << " -B <fraction>  binary collapse fraction\n";
  cerr << " -e <fraction>  reject if median more than <fraction>\n";
  cerr << endl;
  cerr << " -V <minvar>    reject any column with variance less than <minvar>\n";
  cerr << " -V <keep=n>    keep the N most desirable columns, desirability from HV\n";
  cerr << " -V DS=N        most desirable columns are the ones with most different values\n";
  cerr << " -V DS=ACV      most desirable columns are the ones with highest adjusted CV\n";
  cerr << " -V DS=ACV=x    most desirable columns are the ones with highest adjusted CV\n";
  cerr << "                  use X as the adjusted median threshold\n";
  cerr << " -V DS=SCA      most desirable columns are the ones with highest adjusted scaleA\n";
  cerr << endl;
  cerr << " -W             remove columns for which every row is non-zero\n";
  cerr << " -c             once buffer is full, count new values with their nearest neighbour\n";
  cerr << " -N             normalise output (mean 0.0, variance 1.0)\n";
  cerr << " -a             print column statistics (min, max, ave) and missing value stats\n";
  cerr << " -t             print distributions (not recommended for large buffers)\n";
  cerr << " -j             is descriptor file (implies -s 1 -i 1)\n";
  cerr << " -s <skip>      skip the first <skip> records of the file\n";
  cerr << " -i <column>    specify columns to ignore (starting @ 1)\n";
  cerr << " -P <regexp>    only process columns which match <regexp> (implies -j)\n";
  cerr << " -q <number>    only process the first <number> records in the file\n";
  cerr << " -g <fraction>  skip column of <fraction> or more missing values\n";
  cerr << " -g 0.0         skip column if ANY missing values\n";
  cerr << " -w             write file even if no changes needed\n";
  cerr << " -T             translate tabs to spaces\n";
  cerr << " -G             ignore columns containing non-numeric data (normally a fatal error)\n";
  cerr << " -r <range>     remove columns where the range is less than <range>\n";
  cerr << " -u             buffer all the input (avoids seeks, enables stdin)\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

#include "cmdline.h"

/*
  The -V option is so complex that we put it in its own function
*/

static int
parse_dash_V_option(Command_Line & cl, int verbose)
{
  const_IWSubstring v;
  int i = 0;
  while (cl.value('V', v, i++))
  {
    if (v.starts_with("keep="))
    {
      v.remove_leading_chars(5);
      if (! v.numeric_value(highest_variance_nkeep) || highest_variance_nkeep < 1)
      {
        cerr << "The '-V nkeep=' combination must be followed by a positive whole number\n";
        usage(18);
      }

      if (verbose)
        cerr << "The " << highest_variance_nkeep << " columns with the greatest variance will be retained\n";
    }
    else if (v.starts_with("DS="))
    {
      v.remove_leading_chars(3);     // get rid of DS=
  
      if ('N' == v)
      {
        chosen_highest_variance_measure = HV_MOST_VALUES;
        if (verbose)
          cerr << "Columns with the most values will be retained\n";
      }
      else if (v.starts_with("ACV"))
      {
        chosen_highest_variance_measure = HV_CV;
        if (verbose)
          cerr << "Will keep columns with highest adjusted CV values\n";
        if (v.starts_with("ACV="))
        {
          v.remove_leading_chars(4);
          if (! v.numeric_value(adjusted_median_threshold) || adjusted_median_threshold <= 0.0)
          {
            cerr << "Invalid numeric specifier for adjusted median threshold 'DS=ACV=" << v << "'\n";
            usage(4);
          }
    
          if (verbose)
            cerr << "  adjusted median threshold " << adjusted_median_threshold << endl;
        }
      }
      else if ("SCA" == v)
      {
        chosen_highest_variance_measure = HV_SCALEA;
        if (verbose)
          cerr << "Columns with highest adjusted scaleA values will be retained\n";
      }
      else
      {
        cerr << "Unrecognised column desirability qualifier 'DS=" << v << "'\n";
        usage(21);
      }
    }
    else if (! v.numeric_value(min_variance) || min_variance <= 0.0 || min_variance >= 1.0)
    {
      cerr << "The -V option (min variance) must be followed by a positive number between 0 and 1\n";
      usage(42);
    }
    else if (verbose)
      cerr << "Columns with a variance < " << min_variance << " will be eliminated\n";
  }

  if (0 == highest_variance_nkeep)
  {
    cerr << "No 'keep=nn' specifier entered\n";
    usage(6);
  }

  return 1;
}

static int 
notenoughvariance(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "g:q:jae:Ncf:b:B:tTvn:m:s:i:p:P:d:x:wWV:Gr:u");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('n'))
  {
    cerr << "Must specify value for buffer size via -n option\n";
    usage(2);
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', nkeep) || nkeep <= 0)
    {
      cerr << "The buffer size must be a whole positive number\n";
      usage(3);
    }
  }

  if (cl.option_present('t'))
  {
    print_distributions = 1;
    if (verbose)
      cerr << "Will print full distributions !! - this could be large!\n";
  }

  if (cl.option_present('T'))
  {
    translate_tabs = 1;
    if (verbose)
      cerr << "Tabs translated to spaces\n";
  }

  if (cl.option_present('a'))
  {
    print_statistics = 1;
    if (verbose)
      cerr << "Will print column statistics\n";
  }

  if (cl.option_present('V'))
  {
    if (! parse_dash_V_option(cl, verbose))
    {
      return 8;
    }
  }

  if (cl.option_present('G'))
  {
    ignore_columns_with_non_numeric_data = 1;

    if (verbose)
      cerr << "Columns with non-numeric data will be discarded rather than causing a fatal error\n";
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', nskip) || nskip < 1)
    {
      cerr << "the skip option requires a whole positive number\n";
      usage(43);
    }

    if (verbose)
      cerr << "Will skip first " << nskip << " records of the file\n";
  }

  if (cl.option_present('j'))
  {
    is_descriptor_file = 1;
    nskip = 1;
    ignore_columns_from_cl.add(0);
    if (verbose)
      cerr << "Will process as a descriptor file\n";
  }

  if (cl.option_present('m'))
  {
    cl.value('m', missing_value);
    if (verbose)
      cerr << "Missing value string '" << missing_value << "'\n";
  }

  if (cl.option_present('p'))
  {
    if (! cl.value('p', probability) || probability <= 0.0 || probability >= 1.0)
    {
      cerr << "The probability value must be between 0.0 and 1.0\n";
      usage(14);
    }

    (void) iw_random_seed();

    if (verbose)
      cerr << "Will examine records with a probability of " << probability << endl;
  }

  if (cl.option_present('d'))
  {
    if (! cl.value('d', min_different_values) || min_different_values < 2)
    {
      cerr << "The -d switch requires a value >= 2\n";
      usage(18);
    }

    if (verbose)
      cerr << "Columns with fewer than " << min_different_values << " different values will be suppressed\n";
  }

  if (cl.option_present('x'))
  {
    int i = 0;
    const_IWSubstring x;
    while (cl.value('x', x, i++))
    {
      if (x.starts_with("keep="))
      {
        x.remove_leading_chars(5);
        if (! x.numeric_value(lowest_x_nkeep) || lowest_x_nkeep < 1)
        {
          cerr << "The '-x keep=' combination must be followed by a whole number\n";
          usage(5);
        }

        if (verbose)
          cerr << "Will keep the " << lowest_x_nkeep << " columns with the lowest -x value\n";
      }
      else if (! cl.value('x', most_common_max_pop) || most_common_max_pop <= 0.0 || most_common_max_pop >= 1.0)
      {
        cerr << "The most common value's max fraction (-x) must be between 0.0 and 1.0\n";
        usage(16);
      }
      else if (verbose)
        cerr << "Suppress column if the most common value is " << most_common_max_pop << " or more of the values\n";
    }
  }

  if (cl.option_present('e'))
  {
    if (! cl.option_present('o'))
    {
      cerr << "Sorry, the median requires sorted values, use -o\n";
      usage(32);
    }

    if (! cl.value('e', max_fraction_for_median) || max_fraction_for_median < 0.5)
    {
      cerr << "The max fraction for median (-e) value must be > 0.5 and < 1.0\n";
      usage(19);
    }
    if (verbose)
      cerr << "Suppress column if median is " << max_fraction_for_median << " or more of the values\n";
  }

  if (cl.option_present('g'))
  {
    if (! cl.value('g', maximum_missing_value_fraction) || 
        maximum_missing_value_fraction > 1.0 || maximum_missing_value_fraction < 0.0)
    {
      cerr << "The -g option requires a fraction\n";
      usage(47);
    }

    if (0 == verbose)
      ;
    else if (0.0 == maximum_missing_value_fraction)
      cerr << "Will reject any column containing any missing values\n";
    else
      cerr << "Will reject columns having " << maximum_missing_value_fraction << " or more missing values\n";
  }

  if (cl.option_present('q'))
  {
    if (! cl.value('q', records_to_process) || records_to_process < 0)
    {
      cerr << "The -q option requires a whole positive number\n";
      usage(39);
    }

    if (verbose)
      cerr << "Will only process the first " << records_to_process << " records of the input\n";
  }

  if (cl.option_present('i'))
  {
    ignore_columns_from_cl.resize(cl.option_count('i'));

    int i = 0;
    int j;
    while (cl.value('i', j, i++))    // does not handle errors well
    {
      if (j < 1)
      {
        cerr << "Column numbers start with column 1, value " << j << " is invalid\n";
        return 73;
      }

      ignore_columns_from_cl.add(j - 1);
    }
  }

  if (cl.option_present('P'))
  {
    const_IWSubstring p;

    (void) cl.value('P', p);

    if (! process_column_regexp.set_pattern(p))
    {
      cerr << "Cannot parse process column regexp '" << p << "'\n";
      usage(4);
    }

    is_descriptor_file = 1;

    if (verbose)
      cerr << "Only columns matching '" << process_column_regexp.source() << "' will be processed\n";
  }

  if (cl.option_present('b'))
  {
    if (! cl.option_present('f'))
    {
      cerr << "the -b and -f options must be used together\n";
      usage(91);
    }

    if (! cl.value('b', min_values_above_population_threshold) ||
        min_values_above_population_threshold < 1)
    {
      cerr << "The -b option requires a whole number > 1\n";
      usage(6);
    }
  }

  if (cl.option_present('B'))
  {
    if (! cl.value('B', binary_collapse_threshold) || binary_collapse_threshold <= 0.0 || binary_collapse_threshold >= 1.0)
    {
      cerr << "The binary denisity collapse threshold (-B) option must be followed by a valid ratio\n";
      usage(14);
    }

    if (verbose)
      cerr << "Binary collapse denisty set to " << binary_collapse_threshold << endl;
  }

  if (cl.option_present('f'))
  {
    if (! cl.option_present('b'))
    {
      cerr << "the -f and -b options must be used together\n";
      usage(92);
    }

    if (! cl.value('f', population_threshold) || population_threshold < 0.0 || population_threshold >= 1.0)
    {
      cerr << "The population threshold (-f) must be > 0.0 and < 1.0\n";
      usage(24);
    }

    if (verbose)
      cerr << "At least " << min_values_above_population_threshold << " values must be above " <<
               population_threshold << endl;
  }

  if (cl.option_present('c'))
  {
    if (! sort_stored_values)
    {
      cerr << "Sorry, the store with closest (-c) option only works with sorted values\n";
      usage(61);
    }

    store_with_closest_neighbour = 1;
    if (verbose)
      cerr << "Extra values will be stored with their closest neighbour\n";
  }

  if (cl.option_present('N'))
  {
    normalise_output = 1;
    if (verbose)
      cerr << "Output will be normalised (mean 0.0, variance 1.0)\n";
  }

  if (cl.option_present('w'))
  {
    write_file_even_if_no_changes = 1;
    if (verbose)
      cerr << "File will be written even if no changes needed\n";
  }

  if (cl.option_present('W'))
  {
    remove_if_all_greater_than_zero = 1;

    if (verbose)
      cerr << "Will remove columns that have all non-zero values\n";
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', remove_if_all_values_less_than) || remove_if_all_values_less_than <= 0.0)
    {
      cerr << "The remove small range (-r) option must be followed by a positive number\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will remove columns with ranges less than " << remove_if_all_values_less_than << endl;
  }

  if (cl.option_present('u'))
  {
    buffer_file = 1;

    if (verbose)
      cerr << "Will buffer all input\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(3);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Extra arguments '" << cl[1] << "...' ignored\n";
  }

  int rc = notenoughvariance(cl[0], std::cout);

  if (verbose && print_statistics)
  {
    for (int i = 0; i < missing_value_counter.number_elements(); i++)
    {
      if (missing_value_counter[i])
        cerr << missing_value_counter[i] << " records had " << i << " missing values\n";
    }
  }

  if (NULL != ignore_column)
    delete [] ignore_column;

  if (NULL != column_title)
    delete [] column_title;

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = notenoughvariance(argc, argv);

  return rc;
}
