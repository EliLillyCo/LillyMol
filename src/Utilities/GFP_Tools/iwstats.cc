/*
  Compute summary statistics between experimental and predicted numeric values.
*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <random>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Metric.h"
#include "bsquared.h"

using std::cerr;
using std::cout;
using std::endl;

static Enrichment enrichment;

static IWString missing_value;

static int report_missing_values = 1;

static int missing_values_encountered = 0;

static int randomise_ties = 0;

/*
  Nov 2002. Why not allow experimental values to come from another file
*/

static IW_STL_Hash_Map<IWString, IWString> id_activity_hash;

static int strip_leading_zeros_from_identifiers = 0;

static int truncate_predicted_values_to_experimental_range = 0;

static int values_lower_than_experimental_range = 0;
static int values_higher_than_experimental_range = 0;

static int traditional_q2_definition = 0;

static float relative_error_threshold = static_cast<float>(0.0);

static int range_normal_relative_error = 0;

static char input_separator = ' ';

static IWString activity_name;

/*
  When missing values are present, we can impose a minimum number of values
  present for performing correlations
*/

static unsigned int values_needed_for_reporting_correlations = 20;

/*
  Sometimes it is just easier to skip items with missing experimental data
*/

static int skip_items_with_no_activity_data = 0;

static int items_with_no_activity_data = 0;

static float active_inactive_cutoff = 0.5;
static int calculate_enrichment_metrics = 0;

static double bedroc_alpha = 20.0;

static double enrichment_factor_fraction = 0.5;

static double discard_measured_values_below = -std::numeric_limits<double>::max();

static int proper_median = 1;

static int do_cMSD = 0;

/*
  Previous R2 definition could be outside the range of 0 to 1.
*/

static int wikipedia_r2_definition = 1;

/*
  ADME group interested in number of predictions outside N fold difference
  from obs
*/

static resizable_array<float> fold_difference_additive;
static resizable_array<float> fold_difference_multiplicative;

/*
  Useful to be able to collect residuals
*/

static IWString_and_File_Descriptor stream_for_residuals;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(Computes Bsquared and other summary statistics - allows missing values
Compares an experimental set of values (-E,-e) with one more predicted sets of values (-p)
generating summary statistics between the sets of numbers.
Can have the experimental data in a separate file, or as a column in the predicted file.
 -E <file>      activities are in a different file <file>
 -e <col>       column for experimental (measured) values - either
 -z             strip leading 0's from identifiers when using -E
 -p <col>       column for predicted values - can specify any number
 -s <number>    skip <number> records at the top of each file
 -j             treat as a descriptor file
 -n <number>    process only the first/best <number> records
 -t <number>    compute Bsquared for sections of the data. For example,
                -t 10 would report values for the first/best 10, 20, 30 values
 -P <number>    sample the first/best <number> percent of the data
 -M <string>    missing value string
 -q             quietly ignore missing values
 -R .           randomise the sort when duplicate predicted values present
 -c <number>    number of valid pairs needed for producing correlations (default " << values_needed_for_reporting_correlations << ")
 -w             when duplicate values present, suppress computation of best
                and worst BSquared values
 -r <float>     max relative error allowed computations
 -T             truncate predicted values to the experimental range
 -h             use traditional (Wikipedia) Q2 formula
 -b <nbucket>   compute distribution functions across <nbucket> buckets
 -F <f>         calculate number of predictions outside <f> fold of experimental (multiplicative)
 -D <f>         calculate number of predictions differing by <d> from experimental (additive)
 -k             just skip predicted values that have no experimental value (repeat for quiet)
 -d             compute Dave's cMSD rather than cMSR
 -o <float>     cutoff for active/inactive (for enrichment metrics BEDROC EF ...)
 -a <float>     BEDROC alpha value (default " << bedroc_alpha << ")
 -f <float>     Enrichment Factor default fraction (default " << enrichment_factor_fraction << ")
 -u <float>     discard experimental values below <float> - useful for studying actives only 
 -U <fraction>  keep only the <fraction> most active experimental values
 -m <col/name>  do analysis by data in column <col> or descriptor name <name>
 -L <fname>     write residuals to <fname>
 -v             verbose output
)";
  // clang-format on

  exit(rc);
}

static int verbose = 0;

static int records_to_skip = 0;

/*
  We can chop up the data in several ways.
  With step, we can do every n'th group (first 10, first 20, first 30...)
  With percents, we can specify any number of percents..
  With N_to_check we check specific numbers
*/

static int step = 0;
static resizable_array<int> percents_to_check;
static resizable_array<int> N_to_check;
static resizable_array<float> activities_to_check;

/*
  If someone asks for percents, they should have the resulting output show nn%
*/

static resizable_array_p<IWString> chunk_titles;

static int is_descriptor_file = 0;

static int columns_in_input = 0;

static IWString* column_titles = nullptr;

/*
  Since Bsquared can vary depending on the input order, we can do various things...
*/

static int compute_range_of_bsquared_values_when_duplicates_present = 1;

static int number_distribution_buckets = 0;

template <typename T>
T
median(const T* v, int n) {
  if (n == (n / 2 * 2)) {  // N is even, take ave of two in middle
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
  } else {
    return v[n / 2];
  }
}

/*
  Write the descriptor name if we are processing a descriptor file

  When processing the marker_column things, I found a need to have
  some arbitrary text written out before every item. Bit of a kludge
  but the alternative would be to add another argument.
*/

static IWString global_prefix;

template <typename T>
void
write_something_identifying_the_column(int zcol, T& output) {
  if (global_prefix.length() > 0) {
    output << global_prefix;
  }

  if (is_descriptor_file) {
    output << ' ' << column_titles[zcol];
  } else if (zcol >= 0) {
    output << " column " << (zcol + 1);
  }

  return;
}

/*
  We are only really interested in the columns that contain
  predicted values
*/

static int
establish_column_titles(const const_IWSubstring& buffer, int experimental_column,
                        const resizable_array<int>& predicted_column, int& marker_column,
                        const IWString& marker_column_name) {
  if (nullptr != column_titles) {
    delete[] column_titles;
  }

  columns_in_input = buffer.nwords(input_separator);

  if (0 == columns_in_input) {
    cerr << "Huh, no tokens in descriptor file header\n";
    return 0;
  }

  if (experimental_column >= columns_in_input) {
    cerr << "Invalid experimental column " << experimental_column << " only "
         << columns_in_input << " columns in input\n";
    return 0;
  }

  for (int i = 0; i < predicted_column.number_elements(); i++) {
    int p = predicted_column[i];
    if (p >= columns_in_input) {
      cerr << "Invalid predicted value column " << p << " only " << columns_in_input
           << " columns in input\n";
      return 0;
    }
  }

  if (!is_descriptor_file) {
    return 1;
  }

  column_titles = new IWString[columns_in_input];

  int col = 0;
  int i = 0;
  const_IWSubstring token;
  while (buffer.nextword(token, i, input_separator)) {
    column_titles[col] = token;

    // expt column doesn't really mean the same thing when just a name cross reference.
    if (col == experimental_column && id_activity_hash.empty()) {
      cerr << "Experimental column is '" << column_titles[col] << "'\n";
    }

    if (marker_column_name.length() > 0 && marker_column < 0 &&
        token == marker_column_name) {
      marker_column = col;
      if (verbose) {
        cerr << "Found marker descriptor '" << marker_column_name << "' in column "
             << (col + 1) << '\n';
      }
    }

    col++;
  }

  if (experimental_column >= 0 && 0 == activity_name.length()) {
    activity_name = column_titles[experimental_column];
  }

  return 1;
}

/*
  We need an object to hold the raw data.
*/

typedef float experimental_value_t;

class Predicted_Values {
 private:
  experimental_value_t _experimental;

  //  to account for missing values, we mark each value as valid or not

  int* _valid;
  experimental_value_t* _predicted;

  //  when doing randomised breaking of ties, we need a random number with each prediction

  experimental_value_t _random;

  IWString _marker;

 public:
  Predicted_Values();
  ~Predicted_Values();

  int debug_print(std::ostream&) const;

  experimental_value_t
  obs() const {
    return _experimental;
  }

  experimental_value_t
  pred(int c) const {
    return _predicted[c];
  }

  int
  valid(int i) const {
    return _valid[i];
  }

  experimental_value_t
  random_value() const {
    return _random;
  }

  void
  assign_random_value(experimental_value_t r) {
    _random = r;
  }

  const IWString
  marker() const {
    return _marker;
  }

  int parse_buffer(const const_IWSubstring& buffer, int experimental_column,
                   const resizable_array<int>& predicted_column, int marker_column,
                   int& fatal);
};

Predicted_Values::Predicted_Values() {
  assert(columns_in_input > 0);

  _experimental = 0.0;
  _valid = new_int(columns_in_input);
  _predicted = new experimental_value_t[columns_in_input];

  _random = static_cast<experimental_value_t>(0.0);

  return;
}

Predicted_Values::~Predicted_Values() {
  if (nullptr != _valid) {
    delete[] _valid;
  }

  if (nullptr != _predicted) {
    delete[] _predicted;
  }

  return;
}

static int
fetch_activity(const IW_STL_Hash_Map<IWString, IWString>& id_activity_hash, IWString& s) {
  if (strip_leading_zeros_from_identifiers) {
    s.remove_leading_chars('0');
  }

  IW_STL_Hash_Map<IWString, IWString>::const_iterator f = id_activity_hash.find(s);

  if (f == id_activity_hash.end()) {
    return 0;
  }

  s = (*f).second;

  return 1;
}

int
Predicted_Values::parse_buffer(const const_IWSubstring& buffer, int experimental_column,
                               const resizable_array<int>& predicted_column,
                               int marker_column, int& fatal) {
  int i = 0;
  IWString token;
  int col = 0;

  while (buffer.nextword(token, i, input_separator)) {
    int prediction_index = -1;

    if (col == experimental_column) {
      if (id_activity_hash.size()) {
        if (!fetch_activity(id_activity_hash, token)) {
          items_with_no_activity_data++;
          if (skip_items_with_no_activity_data) {
            fatal = 0;
          } else {
            fatal = 1;
          }

          if (skip_items_with_no_activity_data < 2) {
            cerr << "Cannot find activity data for '" << token << "'\n";
          }

          return 0;
        }
      }
    } else if (col == marker_column) {
      _marker = token;
      col++;
      continue;
    } else if ((prediction_index = predicted_column.index(col)) >= 0) {
      ;
    } else {
      col++;
      continue;
    }

    if (missing_value.length() && missing_value == token) {
      if (report_missing_values) {
        cerr << "Ignoring missing value in column " << (col + 1) << '\n';
      }
      missing_values_encountered++;
      col++;
      continue;
    }

    double v;
    if (!token.numeric_value(v))  // great
    {
      cerr << "Invalid numeric value '" << token << "'\n";
      if (col == experimental_column && id_activity_hash.size() > 0) {
        cerr << "Error is in activity file\n";
      } else {
        cerr << buffer << '\n';
      }

      fatal = 1;
      return 0;
    }

    if (col == experimental_column) {
      _experimental = v;
    } else {
      _valid[prediction_index] = 1;
      _predicted[prediction_index] = v;
    }

    col++;
  }

  if (col != columns_in_input) {
    cerr << "Column count mismatch, found " << col << " columns, expected "
         << columns_in_input << '\n';
    fatal = 1;
    return 0;
  }

  return 1;
}

int
Predicted_Values::debug_print(std::ostream& os) const {
  os << "expt: " << _experimental;
  if (nullptr != _predicted) {
    os << " pred:";
    for (int i = 0; i < columns_in_input; i++) {
      if (_valid[i]) {
        os << ' ' << _predicted[i];
      } else {
        os << " *";
      }
    }
  }

  os << '\n';

  return os.good();
}

template class resizable_array_p<Predicted_Values>;
template class resizable_array_base<Predicted_Values*>;

/*
  Are there any duplicate predicted values?
*/

static int
duplicate_predicted_values(const resizable_array_p<Predicted_Values>& zdata,
                           int predicted_column) {
  experimental_value_t predprev = zdata[0]->pred(predicted_column);

  int n = zdata.number_elements();

  for (int i = 1; i < n; i++) {
    experimental_value_t pred = zdata[i]->pred(predicted_column);

    if (pred == predprev) {
      return 1;
    }

    predprev = pred;
  }

  return 0;
}

/*
  The -A option has been used to specify a number of activities to report.
  Convert those activity thresholds to numbers of records to check
*/

static int
determine_numbers_to_check_by_activity(const resizable_array_p<Predicted_Values>& zdata,
                                       resizable_array<int>& numbers_to_check) {
#ifdef NOT_IMPLEMENTED
  int n = zdata.number_elements();

  for (int i = 0; i < activities_to_check.number_elements(); i++) {
    float a = activities_to_check[i];

    int nsample = 0;

    for (int i = 0; i < n; i++) {
    }
  }
#endif

  return 1;
}

/*
  We need to determine how many points to check - the -t option.
  Remember, our convention is that if the number is positive, we
  just return that number. If it is negative, it is assumed to be
  a percentage
*/

static int
determine_numbers_to_check(resizable_array<int>& numbers_to_check, int records_read) {
  numbers_to_check.resize_keep_storage(0);
  chunk_titles.resize_keep_storage(0);

  if (step > 0) {
    if (step > records_read) {
      cerr << "Only read " << records_read << " records, step " << step << " too large\n";
    } else {
      for (int i = step; i <= records_read; i += step) {
        numbers_to_check.add(i);
        IWString* tmp = new IWString;
        tmp->append_number(i);
        tmp->operator+=(" observations");

        chunk_titles.add(tmp);
      }
    }
  }

  if (percents_to_check.number_elements()) {
    for (int i = 0; i < percents_to_check.number_elements(); i++) {
      int p = percents_to_check[i];

      // use a number < 100 to give a result just larger than the likely integer
      int n = int(static_cast<float>(p * records_read) / 99.999);

      if (n > 1 && n <= records_read) {
        if (verbose > 1) {
          cerr << p << " percent is " << n << " observations\n";
        }

        numbers_to_check.add(n);

        IWString* tmp = new IWString;
        tmp->append_number(p);
        tmp->operator+=("% observations");

        chunk_titles.add(tmp);
      } else if (n <= 1) {
        cerr << p << " percent of " << records_read << " is " << n
             << " which is ignored\n";
      } else {
        cerr << "What am I supposed to do with " << n << " records\n";
      }
    }
  }

  if (N_to_check.number_elements()) {
    for (int i = 0; i < N_to_check.number_elements(); i++) {
      int n = N_to_check[i];
      if (n < 2 || n > records_read) {
        cerr << "Skipping invalid number of records to check " << n << '\n';
        continue;
      }

      numbers_to_check.add(n);
      IWString* tmp = new IWString;
      tmp->append_number(n);
      tmp->operator+=(" observations");

      chunk_titles.add(tmp);
    }
  }

  return numbers_to_check.number_elements();
}

class Predicted_Value_Comparitor {
 protected:
  int _predicted_column;

 public:
  Predicted_Value_Comparitor(int p) : _predicted_column(p) {
  }

  int operator()(const Predicted_Values*, const Predicted_Values*);
};

int
Predicted_Value_Comparitor::operator()(const Predicted_Values* pv1,
                                       const Predicted_Values* pv2) {
  experimental_value_t p1 = pv1->pred(_predicted_column);
  experimental_value_t p2 = pv2->pred(_predicted_column);

  // cerr << "Comparing " << p1 << " with " << p2 << '\n';

  if (p1 > p2) {
    return -1;
  }

  if (p1 < p2) {
    return 1;
  }

  if (randomise_ties) {
    p1 = pv1->random_value();
    p2 = pv2->random_value();
    if (p1 > p2) {
      return -1;
    }
    if (p1 < p2) {
      return 1;
    }
  }

  return 0;
}

class Predicted_Value_Comparitor_Best : public Predicted_Value_Comparitor {
 private:
 public:
  Predicted_Value_Comparitor_Best(int p) : Predicted_Value_Comparitor(p){};

  int operator()(const Predicted_Values*, const Predicted_Values*);
};

int
Predicted_Value_Comparitor_Best::operator()(const Predicted_Values* pv1,
                                            const Predicted_Values* pv2) {
  experimental_value_t p1 = pv1->pred(_predicted_column);
  experimental_value_t p2 = pv2->pred(_predicted_column);

  if (p1 > p2) {
    return -1;
  }

  if (p1 < p2) {
    return 1;
  }

  experimental_value_t e1 = pv1->obs();
  experimental_value_t e2 = pv2->obs();

  if (e1 > e2) {
    return -1;
  }

  if (e1 < e2) {
    return 1;
  }

  return 0;
}

class Predicted_Value_Comparitor_Worst : public Predicted_Value_Comparitor {
 private:
 public:
  Predicted_Value_Comparitor_Worst(int p) : Predicted_Value_Comparitor(p){};

  int operator()(const Predicted_Values*, const Predicted_Values*);
};

int
Predicted_Value_Comparitor_Worst::operator()(const Predicted_Values* pv1,
                                             const Predicted_Values* pv2) {
  experimental_value_t p1 = pv1->pred(_predicted_column);
  experimental_value_t p2 = pv2->pred(_predicted_column);

  if (p1 > p2) {
    return -1;
  }

  if (p1 < p2) {
    return 1;
  }

  experimental_value_t e1 = pv1->obs();
  experimental_value_t e2 = pv2->obs();

  if (e1 > e2) {
    return 1;
  }

  if (e1 < e2) {
    return -1;
  }

  return 0;
}

class Float_Comparator {
 private:
 public:
  int operator()(const float& f1, const float& f2) const;
};

int
Float_Comparator::operator()(const float& f1, const float& f2) const {
  if (f1 < f2) {
    return -1;
  }

  if (f1 > f2) {
    return 1;
  }

  return 0;
}

static int
compute_distribution_function(int number_records, const float* v, const float minval,
                              const float maxval, int* d) {
  if (minval >= maxval) {
    cerr << "compute_distribution_function:unsorted or no variability: minval " << minval
         << ", maxval " << maxval << '\n';
    return 0;
  }

  float delta = (maxval - minval) / static_cast<float>(number_distribution_buckets);

  for (int i = 0; i < number_records; i++) {
    int j = static_cast<int>((v[i] - minval) / delta + 0.49999);
    //  cerr << "Bucket " << j << " max " << number_distribution_buckets << '\n';
    if (j >= number_distribution_buckets) {
      j = number_distribution_buckets - 1;
    } else if (j < 0) {
      j = 0;
    }

    d[j]++;
  }

  return 1;
}

static int
write_distribution_function(const char* prefix, int number_records, const int* d,
                            std::ostream& output) {
  output << prefix;

  for (int i = 0; i < number_distribution_buckets; i++) {
    float f = d[i] / static_cast<float>(number_records);
    output << ' ' << f;
  }

  output << '\n';

  return output.good();
}

static int
do_compute_distribution_functions(int number_records, int which_predicted_set,
                                  int predicted_column, int experimental_column,
                                  float* tmp1, float* tmp2,
                                  const resizable_array_p<Predicted_Values>& zdata,
                                  std::ostream& output) {
  assert(experimental_column >= 0);
  assert(number_distribution_buckets > 1);

  float min_expt, max_expt, min_pred, max_pred;

  max_pred = 0.0f;  // keep compiler quiet
  min_pred = 0.0f;  // keep compiler quiet
  min_expt = 0.0f;  // keep compiler quiet
  max_expt = 0.0f;  // keep compiler quiet

  float obs, pred;  // scope here for efficiency

  for (int i = 0; i < number_records; i++) {
    const Predicted_Values& pvi = *(zdata[i]);

    obs = pvi.obs();
    pred = pvi.pred(which_predicted_set);

    if (0 == i) {
      min_expt = obs;
      max_expt = obs;

      min_pred = pred;
      max_pred = pred;
    } else {
      if (obs < min_expt) {
        min_expt = obs;
      } else if (obs > max_expt) {
        max_expt = obs;
      }

      if (pred < min_pred) {
        min_pred = pred;
      } else if (pred > max_pred) {
        max_pred = pred;
      }
    }

    tmp2[i] = pred;
    tmp1[i] = obs;

    //  cerr << " i = " << i << " obs " << obs << " pred " << pred << '\n';
  }

  // cerr << "min obs " << min_expt << " max obs " << max_expt << " min pred " << min_pred
  // << " max pred " << max_pred << '\n';

  int* d1 = new_int(number_distribution_buckets);
  std::unique_ptr<int> free_d1(d1);
  int* d2 = new_int(number_distribution_buckets);
  std::unique_ptr<int> free_d2(d2);

  // cerr << "Computing experimental distribution function\n";
  if (!compute_distribution_function(number_records, tmp1, min_expt, max_expt, d1)) {
    return 0;
  }

  // cerr << "Computing predicted distribution function\n";
  if (!compute_distribution_function(number_records, tmp2, min_pred, max_pred, d2)) {
    return 0;
  }

  write_distribution_function("Dexpt:", number_records, d1, output);
  write_distribution_function("Dpred:", number_records, d2, output);

  output << "Compare:";
  float sum = static_cast<float>(0.0);

  for (int i = 0; i < number_distribution_buckets; i++) {
    if (0 == d1[i]) {
      output << " 0";
      continue;
    }

    float t = static_cast<float>(d2[i]) / static_cast<float>(d1[i]);
    output << ' ' << t;
    if (d1[i] < d2[i]) {
      sum += 1.0 / t;
    } else {
      sum += t;
    }
  }

  output << '\n';

  output << "Distributional similarity "
         << sum / static_cast<float>(number_distribution_buckets) << '\n';

  return output.good();
}

static int
initilize_enrichment(int number_records, int experimental_column,
                     const resizable_array_p<Predicted_Values>& zdata, float* tmp1,
                     int which_predicted_set) {
  int n = 0;

  if (experimental_column >= 0) {
    for (int i = 0; i < number_records; i++) {
      const Predicted_Values& pvi = *(zdata[i]);

      if (!pvi.valid(which_predicted_set)) {
        continue;
      }

      tmp1[n] = pvi.obs();
      n++;
    }
    return enrichment.buildFromSortedArray(tmp1, n, active_inactive_cutoff);
  } else {
    return 0;
  }
}

static int
compute_b_squared(int number_records, int experimental_column,
                  const resizable_array_p<Predicted_Values>& zdata, float* tmp1,
                  float* tmp2, double& bsquared, int which_predicted_set) {
  int n = 0;

  if (experimental_column >= 0) {
    for (int i = 0; i < number_records; i++) {
      const Predicted_Values& pvi = *(zdata[i]);

      if (!pvi.valid(which_predicted_set)) {
        continue;
      }

      if (pvi.obs() < discard_measured_values_below) {
        continue;
      }

      tmp1[n] = pvi.obs();
      n++;
    }
  } else {
    set_vector(tmp1, number_records, static_cast<float>(0.0));
    n = number_records;
  }

// #define CHECK_BSQUARED_COMPUTATION
#ifdef CHECK_BSQUARED_COMPUTATION
  double check;
  compute_b_squared(tmp1, tmp2, n, check);
  cerr << "N = " << n << ", whole set Bsquared " << check << '\n';
#endif

  // cerr << "From " << number_records << " found " << n << " values for B2\n";

  return compute_b_squared(tmp1, tmp2, n, bsquared);
}

static int
gather_markers(const resizable_array_p<Predicted_Values>& zdata,
               IW_STL_Hash_Map_int& markers_present) {
  const auto n = zdata.number_elements();

  for (auto i = 0; i < n; ++i) {
    const IWString& s = zdata[i]->marker();

    auto f = markers_present.find(s);

    if (f == markers_present.end()) {
      markers_present[s] = 1;
    } else {
      (*f).second++;
    }
  }

  if (verbose) {
    cerr << "Across " << n << " records found " << markers_present.size()
         << " unique markers\n";
  }

  return markers_present.size();
}

static int
read_the_data(iwstring_data_source& input, int experimental_column,
              const resizable_array<int>& predicted_column, int& marker_column,
              const IWString& marker_column_name,
              resizable_array_p<Predicted_Values>& zdata) {
  input.set_dos(1);

  const_IWSubstring buffer;
  for (int i = 0; i < records_to_skip; i++) {
    input.next_record(buffer);
    if (1 == input.lines_read()) {
      if (!establish_column_titles(buffer, experimental_column, predicted_column,
                                   marker_column, marker_column_name)) {
        return 0;
      }
    }
  }

  int records_read = 0;
  while (input.next_record(buffer)) {
    if (1 == input.lines_read()) {
      if (!establish_column_titles(buffer, experimental_column, predicted_column,
                                   marker_column, marker_column_name)) {
        return 0;
      }
    }

    Predicted_Values* p = new Predicted_Values();

    int fatal;
    if (!p->parse_buffer(buffer, experimental_column, predicted_column, marker_column,
                         fatal)) {
      delete p;

      if (fatal) {
        cerr << "Fatal error processing '" << buffer << "', line " << input.lines_read()
             << '\n';
      }

      if (fatal) {
        return 0;
      } else {
        continue;
      }
    }

    zdata.add(p);

    records_read++;
  }

  if (verbose) {
    cerr << "Read " << records_read << " records\n";
  }

  return records_read;
}

static int
read_the_data(const char* fname, int experimental_column,
              const resizable_array<int>& predicted_column, int& marker_column,
              const IWString& marker_column_name,
              resizable_array_p<Predicted_Values>& zdata) {
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (0 == zdata.elements_allocated()) {
    zdata.resize(1000);  // just a guess
  } else {
    zdata.resize_keep_storage(0);
  }

  return read_the_data(input, experimental_column, predicted_column, marker_column,
                       marker_column_name, zdata);
}

static int
echo_the_data(const resizable_array_p<Predicted_Values>& zdata, std::ostream& output) {
  for (int i = 0; i < zdata.number_elements(); i++) {
    zdata[i]->debug_print(output);
  }

  return output.good();
}

static void
do_sort_by_predicted_value(resizable_array_p<Predicted_Values>& zdata,
                           int which_predicted_set, int predicted_column) {
  if (randomise_ties) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (auto* item : zdata) {
      item->assign_random_value(u(generator));
    }
  }

  Predicted_Value_Comparitor pvc(which_predicted_set);

  zdata.iwqsort(pvc);

  if (verbose > 2) {
    cerr << "After sorting on column " << predicted_column << '\n';
    echo_the_data(zdata, cerr);
  }

  return;
}

static void
establish_range_of_experimental_data(int number_records, int experimental_column,
                                     const resizable_array_p<Predicted_Values>& zdata,
                                     experimental_value_t& lower_experimental_value,
                                     experimental_value_t& upper_experimental_value) {
  assert(experimental_column >= 0);

  lower_experimental_value = upper_experimental_value = zdata[0]->obs();

  for (int i = 1; i < number_records; i++) {
    const Predicted_Values& pvi = *(zdata[i]);

    experimental_value_t o = pvi.obs();

    if (o < lower_experimental_value) {
      lower_experimental_value = o;
    } else if (o > upper_experimental_value) {
      upper_experimental_value = o;
    }
  }

  return;
}

static void
do_any_truncations_needed(experimental_value_t& pred,
                          const experimental_value_t lower_experimental_value,
                          const experimental_value_t upper_experimental_value) {
  if (pred > upper_experimental_value) {
    values_higher_than_experimental_range++;
    pred = upper_experimental_value;
  } else if (pred < lower_experimental_value) {
    values_lower_than_experimental_range++;
    pred = lower_experimental_value;
  }

  return;
}

static float
compute_observed_range(const resizable_array_p<Predicted_Values>& zdata,
                       int number_records) {
  experimental_value_t minval = zdata[0]->obs();
  experimental_value_t maxval = minval;

  for (int i = 1; i < number_records; i++) {
    experimental_value_t o = zdata[i]->obs();

    if (o < minval) {
      minval = o;
    } else if (o > maxval) {
      maxval = o;
    }
  }

  return maxval - minval;
}

static double
compute_relative_error(float obs, float pred) {
  if (obs == pred) {  // no error
    return 0.0;
  }

  if (fabs(obs - pred) < 1.0e-05) {  // numerically close enough
    return 0.0;
  }

  // both are very small - probably covered above
  if (fabs(obs) < 1.0e-05 && fabs(pred) < 1.0e-05) {
    return 0.0;
  }

  float rc = fabs((obs - pred) / obs);

  if (rc > relative_error_threshold) {
    return relative_error_threshold;
  }

  return rc;
}

static double
compute_range_normal_relative_error(experimental_value_t obs, experimental_value_t pred,
                                    double obs_range) {
  return fabs(obs - pred) / obs_range;
}

class CMSR_Parameters {
 public:
  float _v1;
  float _v2;
  float _v3;
  int _pct;
};

static int
write_cmsr(const CMSR_Parameters& cmsrp, const int which_predicted_set,
           const int predicted_column, const Accumulator<double>& e, const float* tmp1,
           const float* tmp2, float* tmpr, const int ndx, std::ostream& output) {
  // Note that we are using the tmp2 array, so we need to do this computation
  // before we call compute_b_squared which zorches that array.

  Accumulator<double> d;
  for (auto i = 0; i < ndx; ++i) {
    d.extra(tmp2[i]);
  }

  double bias = d.average();
  double sd = sqrt(d.variance());

  if (0.0 == sd) {
    cerr << "write_cmsr:variance is zero, cannot process\n";
    return 0;
  }

  bias = fabs(bias);
  double k = cmsrp._v1 + cmsrp._v2 * exp(-cmsrp._v3 * bias / sd);

// #define DEBUG_CMSR 20
#ifdef DEBUG_CMSR
  cerr << "Bias " << bias << " sd " << sd << ", k = " << k << '\n';
#endif

  float cmsr_sd = pow(10.0, bias + k * sd);

  if (which_predicted_set >= 0 && predicted_column >= 0) {
    write_something_identifying_the_column(predicted_column, output);
  }

  output << " cMSR" << cmsrp._pct << "_sd = " << cmsr_sd << '\n';

  std::copy(tmp2, tmp2 + ndx, tmpr);

#ifdef DEBUG_CMSR
  for (auto i = 0; i < DEBUG_CMSR; ++i) {
    cerr << "tmpr after copy from tmp2 " << tmpr[i] << '\n';
  }
#endif

  float median_difference;

  if (proper_median && ndx == (ndx / 2 * 2)) {
    std::partial_sort(tmpr, tmpr + ndx / 2 + 1, tmpr + ndx);
    median_difference = 0.50 * (tmpr[ndx / 2 - 1] + tmpr[ndx / 2]);
  } else {
    std::nth_element(tmpr, tmpr + ndx / 2 + 1, tmpr + ndx);
    median_difference = tmpr[ndx / 2];
  }
#ifdef DEBUG_CMSR
  cerr << "Median diff " << median_difference << '\n';
  for (auto i = 0; i < DEBUG_CMSR; ++i) {
    cerr << " partially sorted " << i << ' ' << tmpr[i] << '\n';
  }
#endif

  for (int i = 0; i < ndx; ++i) {
    tmpr[i] = fabs(tmp2[i] - median_difference);
  }

#ifdef DEBUG_CMSR
  for (auto i = 0; i < DEBUG_CMSR; ++i) {
    cerr << "abs diff from median " << tmpr[i] << '\n';
  }
#endif

  float mad;
  if (proper_median && ndx == (ndx / 2 * 2)) {
    std::partial_sort(tmpr, tmpr + ndx / 2 + 1, tmpr + ndx);
    mad = 0.50 * (tmpr[ndx / 2 - 1] + tmpr[ndx / 2]);
  } else {
    std::nth_element(tmpr, tmpr + ndx / 2, tmpr + ndx);
    mad = tmpr[ndx / 2];
  }

  mad = mad * 1.4826;

#ifdef DEBUG_CMSR
  cerr << "median diff " << tmpr[ndx / 2] << " MAD " << mad << '\n';
#endif

  if (which_predicted_set >= 0 && predicted_column >= 0) {
    write_something_identifying_the_column(predicted_column, output);
  }

  k = cmsrp._v1 + cmsrp._v2 * exp(-cmsrp._v3 * bias / mad);

  if (do_cMSD) {
    float cmsd = bias + k * mad;

    output << " cMSD" << cmsrp._pct << " = " << cmsd << '\n';
  } else {
    float cmsr_mad = pow(10.0, bias + k * mad);

    output << " cMSR" << cmsrp._pct << "_mad = " << cmsr_mad << '\n';
  }

  return 1;
}

static int
write_fold_data(const int predicted_column, const resizable_array<float>& fold_difference,
                const int* n_fold_difference, const int N, const char star_or_plus,
                std::ostream& output) {
  write_something_identifying_the_column(predicted_column, output);
  output << ' ' << star_or_plus << " fold";
  for (int i = 0; i < fold_difference.number_elements(); ++i) {
    float r = n_fold_difference[i] / static_cast<float>(N);

    output << ' ' << fold_difference[i] << ' ' << n_fold_difference[i] << ' ' << r;
  }

  output << '\n';

  return output.good();
}

static void
compute_fold_difference_values_additive(const float obs, const float pred,
                                        const resizable_array<float>& fold_difference,
                                        int* n_fold_difference) {
  float d = fabs(obs - pred);

  for (int i = 0; i < fold_difference.number_elements(); ++i) {
    if (d < fold_difference[i]) {
      n_fold_difference[i]++;
    }
  }

  return;
}

static void
compute_fold_difference_values_multiplicative(
    const float obs, const float pred, const resizable_array<float>& fold_difference,
    int* n_fold_difference) {
  if (fabs(obs) < 1.0e-04 || fabs(pred) < 1.0e-04) {
    return;
  }

  if (obs * pred < 0.0f) {
    return;
  }

  float ratio;
  if (obs < pred) {
    ratio = pred / obs;
  } else {
    ratio = obs / pred;
  }

  for (int i = 0; i < fold_difference.number_elements(); ++i) {
    if (ratio < fold_difference[i]) {
      n_fold_difference[i]++;
    }
  }

  return;
}

/*
  which_predicted_set is an index into the _predicted array for each prediction type
  predicted_column is the column number in the original file
*/

static int
iwstats(unsigned int number_records, const IWString* chunk_title, int which_predicted_set,
        int predicted_column, int experimental_column, float* tmp1, float* tmp2,
        float* tmpr, resizable_array_p<Predicted_Values>& zdata, std::ostream& output) {
  experimental_value_t lower_experimental_value = static_cast<experimental_value_t>(0.0);
  experimental_value_t upper_experimental_value = static_cast<experimental_value_t>(0.0);

  if (truncate_predicted_values_to_experimental_range) {
    values_lower_than_experimental_range = 0;
    values_higher_than_experimental_range = 0;

    establish_range_of_experimental_data(number_records, experimental_column, zdata,
                                         lower_experimental_value,
                                         upper_experimental_value);
  }

  if (predicted_column >= 0) {
    do_sort_by_predicted_value(zdata, which_predicted_set, predicted_column);
  }

#ifdef ECHO_THE_DATA
  cerr << "Processing " << number_records << " of " << which_predicted_set
       << " in column " << predicted_column << '\n';

  cerr << "N = " << zdata.number_elements() << '\n';
  for (auto i = 0; i < zdata.number_elements(); ++i) {
    cerr << zdata[i]->obs() << ' ' << zdata[i]->pred(0) << '\n';
  }
#endif

  double obs_range;

  if (range_normal_relative_error) {
    obs_range = compute_observed_range(zdata, number_records);
  } else {
    obs_range = 0.0;
  }

  // Do means and R2

  Accumulator<double> o, p;   // obs and pred
  Accumulator<double> e, ae;  // the error and absolute error
  Accumulator<double> are;    // average relative error
  Accumulator<double> rnre;   // range normalised average relative error
  double r = 0.0;
  int ndx = 0;

  int* n_fold_difference_additive = nullptr;
  if (fold_difference_additive.number_elements()) {
    n_fold_difference_additive = new_int(fold_difference_additive.number_elements());
  }

  std::unique_ptr<int[]> free_n_fold_difference_additive(n_fold_difference_additive);

  int* n_fold_difference_multiplicative = nullptr;
  if (fold_difference_multiplicative.number_elements()) {
    n_fold_difference_multiplicative =
        new_int(fold_difference_multiplicative.number_elements());
  }

  std::unique_ptr<int[]> free_n_fold_difference_multiplicative(
      n_fold_difference_multiplicative);

  for (unsigned int i = 0; i < number_records; i++) {
    const Predicted_Values& pvi = *(zdata[i]);

    if (!pvi.valid(which_predicted_set)) {
      continue;
    }

    experimental_value_t pred = pvi.pred(which_predicted_set);

    p.extra(pred);
    //  cerr << " i = " << i << " processing predicted value " << pred << '\n';

    if (truncate_predicted_values_to_experimental_range) {
      do_any_truncations_needed(pred, lower_experimental_value, upper_experimental_value);
    }

    if (experimental_column >= 0) {
      experimental_value_t obs = pvi.obs();

      if (obs < discard_measured_values_below) {
        continue;
      }

      o.extra(obs);

      r += obs * pred;

      e.extra(obs - pred);
      ae.extra(fabs(obs - pred));

      //    tmp1[ndx] = fabs(obs - pred);
      //    tmp2[ndx] =     (obs - pred);
      tmp1[ndx] = obs;
      tmp2[ndx] = pred;

      ndx++;

      if (fold_difference_additive.number_elements()) {
        compute_fold_difference_values_additive(obs, pred, fold_difference_additive,
                                                n_fold_difference_additive);
      }
      if (fold_difference_multiplicative.number_elements()) {
        compute_fold_difference_values_multiplicative(
            obs, pred, fold_difference_multiplicative, n_fold_difference_multiplicative);
      }

      if (relative_error_threshold > static_cast<float>(0.0)) {
        are.extra(compute_relative_error(obs, pred));
      }
      if (range_normal_relative_error) {
        rnre.extra(compute_range_normal_relative_error(obs, pred, obs_range));
      }
    }
  }

  if (0 == p.n()) {
    output << "no values for ";
    write_something_identifying_the_column(predicted_column, output);
    output << '\n';

    return 1;
  }

  if (p.n() < values_needed_for_reporting_correlations) {
    output << "Too few pairs " << p.n() << " for ";
    write_something_identifying_the_column(predicted_column, output);
    output << '\n';

    return 1;
  }

  // cerr << "Read " << p.n() << " predicted and " << o.n() << " obs values\n";

  if (verbose > 1 && p.n() < number_records) {
    cerr << "Only " << p.n() << " value pairs for ";
    write_something_identifying_the_column(predicted_column, cerr);
    cerr << '\n';
  }

  if (experimental_column >= 0) {
    output << "Experimental/Measured values in column " << (experimental_column + 1);
    if (activity_name.length() > 0) {
      output << ' ' << activity_name;
    }
    output << ". Range " << o.minval() << " and " << o.maxval();
    if (o.n() < 2) {
      output << " hmmm, only " << o.n() << " values available\n";
    } else {
      output << " ave " << o.average() << '\n';
    }
  }

  if (predicted_column >= 0) {
    output << "Predicted values in column " << (predicted_column + 1);
    write_something_identifying_the_column(predicted_column, output);
    output << ". Range " << p.minval() << " and " << p.maxval() << " ave " << p.average()
           << " N = " << p.n() << '\n';

    double obar = o.average();
    double pbar = p.average();

    if (p.minval() == p.maxval()) {
      cerr << "All predictions constant " << p.minval() << '\n';
      write_something_identifying_the_column(predicted_column, output);
      output << " R2 " << 0.0 << '\n';
      write_something_identifying_the_column(predicted_column, output);
      output << " Q2 " << 0.0 << '\n';
      output << " Average absolute error " << ae.average() << '\n';
      output << " Predicted Unbiased Bsquared 0\n";
      if (static_cast<float>(0.0) != relative_error_threshold) {
        output << " Average Relative Error 0\n";
      }
      if (range_normal_relative_error) {
        output << "Range Normal Relative Error 0\n";
      }
      return output.good();
    }

    double nx1bx2b = static_cast<double>(o.n()) * obar * pbar;

    //  cerr << "r = " << r << " nx1bx2b = " << nx1bx2b << '\n';

    double v1 = o.sum_of_squares() - o.n() * obar * obar;
    double v2 = p.sum_of_squares() - o.n() * pbar * pbar;

    if (0.0 == v1 || 0.0 == v2) {
      cerr << "Yipes, zero qzyq term " << v1 << " and v2 = " << v2 << '\n';
      return 0;
    }

    //  cerr << " v1 " << v1 << " v2 " << v2 << '\n';

    double rho = (r - nx1bx2b) / sqrt(v1 * v2);

    if (values_lower_than_experimental_range || values_higher_than_experimental_range) {
      output << values_lower_than_experimental_range << " values lower than expt.min, "
             << values_higher_than_experimental_range << " greater than max\n";
    }

    output << "Errors between " << e.minval() << " and " << e.maxval() << " ave "
           << e.average() << '\n';
    output << "Average absolute error " << ae.average() << '\n';
    output << "RMS error " << sqrt(ae.sum_of_squares() / static_cast<double>(ae.n()))
           << '\n';

    if (static_cast<float>(0.0) != relative_error_threshold) {
      output << "Average Relative Error " << static_cast<float>(are.average()) << '\n';
    }

    if (range_normal_relative_error) {
      output << "Average Range Normal Relative Error "
             << static_cast<float>(rnre.average()) << '\n';
    }

    write_something_identifying_the_column(predicted_column, output);
    output << " R2 ";

    if (wikipedia_r2_definition) {
      double n = 0.0;    // numerator
      double dno = 0.0;  // denominator for observed
      double dnp = 0.0;  // denominator for predicted
      for (int i = 0; i < ndx; ++i) {
        n += (tmp1[i] - obar) * (tmp2[i] - pbar);
        dno += (tmp1[i] - obar) * (tmp1[i] - obar);
        dnp += (tmp2[i] - pbar) * (tmp2[i] - pbar);
      }

      if (n < 0.0) {
        output << '-';
      }
      double r2 = n * n / (dno * dnp);
      output << static_cast<float>(r2) << '\n';
    } else {
      if (rho < 0.0) {
        output << '-';
      }
      output << static_cast<float>(rho * rho) << '\n';
    }

    //  Q squared

    double sum1 = 0.0;
    double sum2 = 0.0;
    for (unsigned int i = 0; i < number_records; i++) {
      const Predicted_Values& pvi = *(zdata[i]);

      if (!pvi.valid(which_predicted_set)) {
        continue;
      }

      experimental_value_t pred = pvi.pred(which_predicted_set);

      sum1 += (pvi.obs() - pred) * (pvi.obs() - pred);
      sum2 += (pvi.obs() - obar) * (pvi.obs() - obar);
    }

    double q2;
    if (traditional_q2_definition) {
      q2 = 1.0 - sum1 / sum2;
    } else if (sum1 > sum2) {
      q2 = -(1.0 - sum2 / sum1);
    } else {
      q2 = 1.0 - sum1 / sum2;
    }

    write_something_identifying_the_column(predicted_column, output);

    if (traditional_q2_definition) {
      output << " trQ2 ";
    } else {
      output << " Q2 ";
    }
    output << static_cast<float>(q2) << '\n';

    write_something_identifying_the_column(predicted_column, output);
    output << " Bias " << (o.average() - p.average()) << '\n';

    sum1 = 0.0;
    for (unsigned int i = 0; i < number_records; i++) {
      const Predicted_Values& pvi = *(zdata[i]);

      if (!pvi.valid(which_predicted_set)) {
        continue;
      }

      experimental_value_t obs = pvi.obs();
      experimental_value_t pred = pvi.pred(which_predicted_set);

      sum1 += (obs - obar) * (pred - pbar);
    }

    //  write_something_identifying_the_column (predicted_column, output);
    //  output << " R2 " << (sum1 / (p.variance() * o.variance()) / static_cast<double>
    //  (number_records)) << '\n';
  }

  if (fold_difference_additive.number_elements()) {
    write_fold_data(predicted_column, fold_difference_additive,
                    n_fold_difference_additive, number_records, '+', output);
  }

  if (fold_difference_multiplicative.number_elements()) {
    write_fold_data(predicted_column, fold_difference_multiplicative,
                    n_fold_difference_multiplicative, number_records, '*', output);
  }

  int duplicate_predicted_values_present;
  if (predicted_column >= 0) {
    duplicate_predicted_values_present =
        duplicate_predicted_values(zdata, which_predicted_set);
  } else {
    duplicate_predicted_values_present = 0;
  }

  if (0 == verbose || predicted_column < 0) {  // don't write anything
    ;
  } else if (0 == duplicate_predicted_values_present) {
    cerr << "No duplicate predicted values in column " << (predicted_column + 1) << '\n';
  } else {
    cerr << "Duplicate prediced values present in column " << (predicted_column + 1)
         << '\n';
  }

  if (nullptr != chunk_title) {  // we are dealing with a subset
    output << " first " << *(chunk_title) << " (" << number_records << " records)\n";
  }

  for (int i = 0; i < ndx; ++i) {
    float o = tmp1[i];
    float p = tmp2[i];
    tmp1[i] = fabs(o - p);
    tmp2[i] = (o - p);
  }

  Float_Comparator fc;
  iwqsort(tmp1, ndx, fc);

  write_something_identifying_the_column(predicted_column, output);
  output << " AE50 = " << tmp1[ndx / 2] << "\n";

  write_something_identifying_the_column(predicted_column, output);
  output << " AE75 = " << tmp1[ndx * 75 / 100] << "\n";

  write_something_identifying_the_column(predicted_column, output);
  output << " AE95 = " << tmp1[ndx * 95 / 100] << "\n";

#ifdef DEBUG_CMSR
  for (auto i = 0; i < DEBUG_CMSR; ++i) {
    cerr << "tmp2 " << tmp2[i] << '\n';
  }
#endif

  if (stream_for_residuals.is_open()) {
    write_something_identifying_the_column(predicted_column, stream_for_residuals);
    stream_for_residuals << '\n';
    for (int i = 0; i < ndx; ++i) {
      stream_for_residuals << tmp2[i] << '\n';
      stream_for_residuals.write_if_buffer_holds_more_than(4096);
    }
  }

  CMSR_Parameters cmsrp;

  cmsrp._v1 = 1.645;
  cmsrp._v2 = 0.325;
  cmsrp._v3 = 4.134;
  cmsrp._pct = 95;

  write_cmsr(cmsrp, which_predicted_set, predicted_column, e, tmp1, tmp2, tmpr, ndx,
             output);

  cmsrp._v1 = 1.274;
  cmsrp._v2 = 0.382;
  cmsrp._v3 = 3.448;
  cmsrp._pct = 90;

  write_cmsr(cmsrp, which_predicted_set, predicted_column, e, tmp1, tmp2, tmpr, ndx,
             output);

  if (which_predicted_set >= 0 && predicted_column >= 0) {
    write_something_identifying_the_column(predicted_column, output);
  }

  if (duplicate_predicted_values_present) {
    output << " Input";
  } else {
    output << " Unbiased";
  }

  double Bsquared;
  compute_b_squared(number_records, experimental_column, zdata, tmp1, tmp2, Bsquared,
                    which_predicted_set);

  output << " Bsquared " << static_cast<float>(Bsquared) << '\n';

  if (calculate_enrichment_metrics &&
      initilize_enrichment(number_records, experimental_column, zdata, tmp1,
                           which_predicted_set)) {
    // Calculate the BEDROC metric
    if (which_predicted_set >= 0 && predicted_column >= 0) {
      write_something_identifying_the_column(predicted_column, output);
    }

    calculateBEDROC bedroc(bedroc_alpha);
    double bedroc_value = bedroc(enrichment);
    output << " BEDROC " << static_cast<float>(bedroc_value) << '\n';

    calculateROC roc;

    double roc_value = roc(enrichment);

    if (which_predicted_set >= 0 && predicted_column >= 0) {
      write_something_identifying_the_column(predicted_column, output);
    }

    output << " ROC " << static_cast<float>(roc_value) << '\n';

    calculateEF ef(enrichment_factor_fraction);
    double ef_value = ef(enrichment);

    if (which_predicted_set >= 0 && predicted_column >= 0) {
      write_something_identifying_the_column(predicted_column, output);
    }

    output << " EF " << static_cast<float>(ef_value) << '\n';
  }

  // To avoid horrible problems with iwqsort, which really doesn't do well sorting
  // a nearly sorted list, we sort by a random value before each sort!

  if (which_predicted_set >= 0 && duplicate_predicted_values_present &&
      compute_range_of_bsquared_values_when_duplicates_present) {
    output << "Contains duplicate predicted values in column " << (predicted_column + 1)
           << ", computing best and worst Bsquared\n";

    std::random_shuffle(zdata.begin(), zdata.end());  // randomise before sorting

    Predicted_Value_Comparitor_Best pvcb(which_predicted_set);

    zdata.iwqsort(pvcb);

    double best_bsquared;
    compute_b_squared(number_records, experimental_column, zdata, tmp1, tmp2,
                      best_bsquared, which_predicted_set);
    write_something_identifying_the_column(predicted_column, output);
    output << " Best Bsquared " << best_bsquared << '\n';

    std::random_shuffle(zdata.begin(), zdata.end());  // randomise before sorting

    Predicted_Value_Comparitor_Worst pvcw(which_predicted_set);
    zdata.iwqsort(pvcw);

    double worst_bsquared;
    compute_b_squared(number_records, experimental_column, zdata, tmp1, tmp2,
                      worst_bsquared, which_predicted_set);
    write_something_identifying_the_column(predicted_column, output);
    output << " Worst Bsquared " << worst_bsquared << '\n';

    float b2 = static_cast<float>((best_bsquared + worst_bsquared) *
                                  0.5);  // cut precision to suppress significant figures

    write_something_identifying_the_column(predicted_column, output);
    output << " Unbiased Bsquared " << b2 << '\n';
  }

  if (number_distribution_buckets > 0) {
    do_compute_distribution_functions(number_records, which_predicted_set,
                                      predicted_column, experimental_column, tmp1, tmp2,
                                      zdata, output);
  }

  return output.good();
}

static int
iwstats(int which_predicted_set, int predicted_column, int experimental_column,
        float* tmp1, float* tmp2, float* tmpr,
        const resizable_array<int>& numbers_to_check,
        resizable_array_p<Predicted_Values>& zdata, std::ostream& output) {
  if (numbers_to_check.empty()) {  // no subsetting, doing the whole file at once
    return iwstats(zdata.number_elements(), nullptr, which_predicted_set,
                   predicted_column, experimental_column, tmp1, tmp2, tmpr, zdata,
                   output);
  }

  for (int i = 0; i < numbers_to_check.number_elements(); i++) {
    int n = numbers_to_check[i];

    if (!iwstats(n, chunk_titles[i], which_predicted_set, predicted_column,
                 experimental_column, tmp1, tmp2, tmpr, zdata, output)) {
      return 0;
    }
  }

  return output.good();
}

static int
iwstats(const resizable_array<int>& predicted_column, int experimental_column,
        resizable_array_p<Predicted_Values>& zdata, std::ostream& output) {
  int records_in_file = zdata.number_elements();

  assert(records_in_file > 1);

  // Do we want data on the first n or the first p% of the file

  resizable_array<int> numbers_to_check;

  if (step || percents_to_check.number_elements() || N_to_check.number_elements()) {
    determine_numbers_to_check(numbers_to_check, records_in_file);
  } else if (activities_to_check.number_elements() && experimental_column >= 0) {
    determine_numbers_to_check_by_activity(zdata, numbers_to_check);
  }

  // cerr << "Will check " << numbers_to_check.number_elements() << " items\n";

  float* tmp1 = new float[records_in_file];
  std::unique_ptr<float[]> free_tmp1(tmp1);
  float* tmp2 = new float[records_in_file];
  std::unique_ptr<float[]> free_tmp2(tmp2);
  float* tmpr = new float[records_in_file];
  std::unique_ptr<float[]> free_tmpr(tmpr);

  if (nullptr == tmp1 || nullptr == tmp2 || nullptr == tmpr) {
    cerr << "Memory failure, cannot allocate data for " << records_in_file
         << " records\n";
    return 9;
  }

  for (int i = 0; i < predicted_column.number_elements(); i++) {
    int col = predicted_column[i];

    if (!iwstats(i, col, experimental_column, tmp1, tmp2, tmpr, numbers_to_check, zdata,
                 output)) {
      cerr << "Fatal error processing column " << (col + 1) << '\n';
      return 0;
    }
  }

  return output.good();
}

static int
iwstats(const IW_STL_Hash_Map_int::const_iterator& f,
        const resizable_array<int>& predicted_column, int experimental_column,
        resizable_array_p<Predicted_Values>& zdata, std::ostream& output) {
  resizable_array_p<Predicted_Values> tmp;

  const auto n = zdata.number_elements();

  tmp.resize(n);

  for (auto i = 0; i < n; ++i) {
    if ((*f).first == zdata[i]->marker()) {
      tmp.add(zdata[i]);
    }
  }

  if (tmp.empty()) {
    cerr << "Very strange, no records retrieved for marker '" << (*f).first << "'\n";
    return 0;
  } else if (tmp.number_elements() < 3) {
    tmp.resize_no_delete(0);
    cerr << "Only one instance of '" << (*f).first << "', ignored\n";
    return 1;
  }

  global_prefix.resize_keep_storage(0);

  global_prefix << ' ' << (*f).first;

  auto rc = iwstats(predicted_column, experimental_column, tmp, output);

  global_prefix.resize_keep_storage(0);

  tmp.resize_no_delete(0);

  return rc;
}

static int
read_id_activity_hash_record(const const_IWSubstring& buffer,
                             IW_STL_Hash_Map<IWString, IWString>& id_activity_hash,
                             int experimental_column) {
  int i = 0;
  IWString id;

  if (!buffer.nextword(id, i, input_separator)) {
    cerr << "Cannot extract identifier from input '" << buffer << "'\n";
    return 0;
  }

  if (strip_leading_zeros_from_identifiers) {
    id.remove_leading_chars('0');
  }

  IWString act;

  if (!buffer.nextword(act, i, input_separator)) {
    cerr << "Cannot extract second token from record '" << buffer << "'\n";
    return 0;
  }

  if (experimental_column > 1) {
    for (auto col = 1; col < experimental_column; ++col) {
      if (!buffer.nextword(act, i, input_separator)) {
        cerr << "Cannot find experimental data column '" << buffer << "'\n";
        return 0;
      }
    }
  }

  if (0 == id_activity_hash.size()) {
    activity_name = act;
  }

  if ("." == act) {
    cerr << "Missing value encountered in experimental file '" << buffer
         << "', impossible\n";
    return 0;
  }

  double notused;

  if (act.numeric_value(notused)) {  // great
    ;
  } else if (0 == id_activity_hash.size()) {  // header
    ;
  } else {
    cerr << "Invalid numeric value in experimental file '" << buffer << "'\n";
    return 0;
  }

  id_activity_hash[id] = act;

  return 1;
}

static int
read_id_activity_hash(iwstring_data_source& input,
                      IW_STL_Hash_Map<IWString, IWString>& id_activity_hash,
                      int experimental_column) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    buffer.strip_trailing_blanks();

    if (!read_id_activity_hash_record(buffer, id_activity_hash, experimental_column)) {
      cerr << "Invalid id/activity record, line " << input.lines_read() << " '" << buffer
           << "'\n";
      return 0;
    }
  }

  return id_activity_hash.size();
}

static int
read_id_activity_hash(const const_IWSubstring& fname,
                      IW_STL_Hash_Map<IWString, IWString>& id_activity_hash,
                      int experimental_column) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_translate_tabs(1);

  return read_id_activity_hash(input, id_activity_hash, experimental_column);
}

template <typename T>
int
get_comma_separated_values(Command_Line& cl, const char flag, resizable_array<T>& v,
                           const float minval) {
  const_IWSubstring f;
  for (int i = 0; cl.value(flag, f, i); ++i) {
    int j = 0;
    const_IWSubstring token;
    while (f.nextword(token, j, ',')) {
      float x;
      if (!token.numeric_value(x) || x <= minval) {
        cerr << "The -" << flag << " must be greater than " << minval << " '" << token
             << "' invalid\n";
        return 0;
      }

      v.add(x);

      if (verbose) {
        cerr << "Option " << flag << " got value " << x << '\n';
      }
    }
  }

  return v.number_elements();
}

/*
  Wow, this is really hard because at this point we only have the activity
  as a string. Fix this sometime....
*/

static int
identify_most_active_fraction(const IW_STL_Hash_Map<IWString, IWString>& activity,
                              const double fraction_to_keep,
                              double& discard_measured_values_below) {
  const auto s = activity.size();

  float* v = new float[s];
  std::unique_ptr<float[]> free_v(v);

  int ndx = 0;
  for (auto f : activity) {
    float a;
    if (!f.second.numeric_value(a)) {
      continue;
    }

    v[ndx] = a;
    ndx++;
  }

  std::sort(v, v + ndx);  // ascending. Could do a partial sort if we wanted to

  const int i = static_cast<int>((1.0 - fraction_to_keep) * ndx + 0.4999);

  discard_measured_values_below = v[i];

  if (verbose) {
    cerr << "Keeping the " << fraction_to_keep
         << " most active values, implies an activity cutoff of "
         << discard_measured_values_below << " (range " << v[0] << " to " << v[ndx - 1]
         << ")\n";
  }

  return 1;
}

int
parse_digits(const const_IWSubstring& token, int istart, int& zresult) {
  zresult = 0;

  int rc = 0;  // the number of characters we consume

  for (int i = istart; i < token.length(); i++) {
    int j = token[i] - '0';

    if (j > 9) {
      return rc;
    }

    if (j < 0) {
      return rc;
    }

    zresult = zresult * 10 + j;
    rc++;
  }

  return rc;
}

static int
parse_dash_p(const const_IWSubstring& p, resizable_array<int>& predicted_column) {
  int i = 0;
  const_IWSubstring token;
  while (p.nextword(token, i, ',')) {
    int c1;

    int chars_consumed = parse_digits(token, 0, c1);

    if (chars_consumed == token.length()) {
      predicted_column.add(c1 - 1);
      continue;
    }

    if ('-' != token[chars_consumed]) {
      cerr << "Ranges must be of the form 'n-m'\n";
      return 0;
    }

    int c2;

    int chars_consumed_2 = parse_digits(token, chars_consumed + 1, c2);

    if (chars_consumed + 1 + chars_consumed_2 != token.length()) {
      cerr << "Invalid range specifier '" << token << "'\n";
      usage(6);
    }

    for (int j = c1; j <= c2; j++) {
      predicted_column.add(j - 1);
    }
  }

  return 1;
}

static int
iwstats(int argc, char** argv) {
  Command_Line cl(argc, argv, "ve:E:p:s:jwn:t:P:M:qR:zc:TA:b:hr:ko:a:f:u:U:i:Km:dF:D:L:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(3);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('i')) {
    IWString i = cl.string_value('i');
    if (!char_name_to_char(i)) {
      cerr << "Unrecognised input separator specification '" << i << "'\n";
      return 2;
    }

    input_separator = i[0];
  }

  if (cl.option_present('j')) {
    is_descriptor_file = 1;

    if (verbose) {
      cerr << "Will process as a descriptor file\n";
    }
  }

  int experimental_column = -1;

  if (cl.option_present('e')) {
    if (!cl.value('e', experimental_column) || experimental_column < 1) {
      cerr << "The experimental/measured value column (-e) must be a whole positive "
              "number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Experimental/Measured values in column " << experimental_column << '\n';
    }

    experimental_column--;
  }

  if (cl.option_present('E')) {
    const_IWSubstring e = cl.string_value('E');

    if (cl.option_present('z')) {
      strip_leading_zeros_from_identifiers = 1;

      if (verbose) {
        cerr << "Will strip leading 0's from identifiers\n";
      }
    }

    if (!read_id_activity_hash(e, id_activity_hash, experimental_column) ||
        0 == id_activity_hash.size()) {
      cerr << "Cannot read id/activity hash from '" << e << "'\n";
      return 5;
    }

    if (verbose) {
      cerr << "Read " << id_activity_hash.size() << " activity values from '" << e
           << "'\n";
    }

    is_descriptor_file = 1;
    experimental_column = 0;  // activity data has come from a separate file, don't go
                              // looking for it in the input file

    if (!cl.option_present('p')) {
      cerr << "Must specify one or more predicted columns via the -p option\n";
      usage(5);
    }
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', values_needed_for_reporting_correlations) ||
        values_needed_for_reporting_correlations < 3) {
      cerr << "The number of pairs needed for correlation compuations (-c) must be 3 or "
              "larger\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Correlations only done when at least "
           << values_needed_for_reporting_correlations << " data points available\n";
    }
  }

  resizable_array<int> predicted_column;

  if (cl.option_present('p')) {
    int i = 0;
    const_IWSubstring p;
    while (cl.value('p', p, i++)) {
      if (!parse_dash_p(p, predicted_column)) {
        cerr << "Invalid predicted column(s) specifier '" << p << "'\n";
        usage(6);
      }
    }

    if (predicted_column.contains(experimental_column)) {
      cerr << "Column " << experimental_column
           << " cannot be both predicted and experimental\n";
      return 8;
    }
  }

  if (experimental_column < 0 && predicted_column.empty()) {
    if (is_descriptor_file) {
      experimental_column = -1;
    } else {
      experimental_column = 0;
    }
  }

  if (verbose) {
    if (experimental_column >= 0) {
      cerr << "Measured/Observed values in column " << (experimental_column + 1) << '\n';
    }

    if (predicted_column.number_elements()) {
      cerr << "Predicted values in columns";
      for (int i = 0; i < predicted_column.number_elements(); i++) {
        cerr << ' ' << (predicted_column[i] + 1);
      }
      cerr << '\n';
    }
  }

  if (predicted_column.empty()) {
    predicted_column.add(-1);
  }

  if (cl.option_present('h')) {
    traditional_q2_definition = 1;
    if (verbose) {
      cerr << "Will report traditional Q2 definitions\n";
    }
  }

  if (cl.option_present('q') && !cl.option_present('M')) {
    cerr << "Sorry, the no-report missing values option (-q) only makes sense with the "
            "-M option\n";
    usage(2);
  }

  if (cl.option_present('M')) {
    missing_value = cl.string_value('M');

    if (verbose) {
      cerr << "Missing value set to '" << missing_value << "'\n";
    }

    if (cl.option_present('q')) {
      report_missing_values = 0;

      if (verbose) {
        cerr << "Will quietly ignore missing values\n";
      }
    }
  }

  if (cl.option_present('r')) {
    if (!cl.value('r', relative_error_threshold) ||
        relative_error_threshold <= static_cast<float>(0.0)) {
      cerr << "The max relative error (-r) option must be a non zero, positive floating "
              "point value\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Relative error computation threshold " << relative_error_threshold << '\n';
    }

    range_normal_relative_error = 1;  // turn on by default
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', records_to_skip) || records_to_skip < 1) {
      cerr << "The number of records to skip (-s option) must be a whole positive "
              "number\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Will skip the first " << records_to_skip << " records\n";
    }
  }

  if (is_descriptor_file) {
    if (records_to_skip < 1) {
      records_to_skip = 1;
    }

    if (0 == experimental_column && !cl.option_present('E')) {
      cerr << "Descriptor file, but you gave the experimental column as 1, sorry...\n";
      usage(7);
    }

    if (predicted_column.contains(0)) {
      cerr << "Descriptor file but you said predicted values in column 1. Sorry\n";
      usage(8);
    }
  }

  if (cl.option_present('k')) {
    skip_items_with_no_activity_data = cl.option_count('k');
    if (verbose) {
      cerr << "Will skip items where there is no experimental data\n";
    }
  }

  if (cl.option_present('K')) {
    wikipedia_r2_definition = 0;

    if (verbose) {
      cerr << "Will use R2 definition\n";
    }
  }

  if (cl.option_present('w')) {
    compute_range_of_bsquared_values_when_duplicates_present = 0;

    if (verbose) {
      cerr << "Will suppress computation of best and worst BSquared values\n";
    }
  }

  int n = (cl.option_present('n') > 0) + (cl.option_present('t') > 0) +
          (cl.option_present('P') > 0);

  if (n > 1) {
    cerr << "Hmmm, you have specified more than one of: a constant step (-t) percent "
            "sample (-P) or number (-n)\n";
    cerr << "options. Might be a little strange....\n";
    //  cerr << cl.option_present('n') << ' ' << cl.option_present('t') << ' ' <<
    //  cl.option_present('P') << ' ' << n << '\n';
  }

  if (cl.option_present('n')) {
    int i = 0;
    const_IWSubstring n;
    while (cl.value('n', n, i++)) {
      int j;
      if (!n.numeric_value(j) || j < 2) {
        cerr << "Invalid number of records to process '" << n << "'\n";
        usage(4);
      }

      N_to_check.add_if_not_already_present(j);

      if (verbose) {
        cerr << "Will examine highest predicted " << j << " records in the file\n";
      }
    }
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', step) || step < 1) {
      cerr << "The step every option (-t) must be a whole positive number\n";
      usage(11);
    }

    if (verbose) {
      cerr << "Will compute Bsquared values for " << step
           << " sized sections of the data\n";
    }
  }

  if (cl.option_present('P')) {
    int i = 0;
    int p;
    while (cl.value('P', p, i++)) {
      if (p <= 0 || p > 100) {
        cerr << "Percentages must be between 1 and 100, '" << p << "' is invalid\n";
        usage(7);
      }

      if (verbose) {
        cerr << "Will report the first " << p << " percent of the data\n";
      }

      percents_to_check.add_if_not_already_present(p);
    }
  }

  if (cl.option_present('A')) {
    int i = 0;
    float a;
    while (cl.value('A', a, i++)) {
      if (verbose) {
        cerr << "Will report activities above " << a << "\n";
      }

      activities_to_check.add_if_not_already_present(a);
    }
  }

  if (cl.option_present('T')) {
    truncate_predicted_values_to_experimental_range = 1;

    if (verbose) {
      cerr << "Predicted values truncated to experimental range\n";
    }
  }

  if (cl.option_present('b')) {
    if (!cl.value('b', number_distribution_buckets) || number_distribution_buckets < 2) {
      cerr << "The number of buckets for distribution studies (-b) must be a whole +ve "
              "number\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Will report distribution in " << number_distribution_buckets
           << " buckets\n";
    }
  }

  if (cl.option_present('o')) {
    if (!cl.value('o', active_inactive_cutoff)) {
      cerr << "Invalid cutoff for active/inactive \n";
      usage(5);
    }

    calculate_enrichment_metrics = 1;

    if (verbose) {
      cerr << "Will calculate enrichment metrics, such as BEDROC and EF \n";
    }

    if (cl.option_present('a')) {
      if (!cl.value('a', bedroc_alpha)) {
        cerr << "Invalid value for BEDROC alpha \n";
        usage(5);
      }
    }

    if (cl.option_present('f')) {
      if (!cl.value('f', enrichment_factor_fraction) ||
          enrichment_factor_fraction <= 0.0 || enrichment_factor_fraction >= 1.0) {
        cerr << "The enrichment factor fraction (-f) must be a valid fraction\n";
        usage(2);
      }
    }
  }

  if (cl.option_present('u') && cl.option_present('U')) {
    cerr << "Cannot use both the -u and -U options\n";
    usage(1);
  }

  if (cl.option_present('u')) {
    if (!cl.value('u', discard_measured_values_below)) {
      cerr << "The discard measured values below option (-u) must be a valid float\n";
      usage(2);
    }
    if (verbose) {
      cerr << "Will discard measured values below " << discard_measured_values_below
           << '\n';
    }
  }

  /*
    We use the -u mechanism in order to implement the -U functionality. Just identify the
    cut point of activity that holds the -U fraction of molecules
  */

  if (cl.option_present('U')) {
    double u;
    if (!cl.value('U', u) || u <= 0.0 || u >= 1.0) {
      cerr << "The fraction of most active to retain (-U) option must be a valid "
              "fraction\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will only keep the most " << u << " active molecules\n";
    }

    identify_most_active_fraction(id_activity_hash, u, discard_measured_values_below);
  }

  if (cl.option_present('R')) {
    const_IWSubstring r = cl.string_value('R');

    if ('.' == r) {
    } else {
      cerr << "Setting random seed no longer supported, ignored\n";
    }

    randomise_ties = 1;

    if (verbose) {
      cerr << "Randomise tied predicted values\n";
    }
  }

  if (cl.option_present('F')) {
    if (!get_comma_separated_values(cl, 'F', fold_difference_multiplicative, 1.0f)) {
      cerr << "Cannot parse multiplicative -Nfold option (-F)\n";
      return 1;
    }
  }

  if (cl.option_present('D')) {
    if (!get_comma_separated_values(cl, 'D', fold_difference_additive, 0.0f)) {
      cerr << "Cannot parse additive -Nfold option (-D)\n";
      return 1;
    }
  }

  if (experimental_column < 0 && predicted_column.number_elements() >= 0) {
    cerr << "No experimental column specified, cannot proceed\n";
    usage(4);
  }

  int marker_column = -1;
  IWString marker_column_name;

  if (cl.option_present('m')) {
    IWString m = cl.string_value('m');

    if (m.numeric_value(marker_column)) {
      if (marker_column < 1) {
        cerr << "The marker column option must be a valid column number\n";
        usage(1);
      }

      if (verbose) {
        cerr << "Marker column set to " << marker_column << '\n';
      }

      marker_column--;
    } else {
      marker_column_name = m;
      is_descriptor_file = 1;
      if (records_to_skip <= 0) {
        records_to_skip = 1;
      }

      if (verbose) {
        cerr << "The marker column will be descriptor '" << marker_column_name << "'\n";
      }
    }
  }

  if (cl.option_present('d')) {
    do_cMSD = 1;

    if (verbose) {
      cerr << "Compute cMSD rather than cMSR\n";
    }
  }

  if (cl.empty()) {
    cerr << "INsufficient arguments\n";
    usage(1);
  }

  if (cl.option_present('L')) {
    const char* fname = cl.option_value('L');

    if (!stream_for_residuals.open(fname)) {
      cerr << "Cannot open stream for residuals '" << fname << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Residuals written to '" << fname << "'\n";
    }
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    resizable_array_p<Predicted_Values> zdata(2000);

    if (!read_the_data(cl[i], experimental_column, predicted_column, marker_column,
                       marker_column_name, zdata)) {
      cerr << "Error reading data file '" << cl[i] << "'\n";
      rc = i + 1;
      break;
    }

    if (verbose) {
      cerr << "Read " << zdata.number_elements() << " data items from '" << cl[i]
           << "'\n";
    }

    if (verbose > 2) {
      echo_the_data(zdata, cerr);
    }

    if (1 == zdata.number_elements()) {
      cerr << "Very strange, only one record in '" << cl[i] << "', ignored\n";
      continue;
    }

    if (!iwstats(predicted_column, experimental_column, zdata, cout)) {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      rc = i + 1;
      break;
    }

    if (marker_column >= 0) {
      IW_STL_Hash_Map_int markers_present;
      gather_markers(zdata, markers_present);

      if (verbose) {
        cerr << "Found " << markers_present.size() << " different marker values\n";
      }

      for (auto j = markers_present.begin(); j != markers_present.end(); ++j) {
        if ((*j).second < 3) {
          continue;
        }

        cout << "By " << (*j).first << "\n";

        if (!iwstats(j, predicted_column, experimental_column, zdata, cout)) {
          cerr << "Fatal error processing marker '" << (*j).first << " in '" << cl[i]
               << "'\n";
          rc = i + 1;
          break;
        }
      }
    }
  }

  if (missing_values_encountered) {
    cerr << "Encountered " << missing_values_encountered << " missing values\n";
  }

  if (items_with_no_activity_data) {
    cerr << "Encountered " << items_with_no_activity_data
         << " items with no activity data\n";
  }

  if (nullptr != column_titles) {
    delete[] column_titles;
  }

  return rc;
}

int
main(int argc, char** argv) {
  int rc = iwstats(argc, argv);

  return rc;
}
