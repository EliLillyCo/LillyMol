/*
  Remove correlated descriptors
  This variant computes actual correlation coefficients between columns
*/

#include <math.h>
#include <stdlib.h>

#include <fstream>
#include <limits>
#include <memory>
#include <random>
#include <vector>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iwstring.h"

using std::cerr;

const char* prog_name = NULL;

static int verbose = 0;

static Report_Progress report_progress;

static int write_all_correlated_pairs = 0;

static int print_column_statistics = 0;

static std::ofstream stream_for_column_statistics;

static int sample_records_randomly = 1;

static int default_records_to_sample = 2000;

/*
  Nov 2000. Dave wanted to be able to just ignore columns with invalid data....
*/

static int ignore_columns_with_non_numeric_data = 0;

/*
  Dec 99. Dave wanted the programme to collect all values. Previously
  I stopped accumulating things once we crossed continuous_threshold
*/

static int collect_all_values = 1;

/*
  When collapsing near-binary descriptors down to binary forms, we
  need a threshold
*/

static float binary_collapse_threshold = 0.0;

/*
  Columns with more than continuous_threshold different values are
  considered to be continuous and will never be compared with columns
  which are "binary";
*/

static unsigned int continuous_threshold = 5;

/*
  This will be either 0, 1 or greater than 1
  1 means suppress the column if all values are constant AND THERE ARE
  NO MISSING VALUES
  A value greater than 1 means suppress the column if all values are
  either a constant or missing
*/

static int suppress_constant_columns = 0;

/*
  Sometimes we want to trim down until we have a certain number of columns
*/

static int columns_to_keep = 0;

/*
  October 2010, programme was not dealing with strongly negatively
  correlated values. Make that a default behaviour
*/

static int take_absolute_value_of_correlation_coefficient = 1;

static char input_delimitor(' ');

/*
  Stuff needed for the [float] int hash map
*/

struct eqflt {
  bool
  operator()(const float& d1, const float& d2) const
  {
    return d1 == d2;
  }
};

struct float_hash : public std::hash<float> {
  size_t
  operator()(const float& d) const
  {
    const int* i = (const int*)&d;

    return *i;
  }
};

/*
  Each column keeps a list of other columns that are correlated, together with the
  correlation coef
*/

typedef float correlation_coeff_t;

class Correlated_Column : public std::pair<int, correlation_coeff_t>
{
 private:
 public:
  Correlated_Column();
  Correlated_Column(int, correlation_coeff_t);

  int
  column() const
  {
    return first;
  }

  correlation_coeff_t
  r2() const
  {
    return second;
  }
};

Correlated_Column::Correlated_Column()
{
  first = -1;
  second = 0.0;

  return;
}

Correlated_Column::Correlated_Column(int c, correlation_coeff_t r)
{
  first = c;
  second = r;
}

#ifdef __GNUG__
// template class resizable_array_p<Correlated_Column>;
// template class resizable_array_base<Correlated_Column *>;
#endif

/*
  There are lots of things we need to keep track of about each column

  When deciding which columns to discard we keep track of the number of
  columns correlated to each column above the threshold. That is the
  resizable_array_p<Correlated_Column> part
*/

typedef double stored_value_t;

class SColumn : public resizable_array_p<Correlated_Column>, public Accumulator<double>
{
 private:
  int _column_number;

  //  When deciding whether a column is binary, near binary or continuous,
  //  we need a count of the number of different values present. Note that
  //  we use a float hash deliberately, so close values may be coalesced

  //  IW_Hash_Map<float, int, std::hash<float>, std::equal_to<float>> _different_values;
  //  hash_map<float, int, float_hash, eqflt> _different_values;
  std::unordered_map<float, int, float_hash, eqflt> _different_values;

  int _ignore;
  int _suppress;
  int _missing_value_count;

  double _min_rho;
  int _min_rho_column;
  double _max_rho;
  int _max_rho_column;

  resizable_array<stored_value_t> _values_sampled;

  //  If this is a descriptor file, each column will have a title

  IWString _title;

  //  If we are a near-binary descriptor, we may have some collapsed values

  int _collapsed_values;

  //  But when we collapse values, we need to know the initial number of
  //  different values

  int _initial_different_values;

  //  We decided that when comparing columns for who has the most correlated
  //  columns, we would give extra weight to columns which have been collapsed.
  //  So, if two columns each have 10 other columns correlated, and one was
  //  binary and the other near binary, the near binary column is the more desirable

  int _adjusted_correlated_columns;

  //  When adaptively trimming, we need a means of temporarily storing corrlated columns
  //  that have been eliminated.

  resizable_array<Correlated_Column*> _storage;

  //  Can speed up calculations with some precomputed values

  double _sumx;
  double _sumx2;

  //  private functions

  void
  _do_Q_test(SColumn&);

  void
  _compute_sum_and_sum2();

  void
  _pearsons_test(SColumn&);

  void
  _collapse_values();

  void
  _do_collapse(float f1, float f2, float f, int n);
  void
  _do_collapse(float fto, float f, int n);

  void
  _write_correlation(std::ostream&, const SColumn&, const char*, double) const;

  void
  _compute_adjusted_correlated_columns();

 public:
  SColumn();

  int
  column_number() const
  {
    return _column_number;
  }

  void
  set_column_number(int c)
  {
    _column_number = c;
  }

  int
  write(std::ostream&) const;
  int
  print_column_statistics(std::ostream& os) const;

  int
  missing_value_count() const
  {
    return _missing_value_count;
  }

  int
  write_min_max_rho(std::ostream&) const;

  int
  set_sample_size(int s);

  void
  column_is_correlated(int c, correlation_coeff_t r);

  void
  column_is_suppressed(int col);
  void
  column_is_possibly_suppressed(int col);

  void
  retrieve_suppressed_columns_from_storage();

  int
  adjusted_correlated_columns();

  void
  set_adjusted_correlated_columns(int a)
  {
    _adjusted_correlated_columns = a;
  }

  //  The resizable_array<Correlated_Column> part holds the columns above threshold

  int
  above_threshold() const
  {
    return _number_elements;
  }

  int
  different_values_present() const;

  int
  initial_different_values_present() const
  {
    return _initial_different_values;
  }

  int columns_correlated_above(correlation_coeff_t) const;

  void
  set_ignore(int i)
  {
    _ignore = i;
  }

  int
  ignore() const
  {
    return _ignore;
  }

  void
  set_suppress(int i)
  {
    _suppress = i;
  }

  int
  suppress() const
  {
    return _suppress;
  }

  void
  set_title(const_IWSubstring& t)
  {
    _title = t;
  }

  const IWString&
  title() const
  {
    return _title;
  }

  void
  another_rho_value(double, int);

  void
  extra(double);

  void
  extra_missing_value();

  void
  look_for_correlation(SColumn& rhs);
};

#define STORED_MISSING_VALUE std::numeric_limits<double>::max()

SColumn::SColumn()
{
  _ignore = 0;
  _suppress = 0;

  _min_rho = 1.0;
  _min_rho_column = -1;

  _max_rho = -1.0;
  _max_rho_column = -1;

  _missing_value_count = 0;

  _collapsed_values = 0;

  _initial_different_values = 0;

  _adjusted_correlated_columns = -1;  // shows it is not initialised

  _sumx = 0.0;
  _sumx2 = -1.0;  // very important to be a negative number

  return;
}

void
SColumn::_compute_sum_and_sum2()
{
  _sumx = 0.0;
  _sumx2 = 0.0;

  const stored_value_t* x = _values_sampled.rawdata();

  const int n = _values_sampled.number_elements();

  for (int i = 0; i < n; ++i) {
    _sumx += x[i];
    _sumx2 += x[i] * x[i];
  }

  return;
}

void
SColumn::another_rho_value(double rho, int zcol)
{
  if (rho < _min_rho) {
    _min_rho = rho;
    _min_rho_column = zcol;
  }

  if (rho > _max_rho) {
    _max_rho = rho;
    _max_rho_column = zcol;
  }

  return;
}

/*
 */

void
SColumn::extra(double f)
{
  _values_sampled.add(static_cast<stored_value_t>(f));  // maybe no cast at all

  Accumulator<double>::extra(f);

  if (_different_values.size() < continuous_threshold || collect_all_values) {
    _different_values[static_cast<float>(f)]++;
  }

  return;
}

void
SColumn::extra_missing_value()
{
  _values_sampled.add(STORED_MISSING_VALUE);

  return;
}

int
SColumn::set_sample_size(int s)
{
  return _values_sampled.resize(s);
}

int
SColumn::different_values_present() const
{
  int n = Accumulator<double>::n();

  if (0 == n) {
    return 0;
  }

  return _different_values.size();
}

static SColumn* column = NULL;

/*
  Generic routine for writing column info.
  Note that the invoker considers columns to start with column 1, so
  we increment the value written by 1 to accommodate that.
*/

int
SColumn::write(std::ostream& os) const
{
  os << "column " << (_column_number + 1);  // convert to human numbering
  if (_title.length()) {
    os << " '" << _title << "' ";
  } else {
    os << ' ';
  }

  return os.good();
}

int
SColumn::print_column_statistics(std::ostream& os) const
{
  write(os);

  os << n() << " values";
  if (n()) {
    os << " between " << minval() << " and " << maxval();
  }
  if (_missing_value_count) {
    os << ", " << _missing_value_count << " missing values";
  }

  if (n() > 1) {
    os << " ave = " << average() << " var " << variance();
  }

  os << '\n';

  return os.good();
}

int
SColumn::write_min_max_rho(std::ostream& os) const
{
  write(os);
  os << ' ' << _values_sampled.number_elements() << " samples ";
  if (_collapsed_values) {
    os << _initial_different_values << " initial values, " << _collapsed_values
       << " collapsed values";
  } else {
    os << _different_values.size() << " different values";
  }

  if (_max_rho_column >= 0) {
    os << "\n    max ";
    if (_different_values.size() >= continuous_threshold) {
      os << "rho ";
    } else {
      os << "Q ";
    }
    os << _max_rho << ' ';

    column[_max_rho_column].write(os);
  }

  if (_min_rho_column >= 0) {
    os << "\n    min ";
    if (_different_values.size() >= continuous_threshold) {
      os << "rho ";
    } else {
      os << "Q ";
    }
    os << _min_rho << ' ';

    column[_min_rho_column].write(os);
  }

  os << '\n';

  return os.good();
}

int
SColumn::adjusted_correlated_columns()
{
  if (_adjusted_correlated_columns < 0) {
    //  cerr << "Column " << _column_number << " current value " <<
    //  _adjusted_correlated_columns << " initially " << _initial_different_values  << "
    //  values, " << _number_elements << " correlated\n";

    _adjusted_correlated_columns = _number_elements;  // the default value

    if (2 == _initial_different_values) {  // pure binary descriptor
      ;
    } else if (static_cast<unsigned int>(_initial_different_values) >=
               continuous_threshold) {  // continuous descriptor  - cast to keep the
                                        // compiler quiet
      ;
    } else if (0 == _number_elements) {  // nothing correlated with this column, can't be
                                         // adjusted
      ;
    } else {
      _compute_adjusted_correlated_columns();
    }
  }

  return _adjusted_correlated_columns;
}

/*
  The values from the -y option
*/

static int near_binary_addition_factor = 0;
static int near_binary_percent_factor = 0;

void
SColumn::_compute_adjusted_correlated_columns()
{
  int nva = _number_elements;  // the new value by adding
  int nvp = _number_elements;  // the new value by percent

  if (near_binary_addition_factor) {
    nva += near_binary_addition_factor;
  }

  if (near_binary_percent_factor) {
    nvp += _number_elements * near_binary_percent_factor / 100;
    if (nvp == _number_elements) {
      nvp++;
    }
  }

  if (nva >= nvp) {
    _adjusted_correlated_columns = nva;
  } else {
    _adjusted_correlated_columns = nvp;
  }

  if (verbose > 1 && _adjusted_correlated_columns > _number_elements) {
    cerr << "Column " << _column_number << " adjusted correlated columns "
         << _adjusted_correlated_columns << " initial " << _number_elements << '\n';
  }

  return;
}

static int
parse_number_or_percent(Command_Line& cl, char flag,
                        int& znumber,   // assumed initialised by someone else
                        int& zpercent)  // assumed initialised by someone else
{
  int i = 0;
  const_IWSubstring f;
  while (cl.value(flag, f, i++)) {
    if (f.ends_with('%')) {
      f.chop();
      if (!f.numeric_value(zpercent)) {
        cerr << "Invalid percentage specifier '" << f << "%'\n";
        return 0;
      }
    } else if (!f.numeric_value(znumber)) {
      cerr << "Invalid numeric specifier '" << f << "'\n";
      return 0;
    }
  }

  return 1;  // even if option 'flag' not present
}

/*
  Dealing with the ignore column is a little difficult, because
  the values come from the command line, before we know how many
  records there are in the file. So, we collect them in a
  resizable_array<int> and once the file size is known, we fill
  the array
*/

static resizable_array<int> ignore_columns_from_cl;

/*
  In descriptor files, we can also ignore certain columns based
  on their names
*/

static std::vector<std::unique_ptr<re2::RE2>> column_patterns_to_ignore;

static IWString missing_value('.');

static int is_descriptor_file = 0;

static int nskip = 0;

static int columns_in_input = 0;

static int columns_suppressed = 0;

static int test_mode = 0;

/*
  By default, we cease processing if we find that the file needs no changes.
  We can optionally write the file
  Jan 2007. Change this to being the default
*/

static int write_file_even_if_no_changes = 1;

static int
get_next_token(const const_IWSubstring& buffer, int& i, char sep,
               const_IWSubstring& token)
{
  if (' ' == sep) {
    return buffer.nextword(token, i);
  } else {
    return buffer.nextword_single_delimiter(token, i, sep);
  }
}

static int
establish_column_titles(const const_IWSubstring& buffer)
{
  if (' ' == input_delimitor) {
    columns_in_input = buffer.nwords();
  } else {
    columns_in_input = buffer.nwords_single_delimiter(input_delimitor);
  }

  if (columns_in_input < 2) {
    cerr << "Only " << columns_in_input << " columns in input\n";
    return 0;
  }

  int j = 0;
  const_IWSubstring token;
  for (int i = 0; i < columns_in_input; i++) {
    if (!get_next_token(buffer, j, input_delimitor, token)) {
      cerr << "Cannot extract column " << (i + 1) << " from header\n";
      return 0;
    }

    column[i].set_title(token);
  }

  return 1;
}

static int
echo_file(iwstring_data_source& input, IWString_and_File_Descriptor& output)
{
  input.seekg(0);

  const_IWSubstring buffer;
  while (input.next_record(buffer) && output.good()) {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  output.flush();

  return output.good();
}

/*
  Convert the old header to the new header by taking into account
  values in the suppress_columns array
*/

static void
suppress_columns(const char* s, int nchars, IWString& output_buffer)
{
  int in_word = 0;
  int zcolumn = -1;

  for (int i = 0; i < nchars; i++) {
    if (' ' == *s) {
      if (in_word && !output_buffer.ends_with(' ')) {  // don't want multiple spaces
        output_buffer += ' ';
      }
      in_word = 0;
    } else if (in_word) {
      if (!column[zcolumn].suppress()) {
        output_buffer += *s;
      }
    } else {
      zcolumn++;
      in_word = 1;
      if (!column[zcolumn].suppress()) {
        output_buffer += *s;
      }
    }

    s++;
  }

  output_buffer += '\n';

  return;
}

// A threshold for correlation.

static double threshold = 0.95;

/*
  We can specify the number of records to examine in determining correlations
*/

static int records_to_sample = 2000;

/*
  But because of missing values, we need to specify a minimum number of
  valid values
*/

static unsigned int min_records_needed = 0;

void
SColumn::look_for_correlation(SColumn& rhs)
{
  unsigned int n1 = different_values_present();
  unsigned int n2 = rhs.different_values_present();

  if (0 == _initial_different_values) {
    _initial_different_values = n1;
  }

  if (0 == rhs._initial_different_values) {
    rhs._initial_different_values = n2;
  }

#ifdef DEBUG_LOOK_FOR_CORRELATION
  cerr << "Looking for correlation between '" << _title << "' (" << n1 << " values) and '"
       << rhs._title << "' (" << n2 << " values)\n";
#endif

  if (2 == n1 && 2 == n2)  // compare binary descriptors with the Q test
  {
    _do_Q_test(rhs);
    return;
  }

  // Both continuous

  if (n1 >= continuous_threshold && n2 >= continuous_threshold) {
    _pearsons_test(rhs);
    return;
  }

  // If just one is continuous, we can't compare these

  if ((n2 >= continuous_threshold) || (n1 >= continuous_threshold)) {
    return;
  }

  // We must have one or both of us near binary

  if (0.0 == binary_collapse_threshold) {  // not doing collapsing
    return;
  }

  // Check to see whether collapsing has already been done,

  if (rhs._collapsed_values > 2 || _collapsed_values > 2) {
    return;
  }

  if (n1 > 2 && 0 == _collapsed_values) {
    _collapse_values();
    if (2 != _collapsed_values) {
      return;
    }
  }

  if (n2 > 2 && 0 == rhs._collapsed_values) {
    rhs._collapse_values();
    if (2 != rhs._collapsed_values) {
      return;
    }
  }

  _do_Q_test(rhs);

  return;
}

void
SColumn::_pearsons_test(SColumn& rhs)
{
  if (_sumx2 < 0.0) {
    _compute_sum_and_sum2();
  }

  if (rhs._sumx2 < 0.0) {
    rhs._compute_sum_and_sum2();
  }

  unsigned int valid_pairs = 0;

  int n = _values_sampled.number_elements();

  const stored_value_t* sv1 = _values_sampled.rawdata();
  const stored_value_t* sv2 = rhs._values_sampled.rawdata();

  double r = 0.0;
  double sum_x1 = 0.0;
  double sum_x2 = 0.0;
  double sum_squares_x1 = 0.0;
  double sum_squares_x2 = 0.0;

  if (0 == _missing_value_count && 0 == rhs._missing_value_count) {
    for (int i = 0; i < n; i++) {
      const stored_value_t x1 = sv1[i];
      const stored_value_t x2 = sv2[i];

      r += (x1 * x2);
    }

    valid_pairs = n;
    sum_x1 = _sumx;
    sum_x2 = rhs._sumx;
    sum_squares_x1 = _sumx2;
    sum_squares_x2 = rhs._sumx2;
  } else {
    for (int i = 0; i < n; i++) {
      stored_value_t x1 = sv1[i];
      if (STORED_MISSING_VALUE == x1) {
        continue;
      }
      stored_value_t x2 = sv2[i];
      if (STORED_MISSING_VALUE == x2) {
        continue;
      }

      valid_pairs++;
      r += (x1 * x2);

      sum_x1 += x1;
      sum_x2 += x2;
      sum_squares_x1 += x1 * x1;
      sum_squares_x2 += x2 * x2;

      //  cerr << "Valid pair " << x1 << " and " << x2 << '\n';
    }
  }

  if (valid_pairs < min_records_needed) {
    cerr << "Between ";
    write(cerr);
    cerr << "and ";
    rhs.write(cerr);
    cerr << "only " << valid_pairs << " valid pairs\n";
    return;
  }

  double x1bar = sum_x1 / static_cast<double>(valid_pairs);
  double x2bar = sum_x2 / static_cast<double>(valid_pairs);

  double nx1bx2b = static_cast<double>(valid_pairs) * x1bar * x2bar;

  // cerr << "r = " << r << " nx1bx2b = " << nx1bx2b << " x12bar " << x1bar << " x2bar "
  // << x2bar << '\n'; cerr << "Sum of squares " << sum_squares_x2 << " and " <<
  // sum_squares_x2 << '\n';

  double v1 = sum_squares_x1 - valid_pairs * x1bar * x1bar;
  double v2 = sum_squares_x2 - valid_pairs * x2bar * x2bar;

  // cerr << " v1 " << v2 << " v2 " << v2 << '\n';

  if (0.0 == v1 || 0.0 == v2) {
    cerr << "Yipes, zero qzyq term " << v1 << " ";
    write(cerr);
    cerr << " and " << v2 << " ";
    rhs.write(cerr);
    cerr << "\n";
    return;
  }

  double rho = (r - nx1bx2b) / sqrt(v1 * v2);

  if (rho >= threshold ||
      (take_absolute_value_of_correlation_coefficient && fabs(rho) >= threshold)) {
    column_is_correlated(rhs._column_number, rho);
    rhs.column_is_correlated(_column_number, rho);

    if (write_all_correlated_pairs) {
      _write_correlation(cerr, rhs, "Correlation", rho);
    }
  }

  another_rho_value(rho, rhs._column_number);
  rhs.another_rho_value(rho, _column_number);

  return;
}

// #define DEBUG_COLLAPSE_VALUES

/*
  This is hard coded for continuous_threshold being 5
*/

void
SColumn::_collapse_values()
{
  assert(5 == continuous_threshold);

#ifdef DEBUG_COLLAPSE_VALUES
  cerr << "Collapsing ";
  write(cerr);
  cerr << " ndiff = " << _different_values.size() << '\n';
#endif

  // First work out the binary density

  float f1, f2, f3, f4;
  f1 = 0.0f;  // Keep the compiler quiet.
  f2 = 0.0f;  // Keep the compiler quiet.
  f3 = 0.0f;  // Keep the compiler quiet.
  f4 = 0.0f;  // Keep the compiler quiet.
  int n1 = 0;
  int n2 = 0;
  int n3 = 0;
  int n4 = 0;

  int n = 0;

  // IW_Hash_Map<float, int, std::hash<float>, std::equal_to<float>>::const_iterator f;
  IW_Hash_Map<float, int, float_hash, eqflt>::const_iterator f;
  for (f = _different_values.begin(); f != _different_values.end(); f++) {
    int fn = (*f).second;
    n += fn;

#ifdef DEBUG_COLLAPSE_VALUES
    cerr << "  fn = " << fn << " value " << (*f).first << '\n';
#endif

    if (fn > n1) {
      n4 = n3;
      f4 = f3;
      n3 = n2;
      f3 = f2;
      n2 = n1;
      f2 = f1;

      n1 = fn;
      f1 = (*f).first;
    } else if (fn > n2) {
      n4 = n3;
      f4 = f3;
      n3 = n2;
      f3 = f2;

      n2 = fn;
      f2 = (*f).first;
    } else if (fn > n3) {
      n4 = n3;
      f4 = f3;

      n3 = fn;
      f3 = (*f).first;
    } else if (fn > n4) {
      n4 = fn;
      f4 = (*f).first;
    }
  }

#ifdef DEBUG_COLLAPSE_VALUES
  cerr << "f1 = " << f1 << " n1 = " << n1 << '\n';
  cerr << "f2 = " << f2 << " n2 = " << n2 << '\n';
  cerr << "f3 = " << f3 << " n3 = " << n3 << '\n';
  if (0 != n4) {
    cerr << "f4 = " << f4 << " n4 = " << n4 << '\n';
  }
#endif

  assert(n1 >= n2 && n2 >= n3 && n3 >= n4);

  float binary_density = static_cast<float>(n1 + n2) / static_cast<float>(n);

#ifdef DEBUG_COLLAPSE_VALUES
  cerr << "Binary density is " << binary_density << '\n';
#endif

  assert(binary_density > 0.0 && binary_density < 1.0);

  if (binary_density < binary_collapse_threshold)  // cannot collapse
  {
    _collapsed_values = _different_values.size();
    return;
  }

  // merge the two lower ones in with one of the large ones. n4 will always be set

  int ndiff = _different_values.size();

  _do_collapse(f1, f2, f3, n3);

  if (4 == ndiff) {
    _do_collapse(f1, f2, f4, n4);
  }

  _collapsed_values = _different_values.size();

  assert(2 == _collapsed_values);

  if (verbose > 1) {
    cerr << "Column " << _column_number << " collapsed from " << _initial_different_values
         << " to " << _collapsed_values << " values, binary density " << binary_density
         << '\n';
  }

  return;
}

/*
  F1 and F2 are the highest density items. We need to coalesce F with either F1 or F2
  Our first attempt will be to coalesce it with the one which is closer.
  If F1 and F2 are equidistant from F, then add F to whichever of F1 and F2
  is less populated
*/

void
SColumn::_do_collapse(float f1, float f2, float f, int n)
{
  float d1 = fabs(f - f1);
  float d2 = fabs(f - f2);

  if (d1 < d2) {  // closer to f1
    _do_collapse(f1, f, n);
  } else if (d1 > d2) {  // closer to f2
    _do_collapse(f2, f, n);
  } else if (_different_values[f1] <
             _different_values[f2]) {  // equal, add to lower density item
    _do_collapse(f1, f, n);
  } else {
    _do_collapse(f2, f, n);
  }

  return;
}

/*
  Coalesce F with FTO
*/

void
SColumn::_do_collapse(float fto, float f, int n)
{
  // cerr << "Collapsing " << f << " (" << n << " values) with " << fto << " (" <<
  // _different_values[fto] << " values)\n";

  _different_values.erase(f);

  _different_values[fto] += n;

  // Now update the stored values. Remember that there are N values to change

  int nv = _values_sampled.number_elements();

  const stored_value_t* sv = _values_sampled.rawdata();

  for (int i = 0; i < nv; i++) {
    stored_value_t x = sv[i];
    if (STORED_MISSING_VALUE == x) {
      continue;
    }

    if (f != static_cast<float>(x)) {
      continue;
    }

    _values_sampled[i] = static_cast<stored_value_t>(fto);

    n--;
    if (0 == n) {
      break;
    }
  }

  return;
}

// #define DEBUG_Q_TEST

void
SColumn::_do_Q_test(SColumn& rhs)
{
  assert(2 == _different_values.size());
  assert(2 == rhs._different_values.size());

  // IW_Hash_Map<float, int, std::hash<float>, std::equal_to<float>>::const_iterator f;
  // hash_map<float, int, float_hash, eqflt>::const_iterator f;
  std::unordered_map<float, int, float_hash, eqflt>::const_iterator f;

  // First get the two unique values for each column

  float a11;
  //float a12;

  f = _different_values.begin();
  a11 = (*f).first;
  f++;
  //a12 = (*f).first;

  float a21;
  //float a22;
  f = rhs._different_values.begin();
  a21 = (*f).first;
  f++;
  //a22 = (*f).first;

  int n = _values_sampled.number_elements();

  int valid_pairs = 0;

  // A count of the matrix elements

  int count_11 = 0;
  int count_12 = 0;
  int count_21 = 0;
  int count_22 = 0;

  const stored_value_t* sv1 = _values_sampled.rawdata();
  const stored_value_t* sv2 = rhs._values_sampled.rawdata();

  for (int i = 0; i < n; i++) {
    stored_value_t x1 = sv1[i];
    if (STORED_MISSING_VALUE == x1) {
      continue;
    }

    stored_value_t x2 = sv2[i];
    if (STORED_MISSING_VALUE == x2) {
      continue;
    }

    valid_pairs++;

    const float fx1 = static_cast<float>(x1);
    const float fx2 = static_cast<float>(x2);

    if (a11 == fx1) {
      if (a21 == fx2) {
        count_11++;
      } else {
        count_12++;
      }
    } else {  // a12 == fx1
      if (a21 == fx2) {
        count_21++;
      } else {
        count_22++;
      }
    }
  }

  // Dave Cummins email about problems with zero counts. Correct for this

  if (0 == count_11) {
    count_11 = 1;
  }
  if (0 == count_12) {
    count_12 = 1;
  }
  if (0 == count_21) {
    count_21 = 1;
  }
  if (0 == count_22) {
    count_22 = 1;
  }

  double q = static_cast<double>(count_11 * count_22 - count_12 * count_21) /
             static_cast<double>(count_11 * count_22 + count_12 * count_21);

  if (q * q >= threshold) {
    column_is_correlated(rhs._column_number, q * q);
    rhs.column_is_correlated(_column_number, q * q);

    if (write_all_correlated_pairs) {
      _write_correlation(cerr, rhs, "Q**2", q);
    }
  }

  another_rho_value(q, rhs._column_number);
  rhs.another_rho_value(q, _column_number);

#ifdef DEBUG_Q_TEST
  cerr << "Q test ";
  write(cerr);
  cerr << " and ";
  rhs.write(cerr);
  cerr << " Q = " << q << '\n';
#endif

  return;
}

/*
  Two columns have been found to be correlated. Write that info

  STATISTIC will be either "Correlation" or "Q" depending on the statistic
*/

void
SColumn::_write_correlation(std::ostream& os, const SColumn& rhs, const char* statistic,
                            double d) const
{
  os << statistic << " between ";
  write(os);
  os << "and ";
  rhs.write(os);
  os << "is " << d << '\n';

  return;
}

void
SColumn::column_is_correlated(int c, correlation_coeff_t r)
{
  assert(c >= 0 && c < columns_in_input);

  Correlated_Column* t = new Correlated_Column(c, r);

  resizable_array_p<Correlated_Column>::add(t);
}

/*
  Because it is correlated with some other column being retained, column COL
  is now being removed. If it is on our list of correlated columns, we must remove it
*/

void
SColumn::column_is_suppressed(int col)
{
  int removed = 0;

  for (int i = 0; i < _number_elements; i++) {
    if (col == _things[i]->column()) {
      remove_item(i);
      removed = 1;
      break;
    }
  }

  if (!removed) {  // not correlated with us
    return;
  }

  if (_adjusted_correlated_columns > 0) {
    _adjusted_correlated_columns--;
  }

  return;
}

/*
  During iterative work, we need to mark a column as having been suppressed. But, we need
  to store it somewhere so that when we start another iteration, we can restore the column
*/

void
SColumn::column_is_possibly_suppressed(int col)
{
  for (int i = 0; i < _number_elements; i++) {
    if (col == _things[i]->column()) {
      Correlated_Column* c = _things[i];
      _storage.add(c);
      remove_no_delete(i);

      return;
    }
  }

  return;
}

void
SColumn::retrieve_suppressed_columns_from_storage()
{
  int ns = _storage.number_elements();

  for (int i = 0; i < ns; i++) {
    add(_storage[i]);
  }

  _storage.resize_keep_storage(0);

  return;
}

int
SColumn::columns_correlated_above(correlation_coeff_t r) const
{
  correlation_coeff_t fr = static_cast<correlation_coeff_t>(r);

  int rc = 0;
  for (int i = 0; i < _number_elements; i++) {
    const Correlated_Column* c = _things[i];

    if (c->second >= fr) {
      rc++;
    }
  }

  return rc;
}

/*
  We've read a record. Continue to fill our buffer
*/

static int
update_sample(const_IWSubstring& buffer)
{
  int j = 0;
  const_IWSubstring token;
  for (int i = 0; i < columns_in_input; i++) {
    if (!get_next_token(buffer, j, input_delimitor, token)) {
      cerr << "Yipes, cannot extract column " << i << " from '" << buffer << "'\n";
      return 0;
    }

    SColumn& ci = column[i];

    if (ci.ignore()) {
      continue;
    }

    if (missing_value == token) {
      ci.extra_missing_value();
      continue;
    }

    double x;
    if (!token.numeric_value(x)) {
      cerr << "Yipes, invalid token '" << token << "', column " << (i + 1) << "\n";
      if (ignore_columns_with_non_numeric_data) {
        cerr << "Now ignoring column '" << ci.title() << "'\n";
        ci.set_ignore(1);

        continue;
      }

      return 0;
    }

    ci.extra(x);
  }

  return 1;
}

static int
allocate_values_sampled(int nr)
{
  assert(nr > 1);

  for (int i = 0; i < columns_in_input; i++) {
    if (!column[i].set_sample_size(nr)) {
      cerr << "YIpes, memory failure allocating space for " << nr << " samples\n";
      return 0;
    }
  }

  return 1;
}

static int
fill_sample_buffer_randomly(iwstring_data_source& input)
{
  off_t start_pos = input.tellg();
  off_t file_size = input.file_size();

  // Let's estimate the number of records in the file

  Accumulator_Int<int> record_length;

  if (!allocate_values_sampled(records_to_sample)) {
    return 0;
  }

  const_IWSubstring buffer;  // could be inside the loop

  // Probably should do something to allow setting the random number seed.
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<off_t> u(start_pos - 1,
                                         file_size - 17);  // 17 is just a number

  int rc = 0;
  while (1)  // until we get a valid record
  {
    const off_t p = u(rng);
    //  off_t p = intbtwij(static_cast<off_t>(start_pos) - 1, file_size - 17);    // 17 is
    //  just a number
    if (verbose > 1) {
      cerr << "Selecting random record " << rc << " near " << p << '\n';
    }

    (void)input.seekg(p);

    //  We have probably seek'd to somewhere in the middle of a record.

    (void)input.next_record(buffer);

    if (!input.next_record(buffer)) {  // musthave seeked into the last record
      continue;
    }

    if (!update_sample(buffer)) {
      cerr << "Error on line " << input.lines_read() << '\n';
      return 0;
    }

    rc++;

    record_length.extra(buffer.length());

    if (rc >= records_to_sample) {
      break;
    }
  }

  if (rc < 3) {  // hard to imagine this being meaningful
    return rc;
  }

  // How many records in this file

  int records_in_file = int((file_size - start_pos) / record_length.average()) + 1;

  if (records_in_file <= records_to_sample) {
    cerr << "***WARNING *** potential for serious problems\n";
    cerr << "Estimate only " << records_in_file << " records in file but asking for "
         << records_to_sample << " samples\n";
    return rc;
  }

  return rc;
}

static int
fill_sample_buffer_sequentially(iwstring_data_source& input)
{
  int nr = input.records_remaining();

  if (nr < 3) {
    cerr << "Too few records in file " << nr << '\n';
    return 0;
  }

  if (!allocate_values_sampled(nr)) {
    return 0;
  }

  const_IWSubstring buffer;

  int rc = 0;
  while (input.next_record(buffer)) {
    if (!update_sample(buffer)) {
      cerr << "Error on line " << input.lines_read() << '\n';
      return 0;
    }

    rc++;
  }

  return rc;
}

/*
  Let all the still active columns know that column COL is being suppressed

  Parent is the column number that is being processed. Make sure we don't change it
*/

static void
column_is_suppressed(int col, int parent)
{
  for (int i = 0; i < columns_in_input; i++) {
    SColumn& ci = column[i];

    if (ci.ignore() || ci.suppress()) {
      continue;
    }

    if (0 == ci.number_elements()) {
      continue;
    }

    if (i == col || parent == i) {
      continue;
    }

    ci.column_is_suppressed(col);
  }

  return;
}

static void
column_is_possibly_suppressed(int col, int parent)
{
  for (int i = 0; i < columns_in_input; i++) {
    SColumn& ci = column[i];

    if (ci.ignore() || ci.suppress()) {
      continue;
    }

    if (0 == ci.number_elements()) {
      continue;
    }

    if (i == col || parent == i) {
      continue;
    }

    ci.column_is_possibly_suppressed(col);
  }

  return;
}

static int
column_with_most_correlations()
{
  int rc = -1;
  int cmax = 0;

  for (int i = 0; i < columns_in_input; i++) {
    SColumn& ci = column[i];

    if (ci.ignore() || ci.suppress()) {
      continue;
    }

    int c = ci.adjusted_correlated_columns();

    if (c > cmax) {
      cmax = c;
      rc = i;
    }
  }

  return rc;
}

static int
column_with_most_correlations(double r, const int* column_already_done)
{
  int rc = -1;
  int cmax = 0;

  for (int i = 0; i < columns_in_input; i++) {
    if (column_already_done[i]) {
      continue;
    }

    SColumn& ci = column[i];

    if (ci.ignore() || ci.suppress()) {
      continue;
    }

    int c = ci.columns_correlated_above(r);

    if (c > cmax) {
      cmax = c;
      rc = i;
    }
  }

  return rc;
}

/*
  If the user asks for a random selection of records, but the file has fewer than
  that number of records, we automatically switch them to doing each record
*/

static void
_check_for_too_few_records_in_file(iwstring_data_source& input, int records_needed)
{
  const_IWSubstring buffer;

  int records_read = nskip;

  while (input.next_record(buffer)) {
    records_read++;
    if (records_read >= records_needed) {
      return;
    }
  }

  // If we come to here, there aren't enough records in the file

  cerr << "You asked for a random sample of " << records_needed << " records\n";
  cerr << "Input file contains " << records_read << " records (nskip = " << nskip
       << ")\n";

  if (records_read - nskip <= 2) {
    cerr << "Cannot do correlations with 2 columns\n";
  } else {
    cerr << "All records will be sampled\n";
  }

  sample_records_randomly = 0;

  return;
}

static void
check_for_too_few_records_in_file(iwstring_data_source& input, int records_needed)
{
  off_t o = input.tellg();

  _check_for_too_few_records_in_file(input, records_needed);

  if (!input.seekg(o)) {
    cerr << "Yipes, cannot seek back to beginning of file, catastrophe looming...\n";
  }

  return;
}

/*
  When doing the iterative pruning, we need a means of temporarily
  marking a column as having been suppressed
*/

#define TEMP_SUPPRESS -1
#define PERMANENTLY_SUPPRESS 1

/*
  A column is being suppressed for being correlated to another
*/

static void
process_correlated_column(int zcol, int suppress_value)
{
  if (!column[zcol].suppress())  // may have already been suppressed by another column
  {
    column[zcol].set_suppress(suppress_value);
    if (PERMANENTLY_SUPPRESS == suppress_value) {
      columns_suppressed++;
    }
  }

  if (verbose > 1) {
    cerr << "   ";
    column[zcol].write(cerr);
    cerr << " is correlated\n";
  }

  return;
}

/*
  the user has specified how many columns they want to keep
*/

static int
do_columns_to_keep(double r, int* column_already_done, int& columns_remaining)
{
  for (int i = 0; i < columns_in_input; i++) {
    SColumn& ci = column[i];

    cerr << "TColumn " << (i + 1) << "'" << ci.title() << "' has " << ci.number_elements()
         << " correlated columns\n";

    if (ci.ignore()) {
      column_already_done[i] = 1;
    } else if (PERMANENTLY_SUPPRESS == ci.suppress()) {  // maybe was constant
      column_already_done[i] = 1;
    } else if (TEMP_SUPPRESS == ci.suppress()) {
      ci.set_suppress(0);

      ci.retrieve_suppressed_columns_from_storage();

      column_already_done[i] = 0;
    }
  }

  cerr << "DEBUG:Starting elimination at " << r << '\n';

  int keep;
  while ((keep = column_with_most_correlations(r, column_already_done)) >= 0) {
    SColumn& ck = column[keep];

    if (verbose) {
      ck.write(cerr);
      cerr << " has " << ck.number_elements() << " correlated columns\n";
    }

    for (int i = 0; i < ck.number_elements(); i++) {
      const Correlated_Column* c = ck[i];

      if (c->r2() <= r) {  // only eliminate columns > current threshold
        continue;
      }

      process_correlated_column(c->column(), TEMP_SUPPRESS);

      column_is_possibly_suppressed(c->column(), keep);

      columns_remaining--;
    }

    if (columns_remaining <= columns_to_keep) {
      return 1;
    }

    //  Make sure this column doesn't get chosen next time

    column_already_done[keep] = 1;
  }

  if (verbose > 1) {
    cerr << "At R = " << r << ", " << columns_remaining << " columns remain\n";
  }

  return 0;  // never got below columns_to_keep
}

static int
do_columns_to_keep(int* keep_column, int* column_already_done)
{
  int columns_remaining = 0;  // an initial count

  for (int i = 0; i < columns_in_input; i++) {
    const SColumn& ck = column[i];
    if (ck.suppress() || ck.ignore()) {
      keep_column[i] = 0;
    } else {
      columns_remaining++;
      keep_column[i] = 1;
    }
  }

  if (columns_remaining <= columns_to_keep) {
    cerr << "At start of reduction, only " << columns_remaining
         << " columns active, cannot reduce to " << columns_to_keep << '\n';
    return 1;
  }

  if (verbose > 1) {
    cerr << "Before paring, " << columns_remaining << " columns remaining\n";
  }

  int fewest_columns_remaining = columns_remaining;
  float threshold_with_fewest_columns_remaining = 0.99;

  for (double r = 0.99; r >= threshold; r -= 0.02) {
    int tmp = columns_remaining;
    if (do_columns_to_keep(r, column_already_done, tmp)) {
      if (verbose) {
        cerr << "Trimmed to " << columns_to_keep << " at r = " << r << '\n';
      }

      return 1;
    }

    if (tmp < fewest_columns_remaining) {
      fewest_columns_remaining = tmp;
      threshold_with_fewest_columns_remaining = r;
      for (int i = 0; i < columns_in_input; i++) {
        const SColumn& ck = column[i];
        if (ck.ignore()) {
          continue;
        }

        if (TEMP_SUPPRESS == ck.suppress()) {
          keep_column[i] = 0;
        } else {
          keep_column[i] = 1;
        }
      }
    }
    cerr << "DEBUG:Trimmed to " << tmp << " at r = " << r << '\n';
  }

  cerr << "No threshold allowed reduction to " << columns_to_keep << " columns\n";
  cerr << "Fewest columns remaining " << fewest_columns_remaining
       << " ar R = " << threshold_with_fewest_columns_remaining << '\n';

  return 1;
}

static int
do_columns_to_keep()
{
  int* keep_column = new int[columns_in_input];
  std::unique_ptr<int[]> free_keep_column(keep_column);
  // int * keep_column = new int[columns_in_input]; std::unique_ptr<int>
  // free_keep_column(keep_column);

  int* column_already_done = new int[columns_in_input];
  std::unique_ptr<int[]> free_column_already_done(column_already_done);
  // int * column_already_done = new int[columns_in_input]; std::unique_ptr<int>
  // free_column_already_done(column_already_done);

  int rc = do_columns_to_keep(keep_column, column_already_done);

  for (int i = 0; i < columns_in_input; i++) {
    SColumn& sc = column[i];

    if (sc.ignore()) {
      continue;
    }

    if (0 == keep_column[i]) {
      sc.suppress();
    }
  }

  return rc;
}

/*
  We discard all columns correlated above the threshold
*/

static int
do_discard_correlated_columns()
{
  int keep;
  while ((keep = column_with_most_correlations()) >= 0) {
    SColumn& ck = column[keep];

    if (verbose) {
      ck.write(cerr);
      cerr << " has " << ck.number_elements() << " correlated columns\n";
    }

    for (int i = 0; i < ck.number_elements(); i++) {
      const Correlated_Column* c = ck[i];

      process_correlated_column(c->column(), PERMANENTLY_SUPPRESS);

      column_is_suppressed(c->column(), keep);
    }

    //  Make sure this column doesn't get chosen next time

    ck.resize(0);
    ck.set_adjusted_correlated_columns(0);
  }

  return 1;
}

static int
correlated(iwstring_data_source& input, IWString_and_File_Descriptor& output)
{
  off_t psave = input.tellg();  // need to be able to start @ the same point

  if (sample_records_randomly) {
    check_for_too_few_records_in_file(input, records_to_sample);
  }

  int number_records_sampled;

  if (sample_records_randomly) {
    number_records_sampled = fill_sample_buffer_randomly(input);
  } else {
    number_records_sampled = fill_sample_buffer_sequentially(input);
  }

  // cerr << "Sampled " << number_records_sampled << " records\n";

  if (0 == number_records_sampled) {  // something horrible happened
    return 0;
  }

  if (number_records_sampled <= 2) {
    cerr << "Cannot do correlations on " << number_records_sampled << " rows\n";
    return 1;
  }

  if (sample_records_randomly && number_records_sampled < records_to_sample) {
    cerr << "Must sample " << records_to_sample << " records, but got only "
         << number_records_sampled << '\n';
    return 0;
  }

  if (0 == min_records_needed) {
    min_records_needed = number_records_sampled;
  }

  input.seekg(psave);  // back to where we started

  // Deal with any constant columns

  for (int i = 0; i < columns_in_input; i++) {
    SColumn& ci = column[i];

    if (ci.ignore()) {
      continue;
    }

    if (!print_column_statistics) {
      ;
    } else if (stream_for_column_statistics.rdbuf()->is_open()) {
      ci.print_column_statistics(stream_for_column_statistics);
    } else {
      ci.print_column_statistics(cerr);
    }

    if (0 == ci.n()) {
      ci.set_ignore(1);
      continue;
    }

    if (ci.minval() == ci.maxval()) {
      columns_suppressed++;

      if (verbose) {
        ci.write(cerr);
        cerr << "is constant for " << ci.n() << " samples\n";
      }

      if (suppress_constant_columns) {
        if (0 == ci.missing_value_count()) {
          ci.set_suppress(PERMANENTLY_SUPPRESS);
        } else if (ci.n() >= min_records_needed) {
          ci.set_suppress(PERMANENTLY_SUPPRESS);
        } else if (suppress_constant_columns > 1) {
          ci.set_suppress(PERMANENTLY_SUPPRESS);
        } else {
          ci.set_ignore(PERMANENTLY_SUPPRESS);
        }
      } else {
        ci.set_ignore(1);  // cannot do correlations against these
      }
    }
  }

  if (verbose) {
    int ignored = 0;
    int suppressed = 0;

    for (int i = 0; i < columns_in_input; i++) {
      const SColumn& c = column[i];

      if (c.ignore()) {
        ignored++;
      } else if (c.suppress()) {
        suppressed++;
      }
    }

    cerr << "After initial analysis, of " << columns_in_input << " columns, " << ignored
         << " are ignored and " << suppressed << " are suppressed\n";
  }

  // Check for correlated columns

  for (int i = 0; i < columns_in_input; i++) {
    if (column[i].ignore() || column[i].suppress()) {
      continue;
    }

    SColumn& ci = column[i];

    if (report_progress()) {
      cerr << "Processing ";
      ci.write(cerr);
      cerr << '\n';
    }

    if (ci.n() < min_records_needed) {
      cerr << "Yipes, column ";
      ci.write(cerr);
      cerr << "only " << ci.n() << " valid values found in sampling. Ignored\n";

      continue;
    }

    for (int j = i + 1; j < columns_in_input; j++) {
      if (column[j].ignore() || column[j].suppress()) {
        continue;
      }

      if (column[j].n() < min_records_needed) {
        continue;
      }

      ci.look_for_correlation(column[j]);
    }
  }

  if (verbose) {
    for (int i = 0; i < columns_in_input; i++) {
      const SColumn& ci = column[i];

      if (ci.ignore() || ci.suppress()) {
        continue;
      }

      ci.write_min_max_rho(cerr);

      if (ci.above_threshold()) {
        cerr << "  " << ci.above_threshold()
             << " other columns above correlation threshold\n";
      }
    }
  }

  // Identify the column with the most correlations and remove all that are
  // correlated with it

  if (verbose > 1) {
    cerr << "Beginning column elimination\n";
  }

  if (columns_to_keep > 0) {
    do_columns_to_keep();
  } else {
    do_discard_correlated_columns();
  }

  if (test_mode) {
    return 1;
  }

  if (0 == columns_suppressed) {
    if (verbose) {
      cerr << "No columns suppressed\n";
    }

    if (write_file_even_if_no_changes) {
      return echo_file(input, output);
    } else {
      return 1;
    }
  }

  if (verbose) {
    cerr << "Suppressed " << columns_suppressed << " columns, output will contain "
         << (columns_in_input - columns_suppressed) << " columns\n";
  }

  // write the file

  input.seekg(0);

  IWString buffer;

  while (input.next_record(buffer)) {
    suppress_columns(buffer.rawchars(), buffer.length(), output);

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
correlated(const char* fname, IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  // Determine the number of tokens per record

  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    cerr << "Cannot read header from '" << fname << "'\n";
    return 0;
  }

  if (' ' == input_delimitor) {
    columns_in_input = buffer.nwords();
  } else {
    columns_in_input = buffer.nwords_single_delimiter(input_delimitor);
  }

  if (verbose) {
    cerr << "Input contains " << columns_in_input << " columns\n";
  }

  if (0 == columns_in_input) {
    cerr << "Very bad news, header record contains zero columns\n";
    return 0;
  }

  if (1 == columns_in_input) {
    cerr << "Only one column in input, correlation studies not possible\n";
    return 0;
  }

  column = new SColumn[columns_in_input];
  assert(NULL != column);

  if (is_descriptor_file) {
    if (!establish_column_titles(buffer)) {
      return 0;
    }
  } else if (nskip) {
    for (int i = 0; i < nskip; i++) {
      if (!input.next_record(buffer)) {
        cerr << "Premature EOF, record " << input.lines_read() << '\n';
        return 0;
      }
    }
  } else {
    if (!input.seekg(0)) {
      cerr << "Yipes, cannot seek back to beginning of file\n";
      return 0;
    }
  }

  for (int i = 0; i < columns_in_input; i++) {
    column[i].set_column_number(i);
    column[i].resize(100);  // just a guess
  }

  for (int i = 0; i < ignore_columns_from_cl.number_elements(); i++) {
    int j = ignore_columns_from_cl[i];
    column[j].set_ignore(1);
  }

  for (std::unique_ptr<re2::RE2>& rx : column_patterns_to_ignore) {
    for (int j = 0; j < columns_in_input; j++) {
      if (column[j].ignore()) {
        continue;
      }

      if (iwre2::RE2PartialMatch(column[j].title(), *rx)) {
        column[j].set_ignore(1);
        if (verbose) {
          cerr << "Column " << (j + 1) << " '" << column[j].title()
               << "' ignored by matching '" << rx->pattern() << "'\n";
        }
      }
    }
  }

  input.reset_record_counter();

  return correlated(input, output);
}

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Samples a data file looking for correlated columns\n";
  cerr << "Usage: " << prog_name << " <options> <input_file>\n";
  cerr << "  -s <skip>      skip the first <skip> records of the file\n";
  cerr << "  -i <column>    specify columns to ignore (starting @ 1)\n";
  cerr << "  -I <pattern>   specify columns to ignore based on headers\n";
  cerr << "  -j             is descriptor file (implies -s 1 -i 1)\n";
  cerr << "  -m <string>    specify missing value string (default '" << missing_value << "')\n";
  cerr << "  -N <records>   specify the number of records to sample (default " << default_records_to_sample << ")\n";
  cerr << "  -a             sample the whole file\n";
  cerr << "  -n <records>   min number of pairs needed for correlation determinations\n";
  cerr << "  -t <float>     specify correlation threshold between 0.0 and 1.0 (default " << threshold << ")\n";
  cerr << "  -T             test mode (does no output)\n";
  cerr << "  -c             suppress constant columns\n";
  cerr << "  -c             suppress constant columns (even if missing values present)\n";
  cerr << "  -w             write file even if no changes needed\n";
  cerr << "  -p             print column statistics\n";
  cerr << "  -P <fname>     print column statistics to <fname>\n";
  cerr << "  -e             write all correlated pairs (awful on large files)\n";
  cerr << "  -C <ct>        threshold for continuous descriptors (default " << continuous_threshold << ")\n";
  cerr << "  -b <fraction>  binary density for collapsing near binary descriptors\n";
  cerr << "  -y <number>    add <number> to the correlated columns count for collapsed descriptors\n";
  cerr << "  -y <percent>   add <percent> to the correlated columns count for collapsed descriptors\n";
  cerr << "  -k <number>    trim down to <number> or fewer columns\n";
  cerr << "  -G             ignore columns containing non-numeric data (normally a fatal error)\n";
  cerr << "  -h <seed>      seed for random number generator when sampling records\n";
//cerr << "  -o             take absolute values of correlation coefficients\n";
  cerr << "  -d <char>      specify input token delimitor\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
correlated(int argc, char** argv)
{
  Command_Line cl(argc, argv, "Tt:jvm:s:i:I:n:N:r:wacC:pP:eb:y:k:Goh:d:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('t')) {
    if (!cl.value('t', threshold) || threshold < 0.0 || threshold > 1.0) {
      cerr << "The -t option must be followed by a value possible R^2\n";
      usage(3);
    }
    if (verbose) {
      cerr << "Correlation threshold " << threshold << '\n';
    }
  } else {
    if (cl.option_present('k')) {  // lower default threshold for columns_to_keep option
      threshold = 0.60;
    } else {
      threshold = 0.95;
    }
    if (verbose) {
      cerr << "Using default threshold " << threshold << '\n';
    }
  }

  if (cl.option_present('T')) {
    test_mode = 1;
    cerr << "Working in test mode only\n";
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', nskip) || nskip < 1) {
      cerr << "the skip option requires a whole positive number\n";
      usage(43);
    }

    if (verbose) {
      cerr << "Will skip first " << nskip << " records of the file\n";
    }
  }

  if (cl.option_present('j')) {
    is_descriptor_file = 1;
    nskip = 1;
    ignore_columns_from_cl.add(0);
    if (verbose) {
      cerr << "Will process as a descriptor file\n";
    }
  }

  if (cl.option_present('m')) {
    cl.value('m', missing_value);
    if (verbose) {
      cerr << "Missing value string '" << missing_value << "'\n";
    }
  }

  if (cl.option_present('i')) {
    ignore_columns_from_cl.resize(cl.option_count('i'));

    int i = 0;
    int j;
    while (cl.value('i', j, i++))  // does not handle errors well
    {
      if (j < 1) {
        cerr << "Column numbers start with column 1, value " << j << " is invalid\n";
        return 73;
      }

      ignore_columns_from_cl.add(j - 1);
    }
  }

  if (cl.option_present('I')) {
    is_descriptor_file = 1;
    nskip = 1;
    ignore_columns_from_cl.add_if_not_already_present(0);
    if (verbose) {
      cerr << "Will process as a descriptor file\n";
    }

    int i = 0;
    const_IWSubstring p;
    while (cl.value('I', p, i++)) {
      std::unique_ptr<re2::RE2> rx;
      if (!iwre2::RE2Reset(rx, p)) {
        cerr << "Cannot parse ignore column regular expression '" << p << "'\n";
        return i;
      }

      column_patterns_to_ignore.push_back(std::move(rx));
    }
  }

  if (cl.option_present('N') && cl.option_present('a')) {
    cerr << "The -N and -a options are mutually exclusive\n";
    usage(7);
  }

  if (cl.option_present('a')) {
    sample_records_randomly = 0;
    if (verbose) {
      cerr << "Will sample the whole file\n";
    }
  } else if (cl.option_present('N')) {
    if (!cl.value('N', records_to_sample) || records_to_sample < 3) {
      cerr << "The -N switch requires a whole positive number > 2\n";
      usage(14);
    }

    if (verbose) {
      cerr << "Will sample " << records_to_sample
           << " records in determining correlations\n";
    }
  } else {
    records_to_sample = default_records_to_sample;
    if (verbose) {
      cerr << "Will sample " << records_to_sample << " records\n";
    }
  }

#ifdef IMPLEMENT_AND_UPDATE_THIS_IF_IT_IS_EVER_NEEDED
  if (sample_records_randomly) {
    if (cl.option_present('h')) {
      unsigned long seed;
      if (!cl.value('h', seed)) {
        cerr << "Cannot discern random number seed (-h option)\n";
        usage(2);
      }

      iw_set_rnum_seed(seed);
    } else {
      (void)iw_random_seed();
    }
  }
#endif

  if (cl.option_present('n')) {
    if (!cl.value('n', min_records_needed) || min_records_needed < 3 ||
        min_records_needed > static_cast<unsigned int>(records_to_sample)) {
      cerr << "The -n option must be followed to a number between 3 and "
           << records_to_sample << '\n';
      usage(8);
    }

    if (verbose) {
      cerr << "Must find a minimum of " << min_records_needed
           << " pairs between columns for correlation determination\n";
    }
  }

  if (cl.option_present('w')) {
    write_file_even_if_no_changes = 1;
    if (verbose) {
      cerr << "File will be written even if no changes needed\n";
    }
  }

  if (cl.option_present('c')) {
    suppress_constant_columns = cl.option_count('c');
    if (0 == verbose) {
      ;
    } else if (1 == suppress_constant_columns) {
      cerr << "Will suppress constant columns\n";
    } else {
      cerr << "Will suppress constant columns even with some missing values\n";
    }
  }

  if (cl.option_present('p')) {
    print_column_statistics = 1;
    if (verbose) {
      cerr << "Will print column statistics\n";
    }
  }

  if (cl.option_present('d')) {
    IWString d = cl.string_value('d');
    if (!char_name_to_char(d)) {
      cerr << "Cannot recognise input delimintor name '" << d << "'\n";
      return 1;
    }

    input_delimitor = d[0];
  }

  if (cl.option_present('P')) {
    const char* p = cl.option_value('P');

    stream_for_column_statistics.open(p, std::ios::out);

    if (!stream_for_column_statistics.good()) {
      cerr << "Cannot open column statistics file '" << p << "'\n";
      return 5;
    }

    if (verbose) {
      cerr << "Column statistics written to '" << p << "'\n";
    }

    print_column_statistics = 1;
  }

  if (cl.option_present('G')) {
    ignore_columns_with_non_numeric_data = 1;

    if (verbose) {
      cerr << "Columns with non-numeric data will be discarded rather than causing a "
              "fatal error\n";
    }
  }

  if (cl.option_present('o')) {
    take_absolute_value_of_correlation_coefficient = 0;

    if (verbose) {
      cerr << "Will take absolute value of correlation coefficients\n";
    }
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', continuous_threshold) || continuous_threshold <= 2) {
      cerr << "The continuous threshold option (-C) must be a number > 2\n";
      usage(19);
    }

    if (verbose) {
      cerr << "descriptors having " << continuous_threshold
           << " or more different values will be considered continuous\n";
    }
  }

  if (cl.option_present('e')) {
    write_all_correlated_pairs = 1;
    if (verbose) {
      cerr << "All correlated pairs will be reported\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise progress reporting (-r)\n";
      return 1;
    }
  }

  if (cl.option_present('b')) {
    if (!cl.value('b', binary_collapse_threshold) || binary_collapse_threshold <= 0.0 ||
        binary_collapse_threshold >= 1.0) {
      cerr << "The binary denisity collapse threshold (-b) option must be followed by a "
              "valid ratio\n";
      usage(14);
    }

    if (verbose) {
      cerr << "Binary collapse denisty set to " << binary_collapse_threshold << '\n';
    }
  }

  if (cl.option_present('y')) {
    if (!parse_number_or_percent(cl, 'y', near_binary_addition_factor,
                                 near_binary_percent_factor)) {
      cerr << "Cannot parse the -y option(s)\n";
      usage(21);
    }

    if (verbose && near_binary_addition_factor) {
      cerr << "Will add " << near_binary_addition_factor
           << " to each near-binary descriptor\n";
    }
    if (verbose && near_binary_percent_factor) {
      cerr << "Will increase the correlated column count for near binary descriptors by "
           << near_binary_percent_factor << " percent\n";
    }
  }

  if (cl.option_present('k')) {
    if (!cl.value('k', columns_to_keep) || columns_to_keep < 1) {
      cerr << "The columns to keep option (-k) must be a whole positive number\n";
      usage(11);
    }

    if (verbose) {
      cerr << "Will trim to " << columns_to_keep << " or fewer columns\n";
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(3);
  }

  if (cl.number_elements() > 1) {
    cerr << "Extra arguments '" << cl[1] << "...' ignored\n";
  }

  IWString_and_File_Descriptor output(1);

  int rc = correlated(cl[0], output);

  return !rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = correlated(argc, argv);

  return rc;
}
