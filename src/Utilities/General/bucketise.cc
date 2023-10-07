#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <random>
#include <unordered_set>

#include "Foundational/iwstring/iw_stl_hash_set.h"

/*
  Convert continuous values to bucket form.
  Jun 2002. Implement Glaxo scheme for dealing with outliers.
  The top and bottom N% of the data go in their own buckets,
  the rest go in equally sized buckets
*/

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iwstring.h"

const char* prog_name = nullptr;

using std::cerr;

static int verbose = 0;

static int print_statistics = 0;

static int sample_records_randomly = 0;

static int translate_tabs = 0;

static int only_write_non_ignored_columns = 0;

typedef float dscvalue;

static dscvalue global_min = static_cast<dscvalue>(0.0);
static dscvalue global_max = static_cast<dscvalue>(0.0);
static dscvalue global_dx = static_cast<dscvalue>(0.0);

static int global_range_active = 0;

static int items_below_global_range = 0;
static int items_above_global_range = 0;

/*
  If we are just interested in profiling, we may not want
  individual items bucketised
*/

static int suppress_output = 0;

/*
  Dealing with the ignore column is a little difficult, because
  the values come from the command line, before we know how many
  records there are in the file. So, we collect them in a
  resizable_array<int> and once the file size is known, we fill
  the array
*/

static resizable_array<int> ignore_columns_from_cl;

/*
  Alternatively, we can process only columns which match a regexp (descriptor files only)
*/

static std::unique_ptr<re2::RE2> process_column_regexp;

static resizable_array<int> columns_to_process_from_cl;

static IWString missing_value('.');

static extending_resizable_array<int> missing_value_counter;

static int is_descriptor_file = 0;

static int nskip = 0;

static int columns_in_input = 0;

static int nbuckets = 0;

/*
  It is much more efficient to write text than numbers
*/

static IWString* text_bucket = nullptr;

static int* bucket_count = nullptr;

static int global_bucketise = 0;

#include "Foundational/data_source/iwstring_data_source.h"

/*
  For each column, we need to keep track of the data in the column.
  We work in two modes
    just an accumulator
    The Glaxo method
*/

static double fraction_items_in_extremity_buckets = 0.0;

/*
  When the first column's compute_dx() method is called, it can
  work out number_items_in_extremity_buckets
*/

static int number_items_in_extremity_buckets = -1;

class Column_Data_Base : public Accumulator<dscvalue>
{
 protected:
  IWString _column_id;

  //  Once the extremeties are known, we get minval and dx

  dscvalue _mv;
  dscvalue _dx;

  int _ignore;

  int _is_fingerprint;

 public:
  Column_Data_Base();
  virtual ~Column_Data_Base();

  virtual int debug_print(std::ostream&) const;

  int write_column_title(int, std::ostream&) const;

  const IWString& id() const {
    return _column_id;
  }

  void set_id(const const_IWSubstring& t) {
    _column_id = t;
  }

  void set_ignore(int i) {
    _ignore = i;
  }

  int ignore() const {
    return _ignore;
  }

  int is_fingerprint() const {
    return _is_fingerprint;
  }

  void set_is_fingerprint(int s) {
    _is_fingerprint = s;
  }

  dscvalue dx() const {
    return _dx;
  }

  virtual int extra(const const_IWSubstring&) = 0;

  dscvalue minval() const {
    return Accumulator<dscvalue>::minval();
  }

  dscvalue maxval() const {
    return Accumulator<dscvalue>::maxval();
  }

  virtual int compute_dx() = 0;
  virtual int bucket(dscvalue) const = 0;

  void set_output_specifications(dscvalue m, dscvalue dx);
};

Column_Data_Base::Column_Data_Base()
{
  _dx = -1.0;
  _ignore = 0;
  _is_fingerprint = 1;

  return;
}

Column_Data_Base::~Column_Data_Base()
{
  _dx = -1;
  _ignore = -2;
  _is_fingerprint = -5;

  return;
}

int
Column_Data_Base::debug_print(std::ostream& os) const
{
  os << "Column_Data_Base::debug_print: for '" << _column_id << "'\n";

  os << Accumulator<dscvalue>::n() << " values";
  if (Accumulator<dscvalue>::n() > 0) {
    os << " between " << minval() << " and " << maxval();
    if (Accumulator<dscvalue>::n() > 1) {
      os << " ave " << Accumulator<dscvalue>::average();
    }
  }

  os << '\n';

  return os.good();
}

void
Column_Data_Base::set_output_specifications(dscvalue mv, dscvalue dx)
{
  assert(dx >= 0.0);

  _mv = mv;
  _dx = dx;

  if (0.0 == _dx) {
    _ignore = 1;
  }

  return;
}

class Column_Data_Accumulator : public Column_Data_Base
{
 private:
 public:
  Column_Data_Accumulator();

  int allocate_space_for_this_many_rows(int) {
    return 1;
  }  // not needed

  int compute_dx();

  int debug_print(std::ostream&) const;

  int extra(const const_IWSubstring&);

  int bucket(dscvalue) const;
};

Column_Data_Accumulator::Column_Data_Accumulator()
{
  return;
}

int
Column_Data_Accumulator::debug_print(std::ostream& os) const
{
  Column_Data_Base::debug_print(os);

  return os.good();
}

int
Column_Data_Accumulator::compute_dx()
{
  _mv = Accumulator<dscvalue>::minval();

  _dx = (Accumulator<dscvalue>::maxval() - _mv) / static_cast<dscvalue>(nbuckets);

  if (0.0 == _dx) {
    _ignore = 1;
  }

  return 1;
}

int
Column_Data_Accumulator::extra(const const_IWSubstring& token)
{
  dscvalue d;

  if (_is_fingerprint && '0' == token) {
    d = static_cast<dscvalue>(0.0);
  } else if (_is_fingerprint && '1' == token) {
    d = static_cast<dscvalue>(1.0);
  } else if (!token.numeric_value(d)) {
    cerr << "Column_Data_Accumulator::extra: invalid number '" << token << "'\n";
    return 0;
  } else {
    _is_fingerprint = 0;
  }

  Accumulator<dscvalue>::extra(d);

  return 1;
}

static int
compute_bucket(dscvalue minval, dscvalue dx, dscvalue d)
{
  int rc = static_cast<int>(static_cast<dscvalue>(0.5) + (d - minval) / dx);

  return rc;
}

static int
global_bucket_value(dscvalue d)
{
  if (d <= global_min) {
    if (d < global_min) {
      items_below_global_range++;
    }
    return 0;
  }

  if (d >= global_max) {
    if (d > global_max) {
      items_above_global_range++;
    }
    return nbuckets - 1;
  }

  return static_cast<int>(static_cast<dscvalue>(0.5) + (d - global_min) / global_dx);
}

int
Column_Data_Accumulator::bucket(dscvalue d) const
{
  if (global_range_active) {
    return global_bucket_value(d);
  }

  return compute_bucket(_mv, _dx, d);
}

class Column_Data_Glaxo : public Column_Data_Base
{
 private:
  resizable_array<dscvalue> _values;

  //  We need the locations of the two buckets

  dscvalue _lower_bucket_boundary;
  dscvalue _upper_bucket_boundary;
  dscvalue _dx;

 public:
  Column_Data_Glaxo();
  ~Column_Data_Glaxo();

  int allocate_space_for_this_many_rows(int);
  int extra(const const_IWSubstring&);
  int compute_dx();
  int debug_print(std::ostream&) const;
  int bucket(dscvalue) const;
};

Column_Data_Glaxo::Column_Data_Glaxo()
{
  return;
}

Column_Data_Glaxo::~Column_Data_Glaxo()
{
  return;
}

int
Column_Data_Glaxo::debug_print(std::ostream& os) const
{
  Column_Data_Base::debug_print(os);

  return os.good();
}

int
Column_Data_Glaxo::extra(const const_IWSubstring& token)
{
  dscvalue d;

  if (_is_fingerprint && '0' == token) {
    d = static_cast<dscvalue>(0.0);
  } else if (_is_fingerprint && '1' == token) {
    d = static_cast<dscvalue>(1.0);
  } else if (!token.numeric_value(d)) {
    cerr << "Column_Data_Accumulator::extra: invalid number '" << token << "'\n";
    return 0;
  } else {
    _is_fingerprint = 0;
  }

  _values.add(d);

  Accumulator<dscvalue>::extra(d);

  return 1;
}

int
Column_Data_Glaxo::allocate_space_for_this_many_rows(int nr)
{
  assert(nr > 0);

  return _values.resize(nr);
}

static int
dscvalue_comparitor_larger(const dscvalue* f1, const dscvalue* f2)
{
  if (*f1 > *f2) {
    return -1;
  }

  if (*f1 < *f2) {
    return 1;
  }

  return 0;
}

int
Column_Data_Glaxo::compute_dx()
{
  int n = _values.number_elements();
  if (n < 2) {
    cerr << "Column_Data_Glaxo::compute_dx: only " << n << " values, impossible\n";
    return 0;
  }

  _values.sort(dscvalue_comparitor_larger);

  // do we need to figure out where are the percentiles

  if (number_items_in_extremity_buckets < 0) {
    number_items_in_extremity_buckets =
        static_cast<int>(n * fraction_items_in_extremity_buckets);

    if (verbose) {
      cerr << number_items_in_extremity_buckets << " of " << n
           << " values will be in each extremity bucket\n";
    }

    if (0 == number_items_in_extremity_buckets) {
      cerr << "Column_Data_Glaxo::compute_dx: zero values in extremety bucket\n";
    }
  }

  if (0 == number_items_in_extremity_buckets)  // computed, but not active
  {
    _dx = (Accumulator<dscvalue>::maxval() - Accumulator<dscvalue>::minval()) /
          static_cast<dscvalue>(nbuckets);

    if (0.0 == _dx) {
      _ignore = 1;
    }

    return 1;
  }

  // Now we need to identify how far through the array we need to go in order to get the
  // appropriate number of items
  // TODO:ianwatson finish this sometime.

#ifdef IMPLEMENT_THIS
  int nfound = 1;

  dscvalue dprev = _values[0];

  for (int i = 0; i < n; i++) {
  }
#endif

  return 1;
}

int
Column_Data_Glaxo::bucket(dscvalue f) const
{
  if (global_range_active) {
    return global_bucket_value(f);
  }

  if (f <= _lower_bucket_boundary) {
    return 0;
  }

  if (f >= _upper_bucket_boundary) {
    return nbuckets;
  }

  return compute_bucket(_lower_bucket_boundary, _dx, f);
}

/*
  Generic routine for writing column info.
  Note that the invoker considers columns to start with column 1, so
  we increment the value written by 1 to accommodate that.
*/

int
Column_Data_Base::write_column_title(int c, std::ostream& os) const
{
  os << "column " << (c + 1);
  if (is_descriptor_file) {
    os << " '" << _column_id << "' ";
  }

  return os.good();
}

template <typename T>
int
establish_column_titles(const const_IWSubstring& buffer, T* column_data)
{
  assert(columns_in_input > 0);

  int columns_matched = 0;

  const_IWSubstring token;
  for (int i = 0, column = 0; buffer.nextword(token, i); column++) {
    T& c = column_data[column];

    c.set_id(token);

    if (c.ignore()) {
      continue;
    }

    if (global_range_active)  // if we are doing a global range, we can do it on only one
                              // descriptor, so do NOT do a regular expression match
    {
      if (!process_column_regexp) {
        columns_matched++;
      }
      if (iwstring::Equals(token, process_column_regexp->pattern())) {
        columns_matched++;
      } else {
        c.set_ignore(1);
      }
    } else if (!process_column_regexp) {
      ;
    } else if (!iwre2::RE2PartialMatch(token, *process_column_regexp)) {
      c.set_ignore(1);
    } else {
      columns_matched++;
      if (verbose) {
        c.write_column_title(column, cerr);
        cerr << " matches the process column regular expression\n";
      }
    }
  }

#ifdef ECHO_COLUMN_TITLES
  for (int i = 0; i < columns_in_input; i++) {
    T& c = column_data[i];
    c.write_column_title(i, cerr);
    cerr << " ignore? " << c.ignore() << '\n';
  }
#endif

  if (process_column_regexp && 0 == columns_matched) {
    cerr << "No columns match regular expression '" << process_column_regexp->pattern()
         << "'\n";
    return 0;
  }

  if (global_range_active && columns_matched > 1)  // must be duplicate column titles...
  {
    cerr << "Global range active, but " << columns_matched
         << " columns matched the regular expression\n";
    return 0;
  }

  return 1;
}

template <typename T>
int
_do_bucketise(iwstring_data_source& input, const T* column_data, std::ostream& output)
{
  IWString output_buffer;
  output_buffer.resize(input.longest_record() + 100);

  const_IWSubstring buffer;
  while (input.next_record(buffer) && output.good()) {
    //  cerr << "Read '" << buffer << "'\n";
    int i = 0;
    const_IWSubstring token;
    int column = 0;
    while (buffer.nextword(token, i))  // loop over each token in the record
    {
      const T& cd = column_data[column];

      //    cerr << "Column " << column << " IGNORE " << cd.ignore() << " FP? " <<
      //    cd.is_fingerprint() << '\n';

      if (cd.ignore() || cd.is_fingerprint() || missing_value == token) {
        if (suppress_output) {
          ;
        } else if (is_descriptor_file && 0 == column) {
          output_buffer << ' ' << token;
        } else if (cd.ignore() && only_write_non_ignored_columns) {
          ;
        } else {
          output_buffer << ' ' << token;
        }

        column++;
        continue;
      }

      column++;

      dscvalue d;
      if (!token.numeric_value(d)) {
        cerr << "Ignoring non numeric value '" << token << "' line " << input.lines_read()
             << " column " << column << '\n';
        continue;
      }

      int id = cd.bucket(d);
      //    cerr << "Value " << d << " goes in bucket " << id << '\n';

      if (id >= nbuckets) {
        id = nbuckets - 1;
      }

      bucket_count[id]++;

      if (!suppress_output) {
        output_buffer += text_bucket[id];
      }
    }

    if (!suppress_output) {
      output << output_buffer << '\n';
      output_buffer.resize_keep_storage(0);
    }
  }

  return output.good();
}

template <typename T>
int
do_global_bucketise(iwstring_data_source& input, T* column_data, std::ostream& output)
{
  Accumulator<dscvalue> g;  // over all selected columns

  for (int i = 0; i < columns_in_input; i++) {
    const T& a = column_data[i];

    if (a.ignore()) {
      continue;
    }

    g.extra(a.minval());
    g.extra(a.maxval());
  }

  if (g.minval() == g.maxval()) {
    cerr << "No variance in selected columns\n";
    return 0;
  }

  dscvalue dx = (g.maxval() - g.minval()) / static_cast<dscvalue>(nbuckets);

  for (int i = 0; i < columns_in_input; i++) {
    T& c = column_data[i];

    if (c.ignore() || c.is_fingerprint()) {
      continue;
    }

    c.set_output_specifications(g.minval(), dx);
  }

  return _do_bucketise(input, column_data, output);
}

template <typename T>
int
do_bucketise(iwstring_data_source& input, T* column_data, std::ostream& output)
{
  if (!global_range_active) {
    int rc = 1;

    for (int i = 0; i < columns_in_input; i++) {
      if (!column_data[i].compute_dx()) {
        cerr << "Cannot compute dx for column " << (i + 1) << '\n';
        rc = 0;
      }
    }

    if (0 == rc) {
      return 0;
    }
  }

  return _do_bucketise(input, column_data, output);
}

template <typename T>
int
report_column_statistics(T* column_data, int records_read, std::ostream& output)
{
  output << "Read " << records_read << " records, each with " << columns_in_input
         << " columns\n";
  if (verbose > 1) {
    for (int i = 0; i < columns_in_input; i++) {
      const T& c = column_data[i];
      if (c.ignore()) {
        continue;
      }

      c.debug_print(output);
    }
  }

  return output.good();
}

template <typename T>
int
process_record(const const_IWSubstring& buffer, T* column_data)
{
  int i = 0;
  int nwords_this_record = 0;
  int missing_values_this_record = 0;
  const_IWSubstring token;
  while (buffer.nextword(token, i))  // loop over each token in the record
  {
    if (nwords_this_record >= columns_in_input) {
      cerr << "Column count mismatch, should be only " << columns_in_input
           << " columns\n";
      return 0;
    }

    T& c = column_data[nwords_this_record];

    if (c.ignore()) {
      nwords_this_record++;
      continue;
    }

    if (missing_value == token) {
      missing_values_this_record++;
      nwords_this_record++;
      continue;
    }

    if (!c.extra(token)) {
      return 0;
    }

    nwords_this_record++;
  }

  if (nwords_this_record != columns_in_input) {
    cerr << "Column count mismatch, " << columns_in_input << " vs " << nwords_this_record
         << '\n';
    return 0;
  }

  missing_value_counter[missing_values_this_record]++;

  return 1;
}

template <typename T>
int
read_whole_file(iwstring_data_source& input, T* column_data)
{
  assert(columns_in_input > 0);

  const_IWSubstring buffer;

  cerr << "Reading the whole file\n";

  while (input.next_record(buffer)) {
    if (!process_record(buffer, column_data)) {
      cerr << "Fatal error processing line " << input.lines_read() << '\n';
      return 0;
    }
  }

  if (verbose) {
    report_column_statistics(column_data, input.lines_read(), cerr);
  }

  return 1;
}

template <typename T>
class IW_streampos_hash
{
 private:
 public:
  size_t operator()(const T&) const;
};

template <typename T>
size_t
IW_streampos_hash<T>::operator()(const T& s) const
{
  if (4 == sizeof(s)) {
    return static_cast<size_t>(s);
  }

  if (8 == sizeof(s)) {
    unsigned long i = static_cast<unsigned long>(s);
    return i;
  }

  return static_cast<size_t>(s);
}

template <typename T>
int
read_data_randomly(iwstring_data_source& input, T* column_data)
{
  assert(columns_in_input > 0);

  int nr = input.records_remaining();

  for (int i = 0; i < columns_in_input; i++) {
    if (!column_data[i].allocate_space_for_this_many_rows(nr)) {
      cerr << "Memory failure at column " << i << " of " << columns_in_input
           << " columns\n";
      return 0;
    }
  }

  std::streampos startpos = input.tellg();
  std::streampos filesize = input.file_size();

  int records_sampled = 0;

  int max_attempts = 5 * sample_records_randomly;
  int number_attempts = 0;

  std::unordered_set<off_t> visited;

  const_IWSubstring buffer;  // declared out here just for efficiency

  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<off_t> u(startpos, filesize - 11);  // 11 is just a number

  while (records_sampled < sample_records_randomly) {
    number_attempts++;
    if (number_attempts > max_attempts) {
      cerr << "Tried " << number_attempts << " to sample " << sample_records_randomly
           << " records, giving up...\n";
      break;
    }

    off_t pos = u(rng);

    input.seekg(pos);

    input.next_record(buffer);  // read what is likely to be a partial line

    pos = input.tellg();

    if (visited.find(pos) != visited.end()) {
      continue;
    }

    if (!input.next_record(buffer)) {
      continue;
    }

    visited.insert(pos);

    if (!process_record(buffer, column_data)) {
      cerr << "Fatal error processing record at " << pos << '\n';
      cerr << buffer << '\n';
      return 0;
    }
  }

  if (verbose) {
    report_column_statistics(column_data, records_sampled, cerr);
  }

  return 1;
}

template <typename T>
int
write_non_ignored_columns(const const_IWSubstring& buffer, const T* column_data,
                          std::ostream& output)
{
  int column = 0;
  const_IWSubstring token;

  IWString output_buffer;
  output_buffer.resize(buffer.length() + 1);

  int i = 0;
  while (buffer.nextword(token, i)) {
    if (0 == column && is_descriptor_file) {
      output_buffer << token;
    } else if (column_data[column].ignore()) {
      ;
    } else {
      output_buffer.append_with_spacer(token);
    }

    column++;
  }

  output_buffer << '\n';

  output << output_buffer;

  return output.good();
}

/*
  Note that if we have too many columns in a record, we will spill off
  the end of the counters array. they should have used tcount
*/

template <typename T>
int
bucketise(iwstring_data_source& input, T* column_data, std::ostream& output)
{
  for (int i = 0; i < ignore_columns_from_cl.number_elements(); i++) {
    int j = ignore_columns_from_cl[i];
    //  cerr << "Will ignore column " << j << '\n';
    column_data[j].set_ignore(1);
  }

  if (columns_to_process_from_cl.number_elements()) {
    for (int i = 0; i < columns_in_input; i++) {
      if (!columns_to_process_from_cl.contains(i)) {
        column_data[i].set_ignore(1);
      }
    }
  }

  const_IWSubstring buffer;
  for (int i = 0; i < nskip; i++) {
    input.next_record(buffer);

    if (is_descriptor_file && 0 == i) {
      if (!establish_column_titles(buffer, column_data)) {
        cerr << "Cannot establish column titles\n";
        return 0;
      }
    }
  }

  int rc;
  if (global_range_active)  // don't bother scanning the file
  {
    for (int i = 0; i < columns_in_input; i++) {
      column_data[i].set_is_fingerprint(0);
    }
    rc = 1;
  } else if (sample_records_randomly) {
    rc = read_data_randomly(input, column_data);
  } else {
    rc = read_whole_file(input, column_data);
  }

  if (0 == rc) {
    return 0;
  }

  (void)input.seekg(0);

  for (int i = 0; i < nskip; i++) {
    input.next_record(buffer);

    if (suppress_output) {
      ;
    } else if (only_write_non_ignored_columns) {
      write_non_ignored_columns(buffer, column_data, output);
    } else {
      output << buffer << '\n';
    }
  }

  if (global_bucketise) {
    rc = do_global_bucketise(input, column_data, output);
  } else {
    rc = do_bucketise(input, column_data, output);
  }

  return rc;
}

static int
bucketise(const char* fname, std::ostream& output)
{
  iwstring_data_source input(fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 99;
  }

  if (translate_tabs) {
    input.set_translate_tabs(1);
  }

  // Determine the number of tokens per record

  IWString buffer;
  if (!input.next_record(buffer)) {
    cerr << "Cannot read header from '" << fname << "'\n";
    return 61;
  }

  columns_in_input = buffer.nwords();
  if (verbose) {
    cerr << "Input contains " << columns_in_input << " columns\n";
  }

  if (0 == columns_in_input) {
    cerr << "Very bad news, header record contains zero columns\n";
    return 51;
  }

  if (!input.seekg(0)) {
    cerr << "Yipes, cannot seek back to beginning of file\n";
    return 13;
  }

  input.reset_record_counter();

  int rc;

  if (0.0 != fraction_items_in_extremity_buckets) {
    Column_Data_Glaxo* c = new Column_Data_Glaxo[columns_in_input];

    rc = bucketise(input, c, output);

    delete[] c;
  } else {
    Column_Data_Accumulator* c = new Column_Data_Accumulator[columns_in_input];

    rc = bucketise(input, c, output);

    delete[] c;
  }

  return rc;
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
  cerr << "Bucketises tabular data files\n";
  cerr << "By default, each column is bucketised into <nbuckets> buckets\n";
  cerr << "Usage: " << prog_name << " <options> <input_file>\n";
  cerr << "  -b <nbuckets>    number of buckets\n";
  cerr << "  -g               use dynamic range of all selected columns\n";
  cerr << "  -m <string>      specify missing value string (default '" << missing_value << "')\n";
  cerr << "  -N <records>     randomly sample records\n";
  cerr << '\n';
  cerr << "  -a               print column statistics (min, max, ave) and missing value stats\n";
  cerr << "  -j               is descriptor file (implies -s 1 -i 1)\n";
  cerr << "  -s <skip>        skip the first <skip> records of the file\n";
  cerr << "  -i <column>      specify columns to ignore (starting @ 1)\n";
  cerr << "  -c <col>         specify columns to process (starting @ 1)\n";
  cerr << "  -w               only write columns that are NOT ignored\n";
  cerr << "  -P <regexp>      only process columns which match <regexp>(implies -j)\n";
  cerr << "  -R min,max,dx    specify global range for binning\n";
  cerr << "  -R write=<fname> file for global distribution within global range\n";
  cerr << "  -R norm          normalise the data in the file to probabilities\n";
  cerr << "  -R nob           suppress normal output to stdout\n";
  cerr << "  -T               translate tabs to spaces\n";
  cerr << "  -v               verbose output\n";
  // clang-format on

  exit(rc);
}

#include "Foundational/cmdline/cmdline.h"

static int
allocate_buckets(int nbuckets)
{
  text_bucket = new IWString[nbuckets];
  bucket_count = new_int(nbuckets);

  if (nullptr == text_bucket || nullptr == bucket_count) {
    return 0;
  }

  for (int i = 0; i < nbuckets; i++) {
    text_bucket[i] << ' ' << i;
  }

  return 1;
}

static int
initialise_global_range(const const_IWSubstring& r)
{
  const_IWSubstring token;
  int i = 0;

  if (!r.nextword(token, i, ',')) {
    cerr << "Cannot extract first part of global range\n";
    return 0;
  }

  if (!token.numeric_value(global_min)) {
    cerr << "Invalid global minimum value '" << token << "'\n";
    return 0;
  }

  if (!r.nextword(token, i, ',')) {
    cerr << "Cannot extract second part of global range\n";
    return 0;
  }

  if (!token.numeric_value(global_max) || global_max <= global_min) {
    cerr << "Invalid global minimum value '" << token << "', min is " << global_min
         << '\n';
    return 0;
  }

  if (!r.nextword(token, i, ','))  // they just specified min,max - we compute global_dx
  {
    if (0 == nbuckets) {
      cerr << "Must specify number of buckets via the -b option, cannot continue\n";
      return 0;
    }

    global_dx = (global_max - global_min) / static_cast<dscvalue>(nbuckets);
    if (verbose) {
      cerr << "Range '" << r << "' with " << nbuckets
           << " buckets implies dx = " << global_dx << '\n';
    }
  } else if (nbuckets > 0) {
    cerr << "You have already specified " << nbuckets
         << " buckets, so cannot specify a delta in global range\n";
    return 0;
  } else if (!token.numeric_value(global_dx) || global_dx <= static_cast<dscvalue>(0.0) ||
             global_dx >= (global_max - global_min)) {
    cerr << "Invalid global dx value '" << token << "', range is "
         << (global_max - global_min) << '\n';
    return 0;
  } else {
    nbuckets = static_cast<int>((global_max - global_min) / global_dx + 0.001);
    if (verbose) {
      cerr << "Range '" << r << "' implies " << nbuckets << " buckets\n";
    }
  }

  global_range_active = 1;

  if (!allocate_buckets(nbuckets)) {
    return 0;
  }

  return 1;
}

static int
bucketise(int argc, char** argv)
{
  Command_Line cl(argc, argv, "b:gm:ajs:i:P:TvN:G:wR:c:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!cl.option_present('b') && !cl.option_present('R')) {
    cerr << "Must specify a bucket count via the -b option and/or -R options\n";
    //  usage (2);
  }

  if (cl.option_present('b')) {
    if (!cl.value('b', nbuckets) || nbuckets < 2) {
      cerr << "The number of buckets must be > 1\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will bucketise into " << nbuckets << " buckets\n";
    }

    if (!allocate_buckets(nbuckets)) {
      cerr << "Memory failure allocating " << nbuckets << " buckets\n";
      return 4;
    }
  }

  if (cl.option_present('G')) {
    int percentile_in_each_extremity_bucket;

    if (!cl.value('G', percentile_in_each_extremity_bucket) ||
        percentile_in_each_extremity_bucket < 1 ||
        percentile_in_each_extremity_bucket > 99) {
      cerr << "Must specify a valid percentile for the -G option\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will use the Glaxo binning, with the top and bottom "
           << percentile_in_each_extremity_bucket
           << " percentile in the extremety buckets\n";
    }

    fraction_items_in_extremity_buckets =
        static_cast<double>(percentile_in_each_extremity_bucket) / 100.0;
  }

  if (cl.option_present('N')) {
    if (!cl.value('N', sample_records_randomly) || sample_records_randomly < 2) {
      cerr << "The sample records randomly option (-N) must be followed by a value whole "
              "number\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Will determine ranges by randomly sampling " << sample_records_randomly
           << " records\n";
    }
  }

  if (cl.option_present('T')) {
    translate_tabs = 1;
    if (verbose) {
      cerr << "Tabs translated to spaces\n";
    }
  }

  if (cl.option_present('a')) {
    print_statistics = 1;
    if (verbose) {
      cerr << "Will print column statistics\n";
    }
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

  if (cl.option_present('R') && cl.option_present('g')) {
    cerr << "The global range (-R) and global bucketise (-g) options are mutually "
            "exclusive\n";
    usage(4);
  }

  if (cl.option_present('g')) {
    global_bucketise = 1;
    if (verbose) {
      cerr << "Will bucketise according to the dynamic range of all selected columns\n";
    }
  }

  std::ofstream stream_for_global_profile;
  int normalise_global_profile = 0;

  if (cl.option_present('R')) {
    IWString fname;

    const_IWSubstring r;
    int i = 0;
    while (cl.value('R', r, i++)) {
      if (r.starts_with("write=")) {
        fname = r;
        fname.remove_leading_chars(6);
      } else if (r.starts_with("norm")) {
        normalise_global_profile = 1;
        if (verbose) {
          cerr << "Will normalise the global profile to probabilities\n";
        }
      } else if (r.starts_with("nob")) {
        suppress_output = 1;
      } else if (initialise_global_range(r)) {
        ;
      } else {
        cerr << "Invalid or unrecognised -R qualifier '" << r << "'\n";
        usage(4);
      }
    }

    if (0 == fname.length()) {
      cerr << "No file specified for global profile\n";
      usage(6);
    }

    stream_for_global_profile.open(fname.null_terminated_chars(), std::ios::out);
    if (!stream_for_global_profile.good()) {
      cerr << "Could not open global profile file '" << fname << "'\n";
      return 5;
    }

    if (verbose) {
      cerr << "Global profile written to '" << fname << "'\n";
    }
  }

  if (cl.option_present('w')) {
    only_write_non_ignored_columns = 1;

    if (verbose) {
      cerr << "Will only write non-ignored columns\n";
    }
  }

  if (cl.option_present('c') && cl.option_present('i')) {
    cerr << "The -c and -i options are mutually exclusive\n";
    usage(5);
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

  if (cl.option_present('c')) {
    int i = 0;
    const_IWSubstring c;
    while (cl.value('c', c, i++)) {
      int col;
      if (!c.numeric_value(col) || col < 1) {
        cerr << "Invalid column '" << c << "'\n";
        return 8;
      }

      columns_to_process_from_cl.add(col - 1);
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');

    if (!iwre2::RE2Reset(process_column_regexp, p)) {
      cerr << "Cannot parse process column regexp '" << p << "'\n";
      usage(4);
    }

    is_descriptor_file = 1;
    if (0 == nskip) {
      nskip = 1;
    }
    if (0 == ignore_columns_from_cl.number_elements()) {
      ignore_columns_from_cl.add(0);
    }

    if (verbose) {
      cerr << "Only columns matching '" << process_column_regexp->pattern()
           << "' will be processed\n";
    }
  }

  if (sample_records_randomly && cl.number_elements() > 1) {
    cerr << "Sorry, can only sample records randomly with one input file\n";
    return 4;
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(3);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!bucketise(cl[i], std::cout)) {
      rc = i + 1;
      break;
    }
  }

  if (verbose && print_statistics) {
    for (int i = 0; i < missing_value_counter.number_elements(); i++) {
      if (missing_value_counter[i]) {
        cerr << missing_value_counter[i] << " records had " << i << " missing values\n";
      }
    }
  }

  int items_bucketised = sum_vector(bucket_count, nbuckets);

  if (0 == items_bucketised) {
    cerr << "NO items bucketised, cannot continue\n";
    return 4;
  }

  if (verbose) {
    cerr << "Bucket populations\n";

    for (int i = 0; i < nbuckets; i++) {
      cerr << bucket_count[i] << " items in bucket " << i;
      if (global_range_active) {
        cerr << ' ' << (global_min + i * global_dx);
      }
      cerr << '\n';
      ;
    }

    cerr << items_bucketised << " items bucketised\n";

    if (items_below_global_range) {
      cerr << "Warning, " << items_below_global_range << " were below the global min "
           << global_min << '\n';
    }
    if (items_above_global_range) {
      cerr << "Warning, " << items_above_global_range << " were above the global max "
           << global_max << '\n';
    }
  }

  if (stream_for_global_profile.rdbuf()->is_open()) {
    float divide_by;
    if (normalise_global_profile) {
      divide_by = static_cast<float>(items_bucketised);
    } else {
      divide_by = static_cast<float>(1.0);
    }

    if (verbose) {
      cerr << "Writing " << nbuckets << " buckets of data\n";
    }

    for (int i = 0; i < nbuckets; i++) {
      stream_for_global_profile << (global_min + i * global_dx) << ' ';

      if (normalise_global_profile) {
        stream_for_global_profile << (bucket_count[i] / divide_by);
      } else {
        stream_for_global_profile << bucket_count[i];
      }

      stream_for_global_profile << '\n';
    }

    stream_for_global_profile.close();
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = bucketise(argc, argv);

  return rc;
}
