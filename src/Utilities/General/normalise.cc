#include <stdlib.h>
#include <iostream>
#include <stdint.h>

//using namespace std;

using std::cerr;
using std::endl;

#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "cmdline.h"
#include "accumulator.h"
#include "iwstring_data_source.h"
#include "misc.h"
#include "minmaxspc.h"
#include "iw_stl_hash_map.h"

const char * prog_name = NULL;

static int verbose = 0;
static IWString missing_value('.');
static int header_records_to_skip = 0;
static int is_descriptor_file = 0;
static int initial_columns_to_skip = 0;

#define TAKE_LOG 1

/*
  How can we "normalise" the input.
*/

#define NRML_MIN_TO_MAX 1
#define NRML_UNIT_VARIANCE 2
#define NRML_SPREAD_ZERO 3
#define NRML_MEAN_CENTER 4
#define NRML_0255 5

static int scaling_type = NRML_MIN_TO_MAX;

static int perform_unscaling_operation = 0;

static int allow_out_of_range_unscalings = 0;

static int truncate_low_values = 0;
static int truncate_high_values = 0;

static int allow_unscaled_values_below_zero = 1;

static int operator_to_apply = 0;

static int columns_in_input = 0;

static int write_results_as_float = 0;
static int output_precision = 0;

static double nearest_integer = 0.0;

static double range_min = 0.0;
static double range_max = 1.0;

static double constant_offset = 0.0;

static bool user_specified_range_max = false;
static double user_specified_range_max_value = 0.0;

static int report_out_of_range_values = 1;

static int values_out_of_range = 0;

/*
  Aug 2006
  I have scaling files with a single column, but the header record does
  not match. We have a kludge for that case
*/

static int single_column_kludge = 0;

static IWString name_of_existing_scaling_file_to_use;

static int translate_descriptor_names_to_lowercase = 0;

static char input_separator = ' ';
static char output_separator = ' ';

static int suppress_zero_variance_columns = 1;

static int categorical_variables_possibly_present = 0;
static extending_resizable_array<int> categorical_variables_found_in_column;

/*
  Nov 2006
  I have descriptor files where the first column is the ID, 
  the 2nd column is the activity and the remaining columns are
  data.

  We apply one kind of scaling to the response (column 2) and
  another to the rest.
*/

class NColumn : private Accumulator<double>
{
  private:
    IWString _descriptor_name;

    int _skip;

    int _is_constant;        // useful if we are not suppressing zero variance columns

    int _missing_values_encountered;

//  To avoid recomputing averages and variances, we store them

    double _average;
    double _variance;
    double _range;

//  Since we may read this info from a file, we need to replicate all
//  the values in the accumulator

    double _minval;
    double _maxval;
    int    _n;

  public:
    NColumn();

    void set_descriptor_name (const const_IWSubstring & d) { _descriptor_name = d;}
    const IWString & descriptor_name() const { return _descriptor_name;}

    int skip() const { return _skip;}
    void set_skip (int s) { _skip = s;}

    int n() const { return Accumulator<double>::n();}

    void extra (double f);

    void extra_missing_value() { _missing_values_encountered++;}

    int establish_ranges();
    int establish_range_from_pre_existing_data(const IWString & buffer);

    int scale (double, double &) const;
    int unscale (double, double &) const;

    int report (int, std::ostream &) const;

    int write_scaling_information(int col, IWString_and_File_Descriptor & output) const;
};

static NColumn * column = NULL;

NColumn::NColumn()
{
  _skip = 0;
  _range = 0.0;
  _average = 0.0;
  _variance = 0.0;

  _is_constant = 0;

  _missing_values_encountered = 0;

  return;
}

int
NColumn::report (int col,
                 std::ostream & output) const
{
  output << "Column " << col;

  if (_descriptor_name.length())
    output << ' ' << _descriptor_name;

  output << ' ';

  if (_missing_values_encountered)
    output << _missing_values_encountered << " missing values";

  const int n = Accumulator<double>::n();

  output << n << " values";
  if (0 == n)
  {
    output << '\n';
    return output.good();
  }

  if (_is_constant)
    output << " constant " << Accumulator<double>::minval();
  else
    output << " between " << Accumulator<double>::minval() << " and " << Accumulator<double>::maxval() << " ave " << Accumulator<double>::average() << " var " << Accumulator<double>::variance();

  output << '\n';

  return output.good();
}

int 
NColumn::write_scaling_information(int col,
                                   IWString_and_File_Descriptor & output) const
{
  if (_skip)
    return 1;

  int n = Accumulator<double>::n();

  if (0 == n)
    return 1;

  if (_descriptor_name.length())
    output << _descriptor_name;
  else
    output << "Column " << col;

  output << ' ';

  assert (n > 1);

  output << Accumulator<double>::minval() << ' ' << Accumulator<double>::maxval() << ' ' << Accumulator<double>::average() << ' ';

  if (_is_constant)
    output << '0';
  else
    output << sqrt(Accumulator<double>::variance());
    
  output << " # min max average variance\n";

  return 1;
}

void
NColumn::extra (double f)
{
  Accumulator<double>::extra(f);

  return;
}

int
NColumn::establish_ranges()
{
  if (0 == _n)
    _n = Accumulator<double>::n();

  assert (_n > 1);

  _range = Accumulator<double>::maxval() - Accumulator<double>::minval();

  if (0.0 == _range)
  {
    cerr << "NColumn::establish_ranges: '" << _descriptor_name << "' no variability in column\n";
    _is_constant = 1;
    if (suppress_zero_variance_columns)
    {
      _skip = 1;
      return 0;
    }
  }

  _minval = Accumulator<double>::minval();
  _maxval = Accumulator<double>::maxval();
  _n = Accumulator<double>::n();

  if (NRML_UNIT_VARIANCE == scaling_type)
  {
    _average = Accumulator<double>::average();
    _variance = sqrt(Accumulator<double>::variance());
  }
  else if (NRML_MEAN_CENTER == scaling_type)
    _average = Accumulator<double>::average();

  _skip = 0;

  return 1;
}

template <typename S1, typename S2>
int
do_nextword (const S1 & buffer,
             int & i,
             S2 & token)
{
  if (' ' != input_separator)
    return buffer.nextword_single_delimiter(token, i, input_separator);
  else
    return buffer.nextword(token, i);
}

int
NColumn::establish_range_from_pre_existing_data(const IWString & buffer)
{
  if (4 != buffer.nwords())
  {
    cerr << "NColumn::establish_ranges:must have exactly 4 words\n";
    cerr << buffer << "\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  do_nextword(buffer, i, token);

  if (! token.numeric_value(_minval))
  {
    cerr << "NColumn::establish_ranges:invalid _minval\n";
    return 0;
  }

  do_nextword(buffer, i, token);

  if (! token.numeric_value(_maxval))
  {
    cerr << "NColumn::establish_ranges:invalid _maxval\n";
    return 0;
  }

  do_nextword(buffer, i, token);

  if (! token.numeric_value(_average))
  {
    cerr << "NColumn::establish_ranges:invalid _ave\n";
    return 0;
  }

  do_nextword(buffer, i, token);

  if (! token.numeric_value(_variance))
  {
    cerr << "NColumn::establish_ranges:invalid _variance\n";
    return 0;
  }

  if (_minval < _maxval)
    ;
  else if (_minval > _maxval)
  {
    cerr << "NColumn::establish_ranges:invalid range specification " << _minval << " to " << _maxval << '\n';
    return 0;
  }
  else
  {
    _is_constant = 1;
  }

  _range = _maxval - _minval;

  return 1;
}

int
NColumn::scale (double f,
                double & rc) const
{
  if (_is_constant)
    rc = f;
  else if (NRML_MIN_TO_MAX == scaling_type || NRML_SPREAD_ZERO == scaling_type)
  {
    if (truncate_low_values && f < _minval)
      f = _minval;
    else if (truncate_high_values && f > _maxval)
      f = _maxval;

    rc = range_min + (f - _minval) / _range * (range_max - range_min);

//  if (rc < 0.0)
//    cerr <<"Negative value from " << f << " _minval " << _minval << " _maxval " << _maxval << ", got " << rc << " truncate_low_values " << truncate_low_values << endl;
  }
  else if (NRML_UNIT_VARIANCE == scaling_type)
    rc = range_min + (f - _average) / _variance * (range_max - range_min);
  else if (NRML_MEAN_CENTER == scaling_type)
    rc = f - _average;
  else if (NRML_0255 == scaling_type)
  {
    if (truncate_low_values && f < _minval)
      f = _minval;
    else if (truncate_high_values && f > _maxval)
      f = _maxval;
    rc = static_cast<int>((f - _minval) / (_maxval - _minval) * 255.0 + 0.49999);
  }
  else
  {
    cerr << "NColumn::scale: what kind of scaling " << scaling_type << endl;
    abort();
  }

  rc += constant_offset;

  return 1;
}

int
NColumn::unscale (double f,
                  double & rc) const
{
  if (NRML_MIN_TO_MAX == scaling_type || NRML_SPREAD_ZERO == scaling_type)
  {
    if (_is_constant)   // just echo what came in
      rc = f;
    else if (f >= range_min && f <= range_max)     // in range, great
      rc = _minval + (f - range_min) / (range_max - range_min) * (_maxval - _minval);
    else
    {
      values_out_of_range++;
      if (report_out_of_range_values)
      {
        cerr << "NColumn::unscale:value out of range " << f << " must be " << range_min << " to " << range_max;
        if (_descriptor_name.length() > 0)
          cerr << ' ' << _descriptor_name;
        cerr << endl;
      }

//    cerr << "low " << truncate_low_values << " high " << truncate_high_values << " " << (truncate_low_values && f < range_min) << endl;

      if (truncate_low_values && f < range_min) 
        f = range_min;
      else if (truncate_high_values && f > range_max)
        f = range_max;
      else if (allow_out_of_range_unscalings)
        ;
      else
        return 0;

//    cerr << "Truncated to " << f << endl;

      rc = _minval + (f - range_min) / (range_max - range_min) * (_maxval - _minval);
    }
  }
  else if (NRML_UNIT_VARIANCE == scaling_type)
    rc = _average + (f - range_min) / (range_max - range_min) * _variance;
  else if (NRML_MEAN_CENTER == scaling_type)
    rc += _average;
  else if (NRML_0255 == scaling_type)
    rc = _minval + (f - range_min) / (range_max - range_min) * (_maxval - _minval);
  else
  {
    cerr << "NColumn::unscale: what kind of scaling " << scaling_type << endl;
    abort();
  }

  if (constant_offset != 0.0)
    rc += constant_offset;

  return 1;
}

/*
  Tokens are read in
*/

class Data_Item : public Set_or_Unset<double>
{
  private:
    const_IWSubstring _text_rep;    // we rely on the underlying buffer staying in scope

  public:
    Data_Item();

    int build (const const_IWSubstring &, const int col);

    void set_text_rep (const const_IWSubstring & t) { _text_rep = t;}

    const const_IWSubstring & text_representation() const { return _text_rep;}
};

static Data_Item * d = NULL;

Data_Item::Data_Item()
{
}

int
Data_Item::build (const const_IWSubstring & token, const int col)
{
  _text_rep = token;

  if (missing_value == token)
  {
    unset();
    return 1;
  }

  double xx;
  if (! token.numeric_value(xx))
  {
    cerr << "Data_Item::build:non numeric value ('" << token << "') column " << col << "\n";
    unset();
    return 0;
  }

  set(xx);

  return 1;
}

/*
  For a given value of scaling_type, we need to set the corresponding
  values for range_min and range_max
*/

static int
initialise_ranges_for_scaling()
{
  if (NRML_UNIT_VARIANCE == scaling_type)
  {
    range_min =  0.0;
    range_max =  1.0;
  }
  else if (NRML_SPREAD_ZERO == scaling_type)
  {
    range_min = -1.0;
    range_max =  1.0;
  }
  else if (NRML_MIN_TO_MAX == scaling_type)
  {
    range_min = 0.0;
    if (user_specified_range_max)
      range_max = user_specified_range_max_value;
    else
      range_max = 1.0;
  }
  else if (NRML_MEAN_CENTER == scaling_type)   // min and max don't matter
    ;
  else if (NRML_0255 == scaling_type)
  {
    range_min = 0.0;
    range_max = 255.0;
  }
  else
  {
    cerr << "initialise_ranges_for_scaling:what kind of scaling is this " << scaling_type << endl;
    return 0;
  }

  return 1;
}

static int
determine_tokens_per_line_and_allocate_arrays (iwstring_data_source & input)
{
  const_IWSubstring buffer;
  if (! input.next_record(buffer))
  {
    cerr << "Cannot read header record\n";
    return 0;
  }

  if (' ' == input_separator)
    columns_in_input = buffer.nwords();
  else 
    columns_in_input = buffer.nwords_single_delimiter(input_separator);

  if (0 == columns_in_input)
  {
    cerr << "Fatal error, no columns in input\n";
    return 0;
  }

  if (verbose)
    cerr << "Input contains " << columns_in_input << " columns\n";

  d = new Data_Item[columns_in_input];

  column = new NColumn[columns_in_input];

  if (NULL == column)
  {
    cerr << "Cannot allocate " << columns_in_input << " columns\n";
    return 0;
  }

  if (initial_columns_to_skip >= columns_in_input)
  {
    cerr << "Cannot skip " << initial_columns_to_skip << " columns, only " << columns_in_input << " in file\n";
    return 0;
  }

  int ndx = 0;
  IWString token;

  for (int i = 0; i < initial_columns_to_skip; i++)
  {
    column[i].set_skip(1);
    do_nextword(buffer, ndx, token);
  }

  for (int i = initial_columns_to_skip; i < columns_in_input; i++)
  {
    do_nextword(buffer, ndx, token);
    if (translate_descriptor_names_to_lowercase)
      token.to_lowercase();

    column[i].set_descriptor_name(token);
  }

  if (! input.push_record())
  {
    cerr << "Gack, cannot push record back\n";
    return 0;
  }

  return columns_in_input;
}

int
read_record (iwstring_data_source & input,
             Data_Item * d, 
             int & fatal)
{
  const_IWSubstring buffer;
  if (! input.next_record(buffer))
  {
    fatal = 0;
    return 0;
  }

  int nw;
  if (' ' == input_separator)
    nw = buffer.nwords();
  else
    nw = buffer.nwords_single_delimiter(input_separator);

  if (nw != columns_in_input)
  {
    cerr << "Column count mismatch, found " << nw << " expected " << columns_in_input << endl;
    fatal = 1;
    return 0;
  }

  int i = 0;
  const_IWSubstring token;
  int col = 0;

  while (do_nextword(buffer, i, token))
  {
    Data_Item & di = d[col];

//  cerr << "Reading column " << col << " skip = " << column[col].skip() << " i = " << i << endl;

    if (column[col].skip())
      di.set_text_rep(token);
    else if (categorical_variables_found_in_column[col] > 10)
      di.set_text_rep(token);
    else if (di.build(token, col))
      ;
    else if (categorical_variables_possibly_present)
      categorical_variables_found_in_column[col]++;
    else
    {
      cerr << "Invalid numeric '" << token << "'\n";
      fatal = 1;
      return 0;
    }

    col++;
  }

  return 1;
}

static int
do_write_record(const const_IWSubstring & buffer,
                IWString_and_File_Descriptor & output)
{
  int ndx = 0;
  const_IWSubstring token;
  int col = 0;

  for (int i = 0; i < initial_columns_to_skip; i++, col++)
  {
    do_nextword(buffer, ndx, token);

    output.append_with_spacer(token, output_separator);
  }

  while (do_nextword(buffer, ndx, token))
  {
    if (! column[col].skip())
      output.append_with_spacer(token);

    col++;
  }

  output << '\n';

  return 1;
}

static int
retrieve_scaled_value (const NColumn & col,
                       const Data_Item & d,
                       double & scaled)
{
  double x;

  if (! d.value(x))
    return 0;

  if (perform_unscaling_operation)
    return col.unscale(x, scaled);
  else
    return col.scale(x, scaled);
} 

static int
parse_scaling_info (const const_IWSubstring & buffer,
                    IW_STL_Hash_Map_String & xref)
{
  if (buffer.starts_with("Column "))
  {
    cerr << "Implement Column handling sometime\n";
    abort();
    return 0;
  }

// First token is descriptor name, rest is info

  if (5 != buffer.nwords())
  {
    cerr << "Descriptor scaling info must have exactly 5 tokens, min, max, ave, variance\n";
    cerr << "'" << buffer << "'\n";
    return 0;
  }

  IWString dname, zrest;

  buffer.split(dname, ' ', zrest);

  if (translate_descriptor_names_to_lowercase)
    dname.to_lowercase();

  xref[dname] = zrest;

  return 1;
}

static int
read_scaling_file (iwstring_data_source & input,
                   IW_STL_Hash_Map_String & xref)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (buffer.starts_with("#NTYPE "))
    {
      buffer.remove_leading_words(1);
      if ("UV" == buffer)
      {
        scaling_type = NRML_UNIT_VARIANCE;
        if (verbose)
          cerr << "Scaling information file scaling type UV\n";
      }
      else if ("11" == buffer)
      {
        scaling_type = NRML_SPREAD_ZERO;
        if (verbose)
          cerr << "Scaling information file scaling type 11\n";
      }
      else if ("01" == buffer)
      {
        scaling_type = NRML_MIN_TO_MAX;
        if (verbose)
          cerr << "Scaling information file scaling type 01\n";
      }
      else if ("MC" == buffer)
      {
        scaling_type = NRML_MEAN_CENTER;
        if (verbose)
          cerr << "Scaling information file scaling type mean center\n";
      }
      else if ("0255" == buffer)
      {
        scaling_type = NRML_0255;
        if (verbose)
          cerr << "Scaling will be from 0 to 255\n";
      }
      else
      {
        cerr << "Invalid #NTYPE directive in scaling file '" << buffer << "'\n";
        return 0;
      }

      continue;
    }
    else if (buffer.starts_with("#UMAX"))
    {
      buffer.remove_leading_words(1);
      if (! buffer.numeric_value(user_specified_range_max_value))
      {
        cerr << "read_scaling_file:invalid UMAX directive " << buffer << endl;
        return 0;
      }
      user_specified_range_max = 1;
      continue;
    }
    else if (buffer.starts_with("#OFFSET"))
    {
      buffer.remove_leading_words(1);
      if (! buffer.numeric_value(constant_offset))
      {
        cerr << "read_scaling_file:invalid OFFSET directive " << buffer << endl;
        return 0;
      }
      continue;
    }

    if (buffer.starts_with('#'))
      continue;

    buffer.truncate_at_first('#');
    buffer.strip_trailing_blanks();

    if (! parse_scaling_info(buffer, xref))
    {
      cerr << "Invalid scaling information record '" << buffer << "'\n";
      return 0;
    }
  }

  return xref.size();
}

static int
read_scaling_file (const const_IWSubstring & fname,
                   IW_STL_Hash_Map_String & xref)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open scaling file '" << fname << "'\n";
    return 0;
  }

  return read_scaling_file(input, xref);
}

static int
establish_scaling_from_pre_existing_file (IWString & name_of_existing_scaling_file_to_use)
{
  assert (columns_in_input > 0);

  IW_STL_Hash_Map_String xref;

  if (! read_scaling_file(name_of_existing_scaling_file_to_use, xref))
  {
    cerr << "Cannot read existing scaling file '" << name_of_existing_scaling_file_to_use << "'\n";
    return 0;
  }

  if (single_column_kludge && 1 == xref.size() && 1 == (columns_in_input - initial_columns_to_skip))
  {
    IW_STL_Hash_Map_String::const_iterator f = xref.begin();

    int col = columns_in_input - initial_columns_to_skip;

    if (! column[col].establish_range_from_pre_existing_data((*f).second))
    {
      cerr << "Invalid range specification for descriptor '" << column[col].descriptor_name() << "'\n";
      return 0;
    }

    return 1;
  }

  int columns_being_written = 0;

  for (int i = initial_columns_to_skip; i < columns_in_input; i++)
  {
    const IWString & dname = column[i].descriptor_name();

    IW_STL_Hash_Map_String::const_iterator f = xref.find(dname);

    if (f == xref.end())
    {
      column[i].set_skip(1);
      continue;
    }

    if (! column[i].establish_range_from_pre_existing_data((*f).second))
    {
      cerr << "Invalid range specification for descriptor '" << dname << "'\n";
      return 0;
    }

    columns_being_written++;
  }

  if (0 == columns_being_written)
  {
    cerr << "normalise:no columns selected for output\n";
    return 0;
  }
  
  return 1;
}

/*
  Fills in the array of accumulators according to the columns in INPUT.
  INPUT is assumed to have DATA_ITEMS tokens on each record.
  We skip the first token in each record.
  We return 1 if successful.
*/

static int
determine_profile (iwstring_data_source & input)
{
  assert (input.good());

  if (! determine_tokens_per_line_and_allocate_arrays(input))
  {
    cerr << "Cannot determine number of columns in input\n";
    return 0;
  }

  assert (columns_in_input > 0);

  const_IWSubstring buffer;

  for (int i = 0; i < header_records_to_skip; i++)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Cannot skip " << header_records_to_skip << " header records\n";
      return 0;
    }
  }

  int lines_read = 0;
  int fatal;
  while (read_record(input, d, fatal))
  {
    lines_read++;

    for (int i = initial_columns_to_skip; i < columns_in_input; i++)
    {
      double v;
      if (d[i].value(v))
        column[i].extra(v);
      else
        column[i].extra_missing_value();
    }
  }

  if (fatal)
  {
    cerr << "Error reading data, read " << input.lines_read() << '\n';
    return 0;
  }

  if (verbose)
    cerr << "Read " << lines_read << " records\n";

  if (0 == lines_read)
    return 0;

  assert (columns_in_input > 0);

  int columns_being_written = 0;

  for (int i = 0; i < columns_in_input; i++)
  {
    if (column[i].skip())
      continue;

    if (categorical_variables_found_in_column[i])
      columns_being_written++;
    else if (column[i].establish_ranges())
      columns_being_written++;
  }

  if (0 == columns_being_written)
  {
    cerr << "No columns can be output\n";
    return 0;
  }

  if (verbose)
    cerr << "Writing " << columns_being_written << " active columns\n";

  return 1;
}

static void
maybe_do_nearest_integer (double v, 
                          IWString_and_File_Descriptor & output,
                          const int output_precision)
{
  if (fabs(v) < nearest_integer)
  {
    output << '0';
    return;
  }

  if (v > 0.0)
  {
    int64_t i = static_cast<int64_t>(v + 0.4999);
    if (fabs(v - i) < nearest_integer)
    {
      output << i;
      return;
    }
  }
  else
  {
    int64_t i = static_cast<int64_t>(v - 0.4999);
    if (fabs(v - i) < nearest_integer)
    {
      output << i;
      return;
    }
  }

// just regular output

  if (output_precision > 0)
    output.append_number(v, output_precision);
  else if (write_results_as_float)
    output << static_cast<float>(v);
  else
    output << v;

  return;
}

/*
  The input file has been profiled, actually do it
*/

static int
do_normalise (iwstring_data_source & input, 
              IWString_and_File_Descriptor & output)
{
  assert (input.ok());
  assert (output.good());

  if (name_of_existing_scaling_file_to_use.length())
  {
    if (! determine_tokens_per_line_and_allocate_arrays(input))
      return 0;

    if (! establish_scaling_from_pre_existing_file(name_of_existing_scaling_file_to_use))
    {
      cerr << "Cannot establish scaling from existing scaling file\n";
      return 0;
    }
  }

// Each different scaling type has a min and max associated with it

  initialise_ranges_for_scaling();

  const_IWSubstring buffer;

  for (int i = 0; i < header_records_to_skip; i++)
  {
    input.next_record(buffer);

    do_write_record(buffer, output);
  }

  int items_read = 0;
  int missing_values = 0;
  int records = 0;
  int fatal;

  while (input.good() && output.good())
  {
    if (! read_record(input, d, fatal))
    {
      if (fatal)
        return 0;
      break;
    }

    records++;

    int need_space = 0;

    for (int i = 0; i < initial_columns_to_skip; i++)
    {
      if (need_space)
        output << output_separator;

      output << d[i].text_representation();
      need_space = 1;
    }

    for (int i = initial_columns_to_skip; i < columns_in_input; i++)
    {
      if (column[i].skip())
        continue;

      if (need_space)
        output << output_separator;

      need_space = 1;

      double scaled;

      if (categorical_variables_found_in_column[i])
      {
        output << d[i].text_representation();
        continue;
      }

      if (! retrieve_scaled_value(column[i], d[i], scaled))
      {
        missing_values++;
        output << missing_value;
        continue;
      }

      items_read++;
      if (nearest_integer > 0.0)
        maybe_do_nearest_integer(scaled, output, output_precision);
      else if (write_results_as_float)
        output << static_cast<float>(scaled);
      else if (output_precision)
        output.append_number(scaled, output_precision);
      else
        output << scaled;
    }

    output << "\n";

    output.write_if_buffer_holds_more_than(8192);
  }

  if (verbose)
  {
    cerr << items_read << " items read, " << records << " records written\n";
    if (missing_values)
      cerr << missing_values << " missing values\n";
  }

  return 1;
}
 
static int
normalise (iwstring_data_source & input,
           IWString_and_File_Descriptor & output)
{
  assert (input.ok());

  assert (output.good());

  if (name_of_existing_scaling_file_to_use.length() > 0)   // no need to profile
    ;
  else if (! determine_profile(input))
  {
    cerr << "Cannot profile input file\n";
    return 0;
  }
  else if (! input.seekg(0))   // reset file pointer
  {
    cerr << "Bad news, cannot seek back to start of input\n";
    return 0;
  }
  else
    input.reset_record_counter();

  cerr.setf(std::ios::showpoint);

  if (verbose && 0 == name_of_existing_scaling_file_to_use.length())    // we have profiled the input file
  {
    for (int i = 0; i < columns_in_input; i++)
    {
      NColumn & ci = column[i];

      if (! ci.skip())
        ci.report(i, cerr);
    }
  }

  return do_normalise(input, output);
}

static int
normalise (const char * input_fname,
           IWString_and_File_Descriptor & output)
{
  assert (NULL != input_fname);

  iwstring_data_source input(input_fname);
  if (! input.good())
  {
    cerr << prog_name << " cannot open '" << input_fname << "'\n";
    return 0;
  }

  input.set_dos(1);

  if ('\t' != input_separator)    // preserve historical behaviour which always translated tabs
    input.set_translate_tabs(1);

  return  normalise(input, output);
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Normalises columns\n";
  cerr << " -m <string>        set missing value specifier to <string>\n";
  cerr << " -s <records>       discard <records> from the top of the file\n";
  cerr << " -i <number>        ignore the first <number> columns in each record\n";
  cerr << " -j                 treat as a descriptor file\n";
//cerr << " -O <operator>      apply operator <operator> to all values\n";
  cerr << " -N 01              scale so data in each column is in the range [0,1] (default)\n";
  cerr << " -N 11              scale so data in each column is in the range [-1,1]\n";
  cerr << " -N uv              scale so columns have unit variance\n";
  cerr << " -N mc              scale so columns mean centered\n";
  cerr << " -x <x>             add <x> to all values before output (not saved to -C file)\n";
  cerr << " -X <x>             with 01 scaling, use <x> as max value rather than 1.0\n";
  cerr << " -p f               write output as float rather than double\n";
  cerr << " -p <n>             output precision\n";
  cerr << " -C <fname>         write scaling information to <fname>\n";
  cerr << " -U <fname>         use scaling information to previously created by a -C invocation\n";
  cerr << " -u                 perform a reverse normalisation - convert normalised back to un-normalised\n";
  cerr << " -o allow           allow out of range reverse normalisations\n";
  cerr << " -o trunc           truncate all out of range reverse normalisations\n";
  cerr << " -o trunclow,trunchigh   truncate only low or high reverse normalisations\n";
  cerr << " -b <sep>           input separator (tab, comma, space)";
  cerr << " -q                 do not report out of range errors\n";
  cerr << " -c                 translate all descriptor names to lowercase\n";
  cerr << " -z                 retain zero variance columns\n";
  cerr << " -R <float>         round output values to the nearest integer if within <float> of that nearest int\n";
  cerr << " -A                 categorical variables may be present - beware, may hide erroneous input\n";
  cerr << " -v                 verbose output\n";

  exit(rc);
}

static int
write_scaling_file(IWString_and_File_Descriptor & output,
                   const char * input_fname)
{
  output << "#Scaling file written by normalise, from '" << input_fname << "'\n";
  output << "#NTYPE ";
  if (scaling_type == NRML_UNIT_VARIANCE)
    output << "UV\n";
  else if (scaling_type == NRML_SPREAD_ZERO)
    output << "11\n";
  else if (scaling_type == NRML_MIN_TO_MAX)
    output << "01\n";
  else if (scaling_type == NRML_MEAN_CENTER)
    output << "MC\n";
  else if (NRML_0255 == scaling_type)
    output << "0255\n";
  else
  {
    cerr << "Need to know scaling for -C file\n";
    output << "??\n";
  }

  if (user_specified_range_max_value > 0.0)
    output << "#UMAX " << user_specified_range_max_value << '\n';

  if (constant_offset != 0.0)
    output << "#OFFSET " << constant_offset << '\n';

  if (verbose > 1)
    cerr << "Writing scaling data for " << columns_in_input << " columns\n";

  output << "#N = " << column[initial_columns_to_skip].n() << '\n';

  for (int i = initial_columns_to_skip; i < columns_in_input; i++)
  {
    column[i].write_scaling_information(i, output);
  }

  return 1;
}

static int
write_scaling_file(const char * fname,
                   const char * input_fname)    // the file from which our distribution is built
{
  IWString_and_File_Descriptor output;

  if (! output.open(fname))
  {
    cerr << "Cannot open output scaling file '" << fname << "'\n";
    return 0;
  }

  return write_scaling_file(output, input_fname);
}

static int
normalise (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "O:i:s:m:vjN:p:C:U:uKo:qcb:zR:Ax:X:");

  verbose = cl.option_count('v');

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options\n";
    usage(4);
  }

  if (cl.option_present('m'))
  {
    cl.value('m', missing_value);

    if (verbose)
      cerr << "Missing value '" << missing_value << "'\n";
  }

  if (cl.option_present('j'))
  {
    header_records_to_skip = 1;
    initial_columns_to_skip = 1;
    is_descriptor_file = 1;

    if (verbose)
      cerr << "Will treat as a descriptor file\n";
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', header_records_to_skip) || header_records_to_skip < 0)
    {
      cerr << "Bad or missing value for header records to discard\n";
      usage(9);
    }
    if (verbose)
      cerr << header_records_to_skip << " records skipped at top of file\n";
  }

  if (cl.option_present('i'))
  {
    if (! cl.value('i', initial_columns_to_skip) || initial_columns_to_skip < 0)
    {
      cerr << "Bad or missing value for skip columns\n";
      usage(10);
    }

    if (verbose)
      cerr << initial_columns_to_skip << " columns will be ignored in each record\n";
  }

  if (cl.option_present('K'))
  {
    single_column_kludge = 1;
    if (verbose)
      cerr << "Single column kludge - will not properly match descriptor name\n";
  }

  if (cl.option_present('b'))
  {
    IWString b = cl.string_value('b');
    if (! char_name_to_char(b))
    {
      cerr << "Unrecognised input separator specification '" << b << "'\n";
      usage(1);
    }
    input_separator = b[0];
  }

  if (cl.option_present('z'))
  {
    suppress_zero_variance_columns = 0;

    if (verbose)
      cerr << "Will retain zero variance columns\n";
  }

  if (cl.option_present('R'))
  {
    if (! cl.value('R', nearest_integer) || nearest_integer <= 0.0 || nearest_integer >= 1.0)
    {
      cerr << "The round nearest integer option (-R) must be a number between 0 and 1\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will round output values to nearest integer if within " << nearest_integer << endl;
  }

  if (cl.option_present('A'))
  {
    categorical_variables_possibly_present = 1;

    if (verbose)
      cerr << "Categorial variables may be present\n";
  }

  if (cl.option_present('O'))
  {
    assert (1 == cl.option_count('O'));

    IWString tmp;
    cl.value('O', tmp);
    if ("log" == tmp)
    {
      operator_to_apply = TAKE_LOG;
      if (verbose)
        cerr << "Take the log of all values\n";
    }
    else
    {
      cerr << "Unrecognised operator '" << tmp << "'\n";
      usage(7);
    }
  }

  if (cl.option_present('u'))
  {
    if (cl.option_present('C'))
    {
      cerr << "The -u and -C options are mutually exclusive\n";
      usage(4);
    }

    if (! cl.option_present('U'))
    {
      cerr << "To perform an unscaling (-u) you must specify a Use file (-U)\n";
      usage(8);
    }

    perform_unscaling_operation = 1;
    if (verbose)
      cerr << "Will perform unscaling rather than scaling\n";
  }

  if (cl.option_present('o'))
  {
    const_IWSubstring o;
    for (int i = 0; cl.value('o', o, i); ++i)
    {
      if ("allow" == o)
      {
        allow_out_of_range_unscalings = 1;
        if (verbose)
          cerr << "Out of range scalings will be allowed\n";
      }
      else if ("trunc" == o)
      {
        truncate_low_values = 1;
        truncate_high_values = 1;
        if (verbose)
          cerr << "Out of range scalings will be truncated\n";
      }
      else if ("trunclow" == o)
      {
        truncate_low_values = 1;
        if (verbose)
          cerr << "Will truncate values below range\n";
      }
      else if ("trunchigh" == o)
      {
        truncate_high_values = 1;
        if (verbose)
          cerr << "Will truncate values above range\n";
      }
      else if ("quiet" == o)
      {
        report_out_of_range_values = 0;
        if (verbose)
          cerr << "Will not report out of range errors\n";
      }
      else
      {
        cerr << "Unrecoginsed -o directive '" << o << "'\n";
        usage(6);
      }
    }
  }

  if (cl.option_present('q'))
  {
    report_out_of_range_values = 0;
    if (verbose)
      cerr << "Will not report out of range errors\n";
  }

  if (cl.option_present('c'))
  {
    translate_descriptor_names_to_lowercase = 1;

    if (verbose)
      cerr << "Descriptor names will be translated to lowercase\n";
  }
  

  if (cl.option_present('N'))
  {
    if (1 != cl.option_count('N'))
    {
      cerr << "Only one -N option is allowed\n";
      usage(4);
    }

    const_IWSubstring n(cl.string_value('N'));

    if ("uv" == n || "UV" == n)
    {
      scaling_type = NRML_UNIT_VARIANCE;

      if (verbose)
        cerr << "Data scaled for unit variance\n";
    }
    else if ("01" == n)
    {
      scaling_type = NRML_MIN_TO_MAX;

      if (verbose)
        cerr << "Data scaled for [0,1] range\n";
    }
    else if ("11" == n)
    {
      scaling_type = NRML_SPREAD_ZERO;
      range_min = -1.0;
      range_max = 1.0;

      if (verbose)
        cerr << "Data scaled for [-1,1] range\n";
    }
    else if ("mc" == n)
    {
      scaling_type = NRML_MEAN_CENTER;
//    range_min = -1.0;
//    range_max = 1.0;

      if (verbose)
        cerr << "Data mean centered\n";
    }
    else if ("0255" == n)
    {
      scaling_type = NRML_0255;
    }
    else
    {
      cerr << "Unrecognised scaling type (-N option) '" << n << "'\n";
      usage(5);
    }
  }

  if (cl.option_count('X'))
  {
    if (! cl.value('X', user_specified_range_max_value) || user_specified_range_max_value <= 0.0)
    {
      cerr << "The range maximum value (-x) must be a non negative value\n";
      usage(2);
    }

    if (verbose)
      cerr << "Range max set to " << user_specified_range_max_value << endl;

    user_specified_range_max = true;
  }

  if (cl.option_present('x'))
  {
    if (! cl.value('x', constant_offset))
    {
      cerr << "Invalid constant offset specification (-x)\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will add " << constant_offset << " to each output value\n";
  }

  if (cl.option_present('p'))
  {
    const_IWSubstring p(cl.string_value('p'));

    if ('f' == p)
    {
      write_results_as_float = 1;

      if (verbose)
        cerr << "Will write results as float numbers\n";
    }
    else if (! p.numeric_value(output_precision) || output_precision < 1)
    {
      cerr << "The numeric precision must be a value whole positive number '" << p << "'\n";
      usage(5);
    }
    else if (verbose)
    {
      cerr << "Will show " << output_precision << " digits in output\n";
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments " << argc << "\n";
    usage(3);
  }

  if (cl.option_present('C') && cl.option_present('U'))
  {
    cerr << "Cannot specify both -C and -U options\n";
    usage(6);
  }

  if (cl.option_present('N') && cl.option_present('U'))
  {
    cerr << "The -N and -U options are mutually exclusive\n";
    usage(4);
  }

  if (cl.option_present('U'))
  {
    if (! cl.option_present('j'))
    {
      cerr << "Sorry, only works with descriptor files now (-j option) see Ian...\n";
      return 5;
    }

    cl.value('U', name_of_existing_scaling_file_to_use);

    if (verbose)
      cerr << "Existing distribution data in '" << name_of_existing_scaling_file_to_use << "'\n";
  }

  IWString_and_File_Descriptor output(1);

  if (! normalise(cl[0], output))
  {
    cerr << "Normalise failed\n";
    return 4;
  }

  if (cl.option_present('C'))
  {
    int data_found = 0;

    for (int i = 0; i < columns_in_input; i++)
    {
      if (column[i].n() > 0)
      {
        data_found = 1;
        break;
      }
    }

    if (! data_found)
    {
      cerr << "No data, cannot write the -C file\n";
      return 2;
    }

    const char * c = cl.option_value('C');

    if (! write_scaling_file(c, cl[0]))
    {
      cerr << "Cannot create scaling information file '" << c << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Scaling information written to '" << c << "'\n";
  }

  if (! report_out_of_range_values && values_out_of_range)
    cerr << values_out_of_range << " un-normalise values out of range\n";

  if (categorical_variables_possibly_present)
  {
    for (int i = 0; i < categorical_variables_found_in_column.number_elements(); ++i)
    {
      if (categorical_variables_found_in_column[i])
        cerr << categorical_variables_found_in_column[i] << " non numeric (assumed categorical) variables in column " << i << endl;
    }
  }

  delete [] column;
  delete [] d;

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = normalise(argc, argv);

  return rc;
}
