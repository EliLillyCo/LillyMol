/*
  Finds the average of columns of numbers in a file
*/

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <random>
#include <algorithm>

#include "cmdline.h"
#include "accumulator.h"
#include "set_or_unset.h"
#include "iwstring_data_source.h"

#include "iwtokeniser.h"

#ifdef __INTEL_COMPILERdoesnotwork
#include "istrconv.h"
#endif

const char * prog_name = NULL;

static int verbose = 0;

static IW_Regular_Expression column_rx;

static IW_Regular_Expression descriptor_rx;

static char word_delimeter = ' ';

static int write_as_xydy1dy2 = 0;

static int brief_output = 0;

static int write_file_name_as_first_token = 0;
/*
  What do we do with a record where there are no matching columns
*/

static int ok_no_matching_columns = 0;

static Set_or_Unset<float> special_number;

static int write_column_number = 1;

static int take_absolute_value = 0;

/*
  If we are computing medians and percentiles, we need to store the values. We can
  impose a limit on the number of values stored
*/

static int value_buffer_size = std::numeric_limits<int>::max();

static std::random_device rd;

typedef float percentile_t;

static resizable_array<percentile_t> percentiles;

static int quoted_tokens = 0;

/*
  We may have two columns, where one is the raw values, and the other is the number
  of times that value occurs.
*/

static int prevalence_column = -1;
static IWString prevalence_column_name;
static double prevalence_multiplier = 1.0;

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
  cerr << "Computes statistics on columns in a file\n";
  cerr << " -c <col>       only consider column(s) <col>\n";
  cerr << " -R <rx>        column is the one following a match to <rx>\n";
  cerr << " -d <rx>        descriptor(s) to process\n";
  cerr << " -g             interpret the descriptor names as regular expressions\n";
  cerr << " -p <number>    report the number of values equal to <number>\n";
  cerr << " -M <string>    missing value string\n";
  cerr << " -k             ignore invalid numeric fields\n";
  cerr << " -s <number>    skip the first <number> records in the file\n";
  cerr << " -i <char>      word delimiter, default <space>\n";
  cerr << " -t             input is tab delimited\n";
  cerr << " -w             ignore records where no columns match the specifications\n";
  cerr << " -m             write as 'ave dy1 dy2' - works with xmgrace\n";
  cerr << " -j             treat as a descriptor file\n";
  cerr << " -b             brief output\n";
  cerr << " -f             write the file name as the first token of output\n";
  cerr << " -u             suppress writing the column number - cleaner output\n";
  cerr << " -k             ignore otherwise bad data\n";
  cerr << " -a             take absolute value of all input values\n";
  cerr << " -e <pct>       percentile values to compute\n";
  cerr << " -r <number>    if using percentiles, only keep <number> values\n";
  cerr << " -h <col>       column <col> (number or id) contains prevalence values for values in -c col\n";
  cerr << " -h mult=<x>    prevalence values should be integer. If not, specify multiplier to convert to int\n";
  cerr << " -q             header record may contain quoted variable names\n";
  cerr << " -y <sting>     used delimiter in the file ie. tab, space\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static extending_resizable_array<int> columns_to_process;

/*
  Once we have done the highest active column, we are done
*/

static int highest_column_number_to_check;

static int lines_read = 0;

static int ignore_bad_data = 0;

static int invalid_data_records_skipped = 0;

static IWString missing_value(".");

static int skip_first_records = 0;

static int is_descriptor_file = 1;

class AColumn : public Accumulator<double>
{
  private:
    IWString _name;

    int _number_missing_values;
    int _non_numeric_data;
    int _number_instances_special_value;

    resizable_array<float>_raw_values;

    std::mt19937_64 _rng;
    std::uniform_real_distribution<double> _uniform_01;
    std::uniform_int_distribution<int> _uniform_nraw;

//  private functions

    void _add_to_raw_values (float f);

  public:
    AColumn();
    ~AColumn();

    int extra(const_IWSubstring &, const int prevalence);

    int report(std::ostream &) const;

    void set_name(const const_IWSubstring & s) { _name = s;}

    const IWString & name () const { return _name;}

    int non_numeric_data () const { return _non_numeric_data;}

    void update_global_accumulator (Accumulator<double> &) const;
};

AColumn::AColumn() : _rng(rd()), _uniform_01(0.0, 1.0), _uniform_nraw(0, value_buffer_size)
{
  _number_missing_values = 0;
  _non_numeric_data = 0;
  _number_instances_special_value = 0;

  if (value_buffer_size > 0)
    _raw_values.resize(value_buffer_size);

  return;
}

AColumn::~AColumn ()
{
}

int
AColumn::extra(const_IWSubstring & token,
               const int prevalence)
{
  if (' ' != word_delimeter)
    token.strip_leading_blanks();

  float f;
#ifdef __INTEL_COMPILERdoesnotwork
  char * endptr;
  cerr << "Examining '" << token << "'\n";
  f = __IML_str_to_f(token.rawchars(), token.length(), 0, &endptr);
  cerr << "Value returned " << f << ", errno " << errno << " endptr '" << *endptr << "', diff " << (endptr - token.rawchars()) << endl;
  if ((EINVAL == errno && isspace(*endptr)) || 0 == errno)
#else
  if (token.numeric_value(f))
#endif
  {
    if (special_number.matches(f))
      _number_instances_special_value++;

    if (take_absolute_value)
      f = fabs(f);

    if (percentiles.number_elements())
      _add_to_raw_values(f);

    if (prevalence_column >= 0)
      return Accumulator<double>::extra(f, prevalence);
    else
      return Accumulator<double>::extra(f);
  }

  if (missing_value == token || 0 == token.length())
  {
    _number_missing_values++;
    return 1;
  }

  if (ignore_bad_data)
  {
    _non_numeric_data++;
    invalid_data_records_skipped++;
    return 1;
  }

  cerr << "AColumn::extra:invalid data '" << token << "'\n";
  return 0;
}

void
AColumn::_add_to_raw_values (float f)
{
  if (_raw_values.number_elements() < value_buffer_size)
  {
    _raw_values.add(f);
    return;
  }

  if (_uniform_01(_rng) > 0.01)    // don't replace an existing value
    return;

  int r = _uniform_nraw(_rng);

  _raw_values[r] = f;

  return;
}

int
AColumn::report(std::ostream & output) const
{
  if (_name.length() > 0)
    output << _name << ' ';

  const int n = Accumulator<double>::n();

  if (0 == n)
  {
    output << "AColumn::report: no data\n";
    return output.good();
  }

  const float ave = static_cast<float>(Accumulator<double>::average());
  const float std = static_cast<float>(sqrt(Accumulator<double>::variance()));

  output << n << " values";

  if (1 == n)
  {
    if (brief_output)
      output << " between " << Accumulator<double>::minval() << " and " << Accumulator<double>::maxval() << ", ave " << ave;
    else
      output << " between " << Accumulator<double>::minval() << " and " << Accumulator<double>::maxval() << ", tot " << Accumulator<double>::sum() << " ave " << ave << " std " << 0.0;
  }
  else if (write_as_xydy1dy2)
  {
    output << " ave " << ave << " dy1 " << (Accumulator<double>::maxval() - ave) << " dy2 " << (ave - Accumulator<double>::minval());
  }
  else if (brief_output)
    output << " between " << Accumulator<double>::minval() << " and " << Accumulator<double>::maxval() << " ave " << ave;
  else
    output << " between " << Accumulator<double>::minval() << " and " << Accumulator<double>::maxval() << ", tot " << Accumulator<double>::sum() << " ave " << ave << " std " << std;

  if (special_number.is_set() && n > 0)
  {
    float s;
    (void) special_number.value(s);

    float fraction = static_cast<float>(_number_instances_special_value) / static_cast<float>(n);

    output << ' ' << _number_instances_special_value << " instances of " << s << " fraction " << fraction;

    if (0.0f == s)
      output << " nz " << (1.0f - fraction);
  }

  if (_raw_values.number_elements() > 1)
  {
    float * r = const_cast<float *>(_raw_values.rawdata());

    const int nraw = _raw_values.number_elements();

    std::random_shuffle(r, r + nraw);
    std::sort(r, r + nraw);
    for (auto i = 0; i < percentiles.number_elements(); ++i)
    {
      const auto p = percentiles[i];
      output << ' ' << p << "% " << r[static_cast<int>(p * nraw / 100 + 0.4999)];
    }
  }

  output << '\n';

  return output.good();
}

void
AColumn::update_global_accumulator (Accumulator<double> & acc) const
{
  acc.extra(*this);

  return;
}

static AColumn * acolumn = NULL;

static int 
get_next_token (const const_IWSubstring & buffer,
                const_IWSubstring & token,
                int & i)
{
  if (' ' == word_delimeter)
    return buffer.nextword(token, i);
  else
    return buffer.nextword_single_delimiter(token, i, word_delimeter);
}

#ifdef NOT_USEEEEED_ANY_MORE
static int
get_next_token_csv (const const_IWSubstring & buffer,
                    const_IWSubstring & token,
                    int & i)
{
  int in_quote = 0;
  const int nstart = i;

  for (;i < buffer.length(); ++i)
  {
    const char c = buffer[i];
    if ('"' == c)
      in_quote = ! in_quote;
    else if (in_quote)
      ;
    else if (',' == c)
    {
      buffer.from_to(nstart, i-1, token);
      ++i;
      return 1;
    }
  }

  return 0;
}
#endif

static int
fetch_prevalence_value (const const_IWSubstring & buffer,
                        const int prevalence_column,
                        int & prevalence)
{
  int i = 0;
  const_IWSubstring token;

  for (int col = 0; get_next_token(buffer, token, i); ++col)
  {
    if (col != prevalence_column)
      continue;

    if (1.0 == prevalence_multiplier)
    {
      if (! token.numeric_value(prevalence) || prevalence < 0)
      {
        cerr << "Invalid prevalence value '" << token << "'\n";
        return 0;
      }
    }
    else
    {
      double tmp;
      if (! token.numeric_value(tmp) || tmp < 0.0)
      {
        cerr << "Invalid preference value '" << token << "'\n";
        return 0;
      }

      prevalence = static_cast<int>(tmp * prevalence_multiplier + 0.4999);
    }

    return 1;
  }

  cerr << "Never found prevalence column " << prevalence_column << endl;

  return 0;
}

static int
average (const const_IWSubstring & buffer)
{
  int i = 0;
  const_IWSubstring token;

  int col = 0;

  int found_match_this_record = 0;

  int prevalence = 1.0;

  if (prevalence_column >= 0)
  {
    if (! fetch_prevalence_value(buffer, prevalence_column, prevalence))
      return 0;
  }

  while (get_next_token(buffer, token, i))
  {
    if (column_rx.active())
    {
      if (column_rx.matches(token))
      {
        if (get_next_token(buffer, token, i))
        {
          if (! acolumn[0].extra(token, prevalence))
            return 0;

          found_match_this_record = 1;
          break;
        }
        else
        {
          cerr << "No token following regular expression match '" << buffer << "'\n";
          return 0;
        }
      }
    }
    else if (columns_to_process[col])
    {
      if (! acolumn[col].extra(token, prevalence))
        return 0;

      found_match_this_record = 1;

      if (col >= highest_column_number_to_check)
        break;
    }

    col++;
  }

  if (! found_match_this_record)
  {
    cerr << "Did not find any matching columns!!\n";
    return ok_no_matching_columns;
  }

  return 1;
}

static int
determine_descriptors_to_process (const const_IWSubstring & buffer)
{
  const_IWSubstring token;
  int col = 0;

  int nw;
  if (' ' == word_delimeter)
    nw = buffer.nwords();
  else
    nw = buffer.nwords_single_delimiter(word_delimeter);

  if (nw < 2)
  {
    cerr << "Must be at least two columns in a descriptor file\n";
    return 0;
  }
  
  acolumn = new AColumn[nw];

  if (NULL == acolumn)
  {
    cerr << "Cannot allocate " << col << " column data\n";
    return 0;
  }

  int matches_found = 0;

  resizable_array_p<IWString> d;

  IWTokeniser iwtokeniser(buffer);
  iwtokeniser.set_sep(word_delimeter);
  if (quoted_tokens)
    iwtokeniser.set_quoted_tokens(1);

  while (iwtokeniser.next_token(token))
  {
//  cerr << "Examining header token '" << token << "'\n";
    if (descriptor_rx.matches(token))
    {
      columns_to_process[col] = 1;
      matches_found++;

      highest_column_number_to_check = col;

      acolumn[col].set_name(token);

      if (verbose)
        cerr << "Will process column " << (col + 1) << " descriptor '" << token << "'\n";
    }

    if (token == prevalence_column_name)
      prevalence_column = col;

    col++;
  }

  if (0 == matches_found)
  {
    cerr << "No descriptor names match '" << descriptor_rx.source() << "'\n";
    return 0;
  }

  return 1;
}

static int
assign_column_names_csv (extending_resizable_array<int> & columns_to_process,
                         const const_IWSubstring & buffer)
{
  bool in_quote = false;
  int col = 0;
  IWString token;

  const int n = buffer.length();

  for (int i = 0; i < n; ++i)
  {
    const char c = buffer[i];

    if ('"' == c)
    {
      in_quote = ! in_quote;
      continue;
    }

    if (in_quote)
    {
      token += c;
      continue;
    }

    if (',' == c)
    {
      if (columns_to_process[col])
      {
        acolumn[col].set_name(token);
      }
      token.resize_keep_storage(0);
      col++;
      continue;
    }

  }

  return 1;
}

static int
assign_column_names (extending_resizable_array<int> & columns_to_process,
                     const const_IWSubstring & buffer)
{
  if (',' == word_delimeter && buffer.contains('"'))
    return assign_column_names_csv(columns_to_process, buffer);

  int col = 0;
  const_IWSubstring token;

  for (int i = 0; get_next_token(buffer, token, i) && col < columns_to_process.number_elements(); col++)
  {
    if (columns_to_process[col])
      acolumn[col].set_name(token);
  }

  return 1;
}

static int
average (const char * fname,
         iwstring_data_source & input,
         std::ostream & output)
{
  const_IWSubstring buffer;

  for (int i = 0; i < skip_first_records; i++)
  {
    input.next_record(buffer);

    if (0 > 0)
      ;
    else if (descriptor_rx.active())
    {
      if (! determine_descriptors_to_process(buffer))
      {
        cerr << "Invalid descriptor file header record\n";
        return 0;
      }
    }
    else if (is_descriptor_file)
      assign_column_names(columns_to_process, buffer);
  }

#ifdef SHOW_COLUMNS_TO_PROCESS
  for (int i = 0; i < columns_to_process.number_elements(); i++)
  {
    if (columns_to_process[i])
      cerr << "Will process column " << i << '\n';
  }
#endif

  int rc = 1;

  while (input.next_record(buffer))
  {
    lines_read++;

    if (! average(buffer))
    {
      cerr << "Fatal error on line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      if (lines_read > 2)
      {
        rc = 0;
        break;
      }
      else
        return 0;
    }
  }

  Accumulator<double> global_accumulator;

  int n = columns_to_process.number_elements();

  if (n > 0)
  {
    int columns_processed = 0;

    for (int i = 0; i < n; i++)
    {
      if (0 == columns_to_process[i])
        continue;

      const AColumn & ci = acolumn[i];
  
      if (0 == ci.n())
        continue;
  
      if (write_file_name_as_first_token)
        output << fname << ' ';

      if (! write_column_number)
        ;
      else if (0 == ci.name().length())
        output << "column " << (i + 1) << '\n';
      ci.report(output);

      ci.update_global_accumulator(global_accumulator);
      columns_processed++;
    }

    if (columns_processed > 1)
    {
      if (write_file_name_as_first_token)
        output << fname << ' ';
      output << "Global: " << global_accumulator.n() << " values between " << static_cast<float>(global_accumulator.minval()) << " and " << static_cast<float>(global_accumulator.maxval()) << " ave " << static_cast<float>(global_accumulator.average()) << "\n";
    }
  }
  else
    acolumn[0].report(output);

  return rc;
}

static int
average (const char * fname,
             std::ostream & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return average(fname, input, output);
}

static int
build_regular_expression_from_components (Command_Line & cl,
                                          char flag,
                                          IW_Regular_Expression & rx)
{
  IWString tmp;

  int i = 0;
  const_IWSubstring d;
  while (cl.value(flag, d, i++))
  {
    if (0 == tmp.length())
      tmp << "^(" << d;
    else
      tmp << "|" << d;
  }

  tmp << ")$";

  return rx.set_pattern(tmp);
}

static int
get_range_of_columns (const const_IWSubstring & cs,
                      extending_resizable_array<int> & col)
{
  const_IWSubstring r1, r2;
  if (! cs.split(r1, '-', r2) || 0 == r1.length() || 0 == r2.length())
    return 0;
  
  int c1;
  if (! r1.numeric_value(c1) || c1 < 1)
  {
    cerr << "Invalid start of range '" << r1 << "'\n";
    return 0;
  }
  
  int c2;
  if (! r2.numeric_value(c2) || c2 < 1 || c2 < c1)
  {
    cerr << "Invalid end of range '" << r2 << "'\n";
    return 0;
  }

  c1--;
  c2--;

  for (int i = c1; i <= c2; i++)
  {
    col[i] = 1;
  }

  return 1;
}

static int
get_comma_separated_columns (const const_IWSubstring & cs,
                             extending_resizable_array<int> & col)

{
  const_IWSubstring token;
  int i = 0;

  while (cs.nextword(token, i, ','))
  {
//  cerr << "Token is '" << token << "'\n";
    if (token.contains('-'))
    {
      if (! get_range_of_columns(token, col))
      {
        cerr << "Invalid range '" << token << "'\n";
        return 0;
      }
    }
    else
    {
      int c;
      if (! token.numeric_value(c) ||  c < 1)
      {
        cerr << "Invalid column specification '" << token << "'\n";
        return 0;
      }

      col[c-1] = 1;
    }
  }

  return 1;
}

static int
fetch_columns (Command_Line & cl,
               char flag,
               extending_resizable_array<int> & col)
{
  int i = 0;
  const_IWSubstring cs;
  while (cl.value(flag, cs, i++))
  {
    if (cs.contains('-'))
    {
      if (! get_range_of_columns(cs, col))
      {
        cerr << "Invalid range specification '" << cs << "'\n";
        return 0;
      }

      continue;
    }

    if (cs.contains(','))
    {
      if (! get_comma_separated_columns(cs, col))
      {
        cerr << "Invalid column specification(s) '" << cs << "'\n";
        return 0;
      }

      continue;
    }

    int c;
    if (! cs.numeric_value(c) || c < 1)
    {
       cerr << "Invalid column '" << cs << "'\n";
       return 0;
    }

    c--;

    col[c] = 1;
  }

  highest_column_number_to_check = col.number_elements();

  return col.number_elements();
}

static int
interpret_as_int_or_float (const const_IWSubstring p,
                           resizable_array<percentile_t> & percentiles)
{
  int i;
  if (p.numeric_value(i) && i > 0 && i < 100)
  {
    percentiles.add(static_cast<percentile_t>(i));
    return 1;
  }

  percentile_t f;
  if (! p.numeric_value(f) || f < static_cast<percentile_t>(0) || f > static_cast<percentile_t>(100))
  {
    cerr << "Invalid percentile '" << p << "'\n";
    return 0;
  }

  percentiles.add(f);

  return 1;
}

static int
get_percentiles (resizable_array<percentile_t> & percentiles,
                 const const_IWSubstring & p)
{
  if (! p.contains(','))
    return interpret_as_int_or_float(p, percentiles);

  const_IWSubstring token;
  int i = 0;
  while (p.nextword_single_delimiter(token, i, ','))
  {
    if (0 == token.length())
      continue;

    if (! interpret_as_int_or_float(token, percentiles))
      return 0;
  }

  return 1;
}

static int
average (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vc:R:M:s:kd:y:wp:mjbutafr:e:gh:i:q");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  int specifications = 0;

  if (cl.option_present('c'))
    specifications++;
    
  if (cl.option_present('R'))
    specifications++;
    
  if (cl.option_present('d'))
    specifications++;

  if (0 == specifications)
  {
    columns_to_process[0] = 1;
    acolumn = new AColumn[2];
  }
  else if (1 != specifications)
  {
    cerr << "Must have just one of -c, -d or -R options\n";
    usage(5);
  }

  if (cl.option_present('r'))    // must do before columns are allocated
  {
    if (! cl.value('r', value_buffer_size) || value_buffer_size < 10)
    {
      cerr << "The raw value buffer size (-r) option must be a sensible number of values to store\n";
      usage(2);
    }

    if (verbose)
      cerr << "Median computed based on first " << value_buffer_size << " values in the file\n";
  }

  if (cl.option_present('e'))
  {
    const_IWSubstring t;
    for (auto i = 0; cl.value('e', t, i); ++i)
    {
      if (! get_percentiles(percentiles, t))
      {
        cerr << "Invalid percentile specification '" << t << "'\n";
        return 2;
      }
    }
  }

  if (cl.option_present('q'))
  {
    quoted_tokens = 1;

    if (verbose)
      cerr << "Column names may have quoted tokens\n";
  }

  if (cl.option_present('c'))
  {
    if (! fetch_columns(cl, 'c', columns_to_process))
      return 4;

    acolumn = new AColumn[columns_to_process.number_elements() + 1];
  }
  else if (cl.option_present('R'))
  {
    const_IWSubstring r = cl.string_value('R');

    if (! column_rx.set_pattern(r))
    {
      cerr << "Invalid column regexp '" << r << "'\n";
      return 6;
    }

    if (verbose)
      cerr << "Will check columns that match '" << column_rx.source() << "'\n";

    acolumn = new AColumn[1];
  }
  else if (cl.option_present('d'))
  {
    if (cl.option_count('d') > 1)
    {
      if (! build_regular_expression_from_components(cl, 'd', descriptor_rx))
      {
        cerr << "Invalid descriptor regular expression specifications\n";
        usage(5);
      }
    }
    else
    {
      const_IWSubstring d(cl.string_value ('d'));

      int tmp;

      if (cl.option_present('g'))
        tmp = descriptor_rx.set_pattern(d);
      else
      {
        IWString drx;
        drx << '^' << d << '$';
        tmp = descriptor_rx.set_pattern(drx);
      }

      if (! tmp)
      {
        cerr << "Invalid descriptor regexp '" << d << "'\n";
        return 2;
      }

    }
    if (verbose)
      cerr << "Will check descriptors that match '" << descriptor_rx.source() << "'\n";

    skip_first_records = 1;
  }
  else if (0 == specifications)
    ;
  else
  {
    cerr << "Must specify either the -c or -R options\n";
    usage(6);
  }

  if (cl.option_present('h'))
  {
    const_IWSubstring h;
    for (int i = 0; cl.value('h', h, i); ++i)
    {
      if (h.starts_with("mult="))
      {
        h.remove_leading_chars(5);
        if (! h.numeric_value(prevalence_multiplier) || prevalence_multiplier < 1.0)
        {
          cerr << "The prevalence column multiplier must be a number > 1.0, '" << h << "' invalid\n";
          return 2;
        }

        if (verbose)
          cerr << "Prevalence values multiplied by " << prevalence_multiplier << " to conver to counts\n";
      }
      else if (h.numeric_value(prevalence_column) && prevalence_column > 0)
        prevalence_column = prevalence_column - 1;
      else
        prevalence_column_name = h;
    }

  }

  if (cl.option_present('m'))
  {
    write_as_xydy1dy2 = 1;
    if (verbose)
      cerr << "Will write as x ave dy1 dy2 - for xmgrace\n";
  }

  if (cl.option_present('b'))
  {
    brief_output = 1;

    if (verbose)
      cerr << "Brief output only\n";
  }

  if (cl.option_present('y'))
  {
    IWString tmp;
    cl.value('y', tmp);

    if (! char_name_to_char(tmp))
    {
      cerr << "The word delimiter must be a single character, '" << tmp << "' is invalid\n";
      usage(4);
    }

    word_delimeter = tmp[0];

    if (verbose)
      cerr << "Word delimiter '" << word_delimeter << "'\n";
  }
  else if (cl.option_present('i'))
  {
    IWString tmp;
    cl.value('i', tmp);

    if (! char_name_to_char(tmp))
    {
      cerr << "The word delimiter must be a single character, '" << tmp << "' is invalid\n";
      usage(4);
    }

    word_delimeter = tmp[0];

    if (verbose)
      cerr << "Word delimiter '" << word_delimeter << "'\n";
  }
  else if (cl.option_present('t'))
  {
    word_delimeter = '\t';
    if (verbose)
      cerr << "Input assumed tab delimited\n";
  }

  if (cl.option_present('a'))
  {
    take_absolute_value = 1;

    if (verbose)
      cerr << "Will take average values\n";
  }

  if (cl.option_present('w'))
  {
    ok_no_matching_columns = 1;

    if (verbose)
      cerr << "Will ignore records where nothing matches the column specification(s)\n";
  }

  if (cl.option_present('p'))
  {
    float p;
    if (! cl.value('p', p))
    {
      cerr << "Invalid special number to monitor (-p option)\n";
      usage(4);
    }

    special_number.set(p);

    if (verbose)
      cerr << "Will separately report the number of instances of " << p << '\n';
  }

  if (cl.option_present('u'))
  {
    write_column_number = 0;

    if (verbose)
      cerr << "Will not write column numbers\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('M'))
  {
    missing_value = cl.string_value('M');

    if (verbose)
      cerr << "Missing value string '" << missing_value << "'\n";
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', skip_first_records) || skip_first_records < 0)
    {
      cerr << "Invalid value of the -s option\n";
      usage(8);
    }

    if (verbose)
      cerr << "Will skip the first " << skip_first_records << " records in the file\n";
  }

  if (cl.option_present('j'))
  {
    is_descriptor_file = 1;
    skip_first_records = 1;

    if (verbose)
      cerr << "Will treat as a descriptor file\n";
  }

  if (cl.option_present('k'))
  {
    ignore_bad_data = 1;

    if (verbose)
      cerr << "Will ignore non-numeric input\n";
  }

  if (cl.option_present('f'))
  {
    write_file_name_as_first_token = 1;

    if (verbose)
      cerr << "First token of output is file name\n";
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, only processes one file at a time\n";
    usage(5);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (verbose)
      cerr << "Begin processing '" << cl[i] << "'\n";

    if (! average(cl[i], std::cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (ignore_bad_data && invalid_data_records_skipped)
    cerr << "Skipped " << invalid_data_records_skipped << " invalid data records\n";

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = average(argc, argv);

  return rc;
}
