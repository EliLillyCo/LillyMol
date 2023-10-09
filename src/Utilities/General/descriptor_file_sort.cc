/*
  Sorts a descriptor file
*/

#include <iostream>
#include <memory>
#include <math.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/data_source/iwstring_data_source.h"

using std::endl;
using std::cerr;

const char * prog_name = NULL;

static int verbose = 0;

static int header_records_to_skip = 1;

static resizable_array<int> sort_column;

static resizable_array_p<IWString> sort_by_descriptor;

static int columns_in_input = 0;

static IWString missing_value;

static float missing_value_replacement_value = static_cast<float>(0.0);

static int discard_records_with_missing_data = 0;

static int greater_than = 1;

static int take_absolute_values = 0;

static int less_than = -1;

static char column_separator = ' ';

static int translate_non_numeric_to_number = 0;

static IW_STL_Hash_Map_float string_to_number;

/*
  We can sort by the difference between a pair of columns
*/

class Pair_of_Columns
{
  private:
    int _c1;
    IWString _s1;
    int _c2;
    IWString _s2;

  public:
    Pair_of_Columns();

    int build(const const_IWSubstring &);

    int debug_print(std::ostream &) const;

    int column_is(int col, const const_IWSubstring & token);

    int diff(const const_IWSubstring &, float &) const;

    int fully_specified() const;
};

static resizable_array_p<Pair_of_Columns> column_differences;

Pair_of_Columns::Pair_of_Columns()
{
  _c1 = -1;
  _c2 = -1;

  return;
}

int
Pair_of_Columns::debug_print(std::ostream & os) const
{
  os << "Pair_of_Columns:";
  if (_s1.length())
    os << _s1;
  if (_c1 >= 0)
    os << _c1;

  if (_s2.length())
    os << _s2;
  if (_c2 >= 0)
    os << _c2;

  os << "\n";

  return 1;
}

int
Pair_of_Columns::column_is(int col, const const_IWSubstring & token)
{
  if (_c1 >= 0)
    ;
  else if (_s1 == token)
  {
    _c1 = col;
    return 1;
  }

  if (_c2 >= 0)
    ;
  else if (_s2 == token)
  {
    _c2 = col;
    return 1;
  }

  return 0;
}

int
Pair_of_Columns::build(const const_IWSubstring & token)
{
  const_IWSubstring c1, c2;

  if (! token.split(_s1, '-', _s2) || 0 == _s1.length() || 0 == _s2.length())
  {
    cerr << "Pair_of_Columns::build: invalid specification '" << token << "'\n";
    return 0;
  }

  if (_s1 == _s2)
  {
    cerr << "Pair_of_Columns::build:both columns the same '" << token << "'\n";
    return 0;
  }

// If these are numeric already, we are done

  if (_s1.numeric_value(_c1) && _c1 > 0 &&
      _s2.numeric_value(_c2) && _c2 > 0 &&
      _c1 != _c2)
    return 1;

  return 1;
}

int
Pair_of_Columns::fully_specified() const
{
  if (_c1 < 0)
    return 0;

  if (_c1 < 0)
    return 0;

  return 1;
}

template <typename T>
int
parse_token_as_numeric(const const_IWSubstring & token,
                        T & v,
                        int col)
{
  if (token.numeric_value(v))
    return 1;

  if (token == missing_value)
  {
    v = missing_value_replacement_value;
    return 1;
  }

  cerr << "Illegal numeric '" << token << "' in column " << (col + 1) << endl;
  return 0;
}

template <typename T>
int
get_next_token(const T & buffer,
          const_IWSubstring & token,
          int & i)
{
  if (' ' == column_separator)
    return buffer.nextword(token, i);
  else
    return buffer.nextword_single_delimiter(token, i, column_separator);
}

int
Pair_of_Columns::diff(const const_IWSubstring & buffer,
                       float & d) const
{
  const_IWSubstring token;

  float v1, v2;
  int number_completed = 0;

  for (int i = 0, col = 0; get_next_token(buffer, token, i); col++)
  {
    if (col == _c1)
    {
      if (! parse_token_as_numeric(buffer, v1, col))
        return 0;

      number_completed++;
      if (2 == number_completed)
        break;
    }
    else if (col == _c2)
    {
      if (! parse_token_as_numeric(buffer, v2, col))
        return 0;

      number_completed++;
      if (2 == number_completed)
        break;
    }
  }

  if (2 != number_completed)
  {
    cerr << "Pair_of_Columns::diff:incomplete determination '" << buffer << "'\n";
    debug_print(cerr);
    return 0;
  }

  d = v1 - v2;

  return 1;
}

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "Sorts descriptor files\n";
  cerr << " -d <desc>      sort by the value of descriptor <desc>\n";
  cerr << " -c <col>       sort by the value in column <col>\n";
  cerr << " -s <number>    header records to skip (default 1)\n";
  cerr << " -M <string>    missing value string\n";
  cerr << " -m <float>     replace missing values with <float> (default " << missing_value_replacement_value << ")\n";
  cerr << " -z             discard records with missing data\n";
  cerr << " -i <char>      column separator (space by default)\n";
//cerr << " -t             input is tab separated\n";
  cerr << " -r             reverse sort\n";
  cerr << " -a             take asbolute values of input values\n";
  cerr << " -y             enable sorting by string values. Does not sort, but groups by common values\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

class Descriptor_File_Record : public IWString
{
  private:
    float * _sort_key;
    int _missing_data_count;
  public:
    Descriptor_File_Record (const const_IWSubstring &);
    ~Descriptor_File_Record ();

    int initialise_sort_keys (const int *, int, int, int & missing_values, float *);
    int initialise_sort_keys (const int *, int, const resizable_array_p<Pair_of_Columns>& , int & missing_values);

    const float * sort_keys () const { return _sort_key;}
};

Descriptor_File_Record::Descriptor_File_Record (const const_IWSubstring & d) : IWString(d)
{
  _missing_data_count = 0;
  _sort_key = NULL;
  return;
}

Descriptor_File_Record::~Descriptor_File_Record()
{
//if (NULL != _sort_key)
//  delete [] _sort_key;
}

static float
convert_to_numeric(const const_IWSubstring & token,
                   IW_STL_Hash_Map_float & string_to_number)
{
  if (0 == token.length())
    return -1.0f;

  const auto f = string_to_number.find(token);

  if (f != string_to_number.end())
    return f->second;

  int istop = 4;
  if (istop > token.length())
    istop = token.length();

  int tot = 0;

  for (int i = 0; i < istop; ++i)
  {
    tot = tot * 255 + token[i];
  }

  const int s = 1000000 * tot + string_to_number.size();
  string_to_number[token] = static_cast<float>(s);

  return static_cast<float>(s);
}

int
Descriptor_File_Record::initialise_sort_keys (const int * sort_column,
                                              int columns_in_input,
                                              int n,
                                              int & missing_values,
                                              float * s)
{
  missing_values = 0;

  int i = 0;
  const_IWSubstring token;

  _sort_key = s;

  for (int col = 0; col < columns_in_input && get_next_token(*this, token, i); col++)
  {
    int j = sort_column[col];

    if (j < 0)
      continue;

    token.strip_leading_blanks();
    token.strip_trailing_blanks();

    if (translate_non_numeric_to_number)
    {
      _sort_key[j] = convert_to_numeric(token, string_to_number);
      continue;
    }

    float f;
    if (token.numeric_value(f))
    {
      if (take_absolute_values)
        _sort_key[j] = fabs(f);
      else
        _sort_key[j] = f;
      continue;
    }

    if (missing_value == token)
    {
      missing_values++;
      _sort_key[j] = missing_value_replacement_value;
      continue;
    }

    cerr << "Illegal numeric '" << token << "' in column " << (col + 1) << endl;
    return 0;
  }

#ifdef DEBUG_INITIALISE_SORT_KEYS
  for (int i = 0; i < _sort_key.number_elements(); i++)
  {
    cerr << " i = " << i << " sort key " << _sort_key[i] << endl;
  }
#endif

  return 1;
}

#ifdef INITIALISE_SORT_KEYS
Never finished implemention of the difference idea

int
Descriptor_File_Record::initialise_sort_keys (const int * sort_column,
                                              int n,
                                              const resizable_array_p<Pair_of_Columns> & poc,
                                              int & missing_values)
{
  missing_values = 0;

  int npair = poc.number_elements();

  _sort_key.extend(n + npair, -1);

  int i = 0;
  const_IWSubstring token;

  for (int col = 0; get_next_token(*this, token, i); col++)
  {
    int j = sort_column[col];

    if (j < 0)
      continue;

    token.strip_leading_blanks();
    token.strip_trailing_blanks();

    float f;
    if (token.numeric_value(f))
    {
      if (take_absolute_values)
        _sort_key[j] = fabs(f);
      else
        _sort_key[j] = f;
      continue;
    }

    if (missing_value == token)
    {
      missing_values++;

      _sort_key[j] = missing_value_replacement_value;
      continue;
    }

    cerr << "Illegal numeric '" << token << "' in column " << (col + 1) << endl;
    return 0;
  }

  for (int i = 0; i < npair; i++)
  {
    if (! poc[i]->diff(*this, _sort_key[n + i]))
    {
      cerr << "Descriptor_File_Record::initialise_sort_keys:cannot initialise difference specification\n";
      return 0;
    }
  }

#ifdef DEBUG_INITIALISE_SORT_KEYS
  for (int i = 0; i < _sort_key.number_elements(); i++)
  {
    cerr << " i = " << i << " sort key " << _sort_key[i] << endl;
  }
#endif

  return 1;
}
#endif

static int
identify_difference_specifications (const const_IWSubstring & header,
                                    int columns_in_input,
                                    resizable_array<int> & sort_column,
                                    const resizable_array_p<Pair_of_Columns> & poc)
{
  int n = poc.number_elements();

  int need_to_check_column_names = 0;

  for (int i = 0; i < n; i++)
  {
    if (! poc[i]->fully_specified())
    {
      need_to_check_column_names = 1;
      break;
    }
  }

  if (! need_to_check_column_names)
    return 1;

// Tell each of the Pair_of_Columns objects which column is which

  int i = 0;
  const_IWSubstring token;
  get_next_token(header, token, i);

  for (int col = 0; get_next_token(header, token, i); col++)
  {
    for (int j = 0; j < n; j++)
    {
      poc[j]->column_is(col, token);
    }
  }

  for (int i = 0; i < n; i++)
  {
    if (poc[i]->fully_specified())
      continue;

    cerr << "Column difference specification incomplete\n";
    poc[i]->debug_print(cerr);

    return 0;
  }

  return 1;
}

static int
matches_except_for_quotes(const IWString & descriptor,
                          const const_IWSubstring & d)
{
  if (! d.starts_with('"'))
    return 0;

  if (! d.ends_with('"'))
    return 0;

  const int istop = d.length() - 1;
  for (int i = 1; i < istop; ++i)
  {
    if (descriptor[i-1] != d[i])
      return 0;
  }

  return 1;
}

static int
initialise_sort_conditions (const const_IWSubstring & header)
{
  if (' ' == column_separator)
    columns_in_input = header.nwords(column_separator);
  else
    columns_in_input = header.nwords_single_delimiter(column_separator);

  if (verbose)
    cerr << "Input contains " << columns_in_input << " columns\n";

  for (int i = 0; i < sort_column.number_elements(); i++)
  {
    if (sort_column[i] >= columns_in_input)
    {
      cerr << "Only " << columns_in_input << " columns in input, sorting on column " << (sort_column[i] + 1) << " is impossible\n";
      return 0;
    }
  }

  if (0 == column_differences.number_elements())
    ;
  else if (! identify_difference_specifications(header, columns_in_input, sort_column, column_differences))
  {
    cerr << "Cannot identify column difference specifications '" << header << "'\n";
    return 0;
  }

  int n = sort_by_descriptor.number_elements();
  if (0 == n)
    return 1;

  int i = 0;
  const_IWSubstring d;
  int col = 0;
  int number_matched = 0;

// We need to preserve the order on the command line, so extend the istart array

  int istart = sort_column.number_elements();

  sort_column.extend(n, -1);

  while (get_next_token(header, d, i))
  {
    for (int j = 0; j < n; j++)
    {
      if (sort_column[j + istart] >= 0)
        continue;

      if (*(sort_by_descriptor[j]) == d)
        ;
      else if (matches_except_for_quotes(*sort_by_descriptor[j], d))
        ;
      else
        continue;

      sort_column[j + istart] = col;
      number_matched++;

      if (verbose)
        cerr << "Descriptor '" << d << "' found in column " << (col + 1) << endl;

      if (number_matched == n)
        break;
    }

    if (number_matched == n)
      break;

    col++;
  }

  if (number_matched == n)
    return 1;

  cerr << "Only found " << number_matched << " matches for " << n << " descriptors\n";

  for (int j = 0; j < n; j++)
  {
    if (sort_column[istart + j] < 0)
      cerr << "No match for descriptor '" << *(sort_by_descriptor[j]) << "'\n";
  }

  return 0;
}

static int
read_data (iwstring_data_source & input,
           resizable_array_p<IWString> & header,
           resizable_array_p<Descriptor_File_Record> & zdata)
{
  int reading_header = (0 == header.number_elements());

  const_IWSubstring buffer;

  for (int i = 0; i < header_records_to_skip; i++)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Cannot read header record\n";
      return 0;
    }

    if (! reading_header)
      continue;

    if (0 == header.number_elements())
    {
      if (! initialise_sort_conditions(buffer))
      {
        cerr << "Cannot initialise sort conditions based on header\n";
        return 0;
      }
    }

    IWString * tmp = new IWString(buffer);
    header.add(tmp);
  }

  while (input.next_record(buffer))
  {
    Descriptor_File_Record * d = new Descriptor_File_Record(buffer);
    zdata.add(d);
  }

  return 1;
}

static int
read_data (const char * fname,
           resizable_array_p<IWString> & header,
           resizable_array_p<Descriptor_File_Record> & zdata)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "read_data:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_data(input, header, zdata);
}

class Descriptor_File_Record_Comparitor
{
  private:
    const int _ncols;

  public:
    Descriptor_File_Record_Comparitor(const int c) : _ncols(c) {};

    int operator () (const Descriptor_File_Record *, const Descriptor_File_Record *);
};

int
Descriptor_File_Record_Comparitor::operator () (const Descriptor_File_Record * d1,
                                                const Descriptor_File_Record * d2)
{
  const float * f1 = d1->sort_keys();
  const float * f2 = d2->sort_keys();

  for (int i = 0; i < _ncols; i++)
  {
    const float v1 = f1[i];
    const float v2 = f2[i];

    if (v1 < v2)
      return less_than;

    if (v1 > v2)
      return greater_than;
  }

  return 0;
}

static int
descriptor_file_sort (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vs:c:d:D:rM:m:zrai:qty");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('s'))
  {
    if (! cl.value('s', header_records_to_skip) || header_records_to_skip < 1)
    {
      cerr << "Invalid header records to skip (-c option)\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will skip " << header_records_to_skip << " header records\n";
  }

  if (cl.option_present('i'))
  {
    IWString i = cl.string_value('i');
    if (! char_name_to_char(i))
    {
      cerr << "The -i option must be followed by a single character, '" << i << "' is invalid\n";
      usage(4);
    }

    column_separator = i[0];

    if (verbose)
      cerr << "Column separator set to '" << column_separator << "'\n";
  }
  else if (cl.option_present('t'))
  {
    column_separator = '\t';
    if (verbose)
      cerr << "Input tab separated\n";
  }

  if (cl.option_present('y'))
  {
    translate_non_numeric_to_number = 1;
    if (verbose)
      cerr << "Will sort string values by converting to numeric forms\n";
  }

  if (cl.option_present('M'))
  {
    cl.value('M', missing_value);

    if (verbose)
      cerr << "Missing values represented as '" << missing_value << "'\n";
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', missing_value_replacement_value))
    {
      cerr << "Invalid value for missing value replacement (-m option)\n";
      return 5;
    }

    if (verbose)
      cerr << "Missing values replaced by " << missing_value_replacement_value << endl;
  }

  if (cl.option_present('z'))
  {
    if (cl.option_present('m'))
    {
      cerr << "The -m and -z options are mutually incompatible\n";
      usage(1);
    }

    discard_records_with_missing_data = 1;

    if (verbose)
      cerr << "Will discard records containing missing values\n";
  }

  if (! cl.option_present('c') && ! cl.option_present('d') && ! cl.option_present('D'))
  {
    cerr << "Must specify descriptors (-d), columns (-c) or differences (-D) on which to sort\n";
    usage(4);
  }

  if (cl.option_present('c'))
  {
    int i = 0;
    const_IWSubstring c;
    while (cl.value('c', c, i++))
    {
      int col;
      if (! c.numeric_value(col) || col < 1)
      {
        cerr << "Invalid sort column '" << c << "'\n";
        return 4;
      }

      sort_column.add(col - 1);

      if (verbose)
        cerr << "Will sort on column " << col << "\n";
    }
  }

  if (cl.option_present('d'))
  {
    int i = 0;
    const_IWSubstring d;
    while (cl.value('d', d, i++))
    {
      IWString * tmp = new IWString(d);

      sort_by_descriptor.add(tmp);

      if (verbose)
        cerr << "Will sort by descriptor '" << d << "'\n";
    }
  }

  if (cl.option_present('D'))
  {
    int i = 0;
    const_IWSubstring d;
    while (cl.value('D', d, i++))
    {
      Pair_of_Columns * p = new Pair_of_Columns;

      if (! p->build(d))
      {
        cerr << "invalid descriptor difference specification '" << d << "'\n";
        return i + 1;
      }
      column_differences.add(p);
    }

    if (verbose)
      cerr << "Will sort by " << column_differences.number_elements() << " column difference specifications\n";
  }

  if (cl.option_present('r'))
  {
    greater_than = -1;
    less_than = 1;

    if (verbose)
      cerr << "Reverse sort\n";
  }

  if (cl.option_present('a'))
  {
    take_absolute_values = 1;

    if (verbose)
      cerr << "Will take absolute values\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  resizable_array_p<IWString> header;
  resizable_array_p<Descriptor_File_Record> zdata;

  zdata.resize(20000);

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (verbose > 1)
      cerr << "Reading data from '" << cl[i] << "'\n";

    if (! read_data(cl[i], header, zdata))
    {
      cerr << "Cannot read data from '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  int n = zdata.number_elements();

  if (verbose)
    cerr << "Read " << header.number_elements() << " header records and " << n << " data records\n";

  int * tmp = new_int(columns_in_input, -1); std::unique_ptr<int[]> free_tmp(tmp);

  int nsort = sort_column.number_elements();

  for (int i = 0; i < nsort; i++)
  {
    int c = sort_column[i];
    tmp[c] = i;
  }

  float * sort_keys = new float[nsort * zdata.number_elements()]; std::unique_ptr<float[]> free_sort_keys(sort_keys);

  int records_removed_for_missing_data = 0;

  for (int i = n - 1; i >= 0; --i)
  {
    int missing_values = 0;

    if (! zdata[i]->initialise_sort_keys(tmp, columns_in_input, nsort, missing_values, sort_keys + nsort * i))
    {
      cerr << "Cannot initialise record " << i << endl;
      cerr << *(zdata[i]) << endl;
      return 5;
    }

    if (missing_values && discard_records_with_missing_data)
    {
      zdata.remove_item(i);
      records_removed_for_missing_data++;
    }
  }

  if (verbose && records_removed_for_missing_data)
    cerr << "Discarded " << records_removed_for_missing_data << " records for missing data\n";

  Descriptor_File_Record_Comparitor dfrc(nsort);

  zdata.iwqsort(dfrc);

  IWString_and_File_Descriptor output(1);

  for (int i = 0; i < header.number_elements(); i++)
  {
    output << *(header[i]) << '\n';
  }

  output.write_if_buffer_holds_more_than(32768);

  for (int i = 0; i < zdata.number_elements(); i++)
  {
    output << *(zdata[i]) << '\n';
    output.write_if_buffer_holds_more_than(32768);
  }

  if (output.length())
    output.flush();

  if (verbose)
    cerr << "Sorted " << zdata.number_elements() << " records\n";

  if (cl.option_present('q'))
    _exit(0);   // avoid delays while large arrays get deleted

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = descriptor_file_sort(argc, argv);

  return rc;
}
