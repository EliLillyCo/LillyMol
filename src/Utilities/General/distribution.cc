/*
  We want to get a distribution of a set of numbers
*/

#include <stdlib.h>
#include <limits>
#include <algorithm>

#define RESIZABLE_ARRAY_IMPLEMENTATION 1
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION 1
//#define IWQSORT_FO_IMPLEMENTATION 1

#include "cmdline.h"
#include "accumulator.h"
#include "iwqsort.h"
#include "misc.h"
#include "iwstring_data_source.h"
#include "iw_stl_hash_map.h"

const char * prog_name = NULL;

static int verbose = 0;

//static int normalise_distribution_to_highest = 0;

//static int normalise_distribution_to_unit_area = 0;

#define NORMALISE_DISTRIBUTION_TO_UNIT_AREA 1
#define NORMALISE_DISTRIBUTION_TO_HIGHEST 2

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Create distribution based on a larger input set\n";
  cerr << " -c <col>       column to process\n";
  cerr << " -d <dname>     descriptor to process\n";
  cerr << " -a             normalise output to unit area\n";
  cerr << " -h             normalise output to highest value\n";
  cerr << " -n <int>       round to the nearest int of type <int>\n";
  cerr << " -r <int>       round to this accuracy - 10 means accuracy of 1 in 10\n";
  cerr << " -D <fname>     read per-descriptor processing info from <fname>\n";
  cerr << "                dname <nhv> <nua> acc=<n> nint=<n>\n";
  cerr << " -s <n>         header records to skip\n";
  cerr << " -H <string>    write header record, column name <string>\n";
  cerr << " -S <fname>     create stats file for molecular_property_profile\n";
  cerr << " -R ...         directives for creating R programme fragments\n";
  cerr << " -R fname=<fname>  output file created\n";
  cerr << " -R cname=<fname>  name of the collection (likely PBCHM...)\n";
  cerr << " -i <char>      input delimiter (default ' ')\n";
  cerr << " -z             when rounding is in effect, insert zero values\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

class IWFloatHash
{
  private:
  public:

#if defined (__INTEL_COMPILER)

    static const size_t bucket_size = 4;
    static const size_t min_buckets = 8;
    bool  operator () (float &, float) const;

#endif

    size_t operator () (float) const;
};

size_t
IWFloatHash::operator () (float f) const
{
  void * v = &f;

  size_t * rc = (size_t *) v;

  return *rc;
}

// Basically a scatter operation

static void
transfer_columns_from_cmdline (const resizable_array<int> & columns_specified_on_command_line,
                               int * process_column)
{
  for (int i = 0; i < columns_specified_on_command_line.number_elements(); ++i)
  {
    process_column[columns_specified_on_command_line[i]] = 1;
  }

  return;
}

template <typename T>
int
transfer_rounding_from_cmdline (const int columns_in_input,
                                const int * process_column,
                                const resizable_array<T> & from_cmdline,
                                T * v)
{
  const int cbp = std::count(process_column, process_column + columns_in_input, 1);

  if (cbp == from_cmdline.number_elements())
  {
    int ndx = 0;
    for (int i = 0; i < columns_in_input; ++i)
    {
      if (0 == process_column[i])
        continue;

      v[i] = from_cmdline[ndx];
      ndx++;
    }
  }
  else if (1 == from_cmdline.number_elements())
  {
    for (int i = 0; i < columns_in_input; ++i)
    {
      if (process_column[i])
        v[i] = from_cmdline[0];
    }
  }
  else
  {
    cerr << "transfer_rounding_from_cmdline:processing " << cbp << " columns, but got " << from_cmdline.number_elements() << " rounding specifiers\n";
    return 0;
  }

  return 1;
}

template <typename T>
class Value_Count
{
  private:
    T _v;
    int _count;

  public:
    Value_Count (T v, int n) : _v(v) { _count = n;}

    T v () const { return _v;}

    int count () const { return _count;}
};

typedef Value_Count<float> VCF;

template class resizable_array_p<VCF>;
template class resizable_array_base<VCF *>;

class Value_Count_Comparator
{
  private:
  public:
    int operator () (const VCF * vc1, const VCF * vc2) const;
};

int
Value_Count_Comparator::operator() (const VCF * vc1, const VCF * vc2) const
{
  if (vc1->v() < vc2->v())
    return -1;

  if (vc1->v() > vc2->v())
    return 1;

  cerr << "Yipes, equal values found in comparison, should not happen\n";
  return 0;
}

class Rfile
{
  private:
    IWString _collection;      // PBCHM ....

    IWString_and_File_Descriptor _output;

  public:
    int build (Command_Line & cl, const char flag, const int verbose);

    int active() const { return _output.is_open();}

    int statistics(const IWString & dname, const Accumulator<double> &);
    int data(const IWString & dname, const resizable_array_p<VCF> & vcf);

    template <typename X, typename Y> int data (const IWString & dname, const resizable_array<X> & x, const resizable_array<Y> & y);
};

int
Rfile::build (Command_Line & cl,
              const char flag,
              const int verbose)
{
  IWString fname;

  const_IWSubstring r;
  for (int i = 0; cl.value(flag, r, i); ++i)
  {
    if (r.starts_with("fname="))
    {
      r.remove_leading_chars(6);
      fname = r;
    }
    else if (r.starts_with("cname="))
    {
      r.remove_leading_chars(6);
      _collection = r;
    }
    else if ("help" == r)
    {
    }
    else
    {
      cerr << "Rfile::build:unrecognised directive '" << r << "'\n";
      return 0;
    }
  }

  if (0 == fname.length())
  {
    cerr << "Rfile::build:must specify file name via the 'fname=' directive\n";
    return 0;
  }

  if (0 == _collection.length())
  {
    cerr << "Rfile::build:must specify collection name via the 'cname=' directive\n";
    return 0;
  }

  if (! _output.open(fname.c_str()))
  {
    cerr << "Rfile::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return 1;
}

int
Rfile::statistics (const IWString & dname, const Accumulator<double> & acc)
{
  _output << _collection << '.' << dname << ".min <- " << acc.minval() << '\n';
  _output << _collection << '.' << dname << ".max <- " << acc.maxval() << '\n';
  _output << _collection << '.' << dname << ".ave <- " << static_cast<float>(acc.average()) << '\n';

  return 1;
}

int
Rfile::data (const IWString & dname,
             const resizable_array_p<VCF> & vcf)
{
  const int n = vcf.number_elements();

  _output << _collection << '.' << dname << ".x = c(";

  for (int i = 0; i < n; ++i)
  {
    if (i > 0)
      _output << ',';
    _output << vcf[i]->v();
  }
  _output << ")\n";

  _output << _collection << '.' << dname << ".y = c(";

  for (int i = 0; i < n; ++i)
  {
    if (i > 0)
      _output << ',';
    _output << vcf[i]->count();
  }
  _output << ")\n";


  _output.write_if_buffer_holds_more_than(4096);

  return 1;
}

template <typename X, typename Y>
int
Rfile::data (const IWString & dname,
             const resizable_array<X> & x,
             const resizable_array<Y> & y)
{
  const int n = x.number_elements();
  assert(n == y.number_elements());

  _output << _collection << '.' << dname << ".x = c(";

  for (int i = 0; i < n; ++i)
  {
    if (i > 0)
      _output << ',';
    _output << x[i];
  }

  _output << ")\n";

  _output.write_if_buffer_holds_more_than(4096);

  _output << _collection << '.' << dname << ".y = c(";

  for (int i = 0; i < n; ++i)
  {
    if (i > 0)
      _output << ',';
    _output << y[i];
  }

  _output << ")\n";

  _output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static void
do_nearest_int_rounding (float & v, const int round_to_nearest_int)
{
  if (v >=  static_cast<float>(round_to_nearest_int) / static_cast<float>(2.0))
  {
    int tmp = static_cast<int>(v + 0.4999);
    int i = (tmp / round_to_nearest_int) * round_to_nearest_int;

    float d1 = i + round_to_nearest_int - v;
    float d2 = v - i;

    if (d1 < d2)
      v = static_cast<float>(i + round_to_nearest_int);
    else
      v = static_cast<float>(i);
  }
  else if (v <= static_cast<float>(- round_to_nearest_int) / static_cast<float>(2.0))
  {
    int tmp = static_cast<int>(v - 0.49999);

    int i = (tmp / round_to_nearest_int) * round_to_nearest_int;

    float d1 = i - round_to_nearest_int - v;
    float d2 = v - i;

    if (d1 < d2)
      v = static_cast<float>(i - round_to_nearest_int);
    else
      v = static_cast<float>(i);
  }
  else
  {
    v = 0.0;
    return;
  }

  return;
}

typedef IW_Hash_Map<float, int, IWFloatHash> Hash_Map_Float_Int;

/*
  When reading multiple descriptors, we can read the handling of each from a file
*/

class Descriptor_Processing
{
  private:
    bool _normalisation;

    int  _round_to_nearest_int;
    float  _round_to_significant_figures;

  public:
    Descriptor_Processing();

    int build (const const_IWSubstring & buffer);

    void initialise (float &, int &, int &) const;
};

Descriptor_Processing::Descriptor_Processing()
{
  _normalisation = 0;

  _round_to_nearest_int = 0;
  _round_to_significant_figures = 0.0f;

  return;
}

int
Descriptor_Processing::build (const const_IWSubstring & buffer)
{
  IWString mybuffer(buffer);
  mybuffer.to_lowercase();

  const_IWSubstring token;

  for (int i = 0; mybuffer.nextword(token, i); )
  {
    if ("nhv" == token)
    {
      _normalisation = NORMALISE_DISTRIBUTION_TO_UNIT_AREA;
    }
    else if ("nua" == token)
    {
      _normalisation = NORMALISE_DISTRIBUTION_TO_HIGHEST;
    }
    else if (token.starts_with("acc="))
    {
      token.remove_leading_chars(4);
      if (! token.numeric_value(_round_to_significant_figures) || _round_to_significant_figures < 1.0f)
      {
        cerr << "Descriptor_Processing::build:invalid significant figures directive '" << buffer << "'\n";
        return 0;
      }
    }
    else if (token.starts_with("nint="))
    {
      token.remove_leading_chars(5);
      if (! token.numeric_value(_round_to_nearest_int) || _round_to_nearest_int < 1)
      {
        cerr << "Descriptor_Processing::build:invalid nearest int directive '" << buffer << "'\n";
        return 0;
      }
    }
    else 
    {
      cerr << "Descriptor_Processing::build:unrecognised directive '" << token << "', ignored\n";
    }
  }

//cerr << "Descriptor_Processing::build:values " << _round_to_significant_figures << " and " << _round_to_nearest_int << " from " << buffer << endl;

  return 1;
}

void Descriptor_Processing::initialise (float & rounding, int & round_to_nearest, int & normalisation) const
{
//cerr << "Descriptor_Processing::initialise:input " << rounding << " me " << _round_to_significant_figures << " " << round_to_nearest << " me " << _round_to_nearest_int << endl;
  if (0.0f == rounding)
    rounding = _round_to_significant_figures;
  if (0 == round_to_nearest)
    round_to_nearest = _round_to_nearest_int;

  if (0 == normalisation)
    normalisation = _normalisation;

  return;
}

class Data_From_File
{
  private:
    int _columns_in_input;

    int _header_records_to_skip;

    resizable_array<int> _columns_specified_on_command_line;

    IW_STL_Hash_Map_int _descriptor;

    int * _process_column;

    Accumulator<double> * _acc;
    std::unordered_map<float, int> * _zdata;

    float * _rounding;

    int * _round_to_nearest_int;

    int * _y_normalisation;

    int _insert_zero_values_for_missing;

//  Associated data that helps set things up

    IWString * _descriptor_name;

//  We will get in a lot of descriptors looking like w_
//  We can remove that prefix

    IWString _common_descriptor_prefix;

    IWString _collection_name;

//  If someone wants different rounding in different columns, they either enter each factor
//  in order, or specify a column name  'dname=xx'

    resizable_array<float> _rounding_from_cmdline;    // must contain 1 value or _columns_in_input

    IW_STL_Hash_Map_float _by_descriptor_rounding_from_cmdline;

    resizable_array<int> _round_to_nearest_int_from_cmdline;   // must contain 1 value of _columns_in_input

    IW_STL_Hash_Map_float _by_descriptor_round_to_nearest_int_from_cmdline;

    IW_STL_Hash_Map<IWString, Descriptor_Processing *> _descriptor_processing;

    char _input_delimiter;

//  private functions

    int _initialise_arrays (const const_IWSubstring & buffer);
    int _determine_columns_to_process (const const_IWSubstring & buffer);
    int _extract_data (const const_IWSubstring & buffer);

    int _remove_common_descriptor_prefix();

    int _final_processing (const int col, Rfile & rfile, const IWString & stats_file_stem, const IWString & output_file_stem);
    int _write_sorted_values (const resizable_array_p<VCF> & vcf, const int col, const int n, const int highest, IWString_and_File_Descriptor & output) const;

    int _write_stats_file (const IWString & stats_file_stem, const int col) const;

    int _write_sorted_values_insert_zero_float (const resizable_array_p<VCF> & vcf,
                                                        const int col,
                                                        const int n,
                                                        const int highest,
                                                        IWString_and_File_Descriptor & output) const;
    int _write_sorted_values_insert_zero_int (const resizable_array_p<VCF> & vcf,
                                                        const int col,
                                                        const int n,
                                                        const int highest,
                                                        IWString_and_File_Descriptor & output) const;

    int _write_normalised_value (const int c, const int col, const int n,
                                         const int highest, IWString_and_File_Descriptor & output) const;

    int _transfer_data_to_rfile (Rfile & rfile, const int col, const int n, const int highest, const resizable_array_p<VCF> & vcf) const;
    int _transfer_data_to_rfile_float (Rfile & rfile, const int col, const int n, const int highest, const resizable_array_p<VCF> & vcf) const;
    int _transfer_data_to_rfile_int   (Rfile & rfile, const int col, const int n, const int highest, const resizable_array_p<VCF> & vcf) const;

    template <typename T> int _append_normalised_value (const int c, const int col, const int n, const int highest, resizable_array<T> & x) const;

  public:
    Data_From_File();
    ~Data_From_File();

    void set_header_records_to_skip (int s) { _header_records_to_skip = s;}
    void set_input_delimiter (char s) { _input_delimiter = s;}

    void set_normalisation(int s);

    void set_insert_zero_values_for_missing (int s) { _insert_zero_values_for_missing = s;}

    int set_descriptors_to_process (Command_Line & cl, char flag);
    int set_columns_to_process     (Command_Line & cl, char flag);

    int set_per_column_rounding         (Command_Line & cl, char flag);
    int set_per_column_integer_rounding (Command_Line & cl, char flag);

    void set_common_descriptor_prefix (const char * s) { _common_descriptor_prefix = s;}
    void set_collection_name (const char * s) { _collection_name = s;}

    int read_descriptor_processing (const char * fname);
    int read_descriptor_processing (iwstring_data_source &);

    int columns_being_processed () const { return std::count(_process_column, _process_column + _columns_in_input, 1);}

    int report (std::ostream &) const;

    int build (const char * fname);
    int build (iwstring_data_source &);

    int doit (const IWString & stats_file_stem, Rfile & rfile, const IWString & output_file_stem);
};

Data_From_File::Data_From_File ()
{
  _columns_in_input = -1;
  _header_records_to_skip = 0;

  _process_column = nullptr;
  _acc = nullptr;
  _zdata = nullptr;

  _descriptor_name = nullptr;

  _rounding = nullptr;

  _round_to_nearest_int = nullptr;

  _y_normalisation = nullptr;

  _insert_zero_values_for_missing = 0;

  _input_delimiter = ' ';

  return;
}

Data_From_File::~Data_From_File ()
{
  if (nullptr != _process_column)
  {
    delete [] _process_column;
    delete [] _acc;
    delete [] _zdata;
  }

  if (nullptr != _descriptor_name)
    delete [] _descriptor_name;

  if (nullptr != _rounding)
    delete [] _rounding;

  if (nullptr != _round_to_nearest_int)
    delete [] _round_to_nearest_int;

  if (nullptr != _y_normalisation)
    delete [] _y_normalisation;

  for (auto i : _descriptor_processing)
  {
    delete i.second;
  }

  return;
}

int
Data_From_File::read_descriptor_processing (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Data_From_File::read_descriptor_processing:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_descriptor_processing(input);
}

int
Data_From_File::read_descriptor_processing(iwstring_data_source & input)
{
  input.set_strip_trailing_blanks(1);
  input.set_strip_leading_blanks(1);
  input.set_dos(1);

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (0 == buffer.length() || buffer.starts_with('#'))
      continue;

    IWString dname;
    int i = 0;
    buffer.nextword(dname, i);

    buffer += i;

    if (_descriptor_processing.contains(dname))
    {
      cerr << "Data_From_File::read_descriptor_processing:duplicate '" << dname << "', ignored\n";
      continue;
    }

    Descriptor_Processing * d = new Descriptor_Processing();

    if (! d->build(buffer))
    {
      cerr << "Data_From_File::read_descriptor_processing:invalid directives for '" << dname << "' '" << buffer << "'\n";
      delete d;
      return 0;
    }

    _descriptor_processing[dname] = d;
  }

  return _descriptor_processing.size();
}

int
Data_From_File::set_descriptors_to_process (Command_Line & cl, char flag)
{
  IWString d;
  for (int i = 0; cl.value(flag, d, i); ++i)
  {
    _descriptor[d] = 1;

    if (verbose)
      cerr << "Will read values associated with descriptor '" << d << "'\n";
  }

  return _descriptor.size();
}

int
Data_From_File::set_columns_to_process (Command_Line & cl, char flag)
{
  const_IWSubstring c;

  for (int i = 0; cl.value(flag, c, i); ++i)
  {
    int col;
    if (! c.numeric_value(col) || col < 1)
    {
      cerr << "Data_From_File::set_columns_to_process:invalid column number '" << c << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Will process data in column " << col << endl;

    _columns_specified_on_command_line.add_if_not_already_present(col - 1);
  }

  return _columns_specified_on_command_line.number_elements();
}

template <typename T>
int
add_column_name_directive (const const_IWSubstring & token,
                           IW_STL_Hash_Map<IWString, T> & by_name)
{
  IWString dname, s;

  if (! token.split(dname, '=', s) || 0 == dname.length() || 0 == s.length())
  {
    cerr << "add_column_name_directive:cannot split into directive=value form\n";
    return 0;
  }

  T v;

  if (! s.numeric_value(v) || v < static_cast<T>(0))
  {
    cerr << "Invalid numeric '" << token << "'\n";
    return 0;
  }

  by_name[dname] = v;

  return 1;
}

template <typename T>
int
add_column_number_directive (const const_IWSubstring & token,
                             resizable_array<T> & by_col)
{
  T v;

  if (! token.numeric_value(v) || v <= static_cast<T>(0))
  {
    cerr << "add_column_number_directive:invalid numeric '" << token << "'\n";
    return 0;
  }

  by_col.add(v);

  return 1;
}

template <typename T>
int
fill_by_column_or_by_name_arrays (Command_Line & cl, 
                                  const char flag,
                                  resizable_array<T> & by_col,
                                  IW_STL_Hash_Map<IWString, T> & by_name)
{
  const_IWSubstring token;

  for (int i = 0; cl.value(flag, token, i); ++i)
  {
    if (token.contains('='))
    {
      if (! add_column_name_directive(token, by_name))
      {
        cerr << "fill_by_column_or_by_name_arrays:invalid column name directive '" << token << "'\n";
        return 0;
      }
    }
    else
    {
      if (! add_column_number_directive(token, by_col))
      {
        cerr << "fill_by_column_or_by_name_arrays::invalid numeric '" << token << "'\n";
        return 0;
      }
    }
  }

  if (0 == by_col.size() && 0 == by_name.size())
  {
    cerr << "fill_by_column_or_by_name_arrays::no data for option '" << flag << "'\n";
    return 0;
  }

  if (by_name.size() > 0 && by_col.size() > 0)
  {
    cerr << "fill_by_column_or_by_name_arrays::cannot mix named (" << by_name.size() << ") and positional (" << by_col.size() << ") numeric specifiers\n";
    return 0;
  }

  return 1;
}

int
Data_From_File::set_per_column_rounding (Command_Line & cl, char flag)
{
  if (! fill_by_column_or_by_name_arrays(cl, flag, _rounding_from_cmdline, _by_descriptor_round_to_nearest_int_from_cmdline))
  {
    return 0;
  }

  return 1;
}

int
Data_From_File::set_per_column_integer_rounding (Command_Line & cl, char flag)
{
  const_IWSubstring r;
  for (int i = 0; cl.value(flag, r, i); ++i)
  {
    float x;

    if (! r.numeric_value(x) || x <= 0.0f)
    {
      cerr << "Data_From_File::set_per_column_integer_rounding:rounding (" << flag << ") must be non negative numbers, '" << r << "' invalid\n";
      return 0;
    }
    _round_to_nearest_int_from_cmdline.add(x);
  }

  return 1;
}

void
Data_From_File::set_normalisation (int s)
{
  assert (_columns_in_input > 0);

  std::fill_n(_y_normalisation, _columns_in_input, s);

  return;
}

int
Data_From_File::report (std::ostream & output) const
{
  output << "Data_From_File::report:input file has " << _columns_in_input << " columns\n";

  for (int i = 0; i < _columns_in_input; ++i)
  {
    if (0 == _process_column[i])
      continue;

    output << i;
    if (nullptr != _descriptor_name && _descriptor_name[i].length())
      output << ' ' << _descriptor_name[i];

    output << ' ' << _acc[i].n() << " values btw " << _acc[i].minval() << " and " << _acc[i].maxval() << " ave " << static_cast<float>(_acc[i].average()) << '\n';
  }

  return 1;
}

int
Data_From_File::_extract_data (const const_IWSubstring & buffer)
{
  const_IWSubstring mybuffer(buffer);

  const_IWSubstring token;

  for (int col = 0, i = 0; buffer.nextword(token, i, _input_delimiter); ++col) 
  {
    if (0 == _process_column[col])
      continue;

    float v;
    if (! token.numeric_value(v))
    {
      cerr << "Invalid numeric '" << token << "'\n";
      return 0;
    }

    _acc[col].extra(v);

    if (_rounding[col] > 0.0)
    {
      int tmp;
      if (v < 0.0)
        tmp = static_cast<int>(v * _rounding[col] - 0.49999999);
      else
        tmp = static_cast<int>(v * _rounding[col] + 0.49999999);
      v = static_cast<float>(tmp) / _rounding[col];
    }
    else if (_round_to_nearest_int[col] > 0)
    {
//  #define DEBUG_NEAREST_INT_ROUNDING
#ifdef DEBUG_NEAREST_INT_ROUNDING
      cerr << "Converted " << v;
#endif
      do_nearest_int_rounding(v, _round_to_nearest_int[col]);
#ifdef DEBUG_NEAREST_INT_ROUNDING
      cerr << " to " << v << endl;
#endif
    }
  
    auto f = _zdata[col].find(v);
    if (f == _zdata[col].end())
    {
//    cerr << "col " << col << " first encounter " << v << endl;
//    _zdata[col][v] = 1;
      _zdata[col].insert({v,1});
//    if (_zdata[col].find(v) == _zdata[col].end())
//    {
//      cerr << "hash insertion failed, col " << col << ", size " << _zdata[col].size() << " value " << v << endl;
//    }
    }
    else
    {
      (*f).second++;
//    cerr << "col " << col << " updated value " << _zdata[col][v] << " for " << v << endl;
    }
  }

  return 1;
}

int
Data_From_File::_determine_columns_to_process (const const_IWSubstring & buffer)
{
  int i = 0;
  IWString token;

  int rc = 0;
  for (int col = 0; buffer.nextword(token, i); col++)
  {
    _descriptor_name[col] = token;

    if (_descriptor.contains(token))
    {
      _process_column[col] = 1;
      if (verbose)
        cerr << "descriptor '" << token << " found in column " << (col + 1) << endl;

      rc++;
    }
  }

  if (0 == _descriptor.size() && _descriptor_processing.size())
  {
    for (auto i : _descriptor_processing)
    {
      const auto f = std::find(_descriptor_name, _descriptor_name + _columns_in_input, i.first);

//    cerr << "Looking for '" << i.first << "'\n";

      if (f == _descriptor_name + _columns_in_input)
        continue;

      _process_column[f - _descriptor_name] = 1;

//    cerr << "FROM FILE found in col " << (f - _descriptor_name) << endl;

      rc++;
    }
  }

// First set whatever information we may have read from a file

  for (auto i : _descriptor_processing)
  {
    const auto f = std::find(_descriptor_name, _descriptor_name + _columns_in_input, i.first);

    if (f == _descriptor_name + _columns_in_input)
      continue;

    const auto col = f - _descriptor_name;
//  cerr << "Initialising rounding for '" << i.first << "' found in column " << col << endl;

    i.second->initialise(_rounding[col], _round_to_nearest_int[col], _y_normalisation[col]);
  }

// then the command line

  if (_rounding_from_cmdline.number_elements() > 0)
  {
    if (! transfer_rounding_from_cmdline(_columns_in_input, _process_column, _rounding_from_cmdline, _rounding))
    {
      return 0;
    }
  }
  else if (_round_to_nearest_int_from_cmdline.number_elements() > 0)
  {
    if (! transfer_rounding_from_cmdline(_columns_in_input, _process_column, _round_to_nearest_int_from_cmdline, _round_to_nearest_int))
    {
      return 0;
    }
  }

  if (verbose > 1)
  {
    for (int i = 0; i < _columns_in_input; ++i)
    {
      if (0 == _process_column[i])
        continue;
  
      cerr << "Column " << i << " " << _descriptor_name[i] << " rounding " << _rounding[i] << " nearest int " << _round_to_nearest_int[i] << endl;
    }
  }

  if (static_cast<size_t>(rc) < _descriptor.size())
  {
    cerr << "Data_From_File::_determine_columns_to_process:warning only found " << rc << " of " << _descriptor.size() << " descriptors\n";
    return 0;
  }

  return rc;
}

int
Data_From_File::_initialise_arrays (const const_IWSubstring & buffer)
{
  _columns_in_input = buffer.nwords_single_delimiter(_input_delimiter);
  if (_columns_in_input < 1)
  {
    cerr << "Data_From_File::_initialise_arrays:no tokens in input\n";
    return 0;
  }

  _process_column = new_int(_columns_in_input);
  _acc = new Accumulator<double>[_columns_in_input];
  _zdata = new std::unordered_map<float, int>[_columns_in_input];

  _rounding = new float[_columns_in_input];
  std::fill_n(_rounding, _columns_in_input, 0.0f);
  _round_to_nearest_int = new_int(_columns_in_input);

  _y_normalisation = new_int(_columns_in_input);

  transfer_columns_from_cmdline(_columns_specified_on_command_line, _process_column);

  if (_descriptor.size() || _descriptor_processing.size())
    _descriptor_name = new IWString[_columns_in_input];

  return 1;
}

int
Data_From_File::build (iwstring_data_source & input)
{
  const_IWSubstring buffer;
  
  for (int i = 0; i < _header_records_to_skip; i++)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Cannot read header record\n";
      return 0;
    }

    if (i > 0)
      continue;

    if (! _initialise_arrays(buffer))
      return 0;

    if (_descriptor.size() > 0 || _descriptor_processing.size() > 0)
    {
      if (! _determine_columns_to_process(buffer))
      {
        cerr << "Cannot determine column for one or more descriptors\n";
        return 0;
      }
    }
  }

  while (input.next_record(buffer))
  {
    if (_columns_in_input <= 0)
    {
      if (! _initialise_arrays(buffer))
        return 0;
    }

    if (! _extract_data(buffer))
    {
      cerr << "Cannot extract data from '" << buffer << "', line " << input.lines_read() << endl;
      return 0;
    }
  }

  return input.lines_read();
}

int
Data_From_File::build(const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input);
}

int
Data_From_File::_write_stats_file (const IWString & stats_file_stem,
                                   const int col) const
{

  IWString fname;

  fname << stats_file_stem;
  // Only decorate the file name, if the user had not already come up with something interesting
  if (!fname.find(".stat")) fname << _descriptor_name[col] << ".stats";

  IWString_and_File_Descriptor statfile;
  
  if (! statfile.open(fname.c_str()))
  {
    cerr << "Yipes, cannot open stats file '" << fname << "'\n";
    return 3;
  }

  statfile << "#MPP min " << _acc[col].minval() << " max " << _acc[col].maxval() << " ave " << static_cast<float>(_acc[col].average()) << "\n";

  return 1;
}

int
Data_From_File::_remove_common_descriptor_prefix()
{
  for (int i = 0; i < _columns_in_input; ++i)
  {
    if (0 == _process_column[i])
      continue;

    if (! _descriptor_name[i].starts_with(_common_descriptor_prefix))
      continue;

    _descriptor_name[i].remove_leading_chars(_common_descriptor_prefix.length());
  }

  return 1;
}

int
Data_From_File::doit (const IWString & stats_file_stem,
                      Rfile & rfile,
                      const IWString & output_file_stem)
{
  int cbp = columns_being_processed();

  if (0 == cbp)
  {
    cerr << "Data_From_File::doit:no columns being processed\n";
    return 0;
  }

  for (int i = 0; i < _columns_in_input; ++i)
  {
    if (0 == _process_column[i])
      continue;

    if (! _final_processing(i, rfile, stats_file_stem, output_file_stem))
    {
      cerr << "Data_From_File::doit:fatal error processing column " << i;
      if (nullptr != _descriptor_name)
        cerr << ' ' << _descriptor_name[i];
      cerr << endl;
      return 0;
    }
  }

  return 1;
}

int
Data_From_File::_final_processing (const int col,
                                   Rfile & rfile,
                                   const IWString & stats_file_stem,
                                   const IWString & output_file_stem)
{
  if (_common_descriptor_prefix.length() > 0)
    _remove_common_descriptor_prefix();

  resizable_array_p<VCF> vcf;

  const auto & d = _zdata[col];

  vcf.resize(d.size());

  int n = 0;
  int highest = 0;

  for (const auto i : d)
  {
    float v = i.first;
    int c = i.second;

//  cerr << " gather raw " << v << " " << c << endl;

    n += c;
    if (c > highest)
      highest = c;

    VCF * tmp = new VCF(v, c);

    vcf.add(tmp);
  }

  Value_Count_Comparator vcc;

  vcf.iwqsort(vcc);

  if (verbose)
    cerr << " col " << col << ' ' <<  n << " samples, largest count " << highest << endl;

#ifdef ECHO_SORTED_VALUES_QWE
  for (int i = 0; i < vcf.number_elements(); ++i)
  {
    cerr << " value " << vcf[i]->v() << " count " << vcf[i]->count() << endl;
  }
#endif

  if (stats_file_stem.length() && ! _write_stats_file(stats_file_stem, col))
  {
    cerr << "Data_From_File::_final_processing:cannot write stats file for '" << _descriptor_name[col] << "'\n";
    return 0;
  }

  if (rfile.active())
  {
    rfile.statistics(_descriptor_name[col], _acc[col]);
    _transfer_data_to_rfile (rfile, col, n, highest, vcf);
  }

  if (1 == columns_being_processed() && 0 == output_file_stem.length()) 
  {
    IWString_and_File_Descriptor output(1);
    return _write_sorted_values(vcf, col, n, highest, output);
  }

  IWString fname;

  if (nullptr != _descriptor_name)
    fname << output_file_stem << _descriptor_name[col] << ".dat";
  else
    fname << output_file_stem << col << ".dat";

  IWString_and_File_Descriptor output;

  if (! output.open(fname.c_str()))
  {
    cerr << "Data_From_File::doit:cannot open '" << fname << "'\n";
    return 0;
  }

  return _write_sorted_values(vcf, col, n, highest, output);
}

int
Data_From_File::_write_normalised_value (const int c,
                                         const int col,
                                         const int n,
                                         const int highest,
                                         IWString_and_File_Descriptor & output) const
{
  typedef float ftype;

  if (NORMALISE_DISTRIBUTION_TO_UNIT_AREA == _y_normalisation[col])
    output << static_cast<ftype>(c) / static_cast<ftype>(n);
  else if (NORMALISE_DISTRIBUTION_TO_HIGHEST == _y_normalisation[col])
    output << static_cast<ftype>(c) / static_cast<ftype>(highest);
  else
    output << c;

  return 1;
}

template <typename T>
int
Data_From_File::_append_normalised_value (const int c,
                                          const int col,
                                          const int n,
                                          const int highest,
                                          resizable_array<T> & y) const
{
  typedef float ftype;

  if (NORMALISE_DISTRIBUTION_TO_UNIT_AREA == _y_normalisation[col])
    y.add(static_cast<ftype>(c) / static_cast<ftype>(n));
  else if (NORMALISE_DISTRIBUTION_TO_HIGHEST == _y_normalisation[col])
    y.add(static_cast<ftype>(c) / static_cast<ftype>(highest));
  else
    y.add(c);

  return 1;
}

int
Data_From_File::_write_sorted_values (const resizable_array_p<VCF> & vcf, const int col,
                                      const int n,
                                      const int highest,
                                      IWString_and_File_Descriptor & output) const
{
  IWString colstr = "col" + std::to_string(col);
  if (nullptr != _descriptor_name) {
    colstr = _descriptor_name[col];
    IWString desc = _descriptor_name[col].substr(_descriptor_name[col].find("_") + 1);
    if (desc.length() > 0 && _collection_name.find(desc))
      // We already have an interesting column name with the descriptor embedded
      colstr = "";
  }
  output << "D " << _collection_name << colstr << '\n';

  if (0 == _insert_zero_values_for_missing)
    ;
  else if (_rounding[col] > 0.0f)
    return _write_sorted_values_insert_zero_float(vcf, col, n, highest, output);
  else if (_round_to_nearest_int[col] > 0)
    return _write_sorted_values_insert_zero_int(vcf, col, n, highest, output);

  for (int i = 0; i < vcf.number_elements(); i++)
  {
    const VCF * vi = vcf[i];

    output << vi->v() << ' ';

    _write_normalised_value(vi->count(), col, n, highest, output);
    
    output << '\n';

    output.write_if_buffer_holds_more_than(8192);
  }

  output.flush();

  return 1;
}

/*
  We have rounded data, and need to 'fill in the blanks'.
  To do that, we need to figure out which 'whole number' each
  of our values correspond to
*/

int
Data_From_File::_write_sorted_values_insert_zero_float (const resizable_array_p<VCF> & vcf,
                                                        const int col,
                                                        const int n,
                                                        const int highest,
                                                        IWString_and_File_Descriptor & output) const
{
  int iprev = std::numeric_limits<int>::max() - 1;

  for (int i = 0; i < vcf.number_elements(); ++i)
  {
    const VCF * vi = vcf[i];

    auto v = vi->v();

    int j = static_cast<int>(v * _rounding[col] +  (v < 0.0f ? -0.01 : 0.01));
//  cerr << "Value " << v << " rounding " << _rounding[col] << " int " << j << endl;

    for (int k = iprev + 1; k < j; ++k)    // may not execute at all
    {
      output << (static_cast<float>(k) / _rounding[col]) << " 0\n";
    }

    output << vi->v() << ' ';

//  cerr << _descriptor_name[col] << " value " <<  v << " count " << vi->count() << endl;

    _write_normalised_value(vi->count(), col, n, highest, output);
    
    output << '\n';

    output.write_if_buffer_holds_more_than(8192);

    iprev = j;
  }

  output.flush();

  return 1;
}


int
Data_From_File::_write_sorted_values_insert_zero_int (const resizable_array_p<VCF> & vcf,
                                                        const int col,
                                                        const int n,
                                                        const int highest,
                                                        IWString_and_File_Descriptor & output) const
{
  int iprev = std::numeric_limits<int>::max() - _round_to_nearest_int[col];

  for (int i = 0; i < vcf.number_elements(); ++i)
  {
    const VCF * vi = vcf[i];

    int j = static_cast<int>(vi->v() + (vi->v() < 0.0f ? -0.01 : 0.01));

//  cerr << "INT value " << vi->v() << " int " << j << " prev " << iprev << endl;

    for (int k = iprev + _round_to_nearest_int[col]; k < j; k += _round_to_nearest_int[col])
    {
      output << k << " 0\n";
    }

    output << vi->v() << ' ';
    
    _write_normalised_value(vi->count(), col, n, highest, output);
      
    output << '\n';

    output.write_if_buffer_holds_more_than(8192);

    iprev = j;
  }

  output.flush();

  return 1;
}

int
Data_From_File::_transfer_data_to_rfile (Rfile & rfile,
                                         const int col,
                                         const int n,
                                         const int highest,
                                         const resizable_array_p<VCF> & vcf) const
{
  if (0 == _insert_zero_values_for_missing)
    return rfile.data(_descriptor_name[col], vcf);

  if (_rounding[col] > 0.0f)
    return _transfer_data_to_rfile_float(rfile, col, n, highest, vcf);
  if (_round_to_nearest_int[col] > 0)
    return _transfer_data_to_rfile_int(rfile, col, n, highest, vcf);

  cerr << "Data_From_File::_transfer_data_to_rfile_int:huh, what am I supposed to do?\n";
  return 0;
}

int
Data_From_File::_transfer_data_to_rfile_float (Rfile & rfile,
                                               const int col,
                                               const int n,
                                               const int highest,
                                               const resizable_array_p<VCF> & vcf) const
{
  resizable_array<float> x;
  resizable_array<float> y;

  int iprev = std::numeric_limits<int>::max() - 1;

  for (int i = 0; i < vcf.number_elements(); ++i)
  {
    const VCF * vi = vcf[i];

    auto v = vi->v();

    int j = static_cast<int>(v * _rounding[col] +  (v < 0.0f ? -0.01 : 0.01));
//  cerr << "Value " << v << " rounding " << _rounding[col] << " int " << j << endl;

    for (int k = iprev + 1; k < j; ++k)    // may not execute at all
    {
      x.add(static_cast<float>(k) / _rounding[col]);
      y.add(0.0f);
    }

    x.add(vi->v());

    _append_normalised_value(vi->count(), col, n, highest, y);

    iprev = j;
  }

  return rfile.data(_descriptor_name[col], x, y);
}

int
Data_From_File::_transfer_data_to_rfile_int (Rfile & rfile,
                                             const int col,
                                             const int n,
                                             const int highest,
                                             const resizable_array_p<VCF> & vcf) const
{
  resizable_array<int> x;
  resizable_array<float> y;

  int iprev = std::numeric_limits<int>::max() - _round_to_nearest_int[col];

  for (int i = 0; i < vcf.number_elements(); ++i)
  {
    const VCF * vi = vcf[i];

    int j = static_cast<int>(vi->v() + (vi->v() < 0.0f ? -0.01 : 0.01));

//  cerr << "INT value " << vi->v() << " int " << j << " prev " << iprev << endl;

    for (int k = iprev + _round_to_nearest_int[col]; k < j; k += _round_to_nearest_int[col])
    {
      x.add(k);
      y.add(0.0f);
    }

    x.add(vi->v());
    
    _append_normalised_value(vi->count(), col, n, highest, y);

    iprev = j;
  }

  return rfile.data(_descriptor_name[col], x, y);
}

static int
distribution (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:c:d:r:ahn:H:S:i:O:D:R:zP:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  Data_From_File mydata;

  if (cl.option_present('c') && cl.option_present('d'))
  {
    cerr << "Sorry, cannot use both -c and -d options\n";
    usage(3);
  }

  if (cl.option_present('c') && ! mydata.set_columns_to_process(cl, 'c'))
  {
    cerr << "Cannot determine columns to process (-c)\n";
    usage(1);
  }

  if (cl.option_present('D'))
  {
    const char * d = cl.option_value('D');
    
    if (! mydata.read_descriptor_processing(d))
    {
      cerr << "Cannot initialise per descriptor handling '" << d << "'\n";
      return 1;
    }
  }

  if (cl.option_present('d') && ! mydata.set_descriptors_to_process(cl, 'd') && ! cl.option_present('D'))
  {
    cerr << "Cannot determine descriptors to process (-d)\n";
    usage(1);
  }

  if (cl.option_present('s'))
  {
    int s;
    if (! cl.value('s', s) || s < 0)
    {
      cerr << "The header records to skip option (-s) must be a whole +ve number\n";
      usage(3);
    }

    mydata.set_header_records_to_skip(s);

    if (verbose)
      cerr << "Will skip " << s << " header records\n";
  }
  else if (cl.option_present('d') || cl.option_present('D'))
    mydata.set_header_records_to_skip(1);

  if (cl.option_present('r') && cl.option_present('n'))
  {
    cerr << "Sorry, the -r and -n options are mutually exclusive\n";
    usage(3);
  }

  if (cl.option_present('r') && ! mydata.set_per_column_rounding(cl, 'r'))
  {
    cerr << "Cannot establish per column rounding (-r)\n";
    usage(1);
  }

  if (cl.option_present('n') && ! mydata.set_per_column_integer_rounding(cl, 'n'))
  {
    cerr << "Cannot establish per column integer rounding (-n)\n";
    usage(1);
  }

  if (cl.option_present('h') && cl.option_present('a'))
  {
    cerr << "Sorry, the -a and -h options are mutually exclusive\n";
  }

  if (cl.option_present('h'))
  {
    if (verbose)
      cerr << "Will normalise distribution to highest value\n";
  }

  if (cl.option_present('a'))
  {
    if (verbose)
      cerr << "Will normalise distribution to unit area\n";
  }

  if (cl.option_present('z'))
  {
    mydata.set_insert_zero_values_for_missing(1);

    if (verbose)
      cerr << "With rounded data, will insert 0 at points with no data\n";
  }

  if (cl.option_present('P'))
  {
    const char * p = cl.option_value('P');

    mydata.set_common_descriptor_prefix(p);

    if (verbose)
      cerr << "Common prefix '" << p << "' removed from descriptor names\n";
  }

  if (cl.option_present('H'))
  {
    const char * h = cl.option_value('H');

    mydata.set_collection_name(h);

    if (verbose)
      cerr << "Collection prefix '" << h << "' added to descriptor names\n";
  }

  IWString stats_file_stem;

  if (cl.option_present('S'))
  {
    stats_file_stem = cl.option_value('S');

    if (verbose)
      cerr << "Will create stats file for molecular_property_profile\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  Rfile rfile;

  if (cl.option_present('R'))
  {
    if (! rfile.build(cl, 'R', verbose))
    {
      cerr << "Cannot initialise R programme setup (-R)\n";
      return 1;
    }
  }

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! mydata.build(cl[i]))
    {
      cerr << "Cannot read data from '" << cl[i] << "'\n";
      return (i+1);
    }
  }

  if (cl.option_present('h'))
  {
    mydata.set_normalisation(NORMALISE_DISTRIBUTION_TO_HIGHEST);

    if (verbose)
      cerr << "Output normalised to highest value in each column\n";
  }
  else if (cl.option_present('a'))
  {
    mydata.set_normalisation(NORMALISE_DISTRIBUTION_TO_UNIT_AREA);

    if (verbose)
      cerr << "Output normalised to unit area in each column\n";
  }

  if (verbose)
    mydata.report(cerr);

  if (mydata.columns_being_processed() > 1 && ! cl.option_present('O'))
  {
    cerr << "Multiple columns being processed, must specify output output file name stem via the -O option\n";
    usage(1);
  }

  IWString output_file_stem = cl.string_value('O');

  mydata.doit(stats_file_stem, rfile, output_file_stem);

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = distribution(argc, argv);

  return rc;
}
