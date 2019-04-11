/*
  We have a fingerprint file with sparse fingerprints. Can we convert them into either
  fixed counted, or binary fingerprints
*/

#include <limits>
#include <memory>
#include <ostream>

#include "vector"
#include <algorithm>
#include <functional>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "misc.h"
#include "iw_stl_hash_map.h"
#include "iw_auto_array.h"

#include "gfp.h"

const char * prog_name = NULL;

static int verbose = 0;

static int produce_fixed_binary = 1;
static int produce_fixed_counted = 0;

static int tdts_read = 0;

static int echo_changed_fingerprint = 0;

static vector<int> numeric_threshold;

static int work_as_filter = 0;

static resizable_array_p<IWString> file_records;

static int need_to_determine_fingerprints_present = 1;

static int user_specified_fingerprint_width = std::numeric_limits<int>::max();

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Converts non-colliding fingerprints to fixed counted, or binary forms\n";
  cerr << " -F <tag>       input tag(s) - probably NCXX\n";
  cerr << " -b             produce fixed binary  fingerprints\n";
  cerr << " -f             produce fixed counted fingerprints\n";
  cerr << " -B <nbits>     produce fixed width fingerprints (achieved by hashing)\n";
  cerr << " -p <num>       omit bits set fewer than <num> times. Use xx% for percent\n";
  cerr << " -P <num>       omit bits set more  than <num> times. Use xx% for percent\n";
  cerr << " -e             echo changed fingerprint\n";
  cerr << " -f             work as filter - must hold entire input in RAM\n";
  cerr << " -k n,n,n       numeric thresholds for producing different fingerprints\n";
  cerr << " -S <fname>     write bit translation table to <fname>\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

template <typename T>
class Number_or_Percentage
{
  private:
    T _value;
    int _percentage;
    
    static const int NOP_NOT_INITIALISED;
    static const int NOP_NUMERIC_VALID;

  public:
    Number_or_Percentage();

    int build(const_IWSubstring);   // our local copy

    T numeric_equivalent(T) const;
};

template<> const int Number_or_Percentage<int>::NOP_NOT_INITIALISED = -1;
template<> const int Number_or_Percentage<int>::NOP_NUMERIC_VALID = 0;

template <typename T>
Number_or_Percentage<T>::Number_or_Percentage()
{
  _percentage = NOP_NOT_INITIALISED;

  return;
}

template <typename T>
int
Number_or_Percentage<T>::build(const_IWSubstring buffer)   // note not a reference
{
  if (buffer.ends_with('%'))
  {
    buffer.chop();

    if (! buffer.numeric_value(_percentage) || _percentage <= 0 || _percentage >= 100)
    {
      cerr << "Number_or_Percentage::build: invalid percentage '" << buffer << "'\n";
      return 0;
    }

    return 1;
  }

  if (! buffer.numeric_value(_value))
  {
    cerr << "Number_or_Percentage::build:invalid numeric '" << buffer << "'\n";
    return 0;
  }

  _percentage = NOP_NUMERIC_VALID;

  return 1;
}

template <typename T>
T
Number_or_Percentage<T>::numeric_equivalent(T v) const
{
  if (NOP_NUMERIC_VALID == _percentage)
    return _value;

  assert (NOP_NOT_INITIALISED != _percentage);

  return static_cast<T>(static_cast<float>(v) * static_cast<float>(_percentage) / 100.0f + 0.499f);
}

typedef std::unordered_map<unsigned int, unsigned int> Bit_Data;

/*
  For each non colliding fingerprint we need a mapping of bit numbers
*/

class Sparse_to_Fixed_Bit_Mapping
{
  private:
    int _number_thresholds;

//  For each threshold, we will have a count of how many molecules have that bit set enough times

    Bit_Data * _molecules_containing_bit;

//  An array used for the bit values

    int * _btmp;

//  For every threshold used, we need a mapping

    Bit_Data * _sparse_bit_number_to_fixed_bit_number;
    int * _nbits;

//  When we produce multiple fingerprints, we need a temporary array to
//  form those fingerprints

    int * _above_threshold;

    IWString _fpname;

//  private functions

    int _do_write_fixed_binary (const int * bits, int nb, const IWString & tag, int ndx, IWString_and_File_Descriptor & output) const;
    int _do_write_fixed_counted(const int * bits, int nb, const IWString & tag, int ndx, IWString_and_File_Descriptor & output) const;

  public:
    Sparse_to_Fixed_Bit_Mapping(const IWString &);
    ~Sparse_to_Fixed_Bit_Mapping();

    void set_fpname(const IWString &);
    const IWString & fpname () const { return _fpname;}

    int extra(const const_IWSubstring & buffer);
    int extra(const Sparse_Fingerprint & sfp);

    template <typename T> void do_trim(const Number_or_Percentage<int> & cutoff, T & comparitor);

    int setup_bit_mapping();
    
    int do_write(const Sparse_Fingerprint & sfp, const IWString & tag, IWString_and_File_Descriptor & output);
    int do_write(const const_IWSubstring & buffer, const IWString & tag, IWString_and_File_Descriptor & output);

    template <typename O> int write_bit_mapping(const char sep, O & output) const;
    int write_bit_mapping(const char sep, const char * fname) const;
};

typedef resizable_array_p<Sparse_to_Fixed_Bit_Mapping> Set_of_Bit_Mappings;

Sparse_to_Fixed_Bit_Mapping::Sparse_to_Fixed_Bit_Mapping(const IWString & s)
{
  _number_thresholds = numeric_threshold.size();

  assert (_number_thresholds > 0);

  _btmp = NULL;

  _sparse_bit_number_to_fixed_bit_number = new Bit_Data[_number_thresholds];

  _molecules_containing_bit = new Bit_Data[_number_thresholds];

  _nbits = new_int(_number_thresholds);

  _above_threshold = NULL;

  set_fpname(s);

  return;
}

Sparse_to_Fixed_Bit_Mapping::~Sparse_to_Fixed_Bit_Mapping()
{
  if (NULL != _btmp)
    delete [] _btmp;

  if (NULL != _above_threshold)
    delete [] _above_threshold;

  if (NULL != _sparse_bit_number_to_fixed_bit_number)
    delete [] _sparse_bit_number_to_fixed_bit_number;

  if (NULL != _nbits)
    delete [] _nbits;

  if (NULL != _molecules_containing_bit)
    delete [] _molecules_containing_bit;

  return;
}

void
Sparse_to_Fixed_Bit_Mapping::set_fpname(const IWString & s)
{
  _fpname = s;

  if (_fpname.starts_with("NC"))
    _fpname.remove_leading_chars(2);

  if (_fpname.ends_with('<'))
    _fpname.chop();
}

int
Sparse_to_Fixed_Bit_Mapping::extra(const const_IWSubstring & buffer)
{
  Sparse_Fingerprint sfp;

  (void) sfp.construct_from_tdt_record(buffer);

  return extra(sfp);
}

int
Sparse_to_Fixed_Bit_Mapping::extra(const Sparse_Fingerprint & sfp)
{
  int i = 0;
  unsigned int ibit;
  int icount;

//cerr << "Considering fingerprint with " << sfp.nbits() << " bits\n";

  while (sfp.next_bit_set(i, ibit, icount))
  {
//  cerr << " i = " << i << " bit " << ibit << " icount " << icount << endl;

    for (int j = 0; j < _number_thresholds; j++)
    {
      const int t = numeric_threshold[j];

      if (icount < t)
        continue;

      Bit_Data::iterator f = _molecules_containing_bit[j].find(ibit);

      if (_molecules_containing_bit[j].end() == f)
        _molecules_containing_bit[j][ibit] = 1;
      else 
        _molecules_containing_bit[j][ibit]++;
    }
  }

  return 1;
}

template <typename T>
void
Sparse_to_Fixed_Bit_Mapping::do_trim(const Number_or_Percentage<int> & cutoff,
                                      T & comparitor)
{

  for (int i = 0; i < _number_thresholds; i++)
  {
    int threshold = cutoff.numeric_equivalent(_molecules_containing_bit[i].size());
//  cerr << "For " << _molecules_containing_bit[i].size() << " bits, threshold is " << threshold << endl;

    int trimmed = 0;

    for (Bit_Data::iterator s = _molecules_containing_bit[i].begin(); s != _molecules_containing_bit[i].end(); ++s)
//  for (auto s : _molecules_containing_bit[i])
    {
      unsigned int & c = (*s).second;

      if (comparitor(c, threshold))
      {
        c = 0;
        (*s).second = 0;
        trimmed++;
      }
    }

    if (verbose)
      cerr << _fpname << " threshold " << numeric_threshold[i] << " trimmed " << trimmed << " bits\n";
  }

  return;
}

int
Sparse_to_Fixed_Bit_Mapping::setup_bit_mapping()
{
//cerr << "Sparse_to_Fixed_Bit_Mapping::setup_bit_mapping:_number_thresholds " << _number_thresholds << endl;

  for (int i = 0; i < _number_thresholds; i++)
  {
    for (Bit_Data::const_iterator s = _molecules_containing_bit[i].begin(); s != _molecules_containing_bit[i].end(); ++s)
    {
      unsigned int c = (*s).second;
//    cerr << " bit " << s->first << " count " << s->second << endl;

      if (0 == c)
        continue;

      int x = _sparse_bit_number_to_fixed_bit_number[i].size();

      _sparse_bit_number_to_fixed_bit_number[i][(*s).first] = x;
    }

    _nbits[i] = _sparse_bit_number_to_fixed_bit_number[i].size();

    if (verbose)
      cerr << "Fingerprint '" << _fpname << "' threshold " << numeric_threshold[i] << " initial bits " << _molecules_containing_bit[i].size() << " number trimmed " << (_molecules_containing_bit[i].size() - _nbits[i]) << " final nbits " << _nbits[i] << endl;

    if (produce_fixed_binary)
    {
      if (0 != _nbits[i] % IW_BITS_PER_BYTE)
        _nbits[i] = (_nbits[i] / IW_BITS_PER_BYTE + 1) * IW_BITS_PER_BYTE;
    }
  }

  int highest_nbits = iwmax_of_array(_nbits, _number_thresholds);

  _btmp = new int[highest_nbits];

  return 1;
}

int
Sparse_to_Fixed_Bit_Mapping::do_write(const Sparse_Fingerprint & sfp,
                                       const IWString & tag,
                                       IWString_and_File_Descriptor & output)
{
  unsigned int b;
  int icount;

  for (int i = 0; i < _number_thresholds; i++)
  {
    set_vector(_btmp, _nbits[i], 0);

    int j = 0;
    while (sfp.next_bit_set(j, b, icount))
    {
      Bit_Data::const_iterator fc = _molecules_containing_bit[i].find(b);
      if (fc == _molecules_containing_bit[i].end())
        continue;

      if (0 == (*fc).second)
        continue;

      Bit_Data::const_iterator fi = _sparse_bit_number_to_fixed_bit_number[i].find(b);

      int ndx = (*fi).second;

      if (ndx > user_specified_fingerprint_width)
        ndx = ndx % user_specified_fingerprint_width;

      assert (ndx >= 0 && ndx <= _nbits[i]);

      _btmp[ndx] = 1;
    }

    int t = numeric_threshold[i];
    if (1 == _number_thresholds)
      t = -1;

    int nb = _nbits[i];
    if (nb > user_specified_fingerprint_width)
      nb = user_specified_fingerprint_width;

    if (produce_fixed_binary)
      _do_write_fixed_binary(_btmp, nb, tag, t, output);
    else if (produce_fixed_counted)
      _do_write_fixed_counted(_btmp, nb, tag, t, output);
  }

  return 1;
}

int
Sparse_to_Fixed_Bit_Mapping::do_write(const const_IWSubstring & buffer,
                                       const IWString & tag,
                                       IWString_and_File_Descriptor & output)
{
  Sparse_Fingerprint sfp;

  (void) sfp.construct_from_tdt_record(buffer);

  return do_write(sfp, tag, output);
}

int
Sparse_to_Fixed_Bit_Mapping::_do_write_fixed_binary(const int * bits,
                                        int nb,
                                        const IWString & tag,
                                        int ndx,
                                        IWString_and_File_Descriptor & output) const
{
  IW_Bits_Base b;

  b.construct_from_array_of_ints(bits, nb);

  IWString tmp;
  b.daylight_ascii_representation_including_nset_info(tmp);

  IWString output_tag(tag);
  output_tag[0] = 'F';
  output_tag[1] = 'P';

  if (ndx >= 0)
  {
    if (output_tag.ends_with('<'))
      output_tag.chop();

    output_tag << ndx << '<';
  }

  output << output_tag << tmp << ">\n";

  return 1;
}

int
Sparse_to_Fixed_Bit_Mapping::_do_write_fixed_counted(const int * bits,
                                 int nb,
                                 const IWString & tag,
                                 int ndx,
                                 IWString_and_File_Descriptor & output) const
{
  Fixed_Size_Counted_Fingerprint_uchar fp;

  fp.construct_from_array_of_ints(bits, nb);

  IWString output_tag(tag);
  output_tag[0] = 'F';
  output_tag[1] = 'C';

  if (ndx >= 0)
  {
    if (output_tag.ends_with('<'))
      output_tag.chop();

    output_tag << ndx << '<';
  }

  assert (NULL == "Implment fixed size counted fingerprints sometime, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)");

  return 0;
}

template <typename O>
int
Sparse_to_Fixed_Bit_Mapping::write_bit_mapping(const char sep, O & output) const
{
  for (int i = 0; i < _number_thresholds; ++i)
  {
    output << "FP " << _fpname << " thr " << i << '\n';

    auto & mcb = _molecules_containing_bit[i];     // is actually const, but we use operator [] below

    for (auto j : _sparse_bit_number_to_fixed_bit_number[i])
    {
      output << j.first << sep << j.second << sep << mcb[j.first] << '\n';
      output.write_if_buffer_holds_more_than(4096);
    }
    output << "|\n";
  }

  return 1;
}

int
Sparse_to_Fixed_Bit_Mapping::write_bit_mapping(const char sep, const char * fname) const
{
  IWString_and_File_Descriptor output;
  
  if (! output.open(fname))
  {
    cerr << "Sparse_to_Fixed_Bit_Mapping:;write_bit_mapping:cannot open '" << fname << "'\n";
    return 0;
  }

  return write_bit_mapping(sep, output);
}

static int
identify_tag(const const_IWSubstring & buffer,
              IWString & tag)
{
  int open_angle_bracket = buffer.index('<');

  if (open_angle_bracket < 0)   // should be fatal error
    return 0;

  buffer.from_to(0, open_angle_bracket, tag);

  return 1;
}

static int
sparse_to_fixed_size_buffer(const const_IWSubstring & buffer,
                             const IW_STL_Hash_Map_int & tag_to_ndx,
                             const Set_of_Bit_Mappings & stfbm,
                             IWString_and_File_Descriptor & output)
{
  IWString tag;

  if (! identify_tag(buffer, tag))
  {
    output << buffer << '\n';
    return 1;
  }

  IW_STL_Hash_Map_int::const_iterator f = tag_to_ndx.find(tag);

  if (f == tag_to_ndx.end())   // we are not altering this fingerprint
  {
    output << buffer << '\n';
    return 1;
  }

  if (echo_changed_fingerprint)
    output << buffer << '\n';

  int ndx = (*f).second;

  stfbm[ndx]->do_write(buffer, tag, output);

  output.write_if_buffer_holds_more_than(4096);

  return output.good();
}

static int
sparse_to_fixed_size(const resizable_array_p<IWString> & file_records,
                      const IW_STL_Hash_Map_int & tag_to_ndx,
                      const Set_of_Bit_Mappings & stfbm,
                      IWString_and_File_Descriptor & output)
{
  int n = file_records.number_elements();

  for (int i = 0; i < n; i++)
  {
    const_IWSubstring s = *(file_records[i]);

    sparse_to_fixed_size_buffer(s, tag_to_ndx, stfbm, output);

    output.write_if_buffer_holds_more_than(4096);
  }

  return output.good();
}

static int
sparse_to_fixed_size(iwstring_data_source & input,
                      const IW_STL_Hash_Map_int & tag_to_ndx,
                      const Set_of_Bit_Mappings & stfbm,
                      IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer) && output.good())
  {
    sparse_to_fixed_size_buffer(buffer, tag_to_ndx, stfbm, output);

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
sparse_to_fixed_size(const char * fname,
                      const IW_STL_Hash_Map_int & tag_to_ndx,
                      const Set_of_Bit_Mappings & stfbm,
                      IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return sparse_to_fixed_size(input, tag_to_ndx, stfbm, output);
}

static int
profile_fingerprints(iwstring_data_source & input,
                      IW_STL_Hash_Map_int & tag_to_ndx,
                      Set_of_Bit_Mappings & stfbm)
{
  const_IWSubstring buffer;

  IWString tag;

  while (input.next_record(buffer))
  {
//  cerr << "Profiling read '" << buffer << "', need " << need_to_determine_fingerprints_present << endl;
    if (work_as_filter)
      file_records.add(new IWString(buffer));

    if ('|' == buffer)
    {
      tdts_read++;
      need_to_determine_fingerprints_present = 0;
      continue;
    }

    if (! identify_tag(buffer, tag))
      continue;

    if (! tag.starts_with("NC"))
      continue;

    if (need_to_determine_fingerprints_present)
    {
      int s = tag_to_ndx.size();
      tag_to_ndx[tag] = s;
      stfbm.add(new Sparse_to_Fixed_Bit_Mapping(tag));
      if (verbose)
        cerr << "Added fingerprint '" << tag << "'\n";
    }

    IW_STL_Hash_Map_int::const_iterator f = tag_to_ndx.find(tag);

//  cerr << "TAG " << tag << " iterator " << (f == tag_to_ndx.end()) << endl;

    if (f == tag_to_ndx.end())
      continue;

    int ndx = (*f).second;

//  cerr << "Index is " << ndx << endl;

    stfbm[ndx]->extra(buffer);
  }

  return 1;
}

//#define WRITING_CROSS_REFERENCE_FILE_DEFINED
#ifdef WRITING_CROSS_REFERENCE_FILE_DEFINED

static int
write_cross_reference_file(IWString_and_File_Descriptor & output,
                            const Bit_Data & fp_to_index,
                            const Bit_Data & fp_count)
{
  for (Bit_Data::const_iterator i = fp_count.begin(); i != fp_count.end(); ++i)
  {
    unsigned int c = (*i).second;

    if (0 == c)
      continue;

    unsigned int b = (*i).first;

    Bit_Data::const_iterator f = fp_to_index.find(b);

    output << b << ' ' << (*f).second << ' ' << c << '\n';

    output.write_if_buffer_holds_more_than(4096);
  }

  return output.good();
}

static int
write_cross_reference_file(IWString & fname,
                            const Bit_Data & fp_to_index,
                            const Bit_Data & fp_count)
{
  IWString_and_File_Descriptor output;
  
  if (! output.open(fname))
  {
    cerr << "Cannot open cross reference file '" << fname << "'\n";
    return 0;
  }

  return write_cross_reference_file(output, fp_to_index, fp_count);
}
#endif

template <typename T>
void
do_trim(int nfp,
         Set_of_Bit_Mappings & stfbm,
         const Number_or_Percentage<int> & cutoff,
         T & comparitor)
{
  for (int i = 0; i < nfp; i++)
  {
    stfbm[i]->do_trim(cutoff, comparitor);
  }

  return;
}

static int
profile_fingerprints(const char * fname,
                      IW_STL_Hash_Map_int & tag_to_ndx,
                      Set_of_Bit_Mappings & stfbm)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return profile_fingerprints(input, tag_to_ndx, stfbm);
}


static int
sparse_to_fixed_size(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vF:bfp:P:ek:B:S:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('e'))
  {
    echo_changed_fingerprint = 1;

    if (verbose)
      cerr << "Will echo changed fingerprint\n";
  }

  if (cl.option_present('b'))
  {
    produce_fixed_binary = 1;
    produce_fixed_counted = 0;

    if (verbose)
      cerr << "Will produce fixed binary fingerprints\n";
  }

  if (cl.option_present('c'))
  {
    produce_fixed_counted = 1;
    produce_fixed_binary = 0;

    if (verbose)
      cerr << "Will produce fixed counted fingerprints\n";

    cerr << "Sorry, the -c option is not implemented, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
    return 3;
  }

  if (cl.option_present('B'))
  {
    if (! cl.value('B', user_specified_fingerprint_width) || user_specified_fingerprint_width < 8)
    {
      cerr << "The user specified fixed width fingerprint size (-B) must be a valid bit count\n";
      usage(1);
    }

    if (0 != user_specified_fingerprint_width % 8)
      user_specified_fingerprint_width = (user_specified_fingerprint_width / 8 + 1) * 8;

    if (verbose)
      cerr << "Will produce fixed width fingerprints of width " << user_specified_fingerprint_width << endl;
  }

  if (cl.option_present('f') || (1 == cl.number_elements() && 0 == strcmp("-", cl[0])))
  {
    work_as_filter = 1;

    if (verbose)
      cerr << "Will work as a filter\n";
  }

  if (cl.option_present('k'))   // important to do this near the start
  {
    const_IWSubstring k;
    for (int i = 0; cl.value('k', k, i); i++)
    {
      int v;
      if (k.numeric_value(v) && v > 0)
      {
        numeric_threshold.push_back(v);
        continue;
      }

      int j = 0;
      const_IWSubstring s;
      while (k.nextword(s, j, ','))
      {
        if (! s.numeric_value(v) || v <= 0)
        {
          cerr << "Invalid threshold specification '" << s << "'\n";
          return 3;
        }

        numeric_threshold.push_back(v);
      }
    }
  }
  else
    numeric_threshold.push_back(1);

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

// Generally, we process multiple tags. We need a hash to convert the tag to the specific
// Sparse_to_Fixed_Bit_Mapping object

  IW_STL_Hash_Map_int tag_to_ndx;

  Set_of_Bit_Mappings stfbm;

  if (cl.option_present('F'))
  {
    IWString f;
    for (int i = 0; cl.value('F', f, i); i++)
    {
      if (! f.ends_with('<'))
        f << '<';

      tag_to_ndx[f] = i;

      if (verbose)
        cerr << "non-colliding fingerprint '" << f << "' is fingerprint " << i << "\n";

      stfbm.add(new Sparse_to_Fixed_Bit_Mapping(f));
    }

    need_to_determine_fingerprints_present = 0;
  }

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! profile_fingerprints(cl[i], tag_to_ndx, stfbm))
    {
      cerr << "Cannot read input file '" << cl[i] << "\n";
      return i + 1;
    }
  }

  int nfp = tag_to_ndx.size();

  if (0 == nfp)
  {
    cerr << "No sparse fingerprints in input\n";
    return 2;
  }

  assert (tag_to_ndx.size() == stfbm.size());

  if (verbose)
    cerr << "Profiled " << nfp << " fingerprints across " << tdts_read << " fingerprints\n";

  if (cl.option_present('p'))
  {
    const_IWSubstring c;
    cl.value('p', c);

    Number_or_Percentage<int> cutoff;

    if (! cutoff.build(c))
    {
      cerr << "Invalid lower support level '" << c << "'\n";
      return 2;
    }

    less<unsigned int> intless;

    do_trim(nfp, stfbm, cutoff, intless);
  }

  if (cl.option_present('P'))
  {
    const_IWSubstring c;
    cl.value('P', c);

    Number_or_Percentage<int> cutoff;

    if (! cutoff.build(c))
    {
      cerr << "Invalid upper support level '" << c << "'\n";
      return 2;
    }

    greater<unsigned int> intgreater;

    do_trim(nfp, stfbm, cutoff, intgreater);
  }

  for (int i = 0; i < nfp; i++)
  {
    if (! stfbm[i]->setup_bit_mapping())
    {
      cerr << "Cannot initialise fingerprint " << i << endl;
      return 3;
    }
  }

  if (cl.option_present('S'))
  {
    IWString s = cl.string_value('S');
    for (int i = 0; i < nfp; ++i)
    {
      IWString fname = s;
      fname << stfbm[i]->fpname();
      stfbm[i]->write_bit_mapping(' ', fname.null_terminated_chars());
    }
  }

#ifdef WRITING_CROSS_REFERENCE_FILE_DEFINED
  if (cl.option_present('X'))
  {
    IWString x = cl.string_value('X');
  
    write_cross_reference_file(x, fp_to_index, fp_count);
  }
#endif

  IWString_and_File_Descriptor output(1);

  if (work_as_filter)
  {
    sparse_to_fixed_size(file_records, tag_to_ndx, stfbm, output);
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! sparse_to_fixed_size(cl[i], tag_to_ndx, stfbm, output))
      {
        cerr << "Fatal error processing '" << cl[i] << "'\n";
        return i + 1;
      }
    }
  }

  output.flush();

  return 0;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = sparse_to_fixed_size(argc, argv);

  return rc;
}
