/*
  Converts descriptor(s) to fingerprint form
  The main complexity is the specification of how to do the
  conversion

  dname:B<n>:R,min,max,dx
*/

#include <stdlib.h>
#include <limits>

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "iw_stl_hash_map.h"
#include "sparse_fp_creator.h"
#include "iw_auto_array.h"
#include "misc.h"

using std::numeric_limits;

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static IWString tag("NCFDSC<");

static IW_STL_Hash_Map_String smiles;

static int function_as_filter = 0;

static int default_nbits = 1;

static float default_number_divisions = 10.0f;

/*
  There are two ways of changing a floating point value into an integer.
  Both involve converting the float into a number in the range.
  But then, one uses just one bit, and the count is the converted value.
  The other way is to set different bits depending on which bucket is hit.

  IN the first scenario, floating point values of 3.1 and 3.5 will be
  similar to each other if there is a bucket divide at 3.3, but in 
  the second scenari, they would hit different bits and would 
  have nothing in common.
*/

static int range_sets_bit_count = 1;

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Converts descriptors to a sparse fingerprint\n";
  cerr << " -D ...         specify descriptor(s) to be converted;\n";
  cerr << "                descriptor name(s) found in descriptor file header line\n";
  cerr << " -T <tag>       tag for fingerprints\n";
  cerr << " -S <fname>     smiles file - to get smiles in output\n";
  cerr << " -f             function as a TDT filter\n";
  cerr << " -P <fname>     pre-computed descriptor file - needed with -f option\n";
  cerr << " -P <fname>     test set pre-computed descriptors (2nd -P option)\n";
  cerr << " -r <n>         default bit replicates\n";
  cerr << " -d <n>         default number divisions in range (10)\n";
  cerr << " -h             set different bits depending on bucket value\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

/*
  When reading a pre-computed descriptor file, we need an array of the descriptors
  being processed for every identifier
*/

typedef IW_STL_Hash_Map<IWString, float *>  Descriptor_Values ;

class Descriptor
{
  private:
    IWString _name;
    int _column;

    int _nbits;

    float _min, _max, _dx;

    unsigned int _bstart;

    int * _bucket_count;

//  We need some means of identifying those descriptors that have had their
//  ranges specified, and those for which we need to scan the input

    int _range_specified;

//  private functions

    int _parse_bit_replicats(const_IWSubstring s);
    int _parse_range_specification(const_IWSubstring s);

  public:
    Descriptor();
    ~Descriptor();

    int debug_print(std::ostream &) const;
    int report     (std::ostream &) const;

    const IWString & name() const { return _name;}

    int allocate_bucket_counter();

    void set_name(const const_IWSubstring & s) { _name = s;}
    void set_column(int s) { _column = s;}
    int  column() const { return _column;}

    int build(const const_IWSubstring &);

    int nbits() const { return _nbits;}

    void set_bstart(unsigned int s) { _bstart = s;}

    int set_bits(const const_IWSubstring & buffer, Sparse_Fingerprint_Creator & sfc);

    int set_bits(float v, Sparse_Fingerprint_Creator & sfc);

    int range_specified() const { return _range_specified;}

    void scan_input_values(float);
    int  compute_dx();
};

Descriptor::Descriptor()
{
  _column = -1;

  _nbits = default_nbits;

  _bstart = 0;

  _min = numeric_limits<float>::max();
  _max = -numeric_limits<float>::max();
  _dx = 0.0f;

  _bucket_count = NULL;

  _range_specified = 0;

  return;
}

Descriptor::~Descriptor()
{
  if (NULL != _bucket_count)
    delete [] _bucket_count;

  return;
}

int
Descriptor::debug_print(std::ostream & os) const
{
  os << "Descriptor '" <<_name << "'\n";
  if (_column > 0)
    os << " column " << _column;
  os << " range " << _min << ',' << _max << ',' << _dx;
  os << ", " << _nbits << " replicates.";

  os << " bstart " << _bstart << '\n';

  return 1;
}

int
Descriptor::report(std::ostream & os) const
{
  os << "Descriptor '" << _name << "' range " << _min << ',' << _max << ',' << _dx << endl;
  if (NULL == _bucket_count)
  {
    cerr << "constant\n";
    return 1;
  }

  os << "Bucket occupancy\n";

  int nb = static_cast<int>((_max - _min) / _dx + 0.49999F);

  cerr << "Reporting for " << nb << " buckets\n";

  for (int i = 0; i <= nb; i++)
  {
    os << " b = " << i << " " << _bucket_count[i] << '\n';
  }

  return 1;
}

int
Descriptor::allocate_bucket_counter()
{
  assert (NULL == _bucket_count);

  if (_min == _max)
    return 0;

  int nb = static_cast<int>((_max - _min) / _dx + 0.4999F);

  if (nb <= 0)
  {
    cerr << "Descriptor::allocate_bucket_counter:invalid bucket count, min " << _min << " max " << _max << " dx " << _dx << " nb " << nb << endl;
    return 0;
  }

  _bucket_count = new_int(nb+1);

  return 1;
}

int
Descriptor::_parse_bit_replicats(const_IWSubstring s)    // note pass by value
{
  assert ('B' == s[0]);
  
  s.remove_leading_chars(1);

  if (0 == s.length())
  {
    cerr << "Descriptor::_parse_bit_replicats:zero length bit replicate specification\n";
    return 0;
  }

  if (! s.numeric_value(_nbits) || _nbits <= 0)
  {
    cerr << "Descriptor::_parse_bit_replicats:invalid bit replicate specification\n";
    return 0;
  }

  return 1;
}

int
Descriptor::_parse_range_specification(const_IWSubstring s)    // note local copy
{
  assert ('R' == s[0]);

  s.remove_leading_chars(1);

  if (2 != s.ccount(','))
  {
    cerr << "Descriptor::_parse_range_specification:range must have three components\n";
    return 0;
  }

  const_IWSubstring token;
  int i = 0;

  s.nextword(token, i, ',');

  if (0 == token.length())
  {
    cerr << "Descriptor::_parse_range_specification:zero length min\n";
    return 0;
  }

  if (! token.numeric_value(_min))
  {
    cerr << "Descriptor::_parse_range_specification:invalid min\n";
    return 0;
  }

  s.nextword(token, i, ',');

  if (0 == token.length())
  {
    cerr << "Descriptor::_parse_range_specification:zero length max\n";
    return 0;
  }

  if (! token.numeric_value(_max))
  {
    cerr << "Descriptor::_parse_range_specification:invalid max\n";
    return 0;
  }

  if (_max <= _min)
  {
    cerr << "Descriptor::_parse_range_specification:max inconsistent with min\n";
    return 0;
  }

  s.nextword(token, i, ',');

  if (0 == token.length())
  {
    cerr << "Descriptor::_parse_range_specification:zero length dx\n";
    return 0;
  }

  if (! token.numeric_value(_dx) || _dx <= 0.0)
  {
    cerr << "Descriptor::_parse_range_specification:invalid dx\n";
    return 0;
  }

  if (_min + _dx >= _max)
  {
    cerr << "Descriptor::_parse_range_specification:dx too large\n";
    return 0;
  }

  return 1;
}

/*
  This is kind of messy.
  
  The most complete specification would look like

    dname:min,max,dx,replicates

  The replicates token is optional

  But just 'dname' is a valid specification.
*/

int
Descriptor::build(const const_IWSubstring & s)
{
  if (1 == s.ccount(':') && s.ccount(',') >= 2)
    ;
  else if (0 == s.ccount(':') && 0 == s.ccount('<'))
  {
    _name = s;
    return 1;
  }
  else
  {
    cerr << "Descriptor::build:descriptor bit specifications must contain ':' and ',' separators\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  if (! s.nextword(_name, i, ':'))
  {
    cerr << "Descriptor::build:empty string\n";
    return 0;
  }

  i++;

  for (int ndx = 0; s.nextword(token, i, ','); ndx++)
  {
    if (0 == ndx)
    {
      if (! token.numeric_value(_min))
      {
        cerr << "Descriptor::build:invalid min '" << token << "'\n";
        return 0;
      }
    }
    else if (1 == ndx)
    {
      if (! token.numeric_value(_max))
      {
        cerr << "Descriptor::build:invalid max '" << token << "'\n";
        return 0;
      }

      if (_min >= _max)
      {
        cerr << "Descriptor::build:max is invalid, min = " << _min << " max = " << _max << endl;
        return 0;
      }
    }
    else if (2 == ndx)
    {
      if (! token.numeric_value(_dx) || _dx < 0.0f)
      {
        cerr << "Descriptor::build:invalid dx '" << token << "'\n";
        return 0;
      }

      if (_min + _dx > _max)
      {
        cerr << "Descriptor::build:inconsistent dx, min = " << _min << ", _max = " << _max << " dx = " << _dx << ", impossible\n";
        return 0;
      }
    }
    else if (3 == ndx)
    {
      if (! token.numeric_value(_nbits) || _nbits < 1)
      {
        cerr << "Descriptor::build:invalid nbits specification '" << token << "'\n";
        return 0;
      }
    }
    else
    {
      cerr << "Descriptor::build:too many tokens '" << token << "'\n";
      return 0;
    }
  }

  _range_specified = 1;

  return 1;
}

int
Descriptor::set_bits(float v,
                      Sparse_Fingerprint_Creator & sfc)
{
  int b;
  if (v <= _min)
    b = 0;
  else
  {
    if (v >= _max)
      v = _max;

    b = static_cast<int>((v - _min) / _dx + 0.4999F);
  }

  if (NULL != _bucket_count)
    _bucket_count[b]++;

  if (range_sets_bit_count)
  {
    b++;    // cannot have any zero counts
    for (int i = 0; i < _nbits; i++)
    {
      sfc.hit_bit(_bstart + i, b);
    }
  }
  else
    sfc.hit_bit(_bstart + b, _nbits);

  return 1;
}

int
Descriptor::set_bits(const const_IWSubstring & buffer,
                      Sparse_Fingerprint_Creator & sfc)
{
  const_IWSubstring token;

  if (! buffer.word(_column, token))
  {
    cerr << "Descriptor::set_bits:cannot extract column " << _column << " for '" << _name << "'\n";
    return 0;
  }

  float v;

  if (token.numeric_value(v))
    ;
  else if ('.' == token)
    return 1;
  else
  {
    cerr << "Descriptor::set_bits:invalid floating point value '" << token << "'\n";
    return 0;
  }

//cerr << "Descriptor processing '" << buffer << "' got value " << v << endl;

  return set_bits(v, sfc);
}

void
Descriptor::scan_input_values(float v)
{
  if (v < _min)
    _min = v;

  if (v > _max)
    _max = v;

  return;
}

int
Descriptor::compute_dx()
{
  assert (_min <= _max);

  _dx = (_max - _min) / default_number_divisions;

  return 1;
}

class Set_of_Descriptors
{
  private:
    int _nd;
    Descriptor * _d;

//  private functions

    void _set_bstart();
    int _initialise_all_descriptors(const const_IWSubstring & buffer);
    int _initialise_already_specified(const const_IWSubstring & buffer);

  public:
    Set_of_Descriptors();
    ~Set_of_Descriptors();

    int debug_print (std::ostream &) const;
    int report (std::ostream & os) const;

    int initialise (Command_Line &, char, int);
    int initialise (const const_IWSubstring &);

    int already_initialised() const;

    int number_descriptors () const { return _nd;}

    int identify_descriptors_being_processed (const const_IWSubstring & buffer,
                                      IWString & id,
                                      float * dvalues) const;
    int initialise_any_descriptors_not_already_specified (const Descriptor_Values & descriptor_values);

    void allocate_bucket_counter();

    int all_descriptor_ranges_specified () const;

    void set_bits (const const_IWSubstring &, Sparse_Fingerprint_Creator & sfc);
    void set_bits (const float *, Sparse_Fingerprint_Creator & sfc);
};

Set_of_Descriptors::Set_of_Descriptors()
{
  _nd = 0;
  _d = NULL;

  return;
}

Set_of_Descriptors::~Set_of_Descriptors()
{
  if (NULL != _d)
    delete [] _d;

  return;
}

int
Set_of_Descriptors::already_initialised() const
{
  if (0 == _nd)
    return 0;

  for (int i = 0; i < _nd; i++)
  {
    if (_d[i].column() < 0)
      return 0;
  }

  return 1;
}

int
Set_of_Descriptors::all_descriptor_ranges_specified() const
{
  for (auto i = 0; i < _nd; ++i)
  {
    if (! _d[i].range_specified())
      return 0;
  }

  return 1;
}

int
Set_of_Descriptors::initialise(Command_Line & cl,
                                char flag,
                                int verbose)
{
  _nd = cl.option_count(flag);

  if (0 == _nd)
  {
    cerr << "Set_of_Descriptors::initialise:no descriptors -" << flag << " specified\n";
    return 0;
  }

  _d = new Descriptor[_nd];

  for (int i = 0; i < _nd; i++)
  {
    const_IWSubstring s = cl.string_value(flag, i);

    if (! _d[i].build(s))
    {
      cerr << "Set_of_Descriptors::initialise:invalid descriptor specification '" << s << "'\n";
      return 0;
    }
  }

  return _nd;
}

int
Set_of_Descriptors::debug_print(std::ostream & os) const
{
  os << "Set_of_Descriptors with " << _nd << " descriptors\n";
  for (int i = 0; i < _nd; i++)
  {
    _d[i].debug_print(os);
  }

  return 1;
}

int
Set_of_Descriptors::report(std::ostream & os) const
{
  os << "Set_of_Descriptors:report on " << _nd << " descriptors\n";

  for (int i = 0; i < _nd; i++)
  {
    _d[i].report(os);
  }

  return 1;
}

static int
descriptors_to_fingerprint(IWString buffer,    // note local copy
                            const Descriptor_Values & descriptor_values,
                            Set_of_Descriptors & sod,
                            IWString_and_File_Descriptor & output)
{
  if (! buffer.ends_with('>'))
  {
    cerr << "Invalid TDT form\n";
    return 0;
  }

  buffer.chop();

  buffer.remove_leading_chars(identifier_tag.length());

  buffer.truncate_at_first(' ');

//cerr << "Identifier is '" << buffer << "'\n";

  Descriptor_Values::const_iterator f = descriptor_values.find(buffer);
//auto f = descriptor_values.find(buffer);

  if (f == descriptor_values.end())
  {
    cerr << "Yipes, no descriptors for '" << buffer << "'\n";
    return 0;
  }

  molecules_read++;

  const float * v = (*f).second;

  Sparse_Fingerprint_Creator sfc;
  
  sod.set_bits(v, sfc);

  IWString tmp;

  sfc.daylight_ascii_form_with_counts_encoded(tag, tmp);

  output << tmp << '\n';

  return 1;
}

static int
descriptors_to_fingerprint (iwstring_data_source & input,
                            const Descriptor_Values & descriptor_values,
                            Set_of_Descriptors & sod,
                            IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(32768);

    if (! buffer.starts_with(identifier_tag))
      continue;

    if (! descriptors_to_fingerprint (buffer, descriptor_values, sod, output))
    {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

void
Set_of_Descriptors::set_bits (const const_IWSubstring & buffer,
                              Sparse_Fingerprint_Creator & sfc)
{
  for (int i = 0; i < _nd; i++)
  {
    _d[i].set_bits (buffer, sfc);
  }

  return;
}

void
Set_of_Descriptors::set_bits (const float * v,
                              Sparse_Fingerprint_Creator & sfc)
{
  for (int i = 0; i < _nd; i++)
  {
    _d[i].set_bits(v[i], sfc);
  }

  return;
}

static int
descriptors_to_fingerprint_record (const_IWSubstring & buffer,
                            Set_of_Descriptors & sod,
                            IWString_and_File_Descriptor & output)
{
  Sparse_Fingerprint_Creator sfc;
  
  sod.set_bits(buffer, sfc);

  IWString id;
  buffer.word(0, id);

  if (smiles.size())
  {
    IW_STL_Hash_Map_String::const_iterator f = smiles.find(id);
//  const auto f = smiles.find(id);

    if (f == smiles.end())
    {
      cerr << "Yipes, no smiles for '" << id << "'\n";
      return 0;
    }

    output << smiles_tag << (*f).second << ">\n";
  }

  output << identifier_tag << id << ">\n";

  IWString tmp;

  sfc.daylight_ascii_form_with_counts_encoded(tag, tmp);

  output << tmp << '\n';

  output << "|\n";

  return 1;
}

int
Set_of_Descriptors::initialise (const const_IWSubstring & buffer)
{
  if (_nd > 0)
    return _initialise_already_specified(buffer);
  else
    return _initialise_all_descriptors(buffer);
}

int
Set_of_Descriptors::_initialise_already_specified (const const_IWSubstring & buffer)
{
  assert(_nd > 0);

  int i = 0;
  const_IWSubstring token;

  int nfound = 0;

  for (int col = 0; nfound < _nd ; col++)
  {
    if (! buffer.nextword(token, i))
    {
      cerr << "Did not find one or more descriptors, col = " << col << " nfound " << nfound << "\n";

      for (int j = 0; j < _nd; j++)
      {
        if (_d[j].column() <= 0)
          cerr << "No match for '" << _d[j].name() << endl;
      }

      return 0;
    }

    if (0 == col)
      continue;

    for (int j = 0; j < _nd; j++)
    {
      if (_d[j].name() != token)
        continue;

      _d[j].set_column(col);
      nfound++;
      break;
    }
  }

  if (nfound != _nd)
  {
    cerr << "Set_of_Descriptors::initialise::problems initialising columns, have " << _nd << " descriptors, matched " << nfound <<endl;
    return 0;
  }

  _set_bstart();

  return nfound;
}

void
Set_of_Descriptors::_set_bstart ()
{
  unsigned int bstart = 0;

  for (int i = 0; i < _nd; i++)
  {
    _d[i].set_bstart(bstart);

    bstart += _d[i].nbits();
  }

  return;
}

int
Set_of_Descriptors::_initialise_all_descriptors (const const_IWSubstring & buffer)
{
  assert (0 == _nd);

  int nw = buffer.nwords();

  if (nw < 2)
  {
    cerr << "Set_of_Descriptors::_initialise_all_descriptors:too few tokens in header\n";
    return 0;
  }

  _nd = nw - 1;

  _d = new Descriptor[_nd];

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);   // skip over first column

  for (int col = 0; buffer.nextword(token, i); col++)
  {
    if (! _d[col].build(token))
      return 0;

    _d[col].set_column(col + 1);
  }

  _set_bstart();

  return _nd;
}

static int
descriptors_to_fingerprint (iwstring_data_source & input,
                            Set_of_Descriptors & sod,
                            IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read header record\n";
    return 0;
  }

  if (! sod.initialise (buffer))
  {
    cerr << "Cannot initialise descriptors based on header\n";
    return 0;
  }

  while (input.next_record (buffer))
  {
    molecules_read++;

//  cerr << "Read " << molecules_read << "'th molecule\n";

    if (! descriptors_to_fingerprint_record (buffer, sod, output))
    {
      cerr << "Fatal error processing line " << input.lines_read() << endl;
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

int
Set_of_Descriptors::identify_descriptors_being_processed (const const_IWSubstring & buffer,
                                      IWString & id,
                                      float * dvalues) const
{
  int i = 0;

  if (! buffer.nextword(id, i))
  {
    cerr << "Empty record, cannot extract identifier\n";
    return 0;
  }

  const_IWSubstring token;
  int nfound = 0;

  for (int col = 1; buffer.nextword(token, i); col++)
  {
    for (int j = 0; j < _nd; j++)
    {
      if (col != _d[j].column())
        continue;

      if (! token.numeric_value(dvalues[j]))
      {
        cerr << "Invalid numeric '" << token << "'\n";
        return 0;
      }

      nfound++;

      if (nfound >= _nd)
        return 1;
    }
  }

  cerr << "Yipes, only identified " << nfound << " of " << _nd << " descriptors in record\n";
  cerr << buffer << endl;

  for (int i = 0; i < _nd; i++)
  {
    cerr << " descriptor in column " << _d[i].column() << endl;
  }

  return 0;
}

static int
read_descriptors (iwstring_data_source & input,
                  Descriptor_Values & descriptor_values,
                  Set_of_Descriptors & sod)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read header record from precomputed descriptor file\n";
    return 0;
  }

  if (sod.already_initialised())
    ;
  else if (! sod.initialise (buffer))
  {
    cerr << "Cannot initialise descriptors based on header\n";
    return 0;
  }

  while (input.next_record(buffer))
  {
    float * tmp = new float[sod.number_descriptors()];    // never deallocated, OK
    IWString id;

    if (! sod.identify_descriptors_being_processed(buffer, id, tmp))
    {
      cerr << "Cannot process '" << buffer << "'\n";
      return 0;
    }

    if (descriptor_values.contains(id))
    {
      cerr << "Yipes, duplicate identifier '" << id << "' in descriptor file\n";
      return 0;
    }

    descriptor_values[id] = tmp;
  }

  return descriptor_values.size();
}

static int
read_descriptors (const char * fname, 
                  Descriptor_Values & descriptor_values,
                  Set_of_Descriptors & sod)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open file for precomputed descriptors '" << fname << "'\n";
    return 0;
  }

  return read_descriptors(input, descriptor_values, sod);
}

static int
read_smiles_record (const const_IWSubstring & buffer,
                    IW_STL_Hash_Map_String & smiles)
{
  IWString s, id;

  int i = 0;

  if (! buffer.nextword(s, i) || ! buffer.nextword(id, i))
  {
    cerr << "Cannot extract smiles and id\n";
    return 0;
  }

  smiles[id] = s;

  return smiles.size();
}

static int
read_smiles (iwstring_data_source & input,
             IW_STL_Hash_Map_String & smiles)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (! read_smiles_record (buffer, smiles))
    {
      cerr << "Cannot read smiles record '" << buffer << "'\n";
      return 0;
    }
  }

  return smiles.size();
}

static int
read_smiles (const char * fname,
             IW_STL_Hash_Map_String & smiles)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open smiles file '" << fname << "'\n";
    return 0;
  }

  return read_smiles (input, smiles);
}

static int
descriptors_to_fingerprint (const char * fname,
                            Set_of_Descriptors & sod,
                            IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return descriptors_to_fingerprint(input, sod, output);
}

/*
  Note that this seems not to work if the user specifies some descriptors on
  the command line, because those would get entered into the _d array out of
  order. This is all kind of messy, should undergo some fixing sometime
*/

int
Set_of_Descriptors::initialise_any_descriptors_not_already_specified (const Descriptor_Values & descriptor_values)
{
  int need_initialisation = 0;

  for (int i = 0; i < _nd; i++)
  {
    if (! _d[i].range_specified())
    {
      need_initialisation = 1;
      break;
    }
  }

  if (0 == need_initialisation)
    return 1;

  for (Descriptor_Values::const_iterator i = descriptor_values.begin(); i != descriptor_values.end(); ++i)
//for (auto i = descriptor_values.begin(); i != descriptor_values.end(); ++i)
  {
    const float * v = (*i).second;

    for (int j = 0; j < _nd; j++)
    {
      if (! _d[j].range_specified())
        _d[j].scan_input_values(v[j]);
    }
  }

  for (int i = 0; i < _nd; i++)
  {
    if (! _d[i].range_specified())
      _d[i].compute_dx();
  }

  return 1;
}

void
Set_of_Descriptors::allocate_bucket_counter()
{
  for (int i = 0; i < _nd; i++)
  {
    _d[i].allocate_bucket_counter();
  }

  return;
}

static int
descriptors_to_fingerprint (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vT:D:S:fP:r:d:h");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('T'))
  {
    cl.value('T', tag);

    if (verbose)
      cerr << "Fingerprints written with tag '" << tag << "'\n";

    if (! tag.ends_with('<'))
      tag << '<';
  }

  if (cl.option_present('r'))   // must do before instantiating any Descriptor objects
  {
    if (! cl.value('r', default_nbits) || default_nbits < 1)
    {
      cerr << "The default bit replicates option (-r) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Default bit replicate count set to " << default_nbits << endl;
  }

  if (cl.option_present('d'))
  {
    int tmp;
    if (! cl.value('r', tmp) || tmp < 1)
    {
      cerr << "The default number divisions option (-d) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Default range division count set to " << tmp << endl;

    default_number_divisions = static_cast<float>(tmp);
  }

  if (cl.option_present('h'))
  {
    range_sets_bit_count = 0;

    if (verbose)
      cerr << "The bits set will be a function of bucket occupancy\n";
  }

  Set_of_Descriptors sod;

  if (cl.option_present('D'))
  {
    if (! sod.initialise(cl, 'D', verbose))
    {
      cerr << "Cannot initialise descriptors\n";
      return 2;
    }
  }

  if (cl.option_present('S'))
  {
    const char * s = cl.option_value('S');

    if (! read_smiles (s, smiles))
    {
      cerr << "Cannot read smiles from '" << s << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Read " << smiles.size() << " smiles from '" << s << "'\n";
  }

  Descriptor_Values descriptor_values;

// TODO. If all information about the descriptors has been specified, then we do not really need the -P file(s).

  if (cl.option_present('f'))
  {
    function_as_filter = 1;

    int np = cl.option_count('P');

    if (0 == np)
    {
      cerr << "When working as a TDT filter, must specify precomputed descriptor file via the -P option\n";
      usage(2);
    }

    const char * p = cl.option_value('P');

    if (! read_descriptors(p, descriptor_values, sod))
    {
      cerr << "Cannot read precomputed descriptors from '" << p << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Read " << descriptor_values.size() << " precomputed descriptor values from '" << p << "'\n";

    sod.initialise_any_descriptors_not_already_specified(descriptor_values);

    if (np > 1)
    {
      for (int i = 1; i < cl.option_count('P'); i++)
      {
        p = cl.option_value('P', i);

        if (! read_descriptors (p, descriptor_values, sod))
        {
          cerr << "Cannot read descriptors from '" << p << "'\n";
          return 3;
        }
      }

      if (np > 1)
        cerr << "Read " << descriptor_values.size() << " precomputed descriptor values from " << np <<" files\n";
    }
  }

  if (verbose)
  {
    sod.allocate_bucket_counter();
    sod.debug_print(cerr);
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;

  if (function_as_filter)
  {
    iwstring_data_source input(0);

    if (! descriptors_to_fingerprint(input, descriptor_values, sod, output))
    {
      cerr << "Fatal error processing precomputed descriptors in filter\n";
      rc = 3;
    }
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! descriptors_to_fingerprint(cl[i], sod, output))
      {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    if (molecules_read > 0)
      sod.report(cerr);
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = descriptors_to_fingerprint(argc, argv);

  return rc;
}
