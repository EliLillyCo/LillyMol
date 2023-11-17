/*
  We have run some kind of scoring thing which has associated each atom with one or more
  contributions. We want to produce a coloured molecule that reflects those weights
  Builds an input file for cactvs.
*/

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/mdl.h"

using std::ostream;
using std::numeric_limits;
using std::unique_ptr;
using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int remove_all_chiral_centres = 0;

static int min_radius = 0;
static int max_radius = 3;

static double min_score = 0.0;
static double max_score = 0.0;

static int linear_between_min_and_max = 0;

static IWString_and_File_Descriptor stream_for_isotopically_labelled;

static IWString_and_File_Descriptor stream_for_atomic_contributions;

/*
  We can write a coloured structure for each class, or we can
  aggregate the values across all classes

static int per_class_contribution = 1;
*/

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Produces input for csib based on bit related atomic contributions.\n";
  cerr << "For example the -O file produced by gfp_naive_bayesian is a possible input\n";
  cerr << "  -C <fname>    set of colours - #FF4800FF, one per line (from R's rainbow function)\n";
  cerr << "  -O <fname>    bit contribution file - probably from the -O option of gfp_naive_bayesian\n";
  cerr << "  -B <fname>    bit description  file - probably from the -B option of iwecfp\n";
  cerr << "  -r <rad>      minimum radius\n";
  cerr << "  -R <rad>      maximum radius\n";
  cerr << "  -x <min>      truncate scores below <min>\n";
  cerr << "  -X <max>      truncate scores above <max>\n";
  cerr << "  -I <fname>    write isotopically labelled molecules to <fname>\n";
  cerr << "  -D <dcf>      directcolorfile - produced, feed to csib\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -c            remove all chirality on input\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
preprocess(Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (remove_all_chiral_centres)
    m.remove_all_chiral_centres();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

/*
  Class designed to hold a record from an Ofile

  Bit       Class 1     Class 2

  709623808 -0.04330375 0.2773229
*/

class Bit_and_Class_Weights
{
  private:
    unsigned int _bit;
    const int _nclasses;
    float * _class_weight;

  public:
    Bit_and_Class_Weights(int nc);
    ~Bit_and_Class_Weights();

    unsigned int bit () const { return _bit;}

    int build (const const_IWSubstring & buffer);

    void copy_weights(double *) const;
};

Bit_and_Class_Weights::Bit_and_Class_Weights(int nc) : _nclasses (nc)
{
  _bit = 0;
  _class_weight = new float[_nclasses];

  return;
}

Bit_and_Class_Weights::~Bit_and_Class_Weights ()
{
  if (nullptr != _class_weight)
    delete [] _class_weight;

  return;
}

int
Bit_and_Class_Weights::build(const const_IWSubstring & buffer)
{
  const auto nw = buffer.nwords();

  if (nw < 3)
  {
    cerr << "Bit_and_Class_Weights::build:input buffer must contain at least 3 tokens '" << buffer << "' invalid\n";
    return 0;
  }

  if (_nclasses != nw - 1)
  {
    cerr << "Bit_and_Class_Weights::build:expected " << _nclasses << " classes, so input '" << buffer << "' is incorrect\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);
  if (! token.numeric_value(_bit))
  {
    cerr << "Bit_and_Class_Weights::build:invalid bit specification '" << buffer << "'\n";
    return 0;
  }

  for (auto j = 0; j < _nclasses; ++j)
  {
    buffer.nextword(token, i);
    if (! token.numeric_value(_class_weight[j]))
    {
      cerr << "Bit_and_Class_Weights::build:invalid class weight '" << buffer << "'\n";
      return 0;
    }
  }

  return _nclasses;
}

void
Bit_and_Class_Weights::copy_weights(double * c) const
{
  std::copy_n(_class_weight, _nclasses, c);

  return;
}

class Set_of_Bit_and_Class_Weights : public resizable_array_p<Bit_and_Class_Weights>
{
  private:
  public:
    int build(const int nclasses, iwstring_data_source & input);

    int get_contribution(const unsigned int b, double *) const;
};

int
Set_of_Bit_and_Class_Weights::build(const int nclasses,
                                     iwstring_data_source & input)
{
//cerr << "Set_of_Bit_and_Class_Weights::build:expecting " << nclasses << " classes\n";

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if ('|' == buffer)
      return _number_elements;

    Bit_and_Class_Weights * bcw = new Bit_and_Class_Weights(nclasses);

    if (! bcw->build(buffer))
    {
      cerr << "Set_of_Bit_and_Class_Weights::build:invalid input '" << buffer << "'\n";
      delete bcw;
      return 0;
    }

    add(bcw);
  }

  return _number_elements;
}

int
Set_of_Bit_and_Class_Weights::get_contribution(const unsigned int b,
                                                double * c) const
{
//#define DEBUG_GET_CONTRIBUTION
#ifdef DEBUG_GET_CONTRIBUTION
  cerr << "Set_of_Bit_and_Class_Weights::get_contribution:looking for bit " << b << endl;
  int found_match = 0;
  for (auto i = 0; i < _number_elements; ++i)
  {
    if (_things[i]->bit() == b)
    {
      cerr << "Should match at ndx " << i << endl;
      found_match++;
    }
  }
  if (! found_match)
    cerr << "Set_of_Bit_and_Class_Weights::get_contribution:no matches anticipated\n";
#endif

  const auto f = std::lower_bound(this->cbegin(), this->cend(), b, [] (const Bit_and_Class_Weights * bacw, const unsigned int b) { return b > bacw->bit();});

  if (f == this->cend())    // should not happen
    return 0;

  (*f)->copy_weights(c);

  return 1;
}

class Ofile_Contents
{
  private:
    IW_STL_Hash_Map_off_t _id_to_offset;

    int _nclasses;
    IWString * _class_name;

    IWString _molecule_name;

    resizable_array_p<Set_of_Bit_and_Class_Weights> _fixed_fingerprints;
    resizable_array_p<Set_of_Bit_and_Class_Weights> _sparse_fingerprints;

    float * _score;

    iwstring_data_source _input;

    double * _global_min;
    double * _global_max;

//  private functions

    int _parse_score_record(const const_IWSubstring & buffer);
    int _another_fingerprint(const_IWSubstring buffer, const int nchars,
                              resizable_array_p<Set_of_Bit_and_Class_Weights> & fp);

    int _read_header_record();
    int _determine_min_max_contributions();
    int _determine_min_max_contributions2();
    int _check_global_min_max(const const_IWSubstring & buffer);

  public:
    Ofile_Contents();
    ~Ofile_Contents();

    int nclasses() const { return _nclasses;}
    const IWString & class_name(int ndx) const { return _class_name[ndx];}

    int open(const char * fname);

    const IWString & molecule_name() const { return _molecule_name;}

    int get_next_molecule();
    int get_data_for(const IWString & mname);

    int debug_print(ostream &) const;

    int number_fingerprints() const { return _fixed_fingerprints.size() + _sparse_fingerprints.size();}

    int get_contribution(const unsigned int b, double * c) const;

    float score(int ndx) const { return _score[ndx];}
};

Ofile_Contents::Ofile_Contents()
{
  _nclasses = 0;
  _class_name = nullptr;
  _score = nullptr;

  _global_min = nullptr;
  _global_max = nullptr;

  return;
}

Ofile_Contents::~Ofile_Contents()
{
  if (nullptr != _class_name)
    delete [] _class_name;

  if (nullptr != _score)
    delete [] _score;

  if (nullptr != _global_min)
    delete [] _global_min;

  if (nullptr != _global_max)
    delete [] _global_max;

  return;
}

int
Ofile_Contents::open(const char * fname)
{
  if (! _input.open(fname))
  {
    cerr << "Ofile_Contents::open:cannot open '" << fname << "'\n";
    return 0;
  }

  _input.set_skip_blank_lines(1);

  if (! _read_header_record())
    return 0;

  return _determine_min_max_contributions();
}

/*
  We need to know the range of bit contributions in our file
*/

int
Ofile_Contents::_determine_min_max_contributions()
{
  _global_min = new double[_nclasses];
  _global_max = new double[_nclasses];

  if (0.0 != min_score && 0.0 != max_score)
  {
    for (auto i = 0; i < _nclasses; ++i)
    {
      _global_min[i] = min_score;
      _global_max[i] = max_score;
    }

    return 1;
  }

// kind of a bug here. If just one of min_score and max_score were set, we ignore them

  const auto o = _input.tellg();

  if (! _input.seekg(o))
  {
    cerr << "Ofile_Contents::_determine_min_max_contributions:yipes, cannot seek back to " << o << endl;
    return 0;
  }

  return 1;
}

int
Ofile_Contents::_determine_min_max_contributions2()
{
  std::fill_n(_global_min, _nclasses, numeric_limits<float>::max());
  std::fill_n(_global_max, _nclasses, numeric_limits<float>::min());

  const_IWSubstring buffer;


  for (auto o = _input.tellg(); _input.next_record(buffer); o = _input.tellg())
  {
    if (buffer.starts_with("ID:"))
    {
      buffer.remove_leading_chars(3);
      buffer.strip_leading_blanks();
      buffer.strip_trailing_blanks();
      _id_to_offset[buffer] = o;
      continue;
    }

    if (buffer.starts_with("ID:") || buffer.starts_with("SFP") || buffer.starts_with("FP") || buffer.starts_with("Score"))
      continue;

    if (_nclasses + 1 != buffer.nwords())
      continue;

    if (! _check_global_min_max (buffer))
    {
      cerr << "Ofile_Contents::_determine_min_max_contributions2:invalid input '" << buffer << "'\n";
      return 0;
    }
  }

  if (verbose)
  {
    for (auto c = 0; c < _nclasses; ++c)
    {
      cerr << "Global range of bit contributions: class " << c << " << min << " << _global_min[c] << " max " << _global_max[c] << endl;
    }
  }

  if (0.0 == max_score && 0.0 == min_score)    // then we need to set them
  {
    min_score = _global_min[0];
    max_score = _global_max[0];

    for (auto i = 0; i < _nclasses; ++i)
    {
      if (_global_min[i] < min_score)
        min_score = _global_min[i];
      if (_global_max[i] > max_score)
        max_score = _global_max[i];
    }
  }

  return 1;
}

int
Ofile_Contents::_check_global_min_max (const const_IWSubstring & buffer)
{
  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);

  for (auto c = 0; buffer.nextword(token, i); ++c)
  {
    float v;
    if (! token.numeric_value(v))
      return 0;

    if (v < _global_min[c])
      _global_min[c] = v;

    if (v > _global_max[c])
      _global_max[c] = v;
  }

  return 1;
}

int
Ofile_Contents::debug_print (ostream & output) const
{
  if (0 == _molecule_name.length())
  {
    output << "Ofile_Contents::debug_print:no info\n";
    return 0;
  }

  output << "Ofile details for '" << _molecule_name << "'\n";
  output << _fixed_fingerprints.size() << " fixed fingerprints and " << _sparse_fingerprints.size() << " sparse fingerprints\n";

  return 1;
}

int
Ofile_Contents::_read_header_record ()
{
  const_IWSubstring buffer;

  if (! _input.next_record(buffer))
  {
    cerr << "Bfile_Contents::read_header_record:cannot read header record\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  if (! buffer.nextword(token, i) || '#' != token || ! buffer.nextword(token, i) || "Classes" != token)
  {
    cerr << "Ofile_Contents::read_header_record:invalid header '" << buffer << "'\n";
    return 0;
  }

  const auto nw = buffer.nwords();

  if (nw < 4)
  {
    cerr << "Ofile_Contents:::read_header_record:invalid header record, must be at least four columns '" << buffer << "'\n";
    return 0;
  }

  _nclasses = nw - 2;

  _class_name = new IWString[_nclasses];

  for (auto ndx = 0; buffer.nextword(token, i); ndx++)
  {
    _class_name[ndx] = token;
  }

  return 1;
}

int
Ofile_Contents::get_next_molecule ()
{
  _molecule_name.resize_keep_storage (0);

  _fixed_fingerprints.resize_keep_storage(0);
  _sparse_fingerprints.resize_keep_storage(0);

  const_IWSubstring buffer;

  if (! _input.next_record(buffer))
  {
    cerr << "Ofile_Contents::get_next_molecule:EOF\n";
    return 0;
  }

  if (! buffer.starts_with("ID: "))
  {
    cerr << "Ofile_Contents::get_next_molecule:invalid identifier record '" << buffer << "'\n";
    return 0;
  }

  _molecule_name = buffer;
  _molecule_name.remove_leading_chars(4);

  while (_input.next_record(buffer))
  {
    if (buffer.starts_with("Score"))
      return _parse_score_record(buffer);

    if (buffer.starts_with("FP"))
    {
      if (! _another_fingerprint(buffer, 2, _fixed_fingerprints))
        return 0;
    }
    else if (buffer.starts_with("SFP"))
    {
      if (! _another_fingerprint(buffer, 3, _sparse_fingerprints))
        return 0;
    }
    else
    {
      cerr << "Ofile_Contents::get_next_molecule:unrecognised input '" << buffer << "'\n";
      return 0;
    }
  }

  cerr << "Ofile_Contents::get_next_molecule:premature EOF\n";
  return 0;
}

int
Ofile_Contents::get_data_for (const IWString & mname)
{
  const auto f = _id_to_offset.find(mname);

  if (f == _id_to_offset.end())
  {
    cerr << "Ofile_Contents::get_data_for:no data for '" << mname << "'\n";
    return 0;
  } 

  if (! _input.seekg(f->second))
  {
    cerr << "Ofile_Contents::get_data_for:cannot seek to " << f->second << " for '" << mname << "'\n";
    return 0;
  }

  return get_next_molecule();
}

int
Ofile_Contents::_another_fingerprint (const_IWSubstring buffer,     // note local copy
                                      const int nchars,
                                      resizable_array_p<Set_of_Bit_and_Class_Weights> & fp)
{
  assert (_nclasses > 1);

  buffer.remove_leading_chars(nchars);

  unsigned int ndx;

  if (! buffer.numeric_value(ndx) || ndx < 0)
  {
    cerr << "Ofile_Contents::_another_fingerprint:invalid index '" << buffer << "'\n";
    return 0;
  }

  if (ndx != fp.size())
    cerr << "Ofile_Contents::_another_fixed_width_fingerprint:ordering inconsistency within input, possible serious problems\n";

  auto f = new Set_of_Bit_and_Class_Weights();

  if (! f->build (_nclasses, _input))
  {
    cerr << "Ofile_Contents::_another_fixed_width_fingerprint:cannot read fingerprint data\n";
    delete f;
    return 0;
  }

  fp.add(f);

  return fp.size();
}

int
Ofile_Contents::_parse_score_record (const const_IWSubstring & buffer)
{
  if (_nclasses != buffer.nwords() - 1)
  {
    cerr << "Ofile_Contents::_parse_score_record:invalid score record '" << buffer << "', have " << _nclasses << " classes\n";
    return 0;
  }

  if (nullptr == _score)
    _score = new float[_nclasses];

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);
  assert("Score" == token);

  for (auto c = 0; c < _nclasses; ++c)
  {
    buffer.nextword(token, i);
    if (! token.numeric_value(_score[c]))
    {
      cerr << "Ofile_Contents::_parse_score_record:invalid score value '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Ofile_Contents::get_contribution (const unsigned int b,
                                  double * c) const
{
  for (auto i = 0; i < _sparse_fingerprints; ++i)
  {
    if (_sparse_fingerprints[i]->get_contribution(b, c))
      return 1;
  }

  for (auto i = 0; i < _fixed_fingerprints; ++i)
  {
    if (_fixed_fingerprints[i]->get_contribution(b, c))
      return 1;
  }

  return 0;
}

class Bfile_Contents
{
  private:
    IW_STL_Hash_Map_off_t _id_to_offset;

    int _nclasses;
    IWString _name;

    resizable_array<atom_number_t> _atom_number;
    resizable_array<int> _radius;
    resizable_array<unsigned int> _bit;
    resizable_array_p<IWString> _smarts;

    iwstring_data_source _input;

//  private functions

    int _add_bfile_record (const const_IWSubstring & buffer, const int matoms);
    int _determine_id_offsets ();

  public:
    int open (const char *);

    const IWString & name () const { return _name;}

    unsigned int get_next_molecule (const int matoms);
    unsigned int get_data_for (const IWString &, const int matoms);

    int get_bit (const atom_number_t, const int r, unsigned int & b) const;
};

int
Bfile_Contents::open (const char * fname)
{
  if (! _input.open(fname))
  {
    cerr << "Bfile_Contents::open:cannot open '" << fname << "\n";
    return 0;
  } 

  if (! _determine_id_offsets ())
    return 0;

  if (! _input.seekg(0))
  {
    cerr << "Bfile_Contents::open:cannot seek back to start in '" << fname << "'\n";
    return 0;
  }

  return _id_to_offset.size();
}

int
Bfile_Contents::_determine_id_offsets ()
{
  const_IWSubstring buffer;

  int next_is_id = 1;

  for (auto o = _input.tellg(); _input.next_record(buffer); o = _input.tellg())
  {
    if (next_is_id)
      _id_to_offset[buffer] = o;
    else if ('|' == buffer)
      next_is_id = 1;
  }

  return _id_to_offset.size();
}

/*
  Records look like
    10672
    0 0 168 [cD3H0v4;r5;r6]
    0 1 3704893440 [cD3H0v4;r5;r6]
*/

unsigned int
Bfile_Contents::get_next_molecule (const int matoms)
{
  if (! _input.next_record(_name))
    return 0;

  const_IWSubstring buffer;

  _atom_number.resize_keep_storage(0);
  _radius.resize_keep_storage(0);
  _bit.resize_keep_storage(0);
  _smarts.resize_keep_storage(0);

  while (_input.next_record(buffer))
  {
    if ("|" == buffer)
      break;

    if (! _add_bfile_record (buffer, matoms))
    {
      cerr << "Bfile_Contents::get_next_molecule:invalid record '" << buffer << "', molecule " << _name << endl;
      return 0;
    }
  }

  if (0 == _atom_number.size())
  {
    cerr << "Bfile_Contents::get_next_molecule:no data\n";
    return 0;
  }

  return _atom_number.size();
}

unsigned int
Bfile_Contents::get_data_for (const IWString & mname,
                              const int matoms)
{
  const auto f = _id_to_offset.find(mname);

  if (f == _id_to_offset.end())
  {
    cerr << "Bfile_Contents::get_data_for:no data for '" << mname << "'\n";
    return 0;
  } 

  if (! _input.seekg(f->second))
  {
    cerr << "Bfile_Contents::get_data_for:cannot seek to " << f->second << " for '" << mname << "'\n";
    return 0;
  }

  return get_next_molecule(matoms);
}

int
Bfile_Contents::_add_bfile_record (const const_IWSubstring & buffer,
                                   const int matoms)
{
  int i = 0;
  const_IWSubstring token;

  if (buffer.nwords() < 3)
    return 0;

  buffer.nextword(token, i);
  atom_number_t a;
  if (! token.numeric_value(a) || a < 0 || a >= matoms)
  {
    cerr << "Bfile_Contents::_add_bfile_record:invalid atom number, molecule contains " << matoms << " atoms\n";
    return 0;
  }

  buffer.nextword(token, i);
  int r;
  if (! token.numeric_value(r) || r < 0)
  {
    cerr << "Bfile_Contents::_add_bfile_record:invalid radius\n";
    return 0;
  }

  buffer.nextword(token, i);
  unsigned int b;
  if (! token.numeric_value(b))
  {
    cerr << "Bfile_Contents::_add_bfile_record:invalid bit number\n";
    return 0;
  }

  _atom_number.add(a);
  _radius.add(r);
  _bit.add(b);

  return 1;
}

int
Bfile_Contents::get_bit (const atom_number_t zatom,
                         const int rad,
                         unsigned int & b) const
{
  const auto n = _atom_number.size();
//cerr << "Bfile_Contents::get_bit:searching " << n << " atoms for " << zatom << " rad " << rad << endl;

  const auto f = std::lower_bound(_atom_number.cbegin(), _atom_number.cend(), zatom);

  if (f == _atom_number.cend())    // should not happen
    return 0;

  auto i = f - _atom_number.cbegin();    // we need an index

//cerr << "Found atom " << zatom << " radius " << rad << " at i = " << i << endl;

  while (1)
  {
    if (rad == _radius[i])
    {
      b = _bit[i];
      return 1;
    }

    i++;
    if (i >= n)
      return 0;

    if (_atom_number[i] != zatom)
      return 0;
  }
}

static int
score_to_isotope (double v)
{
  if (linear_between_min_and_max)
    return static_cast<int>((v - min_score) / (max_score - min_score) * 100.0 + 0.4999);

  assert (max_score > 0.0 && min_score < 0.0);

  double range;             // we scale to the largest value, positive or negative
  if (max_score > - min_score)
    range = max_score;
  else
    range = - min_score;

//cerr << "Range is " << range << endl;

  if (v >= 0.0)
    return 50 + static_cast<int>(v / range * 49.0 + 0.4999);

  return 50 - static_cast<int>((v - min_score) / range * 49.0 + 0.49999);
}

class Set_of_Colours
{
  private:
    int _n;
    IWString * _colour;

  public:
    Set_of_Colours ();
    ~Set_of_Colours ();

    int build (const char * fname);
    int build (iwstring_data_source & input);

    const IWString & colour_for_value (double v) const;
};

Set_of_Colours::Set_of_Colours ()
{
  _n = 0;
  _colour = nullptr;

  return;
}

Set_of_Colours::~Set_of_Colours ()
{
  if (nullptr != _colour)
    delete [] _colour;

  _n = -11;

  return;
}

int
Set_of_Colours::build (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Set_of_Colours::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build (input);
}

int
Set_of_Colours::build(iwstring_data_source & input)
{
  input.set_dos(1);

  int n = input.records_remaining();

//cerr << "Colour file has " << n << " records\n";

  _colour = new IWString[n];

  int ndx = 0;

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
//  if (buffer.starts_with("#"))
//    continue

    if (buffer.starts_with('#') && buffer.length() > 7)
      buffer.iwtruncate(7);

    _colour[ndx] = buffer;
//  cerr << "Just got colour " << _colour[ndx] << ", ndx " << ndx << endl;
    ndx++;
  }

  _n = ndx;

  return _n;
}

const IWString &
Set_of_Colours::colour_for_value(double v) const
{
  if (linear_between_min_and_max)
  {
    int ndx = static_cast<int>((v - min_score) / (max_score - min_score) * _n + 0.4999);

    assert (ndx >= 0 && ndx < _n);

    return _colour[ndx];
  }

  double range;

  if (max_score > - min_score)
    range = max_score;
  else
    range = - min_score;

  int ndx;
  if (v >= 0.0)
    ndx = _n/2 + static_cast<int>(v / range * (_n / 2) + 0.4999);
  else
    ndx = static_cast<int>((v - min_score) / range * (_n/2) + 0.4999);

#define DEBUG_COLOUR_FOR_VALUE
#ifdef DEBUG_COLOUR_FOR_VALUE
  cerr << "Converted value " << v << " to index " << ndx << endl;
#endif

  return _colour[ndx];
}

static Set_of_Colours set_of_colours;

static int
do_write_isotopically_labelled (Molecule & m,
                                const int nclasses,
                                int c,
                                const double * atom_contribution,
                                IWString_and_File_Descriptor & output)
{
  const auto matoms = m.natoms();

//cerr << min_score << " min_score " << max_score << " max_score\n";
  for (auto i = 0; i < matoms; ++i)
  {
    double v = atom_contribution[i * nclasses + c];

    int iso = score_to_isotope(v);

//  cerr << "Atom " << i << " class " <<  c << " contribution " << v << " iso " << iso << endl;

    m.set_isotope(i, iso);
  }

  output << m.smiles() << ' ' << m.name() << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
do_write_isotopically_labelled (Molecule & m,
                                const int nclasses,
                                const double * atom_contribution,
                                IWString_and_File_Descriptor & output)
{
  IWString mname(m.name());
  for (auto i = 0; i < nclasses; ++i)
  {
    IWString tmp;
    tmp << mname << '.' << i;

    m.set_name(tmp);
    do_write_isotopically_labelled(m, nclasses, i, atom_contribution, output);
  }

  m.set_name(mname);
  m.transform_to_non_isotopic_form();

  return 1;
}

static int
fingerprint_weights_to_atom_colours (Molecule & m,
                                     const int nclasses,
                                     const Ofile_Contents & ofile,
                                     const double * atom_contribution,
                                     ostream & structure_output,
                                     IWString_and_File_Descriptor & direct_color_file)
{
  if (stream_for_isotopically_labelled.is_open())
    do_write_isotopically_labelled(m, nclasses, atom_contribution, stream_for_isotopically_labelled);

//MDL_File_Supporting_Material * mdlfos = global_default_MDL_File_Supporting_Material();
  
  if (stream_for_atomic_contributions.is_open())
    stream_for_atomic_contributions << m.smiles() << ' ' << m.name() << '\n';

  for (auto c = 0; c < nclasses; ++c)
  {
    const IWString nsave(m.name());
    IWString tmp;
    tmp << m.name() << '.' << ofile.class_name(c) << ' ' << ofile.score(c);

    m.set_name(tmp);
    m.write_molecule_mdl(structure_output, "");
    m.set_name(nsave);

    const auto matoms = m.natoms();

    for (auto i = 0; i < matoms; ++i)
    {
      double v = atom_contribution[i * nclasses + c];
      direct_color_file << i << ' ' << set_of_colours.colour_for_value(v) << '\n';

      if (stream_for_atomic_contributions.is_open())
        stream_for_atomic_contributions << i << ' ' << static_cast<float>(v) << '\n';

      const Atom * a = m.atomi(i);

      const auto acon = a->ncon();

      for (auto j = 0; j < acon; ++j)
      {
        const auto k = a->other(i, j);

        if (k < i)
          continue;

        const double v2 = (v + atom_contribution[k * nclasses + c]) * 0.5;

        direct_color_file << i << ' ' << k << ' ' << set_of_colours.colour_for_value(v2) << '\n';
      }
    }

    direct_color_file << "|\n";
  }

  direct_color_file.write_if_buffer_holds_more_than(8192);

  stream_for_atomic_contributions << "|\n";
  stream_for_atomic_contributions.write_if_buffer_holds_more_than(8192);

  return structure_output.good();
}

static int
fingerprint_weights_to_atom_colours_2 (Molecule & m,
                                       Bfile_Contents & bfile,
                                       Ofile_Contents & ofile,
                                       ostream & structure_output,
                                       IWString_and_File_Descriptor & direct_color_file)
{
  const auto matoms = m.natoms();

  const auto nclasses = ofile.nclasses();

  double * atom_contribution = new double[matoms * nclasses]; unique_ptr<double[]> free_atom_contribution(atom_contribution);
  std::fill_n(atom_contribution, matoms * nclasses, 0.0);

  std::fill_n(atom_contribution, matoms, 0.0);    // not necessary

  for (auto i = 0; i < matoms; ++i)
  {
    for (auto r = min_radius; r <= max_radius; ++r)
    {
      unsigned int b;
      if (! bfile.get_bit(i, r, b))
      {
        cerr << "fingerprint_weights_to_atom_colours_2:no bit in Bfile: " << m.name() << " atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " radius " << r << endl;
        continue;
      }

      if (! ofile.get_contribution(b, atom_contribution + i * nclasses))   
      {
        cerr << "fingerprint_weights_to_atom_colours_2:no bit " << b << " in Ofile: " << m.name() << " atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " radius " << r << endl;
        continue;
      }

    }
  }

  for (auto i = 0; i < matoms * nclasses; ++i)
  {
    if (atom_contribution[i] < min_score)
      atom_contribution[i] = min_score;
    else if (atom_contribution[i] > max_score)
      atom_contribution[i] = max_score;
  }

  return fingerprint_weights_to_atom_colours(m, nclasses, ofile, atom_contribution, structure_output, direct_color_file);
}

static int
fingerprint_weights_to_atom_colours (Molecule & m,
                                     Bfile_Contents & bfile,
                                     Ofile_Contents & ofile,
                                     ostream & structure_output,
                                     IWString_and_File_Descriptor & direct_color_file)
{
  const auto matoms = m.natoms();

#ifdef DEBUG_FINGERPRINT_WEIGHTS_TO_ATOM_COLOURS
  cerr << "Read " << m.name() << " with " << matoms << " atoms\n";
#endif

  if (! bfile.get_data_for(m.name(), matoms))
  {
    cerr << "Premature EOF on Bfile\n";
    return 0;
  }

  if (! ofile.get_data_for(m.name()))
  {
    cerr << "Premature EOF on Ofile\n";
    return 0;
  }

  if (bfile.name() != ofile.molecule_name())
  {
    cerr << "fingerprint_weights_to_atom_colours::name mismatch between bfile '" << bfile.name() << "' and ofile '" << ofile.molecule_name() << "'\n";
    return 0;
  }

  if (m.name() != bfile.name())
  {
    cerr << "fingerprint_weights_to_atom_colours::name mismatch between bfile '" << bfile.name() << "' and molecule '" << m.name() << "'\n";
    return 0;
  }

  if (ofile.number_fingerprints() > 1)
  {
    cerr << "Sorry, do not know how to handle more than one fingerprint, see Ian\n";
    return 0;
  }

  return fingerprint_weights_to_atom_colours_2(m, bfile, ofile, structure_output, direct_color_file);
}

static int
fingerprint_weights_to_atom_colours (data_source_and_type<Molecule> & input,
                                     Bfile_Contents & bfile,
                                     Ofile_Contents & ofile,
                                     ostream & structure_output,
                                     IWString_and_File_Descriptor & direct_color_file)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! fingerprint_weights_to_atom_colours(*m, bfile, ofile, structure_output, direct_color_file))
      return 0;
  }

  return 1;
}

static int
fingerprint_weights_to_atom_colours (const char * fname, FileType input_type, 
                                     Bfile_Contents & bfile,
                                     Ofile_Contents & ofile,
                                     ostream & structure_output,
                                     IWString_and_File_Descriptor & direct_color_file)
{
  assert (NULL != fname);

  if (input_type == FILE_TYPE_INVALID)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return fingerprint_weights_to_atom_colours(input, bfile, ofile, structure_output, direct_color_file);
}

static int
fingerprint_weights_to_atom_colours(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lO:B:r:R:x:X:I:C:D:c");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }
  else
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('c'))
  {
    remove_all_chiral_centres = 1;

    if (verbose)
      cerr << "Will discard chirality\n";
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = FILE_TYPE_SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (! cl.option_present('O'))
  {
    cerr << "Must specify file of bit contributions via the -O file\n";
    usage(2);
  }

  if (! cl.option_present('B'))
  {
    cerr << "Must specify file of bit descriptions via the -B option\n";
    usage(2);
  }

  Bfile_Contents bfile;

  if (cl.option_present('B'))
  {
    const char * b = cl.option_value('B');
    if (! bfile.open(b))
    {
      cerr << "Cannot open -B file '" << b << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "File for bit meanings initialised '" << b << "'\n";
  }

  Ofile_Contents ofile;

  if (cl.option_present('O'))
  {
    const char * o = cl.option_value('O');

    if (! ofile.open(o))
    {
      cerr << "Cannot open file with bit weights (-O) '" << o << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Opened file with bit weights '" << o << "'\n";
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', min_radius) || min_radius < 0)
    {
      cerr << "The minimum radius (-r) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Min radius set to " << min_radius << endl;
  }

  if (cl.option_present('R'))
  {
    if (! cl.value('R', max_radius) || max_radius < min_radius)
    {
      cerr << "The maximum radius (-R) must be a whole +ve number, greater than min radius " << min_radius << "\n";
      usage(1);
    }

    if (verbose)
      cerr << "Max radius set to " << max_radius << endl;
  }

  if (! cl.option_present('C'))
  {
    cerr << "Must specify file of colours via the -C option\n";
    usage(1);
  }

  if (cl.option_present('C'))
  {
    const char * c = cl.option_value('C');

    if (! set_of_colours.build(c))
    {
      cerr << "Cannot initialise colours from '" << c << "'\n";
      return 2;
    }
  }

  int xopts = 0;

  if (cl.option_present('x'))
  {
    if (! cl.value('x', min_score))
    {
      cerr << "The minimum score (-x) must be a valid floating point number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Min score set to " << min_score << endl;

    xopts++;
  }

  if (cl.option_present('X'))
  {
    if (! cl.value('X', max_score) || max_score < min_score)
    {
      cerr << "The maximum score (-X) must be a valid float number, greater than min min_score " << min_score << "\n";
      usage(1);
    }

    if (verbose)
      cerr << "Max score set to " << max_score << endl;

    xopts++;
  }

  if (0 == xopts)
    ;
  else if (2 != xopts)
  {
    cerr << "Sorry, if you are going to specify a min (-x) or max (-X) value, then you must specify both\n";
    return 1;
  }
  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor dcf;

  if (! cl.option_present('D'))
  {
    cerr << "Must specify name of the direct colour file for cactvs\n";
    usage(1);
  }

  if (cl.option_present('D'))
  {
    const char * d = cl.option_value('D');

    if (! dcf.open(d))
    {
      cerr << "Cannot open direct colour file '" << d << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "csib direct color file '" << d << "' opened\n";
  }

  if (cl.option_present('I'))
  {
    IWString i = cl.string_value('I');
    if (! i.ends_with(".smi"))
      i << ".smi";

    if (! stream_for_isotopically_labelled.open(i.null_terminated_chars()))
    {
      cerr << "Cannot open stream for isotopically labelled molecules '" << i << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Isotopically labelled molecules written to '" << i << "'\n";
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! fingerprint_weights_to_atom_colours(cl[i], input_type, bfile, ofile, std::cout, dcf))
    {
      rc = i + 1;
      break;
    }
  }

  std::cout.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = fingerprint_weights_to_atom_colours (argc, argv);

  return rc;
}
