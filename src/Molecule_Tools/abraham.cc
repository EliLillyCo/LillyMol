/*
  Implementation of Abraham Descriptors.
  Initial (wrong) implementation from Richard Lewis

  James A Platts, Darko Butina, Michael Abraham Anne Hersey. 
  J Chem Inf Comput Sci 1999, 39, 835-845
*/

#include <iostream>
#include <memory>
#include <algorithm>
using std::cerr;
using std::endl;

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "cmdline.h"
#include "iwdigits.h"
#include "misc.h"
#include "iw_stl_hash_map.h"
#include "accumulator.h"
#include "sparse_fp_creator.h"

#include "istream_and_type.h"
#include "molecule.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "qry_wstats.h"
#include "aromatic.h"
#include "output.h"
#include "target.h"

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

/*
  We can produce fingerprints. Right now, we only produce some bits
*/

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString fingerprint_tag;

static int bit_replicates = 9;
static double bit_count_dynamic_range = 9.0;

static int function_as_tdt_filter = 0;

static IWString descriptor_prefix("abr_");

// With the -H option, we make implicit hydrogens explicit

static int make_implicit_hydrogens_implicit = 0;

/*
  Since these computations do not come from Queries_and_Additive_Models objects, we
  need to keep their info here
*/

static double min_Vx = 200.0;
static double max_Vx = 700.0;

static double min_logp = -3.00;
static double max_logp = 7.00;

static char output_separator = ' ';

static Accumulator<double> * acc = nullptr;

/*
  A debugging tool to see how queries overlap
*/

static int discern_query_overlaps = 0;

/*
  Sometimes it will make sense to not bother checking
*/

static int check_for_previously_matched_atoms = 1;

/*
  Or we can just run things normally, but give precedence to the last query matched
*/

static int process_queries_first_to_last = 1;

/*
  Indices into the models array of the various types of models
*/

#define ABRAHAM_INDEX 0
#define ALPHA2A_INDEX 1
#define CG1_INDEX 2
#define CG2_INDEX 3
#define CG3_INDEX 4
#define CG4_INDEX 5

static int constantinou_gani_models_present = 0;

/*
  but each Queries_and_Additive_Models creates a different number of model values
  so the results array will contain contributions from various models
  along it. We will interrogate the models to figure out where
*/

static int ap_r2 = -1;
static int ap_pi2h = -1;
static int ap_beta2h = -1;
static int ap_beta2o = -1;
static int ap_alpha2h = -1;

// These two we just have to place manually

static int ap_vx      = -1;
static int ap_logp    = -1;

static int ap_cg[5];

class Abraham_Parameters
{
  protected:
    double * _p;
    const int _n;

  public:
    Abraham_Parameters(int);
    ~Abraham_Parameters();

    void operator += (const Abraham_Parameters & rhs);

    void increment_results (double *) const;

    int write_header (IWString_and_File_Descriptor & output) const;

//  int write_fingerprint (Molecule & m, const IWString & tag, IWString_and_File_Descriptor & output) const;

    double operator [] (int ndx) const { return _p[ndx];}
};

Abraham_Parameters::Abraham_Parameters(int n) : _n(n)
{
  _p = new double[n];

  std::fill_n(_p, n, 0.0);

  return;
}

Abraham_Parameters::~Abraham_Parameters ()
{
  if (nullptr != _p)
    delete [] _p;

  return;
}

/*float
Abraham_Parameters::logp () const
{
  return 0.562 * (_p[AP_R2] + intercept_r2)
         - 1.054 * (_p[AP_PI2H] + intercept_pi2h)
         + 0.034 * (_p[AP_ALPHA2H] + intercept_alpha2h)
         - 3.460 * (_p[AP_BETA2O] + intercept_beta2o)
         + 0.03814 * _p[AP_VX] + 0.088;
}*/

/*static float 
abraham_logp (const double * v)
{
  return    0.562   * (v[AP_R2] + intercept_r2)
          - 1.054   * (v[AP_PI2H] + intercept_pi2h)
          + 0.034   * (v[AP_ALPHA2H] + intercept_alpha2h)
          - 3.460   * (v[AP_BETA2O] + intercept_beta2o)
          + 0.03814 * v[AP_VX] + 0.088;
}*/

static float 
abraham_logp (const double r2,
              const double pi2h,
              const double alpha2h,
              const double beta2o,
              const double vx)
{
//cerr << r2 << ' ' << pi2h << ' ' << alpha2h << ' ' << beta2o << ' ' << vx << endl;

  return    0.562   * r2
          - 1.054   * pi2h
          + 0.034   * alpha2h
          - 3.460   * beta2o
          + 0.03814 * vx + 0.088;
}

static float mcgowan_coefficient[HIGHEST_ATOMIC_NUMBER + 1];

static float
mcgowan(Molecule & m)
{
  float rc = 0.0;

  const auto matoms = m.natoms();

  int hcount = 0;

  for (auto i = 0; i < matoms; ++i)
  {
    const auto z = m.atomic_number(i);

    rc += mcgowan_coefficient[z];

//  cerr << "z = " << z << " mcgowan_coefficient " << mcgowan_coefficient[z] << " sum " << rc << endl;

    hcount += m.hcount(i);
  }

//cerr << rc << " sum so far, hcount " << hcount << " edged " << m.nedges() << endl;

  return rc + hcount * mcgowan_coefficient[1] - m.nedges() * 6.56;
}

static void
convert_to_int_count (const double v,
                      int & i,
                      const double vmin,
                      const double vmax)
{
  if (v <= vmin)
    i = 1;
  else if (v > vmax)
    i = bit_count_dynamic_range + 1;
  else
    i = static_cast<int>((v - vmin) / (vmax - vmin) * bit_count_dynamic_range + 1.4999);

  return;
}

static void
set_bit_replicates (int bstart,
                    int c,
                    int n,
                    Sparse_Fingerprint_Creator & sfc)
{
  if (c < 1)
    c = 1;

  for (int i = 0; i < n; ++i)
  {
    sfc.hit_bit(bstart + i, c);
  }

  return;
}


static int
write_fingerprint (Molecule & m,
                   const IWString & tag,
                   const int * b,
                   const int nb,
                   IWString_and_File_Descriptor & output) 
{
  if (! function_as_tdt_filter)
  {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  Sparse_Fingerprint_Creator sfc;

/*int v = convert_to_int_count(_r2 + intercept_r2);
  for (int j = 0; j < bit_replicates; j++)
  {
    sfc.hit_bit(j, v);
  }*/

  int ndx = 0;    // the actual bit number that gets set

  for (auto i = 0; i < nb; ++i)
  {
    for (auto r = 0; r < bit_replicates; ++r)
    {
      sfc.hit_bit(ndx, b[i]);
      ndx++;
    }
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tag, output);

  output << tmp << '\n';

  if (! function_as_tdt_filter)
    output << "|\n";

  return 1;
}

void
Abraham_Parameters::increment_results (double * v) const
{
  for (auto i = 0; i < _n; ++i)
  {
    v[i] += _p[i];
  }

  return;
}

class Abraham_Substructure_Query : public Substructure_Hit_Statistics,
                                   public Abraham_Parameters
{
  private:
    int _isotope;

//  There are two parameters that control matching. The number of
//  matched atoms to check as being unmatched, and the number of
//  matched atoms that we claim.

    int _number_atoms_to_check;

    int _number_atoms_to_claim;

//  private functions

    template <typename T> int _parse_specifier(const const_IWSubstring &, T &);

  public:
    Abraham_Substructure_Query(int);
    
    int build (const const_IWSubstring &);
    int build (const const_IWSubstring & buffer, const IW_STL_Hash_Map_int & mname_to_ndx, int * recognised);

//  int isotope() const { return _isotope;}

    int place_isotopes(const Set_of_Atoms & e, int * isotope) const;

    int mark_matched_atoms(const Set_of_Atoms & e,
                                               int * already_hit,
                                               int flag) const;

    int all_atoms_unmatched (const Set_of_Atoms & e, const int * already_hit) const;
};

Abraham_Substructure_Query::Abraham_Substructure_Query(int n) : Abraham_Parameters (n)
{
  _isotope = 0;

  _number_atoms_to_check = -1;   // negative means check all

  _number_atoms_to_claim = -1;   // negative means check all

  return;
}

template <typename T>
int
Abraham_Substructure_Query::_parse_specifier(const const_IWSubstring & token,
                                             T & v)
{
  const_IWSubstring directive, vvvvv;

  token.split(directive, '=', vvvvv);

  if (! vvvvv.numeric_value(v))
  {
    cerr << "Abraham_Substructure_Query::_parse_specifier:invalid numeric '" << vvvvv << "'\n";
    return 0;
  }

  return 1;
}

template int Abraham_Substructure_Query::_parse_specifier(const const_IWSubstring &, int &);
template int Abraham_Substructure_Query::_parse_specifier(const const_IWSubstring &, double &);

int
Abraham_Substructure_Query::build (const const_IWSubstring & buffer,
                                   const IW_STL_Hash_Map_int & mname_to_ndx,
                                   int * recognised)
{
  if (buffer.nwords() < 3)   // smiles name directive
  {
    cerr << "Abraham_Substructure_Query::build:must have at least 3 tokens\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);

  if (! create_from_smarts(token))
  {
    cerr << "Abraham_Substructure_Query::build:invalid smarts '" << token << "'\n";
    return 0;
  }

  buffer.nextword(_comment, i);

  std::fill_n(recognised, _n, 0);

  int ncon = -1;

  while (buffer.nextword(token, i))
  {
    if (token.starts_with('#'))    // do not use these any more
      continue;

    if (token == "//")    // inline comment
      break;

    IWString directive, s;
    if (! token.split(directive, '=', s) || 0 == directive.length() || 0 == s.length())
    {
      cerr << "Abraham_Substructure_Query::build:invalid specification token '" << buffer << "' at '" << token << "'\n";
      return 0;
    }

    const auto f = mname_to_ndx.find(directive);

    if (f != mname_to_ndx.end())    // recognised as a model
    {
      const auto ndx = (*f).second;

      if ('*' == s)
        _p[ndx] = 0.0;
      else
      {
        double v;
        if (! s.numeric_value(v))
        {
          cerr << "Abraham_Substructure_Query::build:invalid numeric specification '" << token << "'\n";
          return 0;
        }

        _p[ndx] = v;
      }

      recognised[ndx] = 1;

      continue;
    }

    if (token.starts_with("ISO="))
    {
      if (! _parse_specifier(token, _isotope))
        return 0;
    }
    else if (token.starts_with("NCLAIM="))
    {
      token.remove_leading_chars(7);
      if (! token.numeric_value(_number_atoms_to_claim))
      {
        cerr << "Invalid nclaim '" << buffer << "'\n";
        return 0;
      }
    }
    else if (token.starts_with("NCHECK="))
    {
      token.remove_leading_chars(7);
      if (! token.numeric_value(_number_atoms_to_check))
      {
        cerr << "Invalid ncheck '" << buffer << "'\n";
        return 0;
      }
    }
    else if (token.starts_with("NCON="))
    {
      token.remove_leading_chars(5);
      if (! token.numeric_value(ncon) || ncon < 0)
      {
        cerr << "Invalid ncon '" << buffer << "'\n";
        return 0;
      }
    }
    else if (token.starts_with('#'))   // numeric identifiers
      _comment << ' ' << token;
    else
      cerr << "Abraham_Substructure_Query::build:ignoring unrecognised token '" << token << "'\n";
  }

  set_find_unique_embeddings_only(1);

  if (ncon >= 0)
    set_ncon(ncon);

  if (recognised + _n == std::find(recognised, recognised + _n, 0))    // good, no zero's in the range
    ;
  else
  {
    cerr << "Abraham_Substructure_Query::build:one or more models not recognised\n";
    for (auto i : mname_to_ndx)
    {
      cerr << i.first << ' ' << recognised[i.second] << endl;
    }

    return 0;
  }

  return 1;
}

static IWDigits digits;

static int use_query_name_as_descriptor_name = 0;

static extending_resizable_array<int> unclassified_atom_count;

static IWString descriptor_name;

static Molecule_Output_Object stream_for_labeled_atoms;

/*
  There are 3 types of fingerprints we can write.
    1. Bit for each query
    2. Bit for each ISO= type (the label [] array)
*/

static IWString bit_for_each_query_tag;
static IWString atom_label_tag;

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __TIME__ << " " << __DATE__ << endl;
  cerr << "Usage: " << prog_name << " <options> <file1> <file2> ...\n";
  cerr << " Abraham Descriptors\n";
  cerr << "  -F <file>      specify control file of queries and values\n";
  cerr << "  -P <file>      control file for alpha2H queries\n";
  cerr << "  -C <fname>     control file for ConstantinouGani queries\n";
  cerr << "  -h             make implicit hydrogens explicit\n";
  cerr << "  -L <fname>     write molecules with labeled atoms to <fname>\n";
  cerr << "  -o <type>      specify output file type(s) for -L file\n";
  cerr << "  -G ...         miscellaneous options, enter '-G help' for info\n";
  cerr << "  -S <fname>     specify output stream (default is stdout)\n";
  cerr << "  -J <tag>       fingerprint tag\n";
  cerr << "  -p <int>       bit replicates when producing fingerprints\n";
  cerr << "  -y m,n         display query matches (m) and/or non matches (n)\n";
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  (void) display_standard_aromaticity_options(cerr);
  cerr << "  -f             work as a tdt filter\n";
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static void
display_dash_G_qualifiers (std::ostream & os)
{
  os << " -G blki         aromatic bonds lose their Kekule identity\n";
  os << " -G qndn         use query names as descriptor names\n";
  os << " -G tovl         run in a mode to test query overlaps\n";
  os << " -G ncpm         do NOT check previously matched atoms\n";
  os << " -G ql2f         process queries last to first\n";

  return;
}

#ifdef NOT_USED_ANY_MORE_QWEQWEQWE
static int
write_the_identifier (const IWString & mname,
                      IWString_and_File_Descriptor & output)
{
  if (1 == mname.nwords())
  {
    output << mname;
    return 1;
  }

  output << mname.word(0);

  return output.good();
}
#endif

/*
  Do any of the atoms in this embedding hit atoms already hit.
  we need to do this carefully because there may be -1's in E
*/

int
Abraham_Substructure_Query::all_atoms_unmatched(const Set_of_Atoms & e,
                     const int * already_hit) const
{
  if (! check_for_previously_matched_atoms)
    return 1;

  if (0 == _number_atoms_to_check)   // don't check anything
    return 1;

  auto n = e.number_elements();

  if (_number_atoms_to_check < 0)
    ;
  else if (n > _number_atoms_to_check)
    n = _number_atoms_to_check;

#ifdef DEBUG_ALL_ATOMS_UNMATCHED
  cerr << "all_atoms_unmatched:checking " << n << " atoms\n";
#endif

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = e[i];

#ifdef DEBUG_ALL_ATOMS_UNMATCHED
    cerr << "Matched atom i = " << i << " j = " << j << endl;
#endif

    if (INVALID_ATOM_NUMBER == j)
      continue;

    if (already_hit[j])
      return 0;
  }

  return 1;   // all unmatched
}

template class resizable_array_p<Abraham_Substructure_Query>;
template class resizable_array_base<Abraham_Substructure_Query *>;

class Queries_and_Additive_Models
{
  private:
    int _nmodels;

    IWString * _model_name;

    double * _intercept;

    int * _result_column;    // where to put our results

    resizable_array_p<Abraham_Substructure_Query> _queries;

//  when producing fingerprints, we need to know the min and max possible values for each model

    double * _global_min;
    double * _global_max;

    IW_STL_Hash_Map_int _model_name_to_number;

    IWString _name;

    int _show_all_query_matches, _show_all_query_non_matches;

//  While not a great idea to make this part of the object, it is already not thread safe
//  due to the substructure search objects

    int _molecules_processed;

    extending_resizable_array<int> _unclassified_atoms;

    int * _hit_matrix;

//  private functions

    int _parse_intercept_record (const const_IWSubstring &);
    int _parse_global_min_max (const IWString & buffer, double * v);

    int _process(Molecule_to_Match & target, 
                 int * already_hit,
                 int * isotope,
                 double * results);

  public:
    Queries_and_Additive_Models ();
    ~Queries_and_Additive_Models ();

    int build (const char * fname);
    int build (iwstring_data_source & input);

    void set_name (const char * s) { _name = s;}

    void set_result_column(int &);

    void set_show_all_query_matches     (int s) {_show_all_query_matches = s;}
    void set_show_all_query_non_matches (int s) {_show_all_query_non_matches = s;}

    int examine_query_overlap_behaviour(Molecule & m);

    template <typename T> int append_model_names (const IWString & prefix, T & output) const;

    int nmodels () const { return _nmodels;}

    int result_column (const char *) const;

    int number_queries () const { return _queries.number_elements();}

    int process (Molecule & m, Molecule_to_Match & target, int *, int *, double *);

    void convert_to_int_count (const double * v, int *) const;

    int report (const Accumulator<double> * acc, std::ostream & os) const;

    int report_overlaps(std::ostream & output) const;
};

Queries_and_Additive_Models::Queries_and_Additive_Models ()
{
  _nmodels = 0;
  _model_name = nullptr;
  _intercept = nullptr;
  _result_column = nullptr;
  _global_min = nullptr;
  _global_max = nullptr;

  _show_all_query_matches = 0;
  _show_all_query_non_matches = 0;

  _molecules_processed = 0;

  _hit_matrix = nullptr;

  return;
}

Queries_and_Additive_Models::~Queries_and_Additive_Models ()
{
  if (nullptr != _model_name)
    delete [] _model_name;

  if (nullptr != _intercept)
    delete [] _intercept;

  if (nullptr != _result_column)
    delete [] _result_column;

  if (nullptr != _global_min)
    delete [] _global_min;

  if (nullptr != _global_max)
    delete [] _global_max;

  if (nullptr != _hit_matrix)
    delete [] _hit_matrix;

  return;
}

int
Queries_and_Additive_Models::result_column (const char * s) const
{
  IWString mname(s);

  const auto f = _model_name_to_number.find(mname);

  if (f == _model_name_to_number.end())
    return -1;

  return _result_column[f->second];
}

void
Queries_and_Additive_Models::set_result_column (int & ndx)
{
  for (int i = 0; i < _nmodels; ++i)
  {
    _result_column[i] = ndx;
    ndx++;
  }

  return;
}

template <typename T>
int
Queries_and_Additive_Models::append_model_names (const IWString & prefix, T & output) const
{
  for (auto i = 0; i < _nmodels; ++i)
  {
    output << output_separator << prefix << _model_name[i];
  }

  return 1;
}

int
Queries_and_Additive_Models::build (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Queries_and_Additive_Models::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input);
}

/*
  First we must have the header records

# INTERCEPT model1=x model2=x

*/

int
Queries_and_Additive_Models::_parse_intercept_record (const const_IWSubstring & buffer)
{
  const auto nw = buffer.nwords();

  if (nw <= 2)
  {
    cerr << "Queries_and_Additive_Models::_parse_intercept_record:invalid intercept specification '" << buffer << "'\n";
    return 0;
  }

  _nmodels = nw - 2;

  _intercept = new double[_nmodels];
  _global_min = new double[_nmodels];
  _result_column = new int[_nmodels];
  _global_max = new double[_nmodels];
  _model_name = new IWString[_nmodels];

  const_IWSubstring token;
  int i = 0;

  buffer.nextword(token, i);     // '#' token

  buffer.nextword(token, i);
  if ("INTERCEPT" != token)
  {
    cerr << "Queries_and_Additive_Models::_parse_intercept_record:invalid intercetp record '" << buffer << "'\n";
    return 0;
  }

  while (buffer.nextword(token, i))
  {
    IWString mname, s;
    double v;
    if (! token.split(mname, '=', s) || 0 == mname.length() || 0 == s.length() || ! s.numeric_value(v))
    {
      cerr << "Queries_and_Additive_Models::_parse_intercept_record:invalid model specification '" << buffer << "'\n";
      return 0;
    }

    if (_model_name_to_number.contains(mname))
    {
      cerr << "Queries_and_Additive_Models::_parse_intercept_record:duplicate model specification '" << buffer << "'\n";
      return 0;
    }

    const int ndx = _model_name_to_number.size();

    assert (ndx >= 0 && ndx < _nmodels);

    _model_name_to_number[mname] = ndx;
    _model_name[ndx] = mname;

    _intercept[ndx] = v;
  }

#ifdef ECHO_INTERCEPTS
  for (auto i = 0; i < _nmodels; ++i)
  {
    cerr << "Model " << i << " intercept " << _intercept[i] << endl;
  }
#endif

  return _nmodels;
}

int
Queries_and_Additive_Models::build (iwstring_data_source & input)
{
  const_IWSubstring buffer;

  IWString intercepts, global_min, global_max;

  while (input.next_record(buffer))
  {
    if (0 == buffer.length())
      continue;

    if (! buffer.starts_with('#'))
    {
      input.push_record();
      break;
    }

    if (buffer.starts_with("# INTERCEPT"))
      intercepts = buffer;
    else if (buffer.starts_with("# MIN"))
      global_min = buffer;
    else if (buffer.starts_with("# MAX"))
      global_max = buffer;
    else if (buffer.starts_with("# DM"))
      _show_all_query_matches = 1;
    else if (buffer.starts_with("# DN"))
      _show_all_query_non_matches = 1;
    else if (buffer.starts_with("# NAME"))
    {
      _name = buffer;
      _name.remove_leading_chars(6);
    }
  }

  if (0 == intercepts.length())
  {
    cerr << "Queries_and_Additive_Models::build:did not find intercepts record\n";
    return 0;
  }

  if (! _parse_intercept_record(intercepts))
  {
    cerr << "Queries_and_Additive_Models::build:cannot parse intercepts record '" << intercepts << "'\n";
    return 0;
  }

  if (0 == global_max.length() && 0 == global_max.length())
    ;
  else if (! _parse_global_min_max(global_max, _global_max) || ! _parse_global_min_max(global_min, _global_min))
  {
    cerr << "Queries_and_Additive_Models::build:cannot determine model global min/max\n";
    cerr << "MIN: " << global_min << endl;
    cerr << "MAX: " << global_max << endl;
    return 0;
  }

  int * recognised = new int[_nmodels]; std::unique_ptr<int[]> free_recognised(recognised);

  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;

    if (buffer.starts_with("//"))
      continue;

    if (0 == buffer.length())
      continue;

    Abraham_Substructure_Query * q = new Abraham_Substructure_Query(_nmodels);

    if (! q->build(buffer, _model_name_to_number, recognised))
    {
      cerr << "Queries_and_Additive_Models::build:invalid query '" << buffer << "', line " << input.lines_read() << endl;
      delete q;

      return 0;
    }

    _queries.add(q);
  }

  return 1;
}

int
Queries_and_Additive_Models::_parse_global_min_max (const IWString & buffer,
                                               double * v)
{
  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);    // skip #
  buffer.nextword(token, i);    // skip MIN/MAX

  for (auto col = 0; buffer.nextword(token, i); ++col)
  {
    IWString mname, s;

    if (! token.split(mname, '=', s) || 0 == mname.length() || 0 == s.length())
    {
      cerr << "Queries_and_Additive_Models::_parse_global_min_max:invalid token form '" << buffer << "'\n";
      return 0;
    }

    double x;
    if (! s.numeric_value(x))
    {
      cerr << "Queries_and_Additive_Models::_parse_global_min_max:invalid numeric '" << buffer << "'\n";
      return 0;
    }

    const auto f = _model_name_to_number.find(mname);

    if (f == _model_name_to_number.end())
    {
      cerr << "Queries_and_Additive_Models::_parse_global_min_max:unrecognised model '" << mname << "' from '" << buffer << "'\n";
      return 0;
    }

    v[(*f).second] = x;
  }

  return 1;
}

//#define DEBUG_PROCESS

int
Queries_and_Additive_Models::process(Molecule & m,
                                 Molecule_to_Match & target,
                                 int * already_hit,
                                 int * isotope,
                                 double * v)
{
  _molecules_processed++;

  std::copy_n(_intercept, _nmodels, v + _result_column[0]);

  const auto matoms = m.natoms();

  std::fill_n(already_hit, matoms, 0);
  std::fill_n(isotope,     matoms, 0);

  const auto rc = _process(target, already_hit, isotope, v);

  int unclassified = 0;

  if (verbose || stream_for_labeled_atoms.active())
  {
    unclassified = std::count_if(already_hit, already_hit + matoms, [] (const int s) { return 0 == s;});
/*  for (int i = 0; i < matoms; ++i)
    {
      cerr << " atom " << i << " already_hit " << already_hit[i] << ' ' << m.smarts_equivalent_for_atom(i) << endl;
    }
    cerr << unclassified << " unclassified\n"; */

    _unclassified_atoms[unclassified]++;
  }

  if (stream_for_labeled_atoms.active())
  {
    Molecule mcopy(m);
    IWString tmp;
    tmp << m.name() << ' ' << _name << " UNC " << unclassified;
    mcopy.set_name(tmp);
    mcopy.set_isotopes(isotope);
    stream_for_labeled_atoms.write(mcopy);
  }

#ifdef DEBUG_PROCESS
  cerr << "Queries_and_Additive_Models::process: finished\n";
  for (auto i = 0; i < _nmodels; ++i)
  {
    cerr << " model " << i << " value " << v[i] << endl;
  }
#endif

  return rc;
}

int
Queries_and_Additive_Models::report (const Accumulator<double> * acc,
                                     std::ostream & output) const
{
  output << "Set of Queries " << _name << " processed " << _molecules_processed << " molecules\n";

  for (int i = 0; i < unclassified_atom_count.number_elements(); i++)
  {
    if (_unclassified_atoms[i])
      output << _unclassified_atoms[i] << " molecules had " << i << " unclassified atoms\n";
  }

  if (0 == _molecules_processed)
    return 0;

  for (int i = 0; i < _nmodels; ++i)
  {
    const auto & ai = acc[_result_column[i]];

    if (0 == ai.n())
      continue;

    output << _model_name[i] << ' ' << ai.n() << " results bwt " << ai.minval() << " and " << ai.maxval() << " ave " << static_cast<float>(ai.average()) << '\n';
  }

  return 1;
}

int
Abraham_Substructure_Query::mark_matched_atoms(const Set_of_Atoms & e,
                                               int * already_hit,
                                               int flag) const
{
  if (0 == _number_atoms_to_claim)   // don't mark any
    return 1;

  int n = e.number_elements();

  if (_number_atoms_to_claim < 0)
    ;
  else if (n > _number_atoms_to_claim)
    n = _number_atoms_to_claim;

#ifdef DEBUG_MARK_MATCHED_ATOMS
  cerr << "Marking " << n << " matched atoms\n";
#endif

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = e[i];

#ifdef DEBUG_MARK_MATCHED_ATOMS
    cerr << "How about i = " << i << " j = " << j << endl;
#endif

    if (INVALID_ATOM_NUMBER == j)
      continue;

    already_hit[j] = flag;
  }

  return 1;
}

int
Abraham_Substructure_Query::place_isotopes(const Set_of_Atoms & e, int * isotope) const
{
  if (0 == _isotope)
    return 1;

  return mark_matched_atoms(e, isotope, _isotope);
}

int
Queries_and_Additive_Models::examine_query_overlap_behaviour(Molecule & m)
{
  const int nq = _queries.size();

  if (nullptr == _hit_matrix)
    _hit_matrix = new_int(nq * nq);

  Molecule_to_Match target(&m);

  const int matoms = m.natoms();

  resizable_array<int> * queries_matching_atom = new resizable_array<int>[matoms]; std::unique_ptr<resizable_array<int>[]> free_queries_matching_atom(queries_matching_atom);

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = _queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    for (int j = 0; j < nhits; ++j)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      for (auto k : *e)
      {
        if (k < 0)
          continue;

        queries_matching_atom[k].add_if_not_already_present(i);     // query i has matched atom k
      }
    }
  }

  for (int i = 0; i < matoms; ++i)
  {
    const resizable_array<int> & qmai = queries_matching_atom[i];
    
    const int n = qmai.number_elements();

 // cerr << n << " of " << nq << " queries matched atom " << i << '\n';
    if (n <= 1)
      continue;

    for (int j = 0; j < n; ++j)
    {
      const int q1 = qmai[j];

      for (int k = (j+1); k < n; ++k)
      {
        const int q2 = qmai[k];
        _hit_matrix[q1 * nq + q2]++;
        _hit_matrix[q2 * nq + q1]++;

//      cerr << "Queries " << q1 << " and " << q2 << " both hit atom " << i << endl;
      }
    }
  }

  return 1;
}

int
Queries_and_Additive_Models::report_overlaps(std::ostream & output) const
{
  if (nullptr == _hit_matrix)
  {
    cerr << "Queries_and_Additive_Models::report_overlaps:no data\n";
    return 0;
  }

  output << "Queries_and_Additive_Models::report_overlaps: " << _model_name << '\n';

  const int nq = _queries.number_elements();

  for (int i = 0; i < nq; ++i)
  {
    output << i << ':';
    for (int j = 0; j < nq; ++j)
    {
      output << ' ' << _hit_matrix[i * nq + j];
    }
    output << '\n';
  }

  return 1;
}

int
Queries_and_Additive_Models::_process(Molecule_to_Match & target,
                                 int * already_hit,
                                 int * isotope,
                                 double * results)
{
  const int nq = _queries.size();

#ifdef DEBUG_PROCESS
  cerr << "On entry, processing " << nq << " queries " << target.molecule()->smiles() << " mih " << make_implicit_hydrogens_implicit << endl;
  for (auto k = 0; k < _nmodels; ++k)
  {
    cerr << " " << k << " " << results[k] << endl;
  }
#endif

  int istart, istep, istop;

  if (process_queries_first_to_last)
  {
    istart = 0;
    istop = nq;
    istep = 1;
  }
  else
  {
    istart = nq - 1;
    istop = -1;
    istep = -1;
  }

  for (int i = istart; i != istop; i += istep)
  {
    Substructure_Results sresults;

    int nhits = _queries[i]->substructure_search(target, sresults);

    if (_show_all_query_matches && nhits)
      cerr << "Query " << i << ' ' << _queries[i]->comment() << " " << nhits << " matches\n";
    else if (_show_all_query_non_matches && 0 == nhits)
      cerr << "Query " << i << ' ' << _queries[i]->comment() << " only matched " << sresults.max_query_atoms_matched_in_search() << endl;

    for (int j = 0; j < nhits; j++)
    {
      Set_of_Atoms * e =  const_cast<Set_of_Atoms *>(sresults.embedding(j));   // loss of const OK

      if (verbose > 1)
        cerr << " hit " << j << ' ' << *e << endl;

      if (! _queries[i]->all_atoms_unmatched(*e, already_hit))
        continue;

      _queries[i]->mark_matched_atoms(*e, already_hit, i + 1);
      _queries[i]->place_isotopes(*e, isotope);

      _queries[i]->increment_results(results + _result_column[0]);
#ifdef DEBUG_PROCESS
      cerr << i << " after nit " << j << " to '" << _queries[i]->comment() << "'\n";
      for (auto k = 0; k < _nmodels; ++k)
      {
        cerr << " " << k << " " << results[k] << endl;
      }
#endif
    }
  }

  return 1;
}

void
Queries_and_Additive_Models::convert_to_int_count(const double * v, 
                                                  int * b) const
{
//cerr << "Queries_and_Additive_Models:::convert_to_int_count:writing " << _nmodels << " models\n";

  b += _result_column[0];    // shift to where our results start

  for (auto i = 0; i < _nmodels; ++i)
  {
//  cerr << "line " << __LINE__ << " i = " << i << " cmp " << _global_min[i] << " max " << _global_max[i] << endl;

    if (v[i] <= _global_min[i])
      b[i] = 1;
    else if (v[i] >= _global_max[i])
      b[i] = bit_count_dynamic_range + 1;
    else
    {
      b[i] = static_cast<int>((v[i] - _global_min[i]) / (_global_max[i] - _global_min[i]) * bit_count_dynamic_range + 1.4999);
    }
  }

  return;
}

/*static int
run_a_set_of_queries(const resizable_array_p<Abraham_Substructure_Query> & q,
                     Molecule_to_Match & target,
                     Abraham_Parameters & ap,
                     int * already_hit,
                     int * isotope)
{
  int nq = q.number_elements();

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    const Abraham_Parameters & api = *(q[i]);

    int nhits = q[i]->substructure_search(target, sresults);

    if (0 == nhits && verbose > 2)
      cerr << "Query " << q[i]->comment() << " only matched " << sresults.max_query_atoms_matched_in_search() << endl;

    if (verbose > 1 && nhits)
      cerr << ' ' << nhits << " hits to query " << i << " '" << q[i]->comment() << "'\n";

    for (int j = 0; j < nhits; j++)
    {
      Set_of_Atoms * e =  const_cast<Set_of_Atoms *>(sresults.embedding(j));   // loss of const OK

      if (verbose > 1)
        cerr << " hit " << j << ' ' << *e << endl;

      if (! q[i]->all_atoms_unmatched(*e, already_hit))
        continue;

      q[i]->mark_matched_atoms(*e, already_hit, i + 1);
      q[i]->place_isotopes(*e, isotope);

      ap += api;
//    cerr << "After '" << q[i]->comment() << "' value R2 " << ap.r2() << " alpha2H " << ap.alpha2h() << ", pi2H " << ap.pi2h() << endl;
    }
  }

  return 1;
}*/

static int
compute_number_results (const Queries_and_Additive_Models * models,
                        const int nmodels)
{
  int rc = 0;
  for (auto i = 0; i < nmodels; ++i)
  {
    rc += models[i].nmodels();
  }

  return rc;
}

/*
  The results arrays used each time. They are fixed
  size and do not need to be allocated/deallocated
*/

static int * int_output_for_fingerprints = nullptr;
static double * v = nullptr;

static int
abraham (Molecule & m,
         int * already_hit,
         int * isotope,
         Queries_and_Additive_Models * models,
         const int number_query_sets,
         const int nmodels,
         IWString_and_File_Descriptor & output)
{
  if (verbose > 1)
    cerr << "Processing '" << m.name() << "'\n";

  Molecule_to_Match target(&m);

  int failures = 0;

  for (auto i = 0; i < number_query_sets; ++i)
  {
    if (! models[i].process(m, target, already_hit, isotope, v))
      failures++;

    if (nullptr != int_output_for_fingerprints)
      models[i].convert_to_int_count(v, int_output_for_fingerprints);
  }

  v[ap_vx] = mcgowan(m);

  if (nullptr != int_output_for_fingerprints)
    convert_to_int_count(v[ap_vx], int_output_for_fingerprints[ap_vx], min_Vx, max_Vx);

//cerr << "Indices " << ap_r2 << ' ' << ap_pi2h << ' ' << ap_alpha2h << ' ' << ap_beta2o << " ndx " << ndx << endl;

//cerr << "Setting logp at " << ap_logp << endl;
  v[ap_logp] = abraham_logp(v[ap_r2], v[ap_pi2h], v[ap_alpha2h], v[ap_beta2o], v[ap_vx]);

  if (nullptr != int_output_for_fingerprints)
    convert_to_int_count(v[ap_logp], int_output_for_fingerprints[ap_logp], min_logp, max_logp);

// Now we need to combine the first order Constantinou Gani models

  if (constantinou_gani_models_present)
  {
    v[ap_cg[1]] += v[ap_cg[3]];
    v[ap_cg[2]] += v[ap_cg[4]];

//  and shift down the other models

    v[ap_cg[3]] = v[ap_cg[4] + 1];
    v[ap_cg[4] + 1] = v[ap_cg[4] + 2];
  }

  if (verbose)
  {
    for (int i = 0; i < nmodels; ++i)
    {
      acc[i].extra(v[i]);
    }
  }

  if (nullptr != int_output_for_fingerprints)
    return write_fingerprint(m, fingerprint_tag, int_output_for_fingerprints, nmodels, output);

  append_first_token_of_name(m.name(), output);

  for (auto i = 0; i < nmodels; ++i)
  {
    output << output_separator << static_cast<float>(v[i]);
  }

  output << '\n';

  return output.good();
}

static void
preprocess(Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

static int
abraham (Molecule & m, 
         Queries_and_Additive_Models * models,
         const int number_query_sets,
         const int nmodels,
         IWString_and_File_Descriptor & output)
{
  molecules_read++;

  preprocess(m);

  if (make_implicit_hydrogens_implicit)
    m.make_implicit_hydrogens_explicit();

  const auto matoms = m.natoms();

  if (0 == matoms)
  {
    cerr << "Ignoring empty molecule\n";
    return 1;
  }

  int * tmp = new_int(matoms + matoms); std::unique_ptr<int[]> free_tmp(tmp);   // one for already_hit and one for isotopes

  const auto rc = abraham(m, tmp, tmp + matoms, models, number_query_sets, nmodels, output);

  return rc;
}

static int
do_discern_query_overlaps(Molecule & m,
                          Queries_and_Additive_Models * models, 
                          const int number_query_sets)
{
  for (int i = 0; i < number_query_sets; ++i)
  {
    models[i].examine_query_overlap_behaviour(m);
  }

  return 1;
}

static int
abraham (data_source_and_type<Molecule> & input,
         Queries_and_Additive_Models * models,
         const int number_query_sets,
         const int nmodels,
         IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()) && output.good())
  {
    std::unique_ptr<Molecule> free_m(m);

#ifdef TEST_MCGOWAN
    output << mcgowan(*m) << "\n";
    output.write_if_buffer_holds_more_than(8192);
    continue;
#endif

    if (discern_query_overlaps)
    {
      do_discern_query_overlaps(*m, models, number_query_sets);
      continue;
    }

    if (! abraham(*m, models, number_query_sets, nmodels, output))
      return 0;

    if (verbose)
      output.flush();
    else 
      output.write_if_buffer_holds_more_than(8192);
  }

  return output.good();
}

static int
abraham_record (const_IWSubstring buffer,     // local copy
                Queries_and_Additive_Models * models,
                const int number_query_sets,
                const int nmodels,
                IWString_and_File_Descriptor & output)
{
  buffer.remove_leading_chars(smiles_tag.length());
  assert (buffer.ends_with('>'));
  buffer.chop();

  Molecule m;

  if (! m.build_from_smiles(buffer))
  {
    cerr << "Invalid smiles '" << buffer << "'\n";
    return 0;
  }

  return abraham(m, models, number_query_sets, nmodels, output);
}

static int
abraham (iwstring_data_source & input,
         Queries_and_Additive_Models * models,
         const int number_query_sets,
         const int nmodels,
         IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! buffer.starts_with(smiles_tag))
      output << buffer << '\n';
    else if (! abraham_record(buffer, models, number_query_sets, nmodels, output))
    {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

static int
abraham (const char * fname,
         Queries_and_Additive_Models * models,
         const int number_query_sets,
         const int nmodels,
         IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return abraham(input, models, number_query_sets, nmodels, output);
}


static int
abraham (const char * fname, int input_type, 
         Queries_and_Additive_Models * models,
         const int number_query_sets,
         const int nmodels,
         IWString_and_File_Descriptor & output)
{
  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return abraham(input, models, number_query_sets, nmodels, output);
}

static int
write_header_records (const Queries_and_Additive_Models * models,
                      const int number_query_sets,
                      const int nmodels,
                      IWString_and_File_Descriptor & output)
{
  output << "Name";

  for (auto i = 0; i < number_query_sets; ++i)
  {
    models[i].append_model_names(descriptor_prefix, output);
  }

  output << output_separator << descriptor_prefix << "Vx " << descriptor_prefix << "Logp";

  output << '\n';

  return output.good();
}

/*
  Once upon a time, the make all hydrogens explicit was -H. I changed it so -h does
  that. Leave the -H option there, but it is now ignored
*/

static int
abraham (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:o:g:F:L:lhJ:rP:G:fp:d:C:y:");

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
      return 17;
    }
  }
  else
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot process element specification option(s), -E option\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      return 62;
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;
    if (verbose)
      cerr << "Will reduce multi-fragment molecules to largest fragment\n";
  }

  if (cl.option_present('G'))
  {
    int i = 0;
    const_IWSubstring g;
    while (cl.value('G', g, i++))
    {
      if ("blki" == g)
      {
        set_aromatic_bonds_lose_kekule_identity(1);
        if (verbose)
          cerr << "Will use strict Daylight rules for aliphatic atoms\n";
      }
      else if ("qndn" == g)
      {
        use_query_name_as_descriptor_name = 1;
        if (verbose)
          cerr << "Will use query names as descriptor names\n";
      }
      else if ("tovl" == g)
      {
        discern_query_overlaps = 1;
        if (verbose)
          cerr << "Will study query overlaps\n";
      }
      else if ("ncpm" == g)
      {
        check_for_previously_matched_atoms = 0;
        if (verbose)
          cerr << "Will NOT check for previously matched atoms\n";
      }
      else if ("ql2f" == g)
      {
        process_queries_first_to_last = 0;
        if (verbose)
          cerr << "Will process queries last to first\n";
      }
      else if ("help" == g)
      {
        display_dash_G_qualifiers(cerr);
        return 0;
      }
      else
      {
        cerr << "Unrecognised -G qualifier '" << g << "'\n";
        return 5;
      }
    }
  }

  if (cl.option_present('r'))
  {
    set_aromatic_bonds_lose_kekule_identity(1);
    if (verbose)
      cerr << "Will use strict Daylight rules for aliphatic atoms - use '-G blki' instead\n";
  }

  if (cl.option_present('h'))
  {
    make_implicit_hydrogens_implicit = 1;
    if (verbose)
      cerr << "Implicit hydrogens will be made explicit\n";
  }

  std::fill_n(mcgowan_coefficient, HIGHEST_ATOMIC_NUMBER + 1, 0.0f);

  mcgowan_coefficient[1] = 2.15;
  mcgowan_coefficient[5] = 18.32;
  mcgowan_coefficient[6] = 16.35;
  mcgowan_coefficient[7] = 14.39;
  mcgowan_coefficient[8] = 12.43;
  mcgowan_coefficient[9] = 10.48;
  mcgowan_coefficient[14] = 26.83;     // Si
  mcgowan_coefficient[15] = 24.87;     // P
  mcgowan_coefficient[16] = 22.91;     // S
  mcgowan_coefficient[17] = 20.95;     // Cl
  mcgowan_coefficient[32] = 31.02;     // Ge
  mcgowan_coefficient[33] = 29.42;     // As
  mcgowan_coefficient[34] = 27.81;     // Se
  mcgowan_coefficient[35] = 26.21;     // Br
  mcgowan_coefficient[50] = 39.35;     // Sn
  mcgowan_coefficient[51] = 37.73;     // Sb
  mcgowan_coefficient[52] = 36.14;     // Te
  mcgowan_coefficient[53] = 34.53;     // I

  if (! cl.option_present('F') || ! cl.option_present('P'))
  {
    cerr << "Needs at least two files, -F and -P\n";
    usage(19);
  }

// Each set of queries can contain multiple models.

  int number_query_sets = 2;
  if (cl.option_present('C'))
    number_query_sets += 4;

  Queries_and_Additive_Models * models = new Queries_and_Additive_Models[number_query_sets]; std::unique_ptr<Queries_and_Additive_Models[]> free_models(models);

  int next_result_column = 0;
  int next_query_set = 0;

  if (cl.option_present('F'))
  {
    const char * f = cl.option_value('F');

    if (! models[ABRAHAM_INDEX].build(f))
    {
      cerr << "Cannot read Abraham models control file '" << f << "'\n";
      return 73;
    }

    if (verbose)
      cerr << "Read " << models[ABRAHAM_INDEX].number_queries() << " queries for Abraham models\n";

    models[ABRAHAM_INDEX].set_result_column(next_result_column);

    ap_r2 = models[ABRAHAM_INDEX].result_column("r2");
    ap_pi2h = models[ABRAHAM_INDEX].result_column("pi2h");
    ap_beta2h = models[ABRAHAM_INDEX].result_column("beta2h");
    ap_beta2o = models[ABRAHAM_INDEX].result_column("beta2o");

    next_query_set++;
  }

  if (cl.option_present('P'))
  {
    const char * p = cl.option_value('P');

    if (! models[ALPHA2A_INDEX].build(p))
    {
      cerr << "Cannot read alpha2H queries '" << p << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Read " << models[ALPHA2A_INDEX].number_queries() << " Alpha2A queries from '" << p << "'\n";

    models[ALPHA2A_INDEX].set_result_column(next_result_column);

    ap_alpha2h = models[ALPHA2A_INDEX].result_column("alpha2h");

    next_query_set++;
  }

  if (cl.option_present('C'))
  {
    if (cl.option_present('J'))
    {
      cerr << "Sorry, cannot run Constantinou Gani models and generate fingerprints, see Ian\n";
      return 1;
    }

    IWString table[5];

    const_IWSubstring c;
    for (auto i = 0; cl.value('C', c, i); ++i)
    {
      if (c.starts_with("table1="))
        table[1]=c;
      else if (c.starts_with("table2="))
        table[2]=c;
      else if (c.starts_with("table3="))
        table[3]=c;
      else if (c.starts_with("table4="))
        table[4]=c;
      else
      {
        cerr << "Unrecognised Constantinou Gani file directive '" << c << "'\n";
        return 1;
      }
    }

    for (auto i = 1; i <= 4; ++i)
    {
      table[i].remove_leading_chars(7);

      if (0 == table[i].length())
      {
        cerr << "No specification for Constantinou Gani table " << i << endl;
        return 2;
      }

      ap_cg[i] = next_query_set;

      next_query_set++;

      if (! models[ap_cg[i]].build(table[i]))
      {
        cerr << "Cannot read Constantinou Gani Table 1 queries '" << table[i] << "'\n";
        return 0;
      }

      if (verbose)
        cerr << "Constantinou Gani Table " << i << " read " << models[CG1_INDEX + i - 1].number_queries() << " from '" << table[i] << "'\n";

      models[ap_cg[i]].set_result_column(next_result_column);
    }

    constantinou_gani_models_present = 1;
  }

  ap_vx = next_result_column;    // McGowan
  next_result_column++;
  ap_logp = next_result_column;    // clogp
  next_result_column++; 

  const int nmodels = next_result_column;

// We compute all the results that come from Queries_and_Additive_Models's, then the extra ones

  int tmp = compute_number_results(models, number_query_sets);   // do not change till the 'y' option is done

  if (cl.option_present('y'))
  {
    int show_all_query_matches = 0;
    int show_all_query_non_matches = 0;

    const_IWSubstring y;
    for (auto i = 0; cl.value('y', y, i); ++i)
    {
      if ('m' == y)
        show_all_query_matches = 1;
      else if ('n' == y)
        show_all_query_non_matches = 1;
      else
      {
        cerr << "Unrecognised -y qualifier '" << y << "'\n";
        usage(2);
      }
    }

    for (auto i = 0; i < tmp; ++i)
    {
      if (show_all_query_matches)
        models[i].set_show_all_query_matches(1);
      if (show_all_query_non_matches)
        models[i].set_show_all_query_non_matches(1);
    }
  }

  ap_vx = tmp;

  tmp++;

  ap_logp = tmp;

  tmp++;

  const auto nresults = tmp;

//cerr << "Have " << nresults << " results, nmodels " << nmodels << endl;

  v = new double[nresults]; std::unique_ptr<double[]> free_v(v);

  if (verbose)
    cerr << "Output will have " << nresults << " results\n";

  if (cl.option_present('J'))
  {
    cl.value('J', fingerprint_tag);

    if (! fingerprint_tag.ends_with('<'))
      fingerprint_tag << '<';

    if (verbose)
      cerr << "Results written as fingerprint with tag '" << fingerprint_tag << "'\n";

    if (cl.option_present('f'))
    {
      function_as_tdt_filter = 1;

      if (verbose)
        cerr << "Will function as a TDT filter\n";
    }

    if (cl.option_present('p'))
    {
      if (! cl.value('p', bit_replicates) || bit_replicates < 1)
      {
        cerr << "The number of bit replicates (-p) option must be a whole +ve number\n";
        usage(2);
      }

      if (verbose)
        cerr << "Will produce " << bit_replicates << " bit replicates\n";
    }

    if (cl.option_present('d'))
    {
      if (! cl.value('d', bit_count_dynamic_range) || bit_count_dynamic_range <= 1.0f)
      {
        cerr << "The bit count dynamic range (-d) must be a +ve value\n";
        usage(2);
      }

      if (verbose)
        cerr << "Each replicate will scale to " << bit_count_dynamic_range << endl;
    }

    int_output_for_fingerprints = new int[nresults];
  }

  int input_type = 0;
  if (function_as_tdt_filter)
    ;
  else if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('L'))
  {
    if (! cl.option_present('o'))
      stream_for_labeled_atoms.add_output_type(SMI);
    else if (! stream_for_labeled_atoms.determine_output_types(cl))
    {
      cerr << "Cannot determine output types for labeled atom molecules\n";
      return 0;
    }

    const_IWSubstring l;
    cl.value('L', l);

    if (! stream_for_labeled_atoms.new_stem(l, 1))
    {
      cerr << "Cannot open -L file(s) with stem '" << l << "'\n";
      return 53;
    }

    if (verbose)
      cerr << "Molecules with labeled atoms written to '" << l << "'\n";
  }

  if (verbose)
    acc = new Accumulator<double>[nmodels];

  IWString_and_File_Descriptor output(1);

  if (0 == fingerprint_tag.length())
    write_header_records(models, number_query_sets, nmodels, output);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (function_as_tdt_filter)
    {
      if (! abraham(cl[i], models, number_query_sets, nmodels, output))
      {
        rc = i + 1;
        break;
      }
    }
    else if (! abraham(cl[i], input_type, models, number_query_sets, nmodels, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    for (auto i = 0; i < number_query_sets; ++i)
    {
      models[i].report(acc, cerr);
    }
  }

  if (discern_query_overlaps)
  {
    for (int i = 0; i < number_query_sets; ++i)
    {
      models[i].report_overlaps(std::cerr);
    }
  }

  if (nullptr != int_output_for_fingerprints)
    delete [] int_output_for_fingerprints;

  if (nullptr != acc)
    delete [] acc;

  return rc;
}


int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = abraham(argc, argv);

  return rc;
}
