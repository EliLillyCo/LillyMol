#include <stdlib.h>
#include <math.h>
#include <unordered_map>
#include <memory>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <limits>

#define IWQSORT_FO_IMPLEMENTATION
#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "iw_tdt.h"
#include "accumulator.h"
#include "report_progress.h"
#include "iw_stl_hash_map.h"
#include "iw_stl_hash_set.h"
#include "iwqsort.h"

#include "gfp.h"
#include "bit_subset_v2.h"

#define HEADER_RECORD "# written by gfp_naive_bayesian"
#define COUNT_FIXED "count_fixed"
#define COUNT_SPARSE "count_sparse"
#define COUNT_ACTIVE "count_active"
#define COUNT_INACTIVE "count_inactive"
#define BAYES_METHOD  "bayes_method"
#define SUM_FREQUENCY_ACTIVE "sum_frequency_active"
#define SUM_FREQUENCY_INACTIVE "sum_frequency_inactive"
#define TOTAL_BITS "total_bits"

#define GFP_NB_HEADER "GFP_NB"
#define GFP_NB_ACT "GFP_NB_ACTIVITY"
#define GFP_NB_CLASSES "CLASSES"
#define GFP_NB_CLASS "CLASS"
#define GFP_NB_BP "BP"
#define GFP_NB_FX "FX"
#define GFP_NB_SP "SP"

static char output_separator = ' ';

const char * prog_name = NULL;
static int verbose = 0;

static int strip_leading_zeros = 0;

static int activity_cutoff_set = 0;
static float activity_cutoff = 0.0;

static int activity_column = 1;   // by default identifier in column 0

#define BAYES_METHOD_BINARY 1
#define BAYES_METHOD_MULTINOMIAL 2

static int bayes_method = BAYES_METHOD_MULTINOMIAL;

static int min_support_level = 0;

static int remove_correlated_bits = 0;

static int sum_starts_with_log_class_probability = 1;

static int write_probabilities_in_11_range = 0;

static IWString identifier_tag ("PCN<");

static int brief_output = 0;

static double all_scores_above_threshold = std::numeric_limits<double>::max();

static IWString_and_File_Descriptor stream_for_bit_contributions;

/*
  Let's do a very kludgey implementation of that
*/

static int distinguish_bits_by_count = 0;

/*
  April 2015. 
  Make the implementation of the formula that depends on the number of molecules
  with a feature a switchable thing
*/

static int denominator_is_molecules_with_feature = 0;

/*
  Want to enable some debugging
*/

#define BIT_EXAMINATION_ENABLED
#ifdef BIT_EXAMINATION_ENABLED
static std::unordered_set<unsigned int> bits_to_examine;

static int examine_all_bits = 0;
#endif

//#define DEBUG_IWNB

static unsigned int
bit_number_adjusted_for_count (unsigned int b,
                               int & c)
{
  if (1 == c)
    return b;

  if (distinguish_bits_by_count > 1)
    c = c / distinguish_bits_by_count + 1;

  unsigned int rc;

  if (b < 16777216)
    rc = b + c * 97;
  else
    rc = b - c * 103;

#ifdef DEBUG_IWNB
  cerr << "Bit " << b << " count " << c <<  " converted to " << rc << endl;
#endif

  c = 1;

  return rc;
}

static int
read_activity_data_record (const const_IWSubstring & buffer,
                           IW_STL_Hash_Map_int & activity,
                           IW_STL_Hash_Map_int & string_to_class_number,
                           resizable_array_p<IWString> & class_names)
{
  IWString id, token;
  int i = 0;

  char separator = ' ';

  if (buffer.nextword(id, i) && buffer.nextword(token, i))
    ;
  else if (buffer.nextword(id, i, '\t') && buffer.nextword(token, i, '\t'))
    separator = '\t';
  else if (buffer.nextword(id, i, ',') && buffer.nextword(token, i, ','))
    separator = ',';
  else
  {
    cerr << "Cannot separate into identifier and activity '" << buffer << "'\n";
    return 0;
  }

  if (strip_leading_zeros)
    id.remove_leading_chars('0');

  if (1 != activity_column)
  {
    if (! buffer.word(activity_column, token, separator))
    {
      cerr << "Cannot extract column '" << (activity_column + 1) << " from record\n";
      return 0;
    }
  }


  if (activity_cutoff_set)
  {
    float a;
    if (! token.numeric_value(a))
    {
      cerr << "Invalid activity value '" << token << "'\n";
      return 0;
    }
    if (a < activity_cutoff)
      activity[id] = 0;
    else
      activity[id] = 1;
  }
  else       // textual interpretation
  {
    auto f = string_to_class_number.find(token);
    if (f == string_to_class_number.end())
    {
      int c = string_to_class_number.size();
      string_to_class_number[token] = c;
      class_names.add (new IWString(token));
      activity[id] = c;
    }
    else
      activity[id] = (*f).second;
  }

//cerr << "for id '" << id << "' value '" << activity[id] << "', token '" << token << "'\n";

  return 1;
}

static int
read_activity_data (iwstring_data_source & input,
                    IW_STL_Hash_Map_int & activity,
                    resizable_array_p<IWString> & class_names,
                    IWString & activity_name)
{
  IW_STL_Hash_Map_int string_to_class_number;

  if (! input.next_record(activity_name))
  {
    cerr << "read_activity_data:cannot read header record\n";
    return 0;
  }

  activity_name.remove_leading_words(1);

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;

    if (! read_activity_data_record(buffer, activity, string_to_class_number, class_names))
    {
      cerr << "Cannot read activity data, line " << input.lines_read() << endl;
      cerr << "'" << buffer << "'\n";
      return 0;
    }
  }

  return activity.size();
}

static int
read_activity_data (const const_IWSubstring & fname,
                    IW_STL_Hash_Map_int & activity,
                    resizable_array_p<IWString> & class_names,
                    IWString & activity_name)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);
  input.set_translate_tabs(1);

  return read_activity_data (input, activity, class_names, activity_name);
}

typedef std::pair<int,unsigned int> FP_And_Bit;

class FP_And_Bit_Hasher
{
  public :
    size_t operator () (const FP_And_Bit & f) const { return f.second;}
};

/*
  when we read fingerprints from a file, we need to store that info somewhere so
  we can check to see if the fingerprints read are the same
*/

class Fingerprint_Info_Read_From_Model
{
  private:
    std::unordered_map<int, IWString> _fx;
    std::unordered_map<int, IWString> _sp;

//  private functions

    int _add (const int, const const_IWSubstring &, std::unordered_map<int, IWString> &);
    int _matches (const int ndx, const IWString & tag, const std::unordered_map<int, IWString> & ndx_to_tag) const;

  public:
    int add (const const_IWSubstring &);

    int number_fingerprints () const { return _fx.size() + _sp.size();}

    int check_fingerprint_validity () const;
};

int
Fingerprint_Info_Read_From_Model::add (const const_IWSubstring & buffer)
{
  if (3 != buffer.nwords())
  {
    cerr << "Fingerprint_Info_Read_From_Model::add:invalid fingerprint specification '" << buffer << "'\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring fptype;

  buffer.nextword(fptype, i);

  const_IWSubstring token;
  buffer.nextword(token, i);

  int ndx;
  if (! token.numeric_value(ndx) || ndx < 0)
  {
    cerr << "Fingerprint_Info_Read_From_Model::add:invalid fingerprint number '" << buffer << "'\n";
    return 0;
  }

  buffer.nextword(token, i);

  int rc = 1;

  if (GFP_NB_SP == fptype)
    rc = _add(ndx, token, _sp);
  else if (GFP_NB_FX == fptype)
    rc = _add(ndx, token, _fx);
  else
  {
    cerr << "Fingerprint_Info_Read_From_Model::add:unrecognised fingerprint type '" << buffer << "'\n";
    return 0;
  }

  if (rc)
    return rc;

  cerr << "Fingerprint_Info_Read_From_Model::add:cannot process '" << buffer << "'\n";
  return 0;
}

int
Fingerprint_Info_Read_From_Model::_add (const int ndx,
                                        const const_IWSubstring & tag,
                                        std::unordered_map<int, IWString> & ndx_to_tag)
{
  const auto f = ndx_to_tag.find(ndx);
  if (f != ndx_to_tag.end())
  {
    cerr << "Fingerprint_Info_Read_From_Model::_add:duplicate index " << ndx << endl;
    return 0;
  }

  ndx_to_tag[ndx] = tag;

  return 1;
}

int
Fingerprint_Info_Read_From_Model::_matches (const int ndx,
                                            const IWString & tag,
                                            const std::unordered_map<int, IWString> & ndx_to_tag) const
{
  const auto f = ndx_to_tag.find(ndx);
  if (f == ndx_to_tag.end())
  {
    cerr << "Fingerprint_Info_Read_From_Model::_matches:no fingerprint " << ndx << endl;
    return 0;
  }

  return tag == f->second;
}

int
Fingerprint_Info_Read_From_Model::check_fingerprint_validity () const
{
  int rc = 1;

  for (int i = 0; i < ::number_fingerprints(); ++i)
  {
    if (! _matches(i, fixed_fingerprint_tag(i), _fx))
    {
      cerr << "Fingerprint_Info_Read_From_Model::check_fingerprint_validity:fixed size fingerprint mismatch, i = " << i << " tag " << fixed_fingerprint_tag(i) << endl;
      rc = 0;
    }
  }

  for (int i = 0; i < ::number_sparse_fingerprints(); ++i)
  {
    if (! _matches(i, sparse_fingerprint_tag(i), _sp))
    {
      cerr << "Fingerprint_Info_Read_From_Model::check_fingerprint_validity:sparse fingerprint mismatch, i = " << i << " tag " << sparse_fingerprint_tag(i) << endl;
      rc = 0;
    }
  }

  if (rc)
    return rc;

  cerr << "Fingerprint_Info_Read_From_Model::check_fingerprint_validity:inconsistent fingerprints\n";

  for (auto i : _fx)
  {
    cerr << i.first << " is " << fixed_fingerprint_tag(i.first) << " model expected " << i.second << endl;
  }

  for (auto i : _sp)
  {
    cerr << i.first << " is " << sparse_fingerprint_tag(i.first) << " model expected " << i.second << endl;
  }

  return 0;
}

//typedef unordered_map<FP_And_Bit, Frequency, FP_And_Bit_Hasher> Bit_Frequency;

//typedef unordered_map<FP_And_Bit, double, FP_And_Bit_Hasher> Bit_Weight;

/*
  For each class, we need a list of the bits and the associated weights.
  Within each fingerprint there are a given number of fingerprints;
*/


static double laplace_correction_factor=1.0;

typedef std::unordered_map<FP_And_Bit, int,    FP_And_Bit_Hasher> Bit_Frequency;
typedef std::unordered_map<FP_And_Bit, double, FP_And_Bit_Hasher> Bit_Weight;

class Bayesian_Classification
{
  private:
    int _nclasses;
    int _n;
    int * _items_in_class;
    double * _baseline_class_probability;

    Bit_Frequency * _bit_frequency;   // for each class
    Bit_Weight    * _bit_weight;   // for each class

//  private functions

    int _allocate_arrays ();
    int _score_while_writing_bit_contributions (const IW_General_Fingerprint & fp, double * rc, IWString_and_File_Descriptor & output) const;
    int _write_bit_weight(const int ndx, const IWString * class_label, IWString_and_File_Descriptor & output) const;
    int _read_bit_weight(Bit_Weight &, iwstring_data_source & input);

  public:
    Bayesian_Classification();
    ~Bayesian_Classification();

    int initialise (const IW_STL_Hash_Map_int & activity);
    int initialise (int nc);

    int nclasses () const { return _nclasses;}

    int determine_bit_counts (IW_General_Fingerprint & fp, const int c);
    int determine_bit_counts (IW_General_Fingerprint & fp, const IW_STL_Hash_Map_int & activity);
    int determine_bit_counts (IW_TDT & tdt, const IW_STL_Hash_Map_int & activity);
    int determine_bit_counts (iwstring_data_source & input, const IW_STL_Hash_Map_int & activity);
    int determine_bit_counts (const char * fname, const IW_STL_Hash_Map_int & activity);

    int determine_bit_counts (const char * fname, const IW_STL_Hash_Map_int & activity,
                                               const IW_STL_Hash_Set & training_set_identifiers,
                                               resizable_array_p<IW_General_Fingerprint> & test_set_fingerprints);
    int determine_bit_counts (iwstring_data_source &, const IW_STL_Hash_Map_int & activity,
                                               const IW_STL_Hash_Set & training_set_identifiers,
                                               resizable_array_p<IW_General_Fingerprint> & test_set_fingerprints);
    int determine_bit_counts (IW_TDT & tdt, const IW_STL_Hash_Map_int & activity,
                                               const IW_STL_Hash_Set & training_set_identifiers,
                                               resizable_array_p<IW_General_Fingerprint> & test_set_fingerprints);

    int determine_bit_counts (const IW_General_Fingerprint * fp, const int nfingerprints,
                                               const IW_STL_Hash_Map_int & activity,
                                               const int * training_set_identifiers);
    int initialise_bit_weights ();

    int remove_bits_below_support_level (int support);

    int score (const IW_General_Fingerprint & fp, double *) const;
    int score (IW_TDT & tdt, double * rc) const;
    int score (iwstring_data_source &, double * rc) const;

    int do_write (const IWString * class_label, const IWString & activity_name, const char * fname) const;
    int do_write (const IWString * class_label, const IWString & activity_name, IWString_and_File_Descriptor & output) const;

    int build (const char * fname, IWString * & cname, IWString & activity_name, Fingerprint_Info_Read_From_Model & firfm);
    int build (iwstring_data_source & input, IWString * & cname, IWString & activity_name, Fingerprint_Info_Read_From_Model & firfm);
};

Bayesian_Classification::Bayesian_Classification()
{
  _nclasses = 0;
  _n = 0;
  _items_in_class = nullptr;
  _baseline_class_probability = nullptr;
  _bit_frequency = nullptr;
  
  return;
}

Bayesian_Classification::~Bayesian_Classification ()
{
  if (nullptr != _items_in_class)
  {
    delete [] _items_in_class;
    delete [] _baseline_class_probability;
    delete [] _bit_frequency;
  }

  return;
}

int
Bayesian_Classification::_allocate_arrays ()
{
  _items_in_class = new_int(_nclasses);
  _baseline_class_probability = new double[_nclasses];
  _bit_frequency = new Bit_Frequency[_nclasses];
  _bit_weight = new Bit_Weight[_nclasses];

  return 1;
}

int
Bayesian_Classification::initialise (const IW_STL_Hash_Map_int & activity)
{
  if (0 == activity.size())
  {
    cerr << "Bayesian_Classification::initialise:no data\n";
    return 0;
  }

  _nclasses = 0;

  for (auto i = activity.begin(); i != activity.end(); ++i)
  {
    if ((*i).second > _nclasses)
      _nclasses = (*i).second;
  }

  if (0 == _nclasses)
  {
    cerr << "Bayesian_Classification::initialise:no classes in " << activity.size() << " data values\n";
    return 0;
  }

  _nclasses++;

  return _allocate_arrays();
}

int
Bayesian_Classification::initialise (int nc)
{
  _nclasses = nc;

  return _allocate_arrays ();
}

static int
increment_bit (const int fpi,
               const unsigned int b,
               Bit_Frequency & bit_frequency,
               int c=1)
{
  FP_And_Bit fab(fpi,b);
  
  Bit_Frequency::iterator f = bit_frequency.find(fab);

  if (f == bit_frequency.end())
    bit_frequency[fab] = c;
  else
    (*f).second += c;
  
  return 1;
}

static int
count_occurrences (IW_General_Fingerprint & fp,
                   Bit_Frequency & bit_frequency)
{
  const int number_fixed = number_fingerprints();
  const int number_sparse = number_sparse_fingerprints();

  int fpi=0;
  
  while (fpi<number_fixed)
  {
    const IWDYFP & dyfp = fp[fpi];
    int b;
    int i = 0;

    while ((b = dyfp.next_on_bit(i)) >= 0)
    {
      increment_bit (fpi, b, bit_frequency);
    }
    fpi++;
  }
  
  while (fpi<(number_fixed+number_sparse))
  {
    const Sparse_Fingerprint & sfp = fp.sparse_fingerprint(fpi-number_fixed);
          
    int i = 0;
    unsigned int b;
    int c;
  
    while (sfp.next_bit_set(i, b, c))
    {
      if (distinguish_bits_by_count)
        b = bit_number_adjusted_for_count(b, c);

      if (BAYES_METHOD_MULTINOMIAL == bayes_method)
        increment_bit(fpi, b, bit_frequency, c);
      else
        increment_bit(fpi, b, bit_frequency, c);
    }
    fpi++;
  }

  return 1;
}

int
Bayesian_Classification::determine_bit_counts (IW_TDT & tdt,
                                               const IW_STL_Hash_Map_int & activity)
{
  IW_General_Fingerprint fp;

  int fatal;

  if (! fp.construct_from_tdt (tdt, fatal))
    return 0;

  return determine_bit_counts(fp, activity);
}

static void
standardise_identifier (IWString & id)
{
  if (strip_leading_zeros && id.starts_with('0'))
    id.remove_leading_chars('0');

  id.truncate_at_first(' ');

  return;
}

int 
Bayesian_Classification::determine_bit_counts (IW_General_Fingerprint & fp,
                                               const IW_STL_Hash_Map_int & activity)
{
  IWString id = fp.id();

  standardise_identifier(id);

  const auto f = activity.find(id);

  if (f == activity.end())
  {
    cerr << "No activity data for '" << fp.id() << "'\n";
    return 0;
  }

  _items_in_class[f->second]++;

  _n++;

  return count_occurrences(fp, _bit_frequency[f->second]);
}

int
Bayesian_Classification::determine_bit_counts (IW_General_Fingerprint & fp,
                                               const int c)
{
  _items_in_class[c]++;
  _n++;
  return count_occurrences(fp, _bit_frequency[c]);
}

int
Bayesian_Classification::determine_bit_counts (iwstring_data_source & input,
                                               const IW_STL_Hash_Map_int & activity)
{
  IW_TDT tdt;
  while (tdt.next(input))
  {
    if (! determine_bit_counts (tdt, activity))
    {
      cerr << "Fatal error processing '" << tdt << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Bayesian_Classification::determine_bit_counts (const char * fname,
                                               const IW_STL_Hash_Map_int & activity)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Bayesian_Classification::determine_bits_set:cannot open fingerprint file '" << fname << "'\n";
    return 0;
  }

  return determine_bit_counts (input, activity);
}

int
Bayesian_Classification::determine_bit_counts (IW_TDT & tdt,
                                               const IW_STL_Hash_Map_int & activity,
                                               const IW_STL_Hash_Set & training_set_identifiers,
                                               resizable_array_p<IW_General_Fingerprint> & test_set_fingerprints)
{
  IW_General_Fingerprint * fp = new IW_General_Fingerprint;
  int fatal;

  if (! fp->construct_from_tdt(tdt, fatal))
  {
    delete fp;
    return 0;
  }

  IWString id = fp->id();
  standardise_identifier(id);

  if (! training_set_identifiers.contains(id))
  {
    fp->set_id(id);

    test_set_fingerprints.add(fp);
    return 1;
  }

  std::unique_ptr<IW_General_Fingerprint> free_fp(fp);

  return determine_bit_counts(*fp, activity);
}

int
Bayesian_Classification::determine_bit_counts (iwstring_data_source & input,
                                               const IW_STL_Hash_Map_int & activity,
                                               const IW_STL_Hash_Set & training_set_identifiers,
                                               resizable_array_p<IW_General_Fingerprint> & test_set_fingerprints)
{
  IW_TDT tdt;
  while (tdt.next(input))
  {
    if (! determine_bit_counts (tdt, activity, training_set_identifiers, test_set_fingerprints))
    {
      cerr << "Cannot process '" << tdt << "'\n";
      return 0;
    }
  }

  return test_set_fingerprints.size();
}

int
Bayesian_Classification::determine_bit_counts (const char * fname,
                                               const IW_STL_Hash_Map_int & activity,
                                               const IW_STL_Hash_Set & training_set_identifiers,
                                               resizable_array_p<IW_General_Fingerprint> & test_set_fingerprints)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Bayesian_Classification::determine_bits_set:cannot open fingerprint file '" << fname << "'\n";
    return 0;
  }

  return determine_bit_counts (input, activity, training_set_identifiers, test_set_fingerprints);
}

int
Bayesian_Classification::determine_bit_counts (const IW_General_Fingerprint * fp,
                                               const int nfingerprints,
                                               const IW_STL_Hash_Map_int & activity,
                                               const int * training_set_identifiers)
{
  for (auto i = 0; i < nfingerprints; ++i)
  {
    if (0 == training_set_identifiers[i])
      continue;
  }

  return 1;
}

static int
need_to_examine_bit (const std::unordered_set<unsigned int> & bits_to_examine, 
                     const unsigned int b)
{
  const auto f = bits_to_examine.find(b);

  if (f == bits_to_examine.end())
    return 0;

  return 1;
}

/*
  After we have read all members of the training set, we can determine the
  weights associated with each bit.
  As we read the training set, we have put the counts into per-class specific
  containers.
*/

int
Bayesian_Classification::initialise_bit_weights ()
{
  if (0 == _n)
  {
    cerr << "Bayesian_Classification::initialise_bit_weights:no training data\n";
    return 0;
  }

  Bit_Frequency global_bit_frequency;

  if (verbose)
    cerr << "Training data has " << _n << " items\n";

  for (auto i = 0; i < _nclasses; ++i)
  {
    _baseline_class_probability[i] = static_cast<double>(_items_in_class[i]) / static_cast<double>(_n);

    if (verbose)
      cerr << "Class contains " << _items_in_class[i] << " items, p " << _baseline_class_probability[i] << endl;

    for (auto j = _bit_frequency[i].begin(); j != _bit_frequency[i].end(); ++j)
    {
      auto f = global_bit_frequency.find(j->first);

      int increment;
      if (denominator_is_molecules_with_feature)
        increment = 1;
      else
        increment = j->second;

      if (f == global_bit_frequency.end())
        global_bit_frequency[j->first] = increment;
      else
        f->second += increment;
    }
  }

  if (verbose > 1)
    cerr << "Global bit frequency contains " << global_bit_frequency.size() << " features\n";

  for (auto i = global_bit_frequency.begin(); i != global_bit_frequency.end(); ++i)
  {
    for (auto cls = 0; cls < _nclasses; ++cls)
    {
      int items_with_bit;    // in the class

      auto f = _bit_frequency[cls].find(i->first);

      if (f == _bit_frequency[cls].end())    // no instances of this bit in class 
        items_with_bit = 0;
      else
        items_with_bit = f->second;

      double w = (static_cast<double>(items_with_bit) + laplace_correction_factor) / (static_cast<double>(i->second) * _baseline_class_probability[cls] + laplace_correction_factor);
//    double w = (static_cast<double>(items_with_bit) + laplace_correction_factor) / ((static_cast<double>(i->second)+laplace_correction_factor) * _baseline_class_probability[cls]);
//    double w = (static_cast<double>(items_with_bit) + laplace_correction_factor) / ((static_cast<double>(i->second)) * _baseline_class_probability[cls]);

#ifdef BIT_EXAMINATION_ENABLED
      if (examine_all_bits || bits_to_examine.size())
      {
        if (examine_all_bits || need_to_examine_bit(bits_to_examine, i->first.second))
        {
          cerr << "fp " << i->first.first << " bit " << i->first.second << " in class " << cls << " items_with_bit " << items_with_bit << " total bit count " << i->second << " baseline " << _baseline_class_probability[cls] << ", w " << w << " logw " << log(w) << endl;
//        cerr << i->first.second << ' ' << items_with_bit << ' ' << i->second << ' ' << _baseline_class_probability[cls] << ' ' << w << "\n";
        }
      }
#endif

#ifdef DEBUG_IWNB
      cerr << "fp " << i->first.first << " bit " << i->first.second << " in class " << cls << " " << items_with_bit << " times, total " << i->second << ", w " << w << " logw " << log(w) << endl;
      cerr << i->first.second << ' ' << items_with_bit << ' ' << i->second << ' ' << _baseline_class_probability[cls] << ' ' << w << "\n";
#endif

      _bit_weight[cls][i->first] = log(w);
    }
  }

  return 1;
}

/*
  This is complicated by the fact that we must look across all classes to figure out the number
  of instances of each fingerprint
*/

/*int
Bayesian_Classification::remove_bits_below_support_level (int support)
{
  int nfp = number_fingerprints () + number_sparse_fingerprints();

  for (auto i = 0; i < nfp; ++i)
  {
    unordered_map<unsigned int, int> count;

    for (auto j = 0; j < _nclasses; ++j)
    {
      _bit_frequency[j].increment_count(i, count);
    }

    for (auto j = 0; j < _nclasses; ++j)
    {
      _bit_frequency[j].remove_bits_below_support_level(i, count, support);
    }
  }

  return 1;
}*/

int
Bayesian_Classification::score (const IW_General_Fingerprint & fp, double * rc) const
{
  if (sum_starts_with_log_class_probability)
  {
    for (auto i = 0; i < _nclasses; ++i)
    {
      rc[i] = log(_baseline_class_probability[i]);
    }
  }
  else
    std::fill_n(rc, _nclasses, 0.0);

  if (stream_for_bit_contributions.is_open())
    return _score_while_writing_bit_contributions (fp, rc, stream_for_bit_contributions);

  const int number_fixed = number_fingerprints();
  const int number_sparse = number_sparse_fingerprints();

  int fpi=0;
  
  while (fpi<number_fixed)
  {
    const IWDYFP & dyfp = fp[fpi];
    int b;
    int i = 0;

    while ((b = dyfp.next_on_bit(i)) >= 0)
    {
      FP_And_Bit fpb(fpi, b);

      for (auto cls = 0; cls < _nclasses; ++cls)    // look for this bit in every class
      {
        const auto f = _bit_weight[cls].find(fpb);

        if (f != _bit_weight[cls].end())
          rc[cls] += f->second;
      }
    }
    fpi++;
  }
  
  while (fpi<(number_fixed+number_sparse))
  {
    const Sparse_Fingerprint & sfp = fp.sparse_fingerprint(fpi-number_fixed);
          
    int i = 0;
    unsigned int b;
    int c;

    while (sfp.next_bit_set(i, b, c))
    {
//    const auto ob = b;
      if (distinguish_bits_by_count)
        b = bit_number_adjusted_for_count(b, c);

      FP_And_Bit fpb(fpi, b);

      for (auto cls = 0; cls < _nclasses; ++cls)
      {
        const auto f = _bit_weight[cls].find(fpb);

        if (f != _bit_weight[cls].end())
        {
#ifdef DEBUG_IWNB
          cerr << "Bit " << ob << " cls " << cls << " contribution " << f->second << endl;
#endif
          rc[cls] += f->second;
        }
      }
    }

    fpi++;
  }

  return 1;
}

int
Bayesian_Classification::_score_while_writing_bit_contributions (const IW_General_Fingerprint & fp, double * rc,
                                        IWString_and_File_Descriptor & output) const
{
  output << "ID: " << fp.id() << '\n';

  const int number_fixed = number_fingerprints();
  const int number_sparse = number_sparse_fingerprints();

  int fpi=0;
  
  while (fpi<number_fixed)
  {
    output << "FP" << fpi << '\n';

    const IWDYFP & dyfp = fp[fpi];
    int b;
    int i = 0;

    while ((b = dyfp.next_on_bit(i)) >= 0)
    {
      FP_And_Bit fpb(fpi, b);

      bool bit_written = false;

      for (auto cls = 0; cls < _nclasses; ++cls)    // look for this bit in every class
      {
        const auto f = _bit_weight[cls].find(fpb);

        if (f != _bit_weight[cls].end())
        {
          const auto c = f->second;
          rc[cls] += c;
          if (! bit_written)
          {
            output << b;
            bit_written = true;
          }
          output << ' ' << static_cast<float>(c);
        }
      }
      if (bit_written)
        output << '\n';
    }
    output << "|\n";
    output.write_if_buffer_holds_more_than(4096);
    fpi++;
  }
  
  while (fpi<(number_fixed+number_sparse))
  {
    output << "SFP" << (fpi - number_fixed) << '\n';

    const Sparse_Fingerprint & sfp = fp.sparse_fingerprint(fpi-number_fixed);
          
    int i = 0;
    unsigned int b;
    int c;

    while (sfp.next_bit_set(i, b, c))
    {
//    const auto ob = b;
      if (distinguish_bits_by_count)
        b = bit_number_adjusted_for_count(b, c);

      FP_And_Bit fpb(fpi, b);

      bool bit_written = false;

      for (auto cls = 0; cls < _nclasses; ++cls)
      {
        const auto f = _bit_weight[cls].find(fpb);

        if (f != _bit_weight[cls].end())
        {
#ifdef DEBUG_IWNB
          cerr << "Bit " << ob << " cls " << cls << " contribution " << f->second << endl;
#endif
          const auto c = f->second;
          rc[cls] += c;
          if (! bit_written)
          {
            output << b;
            bit_written = true;
          }

          output << ' ' << static_cast<float>(c);
        }
      }

      if (bit_written)
        output << '\n';
    }

    output << "|\n";
    output.write_if_buffer_holds_more_than(4096);
    fpi++;
  }

  output << "Score";
  for (auto i = 0; i < _nclasses; ++i)
  {
    output << ' ' << rc[i];
  }

  output << "\n";

  output.flush();

  return 1;
}


int
Bayesian_Classification::score (IW_TDT & tdt,
                                double * rc) const
{
  int fatal;

  IW_General_Fingerprint fp;

  if (! fp.construct_from_tdt (tdt, fatal))
  {
    cerr << "Invalid tdt '" << tdt << "', cannot build fingerprint\n";
    return 0;
  }

  return score (fp, rc);
}

int
Bayesian_Classification::score (iwstring_data_source & input,
                                double * rc) const
{
  IW_TDT tdt;

  while (tdt.next(input))
  {
    if (! score (tdt, rc))
      return 0;
  }

  return 1;
}

int
Bayesian_Classification::build (const char * fname,
                                IWString * & cname,
                                IWString & activity_name,
                                Fingerprint_Info_Read_From_Model & firfm)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Bayesian_Classification::build:cannot open input file '" << fname << "'\n";
    return 0;
  }

  return build (input, cname, activity_name, firfm);
}

int
Bayesian_Classification::_write_bit_weight(const int ndx,
                                           const IWString * class_label,
                                           IWString_and_File_Descriptor & output) const
{
  output << GFP_NB_CLASS << class_label[ndx] << '\n';

  for (auto i : _bit_weight[ndx])
  {
    output << output_separator << i.first.first << output_separator << i.first.second << output_separator << i.second << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }
  output << "|\n";

  return 1;
}

int
Bayesian_Classification::_read_bit_weight (Bit_Weight & bw,
                                           iwstring_data_source & input)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if ('|' == buffer)
      break;

    if (3 != buffer.nwords())
    {
      cerr << "Bayesian_Classification::_read_bit_weight:invalid record '" << buffer << "'\n";
      return 0;
    }

    int i = 0;
    const_IWSubstring token;
    buffer.nextword(token, i);

    int fp;
    if (! token.numeric_value(fp) || fp < 0)
    {
      cerr << "Bayesian_Classification::_read_bit_weight:invalid fp specification '" << buffer << "'\n";
      return 0;
    }

    buffer.nextword(token, i);
    unsigned int b;
    if (! token.numeric_value(b))
    {
      cerr << "Bayesian_Classification::_read_bit_weight:invalid bit '" << buffer << "'\n";
      return 0;
    }

    buffer.nextword(token, i);
    double w;
    if (! token.numeric_value(w))
    {
      cerr << "Bayesian_Classification::_read_bit_weight:invalid weight '" << buffer << "'\n";
      return 0;
    }

    FP_And_Bit fpb(fp, b);

    const auto f = bw.find(fpb);
    if (f != bw.end())
    {
      cerr << "Bayesian_Classification::_read_bit_weight:duplicate bit information '" << buffer << "'\n";
      return 0;
    }

    bw[fpb] = w;
  }

  return 1;
}

int
Bayesian_Classification::build (iwstring_data_source & input,
                                IWString * & class_label,
                                IWString & activity_name,
                                Fingerprint_Info_Read_From_Model & firfm)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Bayesian_Classification::build:empty file\n";
    return 0;
  }

  if (GFP_NB_HEADER != buffer)
  {
    cerr << "Bayesian_Classification::build:header record not correct, expected '" << GFP_NB_HEADER << "', got '" << buffer << "'\n";
    return 0;
  }

  if (! input.next_record(buffer))
  {
    cerr << "Bayesian_Classification::build:no data\n";
    return 0;
  }

  if (! buffer.starts_with(GFP_NB_ACT))
  {
    cerr << "Bayesian_Classification::build:expected '" << GFP_NB_ACT << "' got '" << buffer << "'\n";
    return 0;
  }

  activity_name = buffer;
  activity_name.remove_leading_words(1);

  while (input.next_record(buffer))
  {
    if (buffer.starts_with(GFP_NB_FX))
    {
      if (! firfm.add(buffer))
      {
        cerr << "Bayesian_Classification::build:invalid fixed size fingerprint specification '" << buffer << "'\n";
        return 0;
      }
    }
    else if (buffer.starts_with(GFP_NB_SP))
    {
      if (! firfm.add(buffer))
      {
        cerr << "Bayesian_Classification::build:invalid sparse fingerprint specification '" << buffer << "'\n";
        return 0;
      }
    }
    else
      break;
  }

  if (0 == firfm.number_fingerprints())
  {
    cerr << "Bayesian_Classification::build:no fingerprint tags specified\n";
    return 0;
  }

  if (! buffer.starts_with(GFP_NB_CLASSES))
  {
    cerr << "Bayesian_Classification::build:class labels, not '" << GFP_NB_CLASSES << "', got '" << buffer << "'\n"; 
    return 0;
  }

  _nclasses = buffer.nwords() - 1;           // must be class record indicator and at least two classes

  if (_nclasses < 2)         
  {
    cerr << "Bayesian_Classification::build:invalid class label record - must be at least 3 tokens '" << buffer << "'\n";
    return 0;
  }

  class_label = new IWString[_nclasses]; 

  IW_STL_Hash_Map_int label_to_ndx;

  int i = 0;
  const_IWSubstring token;
  buffer.nextword(token, i);

  IWString tmpc;

  for (int j = 0; buffer.nextword(tmpc, i); ++j)
  {
    class_label[j] << output_separator << tmpc;
    label_to_ndx[tmpc] = j;
  }

  if (! input.next_record(buffer))
  {
    cerr << "Bayesian_Classification::build:no baseline probability record\n";
    return 0;
  }

  if (! buffer.starts_with(GFP_NB_BP))
  {
    cerr << "Bayesian_Classification::build:expected '" << GFP_NB_BP << "' got '" << buffer << "', confused\n";
    return 0;
  }

  _baseline_class_probability = new double[_nclasses];

  IWString * tmps = new IWString[_nclasses]; std::unique_ptr<IWString[]> free_tmps(tmps);

  i = 0;
  buffer.nextword(token, i);

  for (int j = 0; buffer.nextword(token, i); ++j)
  {
    if (! token.numeric_value(_baseline_class_probability[j]) || _baseline_class_probability[j] < 0.0 || _baseline_class_probability[j] > 1.0)
    {
      cerr << "Bayesian_Classification::build:invalid baseline probability '" << buffer << "'\n";
      return 0;
    }
  }

  _bit_weight = new Bit_Weight[_nclasses];

  for (int i = 0; i < _nclasses; ++i)
  {
    if (! input.next_record(buffer))
      break;

    if (! buffer.starts_with(GFP_NB_CLASS))
    {
      cerr << "Bayesian_Classification::build:expecting class identifier '" << GFP_NB_CLASS << "', got '" << buffer << "'\n";
      return 0;
    }

    if (2 != buffer.nwords())
    {
      cerr << "Bayesian_Classification::build:class information must contain " << GFP_NB_CLASS << " and class label, '" << buffer << "' not recognised\n";
      return 0;
    }

    IWString token;

    int j = 0;
    buffer.nextword(token, j);
    buffer.nextword(token, j);

    const auto f = label_to_ndx.find(token);
    if (f == label_to_ndx.end())
    {
      cerr << "Bayesian_Classification::build:huh, where did class label '" << token << "' from '" << buffer << "' come from?\n";
      return 0;
    }

    if (! _read_bit_weight(_bit_weight[f->second], input))
    {
      cerr << "Bayesian_Classification::build:cannot read bit weight data, after line " << input.lines_read() << endl;
      return 0;
    }

    label_to_ndx.erase(token);
  }

  if (label_to_ndx.size() > 0)
  {
    cerr << "Bayesian_Classification::build:not all classes read, possibly missing";
    for (auto i : label_to_ndx)
    {
      cerr << ' ' << i.first;
    }
    cerr << endl;
  }

  return 1;
}

int
Bayesian_Classification::do_write (const IWString * class_label,
                                   const IWString & activity_name,
                                   const char * fname) const
{
  IWString_and_File_Descriptor output;

  if (! output.open(fname))
  {
    cerr << "Bayesian_Classification::do_write:cannot open '" << fname << "'\n";
    return 0;
  }

  return do_write(class_label, activity_name, output);
}

int
Bayesian_Classification::do_write (const IWString * class_label,
                                   const IWString & activity_name,
                                   IWString_and_File_Descriptor & output) const
{
  output << GFP_NB_HEADER << '\n';

  output << GFP_NB_ACT << output_separator << activity_name << '\n';

  for (int i = 0; i < number_fingerprints(); ++i)
  {
    output << GFP_NB_FX << output_separator << i << output_separator << fixed_fingerprint_tag(i) << '\n';
  }

  for (int i = 0; i < number_sparse_fingerprints(); ++i)
  {
    output << GFP_NB_SP << output_separator << i << output_separator << sparse_fingerprint_tag(i) << '\n';
  }

  output << GFP_NB_CLASSES;

  for (int i = 0; i < _nclasses; ++i)
  {
    output << class_label[i];
  }
  output << '\n';

  output << GFP_NB_BP;

  for (int i = 0; i < _nclasses; ++i)
  {
    output << output_separator << _baseline_class_probability[i];
  }
  output << '\n';

  for (int i = 0; i < _nclasses; ++i)
  {
    const auto & f = _bit_weight[i];
    _write_bit_weight(i, class_label, output);
  }

  return 1;
}


static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Naive Bayesian training and prediction\n";
  cerr << " Use one of : -A file.activity train.gfp test.gfp\n";
  cerr << "            : -A file.activity -S xxx -R test1.id -R test2.id ... all.gfp\n";
  cerr << "            : -A file.activity -S xxx all.gfp id1 id2 id3 ...\n";
  cerr << " -f <train,test> function of the extra files on the command line, training or test set\n";
  cerr << " -A <fname>     activity file for training\n";
  cerr << " -z             strip leading zero's from identifiers to aid matches\n";
  cerr << " -c <col>       activity in column <col> of activity file\n";
  cerr << " -t <float>     activity cutoff between active and inactive\n";
  cerr << " -H <hdr>       predicted activity output header (default from activity file)\n";
  cerr << " -F <tag>       tag to process - standard gfp_* syntax\n";
  cerr << " -l <float>     Laplace correction factor (default 1.0)\n";
//cerr << " -p <number>    min support level for output\n";
  cerr << " -S <stem>      if multiple -R options, file name stem for outputs\n";
  cerr << " -e <bucket>    adjust bit numbers so different counts become different bits\n";
  cerr << " -C N/P         run N cross validations using P percent as train\n";
  cerr << " -b             brief output, only output most likely class\n";
  cerr << " -B <thr>       only write predicted class likelihood when value above <thr>\n";
  cerr << " -O <fname>     write individual bit contributions to <fname>\n";
  cerr << " -Y  <fname>    write model information to <fname>\n";
  cerr << " -U  <fname>    read model information to <fname>\n";
  cerr << " -d             denominator in NB formula is number of molecules with feature\n";
  cerr << " -D <b1,b2...>  write out detailed information on bits b1, b2, ... Use '-D all' for all bits\n";
  cerr << " -v             verbose output\n";
  
  exit (rc);
}

class Output_Creator
{
  private:
    int _nclasses;
    const double ** _d;
    IWString * _class_names;

  public:
    Output_Creator ();
    ~Output_Creator ();

    int initialise (const IWString * cname, int nc);

    int do_output (const double *, IWString_and_File_Descriptor &);
    int do_output (const IW_TDT & tdt, const double * c, IWString_and_File_Descriptor & output);
};

Output_Creator::Output_Creator ()
{
  _nclasses = 0;
  _d = nullptr;
  _class_names = nullptr;

  return;
}

Output_Creator::~Output_Creator ()
{
  if (nullptr != _d)
    delete [] _d;

  if (nullptr != _class_names)
    delete [] _class_names;

  return;
}

int
Output_Creator::initialise (const IWString * cname,
                            int nc)
{
  _nclasses = nc;

  _class_names = new IWString[_nclasses];

  for (auto i = 0; i < nc; ++i)
  {
    _class_names[i] = cname[i];
  }

  _d = new const double *[nc];

  return 1;
}

int
Output_Creator::do_output (const double * c, IWString_and_File_Descriptor & output)
{
  for (auto i = 0; i < _nclasses; ++i)
  {
    _d[i] = c + i;
  }

  std::sort (_d, _d+_nclasses, [] (const double * d1, const double * d2) { return *d1 > *d2;});

  if (brief_output)
  {
    const auto cls = _d[0] - c;
    output << _class_names[cls];
  }
  else if (std::numeric_limits<double>::max() != all_scores_above_threshold)
  {
    int written = 0;

    for (auto i = 0; i < _nclasses; ++i)
    {
      if (*(_d[i]) < all_scores_above_threshold)
        break;
      const auto cls = _d[i] - c;
      output << _class_names[cls] << output_separator << static_cast<float>(*(_d[i]));
      written++;
    }

    if (0 == written)
    {
      const auto cls = _d[0] - c;
      output << _class_names[cls] << output_separator << static_cast<float>(*(_d[0]));
    }
  }
  else
  {
    for (auto i = 0; i < _nclasses; ++i)
    {
      const auto cls = _d[i] - c;
      output << _class_names[cls] << output_separator << static_cast<float>(*(_d[i]));
    }
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

int
Output_Creator::do_output (const IW_TDT & tdt,
                           const double * c,
                           IWString_and_File_Descriptor & output)
{
  IWString id;
  
  if (! tdt.dataitem_value(identifier_tag, id))
  {
    cerr << "Cannot extract identifier from tdt '" << tdt << "'\n";
    return 0;
  }

  if (id.contains(' '))
    id.truncate_at_first(' ');

  output << id;

  return do_output(c, output);
}

static int
gfp_naive_bayesian (iwstring_data_source & input,
                    Bayesian_Classification & bayesian_classification,
                    const IWString * class_names,
                    IWString_and_File_Descriptor & output)
{
  int nclasses = bayesian_classification.nclasses ();

  double * c = new double[nclasses]; std::unique_ptr<double[]> free_c(c);

  Output_Creator output_creator;

  output_creator.initialise(class_names, nclasses);

  IW_TDT tdt;
  while (tdt.next(input))
  {
    if (! bayesian_classification.score(tdt, c))
    {
      cerr << "Fatal error processing tdt '" << tdt << "'\n";
      return 0;
    }

    output_creator.do_output(tdt, c, output);
  }

  return 1;
}

static int
gfp_naive_bayesian (const char * fname,
                    Bayesian_Classification & bayesian_classification,
                    const IWString * class_names,
                    IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good ())
  {
    cerr << "Cannot open NB test set file '" << fname << "'\n";
    return 0;
  }

  return gfp_naive_bayesian (input, bayesian_classification, class_names, output);
}

static void
write_header (const IWString activity_name,
              const int nclasses,
              IWString_and_File_Descriptor & output)
{
  output << "ID" << output_separator << activity_name << " logPr";
  for (auto i = 1; i < nclasses; ++i)
  {
    output << output_separator << "C" << i << output_separator << "lPr";
  }
  output << '\n';

  return;
}

class NB_Cross_Validation
{
  private:
    IW_General_Fingerprint * _fp;
    int _nfingerprints;

    int * _activity;
    int * _training_set;

    double * _score;

    int _nclasses;

    int _predictions_made;
    int * _correct;

    IWString * _class_names;

    IW_STL_Hash_Map_int _id_to_position;

    IWString _activity_name;

    bool _extra_files_are_test_sets;

    Accumulator<float> _accuracy_accumulator;

//  when working out accuracy, we need a couple of temporary arrays

    int * _iic;

//  We can do cross validation

    int _nsplit;
    double _training_set_fraction;

    std::mt19937_64 _rng;

//  private functions

    int _assign_train_test (iwstring_data_source & input);
    int _assign_train_test (const char * fname);
    
    int _assign_training_set ();      // internal random division

    void _reset_performance_measures ();
    void _update_performance_measures (int ndx);

    void _train_classifier (Bayesian_Classification & bayesian_classification);

    int _doit (const char *);
    int _doit (IWString_and_File_Descriptor &);

  public:
    NB_Cross_Validation ();
    ~NB_Cross_Validation ();

    void set_extra_files_are_test_sets (bool s) { _extra_files_are_test_sets = s;}

    int initialise_cross_validation (const const_IWSubstring &);

    int read_all_fingerprints (const char * fname);
    int read_all_fingerprints (iwstring_data_source & input);

    void set_activity_name (const IWString & s) { _activity_name = s;}

    int assign_activity (const IW_STL_Hash_Map_int & a, const IWString *);

    int doit_cl (const Command_Line & cl, const IWString & output_file_name_stem);
    int doit_R  (const Command_Line & cl, const IWString & output_file_name_stem);

    int do_cross_validation (const IWString &);

    int report_performance_measures (std::ostream &);
    int report_overall_performance_measures (std::ostream & os);
};

static  std::random_device rd;

NB_Cross_Validation::NB_Cross_Validation () : _rng(rd())
{
  _fp = nullptr;
  _nfingerprints = 0;

  _nclasses = 0;

  _activity = nullptr;
  _training_set = nullptr;

  _class_names = nullptr;

  _score = nullptr;

  _correct = nullptr;

  _extra_files_are_test_sets = true;

  _iic = nullptr;

  _nsplit = 0;
  _training_set_fraction = 0.0;

  return;
}

NB_Cross_Validation::~NB_Cross_Validation ()
{
  if (nullptr != _fp)
    delete [] _fp;

  _nfingerprints = -17;

  if (nullptr != _activity)
    delete [] _activity;

  if (nullptr != _training_set)
    delete [] _training_set;

  if (nullptr != _score)
    delete [] _score;

  if (nullptr != _class_names)
    delete [] _class_names;

  _predictions_made = 0;

  if (nullptr != _correct)
    delete [] _correct;

  if (nullptr != _iic)
    delete [] _iic;

  return;
}

void
NB_Cross_Validation::_reset_performance_measures ()
{
  _predictions_made = 0;

  std::fill(_correct, _correct + _nclasses, 0);
}

int
NB_Cross_Validation::report_performance_measures (std::ostream & os)
{
  if (0 == _predictions_made)
  {
    os << "No predictions\n";
    return 0;
  }

  int ntrain = std::count(_training_set, _training_set + _nfingerprints, 1);

  std::fill(_iic, _iic + _nclasses, 0);     // items in class within the test set

  for (auto i = 0; i < _nfingerprints; ++i)
  {
    if (_training_set[i])
      continue;

    _iic[_activity[i]]++;
  }

  os << "Trained on " << ntrain << " of " << _nfingerprints << " fingerprints, test " << (_nfingerprints - ntrain) << " (";
  for (auto i = 0; i < _nclasses; ++i)
  {
    if (i > 0)
      os << ',';

    os << _iic[i];
  }
  os << ") :";

  float sum_accuracy = 0.0f;

  for (auto i = 0; i < _nclasses; ++i)
  {
    float f = static_cast<float>(_correct[i]) / static_cast<float>(_iic[i]);

    os << _class_names[i] << output_separator << _correct[i] << output_separator << f;

    sum_accuracy += f;
  }

  sum_accuracy /= float(_nclasses);

  os << output_separator << sum_accuracy << endl;

  _accuracy_accumulator.extra(sum_accuracy);

  return 1;
}

int
NB_Cross_Validation::report_overall_performance_measures (std::ostream & os)
{
  os << "Performed " << _accuracy_accumulator.n() << " cross validations. Accuracy between " << _accuracy_accumulator.minval() << " and " << _accuracy_accumulator.maxval() << " ave " << static_cast<float>(_accuracy_accumulator.average()) << endl;

  return 1;
}

int
NB_Cross_Validation::initialise_cross_validation (const const_IWSubstring & c)
{
  const_IWSubstring c1, c2;

  if (! c.split(c1, '/', c2) || 0 == c1.length() || 0 == c2.length())
  {
    cerr << "NB_Cross_Validation::initialise_cross_validation:cross validation must be of the form N/P, '" << c << "' invalid\n";
    return 0;
  }

  if (! c1.numeric_value(_nsplit) || _nsplit < 2)
  {
    cerr << "NB_Cross_Validation::initialise_cross_validation:invalid number of splits '" << c << "'\n";
    return 0;
  }

  int p;

  if (! c2.numeric_value(p) || p < 1 || p >= 100)
  {
    cerr << "NB_Cross_Validation::initialise_cross_validation:invalid training set percent '" << c << "'\n";
    return 0;
  }

  _training_set_fraction = static_cast<double>(p) / 100.0;

  return 1;
}

int
NB_Cross_Validation::read_all_fingerprints (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "NB_Cross_Validation::read_all_fingerprints:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_all_fingerprints (input);
}

int
NB_Cross_Validation::read_all_fingerprints (iwstring_data_source & input)
{
  _nfingerprints = input.count_records_starting_with(identifier_tag);

  if (0 == _nfingerprints)
  {
    cerr << "Cannot determine number of fingerprints in input\n";
    return 0;
  }

  _activity = new int[_nfingerprints];
  _training_set = new int[_nfingerprints];

  _fp = new IW_General_Fingerprint[_nfingerprints];

  for (auto i = 0; i < _nfingerprints; ++i)
  {
    IW_TDT tdt;

    if (! tdt.next(input))    // should not happen
    {
      cerr << "HUH, cannot read TDT " << i << endl;
      return 0;
    }

    int fatal;
    if (! _fp[i].construct_from_tdt(tdt, fatal))
    {
      cerr << "Cannot build fingerprint from tdt '" << tdt << "'\n";
      return 0;
    }

    IWString id = _fp[i].id();
    standardise_identifier(id);
    _fp[i].set_id(id);

    _id_to_position[id] = i;
  }

  return _nfingerprints;
}

int
NB_Cross_Validation::assign_activity (const IW_STL_Hash_Map_int & a,
                                      const IWString * cname)
{
  _nclasses = 0;

  for (auto i = 0; i < _nfingerprints; ++i)
  {
    IWString id = _fp[i].id();
    standardise_identifier(id);

    const auto f = a.find(id);

    if (f == a.end())
    {
      cerr << "NB_Cross_Validation::assign_activity:no activity data for '" << id << "'\n";
      return 0;
    }

    _activity[i] = f->second;

    if (_activity[i] > _nclasses)
      _nclasses = _activity[i];
  }

  _nclasses++;

  _score = new double[_nclasses];

  _correct = new int[_nclasses];

  _class_names = new IWString[_nclasses];

  for (auto i = 0; i < _nclasses; ++i)
  {
    _class_names[i] = cname[i];
  }

  _iic = new int[_nclasses];

  return 1;
}

int
NB_Cross_Validation::_assign_train_test (iwstring_data_source & input)
{
  IWString buffer;

  int s;    // set to zero or 1

  if (_extra_files_are_test_sets)
  {
    std::fill(_training_set, _training_set + _nfingerprints, 1);    // all are excluded from training set to start
    s = 0;
  }
  else
  {
    std::fill(_training_set, _training_set + _nfingerprints, 0);    // nothing in training set
    s = 1;
  }

  int rc = 0;

  while (input.next_record(buffer))
  {
    standardise_identifier(buffer);

    const auto f = _id_to_position.find(buffer);

    if (f != _id_to_position.end())
    {
      _training_set[f->second] = s;
      rc++;
    }
  }

  return rc;
}

int
NB_Cross_Validation::_assign_train_test (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "NB_Cross_Validation::_assign_train_test:cannot open '" << fname << "'\n";
    return 0;
  }

  return _assign_train_test(input);
}

int
NB_Cross_Validation::_assign_training_set ()
{
  std::fill(_training_set, _training_set + _nfingerprints, 0);

  int ntrain = static_cast<int>(_training_set_fraction * _nfingerprints + 0.4999);

  std::uniform_int_distribution<int> u(0, _nfingerprints - 1);

  int selected = 0;
  while (selected < ntrain)
  {
    int i = u(_rng);

    if (1 == _training_set[i])
      continue;

    _training_set[i] = 1;
    selected++;
  }

  return selected;
}

void
NB_Cross_Validation::_train_classifier (Bayesian_Classification & bayesian_classification)
{
  for (auto i = 0; i < _nfingerprints; ++i)
  {
    if (_training_set[i])
      bayesian_classification.determine_bit_counts (_fp[i], _activity[i]);
  }

  bayesian_classification.initialise_bit_weights ();

  return;
}

int
NB_Cross_Validation::_doit (IWString_and_File_Descriptor & output)
{
  Bayesian_Classification bayesian_classification;

  bayesian_classification.initialise(_nclasses);

  _train_classifier (bayesian_classification);

  _reset_performance_measures();

  Output_Creator output_creator;

  output_creator.initialise(_class_names, _nclasses);

  write_header(_activity_name, _nclasses, output);

  for (auto i = 0; i < _nfingerprints; ++i)
  {
    if (_training_set[i])
      continue;

    bayesian_classification.score(_fp[i], _score);

    output << _fp[i].id();

    output_creator.do_output(_score, output);

    _update_performance_measures(i);
  }

  output.flush();

  report_performance_measures(std::cout);

  return 1;
}

int
NB_Cross_Validation::_doit (const char * ofile)
{
  IWString_and_File_Descriptor output;
  
  if (! output.open(ofile))
  {
    cerr << "NB_Cross_Validation::_doit:cannot open output file '" << ofile << "'\n";
    return 0;
  }
  
  return _doit (output);
}

int
NB_Cross_Validation::doit_cl (const Command_Line & cl,
                              const IWString & output_file_name_stem)
{

  for (int i = 1; i < cl.number_elements(); ++i)
  {
    if (! _assign_train_test(cl[i]))
    {
      cerr << "Cannot determine train/test identifiers from '" << cl[i] << "'\n";
      return 0;
    }

    IWString ofname;
    ofname << output_file_name_stem << (i-1) << ".dat";

    if (! _doit (ofname.null_terminated_chars()))
    {
      cerr << "Cross validation based on '" << cl[i] << "' failed\n";
      return 0;
    }
  }

  report_overall_performance_measures(std::cout);

  return 1;
}

int
NB_Cross_Validation::doit_R (const Command_Line & cl,
                             const IWString & output_file_name_stem)
{
  IWString fname;
  for (int i = 0; cl.value('R', fname, i); ++i)
  {
    if (! _assign_train_test(fname.null_terminated_chars()))
    {
      cerr << "Cannot determine train/test identifiers from '" << fname << "'\n";
      return 0;
    }

    IWString ofname;
    ofname << output_file_name_stem << i << ".dat";

    if (! _doit (ofname.null_terminated_chars()))
    {
      cerr << "Cross validation based on '" << fname << "' failed\n";
      return 0;
    }
  }

  report_overall_performance_measures(std::cout);

  return 1;
}

void
NB_Cross_Validation::_update_performance_measures (int ndx)
{
  auto h = std::max_element(_score, _score + _nclasses) - _score;

//cerr << "From " << _score[0] << ' ' << _score[1] << " highest " << h << endl;

  _predictions_made++;

  if (h == _activity[ndx])
    _correct[h]++;

  return;
}

int
NB_Cross_Validation::do_cross_validation (const IWString & output_file_name_stem)
{
  for (auto i = 0; i < _nsplit; ++i)
  {
    _reset_performance_measures();

    _assign_training_set ();

    Bayesian_Classification bayesian_classification;

    bayesian_classification.initialise (_nclasses);

    _train_classifier (bayesian_classification);

    for (auto j = 0; j < _nfingerprints; ++j)
    {
      if (_training_set[j])
        continue;

      bayesian_classification.score(_fp[j], _score);

      _update_performance_measures(j);
    }

    report_performance_measures(std::cout);
  }

  report_overall_performance_measures(std::cout);

  return 1;
}

static int
gfp_naive_bayesian (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:c:t:F:l:p:zH:e:xbS:R:f:C:s:B:O:Y:U:dD:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! verbose)
    set_report_fingerprint_status(0);

  if (need_to_call_initialise_fingerprints(cl))
  {
    if (! initialise_fingerprints(cl, verbose))
    {
      cerr << "Cannot initialise fingerprint option(s)\n";
      usage(4);
    }
  }
  
  if (cl.option_present('U'))
  {
    if (cl.option_present('A'))
    {
      cerr << "When using a pre-built model (-U) cannot specify activity file (-A)\n";
      usage(1);
    }
  }
  else if ( ! cl.option_present('A'))
  {
    cerr << "Must specify activity file via the -A option\n";
    usage(4);
  }

  if (cl.option_present('d'))
  {
    denominator_is_molecules_with_feature = 1;

    if (verbose)
      cerr << "In the NB formula, denominator is molecules containing feature\n";
  }

#ifdef BIT_EXAMINATION_ENABLED
  if (cl.option_present('D'))
  {
    const_IWSubstring d;
    for (int i = 0; cl.value('D', d, i); ++i)
    {
      if ("all" == d)
      {
        examine_all_bits = 1;
        break;
      }

      int j = 0;
      const_IWSubstring token;
      while (d.nextword(token, j, ','))
      {
        unsigned int b;
        if (! d.numeric_value(b))
        {
          cerr << "invalid bit specification '" << d << "'\n";
          return 0;
        }

        bits_to_examine.insert(b);
      }
    }

    if (verbose && bits_to_examine.size() > 0)
      cerr << "Will write detailed information on " << bits_to_examine.size() << " bits\n";
  }
#endif

  IW_STL_Hash_Map_int activity;
  resizable_array_p<IWString> class_names;
  IWString activity_name;

  if (cl.option_present('A'))
  {
    if (cl.option_present('c'))
    {
      if (! cl.value('c', activity_column) || activity_column < 1)
      {
        cerr << "The activity column (-c) must be a whole +ve number\n";
        usage(4);
      }

      if (verbose)
        cerr << "Activity data in column " << activity_column << " of activity file\n";

      activity_column--;
    }

    if (cl.option_present('z'))
    {
      strip_leading_zeros = 1;

      if (verbose)
        cerr << "Will strip leading zero's from identifiers\n";
    }
  
    if (cl.option_present('t'))
    {
      if (! cl.value('t', activity_cutoff) )
      {
        cerr << "The activity cutoff (-t) must be a float number\n";
        usage(4);
      }
  
      if (verbose)
        cerr << "Activity cutoff set as " << activity_cutoff << " \n";
  
      activity_cutoff_set = 1;
    }

    const_IWSubstring a = cl.string_value('A');

    if (! read_activity_data(a, activity, class_names, activity_name))
    {
      cerr << "Cannot read activity data from '" << a << "'\n";
      return 5;
    }
                
    if (verbose)
      cerr << "Read " << activity.size() << " activity values from '" << a << "'\n";
  }

  IWString * cname = nullptr;

/*
  As we use the class names, they will always have a leading space.
  Also, we make them an array to avoid pointer dereferencing
*/

  if (! cl.option_present('U'))
  {
    if (0 == class_names.number_elements())
    {
      class_names.add(new IWString("-1"));
      class_names.add(new IWString("1"));
    }

    const auto nc = class_names.number_elements();

    if (verbose)
      cerr << "Found " << nc << " classes\n";

    if (nc < 2)
    {
      cerr << "All molecules same class, or no activity data, cannot continue\n";
      return 2;
    }

    cname = new IWString[nc];

    for (auto i = 0; i < nc; ++i)
    {
      cname[i] << output_separator << *(class_names[i]);
    }
  }

  if (cl.option_present('H'))
  {
    cl.value('H', activity_name);

    if (verbose)
      cerr << "Activity header '" << activity_name << "'\n";
  }

  if (cl.option_present('e'))
  {
    if (! cl.value('e', distinguish_bits_by_count) || distinguish_bits_by_count < 1)
    {
      cerr << "The distinguish bits by count option (-e) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will adjust bit numbers of count is > 1, dividing by " << distinguish_bits_by_count << endl;
  }

  if (cl.option_present('x'))
  {
    sum_starts_with_log_class_probability = 0;

    if (verbose)
      cerr << "Class sums omit log class probability\n";
  }
  
  if (cl.option_present('l'))
  {
    if (! cl.value('l', laplace_correction_factor) || laplace_correction_factor <=0)
    {
      cerr << "The Laplace correction factor (-l) must be a + float number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Laplace correction factor set as " << laplace_correction_factor << " \n";
  }

  if (cl.option_present('b') && cl.option_present('B'))
  {
    cerr << "The -b and -B options are mutually exclusive\n";
    usage(1);
  }

  if (cl.option_present('b'))
  {
    brief_output = 1;

    if (verbose)
      cerr << "Will only write id and prediction\n";
  }

  if (cl.option_present('B'))
  {
    if (! cl.value('B', all_scores_above_threshold))
    {
      cerr << "The output threshold (-B) option must be a valid floating point number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will only write class predictions if the value is >= " << all_scores_above_threshold << endl;
  }

  if (cl.option_present('p'))
  {
    if (! cl.value('p', min_support_level) || min_support_level < 0)
    {
      cerr << "The minimum support level (-p) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will drop bits present in " << min_support_level << " molecules or fewer\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

// Distinguish how we are to be invoked

  if (cl.option_present('U'))
    ;
  else if (cl.option_present('R'))
  {
    if (1 != cl.number_elements())
    {
      cerr << "When training set identifiers come via -R option, can be only one argument\n";
      usage(1);
    }
  
    if (! cl.option_present('S'))
    {
      cerr << "If doing multiple training set splits, must specify an output file name stem via the -S option\n";
      usage(1);
    }
  }
  else if (2 == cl.number_elements())
    ;
  else if (cl.option_present('C'))
  {
    if (1 != cl.number_elements())
    {
      cerr << "When doing cross validation (-C) only one file can be specified\n";
      return 2;
    }
  }
  else
  {
    if (! cl.option_present('S'))
    {
      cerr << "If doing multiple training set splits, must specify an output file name stem via the -S option\n";
      usage(1);
    }
  }

  if (cl.option_present('O'))
  {
    const char * o = cl.option_value('O');

    if (! stream_for_bit_contributions.open(o))
    {
      cerr << "Cannot open stream_for_bit_contributions '" << o << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Bit contributions written to '" << o << "'\n";

    stream_for_bit_contributions << "# Classes";
    for (unsigned int i = 0; i < class_names.size(); ++i)
    {
      stream_for_bit_contributions << ' ' << *(class_names[i]);

    }
    stream_for_bit_contributions << '\n';
  }

  set_sparsefp_warn_empty_data(0);

  Fingerprint_Info_Read_From_Model firfm;

  if (cl.option_present('U'))
  {
    const char * fname = cl.option_value('U');

    Bayesian_Classification bayesian_classification;

    if (! bayesian_classification.build(fname, cname, activity_name, firfm))
    {
      cerr << "Cannot read NB model from '" << fname << "'\n";
      return 1;
    }

    IWString_and_File_Descriptor output(1);

    write_header(activity_name, bayesian_classification.nclasses(), output);

    for (int i = 0; i < cl.number_elements(); ++i)
    {
      if (! gfp_naive_bayesian (cl[i], bayesian_classification, cname, output))
      {
        cerr << "NB model scoring failed '" << cl[i] << "'\n";
        return 2;
      }
    }

    output.flush();

    if (! firfm.check_fingerprint_validity())
    {
      cerr << "NB:YIPES:run just completed with invalid fingerprints!!!!!\n";
      return 1;
    }
  }
  else if (2 == cl.number_elements() && ! cl.option_present('f'))    // train and test gfp on command line
  {
    Bayesian_Classification bayesian_classification;

    bayesian_classification.initialise (activity);

    if (! bayesian_classification.determine_bit_counts (cl[0], activity))
    {
      cerr << "Cannot train NB model from '" << cl[0] << "'\n";
      return 2;
    }

    if (! bayesian_classification.initialise_bit_weights ())
      return 1;

    IWString_and_File_Descriptor output(1);

    write_header(activity_name, bayesian_classification.nclasses(), output);

    if (! gfp_naive_bayesian (cl[1], bayesian_classification, cname, output))
    {
      cerr << "NB model scoring failed '" << cl[1] << "'\n";
      return 2;
    }

    output.flush();

    if (cl.option_present('Y'))
    {
      const char * fname = cl.option_value('Y');

      if (! bayesian_classification.do_write(cname, activity_name, fname))
      {
        cerr << "Cannot write model to '" << fname << "'\n";
        return 1;
      }
    }
  }
  else if (cl.option_present('C'))
  {
    const_IWSubstring c = cl.string_value('C');

    NB_Cross_Validation nbcv;

    nbcv.set_activity_name(activity_name);

    if (! nbcv.initialise_cross_validation(c))
    {
      cerr << "Invalid cross validation specification '" << c << "'\n";
      return 2;
    }

    if (! nbcv.read_all_fingerprints (cl[0]))
    {
      cerr << "Cannot read fingerprints '" << cl[0] << "'\n";
      return 2;
    }

    if (! nbcv.assign_activity (activity, cname))
    {
      cerr << "Cannot associate activity data with fingerprints\n";
      return 2;
    }

    IWString output_file_name_stem;
    if (cl.option_present('S'))
      cl.value('S', output_file_name_stem);

    if (! nbcv.do_cross_validation (output_file_name_stem))
    {
      cerr << "Cross validation failed\n";
      return 2;
    }
  }
  else    // one or more cross validation runs within one fingerprint file
  {
    IWString output_file_name_stem;
    cl.value('S', output_file_name_stem);

    if (verbose)
      cerr << "Output files created with stem '" << output_file_name_stem << "'\n";

    NB_Cross_Validation nbcv;

    if (! cl.option_present('f'))
    {
      cerr << "Must specify whether the extra files are TRAIN or TEST set identifiers (-f)\n";
      usage(2);
    }
    else
    {
      IWString f = cl.string_value('f');
      f.to_lowercase();

      if (f.starts_with("tr"))
      {
        nbcv.set_extra_files_are_test_sets (false);
      }
      else if (f.starts_with("te"))
      {
        nbcv.set_extra_files_are_test_sets (true);
      }
      else
      {
        cerr << "Unrecognised file function directive '" << f << "'\n";
        usage(1);
      }
    }

    nbcv.set_activity_name(activity_name);

    if (! nbcv.read_all_fingerprints (cl[0]))
    {
      cerr << "Cannot read fingerprints '" << cl[0] << "'\n";
      return 2;
    }

    if (! nbcv.assign_activity (activity, cname))
    {
      cerr << "Cannot associate activity data with fingerprints\n";
      return 2;
    }

    int rc;
    if (cl.option_present('R'))
      rc = nbcv.doit_R (cl, output_file_name_stem);
    else
      rc = nbcv.doit_cl(cl, output_file_name_stem);

    if (0 == rc)
    {
      cerr << "NB cross validation failed\n";
      return 2;
    }
  }

  if (nullptr != cname)
    delete [] cname;

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_naive_bayesian (argc, argv);
  return rc;
}
