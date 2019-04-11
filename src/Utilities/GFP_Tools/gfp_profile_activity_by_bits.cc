/*
  Scans fingerprints and computes average activity associated with bits
*/

#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <memory>
#include <algorithm>
#include <random>
using namespace std;

#include "timsort.hpp"

#include "cmdline.h"
#include "iwstring_data_source.h"
#define IWQSORT_FO_IMPLEMENTATION
#include "iwqsort.h"
#include "iw_auto_array.h"
#include "iw_tdt.h"
#include "iw_stl_hash_map.h"
#include "iw_stl_hash_set.h"
#include "misc.h"

#include "gfp.h"

const char * prog_name = NULL;

static int verbose = 0;

static int show_constant_bits = 0;

static int show_absolute_differences = 1;

static int min_support_level = 0;

static int offset_bit_numbers = 0;

static float average_activity_of_whole_collection = 0.0;

static int skip_header_record_in_activity_file = 0;

static int write_most_informative_bits = 0;

static IWString_and_File_Descriptor stream_for_table_output;

static int output_precision = 0;

static int write_full_data = 1;

static double lower_pcorr_report_threshold = -3.0;
static double upper_pcorr_report_threshold =  3.0;

static IWString bit_prefix;
static int classification = 0;

typedef float activity_type_t;

static int records_written = 0;

static IW_STL_Hash_Map_String smiles;

static int strip_leading_zeros = 0;

/*
  Sept 2011
  The initial version of the programme did not report average
  activity for molecules without bits set in sparse fingerprints.
*/

static int show_set_and_nset_averages = 0;

static int write_data_for_gfp_adjust = 0;

static std::random_device rd;

static int nsubset = 100;

static int max_recursion_depth = 3;

static int min_bits_remaining_in_a_subset = max_recursion_depth + 1;

static int search_duration = 20;

/*
  I want to be able to handle arbitrary class labels
*/

//static IW_STL_Hash_Map_int class_to_number;
//static resizable_array_p<IWString> number_to_class;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Scans a fingerprint file and computes average activity associated with each bit\n";
  cerr << " -E fname       activity in file <fname>\n";
  cerr << " -E col=nn      activity is in column <nn> of the identifier\n";
  cerr << " -E skiphdr     skip header record in activity file\n";
  cerr << " -C             classification data\n";
  cerr << " -w             show constant bits\n";
  cerr << " -g             show signed rather than absolute differences\n";
  cerr << " -j             compute set and nset values for sparse fingerprints (slow)\n";
  cerr << " -p <num>       min number of times set or not-set in order to be output\n";
  cerr << " -T <fname>     file for statistical analysis\n";
  cerr << " -o             offset bit numbers by one - helps with tsubstructure files\n";
  cerr << " -n <num>       only write the <num> most significant bits (classification)\n";
  cerr << " -S <fname>     smiles file for bits - probably from dicer\n";
  cerr << " -r <rng>       report pcorr values outside <rng>\n";
  cerr << " -M ...         miscellaneous options, enter -M help for info\n";
  cerr << " -z             strip leading zeros when matching identifiers\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static IW_General_Fingerprint * pool = NULL;

static int pool_size = 0;

/*
  We need some function objects that convert the experimental data string
  representation into the form we need.
  For classification data, we maintain a cross reference between
  external class membership labels and the numbers we use
*/

class Classification_Data
{
  private:
    IW_STL_Hash_Map_int _class_label_to_class_number;
    resizable_array_p<IWString> _class_label_for_class_number;
    extending_resizable_array<int> _items_in_class;
  public:
    int operator()(const const_IWSubstring &, int &);

    const IWString & class_label(int c) const;

    int class_number(const const_IWSubstring &) const;

    int number_classes() const { return _class_label_for_class_number.number_elements();}

    const int * items_in_class() const { return _items_in_class.rawdata();}
};

int
Classification_Data::operator() (const const_IWSubstring & token,
                                 int & c)
{
  IW_STL_Hash_Map_int::const_iterator f = _class_label_to_class_number.find(token);

  if (f != _class_label_to_class_number.end())
  {
    c = (*f).second;
    _items_in_class[c]++;
    return 1;
  }

  c = _class_label_to_class_number.size();

  _class_label_to_class_number[token] = c;

  IWString * tmp = new IWString(token);

  _class_label_for_class_number.add(tmp);

  _items_in_class[c] = 1;

  return 1;
}

const IWString &
Classification_Data::class_label(int c) const
{
  assert (_class_label_for_class_number.ok_index(c));

  return *(_class_label_for_class_number[c]);
}

class Continuous_Response
{
  private:
  public:
    int operator()(const const_IWSubstring &, activity_type_t &) const;
};

int
Continuous_Response::operator() (const const_IWSubstring & token,
                                 activity_type_t & v) const
{
  int success;

  if (token.starts_with(' '))
  {
    const_IWSubstring tmp(token);
    tmp.strip_leading_blanks();

    success = tmp.numeric_value(v);
  }
  else
    success = token.numeric_value(v);

  if (! success)
  {
    cerr << "Invalid activity data '" << token << "'\n";
    return 0;
  }

  return 1;
}

static Classification_Data classification_data;
static Continuous_Response continuous_response;


/*
  For each bit, we need information about how many times it was assigned to the various classes
*/

class Class_Membership
{
  private:
    int _n;
    int * _in_class;
    int _number_classes;

  public:
    Class_Membership();
    ~Class_Membership();

    void extra(int c) { _in_class[c]++; _n++;}

    int n() const { return _n;}

    int report(int n, IWString_and_File_Descriptor &) const;

    void class_membership_ratios (const Classification_Data & classification_data, int n, float * r) const;
};

Class_Membership::Class_Membership()
{
  _n = 0;
  _number_classes = classification_data.number_classes();
  _in_class = new_int(_number_classes);

  return;
}

Class_Membership::~Class_Membership()
{
  delete [] _in_class;

  return;
}

void 
Class_Membership::class_membership_ratios (const Classification_Data & classification_data,
                                           int n,
                                           float * r) const
{
  const int * iic = classification_data.items_in_class();

  for (int i = 0; i < _number_classes; i++)
  {
    double probability_this_class = static_cast<double>(iic[i]) / static_cast<double>(pool_size);

    double pcorr = static_cast<double>(_in_class[i] + 1) / (probability_this_class * _n + 1.0);

    r[i] = log2(pcorr);
  }

  return;
}

int
Class_Membership::report(int n,
                         IWString_and_File_Descriptor & output) const
{
  if (_n == n)
  {
    output << " constant\n";
    return 1;
  }

  output << " N = " << _n;

  const int * iic = classification_data.items_in_class();

  for (int i = 0; i < _number_classes; i++)
  {
    const IWString & cl = classification_data.class_label(i);

    output << " in class " << cl << " " << _in_class[i];

    double probability_this_class = static_cast<double>(iic[i]) / static_cast<double>(pool_size);

//  pcorr = static_cast<double>(_in_class[i] + 1) / static_cast<double>(_n + k) / probability_active;
    double pcorr = static_cast<double>(_in_class[i] + 1) / (probability_this_class * _n + 1.0);

//  Aug 2009. Take log2 of this to keep symmetric around 1.0

    pcorr = log2(pcorr);
//  pcorr = static_cast<double>(_in_class[i] + probability_active * 10) / static_cast<double>(_n + 10) / probability_active;

    if (write_full_data)
      output << " NB " << static_cast<float>(pcorr);

    if (pcorr < lower_pcorr_report_threshold || pcorr > upper_pcorr_report_threshold)
      cerr << "High " << pcorr << " _in_class " << _in_class[i] << " _n " << _n << endl;
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
write_table_output (IWString_and_File_Descriptor & output,
                    unsigned int b,
                    int times_hit,
                    float ave)
{
  output << bit_prefix << b << " hit " << times_hit << " ave " << ave << " diff " << (ave - average_activity_of_whole_collection) << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static IWString_and_File_Descriptor &
write_pair_of_activities (float act1,
                          float act2,
                          IWString_and_File_Descriptor & output)
{
  output << " diff ";

  if (show_absolute_differences)
    output << static_cast<float>(fabs(act1 - act2));
  else 
    output << static_cast<float>(act2 - act1);

  return output;
}

static int
meets_support_requirement (int n, int pool_size)
{
  if (0 == min_support_level)
    return 1;

  if (n < min_support_level)
    return 0;

  if (pool_size - n < min_support_level)
    return 0;

  return 1;
}


static int
build_bit_cross_reference_meeting_support (int f,
                          IW_Hash_Map<unsigned int, unsigned int> & bit_to_sequence,
                          resizable_array<unsigned int> & bit_number)
{
  IW_Hash_Map<unsigned int, unsigned int> molecules_containing;

  for (int i = 0; i < pool_size; i++)
  {
    const IW_General_Fingerprint & fpi = pool[i];

    const Sparse_Fingerprint & sfp = fpi.sparse_fingerprint(f);

    int j = 0;
    unsigned int b;
    int notused;
    while (sfp.next_bit_set(j, b, notused))
    {
      molecules_containing[b]++;
    }
  }

  int suppressed_by_support_requirement = 0;

  for (IW_Hash_Map<unsigned int, unsigned int>::iterator i = molecules_containing.begin(); i != molecules_containing.end(); ++i)
  {
    if (! meets_support_requirement((*i).second, pool_size))
    {
      suppressed_by_support_requirement++;
      (*i).second = 0;
    }
  }
  if (verbose)
    cerr << "Suppressed " << suppressed_by_support_requirement << " of " << molecules_containing.size() << " bits due to support requirements\n";

  for (int i = 0; i < pool_size; i++)
  {
    const IW_General_Fingerprint & fpi = pool[i];

    const Sparse_Fingerprint & sfp = fpi.sparse_fingerprint(f);

    int j = 0;
    unsigned int b;
    int notused;
    while (sfp.next_bit_set(j, b, notused))
    {
      if (0 == molecules_containing[b])    // failed to meet support requirement
        continue;

      IW_Hash_Map<unsigned int, unsigned int>::const_iterator f = bit_to_sequence.find(b);
      if (f == bit_to_sequence.end())
      {
        unsigned int s = bit_to_sequence.size();
        bit_to_sequence[b] = s;
        bit_number.add(b);
      }
    }
  }

  int nbits = bit_to_sequence.size();

  assert (bit_number.number_elements() == nbits);

  if (verbose)
    cerr << "Sparse fingerprint " << f << " contains " << nbits << " bits\n";

  return nbits;
}

static int
build_bit_cross_reference(int f,
                          IW_Hash_Map<unsigned int, unsigned int> & bit_to_sequence,
                          resizable_array<unsigned int> & bit_number)
{
  for (int i = 0; i < pool_size; i++)
  {
    const IW_General_Fingerprint & fpi = pool[i];

    const Sparse_Fingerprint & sfp = fpi.sparse_fingerprint(f);

    int j = 0;
    unsigned int b;
    int notused;
    while (sfp.next_bit_set(j, b, notused))
    {
      IW_Hash_Map<unsigned int, unsigned int>::const_iterator f = bit_to_sequence.find(b);
      if (f == bit_to_sequence.end())
      {
        unsigned int s = bit_to_sequence.size();
        bit_to_sequence[b] = s;
        bit_number.add(b);
      }
    }
  }

  int nbits = bit_to_sequence.size();

  assert (bit_number.number_elements() == nbits);

  if (verbose)
    cerr << "Sparse fingerprint " << f << " contains " << nbits << " bits\n";

  return nbits;
}

static int
write_smiles_if_present_and_bit_number (unsigned int b,
                                        IWString_and_File_Descriptor & output)
{
  IWString bname;
  bname << b;

  IW_STL_Hash_Map_String::const_iterator f = smiles.find(bname);
  
  if (f != smiles.end())
    output << (*f).second << ' ';

  if (write_full_data)
    output << "bit " << bname;

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

/*
  In order to sort by the ratio, we need to keep track of the index of the
  bit and the associated ratio
*/

class Ndx_and_Ratio
{
  private:
    int _ndx;
    float _r;

  public:
    float ratio () const { return _r;}
    void set_ratio (float r) {_r = r;}

    int ndx () const { return _ndx;}
    void set_ndx(int n){ _ndx = n;}
};

class Ndx_and_Ratio_Comparator
{
  private:
  public:
    int operator () (const Ndx_and_Ratio &, const Ndx_and_Ratio &) const;
};

int
Ndx_and_Ratio_Comparator::operator () (const Ndx_and_Ratio & nr1, 
                                       const Ndx_and_Ratio & nr2) const
{
  float r1 = nr1.ratio();
  float r2 = nr2.ratio();

  if (r1 < r2)
    return -1;

  if (r1 > r2)
    return 1;

  return 0;
}

static void
sort_bits_by_activity (const activity_type_t * activity_with_bit,
                       const resizable_array<int> * items_with_bit,
                       const int nbits,
                       int * bits_ordered_by_activity)
{
  typedef std::pair<int, activity_type_t> IA;

  IA * p = new IA[nbits]; unique_ptr<IA[]> free_p(p);

  for (int i = 0; i < nbits; ++i)
  {
    p[i].first = i;
    p[i].second = activity_with_bit[i] / items_with_bit[i].number_elements();
  }

  std::sort(p, p + nbits, [](const IA & lhs, const IA & rhs) { return lhs.second > rhs.second;});

  for (int i = 0; i < nbits; ++i)
  {
    cerr << "Bit " << p[i].first << " ave activity " << p[i].second << endl;
    bits_ordered_by_activity[i] = p[i].first;
  }

  return;
}

static const activity_type_t needs_to_be_recomputed = -2774604.01;

class Subset_of_Bits
{
  private:
    resizable_array<int> _bit;

    float _average_activity;

    std::mt19937_64 _rng;

    int _age;

  public:
    Subset_of_Bits();

    int debug_print (std::ostream & output) const;

    int age () const { return _age;}

    int nbits () const { return _bit.number_elements();}

    int compare_ad_swap_if_better (const Subset_of_Bits & rhs);

    const resizable_array<int> & bits () const { return _bit;}

    void copy_from (const Subset_of_Bits & rhs);
    void mix_with  (Subset_of_Bits & rhs, const int nbits);

    float average_activity (const float *) const;

    activity_type_t average_activity () const { return _average_activity;}

    int initialise (const int pool_size, const int nsel);
    int initialise (const int * bits_ordered_by_activity, const int pool_size, const int nsel);

    activity_type_t compute_activity (const resizable_array<int> * items_with_bit, const activity_type_t * activity_with_bit);
    activity_type_t compute_activity_slow (const resizable_array<int> * items_with_bit, const activity_type_t * activity);

    activity_type_t recompute_average_activity_if_needed (const resizable_array<int> * items_with_bit, const activity_type_t * activity);

    void lose_bit ();
    void gain_bit (const int nbits);

    void add_bit (const int b) { _bit.add(b);}
    void pop_bit () { _bit.pop();}
};

Subset_of_Bits::Subset_of_Bits() : _rng(rd())
{
  _average_activity = needs_to_be_recomputed;

  _age = 0;

  return;
}

int
Subset_of_Bits::debug_print (std::ostream & output) const
{
  output << "Subset_of_Bits::debug_print: " << _bit.number_elements() << " bits";
  for (int i = 0; i < _bit.number_elements(); ++i)
  {
    output << ' ' << _bit[i];
  }

  output << " ave activity " << _average_activity << " age " << _age << '\n';

  return 1;
}

/*float
Subset_of_Bits::average_activity (const float * a) const
{
  const int n = _bit.number_elements();

  if (0 == n)
    return 0.0f;

  float rc = 0.0f;

  const int * b = _bit.rawdata();

  for (int i = 0; i < n; ++i)
  {
    rc += a[b[i]];
  }

  return rc / static_cast<float>(n);
}*/

static activity_type_t
sum_activities (const resizable_array<int> & items_with_bit,
                const activity_type_t * activity)
{
  const int n = items_with_bit.number_elements();

  if (0 == n)
    return 0.0f;

  activity_type_t rc = 0.0f;

  const int * iwb = items_with_bit.rawdata();

  for (int i = 0; i < n; ++i)
  {
    rc += activity[iwb[i]];
  }

  return rc;
}

activity_type_t
Subset_of_Bits::compute_activity_slow (const resizable_array<int> * items_with_bit,
                                       const activity_type_t * activity)
{
  activity_type_t tot = 0.0;
  int nb = 0;

  const int n = _bit.number_elements();

  const int * bb = _bit.rawdata();

  for (int i = 0; i < n; ++i)
  {
    const int b = bb[i];

    const resizable_array<int> & wb = items_with_bit[b];

    tot += sum_activities(wb, activity);
    nb += wb.number_elements();

//  cerr << " bit " << b << " had " << wb.number_elements() << " molecules hitting the bit, ave act " << (sum_activities(wb, activity)/nb) << endl;
  }

  _average_activity = tot / static_cast<activity_type_t>(nb);

  return _average_activity;
}

#ifdef ACCSAVE_STUFF
static activity_type_t * accsave = nullptr;
static int file_scope_nbits;

static int
check_accsave (const activity_type_t * activity_with_bit,
               const int nbits)
{
  int rc = 1;

  for (int i = 0; i < nbits; ++i)
  {
    if (accsave[i] == activity_with_bit[i])
      continue;

    cerr << "check_accsave:bit " << i << " changed, was " << accsave[i] << " now " << activity_with_bit[i] << endl;
    abort();
    rc = 0;
  }

  return rc;
}
#endif

activity_type_t
Subset_of_Bits::compute_activity (const resizable_array<int> * items_with_bit,
                                  const activity_type_t * activity_with_bit)
{
#ifdef DEBUG_COMPUTE_ACTIVITY
  cerr << "Computing activity with subset with " << _bit.number_elements() << " bits set\n";

  check_accsave(activity_with_bit, file_scope_nbits);
#endif

  _average_activity = 0.0f;

  const int n = _bit.number_elements();

  const int * bb = _bit.rawdata();

  activity_type_t tot = 0.0;
  int count = 0;

  for (int i = 0; i < n; ++i)
  {
    if (bb[i] < 0 || bb[i] >= 6271)
    {
      cerr << "Subset_of_Bits::compute_activity:invalid bit number " << bb[i] << endl;
      exit(1);
    }

    const auto xx = activity_with_bit[bb[i]];

    tot += xx;    // activity_with_bit[bb[i]];
    count += items_with_bit[bb[i]].number_elements();
//  cerr << "Fetched activity for items with bit " << bb[i] << " extra " << xx << " sum " << tot << " count " << count << endl;
  }

  if (0 == count)
  {
    cerr << "Subset_of_Bits::compute_activity:huh, no items with these bits\n";
    debug_print(cerr);
    exit(1);
  }

#ifdef DEBUG_COMPUTE_ACTIVITY
  cerr << "Sum " << tot << " count " << count << endl;
#endif

  _average_activity = tot / static_cast<activity_type_t>(count);

  return _average_activity;
}

int
Subset_of_Bits::initialise (const int pool_size,
                            const int nsel)
{
  int rc = 0;
  std::uniform_int_distribution<int> u(min_bits_remaining_in_a_subset, pool_size - 1);

  for (int i = 0; i < 100 && rc < nsel; ++i)
  {
    const int j = u(_rng);
    if (_bit.contains(j))
      continue;

    _bit.insert_in_order(j);
    rc++;
  }

  _age = 0;

  return nsel;
}

int
Subset_of_Bits::initialise (const int * bits_ordered_by_activity, 
                            const int nbits,
                            const int nsel)
{
  int rc = 0;
  std::array<double, 2> intervals{0.0, nbits - 2.0};
  std::array<double, 2> weights {1.0, 0.0};

  std::piecewise_linear_distribution<double> u(intervals.begin(), intervals.end(), weights.begin());

  for (int i = 0; i < 100 && rc < nsel; ++i)
  {
    const int j = static_cast<int>(u(_rng));
    if (_bit.contains(j))
      continue;

    _bit.insert_in_order(bits_ordered_by_activity[j]);
    rc++;
  }

  _age = 0;

  return nsel;
}

void
Subset_of_Bits::copy_from (const Subset_of_Bits & rhs)
{
  _bit = rhs._bit;

  _average_activity = rhs._average_activity;

  _age = 0;

  return;
}

template <typename T, typename O>
void
write_vector (const resizable_array<T> & v, O & output)
{
  for (int i = 0; i < v.number_elements(); ++i)
  {
    output << ' ' << v[i];
  }

  return;
}

static void
remove_dups (resizable_array<int> & bits,
             std::mt19937_64 & rng,
             const int min_bits_remaining_in_a_subset,
             const int nbits)
{
  int * b = bits.rawdata();

  const int n = bits.number_elements();

  int ndx = 1;
  for (int i = 1; i < n; ++i)
  {
    if (b[i-1] < b[i])
    {
      b[ndx] = b[i];
      ndx++;
    }
  }

  bits.resize_keep_storage(ndx);

  if (bits.number_elements() >= min_bits_remaining_in_a_subset)
    return;

  std::uniform_int_distribution<int> u(0, nbits - 1);

  while (bits.number_elements() < min_bits_remaining_in_a_subset)
  {
    const int b = u(rng);

    bits.insert_in_order_if_not_already_present(b);
  }

  return;
}

void
Subset_of_Bits::mix_with (Subset_of_Bits & rhs, const int nbits)
{
  const int nl = _bit.number_elements();
  const int nr = rhs._bit.number_elements();

  if (1 == nl && 1 == nr)
    return;

  int * bl = _bit.rawdata();
  int * br = rhs._bit.rawdata();

  std::shuffle(bl, bl + nl, _rng);
  std::shuffle(br, br + nr, _rng);

  int nswap = nl / 3;
  if (0 == nswap)
    nswap = 1;
  else if (nswap > nr / 2)
    nswap = nr/2;

  for (int i = 0; i < nswap; ++i)
  {
    std::swap(bl[i], br[i]);
  }

  std::sort(bl, bl + nl);
  std::sort(br, br + nr);

  remove_dups(_bit, _rng, min_bits_remaining_in_a_subset, nbits);
  remove_dups(rhs._bit, _rng, min_bits_remaining_in_a_subset, nbits);

  _average_activity = needs_to_be_recomputed;
  rhs._average_activity = needs_to_be_recomputed;

  _age++;

  return;
}

void
Subset_of_Bits::lose_bit ()
{
  if (_bit.number_elements() <= min_bits_remaining_in_a_subset)
    return;

  std::uniform_int_distribution<int> u(0, _bit.number_elements() - 1);

  _bit.remove_item(u(_rng));

  _average_activity = needs_to_be_recomputed;

  _age++;

  return;
}

void 
Subset_of_Bits::gain_bit (const int nbits)
{
  std::uniform_int_distribution<int> u(0, nbits - 1);

  for (int i = 0; i < 100; ++i)    // really should be a while(1) loop, but just being careful
  {
    int b = u(_rng);
    if (_bit.contains(b))
      continue;

    _bit.insert_in_order(b);
    break;
  }

  _average_activity = needs_to_be_recomputed;

  _age++;

  return;
}

activity_type_t
Subset_of_Bits::recompute_average_activity_if_needed (const resizable_array<int> * items_with_bit,
                                                      const activity_type_t * activity_with_bit)
{
  if (needs_to_be_recomputed != _average_activity)
    return _average_activity;

  return compute_activity(items_with_bit, activity_with_bit);
}

int
Subset_of_Bits::compare_ad_swap_if_better (const Subset_of_Bits & rhs)
{
  if (rhs._average_activity < _average_activity)
    return 0;

  copy_from(rhs);

  return 1;
}

bool
operator == (const Subset_of_Bits & lhs, const Subset_of_Bits & rhs)
{
  const auto & l = lhs.bits();
  const auto & r = rhs.bits();

  return l == r;
}

static activity_type_t
sum_activity_with_bit (const resizable_array<int> & items_with_bit,
                        const activity_type_t * activity)
{
  activity_type_t rc = static_cast<activity_type_t>(0.0);

  const auto n = items_with_bit.number_elements();

  for (int i = 0; i < n; ++i)
  {
//  cerr << "sum_activity_with_bit i = " << i << " bit " << items_with_bit[i] << " act " << activity[items_with_bit[i]] << endl;
    rc += activity[items_with_bit[i]];
  }

  return rc;
}

template <typename T>
int
identify_molecules_containing_bits_fixed (const IW_General_Fingerprint * pool,
                                    const int pool_size,
                                    const int fpnum,
                                    const resizable_array<T> & bits,
                                    const unordered_map<int, int> & local_bit_to_global_bit,
                                    resizable_array<int> * items_with_bit)
{
  const int n = bits.number_elements();

  for (int i = 0; i < n; ++i)
  {
    int b = bits[i];

    const auto f = local_bit_to_global_bit.find(b);
    if (f == local_bit_to_global_bit.end())
    {
      cerr << "identify_molecules_containing_bits:yipes, did not find bit " << b << " in local_bit_to_global_bit hash\n";
      return 0;
    }

    resizable_array<int> & c = items_with_bit[f->second];
    c.resize(min_support_level);

    for (int j = 0; j < pool_size; ++j)
    {
      if (pool[j][fpnum].is_set(b))
        c.add(j);
    }
  }

  return 1;
}

static int
identify_molecules_containing_bits_sparse (const IW_General_Fingerprint * pool,
                                    const int pool_size,
                                    const int fpnum,
                                    const resizable_array<unsigned int> & bits,
                                    const unordered_map<unsigned int, int> & local_bit_to_global_bit,
                                    resizable_array<int> * items_with_bit)
{
  const int n = bits.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const auto b = bits[i];

    const auto f = local_bit_to_global_bit.find(b);
    if (f == local_bit_to_global_bit.end())
    {
      cerr << "identify_molecules_containing_bits:yipes, did not find bit " << b << " in local_bit_to_global_bit\n";
      return 0;
    }
//  cerr << "Bit " << b << " is bit number " << f->second << endl;

    resizable_array<int> & c = items_with_bit[f->second];
    if (c.number_elements() > 0)
    {
      cerr << "Yipes, processing sparse bit " << b << " in fingerprint " << fpnum << " alread have " << c.number_elements() << " items counted\n";
      return 0;
    }
    assert (0 == c.number_elements());
    c.resize(min_support_level);

    for (int j = 0; j < pool_size; ++j)
    {
      if (pool[j].sparse_fingerprint(fpnum).is_set(b))
        c.add(j);
    }

//  cerr << "In sparse fp " << fpnum << " found " << c.number_elements() << " molecules with bit " << b << endl;
  }

  return 1;
}

static int
sparse_bits_meeting_support_levels (const IW_General_Fingerprint * pool,
                                    const int pool_size,
                                    const int fpnum,
                                    resizable_array<unsigned int> & bits)
{
  unordered_map<unsigned int, int> bit_count;

  Accumulator<activity_type_t> acc_with;

  for (int i = 0; i < pool_size; ++i)
  {
    const auto & sfp = pool[i].sparse_fingerprint(fpnum);

    int j = 0;
    unsigned int b;
    int c;
    while (sfp.next_bit_set(j, b, c))
    {
      auto f = bit_count.find(b);
      if (f == bit_count.end())
        bit_count[b] = 1;
      else
        f->second++;
    }
  }

  for (auto f : bit_count)
  {
    if (meets_support_requirement(f.second, pool_size))
      bits.add(f.first);
  }

  return bits.number_elements();
}

static int
fixed_bits_meeting_support_levels (const IW_General_Fingerprint * pool,
                                   const int pool_size,
                                   const int fpnum,
                                   resizable_array<int> & bits)
{
  const int nbits = pool[0][fpnum].nbits();

  int * is_set = new_int(nbits); unique_ptr<int[]> free_is_set(is_set);

  for (int i = 0; i < pool_size; ++i)
  {
    const auto & fp = pool[i][fpnum];

    for (int j = 0; j < nbits; ++j)
    {
      if (fp.is_set(j))
        is_set[j]++;
    }
  }

  for (int i = 0; i < nbits; ++i)
  {
    if (meets_support_requirement(is_set[i], pool_size))
      bits.add(i);
  }

  return bits.number_elements();
}


template <typename T>
int
setup_index_arrays (const int fpnum,
                    const resizable_array<T> & bits,     // the bits already found to meet the support requirements
                    unsigned int * global_bit_to_local_bit,
                    int * global_bit_to_fpnum,
                    unordered_map<T, int> & local_bit_to_global_bit,
                    int & ndx)
{
  const int n = bits.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const T b = bits[i];

    global_bit_to_local_bit[ndx] = b;
    global_bit_to_fpnum[ndx] = fpnum;
    local_bit_to_global_bit[b] = ndx;
    ndx++;
  }

  return 1;
}

static int
bits_meeting_support_level (const IW_General_Fingerprint * pool,
                            const int pool_size,
                            resizable_array<int> * fixed,
                            resizable_array<unsigned int> * sparse)
{
  for (int i = 0; i < number_fingerprints(); ++i)
  {
    fixed_bits_meeting_support_levels(pool, pool_size, i, fixed[i]);
  }

  for (int i = 0; i < number_sparse_fingerprints(); ++i)
  {
    sparse_bits_meeting_support_levels(pool, pool_size, i, sparse[i]);
  }

  return 1;
}

static void
recursive_enumeratin (const int b,
                      const int nbits,
                      Subset_of_Bits & sb,
                      const activity_type_t * activity_with_bit,
                      const resizable_array<int> * items_with_bit,
                      const int depth, 
                      const int max_depth, 
                      Accumulator<activity_type_t> * acc,
                      Subset_of_Bits * best,
                      IWString_and_File_Descriptor & output)
{
  for (int i = b; i < nbits; ++i)
  {
    sb.add_bit(i);
    sb.compute_activity(items_with_bit, activity_with_bit);
    acc[depth].extra(sb.average_activity());
    best[depth].compare_ad_swap_if_better(sb);
    if (depth < max_depth)
      recursive_enumeratin(i+1, nbits, sb, activity_with_bit, items_with_bit, depth + 1, max_depth, acc, best, output);
    sb.pop_bit();
    if (0 == depth)
    {
      cerr << "Top level completed bit " << i << " of " << nbits << " bits\n";
      for (int i = 0; i <= max_depth; ++i)
      {
        cerr << " depth " << i << " max " << acc[i].maxval() << " N = " << acc[i].n() << endl;
      }

//    if (i > 20)
//      return;
    }
  }

  return;
}

static int
remove_some_dups (Subset_of_Bits ** sbp, 
                  const int nsubset,
                  const int nbits,
                  std::uniform_int_distribution<int> & initial_bits,
                  std::mt19937_64 & rng)
{
  int rc = 0;

  int nsame = 0;

  for (int i = 1; i < nsubset; ++i)
  {
    if (*sbp[i-1] == *sbp[i])
    {
      nsame++;
      if (nsame < 3)
        continue;

      sbp[i-1]->initialise(nbits, initial_bits(rng));
      rc++;
    }
    else
      nsame = 0;
  }

  return rc;
}

/*
  We want to look for combinations of bits that seem informative.
  To do this, we want to assign a bit number to each bit.
  Therefore we need to keep a cross reference back to the bit number
  and the fingerprint.
  The fixed size fingerprints are easy, we know how many bits are
  in each.
  For the sparse fingerprints, we need a separate index for each
  of the sparse fingerprints - because the same bit could be set
  in different fingerprints
*/

static int
multi_bit_combinations (const IW_General_Fingerprint * pool,
                        const int pool_size,
                        const activity_type_t * activity,
                        const activity_type_t normalised_activity_cutoff,
                        IWString_and_File_Descriptor & output)
{
  const int nfixed = number_fingerprints();
  const int nsparse = number_sparse_fingerprints();

  resizable_array<int> * fixed = new resizable_array<int>[nfixed]; unique_ptr<resizable_array<int>[]> free_fixed(fixed);
  resizable_array<unsigned int> * sparse = new resizable_array<unsigned int>[nsparse]; unique_ptr<resizable_array<unsigned int>[]> free_sparse(sparse);

  bits_meeting_support_level(pool, pool_size, fixed, sparse);

  int nbits = 0;

  for (int i = 0; i < nfixed; ++i)
  {
    nbits += fixed[i].number_elements();
  }

  for (int i = 0; i < nsparse; ++i)
  {
    nbits += sparse[i].number_elements();
  }

  if (verbose)
    cerr << "Across " << nfixed << " fixed and " << nsparse << " sparse fingerprints, find " << nbits << " bits meeting support requirements\n";

  if (0 == nbits)
  {
    cerr << "no bits meet suppore requirements\n";
    return 0;
  }

// We need data structures that can map a bit in the local fingerprint to where that bit is in the global array

  unordered_map<int, int> * local_bit_to_global_bit_fixed = new unordered_map<int, int>[nfixed]; unique_ptr<unordered_map<int, int>[]> free_local_bit_to_global_bit_fixed(local_bit_to_global_bit_fixed);
  unordered_map<unsigned int, int> * local_bit_to_global_bit_sparse = new unordered_map<unsigned int, int>[nsparse]; unique_ptr<unordered_map<unsigned int, int>[]> free_local_bit_to_global_bit_sparse(local_bit_to_global_bit_sparse);

// for each bit, which molecules have that bit

  resizable_array<int> * items_with_bit = new resizable_array<int>[nbits]; unique_ptr<resizable_array<int>[]> free_items_with_bit(items_with_bit);

// for each bit, what is the bit number in the initial fingerprint

  unsigned int * global_bit_to_local_bit = new unsigned int[nbits]; unique_ptr<unsigned int[]> free_global_bit_to_local_bit(global_bit_to_local_bit);

// for each bit, from which fingerprint did it come

  int * global_bit_to_fpnum = new int[nbits]; unique_ptr<int[]> free_global_bit_to_fpnum(global_bit_to_fpnum);

  int ndx = 0;

  for (int i = 0; i < nfixed; ++i)
  {
    setup_index_arrays(i, fixed[i], global_bit_to_local_bit, global_bit_to_fpnum, local_bit_to_global_bit_fixed[i], ndx);
  }

  for (int i = 0; i < nsparse; ++i)
  {
    setup_index_arrays(i, sparse[i], global_bit_to_local_bit, global_bit_to_fpnum, local_bit_to_global_bit_sparse[i], ndx);
  }

  for (int i = 0; i < nfixed; ++i)
  {
    identify_molecules_containing_bits_fixed(pool, pool_size, i, fixed[i], local_bit_to_global_bit_fixed[i], items_with_bit);
  }

  for (int i = 0; i < nsparse; ++i)
  {
    identify_molecules_containing_bits_sparse(pool, pool_size, i, sparse[i], local_bit_to_global_bit_sparse[i], items_with_bit);
  }

  for (int i = 0; i < nbits; ++i)
  {
    const auto & x = items_with_bit[i];
    for (unsigned int j = 0; j < x.size(); ++j)
    {
      if (x[j] < 0 || x[j] >= pool_size)
      {
        cerr << "Possibly invalid items_with_bit value, i = " << i << " j = " << j << " value " << x[j] << endl;
      }
    }
  }

  activity_type_t * activity_with_bit = new activity_type_t[nbits]; unique_ptr<activity_type_t[]> free_activity_with_bit(activity_with_bit);

#ifdef ACCSAVE_STUFF
  accsave = new activity_type_t[nbits];
  file_scope_nbits = nbits;
#endif

  int * include_bit = new_int(nbits, 1); unique_ptr<int[]> free_include_bit(include_bit);

  const auto y = std::minmax_element(activity, activity + pool_size);

  const activity_type_t min_activity = *(y.first);
  const activity_type_t max_activity = *(y.second);

  cerr << "Activity range " << min_activity << " to " << max_activity << endl;
  cerr << "Line " << __LINE__ << " nbits " << nbits <<endl;

  int bits_included = 0;

  for (int b = 0; b < nbits; ++b)
  {
    activity_with_bit[b] = sum_activity_with_bit(items_with_bit[b], activity);
//  accsave[b] = activity_with_bit[b];
//  cerr << "Bit " << b << " computed activity " << activity_with_bit[b] << endl;

    activity_type_t ave = activity_with_bit[b] / items_with_bit[b].number_elements();

    float normalised = (ave  - min_activity) / (max_activity - min_activity);

    if (fabs(normalised - 0.5f) < normalised_activity_cutoff)
      include_bit[b] = 0;
    else
    {
      bits_included++;
      cerr << "Bit " << b << " ave activity " << activity_with_bit[b] << " norm " << normalised << endl;
    }
  }

  int * bits_ordered_by_activity = new int[nbits]; unique_ptr<int[]> free_bits_ordered_by_activity(bits_ordered_by_activity);
  sort_bits_by_activity(activity_with_bit, items_with_bit, nbits, bits_ordered_by_activity);

#ifdef DEBUG_COMPUTE_ACTIVITY
  for (int b = 0; b < nbits; ++b)
  {
    if (include_bit[b])
      cerr << "q Bit " << b << " activity_with_bit " <<  activity_with_bit[b] << endl;
  }
#endif

  if (verbose)
    cerr << "Based on minimum normalised activity threshold " << normalised_activity_cutoff << " left with " << bits_included << " bits\n";

  if (0 == bits_included)
  {
    cerr << "NO bits to process\n";
    return 0;
  }

  for (int i = 0; i < nbits; ++i)
  {
    Subset_of_Bits sbq;
    sbq.add_bit(i);
    activity_type_t a1 = sbq.compute_activity_slow(items_with_bit, activity);
    activity_type_t a2 = sbq.compute_activity(items_with_bit, activity_with_bit);

    if (fabs(a1 - a2) > 1.0e-05)
      cerr << "Mismatch on activity computations " << a1 << " v2 " << a2 << endl;
  }

  cerr << "Completed different method check\n";

  Subset_of_Bits * sb = new Subset_of_Bits[nsubset]; unique_ptr<Subset_of_Bits[]> free_sb(sb);

  Subset_of_Bits ** sbp = new Subset_of_Bits *[nsubset]; unique_ptr<Subset_of_Bits *[]> free_sbp(sbp);

  for (int i = 0; i < nsubset; ++i)
  {
    sb[i].initialise(nbits, 10);
    sb[i].compute_activity(items_with_bit, activity_with_bit);
    sbp[i] = sb + i;
  }

  std::sort(sbp, sbp + nsubset, [] (const Subset_of_Bits * sb1, const Subset_of_Bits * sb2) { return sb1->average_activity() > sb2->average_activity();});

  if (verbose)
  {
    Accumulator<activity_type_t> acc;

    for (int i = 0; i < nsubset; ++i)
    {
      sb[i].debug_print(cerr);
      acc.extra(sb[i].average_activity());
    }
    cerr << "Across " << nsubset << " random subsets, activities between " << acc.minval() << " and " << acc.maxval() << " ave " << static_cast<float>(acc.average()) << endl;
  }

//max_recursion_depth = -1;
  if (max_recursion_depth > 0)
  {
    Accumulator<activity_type_t> * acc = new Accumulator<activity_type_t>[max_recursion_depth + 1]; unique_ptr<Accumulator<activity_type_t>[]> free_acc(acc);
    Subset_of_Bits * best = new Subset_of_Bits[max_recursion_depth + 1]; unique_ptr<Subset_of_Bits[]> free_best(best);

    Subset_of_Bits s;
    recursive_enumeratin(0, nbits, s, activity_with_bit, items_with_bit, 0, max_recursion_depth, acc, best, output);
    cerr << "Following recursive exhaustive enumeration\n";
    for (int i = 0; i <= max_recursion_depth; ++i)
    {
      cerr << (i+1) << " bits btw " << acc[i].minval() << " and " << acc[i].maxval() << " ave " << acc[i].average() << endl;
      best[i].debug_print(cerr);
    }
  }

  std::mt19937_64 rng(rd());
  std::uniform_int_distribution<int> initial_bits(min_bits_remaining_in_a_subset, 10);
  std::uniform_real_distribution<double> r01(0.0, 1.0);

  int begin_restart = nsubset - nsubset / 10;
  int end_restart = nsubset;

  int begin_copy = nsubset / 2;
  int end_copy = nsubset / 2 + nsubset / 5;

  int begin_lose_bit = nsubset / 10;
  int end_lose_bit = nsubset;

  int begin_gain_bit = nsubset / 10;
  int end_gain_bit = nsubset;

  int begin_mix_with = 0;
  int end_mix_with = nsubset / 2;

  std::uniform_int_distribution<int> another_subset(0, nsubset/2);

  Accumulator<activity_type_t> acc2;

  time_t tzero = static_cast<time_t>(0);

  unsigned int istop;
  if (search_duration > 0)
  {
    tzero = time(NULL);
    istop = std::numeric_limits<unsigned int>::max();
  }
  else
    istop = 20000;

  for (unsigned int i = 0; i < istop; ++i)
  {
    for (int j = begin_restart; j < end_restart; ++j)
    {
      const auto r = r01(rng);

      if (r > 0.6)
        continue;

      if (r < 0.3)
        sbp[j]->initialise(nbits, initial_bits(rng));
      else
        sbp[j]->initialise(bits_ordered_by_activity, nbits, initial_bits(rng));
    }

    for (int j = begin_copy, ndx = 0; j < end_copy; ++j, ++ndx)
    {
      if (r01(rng) > 0.1)
        continue;

      if (sbp[j]->age() < 5)
        continue;

      sbp[j]->copy_from(*sbp[ndx]);
    }

    for (int j = begin_lose_bit; j < end_lose_bit; ++j)
    {
      if (r01(rng) < 0.3)
        sbp[j]->lose_bit();
    }

    for (int j = begin_gain_bit; j < end_gain_bit; ++j)
    {
      if (r01(rng) < 0.3)
        sbp[j]->gain_bit(nbits);
    }

    for (int j = begin_mix_with; j < end_mix_with; ++j)
    {
      if (r01(rng) < 0.6)
        continue;

      int k = another_subset(rng);

      if (sbp[k]->age() < 3)
        continue;

      if (k != j)
        sbp[j]->mix_with(*sbp[k], nbits);
    }

    for (int j = 0; j < nsubset; ++j)
    {
      acc2.extra(sb[j].recompute_average_activity_if_needed(items_with_bit, activity_with_bit));
    }

    std::sort(sbp, sbp + nsubset, [] (const Subset_of_Bits * sb1, const Subset_of_Bits * sb2) { return sb1->average_activity() > sb2->average_activity();});

    cerr << "Checking for time, i = " << i << endl;

    if (0 == i % 1000)
    {
      remove_some_dups(sbp, nsubset, nbits, initial_bits, rng);
      if (search_duration > 0 && i > 0)
      {
        const auto tnow = time(NULL);
        cerr << i << " iterations, now " << tnow << " t0 " << tzero << " diff " << (tnow - tzero) << endl;
        if (tnow - tzero > search_duration)
          break;
      }
    }

#ifdef ACCSAVE_STUFF
    for (int j = 0; j < nbits; ++j)
    {
      if (accsave[j] != activity_with_bit[j])
      {
        cerr << "Cycle " << i << " value " << j << " changed, old " << accsave[j] << " new " << activity_with_bit[j] << endl;
      }
    }
#endif
  }

  cerr << "Found " << acc2.n() << " values, btw " << acc2.minval() << " and " << acc2.maxval() << " ave " << static_cast<float>(acc2.average()) << endl;
  for (int i = 0; i < nsubset; ++i)
  {
    sbp[i]->debug_print(cerr);
  }

  return 1;
}

/*
  A Fixed width fingerprint. The index is the same as the bit number
*/

/*static int
do_write_most_informative_bits (const Classification_Data & classification_data,
                                const Class_Membership * cm,
                                int nbits,
                                IWString_and_File_Descriptor & output)
{
  int number_classes = classification_data.number_classes();

  float * rtmp = new float[number_classes]; iw_auto_array<float> free_rtmp(rtmp);

  Ndx_and_Ratio * nr = new Ndx_and_Ratio[nbits]; iw_auto_array<Ndx_and_Ratio> free_nr(nr);

  Ndx_and_Ratio_Comparator nrc;

  for (int i = 0; i < number_classes; i++)
  {
    int nsort = 0;
    for (int j = 0; j < nbits; j++)
    {
      if (! meets_support_requirement(cm[j].n(), pool_size))
        continue;

      cm[j].class_membership_ratios(classification_data, nbits, rtmp);

      nr[nsort].set_ndx(j);
      nr[nsort].set_ratio(rtmp[i]);
      nsort++;
    }

    if (0 == nsort)
      continue;

    if (nsort > 1)
      iwqsort (nr, nsort, nrc);

    int jstop = write_most_informative_bits;
    if (write_most_informative_bits > nsort)
      jstop = nsort;

    for (int j = 0; j < jstop; j++)
    {
      const Ndx_and_Ratio & nrj = nr[j];
      int ndx = nrj.ndx();
      const Class_Membership & cmi = cm[ndx];
      write_smiles_if_present_and_bit_number(ndx, output);
      cmi.report(nbits, output);
    }

    for (int j = 0; j < jstop; j++)
    {
      int k = nsort - j - 1;
      const Ndx_and_Ratio & nrk = nr[k];
      int ndx = nrk.ndx();
      const Class_Membership & cmi = cm[ndx];
      write_smiles_if_present_and_bit_number(ndx, output);
      cmi.report(nbits, output);
    }
  }

  return 1;
}*/

static int
do_write_data_for_gfp_adjust (const Classification_Data & classification_data,
                              const Class_Membership * cm,
                              const resizable_array<unsigned int> & bit_number,
                              IWString_and_File_Descriptor & output)
{
  int nbits = bit_number.number_elements();

  int number_classes = classification_data.number_classes();

  float * rtmp = new float[number_classes]; iw_auto_array<float> free_rtmp(rtmp);

  Ndx_and_Ratio * nr = new Ndx_and_Ratio[nbits]; iw_auto_array<Ndx_and_Ratio> free_nr(nr);

  Ndx_and_Ratio_Comparator nrc;

  for (int i = 0; i < number_classes; i++)
  {
    int nsort = 0;
    for (int j = 0; j < nbits; j++)
    {
      if (! meets_support_requirement(cm[j].n(), pool_size))
        continue;

      cm[j].class_membership_ratios(classification_data, nbits, rtmp);

      nr[nsort].set_ndx(j);
      nr[nsort].set_ratio(rtmp[i]);
      nsort++;
    }

    iwqsort(nr, nsort, nrc);

    for (int j = 0; j < nsort; j++)
    {
      const Ndx_and_Ratio & nrj = nr[j];
      int ndx = nrj.ndx();
      const Class_Membership & cmi = cm[ndx];
      write_smiles_if_present_and_bit_number(bit_number[ndx], output);
      cmi.report(nbits, output);
    }
  }

  output << "|\n";

  return 1;
}


/*
  When we have a sparse fingerprint and there is a cross reference back
  to bit numbers
*/

static int
do_write_most_informative_bits (const Classification_Data & classification_data,
                                const Class_Membership * cm,
                                const resizable_array<unsigned int> & bit_number,
                                IWString_and_File_Descriptor & output)
{
  int nbits = bit_number.number_elements();

  int number_classes = classification_data.number_classes();

  float * rtmp = new float[number_classes]; iw_auto_array<float> free_rtmp(rtmp);

  Ndx_and_Ratio * nr = new Ndx_and_Ratio[nbits]; iw_auto_array<Ndx_and_Ratio> free_nr(nr);

  Ndx_and_Ratio_Comparator nrc;

  for (int i = 0; i < number_classes; i++)
  {
    int nsort = 0;
    for (int j = 0; j < nbits; j++)
    {
      if (! meets_support_requirement(cm[j].n(), pool_size))
        continue;

      cm[j].class_membership_ratios(classification_data, nbits, rtmp);

      nr[nsort].set_ndx(j);
      nr[nsort].set_ratio(rtmp[i]);
      nsort++;
    }

    iwqsort(nr, nsort, nrc);

    int jstop = std::min(write_most_informative_bits, nsort);

    for (int j = 0; j < jstop; j++)
    {
      const Ndx_and_Ratio & nrj = nr[j];
      int ndx = nrj.ndx();
      const Class_Membership & cmi = cm[ndx];
      write_smiles_if_present_and_bit_number(bit_number[ndx], output);
      cmi.report(nbits, output);
    }

    for (int j = 0; j < jstop; j++)
    {
      int k = nsort - j - 1;
      const Ndx_and_Ratio & nrk = nr[k];
      int ndx = nrk.ndx();
      const Class_Membership & cmi = cm[ndx];
      write_smiles_if_present_and_bit_number(bit_number[ndx], output);
      cmi.report(nbits, output);
    }
  }

  return 1;
}

/*
  integer activities correspond to classification problems
*/

static int
profile_sparse_fingerprint(int f,                   // which sparse fingerprint
                           const int * activity,
                           IWString_and_File_Descriptor & output)
{
  IW_Hash_Map<unsigned int, unsigned int> bit_to_sequence;

  resizable_array<unsigned int> bit_number;

  if (! build_bit_cross_reference(f, bit_to_sequence, bit_number))
    return 0;

  int nbits = bit_number.number_elements();

  if (verbose > 1)
    cerr << "Sparse fingerprint contains " << nbits << " bits\n";

  Class_Membership * cm = new Class_Membership[nbits]; iw_auto_array<Class_Membership> free_cm(cm);

  for (int i = 0; i < pool_size; i++)
  {
    const IW_General_Fingerprint & fpi = pool[i];

    const Sparse_Fingerprint & sfp = fpi.sparse_fingerprint(f);

    int j = 0;
    unsigned int b;
    int notused;
    while (sfp.next_bit_set(j, b, notused))
    {
      IW_Hash_Map<unsigned int, unsigned int>::const_iterator f = bit_to_sequence.find(b);
      unsigned int s = (*f).second;

      cm[s].extra(activity[i]);
    }
  }

  if (write_most_informative_bits)
    return do_write_most_informative_bits(classification_data, cm, bit_number, output);

  for (int i = 0; i < nbits; i++)
  {
    const Class_Membership & cmi = cm[i];

    if (! meets_support_requirement(cmi.n(), pool_size))
      continue;

    if (smiles.size())
    {
      IWString bname;
      bname << (bit_number[i] + offset_bit_numbers);

      IW_STL_Hash_Map_String::const_iterator f = smiles.find(bname);
  
      if (f != smiles.end())
        output << (*f).second;
    }
    else
      output << "SFP " << f;

    output << " bit " << (bit_number[i] + offset_bit_numbers);

    records_written++;

    cm[i].report(nbits, output);
  }

  return 1;
}

static int
profile_sparse_fingerprint (int f,                     // which sparse fingerprint
                            const activity_type_t * activity,
                            IWString_and_File_Descriptor & output)
{
  IW_Hash_Map<unsigned int, unsigned int> bit_to_sequence;

  resizable_array<unsigned int> bit_number;

  if (min_support_level > 0)
  {
    if (! build_bit_cross_reference_meeting_support(f, bit_to_sequence, bit_number))
      return 0;
  }
  else if (! build_bit_cross_reference(f, bit_to_sequence, bit_number))
    return 0;

  int nbits = bit_number.number_elements();

  Accumulator<activity_type_t> * acc_set  = new Accumulator<activity_type_t>[nbits];
  Accumulator<activity_type_t> * acc_nset = new Accumulator<activity_type_t>[nbits];

  int * set_this_molecule = new int[nbits]; iw_auto_array<int> free_set_this_molecule(set_this_molecule);

  if (NULL == acc_set || NULL == acc_nset || NULL == set_this_molecule)
  {
    cerr << "Cannot allocate " << nbits << " accumulators\n";
    return 0;
  }

  for (int i = 0; i < pool_size; i++)
  {
    const IW_General_Fingerprint & fpi = pool[i];

    const Sparse_Fingerprint & sfp = fpi.sparse_fingerprint(f);

    if (show_set_and_nset_averages)
      set_vector(set_this_molecule, nbits, 0);

    unsigned int b;
    int j = 0;
    int notused;

    while (sfp.next_bit_set(j, b, notused))
    {
      IW_Hash_Map<unsigned int, unsigned int>::const_iterator f = bit_to_sequence.find(b);

      if (f == bit_to_sequence.end())   // failed to meet support requirements
        continue;

      unsigned int ndx = (*f).second;

      set_this_molecule[ndx] = 1;

      acc_set[ndx].extra(activity[i]);
    }

    if (show_set_and_nset_averages)
    {
      for (int j = 0; j < nbits; j++)
      {
        if (! set_this_molecule[j])
          acc_nset[j].extra(activity[i]);
      }
    }
  }

  Accumulator<activity_type_t> ave_act_acc;

  for (int i = 0; i < nbits; i++)
  {
    const Accumulator<activity_type_t> & acc = acc_set[i];

//  cerr << bit_number[i] << " n = " << acc.n() << " min " << acc.minval() << " max " << acc.maxval() << endl;

//  if (! meets_support_requirement(acc.n(), pool_size))
//    continue;

    write_smiles_if_present_and_bit_number(bit_number[i], output);
    output << " set " << acc.n() << " btw " << acc.minval() << " and " << acc.maxval();
    if (0 == acc.n())   // no range to show
      ;
    else if (show_set_and_nset_averages)
    {
      output << " ave " << static_cast<float>(acc_set[i].average());
      output << " Nset " << acc_nset[i].n() << " " << acc_nset[i].minval() << " to " << acc_nset[i].maxval() << " ave " << static_cast<float>(acc_nset[i].average());
      write_pair_of_activities(acc_set[i].average(), acc_nset[i].average(), output);
    }
    else
    {
      activity_type_t ave = static_cast<activity_type_t>(acc.average());
      output << " ave " << ave;
      ave_act_acc.extra(ave);
      write_pair_of_activities(average_activity_of_whole_collection, ave, output);
    }
    output << '\n';
    output.write_if_buffer_holds_more_than(32768);

    records_written++;

    if (stream_for_table_output.is_open())
      write_table_output(stream_for_table_output, bit_number[i], acc.n(), acc.average_if_available_minval_if_not());
  }

  if (verbose)
    cerr << "Average activities between " << ave_act_acc.minval() << " and " << ave_act_acc.maxval() << "\n";

  return output.good();
}

static int
profile_fingerprint(int f,
                    const int * c,
                    IWString_and_File_Descriptor & output)
{
  const IWDYFP & fp = pool[0].item(f);
  int nbits = fp.nbits();

  if (verbose)
    cerr << "Fingerprint " << f << " contains " << nbits << " bits\n";

  Class_Membership * cm = new Class_Membership[nbits]; iw_auto_array<Class_Membership> free_cm(cm);
  
  for (int i = 0; i < pool_size; i++)
  {
    const IWDYFP & fp = pool[i].item(f);
    for (int j = 0; j < nbits; j++)
    {
      if (fp.is_set(j))
        cm[j].extra(c[i]);
    }
  }

  if (write_most_informative_bits)
    return do_write_most_informative_bits(classification_data, cm, nbits, output);

  if (write_data_for_gfp_adjust)
  {
    output << fixed_fingerprint_tag(f) << '\n';
    return do_write_data_for_gfp_adjust(classification_data, cm, nbits, output);
  }

  for (int i = 0; i < nbits; i++)
  {
    Class_Membership & cmi = cm[i];

    int constant_bit;
    if (0 == cmi.n() || pool_size == cmi.n())
    {
      constant_bit = 1;
      if (! show_constant_bits)
        continue;
    }
    else
      constant_bit = 0;

    if (! meets_support_requirement(cmi.n(), pool_size))
      continue;

    output << "FP " << f << " bit " << i;

    cmi.report(pool_size, output);
  }

  return output.good();
}

static int
profile_fingerprint (int f,
                     const activity_type_t * activity,
                     IWString_and_File_Descriptor & output)
{
  const IWDYFP & fp = pool[0].item(f);
  int nbits = fp.nbits();

  if (verbose)
    cerr << "Fingerprint " << f << " contains " << nbits << " bits\n";

  Accumulator<activity_type_t> * acc_set = new Accumulator<activity_type_t>[nbits];
  Accumulator<activity_type_t> * acc_nset = new Accumulator<activity_type_t>[nbits];

  if (NULL == acc_set || NULL == acc_nset)
  {
    cerr << "Cannot allocate " << nbits << " accumulators\n";
    return 0;
  }
  
  for (int i = 0; i < pool_size; i++)
  {
    const IWDYFP & fp = pool[i].item(f);
    for (int j = 0; j < nbits; j++)
    {
      if (fp.is_set(j))
        acc_set[j].extra(activity[i]);
      else
        acc_nset[j].extra(activity[i]);
    }
  }

  for (int i = 0; i < nbits; i++)
  {
    Accumulator<activity_type_t> & s = acc_set[i];

    int constant_bit;
    if (0 == s.n() || static_cast<unsigned int>(pool_size) == s.n())
    {
      constant_bit = 1;
      if (! show_constant_bits)
        continue;
    }
    else
      constant_bit = 0;

    Accumulator<activity_type_t> & ns = acc_nset[i];

    if (! meets_support_requirement(ns.n(), pool_size))
      continue;

    double aset, anset;

    output << " bit " << i << " set " << s.n() << " ave ";
    if (constant_bit)
      output << '*';
    else
    {
      aset = s.average_if_available_minval_if_not();
      output << static_cast<float>(aset);
    }
      
    output << " not set " << ns.n() << " ave ";
    if (constant_bit)
      output << '*';
    else
    {
      anset = ns.average_if_available_minval_if_not();
      output << static_cast<float>(anset);
    }
      
    if (constant_bit)
      ;
    else
    {
      write_pair_of_activities(anset, aset, output);
      if (stream_for_table_output.is_open())
        write_table_output(stream_for_table_output, i, s.n(), s.average_if_available_minval_if_not());
    }
      
    output << '\n';

    output.write_if_buffer_holds_more_than(32768);
  }

  return output.good();
}

template <typename T>
int
get_experimental_data_from_identifier (int experimental_column,
                                       T * activity)
{
  for (int i = 0; i < pool_size; i++)
  {
    const IWString & id = pool[i].id();

    const_IWSubstring token;

    if (! id.word(experimental_column, token))
    {
      cerr << "Cannot extract column " << (experimental_column + 1) << " from '" << id << "'\n";
      return 0;
    }

    if (! token.numeric_value(activity[i]))
    {
      cerr << "Invalid numeric '" << id << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
add_to_hash (IW_STL_Hash_Set & identifiers_in_pool,
             const IWString & id)
{
  if (! strip_leading_zeros)
    identifiers_in_pool.insert(id);
  else if (! id.starts_with('0'))
    identifiers_in_pool.insert(id);
  else
  {
    IWString tmp(id);
    tmp.remove_leading_chars('0');
    identifiers_in_pool.insert(tmp);
  }

  return 1;
}

static int
in_hash (const IW_STL_Hash_Set & identifiers_in_pool,
         const IWString & query)
{
  if (identifiers_in_pool.contains(query))
    return 1;

  if (! strip_leading_zeros)
    return 0;

  if (! query.starts_with('0'))
    return 0;

  IWString tmp(query);

  tmp.remove_leading_chars('0');

  return identifiers_in_pool.contains(tmp);
}

template <typename T>
int
get_activity (const IW_STL_Hash_Map<IWString, T> & id_2_act,
              const IWString & id,
              T & act)
{
//typedef IW_STL_Hash_Map<IWString, T> Htype;
//typedef typename IW_STL_Hash_Map<IWString, T>::const_iterator Htype_iterator;

  auto f = id_2_act.find(id);

  if (f != id_2_act.end())
  {
    act = (*f).second;
    return 1;
  }

  if (! strip_leading_zeros)
    return 0;

  IWString tmp(id);
  tmp.remove_leading_chars('0');

  f = id_2_act.find(tmp);

  if (f == id_2_act.end())
    return 0;

  act = (*f).second;

  return 1;
}


template <typename T, typename C>
int
read_experimental_data_from_file (iwstring_data_source & input,
                                  int experimental_column,
                                  T * activity,
                                  C & conv)
{
  input.set_strip_trailing_blanks(1);

  IW_STL_Hash_Set identifiers_in_pool;

  for (int i = 0; i < pool_size; i++)
  {
    const IWString & id = pool[i].id();
    if (1 == id.nwords())
      add_to_hash(identifiers_in_pool, id);
    else
    {
      IWString tmp(id);
      tmp.truncate_at_first(' ');
      add_to_hash(identifiers_in_pool, tmp);
      pool[i].set_id(tmp);
    }
  }

  if (verbose)
    cerr << "Read " << identifiers_in_pool.size() << " identifiers from input fingerprint file\n";

  const_IWSubstring buffer;

  if (skip_header_record_in_activity_file)
    (void) input.next_record(buffer);

  typedef IW_STL_Hash_Map<IWString, T> Htype;

  Htype id_2_act;
//typedef typename IW_STL_Hash_Map<IWString, T>::const_iterator Htype_iterator;

  while (input.next_record(buffer))
  {
    IWString id, token;
    int i = 0;

    buffer.nextword(id, i);    // always in colum 1

    if (0 == id.length())   // huh!
      continue;

    if (experimental_column < 0)
      buffer.nextword(token, i);
    else
    {
      if (! buffer.word(experimental_column, token))
      {
        cerr << "Cannot extract column " << (experimental_column+1) << " from '" << buffer << "'\n";
        return 0;
      }
    }

    if (! in_hash(identifiers_in_pool, id))
      continue;

    T a;
    if (conv(token, a))
      id_2_act[id] = a;
    else if (1 == input.lines_read())
      cerr << "Ignoring non numeric activity on first record '" << buffer << "'\n";
    else
    {
      cerr << "Invalid activity '" << buffer << "'\n";
      return 0;
    }
  }

  if (verbose)
    cerr << "Read " << id_2_act.size() << " id-activity value pairs\n";

  if (id_2_act.size() < static_cast<unsigned int>(pool_size))
  {
    cerr << "Read " << id_2_act.size() << " activity values, but " << pool_size << " fingerprints, cannot continue\n";
    return 0;
  }

  for (int i = 0; i < pool_size; i++)
  {
    const IWString & id = pool[i].id();

    if (! get_activity(id_2_act, id, activity[i]))
    {
      cerr << "Yipes, no experimental data for '" << id << "'\n";
      return 0;
    }

//  cerr << " i = " << i << " activity " << activity[i] << endl;
  }

  return 1;
}

template <typename T, typename C>
int
read_experimental_data_from_file (const const_IWSubstring & fname,
                                  int experimental_column,
                                  T * activity,
                                  C & conv)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_translate_tabs(1);

  return read_experimental_data_from_file(input, experimental_column, activity, conv);
}

static int
report_class_memberships (const int * c,
                          int pool_size,
                          std::ostream & os)
{
  int max_class_number = c[0];
  int min_class_number = c[0];

  for (int i = 1; i < pool_size; i++)
  {
    if (c[i] > max_class_number)
      max_class_number = c[i];
    else if (c[i] < min_class_number)
      min_class_number = c[i];
  }

  os << "Classes " << min_class_number << " to " << max_class_number<< endl;

  int * items_in_class = new_int(max_class_number - min_class_number + 1); iw_auto_array<int> free_iic(items_in_class);

  for (int i = 0; i < pool_size; i++)
  {
    int c2 = c[i] - min_class_number;
    items_in_class[c2]++;
  }

  for (int i = min_class_number; i <= max_class_number; i++)
  {
    int ndx = i - min_class_number;
    float frac = static_cast<float>(items_in_class[i]) / static_cast<float>(pool_size);

    os << items_in_class[ndx] << " items in class " << i << " fraction " << frac << '\n';
  }

  return 1;
}

static int
do_classification (const IWString & experimental_file_name,
                   int experimental_column,
                   IWString_and_File_Descriptor & output)
{
  int * c = new int[pool_size]; std::unique_ptr<int[]> free_c(c);

  int rc = 0;
  if (experimental_file_name.length())
    rc = read_experimental_data_from_file(experimental_file_name, experimental_column, c, classification_data);
  else if (experimental_column > 0)
    rc = get_experimental_data_from_identifier(experimental_column, c);

  if (0 == rc)
  {
    cerr << "Cannot determine experimental data\n";
    return 7;
  }

  if (verbose)
    report_class_memberships(c, pool_size, cerr);

  int n = pool[0].nfingerprints();
  for (int i = 0; i < n; i++)
  {
    if (! profile_fingerprint(i, c, output))
    {
      cerr << "Cannot profile fingerprint " << i << endl;
      return i + 1;
    }
  }

  n = pool[0].number_sparse_fingerprints();
  for (int i = 0; i < n; i++)
  {
    if (! profile_sparse_fingerprint(i, c, output))
    {
      cerr << "Cannot profile sparse fingerprint " << i << endl;
      return i + 1;
    }
  }

  return 1;
}

static int 
build_pool (iwstring_data_source & input)
{
  int items_in_pool = 0;

  int tdts_read = 0;

  IW_TDT tdt;
  while (tdt.next(input))
  {
    tdts_read++;

    int fatal;
    if (! pool[items_in_pool].construct_from_tdt(tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      continue;
    }

    items_in_pool++;

    if (items_in_pool == pool_size)
    {
      if (verbose)
        cerr << "Pool is full, max " << pool_size << endl;
      break;
    }
  }

  pool_size = items_in_pool;

  if (verbose)
    cerr << "Read " << tdts_read << " TDT's, pool contains " << pool_size << " fingerprints\n";

  return 1;
}

static int 
build_pool (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (0 == pool_size)
  {
    pool_size = input.count_records_starting_with("|");
    if (0 == pool_size)
    {
      cerr << "Zero TDT items in input '" << fname << "'\n";
      return 0;
    }

    pool = new IW_General_Fingerprint[pool_size];
    if (NULL == pool)
    {
      cerr << "Bad news, cannot allocate pool " << pool_size << endl;
      return 0;
    }

    if (verbose)
      cerr << "File contains " << pool_size << " fingerprints\n";
  }

  return build_pool(input);
}

/*
  The smiles file will probably have come from dicer, so must have 3 tokens
*/

static int
read_smiles_record (const const_IWSubstring & buffer,
                    IW_STL_Hash_Map_String & smiles)
{
  if (buffer.nwords() < 2)
  {
    cerr << "Smiles data must have at least two records\n";
    return 0;
  }

  IWString smi, id;

  int i = 0;
  buffer.nextword(smi, i);
  buffer.nextword(id, i);

  if ("FRAGID" == id)
    buffer.nextword(id, i);

  if (0 == id.length())
  {
    cerr << "No identifier\n";
    return 0;
  }

  if (smiles.contains(id))
  {
    cerr << "Duplicate identifier in smiles file '" << id << "', ignored\n";
    return 1;
  }

  smiles[id] = smi;

  return 1;
}

static int
read_smiles (iwstring_data_source & input,
             IW_STL_Hash_Map_String & smiles)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! read_smiles_record(buffer, smiles))
    {
      cerr << "Cannot process smiles data '" << buffer << "'\n";
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

  return read_smiles(input, smiles);
}

static void
display_dash_m_options(std::ostream & os)
{
  os << " -M prec=<n>      set output precision\n";
  os << " -M brief         suppress output of bit numbers and other disagnostics\n";

  exit(1);
}

static int
gfp_profile_activity (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vs:E:F:P:wgu:p:T:p:Corn:S:M:zjY:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (need_to_call_initialise_fingerprints(cl))
  {
    if (! initialise_fingerprints(cl, verbose))
    {
      cerr << "Cannot initialise general fingerprint options\n";
      usage(17);
    }
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 2)
    {
      cerr << "Invalid pool size (-s option)\n";
      usage(4);
    }

    pool = new IW_General_Fingerprint[pool_size];
    if (NULL == pool)
    {
      cerr << "memory failure, cannot allocate " << pool_size << " items\n";
      return 4;
    }

    if (verbose)
      cerr << "Pool sized for " << pool_size << " items\n";
  }

  if (cl.option_present('o'))
  {
    offset_bit_numbers = 1;

    if (verbose)
      cerr << "Bit numbers offset by one\n";
  }

  if (cl.option_present('z'))
  {
    strip_leading_zeros = 1;

    if (verbose)
      cerr << "Will strip leading zeros\n";
  }

  if (cl.option_present('w'))
  {
    show_constant_bits = 1;

    if (verbose)
      cerr << "Constant bits will be shown\n";
  }

  if (cl.option_present('g'))
  {
    show_absolute_differences = 0;

    if (verbose)
      cerr << "Will show signed differences\n";
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', upper_pcorr_report_threshold))
    {
      cerr << "The pcorr report threshold (-r) must be a valid float\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will report pcorr values outside range " << upper_pcorr_report_threshold << endl;

    lower_pcorr_report_threshold = - upper_pcorr_report_threshold;
  }

// The -u option is deprecated, switch to -p for consistency with other progs

  if (cl.option_present('u'))
  {
    if (! cl.value('u', min_support_level) || min_support_level < 1)
    {
      cerr << "The minumum support level must be a positive whole number (-p)\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will only output bits with a minimum support level of " << min_support_level << " compounds hit\n";
  }

  if (cl.option_present('p'))
  {
    if (! cl.value('p', min_support_level) || min_support_level < 1)
    {
      cerr << "The minumum support level must be a positive whole number (-p)\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will only output bits with a minimum support level of " << min_support_level << " compounds hit\n";
  }

  if (cl.option_present('C'))
  {
    classification = 1;
    if (verbose)
      cerr << "Data is categorical data\n";
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', write_most_informative_bits) || write_most_informative_bits < 1)
    {
      cerr << "The write most informative bits option (-n) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will write the " << write_most_informative_bits << " most informative bits only\n";
  }

  if (cl.option_present('j'))
  {
    show_set_and_nset_averages = 1;

    if (verbose)
      cerr << "For sparse fingerprint output will show set and nset averages\n";
  }

  set_default_iwstring_float_concatenation_precision(4);

  if (cl.option_present('M'))
  {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++))
    {
      if (m.starts_with("prec="))
      {
        m.remove_leading_chars(5);
        if (! m.numeric_value(output_precision) || output_precision < 1)
        {
          cerr << "INvalid output precision '" << m << "'\n";
          return 5;
        }

        if (verbose)
          cerr << "Output precision set to " << output_precision << endl;

        set_default_iwstring_float_concatenation_precision(output_precision);
        set_default_iwstring_double_concatenation_precision(output_precision);
      }
      else if ("brief" == m)
      {
        write_full_data = 0;

        if (verbose)
          cerr << "Display of bit numbers suppressed\n";
      }
      else if ("help" == m)
      {
        display_dash_m_options(cerr);
      }
      else
      {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        display_dash_m_options(cerr);
      }
    }
  }

  if (! cl.option_present('E'))
  {
    cerr << "Must specify the source of experimental data via the -E option\n";
    usage(5);
  }

  int experimental_column = -1;
  IWString experimental_file_name;

  if (cl.option_present('E'))
  {
    int i = 0;
    const_IWSubstring e;
    while (cl.value('E', e, i++))
    {
      if (e.starts_with("col="))
      {
        e.remove_leading_chars(4);
        if (! e.numeric_value(experimental_column) || experimental_column < 1)
        {
          cerr << "Invalid experimental column '" << e << "'\n";
          return 6;
        }

        if (verbose)
          cerr << "Experimental result is column " << experimental_column << " of the identifier\n";

        experimental_column--;
      }
      else if ("skiphdr" == e)
      {
        skip_header_record_in_activity_file = 1;
        if (verbose)
          cerr << "Will skip the header record in the activity file\n";
      }
      else
      {
        experimental_file_name = e;
        if (verbose)
          cerr << "Experimental data in '" << experimental_file_name << "'\n";
      }
    }

    if (0 == experimental_file_name.length() && experimental_column < 0)
    {
      cerr << "No source of data, must specify '-E fname' or '-E col='\n";
      usage(3);
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('S'))
  {
    const char * s = cl.option_value('S');

    if (0 == read_smiles(s, smiles))
    {
      cerr << "Cannot read smiles from '" << s << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Read " << smiles.size() << " smiles from '" << s << "'\n";
  }

  if (cl.option_present('T'))
  {
    if(cl.option_present('C'))
    {
      cerr << "The table file output (-T) is not produced with classification data\n";
      return 3;
    }

    const char * t = cl.option_value('T');
    if (! stream_for_table_output.open(t))
    {
      cerr << "Cannot open file for tabular output '" << t << "'\n";
      return 3;
    }

    if (verbose)
     cerr << "Tabular output written to '" << t << "'\n";
  }

  if (! build_pool(cl[0]))
  {
    cerr << "Cannot read fingerprints from '" << cl[0] << "'\n";
    return 5;
  }

  IWString_and_File_Descriptor output(1);

  if (classification)
  {
    int rc = do_classification(experimental_file_name, experimental_column, output);
    if (NULL != pool)
      delete [] pool;

    return rc;
  }

  activity_type_t * activity = new activity_type_t[pool_size]; unique_ptr<activity_type_t[]> free_activity(activity);

  int rc;
  if (experimental_file_name.length())
    rc = read_experimental_data_from_file(experimental_file_name, experimental_column, activity, continuous_response);
  else if (experimental_column > 0)
    rc = get_experimental_data_from_identifier(experimental_column, activity);
  else
    return 5;

  if (0 == rc)
  {
    cerr << "Cannot determine experimental data\n";
    return 7;
  }

  Accumulator<activity_type_t> acc;
  acc.extra(activity, pool_size);

  if (verbose)
    cerr << "Activity values between " << acc.minval() << " and " << acc.maxval() << " ave " << acc.average() << '\n';

  average_activity_of_whole_collection = acc.average();

  if (cl.option_present('Y'))
  {
    activity_type_t a;
    if (! cl.value('Y', a) || a < 0.0 || a >= 1.0)
    {
      cerr << "The fractional activity cutoff value (-Y) must be a valid fraction\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will discard bits with average normalised activity below " << a << endl;

    multi_bit_combinations(pool, pool_size, activity, a, output);
  }
  else
  {
    int n = pool[0].nfingerprints();
    for (int i = 0; i < n; i++)
    {
      if (! profile_fingerprint(i, activity, output))
      {
        cerr << "Cannot profile fingerprint " << i << endl;
        return i + 1;
      }
    }

    n = pool[0].number_sparse_fingerprints();
    for (int i = 0; i < n; i++)
    {
      if (! profile_sparse_fingerprint(i, activity, output))
      {
        cerr << "Cannot profile sparse fingerprint " << i << endl;
        return i + 1;
      }
    }
  }

  output.flush();

  if (NULL != pool)
    delete [] pool;

  delete_gfp_file_scope_static_objects();

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_profile_activity(argc, argv);

  return rc;
}
