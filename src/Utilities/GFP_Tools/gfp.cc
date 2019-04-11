#include <stdlib.h>
#include <math.h>
#include <string.h>

//#include "tbb/scalable_allocator.h"

#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "iw_tdt.h"
#include "cmdline.h"
#include "iwcrex.h"

#include "gfp.h"
#include "tversky.h"
#include "various_distance_metrics.h"

#ifdef FB_ENTROPY_WEIGHTED_FPS
#include "fb_bits_and_weights.h"

static FB_Bits_and_Weights_Fixed_Width ffbwfw;

int
initialise_entropy_weighting_for_fingerprints (Command_Line & cl,
                                char flag,
                                int verbose)
{
  return ffbwfw.construct_from_command_line(cl, flag, verbose);
}

#endif

/*
  All fingerprints must have the same number of fingerprints.
  These variables could be made static class members

  When we get into convering sparse fingerprints to denser forms, we get
  great complexity around the number of fingerprints at various times.

  When the fingerprint is read from a file, there will be a given number
  of each kind of fingerprint, but during a computation, there may be
  a different number of each kind - specifically, the sparse fingerprints
  may have been converted to a denser form.
*/

static int bit_fingerprints_in_file = 0;

/*
  Because we may change sparse to dense fingerprints, we have separate variables for
  the number of fingerprints to allocate. Generally, fingerprints_to_allocate will be
  the same as number_fingerprints_to_use_in_computations
*/

static int fingerprints_to_allocate = 0;

static int number_fingerprints_to_use_in_computations = 0;

static int _number_sparse_fingerprints = 0;

static int _number_sparse_fingerprints_to_use_in_computations = 0;

static int number_fixed_size_counted_fingerprints_in_file = 0;

static int number_fixed_size_counted_fingerprints_to_allocate = 0;

static int number_fixed_size_counted_fingerprints_to_use_in_computations = 0;

static int number_multiconformer_01_fingerprints_in_file = 0;

static int number_multiconformer_fixed_size_counted_fingerprints_in_file = 0;

static int number_multiconformer_sparse_fingerprints_in_file = 0;

//  The weight given to each component will be the same for all fingerprints

static float * _fingerprint_weight = NULL;
static float * _sparse_fingerprint_weight = NULL;
static float * _fixed_size_counted_fingerprint_weight = NULL;

static int * _bits_in_fingerprint = NULL;

/*
  When we read fingerprints in Daylight ascii representation, the counts are truncated at 255.
  But if we read the ascii form 'b,c', then counts are not limited
*/

static int sparse_fingerprint_counts_limited = 1;   

static int report_fingerprint_status = 1;

void
set_report_fingerprint_status (int s)
{
  report_fingerprint_status = s;
}

void
set_fixed_width_fingerprint_weight (int ndx, double w)
{
  assert (NULL != _fingerprint_weight);
  assert (ndx >= 0 && ndx < bit_fingerprints_in_file);

  _fingerprint_weight[ndx] = w;

  return;
}

void
set_sparse_fingerprint_weight (int ndx, double w)
{
  assert (NULL != _sparse_fingerprint_weight);
  assert (ndx >= 0 && ndx < _number_sparse_fingerprints);

  _sparse_fingerprint_weight[ndx] = w;

  return;
}

/*
  Aug 2005. In order to keep sparse fingerprints to similar distances, we can scale their
  values. Note that this is different from _sparse_fingerprint_weight
*/

static float * _sparse_fingerprint_scaling_factor = NULL;

/*
  With multiconformer fingerprints, by default we find the highest
  similarity between the two sets of fingerprints. But for some
  applications we want the average similarity
*/

static int multiconformer_use_average_similarity = 0;

/*
  We can use Tanimoto, Dice or Cosine for the distances between sparse fingerprints
*/

static int * _sparse_fingerprint_distance_metric = NULL;

static int multiconformer_fingerprints_present = 0;

static float * _multiconformer_01_weight = NULL;
static int * _multiconformer_01_distance_metric = NULL;
static float * _multiconformer_fixed_size_counted_weight = NULL;
static int * _multiconformer_fixed_size_counted_distance_metric = NULL;
static float * _multiconformer_sparse_weight = NULL;
static int * _multiconformer_sparse_distance_metric = NULL;

/*
  Very important for both Tanimoto types to be 0
*/

#define SPARSE_DISTANCE_METRIC_TANIMOTO 0
#define DENSE_DISTANCE_METRIC_TANIMOTO 0

static int default_sparse_fingerprint_distance_metric = SPARSE_DISTANCE_METRIC_TANIMOTO;

static int default_dense_fingerprint_distance_metric = DENSE_DISTANCE_METRIC_TANIMOTO;

static int default_multiconformer_01_distance_metric = DENSE_DISTANCE_METRIC_TANIMOTO;
static int default_multiconformer_fixed_size_counted_distance_metric = DENSE_DISTANCE_METRIC_TANIMOTO;
static int default_multiconformer_sparse_distance_metric = SPARSE_DISTANCE_METRIC_TANIMOTO;

#define DICE_DISTANCE_METRIC 1
#define SPARSE_DISTANCE_METRIC_COSINE 4

/*
  The modified Tanimoto from Technometrics, May 2002 Vol 44 No. 2 Page 110
  Fligner, Verducci and Blower
*/

#define FVB_MODIFIED_TANIMOTO 7

#define RUSSEL_RAO 8

#define FORBES_SIMILARITY 9

#define SIMPLE_MATCHING 10

/*
  Sparse counted fingerprints can use the Manhattan distance metric
*/

#define MANHATTAN_DISTANCE_METRIC 11

#define SOERGEL_DISTANCE_METRIC 12
#define SOERGEL_VARIANT_DISTANCE_METRIC 13

#define CTAN_DISTANCE_METRIC 14

#define SPARSE_DISTANCE_METRIC_BINARY 15

#define DENSE_DISTNANCE_OVERLAP 16

static int * _dense_fingerprint_distance_metric = NULL;

/*
  With continuous properties, there are also a number of distance metrics we can use
*/

#define CONTINUOUS_PROPERTY_DISTANCE_METRIC_RATIO 0
#define CONTINUOUS_PROPERTY_DISTANCE_METRIC_EXP1  1
#define CONTINUOUS_PROPERTY_DISTANCE_METRIC_EXP2  2
#define CONTINUOUS_PROPERTY_DISTANCE_METRIC_CART  3
#define CONTINUOUS_PROPERTY_DISTANCE_METRIC_DICE  4

static similarity_type_t continuous_property_distance_metric_exp1_value = 0.0;
static similarity_type_t continuous_property_distance_metric_exp2_value = 0.0;

static int _continuous_property_distance_metric = CONTINUOUS_PROPERTY_DISTANCE_METRIC_RATIO;

#define INTEGER_PROPERTY_DISTANCE_METRIC_RATIO 0
#define INTEGER_PROPERTY_DISTANCE_METRIC_DICE 1

static int _integer_property_distance_metric = INTEGER_PROPERTY_DISTANCE_METRIC_RATIO;

/*
  The counted fingerprints can have their own distance metric
*/

#define FIXED_SIZE_COUNTED_FINGERPRINT_DISTANCE_METRIC_TANIMOTO 0

static int default_fixed_size_counted_fingerprint_distance_metric = FIXED_SIZE_COUNTED_FINGERPRINT_DISTANCE_METRIC_TANIMOTO;

static int * _fixed_size_counted_fingerprint_distance_metric = NULL;

//  The tag for each fingerprint
//  Note that these are never freed

static IWString * _fingerprint_tag = NULL;
static IWString * _sparse_fingerprint_tag = NULL;
static IWString * _fixed_size_counted_fingerprint_tag = NULL;
static IWString * _multiconformer_01_fingerprint_tag = NULL;
static IWString * _multiconformer_fixed_size_counted_fingerprint_tag = NULL;
static IWString * _multiconformer_sparse_fingerprint_tag = NULL;

/*
  When processing sparse fingerprints, we can sometimes speed things up
  by converting them to fixed width forms

  Note that we make no attempt to maintain consistency between the
  two possible conversions
*/

static int convert_sparse_fingerprint_to_bits = 0;

void
set_convert_sparse_fingerprint_to_bits (int s)
{
  convert_sparse_fingerprint_to_bits = s;
}

static int convert_sparse_fingerprint_to_fixed_width_counted = 0;

void
set_convert_sparse_fingerprint_to_fixed_width_counted (int s)
{
  convert_sparse_fingerprint_to_fixed_width_counted = s;
}

void
set_sparse_fingerprint_counts_limited(const int s)
{
  sparse_fingerprint_counts_limited = s;
}

/*
  We may have a guard fingerprint. When this fingerprint is evaluated, we
  check the distance at that stage and see whether or not it is in range
*/

static IWString guard_fingerprint_tag;
static similarity_type_t guard_fingerprint_distance = 0.0;

/*
  Fingerprints can be read in as either
    Daylight compressed form,
    Hex form (uppercase or lowercase)
    ASCII 0/1 form
    Sparse form <2,5-8,99-166;1028>
    Non colliding sparse form, Daylight encoded
    Sparse SVM-Lite form  b:c b:c
*/

#define DAYLIGHT_COMPRESSED_FORM 0
#define UPPERCASE_HEX_FORM 1
#define LOWERCASE_HEX_FORM 2
#define ASCII_01_FORM 3
#define SPARSE_FORM 4
#define NON_COLLIDING_FORM 5
#define FIXED_COUNTED_FORM 6
#define CONVERTED_FROM_SPARSE 7
#define MULTICONFORMER_01_FORM 8
#define MULTICONFORMER_FSC_FORM 9
#define MULTICONFORMER_SPARSE_FORM 10
#define SPARSE_ASCII_FORM 11
#define SPARSE_SVMLITE_FORM 12

static int * _fingerprint_type = NULL;

static int * _sparse_fingerprint_type = NULL;

static int
set_nfingerprints (int nfp, int sfp, int cfp,
                   int mc01,     // multiconformer01
                   int mcfsc,    // multiconformer fixed size counted
                   int mcsparse) // multiconformer, sparse
{
  assert (nfp > 0 || sfp > 0 || cfp > 0 || mc01 > 0 || mcfsc > 0 || mcsparse > 0);

  int total_fingerprints = nfp + sfp + cfp + mc01 + mcfsc + mcsparse;

  int nfp_to_allocate = nfp;
  if (convert_sparse_fingerprint_to_bits)
    nfp_to_allocate += sfp;

  int cfp_to_allocate = cfp;
  if (convert_sparse_fingerprint_to_fixed_width_counted)
    cfp_to_allocate += sfp;

  if (nfp_to_allocate > 0)
  {
    _fingerprint_weight = new float[nfp_to_allocate];
    _bits_in_fingerprint = new int[nfp_to_allocate];
    _fingerprint_type = new int[nfp_to_allocate];
    _dense_fingerprint_distance_metric = new int[nfp_to_allocate];
    _fingerprint_tag = new IWString[nfp_to_allocate];

    if (NULL == _fingerprint_tag)
    {
      cerr << "set_nfingerprints: " << nfp_to_allocate << " memory failure\n";
      return 0;
    }

    for (int i = 0; i < nfp_to_allocate; i++)
    {
      _fingerprint_weight[i] = 1.0 / static_cast<float>(total_fingerprints);
      _bits_in_fingerprint[i] = -1;
      _fingerprint_type[i] = DAYLIGHT_COMPRESSED_FORM;
      _dense_fingerprint_distance_metric[i] = default_dense_fingerprint_distance_metric;
    }

    bit_fingerprints_in_file = nfp;

    fingerprints_to_allocate = bit_fingerprints_in_file;

    number_fingerprints_to_use_in_computations = bit_fingerprints_in_file;
  }

  if (cfp_to_allocate > 0)
  {
    _fixed_size_counted_fingerprint_weight = new float[cfp_to_allocate];
    _fixed_size_counted_fingerprint_distance_metric = new int[cfp_to_allocate];
    _fixed_size_counted_fingerprint_tag = new IWString[cfp_to_allocate];
    if (NULL == _fixed_size_counted_fingerprint_tag)
    {
      cerr << "set_nfingerprints: " << cfp_to_allocate << " fixed size counted fingerprints, memory failure\n";
      return 0;
    }

    for (int i = 0; i < cfp_to_allocate; i++)
    {
      _fixed_size_counted_fingerprint_weight[i] = 1.0 / static_cast<float>(total_fingerprints);
      _fixed_size_counted_fingerprint_distance_metric[i] = default_fixed_size_counted_fingerprint_distance_metric;
    }

    number_fixed_size_counted_fingerprints_in_file = cfp;

    number_fixed_size_counted_fingerprints_to_use_in_computations = cfp;   // for now
  }

  if (sfp > 0)
  {
    _sparse_fingerprint_weight = new float[sfp];
    _sparse_fingerprint_scaling_factor = new float[sfp];
    _sparse_fingerprint_distance_metric = new int[sfp];
    _sparse_fingerprint_tag = new IWString[sfp];
    _sparse_fingerprint_type = new int[sfp];
    if (NULL == _sparse_fingerprint_tag)
    {
      cerr << "set_nfingerprints: " << sfp << " sparse fingerprints, memory failure\n";
      return 0;
    }

    for (int i = 0; i < sfp; i++)
    {
      _sparse_fingerprint_weight[i] = 1.0 / static_cast<float>(total_fingerprints);
      _sparse_fingerprint_distance_metric[i] = default_sparse_fingerprint_distance_metric;
      _sparse_fingerprint_scaling_factor[i] = 1.0;
      _sparse_fingerprint_type[i] = DAYLIGHT_COMPRESSED_FORM;
    }

    _number_sparse_fingerprints = sfp;

    _number_sparse_fingerprints_to_use_in_computations = _number_sparse_fingerprints;

    if (convert_sparse_fingerprint_to_bits)
    {
      fingerprints_to_allocate += sfp;
      number_fingerprints_to_use_in_computations += sfp;
      _number_sparse_fingerprints_to_use_in_computations = 0;

      for (int i = bit_fingerprints_in_file; i < fingerprints_to_allocate; i++)
      {
        _fingerprint_type[i] = CONVERTED_FROM_SPARSE;
      }
    }
    else if (convert_sparse_fingerprint_to_fixed_width_counted)
    {
      number_fixed_size_counted_fingerprints_to_allocate += sfp;
      number_fixed_size_counted_fingerprints_to_use_in_computations += sfp;
      _number_sparse_fingerprints_to_use_in_computations = 0;
    }
  }

  if (mc01 > 0)
  {
    _multiconformer_01_weight = new float[mc01];
    _multiconformer_01_distance_metric = new int[mc01];
    _multiconformer_01_fingerprint_tag = new IWString[mc01];
    if (NULL == _multiconformer_01_fingerprint_tag)
    {
      cerr << "set_nfingerprints: " << mc01 << " multiconformer 0/1 fingerprints, memory failure\n";
      return 0;
    }

    for (int i = 0; i < mc01; i++)
    {
      _multiconformer_01_weight[i] = 1.0 / float(total_fingerprints);
      _multiconformer_01_distance_metric[i] = default_multiconformer_01_distance_metric;
    }

    number_multiconformer_01_fingerprints_in_file = mc01;

    multiconformer_fingerprints_present = 1;
  }

  if (mcfsc > 0)
  {
    _multiconformer_fixed_size_counted_weight = new float[mcfsc];
    _multiconformer_fixed_size_counted_distance_metric = new int[mcfsc];
    _multiconformer_fixed_size_counted_fingerprint_tag = new IWString[mcfsc];
    if (NULL == _multiconformer_fixed_size_counted_fingerprint_tag)
    {
      cerr << "set_nfingerprints: " << mcfsc << " multiconformer fixed size counted fingerprints, memory failure\n";
      return 0;
    }

    for (int i = 0; i < mcfsc; i++)
    {
      _multiconformer_fixed_size_counted_weight[i] = 1.0 / float(total_fingerprints);
      _multiconformer_fixed_size_counted_distance_metric[i] = default_multiconformer_fixed_size_counted_distance_metric;
    }

    number_multiconformer_fixed_size_counted_fingerprints_in_file = mcfsc;

    multiconformer_fingerprints_present = 1;
  }

  if (mcsparse > 0)
  {
    _multiconformer_sparse_weight = new float[mcsparse];
    _multiconformer_sparse_distance_metric = new int[mcsparse];
    _multiconformer_sparse_fingerprint_tag = new IWString[mcsparse];
    if (NULL == _multiconformer_sparse_fingerprint_tag)
    {
      cerr << "set_nfingerprints: " << mcsparse << " multiconformer sparse fingerprints, memory failure\n";
      return 0;
    }

    for (int i = 0; i < mcsparse; i++)
    {
      _multiconformer_sparse_weight[i] = 1.0 / float(total_fingerprints);
      _multiconformer_sparse_distance_metric[i] = default_multiconformer_sparse_distance_metric;
    }

    number_multiconformer_sparse_fingerprints_in_file = mcsparse;

    multiconformer_fingerprints_present = 1;
  }

  return 1;
}

int
number_fingerprints()
{
  return bit_fingerprints_in_file;
}

int
number_sparse_fingerprints()
{
  return _number_sparse_fingerprints;
}

static IWString identifier_tag("PCN<");

void
set_identifier_tag (const_IWSubstring const & id)
{
  identifier_tag = id;

  if (! identifier_tag.ends_with('<'))
    identifier_tag += '<';

  return;
}

/*
  This regular expression is a generic form for the basic fingerprint tag
  Many problems trying to initialise it...
*/

static IW_Regular_Expression fingerprint_rx; // ("^FP[A-Z]*<.+[0-9]>");

//static IWString fingerprint_rx_default ("^(FP[A-Z,0-9]+|HX[A-Z,0-9]+|NC[A-Z,0-9]+)<");
static IWString fingerprint_rx_default("^(FP[A-Z,0-9]+|HX[A-Z,0-9]+|NC[A-Z,0-9]+|FC[A-Z,0-9]+|MC[A-Z,0-9]+)<");

/*
  Molecular properties may or may not be present. The default is present
*/

static int dash_P_specified = 0;

static IWString molecular_property_integer_tag("MPR");
static IWString molecular_property_continuous_tag;

static float _property_weight_integer = -1.0;
static float _property_weight_continuous = -1.0;

const IWString &
property_tag ()
{
  return molecular_property_integer_tag;
}


void
set_property_weight_integer (double w)
{
  _property_weight_integer = w;

  return;
}

int
set_fingerprint_regular_expression (const const_IWSubstring & rx)
{
  if (! rx.starts_with('^'))
    cerr << "set_fingerprint_regular_expression: WARNING, rx does not start with '^' :" << rx << endl;

  return fingerprint_rx.set_pattern(rx);
}

static int
print_tags_and_weights (std::ostream & os)
{
  for (int i = 0; i < number_fingerprints_to_use_in_computations; i++)
  {
    os << "Fingerprint " << i << " '" << _fingerprint_tag[i] << "' weight " << _fingerprint_weight[i];

    if (DAYLIGHT_COMPRESSED_FORM == _fingerprint_type[i])
      ;
    else if (LOWERCASE_HEX_FORM == _fingerprint_type[i])
      os << " hex";
    else if (UPPERCASE_HEX_FORM == _fingerprint_type[i])
      os << " HEX";
    else if (ASCII_01_FORM == _fingerprint_type[i])
      os << " ascii";
    else if (SPARSE_FORM == _fingerprint_type[i])
      os << " sparse";
    else if (CONVERTED_FROM_SPARSE == _fingerprint_type[i])
      os << " from sparse";
    else if (0 == _fingerprint_type[i])
      ;
    else
      os << " unknown type " << _fingerprint_type[i];
    os << endl;
  }

  if (convert_sparse_fingerprint_to_fixed_width_counted || convert_sparse_fingerprint_to_bits)
    ;
  else
  {
    for (int i = 0; i < _number_sparse_fingerprints; i++)
    {
      os << "Non colliding fingerprint " << i << " '" << _sparse_fingerprint_tag[i] << "' weight " << _sparse_fingerprint_weight[i];
      if (DAYLIGHT_COMPRESSED_FORM == _sparse_fingerprint_type[i])
        ;
      else if (SPARSE_ASCII_FORM == _sparse_fingerprint_type[i])
        os << " ascii";
      os << '\n';
    }
  }

  for (int i = 0; i < number_fixed_size_counted_fingerprints_to_use_in_computations; i++)
  {
    os << "Fixed counted fingerprint " << i << " '" << _fixed_size_counted_fingerprint_tag[i] << "' weight " << _fixed_size_counted_fingerprint_weight[i] << endl;
  }

  for (int i = 0; i < number_multiconformer_01_fingerprints_in_file; i++)
  {
    os << "Multiconformer 01 " << i << " '" << _multiconformer_01_fingerprint_tag[i] << "' weight " << _multiconformer_01_weight[i] << endl;
  }

  for (int i = 0; i < number_multiconformer_fixed_size_counted_fingerprints_in_file; i++)
  {
    os << "Multiconformer FSC " << i << " '" << _multiconformer_fixed_size_counted_fingerprint_tag[i] << "' weight " << _multiconformer_fixed_size_counted_weight[i] << endl;
  }

  for (int i = 0; i < number_multiconformer_sparse_fingerprints_in_file; i++)
  {
    os << "Multiconformer SPARSE " << i << " '" << _multiconformer_sparse_fingerprint_tag[i] << "' weight " << _multiconformer_sparse_weight[i] << endl;
  }

  if (molecular_property_integer_tag.length())
    os << "Property weight " << _property_weight_integer << endl;

  if (molecular_property_continuous_tag.length())
    os << "Continuous property weight " << _property_weight_continuous << endl;

  return os.good();
}

static float individual_property_scaling_factor = 1.0;

static int file_scope_number_integer_molecular_properties = 0;

/*
  When someone enters 'desc=all', we need a special value for
  number_continuous_molecular_properties that indicates take all
  descriptors
*/

#define ALL_CONTINUOUS_MOLECULAR_PROPERTIES -1

static int number_continuous_molecular_properties = 0;

int
set_number_integer_molecular_properties (int n)
{
  assert (n > 0);

  file_scope_number_integer_molecular_properties = n;

  return 1;
}

int
number_integer_molecular_properties()
{
  return file_scope_number_integer_molecular_properties;
}

int
set_number_continuous_molecular_properties (int n)
{
  assert (n > 0);

  number_continuous_molecular_properties = n;

  return 1;
}

template <typename T>
Molecular_Properties<T>::~Molecular_Properties()
{
  if (NULL != _property)
    delete [] _property;

  return;
}

Molecular_Properties_Continuous::Molecular_Properties_Continuous()
{
  if (number_continuous_molecular_properties <= 0)
  {
    _nproperties = 0;
    _property = NULL;
  }
  else
  {
    _nproperties = number_continuous_molecular_properties;
    _property = new iwproperty_t[_nproperties];

    if (NULL == _property)
      cerr << "Molecular_Properties_Continuous::cannot allocate space for " << _nproperties << " properties\n";
  }

  return;
}

int
Molecular_Properties_Continuous::construct_from_descriptor_record (const const_IWSubstring & buffer)
{
  const_IWSubstring tmp = buffer;

  if (tmp.ends_with('>'))
    tmp.chop();
  tmp.remove_up_to_first('<');

  if (ALL_CONTINUOUS_MOLECULAR_PROPERTIES == number_continuous_molecular_properties)
  {
    number_continuous_molecular_properties = tmp.nwords();
    if (0 == number_continuous_molecular_properties)
    {
      cerr << "molecular_properties::construct_from_descriptor_record: no descriptors!\n";
      return 0;
    }

    cerr << "Molecular_Properties::construct_from_descriptor_record: " << number_continuous_molecular_properties << " descriptors in input\n";
  }

  _nproperties = number_continuous_molecular_properties;

  if (NULL == _property)
    _property = new iwproperty_t[_nproperties];

  int j = 0;      // index into string
  for (int i = 0; i < number_continuous_molecular_properties; i++)
  {
    const_IWSubstring token;
    if (! tmp.nextword(token, j))
    {
      cerr << "Molecular_Properties::construct_from_descriptor_record: insufficient tokens, i = " << i << "\n";
      cerr << tmp << endl;
      return 0;
    }

    if (! token.numeric_value(_property[i]))
    {
      cerr << "Molecular_Properties::construct_from_descriptor_record: invalid token '" << token << "'\n";
      return 0;
    }

    _property[i] = _property[i] * individual_property_scaling_factor;
  }

  return 1;
}

Molecular_Properties_Integer::Molecular_Properties_Integer()
{
  if (0 == file_scope_number_integer_molecular_properties)
  {
    _nproperties = 0;
    _property = NULL;
  }
  else
  {
    _nproperties = file_scope_number_integer_molecular_properties;
    _property = new int[_nproperties];

    if (NULL == _property)
      cerr << "Molecular_Properties_Integer:: cannot allocate " << _nproperties << " properties\n";
  }

  return;
}

int
Molecular_Properties_Integer::construct_from_tdt_fp_record (const const_IWSubstring & buffer)
{
  IWDYFP fp;

  if (! fp.construct_from_tdt_record(buffer))
  {
    cerr << "Molecular_Properties::construct_from_tdt_record: cannot parse '" << buffer << "' as FP\n";
    return 0;
  }

  if (0 == _nproperties)
  {
    _nproperties = fp.nbits() / IW_BITS_PER_BYTE;

    set_number_integer_molecular_properties(_nproperties);
  }
  else
  {
    assert (fp.nbits() / IW_BITS_PER_BYTE == _nproperties);
  }

  if (NULL == _property)
    _property = new int[_nproperties];

//const unsigned char * b = static_cast<const unsigned char *> (fp.bits());
  const unsigned char * b = (const unsigned char *) fp.bits();

  for (int i = 0; i < _nproperties; i++)
  {
    _property[i] = static_cast<int>(b[i]);
  }

  return 1;
}

int
Molecular_Properties_Integer::natoms() const
{
  assert (NULL != _property);

  return static_cast<int>(_property[0]);
}

int
Molecular_Properties_Integer::nrings() const
{
  assert (NULL != _property);

  return static_cast<int>(_property[1]);
}

int
Molecular_Properties_Integer::aromatic_atoms() const
{
  assert (NULL != _property);

  return static_cast<int>(_property[4]);
}

template <typename T>
Molecular_Properties<T> &
Molecular_Properties<T>::operator = (const Molecular_Properties<T> & rhs)
{
  if (NULL == rhs._property)
  {
    if (NULL != _property)
    {
      delete _property;
      _property = NULL;
    }

    return * this;
  }

  if (0 == _nproperties)
    return *this;

  if (NULL == _property)
    _property = new T[_nproperties];

// could just use a for loop here...

  memcpy(_property, rhs._property, _nproperties * sizeof(T));

  return *this;
}

template class Molecular_Properties<int>;
template class Molecular_Properties<float>;

/*
  We can greatly speed things up by pre-computing ratios
  Hmmm, bad idea, we use properties to store descriptors, so they need to
  be floats.

  External linkage
*/

similarity_type_t * precomputed_ratio = NULL;

int
initialise_properties_ratios()
{
  if (NULL != precomputed_ratio)    // must have been set up already
    return 1;

  precomputed_ratio = new similarity_type_t[256*256];

  if (NULL == precomputed_ratio)
  {
    cerr << "Yipes, cannot allocate the properties ratio array\n";
    return 0;
  }

  for (int i = 0; i < 256; i++)
  {
    similarity_type_t si = static_cast<similarity_type_t>(i);
    for (int j = 0; j < 256; j++)
    {
      if (i == j)
      {
        precomputed_ratio[i * 256 + j] = 1.0;
      }
      else if (0 == i || 0 == j)
      {
        precomputed_ratio[i * 256 + j] = 0.5;
        precomputed_ratio[j * 256 + i] = 0.5;
      }
      else
      {
        similarity_type_t sj = static_cast<similarity_type_t>(j);
        if (j > i)
        {
          precomputed_ratio[i * 256 + j] = si / sj;
          precomputed_ratio[j * 256 + i] = si / sj;
        }
        else
        {
          precomputed_ratio[i * 256 + j] = sj / si;
          precomputed_ratio[j * 256 + i] = sj / si;
        }
      }
    }
  }

  for (int i = 0; i < 256 * 256; i++)
  {
    precomputed_ratio[i] = precomputed_ratio[i] / static_cast<similarity_type_t>(8.0);
  }

  return 1;
}

class Ratio_Deallocator
{
  private:
  public:
    ~Ratio_Deallocator();
};

Ratio_Deallocator::~Ratio_Deallocator()
{
  if (NULL == precomputed_ratio)
    return;

  delete [] precomputed_ratio;
  precomputed_ratio = NULL;

  return;
}

static Ratio_Deallocator rd;

similarity_type_t
Molecular_Properties_Integer::similarity (const Molecular_Properties_Integer & rhs) const
{
  if (INTEGER_PROPERTY_DISTANCE_METRIC_RATIO == _integer_property_distance_metric)
    ;
  else if (INTEGER_PROPERTY_DISTANCE_METRIC_DICE == _integer_property_distance_metric)
    return dice_coefficient(rhs);

/*assert (NULL != precomputed_ratio);*/

  similarity_type_t rc = static_cast<similarity_type_t>(0.0);
  
  int j;
  for (int i = 0; i < file_scope_number_integer_molecular_properties; i++)
  {
    j = _property[i] * 256 + rhs._property[i];

    rc += precomputed_ratio[j];
  }

/*assert (rc >= 0.0 && rc <= 1.0);*/

  return rc;
}

#ifdef VERSION_WITHOUT_PRECOMPUTED_RATIOS

similarity_type_t
Molecular_Properties_Integer::similarity (const Molecular_Properties_Integer & rhs) const
{
  if (INTEGER_PROPERTY_DISTANCE_METRIC_RATIO == _integer_property_distance_metric)
    ;
  else if (INTEGER_PROPERTY_DISTANCE_METRIC_DICE == _integer_property_distance_metric)
    return dice_coefficient(rhs);

  similarity_type_t rc = static_cast<similarity_type_t>(0.0);
  
  int p1, p2;

  for (int i = 0; i < file_scope_number_integer_molecular_properties; i++)
  {
    p1 = _property[i];
    p2 = rhs._property[i];

    if (p1 == p2)
      rc += static_cast<similarity_type_t>(1.0);
    else if (0.0 == p1 || 0.0 == p2)
      rc += 0.5;
    else if (p1 > p2)
      rc += static_cast<float>(p2) / static_cast<float>(p1);
    else
      rc += static_cast<float>(p1) / static_cast<float>(p2);
//  cerr << " property " << i << " values " << _property[i] << " and " << rhs._property[i] << " precomputed " << precomputed_ratio[j] << endl;
  }

  return rc / static_cast<float>(file_scope_number_integer_molecular_properties);
}

#endif


similarity_type_t
Molecular_Properties_Integer::dice_coefficient (const Molecular_Properties_Integer & rhs) const
{
  int bits_in_common = 0;
  int na = 0;
  int nb = 0;

  for (int i = 0; i < file_scope_number_integer_molecular_properties; i++)
  {
    if (_property[i] > rhs._property[i])
      bits_in_common += rhs._property[i];
    else
      bits_in_common += _property[i];

    na += _property[i];
    na += rhs._property[i];
  }

  return static_cast<float>(bits_in_common) / static_cast<float> (na + nb - bits_in_common);
}

int
Molecular_Properties_Integer::bits_in_common (const Molecular_Properties_Integer & rhs) const
{
  int bic = 0;

  for (int i = 0; i < _nproperties; i++)
  {
    if (_property[i] > rhs._property[i])
      bic += rhs._property[i];
    else
      bic += _property[i];
  }

  return bic;
}

similarity_type_t
Molecular_Properties_Continuous::similarity (const Molecular_Properties_Continuous & rhs) const
{
  assert (NULL != _property && NULL != rhs._property);

#ifdef ECHO_PROPERTIES
  for (int i = 0; i < _nproperties; i++)
  {
    cerr << " i = " << i << ' ' << _property[i] << " " << rhs._property[i] << endl;
  }
#endif

  if (CONTINUOUS_PROPERTY_DISTANCE_METRIC_CART == _continuous_property_distance_metric)
    return static_cast<similarity_type_t>(1.0) - cartesian_distance(rhs);
  if (CONTINUOUS_PROPERTY_DISTANCE_METRIC_EXP1 == _continuous_property_distance_metric)
    return static_cast<similarity_type_t>(1.0) - exp1_distance(rhs);
  if (CONTINUOUS_PROPERTY_DISTANCE_METRIC_EXP2 == _continuous_property_distance_metric)
    return static_cast<similarity_type_t>(1.0) - exp2_distance(rhs);
  if (CONTINUOUS_PROPERTY_DISTANCE_METRIC_DICE == _continuous_property_distance_metric)
    return static_cast<similarity_type_t>(1.0) - dice_coefficient(rhs);

  similarity_type_t rc;

  if (_property[0] > rhs._property[0])
    rc = rhs._property[0] / _property[0];
  else if (_property[0] < rhs._property[0])
    rc = _property[0] / rhs._property[0];
  else
    rc = static_cast<similarity_type_t>(1.0);

  for (int i = 1; i < _nproperties; i++)
  {
    if (_property[i] == rhs._property[i])
      rc += static_cast<similarity_type_t>(1.0);
    else if (static_cast<float>(0.0) == _property[i] || static_cast<float>(0.0) == rhs._property[i])
      rc += static_cast<similarity_type_t>(0.5);
    else if (_property[i] > rhs._property[i])
      rc += rhs._property[i] / _property[i];
    else
      rc += _property[i] / rhs._property[i];
  }

  rc = rc / static_cast<similarity_type_t>(_nproperties);

  if (rc < static_cast<similarity_type_t>(0.0) || rc > static_cast<similarity_type_t>(1.0))
  {
    cerr << "Molecular_Properties::distance: properties sum out of range " << rc << endl;
    for (int i = 0; i < _nproperties; i++)
    {
      cerr << "i = " << i << " lhs " << _property[i] << " rhs " << rhs._property[i];
      if (0.0 == _property[i] || 0.0 == rhs._property[i])
        cerr << " zero";
      cerr << endl;
    }

    if (rc < static_cast<similarity_type_t>(0.0))
      rc = static_cast<similarity_type_t>(0.0);
    else if (rc > static_cast<similarity_type_t>(1.0))
      rc = static_cast<similarity_type_t>(1.0);
  }

  return rc;
}

similarity_type_t
Molecular_Properties_Continuous::cartesian_distance (const Molecular_Properties_Continuous & rhs) const
{
  similarity_type_t rc = static_cast<similarity_type_t>(0.0);

  for (int i = 0; i < _nproperties; i++)
  {
    rc += (_property[i] - rhs._property[i]) * (_property[i] - rhs._property[i]);
//  cerr << _property[i] << " and " << rhs._property[i] << " rc = " << rc << endl;
  }

//cerr << "Sum over " << _nproperties << " properties, rc = " << rc << endl;

  return sqrt(rc);
}

similarity_type_t
Molecular_Properties_Continuous::exp1_distance (const Molecular_Properties_Continuous & rhs) const
{
  double rc = 0.0;

  for (int i = 0; i < _nproperties; i++)
  {
    rc += exp(continuous_property_distance_metric_exp1_value * fabs(_property[i] - rhs._property[i]));
  }

  return static_cast<similarity_type_t>(1.0 - rc / _nproperties);
}

similarity_type_t
Molecular_Properties_Continuous::exp2_distance (const Molecular_Properties_Continuous & rhs) const
{
  double rc = 0.0;

  for (int i = 0; i < _nproperties; i++)
  {
    rc += exp(continuous_property_distance_metric_exp2_value * fabs((_property[i] - rhs._property[i]) * (_property[i] - rhs._property[i])));
  }

  return static_cast<similarity_type_t>(1.0 - rc / _nproperties);
}

similarity_type_t
Molecular_Properties_Continuous::dice_coefficient (const Molecular_Properties_Continuous & rhs) const
{
  float bits_in_common = static_cast<float>(0.0);
  float na = static_cast<float>(0.0);
  float nb = static_cast<float>(0.0);

  for (int i = 0; i < _nproperties; i++)
  {
    if (_property[i] > rhs._property[i])
      bits_in_common += rhs._property[i];
    else
      bits_in_common += _property[i];

    na += _property[i];
    nb += rhs._property[i];
  }

  return bits_in_common / (na + nb - bits_in_common);
}

IWDYFP &
IW_General_Fingerprint::operator [] (int i) const
{
  assert (i >= 0 && i < bit_fingerprints_in_file);

  return _fingerprint[i];
}

IWDYFP &
IW_General_Fingerprint::item (int i) const
{
  assert (i >= 0 && i < bit_fingerprints_in_file);

  return _fingerprint[i];
}

int
IW_General_Fingerprint::nfingerprints() const
{
  return bit_fingerprints_in_file;
}

IW_General_Fingerprint &
IW_General_Fingerprint::operator = (const IW_General_Fingerprint & rhs)
{
  int need_to_allocate = 0;

  if (number_fingerprints_to_use_in_computations > 0 && NULL == _fingerprint)
    need_to_allocate = 1;
  else if (_number_sparse_fingerprints > 0 && NULL == _sparse_fingerprint)
    need_to_allocate = 1;

// should check other forms...

  if (need_to_allocate)
    _allocate_fingerprint_array();

  for (int i = 0; i < bit_fingerprints_in_file; i++)
  {
    _fingerprint[i] = rhs._fingerprint[i];
  }

  for (int i = 0; i < _number_sparse_fingerprints; i++)
  {
    _sparse_fingerprint[i] = rhs._sparse_fingerprint[i];
  }

  _molecular_properties_integer = rhs._molecular_properties_integer;
  _molecular_properties_continuous = rhs._molecular_properties_continuous;

  for (int i = 0; i < number_multiconformer_fixed_size_counted_fingerprints_in_file; i++)
  {
    _multiconformer_fixed_size_counted[i] = rhs._multiconformer_fixed_size_counted[i];
  }

  for (int i = 0; i < number_multiconformer_01_fingerprints_in_file; i++)
  {
    _multiconformer_01[i] = rhs._multiconformer_01[i];
  }

  for (int i = 0; i < number_multiconformer_sparse_fingerprints_in_file; i++)
  {
    _multiconformer_sparse[i] = rhs._multiconformer_sparse[i];
  }

  _id = rhs._id;

  off_t o;
  if (rhs._offset.value(o))
    _offset.set(o);
  else
    _offset.unset();

  return *this;
}

int
IW_General_Fingerprint::ok() const
{
  return _id.ok();
}

void
IW_General_Fingerprint::flip()
{
  for (int i = 0; i < bit_fingerprints_in_file; i++)
  {
    _fingerprint[i].flip();
  }

  return;
}

//#define DEBUG_TANIMOTO

similarity_type_t
IW_General_Fingerprint::tanimoto(IW_General_Fingerprint & rhs)
{
#ifdef DEBUG_TANIMOTO
  cerr << "Between " << id() << " and " << rhs.id() << endl;
#endif

  similarity_type_t result = 0;

  for (int i = 0; i < number_fingerprints_to_use_in_computations; i++)
  {
    IWDYFP & fp2 = rhs._fingerprint[i];

#ifdef DEBUG_TANIMOTO
    cerr << "Metric " << _dense_fingerprint_distance_metric[i] << " weight " << _fingerprint_weight[i] << endl;
#endif

    if (DENSE_DISTANCE_METRIC_TANIMOTO == _dense_fingerprint_distance_metric[i])
      result += _fingerprint_weight[i] * _fingerprint[i].tanimoto(fp2);
    else if (FVB_MODIFIED_TANIMOTO == _dense_fingerprint_distance_metric[i])
      result += _fingerprint_weight[i] * _fingerprint[i].fvb_modified_tanimoto(fp2);
    else if (RUSSEL_RAO == _dense_fingerprint_distance_metric[i])
      result += _fingerprint_weight[i] * _fingerprint[i].russel_rao(fp2);
    else if (FORBES_SIMILARITY == _dense_fingerprint_distance_metric[i])
      result += _fingerprint_weight[i] * _fingerprint[i].forbes_similarity(fp2);
    else if (SIMPLE_MATCHING == _dense_fingerprint_distance_metric[i])
      result += _fingerprint_weight[i] * _fingerprint[i].simple_matching(fp2);
    else if (DICE_DISTANCE_METRIC == _dense_fingerprint_distance_metric[i])
      result += _fingerprint_weight[i] * _fingerprint[i].sorensendice(fp2);
    else if (DENSE_DISTNANCE_OVERLAP == _dense_fingerprint_distance_metric[i])
      result += _fingerprint_weight[i] * _fingerprint[i].overlap(fp2);

#ifdef DEBUG_TANIMOTO
    cerr << _fingerprint[i].nbits() << " bits, " << _fingerprint[i].nset() << " bits set and " << fp2.nset() << endl;
    _fingerprint[i].printon(cerr);
    cerr << endl;
    fp2.printon(cerr);
    cerr << endl;
    cerr << " FP " << i << " similarity " <<  _fingerprint[i].tanimoto(fp2) << " weight " << _fingerprint_weight[i] << " bic " << _fingerprint[i].bits_in_common (rhs._fingerprint[i]) << endl;
#endif

#ifdef DEBUG_TANIMOTO
    cerr << "Fingerprint component " << i << " " << (_fingerprint[i].tanimoto(fp2) * _fingerprint_weight[i]) << ", result so far " << result << endl;
#endif
  }

  for (int i = 0; i < _number_sparse_fingerprints_to_use_in_computations; i++)
  {
    similarity_type_t tmp;

    if (SPARSE_DISTANCE_METRIC_TANIMOTO == _sparse_fingerprint_distance_metric[i])
    {
      if (sparse_fingerprint_counts_limited)
        tmp = _sparse_fingerprint_weight[i] * _sparse_fingerprint[i].tanimoto(rhs._sparse_fingerprint[i]);
      else
        tmp = _sparse_fingerprint_weight[i] * _sparse_fingerprint[i].tanimoto_with_unlimited_counts(rhs._sparse_fingerprint[i]);
    }
    else if (SPARSE_DISTANCE_METRIC_BINARY == _sparse_fingerprint_distance_metric[i])
      tmp = _sparse_fingerprint_weight[i] * _sparse_fingerprint[i].tanimoto_binary(rhs._sparse_fingerprint[i]);
    else if (FVB_MODIFIED_TANIMOTO == _sparse_fingerprint_distance_metric[i])
      tmp = _sparse_fingerprint_weight[i] * _sparse_fingerprint[i].fvb_modified_tanimoto(rhs._sparse_fingerprint[i]);
    else if (MANHATTAN_DISTANCE_METRIC == _sparse_fingerprint_distance_metric[i])
      tmp = _sparse_fingerprint_weight[i] * _sparse_fingerprint[i].manhattan_distance(rhs._sparse_fingerprint[i]);
    else if (SOERGEL_DISTANCE_METRIC == _sparse_fingerprint_distance_metric[i])
      tmp = _sparse_fingerprint_weight[i] * _sparse_fingerprint[i].soergel_similarity(rhs._sparse_fingerprint[i]);
    else if (SOERGEL_VARIANT_DISTANCE_METRIC == _sparse_fingerprint_distance_metric[i])
      tmp = _sparse_fingerprint_weight[i] * _sparse_fingerprint[i].soergel_variant_similarity(rhs._sparse_fingerprint[i]);
    else if (CTAN_DISTANCE_METRIC == _sparse_fingerprint_distance_metric[i])
      tmp = _sparse_fingerprint_weight[i] * _sparse_fingerprint[i].continuous_tanimoto(rhs._sparse_fingerprint[i]);
    else if (DICE_DISTANCE_METRIC == _sparse_fingerprint_distance_metric[i])
      tmp = _sparse_fingerprint_weight[i] * _sparse_fingerprint[i].cosine_coefficient(rhs._sparse_fingerprint[i]);
    else
      tmp = _sparse_fingerprint_weight[i] * _sparse_fingerprint[i].cosine_measure(rhs._sparse_fingerprint[i]);

    if (1.0 != _sparse_fingerprint_scaling_factor[i])
      tmp = 1.0 - _sparse_fingerprint_scaling_factor[i] * (1.0 - tmp);

    result += tmp;

#ifdef DEBUG_TANIMOTO
    cerr << "sparse fingerprint component " << i << ", weight " << _sparse_fingerprint_weight[i] << ", value " << tmp << " total " << result <<endl;
#endif
  }

  for (int i = 0; i < number_fixed_size_counted_fingerprints_to_use_in_computations; i++)
  {
    if (SPARSE_DISTANCE_METRIC_TANIMOTO == _fixed_size_counted_fingerprint_distance_metric[i])
      result += _fixed_size_counted_fingerprint_weight[i] * _fixed_size_counted_fingerprint[i].tanimoto(rhs._fixed_size_counted_fingerprint[i]);
    else if (FVB_MODIFIED_TANIMOTO == _fixed_size_counted_fingerprint_distance_metric[i])
      result += _fixed_size_counted_fingerprint_weight[i] * _fixed_size_counted_fingerprint[i].fvb_modified_tanimoto(rhs._fixed_size_counted_fingerprint[i]);
    else
    {
      cerr << "IW_General_Fingerprint::unrecognised fixed counted distance metric " << _fixed_size_counted_fingerprint_distance_metric[i] << endl;
      abort();
    }
#ifdef DEBUG_TANIMOTO
    cerr << "After fixed size counted fingerprint " << i << " result " << result <<endl;
#endif
  }

#ifdef DEBUG_TANIMOTO
  cerr << "Any multiconformer_fingerprints_present " << multiconformer_fingerprints_present << ", result so far " << result << endl;
#endif

  if (multiconformer_fingerprints_present)
  {
    for (int i = 0; i < number_multiconformer_01_fingerprints_in_file; i++)
    {
      if (multiconformer_use_average_similarity)
        result += _multiconformer_01[i].average_tanimoto(rhs._multiconformer_01[i]) * _multiconformer_01_weight[i];
      else
        result += _multiconformer_01[i].tanimoto(rhs._multiconformer_01[i]) * _multiconformer_01_weight[i];
#ifdef DEBUG_TANIMOTO
      cerr << "After multiconformer 01 " << result << endl;
#endif
    }

    for (int i = 0; i < number_multiconformer_fixed_size_counted_fingerprints_in_file; i++)
    {
      result += _multiconformer_fixed_size_counted[i].tanimoto(rhs._multiconformer_fixed_size_counted[i]) * _multiconformer_fixed_size_counted_weight[i];
    }

    for (int i = 0; i < number_multiconformer_sparse_fingerprints_in_file; i++)
    {
      result += _multiconformer_sparse[i].tanimoto(rhs._multiconformer_sparse[i]) * _multiconformer_sparse_weight[i];
#ifdef DEBUG_TANIMOTO
      cerr << "After multiconformer sparse " << i << " result " << result << ", w = " << _multiconformer_sparse_weight[i] << endl;
#endif
    }
  }

#ifdef FB_ENTROPY_WEIGHTED_FPS
  const int * p1 = _molecular_properties_integer.rawdata();
  const int * p2 = rhs._molecular_properties_integer.rawdata();
  result = result + _property_weight_integer * ffbwfw.tanimoto(p1, p2);
#else
  if (_molecular_properties_integer.active())
    result = result + _property_weight_integer * _molecular_properties_integer.similarity(rhs._molecular_properties_integer);

#ifdef DEBUG_TANIMOTO
  if (_molecular_properties_integer.active())
    cerr << "Molecular properties integer " << ( _property_weight_integer * _molecular_properties_integer.similarity (rhs._molecular_properties_integer)) << ", result so far " << result << endl;
#endif
#endif

  if (_molecular_properties_continuous.active())
    result = result + _property_weight_continuous * _molecular_properties_continuous.similarity (rhs._molecular_properties_continuous);

#ifdef DEBUG_TANIMOTO
  if (_molecular_properties_continuous.active())
    cerr << "Continuous molecular properties " << (_property_weight_continuous * _molecular_properties_continuous.similarity(rhs._molecular_properties_continuous)) << ", result so far " << result << endl;
#endif

#ifdef DEBUG_TANIMOTO
  cerr << "result " << result << endl;
#endif

// Jan 02. Sometimes we get numerical roundoff issues

  if (result > static_cast<similarity_type_t>(1.0))
  {
    if (result < static_cast<similarity_type_t>(1.01))
    {
//    cerr << "IW_General_Fingerprint::tanimoto: ignoring numeric roundoff " << (result - 1.0) << endl;   happens often, don't bother printing
      return static_cast<similarity_type_t>(1.0);
    }

    cerr << "IW_General_Fingerprint::tanimoto: invalid similarity " << result << " between '" << _id << "' and '" << rhs._id << "'\n";
    cerr << _property_weight_integer << " property weight\n";
    for (int i = 0; i < number_fingerprints_to_use_in_computations; i++)
    {
      cerr << " i = " << i << " fixed weight " << _fingerprint_weight[i] << endl;
    }
    for (int i = 0; i < _number_sparse_fingerprints_to_use_in_computations; i++)
    {
      cerr << " i = " << i << " sparse weight " << _sparse_fingerprint_weight[i] << endl;
    }
    return static_cast<similarity_type_t>(1.0);
  }


  return result;
}

//#define DEBUG_TVERSKY

similarity_type_t
IW_General_Fingerprint::tversky (IW_General_Fingerprint & rhs,
                                 const Tversky & t)
{
  if (t.optimistic_mode())
    return static_cast<similarity_type_t>(1.0) - optimistic_distance(rhs, t);

  similarity_type_t result = static_cast<similarity_type_t>(0.0);

  for (int i = 0; i < number_fingerprints_to_use_in_computations; i++)
  {
    IWDYFP & fp2 = rhs._fingerprint[i];

    result += _fingerprint[i].tversky(fp2, t) * _fingerprint_weight[i];
#ifdef DEBUG_TVERSKY
    cerr << "Tversky to component " << i << " value " <<  (_fingerprint[i].tversky(fp2, t) * _fingerprint_weight[i]) << endl;
#endif
  }

  for (int i = 0; i < number_fixed_size_counted_fingerprints_to_use_in_computations; i++)
  {
    result += _fixed_size_counted_fingerprint[i].tversky(rhs._fixed_size_counted_fingerprint[i], t) * _fixed_size_counted_fingerprint_weight[i];
  }

  for (int i = 0; i < _number_sparse_fingerprints_to_use_in_computations; i++)
  {
    result += _sparse_fingerprint[i].tversky(rhs._sparse_fingerprint[i], t) * _sparse_fingerprint_weight[i];
  }

  if (_molecular_properties_integer.active())
    result += _property_weight_integer * _molecular_properties_integer.similarity(rhs._molecular_properties_integer);

  if (_molecular_properties_continuous.active())
    result += _property_weight_continuous * _molecular_properties_continuous.similarity(rhs._molecular_properties_continuous);

#ifdef DEBUG_TVERSKY
  cerr << "Tversky between " << _id << " and " << rhs.id() << " tversky is " << result << endl;
#endif

  return result;
}

//#define DEBUG_OPTIMISTIC_DISTANCE

similarity_type_t
IW_General_Fingerprint::optimistic_distance (IW_General_Fingerprint & rhs,
                                             const Tversky & tv)
{
  similarity_type_t result = static_cast<similarity_type_t>(0.0);

  for (int i = 0; i < number_fingerprints_to_use_in_computations; i++)
  {
    IWDYFP & fp2 = rhs._fingerprint[i];

    result += _fingerprint[i].optimistic_distance(fp2, tv) * _fingerprint_weight[i];

#ifdef DEBUG_OPTIMISTIC_DISTANCE
    cerr << "To fingerprint component " << i << " value " <<  (_fingerprint[i].optimistic_distance(fp2, tv) * _fingerprint_weight[i]) << ", result so far " << result << endl;
#endif
  }

  for (int i = 0; i < _number_sparse_fingerprints_to_use_in_computations; i++)
  {
    result += _sparse_fingerprint[i].optimistic_distance(rhs._sparse_fingerprint[i], tv) * _sparse_fingerprint_weight[i];

#ifdef DEBUG_OPTIMISTIC_DISTANCE
    cerr << "Sparse component " << i << " dist " <<  (_sparse_fingerprint[i].optimistic_distance(rhs._sparse_fingerprint[i], tv) * _sparse_fingerprint_weight[i]) << ", result so far " << result << endl;
#endif
  }

  if (_molecular_properties_integer.active())
    result += _property_weight_integer * (static_cast<similarity_type_t>(1.0) - _molecular_properties_integer.similarity (rhs._molecular_properties_integer));

#ifdef DEBUG_OPTIMISTIC_DISTANCE
  if (_molecular_properties_integer.active())
    cerr << "Molecular properties integer " << ( _property_weight_integer * (static_cast<similarity_type_t>(1.0) - _molecular_properties_integer.similarity(rhs._molecular_properties_integer))) << ", result so far " << result << endl;
#endif

  if (_molecular_properties_continuous.active())
    result += _property_weight_continuous * (static_cast<similarity_type_t>(1.0) - _molecular_properties_continuous.similarity(rhs._molecular_properties_continuous));

#ifdef DEBUG_OPTIMISTIC_DISTANCE
  if (_molecular_properties_continuous.active())
    cerr << "Continuous molecular properties " << (_property_weight_continuous * (static_cast<similarity_type_t>(1.0) - _molecular_properties_continuous.similarity(rhs._molecular_properties_continuous))) << ", result so far " << result << endl;
#endif

#ifdef DEBUG_OPTIMISTIC_DISTANCE
  cerr << "Between " << _id << " and " << rhs.id() << " optimistic dist is " << result << endl;
#endif

  return result;
}

void
IW_General_Fingerprint::_allocate_fingerprint_array()
{
  if (fingerprints_to_allocate > 0)
  {
    _fingerprint = new IWDYFP [fingerprints_to_allocate];

    assert (NULL != _fingerprint);
  }

  if (_number_sparse_fingerprints > 0)
  {
    _sparse_fingerprint = new Sparse_Fingerprint[_number_sparse_fingerprints];

    assert (NULL != _sparse_fingerprint);
  }

  if (number_fixed_size_counted_fingerprints_to_allocate > 0)
  {
    _fixed_size_counted_fingerprint = new Fixed_Size_Counted_Fingerprint[number_fixed_size_counted_fingerprints_to_allocate];

    assert (NULL != _fixed_size_counted_fingerprint);
  }

  if (number_multiconformer_01_fingerprints_in_file > 0)
  {
    _multiconformer_01 = new Multiconformer_01[number_multiconformer_01_fingerprints_in_file];
    assert (NULL != _multiconformer_01);
  }

  if (number_multiconformer_fixed_size_counted_fingerprints_in_file > 0)
  {
    _multiconformer_fixed_size_counted = new Multiconformer_Fixed_Counted[number_multiconformer_fixed_size_counted_fingerprints_in_file];
    assert (NULL != _multiconformer_fixed_size_counted);
  }

  if (number_multiconformer_sparse_fingerprints_in_file > 0)
  {
    _multiconformer_sparse = new Multiconformer_Sparse[number_multiconformer_sparse_fingerprints_in_file];
    assert (NULL != _multiconformer_sparse);
  }

  return;
}

void
IW_General_Fingerprint::_default_values()
{
  _fingerprint = NULL;
  _sparse_fingerprint = NULL;

  _fixed_size_counted_fingerprint = NULL;

  _multiconformer_01 = NULL;

  _multiconformer_fixed_size_counted = NULL;

  _multiconformer_sparse = NULL;

  if (bit_fingerprints_in_file > 0 || _number_sparse_fingerprints || multiconformer_fingerprints_present)
    _allocate_fingerprint_array();

  return;
}

/*
  Notice that we don't initialise the arrays to NULL in the constructor
  because they are shared by all objects of this type. Hmmm, maybe
  static arrays aren't such a good idea after all
*/

IW_General_Fingerprint::IW_General_Fingerprint()
{
  _default_values();

  return;
}

IW_General_Fingerprint::IW_General_Fingerprint (const IW_General_Fingerprint & rhs)
{
  _default_values();

  for (int i = 0; i < bit_fingerprints_in_file; i++)
  {
    _fingerprint[i] = rhs._fingerprint[i];
  }

  for (int i = 0; i < _number_sparse_fingerprints; i++)
  {
    _sparse_fingerprint[i] = rhs._sparse_fingerprint[i];
  }

  for (int i = 0; i < number_multiconformer_01_fingerprints_in_file; i++)
  {
    _multiconformer_01[i] = rhs._multiconformer_01[i];
  }

  for (int i = 0; i < number_multiconformer_fixed_size_counted_fingerprints_in_file; i++)
  {
    _multiconformer_fixed_size_counted[i] = rhs._multiconformer_fixed_size_counted[i];
  }

  for (int i = 0; i < number_multiconformer_sparse_fingerprints_in_file; i++)
  {
    _multiconformer_sparse[i] = rhs._multiconformer_sparse[i];
  }

  _molecular_properties_integer = rhs._molecular_properties_integer;

  _molecular_properties_continuous = rhs._molecular_properties_continuous;

  return;
}

IW_General_Fingerprint::~IW_General_Fingerprint()
{
  if (_fingerprint)
  {
    delete [] _fingerprint;
    _fingerprint = NULL;
  }

  if (_sparse_fingerprint)
  {
    delete [] _sparse_fingerprint;
    _sparse_fingerprint = NULL;
  }

  if (NULL != _fixed_size_counted_fingerprint)
    delete [] _fixed_size_counted_fingerprint;


  if (NULL !=  _multiconformer_01)
    delete [] _multiconformer_01;

  if (NULL != _multiconformer_fixed_size_counted)
    delete [] _multiconformer_fixed_size_counted;

  if (NULL != _multiconformer_sparse)
    delete [] _multiconformer_sparse;

  return;
}

int
IW_General_Fingerprint::offset (off_t & result) const
{
  return _offset.value(result);
}

/*
  Fingerprints can be initialised by examining the first TDT in the dataset
  We should call this for the first TDT in each run.
*/

int
IW_General_Fingerprint::initialise (const IW_TDT & tdt)
{
  assert (0 == bit_fingerprints_in_file);

  if (! fingerprint_rx.active())
  {
    fingerprint_rx.set_pattern(fingerprint_rx_default);
  }

  int nfp = 0;    // regular fingerprints
  int sfp = 0;    // sparse form fingerprints
  int cfp = 0;    // counted form fingerprints
  int ddat = 0;   // descriptors as continuous properties
  int mc01 = 0;   // multiconformer 01
  int mcfsc = 0;  // multiconformer fixed size counted
  int mcsparse = 0;  // multiconformer sparse
 
  int i = 0;
  const_IWSubstring dataitem;

  while (tdt.next_dataitem(dataitem, i))
  {
    if (dataitem.starts_with("DDAT<"))   // the default descriptor tag
    {
      ddat++;
      continue;
    }

    if (! fingerprint_rx.matches(dataitem))
      continue;

    if (dataitem.starts_with("NC"))
      sfp++;
    else if (dataitem.starts_with("FC"))
      cfp++;
    else if (dataitem.starts_with("MC01"))
      mc01++;
    else if (dataitem.starts_with("MCFSC"))
      mcfsc++;
    else if (dataitem.starts_with("MCSP"))
      mcsparse++;
    else
      nfp++;
  }

  if (nfp <= 0 && sfp <= 0 && cfp <= 0 && mc01 <= 0 && mcfsc <= 0 && mcsparse <= 0)
  {
    cerr << "IW_General_Fingerprint::_initalize_auto: no fingerprints '" << fingerprint_rx.source() << "' in tdt\n";
    return 1;   // should this be return 0
  }

  if (! set_nfingerprints(nfp, sfp, cfp, mc01, mcfsc, mcsparse))
    return 0;

  if (report_fingerprint_status)
  {
    cerr << "Auto sized for\n";
    if (bit_fingerprints_in_file)
      cerr << ' ' << bit_fingerprints_in_file << " fingerprints\n";
    if (_number_sparse_fingerprints)
      cerr << ' ' << _number_sparse_fingerprints << " sparse fingerprints\n";
    if (number_fixed_size_counted_fingerprints_in_file)
      cerr << ' ' << number_fixed_size_counted_fingerprints_in_file << " fixed counted fingerprints\n";
    if (number_multiconformer_01_fingerprints_in_file)
      cerr << ' ' << number_multiconformer_01_fingerprints_in_file << " mc01\n";
    if (number_multiconformer_fixed_size_counted_fingerprints_in_file)
      cerr << ' ' << number_multiconformer_fixed_size_counted_fingerprints_in_file << " mcfsc\n";
    if (number_multiconformer_sparse_fingerprints_in_file)
      cerr << ' ' << number_multiconformer_sparse_fingerprints_in_file << " mcsp\n";
  }

// We need indices into each fingerprint type array

  nfp = 0;
  sfp = 0;
  cfp = 0;
  mc01 = 0;
  mcfsc = 0;
  mcsparse = 0;

  i = 0;
  const_IWSubstring tag;
  while (tdt.next_dataitem(tag, i))
  {
    if (! fingerprint_rx.matches(tag))     // doesn't look like a TDT
      continue;

    tag.truncate_at_first('<');     // just the tag part of the record

    if (tag.starts_with("NC"))
    {
      _sparse_fingerprint_tag[sfp] = tag;
      sfp++;
    }
    else if (tag.starts_with("FC"))
    {
      _fixed_size_counted_fingerprint_tag[cfp] = tag;
      cfp++;
    }
    else if (tag.starts_with("MC01"))
    {
      _multiconformer_01_fingerprint_tag[mc01] = tag;
      mc01++;
    }
    else if (tag.starts_with("MCFSC"))
    {
      _multiconformer_fixed_size_counted_fingerprint_tag[mcfsc] = tag;
      mcfsc++;
    }
    else if (tag.starts_with("MCSP"))
    {
      _multiconformer_sparse_fingerprint_tag[mcsparse] = tag;
      mcsparse++;
    }
    else
    {
      _fingerprint_tag[nfp] = tag;
      if (report_fingerprint_status)
        cerr << "initialise auto, tag " << nfp << " set to '" << _fingerprint_tag[nfp] << "'\n";

      if (tag.starts_with("HX"))
      {
        cerr << "Fingerprint '" << tag << "' automatically hex type\n";
        _fingerprint_type[nfp] = LOWERCASE_HEX_FORM;
      }
      nfp++;
    }
  }

  assert (nfp == bit_fingerprints_in_file);
  assert (sfp == _number_sparse_fingerprints);
  assert (cfp == number_fixed_size_counted_fingerprints_in_file);
  assert (mc01 == number_multiconformer_01_fingerprints_in_file);
  assert (mcfsc == number_multiconformer_fixed_size_counted_fingerprints_in_file);
  assert (mcsparse == number_multiconformer_sparse_fingerprints_in_file);

  if (molecular_property_integer_tag.length())    // look for molecular properties
  {
    const_IWSubstring notused;
    if (tdt.dataitem_value(molecular_property_integer_tag, notused))  // great
    {
      if (report_fingerprint_status)
        cerr << "initialise auto, properties present\n";
    }
    else if (! dash_P_specified)    // OK for properties to be absent
      molecular_property_integer_tag = "";
    else                          // they specified something but it is missing
    {
      cerr << "IW_General_Fingerprint::_initalize_auto: no '" << molecular_property_integer_tag << "' present\n";
      molecular_property_integer_tag = "";
    }
  }

  if (molecular_property_continuous_tag.length())    // look for molecular properties
  {
    const_IWSubstring notused;

    if (! tdt.dataitem_value(molecular_property_continuous_tag, notused))
    {
      cerr << "IW_General_Fingerprint::_initalize_auto: no '" << molecular_property_continuous_tag << "' present\n";
      molecular_property_continuous_tag = "";
    }
    else
      cerr << "initialise auto, continuous properties present\n";
  }
  else if (ddat > 0)
  {
    molecular_property_continuous_tag = "DDAT<";
    number_continuous_molecular_properties = ALL_CONTINUOUS_MOLECULAR_PROPERTIES;
  }

// We need to adjust the default weights. Beware, the properties weight may
// have been set via the -P option.

  if (report_fingerprint_status && molecular_property_integer_tag.length())
    cerr << "Molecular property tag '" << molecular_property_integer_tag << "', weight " << _property_weight_integer << ", nfp = " << nfp << endl;

  int items_specified = nfp + sfp + cfp + mc01 + mcfsc + mcsparse;
  if (molecular_property_integer_tag.length())
    items_specified++;
  if (molecular_property_continuous_tag.length())
    items_specified++;

  float weights_specified = static_cast<float>(0.0);
  int items_needing_default_weight = nfp + sfp + cfp + mc01 + mcfsc + mcsparse;

  if (_property_weight_integer > static_cast<float>(0.0))
    weights_specified += _property_weight_integer;
  else if (molecular_property_integer_tag.length())
    items_needing_default_weight++;

  if (_property_weight_continuous > static_cast<float>(0.0))
    weights_specified += _property_weight_continuous;
  else if (molecular_property_continuous_tag.length())
    items_needing_default_weight++;

  if (0 == items_needing_default_weight)
  {
    print_tags_and_weights(cerr);
    return 1;
  }

// The default weight that will be assigned to everything that doesn't have a weight already

  if (weights_specified >= static_cast<float>(1.0))
  {
    cerr << "IW_General_Fingerprint::initialise: huh, there are " << items_needing_default_weight << " items without weights, but total specified = " << weights_specified << endl;
    return 0;
  }

  float w = (1.0 - weights_specified) / static_cast<float>(items_needing_default_weight);

  if (w > 1.0 || w <= 0.0)    // how could this happen
  {
    cerr << "IW_General_Fingerprint::initialise: bad news default weight out of range " << w << endl;
    cerr << "Total weight specified " << weights_specified << ", " << items_needing_default_weight << " items needing default weight\n";
    return 0;
  }

  if (molecular_property_integer_tag.length())
    (void) initialise_properties_ratios();

  if (report_fingerprint_status)
  {
    cerr << "Setting ";
    if (number_fingerprints_to_use_in_computations)
      cerr << number_fingerprints_to_use_in_computations << " fingerprints ";
    if (_number_sparse_fingerprints_to_use_in_computations)
      cerr << _number_sparse_fingerprints_to_use_in_computations << " sparse fingerprint ";
    if (number_fixed_size_counted_fingerprints_to_use_in_computations)
      cerr << number_fixed_size_counted_fingerprints_to_use_in_computations << " fixed counted ";
    cerr << "weights to " << w << endl;
  }

// Note that we use the variables use_in_computations for filling the weight arrays

  for (int i = 0; i < number_fingerprints_to_use_in_computations; i++)
  {
    _fingerprint_weight[i] = w;
  }
  for (int i = 0; i < number_fixed_size_counted_fingerprints_to_use_in_computations; i++)
  {
    _fixed_size_counted_fingerprint_weight[i] = w;
  }
  for (int i = 0; i < _number_sparse_fingerprints_to_use_in_computations; i++)
  {
    _sparse_fingerprint_weight[i] = w;
  }
  for (int i = 0; i < number_multiconformer_01_fingerprints_in_file; i++)
  {
    _multiconformer_01_weight[i] = w;
  }
  for (int i = 0; i < number_multiconformer_fixed_size_counted_fingerprints_in_file; i++)
  {
    _multiconformer_fixed_size_counted_weight[i] = w;
  }
  for (int i = 0; i < number_multiconformer_sparse_fingerprints_in_file; i++)
  {
    _multiconformer_sparse_weight[i] = w;
  }

  if (molecular_property_integer_tag.length() && _property_weight_integer < 0.0)
    _property_weight_integer = w;
 
  if (molecular_property_continuous_tag.length() && _property_weight_continuous < 0.0)
    _property_weight_continuous = w;
 
  if (report_fingerprint_status)
    print_tags_and_weights(cerr);

  return 1; 
}

int
IW_General_Fingerprint::_read_molecular_properties_integer (IW_TDT & tdt,
                                                            const IWString & tag)
{
  const_IWSubstring mpr;

  if (! tdt.dataitem(tag, mpr))
  {
    cerr << "IW_General_Fingerprint::_read_molecular_properties_integer: no property tag '" << tag << "'\n";
    return 0;
  }

  return _molecular_properties_integer.construct_from_tdt_fp_record(mpr);
}

int
Molecular_Properties_Continuous::construct_from_tdt_record (const const_IWSubstring & buffer)
{
  const_IWSubstring p(buffer);
  p.remove_up_to_first('<');
  p.chop(1);

  return construct_from_descriptor_record(p);
}

int
IW_General_Fingerprint::_read_molecular_properties_continuous (IW_TDT & tdt,
                                                               const IWString & tag)
{
  const_IWSubstring mpr;
  if (! tdt.dataitem(tag, mpr))
  {
    cerr << "IW_General_Fingerprint::_read_molecular_properties_continuous: no property tag '" << tag << "'\n";
    return 0;
  }

  return _molecular_properties_continuous.construct_from_tdt_record(mpr);
}


/*
  Very simplistic. If we find a space in the first few characters of BUFFER,
  we decide that it is a record from a descriptor file
*/

static int
looks_like_descriptors (const const_IWSubstring & buffer)
{
  int i = buffer.index('<');

  assert (i > 0);      // cannot be a valid TDT record otherwise

// How many characters to scan looking for a space. Say 10, but really, 2 should
// be enough

  int istop = i + 10;
  if (istop > buffer.length())
    istop = buffer.length() - 1;

  for ( ; i < istop; i++)
  {
    if (' ' == buffer[i])
      return 1;
  }

// no spaces encountered

  return 0;
}

static void
remove_tdt_stuff (const_IWSubstring & buffer)
{
  assert (buffer.ends_with('>'));

  buffer.remove_up_to_first('<');

  buffer.chop();      // get rid of trailing angle bracket

  return;
}

int
IW_General_Fingerprint::_construct_from_tdt (IW_TDT & tdt,
                                             int & fatal)
{
  if (! tdt.dataitem_value(identifier_tag, _id))
  {
    cerr << "Ignoring TDT with no identifier '" << identifier_tag << "' in tdt\n";
    return 0;
  }

  if (0 == bit_fingerprints_in_file &&
      0 == _number_sparse_fingerprints &&
      0 == number_fixed_size_counted_fingerprints_in_file &&
      0 == number_multiconformer_fixed_size_counted_fingerprints_in_file &&
      0 == number_multiconformer_01_fingerprints_in_file &&
      0 == number_multiconformer_sparse_fingerprints_in_file)
  {
//  cerr << "Calling initialise\n";
    if (! initialise(tdt))
    {
      fatal = 1;
      return 0;
    }
  }

//#define ECHO_FINGERPRINT_DATA
#ifdef ECHO_FINGERPRINT_DATA
  cerr << "_nfingerprints = " << bit_fingerprints_in_file << endl;
  for (int i = 0; i < bit_fingerprints_in_file ; i++)
  {
    cerr << _bits_in_fingerprint[i] << " bits in fingerprint " << i << endl;
  }
#endif

  assert (bit_fingerprints_in_file > 0 || _number_sparse_fingerprints > 0 || number_fixed_size_counted_fingerprints_in_file || multiconformer_fingerprints_present);

  if (NULL == _fingerprint && NULL == _sparse_fingerprint && NULL == _fixed_size_counted_fingerprint)
    _allocate_fingerprint_array();

  for (int i = 0; i < bit_fingerprints_in_file; i++)
  {
    const_IWSubstring fp;
    if (! tdt.dataitem(_fingerprint_tag[i], fp))
    {
      cerr << "IW_General_Fingerprint::_construct_from_tdt: cannot find fixed width fingerprint '" << _fingerprint_tag[i] << "'\n";
      fatal = 1;
      return 0;
    }

    if (fp.ends_with('\n'))
      fp.chop();

    if (LOWERCASE_HEX_FORM == _fingerprint_type[i])
    {
      const_IWSubstring zhex(fp);
      remove_tdt_stuff(zhex);
      if (! _fingerprint[i].construct_from_hex(zhex))
      {
        cerr << "IW_General_Fingerprint::_construct_from_tdt: invalid hex form\n";
        cerr << fp << endl;
        fatal = 1;
        return 0;
      }
    }
    else if (UPPERCASE_HEX_FORM == _fingerprint_type[i])
    {
      IWString tmp(fp);
      tmp.to_lowercase();

      const_IWSubstring zhex(tmp);
      remove_tdt_stuff(zhex);
      if (! _fingerprint[i].construct_from_hex(zhex))
      {
        cerr << "IW_General_Fingerprint::_construct_from_tdt: invalid hex form\n";
        cerr << fp << endl;
        fatal = 1;
        return 0;
      }
    }
    else if (ASCII_01_FORM == _fingerprint_type[i])
    {
      const_IWSubstring zascii(fp);
      remove_tdt_stuff(zascii);

      if (! _fingerprint[i].construct_from_ascii_01_representation(zascii.rawchars(), zascii.length()))
      {
        cerr << "Cannot parse ascii form '" << zascii << "'\n";
        fatal = 1;
        return 0;
      }
//    std::cout << "Fingerprint has " << _fingerprint[i].nbits() << " bits, nset " << _fingerprint[i].nset() << endl;
//    _fingerprint[i].printon(cout, '1', '0');
      std::cout << '\n';
      if (0 != _fingerprint[i].nbits() % 32)
      {
        int b = _fingerprint[i].nbits();

        _fingerprint[i].allocate_space_for_bits( 32 * (b / 32 + 1));
      }
//    cout << "Fingerprint has " << _fingerprint[i].nbits() << " bits, nset " << _fingerprint[i].nset() << endl;
//    _fingerprint[i].printon(cout, '1', '0');
//    cout << endl;
    }
    else if (SPARSE_FORM == _fingerprint_type[i])
    {
      const_IWSubstring zascii(fp);
      remove_tdt_stuff(zascii);

      if (! _fingerprint[i].construct_from_sparse_representation(zascii))
      {
        cerr << "IW_General_Fingerprint::construct_from_tdt: cannot parse sparse bits representation\n";
        cerr << "'" << zascii << "'\n";
        return 0;
      }
    }
    else if (looks_like_descriptors(fp))
    {
      if (! _fingerprint[i].construct_from_descriptor_record(fp))
      {
        cerr << "IW_General_Fingerprint::_construct_from_tdt: invalid descriptor record\n";
        cerr << fp << endl;
        fatal = 1;
        return 0;
      }
    }
    else if (! _fingerprint[i].construct_from_tdt_record(fp))
    {
      cerr << "IW_General_Fingerprint::_construct_from_tdt: cannot parse fp " << _fingerprint_tag[i] << "\n";
      return 0;
    }

    if (_bits_in_fingerprint[i] <= 0)     // must be first time through
      _bits_in_fingerprint[i] = _fingerprint[i].nbits();
    else if (_bits_in_fingerprint[i] != _fingerprint[i].nbits())
    {
      cerr << "Bit count mismatch, stored = " << _bits_in_fingerprint[i] << ", input has " << _fingerprint[i].nbits() << endl;
      fatal = 1;
      return 0;
    }
  }

//cerr << "Reading " << _number_sparse_fingerprints << " sparse fingerprints\n";
  for (int i = 0; i < _number_sparse_fingerprints; i++)
  {
    const_IWSubstring fp;

    if (! tdt.dataitem(_sparse_fingerprint_tag[i], fp))
    {
      cerr << "IW_General_Fingerprint::_construct_from_tdt: cannot find sparse '" << _sparse_fingerprint_tag[i] << "'\n";
      fatal = 1;
      return 0;
    }
//  cerr << "Fetched data for '" << _sparse_fingerprint_tag[i] << "'\n";

    if (fp.ends_with('\n'))
      fp.chop();

    int rc = 0;
 
    if (DAYLIGHT_COMPRESSED_FORM == _sparse_fingerprint_type[i])
      rc = _sparse_fingerprint[i].construct_from_tdt_record(fp);
    else if (SPARSE_ASCII_FORM == _sparse_fingerprint_type[i])
    {
      remove_tdt_stuff(fp);
      rc = _sparse_fingerprint[i].construct_from_sparse_ascii_representation(fp);
    }
    else
    {
      cerr << "IW_General_Fingerprint:_construct_from_tdt:what to do with sparsefptype " << _sparse_fingerprint_type[i] << endl;
      return 0;
    }

    if (0 == rc)
    {
      cerr << "IW_General_Fingerprint::_construct_from_tdt: invalid sparse fingerprint\n";
      cerr << fp << endl;
      fatal = 1;
      return 0;
    }
  }

  for (int i = 0; i < number_fixed_size_counted_fingerprints_in_file; i++)
  {
    const_IWSubstring fp;
    if (! tdt.dataitem(_fixed_size_counted_fingerprint_tag[i], fp))
    {
      cerr << "IW_General_Fingerprint::_construct_from_tdt: cannot find fixed counted '" << _fixed_size_counted_fingerprint_tag[i] << "'\n";
      fatal = 1;
      return 0;
    }

    if (! _fixed_size_counted_fingerprint[i].construct_from_tdt_record(fp))
    {
      cerr << "IW_General_Fingerprint::_construct_from_tdt: invalid fixed counted fingerprint\n";
      cerr << fp << endl;
      fatal = 1;
      return 0;
    }
  }

  for (int i = 0; i < number_multiconformer_01_fingerprints_in_file; i++)
  {
    const_IWSubstring fp;
    if (! tdt.dataitem(_multiconformer_01_fingerprint_tag[i], fp))
    {
      cerr << "IW_General_Fingerprint::_construct_from_tdt: cannot find multiconformer 01 '" << _multiconformer_01_fingerprint_tag[i] << "'\n";
      fatal = 1;
      return 0;
    }

    if (fp.ends_with('\n'))
      fp.chop();

    if (! _multiconformer_01[i].construct_from_tdt_record(fp))
    {
      cerr << "IW_General_Fingerprint::_construct_from_tdt: invalid multiconformer 01 fingerprint\n";
      cerr << fp << endl;
      fatal = 1;
      return 0;
    }
  }

  for (int i = 0; i < number_multiconformer_fixed_size_counted_fingerprints_in_file; i++)
  {
    const_IWSubstring fp;
    if (! tdt.dataitem(_multiconformer_fixed_size_counted_fingerprint_tag[i], fp))
    {
      cerr << "IW_General_Fingerprint::_construct_from_tdt: cannot find multiconformer fsc '" << _multiconformer_fixed_size_counted_fingerprint_tag[i] << "'\n";
      fatal = 1;
      return 0;
    }

    if (fp.ends_with('\n'))
      fp.chop();

    if (! _multiconformer_fixed_size_counted[i].construct_from_tdt_record(fp))
    {
      cerr << "IW_General_Fingerprint::_construct_from_tdt: invalid multiconformer fsc fingerprint\n";
      cerr << fp << endl;
      fatal = 1;
      return 0;
    }
  }

  for (int i = 0; i < number_multiconformer_sparse_fingerprints_in_file; i++)
  {
    const_IWSubstring fp;
    if (! tdt.dataitem(_multiconformer_sparse_fingerprint_tag[i], fp))
    {
      cerr << "IW_General_Fingerprint::_construct_from_tdt: cannot find multiconformer sparse '" << _multiconformer_sparse_fingerprint_tag[i] << "'\n";
      fatal = 1;
      return 0;
    }

    if (fp.ends_with('\n'))
      fp.chop();

    if (! _multiconformer_sparse[i].construct_from_tdt_record(fp))
    {
      cerr << "IW_General_Fingerprint::_construct_from_tdt: invalid multiconformer sparse fingerprint\n";
      cerr << fp << endl;
      fatal = 1;
      return 0;
    }
  }

  if (molecular_property_integer_tag.length())
  {
    if (NULL == precomputed_ratio)
    {
      cerr << "Properties present '" << molecular_property_integer_tag << "', but not specified on command line, suppressed\n";
      molecular_property_integer_tag.resize(0);
    }
    else if (! _read_molecular_properties_integer(tdt, molecular_property_integer_tag))
    {
      fatal = 1;
      return 0;
    }
  }

  if (molecular_property_continuous_tag.length())
  {
    if (! _read_molecular_properties_continuous(tdt, molecular_property_continuous_tag))
    {
      fatal = 1;
      return 0;
    }
  }

  return 1;
}

int
IW_General_Fingerprint::construct_from_tdt (IW_TDT & tdt,
                                           int & fatal)
{

  fatal = 0;    // errors can be ignored unless we say so.

  if (_construct_from_tdt(tdt, fatal))   // success
    return 1;

  cerr << tdt;   // failure, echo the offending TDT

  return 0;
}

void
IW_General_Fingerprint::iwor (const IW_General_Fingerprint & rhs)
{
  if (_molecular_properties_continuous.active() || _molecular_properties_integer.active())
    cerr << "IW_General_Fingerprint::or: properties present, ignored\n";

  for (int i = 0; i < bit_fingerprints_in_file; i++)
  {
    _fingerprint[i].iwor( rhs[i]);
  }

  return;
}

/*
  A possibly strange function.
  We fail if any COMPONENT is out of range
*/

int
IW_General_Fingerprint::tanimoto(IW_General_Fingerprint & rhs,
                                 similarity_type_t must_be_larger_than,
                                 similarity_type_t & zresult)
{
  zresult = static_cast<similarity_type_t>(0.0);

// Check molecular properties first as they are cheapest

  if (_molecular_properties_integer.active())
  {
    const similarity_type_t mpr_integer = _molecular_properties_integer.similarity(rhs._molecular_properties_integer);

    if (mpr_integer < must_be_larger_than)   // sometimes including properties is a good idea, sometimes not
      return 0;

    zresult = zresult + _property_weight_integer * mpr_integer;
  }

  if (_molecular_properties_continuous.active())
  {
    const similarity_type_t mpr_continuous = _molecular_properties_continuous.similarity(rhs._molecular_properties_continuous);

    if (mpr_continuous < must_be_larger_than)
      return 0;

    zresult = zresult + _property_weight_continuous * mpr_continuous;
  }

  for (int i = 0; i < number_fingerprints_to_use_in_computations; i++)
  {
    similarity_type_t tmp = _fingerprint[i].tanimoto(rhs._fingerprint[i]);
    if (tmp < must_be_larger_than)
      return 0;

    zresult += _fingerprint_weight[i] * tmp;
  }

  for (int i = 0; i < number_fixed_size_counted_fingerprints_to_use_in_computations; i++)
  {
    similarity_type_t tmp = _fixed_size_counted_fingerprint[i].tanimoto(rhs._fixed_size_counted_fingerprint[i]);
    if (tmp < must_be_larger_than)
      return 0;

    zresult += _fixed_size_counted_fingerprint_weight[i] * tmp;
  }

  for (int i = 0; i < _number_sparse_fingerprints_to_use_in_computations; i++)
  {
    similarity_type_t tmp = _sparse_fingerprint[i].tanimoto(rhs._sparse_fingerprint[i]);
    if (tmp < must_be_larger_than)
      return 0;

    if (1.0 != _sparse_fingerprint_scaling_factor[i])
      tmp = 1.0 - _sparse_fingerprint_scaling_factor[i] * (1.0 - tmp);

    zresult += _sparse_fingerprint_weight[i] * tmp;
  }

  if (multiconformer_fingerprints_present)
  {
    for (int i = 0; i < number_multiconformer_01_fingerprints_in_file; i++)
    {
      similarity_type_t tmp = _multiconformer_01[i].tanimoto(rhs._multiconformer_01[i]);

      if (tmp < must_be_larger_than)
        return 0;

      zresult += _multiconformer_01_weight[i] * tmp;
    }

    for (int i = 0; i < number_multiconformer_fixed_size_counted_fingerprints_in_file; i++)
    {
      similarity_type_t tmp = _multiconformer_fixed_size_counted[i].tanimoto(rhs._multiconformer_fixed_size_counted[i]);
      if (tmp < must_be_larger_than)
        return 0;
      zresult += _multiconformer_fixed_size_counted_weight[i] * tmp;
    }

    for (int i = 0; i < number_multiconformer_sparse_fingerprints_in_file; i++)
    {
      similarity_type_t tmp = _multiconformer_sparse[i].tanimoto(rhs._multiconformer_sparse[i]);
      if (tmp < must_be_larger_than)
        return 0;
      zresult += _multiconformer_sparse_weight[i] * tmp;
    }
  }

  return 1;
}

/*
  We can establish an atom count window for comparing molecules.
*/

#define MAX_NATOMS 256

static int lower_atom_count_cutoff[MAX_NATOMS];
static int upper_atom_count_cutoff[MAX_NATOMS];

static int atom_count_window_present = 0;
static int ring_count_window_present = 0;

static int lower_aromatic_atom_count_cutoff[MAX_NATOMS];
static int upper_aromatic_atom_count_cutoff[MAX_NATOMS];

static int aromatic_atom_count_window_present = 0;

static int windows_present = 0;

int
property_based_windows_present()
{
  return windows_present;
}

/*static void
no_atom_count_cutoffs()
{
  for (int i = 0; i < MAX_NATOMS; i++)
  {
    lower_atom_count_cutoff[i] = 1;
    upper_atom_count_cutoff[i] = MAX_NATOMS;
  }

  return;
}*/

/*static void
no_ring_count_cutoffs()
{
  for (int i = 0; i < MAX_NRINGS; i++)
  {
    lower_ring_count_cutoff[i] = 0;
    upper_ring_count_cutoff[i] = MAX_NRINGS;
  }

  return;
}*/

int
set_atom_count_window (int w, int verbose)
{
  assert (w > 1 && w <= 100);     // must be a percentage

  for (int i = 1; i < MAX_NATOMS; i++)
  {
    int delta = int (static_cast<float> (i * w) / 100.0);

    if (i - delta < 1)
      lower_atom_count_cutoff[i] = 1;
    else
      lower_atom_count_cutoff[i] = i - delta;

    upper_atom_count_cutoff[i] = i + delta;

    if (verbose)
      cerr << "Molecules with " << i << " atoms compared with " << lower_atom_count_cutoff[i] << " to " << upper_atom_count_cutoff[i] << " atoms\n";
  }

  atom_count_window_present = 1;

  return 1;
}

static int
set_atom_count_delta_window (int min_extra_atoms, int max_extra_atoms)
{
  for (int i = 1; i < MAX_NATOMS; i++)
  {
    if (max_extra_atoms >= 0)
      upper_atom_count_cutoff[i] = i + max_extra_atoms;
    else
      upper_atom_count_cutoff[i] = MAX_NATOMS;

    if (min_extra_atoms >= 0)
      lower_atom_count_cutoff[i] = i + min_extra_atoms;
    else
      lower_atom_count_cutoff[i] = 1;
  }

  atom_count_window_present = 1;

  return 1;
}

static int
set_aromatic_atom_count_window (int w)
{
  for (int i = 0; i < MAX_NATOMS; i++)
  {
    lower_aromatic_atom_count_cutoff[i] = i - w;
    if (lower_aromatic_atom_count_cutoff[i] < 0)
      lower_aromatic_atom_count_cutoff[i] = 0;

    upper_aromatic_atom_count_cutoff[i] = i + w;
  }

  aromatic_atom_count_window_present = 1;

  return 1;
}

/*
  Same thing for rings
*/

#define MAX_NRINGS 30

static int lower_ring_count_cutoff[MAX_NRINGS];
static int upper_ring_count_cutoff[MAX_NRINGS];

int
set_ring_count_window (int w, int verbose)
{
  assert (w > 1 && w <= 100);     // must be a percentage

  for (int i = 1; i < MAX_NRINGS; i++)
  {
    int delta = int (static_cast<float> (i * w) / 100.0);

    if (i - delta < 0)
      lower_ring_count_cutoff[i] = 0;
    else
      lower_ring_count_cutoff[i] = i - delta;

    upper_ring_count_cutoff[i] = i + delta;

    if (verbose)
      cerr << "Molecules with " << i << " rings compared with " << lower_ring_count_cutoff[i] << " to " << upper_ring_count_cutoff[i] << " rings\n";
  }

  ring_count_window_present = 1;

  return 1;
}

//#define DEBUG_CAN_BE_COMPARED

/*
  Programmes can call this function to see if two fingerprints can
  be compared. Typically fp2 is the member of the haystack
*/

int
can_be_compared (const IW_General_Fingerprint & fp1,
                 const IW_General_Fingerprint & fp2)
{
  if (0 == windows_present)
    return 1;

  if (atom_count_window_present)
  {
    int na1 = fp1.natoms();
    int na2 = fp2.natoms();

    if (na1 >= MAX_NATOMS || na2 >= MAX_NATOMS)
      return 1;

#ifdef DEBUG_CAN_BE_COMPARED
    cerr << "Can we compare '" << fp1.id() << " " << na1 << " (" << lower_atom_count_cutoff[na1] << ',' << upper_atom_count_cutoff[na1] << ") and '" << fp2.id() << " " << na2 << endl; //" (" << lower_atom_count_cutoff[na2] << ',' << upper_atom_count_cutoff[na2] << ")\n";
#endif

    if (na2 < lower_atom_count_cutoff[na1])
      return 0;

    if (na2 > upper_atom_count_cutoff[na1])
      return 0;

//  Jul 2016. In general however, we do need to compare both ways. Have not bothered changing the others, but should do it...

    if (na1 < lower_atom_count_cutoff[na2])
      return 0;

    if (na1 > upper_atom_count_cutoff[na2])
      return 0;
  }

  if (ring_count_window_present)
  {
    int nr1 = fp1.nrings();
    int nr2 = fp2.nrings();

//  should check out of range, but seems unlikely

#ifdef DEBUG_CAN_BE_COMPARED
    cerr << "Rings?? " << nr1 << " (" << lower_ring_count_cutoff[nr1] << ',' << upper_ring_count_cutoff[nr1] << ") and " << nr2 << " (" << lower_ring_count_cutoff[nr2] << ',' << upper_ring_count_cutoff[nr2] << ")\n";
#endif
    if (nr2 < lower_ring_count_cutoff[nr1])
      return 0;

    if (nr2 > upper_ring_count_cutoff[nr1])
      return 0;
  }

  if (aromatic_atom_count_window_present)
  {
    int na1 = fp1.aromatic_atoms();
    int na2 = fp2.aromatic_atoms();

    if (na2 < lower_aromatic_atom_count_cutoff[na1])
      return 0;

    if (na2 > upper_aromatic_atom_count_cutoff[na1])
      return 0;
  }

#ifdef DEBUG_CAN_BE_COMPARED
  cerr << "Yes, can be compared\n";
#endif

  return 1;      // these fingerprints can be compared
}

int
display_standard_gfp_options (std::ostream & os)
{
  os << " -F <TAG>         fingerprint tag 'TAG'\n";
  os << " -F <TAG,w=0.2>   fingerprint 'TAG' weight 0.2\n";
  os << " -F <TAG,fold=n>  fingerprint 'TAG' fold 'n' times DNU\n";
  os << " -F <TAG,hex>     fingerprint 'TAG', lowercase hex form DNU\n";
  os << " -F <TAG,HEX>     fingerprint 'TAG', uppercase hex form DNU\n";
  os << " -F <TAG,ascii>   fingerprint 'TAG', ascii 0,1 form DNU\n";
  os << " -F <TAG,sparse>  fingerprint 'TAG', sparse bit representation <1,6-10,22;32> DNU\n";
  os << " -F <TAG,nc>      non-colliding fingerprint 'TAG'\n";
  os << " -F <TAG,cosine>  fingerprint 'TAG', use cosine measure as distance DNU\n";
  os << " -F <TAG,fc>      fixed size counted fingerprint 'TAG' DNU\n";
  os << " -F <TAG,mc01>    multiconformer 01 fingerprint in 'TAG' DNU\n";
  os << " -F <TAG,mcfsc>   multiconformer fixed size counted fingerprint in 'TAG' DNU\n";
  os << " -F <TAG,mcsp>    multiconformer sparse counted fingerprint in 'TAG' DNU\n";
  os << " -P <TAG>         integer molecular properties in 'TAG' (default '" << molecular_property_integer_tag << "')\n";
  os << " -P none          no processing of molecular properties\n";
  os << " -P <TAG,w=0.1>   integer molecular properties in 'TAG' weight 0.1\n";
  os << " -P <TAG,desc=n,w=0.1> continuous (descriptor) properties, <n> of them DNU\n";
  os << "                  'desc=all' to use all DNU\n";
  os << " -P cartesian     use cartesian distances - make sure you use scaling DNU\n";
  os << " -P scale=x.xx    scale factor applied to each continuous property as read in DNU\n";
  os << " -W A:nn          atom count window of nn% (nn between 1 and 100) DNU\n";
  os << " -W A:maxextra=nn at most  nn extra atoms allowed  for comparison DNU\n";
  os << " -W A:minextra=nn at least nn extra atoms required for comparison DNU\n";
  os << " -W A:-nn         at most nn fewer atoms allowed for comparison DNU\n";
  os << " -W R:nn          ring count window of nn% (nn between 1 and 100)\n";
  os << " -Q dddm=xxx      default dense distance metric (tan fvb) DNU\n";
  os << " -Q dsdm=xxx      default sparse distance metric (tan fvb manh soergel soergelv ctan dice binary binary) DNU\n";
  os << " -Q dcdm=xxx      default counted distance metric (tan fvb) DNU\n";

  return os.good();
}

static int
parse_guard_fingerprint_specification (const const_IWSubstring & g)
{
  IWString tmp;

  if (! g.split(guard_fingerprint_tag, ',', tmp))
  {
    cerr << "Guard fingerprint must contain a tag and a distance\n";
    return 0;
  }

  if (! tmp.numeric_value(guard_fingerprint_distance) || guard_fingerprint_distance < 0.0)
  {
    cerr << "Invalid guard fingerprint distance '" << tmp << "'\n";
    return 0;
  }

  return 1;
}

/*
  See if there is some kind of whole number at the front of buffer...
*/

/*static int
fetch_numeric (const const_IWSubstring & buffer,
               int & zresult)
{
  zresult = 0;

  for (int i = 0; i < buffer.length(); i++)
  {
    char c = buffer[i];
    if (c >= '0' && c <= '9')
      zresult = 10 * zresult + c - '0';
    else
      return i + 1;
  }

  return buffer.length();
}*/

static int
display_dash_p_options (std::ostream & output)
{
  output << "There are two kinds of properties in gfp programmes, integer and continuous\n";
  output << "The integer properties are generated by default by gfp_make.sh, and are in the tag MPR\n";
  output << "To suppress properties from your computation, enter '-P none'\n";
  output << "To give a weight to integer molecular properties try something like '-P MPR,wt=0.1'\n";
  output << "To use the Dice similarity metric with integer properties use '-P dice'\n";
  output << "Continuous properties are any floating point numbers. Various means are available for adding\n";
  output << "descriptors and other numbers to gfp files, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol) for details\n";
  output << "Specification of continuous properties must be via one -P option, and must contain 'desc='\n";
  output << "To use properties in the tag FOO, try '-P FOO,desc=all' to take all the numeric descriptors\n";
  output << "Add 'w=0.3' to change the weight of continuous properties\n";
  output << "Add 'dice' or 'cartesian' to change the distance metric used by continuous properties\n";

  return output.good();
}

/*
  Avoid passing lots and lots of arguments
*/

class Parse_Tag_Weight_Args
{
  private:
    int _ftype;
    float _weight;
    int _nfold;
    int _ascii;
    float _ncscale;
    int _distance_metric;

  public:
    Parse_Tag_Weight_Args();

    int ftype() const { return _ftype;}
    void set_ftype(int s) { _ftype = s;}

    float weight() const { return _weight;}
    void set_weight(float s) { _weight = s;}

    int nfold() const { return _nfold;}
    void set_nfold(int s) { _nfold = s;}

    int ascii() const { return _ascii;}
    void set_ascii(int s) { _ascii = s;}

    float ncscale() const { return _ncscale;}
    void set_ncscale(float s) { _ncscale = s;}

    int distance_metric() const { return _distance_metric;}
    void set_distance_metric(int s) { _distance_metric = s;}
};

Parse_Tag_Weight_Args::Parse_Tag_Weight_Args()
{
  _ftype = DAYLIGHT_COMPRESSED_FORM;
  _weight = static_cast<float>(-1.0);
  _nfold = 0;
  _ascii = 0;
  _ncscale = static_cast<float>(1.0);
  _distance_metric = 0;    // same as the two default tanimoto methods

  return;
}

static int
parse_tag_weight (const const_IWSubstring & f,
                  const_IWSubstring & tag,
                  Parse_Tag_Weight_Args & ptwa)
{
  tag = f;

  int i = f.index(',');

  if (i < 0)     // no qualifiers
    return 1;

  tag.iwtruncate(i);
  i++;

  const_IWSubstring token;
  while (f.nextword(token, i, ','))
  {
    if (token.starts_with("w=") || token.starts_with("W="))
    {
      float w;
      token.remove_leading_chars(2);
      if (! token.numeric_value(w) || w < 0.0)
      {
        cerr << "Invalid weight '" << token << "'\n";
        return 0;
      }
      ptwa.set_weight(w);
    }
    else if (token.starts_with("fold=") || token.starts_with("FOLD="))
    {
      int nfold;
      token.remove_leading_chars(5);
      if (! token.numeric_value(nfold) || nfold < 2)
      {
        cerr << "Invalid fold value '" << token << "'\n";
        return 0;
      }
      ptwa.set_nfold(nfold);
    }
    else if ("hex" == token)
    {
      ptwa.set_ftype(LOWERCASE_HEX_FORM);
    }
    else if ("HEX" == token)
    {
      ptwa.set_ftype(UPPERCASE_HEX_FORM);
    }
    else if ("ascii" == token)
    {
      ptwa.set_ascii(1);
    }
    else if ("sparse" == token)
    {
      ptwa.set_ftype(SPARSE_FORM);
    }
    else if ("counted" == token || "fc" == token)
    {
      ptwa.set_ftype(FIXED_COUNTED_FORM);
    }
    else if ("cosine" == token)
    {
      ptwa.set_distance_metric(SPARSE_DISTANCE_METRIC_COSINE);
    }
    else if ("fvb" == token)
    {
      ptwa.set_distance_metric(FVB_MODIFIED_TANIMOTO);
    }
    else if ("rr" == token)
    {
      ptwa.set_distance_metric(RUSSEL_RAO);
    }
    else if ("forbes" == token)
    {
      ptwa.set_distance_metric(FORBES_SIMILARITY);
    }
    else if ("sm" == token)
    {
      ptwa.set_distance_metric(SIMPLE_MATCHING);
    }
    else if ("manh" == token)
    {
      ptwa.set_distance_metric(MANHATTAN_DISTANCE_METRIC);
    }
    else if ("nc" == token || "NC" == token)
    {
      ptwa.set_ftype(NON_COLLIDING_FORM);
    }
    else if ("mc01" == token)
    {
      ptwa.set_ftype(MULTICONFORMER_01_FORM);
    }
    else if ("mcfsc" == token)
    {
      ptwa.set_ftype(MULTICONFORMER_FSC_FORM);
    }
    else if ("mcsp" == token)
    {
      ptwa.set_ftype(MULTICONFORMER_SPARSE_FORM);
    }
    else if (token.starts_with("scale="))
    {
      token.remove_leading_chars(6);
      float ncscale;
      if (! token.numeric_value(ncscale) || ncscale < 0.0)
      {
        cerr << "Invalid sparse fingerprint scaling value '" << token << "'\n";
        return 0;
      }
      ptwa.set_ncscale(ncscale);
    }
    else if ("soergel" == token)
    {
      ptwa.set_distance_metric(SOERGEL_DISTANCE_METRIC);
    }
    else if ("soergelv" == token)
    {
      ptwa.set_distance_metric(SOERGEL_VARIANT_DISTANCE_METRIC);
    }
    else if ("ctan" == token)
      ptwa.set_distance_metric(CTAN_DISTANCE_METRIC);
    else if ("dice" == token)
      ptwa.set_distance_metric(DICE_DISTANCE_METRIC);
    else if ("overlap" == token)
      ptwa.set_distance_metric(DENSE_DISTNANCE_OVERLAP);
    else
    {
      cerr << "Unrecognised token qualifier '" << token << "'\n";
      return 0;
    }
  }

  if (SPARSE_DISTANCE_METRIC_COSINE == ptwa.distance_metric() && NON_COLLIDING_FORM != ptwa.ftype())
  {
    cerr << "Cosine measure only valid with non-colliding counted fingerprints\n";
    return 0;
  }

  return 1;
}


/*
  A tag will look like '^[:alnum:]+$' 
  or
  '^[:alnum:]+w=[:number:]$'
*/

static int
parse_tag_weight (const const_IWSubstring & f,
                  const_IWSubstring & tag,
                  float & weight,
                  int & nfold,
                  int & ztype,
                  double & ncscale,
                  int & distance_metric)
{
  tag = f;
  ztype = DAYLIGHT_COMPRESSED_FORM;
  ncscale = 1.0;

  int i = f.index(',');

  if (i < 0)     // no qualifiers
    return 1;

  tag.iwtruncate(i);
  i++;

  const_IWSubstring token;
  while (f.nextword(token, i, ','))
  {
    if (token.starts_with("w=") || token.starts_with("W="))
    {
      token.remove_leading_chars(2);
      if (! token.numeric_value(weight) || weight < 0.0)
      {
        cerr << "Invalid weight '" << token << "'\n";
        return 0;
      }
    }
    else if (token.starts_with("fold=") || token.starts_with("FOLD="))
    {
      token.remove_leading_chars(5);
      if (! token.numeric_value(nfold) || nfold < 2)
      {
        cerr << "Invalid fold value '" << token << "'\n";
        return 0;
      }
    }
    else if ("hex" == token)
    {
      ztype = LOWERCASE_HEX_FORM;
    }
    else if ("HEX" == token)
    {
      ztype = UPPERCASE_HEX_FORM;
    }
    else if ("ascii" == token)
    {
      ztype = ASCII_01_FORM;
    }
    else if ("sparse" == token)
    {
      ztype = SPARSE_FORM;
    }
    else if ("counted" == token || "fc" == token)
    {
      ztype = FIXED_COUNTED_FORM;
    }
    else if ("cosine" == token)
    {
      distance_metric = SPARSE_DISTANCE_METRIC_COSINE;
    }
    else if ("fvb" == token)
    {
      distance_metric = FVB_MODIFIED_TANIMOTO;
    }
    else if ("rr" == token)
    {
      distance_metric = RUSSEL_RAO;
    }
    else if ("forbes" == token)
    {
      distance_metric = FORBES_SIMILARITY;
    }
    else if ("sm" == token)
    {
      distance_metric = SIMPLE_MATCHING;
    }
    else if ("nc" == token || "NC" == token)
    {
      ztype = NON_COLLIDING_FORM;
    }
    else if ("mc01" == token)
    {
      ztype = MULTICONFORMER_01_FORM;
    }
    else if ("mcfsc" == token)
    {
      ztype = MULTICONFORMER_FSC_FORM;
    }
    else if ("mcsp" == token)
    {
      ztype = MULTICONFORMER_SPARSE_FORM;
    }
    else if ("ctan" == token)
      distance_metric = CTAN_DISTANCE_METRIC;
    else if ("dice" == token)
      distance_metric = DICE_DISTANCE_METRIC;
    else if ("manh" == token)
      distance_metric = MANHATTAN_DISTANCE_METRIC;
    else if (token.starts_with("scale="))
    {
      token.remove_leading_chars(6);
      if (! token.numeric_value(ncscale) || ncscale < 0.0)
      {
        cerr << "Invalid sparse fingerprint scaling value '" << token << "'\n";
        return 0;
      }
    }
    else
    {
      cerr << "Unrecognised token qualifier '" << token << "'\n";
      return 0;
    }
  }

  if (SPARSE_DISTANCE_METRIC_COSINE == distance_metric && NON_COLLIDING_FORM != ztype)
  {
    cerr << "Cosine measure only valid with non-colliding counted fingerprints\n";
    return 0;
  }

  return 1;
}

/*
  Must be of the form

  TAG,desc=X,....
*/

static int
parse_desc_specification (const const_IWSubstring & p,
                          int verbose)
{
  assert (p.contains("desc="));

  int i = 0;
  const_IWSubstring token;

  while (p.nextword(token, i, ','))
  {
    cerr << "Examining token '" << token << "'\n";
    if (token.starts_with("desc="))
    {
      token.remove_leading_chars(5);
      if (0 == token.length())
        number_continuous_molecular_properties = ALL_CONTINUOUS_MOLECULAR_PROPERTIES;
      else if("all" == token || "ALL" == token)
        number_continuous_molecular_properties = ALL_CONTINUOUS_MOLECULAR_PROPERTIES;
      else
      {
        if (! token.numeric_value(number_continuous_molecular_properties) || number_continuous_molecular_properties <= 0)
        {
          cerr << "Invalid number of descriptors '" << p << "'\n";
          return 0;
        }
      }
    }
    else if (token.starts_with("w="))
    {
      token.remove_leading_chars(2);
      if (! token.numeric_value(_property_weight_continuous) || _property_weight_continuous < 0.0)
      {
        cerr << "Invalid continuous property weight '" << token << "' in '" << p << "'\n";
        return 0;
      }
    }
    else if ("cartesian" == token)
    {
      _continuous_property_distance_metric = CONTINUOUS_PROPERTY_DISTANCE_METRIC_CART; 
    }
    else if (token.starts_with("exp1="))
    {
      _continuous_property_distance_metric = CONTINUOUS_PROPERTY_DISTANCE_METRIC_EXP1; 
      token.remove_leading_chars(5);
      if (! token.numeric_value(continuous_property_distance_metric_exp1_value))
      {
        cerr << "Invalid exp1 qualifier '" << token << "'\n";
        return 0;
      }

      if (continuous_property_distance_metric_exp1_value > 0.0)
        continuous_property_distance_metric_exp1_value = - continuous_property_distance_metric_exp1_value;
    }
    else if (token.starts_with("exp2="))
    {
      _continuous_property_distance_metric = CONTINUOUS_PROPERTY_DISTANCE_METRIC_EXP2; 
      token.remove_leading_chars(5);
      if (! token.numeric_value(continuous_property_distance_metric_exp2_value))
      {
        cerr << "Invalid exp2 qualifier '" << token << "'\n";
        return 0;
      }

      if (continuous_property_distance_metric_exp2_value > 0.0)
        continuous_property_distance_metric_exp2_value = - continuous_property_distance_metric_exp2_value;
    }
    else if (token.starts_with("scale="))
    {
      token.remove_leading_chars(6);

      if (! token.numeric_value(individual_property_scaling_factor))
      {
        cerr << "Invalid property scaling factor '" << p << "'\n";
        return 0;
      }

      if (verbose)
        cerr << "Each individual continuous property multipled by '" << individual_property_scaling_factor << " on reading\n";
    }
    else if ("dice" == token)
    {
      _continuous_property_distance_metric = CONTINUOUS_PROPERTY_DISTANCE_METRIC_DICE;
    }
    else if (0 == molecular_property_continuous_tag.length())
    {
      molecular_property_continuous_tag = token;
    }
    else
    {
      cerr << "Unrecognised -P qualifier '" << token << "' in '" << p << "'\n";
      return 0;
    }
  }

  if (verbose)
    cerr << "Continuous properties in '" << molecular_property_continuous_tag << "', ";
  if (number_continuous_molecular_properties > 0)
    cerr << number_continuous_molecular_properties << " descriptors.";
  else
    cerr << "all descriptors.";
  if (_property_weight_continuous >= 0.0)
    cerr << " weight " << _property_weight_continuous;
  if (CONTINUOUS_PROPERTY_DISTANCE_METRIC_RATIO == _continuous_property_distance_metric)
    ;
  else if (CONTINUOUS_PROPERTY_DISTANCE_METRIC_CART == _continuous_property_distance_metric)
    cerr << ", cartesian distances\n";
  else if (CONTINUOUS_PROPERTY_DISTANCE_METRIC_EXP1 == _continuous_property_distance_metric)
    cerr << ", exp1 distances\n";
  else if (CONTINUOUS_PROPERTY_DISTANCE_METRIC_EXP2 == _continuous_property_distance_metric)
    cerr << ", exp2 distances\n";
  cerr << endl;

  return 1;
}

int
initialise_fingerprints (Command_Line & cl, int verbose)
{
  bool options_specified_by_user = 0;
  if (cl.option_present('Q'))
    options_specified_by_user = 1;
  else if (cl.option_present('F'))
    options_specified_by_user = 1;
  else if (cl.option_present('P'))
    options_specified_by_user = 1;

  int i = 0;
  const_IWSubstring q;
  while (cl.value('Q', q, i++))
  {
    if (q.starts_with("dddm="))
    {
      q.remove_leading_chars(5);
      if (q.starts_with("tan"))
        default_dense_fingerprint_distance_metric = DENSE_DISTANCE_METRIC_TANIMOTO;
      else if ("fvb" == q)
        default_dense_fingerprint_distance_metric = FVB_MODIFIED_TANIMOTO;
      else if ("rr" == q)
        default_dense_fingerprint_distance_metric = RUSSEL_RAO;
      else if ("forbes" == q)
        default_dense_fingerprint_distance_metric = FORBES_SIMILARITY;
      else if ("sm" == q)
        default_dense_fingerprint_distance_metric = SIMPLE_MATCHING;
      else
      {
        cerr << "Unrecognised default dense fingerprint metric 'dfdm=" << q << "'\n";
        return 0;
      }
    }
    else if (q.starts_with("dsdm="))
    {
      q.remove_leading_chars(5);
      if (q.starts_with("tan"))
      {
        default_sparse_fingerprint_distance_metric = SPARSE_DISTANCE_METRIC_TANIMOTO;
      }
      else if ("fvb" == q)
      {
        default_sparse_fingerprint_distance_metric = FVB_MODIFIED_TANIMOTO;
      }
      else if ("ctan" == q)
        default_sparse_fingerprint_distance_metric = CTAN_DISTANCE_METRIC;
      else if (q.starts_with("ctan="))
      {
        q.remove_leading_chars(5);
        double x;
        if (! q.numeric_value(x)  || 0.0 == x)
        {
          cerr << "The continuous tanimoto exponent (ctan=" << q << ") must be a valid number\n";
          return 0;
        }
        set_continuous_tanimoto_exponent(x);
      }
      else if ("dice" == q)
        default_sparse_fingerprint_distance_metric = DICE_DISTANCE_METRIC;
      else if ("soergel" == q)
        default_sparse_fingerprint_distance_metric = SOERGEL_DISTANCE_METRIC;
      else if ("soergelv" == q)
        default_sparse_fingerprint_distance_metric = SOERGEL_VARIANT_DISTANCE_METRIC;
      else if ("binary" == q)
        default_sparse_fingerprint_distance_metric = SPARSE_DISTANCE_METRIC_BINARY;
      else
      {
        cerr << "Unrecognised default sparse fingerprint metric 'dfsm=" << q << "'\n";
        return 0;
      }
    }
    else if (q.starts_with("dcdm="))
    {
      q.remove_leading_chars(5);
      if (q.starts_with("tan"))
        _fixed_size_counted_fingerprint_distance_metric = FIXED_SIZE_COUNTED_FINGERPRINT_DISTANCE_METRIC_TANIMOTO;
      else if (q.starts_with("fvb"))
        default_fixed_size_counted_fingerprint_distance_metric = FVB_MODIFIED_TANIMOTO;
      else
      {
        cerr << "Unrecognised default counted distance metric qualifier '" << q << "'\n";
        return 0;
      }
    }
    else if (q.starts_with("FVB="))
    {
      q.remove_leading_chars(4);
      if (! set_fvb_ratios(q))
      {
        cerr << "INvalid FVB ratio specification '" << q << "'\n";
        return 0;
      }
    }
    else if ("mcave" == q)
    {
      multiconformer_use_average_similarity = 1;
    }
    else if ("help" == q)
    {
      cerr << " -Q dddm=...            default dense  distance metric\n";
      cerr << " -Q dsdm=...            default sparse distance metric\n";
      cerr << " -Q FBV=a,b             specify the two FVB parameters\n";
      cerr << " mcave                  use averages between multiconformer fp's\n";

      exit(2);
    }
    else
    {
      cerr << "Unrecognised 'Q qualifier '" << q << "'\n";
      return 0;
    }
  }

  int number_dash_F = cl.option_count('F');

  int nfp = 0;
  int sfp = 0;
  int cfp = 0;
  int mc01 = 0;
  int mcfsc = 0;
  int mcsparse = 0;

  if (0 == number_dash_F)
  {
    fingerprint_rx.set_pattern(fingerprint_rx_default);
    if (verbose > 1)
      cerr << "All fingerprints matching '" << fingerprint_rx.source() << "' will be processed\n";
  } 
  else
  {
    for (int i = 0; i < number_dash_F; i++)
    {
      const_IWSubstring f;
      cl.value('F', f, i);

      if (f.contains(",nc") || f.contains(",NC"))
        sfp++;
      else if (f.contains(",fc") || f.contains(",FC"))
        cfp++;
      else if (f.contains(",mc01"))
        mc01++;
      else if (f.contains(",mcfsc"))
        mcfsc++;
      else if (f.contains(",mcsp"))
        mcsparse++;
      else
        nfp++;
    }

    if (! set_nfingerprints(nfp, sfp, cfp, mc01, mcfsc, mcsparse))
    {
      cerr << "initialise_fingerprints: yipes, cannot initialise for " << nfp << " fingerprints\n";
      return 0;
    }
  }

  if (nfp > 0 || sfp > 0 || cfp > 0 || mc01 > 0 || mcfsc > 0 || mcsparse > 0)
  {
    nfp = 0;
    sfp = 0;
    int cfp = 0;
    mc01 = 0;
    mcfsc = 0;
    mcsparse = 0;
    for (int i = 0; i < number_dash_F; i++)
    {
      const_IWSubstring f;
      cl.value('F', f, i);

      if ("help" == f)
      {
        display_standard_gfp_options(cerr);
        exit (2);
      }

      const_IWSubstring tag;

      Parse_Tag_Weight_Args ptwa;
      if (! parse_tag_weight(f, tag, ptwa))
      {
        cerr << "Invalid -F qualifier '" << f << "'\n";
        return 0;
      }

      int ftype = ptwa.ftype();
      float weight = ptwa.weight();
      double ncscale = ptwa.ncscale();
      int distance_metric = ptwa.distance_metric();

//    cerr << "ftype " << ftype << " ascii? " << ptwa.ascii() << "\n";

      if (NON_COLLIDING_FORM == ftype && 0 == ptwa.ascii())
      {
        _sparse_fingerprint_weight[sfp] = weight;
        _sparse_fingerprint_tag[sfp] = tag;
        _sparse_fingerprint_scaling_factor[sfp] = ncscale;

        if (0 != distance_metric)
          _sparse_fingerprint_distance_metric[sfp] = distance_metric;
//      cerr << "Sparse fingerprint " << sfp << " distance metric " << distance_metric << endl;
        sfp++;
      }
      else if (NON_COLLIDING_FORM == ftype && ptwa.ascii())
      {
        _sparse_fingerprint_weight[sfp] = weight;
        _sparse_fingerprint_tag[sfp] = tag;
        _sparse_fingerprint_scaling_factor[sfp] = ncscale;
        _sparse_fingerprint_type[sfp] = SPARSE_ASCII_FORM;

        if (0 != distance_metric)
          _sparse_fingerprint_distance_metric[sfp] = distance_metric;
//      cerr << "Sparse fingerprint " << sfp << " distance metric " << distance_metric << endl;
        sfp++;
        sparse_fingerprint_counts_limited = 0;
      }
      else if (FIXED_COUNTED_FORM == ftype)
      {
        _fixed_size_counted_fingerprint_weight[cfp] = weight;
        _fixed_size_counted_fingerprint_tag[cfp] = tag;
        if (0 != distance_metric)
          _fixed_size_counted_fingerprint_distance_metric[cfp] = distance_metric;
         cfp++;
      }
      else if (MULTICONFORMER_01_FORM == ftype)
      {
        _multiconformer_01_weight[mc01] = weight;
        _multiconformer_01_fingerprint_tag[mc01] = tag;
        if (0 != distance_metric)
          _multiconformer_01_distance_metric[mc01] = distance_metric;
        cerr << "Multiconformer 01 " << mc01 << " distance metric " << distance_metric << endl;
        mc01++;
      }
      else if (MULTICONFORMER_FSC_FORM == ftype)
      {
        _multiconformer_fixed_size_counted_weight[mcfsc] = weight;
        _multiconformer_fixed_size_counted_fingerprint_tag[mcfsc] = tag;
        if (0 != distance_metric)
          _multiconformer_fixed_size_counted_distance_metric[mcfsc] = distance_metric;
        cerr << "Multiconformer FSC " << mcfsc << " distance metric " << distance_metric << endl;
        mcfsc++;
      }
      else if (MULTICONFORMER_SPARSE_FORM == ftype)
      {
        _multiconformer_sparse_weight[mcsparse] = weight;
        _multiconformer_sparse_fingerprint_tag[mcsparse] = tag;
        if (0 != distance_metric)
          _multiconformer_sparse_distance_metric[mcsparse] = distance_metric;
        cerr << "Multiconformer sparse " << mcsparse << " distance metric " << distance_metric << endl;
        mcsparse++;
      }
      else
      {
        if (ptwa.ascii())
          _fingerprint_type[nfp] = ASCII_01_FORM;
        else
          _fingerprint_type[nfp] = ftype;

        _fingerprint_tag[nfp] = tag;
        _fingerprint_weight[nfp] = weight;
        if (0 != distance_metric)
          _dense_fingerprint_distance_metric[nfp] = distance_metric;

        if (tag.starts_with("HX") && DAYLIGHT_COMPRESSED_FORM == _fingerprint_type[nfp])
        {
          cerr << "Automatic hex conversion for fingerprint '" << tag << "'\n";
          _fingerprint_type[nfp] = UPPERCASE_HEX_FORM;
        }

        nfp++;
      }
    }
  }

  if (cl.option_present('P'))
  {
    dash_P_specified = 1;

    const_IWSubstring p;
    for (int i = 0; cl.value('P', p, i); ++i)
    {
      if ("none" == p)
      {
        molecular_property_integer_tag = "";
        if (verbose)
          cerr << "Properties not processed\n";
      }
      else if (p.contains("desc="))
      {
        if (! parse_desc_specification(p, verbose))
        {
          cerr << "INvalid -P qualifier '" << p << "'\n";
          return 0;
        }
      }
      else if ("dice" == p)
      {
        _integer_property_distance_metric = INTEGER_PROPERTY_DISTANCE_METRIC_DICE;
      }
      else if (0 == p.length())
      {
        molecular_property_integer_tag = "";
        if (verbose)
          cerr << "Properties in '" << molecular_property_integer_tag << "'\n";
      }
      else if ("help" == p)
      {
        display_dash_p_options(cerr);
        exit (0);
      }
      else
      {
        const_IWSubstring tag;
        float weight = -1.0;
  
        int notused1 = 0;
        int notused2 = 0;
        double notused3;
        int notused4 = 0;

        if (! parse_tag_weight(p, tag, weight, notused1, notused2, notused3, notused4))
        {
          cerr << "Invalid -P qualifier '" << molecular_property_integer_tag << "'\n";
          return 0;
        }
  
        molecular_property_integer_tag = tag;
       _property_weight_integer = weight;
  
        if (verbose)
          cerr << "Molecular properties in '" << molecular_property_integer_tag << "' field, weight " << _property_weight_integer << endl;
      }
    }
  }
  else if (options_specified_by_user)
    molecular_property_integer_tag.resize(0);

  if (molecular_property_integer_tag.length())
    (void) initialise_properties_ratios();

// If they entered any fingerprints, normalise the weights.

//if (nfp || sfp || cl.option_present ('P'))      // they entered something
  if (nfp || sfp || cfp || mc01 || mcfsc || mcsparse)      // they entered something
  {
    int nset = 0;
    float wset = 0.0;

    for (int i = 0; i < bit_fingerprints_in_file; i++)
    {
      if (_fingerprint_weight[i] >= 0.0)    // was set by the user
      {
        nset++;
        wset += _fingerprint_weight[i];
        cerr << "Fingerprint " << i << " weight " << _fingerprint_weight[i] << endl;
      }
    }

    for (int i = 0; i < _number_sparse_fingerprints; i++)
    {
      if (_sparse_fingerprint_weight[i] >= 0.0)
      {
        nset++;
        wset += _sparse_fingerprint_weight[i];
        cerr << "Non colliding fingerprint " << i << " weight " << _sparse_fingerprint_weight[i] << endl;
      }
    }

    for (int i = 0; i < number_fixed_size_counted_fingerprints_in_file; i++)
    {
      if (_fixed_size_counted_fingerprint_weight[i] >= 0.0)
      {
        nset++;
        wset += _fixed_size_counted_fingerprint_weight[i];
        cerr << "Fixed counted fingerprint " << i << " weight " << _fixed_size_counted_fingerprint_weight[i] << endl;
      }
    }

    int nitems = nfp + sfp + cfp + mc01 + mcfsc + mcsparse;

    if (molecular_property_integer_tag.length())
    {
      nitems++;
      if (_property_weight_integer >= 0.0)    // was set by the user
      {
        nset++;
        wset += _property_weight_integer;
      }
    }

    if (molecular_property_continuous_tag.length())
    {
      nitems++;
      if (_property_weight_continuous >= 0.0)
      {
        nset++;
        wset += _property_weight_continuous;
      }
    }

    if (nset == nitems)    // all items have a pre-defined weight
    {
      if (wset > 1.0)
        cerr << "Warning, weights sum to " << wset << endl;
    }
    else if (wset > 1.0)
    {
      cerr << "Sorry, the weights entered are > 1.0, sum = " << wset << endl;
      return 0;
    }
    else      // some entered without a weight
    {
      wset = (1.0 - wset) / static_cast<float>(nitems - nset);

      for (int i = 0; i < bit_fingerprints_in_file; i++)
      {
        if (_fingerprint_weight[i] < 0.0)
          _fingerprint_weight[i] = wset;
      }
      for (int i = 0; i < _number_sparse_fingerprints; i++)
      {
        if (_sparse_fingerprint_weight[i] < 0.0)
          _sparse_fingerprint_weight[i] = wset;
      }
      for (int i = 0; i < number_fixed_size_counted_fingerprints_in_file; i++)
      {
        if (_fixed_size_counted_fingerprint_weight[i] < 0.0)
          _fixed_size_counted_fingerprint_weight[i] = wset;
      }
      for (int i = 0; i < number_multiconformer_01_fingerprints_in_file; i++)
      {
        if (_multiconformer_01_weight[i] < 0.0)
          _multiconformer_01_weight[i] = wset;
      }
      for (int i = 0; i < number_multiconformer_fixed_size_counted_fingerprints_in_file; i++)
      {
        if (_multiconformer_fixed_size_counted_weight[i] < 0.0)
          _multiconformer_fixed_size_counted_weight[i] = wset;
      }
      for (int i = 0; i < number_multiconformer_sparse_fingerprints_in_file; i++)
      {
        if (_multiconformer_sparse_weight[i] < 0.0)
          _multiconformer_sparse_weight[i] = wset;
      }
    }

    if (_property_weight_integer < 0.0 && molecular_property_integer_tag.length())
      _property_weight_integer = wset;

    if (_property_weight_continuous < 0.0 && molecular_property_continuous_tag.length())
      _property_weight_continuous = wset;

    if (verbose)
      print_tags_and_weights(cerr);
  }

  if (cl.option_present('W'))
  {
    if (0 == molecular_property_integer_tag.length())
    {
      cerr << "Window specified, but no molecular properties, windows use MPR\n";
      return 0;
    }

    int minxtra = -1;
    int maxxtra = -1;

    windows_present = 0;

    const_IWSubstring w;
    int i = 0;
    while (cl.value('W', w, i++))
    {
      if (w.starts_with("R:"))
      {
        w += 2;
        int window;
        if (! w.numeric_value(window) || window < 0 || window > 100)
        {
          cerr << "The '-w R:' option must be followed by a whole percentage\n";
          display_standard_gfp_options(cerr);
          return 0;
        }

        set_ring_count_window(window, verbose);
        windows_present++;
      }
      else if (w.starts_with("A:minextra="))
      {
        w.remove_leading_chars(11);
        if (! w.numeric_value(minxtra) || minxtra < 0)
        {
          cerr << "Invalid min extra atoms 'A:minextra=" << w << "'\n";
          return 0;
        }
        if (verbose)
          cerr << "Will only compare molecules with at least " << minxtra << " extra atoms\n";
        windows_present++;
      }
      else if (w.starts_with("A:maxextra="))
      {
        w.remove_leading_chars(11);
        if (! w.numeric_value(maxxtra) || maxxtra < 0)
        {
          cerr << "Invalid max extra atoms 'A:maxextra=" << w << "'\n";
          return 0;
        }
        if (verbose)
          cerr << "Will only compare molecules with at most " << maxxtra << " extra atoms\n";
        windows_present++;
      }
      else if (w.starts_with("A:"))
      {
        w += 2;
        int window;
        if (! w.numeric_value(window) || window < 0 || window > 100)
        {
          cerr << "The '-w A:' option must be followed by a whole percentage\n";
          display_standard_gfp_options(cerr);
          return 0;
        }

        set_atom_count_window(window, verbose);
        windows_present++;
      }
      else if (w.starts_with("a:"))
      {
        w += 2;
        int window;

        if (! w.numeric_value(window) || window < 0)
        {
          cerr << "The '-W a:' option must be followed by a valid number\n";
          display_standard_gfp_options(cerr);
          return 0;
        }

        if (verbose)
          cerr << "Aromatic atom count window set to " << window << endl;
        set_aromatic_atom_count_window(window);
        windows_present++;
      }
      else
      {
        cerr << "Unrecognised -W qualifier '" << w << "'\n";
        display_standard_gfp_options(cerr);
        exit(9);
      }
    }

    if (maxxtra >= 0 || minxtra >= 0)
    {
      if (atom_count_window_present)
      {
        cerr << "Sorry, cannot have a percentage based atom count window as well as minextra/maxextra\n";
        return 0;
      }
      if (maxxtra >= 0 && minxtra >= 0 && minxtra > maxxtra)
      {
        cerr << "Inconsistent values for max extra atoms " << maxxtra << " and min extr atoms " << minxtra << endl;
        return 0;
      }

      set_atom_count_delta_window(minxtra, maxxtra);
    }
  }

// The guard fingerprint specification must be of the form 'TAG,dist'

  if (cl.option_present('G'))
  {
    const_IWSubstring g = cl.string_value('G');
    if (! parse_guard_fingerprint_specification(g))
    {
      cerr << "IW_General_Fingerprint::initialise: invalid guard fingerprint specification '" << g << "'\n";
      return 0;
    }
  }

  if (number_multiconformer_01_fingerprints_in_file || number_multiconformer_fixed_size_counted_fingerprints_in_file || number_multiconformer_sparse_fingerprints_in_file)
    multiconformer_fingerprints_present = 1;

  return 1;
}

int
need_to_call_initialise_fingerprints (const Command_Line & cl)
{
  if (cl.option_present('F'))
    return 1;

  if (cl.option_present('P'))
    return 1;

  if (cl.option_present('W'))
    return 1;

  if (cl.option_present('Q'))
    return 1;

  if (cl.option_present('G'))
    return 1;

  return 0;
}

IW_GFP_D::IW_GFP_D()
{
  _distance = static_cast<similarity_type_t>(1.0);

  return;
}

static int
initialise_fingerprints_2 (iwstring_data_source & input,
                           int verbose)
{
  IW_TDT tdt;
  while (tdt.next(input))
  {
    IW_General_Fingerprint fp;
    int fatal;

    if (fp.construct_from_tdt(tdt, fatal))   // great!
      return 1;

    if (fatal)
    {
      cerr << "initialise_fingerprints: bad TDT\n";
      cerr << tdt;
      return 0;
    }
  }

  return 0;
}

/*
  We need to make sure this is thread safe.

  We have two variables. An int, which should be cheap to check,
  and if that it set, then we need to check for serialisation

  Note that this is quite sub-optimal. Depending on the timing,
  potentially each thread is forced into serial execution, whereas
  all we need to do is ensure that the first thread completes
  before anyone else starts. But this is executed once and this
  is nice and simple so should not cause problems...
*/

#ifdef GFP_PARALLEL

static int need_to_check_mutex = 1;

static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

#endif

static int
initialise_fingerprints (iwstring_data_source & input,
                         int verbose)
{
#ifdef GFP_PARALLEL
  int need_to_unlock_mutex = 0;

  if (need_to_check_mutex)
  {
    pthread_mutex_lock(&mutex);
    need_to_unlock_mutex = 1;
  }
#endif

  int rc = initialise_fingerprints_2(input, verbose);
  
#ifdef GFP_PARALLEL
  if (need_to_unlock_mutex)
  {
    need_to_check_mutex = 0;
    pthread_mutex_unlock(&mutex);
  }
#endif

  return rc;
}

int
initialise_fingerprints (const char * fname, int verbose)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "initialise_fingerprints: cannot open '" << fname << "'\n";
    return 0;
  }

  return initialise_fingerprints(input, verbose);
}

int
initialise_fingerprints (const IWString & fname, int verbose)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "initialise_fingerprints: cannot open '" << fname << "'\n";
    return 0;
  }

  return initialise_fingerprints(input, verbose);
}

int
IW_General_Fingerprint::print_all_bits (std::ostream & output) const
{
  output << _id << " has ";
  if (bit_fingerprints_in_file)
    output << bit_fingerprints_in_file << " fingerprints";
  if (_number_sparse_fingerprints)
    output << _number_sparse_fingerprints << " sparse fingerprints";
  output << endl;

  for (int i = 0; i < bit_fingerprints_in_file; i++)
  {
    output << i << " " << _fingerprint_tag[i] << ' ' << _fingerprint[i].nbits() << " bits " << _fingerprint[i].nset() << " bits set\n";

    _fingerprint[i].printon(output);

    output << endl;
  }

  for (int i = 0; i < _number_sparse_fingerprints; i++)
  {
    _sparse_fingerprint[i].debug_print(output);
  }

  for (int i = 0; i < number_multiconformer_01_fingerprints_in_file; i++)
  {
    _multiconformer_01[i].debug_print(output);
  }

  for (int i = 0; i < number_multiconformer_fixed_size_counted_fingerprints_in_file; i++)
  {
    _multiconformer_fixed_size_counted[i].debug_print(output);
  }

  for (int i = 0; i < number_multiconformer_sparse_fingerprints_in_file; i++)
  {
    _multiconformer_sparse[i].debug_print(output);
  }

  return output.good();
}

int
IW_General_Fingerprint::convert_to_non_sparse_forms (const Set_of_Sparse_Fingerprint_Collection_Profile & sfpcp)
{
  if (convert_sparse_fingerprint_to_bits)
    _convert_sparse_fingerprints_to_bits(sfpcp);

  if (convert_sparse_fingerprint_to_fixed_width_counted)
    _convert_sparse_fingerprints_to_fixed_width_counted(sfpcp);

  delete [] _sparse_fingerprint;
  _sparse_fingerprint = NULL;

  return 1;
}

int
IW_General_Fingerprint::_convert_sparse_fingerprints_to_bits (const Set_of_Sparse_Fingerprint_Collection_Profile & sfpcp)
{
  for (int i = 0; i < _number_sparse_fingerprints; i++)
  {
    int j = bit_fingerprints_in_file + i;

    int extra_bits;

    sfpcp.convert_to_fixed_width(i, _sparse_fingerprint[i], _fingerprint[j], extra_bits);

//  if (extra_bits)
//    _fingerprint[j].set_nset (_fingerprint[j].nset() + extra_bits);
//  else
//    _fingerprint[j].compute_nset();
  }

  return 1;
}

int
IW_General_Fingerprint::_convert_sparse_fingerprints_to_fixed_width_counted (const Set_of_Sparse_Fingerprint_Collection_Profile & sfpcp)
{
  for (int i = 0; i < _number_sparse_fingerprints; i++)
  {
    int extra_bits, extra_count;

    int j = number_fixed_size_counted_fingerprints_in_file + i;

    sfpcp.convert_to_fixed_width(i, _sparse_fingerprint[i], _fixed_size_counted_fingerprint[j], extra_bits, extra_count);

//  if (extra_count)
//    _fixed_size_counted_fingerprint[j].set_nset (_fixed_size_counted_fingerprint[j].nset() + extra_count);
  }

  return 1;
}

int 
IW_General_Fingerprint::number_sparse_fingerprints(void) const
{
  return _number_sparse_fingerprints;
}

float
IW_General_Fingerprint::equal_weight_tanimoto(IW_General_Fingerprint & rhs)
{
  int bic = 0;
  int na = 0;
  int nb = 0;

  if (_molecular_properties_integer.active())
  {
    na += _molecular_properties_integer.sum();
    nb += rhs._molecular_properties_integer.sum();

    bic += _molecular_properties_integer.bits_in_common(rhs._molecular_properties_integer);
  }

  for (int i = 0; i < number_fingerprints_to_use_in_computations; i++)
  {
    na += _fingerprint[i].nset();
    nb += rhs._fingerprint[i].nset();

    bic += _fingerprint[i].bits_in_common(rhs._fingerprint[i]);
  }

  for (int i = 0; i < _number_sparse_fingerprints_to_use_in_computations; i++)
  {
    na += _sparse_fingerprint[i].nset();
    nb += rhs._sparse_fingerprint[i].nset();

    bic += _sparse_fingerprint[i].bits_in_common(rhs._sparse_fingerprint[i]);
  }

  if (0 == bic)
  {
    if (0 == na && 0 == nb)      // otherwise we have identical molecules that will have non-zero distances
      return static_cast<similarity_type_t> (1.0);

    return static_cast<similarity_type_t> (0.0);
  }

  return static_cast<double>(bic) / static_cast<double>(na + nb - bic);
}

float
IW_General_Fingerprint::equal_weight_dot_product(IW_General_Fingerprint & rhs)
{
  int rc = 0;

  if (_molecular_properties_integer.active())
  {
//  rc += _molecular_properties_integer.bits_in_common(rhs._molecular_properties_integer);   wrong, but this will never be used, implement if ever needed
  }

  for (int i = 0; i < number_fingerprints_to_use_in_computations; i++)
  {
    rc += _fingerprint[i].bits_in_common(rhs._fingerprint[i]);
  }

  for (int i = 0; i < _number_sparse_fingerprints_to_use_in_computations; i++)
  {
    rc += _sparse_fingerprint[i].dot_product(rhs._sparse_fingerprint[i]);
  }

  return static_cast<float>(rc);
}


const IWString &
fixed_fingerprint_tag (int i)
{
  return _fingerprint_tag[i];
}

const IWString & 
sparse_fingerprint_tag (int i)
{
  return _sparse_fingerprint_tag[i];
}

#ifdef FB_ENTROPY_WEIGHTED_FPS

int
IW_General_Fingerprint::convert_01_fingerprint_to_integer_molecular_properties()
{
  if (1 != bit_fingerprints_in_file)
  {
    cerr << "IW_General_Fingerprint::convert_01_fingerprint_to_integer_molecular_properties:must be one, and only one 01 fingerprint in the input, " << bit_fingerprints_in_file << " impossible\n";
    assert (NULL == "Cannot continue");
  }

  if (_property_weight_integer <= 0.0)
  {
    _property_weight_integer = 1.0;
    _fingerprint_weight[0] = 0.0;
  }

  assert (! _molecular_properties_integer.active());

  if (! _molecular_properties_integer.construct_from_bit_vector (_fingerprint[0]))
  {
    cerr << "IW_General_Fingerprint::convert_01_fingerprint_to_integer_molecular_properties:error '" << _id << "'\n";
    return 0;
  }

  number_fingerprints_to_use_in_computations = 0;

  return 1;
}

int
Molecular_Properties_Integer::construct_from_bit_vector (const IWDYFP & fp)
{
  assert (NULL == _property);

  _nproperties = fp.nbits();

  assert (_nproperties > 0);

  _property = new int[_nproperties];

  fp.set_vector(_property);

  return 1;
}

#endif

int
IW_General_Fingerprint::fixed_fingerprints_as_hex (IWString & s) const
{
  for (int i = 0; i < bit_fingerprints_in_file; i++)
  {
    IWString tmp;
    _fingerprint[i].hex_form(tmp);

    if (i > 0)
      s << ',';

    s << tmp;
  }

  return 1;
}

int
IW_General_Fingerprint::debug_print(std::ostream & output) const
{
  output << "IW_General_Fingerprint::debug_print: " << _id << endl;

  return 1;
}

void
delete_gfp_file_scope_static_objects ()
{
  if (NULL != _fingerprint_weight)
    delete [] _fingerprint_weight;
  if (NULL != _sparse_fingerprint_weight)
    delete [] _sparse_fingerprint_weight;
  if (NULL != _fixed_size_counted_fingerprint_weight)
    delete [] _fixed_size_counted_fingerprint_weight;
  if (NULL != _bits_in_fingerprint)
    delete [] _bits_in_fingerprint;

  if (NULL != _fingerprint_type)
    delete [] _fingerprint_type;
  if (NULL != _sparse_fingerprint_type)
    delete [] _sparse_fingerprint_type;

  if (NULL != _fingerprint_tag)
    delete [] _fingerprint_tag;
  if (NULL != _sparse_fingerprint_tag)
    delete [] _sparse_fingerprint_tag;
  if (NULL != _fixed_size_counted_fingerprint_tag)
    delete [] _fixed_size_counted_fingerprint_tag;
  if (NULL != _multiconformer_01_fingerprint_tag)
    delete [] _multiconformer_01_fingerprint_tag;
  if (NULL != _multiconformer_fixed_size_counted_fingerprint_tag)
    delete [] _multiconformer_fixed_size_counted_fingerprint_tag;
  if (NULL != _multiconformer_sparse_fingerprint_tag)
    delete [] _multiconformer_sparse_fingerprint_tag;

  if (NULL != _sparse_fingerprint_scaling_factor)
    delete [] _sparse_fingerprint_scaling_factor;
  if (NULL != _sparse_fingerprint_distance_metric)
    delete [] _sparse_fingerprint_distance_metric;

  if (NULL != _dense_fingerprint_distance_metric)
    delete [] _dense_fingerprint_distance_metric;


  return;
}
