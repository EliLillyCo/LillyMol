/*
  Converts a gfp file to input for svm lite
*/

#include <stdlib.h>
#ifdef _WIN32
#else
#include <unistd.h>
#endif
#include <math.h>

#include <iostream>
#include <memory>
#include <random>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwmmap.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "gfp_bit_subset.h"

#include "Utilities/GFP_Tools/gfp.h"

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

static IW_STL_Hash_Map_String activity;
static int activity_data_read_from_file = 0;

static int activity_column = 1;  // by default identifier in column 0

// Many times we do not need activity values, and instead we can write
// random values. I did think about a fixed value, but that seems riskier.
static int write_random_activity_values = 0;

static int tdts_read = 0;

static int ignore_missing_activity = 0;

static unsigned int highest_sparse_bit_number = 0;

static unsigned int lowest_sparse_bit_number = 0;

typedef IW_Hash_Map<unsigned int, int> Bit_to_Feature_Map;

static int identifiers_truncated_to_first_token = 0;

static int strip_leading_zeros = 0;

static int produce_normalised_output = 0;

/*
  If descriptors are present, we can write them in either of two ways.
  Standard svm-lite input format with feature_number:value or just
  as descriptors
*/

static int write_descriptor_data_as_descriptors = 0;

/*
  svmLite wants an identifier at the end of each record,
  LIBSVM does not
*/

static int append_identifier_at_end_of_each_record = 1;

static int max_count = 255;

static Accumulator_Int<int> acc_count;

static IWDigits bit_number_iwdigits;
static IWDigits count_iwdigits;

/*
  It can be informative to know the number of bits not present
  in the training set for each molecule
*/

static IWString_and_File_Descriptor bnpts;

static int processing_difference_fingerprint_offset = 0;

/*
  April 2015. A lot of time is spent in the kernel function counting the
  number of features set. We can save time by using a convention that the
  first feature is not really a feature, but is the total for that document
*/

static int first_features_are_nbits_and_weight = 0;

static int warn_non_numeric_activity = 1;
static int non_numeric_activity_encountered = 0;

static int flush_after_every_molecule = 0;

template <typename T>
class Feature_Value_Pair_Template
{
 private:
  unsigned int _feature_number;
  T _v;

 public:
  Feature_Value_Pair_Template(const unsigned int b, const T v) : _feature_number(b), _v(v)
  {
  }

  Feature_Value_Pair_Template(const Feature_Value_Pair_Template<T>& rhs)
  {
    _feature_number = rhs._feature_number;
    _v = rhs._v;
  }

  Feature_Value_Pair_Template()
  {
    _feature_number = 0;
    _v = 0;
  }

  void
  set_feature_number(unsigned int s)
  {
    _feature_number = s;
  }

  void
  set_value(T v)
  {
    _v = v;
  }

  unsigned int
  feature_number() const
  {
    return _feature_number;
  }

  T
  zvalue() const
  {
    return _v;
  }
};

typedef class Feature_Value_Pair_Template<int> Feature_Value_Pair;

/*
  The Intel machine learning tools need a CSR file format
  Compressed Sparse Row
*/

template <typename T>
class Item_For_Intel_Template
{
 private:
  IWString _id;
  unsigned int* _b;
  T* _v;
  int _n;

 public:
  Item_For_Intel_Template();
  ~Item_For_Intel_Template();

  void
  set_id(const const_IWSubstring& s)
  {
    _id = s;
  }

  const IWString&
  id() const
  {
    return _id;
  }

  int
  extra(const Feature_Value_Pair_Template<T>*, const int n);

  int
  nbits() const
  {
    return _n;
  }

  int
  fill_column_array(int* c) const;

  int
  write_values(const int ndx, IWString_and_File_Descriptor& output) const;
  int
  write_column_array(const int ndx, IWString_and_File_Descriptor& output) const;
};

typedef Item_For_Intel_Template<int> Item_For_Intel;

template <typename T>
Item_For_Intel_Template<T>::Item_For_Intel_Template()
{
  _b = nullptr;
  _v = nullptr;
  _n = 0;
}

template <typename T>
Item_For_Intel_Template<T>::~Item_For_Intel_Template()
{
  if (nullptr != _b) {
    delete[] _b;
  }

  if (nullptr != _v) {
    delete[] _v;
  }
}

template <>
int
Item_For_Intel_Template<int>::write_values(const int ndx,
                                           IWString_and_File_Descriptor& output) const
{
  if (0 == _n) {
    cerr << "Item_For_Intel_Template::write_values:no bits\n";
    return 0;
  }

  int istart;
  if (0 == ndx) {
    istart = 1;
    output << _v[0];
  } else {
    istart = 0;
  }

  for (int i = istart; i < _n; ++i) {
    bit_number_iwdigits.append_number(output, _v[i]);
  }

  return 1;
}

template <typename T>
int
Item_For_Intel_Template<T>::write_column_array(const int ndx,
                                               IWString_and_File_Descriptor& output) const
{
  if (0 == _n) {
    cerr << "Item_For_Intel_Template::write_column_array:no bits\n";
    return 0;
  }

  int istart;
  if (0 == ndx) {
    istart = 1;
    output << _b[0];
  } else {
    istart = 0;
  }

  for (int i = istart; i < _n; ++i) {
    bit_number_iwdigits.append_number(output, _b[i]);
  }

  return 1;
}

#ifdef FINISH_THIS_SOMETIME_PERHAPS
template <typename T>
int
Item_For_Intel_Template<T>::fill_column_array(int* c) const
{
  const int n = _bit.number_elements();

  for (int i = 0; i < n; ++i) {
    const auto* fvpi = _bit[i];
    c[i] = fvpi->feature_number();
  }

  return n;
}
#endif

template <typename T>
int
Item_For_Intel_Template<T>::extra(const Feature_Value_Pair_Template<T>* fvp, const int n)
{
  _n = n;
  _b = new unsigned int[n];
  _v = new T[n];

  for (int i = 0; i < n; ++i) {
    _b[i] = fvp[i].feature_number();
    _v[i] = fvp[i].zvalue();
  }

  return n;
}

template <typename T>
class Intel_CSR
{
 private:
  Item_For_Intel_Template<T>* _s;
  int _capacity;
  int _n;

  int _index_base;

 public:
  Intel_CSR();
  ~Intel_CSR();

  int
  active() const
  {
    return _capacity;
  }

  int
  activate(const int n);

  void
  set_index_base(const int s)
  {
    _index_base = s;
  }

  //  int extra (Item_For_Intel_Template<T> * s) { _s.add(s);}
  int
  extra(const Feature_Value_Pair_Template<T>* fvp, const int n);

  int
  do_write(const int highest_feature_number, IWString_and_File_Descriptor& output) const;
};

template <typename T>
Intel_CSR<T>::Intel_CSR()
{
  _capacity = 0;
  _n = -1;
  _s = nullptr;

  _index_base = 1;

  return;
}

template <typename T>
Intel_CSR<T>::~Intel_CSR()
{
  if (nullptr != _s) {
    delete[] _s;
  }

  return;
}

template <typename T>
int
Intel_CSR<T>::activate(const int n)
{
  _s = new Item_For_Intel_Template<T>[n];

  _capacity = n;

  _n = 0;

  return 1;
}

template <typename T>
int
Intel_CSR<T>::extra(const Feature_Value_Pair_Template<T>* fvp, const int n)
{
  Item_For_Intel_Template<T>& x = _s[_n];

  _n++;

  assert(_n <= _capacity);

  x.extra(fvp, n);

  return 1;
}

template <typename T>
int
Intel_CSR<T>::do_write(const int highest_feature_number,
                       IWString_and_File_Descriptor& output) const
{
  int* rowindex = new int[_n];
  std::unique_ptr<int[]> free_rowindex(rowindex);

  int ndx = 0;
  for (int i = 0; i < _n; ++i) {
    rowindex[i] = ndx;
    ndx += _s[i].nbits();
  }

  output << _index_base;
  for (int i = 1; i < _n; ++i) {
    bit_number_iwdigits.append_number(output, rowindex[i] + _index_base);
    output.write_if_buffer_holds_more_than(32768);
  }
  output << ',' << (ndx + _index_base) << '\n';

  for (int i = 0; i < _n; ++i) {
    _s[i].write_column_array(i, output);
    output.write_if_buffer_holds_more_than(32768);
  }
  output << '\n';

  for (int i = 0; i < _n; ++i) {
    _s[i].write_values(i, output);
    output.write_if_buffer_holds_more_than(32768);
  }

  output << '\n';

  return 1;
}

static Intel_CSR<int> intel_csr;
static IWString_and_File_Descriptor stream_for_activity_values;

static void
usage(int rc)
{
  // clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Converts a gfp file to input for svm-lite\n";
  cerr << " -A <fname>     activity file, use 'RANDOM' to write arbitrary values instead\n";
  cerr << " -c <col>       activity in column <col> of activity file\n";
  cerr << " -p <support>   support level for inclusion\n";
  cerr << " -m <max>       truncate counted bits at <max>\n";
  cerr << " -F <tag>       tag to process - standard gfp_* syntax\n";
  cerr << " -C <fname>     create a bit->feature number cross reference file\n";
  cerr << " -U <fname>     use an existing bit->feature number cross reference file\n";
  cerr << " -N <fname>     for each molecule write the number of bits not present in training set\n";
  cerr << " -D <fname>     merge in descriptor file <fname>\n";
  cerr << " -d             when processing descriptors, write them as descriptors\n";
  cerr << " -g             ignore missing descriptors - don't write any data for that id\n";
  cerr << " -V             output is for LIBSVM, suppress addition of identifier\n";
  cerr << " -I <fname>     output is CSR for Intel. <fname> will contain activity values in order\n";
  cerr << " -u             write output as normalised unit vectors\n";
  cerr << " -X <offset>    processing a difference fingerprint, values offset around <offset>\n";
  cerr << " -Y ...         more options, enter '-Y help' for info\n";
  cerr << " -t             first feature written is the total count\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

class Count_and_Column
{
 private:
  unsigned int _feature_number;
  int _count;
  int _column;

 public:
  Count_and_Column(int col) : _count(1), _column(col)
  {
  }

  void in_another_molecule() {
    _count++;
  }

  int column() const {
    return _column;
  }

  int count() const {
    return _count;
  }
};

class Bit_Data
{
 private:
  typedef IW_Hash_Map<unsigned int, Count_and_Column*> B2CC;

  B2CC _bit_to_column_and_count;

 public:
  int bit_found(unsigned int b);
};

int
Bit_Data::bit_found(unsigned int b)
{
  B2CC::iterator f = _bit_to_column_and_count.find(b);

  if (f != _bit_to_column_and_count.end()) {
    Count_and_Column* cc = (*f).second;

    cc->in_another_molecule();

    return 1;
  }

  unsigned int s = _bit_to_column_and_count.size();

  _bit_to_column_and_count[b] = new Count_and_Column(s);

  return 1;
}

class Column_Number_Comparator
{
 public:
  int
  operator()(const Count_and_Column&, const Count_and_Column&) const;
};

class Descriptor_File
{
 private:
  IWString _header;

  iwstring_data_source _input;

  //  We record the offset inthe file of each identifier

  IW_STL_Hash_Map<IWString, off_t> _id_to_offset;

 public:
  int build(const char* fname);
  int build();

  const IWString& header() const {
    return _header;
  }

  int fetch_data_for_id(const IWString& id, const_IWSubstring& zdata);
};

int
Descriptor_File::build(const char* fname)
{
  if (!_input.open(fname)) {
    cerr << "Descriptor_File::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build();
}

int
Descriptor_File::build()
{
  assert(_input.is_open());

  if (!_input.next_record(_header)) {
    cerr << "Descriptor_File::build:cannot read header record\n";
    return 0;
  }

  off_t o = _input.tellg();
  const_IWSubstring buffer;
  while (_input.next_record(buffer)) {
    buffer.truncate_at_first(' ');

    _id_to_offset[buffer] = o;

    o = _input.tellg();
  }

  if (0 == _id_to_offset.size()) {
    cerr << "Descriptor_File::build:no data\n";
    return 0;
  }

  if (verbose) {
    cerr << "Descriptor_File::build:read offsets for " << _id_to_offset.size()
         << " identifiers\n";
  }

  return _id_to_offset.size();
}

int
Descriptor_File::fetch_data_for_id(const IWString& id, const_IWSubstring& zdata)
{
  IW_STL_Hash_Map<IWString, off_t>::const_iterator f = _id_to_offset.find(id);

  if (f == _id_to_offset.end()) {
    cerr << "Descriptor_File::fetch_data_for_id:no offset for '" << id << "'\n";
    return 0;
  }

  if (!_input.seekg((*f).second)) {
    cerr << "Descriptor_File::fetch_data_for_id:cannot seek to " << (*f).second
         << " for '" << id << "'\n";
    return 0;
  }

  if (!_input.next_record(zdata)) {
    cerr << "Descriptor_File::fetch_data_for_id:cannot fetch record for '" << id << "'\n";
    return 0;
  }

  zdata.remove_leading_words(1);

  return 1;
}

static int
write_descriptor_data(const const_IWSubstring& zdata,
                      IWString_and_File_Descriptor& output)
{
  if (write_descriptor_data_as_descriptors) {
    output << ' ' << zdata;
    return 1;
  }

  float float_zero = static_cast<float>(0.0);

  int i = 0;
  const_IWSubstring token;
  for (int feature_number = 1; zdata.nextword(token, i); feature_number++) {
    if ('0' == token) {
      continue;
    }

    float v;
    if (!token.numeric_value(v)) {
      cerr << "write_descriptor_data:invalid numeric '" << token << "'\n";
      return 0;
    }

    if (float_zero == v) {
      continue;
    }

    bit_number_iwdigits.append_number(output, feature_number);
    output << v;
  }

  return 1;
}

static Descriptor_File* dfile = nullptr;

static int ignore_missing_descriptor_data = 0;

static int missing_descriptor_data_ignored = 0;

/*
  This is a two pass process.
  First we read the whole fingerprint file and determine which bits
  are present, and in how many molecules is each bit set.
  Once we have decided which bits match the support criterion, we
  then use the integer as the feature number

  In general, a mapping must consider each of the fixed width fingerprints
  and each of the sparse fingerprints

  We can have a descriptor file as well. In that case, we record the
  header and assume that features are assigned sequentially according
  to the descriptors in the header
*/

typedef IW_Hash_Map<unsigned int, int> Bit_and_Integer;

class GFP_to_Feature_Map
{
 private:
  int _next_offset_to_assign;

  IWString _header;

  Bit_and_Integer _properties;

  int _nproperties;

  int _number_fixed;
  Bit_and_Integer* _bai_fixed;
  int* _nbits_fixed;

  int _number_sparse;
  Bit_and_Integer* _bai_sparse;

  //  We keep track of the largest number of bits set in any given fingerprint

  int _highest_number_features;

  //  When being built from a set of fingerprints, we need to keep track of
  //  the whether this is the first one or not

  int _needs_to_be_initialised;

  //  private functions

  int
  _profile_bits_fixed_width(const IWDYFP& fp, Bit_and_Integer& bai,
                            int& features_this_fingerprint);
  int
  _profile_bits_sparse(const Sparse_Fingerprint& fp, Bit_and_Integer& bai,
                       int& features_this_fingerprint);
  int
  _profile_properties(const Molecular_Properties_Integer& mpr, Bit_and_Integer&,
                      int& features_this_fingerprint);

  int
  _gather_bits_set_properties(const Molecular_Properties_Integer& mpr,
                              Feature_Value_Pair*, int& ndx,
                              int& bits_not_present_in_training_set) const;
  int
  _gather_bits_set_fixed(const IWDYFP& fp, const Bit_and_Integer&,
                         Feature_Value_Pair* fvp, int& ndx,
                         int& bits_not_present_in_training_set) const;
  int
  _gather_bits_set_sparse(const Sparse_Fingerprint& fp, const Bit_and_Integer&,
                          Feature_Value_Pair* fvp, int& ndx,
                          int& bits_not_present_in_training_set) const;

  void
  _assign_offsets(Bit_and_Integer& bai);
  void
  _assign_offsets(Bit_and_Integer& bai, int support,
                  int& bits_suppressed_by_support_requirement);

 public:
  GFP_to_Feature_Map();
  ~GFP_to_Feature_Map();

  int debug_print(std::ostream&) const;

  int initialise();

  int highest_number_features() const {
    return _highest_number_features;
  }

  int number_features_defined() const;

  int extract_information_from_descriptor_file(const Descriptor_File&);

  int profile_bits_set(const IW_General_Fingerprint&);

  //  Offsets will be different if there is a support level specified

  int assign_offsets();
  int assign_offsets(int support, int& bits_suppressed_by_support_requirement);

  int gather_bits_set(const IW_General_Fingerprint&, Feature_Value_Pair* fvp,
                  int& bits_not_present_in_training_set) const;  // second phase

  int write_mapping(const char* fname) const;
  int write_mapping(IWString_and_File_Descriptor& output) const;
  int read_existing_mapping(const char* fname);
  int read_existing_mapping(iwstring_data_source& input);
};

GFP_to_Feature_Map::GFP_to_Feature_Map()
{
  _next_offset_to_assign = 1;

  if (first_features_are_nbits_and_weight) {
    _next_offset_to_assign = 3;
  }

  _number_fixed = 0;
  _number_sparse = 0;

  _nproperties = 0;

  _bai_fixed = nullptr;
  _nbits_fixed = nullptr;

  _bai_sparse = nullptr;

  _highest_number_features = 0;

  _needs_to_be_initialised = 1;

  return;
}

GFP_to_Feature_Map::~GFP_to_Feature_Map()
{
  if (nullptr != _bai_fixed) {
    delete[] _bai_fixed;
  }

  if (nullptr != _bai_sparse) {
    delete[] _bai_sparse;
  }

  if (nullptr != _nbits_fixed) {
    delete[] _nbits_fixed;
  }

  return;
}

int
GFP_to_Feature_Map::debug_print(std::ostream& output) const
{
  output << "GFP_to_Feature_Map::debug_print\n";
  if (_nproperties) {
    output << _nproperties << " integer properties\n";
  }
  output << _number_fixed << " fixed fingerprints and " << _number_sparse
         << " sparse fingerprints\n";

  return output.good();
}

/*
  Calls a bunch of global functions to initialise info on numbers
  of fingerprints active
*/

int
GFP_to_Feature_Map::initialise()
{
  _number_fixed = number_fingerprints();

  if (_number_fixed) {
    _bai_fixed = new Bit_and_Integer[_number_fixed];

    _nbits_fixed = new_int(_number_fixed);
  }

  _number_sparse = number_sparse_fingerprints();

  if (_number_sparse) {
    _bai_sparse = new Bit_and_Integer[_number_sparse];
  }

  return 1;
}

int
GFP_to_Feature_Map::number_features_defined() const
{
  int rc = _properties.size();

  for (int i = 0; i < _number_fixed; i++) {
    rc += _bai_fixed[i].size();
  }

  for (int i = 0; i < _number_sparse; i++) {
    rc += _bai_sparse[i].size();
  }

  return rc;
}

int
GFP_to_Feature_Map::profile_bits_set(const IW_General_Fingerprint& gfp)
{
  if (_needs_to_be_initialised) {
    if (!initialise()) {
      return 0;
    }

    _needs_to_be_initialised = 0;
  }

  int features_this_fingerprint = 0;

  const Molecular_Properties_Integer& mpi = gfp.molecular_properties_integer();
  if (mpi.active()) {
    _profile_properties(mpi, _properties, features_this_fingerprint);
  }

  assert(_number_fixed == gfp.nfingerprints());

  for (int i = 0; i < _number_fixed; i++) {
    const IWDYFP& fpi = gfp[i];

    _nbits_fixed[i] = fpi.nbits();

    _profile_bits_fixed_width(fpi, _bai_fixed[i], features_this_fingerprint);
  }

  assert(_number_sparse == gfp.number_sparse_fingerprints());

  for (int i = 0; i < _number_sparse; i++) {
    const Sparse_Fingerprint& fpi = gfp.sparse_fingerprint(i);

    _profile_bits_sparse(fpi, _bai_sparse[i], features_this_fingerprint);
  }

  if (features_this_fingerprint > _highest_number_features) {
    _highest_number_features = features_this_fingerprint;
  }

  return 1;
}

static int
increment_count(Bit_and_Integer& bai, unsigned int b)
{
  Bit_and_Integer::iterator f = bai.find(b);

  if (f == bai.end()) {
    bai[b] = 1;
  } else {
    (*f).second += 1;
  }

  return 1;
}

int
GFP_to_Feature_Map::_profile_bits_fixed_width(const IWDYFP& fp, Bit_and_Integer& bai,
                                              int& features_this_fingerprint)
{
  int b;
  int i = 0;

  while ((b = fp.next_on_bit(i)) >= 0) {
    //  cerr << "Fixed bit " << b << " turned on\n";
    increment_count(bai, b);
    features_this_fingerprint++;
  }

  return 1;
}

int
GFP_to_Feature_Map::_profile_bits_sparse(const Sparse_Fingerprint& fp,
                                         Bit_and_Integer& bai,
                                         int& features_this_fingerprint)
{
  int i = 0;
  unsigned int b;
  int c;

  while (fp.next_bit_set(i, b, c)) {
    increment_count(bai, b);
    features_this_fingerprint++;
  }

  return 1;
}

int
GFP_to_Feature_Map::_profile_properties(const Molecular_Properties_Integer& mpi,
                                        Bit_and_Integer& bai,
                                        int& features_this_fingerprint)
{
  _nproperties = mpi.nproperties();

  const int* p = mpi.rawdata();

  for (int i = 0; i < _nproperties; i++) {
    if (p[i]) {
      increment_count(bai, i);
      features_this_fingerprint++;
    }
  }

  return 1;
}

void
GFP_to_Feature_Map::_assign_offsets(Bit_and_Integer& bai)
{
  // for (Bit_and_Integer::iterator i = bai.begin(); i != bai.end(); ++i)
  for (auto i = bai.begin(); i != bai.end(); ++i) {
    (*i).second = _next_offset_to_assign;

    //  cerr << "Bit " << (*i).first << " assigned " << offset << endl;
    _next_offset_to_assign++;
  }
}

int
GFP_to_Feature_Map::assign_offsets()
{
  assert(_next_offset_to_assign > 0);

  _assign_offsets(_properties);

  for (int i = 0; i < _number_fixed; i++) {
    _assign_offsets(_bai_fixed[i]);
  }

  for (int i = 0; i < _number_sparse; i++) {
    _assign_offsets(_bai_sparse[i]);
  }

  return _next_offset_to_assign - 1;
}

void
GFP_to_Feature_Map::_assign_offsets(Bit_and_Integer& bai, int support,
                                    int& bits_suppressed_by_support_requirement)
{
  // for (Bit_and_Integer::iterator i = bai.begin(); i != bai.end(); ++i)
  for (auto i = bai.begin(); i != bai.end(); ++i) {
    if ((*i).second >= support) {
      (*i).second = _next_offset_to_assign;
      _next_offset_to_assign++;
    } else {
      (*i).second = -1;
      bits_suppressed_by_support_requirement++;
    }
  }

  return;
}

int
GFP_to_Feature_Map::assign_offsets(int support,
                                   int& bits_suppressed_by_support_requirement)
{
  assert(_next_offset_to_assign > 0);

  _assign_offsets(_properties, 0, bits_suppressed_by_support_requirement);

  for (int i = 0; i < _number_fixed; i++) {
    _assign_offsets(_bai_fixed[i], support, bits_suppressed_by_support_requirement);
  }

  for (int i = 0; i < _number_sparse; i++) {
    _assign_offsets(_bai_sparse[i], support, bits_suppressed_by_support_requirement);
  }

  return _next_offset_to_assign - 1;  // number of features we have defined
}

int
GFP_to_Feature_Map::gather_bits_set(const IW_General_Fingerprint& fp,
                                    Feature_Value_Pair* fvp,
                                    int& bits_not_present_in_training_set) const
{
  bits_not_present_in_training_set = 0;

  int ndx = 0;

  const Molecular_Properties_Integer& mpi = fp.molecular_properties_integer();

  if (mpi.active()) {
    _gather_bits_set_properties(mpi, fvp, ndx, bits_not_present_in_training_set);
  }

  for (int i = 0; i < _number_fixed; i++) {
    _gather_bits_set_fixed(fp[i], _bai_fixed[i], fvp, ndx,
                           bits_not_present_in_training_set);
  }

  for (int i = 0; i < _number_sparse; i++) {
    _gather_bits_set_sparse(fp.sparse_fingerprint(i), _bai_sparse[i], fvp, ndx,
                            bits_not_present_in_training_set);
  }

  return ndx;
}

int
GFP_to_Feature_Map::_gather_bits_set_properties(
    const Molecular_Properties_Integer& mpr, Feature_Value_Pair* fvp, int& ndx,
    int& bits_not_present_in_training_set) const
{
  int n = mpr.nproperties();

  const int* v = mpr.rawdata();

  for (int i = 0; i < n; i++) {
    if (0 == v[i]) {
      continue;
    }

    Bit_and_Integer::const_iterator f = _properties.find(i);

    if (f == _properties.end()) {
      bits_not_present_in_training_set++;
    } else {
      //    cerr << "Property " << i << " value " << (*f).second << endl;
      fvp[ndx].set_feature_number((*f).second);
      fvp[ndx].set_value(v[i]);
      ndx++;
    }
  }

  return 1;
}

int
GFP_to_Feature_Map::_gather_bits_set_fixed(const IWDYFP& fp, const Bit_and_Integer& bai,
                                           Feature_Value_Pair* fvp, int& ndx,
                                           int& bits_not_present_in_training_set) const
{
  int b;
  int i = 0;

  while ((b = fp.next_on_bit(i)) >= 0) {
    Bit_and_Integer::const_iterator f = bai.find(b);

    if (f == bai.end()) {
      bits_not_present_in_training_set++;
    } else if ((*f).second < 0) {
      ;
    } else {
      fvp[ndx].set_feature_number((*f).second);
      fvp[ndx].set_value(1);
      ndx++;
    }
  }

  return 1;
}

int
GFP_to_Feature_Map::_gather_bits_set_sparse(const Sparse_Fingerprint& fp,
                                            const Bit_and_Integer& bai,
                                            Feature_Value_Pair* fvp, int& ndx,
                                            int& bits_not_present_in_training_set) const
{
  int i = 0;
  unsigned int b;
  int c;

  while (fp.next_bit_set(i, b, c)) {
    Bit_and_Integer::const_iterator f = bai.find(b);

    if (f == bai.end()) {
      bits_not_present_in_training_set++;
    } else if ((*f).second < 0) {
      ;
    } else {
      fvp[ndx].set_feature_number((*f).second);
      fvp[ndx].set_value(c);
      ndx++;
    }
  }

  return 1;
}

int
GFP_to_Feature_Map::extract_information_from_descriptor_file(const Descriptor_File& dfile)
{
  const IWString& h = dfile.header();

  _next_offset_to_assign = h.nwords();  // remember the header still contains column 1

  return 1;
}

int
GFP_to_Feature_Map::write_mapping(const char* fname) const
{
  IWString_and_File_Descriptor output;

  if (!output.open(fname)) {
    cerr << "GFP_to_Feature_Map::write_mapping:cannot open '" << fname << "'\n";
    return 0;
  }

  return write_mapping(output);
}

// After each molecule is processed, we might need to flush
// `output`, or write if it gets too full.
void
MaybeFlush(IWString_and_File_Descriptor& output)
{
  if (flush_after_every_molecule) {
    output.flush();
  } else {
    output.write_if_buffer_holds_more_than(32768);
  }
}

static int
_write_mapping(const Bit_and_Integer& bai, IWString_and_File_Descriptor& output)
{
  // for (Bit_and_Integer::const_iterator i = bai.begin(); i != bai.end(); ++i)
  for (auto i = bai.cbegin(); i != bai.cend(); ++i) {
    int v = (*i).second;

    if (v <= 0) {
      continue;
    }

    output << (*i).first << ' ' << v << '\n';
  }

  output << "|\n";

  MaybeFlush(output);

  return 1;
}

int
GFP_to_Feature_Map::write_mapping(IWString_and_File_Descriptor& output) const
{
  output << HEADER_RECORD << '\n';
  output << '\n';

  if (_properties.size()) {
    output << PROPERTIES_TAG << '\n';
  }
  if (_number_fixed) {
    output << COUNT_FIXED << ' ' << _number_fixed << '\n';
  }
  if (_number_sparse) {
    output << COUNT_SPARSE << ' ' << _number_sparse << '\n';
  }

  output << "|\n";

  if (_properties.size()) {
    output << PROPERTIES_TAG << ' ' << _nproperties << '\n';
    _write_mapping(_properties, output);
  }

  for (int i = 0; i < _number_fixed; i++) {
    output << FIXED_TAG << ' ' << _nbits_fixed[i] << '\n';
    _write_mapping(_bai_fixed[i], output);
  }

  for (int i = 0; i < _number_sparse; i++) {
    output << SPARSE_TAG << '\n';

    _write_mapping(_bai_sparse[i], output);
  }

  output.flush();

  return output.good();
}

int
GFP_to_Feature_Map::read_existing_mapping(const char* fname)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "GFP_to_Feature_Map::read_existing_mapping:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_existing_mapping(input);
}

static int
read_mapping_record(const const_IWSubstring& buffer, Bit_and_Integer& bai)
{
  int i = 0;
  const_IWSubstring token;

  if (!buffer.nextword(token, i)) {  // skip blank lines
    return 1;
  }

  unsigned int b;
  if (!token.numeric_value(b)) {
    cerr << "Invalid bit number '" << buffer << "'\n";
    return 0;
  }

  if (!buffer.nextword(token, i)) {
    cerr << "Mapping records must contain at least two tokens '" << buffer << "'\n";
    return 0;
  }

  int f;

  if (!token.numeric_value(f) || f < 1) {
    cerr << "Feature numbers must be while +ve numbers '" << buffer << "'\n";
    return 0;
  }

  Bit_and_Integer::const_iterator c = bai.find(b);
  if (c != bai.end()) {
    cerr << "Duplicate bit number " << b << endl;
    return 0;
  }

  bai[b] = f;

  return 1;
}

static int
read_mapping(iwstring_data_source& input, Bit_and_Integer& bai)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with('#') || 0 == buffer.length()) {
      continue;
    }

    if ('|' == buffer) {
      break;
    }

    if (!read_mapping_record(buffer, bai)) {
      cerr << "Invalid record in bit->feature map '" << buffer << "', line "
           << input.lines_read() << endl;
      return 0;
    }
  }

  return 1;
}

/*static int
read_mapping (iwstring_data_source & input,
              Bit_and_Integer & bai,
              int & ndx,
              int max_ndx)
{
  if (ndx >= max_ndx)
  {
    cerr << "Too many fingerprints of given type " << ndx << ", max is " << max_ndx <<
endl; return 0;
  }

  int rc = read_mapping (input, bai);

  ndx++;

  return rc;
}*/

int
GFP_to_Feature_Map::read_existing_mapping(iwstring_data_source& input)
{
  const_IWSubstring buffer;

  if (!input.next_record(buffer)) {
    cerr << "GFP_to_Feature_Map::read_existing_mapping:cannot read header\n";
    return 0;
  }

  buffer.strip_trailing_blanks();

  if (HEADER_RECORD != buffer) {
    cerr << "GFP_to_Feature_Map::read_existing_mapping:header record mismatch\n";
    cerr << "Expected '" << HEADER_RECORD << "' got '" << buffer << "'\n";
    return 0;
  }

  int properties_present = 0;
  _number_fixed = 0;
  _number_sparse = 0;

  while (input.next_record(buffer)) {
    if (0 == buffer.length() || buffer.starts_with('#')) {
      continue;
    }

    if (buffer.starts_with(COUNT_FIXED)) {
      buffer.remove_leading_words(1);
      if (!buffer.numeric_value(_number_fixed) || _number_fixed < 1) {
        cerr << "The number of fixed fingerprints must be a whole +ve number '" << buffer
             << "'\n";
        return 0;
      }
    } else if (buffer.starts_with(COUNT_SPARSE)) {
      buffer.remove_leading_words(1);
      if (!buffer.numeric_value(_number_sparse) || _number_sparse < 1) {
        cerr << "The number of sparse fingerprints must be a whole +ve number '" << buffer
             << "'\n";
        return 0;
      }
    } else if (PROPERTIES_TAG == buffer) {
      properties_present = 1;
    } else if ('|' == buffer) {  // end of preamble
      break;
    } else {
      cerr << "Unrecognised record in bit->feature cross reference '" << buffer << "'\n";
      return 0;
    }
  }

  // Did we get anything?

  if (properties_present) {
    ;
  } else if (_number_fixed) {
    ;
  } else if (_number_sparse) {
    ;
  } else {
    cerr << "No bits present in preamble\n";
    return 0;
  }

  if (_number_fixed > 0) {
    _bai_fixed = new Bit_and_Integer[_number_fixed];
  }

  if (_number_sparse > 0) {
    _bai_sparse = new Bit_and_Integer[_number_sparse];
  }

  int ndx_fixed = 0;
  int ndx_sparse = 0;
  int ndx_mpr = 0;
  while (input.next_record(buffer)) {
    int rc = 0;

    if (buffer.starts_with(PROPERTIES_TAG)) {
      if (ndx_mpr > 0) {
        cerr << "GFP_to_Feature_Map::read_existing_mapping:only one set of properties "
                "allowed, line "
             << input.lines_read() << endl;
        return 0;
      }
      rc = read_mapping(input, _properties);

      ndx_mpr++;
    } else if (buffer.starts_with(FIXED_TAG))
    //  else if (FIXED_TAG == buffer)
    {
      if (ndx_fixed >= _number_fixed) {
        cerr << "GFP_to_Feature_Map::read_existing_mapping:too many fixed fingerprints\n";
        return 0;
      }

      rc = read_mapping(input, _bai_fixed[ndx_fixed]);
      ndx_fixed++;
    } else if (SPARSE_TAG == buffer) {
      if (ndx_sparse >= _number_sparse) {
        cerr
            << "GFP_to_Feature_Map::read_existing_mapping:too many sparse fingerprints\n";
        return 0;
      }

      rc = read_mapping(input, _bai_sparse[ndx_sparse]);
      ndx_sparse++;
    } else {
      cerr << "Unrecognised group heading in bit xref file '" << buffer << "'\n";
      return 0;
    }

    if (0 == rc) {
      return 0;
    }
  }

  if (properties_present && 0 == ndx_mpr) {
    cerr << "GFP_to_Feature_Map::read_existing_mapping:no properties\n";
    return 0;
  }

  if (ndx_fixed != _number_fixed) {
    cerr << "GFP_to_Feature_Map::read_existing_mapping:missing fixed fingerprints, got "
         << (ndx_fixed + 1) << ", expected " << _number_fixed << "\n";
    return 0;
  }

  if (ndx_sparse != _number_sparse) {
    cerr << "GFP_to_Feature_Map::read_existing_mapping:missing sparse fingerprints, got "
         << (ndx_sparse + 1) << ", expected " << _number_sparse << "\n";
    return 0;
  }

  if (verbose) {
    cerr << "Read mapping";
    if (_number_fixed) {
      cerr << ' ' << _number_fixed << " fixed";
    }
    if (_number_sparse) {
      cerr << ' ' << _number_sparse << " sparse";
    }
    cerr << endl;
  }

  return 1;
}

/*
  Now bit Bit_and_Integer objects are holding bit and feature numbers
*/

class Feature_Number_Comparator
{
 public:
  int operator()(const Feature_Value_Pair&, const Feature_Value_Pair&) const;
};

int
Feature_Number_Comparator::operator()(const Feature_Value_Pair& p1,
                                      const Feature_Value_Pair& p2) const
{
  int f1 = p1.feature_number();
  int f2 = p2.feature_number();

  if (f1 < f2) {
    return -1;
  } else if (f1 > f2) {
    return 1;
  } else {
    return 0;  // should never happen
  }
}

/*
  One time we just need to check existence of the key, another time we
  need to return the associated data
*/

static int
find_identifier_in_activity_hash(const IW_STL_Hash_Map_String& activity,
                                 const IWString& id, IWString* activity_for_id = nullptr)
{
  IW_STL_Hash_Map_String::const_iterator f = activity.find(id);

  if (f != activity.end()) {
    if (nullptr != activity_for_id) {
      *activity_for_id = (*f).second;
    }

    return 1;
  }

  if (id.nwords() > 1) {
    IWString tmp(id);
    tmp.truncate_at_first(' ');
    return find_identifier_in_activity_hash(activity, tmp, activity_for_id);
  }

  if (!strip_leading_zeros) {
    return 0;
  }

  if (!id.starts_with('0')) {
    return 0;
  }

  IWString tmp(id);

  tmp.remove_leading_chars('0');

  f = activity.find(tmp);

  if (f == activity.end()) {
    return 0;
  }

  if (nullptr != activity_for_id) {
    *activity_for_id = (*f).second;
  }

  return 1;
}

static int
handle_descriptor_data(Descriptor_File& dfile, const IWString& id,
                       const IWString& activity_this_molecule,
                       IWString_and_File_Descriptor& output)
{
  const_IWSubstring zdata;

  if (!dfile.fetch_data_for_id(id, zdata)) {
    cerr << "handle_descriptor_data:no descriptor data for '" << id << "'\n";
    return 0;
  }

  output << id << ' ' << activity_this_molecule;

  return write_descriptor_data(zdata, output);
}

static int
do_produce_normalised_output(const Feature_Value_Pair* fvp, int nfeatures,
                             IWString_and_File_Descriptor& output)
{
  int sum2 = 0;

  for (int i = 0; i < nfeatures; i++) {
    int vi = fvp[i].zvalue();

    sum2 += (vi * vi);
  }

  double norm = sqrt(static_cast<double>(sum2));

  for (int i = 0; i < nfeatures; i++) {
    const Feature_Value_Pair& fvpi = fvp[i];

    double v = static_cast<double>(fvp[i].zvalue()) / norm;

    bit_number_iwdigits.append_number(output, fvpi.feature_number());
    output << static_cast<float>(v);
  }

  output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
count_features(const Feature_Value_Pair* fvp, const int nfeatures)
{
  int rc = 0;

  for (int i = 0; i < nfeatures; ++i) {
    rc += fvp[i].zvalue();
  }

  return rc;
}

static int
gfp_to_svm_lite(const IW_General_Fingerprint& fp, GFP_to_Feature_Map& gfp2fm,
                Feature_Value_Pair* fvp, IWString_and_File_Descriptor& output)
{
  const IWString& id = fp.id();

  IWString activity_this_molecule;
  if (write_random_activity_values) {
    static std::random_device rd;
    std::uniform_real_distribution u(0.0, 1.0);
    activity_this_molecule << u(rd);
  } else if (0 == activity.size()) {
    activity_this_molecule = "0";
  } else {
    if (!find_identifier_in_activity_hash(activity, id, &activity_this_molecule)) {
      cerr << "No activity data for '" << id << "'\n";
      if (ignore_missing_activity) {
        return 1;
      }

      return 0;
    }

    //  cerr << "Activity for '" << id << "' is '" << activity_this_molecule << "'\n";
  }

  int bits_not_present_in_training_set = 0;

  int nfeatures = gfp2fm.gather_bits_set(fp, fvp, bits_not_present_in_training_set);

  if (verbose > 2) {
    cerr << fp.id() << " nbits " << nfeatures << ", not in training set "
         << bits_not_present_in_training_set << '\n';
  }

  //  nfeatures is the number of bits that were set.
  //  bits_not_present_in_training_set is number of fp's in initial fp not present
  //  therefore the sum of them is the total number of bits in the initial fingerprint(s)

  if (bnpts.active()) {
    int bits_in_initial_fps = bits_not_present_in_training_set + nfeatures;

    float fraction_missing = static_cast<float>(bits_not_present_in_training_set) /
                             static_cast<float>(bits_in_initial_fps);

    bnpts << fp.id() << ' ' << bits_not_present_in_training_set << ' '
          << bits_in_initial_fps << ' ' << fraction_missing << '\n';

    bnpts.write_if_buffer_holds_more_than(32768);
  }

  Feature_Number_Comparator fnc;

  iwqsort(fvp, nfeatures, fnc);

  if (nullptr != dfile) {
    if (handle_descriptor_data(*dfile, id, activity_this_molecule, output)) {
      ;
    } else if (ignore_missing_descriptor_data) {
      if (verbose > 1) {
        cerr << "Ignoring missing descriptor data for '" << id << "'\n";
      }
      missing_descriptor_data_ignored++;
    } else {
      cerr << "Cannot write descriptor data for '" << id << "'\n";
      return 0;
    }
  } else if (intel_csr.active()) {
    if (stream_for_activity_values.is_open()) {
      stream_for_activity_values << activity_this_molecule << '\n';
      stream_for_activity_values.write_if_buffer_holds_more_than(4096);
    }
    return intel_csr.extra(fvp, nfeatures);
  } else {
    output << activity_this_molecule;
  }

  if (first_features_are_nbits_and_weight) {
    bit_number_iwdigits.append_number(output, 1);
    count_iwdigits.append_number(output, nfeatures);

    bit_number_iwdigits.append_number(output, 2);
    const int nset = count_features(fvp, nfeatures);
    count_iwdigits.append_number(output, nset);
  }

  if (produce_normalised_output) {
    do_produce_normalised_output(fvp, nfeatures, output);
  } else {
    int istart = 0;

    for (int i = istart; i < nfeatures; i++) {
      const Feature_Value_Pair& fvpi = fvp[i];

      //    output << ' ' << b.column() << ':' << b.count();
      bit_number_iwdigits.append_number(output, fvpi.feature_number());
      //    output << fvpi.zvalue();
      count_iwdigits.append_number(
          output, fvpi.zvalue() - processing_difference_fingerprint_offset);
    }
  }

  if (append_identifier_at_end_of_each_record) {
    output << " # " << id;
  }

  output << "\n";

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
gfp_to_svm_lite(IW_TDT& tdt, GFP_to_Feature_Map& gfp2fm, Feature_Value_Pair* fvp,
                IWString_and_File_Descriptor& output)
{
  IW_General_Fingerprint fp;

  int fatal;

  if (!fp.construct_from_tdt(tdt, fatal)) {
    cerr << "Cannot build fingerprint\n";
    return 0;
  }

  if (fp.id().nwords() > 1) {
    IWString tmp = fp.id();

    if (verbose > 1) {
      cerr << "ID '" << tmp << "' being truncated to first token\n";
    }

    tmp.truncate_at_first(' ');
    fp.set_id(tmp);
    identifiers_truncated_to_first_token++;
  }

  return gfp_to_svm_lite(fp, gfp2fm, fvp, output);
}

static int
gfp_to_svm_lite(iwstring_data_source& input, GFP_to_Feature_Map& gfp2fm,
                Feature_Value_Pair* fvp, IWString_and_File_Descriptor& output)
{
  IW_TDT tdt;
  while (tdt.next(input)) {
    tdts_read++;

    if (!gfp_to_svm_lite(tdt, gfp2fm, fvp, output)) {
      cerr << "Fatal error processing '" << tdt << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
gfp_to_svm_lite(const char* fname, GFP_to_Feature_Map& gfp2fm, Feature_Value_Pair* fvp,
                IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return gfp_to_svm_lite(input, gfp2fm, fvp, output);
}

static int
extract_activity_from_identifier_token(const IW_General_Fingerprint& fp,
                                       int activity_column)
{
  const_IWSubstring id = fp.id();

  const_IWSubstring token;

  if (!id.word(activity_column, token)) {
    cerr << "Cannot extract column " << (activity_column + 1) << " from '" << id << "'\n";
    return 0;
  }

  float a;

  if (!token.numeric_value(a)) {
    cerr << "INvalid activity value '" << token << "'\n";
    return 0;
  }

  id.truncate_at_first(' ');

  activity[id] = token;

  // cerr << "Set activity of '" << id << "' to " << a << endl;

  return 1;
}

static int
profile_bits_set(IW_General_Fingerprint& fp, GFP_to_Feature_Map& gfp2fm)
{
  const IWString& id = fp.id();

  // cerr << "Examining '" << id << "', activity_data_read_from_file " <<
  // activity_data_read_from_file << endl;

  if (activity_data_read_from_file) {
    if (!find_identifier_in_activity_hash(activity, id)) {
      cerr << "No activity data for '" << id << "'\n";
      if (ignore_missing_activity) {
        return 1;
      }

      return 0;
    }
  } else if (write_random_activity_values) {
  } else if (!extract_activity_from_identifier_token(fp, activity_column)) {
    return 0;
  }

  if (id.nwords() > 1) {
    const_IWSubstring tmp(id);

    tmp.truncate_at_first(' ');
    fp.set_id(tmp);
  }

  return gfp2fm.profile_bits_set(fp);
}

static int
profile_bits_set(IW_TDT& tdt, GFP_to_Feature_Map& gfp2fm, int& fingerprints_in_input)
{
  IW_General_Fingerprint fp;

  int fatal;

  if (!fp.construct_from_tdt(tdt, fatal)) {
    cerr << "Cannot build fingerprint\n";
    return 0;
  }

  fingerprints_in_input++;

  return profile_bits_set(fp, gfp2fm);
}

static int
profile_bits_set(iwstring_data_source& input, GFP_to_Feature_Map& gfp2fm,
                 int& fingerprints_in_input)
{
  IW_TDT tdt;
  while (tdt.next(input)) {
    if (!profile_bits_set(tdt, gfp2fm, fingerprints_in_input)) {
      cerr << "Fatal error processing '" << tdt << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
profile_bits_set(const char* fname, GFP_to_Feature_Map& gfp2fm,
                 int& fingerprints_in_input)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return profile_bits_set(input, gfp2fm, fingerprints_in_input);
}

static int
read_activity_data_record(const const_IWSubstring& buffer,
                          IW_STL_Hash_Map_String& activity)
{
  IWString id, token;
  int i = 0;

  char separator = ' ';

  if (buffer.nextword(id, i) && buffer.nextword(token, i)) {
    ;
  } else if (buffer.nextword(id, i, '\t') && buffer.nextword(token, i, '\t')) {
    separator = '\t';
  } else if (buffer.nextword(id, i, ',') && buffer.nextword(token, i, ',')) {
    separator = ',';
  } else {
    cerr << "Cannot separate into identifier and activity '" << buffer << "'\n";
    return 0;
  }

  if (1 != activity_column) {
    if (!buffer.word(activity_column, token, separator)) {
      cerr << "Cannot extract column '" << (activity_column + 1) << " from record\n";
      return 0;
    }
  }

  double a;
  if (token.numeric_value(a)) {
    ;
  } else if (0 == activity.size()) {  // header record
    ;
  } else {
    if (warn_non_numeric_activity) {
      cerr << "Warning, non numeric activity value '" << token << "', id '" << id
           << "'\n";
    }
    non_numeric_activity_encountered++;
  }

  if (strip_leading_zeros) {
    id.remove_leading_chars('0');
  }

  activity[id] = token;

  // cerr << "for id '" << id << "' value '" << activity[id] << "', token '" << token <<
  // "'\n";

  return 1;
}

#ifdef TABS_NOT_ALLOWED
static int
read_activity_data_record(const const_IWSubstring& buffer,
                          IW_STL_Hash_Map_String& activity)
{
  int i = 0;
  IWString id;

  if (!buffer.nextword(id, i)) {
    cerr << "Cannot extract identifier\n";
    return 0;
  }

  IWString token;

  if (!buffer.nextword(token, i)) {
    cerr << "Not enough tokens on experimental data record\n";
    return 0;
  }

  if (1 != activity_column) {
    if (!buffer.word(activity_column, token)) {
      cerr << "Cannot extract column '" << (activity_column + 1) << " from record\n";
      return 0;
    }
  }

  double a;
  if (token.numeric_value(a)) {
    ;
  } else if (0 == activity.size()) {  // header record
    ;
  } else {
    if (warn_non_numeric_activity) {
      cerr << "Warning, non numeric activity value '" << token << "', id '" << id
           << "'\n";
    }
    non_numeric_activity_encountered++;
  }

  if (strip_leading_zeros) {
    id.remove_leading_chars('0');
  }

  activity[id] = token;

  cerr << "for id '" << id << "' value '" << activity[id] << "', token '" << token
       << "'\n";

  return 1;
}
#endif

static int
read_activity_data(iwstring_data_source& input, IW_STL_Hash_Map_String& activity)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with('#')) {
      continue;
    }

    if (!read_activity_data_record(buffer, activity)) {
      cerr << "Cannot read activity data, line " << input.lines_read() << endl;
      cerr << "'" << buffer << "'\n";
      return 0;
    }
  }

  return activity.size();
}

static int
read_activity_data(const const_IWSubstring& fname, IW_STL_Hash_Map_String& activity)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);
  input.set_translate_tabs(1);

  return read_activity_data(input, activity);
}

static int
ReadActivityData(const const_IWSubstring& a) {
  if (a == "RANDOM") {
    write_random_activity_values = 1;
    if (verbose) {
      cerr << "Random activity values will be written";
    }

    set_default_iwstring_double_concatenation_precision(4);

    return 1;
  }

  if (!read_activity_data(a, activity)) {
    cerr << "Cannot read activity data from '" << a << "'\n";
    return 5;
  }

  if (verbose) {
    cerr << "Read " << activity.size() << " activity values from '" << a << "'\n";
  }

  activity_data_read_from_file = 1;

  return 1;
}

static int
count_tdts(const char* fname)
{
  IWString_Data_Source_MMAP input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return input.count_records_starting_with("|");
}

static void
DisplayDashYOptions(std::ostream& output)
{
  output << " -Y flush            flush output after each molecule\n";

  ::exit(0);
}

static int
gfp_to_svm_lite(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:c:F:P:Q:p:C:U:m:VN:zD:dguI:0tX:qY:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (need_to_call_initialise_fingerprints(cl)) {
    if (!initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise fingerprint option(s)\n";
      usage(4);
    }
  }

  if (cl.option_present('z')) {
    strip_leading_zeros = 1;

    if (verbose) {
      cerr << "Will strip leading zero's from identifiers\n";
    }
  }

  if (cl.option_present('V')) {
    append_identifier_at_end_of_each_record = 0;

    if (verbose) {
      cerr << "No identifiers appended (SVMLITE)\n";
    }
  }

  if (cl.option_present('I')) {  // activity file not needed with Intel
    ;
  } else if (cl.option_present('C') && !cl.option_present('A') &&
             !cl.option_present('c')) {
    cerr << "Must specify activity file via the -A option or -c option\n";
    usage(4);
  }

  if (!cl.option_present('C') && !cl.option_present('U') && !cl.option_present('A')) {
    cerr << "Must specify activity file via the -A option\n";
    usage(5);
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', activity_column) || activity_column < 1) {
      cerr << "The activity column (-c) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Activity data in column " << activity_column << " of activity file\n";
    }

    activity_column--;
  }

  if (cl.option_present('u')) {
    produce_normalised_output = 1;

    if (verbose) {
      cerr << "Output will be normalised unit vectors\n";
    }
  }

  if (cl.option_present('X')) {
    if (!cl.value('X', processing_difference_fingerprint_offset) ||
        processing_difference_fingerprint_offset < 0) {
      cerr
          << "The difference fingerprint offset (-X) option must be a whole +ve number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Processing a difference fingerprint, offset "
           << processing_difference_fingerprint_offset << endl;
    }
  }

  if (cl.option_present('t')) {
    first_features_are_nbits_and_weight = 1;

    if (verbose) {
      cerr << "First feature is nset\n";
    }
  }

  if (cl.option_present('q')) {
    warn_non_numeric_activity = 0;

    if (verbose) {
      cerr << "Will not warn about non numeric activity values\n";
    }
  }

  if (cl.option_present('A')) {
    const_IWSubstring a = cl.string_value('A');
    if (! ReadActivityData(a)) {
      cerr << "Cannot read activity data '" << a << "'\n";
      return 1;
    }
  }

  if (cl.option_present('D')) {
    const char* fname = cl.option_value('D');

    dfile = new Descriptor_File();

    if (!dfile->build(fname)) {
      cerr << "Cannot initialise descriptor file '" << fname << "'\n";
      return 4;
    }

    if (cl.option_present('d')) {
      write_descriptor_data_as_descriptors = 1;
      append_identifier_at_end_of_each_record = 0;

      if (verbose) {
        cerr << "Descriptors written in descriptor file form\n";
      }
    }

    if (cl.option_present('g')) {
      ignore_missing_descriptor_data = 1;

      if (verbose) {
        cerr << "Will ignore missing descriptor data\n";
      }
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', max_count) || max_count < 1 || max_count > 255) {
      cerr << "The max count must be a +ve whole number between 1 and 255\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Counted fingerprints truncated to max count " << max_count << endl;
    }
  }

  set_sparsefp_warn_empty_data(0);

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('U') && cl.option_present('C')) {
    cerr << "Cannot specify both -U and -C options\n";
    usage(3);
  }

  if (cl.option_present('U') && cl.option_present('p')) {
    cerr << "Cannot use the -p option with the -U option\n";
    usage(4);
  }

  if (cl.option_present('U') && cl.option_present('A')) {
    cerr << "Both -U and -A options specified, very unusual\n";
  }

  // This hash plays multiple roles. Ultimately it will be turned into a mapping
  // between bit number and (numbered from 1) feature numbers.

  GFP_to_Feature_Map gfp2fm;

  int fingerprints_in_input = 0;

  if (cl.option_present('U'))  // use an existing mapping
  {
    const char* u = cl.option_value('U');

    if (!gfp2fm.read_existing_mapping(u)) {
      cerr << "Cannot read existing mapping '" << u << "'\n";
      return 4;
    }

    for (int i = 0; i < cl.number_elements(); ++i) {
      fingerprints_in_input += count_tdts(cl[i]);
    }
  } else  // establish the mapping from the current set of file(s)
  {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!profile_bits_set(cl[i], gfp2fm, fingerprints_in_input)) {
        cerr << "Cannot profile bits set in '" << cl[i] << "'\n";
        return 3;
      }
    }

    if (0 == gfp2fm.number_features_defined()) {
      cerr << "Yipes, no features identified in input file(s), cannot continue\n";
      exit(2);
    }

    if (nullptr != dfile) {
      gfp2fm.extract_information_from_descriptor_file(*dfile);
    }

    if (cl.option_present('p')) {
      int support;
      if (!cl.value('p', support) || support < 1) {
        cerr << "The minimum support level (-p) must be a whole +ve number\n";
        usage(4);
      }

      int bits_suppressed_by_support_requirement = 0;

      //    int initial_size = gfp2fm.number_features_defined();

      int final_size =
          gfp2fm.assign_offsets(support, bits_suppressed_by_support_requirement);

      if (verbose) {
        cerr << bits_suppressed_by_support_requirement
             << " bits suppressed by support requirement " << support << endl;
      }

      //    if (initial_size != final_size)
      //      cerr << "Hash size mismatch " << initial_size << " vs " << final_size <<
      //      endl;

      if (0 == final_size) {
        cerr << "No bits selected for output, try decreasing support level\n";
        return 6;
      }
    } else {
      gfp2fm.assign_offsets();
    }

    if (cl.option_present('C')) {
      const char* x = cl.option_value('C');

      if (!gfp2fm.write_mapping(x)) {
        cerr << "Cannot write bit-to-feature translation to '" << x << "'\n";
        return 6;
      }
    }

    if (verbose > 2) {
      gfp2fm.debug_print(cerr);
    }
  }

  if (activity.size() && cl.option_present('N')) {
    cerr << "The -N option doesn't make sense with the -A option\n";
    usage(3);
  }

  if (cl.option_present('N')) {
    const char* n = cl.option_value('N');

    if (!bnpts.open(n)) {
      cerr << "Cannot open stream for bits not present in training set '" << n << "'\n";
      return 4;
    }

    if (verbose) {
      cerr << "The number of bits not present in the taining set written to '" << n
           << "'\n";
    }

    set_default_iwstring_float_concatenation_precision(2);
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "flush") {
        flush_after_every_molecule = 1;
        if (verbose) {
          cerr << "Will flush output after every molecule\n";
        }
      } else if (y == "help") {
        DisplayDashYOptions(cerr);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  // We need to know the largest size of an output record. When reading a
  // training set file, we must scan the whole dataset first. When reading
  // a test set, the largest output record is limited to the number of
  // features in the cross reference file

  int h;

  if (cl.option_present('U')) {
    h = gfp2fm.number_features_defined();
    if (verbose) {
      cerr << "Mapping contains " << h << " features\n";
    }
  } else {
    h = gfp2fm.highest_number_features();

    if (verbose) {
      cerr << "Input records have as many as " << h << " features\n";
    }
  }

  if (h < 1) {
    cerr << "No features in input\n";
    return 0;
  }

  int s = h;
  if (first_features_are_nbits_and_weight) {
    s += 2;
  }

  Feature_Value_Pair* fvp = new Feature_Value_Pair[s];
  std::unique_ptr<Feature_Value_Pair[]> free_fvp(fvp);

  if (cl.option_present('I')) {
    const char* fname = cl.option_value('I');

    if (!stream_for_activity_values.open(fname)) {
      cerr << "Cannot open stream for ordered activity values '" << fname << "'\n";
      return 1;
    }

    intel_csr.activate(fingerprints_in_input);
    cerr << "active " << intel_csr.active() << endl;

    bit_number_iwdigits.set_leading_string(",");
    append_identifier_at_end_of_each_record = 0;

    if (cl.option_present('0')) {
      intel_csr.set_index_base(0);
    }
  } else {
    bit_number_iwdigits.set_include_leading_space(1);
  }

  bit_number_iwdigits.initialise(10000);

  if (!intel_csr.active()) {
    bit_number_iwdigits.append_to_each_stored_string(":");
  }

  count_iwdigits.initialise(256);

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!gfp_to_svm_lite(cl[i], gfp2fm, fvp, output)) {
      rc = i + 1;
      break;
    }
  }

  if (intel_csr.active()) {
    intel_csr.do_write(gfp2fm.number_features_defined(), output);
  }

  output.flush();

  if (identifiers_truncated_to_first_token) {
    cerr << identifiers_truncated_to_first_token
         << " identifiers truncated to first token\n";
  }

  if (verbose) {
    cerr << "Read " << tdts_read << " tdt's\n";
    if (highest_sparse_bit_number > 0) {
      cerr << "Sparse fingerprint bits numbers between " << lowest_sparse_bit_number
           << " and " << highest_sparse_bit_number << endl;
      if (acc_count.n() > 1) {
        cerr << "Sparse bit counts between " << acc_count.minval() << " and "
             << acc_count.maxval() << " ave " << static_cast<float>(acc_count.average())
             << endl;
      }
    }

    if (ignore_missing_descriptor_data) {
      cerr << missing_descriptor_data_ignored << " items with missing descriptor data\n";
    }
  }

  if (non_numeric_activity_encountered) {
    cerr << non_numeric_activity_encountered
         << " non numeric activity values encountered\n";
  }

  // Destroying the STL hashes is expensive, so take a quick exit

  cerr.flush();

  //_exit(rc);

  if (nullptr != dfile) {  // another expensive thing to be deleted
    delete dfile;
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = gfp_to_svm_lite(argc, argv);

  return rc;
}
