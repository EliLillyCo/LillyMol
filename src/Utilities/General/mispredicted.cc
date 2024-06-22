/*
  We have generated some predictions, which are the molecules
  that are most frequently mipredicted
*/

#include <dirent.h>
#include <stdlib.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <sys/types.h>

#define IWQSORT_FO_IMPLEMENTATION 1

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;

const char* prog_name = NULL;

static int verbose = 0;

static int predicted_column = 1;
static int experimental_column = 1;

static int ignore_missing_experimental_data = 0;

static int predicted_values_read = 0;

static int use_cutoff_to_convert_to_class = 0;
static float cutoff = static_cast<float>(0.0);

static int fault_tolerant_mode = 0;

static int take_absolute_differences = 1;

static int write_worst_predicted = 0;

static float write_differences_threshold = -1.0f;

static int produce_descriptor_file = 0;

static int balance_numeric_filtering = 0;

static IWString_and_File_Descriptor stream_for_summary;

static int data_includes_distance_to_support_vector = 0;

static int sort_best_to_worst = 1;

static IW_STL_Hash_Map_String id_to_smiles;

static int ignore_duplicate_experimental_values = 0;
static uint32_t duplicate_experimental_values_ignored = 0;

class Binning_Data
{
 private:
  std::vector<float> _bin_start;
  std::vector<IWString> _bin_text;

 public:
  int
  active() const
  {
    return _bin_start.size();
  }

  int
  build(const char*);
  int
  build(iwstring_data_source&);

  int
  bucket(float) const;
};

int
Binning_Data::build(const char* fname)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Binning_Data::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input);
}

int
Binning_Data::build(iwstring_data_source& input)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (buffer.starts_with('#')) {
      continue;
    }

    if (0 == buffer.length()) {
      continue;
    }

    const_IWSubstring token;
    int i = 0;

    float v;
    if (!buffer.nextword(token, i) || !token.numeric_value(v) ||
        !buffer.nextword(token, i)) {
      cerr << "Binning_Data::build:cannot extract distance or invalid distance '"
           << buffer << "'\n";
      return 0;
    }

    _bin_start.push_back(v);
    _bin_text.push_back(token);
  }

  if (0 == _bin_start.size()) {
    cerr << "Binning_Data::build:no data\n";
    return 0;
  }

  return 1;
}

int
Binning_Data::bucket(float v) const
{
  for (unsigned int i = 0; i < _bin_start.size(); ++i) {
    if (v <= _bin_start[i]) {
      return i;
    }
  }

  return _bin_start.size();
}

class Entity
{
 protected:
  const IWString _id;
  int _times_predicted;

 public:
  Entity(const IWString&);

  const IWString&
  id() const
  {
    return _id;
  }

  int
  times_predicted() const
  {
    return _times_predicted;
  }

  //  virtual int append_details(IWString_and_File_Descriptor &) const;
};

Entity::Entity(const IWString& s) : _id(s)
{
  _times_predicted = 0;

  return;
}

class Entity_Classification : public Entity
{
 private:
  IWString _correct_class;

  int _correct_predictions;
  int _incorrect_predictions;

 public:
  Entity_Classification(const IWString&);

  void
  set_class(const const_IWSubstring& s)
  {
    _correct_class = s;
  }

  const IWString&
  correct_class() const
  {
    return _correct_class;
  }

  void
  predicted_as(const const_IWSubstring& s);

  float accuracy() const;

  int report(IWString_and_File_Descriptor&) const;

  int write_summary_data(IWString_and_File_Descriptor&) const;
};

Entity_Classification::Entity_Classification(const IWString& s) : Entity(s)
{
  _correct_predictions = 0;
  _incorrect_predictions = 0;
}

void
Entity_Classification::predicted_as(const const_IWSubstring& s)
{
  _times_predicted++;

  if (_correct_class == s) {
    _correct_predictions++;
  } else {
    _incorrect_predictions++;
  }

  return;
}

float
Entity_Classification::accuracy() const
{
  if (0 == _times_predicted) {
    return static_cast<float>(0.0);
  }

  return static_cast<float>(_correct_predictions) / static_cast<float>(_times_predicted);
}

int
Entity_Classification::report(IWString_and_File_Descriptor& output) const
{
  static constexpr char kSep = ' ';

  output << _id << kSep;

  if (produce_descriptor_file) {
    output << _correct_class << kSep << _times_predicted << kSep <<
              _correct_predictions << kSep <<
              iwmisc::Fraction<float>(_correct_predictions, _times_predicted);
  } else {
    output << _correct_class << " predicted " << _times_predicted;
    if (_times_predicted > 0) {
      output << " correct " << _correct_predictions << kSep << accuracy();
    }
  }

  output << '\n';

  return 1;
}

int
Entity_Classification::write_summary_data(IWString_and_File_Descriptor& output) const
{
  output << _id << ' '
         << static_cast<float>(_correct_predictions) /
                static_cast<float>(_times_predicted)
         << " . " << _times_predicted << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

template <typename T>
class Entity_Continuous : public Entity
{
 private:
  T _experimental_value;
  Accumulator<T> _predictions;
  Accumulator<T> _diffs;

  std::vector<T> _diff;
  std::vector<float> _dist;

 public:
  Entity_Continuous(const IWString&);

  void
  set_experimental_value(T s) {
    _experimental_value = s;
  }

  T
  experimental_value() const {
    return _experimental_value;
  }

  int predicted_as(const const_IWSubstring&);
  void predicted_as(T);
  int prediction_and_distance(const const_IWSubstring&, const const_IWSubstring&);

  T accuracy() const;

  int report(IWString_and_File_Descriptor&) const;

  int write_summary_data(IWString_and_File_Descriptor&) const;
};

template <typename T>
Entity_Continuous<T>::Entity_Continuous(const IWString& s) : Entity(s)
{
  _experimental_value = static_cast<T>(-99);
}

template <typename T>
void
Entity_Continuous<T>::predicted_as(T p)
{
  _times_predicted++;

  _predictions.extra(p);

  if (take_absolute_differences) {
    _diffs.extra(fabs(_experimental_value - p));
  } else {
    _diffs.extra(_experimental_value - p);
  }

  return;
}

template <typename T>
int
Entity_Continuous<T>::predicted_as(const const_IWSubstring& p)
{
  T v;
  if (!p.numeric_value(v)) {
    cerr << "Entity_Continuous::predicted_as:invalid numeric '" << p << "'\n";
    return 0;
  }

  predicted_as(v);

  return 1;
}

template <typename T>
int
Entity_Continuous<T>::prediction_and_distance(const const_IWSubstring& sp,
                                              const const_IWSubstring& sd)
{
  T p;
  float d;
  if (!sp.numeric_value(p) || !sd.numeric_value(d) || d < 0.0f || d > 1.0f) {
    cerr << "Entity_Continuous::prediction_and_distance:invalid prediction '" << sp
         << "' and/or distance '" << sd << "'\n";
    return 0;
  }

  predicted_as(p);
  _diff.push_back(fabs(_experimental_value - p));
  _dist.push_back(d);

  return 1;
}

template <typename T>
T
Entity_Continuous<T>::accuracy() const
{
  if (0 == _times_predicted) {
    return static_cast<T>(0);
  }

  return static_cast<T>(_diffs.average_if_available_minval_if_not());
}

template <typename T>
int
Entity_Continuous<T>::report(IWString_and_File_Descriptor& output) const
{
  static constexpr char kSep = ' ';

  if (! id_to_smiles.empty()) {
    const auto iter = id_to_smiles.find(_id);
    if (iter == id_to_smiles.end()) {
      cerr << "NO smiles for '" << _id << "'\n";
      output << '*' << kSep;
    } else {
      output << iter->second << kSep;
    }
  }

  output << _id << kSep << _experimental_value << kSep;

  if (produce_descriptor_file) {
    output << accuracy();
  } else {
    output << "predicted " << _times_predicted;
    if (_times_predicted > 0) {
      output << " btw " << _predictions.minval() << " and " << _predictions.maxval()
             << " diffs btw " << _diffs.minval() << " and " << _diffs.maxval() << " ave "
             << accuracy();
    }
  }

  output << '\n';

  return 1;
}

template <typename T>
int
Entity_Continuous<T>::write_summary_data(IWString_and_File_Descriptor& output) const
{
  static constexpr char kSep = ' ';

  // clang-format off
  output << _id << kSep <<
            _experimental_value << kSep <<
            _diffs.n() << kSep <<
            _predictions.minval() << kSep <<
            _predictions.maxval() << kSep <<
            static_cast<float>(_predictions.average()) << kSep <<
            static_cast<float>(_diffs.average()) <<  kSep <<
            _diffs.maxval() << kSep <<
            static_cast<float>(sqrt(_diffs.variance())) << "\n";
  // clang-format on

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

class Classification_Data_Comparator
{
 private:
 public:
  int
  operator()(const Entity_Classification*, const Entity_Classification*) const;
};

int
Classification_Data_Comparator::operator()(const Entity_Classification* e1,
                                           const Entity_Classification* e2) const
{
  if (!sort_best_to_worst) {
    std::swap(e1, e2);
  }

  float a1 = e1->accuracy();
  float a2 = e2->accuracy();

  if (a1 < a2) {
    return 1;
  }
  if (a1 > a2) {
    return -1;
  }

  if (e1->times_predicted() < e2->times_predicted()) {
    return 1;
  }
  if (e1->times_predicted() > e2->times_predicted()) {
    return -1;
  }

  const IWString& c1 = e1->correct_class();
  const IWString& c2 = e2->correct_class();

  if (c1 == c2) {
    return 0;
  }

  if (c1.length() < c2.length()) {
    return -1;
  }

  if (c1.length() > c2.length()) {
    return 1;
  }

  return c1.strcmp(c2);
}

class Continuous_Data_Comparator
{
 private:
 public:
  int
  operator()(const Entity_Continuous<float>*, const Entity_Continuous<float>*) const;
};

int
Continuous_Data_Comparator::operator()(const Entity_Continuous<float>* e1,
                                       const Entity_Continuous<float>* e2) const
{
  float a1 = e1->accuracy();
  float a2 = e2->accuracy();

  if (a1 < a2) {
    return 1;
  }
  if (a1 > a2) {
    return -1;
  }

  if (e1->times_predicted() < e2->times_predicted()) {
    return 1;
  }
  if (e1->times_predicted() > e2->times_predicted()) {
    return -1;
  }

  return 0;
}

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
  cerr << R"(Identifies molecules that are poorly predicted across multiple prediction files.
 -E <fname>     experimental (Correct) results in <fname>
 -g             ignore missing experimental data
 -h             ignore duplicate entries in the -E file.
 -e <col>       experimental values in column <col>
 -p <col>       predicted    values in column <col>
 -C             classification problem
 -c <float>     cutoff for recomputing class membership using raw score
 -X -1,1        convert class labels to -1,1 based on prevalence (svmfp_calibrate)
 -F             run in a fault tolerant mode - skip seemingly incorrect files
 -s             consider signed differences rather than absolute
 -w <number>    only write the <number> worst predicted values
 -b             'balance' numeric filtering - approximately equal large and small values
 -W <value>     only write items with prediction errors above <value>
 -d             produce a 'descriptor' like file, just id and ave activity diff
 -k             input files were produced by svmfp_evaluate with the -c option (dist to SV)
 -D <fname>     produce file of distance vs abs prediction error
 -Y <fname>     produce file 'id ave_diff N'
 -R <regex>     process all files that match <regex> in current directory
 -S <fname>     file of smiles - output will include the smiles
 -r             reverse sort order, sort output from worst (top of file) to best (bottom)
 -v             verbose output
)";
  // clang-format on


  exit(rc);
}

static int
ReadSmilesLine(const_IWSubstring& buffer, IW_STL_Hash_Map_String& id_to_smiles) {
  IWString smiles, id;
  int i = 0;
  if (! buffer.nextword(smiles, i) || 
      ! buffer.nextword(id, i)) {
    cerr << "Cannot extract smiles and/or id\n";
    return 0;
  }

  id_to_smiles[id] = smiles;

  return 1;
}

static int
ReadSmiles(iwstring_data_source& input, IW_STL_Hash_Map_String& id_to_smiles) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! ReadSmilesLine(buffer, id_to_smiles)) {
      cerr << "ReadSmiles:cannot process '" << buffer << "'\n";
      return 0;
    }
  }

  return id_to_smiles.size();
}

static int
ReadSmiles(const char* fname, IW_STL_Hash_Map_String& id_to_smiles) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "ReadSmiles:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadSmiles(input, id_to_smiles);
}

typedef IW_STL_Hash_Map<IWString, Entity_Classification*> ID_to_Classification_Data;
typedef IW_STL_Hash_Map<IWString, Entity_Continuous<float>*> ID_to_Continuous_Data;

template <typename T>
int
do_summary(const T* ec, int n, IWString_and_File_Descriptor& output)
{
  for (auto i = 0; i < n; ++i) {
    ec[i]->write_summary_data(output);
  }

  return 1;
}

template <typename T>
int
report_data(const T* ec, int n, IWString_and_File_Descriptor& output)
{
  for (int i = 0; i < n; i++) {
    ec[i]->report(output);

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
report_classification_results(const ID_to_Classification_Data& idcd,
                              IWString_and_File_Descriptor& output)
{
  unsigned int n = idcd.size();
  if (0 == n) {
    cerr << "No data to process\n";
    return 0;
  }

  Entity_Classification** ec = new Entity_Classification*[n];
  std::unique_ptr<Entity_Classification*[]> free_ec(ec);

  int ndx = 0;
  for (ID_to_Classification_Data::const_iterator i = idcd.begin(); i != idcd.end(); ++i) {
    ec[ndx] = (*i).second;

    if (0 == ec[ndx]->times_predicted()) {
      continue;
    }

    ndx++;
  }

  n = ndx;

  if (0 == n) {
    return 0;
  }

  Classification_Data_Comparator cdc;
  iwqsort(ec, n, cdc);

  IW_STL_Hash_Map<IWString, Accumulator<double>*> acc;

  for (unsigned int i = 0; i < n; i++) {
    const Entity_Classification* eci = ec[i];

    const IWString& c = eci->correct_class();

    if (!acc.contains(c)) {
      acc[c] = new Accumulator<double>;
    }

    acc[c]->extra(eci->accuracy());
  }

  cerr << "Data on " << acc.size() << " classes\n";

  for (IW_STL_Hash_Map<IWString, Accumulator<double>*>::const_iterator i = acc.begin();
       i != acc.end(); ++i) {
    const IWString& c = (*i).first;

    const Accumulator<double>* a = (*i).second;

    cerr << "Class " << c << " " << a->n() << " items, accuracy between " << a->minval()
         << " and " << a->maxval() << " ave "
         << static_cast<float>(a->average_if_available_minval_if_not()) << "\n";
  }

  return report_data(ec, n, output);
}

static int
filter_to_worst_with_balancing(Entity_Continuous<float>** ec, unsigned int& n,
                               int write_worst_predicted)
{
#ifdef DEBUG_FILTER_TO_WORST
  cerr << "Filtering " << n << " items\n";

  for (int i = 0; i < n; i++) {
    cerr << ec[i]->accuracy() << '\n';
  }
#endif

  double low_sum = ec[0]->accuracy();
  double high_sum = -ec[n - 1]->accuracy();

  // cerr << "Highsum " << high_sum << '\n';

  int low_index = 1;
  int high_index = n - 2;

  int nsel = 2;

  // cerr << "low_index " << low_index << " high_index " << high_index << '\n';

  while (low_index < high_index && nsel < write_worst_predicted) {
    if (low_sum < high_sum) {
      if (ec[low_index]->accuracy() < 0.0F) {
        break;
      }

      low_sum += ec[low_index]->accuracy();
      low_index++;
    } else {
      if (ec[high_index]->accuracy() > 0.0F) {
        break;
      }

      high_sum -= ec[high_index]->accuracy();
      high_index--;
    }

    nsel++;
    //  cerr << "End of loop, high_sum " << high_sum << " low_sum " << low_sum << '\n';
  }

  // balance

  if (verbose > 1) {
    cerr << "Balancing        low " << low_index << " sum " << low_sum << " high "
         << (n - high_index) << " sum " << high_sum << '\n';
  }

  if (low_sum < high_sum) {
    while (low_sum + ec[low_index]->accuracy() < high_sum) {
      low_sum += ec[low_index]->accuracy();
      low_index++;
      if (low_index >= high_index) {
        break;
      }
    }
  } else {
    while (high_sum - ec[high_index]->accuracy() < low_sum) {
      high_sum -= ec[high_index]->accuracy();
      high_index--;
      if (high_index <= low_index) {
        break;
      }
    }
  }

  if (verbose > 1) {
    cerr << "After balancing, low " << low_index << " sum " << low_sum << " high "
         << (n - high_index) << " sum " << high_sum << '\n';
  }

  // now shift things down;

  // low_index++;

  for (unsigned int i = high_index + 1; i < n; i++) {
    ec[low_index] = ec[i];
    low_index++;
  }

  n = low_index;

  return n;
}

static int
filter_to_worst(Entity_Continuous<float>** ec, unsigned int& n, int write_worst_predicted)
{
  int low_index = 0;
  int high_index = n - 1;

  for (int nsel = 0; nsel < write_worst_predicted; nsel++) {
    if (ec[low_index]->accuracy() > fabs(ec[high_index]->accuracy())) {
      low_index++;
    } else {
      high_index--;
    }
  }

  for (unsigned int i = high_index + 1; i < n; i++) {
    ec[low_index] = ec[i];
    low_index++;
  }

  n = low_index;

  return n;
}

static int
filter_to_threshold(Entity_Continuous<float>** ec, unsigned int& n,
                    float write_differences_threshold)
{
  int low_index = 0;
  while (ec[low_index]->accuracy() >= write_differences_threshold) {
    low_index++;
  }

  unsigned int high_index = n - 1;

  while (-ec[high_index]->accuracy() >= write_differences_threshold) {
    high_index--;
  }

  high_index++;

  while (high_index < n) {
    ec[low_index] = ec[high_index];
    low_index++;
    high_index++;
  }

  n = low_index;

  return 1;
}

static int
report_continuous_results(const ID_to_Continuous_Data& idcd,
                          IWString_and_File_Descriptor& output)
{
  unsigned int n = idcd.size();
  if (0 == n) {
    cerr << "No data to process\n";
    return 0;
  }

  // cerr << "Start with " << n << " tiems\n";

  Entity_Continuous<float>** ec = new Entity_Continuous<float>*[n];
  std::unique_ptr<Entity_Continuous<float>*[]> free_ec(ec);

  int ndx = 0;
  for (ID_to_Continuous_Data::const_iterator i = idcd.begin(); i != idcd.end(); ++i) {
    ec[ndx] = (*i).second;

    //  cerr << "Item '" << (*i).first << "' predicted " << ec[ndx]->times_predicted() <<
    //  " times\n";

    if (0 == ec[ndx]->times_predicted()) {
      continue;
    }

    ndx++;
  }

  n = ndx;

  //  cerr << n << " items with non-zero observations\n";

  if (0 == n) {
    return 0;
  }

  Continuous_Data_Comparator cdc;

  iwqsort(ec, n, cdc);

  if (write_worst_predicted > 0) {
    if (balance_numeric_filtering) {
      filter_to_worst_with_balancing(ec, n, write_worst_predicted);
    } else {
      filter_to_worst(ec, n, write_worst_predicted);
    }
  } else if (write_differences_threshold > 0.0F) {
    filter_to_threshold(ec, n, write_differences_threshold);
  }

  if (stream_for_summary.active()) {
    do_summary(ec, n, stream_for_summary);
  }

  return report_data(ec, n, output);
}

static int
get_token_from_column(const const_IWSubstring& buffer, int& i, int desired_column,
                      const_IWSubstring& s)
{
  int col = 1;
  while (buffer.nextword(s, i)) {
    if (col == desired_column) {
      return 1;
    }

    col++;
  }

  cerr << "No column " << (desired_column + 1) << " in '" << buffer << "'\n";
  return 0;
}

static int
parse_experimental_data_record(const const_IWSubstring& buffer,
                               ID_to_Continuous_Data& idd)
{
  if (buffer.nwords() <= experimental_column) {
    cerr << "Must be at least " << (experimental_column + 1)
         << " tokens in an experimental data record\n";
    return 0;
  }

  IWString id;
  int i = 0;
  buffer.nextword(id, i);

  if (! idd.contains(id)) {
    // new value, great
  } else if (ignore_duplicate_experimental_values) {
    ++duplicate_experimental_values_ignored;
  } else {
    cerr << "Duplicate experimental data for '" << id << "', cannot continue\n";
    return 0;
  }

  Entity_Continuous<float>* e = new Entity_Continuous<float>(id);

  const_IWSubstring token;

  if (!get_token_from_column(buffer, i, experimental_column, token)) {
    return 0;
  }

  float v;
  if (!token.numeric_value(v)) {
    cerr << "Invalid experimental value '" << token << "'\n";
    return 0;
  }

  e->set_experimental_value(v);

  idd[id] = e;

  return idd.size();
}

static int
read_experimental_data_continuous(iwstring_data_source& input, ID_to_Continuous_Data& idd)
{
  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    cerr << "Empty experimental data file\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (!parse_experimental_data_record(buffer, idd)) {
      cerr << "Invalid experimental data record '" << buffer << "'\n";
      return 0;
    }
  }

  return idd.size();
}

static int
read_experimental_data_continuous(const const_IWSubstring& fname,
                                  ID_to_Continuous_Data& idd)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open experimental data file '" << fname << "'\n";
    return 0;
  }

  return read_experimental_data_continuous(input, idd);
}

static int
parse_experimental_data_record(const const_IWSubstring& buffer,
                               ID_to_Classification_Data& idd)
{
  if (buffer.nwords() <= experimental_column) {
    cerr << "Must be at least " << (experimental_column + 1)
         << " tokens in an experimental data record\n";
    return 0;
  }

  IWString id;
  int i = 0;
  buffer.nextword(id, i);

  if (idd.contains(id)) {
    cerr << "Duplicate experimental data for '" << id << "', cannot continue\n";
    return 0;
  }

  Entity_Classification* e = new Entity_Classification(id);

  const_IWSubstring token;

  if (!get_token_from_column(buffer, i, experimental_column, token)) {
    return 0;
  }

  e->set_class(token);

  idd[id] = e;

  return idd.size();
}

static int
read_experimental_data_classification(iwstring_data_source& input,
                                      ID_to_Classification_Data& idd)
{
  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    cerr << "Empty experimental data file\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (!parse_experimental_data_record(buffer, idd)) {
      cerr << "Invalid experimental data record '" << buffer << "'\n";
      return 0;
    }
  }

  return idd.size();
}

static int
read_experimental_data_classification(const const_IWSubstring& fname,
                                      ID_to_Classification_Data& idd)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open experimental data file '" << fname << "'\n";
    return 0;
  }

  input.set_translate_tabs(1);

  return read_experimental_data_classification(input, idd);
}

static int
process_classification_prediction(const const_IWSubstring& buffer,
                                  ID_to_Classification_Data& idd,
                                  const IW_STL_Hash_Map_String& xref)
{
  if (buffer.nwords() <= predicted_column) {
    cerr << "Experimental data must contain at least " << (predicted_column + 1)
         << " tokens\n";
    return 0;
  }

  IWString id;
  int i = 0;

  buffer.nextword(id, i);

  ID_to_Classification_Data::const_iterator f = idd.find(id);

  if (f == idd.end()) {
    cerr << "No experimental data for '" << id << "'\n";
    return ignore_missing_experimental_data;
  }

  const_IWSubstring predicted;

  if (!get_token_from_column(buffer, i, predicted_column, predicted)) {
    return 0;
  }

  Entity_Classification* e = (*f).second;

  if (use_cutoff_to_convert_to_class) {
    float p;
    if (!predicted.numeric_value(p)) {
      cerr << "Invalid prediction '" << predicted << "', non numeric\n";
      return 0;
    }
    if (p < cutoff) {
      predicted = "-1";
    } else {
      predicted = "1";
    }
  }

  predicted_values_read++;

  if (0 == xref.size()) {
    e->predicted_as(predicted);
    return 1;
  }

  const auto g = xref.find(predicted);

  if (g == xref.end()) {
    cerr << "No class label translation for '" << predicted << "'\n";
    return 0;
  }

  e->predicted_as(g->second);

  return 1;
}

static int
gather_predicted_data_classification(iwstring_data_source& input,
                                     ID_to_Classification_Data& idd,
                                     const IW_STL_Hash_Map_String& xref)
{
  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    cerr << "Empty prediction file\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (!process_classification_prediction(buffer, idd, xref)) {
      cerr << "Bad predicted data record '" << buffer << "'\n";
      return fault_tolerant_mode;
    }
  }

  return 1;
}

static int
gather_predicted_data_classification(const char* fname, ID_to_Classification_Data& idd,
                                     const IW_STL_Hash_Map_String& xref)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_translate_tabs(1);

  return gather_predicted_data_classification(input, idd, xref);
}

static int
process_continuous_prediction(const const_IWSubstring& buffer, ID_to_Continuous_Data& idd)
{
  if (buffer.nwords() <= predicted_column) {
    cerr << "Experimental data must contain at least " << (predicted_column + 1)
         << " tokens\n";
    return 0;
  }

  IWString id;
  int i = 0;

  buffer.nextword(id, i);

  ID_to_Continuous_Data::const_iterator f = idd.find(id);

  if (f == idd.end()) {
    cerr << "No experimental data for '" << id << "'\n";
    return ignore_missing_experimental_data;
  }

  const_IWSubstring predicted;

  if (!get_token_from_column(buffer, i, predicted_column, predicted)) {
    return 0;
  }

  predicted_values_read++;

  Entity_Continuous<float>* e = (*f).second;

  if (!data_includes_distance_to_support_vector) {
    e->predicted_as(predicted);

    return 1;
  }

  const_IWSubstring dist;

  if (!buffer.nextword(dist, i)) {
    cerr << "Did not find distance to support vector\n";
    return 0;
  }

  return e->prediction_and_distance(predicted, dist);
}

static int
gather_predicted_data_continuous(iwstring_data_source& input, ID_to_Continuous_Data& idd)
{
  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    cerr << "Empty prediction file\n";
    return fault_tolerant_mode;
  }

  while (input.next_record(buffer)) {
    if (!process_continuous_prediction(buffer, idd)) {
      cerr << "Bad predicted data record '" << buffer << "'\n";
      return fault_tolerant_mode;
    }
  }

  return 1;
}

static int
gather_predicted_data_continuous(const char* fname,
                                 ID_to_Continuous_Data& idcd);  // forward declaration

static int
gather_predicted_data_continuous_from_file(const char* fname, ID_to_Continuous_Data& idcd)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return fault_tolerant_mode;
  }

  IWString buffer;

  while (input.next_record(buffer)) {
    buffer.strip_trailing_blanks();
    if (0 == buffer.length()) {
      continue;
    }

    if (buffer.starts_with('#')) {
      continue;
    }

    if (!gather_predicted_data_continuous(buffer.null_terminated_chars(), idcd)) {
      return 0;
    }
  }

  return 1;
}

static int
gather_predicted_data_continuous(const char* fname, ID_to_Continuous_Data& idcd)
{
  if (strlen(fname) > 2 && 0 == strncmp(fname, "F:", 2)) {
    return gather_predicted_data_continuous_from_file(fname + 2, idcd);
  }

  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 2) {
    cerr << "Processing '" << fname << "'\n";
  }

  return gather_predicted_data_continuous(input, idcd);
}

static int
gather_predicted_data_classification(DIR* dir, std::unique_ptr<re2::RE2>& fname_rx,
                                     ID_to_Classification_Data& idcd,
                                     const IW_STL_Hash_Map_String& xref)
{
  struct dirent* d;

  int files_processed = 0;

  while (NULL != (d = readdir(dir))) {
    const_IWSubstring fname(d->d_name);
    if (!iwre2::RE2PartialMatch(fname, *fname_rx)) {
      continue;
    }

    files_processed++;

    if (!gather_predicted_data_classification(d->d_name, idcd, xref)) {
      return 0;
    }
  }

  return files_processed;
}

static int
gather_predicted_data_classification(std::unique_ptr<re2::RE2>& fname_rx,
                                     ID_to_Classification_Data& idcd,
                                     const IW_STL_Hash_Map_String& xref)
{
  DIR* dir = opendir(".");

  if (NULL == dir) {
    cerr << "Cannot open directory\n";
    return 0;
  }

  int rc = gather_predicted_data_classification(dir, fname_rx, idcd, xref);

  closedir(dir);

  return rc;
}

static int
gather_predicted_data_continuous(DIR* dir, std::unique_ptr<re2::RE2>& fname_rx,
                                 ID_to_Continuous_Data& idcd)
{
  struct dirent* d;

  int files_processed = 0;

  while (NULL != (d = readdir(dir))) {
    const_IWSubstring fname(d->d_name);
    if (!iwre2::RE2PartialMatch(fname, *fname_rx)) {
      continue;
    }

    files_processed++;

    if (!gather_predicted_data_continuous(d->d_name, idcd)) {
      return 0;
    }
  }

  if (0 == files_processed) {
    cerr << "Warning no files matched file name regular expression\n";
  }

  return files_processed;
}

static int
gather_predicted_data_continuous(std::unique_ptr<re2::RE2>& fname_rx,
                                 ID_to_Continuous_Data& idcd)
{
  DIR* dir = opendir(".");

  if (NULL == dir) {
    cerr << "Cannot open directory\n";
    return 0;
  }

  int rc = gather_predicted_data_continuous(dir, fname_rx, idcd);

  closedir(dir);

  return rc;
}

/*
  Pretty much hard coded for two classes, but would not be hard to expand
*/

static int
determine_class_label_translation(ID_to_Classification_Data& idcd,
                                  const const_IWSubstring& s,
                                  IW_STL_Hash_Map_String& xref)
{
  IWString l1, l2;

  if (!s.split(l1, ',', l2) || 0 == l1.length() || 0 == l2.length()) {
    cerr << "The class label translation must be of the form 'l1,l2', '" << s
         << "' not valid\n";
    return 0;
  }

  IW_STL_Hash_Map_int count;  // for each class, how many instances

  for (auto i = idcd.cbegin(); i != idcd.cend(); ++i) {
    count[i->second->correct_class()]++;
  }

  const auto nc = count.size();

  if (2 != nc) {
    cerr << "determine_class_label_translation:sorry, only works on two classes, not "
         << count.size() << '\n';
    return 0;
  }

  std::vector<std::pair<IWString, int>> ccount(count.cbegin(), count.cend());

  std::sort(ccount.begin(), ccount.end(),
            [](const std::pair<IWString, int>& c1, const std::pair<IWString, int>& c2) {
              return c1.second > c2.second;
            });

  for (unsigned int i = 0; i < nc; ++i) {
    const auto c = ccount[i];

    if (verbose) {
      cerr << c.second << " instances of " << c.first << '\n';
    }

    if (0 == i) {
      xref[l1] = c.first;
    } else {
      xref[l2] = c.first;
    }
  }

  return 1;
}

static int
report_coverage(const ID_to_Classification_Data& idcd, std::ostream& os)
{
  Accumulator_Int<int> predicted;
  int not_predicted = 0;

  for (ID_to_Classification_Data::const_iterator i = idcd.begin(); i != idcd.end(); ++i) {
    const Entity_Classification* e = (*i).second;
    predicted.extra(e->times_predicted());
    if (0 == e->times_predicted()) {
      not_predicted++;
    }
  }

  if (0 == predicted.n()) {
    return 0;
  }

  os << idcd.size() << " items predicted between " << predicted.minval() << " and "
     << predicted.maxval() << " times. Ave "
     << static_cast<float>(predicted.average_if_available_minval_if_not()) << "\n";

  return os.good();
}

static int
mispredicted(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vgCE:e:p:c:Fsw:W:dbY:R:kX:rS:h");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    ignore_missing_experimental_data = 1;

    if (verbose) {
      cerr << "Will ignore missing experimental data\n";
    }
  }

  if (cl.option_present('p')) {
    if (!cl.value('p', predicted_column) || predicted_column < 2) {
      cerr << "The predicted column must be a valid column number > 1\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Predicted values in column " << predicted_column << '\n';
    }

    predicted_column--;
  }

  if (cl.option_present('e')) {
    if (!cl.value('e', experimental_column) || experimental_column < 2) {
      cerr << "The experimental column must be a valid column number > 1\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Experimental values in column " << experimental_column << '\n';
    }

    experimental_column--;
  }

  if (cl.option_present('F')) {
    fault_tolerant_mode = 1;

    if (verbose) {
      cerr << "WIll ignore most errors\n";
    }
  }

  if (cl.option_present('s')) {
    take_absolute_differences = 0;

    if (verbose) {
      cerr << "Will report signed differences\n";
    }
  }

  if (cl.option_present('k')) {
    data_includes_distance_to_support_vector = 1;

    if (verbose) {
      cerr << "Will look for extra column of distance to nearest support vector\n";
    }
  }

  if (cl.option_present('w')) {
    if (!cl.value('w', write_worst_predicted) || write_worst_predicted < 1) {
      cerr << "The write worst predicted items (-w) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will write the " << write_worst_predicted << " worst predicted items\n";
    }

    take_absolute_differences = 0;

    if (cl.option_present('b')) {
      balance_numeric_filtering = 1;

      if (verbose) {
        cerr << "Will produce approximately balanced samples of high and low "
                "differences\n";
      }
    }
  }

  if (cl.option_present('W')) {
    if (!cl.value('W', write_differences_threshold) ||
        write_differences_threshold <= 0.0) {
      cerr << "The write differences threshold (-W) option must be a +ve real value\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Will write items predicted worse than " << write_differences_threshold
           << " from obs\n";
    }

    take_absolute_differences = 0;
  }

  if (cl.option_present('d')) {
    produce_descriptor_file = 1;

    if (verbose) {
      cerr << "Will produce an activity file\n";
    }
  }

  int classification = cl.option_present('C');

  if (cl.option_present('X')) {
    classification = 1;
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', cutoff)) {
      cerr << "The class cutoff (-c) must be a valid number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "Predictions less than " << cutoff << " will go to the '-1' class\n";
    }

    use_cutoff_to_convert_to_class = 1;
    classification = 1;
  }

  if (!cl.option_present('E')) {
    cerr << "Must specify experimental (Correct) data via the -E option\n";
    usage(3);
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (! ReadSmiles(fname.null_terminated_chars(), id_to_smiles)) {
      cerr << "Cannot read smiles '" << fname << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Read " << id_to_smiles.size() << " id to smiles from " << fname << "n";
    }
  }

  if (cl.option_present('R')) {
    ;
  } else if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  const_IWSubstring experimental_data_fname = cl.string_value('E');

  IWString_and_File_Descriptor output(1);

  if (produce_descriptor_file) {
    static constexpr char kSep = ' ';

    if (id_to_smiles.size()) {
      output << "SMILES" << kSep;
    }
    if (classification) {
      output << "ID EXPT N Correct Accuracy";
    } else {
      output << "ID EXPT AVDIFF";
    }
    output << '\n';
  }

  if (cl.option_present('Y')) {
    if (!take_absolute_differences) {
      cerr << "Warning, summary file being generated, but absolute differences not being "
              "used\n";
    }

    const char* y = cl.option_value('Y');

    if (!stream_for_summary.open(y)) {
      cerr << argv[0] << " cannot open summary file '" << y << "'\n";
      return 2;
    }

    if (verbose) {
      cerr << "Summary data written to '" << y << "'\n";
    }

    stream_for_summary << "ID Expt N minpred maxpred avepred avediff maxdiff stddiff\n";
  }

  if (cl.option_present('r')) {
    sort_best_to_worst = 0;

    if (verbose) {
      cerr << "Will sort worst to best in output file\n";
    }
  }

  if (cl.option_present('h')) {
    ignore_duplicate_experimental_values = 1;
    if (verbose) {
      cerr << "Will ignore duplicate experimental values in the -E file\n";
    }
  }

  std::unique_ptr<re2::RE2> fname_rx;

  bool fname_rx_active = false;

  if (cl.option_present('R')) {
    const const_IWSubstring r = cl.string_value('R');

    if (!iwre2::RE2Reset(fname_rx, r)) {
      fname_rx_active = true;
    }
    if (verbose) {
      cerr << "File name regular expression initialised '" << r << "'\n";
    }

    //  cerr << "match? " << std::regex_search("PREDQ", fname_rx) << '\n';
  }

  int rc = 0;

  if (classification) {
    ID_to_Classification_Data idcd;

    if (!read_experimental_data_classification(experimental_data_fname, idcd)) {
      cerr << "Cannot read experimental data from '" << experimental_data_fname << "'\n";
      return 3;
    }

    if (verbose) {
      cerr << "Read " << idcd.size() << " experimental data values from '" << experimental_data_fname << "'\n";
    }

    IW_STL_Hash_Map_String xref;

    if (cl.option_present('X')) {
      const_IWSubstring x = cl.string_value('X');

      if (!determine_class_label_translation(idcd, x, xref)) {
        cerr << "Cannot translate class labels based on '" << x << "'\n";
        return 1;
      }
    }

    if (fname_rx_active) {
      if (!gather_predicted_data_classification(fname_rx, idcd, xref)) {
        return 1;
      }
    } else {
      for (int i = 0; i < cl.number_elements(); i++) {
        if (!gather_predicted_data_classification(cl[i], idcd, xref)) {
          rc = i + 1;
          break;
        }
      }
    }

    if (verbose) {
      cerr << "Read " << predicted_values_read << " predicted values\n";
    }

    if (0 == rc) {
      report_classification_results(idcd, output);
      output.flush();
      report_coverage(idcd, cerr);
    }
  } else {
    ID_to_Continuous_Data idcd;

    if (!read_experimental_data_continuous(experimental_data_fname, idcd)) {
      cerr << "Cannot read experimental data from '" << experimental_data_fname << "'\n";
      return 3;
    }

    if (verbose) {
      cerr << "Read " << idcd.size() << " experimental data values from '" << experimental_data_fname << "'\n";
    }

    if (fname_rx_active) {
      if (!gather_predicted_data_continuous(fname_rx, idcd)) {
        return 1;
      }
    } else {
      for (int i = 0; i < cl.number_elements(); i++) {
        if (!gather_predicted_data_continuous(cl[i], idcd)) {
          rc = i + 1;
          break;
        }
      }
    }

    if (verbose) {
      cerr << "Read " << predicted_values_read << " predicted values\n";
    }

    if (0 == rc) {
      report_continuous_results(idcd, output);
      output.flush();
    }
  }

  if (verbose) {
  }

  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = mispredicted(argc, argv);

  return rc;
}
