/*
  Performs stratified sampling
*/

#include <stdlib.h>

#include <algorithm>
#include <memory>
#include <random>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

static int activity_column = 1;

static int records_to_select = 2;

static int percent_of_records_to_select = 0;

static int nsplit = 1;

static int number_strata = 10;

static int items_per_stratum = 1;

static int header_records_to_skip = 0;

/*
  We want to correctly echo header records
*/

static resizable_array_p<IWString> header_records;

static IWString_and_File_Descriptor stream_for_test_set;

static int just_write_identifiers = 0;

static IWString stem_for_training_set;
static IWString stem_for_test_set;

static IWString suffix;

static int venetian_stratification = 0;

static int balanced_with_respect_to_first_stratum = 0;
static double balanced_with_respect_to_first_stratum_fraction = 0.0;

static int balanced_with_respect_to_first_stratum_value_valid = 0;
static double balanced_with_respect_to_first_stratum_value = 0.0;

/*
  The default is to take the top N active values, then another N
  stratified from the rest of the dataset.
*/

static float balanced_with_respect_to_first_stratum_expand = 1.0f;

/*
  We can also randomly sample from withing that first stratum
*/

static float sample_fraction_within_first_stratum = 0.0f;

static IW_STL_Hash_Map_String id_to_smiles;

static int seq_start = 0;

static int produce_chronological_split = 0;

/*
  When building subset models, we often have extra values in the activity file, but want
  to restrict the sampling to just those values for which we have the smiles
*/

static int remove_activity_values_not_in_smiles_file = 0;

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
  cerr << "Performs multiple stratified samplings of a file\n";
  cerr << " -E <fname>     write test set records to <fname>\n";
  cerr << " -R <fname>     write training set records to <fname>\n";
  cerr << " -c <col>       activity in column<c>\n";
  cerr << " -N <nsplit>    number of splits to create\n";
  cerr << " -n <nrecs>     number of items in training set(s)\n";
  cerr << " -p <pct>       percentage of items in training set(s)\n";
  cerr << " -s <nrecs>     header records to skip\n";
  cerr << " -a <nrecs>     sample the most active <nrecs> items, equal size stratified sample of the rest\n";
  cerr << " -a <p>%        sample the p% most active items, equal size stratified sample of the rest\n";
  cerr << " -a .gt.<y>     sample all records > y, equal size stratified sample of the rest\n";
  cerr << " -a sample=<x>  within the most active items selected above, randomly sample fraction <x> of them\n";
  cerr << " -a expand=<x>  rather than an equal number of less active items, multiply the number active by <x>\n";
  cerr << " -i             just write the identifiers, not whole records\n";
  cerr << " -S <suffix>    create files with suffix <suffix>\n";
  cerr << " -M <fname>     produce corresponding smiles files, smiles in <fname>\n";
  cerr << " -K             drop input values for which smiles not available\n";
  cerr << " -b <start>     normally files produced start with 0, start with <start> instead\n";
  cerr << " -C             produce an extra chronological split\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

/*
  There are some variables here that are no longer used in the new version.
*/

class ID_Stratum_Selected
{
 private:
  IWString _id;

  IWString _buffer;

  int _stratum;

  float _activity;

  int _selected_this_iteration;

  int _times_selected_across_splits;

  int _id_converted_to_number;

 public:
  ID_Stratum_Selected();

  const IWString&
  id() const
  {
    return _id;
  }

  int
  initialise(const const_IWSubstring& buffer);

  int
  stratum() const
  {
    return _stratum;
  }

  void
  set_stratum(int s)
  {
    _stratum = s;
  }

  float
  activity() const
  {
    return _activity;
  }

  int
  times_selected_across_splits() const
  {
    return _times_selected_across_splits;
  }

  void
  set_selected_this_iteration();

  void
  set_selected_this_iteration(int s)
  {
    _selected_this_iteration = s;
  }

  int
  selected_this_iteration() const
  {
    return _selected_this_iteration;
  }

  //  we need the input file if we are echoing everything

  int
  do_write(IWString_and_File_Descriptor&) const;

  int
  id_converted_to_number() const
  {
    return _id_converted_to_number;
  }
};

template class resizable_array_p<ID_Stratum_Selected>;
template class resizable_array_base<ID_Stratum_Selected*>;

ID_Stratum_Selected::ID_Stratum_Selected()
{
  _stratum = -9999;

  _activity = static_cast<float>(-1234567.0);

  _times_selected_across_splits = 0;

  _selected_this_iteration = 0;

  _id_converted_to_number = -1;

  return;
}

void
ID_Stratum_Selected::set_selected_this_iteration()
{
  _selected_this_iteration++;  // should never be larger than 1

  _times_selected_across_splits++;

  return;
}

int
ID_Stratum_Selected::initialise(const const_IWSubstring& buffer)
{
  _buffer = buffer;

  int i = 0;

  if (!buffer.nextword(_id, i)) {
    cerr << "ID_Stratum_Selected::initialise:cannot extract identifier\n";
    return 0;
  }

  _id_converted_to_number = -1;

  if (_id.numeric_value(_id_converted_to_number)) {  // breat
    ;
  } else {
    const_IWSubstring tmp(_id);
    tmp.remove_leading_chars('0');
    if (!tmp.numeric_value(_id_converted_to_number)) {
      _id_converted_to_number = -1;
    }
  }

  IWString c;

  if (1 == activity_column) {
    buffer.nextword(c, i);
  } else {
    const_IWSubstring token;
    int col = 1;
    while (buffer.nextword(token, i)) {
      col++;

      if (col == activity_column) {
        c = token;
        break;
      }
    }
  }

  if (0 == c.length()) {
    cerr << "ID_Stratum_Selected::initialise:no activity data for '" << _id << "'\n";
    return 0;
  }

  if (!c.numeric_value(_activity)) {
    cerr << "ID_Stratum_Selected::initialise:invalid activity '" << c << "'\n";
    return 0;
  }

  return 1;
}

int
ID_Stratum_Selected::do_write(IWString_and_File_Descriptor& output) const
{
  if (just_write_identifiers) {
    output << _id << '\n';
  } else {
    output << _buffer << '\n';
  }

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
smiles_available_for_all_identifiers(resizable_array_p<ID_Stratum_Selected>& idss,
                                     const IW_STL_Hash_Map_String& id_to_smiles)
{
  int number_removed = 0;

  for (auto i = idss.number_elements() - 1; i >= 0; --i) {
    const IWString& id = idss[i]->id();

    if (id_to_smiles.contains(id)) {
      continue;
    }

    if (remove_activity_values_not_in_smiles_file) {
      idss.remove_item(i);
      number_removed++;
    } else {
      cerr << "NO smiles for '" << id << "'\n";
      return 0;
    }
  }

  if (0 == idss.size()) {
    cerr << "smiles_available_for_all_identifiers:all values removed\n";
    return 0;
  }

  if (verbose) {
    cerr << "Dropped " << number_removed << " values for which smiles not available\n";
  }

  return 1;
}

/*
  Our first sort, in order to assign strata, must be just by activity
*/

class IDSS_Activity_Comparator
{
 private:
 public:
  int
  operator()(const ID_Stratum_Selected*, const ID_Stratum_Selected*) const;
};

int
IDSS_Activity_Comparator::operator()(const ID_Stratum_Selected* pidss1,
                                     const ID_Stratum_Selected* pidss2) const
{
  float s1 = pidss1->activity();
  float s2 = pidss2->activity();

  if (s1 < s2) {
    return 1;
  }

  if (s1 > s2) {
    return -1;
  }

  return 0;
}

// #define DEBUG_DO_BALANCED_WITH_RESPECT_TO_FIRST_STRATUM

static int
do_balanced_with_respect_to_first_stratum(
    resizable_array_p<ID_Stratum_Selected>& idds,
    const resizable_array<int>& stratum_start, const resizable_array<int>& stratum_stop,
    IWString_and_File_Descriptor& training_set_output,
    IWString_and_File_Descriptor& test_set_output)
{
#ifdef DEBUG_DO_BALANCED_WITH_RESPECT_TO_FIRST_STRATUM
  assert(balanced_with_respect_to_first_stratum > 0 &&
         balanced_with_respect_to_first_stratum < n);
#endif

  std::random_device rd;
  std::mt19937_64 rng(rd());

#ifdef DEBUG_DO_BALANCED_WITH_RESPECT_TO_FIRST_STRATUM
  cerr << "From " << n << " items selected " << balanced_with_respect_to_first_stratum
       << " remaining " << (n - balanced_with_respect_to_first_stratum) << endl;
#endif

  if (sample_fraction_within_first_stratum) {
    std::shuffle(idds.rawdata(), idds.rawdata() + stratum_stop[0], rng);
    int istop =
        static_cast<int>(sample_fraction_within_first_stratum * stratum_stop[0] + 0.4999);
    if (0 == istop) {  // should never happen
      istop = 1;
    }

    for (int i = 0; i < istop; ++i) {
      idds[i]->set_selected_this_iteration();
      idds[i]->do_write(training_set_output);
    }
  } else {
    for (int i = 0; i < stratum_stop[0]; ++i)  // all of first stratum gets written
    {
      idds[i]->set_selected_this_iteration();
      idds[i]->do_write(training_set_output);
    }
  }

  int items_per_stratum =
      static_cast<int>((balanced_with_respect_to_first_stratum *
                        balanced_with_respect_to_first_stratum_expand) /
                           (number_strata - 1) +
                       0.499999f);

  if (0 == items_per_stratum) {
    items_per_stratum = 1;
  }

  int nextra = static_cast<int>((balanced_with_respect_to_first_stratum *
                                 balanced_with_respect_to_first_stratum_expand) -
                                ((number_strata - 1) * items_per_stratum) + 0.49999f);

  for (int s = 1; s < number_strata; ++s) {
    std::shuffle(idds.rawdata() + stratum_start[s], idds.rawdata() + stratum_stop[s],
                 rng);

    int istop = stratum_start[s] + items_per_stratum;

    if (nextra) {
      istop++;
      nextra--;
    }

    for (int i = stratum_start[s]; i < istop; ++i) {
      idds[i]->set_selected_this_iteration();
      idds[i]->do_write(training_set_output);
    }

    if (test_set_output.is_open()) {
      for (int i = istop; i < stratum_stop[s]; ++i) {
        idds[i]->do_write(test_set_output);
      }
    }
  }

  return 1;
}

static int
echo_header_records(const resizable_array_p<IWString>& header_records,
                    IWString_and_File_Descriptor& output)
{
  for (int i = 0; i < header_records.number_elements(); i++) {
    const IWString& s = *(header_records[i]);

    output << s << '\n';
  }

  return 1;
}

static int
do_chronological_split(resizable_array_p<ID_Stratum_Selected>& idds,
                       IWString_and_File_Descriptor& training_set_output,
                       IWString_and_File_Descriptor& test_set_output)
{
  const int n = idds.number_elements();

  const int ntrain = n - records_to_select;
  for (int i = 0; i < ntrain; ++i) {
    idds[i]->do_write(training_set_output);
  }

  if (test_set_output.is_open()) {
    for (int i = ntrain; i < n; ++i) {
      idds[i]->do_write(test_set_output);
    }
  }

  return 1;
}

static int
do_chronological_split(resizable_array_p<ID_Stratum_Selected>& idds, const int ndx)
{
  idds.iwqsort_lambda([](const ID_Stratum_Selected* s1, const ID_Stratum_Selected* s2) {
    return s1->id_converted_to_number() < s2->id_converted_to_number();
  });

  IWString fname;
  fname << stem_for_training_set << ndx << suffix;
  IWString_and_File_Descriptor training_set_output;
  if (!training_set_output.open(fname.c_str())) {
    cerr << "Cannot open chronological training set output '" << fname << "'\n";
    return 0;
  }

  IWString_and_File_Descriptor test_set_output;
  if (stem_for_test_set.length()) {
    fname = stem_for_test_set;
    fname << ndx << suffix;

    if (!test_set_output.open(fname.c_str())) {
      cerr << "Cannot open chronological test set output '" << fname << "'\n";
      return 0;
    }
  }

  return do_chronological_split(idds, training_set_output, test_set_output);
}

static int
do_split(resizable_array_p<ID_Stratum_Selected>& idss, int iteration_number,
         const resizable_array<int>& stratum_start,
         const resizable_array<int>& stratum_stop,
         IWString_and_File_Descriptor& training_set_output,
         IWString_and_File_Descriptor& test_set_output)
{
  echo_header_records(header_records, training_set_output);
  if (test_set_output.is_open()) {
    echo_header_records(header_records, test_set_output);
  }

  int n = idss.number_elements();

  for (int i = 0; i < n; ++i) {
    idss[i]->set_selected_this_iteration(0);
  }

  if (balanced_with_respect_to_first_stratum) {
    return do_balanced_with_respect_to_first_stratum(
        idss, stratum_start, stratum_stop, training_set_output, test_set_output);
  }

  std::random_device rd;
  std::mt19937_64 rng(rd());

  for (int i = 0; i < number_strata; i++) {
    std::shuffle(idss.rawdata() + stratum_start[i], idss.rawdata() + stratum_stop[i],
                 rng);

    int train_end = stratum_start[i] + records_to_select / number_strata;
    assert(train_end < stratum_stop[i]);

    // cerr << "Processing stratum " << i << " from " << stratum_start[i] << " to " <<
    // stratum_stop[i] << " items " << (train_end - stratum_start[i]) << endl;;

    for (int j = stratum_start[i]; j < train_end; ++j) {
      idss[j]->set_selected_this_iteration();
      idss[j]->do_write(training_set_output);
    }

    if (test_set_output.is_open()) {
      for (int j = train_end; j < stratum_stop[i]; ++j) {
        idss[j]->do_write(test_set_output);
      }
    }
  }

  return 1;
}

static int
_do_split(resizable_array_p<ID_Stratum_Selected>& idss, int ndx,
          const resizable_array<int>& stratum_start,
          const resizable_array<int>& stratum_stop,
          IWString_and_File_Descriptor& training_set_output)
{
  IWString_and_File_Descriptor test_set_output;

  if (0 == stem_for_test_set.length()) {
    return do_split(idss, ndx, stratum_start, stratum_stop, training_set_output,
                    test_set_output);
  }

  IWString fname;

  if (1 == nsplit) {
    fname = stem_for_test_set;
  } else {
    fname << stem_for_test_set << ndx << suffix;
  }

  if (!test_set_output.open(fname.null_terminated_chars())) {
    cerr << "Cannot open test set file '" << fname << "'\n";
    return 0;
  }

  return do_split(idss, ndx, stratum_start, stratum_stop, training_set_output,
                  test_set_output);
}

template <typename T>
int
write_smiles(const resizable_array_p<ID_Stratum_Selected>& idss, const char* fname,
             const T& c)
{
  IWString_and_File_Descriptor output;

  if (!output.open(fname)) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  const int n = idss.number_elements();

  for (int i = 0; i < n; ++i) {
    const ID_Stratum_Selected* idssi = idss[i];

    if (!c(idssi->selected_this_iteration())) {
      continue;
    }

    const auto f = id_to_smiles.find(idssi->id());

    output << f->second << ' ' << f->first << '\n';

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
do_split(resizable_array_p<ID_Stratum_Selected>& idss, int ndx,
         const resizable_array<int>& stratum_start,
         const resizable_array<int>& stratum_stop,
         IWString_and_File_Descriptor& training_set_output)
{
  if (!_do_split(idss, ndx, stratum_start, stratum_stop, training_set_output)) {
    return 0;
  }

  if (0 == id_to_smiles.size()) {
    return 1;
  }

  IWString fname;

  fname << stem_for_training_set << ndx << ".smi";

  if (!write_smiles(idss, fname, [](const int s) { return s > 0;})) {
    return 0;
  }

  if (stem_for_test_set.length()) {
    IWString fname;

    fname << stem_for_test_set << ndx << ".smi";

    if (!write_smiles(idss, fname, [](const int s) { return s == 0;})) {
      return 0;
    }
  }

  return 1;
}

static int
assign_strata_activity_balanced(const resizable_array_p<ID_Stratum_Selected>& idds,
                                resizable_array<int>& stratum_start,
                                resizable_array<int>& stratum_stop)
{
  const int n = idds.number_elements();

  if (balanced_with_respect_to_first_stratum_value_valid) {
    for (int i = n - 1; i >= 0; --i)  // array is sorted highest to smallest
    {
      //    cerr << " i = " << i << " activity " << idds[i]->activity() << endl;
      if (idds[i]->activity() < balanced_with_respect_to_first_stratum_value) {
        continue;
      }

      balanced_with_respect_to_first_stratum = i;
      break;
    }

    if (balanced_with_respect_to_first_stratum < 0) {
      cerr << "No values found greater than "
           << balanced_with_respect_to_first_stratum_value << endl;
      return 0;
    }

    if (verbose) {
      cerr << "Found " << balanced_with_respect_to_first_stratum << " of " << n
           << " values above " << balanced_with_respect_to_first_stratum_value << endl;
    }
  }

  assert(balanced_with_respect_to_first_stratum > 0);

  const int items_remaining =
      idds.number_elements() - balanced_with_respect_to_first_stratum;

  int items_per_stratum = items_remaining / (number_strata - 1);

  if (items_per_stratum <= 0) {
    cerr << "Starting with " << idds.number_elements() << " and selecting "
         << balanced_with_respect_to_first_stratum
         << " in first stratum, does not allow filling " << number_strata << " strata\n";
    return 0;
  }

  if (1 == items_per_stratum) {
    cerr << "Yipes, one item per stratum, very unlikely this is what you want\n";
  }

  stratum_start.add(0);
  stratum_stop.add(balanced_with_respect_to_first_stratum);

  for (int i = 0; i < (number_strata - 1); ++i) {
    stratum_start.add(balanced_with_respect_to_first_stratum + i * items_per_stratum);
    stratum_stop.add(stratum_start.last_item() + items_per_stratum);
  }

  stratum_stop[stratum_stop.number_elements() - 1] = idds.number_elements();

  if (verbose > 1) {
    for (int i = 0; i < number_strata; ++i) {
      cerr << " stratum " << i << " from " << stratum_start[i] << " to "
           << stratum_stop[i] << ", " << (stratum_stop[i] - stratum_start[i])
           << " items\n";
    }
  }

  return 1;
}

static int
assign_strata(const resizable_array_p<ID_Stratum_Selected>& idds,
              resizable_array<int>& stratum_start, resizable_array<int>& stratum_stop)
{
  stratum_start.resize(number_strata);
  stratum_stop.resize(number_strata);

  if (balanced_with_respect_to_first_stratum ||
      balanced_with_respect_to_first_stratum_value_valid) {
    return assign_strata_activity_balanced(idds, stratum_start, stratum_stop);
  }

  for (int s = 0; s < number_strata; ++s) {
    stratum_start.add(s * items_per_stratum);
    stratum_stop.add((s + 1) * items_per_stratum);
  }

  stratum_stop[number_strata - 1] = idds.number_elements();

  return 1;
}

static int
all_ids_numeric(const resizable_array_p<ID_Stratum_Selected>& idss)
{
  const int n = idss.number_elements();

  for (int i = 0; i < n; ++i) {
    if (idss[i]->id_converted_to_number() < 0) {
      return 0;
    }
  }

  return 1;
}

static int
stratified_sample(iwstring_data_source& input)
{
  const_IWSubstring buffer;

  for (int i = 0; i < header_records_to_skip; i++) {
    if (!input.next_record(buffer)) {
      cerr << "Cannot skip " << header_records_to_skip << " header records\n";
      return 0;
    }

    if (just_write_identifiers) {
      header_records.add(new IWString("ID"));
    } else {
      header_records.add(new IWString(buffer));
    }
  }

  resizable_array_p<ID_Stratum_Selected> idss;

  Accumulator<float> acc;

  while (input.next_record(buffer)) {
    ID_Stratum_Selected* tmp = new ID_Stratum_Selected;

    if (!tmp->initialise(buffer)) {
      cerr << "Fatal error on line " << input.lines_read() << " '" << buffer << "'\n";
      delete tmp;
      return 0;
    }

    acc.extra(tmp->activity());

    idss.add(tmp);
  }

  int n = idss.number_elements();

  if (0 == n) {
    cerr << "No data\n";
    return 0;
  }

  if (n <= 2) {
    cerr << "Too few records " << n << " in this file to do anything sensible\n";
    return 0;
  }

  if (id_to_smiles.size()) {
    if (!smiles_available_for_all_identifiers(idss, id_to_smiles)) {
      return 0;
    }

    n = idss.number_elements();
  }

  if (number_strata >= n) {
    cerr << "Input contains " << n << " items, cannot use " << number_strata
         << " strata\n";
    return 0;
  }

  items_per_stratum = n / number_strata;

  if (verbose) {
    cerr << "Input contains " << n << " items, " << items_per_stratum
         << " items assigned to each stratum\n";
  }

  if (verbose) {
    cerr << "Activities between " << acc.minval() << " and " << acc.maxval() << " ave "
         << static_cast<float>(acc.average()) << '\n';
  }

  if (1 == items_per_stratum) {
    cerr << "Input contains " << n << " items, number_strata " << number_strata
         << " leaves just one item per stratum\n";
  }

  if (produce_chronological_split && !all_ids_numeric(idss)) {
    cerr << "Chronological split requested, but all ids not numeric\n";
    return 0;
  }

  if (percent_of_records_to_select > 0) {
    records_to_select = (n * percent_of_records_to_select) / 100;

    if (0 == records_to_select) {
      records_to_select = 1;
    }

    if (verbose) {
      cerr << "Input contains " << n << " records, will select " << records_to_select
           << " items in each split\n";
    }
  } else if (verbose) {
    cerr << "Input contains " << n << " items\n";
  }

  if (records_to_select >= n) {
    cerr << "Asked for " << records_to_select << " records, but input only contains " << n
         << endl;
    return 0;
  }

  if (balanced_with_respect_to_first_stratum_fraction > 0.0 &&
      0 == balanced_with_respect_to_first_stratum) {
    balanced_with_respect_to_first_stratum = static_cast<int>(
        balanced_with_respect_to_first_stratum_fraction * static_cast<double>(n) +
        0.4999);
    if (balanced_with_respect_to_first_stratum < 1 ||
        balanced_with_respect_to_first_stratum == n) {
      cerr << "Cannot accommodate " << balanced_with_respect_to_first_stratum_fraction
           << " as most active fraction\n";
      return 0;
    }

    if (verbose) {
      cerr << "Balance most active fraction "
           << balanced_with_respect_to_first_stratum_fraction << " translates to "
           << balanced_with_respect_to_first_stratum << " items\n";
    }
  }

  // All methods depend on the initial array being sorted

  IDSS_Activity_Comparator idssac;
  idss.iwqsort(idssac);

  resizable_array<int> stratum_start, stratum_stop;

  if (!assign_strata(idss, stratum_start, stratum_stop)) {
    cerr << "Cannot assign strata to " << idss.number_elements() << " dataitems\n";
    return 0;
  }

// #define DEBUG_STRATA
#ifdef DEBUG_STRATA
  cerr << "Strata for " << idss.size() << " points, into " << number_strata
       << " strata\n";
  for (int i = 0; i < number_strata; ++i) {
    cerr << "Stratum " << i << " start " << stratum_start[i] << " to " << stratum_stop[i]
         << ", " << (stratum_stop[i] - stratum_start[i]) << " values\n";
  }
#endif

  // Handle the case of just one split specially

  if (1 == nsplit) {
    if (0 == stem_for_training_set.length()) {
      IWString_and_File_Descriptor output(1);
      do_split(idss, 0, stratum_start, stratum_stop, output);
    } else {
      IWString_and_File_Descriptor output;

      if (!output.open(stem_for_training_set.null_terminated_chars())) {
        cerr << "Cannot open training set output file '" << stem_for_training_set
             << "'\n";
        return 0;
      }

      do_split(idss, 0, stratum_start, stratum_stop, output);
    }
  } else {
    for (int i = 0; i < nsplit; i++) {
      IWString fname;

      fname << stem_for_training_set << (i + seq_start) << suffix;

      IWString_and_File_Descriptor training_set_output;

      if (!training_set_output.open(fname.null_terminated_chars())) {
        cerr << "Cannot open training set file '" << fname << "'\n";
        return 0;
      }

      if (!do_split(idss, i, stratum_start, stratum_stop, training_set_output)) {
        return 0;
      }
    }
  }

  if (produce_chronological_split) {
    do_chronological_split(idss, nsplit + seq_start);
  }

  return 1;

  if (verbose > 1) {
    IDSS_Activity_Comparator f;
    idss.iwqsort(f);  // should sort within each stratum...

    for (int i = 0; i < n; i++) {
      const ID_Stratum_Selected* idssi = idss[i];

      int c = idssi->stratum();

      cerr << idssi->id() << " stratum " << c << " train "
           << idssi->times_selected_across_splits() << " test "
           << (nsplit - idssi->times_selected_across_splits()) << '\n';
    }
  }

  return 1;
}

static int
stratified_sample(const char* fname)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return stratified_sample(input);
}

static int
read_smiles(iwstring_data_source& input, IW_STL_Hash_Map_String& id_to_smiles)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    IWString smiles, id;

    if (!buffer.split(smiles, ' ', id) || 0 == smiles.length() || 0 == id.length()) {
      cerr << "Cannot split into smiles and id '" << buffer << "'\n";
      return 0;
    }

    id.truncate_at_first(' ');

    id_to_smiles[id] = smiles;
  }

  return id_to_smiles.size();
}

static int
read_smiles(const char* fname, IW_STL_Hash_Map_String& id_to_smiles)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open smiles file '" << fname << "'\n";
    return 0;
  }

  return read_smiles(input, id_to_smiles);
}

static int
stratified_sample(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vE:R:c:N:n:p:s:iS:a:M:b:CK");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    if (!cl.value('c', activity_column) || activity_column < 0) {
      cerr << "The column with the class information must be a whole +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Class data in column " << activity_column << endl;
    }

    activity_column--;
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', header_records_to_skip) || header_records_to_skip < 1) {
      cerr << "The header records to skip (-s) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will skip " << header_records_to_skip << " header records\n";
    }
  }

  if (cl.option_present('i')) {
    just_write_identifiers = 1;
    if (verbose) {
      cerr << "Will just write identifier fields rather than whole records\n";
    }
  }

  if (cl.option_present('n') && cl.option_present('p')) {
    cerr << "The -p and -n options are mutually exclusive\n";
    usage(5);
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', records_to_select) || records_to_select < 2) {
      cerr << "The number of records to select (-n) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Requested " << records_to_select << " records\n";
    }
  } else if (cl.option_present('p')) {
    if (!cl.value('p', percent_of_records_to_select) ||
        percent_of_records_to_select < 1 || percent_of_records_to_select >= 100) {
      cerr << "The percent of records per sample (-p) must be valid percentage\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will choose " << percent_of_records_to_select << "% of the records\n";
    }
  } else if (cl.option_present('a')) {
    ;
  } else {
    cerr << "Must sspecify the size of the training set via the -n, -p or -a options\n";
    usage(2);
  }

  if (cl.option_present('R')) {
    cl.value('R', stem_for_training_set);

    if (verbose) {
      cerr << "Training set files created with stem '" << stem_for_training_set << "'\n";
    }
  }

  if (cl.option_present('N')) {
    if (!cl.value('N', nsplit) || nsplit < 1) {
      cerr << "The nsplit option (-N) must be a whole +ve number\n";
      usage(3);
    }

    if (1 == nsplit) {
      ;
    } else if (0 == stem_for_training_set.length()) {
      cerr << "When producing multiple splits, you must specify the -R option\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will produce " << nsplit << " splits of the data\n";
    }
  }

  if (cl.option_present('V')) {
    if (cl.option_present('t') || cl.option_present('p')) {
      cerr << "The -V and (-t, -p) options are mutually exclusive\n";
      usage(1);
    }

    venetian_stratification = 1;

    if (verbose) {
      cerr << "Will do Venetian blind stratification\n";
    }

    number_strata = nsplit;

    if (number_strata < 2) {
      cerr << "Doing Venetian blind stratification, but only one split, does not make "
              "sense\n";
      return 2;
    }
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', number_strata) || number_strata < 2) {
      cerr << "The number of strata (-t) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will divide the data into " << number_strata << " strata\n";
    }
  }

  if (cl.option_present('M')) {
    const char* fname = cl.option_value('M');

    if (!read_smiles(fname, id_to_smiles)) {
      cerr << "Cannot read smiles from '" << fname << "'\n";
      return 0;
    }

    if (verbose) {
      cerr << "Read " << id_to_smiles.size() << " identifier/smiles values from '"
           << fname << "'\n";
    }

    if (cl.option_present('K')) {
      remove_activity_values_not_in_smiles_file = 1;

      if (verbose) {
        cerr << "Will drop input values for which no smiles available\n";
      }
    }
  }

  if (cl.option_present('E')) {
    cl.value('E', stem_for_test_set);

    if (verbose) {
      cerr << "Test set records written to '" << stem_for_test_set << "'\n";
    }
  }

  if (cl.option_present('S')) {
    cl.value('S', suffix);

    if (verbose) {
      cerr << "Files created with suffix '" << suffix << "'\n";
    }
  }

  if (cl.option_present('a')) {
    const_IWSubstring a;
    for (int i = 0; cl.value('a', a, i); ++i) {
#ifdef DEBUG_DO_BALANCED_WITH_RESPECT_TO_FIRST_STRATUM
      cerr << "Building balanced sampling from '" << a << "'\n";
#endif

      if (a.starts_with(".gt.")) {
        a.remove_leading_chars(4);
        if (!a.numeric_value(balanced_with_respect_to_first_stratum_value)) {
          cerr << "The numeric cutoff for first stratum \"-a .gt.x\" must be a valid "
                  "number\n";
          usage(1);
        }

        if (verbose) {
          cerr << "All values above " << balanced_with_respect_to_first_stratum_value
               << " automatically selected\n";
        }

        balanced_with_respect_to_first_stratum_value_valid = 1;
      } else if (a.starts_with("expand=")) {
        a.remove_leading_chars(7);
        if (!a.numeric_value(balanced_with_respect_to_first_stratum_expand) ||
            balanced_with_respect_to_first_stratum_expand <= 0.0f) {
          cerr << "The expand balanced with respect to first stratum \"-a expand=n\" "
                  "must be a whole +ve number\n";
          usage(1);
        }

        if (verbose) {
          cerr << "Across other strata, will select "
               << balanced_with_respect_to_first_stratum_expand
               << " times as many points as from first stratum\n";
        }
      } else if (a.ends_with('%') || a.ends_with("pct")) {
        if (a.ends_with('%')) {
          a.chop();
        } else {
          a.chop(3);
        }

        if (!a.numeric_value(balanced_with_respect_to_first_stratum_fraction) ||
            balanced_with_respect_to_first_stratum_fraction < 1.0 ||
            balanced_with_respect_to_first_stratum_fraction > 99.0) {
          cerr << "Invalid active stratum selection '-a " << a << "%'\n";
          usage(1);
        }

#ifdef DEBUG_DO_BALANCED_WITH_RESPECT_TO_FIRST_STRATUM
        cerr << "After conversion " << balanced_with_respect_to_first_stratum_fraction
             << endl;
#endif
        balanced_with_respect_to_first_stratum_fraction /= 100.0;
#ifdef DEBUG_DO_BALANCED_WITH_RESPECT_TO_FIRST_STRATUM
        cerr << "Converted to fraction "
             << balanced_with_respect_to_first_stratum_fraction << endl;
#endif
      } else if (a.starts_with("sample=")) {
        a.remove_leading_chars(7);
        if (!a.numeric_value(sample_fraction_within_first_stratum) ||
            sample_fraction_within_first_stratum <= 0.0 ||
            sample_fraction_within_first_stratum >= 1.0) {
          cerr << "The active stratum sample fraction must be a valid fraction, '" << a
               << "' invalid\n";
          return 2;
        }

        if (verbose) {
          cerr << "Will randomly sample " << sample_fraction_within_first_stratum
               << " of the active stratum\n";
        }
      } else if (a.numeric_value(balanced_with_respect_to_first_stratum_fraction) &&
                 balanced_with_respect_to_first_stratum_fraction > 0.0 &&
                 balanced_with_respect_to_first_stratum_fraction < 1.0) {
        ;
      } else if (a.numeric_value(balanced_with_respect_to_first_stratum) &&
                 balanced_with_respect_to_first_stratum > 1) {
      } else {
        cerr << "Invalid most active stratum selection '-a " << a << "'\n";
        usage(1);
      }
    }

    if (verbose) {
      cerr << "Balanced fraction " << balanced_with_respect_to_first_stratum_fraction
           << " or balanced N " << balanced_with_respect_to_first_stratum << endl;
    }
  }

  if (cl.option_present('b')) {
    if (!cl.value('b', seq_start) || seq_start < 0) {
      cerr << "The file sequence start value option (-b) must be a non-negative number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Output files start at " << seq_start << endl;
    }
  }

  if (cl.option_present('C')) {
    produce_chronological_split = 1;

    if (verbose) {
      cerr << "Will produce an extra chronological split\n";
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (1 != cl.number_elements()) {
    cerr << "Sorry, don't know how to handle multiple input files\n";
    usage(3);
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!stratified_sample(cl[i])) {
      rc = i + 1;
      break;
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

  int rc = stratified_sample(argc, argv);

  return rc;
}
