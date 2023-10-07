/*
  Construct test train splits based on leader clustering
*/

#include <stdlib.h>

#include <iostream>
#include <random>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwqsort/iwqsort.h"

using std::cerr;
using std::endl;

const char* prog_name = NULL;

static int verbose = 0;

static int number_splits = 1;

static IWString output_stem = "ttCLsplit";

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static int molecules_in_training_set = 0;

static int max_attempts = 5;

static int duplicate_splits_discarded = 0;

static int nbits = 0;  // multiple of 8 near number clusters

static int write_header_to_split_files = 0;

static int write_smiles = 0;

static Accumulator_Int<int> split_size;

using random_number_t = uint32_t;

static std::random_device rd;
static std::mt19937 rng(rd());
static std::uniform_int_distribution<random_number_t> uniform(
    0, std::numeric_limits<random_number_t>::max());

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
cerr << R"(Creates train/test splits based on leader clustering.
All members of each cluster are assigned to either train or test randomly.
Run gfp_leader and pipe data to this programme\n";
 -S <stem>      prefix for output files <stem> (default 'ttCLsplit')
 -p <pct>       percent of molecules in the training set (default 50)
 -n <number>    number of splits to create
 -i             write smiles files rather than identifier files
 -h             write a header record to each of the split files
 -v             verbose output
)";
 // clang-format on

  exit(rc);
}

class Smiles_and_ID {
 private:
  const IWString _smiles;
  const IWString _id;

 public:
  Smiles_and_ID(const IWString&, const IWString&);

  const IWString&
  smiles() const {
    return _smiles;
  }

  const IWString&
  id() const {
    return _id;
  }
};

Smiles_and_ID::Smiles_and_ID(const IWString& s, const IWString& i) : _smiles(s), _id(i) {
  return;
}

class Cluster {
 private:
  resizable_array_p<Smiles_and_ID> _sid;

  // Clusters are assigned random values and sorted based on that random value.
  // We should also use shuffle.
  random_number_t _random_value;

  int _cluster_id;

  // Is this cluster in train (0) or test (1).
  int _set_membership;

  //  private functions

  int _get_smiles_and_id(iwstring_data_source& input, int& fatal);

 public:
  int build(iwstring_data_source&, int&);

  void
  set_set_membership(int s) {
    _set_membership = s;
  }

  int
  set_membership() const {
    return _set_membership;
  }

  void
  choose_random_number() {
    _random_value = uniform(rng);
  }

  random_number_t
  random_number() const {
    return _random_value;
  }

  int
  cluster_size() const {
    return _sid.number_elements();
  }

  void
  set_cluster_id(int s) {
    _cluster_id = s;
  }

  int
  cluster_id() const {
    return _cluster_id;
  }

  int write_smiles(IWString_and_File_Descriptor& output) const;
  int write_identifiers(IWString_and_File_Descriptor& output) const;
};

/*
  Input looks like

$SMI<smiles1>
PCN<id1>
CLUSTER<0>
CSIZE<4>
$SMI<smiles2>
PCN<id2>
DIST<0.32>
$SMI<smiles3>

*/

int
Cluster::build(iwstring_data_source& input, int& fatal) {
  fatal = 0;

  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    return 0;
  }

  input.push_record();

  const_IWSubstring smiles, id;
  if (!_get_smiles_and_id(input, fatal)) {
    cerr << "Cluster::build:invalid cluster, near line " << input.lines_read() << endl;
    return 0;
  }

  if (!input.next_record(buffer)) {
    cerr << "Cluster::build:premature EOF\n";
    fatal = 1;
    return 0;
  }

  if (!buffer.starts_with("CLUSTER<")) {
    cerr << "Cluster::build:CLUSTER record invalid '" << buffer << "'\n";
    return 0;
  }

  if (!input.next_record(buffer)) {
    cerr << "Cluster::build:premature EOF\n";
    fatal = 1;
    return 0;
  }

  if (!buffer.starts_with("CSIZE<")) {
    cerr << "Cluster::build:CSIZE record invalid '" << buffer << "'\n";
    fatal = 1;
    return 0;
  }

  while (_get_smiles_and_id(input, fatal)) {
  }

  if (fatal) {
    return 0;
  }
  return 1;
}

static int
extract_tdt_value(const const_IWSubstring& buffer, int tag_length, IWString& result) {
  if (!buffer.ends_with('>')) {
    cerr << "TDT items must end in > '" << buffer << "'\n";
    return 0;
  }

  int bstop = buffer.length() - 2;
  if (bstop < tag_length)  // happens with 'PCN<>'
  {
    result = "";
    return 1;
  }

  buffer.from_to(tag_length, buffer.length() - 2, result);

  return 1;
}

int
Cluster::_get_smiles_and_id(iwstring_data_source& input, int& fatal) {
  const_IWSubstring buffer;

  while (1) {
    if (!input.next_record(buffer)) {
      cerr << "Cluster::_get_smiles_and_id:premature EOF\n";
      fatal = 1;
      return 0;
    }

    if ('|' == buffer) {
      fatal = 0;
      return 0;
    }

    if (buffer.starts_with(smiles_tag)) {
      break;
    }
  }

  assert(buffer.starts_with(smiles_tag));

  IWString smiles;
  if (!extract_tdt_value(buffer, smiles_tag.length(), smiles)) {
    cerr << "Cluster::_get_smiles_and_id:cannot extract smiles from '" << buffer << "'\n";
    fatal = 1;
    return 0;
  }

  if (!input.next_record(buffer)) {
    cerr << "Cluster::_get_smiles_and_id:premature EOF\n";
    fatal = 1;
    return 0;
  }

  if (!buffer.starts_with(identifier_tag)) {
    cerr << "Cluster::_get_smiles_and_id:not PCN '" << buffer << "'\n";
    return 0;
  }

  IWString id;
  if (!extract_tdt_value(buffer, identifier_tag.length(), id)) {
    cerr << "Cluster::_get_smiles_and_id:cannot extract ID from '" << buffer << "'\n";
    fatal = 1;
    return 0;
  }

  Smiles_and_ID* s = new Smiles_and_ID(smiles, id);

  _sid.add(s);

  return 1;
}

int
Cluster::write_smiles(IWString_and_File_Descriptor& output) const {
  int n = _sid.number_elements();

  for (int i = 0; i < n; i++) {
    const Smiles_and_ID* si = _sid[i];

    output << si->smiles() << ' ' << si->id() << '\n';

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
Cluster::write_identifiers(IWString_and_File_Descriptor& output) const {
  int n = _sid.number_elements();

  for (int i = 0; i < n; i++) {
    const Smiles_and_ID* si = _sid[i];

    output << si->id() << '\n';

    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

class Random_Number_Comparator {
 private:
 public:
  int operator()(const Cluster*, const Cluster*) const;
};

int
Random_Number_Comparator::operator()(const Cluster* c1, const Cluster* c2) const {
  random_number_t r1 = c1->random_number();
  random_number_t r2 = c2->random_number();

  if (r1 < r2) {
    return -1;
  }
  if (r1 > r2) {
    return 1;
  }
  return 0;
}

static int
write_split(const resizable_array_p<Cluster>& cluster, int flag,
            IWString_and_File_Descriptor& output) {
  if (write_header_to_split_files) {
    output << "ID\n";
  }

  int nc = cluster.number_elements();

  for (int i = 0; i < nc; i++) {
    const Cluster* ci = cluster[i];

    if (flag != ci->set_membership()) {
      continue;
    }

    if (write_smiles) {
      ci->write_smiles(output);
    } else {
      ci->write_identifiers(output);
    }
  }

  return 1;
}

static int
write_split(const resizable_array_p<Cluster>& cluster, IWString& fname, int flag) {
  IWString_and_File_Descriptor output;

  if (!output.open(fname.null_terminated_chars())) {
    cerr << "Cannot open split file '" << fname << "'\n";
    return 0;
  }

  return write_split(cluster, flag, output);
}

static int
write_split(const IWString& output_stem, int ndx,
            const resizable_array_p<Cluster>& cluster) {
  IWString fname;

  if (write_smiles) {
    fname << output_stem << 'R' << ndx << ".smi";
  } else {
    fname << output_stem << ndx << ".train";
  }

  write_split(cluster, fname, 0);

  fname.resize_keep_storage(0);

  if (write_smiles) {
    fname << output_stem << 'E' << ndx << ".smi";
  } else {
    fname << output_stem << ndx << ".test";
  }

  write_split(cluster, fname, 1);

  return 1;
}

static Random_Number_Comparator rnc;

static int
find_unique_split(resizable_array_p<Cluster>& cluster,
                  resizable_array_p<IW_Bits_Base>& splits_formed) {
  for (Cluster* c : cluster) {
    c->choose_random_number();
  }

  cluster.iwqsort(rnc);

  const int nc = cluster.number_elements();

  int ntrain = 0;

  for (int i = 0; i < nc; i++) {
    int n = cluster[i]->cluster_size();

    if ((ntrain + n) <= molecules_in_training_set)  // not there yet
    {
      cluster[i]->set_set_membership(0);
      ntrain += n;
      continue;
    }

    //  OK, we just passed the dividing point. Should this point go in train or test?

    if (molecules_in_training_set - ntrain < ntrain + n - molecules_in_training_set) {
      cluster[i]->set_set_membership(0);
      split_size.extra(ntrain);
    } else {
      cluster[i]->set_set_membership(1);
      split_size.extra(ntrain + n);
    }

    if (verbose > 1) {
      cerr << "At cluster " << i << " ntrain " << ntrain << " got " << n << " items\n";
    }

    for (int j = i + 1; j < nc; j++) {
      cluster[j]->set_set_membership(1);
    }
    break;
  }

  IW_Bits_Base* b = new IW_Bits_Base(nbits);

  for (int i = 0; i < nc; i++) {
    if (0 == cluster[i]->set_membership()) {  // in training set
      continue;
    }

    int j = cluster[i]->cluster_id();

    b->set(j);
  }

  int ns = splits_formed.number_elements();

  if (0 == ns) {
    splits_formed.add(b);
    return 1;
  }

  for (int i = 0; i < ns; i++) {
    const IW_Bits_Base* si = splits_formed[i];

    if ((*b) == (*si)) {
      delete b;
      duplicate_splits_discarded++;
      return 0;
    }
  }

  splits_formed.add(b);

  return 1;
}

static int
test_train_split_cluster(const IWString& output_stem, int ndx,
                         resizable_array_p<IW_Bits_Base>& splits_formed,
                         resizable_array_p<Cluster>& cluster) {
  for (int a = 0; a < max_attempts; a++) {
    if (find_unique_split(cluster, splits_formed)) {
      return write_split(output_stem, ndx, cluster);
    }
  }

  cerr << "Cannot find a unique split\n";
  return 0;  // no unique split found
}

static int
read_cluster_data(iwstring_data_source& input, resizable_array_p<Cluster>& cluster) {
  while (1) {
    Cluster* c = new Cluster;
    int fatal;
    if (c->build(input, fatal)) {
      cluster.add(c);
    } else {
      delete c;
      if (fatal) {
        return 0;
      }
      return 1;
    }
  }

  cerr << "Should not come to here\n";
  return 0;
}

static int
read_cluster_data(const char* fname, resizable_array_p<Cluster>& cluster) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_cluster_data(input, cluster);
}

static int
test_train_split_cluster(int argc, char** argv) {
  Command_Line cl(argc, argv, "vp:S:n:ih");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('n')) {
    if (!cl.value('n', number_splits) || number_splits < 1) {
      cerr << "The number of splits to create option (-n) must be a whole positive "
              "number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will create " << number_splits << " random splits\n";
    }
  }

  if (cl.option_present('i')) {
    write_smiles = 1;

    if (verbose) {
      cerr << "Will write smiles files\n";
    }
  }

  float fraction_in_training_set = 0.5;

  if (cl.option_present('p')) {
    int p;
    if (!cl.value('p', p) || p < 1 || p > 99) {
      cerr << "The percent in training set option (-p) must be a valid percentage\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Will put " << p << " percent of the items in the training set\n";
    }

    fraction_in_training_set = static_cast<float>(p) / static_cast<float>(100.0);
  }

  if (cl.option_present('h')) {
    if (cl.option_present('i')) {
      cerr << "The -i (write smiles) and -h (write header) options are mutually "
              "incompatible\n";
      usage(3);
    }

    write_header_to_split_files = 1;

    if (verbose) {
      cerr << "Will write a header record to each split\n";
    }
  }

  if (cl.option_present('S')) {
    output_stem = cl.string_value('S');

    if (verbose) {
      cerr << "Will create file(s) will stem '" << output_stem << "'\n";
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1) {
    cerr << "Only makes sense with one input file\n";
    return 3;
  }

  resizable_array_p<Cluster> cluster;

  for (int i = 0; i < cl.number_elements(); i++) {
    if (!read_cluster_data(cl[i], cluster)) {
      cerr << "Cannot read cluster data in '" << cl[i] << "'\n";
      return 3;
    }
  }

  int nc = cluster.number_elements();

  if (0 == nc) {
    cerr << "No data\n";
    return 2;
  }

  if (nc < 2) {
    cerr << "Only one cluster, cannot process\n";
    return 3;
  }

  if (0 == nc % 8) {
    nbits = nc;
  } else {
    nbits = (nc / 8 + 1) * 8;
  }

  int nmolecules = 0;
  Accumulator_Int<int> cluster_size;
  for (int i = 0; i < nc; i++) {
    cluster[i]->set_cluster_id(i);
    int n = cluster[i]->cluster_size();
    nmolecules += n;
    if (verbose) {
      cluster_size.extra(n);
    }
  }

  if (verbose) {
    cerr << "Read data on " << nmolecules << " molecules in " << nc << " clusters\n";
    cerr << "Cluster sizes between " << cluster_size.minval() << " and "
         << cluster_size.maxval() << " ave "
         << cluster_size.average_if_available_minval_if_not() << endl;
  }

  molecules_in_training_set =
      static_cast<int>(fraction_in_training_set * nmolecules + 0.4999);

  if (verbose) {
    cerr << "Will try to place " << molecules_in_training_set
         << " molecules in the training set\n";
  }

  // Use the cluster id to keep track of which splits have been formed

  resizable_array_p<IW_Bits_Base> splits_formed;

  for (int i = 0; i < number_splits; i++) {
    if (!test_train_split_cluster(output_stem, i + 1, splits_formed, cluster)) {
      cerr << "Cannot create split " << i << endl;
      return i + 1;
    }
  }

  if (verbose) {
    int n = split_size.n();

    cerr << "Created " << n << " splits\n";
    cerr << "Sizes between " << split_size.minval() << " and " << split_size.maxval()
         << " ave " << split_size.average_if_available_minval_if_not() << endl;
  }

  return 0;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = test_train_split_cluster(argc, argv);

  return rc;
}
