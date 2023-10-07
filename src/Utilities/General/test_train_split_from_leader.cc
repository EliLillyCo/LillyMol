/*
  Suntara's idea of generating splits based on leader clusters.
  We anticipate that the workflow will be to run:
  1.  gfp_leader_multiple_thresholds to produce a bunch of leader
      results.
  2.  We read one or more of those files and decide which one(s) to
      use in order to produce splits of the appropriate size

*/

#include <stdlib.h>

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

using std::cerr;
using std::endl;

const char* prog_name = NULL;

static int verbose = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static int nsplit = 1;

static IWString suffix;

static IWString stem_for_training_set;
static IWString stem_for_test_set;

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
  cerr << "Creates train/test splits based on leader clustering results (possibly from gfp_leader_multiple_thresholds)\n";
  cerr << " -E <fname>     write test set records to <fname>\n";
  cerr << " -R <fname>     write training set records to <fname>\n";
  cerr << " -S <suffix>    create files with <suffix>\n";
  cerr << " -N <nsplit>    number of splits to create\n";
  cerr << "                two ways of specifying number of items in training set:\n";
  cerr << " -n <nitems>    number  of items in the training set\n";
  cerr << " -p <pct>       percent of items in the training set\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

class Leader_Cluster
{
 private:
  int _leader;
  resizable_array<int> _nbr;

 public:
  Leader_Cluster();

  void
  set_leader(const int s)
  {
    _leader = s;
  }

  int
  leader() const
  {
    return _leader;
  }

  int
  build(iwstring_data_source&, const IW_STL_Hash_Map_int& id_to_ndx);
};

Leader_Cluster::Leader_Cluster()
{
  _leader = -1;
}

int
Leader_Cluster::build(iwstring_data_source& input, const IW_STL_Hash_Map_int& id_to_ndx)
{
  IWString buffer;

  _leader = -1;

  while (input.next_record(buffer)) {
    if (!buffer.starts_with(identifier_tag)) {
      continue;
    }

    buffer.remove_leading_chars(identifier_tag.length());
    buffer.chop();
    buffer.truncate_at_first(' ');

    const auto f = id_to_ndx.find(buffer);

    if (f == id_to_ndx.end()) {
      cerr << "Leader_Cluster::build:no info for '" << buffer << "'\n";
      return 0;
    }

    if (_leader < 0) {
      _leader = f->second;
    } else {
      _nbr.add(f->second);
    }
  }

  return _leader >= 0;
}

class Leader_Results
{
 private:
  int _nclusters;

  Leader_Cluster* _cluster;

 public:
  Leader_Results();
  ~Leader_Results();

  int
  build(const char*, const IW_STL_Hash_Map_int& id_to_ndx);
  int
  build(iwstring_data_source&, const IW_STL_Hash_Map_int& id_to_ndx);
};

Leader_Results::Leader_Results()
{
  _nclusters = 0;
  _cluster = nullptr;

  return;
}

Leader_Results::~Leader_Results()
{
  if (nullptr != _cluster) {
    delete[] _cluster;
  }

  return;
}

int
Leader_Results::build(const char* fname, const IW_STL_Hash_Map_int& id_to_ndx)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Leader_Results::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input, id_to_ndx);
}

int
Leader_Results::build(iwstring_data_source& input, const IW_STL_Hash_Map_int& id_to_ndx)
{
  std::unique_ptr<re2::RE2> rx = std::make_unique<re2::RE2>("^|");

  _nclusters = input.grep(*rx);

  if (0 == _nclusters) {
    cerr << "Leader_Results::build:no clusters in input\n";
    return 0;
  }

  _cluster = new Leader_Cluster[_nclusters];

  for (int i = 0; i < _nclusters; ++i) {
    const int cstart = input.lines_read();  // in case of error

    if (_cluster[i].build(input, id_to_ndx)) {
      ;
    } else if (input.eof()) {
      return _nclusters;
    } else {
      cerr << "Leader_Cluster::build:cannot build cluster starting at line " << cstart
           << endl;
      return 0;
    }
  }

  return _nclusters;
}

static int
test_train_split_from_leader(iwstring_data_source& input,
                             const IW_STL_Hash_Map_int& id_to_ndx,
                             const IW_STL_Hash_Map_String& id_to_smiles,
                             IWString_and_File_Descriptor& output)
{
  return 1;
}

static int
test_train_split_from_leader(const char* fname, const IW_STL_Hash_Map_int& id_to_ndx,
                             const IW_STL_Hash_Map_String& id_to_smiles,
                             IWString_and_File_Descriptor& output)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return test_train_split_from_leader(input, id_to_ndx, id_to_smiles, output);
}

/*static int
read_smiles_record (const const_IWSubstring & buffer,
                    IW_STL_Hash_Map_int & id_to_ndx,
                    IW_STL_Hash_Map_String & id_to_smiles)
{
  IWString smiles, id;
  int i = 0;
  buffer.nextword(smiles, i);
  if (0 == smiles.length() || ! buffer.nextword(id, i))
    return 0;

  id_to_smiles.emplace(id, smiles);

  const int s = id_to_ndx.size();
  id_to_ndx.emplace(id, s);

  return 1;
}*/

static int
get_tag_value(const_IWSubstring buffer,  // note local copy
              const int nremove, IWString& s)
{
  buffer.remove_leading_chars(nremove);
  buffer.chop();
  s = buffer;

  return s.length();
}

static int
read_smiles(iwstring_data_source& input, IW_STL_Hash_Map_int& id_to_ndx,
            IW_STL_Hash_Map_String& id_to_smiles)
{
  const_IWSubstring buffer;

  IWString smiles;

  while (input.next_record(buffer)) {
    if (buffer.starts_with(smiles_tag)) {
      get_tag_value(buffer, smiles_tag.length(), smiles);
      continue;
    }

    if (!buffer.starts_with(identifier_tag)) {
      continue;
    }

    IWString id;
    get_tag_value(buffer, identifier_tag.length(), id);

    const auto f = id_to_ndx.find(id);
    if (f != id_to_ndx.end()) {
      continue;
    }

    id_to_smiles.emplace(id, smiles);
  }

  return id_to_smiles.size();
}

static int
read_smiles(const char* fname, IW_STL_Hash_Map_int& id_to_ndx,
            IW_STL_Hash_Map_String& id_to_smiles)
{
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "read_smiles:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_smiles(input, id_to_ndx, id_to_smiles);
}

static int
test_train_split_from_leader(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vS:E:R:N:n:p:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

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

  const int n = cl.number_elements();

  if (0 == n) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IW_STL_Hash_Map_String id_to_smiles;
  IW_STL_Hash_Map_int id_to_ndx;

  if (!read_smiles(cl[0], id_to_ndx, id_to_smiles)) {
    cerr << "Cannot determine smiles and id's from '" << cl[0] << "'\n";
    return 1;
  }

  assert(id_to_smiles.size() == id_to_ndx.size());

  Leader_Results* l = new Leader_Results[n];
  std::unique_ptr<Leader_Results[]> free_l(l);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!l[i].build(cl[i], id_to_ndx)) {
      cerr << "Cannot read leader results from '" << cl[i] << "'\n";
      return i + 1;
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

  int rc = test_train_split_from_leader(argc, argv);

  return rc;
}
