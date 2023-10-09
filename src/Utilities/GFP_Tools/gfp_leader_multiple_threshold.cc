/*
  We have run gfp_nearneighbours and can quickly do several clusterings.
*/

#include <stdlib.h>
#include <values.h>

#include <fstream>
#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/set_or_unset.h"

#include "nn_results_.h"

using std::cerr;
using std::endl;

typedef float similarity_type_t;

/*
  We can have the threshold for each item read from the file
*/

static IWString score_tag;

/*
  Or the score may be embedded in the id
*/

static int score_column = -1;

static int verbose = 0;

/*
  The variables which control the clustering
*/

static int max_clusters_to_find = INT_MAX;
static int max_cluster_size = 0;

/*
  We can process any number of thresholds
*/

static int number_thresholds = 0;
static similarity_type_t* threshold = nullptr;

extending_resizable_array<int> cluster_size;

/*
  The identifier tag used in each TDT
*/

/*static IWString smiles_tag ("$SMI<");

static IWString identifier_tag ("PCN<");
static IWString neighbour_tag ("NBR<");

static IWString distance_tag ("DIST<");*/

static int number_threshold_tags = 0;
static IWString* threshold_from_file_tag = nullptr;

static IWString max_cluster_size_tag;
static int max_cluster_size_column = -1;

static int
fetch_tdt_value(const const_IWSubstring& buffer, IWString& zresult)
{
  assert(buffer.ends_with('>'));

  zresult = buffer;
  zresult.remove_up_to_first('<');
  zresult.chop();

  return 1;
}

template <typename T>
int
fetch_tdt_value(const const_IWSubstring& buffer, T& zresult)
{
  IWString tmp;
  if (!fetch_tdt_value(buffer, tmp)) {
    return 0;
  }

  if (!tmp.numeric_value(zresult)) {
    cerr << "Invalid numeric '" << tmp << "'\n";
    return 0;
  }

  return 1;
}

template <typename T>
int
fetch_tdt_value(const const_IWSubstring& buffer, Set_or_Unset<T>& zresult)
{
  T tmp;
  if (!fetch_tdt_value(buffer, tmp)) {
    return 0;
  }

  zresult.set(tmp);

  return 1;
}

#ifdef __GNUG__
template int
fetch_tdt_value(const const_IWSubstring&, float&);
template int
fetch_tdt_value(const const_IWSubstring&, int&);
template int
fetch_tdt_value(const const_IWSubstring&, Set_or_Unset<similarity_type_t>&);
template int
fetch_tdt_value(const const_IWSubstring&, Set_or_Unset<int>&);
#endif

/*
  When a neighbour is read in, all it knows is its index. Once the whole pool is read in,
  that will be converted to a pointer to a Leader_Item
*/

class Leader_Item;  // forward declaration

class Neighbour
{
 private:
  similarity_type_t _distance;

  int _ndx;

  Leader_Item* _myparent;

 public:
  Neighbour();

  const Leader_Item*
  parent() const
  {
    return _myparent;
  }

  void
  set_parent(Leader_Item* p)
  {
    _myparent = p;
  }

  int
  index_in_pool() const
  {
    return _ndx;
  }

  void
  set_index_in_pool(int s)
  {
    _ndx = s;
  }

  const IWString&
  smiles() const;
  const IWString&
  id() const;

  similarity_type_t
  distance() const
  {
    return _distance;
  }

  void
  set_distance(similarity_type_t s)
  {
    _distance = s;
  }

  int
  selected() const;
  void
  set_selected(int);
  int
  threshold(similarity_type_t&, int) const;

  int
  build(iwstring_data_source&, int&);
};

Neighbour::Neighbour()
{
  _distance = 2.0;

  _ndx = -1;

  _myparent = nullptr;

  return;
}

int
Neighbour::build(iwstring_data_source& input, int& fatal)
{
  const_IWSubstring buffer;
  if (!input.next_record(buffer))  // should never happen, should hit the vertical bar
  {
    fatal = 1;
    return 0;
  }

  if ('|' == buffer) {
    fatal = 0;
    return 0;
  }

  if (!buffer.starts_with(neighbour_tag)) {
    cerr << "Neighbour::build: Not '" << neighbour_tag << "' -> '" << buffer << "'\n";
    fatal = 1;
    return 0;
  }

  if (!fetch_tdt_value(buffer, _ndx) || _ndx < 0) {
    cerr << "Neighbour::build: invalid index '" << buffer << "'\n";
    fatal = 1;
    return 0;
  }

  if (!input.next_record(buffer)) {
    cerr << "Neighbour::build: eof\n";
    fatal = 1;
    return 0;
  }

  if (!buffer.starts_with(distance_tag)) {
    cerr << "Neighbour::build: Not '" << distance_tag << "' -> '" << buffer << "'\n";
    fatal = 1;
    return 0;
  }

  if (!fetch_tdt_value(buffer, _distance) || _distance < 0.0 || _distance > 1.0) {
    cerr << "Neighbour::build: invalid distance '" << buffer << "'\n";
    fatal = 1;
    return 0;
  }

  return 1;
}

typedef float score_t;

class Leader_Item : public NN_Item_Base<Neighbour>
{
 private:
  int _selected;

  //  Each item can limit the number of items in its cluster

  Set_or_Unset<int> _max_cluster_size;

  //  Each item may have any number of thresholds associated with it

  similarity_type_t* _threshold;

  score_t _score;

  //  private functions

  int
  _fetch_score_value_from_id();
  int
  _matches_threshold_from_file_tag(const const_IWSubstring& buffer);

 public:
  Leader_Item();
  ~Leader_Item();

  int
  selected() const
  {
    return _selected;
  }

  void
  set_selected(int s)
  {
    _selected = s;
  }

  int
  build(iwstring_data_source& input);

  score_t
  score() const
  {
    return _score;
  }

  int
  threshold(similarity_type_t& t, int) const;

  int
  max_cluster_size(int& m) const
  {
    return _max_cluster_size.value(m);
  }

  int
  establish_cross_reference_pointers(Leader_Item**);
};

int
Neighbour::selected() const
{
  return _myparent->selected();
}

int
Neighbour::threshold(similarity_type_t& t, int which_threshold) const
{
  return _myparent->threshold(t, which_threshold);
}

void
Neighbour::set_selected(int s)
{
  _myparent->set_selected(s);

  return;
}

const IWString&
Neighbour::smiles() const
{
  return _myparent->smiles();
}

const IWString&
Neighbour::id() const
{
  return _myparent->id();
}

#ifdef __GNUG__
template class resizable_array_p<Neighbour>;
template class resizable_array_base<Neighbour*>;
template class resizable_array<Neighbour*>;

template class NN_Results_Base<Leader_Item>;
template class NN_Item_Base<Neighbour>;
#endif

Leader_Item::Leader_Item()
{
  _selected = 0;
  _score = 0.0;

  if (number_threshold_tags) {
    _threshold = new similarity_type_t[number_threshold_tags];
    set_vector(_threshold, number_threshold_tags, static_cast<similarity_type_t>(-1.0));
  } else {
    _threshold = nullptr;
  }

  return;
}

Leader_Item::~Leader_Item()
{
  if (NULL != _threshold) {
    delete _threshold;
  }

  return;
}

int
Leader_Item::threshold(similarity_type_t& zresult, int which_one) const
{
  assert(which_one >= 0);

  similarity_type_t t = _threshold[which_one];

  if (t >= 0.0)  // valid value, therefore we have a threshold
  {
    zresult = t;
    return 1;
  }

  return 0;
}

int
Leader_Item::_matches_threshold_from_file_tag(const const_IWSubstring& buffer)
{
  for (int i = 0; i < number_threshold_tags; i++) {
    if (!buffer.starts_with(threshold_from_file_tag[i])) {
      continue;
    }

    if (!fetch_tdt_value(buffer, _threshold[i]) || _threshold[i] < 0.0 ||
        _threshold[i] >= 1.0) {
      cerr << "Leader_Item::_matches_threshold_from_file_tag: invalid threshold '"
           << buffer << "'\n";
      return 0;
    }

    return 1;
  }

  return 0;  // not recognised
}

int
Leader_Item::build(iwstring_data_source& input)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with(smiles_tag)) {
      fetch_tdt_value(buffer, _smiles);
    } else if (_matches_threshold_from_file_tag(buffer)) {
    } else if (score_tag.length() && buffer.starts_with(score_tag)) {
      if (!fetch_tdt_value(buffer, _score)) {
        return 0;
      }
    } else if (max_cluster_size_tag.length() &&
               buffer.starts_with(max_cluster_size_tag)) {
      if (!fetch_tdt_value(buffer, _max_cluster_size)) {
        return 0;
      }
    } else if (buffer.starts_with(neighbour_tag)) {
      input.push_record();
      break;
    } else if (buffer.starts_with(identifier_tag)) {
      if (_id.length())  // how could this happen?
      {
        cerr << "More than one '" << identifier_tag << "' found at line "
             << input.lines_read() << endl;
        return 0;
      }

      fetch_tdt_value(buffer, _id);
      if (score_column > 0) {
        if (!_fetch_score_value_from_id()) {
          return 0;
        }
      }
    } else if ('|' == buffer) {
      if (0 == _smiles.length() || 0 == _id.length()) {
        cerr << "Leader_Item::build: missing info at line " << input.lines_read() << endl;
        return 0;
      }

      return 1;
    }
  }

  // We have read the info for the target, not the neighbours

  while (1) {
    Neighbour* tmp = new Neighbour;

    if (tmp->build(input, fatal)) {
      add(tmp);
      continue;
    }

    //  All the various failures

    delete tmp;

    if (input.eof())  // no real problem
    {
      fatal = 0;
      return 0;
    }

    if (!fatal) {  // probably just got to end
      return 1;
    }

    cerr << "Fatal error fetching nbr\n";
    return 0;
  }

  cerr << _id << " has " << number_elements() << " neighbours\n";
  return 1;
}

int
Leader_Item::_fetch_score_value_from_id()
{
  assert(score_column >= 0);

  const_IWSubstring s;
  if (!_id.word(score_column, s)) {
    cerr << "Leader_Item::_fetch_score_value_from_id: cannot access column "
         << score_column << " in '" << _id << "'\n";
    return 0;
  }

  if (!s.numeric_value(_score)) {
    cerr << "Leader_Item::_fetch_score_value_from_id: invalid score value, column "
         << score_column << " in '" << _id << "'\n";
    return 0;
  }

  return 1;
}

static Leader_Item* pool = nullptr;

static int pool_size = 0;

int
Leader_Item::establish_cross_reference_pointers(Leader_Item** xref)
{
  for (int i = 0; i < _number_elements; i++) {
    int ndx = _things[i]->index_in_pool();
    assert(ndx >= 0 && ndx < pool_size);

    Leader_Item* p = xref[ndx];

    _things[i]->set_parent(p);
  }

  return 1;
}

static void
reset_pool()
{
  for (int i = 0; i < pool_size; i++) {
    Leader_Item& l = pool[i];
    l.set_selected(0);
  }

  cluster_size.resize_keep_storage(0);

  return;
}

/*
  A cluster is a set of pointers to such objects
*/

typedef resizable_array<Neighbour*> Cluster;

int
distance_comparitor(Neighbour* const* ppfp1, Neighbour* const* ppfp2)
{
  const Neighbour* pfp1 = *ppfp1;
  const Neighbour* pfp2 = *ppfp2;

  if (pfp1->distance() < pfp2->distance()) {
    return -1;
  } else if (pfp1->distance() > pfp2->distance()) {
    return 1;
  } else {
    return 0;
  }
}

static int
report_clustering(int clusters_found)
{
  cerr << "Clustered " << pool_size << " fingerprints into " << clusters_found
       << " clusters\n";
  int isum = 0;
  for (int i = 0; i < cluster_size.number_elements(); i++) {
    int j = cluster_size[i];
    if (0 == j) {
      continue;
    }

    cerr << j << " clusters were of size " << i << " members\n";

    isum += j * i;
  }

  cerr << "In clusters " << isum << endl;

  return cerr.good();
}

static int
write_cluster(const Leader_Item& ldr, const Cluster& cluster, int threshold_number,
              int clusters_found, int items_selected, ostream& output)
{
  int cs = cluster.number_elements();
  cluster_size[cs + 1]++;  // the leader isn't in the cluster

  if (verbose > 1) {
    cerr << "Cluster " << clusters_found << ' ' << (cs + 1) << " items, centre '"
         << ldr.id() << "', ";
    if (threshold_number >= 0) {
      similarity_type_t threshold;
      (void)ldr.threshold(threshold, threshold_number);
      cerr << "threshold " << threshold << ", ";
    }
    cerr << (items_selected + cs + 1) << " items selected\n";
  }

  output << smiles_tag << ldr.smiles() << ">\n";
  output << identifier_tag << ldr.id() << ">\n";

  output << "CLUSTER<" << clusters_found << ">\n";
  output << "CSIZE<" << (cs + 1) << ">\n";

  for (int i = 0; i < cs && output.good(); i++) {
    const Neighbour* n = cluster[i];

    output << smiles_tag << n->smiles() << ">\n";
    output << identifier_tag << n->id() << ">\n";
    output << distance_tag << n->distance() << ">\n";
  }

  output << "|\n";

  return output.good();
}

static int
choose_next_centre(int& icentre)
{
  icentre = -1;
  if (0 == score_tag.length())  // just grab the first unselected item
  {
    for (int i = 0; i < pool_size; i++) {
      if (!pool[i].selected()) {
        icentre = i;
        return 1;
      }
    }
  } else  // raw scores
  {
    score_t max_score = 0.0;
    for (int i = 0; i < pool_size; i++) {
      if (pool[i].selected()) {
        continue;
      }

      score_t s = pool[i].score();

      if (icentre < 0 || s > max_score) {
        max_score = s;
        icentre = i;
      }
    }
  }

  return icentre >= 0;
}

static int
form_cluster_threshold(int icentre, similarity_type_t my_threshold, Cluster& cluster)
{
  Leader_Item& ldr = pool[icentre];

  int number_neighbours = ldr.number_elements();

  for (int i = 0; i < number_neighbours; i++) {
    Neighbour* n = ldr[i];

    if (n->selected()) {
      continue;
    }

    if (n->distance() > my_threshold) {
      continue;
    }

    cluster.add(n);
    n->set_selected(1);
  }

  return cluster.number_elements();
}

static int
form_cluster_max_cluster_size(int icentre, Cluster& cluster,
                              similarity_type_t my_threshold,
                              int max_cluster_size_this_molecule)
{
  assert(max_cluster_size_this_molecule > 0);

  Leader_Item& ldr = pool[icentre];

  int number_neighbours = ldr.number_elements();

  cluster.resize(number_neighbours);

  for (int i = 0; i < number_neighbours; i++) {
    Neighbour* n = ldr[i];

    if (n->selected()) {
      continue;
    }

    if (n->distance() > my_threshold) {
      continue;
    }

    cluster.add(n);
    n->set_selected(1);

    if (cluster.number_elements() >= max_cluster_size_this_molecule) {
      return 1;
    }
  }

  return 1;
}

int
leader(similarity_type_t global_threshold, int threshold_number, ostream& output)
{
  if (verbose > 1) {
    cerr << "Beginning clustering at " << global_threshold << " threshold "
         << threshold_number << endl;
  }

  int clusters_found = 0;
  int items_selected = 0;

  int icentre;
  if (!choose_next_centre(icentre)) {
    cerr << "Yipes, cannot find initial leader\n";
    return 0;
  }

  Cluster cluster;

  while (items_selected < pool_size) {
    Leader_Item& ldr = pool[icentre];

    similarity_type_t my_threshold;
    if (threshold_number < 0) {
      my_threshold = global_threshold;
    } else if (ldr.threshold(my_threshold, threshold_number)) {  // has come from the file
      ;
    } else {
      my_threshold = global_threshold;
    }

    int max_cluster_size_this_molecule;

    if (ldr.max_cluster_size(max_cluster_size_this_molecule)) {
      ;
    } else {
      max_cluster_size_this_molecule = max_cluster_size;
    }

    if (verbose > 1) {
      cerr << "Start cluster " << clusters_found << ". Threshold = " << my_threshold
           << ", max size " << max_cluster_size_this_molecule << endl;
    }

    ldr.set_selected(1);

    if (max_cluster_size_this_molecule) {
      form_cluster_max_cluster_size(icentre, cluster, my_threshold,
                                    max_cluster_size_this_molecule);
    } else {
      form_cluster_threshold(icentre, my_threshold, cluster);
    }

    (void)write_cluster(ldr, cluster, threshold_number, clusters_found, items_selected,
                        output);

    clusters_found++;
    if (clusters_found >= max_clusters_to_find) {
      break;
    }

    items_selected += cluster.number_elements() + 1;

    if (!choose_next_centre(icentre)) {
      break;
    }

    cluster.resize_keep_storage(0);
  }

  if (verbose) {
    report_clustering(clusters_found);
  }

  return output.good();
}

static int
leader(similarity_type_t global_threshold, int threshold_number,
       const const_IWSubstring& stem, int ndx)
{
  IWString fname(stem);
  fname << ndx;
  fname << ".ldr";

  std::ofstream output(fname.null_terminated_chars(), std::ios::out);
  if (!output.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  static int already_been_here = 0;

  if (already_been_here) {
    reset_pool();
  }

  already_been_here = 1;

  if (verbose) {
    cerr << "Begin clustering at " << global_threshold << endl;
  }

  return leader(global_threshold, threshold_number, output);
}

static int
establish_cross_reference_pointers()
{
  if (verbose > 1) {
    cerr << "Establishing cross reference pointers for " << pool_size << " items\n";
  }

  Leader_Item** tmp = new Leader_Item*[pool_size];
  std::unique_ptr<Leader_Item*[]> free_tmp(tmp);

  for (int i = 0; i < pool_size; i++) {
    tmp[i] = &pool[i];
  }

  for (int i = 0; i < pool_size; i++) {
    pool[i].establish_cross_reference_pointers(tmp);
  }

  return 1;
}

static int
do_read_nn_results(iwstring_data_source& input, int& items_in_pool)
{
  while (1) {
    Leader_Item& pi = pool[items_in_pool];
    int fatal;
    if (pi.build(input, fatal)) {
      ;
    } else if (fatal) {
      cerr << "Cannot read nn results before line " << input.lines_read() << endl;
      return 0;
    } else {
      return items_in_pool;
    }

    items_in_pool++;
    if (items_in_pool >= pool_size) {
      break;
    }
  }

  return items_in_pool;
}

static int
read_nn_results(iwstring_data_source& input, int& items_in_pool)
{
  if (0 == pool_size) {
    IWString tmp("^\\|");

    std::unique_ptr<re2::RE2> pcn = std::make_unique<RE2>(tmp);
    pool_size = input.grep(*pcn);

    if (0 == pool_size) {
      cerr << "No occurrences of " << pcn->pattern() << "' in input\n";
      return 0;
    }

    cerr << "Pool automatically sized to " << pool_size << endl;

    pool = new Leader_Item[pool_size];
    if (NULL == pool) {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 0;
    }
  }

  int rc = do_read_nn_results(input, items_in_pool);

  return rc;
}

static int
read_nn_results(const char* fname, int& items_in_pool)
{
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    cerr << "Reading nn results from '" << fname << "'\n";
  }

  if (!read_nn_results(input, items_in_pool)) {
    cerr << "Yipes, cannot read NN results from '" << fname << "'\n";
    return 0;
  }

  if (verbose) {
    cerr << "Read " << items_in_pool << " NN results\n";
  }

  return 1;
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
  cerr << "Performs leader clustering on a set of fingerprints\n";
  cerr << "Usage <options> <input_file>\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -t <dis>         specify distance threshold\n";
  cerr << " -C <number>      maximum number of clusters to find\n";
  cerr << " -H <TAG>         threshold for each molecule in dataitem <TAG>\n";
  cerr << " -m <number>      maximum cluster size\n";
  cerr << " -M <tag>         max cluster size for each molecule in <TAG>\n";
  cerr << " -M col=<col>     max cluster size for each molecule in column <col> of name\n";
  cerr << " -I <TAG>         specify identifier tag\n";
  cerr << " -S <TAG>         score tag\n";
  cerr << " -S col=<col>     scores are in column <col> of the name\n";
  cerr << " -B <stem>        file name stem for multiple thresholds\n";
  cerr << " -v               verbose output\n";
  // clang-format on

  exit(rc);
}

static int
leader(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vs:I:E:t:F:P:W:H:S:C:m:M:B:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 1) {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }

    pool = new Leader_Item[pool_size];

    if (NULL == pool) {
      cerr << "Cannot allocate space for system of size " << pool_size << endl;
      return 11;
    }

    if (verbose) {
      cerr << "system sized to " << pool_size << endl;
    }
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', max_clusters_to_find) || max_clusters_to_find < 1) {
      cerr << "The -C option (max clusters to find) must be followed by a positive "
              "integer\n";
      usage(41);
    }

    if (verbose) {
      cerr << "Will find a max of " << max_clusters_to_find << " clusters\n";
    }
  }

  // We need to be careful with the -i and -I options. Remember
  // that the pool is built first

  if (cl.option_present('I')) {
    (void)cl.value('I', identifier_tag);

    if (verbose) {
      cerr << "Identifiers in dataitem '" << identifier_tag << "'\n";
    }
  }

  if (cl.option_present('H')) {
    number_threshold_tags = cl.option_count('H');
    threshold_from_file_tag = new IWString[number_threshold_tags];

    for (int i = 0; i < number_threshold_tags; i++) {
      threshold_from_file_tag[i] = cl.string_value('H', i);

      if (verbose) {
        cerr << "Threshold in tag '" << threshold_from_file_tag[i] << "' tag\n";
      }
    }
  }

  if (cl.option_present('S')) {
    const_IWSubstring s = cl.string_value('S');

    if (s.starts_with("col=")) {
      s.remove_leading_chars(4);
      if (!s.numeric_value(score_column) || score_column < 1) {
        cerr << "INvalid score column specifier 'col=" << s << "'\n";
        return 8;
      }

      if (verbose) {
        cerr << "Scores taken from column " << score_column << " of the identifier\n";
      }

      score_column--;
    } else {
      score_tag = cl.string_value('S');
      if (verbose) {
        cerr << "Score tag is " << score_tag << "'\n";
      }
    }
  }

  if (!cl.option_present('t') && !cl.option_present('H') && !cl.option_present('Y')) {
    cerr << "Threshold distance must be specified via -t, -H or -V options\n";
    usage(28);
  }

  number_thresholds = cl.option_count('t');
  if (number_thresholds) {
    threshold = new similarity_type_t[number_thresholds];
    for (int i = 0; i < number_thresholds; i++) {
      if (!cl.value('t', threshold[i], i) || threshold[i] < 0.0 || threshold[i] > 1.0) {
        cerr << "Invalid threshold\n";
        return 6;
      }

      if (verbose) {
        cerr << "Threshold " << i << " set to " << threshold[i] << endl;
      }
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', max_cluster_size) || max_cluster_size < 2) {
      cerr << "The -m (max cluster size) option must be followed by a whole number > 1\n";
      usage(43);
    }

    if (verbose) {
      cerr << "Max cluster size " << max_cluster_size << endl;
    }
  }

  if (cl.option_present('M')) {
    const_IWSubstring m = cl.string_value('M');

    if (m.starts_with("col=")) {
      m.remove_leading_chars(4);
      if (!m.numeric_value(max_cluster_size_column) || max_cluster_size_column < 1) {
        cerr << "The column for the per molecule maximum cluster size must be a whole "
                "positive number\n";
        usage(11);
      }

      if (verbose) {
        cerr << "The maximum cluster size per molecule will be in column "
             << max_cluster_size_column << endl;
      }

      max_cluster_size_column--;
    } else {
      max_cluster_size_tag = m;

      if (verbose) {
        cerr << "Max cluster size in '" << max_cluster_size_tag << "' tag\n";
      }
    }
  }

  int items_in_pool = 0;

  for (int i = 0; i < cl.number_elements(); i++) {
    if (!read_nn_results(cl[i], items_in_pool)) {
      cerr << "Cannot build pool from '" << cl[i] << "'\n";
      return 21;
    }
  }

  pool_size = items_in_pool;

  if (!establish_cross_reference_pointers()) {
    cerr << "Cannot establish cross reference pointers\n";
    return 15;
  }

  int number_clusterings = number_thresholds + number_threshold_tags;
  if (number_clusterings > 1 && !cl.option_present('B')) {
    cerr << "If doing multiple thresholds, must specify file name stem via the -B "
            "option\n";
    usage(12);
  }

  int computation_number = 0;

  if (1 == number_thresholds && 0 == number_threshold_tags) {
    leader(threshold[0], -1, std::cout);
  } else {
    const_IWSubstring b = cl.string_value('B');

    for (int i = 0; i < number_thresholds; i++) {
      if (!leader(threshold[i], -1, b, computation_number)) {
        cerr << "Clustering at " << threshold[i] << " failed\n";
        return i + 1;
      }

      computation_number++;
    }
  }

  if (1 == number_threshold_tags && 0 == number_thresholds) {
    leader(0.0, 0, std::cout);
  } else {
    const_IWSubstring b = cl.string_value('B');

    for (int i = 0; i < number_threshold_tags; i++) {
      if (!leader(0.0, i, b, computation_number)) {
        cerr << "Clustering on threshold tag " << i << " failed\n";
        return i + 1;
      }

      computation_number++;
    }
  }

  return 0;
}

int
main(int argc, char** argv)
{
  int rc = leader(argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status(stderr);
#endif

  return rc;
}
