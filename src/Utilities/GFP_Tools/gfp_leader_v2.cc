/*
  An implementation of the leader algorithm for fingerprints
  This variant produces output that can be processed by nplotnn
*/

#include <stdint.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <limits>

#include "google/protobuf/text_format.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"
#include "leader.h"
#include "sparse_collection.h"
#include "tversky.h"

#include "Utilities/GFP_Tools/nearneighbours.pb.h"

using std::cerr;

static Tversky tversky;

// If there is a per-item score. Can be specified as either
//  TDT dataitem value
//  Numeric value in a column in the name.

static IWString score_tag;
static int score_column = -1;

/*
  May 99. When looking for the next molecule to select, we can
  add a multiple of the distance to the nearest cluster centre
  to the score. That way, we can include some varying function
  of diversity into the next cluster selection
*/

static similarity_type_t cluster_distance_scale_factor = 0.0;

static resizable_array_p<IWString> dataitems_to_echo;

static int verbose = 0;

static similarity_type_t abandon_distance_cutoff = -1.0;

/*
  The variables which control the clustering
*/

static int max_clusters_to_find = std::numeric_limits<int>::max();
static int clusters_found = 0;
static int max_cluster_size = 0;

static similarity_type_t threshold = 0.0;

static int threshold_column = -1;

/*
  when dealing with clusters which are decided by the threshold, we can
  optinally sort the cluster members by their distance from the leader
*/

static int sort_by_distance_from_centre = 0;

static int items_selected = 0;

extending_resizable_array<int> cluster_size;

static Accumulator<float> distance_stats;

static std::ofstream stream_for_discarded_by_previously_selected;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static IWString distance_tag("DIST<");

/*
  Mar 2004. in order to
  sample an SAR within each cluster. Note that the way this is implemented
  now is somewhat strange because rather than marking the discarded molecules
  in some way, they simply do not appear in the output! Change if this ever
  becomes a problem.

  There are two ways of doing the within-cluster clustering.
  A fixed threshold applied to all clusters
  A constant multiple of the threshold used to form the cluster
*/

static similarity_type_t sub_cluster_threshold = static_cast<similarity_type_t>(0.0);

static float sub_cluster_threshold_ratio = static_cast<float>(0.0);

/*
  A single variable that indicates whether or not any form of cluster clustering
  is active
*/

static int leader_on_clusters = 0;

// Sept 2023. Add proto output as an option.
static int write_as_proto = 0;

// When writing textproto values a zero distance will not be written.
// We have the option of writing a small value instead of zero.
static int offset_zero_distances_in_textproto = 0;
// The actual distance used.
static constexpr float kNotZero = 0.00001f;

GFP_L::GFP_L() {
  _selected = 0;

  _score = 0.0;

  _shortest_distance_to_cluster_centre = static_cast<similarity_type_t>(1.0);

  return;
}

static IWString threshold_from_file_tag;

static IWString max_cluster_size_tag;

static int max_cluster_size_column = -1;

int
GFP_L::construct_from_tdt(IW_TDT &tdt, int &fatal) {
  if (!IW_General_Fingerprint::construct_from_tdt(tdt, fatal)) {
    return 0;
  }

  if (threshold_from_file_tag.length()) {
    similarity_type_t tmp;
    if (!tdt.dataitem_value(threshold_from_file_tag, tmp) || tmp < 0.0) {
      cerr << "GFP_L::construct_from_tdt: invalid '" << threshold_from_file_tag
           << "' in tdt\n";
      return 0;
    }

    _threshold.set(tmp);
  } else if (threshold_column >= 0) {
    const_IWSubstring t;
    if (!_id.word(threshold_column, t)) {
      cerr << "GFP_L::construct_from_tdt: no " << threshold_column << " column in '"
           << _id << "'";
      if (::threshold > 0.0) {
        cerr << " using default threshold\n";
      } else {
        cerr << '\n';
        return 0;
      }
    } else {
      similarity_type_t d;
      if (!t.numeric_value(d) || d < 0.0 || d > 1.0) {
        cerr << "Invalid threshold '" << t << "' in '" << _id << "'\n";
        return 0;
      }

      _threshold.set(d);
    }
  }

  if (score_tag.length()) {
    if (!tdt.dataitem_value(score_tag, _score)) {
      cerr << "GFP_L::construct_from_tdt: cannot extract '" << score_tag
           << "' from tdt\n";
      return 0;
    }
  }

  if (score_column >= 0) {
    const_IWSubstring c;
    if (!_id.word(score_column, c)) {
      cerr << "Cannot extract column " << score_column << " from '" << _id << "'\n";
      return 0;
    }

    if (!c.numeric_value(_score)) {
      cerr << "Invalid score, column " << score_column << " in '" << _id << "'\n";
      return 0;
    }

    if (verbose > 2) {
      cerr << _id << " set score to " << _score << '\n';
    }
  }

  if (max_cluster_size_tag.length() > 0) {
    int tmp;
    if (!tdt.dataitem_value(max_cluster_size_tag, tmp) || tmp < 1) {
      cerr << "GFP_L::construct_from_tdt: missing or invalid '" << max_cluster_size_tag
           << "' in tdt\n";
      return 0;
    }

    _max_cluster_size.set(tmp);
  }

  if (max_cluster_size_column >= 0) {
    const_IWSubstring c;
    if (!_id.word(max_cluster_size_column, c)) {
      cerr << "Cannot extract column " << max_cluster_size_column << " from '" << _id
           << "'\n";
      return 0;
    }

    int tmp;
    if (!c.numeric_value(tmp) || tmp < 1) {
      cerr << "Invalid maximum cluster size, column " << max_cluster_size_column
           << " in '" << _id << "'\n";
      return 0;
    }

    //  cerr << "Max cluster size for '" << _id << "' is " << tmp << '\n';
    _max_cluster_size.set(tmp);
  }

  return 1;
}

/*
  Our pool is an array of FP objects
*/

static GFP_L *pool = nullptr;

static int pool_size = 0;

/*
  A cluster is a set of pointers to such objects
*/

typedef resizable_array<GFP_L *> Cluster;

/*
  Within the pool it is convenient to keep track of the first unselected
  item
*/

static int first_unselected = 0;
static int last_unselected = 0;

static int
echo_selected_dataitems(IW_TDT &tdt, const resizable_array_p<IWString> &items_to_echo,
                        IWString_and_File_Descriptor &output) {
  int ne = items_to_echo.number_elements();
  for (int i = 0; i < ne; i++) {
    const IWString &tag = *(items_to_echo[i]);

    if (!tdt.echo_dataitem(tag, 0, output)) {
      cerr << "Cannot echo '" << tag << "' from TDT\n";
      throw "Missing dataitem";
      return 0;
    }
  }

  return output.good();
}

static int
echo_all_dataitems(IW_TDT &tdt, IWString_and_File_Descriptor &output) {
  return tdt.write_all_except_vbar(output);
}

static int
write_cluster_data(IW_TDT &tdt, int clusters_found, int id_within_cluster,
                   similarity_type_t distance_to_centre,
                   IWString_and_File_Descriptor &output) {
  if (dataitems_to_echo.number_elements()) {
    echo_selected_dataitems(tdt, dataitems_to_echo, output);
  } else {
    echo_all_dataitems(tdt, output);
  }

  if (distance_to_centre >= 0.0) {
    output << distance_tag << distance_to_centre << ">\n";
  }

  if (!output.good()) {
    throw "Bad output stream";
  }

  return output.good();
}

static int
get_tdt(IW_TDT &tdt, iwstring_data_source &input, const GFP_L &fp) {
  off_t offset;
  (void)fp.offset(offset);

  if (!input.seekg(offset)) {
    cerr << "Cannot seek to offset '" << offset << '\n';
    return 0;
  }

  return tdt.next(input);
}

// We may write a small floating point number for zero distances just to 
// ensure consistent tabular output in textproto form.
static int
WriteAsProto(const Cluster& cluster,
             int cluster_number,
             const IW_TDT& tdt,
             iwstring_data_source& input,
             IWString_and_File_Descriptor& output) {
  nnbr::NearNeighbours proto;
  IWString value;
  tdt.dataitem_value(smiles_tag, value);
  proto.set_smiles(value.data(), value.length());
  tdt.dataitem_value(identifier_tag, value);
  proto.set_name(value.data(), value.length());
  proto.set_csize(cluster.size());
  proto.set_cluster(cluster_number);

  for (const GFP_L* c : cluster) {
    IW_TDT tdt;
    if (! get_tdt(tdt, input, *c)) {
      continue;
    }

    nnbr::Nbr* nbr = proto.add_nbr();

    tdt.dataitem_value(smiles_tag, value);
    nbr->set_smi(value.data(), value.length());
    tdt.dataitem_value(identifier_tag, value);
    nbr->set_id(value.data(), value.length());
    if (c->distance() > 0.0f) {
      nbr->set_dist(c->distance());
    } else if (offset_zero_distances_in_textproto) {
      nbr->set_dist(kNotZero);
    }
  }

  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;
  if (!printer.PrintToString(proto, &buffer)) {
    cerr << "WriteAsProto:cannot print " << proto.ShortDebugString() << '\n';
    return 0;
  }

  output << buffer << '\n';
  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

static int
do_leader_on_clusters(Cluster &cluster, similarity_type_t threshold) {
  similarity_type_t my_threshold;
  if (sub_cluster_threshold > static_cast<similarity_type_t>(0.0)) {
    my_threshold = sub_cluster_threshold;
  } else if (sub_cluster_threshold_ratio > static_cast<float>(0.0)) {
    my_threshold = threshold * sub_cluster_threshold_ratio;
  } else {
    cerr << "Not sure how to sub cluster\n";
    return 0;
  }

  for (int i = 0; i < cluster.number_elements(); i++) {
    GFP_L *ldr = cluster[i];

    for (int j = i + 1; j < cluster.number_elements(); j++) {
      GFP_L *cj = cluster[j];

      similarity_type_t d = ldr->IW_General_Fingerprint::distance(*cj);
      if (d <= my_threshold) {
        cluster.remove_item(j);
        j--;
      }
    }
  }

  return cluster.number_elements();
}

int
distance_comparitor(GFP_L *const *ppfp1, GFP_L *const *ppfp2) {
  const GFP_L *pfp1 = *ppfp1;
  const GFP_L *pfp2 = *ppfp2;

  if (pfp1->distance() < pfp2->distance()) {
    return -1;
  } else if (pfp1->distance() > pfp2->distance()) {
    return 1;
  } else {
    return 0;
  }
}

class Distance_Comparitor {
 private:
 public:
  int operator()(GFP_L *const, GFP_L *const);
};

int
Distance_Comparitor::operator()(GFP_L *const p1, GFP_L *const p2) {
  if (p1->distance() < p2->distance()) {
    return -1;
  }

  if (p1->distance() > p2->distance()) {
    return 1;
  }

  return 0;
}

template void resizable_array_base<GFP_L *>::iwqsort<Distance_Comparitor>(
    Distance_Comparitor &);
template void iwqsort<GFP_L *, Distance_Comparitor>(GFP_L **, int, Distance_Comparitor &);
template void iwqsort<GFP_L *, Distance_Comparitor>(GFP_L **, int, Distance_Comparitor &,
                                                    void *);
template void compare_two_items<GFP_L *, Distance_Comparitor>(GFP_L **,
                                                              Distance_Comparitor &,
                                                              void *);
template void move_in_from_left<GFP_L *, Distance_Comparitor>(GFP_L **, int &, int &, int,
                                                              Distance_Comparitor &,
                                                              void *);
// template void move_in_from_right<GFP_L, Distance_Comparitor>(GFP_L**, int&, int&,
// Distance_Comparitor&);
template void swap_elements<GFP_L *>(GFP_L *&, GFP_L *&, void *);
template void move_in_from_right<GFP_L *, Distance_Comparitor>(GFP_L **, int &, int &,
                                                               Distance_Comparitor &);

static int
process_cluster(Cluster &cluster, similarity_type_t my_threshold,
                iwstring_data_source &input, IWString_and_File_Descriptor &output) {
  if (sort_by_distance_from_centre) {
    Distance_Comparitor dc;
    cluster.iwqsort(dc);
  }

  if (leader_on_clusters) {
    do_leader_on_clusters(cluster, my_threshold);
  }

  int cs = cluster.number_elements();
  cluster_size[cs]++;  // the leader isn't in the cluster

  GFP_L *centre = cluster[0];

  if (verbose) {
    cerr << "Cluster " << clusters_found << ' ' << cs << " items, centre '"
         << cluster[0]->id() << "', ";
    if (threshold_from_file_tag.length()) {
      similarity_type_t threshold = 0.0f;
      (void)centre->threshold(threshold);
      cerr << "threshold " << threshold << ", ";
    }
    cerr << (items_selected + cs) << " items selected\n";
  }

  IW_TDT tdt;
  if (!get_tdt(tdt, input, *centre)) {
    return 0;
  }

  if (write_as_proto) {
    return WriteAsProto(cluster, clusters_found, tdt, input, output);
  }

  if (dataitems_to_echo.number_elements()) {
    echo_selected_dataitems(tdt, dataitems_to_echo, output);
  } else {
    echo_all_dataitems(tdt, output);
  }

  output << "CLUSTER<" << clusters_found << ">\n";
  output << "CSIZE<" << cs << ">\n";

  // start at 1, we've already done centre above
  for (int i = 1; i < cs && output.good(); i++) {
    GFP_L &fp = *(cluster[i]);

    if (!get_tdt(tdt, input, fp)) {
      return 0;
    }

    if (!write_cluster_data(tdt, clusters_found, i, fp.distance(), output)) {
      return 0;
    }
  }

  output << "|\n";

  return output.good();
}

static int
choose_next_centre(int &icentre) {
  icentre = -1;

  // just grab the first unselected item
  if (0 == score_tag.length() && score_column < 0) {
    for (int i = 0; i < pool_size; i++) {
      if (!pool[i].selected()) {
        icentre = i;
        return 1;
      }
    }
  } else if (cluster_distance_scale_factor > static_cast<similarity_type_t>(0.0)) {
    score_t max_score = static_cast<score_t>(0.0);
    for (int i = 0; i < pool_size; i++) {
      if (pool[i].selected()) {
        continue;
      }

      score_t s = pool[i].score() + cluster_distance_scale_factor *
                                        pool[i].shortest_distance_to_cluster_centre();

      if (icentre < 0 || s > max_score) {
        max_score = s;
        icentre = i;
      }
    }
  } else  // raw scores
  {
    score_t max_score = static_cast<score_t>(0.0);
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
compute_the_distance(IW_General_Fingerprint &fp, IW_General_Fingerprint &p,
                     similarity_type_t &d) {
  if (!can_be_compared(fp, p)) {
    return 0;
  }

  if (tversky.active()) {
    d = static_cast<similarity_type_t>(1.0) -
        fp.IW_General_Fingerprint::tversky(p, tversky);
    return 1;
  }

  if (abandon_distance_cutoff > static_cast<similarity_type_t>(0.0)) {
    if (!fp.IW_General_Fingerprint::tanimoto(p, abandon_distance_cutoff, d)) {
      return 0;
    }

    d = static_cast<similarity_type_t>(1.0) - d;
    return 1;
  }

  d = fp.IW_General_Fingerprint::distance(p);

  return 1;
}

static int
form_cluster_threshold(int icentre, Cluster &cluster,
                       const similarity_type_t my_threshold) {
  GFP_L &fp = pool[icentre];

  if (verbose > 2) {
    cerr << "Leader is " << fp.id() << '\n';
  }

  int istop = last_unselected;

  int next_first_unselected = -1;

  for (int i = first_unselected; i <= istop; i++) {
    GFP_L &p = pool[i];

    if (p.selected()) {
      continue;
    }

    if (next_first_unselected < 0) {
      next_first_unselected = i;
    }

    last_unselected = i;

    similarity_type_t d;
    if (!compute_the_distance(fp, p, d)) {
      continue;
    }

    if (d <= my_threshold) {
      cluster.add(&(pool[i]));
      p.set_selected(1);
      p.set_distance(d);
      if (i == next_first_unselected) {
        next_first_unselected = -1;
      }
    } else if (d < p.shortest_distance_to_cluster_centre()) {
      p.set_shortest_distance_to_cluster_centre(d);
    }
  }

  if (next_first_unselected > first_unselected) {
    first_unselected = next_first_unselected;
  }

  return cluster.number_elements();
}

static int
form_cluster_max_cluster_size(int icentre, Cluster &cluster,
                              similarity_type_t my_threshold,
                              int max_cluster_size_this_molecule) {
  assert(max_cluster_size_this_molecule > 0);

  cluster.resize(pool_size);

  GFP_L &fp = pool[icentre];

  for (int i = 0; i < pool_size; i++) {
    GFP_L &p = pool[i];

    if (p.selected()) {
      continue;
    }

    if (!can_be_compared(fp, p)) {
      continue;
    }

    similarity_type_t d;
    if (!compute_the_distance(fp, p, d)) {
      continue;
    }

    if (my_threshold > static_cast<similarity_type_t>(0.0) && d > my_threshold) {
      continue;
    }

    cluster.add(&(pool[i]));
    p.set_distance(d);
    if (d < p.shortest_distance_to_cluster_centre()) {
      p.set_shortest_distance_to_cluster_centre(d);
    }
  }

  cluster.sort(&distance_comparitor);

  cluster.resize_keep_storage(max_cluster_size_this_molecule);

  int istop;
  if (cluster.number_elements() < max_cluster_size_this_molecule) {
    istop = cluster.number_elements();
  } else {
    istop = max_cluster_size_this_molecule;
  }

  for (int i = 0; i < istop; i++) {
    GFP_L *p = cluster[i];

    p->set_selected(1);
  }

  return 1;
}

/*
  The clustering will be limited either by the maximum number of items which
  can be in a cluster, or a threshold
*/

static int
form_cluster(int icentre, Cluster &cluster, const similarity_type_t my_threshold,
             int max_cluster_size_this_molecule) {
  cluster.resize_keep_storage(0);

  cluster.add(&(pool[icentre]));

  pool[icentre].selected() = 1;
  pool[icentre].set_distance(static_cast<similarity_type_t>(0.0));

  if (max_cluster_size_this_molecule) {
    return form_cluster_max_cluster_size(icentre, cluster, my_threshold,
                                         max_cluster_size_this_molecule);
  } else {
    return form_cluster_threshold(icentre, cluster, my_threshold);
  }
}

int
leader(iwstring_data_source &input, IWString_and_File_Descriptor &output) {
  assert(pool_size > 1);

  assert(0 == items_selected);
  assert(0 == clusters_found);

  first_unselected = 0;
  last_unselected = pool_size - 1;

  int icentre;
  if (!choose_next_centre(icentre)) {
    cerr << "Yipes, cannot find initial leader\n";
    return 0;
  }

  Cluster cluster;
  if (!cluster.resize(pool_size)) {
    cerr << "Yipes, cannot allocate " << pool_size << " elements in pool\n";
    return 0;
  }

  while (items_selected < pool_size) {
    GFP_L &centre = pool[icentre];

    similarity_type_t my_threshold = 0.0f;
    if (centre.threshold(my_threshold)) {  // has come from the file
      ;
    } else {
      my_threshold = threshold;
    }

    int max_cluster_size_this_molecule;

    if (max_cluster_size_tag.length() &&
        pool[icentre].max_cluster_size(max_cluster_size_this_molecule)) {
      ;
    } else if (max_cluster_size_column >= 0 &&
               pool[icentre].max_cluster_size(max_cluster_size_this_molecule)) {
      ;
    } else {
      max_cluster_size_this_molecule = max_cluster_size;
    }

    if (verbose > 1) {
      cerr << "Start cluster " << clusters_found << ". ndx " << icentre
           << ", threshold = " << my_threshold;
      if (max_cluster_size_this_molecule > 0) {
        cerr << ", max size " << max_cluster_size_this_molecule;
      }
      cerr << '\n';
    }

    (void)form_cluster(icentre, cluster, my_threshold, max_cluster_size_this_molecule);

    (void)process_cluster(cluster, my_threshold, input, output);

    clusters_found++;
    if (clusters_found >= max_clusters_to_find) {
      break;
    }

    items_selected += cluster.number_elements();

    if (!choose_next_centre(icentre)) {
      break;
    }
  }

  return 1;
}

static int
build_pool(iwstring_data_source &input) {
  off_t offset = input.tellg();

  int items_in_pool = 0;

  int tdts_read = 0;

  IW_TDT tdt;
  while (tdt.next(input)) {
    tdts_read++;

    int fatal;
    if (!pool[items_in_pool].construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      offset = input.tellg();
      continue;
    }

    pool[items_in_pool].set_offset(offset);

    items_in_pool++;

    if (items_in_pool == pool_size) {
      cerr << "Pool is full, max " << pool_size << '\n';
      break;
    }

    offset = input.tellg();
  }

  pool_size = items_in_pool;

  if (verbose) {
    cerr << "Read " << tdts_read << " TDT's, pool contains " << pool_size
         << " fingerprints\n";
  }

  return 1;
}

static int
build_pool(const const_IWSubstring &fname, iwstring_data_source &input) {
  IWString tmp(fname);

  if (!input.open(tmp))  // method is non-const on its argument!
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size) {
    pool_size = input.count_records_starting_with(identifier_tag);

    if (0 == pool_size) {
      cerr << "No occurrences of " << identifier_tag << "' in input\n";
      return 0;
    }

    pool = new GFP_L[pool_size];
    if (nullptr == pool) {
      cerr << "Yipes, could not allocate pool of size " << pool_size << '\n';
      return 62;
    }

    cerr << "Pool automatically sized to " << pool_size << '\n';
  }

  return build_pool(input);
}

/*
  If we have a previously selected file, we keep track of the number
  of members of the pool that get selected by the previously selected file
*/

static int molecules_selected_by_previously_selected_file = 0;

static int
do_previously_selected_file(IW_General_Fingerprint &fp, similarity_type_t t) {
  for (int i = 0; i < pool_size; i++) {
    GFP_L &p = pool[i];

    if (p.selected()) {  // probably hit during an earlier pass
      continue;
    }

    if (!can_be_compared(fp, p)) {
      continue;
    }

    similarity_type_t d;
    if (!compute_the_distance(p, fp, d)) {
      continue;
    }

    if (d > t) {
      continue;
    }

    p.set_selected(1);

    molecules_selected_by_previously_selected_file++;

    if (verbose > 1) {
      cerr << p.id() << " distance " << d << " from previously selected '" << fp.id()
           << "'\n";
    }

    if (stream_for_discarded_by_previously_selected.rdbuf()->is_open()) {
      stream_for_discarded_by_previously_selected << identifier_tag << p.id() << ">\n";
      stream_for_discarded_by_previously_selected << identifier_tag << fp.id() << ">\n";
      stream_for_discarded_by_previously_selected << distance_tag << d << ">\n";
      stream_for_discarded_by_previously_selected << "|\n";
    }
  }

  return 1;
}

static int
do_previously_selected_file(IW_TDT &tdt, similarity_type_t t) {
  IW_General_Fingerprint fp;

  int fatal;
  if (!fp.construct_from_tdt(tdt, fatal)) {
    cerr << "Cannot construct fingerprint from TDT\n";
    cerr << tdt;

    if (fatal) {
      return 0;
    }

    return 1;
  }

  return do_previously_selected_file(fp, t);
}

static int
do_previously_selected_file(iwstring_data_source &input, similarity_type_t t) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    if (!do_previously_selected_file(tdt, t)) {
      return 0;
    }
  }

  return 1;
}

static int
do_previously_selected_file(const IWString &fname, similarity_type_t t) {
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open previously selected file '" << fname << "'\n";
    return 0;
  }

  return do_previously_selected_file(input, t);
}

static int
set_thresholds_via_factor(score_t mvieth_factor) {
  score_t global_max_score = -std::numeric_limits<float>::max();
  for (int i = 0; i < pool_size; i++) {
    score_t s = pool[i].score();
    if (s > global_max_score) {
      global_max_score = s;
    }
  }

  if (verbose) {
    cerr << "Max score in pool is " << global_max_score << '\n';
  }

  for (int i = 0; i < pool_size; i++) {
    GFP_L &fp = pool[i];

    similarity_type_t t = (global_max_score - fp.score()) / mvieth_factor;
    fp.set_threshold(t);

    if (verbose > 1) {
      cerr << "i = " << i << " '" << fp.id() << " score " << fp.score()
           << " threshold set to " << t << '\n';
    }
  }

  return 1;
}

static int
process_dash_e_option(Command_Line &cl, char e,
                      resizable_array_p<IWString> &items_to_echo) {
  if (!cl.option_present(e)) {
    items_to_echo.resize(2);
    IWString *t = new IWString(smiles_tag);
    items_to_echo.add(t);
    t = new IWString(identifier_tag);
    items_to_echo.add(t);

    return 1;
  }

  int all_found = 0;

  const_IWSubstring evalue;
  int i = 0;
  while (cl.value(e, evalue, i++)) {
    if ("ALL" == evalue) {
      all_found = 1;
      if (verbose) {
        cerr << "Will echo entire tdt on output\n";
      }
    } else {
      IWString *t = new IWString(evalue);
      items_to_echo.add(t);
      if (verbose) {
        cerr << "Will echo item '" << evalue << "'\n";
      }
    }
  }

  if (all_found && items_to_echo.number_elements()) {
    cerr << "Using '-" << e << " ALL' and other -" << e
         << " options doesn't make sense\n";
    return 0;
  }

  return 1;
}

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
  cerr << "Performs leader clustering on a set of fingerprints\n";
  cerr << "Usage <options> <input_file>\n";
  cerr << " -C <ncluster>    maximum number of clusters to find\n";
  cerr << " -t <dis>         specify distance threshold\n";
  cerr << " -t col=<nn>      threshold is column <nn> of the name field\n";
  cerr << " -t tag=<TAG>     threshold for each molecule in dataitem <TAG>\n";
  cerr << " -H <TAG>         threshold for each molecule in dataitem <TAG>\n";
  cerr << " -m <number>      maximum cluster size\n";
  cerr << " -M <tag>         max cluster size for each molecule in <TAG>\n";
  cerr << " -M col=nn        max cluster size is column <nn> of the name field\n";
  cerr << " -S <TAG>         score tag\n";
  cerr << " -S col=nn        score is column <nn> in the name field\n";
  cerr << " -I <TAG>         specify identifier tag\n";
  cerr << " -r               sort clusters by distance from leader\n";
  cerr << " -E <dataitem>    specify pool object dataitems to be echo'd (default $SMI and PCN)\n";
  cerr << " -E ALL           echo all dataitems from the pool file\n";
  // cerr << " -X <distance>    abandon distance computation if any component >
  // distance\n"; cerr << " -Y <factor>      thresold = (max_score - score) / factor\n";
  // cerr << " -R <factor>      score = score + <factor> * distance to cluster\n";
  cerr << " -A <file>        file(s) of previously selected molecules - discard all within threshold\n";
  cerr << " -a <dist>        use <dist> as the threshold when comparing against the -A file\n";
  cerr << " -L <file>        write fingerprints discarded by -A file(s)\n";
  cerr << " -s <size>        specify max pool size\n";
  cerr << " -D ...           miscellaneous options, enter '-D help' for info\n";
  cerr << " -F ...           gfp options, enter '-F help' for details\n";
  cerr << " -V ...           Tversky specification, enter '-V help' for details\n";
  cerr << " -k .             write neighbours as textproto\n";
  cerr << " -v               verbose output\n";
// clang-format on

  exit(rc);
}

static int
display_dash_D_options(std::ostream &output) {
  output << " -D sbc=<dist>    cluster each cluster, constant threshold <dist>\n";
  output << " -D sbcr=<ratio>  cluster each cluster. Use ratio of cluster threshold\n";

  exit(0);

  return 0;
}

static void
DisplayDashKOptions(std::ostream& output) {
  output << R"(-k .             default textproto output
-k nozero               write a small value rather than zero in textproto output
)";
  ::exit(0);
}

static int
leader(int argc, char **argv) {
  Command_Line cl(argc, argv, "vs:I:E:t:F:P:W:X:rH:S:C:Y:V:R:m:M:A:L:a:D:k:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('D')) {
    int i = 0;
    const_IWSubstring d;
    while (cl.value('D', d, i++)) {
      const_IWSubstring directive;
      float dvalue;

      if (!d.split_into_directive_and_value(directive, '=', dvalue)) {
        cerr << "Invalid -D qualifier '" << d << "'\n";
        display_dash_D_options(cerr);
      }

      if ("sbc" == directive) {
        sub_cluster_threshold = dvalue;
        if (verbose) {
          cerr << "Fixed sub-cluster threshold " << sub_cluster_threshold << '\n';
        }

        leader_on_clusters = 1;
      } else if ("sbcr" == directive) {
        sub_cluster_threshold_ratio = dvalue;
        if (verbose) {
          cerr << "Variable sub-cluster threshold ratio " << sub_cluster_threshold_ratio
               << '\n';
        }

        leader_on_clusters = 1;
      } else if ("help" == directive) {
        display_dash_D_options(cerr);
      } else {
        cerr << "Unrecognised -D qualifier '" << d << "'\n";
        display_dash_D_options(cerr);
      }
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 1) {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }

    pool = new GFP_L[pool_size];
    if (nullptr == pool) {
      cerr << "Yipes, could not allocate pool of size " << pool_size << '\n';
      return 62;
    }

    if (verbose) {
      cerr << "system sized to " << pool_size << '\n';
    }
  }

  if (cl.option_present('V')) {
    if (!tversky.parse_command_line(cl, 'V', verbose)) {
      cerr << "Cannot parse Tversky specifications\n";
      usage(38);
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

  if (cl.option_present('r')) {
    sort_by_distance_from_centre = 1;
    if (verbose) {
      cerr << "Clusters sorted by distance from centre\n";
    }
  }

  if (cl.option_present('k')) {
    write_as_proto = 1;
    const_IWSubstring k;
    for (int i = 0; cl.value('k', k, i); ++i) {
      if (k == '.') {
        continue;
      }
      if (k == "nozero") {
        offset_zero_distances_in_textproto = 1;
        if (verbose) {
          cerr << "Will offset zero distance values in textproto output\n";
        }
      } else if (k == "help") {
        DisplayDashKOptions(cerr);
      } else {
        cerr << "Unrecognised -k qualifier '" << k << "'\n";
        DisplayDashKOptions(cerr);
      }
    }
    if (verbose) {
      cerr << "Will write neighbours as nnbr::NearNeighbours textproto form\n";
    }
  }

  // We need to be careful with the -i and -I options. Remember
  // that the pool is built first

  if (cl.option_present('I')) {
    (void)cl.value('I', identifier_tag);

    set_identifier_tag(identifier_tag);

    if (verbose) {
      cerr << "Identifiers in dataitem '" << identifier_tag << "'\n";
    }
  }

  if (cl.option_present('H')) {
    cl.value('H', threshold_from_file_tag);
    if (verbose) {
      cerr << "Each threshold from the '" << threshold_from_file_tag
           << "' dataitem in the input\n";
    }

    if (!threshold_from_file_tag.ends_with('<')) {
      threshold_from_file_tag << '<';
    }
  }

  if (cl.option_present('S')) {
    const_IWSubstring s = cl.string_value('S');

    if (s.starts_with("col=")) {
      s.remove_leading_chars(4);
      if (!s.numeric_value(score_column) || score_column < 1) {
        cerr << "Invalid column for score '" << s << "'\n";
        usage(14);
      }

      if (verbose) {
        cerr << "Score for each item in column " << score_column << '\n';
      }

      score_column--;
    } else {
      score_tag = s;
      if (verbose) {
        cerr << "Score tag is " << score_tag << "'\n";
      }
    }
  }

  if (cl.option_present('R')) {
    if (!cl.value('R', cluster_distance_scale_factor) ||
        cluster_distance_scale_factor <= 0.0) {
      cerr << "The cluster distance scale factor option (-R) must be followed by a "
              "positive number\n";
      usage(19);
    }

    if (verbose) {
      cerr << "Scores adjusted by " << cluster_distance_scale_factor
           << " times distance to nearest cluster centre\n";
    }
  }

  if (cl.option_present('F') || cl.option_present('P') || cl.option_present('W')) {
    if (!initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise general fingerprint options\n";
      usage(17);
    }
  }

  if (cl.option_present('X')) {
    if (!cl.value('X', abandon_distance_cutoff) || abandon_distance_cutoff < 0.0 ||
        abandon_distance_cutoff > 1.0) {
      cerr << "The -X option must be followed by a valid distance (0.0, 1.0)\n";
      usage(13);
    }

    if (verbose) {
      cerr << "Distance compuations abandoned if any component > "
           << abandon_distance_cutoff << '\n';
    }
  }

  if (!process_dash_e_option(cl, 'E', dataitems_to_echo)) {
    cerr << "Cannot process -E option\n";
    usage(15);
  }

  if (!cl.option_present('t') && !cl.option_present('H') && !cl.option_present('Y')) {
    cerr << "Threshold distance must be specified via -t, -H or -V options\n";
    usage(28);
  }

  if (need_to_call_initialise_fingerprints(cl)) {
    if (!initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise GFP options\n";
      usage(23);
    }
  }

  if (cl.option_present('t')) {
    int i = 0;
    const_IWSubstring t;
    while (cl.value('t', t, i++)) {
      if (t.starts_with("col=")) {
        t.remove_leading_chars(4);
        if (!t.numeric_value(threshold_column) || threshold_column < 1) {
          cerr << "Invalid column for threshold '" << t << "'\n";
          usage(14);
        }

        if (verbose) {
          cerr << "Threshold for each item in column " << threshold_column << '\n';
        }

        threshold_column--;
      } else if (t.starts_with("tag=")) {
        threshold_from_file_tag = t;
        threshold_from_file_tag.remove_leading_chars(4);
        if (verbose) {
          cerr << "Threshold in tag " << threshold_from_file_tag << "'\n";
        }

        if (!threshold_from_file_tag.ends_with('<')) {
          threshold_from_file_tag.add('<');
        }
      } else if (!t.numeric_value(threshold) || threshold < 0.0 || threshold > 1.0) {
        cerr << "The -t option must be followed by a valid distance value\n";
        usage(12);
      } else {
        if (verbose) {
          cerr << "Distance threshold set to " << threshold << '\n';
        }
      }
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', max_cluster_size) || max_cluster_size < 2) {
      cerr << "The -m (max cluster size) option must be followed by a whole number > 1\n";
      usage(43);
    }

    if (verbose) {
      cerr << "Max cluster size " << max_cluster_size << '\n';
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
             << max_cluster_size_column << '\n';
      }

      max_cluster_size_column--;
    } else {
      max_cluster_size_tag = m;

      if (verbose) {
        cerr << "Max cluster size in '" << max_cluster_size_tag << "' tag\n";
      }
    }
  }

  if (cl.number_elements() > 1) {
    cerr << "Extra arguments ignored\n";
  }

  iwstring_data_source pool_file;

  if (!build_pool(cl[0], pool_file) || 0 == pool_size) {
    cerr << "Cannot build pool from '" << cl[0] << "'\n";
    return 21;
  }

  if (cl.option_present('A')) {
    if (!cl.option_present('t') && !cl.option_present('a')) {
      cerr << "Must have a threshold available with the -A option (use -t or -a)\n";
      usage(11);
    }

    if (cl.option_present('L')) {
      IWString fname = cl.string_value('L');
      stream_for_discarded_by_previously_selected.open(fname.null_terminated_chars(),
                                                       std::ios::out);
      if (!stream_for_discarded_by_previously_selected.good()) {
        cerr << "Cannot open stream for molecules discarded by previously selected '"
             << fname << "'\n";
        return 4;
      }

      if (verbose) {
        cerr << "Molecules too close to previously selected file(s) written to '" << fname
             << "'\n";
      }
    }

    similarity_type_t t;

    if (cl.option_present('a')) {
      if (!cl.value('a', t) || t < 0.0 || t >= 1.0) {
        cerr << "Invalid value for previously selected threshold (-a option)\n";
        usage(4);
      }

      if (verbose) {
        cerr << "Will use " << t
             << " as the threshold for the previously selected list\n";
      }
    } else {
      t = threshold;
    }

    IWString fname;
    int i = 0;
    while (cl.value('A', fname, i++)) {
      if (!do_previously_selected_file(fname, t)) {
        cerr << "Cannot process previously selected file (-A option)\n";
        return 8;
      }
    }

    if (verbose) {
      cerr << "Rejected " << molecules_selected_by_previously_selected_file
           << " molecules by previously selected file(s)\n";
    }

    if (molecules_selected_by_previously_selected_file == pool_size) {
      cerr << "Yipes, the previously selected file knocked out the whole pool\n";
      return 1;
    }
  }

  /*
    threshold = (MAX_SCORE - score) / factor
    where factor is a user settable number
  */

  if (cl.option_present('Y')) {
    if (!cl.option_present('S')) {
      cerr << "Scores must be present(-S) in order to use -V\n";
      usage(42);
    }

    if (cl.option_present('H')) {
      cerr << "The -H (threshold in file) and -Y options are mutually exclusive\n";
      usage(31);
    }

    score_t mvieth_factor;

    if (!cl.value('Y', mvieth_factor) || mvieth_factor <= 0.0) {
      cerr << "The Vieth factor (-Y) must be followed by a positive number\n";
      usage(12);
    }

    if (verbose) {
      cerr << "Vieth factor " << mvieth_factor << '\n';
    }

    set_thresholds_via_factor(mvieth_factor);
  }

  IWString_and_File_Descriptor output(1);

  try {
    leader(pool_file, output);
  } catch (const char *err) {
    cerr << "Caught '" << err << "' terminated\n";
    return 81;
  }

  std::cout.flush();

  if (verbose) {
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

    cerr << "In clusters " << isum << '\n';
  }

  return 0;
}

int
main(int argc, char **argv) {
  int rc = leader(argc, argv);

  return rc;
}
