/*
  An implementation of the leader algorithm for fingerprints
  Implements an optimisation for determining the new leader
*/

#include <stdlib.h>
#include <values.h>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/misc.h"

#include "Utilities/GFP_Tools/leader.h"
#include "Utilities/GFP_Tools/tversky.h"

using std::cerr;
using std::endl;

static Tversky tversky;

static similarity_type_t abandon_distance_cutoff = -2.0;

/*
  We can have the threshold for each item read from the file
*/

static IWString score_tag;

/*
  The number of optimsations to do
*/

static int nopt = 0;

/*
  We can keep track of how the cluster sizes change during optimisation
*/

static Accumulator_Int<int> initial_cluster_size, final_cluster_size;

static extending_resizable_array<int> optimisation_attempted;

/*
  May 99. When looking for the next molecule to select, we can
  add a multiple of the distance to the nearest cluster centre
  to the score. That way, we can include some varying function
  of diversity into the next cluster selection
*/

static similarity_type_t cluster_distance_scale_factor = 0.0;

static resizable_array_p<IWString> dataitems_to_echo;

static int verbose = 0;

static similarity_type_t abandon_distance_threshold = -1.0;

/*
  The variables which control the clustering
*/

static int max_clusters_to_find = INT_MAX;
static int clusters_found = 0;
static int max_cluster_size = 0;

static similarity_type_t threshold = 0.0;

/*
  when dealing with clusters which are decided by the threshold, we can
  optinally sort the cluster members by their distance from the leader
*/

static int sort_by_distance_from_centre = 0;

static int items_selected = 0;

extending_resizable_array<int> cluster_size;

static Accumulator<float> distance_stats;

/*
  The identifier tag used in each TDT
*/

static IWString identifier_tag ("PCN<");

GFP_L::GFP_L()
{
  _selected = 0;

  _score = 0.0;

  _shortest_distance_to_cluster_centre = 1.0;

  return;
}

static IWString threshold_from_file_tag;

int
GFP_L::construct_from_tdt (IW_TDT & tdt, int & fatal)
{
  if (! IW_General_Fingerprint::construct_from_tdt (tdt, fatal))
  {
    return 0;
  }

  if (threshold_from_file_tag.length())
  {
    similarity_type_t tmp;
    if (! tdt.dataitem_value (threshold_from_file_tag, tmp) || tmp < 0.0)
    {
      cerr << "GFP_L::construct_from_tdt: invalid '" << threshold_from_file_tag << "' in tdt\n";
      return 0;
    }

    _threshold.set (tmp);
  }

  if (0 == score_tag.length())
    return 1;

  if (! tdt.dataitem_value (score_tag, _score))
  {
    cerr << "GFP_L::construct_from_tdt: cannot extract '" << score_tag << "' from tdt\n";
    return 0;
  }

  return 1;
}

/*
  Our pool is an array of FP objects
*/

static GFP_L * pool = nullptr;

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
echo_selected_dataitems (IW_TDT & tdt,
                         const resizable_array_p<IWString> & items_to_echo,
                         IWString_and_File_Descriptor & output)
{
  int ne = items_to_echo.number_elements();
  for (int i = 0; i < ne; i++)
  {
    const IWString & tag = *(items_to_echo[i]);

    const_IWSubstring zdata;

    if (! tdt.dataitem_value (tag, zdata))
    {
      cerr << "Yipes, cannot find tag '" << (*(items_to_echo[i])) << "' in tdt\n";
      throw "Missing dataitem";
      return 0;
    }

    output << zdata << '\n';
  }

  return output.good();
}

static int
echo_all_dataitems (IW_TDT & tdt,
                    IWString_and_File_Descriptor & output)
{
  return tdt.write_all_except_vbar (output);
}

static int
write_cluster_data (IW_TDT & tdt,
                    int clusters_found,
                    int id_within_cluster, 
                    similarity_type_t distance_to_centre,
                    IWString_and_File_Descriptor & output)
{
  if (dataitems_to_echo.number_elements())
    echo_selected_dataitems (tdt, dataitems_to_echo, output);
  else
    echo_all_dataitems (tdt, output);

  output << "DIST<" << distance_to_centre << ">\n";

  return output.good();
}

static int
get_tdt (IW_TDT & tdt, iwstring_data_source & input,
         GFP_L & fp)
{
  off_t offset;
  (void) fp.offset (offset);

  if (! input.seekg (offset))
  {
    cerr << "Cannot seek to offset '" << offset << endl;
    return 0;
  }

  return tdt.next (input);
}

int
distance_comparitor (GFP_L * const * ppfp1, GFP_L * const * ppfp2)
{
  const GFP_L * pfp1 = *ppfp1;
  const GFP_L * pfp2 = *ppfp2;

  if (pfp1->distance() < pfp2->distance())
    return -1;
  else if (pfp1->distance() > pfp2->distance())
    return 1;
  else
    return 0;
}

static int
process_cluster (Cluster & cluster,
                 iwstring_data_source & input,
                 IWString_and_File_Descriptor & output)
{
  if (sort_by_distance_from_centre)
    cluster.sort (&distance_comparitor);

  int cs = cluster.number_elements();
  cluster_size[cs]++;      // the leader isn't in the cluster

  GFP_L * centre = cluster[0];

  if (verbose)
  {
    cerr << "Cluster " << clusters_found << " contains " << cs << " items, centre '" << cluster[0]->id() << "', ";
    if (threshold_from_file_tag.length())
    {
      similarity_type_t threshold = 0.0f;
      (void) centre->threshold (threshold);
      cerr << "threshold " << threshold << ", ";
    }
    cerr << (items_selected + cs) << " items selected\n";
  }

  IW_TDT tdt;
  if (! get_tdt (tdt, input, *centre))
    return 0;

  if (dataitems_to_echo.number_elements())
    echo_selected_dataitems (tdt, dataitems_to_echo, output);
  else
    echo_all_dataitems (tdt, output);

  output << "CLUSTER<" << clusters_found << ">\n";
  output << "CSIZE<" << cs << ">\n";

  for (int i = 1; i < cs && output.good(); i++)    // start at 1, we've already done centre above
  {
    GFP_L & fp = *(cluster[i]);

    assert (fp.selected());

    if (! get_tdt (tdt, input, fp))
      return 0;

    if (! write_cluster_data (tdt, clusters_found, i, fp.distance(), output))
      return 0;
  }

  output << "|\n";

  output.write_if_buffer_holds_more_than(32768);

  return output.good();
}

static int
choose_next_centre (int & icentre)
{
  icentre = -1;
  if (0 == score_tag.length())    // just grab the first unselected item
  {
    for (int i = 0; i < pool_size; i++)
    {
      if (! pool[i].selected())
      {
        icentre = i;
        return 1;
      }
    }
  }
  else if (cluster_distance_scale_factor > 0.0)
  {
    score_t max_score = 0.0;
    for (int i = 0; i < pool_size; i++)
    {
      if (pool[i].selected())
        continue;

      score_t s = pool[i].score() + cluster_distance_scale_factor * pool[i].shortest_distance_to_cluster_centre();

      if (icentre < 0 || s > max_score)
      {
        max_score = s;
        icentre = i;
      }
    }
  }
  else    // raw scores
  {
    score_t max_score = 0.0;
    for (int i = 0; i < pool_size; i++)
    {
      if (pool[i].selected())
        continue;

      score_t s = pool[i].score();

      if (icentre < 0 || s > max_score)
      {
        max_score = s;
        icentre = i;
      }
    }
  }

  return icentre >= 0;
}

static int
form_cluster_threshold (int icentre, Cluster & cluster,
                        const similarity_type_t my_threshold)
{
  GFP_L & leader = pool[icentre];

  int istart = first_unselected;
  int istop = last_unselected;

  first_unselected = -1;
  for (int i = istart; i <= istop; i++)
  {
    GFP_L & p = pool[i];

    if (p.selected())
      continue;

    similarity_type_t d;

    if (! can_be_compared (leader, p))
      d = my_threshold + 1.0;
    else if (tversky.active())
      d = (1.0 - leader.IW_General_Fingerprint::tversky (p, tversky));
    else if (abandon_distance_cutoff > 0.0)
    {
      if (! leader.IW_General_Fingerprint::tanimoto (p, abandon_distance_cutoff, d))
        d = my_threshold + 1.0;
      else
        d = static_cast<similarity_type_t> (1.0) - d;
    }
    else
      d = leader.IW_General_Fingerprint::distance (p);

    if (d > my_threshold)
    {
      if (first_unselected < 0)
        first_unselected = i;

      last_unselected = i;
      if (d < p.shortest_distance_to_cluster_centre())
        p.set_shortest_distance_to_cluster_centre (d);
    }
    else
    {
      cluster.add (&(pool[i]));
      p.set_selected (clusters_found + 1);
      p.set_distance (d);
    }
  }

  return cluster.number_elements();
}

static int
form_cluster_max_cluster_size (int icentre, Cluster & cluster)
{
  assert (max_cluster_size > 0);

  cluster.resize (max_cluster_size);

  GFP_L & fp = pool[icentre];

  for (int i = 0; i < pool_size; i++)
  {
    GFP_L & p = pool[i];

    if (p.selected())
      continue;

    similarity_type_t d;

    if (! can_be_compared (fp, p))
      continue;
    else if (tversky.active())
      d = (1.0 - fp.IW_General_Fingerprint::tversky (p, tversky));
    else if (abandon_distance_cutoff > 0.0)
    {
      if (! fp.IW_General_Fingerprint::tanimoto (p, abandon_distance_cutoff, d))
        continue;
      else
        d = static_cast<similarity_type_t> (1.0) - d;
    }
    else
      d = fp.IW_General_Fingerprint::distance (p);

    cluster.add (&(pool[i]));
    p.set_distance (d);
    if (d < p.shortest_distance_to_cluster_centre())
      p.set_shortest_distance_to_cluster_centre (d);
  }

  cluster.resize_keep_storage (max_cluster_size);

  cluster.sort (&distance_comparitor);

  int istop;
  if (cluster.number_elements() < max_cluster_size)
    istop = cluster.number_elements();
  else
    istop = max_cluster_size;

  for (int i = 0; i < istop; i++)
  {
    GFP_L * p = cluster[i];

    p->set_selected (1);
  }

  return 1;
}

/*
  The clustering will be limited either by the maximum number of items which
  can be in a cluster, or a threshold
*/

static int
form_cluster (int icentre, Cluster & cluster)
{
  cluster.resize_keep_storage (0);

  cluster.add (&(pool[icentre]));

  pool[icentre].selected() = 1;
  pool[icentre].set_distance (0.0);

  similarity_type_t my_threshold;
  if (pool[icentre].threshold (my_threshold))     // has come from the file
    ;
  else
    my_threshold = threshold;

  if (verbose > 1)
    cerr << "Cluster " << clusters_found << " being determined. Threshold = " << my_threshold << endl;

  if (max_cluster_size)
    return form_cluster_max_cluster_size (icentre, cluster);
  else
    return form_cluster_threshold (icentre, cluster, my_threshold);
}

static void
unselect_cluster (Cluster & cluster)
{
  int cluster_size = cluster.number_elements();

  for (int i = 0; i < cluster_size; i++)
  {
    cluster[i]->set_selected (0);
    cluster[i]->set_distance (1.0);
  }

  first_unselected = 0;
  last_unselected = pool_size - 1;

  return;
}

/*
  Work out the sum of distances from each cluster member to all other cluster members
*/

static void
compute_intra_cluster_distances (const Cluster & cluster,
                                 similarity_type_t * tmp)
{
  int cluster_size = cluster.number_elements();

  set_vector (tmp, cluster_size, static_cast<similarity_type_t> (0.0));

  for (int i = 0; i < cluster_size; i++)
  {
    GFP_L & fpi = *(cluster[i]);

    for (int j = i + 1; j < cluster_size; j++)
    {
      GFP_L & fpj = *(cluster[j]);

      if (! can_be_compared (fpi, fpj))    // hmmm, is this correct?
        continue;

      similarity_type_t d;

      if (tversky.active())
        d = (1.0 - fpi.IW_General_Fingerprint::tversky (fpj, tversky));
      else
        d = fpi.IW_General_Fingerprint::distance (fpj);

      tmp[i] += d;    // this distance is associated with both cluster members
      tmp[j] += d;
    }
  }

  return;
}

static int
optimise_cluster (Cluster & cluster, 
                  resizable_array<GFP_L *> & already_tried,
                  int optimisations_completed)
{
  optimisation_attempted[optimisations_completed]++;

  int cluster_size = cluster.number_elements();

  if (verbose > 1)
    cerr << "Optimising level " << optimisations_completed << " with " << cluster_size << " items, leader '" << cluster[0]->id() << "'\n";

  similarity_type_t * tmp = new similarity_type_t[cluster_size];    // distances associated with each member of the cluster
  std::unique_ptr<similarity_type_t[]> free_tmp(tmp);

  compute_intra_cluster_distances (cluster, tmp);

// which one has the shortest total distance

  similarity_type_t shortest = tmp[0];
  GFP_L * zshortest = cluster[0];

  if (verbose > 2)
    cerr << "Cluster member " << 0 << " distances " << tmp[0] << " id '" << cluster[0]->id() << "' (leader)\n";
  for (int i = 1; i < cluster_size; i++)
  {
    if (verbose > 2)
      cerr << "Cluster member " << i << " distances " << tmp[i] << " id '" << cluster[i]->id() << "'\n";

    if (tmp[i] <= shortest)    // Using <= means we can optimise clusters with two members
    {
      shortest = tmp[i];
      zshortest = cluster[i];
    }
  }

  if (verbose > 1)
    cerr << "Shortest distance " << shortest << " for '" << zshortest->id() << "'\n";

  if (already_tried.last_item() == zshortest)    // no change from previous iteration
  {
    if (verbose > 2)
      cerr << "Same as last iteration\n";
    return 1;
  }

  GFP_L * initial_leader = cluster[0];

// We are going to make a change. Unset everything in this cluster

  unselect_cluster (cluster);

// Form cluster needs the pool index of the leader, use pointer arithmetic

  int leader = zshortest - pool;

  form_cluster (leader, cluster);

  assert (cluster[0] == zshortest);

  if (verbose > 2)
    cerr << "New cluster has " << cluster.number_elements() << " members, previous " << cluster_size << ", opt = " << optimisations_completed << endl;

  if (cluster.number_elements() < cluster_size)    // we cannot shrink the cluster size
  {
    if (verbose > 2)
      cerr << "Size reduced, optimisation terminated\n";

    unselect_cluster (cluster);

    leader = initial_leader - pool;

    return form_cluster (leader, cluster);    // re-form the initial cluster
  }

  if (already_tried.contains (zshortest))    // we are cycling
  {
    if (verbose > 2)
      cerr << "Cycling, done\n";

    return 1;
  }

  if (optimisations_completed == nopt)    // cannot do any more optimisation
  {
    if (verbose > 2)
      cerr << nopt << " cycles of optimisation done\n";
    return 1;
  }

  already_tried.add (cluster[0]);

  if (verbose > 2)
    cerr << "Leader optimisation trying " << leader << endl;

  return optimise_cluster (cluster, already_tried, optimisations_completed + 1);
}

static int
optimise_cluster (Cluster & cluster)
{
  resizable_array<GFP_L *> tmp;    // keep track of which ones have been tried as leaders
  tmp.resize (nopt);

  tmp.add (cluster[0]);

  return optimise_cluster (cluster, tmp, 0);
}

int
leader (iwstring_data_source & input,
        IWString_and_File_Descriptor & output)
{
  assert (pool_size > 1);

  assert (0 == items_selected);
  assert (0 == clusters_found);

  first_unselected = 0;
  last_unselected = pool_size - 1;

  int icentre = 0;

  Cluster cluster;
  if (! cluster.resize (pool_size))
  {
    cerr << "Yipes, cannot allocate " << pool_size << " elements in pool\n";
    return 0;
  }

  while (items_selected < pool_size)
  {
    (void) form_cluster (icentre, cluster);

    if (nopt && cluster.number_elements() > 1)
    {
      initial_cluster_size.extra (cluster.number_elements());
      optimise_cluster (cluster);
      final_cluster_size.extra (cluster.number_elements());
    }

    (void) process_cluster (cluster, input, output);

    clusters_found++;
    if (clusters_found >= max_clusters_to_find)
      break;

    items_selected += cluster.number_elements();

    if (! choose_next_centre (icentre))
      break;
  }

  return 1;
}

static int
build_pool (iwstring_data_source & input)
{
  off_t offset = input.tellg();

  int items_in_pool = 0;

  int tdts_read = 0;

  IW_TDT tdt;
  while (tdt.next (input))
  {
    tdts_read++;

    int fatal;
    if (! pool[items_in_pool].construct_from_tdt (tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      offset = input.tellg();
      continue;
    }

    pool[items_in_pool].set_offset (offset);

    items_in_pool++;

    if (items_in_pool == pool_size)
    {
      cerr << "Pool is full, max " << pool_size << endl;
      break;
    }

    offset = input.tellg();
  }

  pool_size = items_in_pool;

  if (verbose)
    cerr << "Read " << tdts_read << " TDT's, pool contains " << pool_size << " fingerprints\n";

  return 1;
}

static int
build_pool (const const_IWSubstring & fname,
            iwstring_data_source & input)
{
  IWString tmp (fname);

  if (! input.open (tmp))    // method is non-const on its argument!
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size)
  {
    IWString tmp;
    tmp << '^' << identifier_tag;

    std::unique_ptr<re2::RE2> pcn = std::make_unique<RE2>(tmp);
    pool_size = input.grep(*pcn);

    if (0 == pool_size) {
      cerr << "No occurrences of " << pcn->pattern() << "' in input\n";
      return 0;
    }

    pool = new GFP_L[pool_size];
    if (pool == nullptr) {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    cerr << "Pool automatically sized to " << pool_size << endl;
  }

  return build_pool (input);
}

static int
set_thresholds_via_factor (score_t mvieth_factor)
{
  score_t global_max_score = -FLT_MAX;
  for (int i = 0; i < pool_size; i++)
  {
    score_t s = pool[i].score();
    if (s > global_max_score)
      global_max_score = s;
  }

  if (verbose)
    cerr << "Max score in pool is " << global_max_score << endl;

  for (int i = 0; i < pool_size; i++)
  {
    GFP_L & fp = pool[i];

    similarity_type_t t = (global_max_score - fp.score()) / mvieth_factor;
    fp.set_threshold (t);

    if (verbose > 1)
      cerr << "i = " << i << " '" << fp.id() << " score " << fp.score() << " threshold set to " << t << endl;
  }

  return 1;
}

static int
process_dash_e_option (Command_Line & cl,
                       char e, 
                       resizable_array_p<IWString> & items_to_echo)
{
  if (! cl.option_present (e))
  {
    items_to_echo.resize (2);
    IWString * t = new IWString ("$SMI<");
    items_to_echo.add (t);
    t = new IWString ("PCN<");
    items_to_echo.add (t);

    return 1;
  }

  int all_found = 0;

  const_IWSubstring evalue;
  int i = 0;
  while (cl.value (e, evalue, i++))
  {
    if ("ALL" == evalue)
    {
      all_found = 1;
      if (verbose)
        cerr << "Will echo entire tdt on output\n";
    }
    else
    {
      IWString * t = new IWString (evalue);
      items_to_echo.add (t);
      if (verbose)
        cerr << "Will echo item '" << evalue << "'\n";
    }
  }

  if (all_found && items_to_echo.number_elements())
  {
    cerr << "Using '-" << e << " ALL' and other -" << e << " options doesn't make sense\n";
    return 0;
  }

  return 1;
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Finds near neighbours of a set of fingerprints\n";
  cerr << "Usage <options> <input_file>\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -t <dis>         specify distance threshold\n";
  cerr << " -C <number>      maximum number of clusters to find\n";
  cerr << " -H <TAG>         threshold for each molecule in dataitem <TAG>\n";
  cerr << " -I <TAG>         specify identifier tag\n";
  cerr << " -S <TAG>         score tag\n";
  cerr << " -R <factor>      score = score + <factor> * distance to cluster\n";
  cerr << " -r               sort clusters by distance from leader\n";
  cerr << " -E <dataitem>    specify pool object dataitems to be echo'd (default $SMI and PCN)\n";
  cerr << " -E ALL           echo all dataitems from the pool file\n";
  cerr << " -X <distance>    abandon distance computation if any component > distance\n";
  cerr << " -Z <factor>      thresold = (max_score - score) / factor\n";
  cerr << " -o <nopt>        the number of optimisation steps to perform\n";
  cerr << " -F -P ...        standard gfp options, enter '-F help' for info\n";
  cerr << " -V ...           standard Tversky options, enter '-V help' for info\n";
  cerr << " -v               verbose output\n";

  exit (rc);
}

static int
leader (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:I:E:t:F:P:W:X:rH:S:C:Z:V:R:o:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (1);
  }

  if (cl.option_present ('s'))
  {
    if (! cl.value ('s', pool_size) || pool_size < 1)
    {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage (3);
    }

    pool = new GFP_L[pool_size];
    if (pool == nullptr) {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    if (verbose)
      cerr << "system sized to " << pool_size << endl;
  }

  if (cl.option_present ('V'))
  {
    if (! tversky.parse_command_line (cl, 'V', verbose))
    {
      cerr << "Cannot parse Tversky specifications (-V option)\n";
      usage (38);
    }
  }

  if (cl.option_present ('C'))
  {
    if (! cl.value ('C', max_clusters_to_find) || max_clusters_to_find < 1)
    {
      cerr << "The -C option (max clusters to find) must be followed by a positive integer\n";
      usage (41);
    }

    if (verbose)
      cerr << "Will find a max of " << max_clusters_to_find << " clusters\n";
  }

  if (cl.option_present ('r'))
  {
    sort_by_distance_from_centre = 1;
    if (verbose)
      cerr << "Clusters sorted by distance from centre\n";
  }

// We need to be careful with the -i and -I options. Remember
// that the pool is built first

  if (cl.option_present ('I'))
  {
    (void) cl.value ('I', identifier_tag);

    set_identifier_tag (identifier_tag);

    if (verbose)
      cerr << "Identifiers in dataitem '" << identifier_tag << "'\n";
  }

  if (cl.option_present ('H'))
  {
    cl.value ('H', threshold_from_file_tag);
    if (verbose)
      cerr << "Each threshold from the '" << threshold_from_file_tag << "' dataitem in the input\n";
  }

  if (cl.option_present ('S'))
  {
    score_tag = cl.string_value ('S');
    if (verbose)
      cerr << "Score tag is " << score_tag << "'\n";
  }

  if (cl.option_present ('R'))
  {
    if (! cl.value ('R', cluster_distance_scale_factor) || cluster_distance_scale_factor <= 0.0)
    {
      cerr << "The cluster distance scale factor option (-R) must be followed by a positive number\n";
      usage (19);
    }

    if (verbose)
      cerr << "Scores adjusted by " << cluster_distance_scale_factor << " times distance to nearest cluster centre\n";
  }

  if (cl.option_present ('F') || cl.option_present ('P') || cl.option_present ('W'))
  {
    if (! initialise_fingerprints (cl, verbose))
    {
      cerr << "Cannot initialise general fingerprint options\n";
      usage (17);
    }
  }

  if (cl.option_present ('X'))
  {
    if (! cl.value ('X', abandon_distance_threshold) || abandon_distance_threshold < 0.0 || abandon_distance_threshold > 1.0)
    {
      cerr << "The -X option must be followed by a valid distance (0.0, 1.0)\n";
      usage (13);
    }

    if (verbose)
      cerr << "Distance compuations abandoned if any component > " << abandon_distance_threshold << endl;
  }

  if (! process_dash_e_option (cl, 'E', dataitems_to_echo))
  {
    cerr << "Cannot process -E option\n";
    usage (15);
  }

  if (! cl.option_present ('t') && ! cl.option_present ('H') && ! cl.option_present ('Z'))
  {
    cerr << "Threshold distance must be specified via -t, -H or -Z options\n";
    usage (28);
  }

  if (cl.option_present ('t'))
  {
    if (! cl.value ('t', threshold) || threshold < 0.0 || threshold > 1.0)
    {
      cerr << "The -t option must be followed by a valid distance value\n";
      usage (12);
    }

    if (verbose)
      cerr << "Distance threshold set to " << threshold << endl;
  }

  if (cl.option_present ('m'))
  {
    if (! cl.value ('m', max_cluster_size) || max_cluster_size < 2)
    {
      cerr << "The -m (max cluster size) option must be followed by a whole number > 1\n";
      usage (43);
    }

    if (verbose)
      cerr << "Max cluster size " << max_cluster_size << endl;
  }

  if (cl.option_present ('o'))
  {
    if (! cl.value ('o', nopt) || nopt <= 0)
    {
      cerr << "The optimisation level (-o flag) must be a whole positive number\n";
      usage (16);
    }

    if (verbose)
      cerr << "Will perform " << nopt << " optimisation steps\n";
  }

  if (cl.number_elements() > 1)
    cerr << "Extra arguments ignored\n";

  iwstring_data_source pool_file;

  if (! build_pool (cl[0], pool_file))
  {
    cerr << "Cannot build pool from '" << cl[0] << "'\n";
    return 21;
  }

/*
  Michael Vieth wants to be able to set the threshold as

    threshold = (MAX_SCORE - score) / factor

  where factor is a user settable number
*/

  if (cl.option_present ('Z'))
  {
    if (! cl.option_present ('S'))
    {
      cerr << "Scores must be present (-S) in order to use -Z\n";
      usage (42);
    }
    
    if (cl.option_present ('H'))
    {
      cerr << "The -H (threshold in file) and -Z options are mutually exclusive\n";
      usage (31);
    }

    score_t mvieth_factor;

    if (! cl.value ('Z', mvieth_factor) || mvieth_factor <= 0.0)
    {
      cerr << "The Vieth factor (-Z) must be followed by a positive number\n";
      usage (12);
    }

    if (verbose)
      cerr << "Vieth factor " << mvieth_factor << endl;

    set_thresholds_via_factor (mvieth_factor);
  }

  IWString_and_File_Descriptor output(1);

  try {
    leader (pool_file, output);
  }
  catch (const char * err)
  {
    cerr << "Caught '" << err << "' terminated\n";
    return 81;
  }

  if (verbose)
  {
    cerr << "Clustered " << pool_size << " fingerprints into " << clusters_found << " clusters\n";
    int isum = 0;
    for (int i = 0; i < cluster_size.number_elements(); i++)
    {
      int j = cluster_size[i];
      if (0 == j)
        continue;

      cerr << j << " clusters were of size " << i << " members\n";

      isum += j * i;
    }

    cerr << "In clusters " << isum << endl;

    int non_singleton_clusters = clusters_found - cluster_size[1];

    cerr << "Average cluster size for " << non_singleton_clusters << " non singleton clusters " << float (isum - cluster_size[1]) / float (non_singleton_clusters) << endl;

    for (int i = 0; i < optimisation_attempted.number_elements(); i++)
    {
      cerr << optimisation_attempted[i] << " optimisations attempted at level " << i << endl;
    }
    cerr << "Initial cluster sizes between " << initial_cluster_size.minval() << " and " << initial_cluster_size.maxval();
    if (initial_cluster_size.n() > 1)
      cerr << ", ave " << initial_cluster_size.average();
    cerr << endl;
    cerr << "Final   cluster sizes between " << final_cluster_size.minval() << " and " << final_cluster_size.maxval();
    if (final_cluster_size.n() > 1)
      cerr << ", ave " << final_cluster_size.average();
    cerr << endl;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = leader (argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status (stderr);
#endif

  return rc;
}
