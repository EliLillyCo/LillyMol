/*
  Someone has done some clustering by some method.
  Report statistics on the clusters
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "cluster_eval.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

static int verbose = 0;

static int highest_cluster_number = 0;

static int pool_size = 0;

static GFP_C * pool = nullptr;

static IWString identifier_tag ("PCN<");

static int do_intra_cluster_stats = 0;

static int do_inter_cluster_stats = 0;

static int cummins_centroid_stuff = 0;

static int remove_leading_zeros = 0;

static int
inter_cluster_distribution (int cluster, 
                            Accumulator<similarity_type_t> & global_inter,
                            ostream & output)
{
  similarity_type_t max_distance = 0.0;    // max for this cluster
  Accumulator<similarity_type_t> acc;      // for this cluster

  int items_in_cluster = 0;
  for (int i = 0; i < pool_size; i++)
  {
    GFP_C & fi = pool[i];

    if (cluster != fi.cluster ())
      continue;

    items_in_cluster++;

    for (int j = i + 1; j < pool_size; j++)
    {
      GFP_C & fj = pool[j];

      if (cluster == fj.cluster ())    // skip items in the same cluster
        continue;

      similarity_type_t d = fi.distance (fj);

      acc.extra (d);
      if (d > max_distance)
        max_distance = d;
    }
  }

  if (acc.n () <= 1)    // not sure this can happen
    return output.good ();

  if (verbose)
  {
    output << "Distances to cluster " << cluster << ", " << items_in_cluster << " members\n";
    output << "  between " << acc.minval () << " and " << acc.maxval ();
    if (acc.n () > 1)
      output << " average " << acc.average () << " variance " << acc.variance ();
    output << endl;
  }

  global_inter.extra (acc);

  return output.good ();
}

/*
*/

static int
intra_cluster_distribution_centroid (int cluster, 
                            Accumulator<similarity_type_t> & global_intra,
                            double * d,
                            double & global_intra_d,
                            ostream & output)
{
  resizable_array<int> cluster_member;
  cluster_member.resize (pool_size);

// First get the members of the cluster

  for (int i = 0; i < pool_size; i++)
  {
    GFP_C & fi = pool[i];

    if (cluster == fi.cluster ())
      cluster_member.add (i);
  }

  int cs = cluster_member.number_elements ();
  if (cs <= 1)    // singleton or no cluster with this ID
    return 1;

  for (int i = 0; i < cs; i++)
  {
    d[i] = 0.0;
  }

  for (int i = 0; i < cs; i++)
  {
    GFP_C & fpi = pool[cluster_member[i]];

    for (int j = i + 1; j < cs; j++)
    {
      GFP_C & fpj = pool[cluster_member[j]];

      if (! can_be_compared (fpi, fpj))
        continue;

      double dij = static_cast<double> (fpi.distance (fpj));

      d[i] += dij;
      d[j] += dij;
    }
  }

// The centroid is the point with the smallest summed distances

  double mind = d[0];
  int centroid = 0;    // will convert to actual pool member number later

  for (int i = 1; i < cs; i++)
  {
    if (d[i] < mind)
    {
      mind = d[i];
      centroid = i;    // will convert to actual pool member number later
    }
  }

  global_intra_d += mind;
  global_intra.extra (mind);

  centroid = cluster_member[centroid];

  if (verbose)
  {
    output << "Cluster " << cluster << ", " << cs << " members, centroid " << centroid << " (" << pool[centroid].id () << ")\n";
    output << "Intra-cluster distance sum " << mind << endl;
  }

  return 1;
}

/*
  Determine the centroid of the entire pool
*/

static int
global_centroid (double * d,
                 double & global_centroid_inter,
                 ostream & output)
{
  global_centroid_inter = 0.0;

  for (int i = 0; i < pool_size; i++)
  {
    d[i] = 0.0;
  }

  for (int i = 0; i < pool_size; i++)
  {
    GFP_C & fi = pool[i];
    for (int j = i + 1; j < pool_size; j++)
    {
//    Thought about putting in a call to can_be_compared here, but that wouldn't
//    work in this case. We want to know the total distance

      double dij = static_cast<double> (fi.distance (pool[j]));

      d[i] += dij;
      d[j] += dij;
    }
  }

  double mind = d[0];
  int centroid = 0;

  for (int i = 1; i < pool_size; i++)
  {
    if (d[i] < mind)
    {
      mind = d[i];
      centroid = i;
    }
  }

  global_centroid_inter = d[centroid];

  if (verbose)
    output << "Global centroid is " << centroid << ", '" << pool[centroid].id () << "', dist " << global_centroid_inter << endl;

  return 1;
}

static int
do_cummins_centroid_stuff (const int * items_in_cluster,
                           ostream & output)
{
  double * d = new double[pool_size];

  double global_intra_d = 0.0;
  Accumulator<similarity_type_t> global_intra;

  for (int i = 0; i <= highest_cluster_number; i++)
  {
    if (items_in_cluster[i] < 2)
      continue;

    (void) intra_cluster_distribution_centroid (i, global_intra, d, global_intra_d, output);
  }

  if (verbose)
  {
    output << global_intra.n () << " clusters. Intra cluster (centroid) distances between " << global_intra.minval () << " and " << global_intra.maxval () << endl;
  }

  double global_centroid_inter = 0.0;

  (void) global_centroid (d, global_centroid_inter, output);

  delete[] d;

  output << "Intra cluster distance " << global_intra_d << " across cluster total distance " << global_centroid_inter << endl;
  output << "Clustering metric " << (1.0 - global_intra_d / global_centroid_inter) << endl;

  return 1;
}

static int
intra_cluster_distribution (int cluster, 
                            Accumulator<similarity_type_t> & global_intra,
                            ostream & output)
{
  similarity_type_t max_distance = 0.0;    // max within this cluster
  Accumulator<similarity_type_t> acc;      // within this cluster

  int items_in_cluster = 0;
  for (int i = 0; i < pool_size; i++)
  {
    GFP_C & fi = pool[i];

    if (cluster != fi.cluster ())
      continue;

    items_in_cluster++;

    for (int j = i + 1; j < pool_size; j++)
    {
      GFP_C & fj = pool[j];

      if (cluster != fj.cluster ())
        continue;

      similarity_type_t d = fi.distance (fj);

      acc.extra (d);
      if (d > max_distance)
        max_distance = d;
    }
  }

  if (acc.n () <= 1)    // must be a singleton or maybe no clusters with this cluster ID
    return output.good ();

  if (verbose)
  {
    output << "Distances within cluster " << cluster << ", " << items_in_cluster << " members\n";
    output << "  between " << acc.minval () << " and " << acc.maxval ();
    if (acc.n () > 1)
      output << " average " << acc.average () << " variance " << acc.variance ();
    output << endl;
  }

  global_intra.extra (acc);

  return output.good ();
}

static int
evaluate_clustering (int * items_in_cluster,
                     ostream & output)
{
  int clusters_found = 0;

  for (int i = 0; i < pool_size; i++)
  {
    int c = pool[i].cluster ();

    if (0 == items_in_cluster[c])
      clusters_found++;

    items_in_cluster[c]++;
  }

  output << "Found " << clusters_found << " clusters\n";

  extending_resizable_array<int> cluster_size;

  cluster_size.resize (pool_size);

  Accumulator_Int<int> cs;
  for (int i = 0; i <= highest_cluster_number; i++)
  {
    int iic = items_in_cluster[i];

    if (0 == iic)
      continue;

    if (verbose > 1)
      output << iic << " items in cluster " << i << endl;

    cluster_size[iic]++;

    if (iic > 1)
      cs.extra (iic);
  }

  if (verbose)
  {
    output << "Non singleton cluster sizes\n";
    output << "  between " << cs.minval () << " and " << cs.maxval ();
    if (cs.n () > 1)
      output << " average " << cs.average () << " variance " << cs.variance ();
    output << endl;
  }

  for (int i = 0; i < cluster_size.number_elements (); i++)
  {
    if (0 == cluster_size[i])
      continue;

    output << cluster_size[i] << " clusters had " << i << " members\n";
  }

  if (cummins_centroid_stuff)
    do_cummins_centroid_stuff (items_in_cluster, output);

  if (do_intra_cluster_stats)
  {
    Accumulator<similarity_type_t> global_intra;

    for (int i = 0; i <= highest_cluster_number; i++)
    {
      if (0 == items_in_cluster[i])
        continue;
  
      intra_cluster_distribution (i, global_intra, output);
    }

    output << "Intra cluster distances, " << global_intra.n () << " distances\n";
    output << "  between " << global_intra.minval () << " and " << global_intra.maxval ();
    if (global_intra.n () > 1)
      output << " average " << global_intra.average () << " variance " << global_intra.variance ();
    output << endl;
  }

  if (do_inter_cluster_stats)
  {
    Accumulator<similarity_type_t> global_inter;

    for (int i = 0; i <= highest_cluster_number; i++)
    {
      if (0 == items_in_cluster[i])
        continue;
  
      inter_cluster_distribution (i, global_inter, output);
    }

    output << "Inter cluster distances, " << global_inter.n () << " distances\n";
    output << "  between " << global_inter.minval () << " and " << global_inter.maxval ();
    if (global_inter.n () > 1)
      output << " average " << global_inter.average () << " variance " << global_inter.variance ();
    output << endl;
  }

  return output.good ();
}

static int
build_pool (iwstring_data_source & input)
{
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

      continue;
    }

    if (pool[items_in_pool].cluster () > highest_cluster_number)
      highest_cluster_number = pool[items_in_pool].cluster ();

    items_in_pool++;

    if (items_in_pool == pool_size)
    {
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
  iwstring_data_source input;

  if (! input.open (fname))    // method is non-const on its argument!
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size)
  {
    IWString tmp;
    tmp << '^' << identifier_tag;

    std::unique_ptr<re2::RE2> pcn;
    iwre2::RE2Reset(pcn, tmp);
    pool_size = input.grep (*pcn);

    if (0 == pool_size)
    {
      cerr << "No occurrences of " << pcn->pattern() << "' in input\n";
      return 0;
    }

    pool = new GFP_C[pool_size];
    if (NULL == pool)
    {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    cerr << "Pool automatically sized to " << pool_size << endl;
  }

  return build_pool (input);
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Evaluate the results of a clustering\n";
  cerr << "Input must be fingerprints with an integer cluster number with each item\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -C <tag>         TAG for cluster number (default \"CLUSTER\")\n";
  display_standard_gfp_options (cerr);
  cerr << " -E               compute intEr cluster distances (can be expensive)\n";
  cerr << " -R               compute intRa cluster distances (the default)\n";
  cerr << " -D               do centroid based metrics (from Dave Cummins)\n";
  cerr << " -z               remove leading 0's when doing cross references\n";
  cerr << " -v               verbose output\n";

  exit (rc);
}

static int
determine_cluster_from_file (iwstring_data_source & input)
{
  IW_STL_Hash_Map_int id_cluster;

  const_IWSubstring buffer;

  int rc = 1;
  while (input.next_record (buffer))
  {
    buffer.strip_leading_blanks ();
    buffer.strip_trailing_blanks ();

    IWString id, cluster;
    int c;

    if (! buffer.split (id, ' ', cluster) || ! cluster.numeric_value (c) || c < 0)
    {
      cerr << "Invalid cluster info record '" << buffer << "'\n";
      rc = 0;
    }

    if (remove_leading_zeros)
      id.remove_leading_chars ('0');

    id_cluster[id] = c;
  }

  if (0 == rc)
    return 0;

  if (verbose)
    cerr << "Read " << id_cluster.size () << " id->cluster relationships\n";

  for (int i = 0; i < pool_size; i++)
  {
    GFP_C & pi = pool[i];

    IWString id = pi.id ();

    if (remove_leading_zeros)
      id.remove_leading_chars ('0');

    IW_STL_Hash_Map_int::const_iterator f = id_cluster.find (id);

    if (f == id_cluster.end ())
    {
      cerr << "No cluster data for '" << id << "'\n";
      rc = 0;
    }
    else
    {
      int c = (*f).second;
      pi.set_cluster (c);
      if (c > highest_cluster_number)
        highest_cluster_number = c;
    }
  }

  return rc;
}

static int
determine_cluster_from_file (const IWString & fname)
{
  iwstring_data_source input (fname);

  if (! input.ok ())
  {
    cerr << "Cannot open id-cluster cross reference file '" << fname << "'\n";
    return 0;
  }

  return determine_cluster_from_file (input);
}

static int
determine_cluster_from_identifier (int cluster_in_column)
{
  int rc = 1;

  for (int i = 0; i < pool_size; i++)
  {
    GFP_C & pi = pool[i];

    const_IWSubstring id = pi.id ();

    const_IWSubstring token;
    int c;
    if (! id.word (cluster_in_column, token) || ! token.numeric_value (c) || c < 0)
    {
      cerr << "Missing or invalid cluster data '" << id << "'\n";
      rc = 0;
    }
    else
    {
      pi.set_cluster (c);
      if (c > highest_cluster_number)
        highest_cluster_number = c;
    }
  }

  return rc;
}

int
evaluate_clustering (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:C:F:P:W:I:ERDz");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('F') || cl.option_present ('P') || cl.option_present ('W'))
  {
    if (! initialise_fingerprints (cl, verbose))
    {
      cerr << "Cannot initialise general fingerprint options\n";
      usage (17);
    }
  }

  if (cl.option_present ('z'))
  {
    remove_leading_zeros = 1;

    if (verbose)
      cerr << "Leading 0's will be discarded when doing cross references\n";
  }

  int cluster_in_column = -1;
  IWString file_with_cluster_info;

  if (cl.option_present ('C'))
  {
    const_IWSubstring c = cl.string_value ('C');

    if (c.starts_with ("FILE="))
    {
      c.remove_leading_chars (5);
      file_with_cluster_info = c;

      if (verbose)
        cerr << "Cluster info in file '" << file_with_cluster_info << "'\n";

      (void) set_cluster_id_tag ("");
    }
    else if (c.starts_with ("col="))
    {
      c.remove_leading_chars (4);
      if (! c.numeric_value (cluster_in_column) || cluster_in_column < 1)
      {
        cerr << "The cluster in column directive must be a whole positive number\n";
        usage (5);
      }

      if (verbose)
        cerr << "Cluster info in column " << cluster_in_column << " of the identifiers\n";

      (void) set_cluster_id_tag ("");
    }
    else
    {
      (void) set_cluster_id_tag (c);

      if (verbose)
        cerr << "Cluster ID in tag '" << c << "'\n";
    }
  }

  if (cl.option_present ('I'))
  {
    (void) cl.value ('I', identifier_tag);

    set_identifier_tag (identifier_tag);

    if (verbose)
      cerr << "Identifiers in dataitem '" << identifier_tag << "'\n";
  }

  if (cl.option_present ('E'))
  {
    do_inter_cluster_stats = 1;
    if (verbose)
      cerr << "Will compute inter cluster distances\n";
  }

  if (cl.option_present ('R'))
  {
    do_intra_cluster_stats = 1;
    if (verbose)
      cerr << "Will compute intra cluster distances\n";
  }

  if (cl.option_present ('D'))
  {
    cummins_centroid_stuff = 1;
    if (verbose)
      cerr << "Will compute Dave Cummins centroid based metrics\n";
  }

  if (0 == do_inter_cluster_stats && 0 == do_intra_cluster_stats && 0 == cummins_centroid_stuff)
  {
    do_intra_cluster_stats = 1;
    if (verbose)
      cerr << "Will compute intra cluster distances by default\n";
  }

  if (0 == cl.number_elements ())
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

    pool = new GFP_C[pool_size];
    if (NULL == pool)
    {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    if (verbose)
      cerr << "system sized to " << pool_size << endl;
  }

  

  if (cl.number_elements () > 1)
    cerr << "Warning, extra arguments ignored\n";

  if (! build_pool (cl[0]))
  {
    cerr << "Cannot build pool from '" << cl[0] << "'\n";
    return 13;
  }

  int rc;
  if (cluster_in_column >= 0)
    rc = determine_cluster_from_identifier (cluster_in_column);
  else if (file_with_cluster_info.length ())
    rc = determine_cluster_from_file (file_with_cluster_info);
  else
    rc = 1;

  if (0 == rc)
  {
    cerr << "Cannot determine cluster info\n";
    return 5;
  }

  if (highest_cluster_number < 0)
  {
    cerr << "No cluster number values found\n";
    return 11;
  }

  int * items_in_cluster = new int[highest_cluster_number + 1];
  for (int i = 0; i <= highest_cluster_number; i++)
  {
    items_in_cluster[i] = 0;
  }

  rc = evaluate_clustering (items_in_cluster, cout);

  delete[] items_in_cluster;

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = evaluate_clustering (argc, argv);

  return rc;
}
