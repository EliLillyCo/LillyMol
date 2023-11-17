/*
  We have divided a set of molecules into two. We generate statistics
  on the nature of the split
*/

#include <stdlib.h>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/histogram/iwhistogram.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "nn_results_.h"

const char * prog_name = nullptr;

static int verbose = 0;

static int input_is_distance_matrix = 0;

static int input_from_gfp_nearneighbours_dash_o = 0;

static float upper_distance_cutoff = static_cast<float> (1.0);

static int write_histogram = 1;

static IWString histogram_initialisation_string;

static float truncate_distances_at = static_cast<float> (1.0);

static int identifier_column = -1;

/*
  We can just consider a given number of neighbours across the divide
*/

static int neighbours_to_consider = 0;

/*
  When we are evaluating multiple sets, we need something to keep track of
  the overall statistics. 
  We determine the number of buckets in the histogram and keep an accumulator
  for each
*/

static int nsplits = 0;

static Accumulator_Int<int> * overall_accumulator_count = nullptr;
static Accumulator<float> * overall_accumulator_fraction = nullptr;

static Accumulator_Int<int> items_in_set_acc;

static resizable_array<double> distances_to_check;

static IWString index_tag ("NBR<");

/*
  When computing similarities between splits, we need an array of fingerprints
*/

static int find_inter_split_similarities = 0;

static IW_Bits_Base * fp = nullptr;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  if (verbose)
    cerr << "$Id: test_train_split_evaluate.cc,v 1.2 2002/10/03 15:25:07 ian Exp $\n";
  cerr << "Uses a distance matrix or output from gfp_nearneighbours to split a dataset\n";
  cerr << " -F <type>      what created the input file\n";
  cerr << "    gfp         gfp_nearneighbours (use the -h and -o options)\n";
  cerr << "    dm          gfp_distance_matrix\n";
  cerr << " -D <fname>     distance matrix or gfp_nearneighbours output\n";
  cerr << " -t <dist>      one or more distances to check\n";
  cerr << " -x <dist>      group all distances beyond <dist> into the same group\n";
  cerr << " -n <number>    only consider the <number> closest neighbours across the split\n";
  cerr << " -b             suppress writing the histogram\n";
  cerr << " -c <col>       identifiers are in column <col> of the input file\n";
  cerr << " -s             compute inter-split similarities - bits in common\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

class ID_Dist
{
  public:
    typedef float distance_t;

  private:

    distance_t _dist;

    int _index_in_pool;

  public:
    ID_Dist ();

    void set_index_in_pool (int i) { _index_in_pool = i;}
    int index_in_pool () const { return _index_in_pool;}

    int id () const { return _index_in_pool;}

    void set_distance (distance_t);

    distance_t distance () const { return _dist;}

    int build (iwstring_data_source &, int &);
};

ID_Dist::ID_Dist ()
{
  _dist = static_cast<distance_t> (-1.0);
  _index_in_pool = -1;

  return;
}

//static similarity_type_t f100 = static_cast<similarity_type_t> (100.0);

void
ID_Dist::set_distance (distance_t d)
{
  if (d < static_cast<distance_t> (0.0))
  {
    cerr << "ID_Dist::set_distance: invalid distance " << d << endl;
    _dist = static_cast<distance_t> (0.0);
    return;
  }

  _dist = d;

  return;
}

int
ID_Dist::build (iwstring_data_source & input,
                int & fatal)

{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (buffer.starts_with (index_tag))
    {
      if (! extract_from_tdt_form (buffer, index_tag, _index_in_pool) || _index_in_pool < 0)
      {
        cerr << "ID_Dist::build: invalid index specification '" << buffer << "'\n";
        return 0;
      }
    }
    else if (buffer.starts_with (distance_tag))
    {
      similarity_type_t tmp;
      if (! extract_from_tdt_form (buffer, distance_tag, tmp) || tmp < static_cast<similarity_type_t> (0.0))
      {
        cerr << "ID_Dist::build: invalid distance specification '" << buffer << "'\n";
        return 0;
      }

      set_distance (tmp);

      if (_index_in_pool >= 0)    // the order must be NBR followed by DIST
        return 1;

      cerr << "ID_Dist::build: got distance, but no index\n";
      fatal = 1;
      return 0;
      break;
    }
    else if ('|' == buffer)
    {
      if (_dist < 0 && _index_in_pool < 0)
      {
        fatal = 0;
        return 0;
      }
      break;
    }
  }

  cerr << index_tag << endl;
  cerr << distance_tag << endl;
  cerr << "ID_Dist::build: build incomplete, dist " << _dist << " ndx " << _index_in_pool << "'\n";
  fatal = 1;
  return 0;
}

template class resizable_array_p<ID_Dist>;
template class resizable_array_base<ID_Dist *>;

class NN_Item : public NN_Item_Base<ID_Dist>
{
  private:
  public:
};

class NN_Results : public NN_Results_Base<NN_Item>
{
  private:
  public:
};

static int
do_find_inter_split_similarities (ostream & output)
{
  assert (NULL != fp);

  Accumulator_Int<int> bic;

  for (int i = 0; i < nsplits; i++)
  {
    const IW_Bits_Base & fpi = fp[i];

    for (int j = i + 1; j < nsplits; j++)
    {
      IW_Bits_Base & fpj = fp[j];

      int b = fpi.bits_in_common(fpj);
      bic.extra(b);
    }
  }

  output << "Between " << nsplits << " splits, compared " << bic.n() << " pairs of bits in common\n";
  if (0 == bic.n())
    return 1;

  if (1 == bic.n())
  {
    output << "Single comparison had " << bic.minval() << " bits in common\n";
    return 1;
  }

  output << "Between " << bic.minval() << " and " << bic.maxval() << " items in common, ave " << static_cast<float> (bic.average()) << endl;

  return 1;
}

static void
do_truncate_distance_matrix (IWDistanceMatrixBase<float> & dm,
                             float t)
{
  int n = dm.number_molecules ();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    for (int j = i + 1; j < n; j++)
    {
      float d = dm.zvalue (i, j);

      if (d > t)
      {
        dm.set (i, j, t);
        rc++;
      }
    }
  }

  if (verbose)
    cerr << "Truncated " << rc << " distances to " << t << endl;
}

static void
do_truncate_neighbour_distances (NN_Results & nn_results,
                                 float t)
{
  int nr = nn_results.number_results ();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    NN_Item & ni = nn_results[i];

    for (int j = 0; j < ni.number_elements (); j++)
    {
      ID_Dist * nj = ni[j];

      if (nj->distance () > t)
      {
        nj->set_distance (t);
        rc++;
      }
    }
  }

  if (verbose)
    cerr << "Truncated " << rc << " distances to " << t << endl;

  return;
}


static int
report_accumulator (const Accumulator<double> & acc,
                    ostream & output)
{
  output << "N = " << acc.n ();
  output << " min " << acc.minval ();
  output << " ave " << static_cast<float> (acc.average_if_available_minval_if_not ());
  output << " max " << acc.maxval ();
  output << " tot " << acc.sum () << endl;

  return output.good ();
}

/*
  Sometimes we don't actually do any output, we just update the overal_accumulator
  objects
*/

static int
write_results (const Accumulator<double> & acc_set0,
               const Accumulator<double> & acc_set1,
               const Accumulator<double> & acc_across,
               const IWHistogram & histogram,
               ostream & output)
{
  int do_writes = 1;

//if (nsplits > 0 && 0 == verbose)
//  do_writes = 0;

  if (do_writes)
  {
    output << "Within set0\n";
    report_accumulator (acc_set0, output);
    output << "Within set\n";
    report_accumulator (acc_set1, output);
    output << "Between sets\n";
    report_accumulator (acc_across, output);
  }

  if (write_histogram)
  {
    const unsigned int * c = histogram.raw_counts ();
    unsigned int csum = 0;
    for (int i = 0; i < histogram.nbuckets () - 1; i++)
    {
      csum += c[i];

      if (NULL != overall_accumulator_count)
      {
        overall_accumulator_count[i].extra (c[i]);
        overall_accumulator_fraction[i].extra (static_cast<float> (c[i]) / static_cast<float> (histogram.nsamples ()));
      }
  
      if (do_writes)
        output << histogram.minval () + (i+1) * histogram.delta () << ' ' << c[i] << ' ' << csum << endl;
    }
  }

  if (distances_to_check.number_elements ())
  {
    output << "Inter set distances within -t distances\n";

    const unsigned int * c = histogram.raw_counts ();
    unsigned int csum = 0;
    int j = 0;
    for (int i = 0; i < histogram.nbuckets () - 1 && j < distances_to_check.number_elements (); i++)
    {
      csum += c[i];

      double d = histogram.minval () + (i+1) * histogram.delta ();

      if (d >= distances_to_check[j])
      {
        output << d << ' ' << csum << endl;
        j++;
      }
    }
  }

  return output.good ();
}

/*
  Special case where we are considering only the N nearest neighbours across the divide
*/

static int
test_train_split_evaluate_N (const IWDistanceMatrixBase<float> & dm,
                             const int * in_set,
                             ostream & output)
                                
{
  int nr = dm.number_molecules ();

  IWHistogram histogram;

  if (histogram_initialisation_string.length ())
    histogram.initialise (histogram_initialisation_string);
  else
    histogram.initialise (0.0, upper_distance_cutoff, upper_distance_cutoff / 10.0);

  Accumulator<double> acc_across;
  Accumulator<double> acc_set0;
  Accumulator<double> acc_set1;

  resizable_array<int> nbr;    // scope here for efficiency

  for (int i = 0; i < nr; i++)
  {
    int jstop = dm.neighbours_sorted_by_distance (i, nbr);

    for (int j = 0; j < jstop; j++)
    {
      int k = nbr[j];
      float d = dm.zvalue (i, k);

      if (0 == in_set[i] && 0 == in_set[k])
        acc_set0.extra (d);
      else if (1 == in_set[i] && 1 == in_set[k])
        acc_set1.extra (d);
      else if (neighbours_to_consider > 0 && acc_across.n () > neighbours_to_consider)   // already got the required number of inter-set distances
        ;
      else
      {
//      cerr << "Items " << i << " and " << k << " at distance " << nj->distance () << endl;
        histogram.extra (d);
        acc_across.extra (d);
      }
    }
  }

  return write_results (acc_set0, acc_set1, acc_across, histogram, output);
}

static int
test_train_split_evaluate (const IWDistanceMatrixBase<float> & dm,
                           const int * in_set,
                           ostream & output)
                                
{
  int nr = dm.number_molecules ();

  IWHistogram histogram;

  if (histogram_initialisation_string.length ())
    histogram.initialise (histogram_initialisation_string);
  else
    histogram.initialise (0.0, upper_distance_cutoff, upper_distance_cutoff / 10.0);

  Accumulator<double> acc_across;
  Accumulator<double> acc_set0;
  Accumulator<double> acc_set1;

  for (int i = 0; i < nr; i++)
  {
    for (int j = i + 1; j < nr; j++)
    {
      float d = dm.zvalue (i, j);

      if (0 == in_set[i] && 0 == in_set[j])
        acc_set0.extra (d);
      else if (1 == in_set[i] && 1 == in_set[j])
        acc_set1.extra (d);
      else if (neighbours_to_consider > 0 && acc_across.n () > neighbours_to_consider)
        ;
      else
      {
//      cerr << "Items " << i << " and " << k << " at distance " << nj->distance () << endl;
        histogram.extra (d);
        acc_across.extra (d);
      }
    }
  }

  return write_results (acc_set0, acc_set1, acc_across, histogram, output);
}

static int
test_train_split_evaluate (const NN_Results & nn_results,
                           const int * in_set,
                           ostream & output)
{
  int nr = nn_results.number_results ();

  IWHistogram histogram;

  if (histogram_initialisation_string.length ())
    histogram.initialise (histogram_initialisation_string);
  else
    histogram.initialise (0.0, upper_distance_cutoff, upper_distance_cutoff / 10);

  Accumulator<double> acc_across;
  Accumulator<double> acc_set0;
  Accumulator<double> acc_set1;

  for (int i = 0; i < nr; i++)
  {
    const NN_Item & ni = nn_results[i];

    for (int j = 0; j < ni.number_elements (); j++)
    {
      const ID_Dist * nj = ni[j];

      int k = nj->index_in_pool ();

      assert (k >= 0 && k < nr);

      ID_Dist::distance_t d = nj->distance ();

      if (0 == in_set[i] && 0 == in_set[k])
        acc_set0.extra (d);
      else if (1 == in_set[i] && 1 == in_set[k])
        acc_set1.extra (d);
      else if (neighbours_to_consider > 0 && acc_across.n () > neighbours_to_consider)
        ;
      else
      {
//      cerr << "Items " << i << " and " << k << " at distance " << nj->distance () << endl;
        histogram.extra (d);
        acc_across.extra (d);
      }
    }
  }

  return write_results (acc_set0, acc_set1, acc_across, histogram, output);
}

static int
identify_items_in_subset (iwstring_data_source & input,
                          const IW_STL_Hash_Map_int & id_to_index,
                          int nr,
                          int * in_set)
{
  set_vector (in_set, nr, 0);

  int items_in_set = 0;

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (buffer.starts_with ('#'))
      continue;

    if (identifier_column > 0)
      buffer.remove_leading_words (identifier_column);

    buffer.truncate_at_first (' ');

    if (0 == buffer.length ())
      continue;

    if (input.lines_read () > 1)
      ;
    else if ("Name" == buffer || "name" == buffer || "id" == buffer)
      continue;

    IW_STL_Hash_Map_int::const_iterator f = id_to_index.find (buffer);

    if (f == id_to_index.end ())
    {
      cerr << "Warning, subset member '" << buffer << "' not part of input\n";
      continue;
    }

    int j = (*f).second;

    if (in_set[j])
    {
      cerr << "Duplicate id in subset '" << buffer << "', ignored\n";
      continue;
    }

    items_in_set++;

    in_set[j] = 1;
  }

  if (0 == items_in_set || nr == items_in_set)
  {
    cerr << "IMpossible, cannot have " << items_in_set << " of " << nr << " items in a split\n";
    return 0;
  }

  if (verbose)
    std::cout << "Read " << items_in_set << " of " << nr << " items to be part of training set\n";

  items_in_set_acc.extra (items_in_set);

  return 1;
}

static int
identify_items_in_subset (const char * fname,
                          const IW_STL_Hash_Map_int & id_to_index,
                          int nr,
                          int * tmp)
{
  iwstring_data_source input (fname);
  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return identify_items_in_subset (input, id_to_index, nr, tmp);
}

static int
read_distance_matrix (IWDistanceMatrixBase<float> & dm,
                      const IWString & fname)
{
  return dm.do_read (fname);
}

static int
read_neighbour_data (NN_Results & nn_results,
                     iwstring_data_source & input)
{
  return nn_results.build_from_gfp_nearneighbours_dash_o (input);
}

static int
read_neighbour_data (NN_Results & nn_results,
                     const const_IWSubstring & fname)
{
  iwstring_data_source input (fname);
  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_neighbour_data (nn_results, input);
}

static int
test_train_split_evaluate (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vF:D:H:bn:t:x:c:s");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('c'))
  {
    if (! cl.value ('c', identifier_column) || identifier_column < 1)
    {
      cerr << "The identifier column specification (-c) must be a whole +ve number\n";
      usage (5);
    }

    if (verbose)
      cerr << "Identifiers in column " << identifier_column << endl;

    identifier_column--;
  }

  if (! cl.option_present ('F'))
  {
    cerr << "Must specify how the input was created with the -F option\n";
    usage (5);
  }

  if (cl.option_present ('F'))
  {
    const_IWSubstring f = cl.string_value ('F');

    if ("dm" == f)
    {
      input_is_distance_matrix = 1;
      if (verbose)
        cerr << "Input is a distance matrix\n";
    }
    else if ("gfp" == f)
    {
      input_from_gfp_nearneighbours_dash_o = 1;
      if (verbose)
        cerr << "Input from gfp_nearneighbours -o\n";
    }
    else
    {
      cerr << "Unrecognised -F qualifier '" << f << "'\n";
    }
  }

  if (cl.option_present('s'))
  {
    find_inter_split_similarities = 1;
    if (verbose)
      cerr << "Will compute inter-split similarities - bits in common\n";
  }

  if (cl.option_present ('x'))
  {
    if (! cl.value ('x', truncate_distances_at) || truncate_distances_at <= 0.0)
    {
      cerr << "The truncate distances at option (-x) must be a valid distance\n";
      usage (4);
    }

    if (verbose)
      cerr << "Will truncate all distances at " << truncate_distances_at << endl;
  }

  nsplits = cl.number_elements ();

  if (0 == nsplits)
  {
    cerr << "Insufficient arguments, must specify one or more subset files\n";
    usage (1);
  }

  if (find_inter_split_similarities)
  {
    fp = new IW_Bits_Base[nsplits];
    if (NULL == fp)
    {
      cerr << "Cannot allocate " << nsplits << " fingerprints\n";
      return 4;
    }
  }

  IWHistogram tmp_histogram;   // used for bucket info

  int overall_nbuckets = 0;

  if (cl.option_present ('H'))
  {
    histogram_initialisation_string = cl.string_value ('H');

    if (verbose)
      cerr << "Histogram initialised from " << histogram_initialisation_string << "'\n";

    if (! tmp_histogram.initialise (histogram_initialisation_string))
    {
      cerr << "Invalid histogram initialisation string '" << histogram_initialisation_string << "'\n";
      return 8;
    }

    overall_nbuckets = tmp_histogram.nbuckets ();

    overall_accumulator_count = new Accumulator_Int<int>[overall_nbuckets];
    overall_accumulator_fraction = new Accumulator<float>[overall_nbuckets];

    if (NULL == overall_accumulator_fraction)
    {
      cerr << "Bad news, could not allocate " << overall_nbuckets << " accumulator objects\n";
      return 5;
    }
  }

  if (cl.option_present ('b'))
  {
    write_histogram = 0;

    if (verbose)
      cerr << "Histogram not written\n";
  }

  if (cl.option_present ('n'))
  {
    if (! cl.value ('n', neighbours_to_consider) || neighbours_to_consider < 1)
    {
      cerr << "The number of neighbours to consider option (-n) must be a whole positive number\n";
      usage (4);
    }

    if (verbose)
      cerr << "Will consider on the " << neighbours_to_consider << " closest neighbours of each molecule\n";
  }

  if (cl.option_present ('t'))
  {
    int i = 0;
    const_IWSubstring t;
    while (cl.value ('t', t, i++))
    {
      int j = 0;
      const_IWSubstring token;
      while (t.nextword (token, j, ','))
      {
        double d;
        if (! token.numeric_value (d) || d < 0.0)
        {
          cerr << "Invalid distance '" << token << "' in -t qualifier\n";
          usage (5);
        }

        distances_to_check.add (d);
      }
    }
  }

  if (! cl.option_present ('D'))
  {
    cerr << "Must specify the distance matrix or near neighbour file via the -D option\n";
    usage (5);
  }

// Only one of these will actually be used

  IWDistanceMatrixBase<float> dm;

  NN_Results nn_results;

  if (cl.option_present ('D'))
  {
    const_IWSubstring fname = cl.string_value ('D');

    int rc;
    if (input_is_distance_matrix)
      rc = read_distance_matrix (dm, fname);
    else
      rc = read_neighbour_data (nn_results, fname);

    if (0 == rc)
    {
      cerr << "Cannot read input data '" << fname << "'\n";
      return 5;
    }
  }

  if (verbose > 2)
    nn_results.debug_print (cerr, verbose);

  IW_STL_Hash_Map_int id_to_index;
  
  int nr;    // the number of molecules in the set

  if (input_is_distance_matrix)
  {
    nr = dm.number_molecules ();
    for (int i = 0; i < nr; i++)
    {
      IWString tmp = dm.id (i);
      tmp.truncate_at_first (' ');

      int j = id_to_index.size ();
  
      id_to_index[tmp] = j;
    }
  }
  else
  {
    nr = nn_results.number_results ();
    for (int i = 0; i < nr; i++)
    {
      IWString tmp = nn_results[i].id ();
      tmp.truncate_at_first (' ');
  
      int j = id_to_index.size ();
  
      id_to_index[tmp] = j;
    }
  }

  if (0.0 == truncate_distances_at)
    ;
  else if (input_is_distance_matrix)
    do_truncate_distance_matrix (dm, truncate_distances_at);
  else
    do_truncate_neighbour_distances (nn_results, truncate_distances_at);

  int * tmp = new int[nr]; std::unique_ptr<int> free_tmp (tmp);

  for (int i = 0; i < nsplits; i++)
  {
    if (! identify_items_in_subset (cl[i], id_to_index, nr, tmp))
    {
      cerr << "Could not find subset members from '" << cl[i] << "'\n";
      return i + 1;
    }

    if (find_inter_split_similarities)
      fp[i].construct_from_array_of_ints(tmp, id_to_index.size());

    int rc;

    if (neighbours_to_consider && input_is_distance_matrix)
      rc = test_train_split_evaluate_N (dm, tmp, std::cout);
    else if (input_is_distance_matrix)
      rc = test_train_split_evaluate (dm, tmp, std::cout);
    else
      rc = test_train_split_evaluate (nn_results, tmp, std::cout);

    if (0 == rc)
    {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      return i;
    }
  }

  if (NULL != overall_accumulator_count)
  {
    int splits_examined = items_in_set_acc.n ();

    cerr << "Collected data on " << splits_examined << " splits\n";

    if (splits_examined)
    {
      cerr << "Between " << items_in_set_acc.minval () << " and " << items_in_set_acc.maxval ();
      if (items_in_set_acc.n () > 1)
        cerr << ", ave " << items_in_set_acc.average ();
      cerr << " items in sets\n";

      for (int i = 0; i < overall_nbuckets - 1; i++)
      {
        const Accumulator_Int<int> & c = overall_accumulator_count[i];
        const Accumulator<float> & f = overall_accumulator_fraction[i];
  
        float d = tmp_histogram.minval () + (i+1) * tmp_histogram.delta ();
  
        cerr << d << " min " << c.minval ();
  
        if (splits_examined > 1)
          cerr << " max " << c.maxval () << " ave " << c.average ();

        cerr << " fraction: minval " << f.minval ();
  
        if (splits_examined > 1)
          cerr << " max " << f.maxval () << " ave " << static_cast<float> (f.average ());

        cerr << endl;
      }
    }

    delete [] overall_accumulator_count;
    delete [] overall_accumulator_fraction;
  }

  if (find_inter_split_similarities)
    do_find_inter_split_similarities(cerr);

  if (verbose)
  {
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_train_split_evaluate (argc, argv);

  return rc;
}
