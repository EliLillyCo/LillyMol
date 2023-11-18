/*
  We have divided a set of molecules into two. We generate statistics
  on the nature of the split

  This variant converts the distance matrix to integers
*/

#include <stdlib.h>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/histogram/iwhistogram.h"
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "nn_results_.h"

const char * prog_name = nullptr;

static int verbose = 0;

static int input_is_distance_matrix = 0;

static int input_from_gfp_nearneighbours_dash_o = 0;

static IWString histogram_initialisation_string;

#define NPARTS 100

static float distance_minval = 0.0;

static float distance_maxval = 1.0;

static float distance_delta = distance_maxval / static_cast<float> (NPARTS);

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

static Accumulator_Int<int>   * overall_accumulator_count    = nullptr;
static Accumulator<float> * overall_accumulator_fraction = nullptr;

static Accumulator_Int<int> items_in_set_acc;

static int indices_to_check[NPARTS];

static IWString index_tag ("NBR<");

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  if (verbose)
    cerr << "$Id: test_train_split_evaluate.cc,v 1.2 2002/10/03 15:25:07 ian Exp $\n";
  cerr << "Analyses one or more test/train splits for inter-set distances\n";
  cerr << " -D <fname>     distance matrix or gfp_nearneighbours output\n";
  cerr << " -F <type>      what created the -D file\n";
  cerr << "    gfp         gfp_nearneighbours (use the -h and -o options)\n";
  cerr << "    dm          distance_matrix\n";
  cerr << " -t <dist>      one or more distances to check\n";
  cerr << " -x <dist>      group all distances beyond <dist> into the same group\n";
  cerr << " -n <number>    only consider the <number> closest neighbours across the split\n";
  cerr << " -T <dist>      specify largest possible distance\n";
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

static float
convert_int_back_to_distance (int i)
{
  return distance_minval + static_cast<float> (i) * distance_delta;
}

static int
convert_distance_to_index (float d)
{
  return static_cast<int> ((d - distance_minval) / distance_delta);
}

static int
report_accumulator (const int * zcount,
                    ostream & output)
{
  Accumulator_Int<int> acc;

  for (int i = 0; i < NPARTS; i++)
  {
    if (zcount[i])
      acc.extra (i, zcount[i]);
  }

  output << "N = " << acc.n ();
  output << " min " << convert_int_back_to_distance (acc.minval ());
  output << " ave " << static_cast<float> (acc.average_if_available_minval_if_not () / NPARTS);
  output << " max " << convert_int_back_to_distance (acc.maxval ());
  output << " tot " << static_cast<float> (acc.sum () / static_cast<float> (NPARTS) ) << endl;

  return output.good ();
}

/*
  Sometimes we don't actually do any output, we just update the overal_accumulator
  objects
*/

static int
write_results (const int * acc_set0,
               const int * acc_set1,
               int nacross,
               const int * acc_across,
               ostream & output)
{
  int do_writes = 1;

  if (nsplits > 0 && 0 == verbose)
    do_writes = 0;

  if (do_writes)
  {
    output << "Within set0\n";
    report_accumulator (acc_set0, output);
    output << "Within set1\n";
    report_accumulator (acc_set1, output);
    output << "Between sets";
    if (neighbours_to_consider)
      output << ", " << neighbours_to_consider << " neighbours";
    output << endl;
    report_accumulator (acc_across, output);
  }

  int j = 0;     // index into the array of accumulators

  unsigned int csum = 0;
  unsigned int sum_to_next_output = 0;

  for (int i = 0; i < NPARTS; i++)
  {
    sum_to_next_output += acc_across[i];

    if (0 == indices_to_check[i])
      continue;

    if (NULL != overall_accumulator_count)
    {
      overall_accumulator_count[j].extra (sum_to_next_output);
      overall_accumulator_fraction[j].extra (static_cast<float> (sum_to_next_output) / static_cast<float> (nacross));
      j++;
    }

    if (do_writes)
    {
      float d = convert_int_back_to_distance (i);

      csum += sum_to_next_output;

      output << d << ' ' << sum_to_next_output << ' ' << csum << endl;
    }

    sum_to_next_output = 0;
  }

  return output.good ();
}

/*
  Found that I got dramatically worse run-times if I didn't have an offset. Probably
  because some invalid floating point numbers were being represented. This arbitrary
  offset solves the problem
*/

#define NFPBAD_OFFSET -121213

static inline float
convert_to_int (float f)
{
  int i = static_cast<int> (f);

  int * iptr = reinterpret_cast<int *> (&f);

  *iptr = i + NFPBAD_OFFSET;

  return f;
}

static inline float
convert_to_scaled_value (float f)
{
  if (distance_maxval - f < distance_delta)
    return static_cast<float> (NPARTS - 1);

  f = (f - distance_minval) / distance_delta;

  assert (f >= 0.0 && f < NPARTS);

  return f;
}

static int
convert_distance_matrix_to_int (const IWDistanceMatrixBase<float> & dfrom,
                                IWDistanceMatrixBase<int> & dto)
{
  int n = dfrom.number_molecules ();

  if (! dto.resize (n))
  {
    cerr << "Memory failure, cannot resize in array to " << n << " size\n";
    return 0;
  }

  for (int i = 0; i < n; i++)
  {
    for (int j = i + 1; j < n; j++)
    {
      float d = dfrom.zvalue (i, j);

      d = convert_to_scaled_value (d);

      dto.set (i, j, static_cast<int> (d));
    }
  }

  return 1;
}

static void
convert_neighbour_distances_to_int (NN_Results & nn_results)
{
  int n = nn_results.number_results ();

  for (int i = 0; i < n; i++)
  {
    NN_Item & ni = nn_results[i];

    for (int j = 0; j < ni.number_elements (); j++)
    {
      ID_Dist * nj = ni[j];

      float d = nj->distance ();

      d = convert_to_scaled_value (d);

      nj->set_distance (convert_to_int (d));
    }
  }

  return;
}

/*
  We are storing int's in float representations!
*/

static inline int
get_int_value_stored_in_float (float f)
{
  const int * i = reinterpret_cast<const int *> (&f);

  return *i - NFPBAD_OFFSET;
}

/*
  Special case where we are considering only the N nearest neighbours across the divide
*/

static int
test_train_split_evaluate_N (const IWDistanceMatrixBase<int> & dm,
                             const int * in_set,
                             ostream & output)
                                
{
  int acc_across[NPARTS];
  set_vector (acc_across, NPARTS, 0);
  int acc_set0[NPARTS];
  set_vector (acc_set0, NPARTS, 0);
  int acc_set1[NPARTS];
  set_vector (acc_set1, NPARTS, 0);

  resizable_array<int> nbr;    // scope here for efficiency

  int nr = dm.number_molecules ();

  int nacross = 0;

  for (int i = 0; i < nr; i++)
  {
    int jstop = dm.neighbours_sorted_by_distance (i, nbr);

    int nacross_for_i = 0;

    for (int j = 0; j < jstop; j++)
    {
      int k = nbr[j];
      int d = dm.zvalue (i, k);

      assert (d >= 0 && d < NPARTS);

      if (0 == in_set[i] && 0 == in_set[k])
        acc_set0[d]++;
      else if (1 == in_set[i] && 1 == in_set[k])
        acc_set1[d]++;
      else
      {
        if (nacross_for_i <= neighbours_to_consider)
        {
//        cerr << "Items " << i << " and " << k << " at distance " << nj->distance () << endl;
          acc_across[d]++;
          nacross_for_i++;
          nacross++;
        }
      }
    }
  }

  return write_results (acc_set0, acc_set1, nacross, acc_across, output);
}

static int
test_train_split_evaluate (const IWDistanceMatrixBase<int> & dm,
                           const int * in_set,
                           ostream & output)
                                
{
  int acc_across[NPARTS];
  set_vector (acc_across, NPARTS, 0);
  int acc_set0[NPARTS];
  set_vector (acc_set0, NPARTS, 0);
  int acc_set1[NPARTS];
  set_vector (acc_set1, NPARTS, 0);

  int nr = dm.number_molecules ();

  int nacross = 0;

  for (int i = 0; i < nr; i++)
  {
    for (int j = i + 1; j < nr; j++)
    {
      int d = dm.zvalue (i, j);

      if (0 == in_set[i] && 0 == in_set[j])
        acc_set0[d]++;
      else if (1 == in_set[i] && 1 == in_set[j])
        acc_set1[d]++;
      else if (neighbours_to_consider > 0 && nacross > neighbours_to_consider)
        ;
      else
      {
//      cerr << "Items " << i << " and " << k << " at distance " << nj->distance () << endl;
        acc_across[d]++;
        nacross++;
      }
    }
  }

  return write_results (acc_set0, acc_set1, nacross, acc_across, output);
}

static int
test_train_split_evaluate_N (const NN_Results & nn_results,
                             const int * in_set,
                             ostream & output)
{
  assert (neighbours_to_consider > 0);

  int acc_across[NPARTS];
  set_vector (acc_across, NPARTS, 0);
  int acc_set0[NPARTS];
  set_vector (acc_set0, NPARTS, 0);
  int acc_set1[NPARTS];
  set_vector (acc_set1, NPARTS, 0);

  int nacross = 0;

  int nr = nn_results.number_results ();

  for (int i = 0; i < nr; i++)
  {
    const NN_Item & ni = nn_results[i];

    int nacross_for_i = 0;

    for (int j = 0; j < ni.number_elements (); j++)
    {
      const ID_Dist * nj = ni[j];

      int k = nj->index_in_pool ();

      assert (k >= 0 && k < nr);

      int d = get_int_value_stored_in_float (nj->distance ());

      if (0 == in_set[i] && 0 == in_set[k])
        acc_set0[d]++;
      else if (1 == in_set[i] && 1 == in_set[k])
        acc_set1[d]++;
      else if (nacross_for_i < neighbours_to_consider)
      {
//      cerr << "Items " << i << " and " << k << " at distance " << nj->distance () << endl;
        acc_across[d]++;
        nacross_for_i++;
        nacross++;
      }
    }
  }

  return write_results (acc_set0, acc_set1, nacross, acc_across, output);
}


static int
test_train_split_evaluate (const NN_Results & nn_results,
                           const int * in_set,
                           ostream & output)
{
  int acc_across[NPARTS];
  set_vector (acc_across, NPARTS, 0);
  int acc_set0[NPARTS];
  set_vector (acc_set0, NPARTS, 0);
  int acc_set1[NPARTS];
  set_vector (acc_set1, NPARTS, 0);

  int nacross = 0;

  int nr = nn_results.number_results ();

  for (int i = 0; i < nr; i++)
  {
    const NN_Item & ni = nn_results[i];

    for (int j = 0; j < ni.number_elements (); j++)
    {
      const ID_Dist * nj = ni[j];

      int k = nj->index_in_pool ();

      assert (k >= 0 && k < nr);

      int d = get_int_value_stored_in_float (nj->distance ());

      if (0 == in_set[i] && 0 == in_set[k])
        acc_set0[d]++;
      else if (1 == in_set[i] && 1 == in_set[k])
        acc_set1[d]++;
      else
      {
//      cerr << "Items " << i << " and " << k << " at distance " << nj->distance () << endl;
        acc_across[d]++;
        nacross++;
      }
    }
  }

  return write_results (acc_set0, acc_set1, nacross, acc_across, output);
}

static float
largest_distance_in_neighbour_lists (const NN_Results & nn_results)
{
  int nr = nn_results.number_results ();

  ID_Dist::distance_t maxd = 0.0;

  for (int i = 0; i < nr; i++)
  {
    const NN_Item & ni = nn_results[i];

    for (int j = 0; j < ni.number_elements (); j++)   // could just check last neighbour, but what if neighbours are out of order
    {
      ID_Dist::distance_t d = ni[j]->distance ();
      if (d > maxd)
        maxd = d;
    }
  }

  return static_cast<float> (maxd);
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
float_comparitor (const float * f1, const float * f2)
{
  if (*f1 < *f2)
    return -1;

  if (*f1 > *f2)
    return 1;

  return 0;
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
  Command_Line cl (argc, argv, "vF:D:H:n:t:T:x:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

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

      cerr << "Note that numbers will be about double what they are compared with a distance matrix\n";
    }
    else
    {
      cerr << "Unrecognised -F qualifier '" << f << "'\n";
    }
  }

  float truncate_distances_at = 0.0;

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

  if (! cl.option_present ('D'))
  {
    cerr << "Must specify the distance matrix or near neighbour file via the -D option\n";
    usage (5);
  }

  if (! cl.option_present ('t'))
  {
    cerr << "Must specify one or more distances to check via the -t option\n";
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
    else if (input_from_gfp_nearneighbours_dash_o)
      rc = read_neighbour_data (nn_results, fname);
    else if (fname.ends_with (".nn"))
    {
      rc = read_neighbour_data (nn_results, fname);
      input_is_distance_matrix = 1;
    }
    else
    {
      cerr << "NOt sure what to do with '" << fname << "' use the -F option to specify what type\n";
      usage (5);
      rc = 99;
    }

    if (0 == rc)
    {
      cerr << "Cannot read input data '" << fname << "'\n";
      return 5;
    }
  }

  if (verbose)
  {
    float maxd;

    if (input_is_distance_matrix)
      maxd = dm.maxval ();
    else
      maxd = largest_distance_in_neighbour_lists (nn_results);

    cerr << "Largest distance in input " << maxd << endl;
  }

  if (cl.option_present ('T'))
  {
    if (! cl.value ('T', distance_maxval) || distance_maxval <= distance_minval)
    {
      cerr << "Invalid maximum distance (-T option) must be a valid distance\n";
      usage (11);
    }

    if (verbose)
      cerr << "Maximum distance set to " << distance_maxval << endl;

    distance_delta = (distance_maxval - distance_minval) / static_cast<float> (NPARTS);
  }

  resizable_array<float> distances_to_check;

  if (cl.option_present ('t'))
  {
    set_vector (indices_to_check, NPARTS, 0);

    int i = 0;
    const_IWSubstring t;
    while (cl.value ('t', t, i++))
    {
      int j = 0;
      const_IWSubstring token;
      while (t.nextword (token, j, ','))
      {
        float d;
        if (! token.numeric_value (d) || d < distance_minval || d > distance_maxval)
        {
          cerr << "Invalid distance '" << token << "' in -t qualifier, must be between " << distance_minval << " and " << distance_maxval << endl;
          usage (5);
        }

        int k = convert_distance_to_index (d);

        indices_to_check[k] = 1;

        if (verbose)
          cerr << "Will report distance " << d << " index " << k << endl;

        distances_to_check.add_if_not_already_present (d);
      }
    }

    distances_to_check.sort (float_comparitor);

    int nd = distances_to_check.number_elements ();

    assert (nd > 0);

    overall_accumulator_count    = new Accumulator_Int<int>[nd];
    overall_accumulator_fraction = new Accumulator<float>[nd];

    if (NULL == overall_accumulator_fraction)
    {
      cerr << "Bad news, could not allocate " << nd << " accumulator objects\n";
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

  IWDistanceMatrixBase<int> intdm;

  if (input_is_distance_matrix)
    convert_distance_matrix_to_int (dm, intdm);
  else
    convert_neighbour_distances_to_int (nn_results);

  int * tmp = new int[nr];

  for (int i = 0; i < nsplits; i++)
  {
    if (! identify_items_in_subset (cl[i], id_to_index, nr, tmp))
    {
      cerr << "Could not find subset members from '" << cl[i] << "'\n";
      return i + 1;
    }

    int rc;

    if (neighbours_to_consider && input_is_distance_matrix)
      rc = test_train_split_evaluate_N (intdm, tmp, std::cout);
    else if (input_is_distance_matrix)
      rc = test_train_split_evaluate (intdm, tmp, std::cout);
    else if (neighbours_to_consider)
      rc = test_train_split_evaluate_N (nn_results, tmp, std::cout);
    else
      rc = test_train_split_evaluate (nn_results, tmp, std::cout);

    if (0 == rc)
    {
      cerr << "Fatal error processing '" << cl[i] << "'\n";
      return i;
    }
  }

  delete [] tmp;

  if (NULL != overall_accumulator_count)
  {
    int splits_examined = items_in_set_acc.n ();

    cerr << "Collected data on " << splits_examined << " splits\n";

    if (splits_examined)
    {
      cerr << "Splits between " << items_in_set_acc.minval () << " and " << items_in_set_acc.maxval ();
      if (items_in_set_acc.n () > 1)
        cerr << ", ave " << items_in_set_acc.average ();
      cerr << " items in set\n";

      for (int i = 0; i < distances_to_check.number_elements (); i++)
      {
        const Accumulator_Int<int> & c = overall_accumulator_count[i];
        const Accumulator<float> & f = overall_accumulator_fraction[i];
  
        cerr << distances_to_check[i] << " min " << c.minval ();
  
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
