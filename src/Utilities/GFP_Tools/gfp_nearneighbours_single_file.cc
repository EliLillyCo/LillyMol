/*
  Finds near neighbours within a single fingerprint file
*/

#include <stdlib.h>
#include <memory>
#include <fstream>
#include <limits>
#include <cmath>

#include "cmdline.h"
#include "accumulator.h"
#include "iwstring_data_source.h"
#include "iw_tdt.h"
#include "iwhistogram.h"
#include "iwdigits.h"
#include "report_progress.h"

#include "gfp.h"
#include "sparse_collection.h"
#include "tversky.h"

static int verbose = 0;

static int neighbours_to_find = 1;

static int write_neighbours_as_index_numbers = 0;

/*
  If we have a filter active, we may not get enough neighbours. We can rescan the
  pool with all the exclusions turned off
*/

static int rescan_if_not_enough_neighbours = 0;

static int molecules_rescanned = 0;

static similarity_type_t abandon_distance_threshold = -1.0;

static Accumulator<float> distance_stats;

static Fraction_as_String fraction_as_string;

/*
  Sometimes it is useful to know the statistics of the nearest neighbour of each molecule
*/

static Accumulator<similarity_type_t> nearest_neighbour_distance_stats;

static IWHistogram histogram_nearnest_neighbour_distances;

static int create_histogram = 0;

static int write_minimal_histogram = 0;

static extending_resizable_array<int> neighbour_count;

static int fingerprints_read = 0;

static Tversky tversky;

static IWString smiles_tag("$SMI<");

/*
  The identifier tag used in each TDT
*/

static IWString identifier_tag("PCN<");

static IWString number_neighbours_tag;

/*
  If we are writing neighbours as index numbers
*/

static IWString neighbour_tag("NBR<");

/*
  May 99. For each input TDT, I need to know the average distance
  of the neighbours within the pool
*/

static IWString tag_for_average_distance;

static const_IWSubstring distance_tag("DIST<");

/*
  We ignore distances longer than DISTANCE_THRESHOLD
*/

static similarity_type_t upper_distance_threshold = std::numeric_limits<float>::max();

static similarity_type_t lower_distance_threshold = -1.0;

static int allow_arbitrary_distances = 0;

/*
  When we have thresholds, we may choose to not write molecules with no neighbours
*/

static int write_molecules_with_no_neighbours = 1;

static int molecules_with_no_neighbours = 0;

static int write_smiles = 1;

static Report_Progress report_progress;

/*
  We often do a run where we find the neighbours of a set of molecules within itself.
  We may or may not want to find the molecule itself as its own nearest neighbour
*/

static int do_not_compare_molecules_with_themselves = 1;   // Oct 2015, now the default

/*
  When doing near neighbour determinations within a single set of
  molecules, we gain great efficiencies by writing the index of the
  item rather than its name.  This helps programmes that read the nn
  file
*/

class IW_GFP_D_ID : public IW_GFP_D
{
  private:
    int _ndx;
    IWString _smiles;

  public:
    IW_GFP_D_ID();

    void set_index (int n) { _ndx = n;}
    int  index_in_pool() const { return _ndx;}

    int write_smiles_and_id (IWString &) const;

    void set_smiles (const const_IWSubstring & s) { _smiles = s;}
    const IWString & smiles() const { return _smiles;}
};

IW_GFP_D_ID::IW_GFP_D_ID()
{
  _ndx = -1;
}

int
IW_GFP_D_ID::write_smiles_and_id (IWString & output) const
{
  if (_smiles.length())
    output << smiles_tag << _smiles << ">\n";

  output << identifier_tag << IW_General_Fingerprint::id() << ">\n";

  return 1;
}

/*
  Our pool is an array of FP objects
*/

static IW_GFP_D_ID * pool = NULL;

static int pool_size = 0;

static int
build_pool (iwstring_data_source & input)
{
  int items_in_pool = 0;

  int tdts_read = 0;

  IW_TDT tdt;
  while (tdt.next(input))
  {
    tdts_read++;

    int fatal;
    if (! pool[items_in_pool].construct_from_tdt(tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      continue;
    }

    pool[items_in_pool].set_index(items_in_pool);

    if (write_smiles)
    {
      const_IWSubstring smi;
      if (! tdt.dataitem_value(smiles_tag, smi))
      {
        cerr << "Cannot extract smiles\n";
        cerr << tdt;
        return 0;
      }

      pool[items_in_pool].set_smiles(smi);
    }

    items_in_pool++;

    if (items_in_pool == pool_size)
    {
      if (verbose)
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
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size)
  {
    pool_size = input.count_records_starting_with(identifier_tag);

    if (0 == pool_size)
    {
      cerr << "No occurrences of " << identifier_tag << "' in input\n";
      return 0;
    }

    pool = new IW_GFP_D_ID[pool_size];
    if (NULL == pool)
    {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    if (verbose)
      cerr << "Pool automatically sized to " << pool_size << endl;
  }

  return build_pool(input);
}

static int
write_average_neighbour_distance(IW_GFP_D_ID ** neighbours,
                      int number_neighbours,
                      IWString & output)
{
  Accumulator<similarity_type_t> acc;

  for (int i = 0; i < number_neighbours; i++)
  {
    IW_GFP_D_ID * n = neighbours[i];
    assert (NULL != neighbours[i]);

    acc.extra(n->distance());
  }

  output << tag_for_average_distance << acc.average_if_available_minval_if_not() << ">\n";

  return output.good();
}


static int
write_neighbour_list (const IW_GFP_D_ID & target,
                      IW_GFP_D_ID ** neighbours,
                      int number_neighbours,
                      IWString & output)
{
  target.write_smiles_and_id(output);

  if (number_neighbours && tag_for_average_distance.length())
    write_average_neighbour_distance(neighbours, number_neighbours, output);

  if (number_neighbours_tag.length())
    output << number_neighbours_tag << number_neighbours << ">\n";

  for (int i = 0; i < number_neighbours; i++)
  {
    IW_GFP_D_ID * n = neighbours[i];
    assert (NULL != neighbours[i]);

    if (write_smiles)
      output << smiles_tag << n->smiles() << ">\n";

    if (write_neighbours_as_index_numbers)
      output << neighbour_tag << n->index_in_pool() << ">\n";
    else
      output << identifier_tag << n->id() << ">\n";

    similarity_type_t d = n->distance();

    output << distance_tag;
    output.append_number(d, 3);
    output << ">\n";

    distance_stats.extra(d);

    if (allow_arbitrary_distances)
      ;
    else if (0 == i && (create_histogram || verbose))
    {
      nearest_neighbour_distance_stats.extra(d);
      if (create_histogram)
        histogram_nearnest_neighbour_distances.extra(d);
    }
  }

  output << "|\n";

  return 1;
}

//#define DEBUG_INSERTION

/*
  Insert NEIGHBOUR into the neighbour list of maximum size NEIGHBOURS_TO_FIND.
  NEIGHBOUR becomes the new item WHERE
*/

static void
insert_in_neighbour_list (IW_GFP_D_ID ** neighbours, int neighbours_to_find, 
                          int & neighbours_found,
                          IW_GFP_D_ID & neighbour, int where)
{
  assert (neighbours_found <= neighbours_to_find);
  assert (where < neighbours_to_find);

#ifdef DEBUG_INSERTION
  cerr << "Inserting " << neighbour.distance() << ". To find = " << neighbours_to_find << " found = " << neighbours_found << " where = " << where << endl;
#endif

  if (0 == neighbours_found)     // list is empty, easy...
  {
    neighbours[0] = &neighbour;
    neighbours_found = 1;
    return;
  }

// items already in the list, shift right to make room for new member

  int istart;
  if (neighbours_found == neighbours_to_find)     // list is full
    istart = neighbours_found - 1;
  else
  {
    istart = neighbours_found;
    neighbours_found++;
  }

  for (int i = istart; i > where; i--)
  {
    neighbours[i] = neighbours[i - 1];
  }

  neighbours[where] = &neighbour;

//#define CHECK_INSERTION
#ifdef CHECK_INSERTION
  int failure = 0;

  for (int i = 1; i < neighbours_found; i++)
  {
    if (neighbours[i - 1]->distance() > neighbours[i]->distance())
    {
      cerr << "Sort/insertion failed, out of order, i = " << i << ' ' << neighbours[i - ]->distance() << " vs " << neighbours[i]->distance() << endl;
      failure++;
    }
  }

#ifdef DEBUG_INSERTION
  cerr << "After insertion, Neighbours found = " << neighbours_found << endl;
  for (int i = 0; i < neighbours_found; i++)
  {
    cerr << "i = " << i << " distance " << neighbours[i]->distance() << endl;
  }
#endif

  if (failure)
  {
    for (int i = 0; i < neighbours_found; i++)
    {
      cerr << "i = " << i << " distance " << neighbours[i]->distance() << endl;
    }

    exit(87);
  }

#endif

  return;
}

/*
  We need to all NEIGHBOUR to a list of neighbours
*/

static void
do_neighbour_list_insertion (IW_GFP_D_ID ** neighbours, int neighbours_to_find_this_fingerprint, 
                             int & neighbours_found,
                             similarity_type_t t,
                             IW_GFP_D_ID & neighbour)
{

// Can we put it at the head of the list

  if (0 == neighbours_found || t <= neighbours[0]->distance())
  {
    insert_in_neighbour_list(neighbours, neighbours_to_find_this_fingerprint, neighbours_found, neighbour, 0);
    return;
  }

// Does it go at the end of the list

  if (t >= neighbours[neighbours_found - 1]->distance())
  {
    if (neighbours_found == neighbours_to_find_this_fingerprint)     // list is full
      return;

    insert_in_neighbour_list(neighbours, neighbours_to_find_this_fingerprint, neighbours_found, neighbour, neighbours_found);
    return;
  }

// It goes somewhere in the list

#ifdef DEBUG_INSERTION
    cerr << "Inserting " << t << " into list of " << neighbours_found << " neighbours\n";

    for (int j = 0; j < neighbours_found; j++)
    {
      cerr << "j = " << j << " dist " << neighbours[j]->distance();
      if (neighbours[j]->distance() < t)
        cerr << endl;
      else
        cerr << " *\n";
    }
#endif

    int left = 0;    // loop needs to set left
    int right = neighbours_found - 1;
    int middle = (left + right) / 2;
    while (middle > left)
    {
      similarity_type_t m = neighbours[middle]->distance();
#ifdef DEBUG_INSERTION
      cerr << "left " << left << " middle " << middle << " right " << right << endl;
      cerr << neighbours[left]->distance() << ',' << neighbours[middle]->distance() << ',' << neighbours[right]->distance() << endl;
#endif

      if (t < m)
        right = middle;
      else if (t > m)
        left = middle;
      else    // they are equal, will insert at middle
        break;

      middle = (left + right) / 2;
    }

    insert_in_neighbour_list(neighbours, neighbours_to_find_this_fingerprint, neighbours_found, neighbour, middle + 1);

  return;
}

/*
  Common function for doing the distance computation
*/

similarity_type_t
compute_the_distance (IW_GFP_D_ID & fp1, IW_GFP_D_ID & fp2,
                      const Tversky & tversky)
{
  if (tversky.active())
    return static_cast<similarity_type_t>(1.0) - fp1.IW_General_Fingerprint::tversky(fp2, tversky);

  return static_cast<similarity_type_t>(1.0) - fp1.tanimoto(fp2);
}

static int
three_column_output_all_pairs (IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < pool_size; i++)
  {
    IW_GFP_D_ID & fpi = pool[i];

    for (int j = i + 1; j < pool_size; j++)
    {
      similarity_type_t d = compute_the_distance(fpi, pool[j], tversky);

      if (d < lower_distance_threshold || d > upper_distance_threshold)
        continue;

      output << fpi.id() << ' ' << pool[j].id() << ' ' << d << '\n';
      output.write_if_buffer_holds_more_than(32768);
    }
  }

  output.flush();

  return 1;
}

/*
  Sometimes we need to get at least some neighbours. due to things like atom count
  windows or other factors, we may end up with no neighbours
*/

static int
do_rescan_for_not_enough_neighbours (IW_GFP_D_ID & fp,
                IW_GFP_D_ID ** neighbours,
                int neighbours_to_find_this_fingerprint)
{
  molecules_rescanned++;

  int neighbours_found = 0;

  for (int i = 0; i < pool_size; i++)
  {
    similarity_type_t t = compute_the_distance(fp, pool[i], tversky);

    if (0.0 == t && do_not_compare_molecules_with_themselves && fp.id() == pool[i].id())   // delay expensive and rare string comparison till here
      continue;

    pool[i].set_distance(t);

    do_neighbour_list_insertion(neighbours, neighbours_to_find_this_fingerprint, neighbours_found, t, pool[i]);
  }

  return neighbours_found;
}

static int
nearneighbours(IW_GFP_D_ID & fp,
               IW_GFP_D_ID ** neighbours,
               int neighbours_to_find_this_fingerprint)
{
  int neighbours_found = 0;
  assert (neighbours_to_find_this_fingerprint > 0);

  for (int i = 0; i < neighbours_to_find_this_fingerprint; i++)
  {
    neighbours[i] = NULL;
  }

  neighbours_found = 0;

  for (int i = 0; i < pool_size; i++)
  {
    if (! can_be_compared(fp, pool[i]))
      continue;

    similarity_type_t t;
    if (abandon_distance_threshold > static_cast<similarity_type_t>(0.0))
    {
      if (! fp.IW_General_Fingerprint::tanimoto(pool[i], abandon_distance_threshold, t))
        continue;

      t = static_cast<similarity_type_t>(1.0) - t;
    }
    else 
      t = compute_the_distance(fp, pool[i], tversky);

//#define DEBUG_NN
#ifdef DEBUG_NN
    cerr << "Distance between '" << fp.id() << " and pool " << i << " '" << pool[i].id() << "' is " << t << endl;
#endif

    if (t < static_cast<similarity_type_t>(0.0))
    {
      cerr << "INvalid similarity " << t << " i = " << i << endl;
      abort();
    }

    if (allow_arbitrary_distances)
      ;
    else if (t > static_cast<similarity_type_t>(1.0))
    {
      cerr << "Distance out of range " << t << " '" << fp.id() << "' to '" << pool[i].id() << "'\n";
      abort();
    }

    if (t < upper_distance_threshold && t > lower_distance_threshold)   // in range
      ;
    else
      continue;

    if (neighbours_found == neighbours_to_find_this_fingerprint && t >= neighbours[neighbours_found - 1]->distance())
      continue;

    if (static_cast<similarity_type_t>(0.0) == t && do_not_compare_molecules_with_themselves && fp.id() == pool[i].id())   // delay expensive and rare string comparison till here
      continue;

    pool[i].set_distance(t);

    do_neighbour_list_insertion(neighbours, neighbours_to_find_this_fingerprint, neighbours_found, t, pool[i]);
  }

  if (neighbours_found < rescan_if_not_enough_neighbours)
    return do_rescan_for_not_enough_neighbours(fp, neighbours, rescan_if_not_enough_neighbours);

  return neighbours_found;
}

static int
nearneighbours(IWString_and_File_Descriptor & output)
{
  int nbrs = neighbours_to_find;
  if (std::numeric_limits<float>::max() != upper_distance_threshold)
    nbrs = pool_size;
  else if (rescan_if_not_enough_neighbours > nbrs)
    nbrs = rescan_if_not_enough_neighbours;

  IW_GFP_D_ID ** neighbours = new IW_GFP_D_ID * [nbrs]; std::unique_ptr<IW_GFP_D_ID *[]> free_neighbours(neighbours);

  for (int i = 0; i < pool_size; i++)
  {
    int nn = nearneighbours(pool[i], neighbours, nbrs);

    neighbour_count[nn]++;

    if (0 == nn && ! write_molecules_with_no_neighbours)
      molecules_with_no_neighbours++;
    else
    {
      (void) write_neighbour_list(pool[i], neighbours, nn, output);
      output.write_if_buffer_holds_more_than(32768);
    }

    if (report_progress())
      cerr << "Processed " << i << " computations\n";
  }

  output.flush();

  return output.good();
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;

  cerr << "Finds near neighbours of a set of fingerprints\n";
  cerr << "Usage <options> <input_file>\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -n <number>      specify how many neighbours to find\n";
  cerr << " -t <dis>         discard distances shorter than <dis>\n";
  cerr << " -T <dis>         discard distances longer than <dis>\n";
  cerr << " -z               don't write molecules with no neighbours\n";
  cerr << " -I <tag>         specify identifier dataitem (default '" << identifier_tag << ")\n";
  cerr << " -A <TAG>         write average neighbour distance to <TAG>\n";
  cerr << " -X <distance>    abandon distance computation if any component > distance\n";
  cerr << " -r <number>      ensure that all molecules have at least <number> neighbours\n";
//cerr << " -h               discard neighbours with zero distance and the same ID as the target\n";
  cerr << " -o               cross referencing a single file. Write neighbours as index numbers\n";
  cerr << " -H <fname>       write histogram of closest distances to <fname>\n";
  cerr << " -b               write minimal histogram data - two columns\n";
  cerr << " -F ...           gfp options, enter '-F help' for details\n";
  cerr << " -V ...           Tversky specification, enter '-V help' for details\n";
  cerr << " -K ...           options for converting sparse fingerprints to fixed\n";
  cerr << " -N <tag>         write number neighbours as <tag>\n";
  cerr << " -p               write all pair-wise distances in 3 column form\n";
  cerr << " -j <precision>   output precision for distances\n";
  cerr << " -x               exclude smiles from the output\n";
  cerr << " -y               allow arbitrary distances\n";
  cerr << " -R <number>      report progress every <number> items processed\n";
  cerr << " -v               verbose output\n";

  exit(rc);
}

static int
do_write_histogram(std::ostream & os)
{
  if (write_minimal_histogram)
    return histogram_nearnest_neighbour_distances.write_terse(os);
  else
    return histogram_nearnest_neighbour_distances.write(os);
}

static int
nearneighbours(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vs:n:I:t:T:P:F:W:Q:X:V:A:zr:hoK:xH:N:bypj:R:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('V'))
  {
    if (! tversky.parse_command_line(cl, 'V', verbose))
    {
      cerr << "Cannot parse Tversky specifications\n";
      usage(4);
    }

    cerr << "Tversky parameters " << tversky.a() << " and " << tversky.b() << endl;
  }

  if (cl.option_present('I'))
  {
    (void) cl.value('I', identifier_tag);

    set_identifier_tag(identifier_tag);

    if (verbose)
      cerr << "Identifiers tagged as '" << identifier_tag << "'\n";
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', rescan_if_not_enough_neighbours) || rescan_if_not_enough_neighbours < 1)
    {
      cerr << "Invalid value for min number of neighbours (-r option)\n";
      usage(13);
    }

    if (verbose)
      cerr << "Molecules will always have at least " << rescan_if_not_enough_neighbours << " neighbours\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  Set_of_Sparse_Fingerprint_Collection_Profile sfcp;

  if (cl.option_present('K'))
  {
    if (! parse_sparse_to_dense_fingerprint_specifications(cl, 'K', verbose))
    {
      cerr << "Invalid sparse->fixed specification(s) (-K option)\n";
      return 4;
    }
  }

  if (need_to_call_initialise_fingerprints(cl))
  {
    if (! initialise_fingerprints(cl, verbose))
    {
      cerr << "Cannot initialise general fingerprint options\n";
      usage(17);
    }
  }
  else if (! initialise_fingerprints(cl[0], verbose))
  {
    cerr << "Cannot initialise fingerprints from '" << cl[0] << "'\n";
    return 11;
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 1)
    {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }

    pool = new IW_GFP_D_ID[pool_size];
    if (NULL == pool)
    {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    if (verbose)
      cerr << "system sized to " << pool_size << endl;
  }

  if (! build_pool(cl[0]))
  {
    cerr << "Cannot build pool from '" << cl[0] << "'\n";
    return 5;
  }

  if (0 == pool_size)
  {
    cerr << "No fingerprints\n";
    return 5;
  }

  if (cl.option_present('y'))
  {
    allow_arbitrary_distances = 1;

    if (verbose)
      cerr << "Distances not constrained to [0-1]\n";
  }

  if (cl.option_present('R'))
  {
    if (! report_progress.initialise(cl, 'R', verbose))
    {
      cerr << "The -R option must be followed by a positive whole number\n";
      usage(18);
    }
  }

  std::ofstream stream_for_nearest_neighbour_histogram;

  if (cl.option_present('H'))
  {
    if (allow_arbitrary_distances)
    {
      cerr << "Sorry the histogram option is inconsistent with the -y option\n";
      usage(4);
    }

    create_histogram = 1;

    if (cl.option_present('b'))
    {
      write_minimal_histogram = 1;
      if (verbose)
        cerr << "Only a minimal histogram will be written\n";
    }

    const char * h = cl.option_value('H');
    stream_for_nearest_neighbour_histogram.open(h, std::ios::out);
    if (! stream_for_nearest_neighbour_histogram.good())
    {
      cerr << "Sorry, cannot open histogram stream '" << h << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Histogram data written to '" << h << "'\n";
  }

  if (cl.option_present('N'))
  {
    cl.value('N', number_neighbours_tag);

    if (verbose)
      cerr << "The number of neighbours written as '" << number_neighbours_tag << "'\n";

    if (! number_neighbours_tag.ends_with('<'))
      number_neighbours_tag << '<';
  }

  if (cl.option_present('K'))
  {
    for (int i = 0; i < pool_size; i++)
    {
      sfcp.build_profile(pool[i]);
    }

    sfcp.finished_profiling(verbose);

    if (verbose)
      sfcp.report(cerr);

    for (int i = 0; i < pool_size; i++)
    {
      pool[i].convert_to_non_sparse_forms(sfcp);
    }
  }

  if (cl.option_present('X'))
  {
    if (! cl.value('X', abandon_distance_threshold) || abandon_distance_threshold < 0.0 || abandon_distance_threshold > 1.0)
    {
      cerr << "The -X option must be followed by a valid distance (0.0, 1.0)\n";
      usage(13);
    }

    if (verbose)
      cerr << "Distance compuations abandoned if any component > " << abandon_distance_threshold << endl;
  }

// process the -o option here, after the -e but before the -E option

  if (cl.option_present('o'))
  {
    write_neighbours_as_index_numbers = 1;
    write_smiles = 0;

    if (verbose)
      cerr << "Neighbours written as index numbers\n";
  }

  if (cl.option_present('x'))
  {
    write_smiles = 0;

    if (verbose)
      cerr << "Smiles excluded from output\n";
  }

  if (cl.option_present('A'))
  {
    tag_for_average_distance = cl.string_value('A');

    if (verbose)
      cerr << "The average neighbour distance will be written to '" << tag_for_average_distance << "'\n";

    if (! tag_for_average_distance.ends_with('<'))
      tag_for_average_distance += '<';
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', lower_distance_threshold) || lower_distance_threshold < 0.0 || lower_distance_threshold > 1.0)
    {
      cerr << "The -t option must be followed by a valid similarity value\n";
      usage(12);
    }

    if (verbose)
      cerr << "Lower distance threshold set to " << lower_distance_threshold << endl;
  }

  if (cl.option_present('T'))
  {
    if (! cl.value('T', upper_distance_threshold) || upper_distance_threshold < 0.0 || upper_distance_threshold > 1.0)
    {
      cerr << "The -T option must be followed by a valid similarity value\n";
      usage(12);
    }

    if (verbose)
      cerr << "Upper distance threshold set to " << upper_distance_threshold << endl;
  }

  if (cl.option_present('z') && ! cl.option_present('t') && ! cl.option_present('T'))
  {
    cerr << "The don't write molecules with no neighbours option (-z) only makes sense with thresholds\n";
    usage(13);
  }

  if (cl.option_present('z'))
  {
    write_molecules_with_no_neighbours = 0;
    if (verbose)
      cerr << "Will not write molecules with no neighbours\n";
  }

  if (cl.option_present('n'))
  {
    const_IWSubstring nvalue = cl.string_value('n');
    if ("all" == nvalue)
    {
      neighbours_to_find = pool_size;

      if (verbose)
        cerr << "May get as many as " << pool_size << " neighbours\n";
    }
    else
    {
      if (! nvalue.numeric_value(neighbours_to_find) || neighbours_to_find < 1)
      {
        cerr << "Invalid neighbours to find specifier '" << nvalue << "'\n";
        usage(19);
      }
    
      if (neighbours_to_find > pool_size)
      {
        cerr << "You asked for " << neighbours_to_find << " neighbours, but pool only contains " << pool_size << ". Shortened\n";
        neighbours_to_find = pool_size;
      }

      if (verbose)
        cerr << "A maximum of " << neighbours_to_find << " neighbours of each molecule will be found\n";
    }

    neighbour_count.resize(neighbours_to_find + 1);
  }

  if (cl.option_present('h'))
  {
    do_not_compare_molecules_with_themselves = 1;

    if (verbose)
      cerr << "Will discard neighbours with zero distance and the same id as the target\n";
  }

// If verbose and a threshold specified, they still need the neighbour characteristics

  if (0 == neighbour_count.elements_allocated())
    neighbour_count.resize(pool_size + 1);

  if (allow_arbitrary_distances)
    ;
  else if (create_histogram)
    histogram_nearnest_neighbour_distances.initialise(0.0, 1.0, 0.01);

  if (cl.option_present('j'))
  {
    int j;
    if (! cl.value('j', j) || j < 2)
    {
      cerr << "The output precision option (-j) must be a whole +ve number\n";
      usage(3);
    }

    set_default_iwstring_float_concatenation_precision(j);

    if (verbose)
      cerr << "Default float concatenation precision " << j << endl;
  }

  int rc = 0;

  IWString_and_File_Descriptor output(1);

  if (cl.option_present('p'))
  {
    if (! three_column_output_all_pairs(output))
    {
      cerr << "Cannot create three column output form\n";
      rc = 4;
    }
  }
  else if (! nearneighbours(output))
  {
    cerr << "Fatal error during neighbour determination\n";
    rc = 3;
  }

  if (verbose)
  {
    cerr << "Read " << fingerprints_read << " fingerprints\n";
    cerr << "Neighbour distances for " << distance_stats.n() << " neighbours between " << distance_stats.minval() << " and " << distance_stats.maxval() << endl;
    if (distance_stats.n() > 1)
      cerr << "Average " << distance_stats.average() << " variance " << distance_stats.variance();
    cerr << endl;

    cerr << "Nearest neighbour distances between " << nearest_neighbour_distance_stats.minval() << " and " << nearest_neighbour_distance_stats.maxval();
    if (nearest_neighbour_distance_stats.n() > 1)
      cerr << " ave " << nearest_neighbour_distance_stats.average() << " std dev " << static_cast<float>(sqrt(nearest_neighbour_distance_stats.variance()));
    cerr << endl;

    if (molecules_rescanned)
      cerr << molecules_rescanned << " molecules rescanned to get at least " << rescan_if_not_enough_neighbours << " neighbours\n";

    for (int i = 0; i < neighbour_count.number_elements(); i++)
    {
      if (neighbour_count[i])
        cerr << neighbour_count[i] << " molecules had " << i << " neighbours\n";
    }
  }

  if (allow_arbitrary_distances)   // don't have a histogram
    ;
  else if (stream_for_nearest_neighbour_histogram.rdbuf()->is_open())
    do_write_histogram(stream_for_nearest_neighbour_histogram);

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = nearneighbours(argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status(stderr);
#endif

  return rc;
}
