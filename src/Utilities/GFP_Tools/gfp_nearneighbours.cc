/*
  Special purpose near neighbour programme
*/

#include <stdlib.h>

#include "cmdline.h"
#include "accumulator.h"
#include "iwstring_data_source.h"
#include "iw_tdt.h"
#include "iwhistogram.h"

#include "gfp.h"
#include "sparse_collection.h"

#include "tversky.h"

#include "iw_tdt_filter.h"
/*
  gfp_nearneighbours.cc,v 1.3 2002/06/21 19:05:32 
*/

/*
  When doing near neighbour determinations within a single set of molecules, we gain
  great efficiencies by writing the index of the item rather than its name. This helps
  programmes that read the nn file
*/

class IW_GFP_D_ID : public IW_GFP_D
{
  private:
    int _ndx;

  public:
    IW_GFP_D_ID ();

    void set_index (int n) { _ndx = n;}
    int  index_in_pool () const { return _ndx;}
};

IW_GFP_D_ID::IW_GFP_D_ID ()
{
  _ndx = -1;
}

/*
  Our pool is an array of FP objects
*/

static IW_GFP_D_ID * pool = NULL;

static int write_neighbours_as_index_numbers = 0;

static int pool_size = 0;

static int verbose = 0;

static int neighbours_to_find = 1;

/*
  Or we can specify the number of neighbours with each item
*/

static IWString neighbours_to_find_tag;

/*
  If we have a filter active, we may not get enough neighbours. We can rescan the
  pool with all the exclusions turned off
*/

static int rescan_if_not_enough_neighbours = 0;

static int molecules_rescanned = 0;

static int flush_output = 0;

static similarity_type_t abandon_distance_threshold = -1.0;

static Accumulator<float> distance_stats;

/*
  Sometimes it is useful to know the statistics of the nearest neighbour of each molecule
*/

static Accumulator<similarity_type_t> nearest_neighbour_distance_stats;

static IWHistogram histogram_nearnest_neighbour_distances;

static extending_resizable_array<int> neighbour_count;

static int fingerprints_read = 0;

static Tversky tversky;

static Set_of_Sparse_Fingerprint_Collection_Profile sfcp;

/*
  The identifier tag used in each TDT
*/

static IWString identifier_tag ("PCN<");

/*
  If we are writing neighbours as index numbers
*/

static IWString neighbour_tag ("NBR<");

/*
  May 99. For each input TDT, I need to know the average distance
  of the neighbours within the pool
*/

static IWString tag_for_average_distance;

/*
  During output we need to specify which items from the pool and
  from the target file to echo
*/

static resizable_array_p<IWString> pool_items_to_echo, target_items_to_echo;

static const_IWSubstring distance_tag ("DIST<");

/*
  We ignore distances longer than DISTANCE_THRESHOLD
*/

static similarity_type_t upper_distance_threshold = 2.0;

static similarity_type_t lower_distance_threshold = -1.0;

/*
  When we have thresholds, we may choose to not write molecules with no neighbours
*/

static int write_molecules_with_no_neighbours = 1;

static int molecules_with_no_neighbours = 0;

/*
  We often do a run where we find the neighbours of a set of molecules within itself.
  We may or may not want to find the molecule itself as its own nearest neighbour
*/

static int do_not_compare_molecules_with_themselves = 0;

/*
  Since we will want to go back and get members of the pool, we must make the 
  file from which they come available throughout programme execution.

  That's why we need our own copy of build pool
*/

static iwstring_data_source pool_data_source;

static int
build_pool (iwstring_data_source & input,
            IW_TDT_Filter & filter)
{
  off_t offset = input.tellg ();

  int items_in_pool = 0;

  int tdts_read = 0;

  IW_TDT tdt;
  while (tdt.next (input))
  {
    tdts_read++;

    if (filter.active ())
    {
      if (! filter.matches (tdt))
      {
        cerr << "failed filter\n";
        offset = input.tellg ();
        continue;
      }
    }

    int fatal;
    if (! pool[items_in_pool].construct_from_tdt (tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      offset = input.tellg ();
      continue;
    }

    pool[items_in_pool].set_offset (offset);
    pool[items_in_pool].set_index (items_in_pool);

    items_in_pool++;

    if (items_in_pool == pool_size)
    {
      if (verbose)
        cerr << "Pool is full, max " << pool_size << endl;
      break;
    }

    offset = input.tellg ();
  }

  pool_size = items_in_pool;

  if (verbose)
    cerr << "Read " << tdts_read << " TDT's, pool contains " << pool_size << " fingerprints\n";

  return 1;
}

static int
build_pool (const const_IWSubstring & fname, 
            IW_TDT_Filter & filter)
{
  IWString tmp (fname);

  if (! pool_data_source.open (tmp))    // method is non-const on its argument!
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size)
  {
    pool_size = pool_data_source.count_records_starting_with(identifier_tag);

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

  return build_pool(pool_data_source, filter);
}

static int
write_average_neighbour_distance (IW_GFP_D_ID ** neighbours,
                      int number_neighbours,
                      std::ostream & output)
{
  similarity_type_t total = static_cast<similarity_type_t>(0.0);

  for (int i = 0; i < number_neighbours; i++)
  {
    IW_GFP_D_ID * n = neighbours[i];
    assert (NULL != neighbours[i]);

    total += n->distance();
  }

  float tmp = float(total) / float(number_neighbours);

  output << tag_for_average_distance << tmp << ">\n";

  return output.good();
}

static int
echo_tdt_items (IW_TDT & tdt,
                resizable_array_p<IWString> & items_to_echo,
                std::ostream & output)
{
  int ne = items_to_echo.number_elements();

  if (0 == ne)      // Echo the whole TDT
    tdt.write_all_except_vbar(output);
  else
  {
    for (int i = 0; i < ne; i++)
    {
      const IWString & t = *(items_to_echo[i]);
      if (! tdt.echo_dataitem(t, 0, output))
        cerr << "Yipes, cannot echo item '" << t << "'\n";
    }
  }

  return output.good();
}

static int
write_neighbour_list (IW_TDT & target,
                      IW_GFP_D_ID ** neighbours,
                      int number_neighbours,
                      std::ostream & output)
{
  echo_tdt_items(target, target_items_to_echo, output);

  if (tag_for_average_distance.length())
    write_average_neighbour_distance(neighbours, number_neighbours, output);

  for (int i = 0; i < number_neighbours; i++)
  {
    IW_GFP_D_ID * n = neighbours[i];
    assert (NULL != neighbours[i]);

    if (write_neighbours_as_index_numbers)
    {
      output << neighbour_tag << n->index_in_pool() << ">\n";
    }
    else
    {
      off_t o;
      (void) n->offset(o);

      if (! pool_data_source.seekg(o))
      {
        cerr << "Yipes, neighbour " << i << " offset " << o << " cannot seek\n";
        return 0;
      }

      IW_TDT tdt;
      (void) tdt.next(pool_data_source);    // should not fail!
    
      echo_tdt_items(tdt, pool_items_to_echo, output);
    }

    similarity_type_t d = n->distance();

    IWString tmp(distance_tag);
    tmp.append_number(d, 3);
    tmp += ">\n";
    output << tmp;

//  output << distance_tag << d << ">\n";

    distance_stats.extra(d);

    if (0 == i && verbose)
    {
      nearest_neighbour_distance_stats.extra(d);
      histogram_nearnest_neighbour_distances.extra(d);
    }
  }

  output << "|\n";

  if (flush_output)
    output.flush();

  return output.good();
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
    neighbours[0] = & neighbour;
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
nearneighbours (IW_GFP_D_ID & fp,
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

    assert (t >= static_cast<similarity_type_t>(0.0) && t <= static_cast<similarity_type_t>(1.0));

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
nearneighbours (iwstring_data_source & input, std::ostream & output)
{
  IW_GFP_D_ID ** neighbours;

  if (neighbours_to_find_tag.length())    // maybe every pool member is a neighbour
    neighbours = new IW_GFP_D_ID * [pool_size];
  else
    neighbours = new IW_GFP_D_ID * [neighbours_to_find];

  IW_TDT tdt;
  while (tdt.next(input) && output.good())
  {
    IW_GFP_D_ID fp;
    int fatal;

    if (! fp.construct_from_tdt(tdt, fatal))
    {
      if (fatal)
        return 0;

      continue;
    }

    int neighbours_to_find_this_fingerprint = neighbours_to_find;   // the default value

    if (neighbours_to_find_tag.length())
    {
      if (tdt.index_of_dataitem(neighbours_to_find_tag) < 0)    // no value specified, use the global default
        cerr << "No " << neighbours_to_find_tag << " tag present\n";
      else if (! tdt.dataitem_value(neighbours_to_find_tag, neighbours_to_find_this_fingerprint))
      {
        cerr << "Invalid neighbours to find value\n";
        cerr << tdt;
        return 0;
      }
      else if (neighbours_to_find_this_fingerprint > pool_size)
      {
        cerr << "Cannot find " << neighbours_to_find_this_fingerprint << " neighbours, only " << pool_size << " available\n";
        neighbours_to_find_this_fingerprint = pool_size;
      }
      else if (verbose > 1)
        cerr << "Will find " << neighbours_to_find_this_fingerprint << " neighbours for '" << fp.id() << "'\n";
    }

    fingerprints_read++;

    if (sfcp.active())
      fp.convert_to_non_sparse_forms(sfcp);

    int nn = nearneighbours(fp, neighbours, neighbours_to_find_this_fingerprint);

    neighbour_count[nn]++;

    if (0 == nn && ! write_molecules_with_no_neighbours)
      molecules_with_no_neighbours++;
    else
      (void) write_neighbour_list(tdt, neighbours, nn, output);
  }

  delete neighbours;

  return output.good();
}

static int
nearneighbours (const char * fname, std::ostream & output)
{
  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << "Cannot open input '" << fname << "'\n";
    return 0;
  }

  return nearneighbours(input, output);
}

/*
  Make sure write_neighbours_as_index_numbers is set at the appropriate time - after
  the -e option and before the -E option.
*/

static int
process_dash_e_option (Command_Line & cl,
                       char e, 
                       resizable_array_p<IWString> & items_to_echo)
{
  if (! cl.option_present(e))
  {
    if (write_neighbours_as_index_numbers)    // order dependency
      return 1;

    items_to_echo.resize(2);
    IWString * t = new IWString("$SMI<");
    items_to_echo.add(t);
    t = new IWString("PCN<");
    items_to_echo.add(t);

    return 1;
  }

  int all_found = 0;

  const_IWSubstring evalue;

  if (1 == cl.option_count(e))
  {
    evalue = cl.string_value(e);

    if ("ALL" == evalue)
    {
      return 1;
    }

    IWString * tmp = new IWString(e);
    items_to_echo.add(tmp);

    return 1;
  }

  int i = 0;
  while (cl.value(e, evalue, i++))
  {
    if ("ALL" == evalue)
    {
      all_found = 1;
      if (verbose)
        cerr << "Will echo entire tdt on output\n";
    }
    else
    {
      IWString * t = new IWString(evalue);
      items_to_echo.add(t);
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
  cerr << endl;

  cerr << "Finds near neighbours of a set of fingerprints\n";
  cerr << "Usage <options> <input_file>\n";
  cerr << " -p <file>        file against which input is to be compared (haystack)\n";
  cerr << " -s <number>      max pool size\n";
  cerr << " -n <number>      how many neighbours to find\n";
  cerr << " -t <dis>         discard distances shorter than <dis>\n";
  cerr << " -T <dis>         discard distances longer than <dis>\n";
  cerr << " -z               don't write molecules with no neighbours\n";
//cerr << " -i <dataitem>    identifier dataitem in pool\n";
//cerr << " -I <dataitem>    identifier dataitem in input file\n";
//cerr << " -e <dataitem>    target object dataitems to be echo'd (default $SMI and PCN)\n";
//cerr << " -e ALL           echo all dataitems in the input in target\n";
//cerr << " -E <dataitem>    pool object dataitems to be echo'd (default $SMI and PCN)\n";
//cerr << " -E ALL           echo all dataitems from the pool file\n";
  cerr << " -A <TAG>         write average neighbour distance to <TAG>\n";
//cerr << " -O <qualifier>   options for TDT filter (enter -O help for details)\n";
  cerr << " -X <distance>    abandon distance computation if any component > distance\n";
  cerr << " -r <number>      ensure all molecules have at least <number> neighbours\n";
  cerr << " -h               discard nbrs with at dist and same ID as target\n";
  cerr << " -o               cross referencing single file. Write nbrs as index numbers\n";
  cerr << " -F ...           gfp options, enter '-F help' for details\n";
  cerr << " -V ...           Tversky specification, enter '-V help' for details\n";
  cerr << " -v               verbose output\n";

  exit(rc);
}

static int
nearneighbours(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vs:n:p:I:i:e:E:t:T:fF:P:W:Q:X:V:O:A:zN:r:hoK:");

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

  if (cl.option_present('f'))
  {
    flush_output = 1;
    if (verbose)
      cerr << "Output file will be flushed\n";
  }

  if (! cl.option_present('p'))
  {
    cerr << "Must specify a pool file via the -p option\n";
    usage(5);
  }

  if (1 != cl.option_count('p'))
  {
    cerr << "Only one pool file may be specified\n";
    usage(6);
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

  if (cl.option_present('K'))
  {
    if (! parse_sparse_to_dense_fingerprint_specifications(cl, 'K', verbose))
    {
      cerr << "Invalid sparse->fixed specification(s) (-K option)\n";
      return 4;
    }
  }

// We need to be very careful if we are reading from stdin

  const char * p = cl.option_value('p');

  if (! need_to_call_initialise_fingerprints(cl))
    ;
  else if (1 == cl.number_elements() && 0 == strcmp("-", cl[0]))
  {
    if (! initialise_fingerprints(p, verbose))  // must initialise from pool file
    {
      cerr << "Cannot initialise GFP options\n";
      usage (23);
    }
    cerr << "Warning, fingerprints initialised from -p file\n";
  }
  else if (! initialise_fingerprints(cl, verbose))
  {
    cerr << "Cannot initialise GFP options\n";
    usage(23);
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

// We need to be careful with the -i and -I options. Remember
// that the pool is built first

  if (cl.option_present('i'))
  {
    (void) cl.value('i', identifier_tag);

    set_identifier_tag(identifier_tag);

    if (verbose)
      cerr << "Identifiers in pool tagged as '" << identifier_tag << "'\n";
  }

  if (cl.option_present('p'))
  {
    IW_TDT_Filter filter;

    if (cl.option_present('O'))
    {
      const_IWSubstring o = cl.string_value('O');
      if ("help" == o)
      {
        display_tdt_filter_syntax(cerr);
        return 2;
      }

      if (! filter.build_from_string(o))
      {
        cerr << "Cannot parse the -O qualifier '" << o << "'\n";
        usage(21);
      }

      if (verbose)
        cerr << "Will filter pool TDT's according to filter\n";
    }

    const_IWSubstring fname;
    cl.value('p', fname);

    if (! build_pool(fname, filter))
    {
      cerr << "Cannot build pool from '" << fname << "'\n";
      return 76;
    }

    if (0 == pool_size)
    {
      cerr << "Yipes, pool is empty\n";
      return 12;
    }

    if (verbose && filter.active())
      filter.report(cerr);
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

// Now that the pool is built, we can switch identifiers if needed

  if (cl.option_present('I'))
  {
    const_IWSubstring id;
    cl.value('I', id);

    set_identifier_tag(id);

    if (verbose)
      cerr << "Identifiers in input tagged as '" << id << "'\n";
  }

  if (! process_dash_e_option(cl, 'e', target_items_to_echo))
  {
    cerr << "Cannot process -e option\n";
    usage(15);
  }

// process the -o option here, after the -e but before the -E option

  if (cl.option_present('o'))
  {
    if (cl.option_present('E'))
    {
      cerr << "The -E option is likely to mess up the -o option, beware!\n";
    }

    write_neighbours_as_index_numbers = 1;

    if (verbose)
      cerr << "Neighbours written as index numbers\n";
  }

  if (! process_dash_e_option(cl, 'E', pool_items_to_echo))
  {
    cerr << "Cannot process -E option\n";
    usage(15);
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

  if (cl.option_present('N'))
  {
    neighbours_to_find_tag = cl.string_value('N');

    if (verbose)
      cerr << "The number of neighbours to find will be from the '" << neighbours_to_find_tag << "' dataitem\n";

    neighbour_count.resize(pool_size);    // worst case scenario
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

  if (verbose)
    histogram_nearnest_neighbour_distances.initialise(0.0, 1.0, 0.01);

  for (int i = 0; i < cl.number_elements(); i++)
  {
    (void) nearneighbours(cl[i], std::cout);
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
      cerr << " ave " << nearest_neighbour_distance_stats.average() << " variance " << nearest_neighbour_distance_stats.variance();
    cerr << endl;

    if (molecules_rescanned)
      cerr << molecules_rescanned << " molecules rescanned to get at least " << rescan_if_not_enough_neighbours << " neighbours\n";

    histogram_nearnest_neighbour_distances.write(cerr);

    for (int i = 0; i < neighbour_count.number_elements(); i++)
    {
      if (neighbour_count[i])
        cerr << neighbour_count[i] << " molecules had " << i << " neighbours\n";
    }
  }

  return 0;
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
