#include <stdlib.h>
#include <iostream>
#include <memory>

/*
  Activity imputation by nearest neighbours
*/

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/tversky.h"

using std::cerr;
using std::endl;
using std::flush;

static Tversky tversky;

static IW_GFP_DA * pool = nullptr;

static int pool_size = 0;

static int verbose = 0;

static int neighbours_to_find = 1;

static int flush_output = 0;

static similarity_type_t abandon_distance_threshold = -1.0;

static Accumulator<float> distance_stats;

static int * neighbour_count = nullptr;

/*
  The identifier tag used in each TDT
*/

static IWString identifier_tag ("PCN<");

/*
  The activity tag
*/

static IWString activity_tag;

/*
  During output we need to specify which items from the pool and
  from the target file to echo
*/

static resizable_array_p<IWString> pool_items_to_echo, target_items_to_echo;

static const_IWSubstring distance_tag ("DIST<");

/*
  We ignore distances longer than DISTANCE_THRESHOLD
*/

static similarity_type_t upper_distance_threshold = 1.0;

static similarity_type_t lower_distance_threshold = -1.0;

/*
  Since we will want to go back and get members of the pool, we must make the 
  file from which they come available throughout programme execution.

  That's why we need our own copy of build pool
*/

static iwstring_data_source pool_data_source;

/*
  Our pool is an array of FP objects
*/

typedef float activity_type_t;

class IW_GFP_DA : public IW_GFP_D
{
  private:
    activity_type_t _activity;

  public:
    IW_GFP_DA ();

    int construct_from_tdt (IW_TDT &, int &);

    activity_type_t activity () const { return _activity;}
    activity_type_t & activity () { return _activity;}
};

IW_GFP_DA::IW_GFP_DA ()
{
  _activity = -999.9;

  return;
}

/*
  We keep track of stats on the activities
*/

static Accumulator<activity_type_t> activity_stats;

int
IW_GFP_DA::construct_from_tdt (IW_TDT & tdt, int & fatal)
{
  if (! IW_GFP_D::construct_from_tdt (tdt, fatal))
    return 0;

  if (! tdt.dataitem_value (activity_tag, _activity))
  {
    cerr << "Cannot extract '" << activity_tag << "' from tdt\n";
    fatal = 1;
    return 0;
  }

  activity_stats.extra (_activity);

  return 1;
}

static int
build_pool (iwstring_data_source & input)
{
  std::streampos offset = input.tellg ();

  int items_in_pool = 0;

  IW_TDT tdt;
  while (tdt.next (input))
  {
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

    items_in_pool++;

    if (items_in_pool >= pool_size)
    {
      cerr << "Pool is full, max " << pool_size << endl;
      return 1;
    }

    offset = input.tellg ();
  }

  pool_size = items_in_pool;

  return 1;
}

static int
build_pool (const const_IWSubstring & fname)
{
  IWString tmp (fname);

  if (! pool_data_source.open (tmp))    // method is non-const on its argument!
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
    pool_size = pool_data_source.grep (*pcn);

    if (0 == pool_size)
    {
      cerr << "No occurrences of " << pcn->pattern() << "' in input\n";
      return 0;
    }

    pool = new IW_GFP_DA[pool_size];
    if (NULL == pool)
    {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    cerr << "Pool automatically sized to " << pool_size << endl;
  }

  return build_pool (pool_data_source);
}

static int
echo_tdt_items (IW_TDT & tdt,
                resizable_array_p<IWString> & items_to_echo,
                std::ostream & output)
{
  int nt = tdt.number_elements ();

  int ne = items_to_echo.number_elements ();

  if (0 == ne)      // Echo the whole TDT
  {
    nt--;      // we don't want to echo the trailing '|'
    for (int i = 0; i < nt; i++)
    {
      const IWString & t = *(tdt[i]);
      output << t << endl;
    }
  }
  else
  {
    for (int i = 0; i < ne; i++)
    {
      const IWString & t = *(items_to_echo[i]);
      if (! tdt.echo_dataitem (t, 0, output))
        cerr << "Yipes, cannot echo item '" << t << "'\n";
    }
  }

  return output.good ();
}

/*
  Impute activity based on near neighbours and associated distance
*/

static int
compute_predicted_activity (IW_GFP_DA ** neighbours,
                            int number_neighbours,
                            activity_type_t & predicted_activity)
{
  predicted_activity = 0.0;

  for (int i = 0; i < number_neighbours; i++)
  {
    const IW_GFP_DA * n = neighbours[i];

    similarity_type_t d = (1.0 - n->distance ());
    activity_type_t   a = n->activity ();

    predicted_activity += a * d;
  }

  predicted_activity = predicted_activity / float (number_neighbours);

  return 1;
}

static int
write_neighbour_list (IW_TDT & target,
                      IW_GFP_D ** neighbours,
                      int number_neighbours,
                      std::ostream & output)
{
  activity_type_t predicted_activity;
  if (! compute_predicted_activity (neighbours, number_neighbours, predicted_activity))
  {
    return 1;
  }

  echo_tdt_items (target, target_items_to_echo, output);

  if (NULL != neighbour_count)
    neighbour_count[number_neighbours]++;

  for (int i = 0; i < number_neighbours; i++)
  {
    IW_GFP_DA * n = neighbours[i];
    assert (NULL != neighbours[i]);

    off_t o;
    (void) n->offset (o);

    if (! pool_data_source.seekg (o))
    {
      cerr << "Yipes, neighbour " << i << " offset " << o << " cannot seek\n";
      return 0;
    }

    IW_TDT tdt;
    (void) tdt.next (pool_data_source);    // should not fail!
    
    echo_tdt_items (tdt, pool_items_to_echo, output);

    output << distance_tag << n->distance () << ">\n";
    distance_stats.extra (n->distance ());
  }

  output << "|\n";

  if (flush_output)
    output << flush;

  return output.good ();
}

//#define DEBUG_INSERTION

/*
  Insert NEIGHBOUR into the neighbour list of maximum size NEIGHBOURS_TO_FIND.
  NEIGHBOUR becomes the new item WHERE
*/

static void
insert_in_neighbour_list (IW_GFP_D ** neighbours, int neighbours_to_find, 
                          int & neighbours_found,
                          IW_GFP_DA & neighbour, int where)
{
  assert (neighbours_found <= neighbours_to_find);
  assert (where < neighbours_to_find);

#ifdef DEBUG_INSERTION
  cerr << "Inserting " << neighbour.distance () << ". To find = " << neighbours_to_find << " found = " << neighbours_found << " where = " << where << endl;
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

#define CHECK_INSERTION
#ifdef CHECK_INSERTION
  int failure = 0;

  for (int i = 1; i < neighbours_found; i++)
  {
    if (neighbours[i - 1]->distance () > neighbours[i]->distance ())
    {
      cerr << "Sort/insertion failed, out of order\n";
      failure = 1;
    }
  }

#ifdef DEBUG_INSERTION
  cerr << "After insertion, Neighbours found = " << neighbours_found << endl;
  for (int i = 0; i < neighbours_found; i++)
  {
    cerr << "i = " << i << " distance " << neighbours[i]->distance () << endl;
  }
#endif

  if (failure)
  {
    for (int i = 0; i < neighbours_found; i++)
    {
      cerr << "i = " << i << " distance " << neighbours[i]->distance () << endl;
    }

    exit (87);
  }

#endif

  return;
}

static int
nearneighbours (IW_General_Fingerprint & fp,
                IW_GFP_D ** neighbours)
{
  int neighbours_found = 0;
  assert (neighbours_to_find > 0);

  for (int i = 0; i < neighbours_to_find; i++)
    neighbours[i] = nullptr;

  neighbours_found = 0;

  for (int i = 0; i < pool_size; i++)
  {
    if (! can_be_compared (fp, pool[i]))
      continue;

    similarity_type_t t;
    if (abandon_distance_threshold > 0.0)
    {
      if (! fp.IW_General_Fingerprint::tanimoto(pool[i], 1.0f - abandon_distance_threshold, t))
        continue;
      t = 1.0f - t;
    }
    else if (tversky.active ())
    {
      t = (1.0 - fp.IW_General_Fingerprint::tversky (pool[i], tversky));
    }
    else
    {
      t = fp.IW_General_Fingerprint::distance (pool[i]);
    }

//#define DEBUG_NN
#ifdef DEBUG_NN
    cerr << "Distance between '" << fp.id () << " and pool " << i << " '" << pool[i].id () << "' is " << t << endl;
#endif

    if (t >= upper_distance_threshold || t <= lower_distance_threshold)
      continue;

    if (neighbours_found == neighbours_to_find && t >= neighbours[neighbours_found - 1]->distance ())
      continue;

    pool[i].set_distance (t);

//  Do we put this new neighbour at the head of the list

    if (0 == neighbours_found || t <= neighbours[0]->distance ())
    {
      insert_in_neighbour_list (neighbours, neighbours_to_find, neighbours_found, pool[i], 0);
      continue;
    }

//  Does it go at the end of the list

    if (t >= neighbours[neighbours_found - 1]->distance ())
    {
      if (neighbours_found == neighbours_to_find)     // list is full
        continue;

      insert_in_neighbour_list (neighbours, neighbours_to_find, neighbours_found, pool[i], neighbours_found);
      continue;
    }

//  It goes somewhere in the list

#ifdef DEBUG_INSERTION
    cerr << "Inserting " << t << " into list of " << neighbours_found << " neighbours\n";

    for (int j = 0; j < neighbours_found; j++)
    {
      cerr << "j = " << j << " dist " << neighbours[j]->distance ();
      if (neighbours[j]->distance () < t)
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
      similarity_type_t m = neighbours[middle]->distance ();
#ifdef DEBUG_INSERTION
      cerr << "left " << left << " middle " << middle << " right " << right << endl;
      cerr << neighbours[left]->distance () << ',' << neighbours[middle]->distance () << ',' << neighbours[right]->distance () << endl;
#endif

      if (t < m)
        right = middle;
      else if (t > m)
        left = middle;
      else    // they are equal, will insert at middle
        break;

      middle = (left + right) / 2;
    }

    insert_in_neighbour_list (neighbours, neighbours_to_find, neighbours_found, pool[i], middle + 1);
  }

  return neighbours_found;
}

static int
nearneighbours (iwstring_data_source & input, std::ostream & output)
{

  IW_GFP_D ** neighbours = new IW_GFP_D * [neighbours_to_find];

  IW_TDT tdt;
  while (tdt.next (input) && output.good ())
  {
    IW_GFP_D fp;
    int fatal;

    if (! fp.construct_from_tdt (tdt, fatal))
    {
      if (fatal)
        return 0;

      continue;
    }

    int nn = nearneighbours (fp, neighbours);

    (void) write_neighbour_list (tdt, neighbours, nn, output);
  }

  delete neighbours;

  return output.good ();
}

static int
nearneighbours (const char * fname, std::ostream & output)
{
  iwstring_data_source input (fname);
  if (! input.ok ())
  {
    cerr << "Cannot open input '" << fname << "'\n";
    return 0;
  }

  return nearneighbours (input, output);
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

  if (all_found && items_to_echo.number_elements ())
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
  cerr << " -p <file>        specify file against which input is to be compared\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -n <number>      specify how many neighbours to find\n";
  cerr << " -t <dis>         specify lower distance threshold\n";
  cerr << " -T <dis>         specify upper distance threshold\n";
  cerr << " -i <dataitem>    specify identifier dataitem in pool\n";
  cerr << " -I <dataitem>    specify identifier dataitem in input file\n";
  cerr << " -e <dataitem>    specify target object dataitems to be echo'd (default $SMI and PCN)\n";
  cerr << " -e ALL           echo all dataitems in the input in target\n";
  cerr << " -E <dataitem>    specify pool object dataitems to be echo'd (default $SMI and PCN)\n";
  cerr << " -E ALL           echo all dataitems from the pool file\n";
  cerr << " -X <distance>    abandon distance computation if any component > distance\n";
  display_standard_tversky_options (cerr);
  display_standard_gfp_options (cerr);
  cerr << " -A <tag>         activity tag\n";
  cerr << " -v               verbose output\n";

  exit (rc);
}

static int
nearneighbours (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:n:p:I:i:e:E:t:T:fF:P:W:X:a:b:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

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

    pool = new IW_GFP_DA[pool_size];
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

  if (cl.option_present ('i'))
  {
    (void) cl.value ('i', identifier_tag);

    set_identifier_tag (identifier_tag);

    if (verbose)
      cerr << "Identifiers in pool tagged as '" << identifier_tag << "'\n";
  }

  if (cl.option_present ('f'))
  {
    flush_output = 1;
    if (verbose)
      cerr << "Output file will be flushed\n";
  }

  if (! cl.option_present ('p'))
  {
    cerr << "Must specify a pool file via the -p option\n";
    usage (5);
  }

  if (1 != cl.option_count ('p'))
  {
    cerr << "Only one pool file may be specified\n";
    usage (6);
  }

  if (cl.option_present ('F') || cl.option_present ('P') || cl.option_present ('W'))
  {
    if (! initialise_fingerprints (cl, verbose))
    {
      cerr << "Cannot initialise general fingerprint options\n";
      usage (17);
    }
  }

  if (! tversky.parse_command_line (cl, 'a', 'b', verbose))
  {
    cerr << "Cannot initialise Tversky conditions\n";
    usage (26);
  }
  
  if (! cl.option_present ('A'))
  {
    cerr << "Must specify activity tag via the -A option\n";
    usage (83);
  }

  if (cl.option_present ('A'))
  {
    activity_tag = cl.string_value ('A');

    if (verbose)
      cerr << "Activity in the '" << activity_tag << "' dataitem\n";
  }

  if (cl.option_present ('p'))
  {
    const_IWSubstring fname;
    cl.value ('p', fname);

    if (! build_pool (fname))
    {
      cerr << "Cannot build pool from '" << fname << "'\n";
      return 76;
    }

    if (verbose)
    {
      cerr << "Activities in pool between " << activity_stats.minval () << " and " << activity_stats.maxval ();
      if (activity_stats.n () > 1)
        cerr << " ave " << activity_stats.average ();
      cerr << endl;
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

// Now that the pool is built, we can switch identifiers if needed

  if (cl.option_present ('I'))
  {
    const_IWSubstring id;
    cl.value ('I', id);

    set_identifier_tag (id);

    if (verbose)
      cerr << "Identifiers in input tagged as '" << id << "'\n";
  }

  if (! process_dash_e_option (cl, 'e', target_items_to_echo))
  {
    cerr << "Cannot process -e option\n";
    usage (15);
  }

  if (! process_dash_e_option (cl, 'E', pool_items_to_echo))
  {
    cerr << "Cannot process -E option\n";
    usage (15);
  }

  if (cl.option_present ('t'))
  {
    if (! cl.value ('t', lower_distance_threshold) || lower_distance_threshold < 0.0 || lower_distance_threshold > 1.0)
    {
      cerr << "The -t option must be followed by a valid similarity value\n";
      usage (12);
    }

    if (verbose)
      cerr << "Lower distance threshold set to " << lower_distance_threshold << endl;
  }

  if (cl.option_present ('T'))
  {
    if (! cl.value ('T', upper_distance_threshold) || upper_distance_threshold < 0.0 || upper_distance_threshold > 1.0)
    {
      cerr << "The -T option must be followed by a valid similarity value\n";
      usage (12);
    }

    if (verbose)
      cerr << "Upper distance threshold set to " << upper_distance_threshold << endl;
  }

  if (cl.option_present ('n'))
  {
    int n;
    if (! cl.value ('n', n) || n < 1)
    {
      cerr << "the -n option must be followed by a whole positive number\n";
      usage (13);
    }
    
    if (n > pool_size)
    {
      cerr << "You asked for " << n << " neighbours, but pool only contains " << pool_size << ". Shortened\n";
      n = pool_size;
    }

    neighbours_to_find = n;
    if (verbose)
      cerr << "A maximum of " << n << " neighbours of each molecule will be found\n";

    neighbour_count = new int[neighbours_to_find + 1];
    for (int i = 0; i < neighbours_to_find + 1; i++)
      neighbour_count[i] = 0;
  }

  for (int i = 0; i < cl.number_elements (); i++)
  {
    (void) nearneighbours (cl[i], std::cout);
  }

  if (verbose)
  {
    cerr << "Neighbour distances for " << distance_stats.n () << " neighbours between " << distance_stats.minval () << " and " << distance_stats.maxval () << endl;
    if (distance_stats.n () > 1)
      cerr << "Average " << distance_stats.average () << " variance " << distance_stats.variance () << endl;

    for (int i = 0; i <= neighbours_to_find; i++)
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
  int rc = nearneighbours (argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status (stderr);
#endif

  return rc;
}
