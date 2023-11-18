/*
  Special purpose near neighbour programme
  This variant can go parallel after building the pool. Each
  process works on a subset of the input file.
  Since the iwstring_data_source object does not support parallel
  access, we write an intermediate file only with offsets into the
  pool file
*/

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fstream>

#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/tversky.h"

using std::cerr;
using std::endl;

/*
  Our pool is an array of FP objects
*/

static IW_General_Fingerprint * pool = nullptr;

static int pool_size = 0;

static int verbose = 0;

static int neighbours_to_find = 1;

/*
  Or we can specify the number of neighbours with each item
*/

static IWString neighbours_to_find_tag;

/*
  If we have a filter active, we may find no neighbours for a molecule. In those cases,
  we can re-scan the pool
*/

static int rescan_if_no_neighbours = 0;

static similarity_type_t abandon_distance_threshold = -1.0;

static Tversky tversky;

/*
  The identifier tag used in each TDT
*/

static IWString identifier_tag("PCN<");

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

static IWString distance_tag("DIST<");

/*
  We ignore distances longer than DISTANCE_THRESHOLD
*/

static similarity_type_t upper_distance_threshold = 1.0;

static similarity_type_t lower_distance_threshold = -1.0;

/*
  When we have thresholds, we may choose to not write molecules with no neighbours
*/

static int write_molecules_with_no_neighbours = 1;

static int molecules_with_no_neighbours = 0;

static int
build_pool(iwstring_data_source & input)
{
  off_t offset = input.tellg();

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

      offset = input.tellg();
      continue;
    }

    pool[items_in_pool].set_offset(offset);

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
build_pool(const const_IWSubstring & fname)
{
  IWString tmp(fname);

  iwstring_data_source input;

  if (! input.open(tmp))    // method is non-const on its argument!
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
    pool_size = input.grep(*pcn);

    if (0 == pool_size)
    {
      cerr << "No occurrences of " << pcn->pattern() << "' in input\n";
      return 0;
    }

    pool = new IW_General_Fingerprint[pool_size];
    if (NULL == pool)
    {
      cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
      return 62;
    }

    cerr << "Pool automatically sized to " << pool_size << endl;
  }

  return build_pool(input);
}

static int
echo_tdt_items (IW_TDT & tdt,
                resizable_array_p<IWString> & items_to_echo,
                std::ostream & output)
{
  int nt = tdt.number_elements();

  int ne = items_to_echo.number_elements();

  if (0 == ne)      // Echo the whole TDT
  {
    tdt.write_all_except_vbar(output);
    nt--;
  }
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

class ID_and_Distance
{
  private:
    IW_General_Fingerprint * _fp;
    similarity_type_t _dist;
  public:
    similarity_type_t distance() const { return _dist;}

    void set_distance(similarity_type_t d) { _dist = d;}

    void set_pool_member(IW_General_Fingerprint * f) { _fp = f;}
    IW_General_Fingerprint * pool_member() const { return _fp;}
};

static int
write_average_neighbour_distance(const ID_and_Distance * neighbours,
                      int number_neighbours,
                      std::ostream & output)
{
  similarity_type_t total = 0.0;

  for (int i = 0; i < number_neighbours; i++)
  {
    const ID_and_Distance & n = neighbours[i];

    total += n.distance();
  }

  float tmp = float(total) / float(number_neighbours);

  output << tag_for_average_distance << tmp << ">\n";

  return output.good();
}

static int
write_neighbour_list(IW_TDT & target,
                      ID_and_Distance * neighbours,
                      int number_neighbours,
                      iwstring_data_source & pool_data_source,
                      std::ostream & output)
{
  echo_tdt_items(target, target_items_to_echo, output);

  if (tag_for_average_distance.length())
    write_average_neighbour_distance(neighbours, number_neighbours, output);

  for (int i = 0; i < number_neighbours; i++)
  {
    ID_and_Distance & n = neighbours[i];

    IW_General_Fingerprint * p = n.pool_member();

    off_t o;
    (void) p->offset(o);

    if (! pool_data_source.seekg(o))
    {
      cerr << "Yipes, neighbour " << i << " offset " << o << " cannot seek\n";
      return 0;
    }

    IW_TDT tdt;
    (void) tdt.next(pool_data_source);    // should not fail!
    
    echo_tdt_items(tdt, pool_items_to_echo, output);

    IWString tmp = distance_tag;
    tmp.append_number(n.distance(), 6);
    tmp << ">\n";

    output << tmp;
  }

  output << "|\n";

  return output.good();
}


static int
id_and_dist_comparitor_function(const void * pid1,
                                 const void * pid2)
{
  const ID_and_Distance * id1 = reinterpret_cast<const ID_and_Distance *>(pid1);
  const ID_and_Distance * id2 = reinterpret_cast<const ID_and_Distance *>(pid2);

//cerr << "Comparing " << pid1 << " with " << pid2 << endl;
//cerr << "Comparing " << id1->pool_member()->id() << '@' << id1->distance() << " with " << id2->pool_member()->id() << '@' << id2->distance() << endl;

  if (id1->distance() < id2->distance())
    return -1;

  if (id1->distance() > id2->distance())
    return 1;

  return 0;
}

/*
  Common code for computing the distance. Note that when an abandon_distance_threshold
  is specified, no Tversky stuff is possible. this is really just for efficiency and
  to avoid makign this function too complicated
*/

static similarity_type_t
compute_the_distance(IW_General_Fingerprint & fp1, IW_General_Fingerprint & fp2,
                      const Tversky & tversky)
{
#ifdef DOES_NOT_WORK
  if (abandon_distance_threshold > 0.0)
  {
    similarity_type_t t;
    if (! fp1.IW_General_Fingerprint::distance(fp2, abandon_distance_threshold, t))
      return 2.0;

    return t;
  }
#endif

  if (! tversky.active())
    return fp1.IW_General_Fingerprint::distance(fp2);

  if (tversky.optimistic_mode())
    return fp1.optimistic_distance(fp2, tversky);

  return 1.0 - fp1.IW_General_Fingerprint::tversky(fp2, tversky);
}

static int
nearneighbours(IW_General_Fingerprint & fp,
                ID_and_Distance * id_and_dist,
                int neighbours_to_find_this_fingerprint)
{
  assert (neighbours_to_find_this_fingerprint > 0);

  for (int i = 0; i < pool_size; i++)
  {
    similarity_type_t t = compute_the_distance(fp, pool[i], tversky);

//#define DEBUG_NN
#ifdef DEBUG_NN
    cerr << "Distance between '" << fp.id() << " and pool " << i << " '" << pool[i].id() << "' is " << t << endl;
#endif

    id_and_dist[i].set_pool_member(&(pool[i]));
    id_and_dist[i].set_distance(t);
  }

  qsort(id_and_dist, pool_size, sizeof(ID_and_Distance), id_and_dist_comparitor_function);

  int neighbours_found = 0;

  for (int i = 0; i < pool_size; i++)
  {
    const ID_and_Distance & idd = id_and_dist[i];

    similarity_type_t d = idd.distance();

    if (d > 1.0)
      break;

    if (d >= upper_distance_threshold || d <= lower_distance_threshold)
      continue;

    neighbours_found++;
  }

  if (neighbours_found > neighbours_to_find_this_fingerprint)
    neighbours_found = neighbours_to_find_this_fingerprint;

  return neighbours_found;
}

static int
nearneighbours(iwstring_data_source & input,
                int molecules_to_process,
                iwstring_data_source & pool_data_source,
                std::ostream & output)
{
  ID_and_Distance * neighbours = new ID_and_Distance[pool_size];

  int fingerprints_read = 0;

  IW_TDT tdt;
  while (tdt.next(input) && output.good())
  {
    IW_General_Fingerprint fp;
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

    int nn = nearneighbours(fp, neighbours, neighbours_to_find_this_fingerprint);

    if (0 == nn && ! write_molecules_with_no_neighbours)
      molecules_with_no_neighbours++;
    else
      (void) write_neighbour_list(tdt, neighbours, nn, pool_data_source, output);

    if (fingerprints_read >= molecules_to_process)
      break;
  }

  delete [] neighbours;

  return output.good();
}

static int
nearneighbours(const char * fname,
                int which_one_am_i,
                int nskip,
                int nprocess,
                iwstring_data_source & pool_data_source,
                std::ostream & output)
{
  iwstring_data_source input(fname);

  if (! input.ok())     // should never happen
  {
    cerr << "Cannot open input '" << fname << "'\n";
    return 0;
  }

  if (nskip > 0)
  {
    std::unique_ptr<re2::RE2> vbar = std::make_unique<re2::RE2>("^\\|");

    if (! input.skip_records(*vbar, nskip))
    {
      cerr << "Yipes, thread " << which_one_am_i << " cannot skip " << nskip << " items\n";
      return 0;
    }
  }

  return nearneighbours(input, nprocess, pool_data_source, output);
}


// Each thread needs to know where to start, how many to do, etc

struct Thread_Data
{
  int _which_one;
  const char * _input_fname;
  int _nskip;
  int _nprocess;
  IWString _pool_file_name;
  IWString _output_file_stem;
  int _rc;
};

static void *
nearneighbour_thread(void * vtd)
{
  Thread_Data * td = reinterpret_cast<Thread_Data *>(vtd);

  if (verbose)
    cerr << "Thread " << td->_which_one << " skip = " << td->_nskip << " do = " << td->_nprocess << endl;

  IWString output_fname = td->_output_file_stem;
  output_fname << td->_which_one;

  std::ofstream output(output_fname.null_terminated_chars(), std::ios::out);

  if (! output.good())
  {
    cerr << "Cannot open output file '" << output_fname << "'\n";
    td->_rc = 0;
    return NULL;
  }

  iwstring_data_source pool_data_source(td->_pool_file_name);

  if (! pool_data_source.good())    // should never happen
  {
    cerr << "Yipes, thread " << td->_which_one << " could not open the pool file '" << td->_pool_file_name << "'\n";
    td->_rc = 0;
    return NULL;
  }

  td->_rc = nearneighbours(td->_input_fname, td->_which_one, td->_nskip, td->_nprocess, pool_data_source, output);

  return NULL;
}

static int
nearneighbours(const char * fname,
                int nfork,
                const IWString & pool_file_name,
                const IWString & output_file_stem)
{
  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << "Cannot open input '" << fname << "'\n";
    return 0;
  }

  IWString tmp;
  tmp << '^' << identifier_tag;

  std::unique_ptr<re2::RE2> pcn;
  iwre2::RE2Reset(pcn, tmp);
  int molecules_in_input = input.grep(*pcn);

  if (0 == molecules_in_input)
  {
    cerr << "Input contains no molecules\n";
    return 0;
  }

  if (molecules_in_input < nfork)
  {
    cerr << "Requested " << nfork << " processed, but only " << molecules_in_input << " molecules in the input\n";
    return 0;
  }

  int molecules_per_thread = molecules_in_input / nfork;

  assert (molecules_per_thread >= 1);

#ifdef RUN_SERIAL
  nfork = 1;
  molecules_per_thread = molecules_in_input;
#endif

  Thread_Data * thread_data = new Thread_Data[nfork];

  pthread_t * neighbour_threads = new pthread_t[nfork];

  pthread_attr_t attr;

  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  for (int i = 0; i < nfork; i++)
  {
    Thread_Data & td = thread_data[i];

    td._which_one = i;
    td._input_fname = fname;
    td._pool_file_name = pool_file_name;
    td._output_file_stem = output_file_stem;

    td._nskip = i * molecules_per_thread;

    if (nfork - 1 == i)    // the last thread to launch
      td._nprocess = molecules_in_input - i * molecules_per_thread;
    else
      td._nprocess = molecules_per_thread;

    if (verbose)
      cerr << "Child process " << i << ", skip = " << td._nskip << " process = " << td._nprocess << endl;

    int err = pthread_create(&(neighbour_threads[i]), &attr, nearneighbour_thread, &td);
    if (0 != err)
    {
      cerr << "Error " << err << " creating thread " << i;
      perror(" thread creation error");
      return 0;
    }
  }

  for (int i = 0; i < nfork; i++)
  {
    pthread_join(neighbour_threads[i], NULL);
  }

  delete [] neighbour_threads;

  delete [] thread_data;

  return 1;
}

static int
process_dash_e_option(Command_Line & cl,
                       char e, 
                       resizable_array_p<IWString> & items_to_echo)
{
  if (! cl.option_present(e))
  {
    items_to_echo.resize(2);
    IWString * t = new IWString("$SMI<");
    items_to_echo.add(t);
    t = new IWString("PCN<");
    items_to_echo.add(t);

    return 1;
  }

  int all_found = 0;

  const_IWSubstring evalue;
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
  cerr << " -p <file>        specify file against which input is to be compared\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -n <number>      specify how many neighbours to find\n";
  cerr << " -N <tag>         the number of neighbours to find is in the TAG dataitem\n";
  cerr << " -t <dis>         specify lower distance threshold\n";
  cerr << " -T <dis>         specify upper distance threshold\n";
  cerr << " -z               don't write molecules with no neighbours\n";
  cerr << " -i <tag>         specify identifier tag in pool\n";
  cerr << " -I <tag>         specify identifier tag in input file\n";
  cerr << " -e <dataitem>    specify target object dataitems to be echo'd (default $SMI and PCN)\n";
  cerr << " -e ALL           echo all dataitems in the input in target\n";
  cerr << " -E <dataitem>    specify pool object dataitems to be echo'd (default $SMI and PCN)\n";
  cerr << " -E ALL           echo all dataitems from the pool file\n";
  cerr << " -F ...           standard gfp options\n";
  cerr << " -V ...           standard Tversky options\n";
  cerr << " -X <distance>    abandon distance computation if any component > distance\n";
  cerr << " -r               re-scan the pool if a molecule doesn't have any neighbours\n";
  cerr << " -f <number>      number of processes to fork\n";
  cerr << " -S <stem>        output file name stem\n";
  cerr << " -v               verbose output\n";

  exit(rc);
}

static int
nearneighbours(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vs:n:p:i:I:e:E:t:T:F:P:W:X:V:zN:rf:S:A:");

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

  if (! cl.option_present('f'))
  {
    cerr << "Must specify number of processes to fork via the -f option\n";
    usage(3);
  }

  int nprocesses;

  if (! cl.value('f', nprocesses) || nprocesses < 2)
  {
    cerr << "The number of processes to fork (-f option) must be a whole number > 1\n";
    usage(3);
  }

  if (cl.option_present('r'))
  {
    rescan_if_no_neighbours = 1;

    if (verbose)
      cerr << "If molecule ends up with no neighbours due to cutoffs, the pool will be rescanned\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, forked version only works on a single input file\n";
    return 9;
  }

  if (cl.option_present('F') || cl.option_present('P') || cl.option_present('W'))
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

  IWString stem;

  if (! cl.option_present('S'))
  {
    cerr << "Must specify stem for output files (-S option)\n";
    usage(6);
  }

  stem = cl.option_value('S');

  if (verbose)
    cerr << "Output files created with stem '" << stem << "'\n";

  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 1)
    {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }

    pool = new IW_General_Fingerprint[pool_size];
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

  IWString pool_file_name;

  if (cl.option_present('p'))
  {
    

    cl.value('p', pool_file_name);

    if (! build_pool(pool_file_name))
    {
      cerr << "Cannot build pool from '" << pool_file_name << "'\n";
      return 76;
    }

    if (0 == pool_size)
    {
      cerr << "Yipes, pool is empty\n";
      return 12;
    }

  }

  if (cl.option_present('I'))
  {
    (void) cl.value('I', identifier_tag);

    set_identifier_tag(identifier_tag);

    if (verbose)
      cerr << "Identifiers in input file tagged as '" << identifier_tag << "'\n";
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

  if (! process_dash_e_option(cl, 'e', target_items_to_echo))
  {
    cerr << "Cannot process -e option\n";
    usage(15);
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
    int n;
    if (! cl.value('n', n) || n < 1)
    {
      cerr << "the -n option must be followed by a whole positive number\n";
      usage(13);
    }
    
    if (n > pool_size)
    {
      cerr << "You asked for " << n << " neighbours, but pool only contains " << pool_size << ". Shortened\n";
      n = pool_size;
    }

    neighbours_to_find = n;
    if (verbose)
      cerr << "A maximum of " << n << " neighbours of each molecule will be found\n";
  }

  if (cl.option_present('N'))
  {
    neighbours_to_find_tag = cl.string_value('N');

    if (verbose)
      cerr << "The number of neighbours to find will be from the '" << neighbours_to_find_tag << "' dataitem\n";
  }

// If verbose and a threshold specified, they still need the neighbour characteristics

  int rc = nearneighbours(cl[0], nprocesses, pool_file_name, stem);

  return ! rc;
}

int
main(int argc, char ** argv)
{
  int rc = nearneighbours(argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status(stderr);
#endif

  return rc;
}
