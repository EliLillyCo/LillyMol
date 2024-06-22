/*
  Spread implementation
*/

#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <iostream>
#include <memory>
#include <random>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/numeric_data_from_file.h"

#include "Utilities/GFP_Tools/tversky.h"
#include "Utilities/GFP_Tools/spread.h"

using std::cerr;
using std::endl;

static int verbose = 0;

/*
  In retrospect, I'm really not sure why anyone would use a Tversky measure
  in this programme
*/

static Tversky tversky;

/*
  Idea from Dave Cummins.
  When doing a run with no pre-selected molecules, start with the
  object which is furthest from the first fingerprint.
*/

static int start_with_object_furthest_from_first = 0;

static int start_with_object_furthest_from_everything = 0;

static int choose_first_item_randomly = 0;

static int first_item_is_one_with_highest_scale_factor = 0;

static int already_selected_molecules_present = 0;

/*
  Our pool is an array of FP objects
*/

static Spread_Object ** pool = nullptr;

static int pool_size = 0;

static int poolptr = 0;    // the next item in the pool to be filled

static int number_to_select = 0;

static int report_establish_initial_distances = 0;

/*
  We keep track of the distances of the nearest selected items
*/

static Accumulator<similarity_type_t> nearest_selected_neighbour_distance;

static int retest_no_neighbours = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");

static IWString previously_computed_nn_distance_tag;

static Numeric_Data_From_File<float> previously_computed_distances;

/*
  The normal output of the near neighbour programmes leaves $SMI and PCN tag in the
  output file. Those would confuse things, so they must be transformed to something
  else. We assume that has been done
*/

static IWString nn_smiles_tag("NNSMI<");
static IWString nn_id_tag    ("NNID<");

static float blurr_distances = static_cast<float>(0.0);

static int squeeze_frequency = 100;

/*
  The threads operate in two phases.
  Wait for time_for_child_to_wake1
  Scan for the furthest unselected molecule
  Post  to child_done
  wait for time_for_child_to_wake2
  update unselected items for newly selected item
  post to child_done
*/

static int nthreads = 0;

static pthread_t * child_thread = nullptr;

static sem_t * time_for_child_to_wake1 = nullptr;
static sem_t * time_for_child_to_wake2 = nullptr;

static sem_t child_done;

static pthread_t * ps_child_thread = nullptr;

static int
get_previously_computed_nearest_neighbour (Spread_Object & p,
                                           const IW_STL_Hash_Map_float & previously_computed_distances)
{
  IW_STL_Hash_Map_float::const_iterator f = previously_computed_distances.find(p.id());

  if (f == previously_computed_distances.end())    // OK if no previously computed distance
    return 1;

  p.set_nearest_previously_selected_neighbour("C", "UNK", (*f).second);

  return 1;
}

static int
get_previously_computed_nearest_neighbour (const IW_TDT & tdt,
                                           Spread_Object & p,
                                           const IWString & previously_computed_nn_distance_tag)
{
  similarity_type_t d;
  if (! tdt.dataitem_value(previously_computed_nn_distance_tag, d))
    return 0;

  IWString nnsmiles;
  if (! tdt.dataitem_value(nn_smiles_tag, nnsmiles))
    return 0;

  IWString nnid;
  if (! tdt.dataitem_value(nn_id_tag, nnid))
    return 0;

  p.set_nearest_previously_selected_neighbour(nnsmiles, nnid, d);

  return 1;
}

static int
build_pool (iwstring_data_source & input,
            const IWString & previously_computed_nn_distance_tag)
{
  assert (pool_size > 0);
//cerr << "Pool ptr " << poolptr << ", pool size " << pool_size << endl;
  assert (poolptr >= 0 && poolptr < pool_size);

  int items_with_previously_computed_distances = 0;

  IW_TDT tdt;
  while (tdt.next(input))
  {
    Spread_Object * tmp = new Spread_Object;

    int fatal;
    if (! tmp->construct_from_tdt(tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot parse tdt " << tdt;
        delete tmp;
        return 0;
      }

      continue;
    }

    if (previously_computed_nn_distance_tag.length())
    {
      get_previously_computed_nearest_neighbour(tdt, *tmp, previously_computed_nn_distance_tag);  // should really check for fatal errors
    }
    else if (previously_computed_distances.size())
    {
      if (! get_previously_computed_nearest_neighbour(*tmp, previously_computed_distances))
        return 0;
    }
      
    if (tmp->has_a_nearest_selected_neighbour())
      items_with_previously_computed_distances++;

    pool[poolptr] = tmp;
    poolptr++;

    if (poolptr >= pool_size)
    {
      cerr << "Pool is full, max " << pool_size << endl;
      break;
    }
  }

  poolptr--;

  if (verbose)
  {
    cerr << "Pool now contains " << (poolptr + 1) << " objects\n";
    if (previously_computed_nn_distance_tag.length())
      cerr << items_with_previously_computed_distances << " items had previously computed distances\n";
  }

  return 1;
}

static int
allocate_pool()
{
  assert (pool_size > 0);
  assert (NULL == pool);

  pool = new Spread_Object *[pool_size];
  if (NULL == pool)
  {
    cerr << "Yipes, could not allocate pool of size " << pool_size << endl;
    return 62;
  }

  if (verbose)
    cerr << "system sized to " << pool_size << endl;

  return 1;
}

static int
build_pool (const char * fname,
            const IWString & previously_computed_nn_distance_tag)
{
  iwstring_data_source input(fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size)
  {
    std::unique_ptr<re2::RE2> pcn_rx = std::make_unique<re2::RE2>("^PCN<");

    pool_size = input.grep(*pcn_rx);

    if (0 == pool_size)
    {
      cerr << "Zero occurrences of '" << pcn_rx->pattern() << "' in '" << fname << "'\n";
      return 0;
    }

    if (! allocate_pool())
      return 0;
  }

  return build_pool(input, previously_computed_nn_distance_tag);
}

/*
  After establishing initial distances, some pool members don't have a nearest selected
  neighbour. re-scan the pool
*/

static int
rescan_for_no_neighbours (iwstring_data_source & input)
{
  resizable_array<Spread_Object *> to_scan;
  to_scan.resize(pool_size);

  for (int i = 0; i < pool_size; i++)
  {
    if (! pool[i]->has_a_nearest_selected_neighbour())
      to_scan.add(pool[i]);
  }

  int nts = to_scan.number_elements();

  if (0 == nts)
    return 1;

  if (verbose)
    cerr << "After reading previously selected, " << nts << " items with no neighbour\n";

  if (! input.seekg(0))
  {
    cerr << "rescan_for_no_neighbours: yipes, cannot seek back to beginning of file\n";
    return 0;
  }

  IW_TDT tdt;

  while (tdt.next(input))
  {
    int fatal;
    Spread_Object fp;
    if (! fp.construct_from_tdt(tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot build fingerpint\n" << tdt << endl;
        return 0;
      }

      continue;
    }

    for (int i = 0; i < to_scan.number_elements(); i++)
    {
      Spread_Object * p = to_scan[i];

      if (tversky.active())
        p->object_has_been_selected(fp, tversky);
      else
        p->object_has_been_selected(fp);
    }
  }

  return 1;
}

static int
post_to_array_of_semaphores (sem_t * sem,
                             int nsem)
{
  for (int i = 0; i < nsem; i++)
  {
    if (0 != sem_post(&(sem[i])))
    {
      cerr << "Cannot do semaphore post, semphore " << i << endl;
      perror("sem_post");
      return 0;
    }
  }

  return 1;
}

struct PS_Thread_Communication
{
  int thread_number;
  int istart;
  int istop;
  Spread_Object * previously_selected_fp;
};

static PS_Thread_Communication * pstc = nullptr;

static void
thread_establish_initial_distances (PS_Thread_Communication * pstc)
{
  Spread_Object * fp = pstc->previously_selected_fp;

  for (int i = pstc->istart; i < pstc->istop; i++)
  {
    if (tversky.active())
      pool[i]->object_has_been_selected(*fp, tversky);
    else
      pool[i]->object_has_been_selected(*fp);
  }

  return;
}

static void *
thread_establish_initial_distances (void * v)
{
  PS_Thread_Communication * pstc = reinterpret_cast<PS_Thread_Communication *>(v);

  int mythread = pstc->thread_number;

  while (1)
  {
    sem_wait(&(time_for_child_to_wake1[mythread]));

    if (pstc->thread_number < 0)
    {
      sem_post(&child_done);
      pthread_exit(NULL);
    }

    thread_establish_initial_distances(pstc);
    sem_post(&child_done);
  }
}

static void
wait_for_children ()
{
  for (int i = 0; i < nthreads; i++)
  {
    sem_wait(&child_done);
  }

  return;
}

static int
establish_initial_distances (iwstring_data_source & input)
{
  int ntdt = 0;

  IW_TDT tdt;

  if (! tdt.next(input))
  {
    cerr << "Cannot read first tdt from previously selected file\n";
    return 0;
  }

  Spread_Object fp1;

  int fatal;
  if (! fp1.construct_from_tdt(tdt, fatal))   // all errors fatal here
    return 0;

  for (int i = 0; i < nthreads; i++)
  {
    pstc[i].previously_selected_fp = &fp1;
  }

  post_to_array_of_semaphores(time_for_child_to_wake1, nthreads);

  Spread_Object fp2;
  int ndx = 0;

  while (tdt.next(input))
  {
    Spread_Object * s;
    if (0 == ndx)
    {
      s = &fp2;
      ndx = 1;
    }
    else
    {
      s = &fp1;
      ndx = 0;
    }

    int fatal;
    if (! s->construct_from_tdt(tdt, fatal))
    {
      cerr << "Cannot build fingerpint\n" << tdt << endl;
      return 0;
    }

    ntdt++;

    wait_for_children();

    for (int i = 0; i < nthreads; i++)
    {
      pstc[i].previously_selected_fp = s;
    }
    
    post_to_array_of_semaphores(time_for_child_to_wake1, nthreads);

    if (report_establish_initial_distances && 0 == ntdt % report_establish_initial_distances)
      cerr << "Established initial distances " << ntdt << endl;
  }

  wait_for_children();

  if (0 == ntdt)
    cerr << "establish_initial_distances:: warning, no TDT's read\n";

  if (verbose > 1)
  {
    cerr << "INitial distances established\n";
    for (int i = 0; i < pool_size; i++)
    {
      cerr << "Pool " << i << " is " << pool[i]->id() << " dist " << pool[i]->distance() << endl;
    }
  }

  if (retest_no_neighbours)
    rescan_for_no_neighbours(input);

  return 1;
}

static int
initialise_ps_threads(int nthreads)
{
  assert (nthreads > 0);
  assert (NULL != pstc);

  int items_per_chunk = pool_size / nthreads;

  if (items_per_chunk < 1)
  {
    cerr << "Too many threads " << nthreads << " for " << pool_size << " fingerprints\n";
    return 0;
  }

  ps_child_thread = new pthread_t[nthreads];

  pthread_attr_t attr;

  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  for (int i = 0; i < nthreads; i++)
  {
    PS_Thread_Communication & pstci = pstc[i];

    pstci.thread_number = i;
    pstci.istart = i * items_per_chunk;
    if (i == nthreads - 1)
      pstci.istop = pool_size;
    else
      pstci.istop = pstci.istart + items_per_chunk;

    if (0 != pthread_create(&(ps_child_thread[i]), &attr, thread_establish_initial_distances, &pstci))
    {
      perror("cannot create child process");
      return 0;
    }
  }

  return 1;
}

static int
initialise_array_of_semaphores (sem_t * sem,
                                int nsem)
{
  for (int i = 0; i < nsem; i++)
  {
    if (0 != sem_init(&(sem[i]), 0, 0))
    {
      perror("sem_init error");
      return 0;
    }
  }

  return nsem;
}

static int
establish_initial_distances (const const_IWSubstring & fname)
{
  static int first_call = 1;

  if (first_call)
  {
    for (int i = 0; i < pool_size; i++)
    {
      pool[i]->set_distance(2.0);
    }

    first_call = 0;

    if (! initialise_array_of_semaphores(time_for_child_to_wake1, nthreads))
      return 0;

    pstc = new PS_Thread_Communication[nthreads];

    if (0 != sem_init(&child_done, 0, 0))
    {
      perror("Cannot init child done");
      return 0;
    }

    if (! initialise_ps_threads(nthreads))
      return 0;
  }

  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << "Cannot open already selected file '" << fname << "'\n";
    return 0;
  }

  if (verbose)
    cerr << "Establishing initial distances wrt '" << fname << "'\n";

  return establish_initial_distances(input);
}

static void
do_object_has_been_selected_tversky (int isel)
{
  Spread_Object & fpsel = *(pool[isel]);

  for (int i = 0; i < pool_size; i++)
  {
    if (pool[i]->selected() || i == isel)
      continue;

    if (static_cast<float>(0.0) != blurr_distances)
      pool[i]->object_has_been_selected(fpsel, tversky, blurr_distances);
    else
      pool[i]->object_has_been_selected(fpsel, tversky);
  }

  return;
}

static void
do_object_has_been_selected_no_blurring (int isel)
{
  Spread_Object & fpsel = *(pool[isel]);

  for (int i = 0; i < pool_size; i++)
  {
    if (pool[i]->selected() || i == isel)
      continue;

    pool[i]->object_has_been_selected(fpsel);
  }

  return;
}

static void
do_object_has_been_selected_with_blurring (int isel)
{
  Spread_Object & fpsel = *(pool[isel]);

  for (int i = 0; i < pool_size; i++)
  {
    if (pool[i]->selected() || i == isel)
      continue;

    pool[i]->object_has_been_selected(fpsel, blurr_distances);
  }

  return;
}

static void
do_object_has_been_selected (int isel)
{
  if (tversky.active())
    do_object_has_been_selected_tversky(isel);
  else if (static_cast<float>(0.0) == blurr_distances)
    do_object_has_been_selected_no_blurring(isel);
  else
    do_object_has_been_selected_with_blurring(isel);

  return;
}

static similarity_type_t
compute_the_distance (Spread_Object & fp1, Spread_Object & fp2)
{
  if (tversky.active())
    return static_cast<similarity_type_t>(1.0) - fp1.IW_General_Fingerprint::tversky(fp2, tversky);

  return fp1.IW_General_Fingerprint::distance(fp2);
}

static int
choose_largest_previously_computed_distance()
{
  int rc = -1;
  similarity_type_t dmax = static_cast<similarity_type_t>(0.0);

  for (int i = 0; i < pool_size; i++)
  {
    similarity_type_t d = pool[i]->distance();

    if (static_cast<similarity_type_t>(1.0) == d)
      continue;

    if (d > dmax)
    {
      rc = i;
      dmax = d;
    }
  }

  if (rc < 0)
  {
    cerr << "Warning, no items have previously computed distances!\n";
    rc = 0;
  }

  return rc;
}

/*
  No checks as to whether things have scaling factors or not
*/

static int
item_with_highest_scale_factor ()
{
  float highest_scale = pool[0]->scale();
  int rc = 0;

  for (int i = 1; i < pool_size; i++)
  {
    if (pool[i]->scale() > highest_scale)
    {
      highest_scale = pool[i]->scale();
      rc = i;
    }
  }

  return rc;
}

/*
  Since Tversky is asymmetric, we should probably do both loops 1-pool_size
*/

static int
do_start_with_object_furthest_from_everything (int & istart)
{
  int id_of_further_distance_encountered = -1;
  similarity_type_t furthest_distance_encountered = 0.0;

  for (int i = 0; i < pool_size; i++)
  {
    Spread_Object & pi = *(pool[i]);

    for (int j = i + 1; j < pool_size; j++)
    {
      similarity_type_t d = compute_the_distance(pi, *(pool[j]));
      if (d > furthest_distance_encountered)
      {
        furthest_distance_encountered = d;
        id_of_further_distance_encountered = i;
      }
    }
  }

  istart = id_of_further_distance_encountered;

  if (verbose)
    cerr << "Starting with '" << pool[id_of_further_distance_encountered]->id() << "' dist " << furthest_distance_encountered << endl;

  return 1;
}

static int
do_start_with_object_furthest_from_first (int & istart)
{
  resizable_array<int> already_done;
  already_done.resize(start_with_object_furthest_from_first);
  similarity_type_t furthest_distance_encountered = 0.0;
  int id_of_further_distance_encountered = 0;

  for (int i = 0; i < start_with_object_furthest_from_first; i++)
  {
    already_done.add(istart);

    Spread_Object & fp0 = *(pool[istart]);

    int furthest_away = -1;

    similarity_type_t d0 = 0.0;
    for (int j = 1; j < pool_size; j++)
    {
      similarity_type_t d = compute_the_distance(*(pool[j]), fp0);

      if (d <= d0)
        continue;

      if (j == istart || already_done.contains(j))
        continue;

      d0 = d;
      furthest_away = j;
    }

    assert (furthest_away > 0);

    if (verbose)
      cerr << "Furthest from first fingerprint is " << furthest_away << " '" << pool[furthest_away]->id() << "', distance " << d0 << endl;

    istart = furthest_away;

    if (d0 > furthest_distance_encountered)
    {
      furthest_distance_encountered = d0;
      id_of_further_distance_encountered = istart;
    }
  }

  istart = id_of_further_distance_encountered;

  if (verbose)
    cerr << "Starting with '" << pool[id_of_further_distance_encountered]->id() << "' dist " << furthest_distance_encountered << endl;

  return 1;
}

static int 
furthest_from_already_selected()
{
  int rc = 0;
  similarity_type_t maxd = pool[0]->distance();

  for (int i = 1; i < pool_size; i++)
  {
    similarity_type_t d = pool[i]->distance();

    if (d <= maxd)
      continue;

    maxd = d;
    rc = i;
  }

  return rc;
}

/*
  We need a structure via which we communicate between threads
*/

struct Thread_Communication
{
  int thread_number;
  int istart;
  int istop;

// In the first phase, we will fill these two values

  int most_distant;
  similarity_type_t furthest_distance;

// In the second phase, we will receive this value from the controlling
// process

  int sel;
};

static void
thread_identify_furthest_away (struct Thread_Communication * tc)
{
  if (tc->thread_number < 0)    // special value for finished
  {
    pthread_exit(NULL);
    return;
  }

  similarity_type_t maxdist = static_cast<similarity_type_t>(-1.0);
  int ichoose = -1;

  for (int i = tc->istart; i < tc->istop; i++)
  {
    if (pool[i]->selected())
      continue;

//    cerr << "Distance to " << i << " is " << pool[i]->distance() << endl;
    if (pool[i]->distance() > maxdist)
    {
      maxdist = pool[i]->distance();
      ichoose = i;
    }
  }

  tc->most_distant = ichoose;
  tc->furthest_distance = maxdist;

  return;
}

static void
thread_object_has_been_selected (Thread_Communication * tc)
{
  if (tc->thread_number < 0)    // special value for finished
  {
    pthread_exit(NULL);
    return;
  }

  Spread_Object & fpsel = *(pool[tc->sel]);

  for (int i = tc->istart; i < tc->istop; i++)
  {
    if (pool[i]->selected())
      continue;

    pool[i]->object_has_been_selected(fpsel);
  }

  return;
}

/*
  The threads that look for the furthest away are very simple.
  Just wait until it is our turn, then scan and post
*/

static void *
thread_process (void * v)
{
  struct Thread_Communication * tc = reinterpret_cast<Thread_Communication *>(v);

  int mythread = tc->thread_number;

  while (1)
  {
    sem_wait(&(time_for_child_to_wake1[mythread]));
    thread_identify_furthest_away(tc);
    sem_post(&child_done);
    sem_wait(&(time_for_child_to_wake2[mythread]));
    thread_object_has_been_selected(tc);
    sem_post(&child_done);
  }
}

static struct Thread_Communication * tc = nullptr;

static int
initialise_threads()
{
  assert (nthreads > 0);
  assert (NULL != tc);

  int items_per_chunk = pool_size / nthreads;

  if (items_per_chunk < 1)
  {
    cerr << "Too many threads " << nthreads << " for " << pool_size << " fingerprints\n";
    return 0;
  }

  child_thread = new pthread_t[nthreads];

  pthread_attr_t attr;

  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  for (int i = 0; i < nthreads; i++)
  {
    Thread_Communication & tci = tc[i];

    tci.thread_number = i;
    tci.istart = i * items_per_chunk;
    if (i == nthreads - 1)
      tci.istop = pool_size;
    else
      tci.istop = tci.istart + items_per_chunk;

    if (0 != pthread_create(&(child_thread[i]), &attr, thread_process, &tci))
    {
      perror("cannot create child process");
      return 0;
    }
  }

  return 1;
}

static int
at_least_this_many_active (int n)
{
  int active = 0;

  for (int i = 0; i < pool_size; i++)
  {
    if (pool[i]->selected())
      continue;

    active++;

    if (active >= n)
      return 1;
  }

  return 0;
}

/*
  We can compress the pool array.
*/

static void
squeeze_selected_items ()
{
  if (! at_least_this_many_active(5 * nthreads))
    return;

  if (pool_size < 1000)
    return;

  int ndx = 0;

  for (int i = 0; i < pool_size; i++)
  {
    if (pool[i]->selected())
      delete pool[i];
    else
    {
      pool[ndx] = pool[i];
      ndx++;
    }
  }

  pool_size = ndx;

  int items_per_chunk = pool_size / nthreads;

  for (int i = 0; i < nthreads; i++)
  {
    Thread_Communication & tci = tc[i];

    tci.thread_number = i;
    tci.istart = i * items_per_chunk;
    if (i == nthreads - 1)
      tci.istop = pool_size;
    else
      tci.istop = tci.istart + items_per_chunk;
  }

  return;
}

static int
RandomNumberBetween(int zmin, int zmax) {
  std::random_device rd;
  std::uniform_int_distribution<int> u(zmin, zmax);
  return u(rd);
}

static int
fpobj_spread (IWString_and_File_Descriptor & output)
{
  int first_selected;

  if (choose_first_item_randomly)
    first_selected = RandomNumberBetween(0, pool_size - 1);
  else if (already_selected_molecules_present)
    first_selected = furthest_from_already_selected();
  else if (previously_computed_nn_distance_tag.length())
    first_selected = choose_largest_previously_computed_distance();
  else if (first_item_is_one_with_highest_scale_factor)
    first_selected = item_with_highest_scale_factor();
  else
    first_selected = 0;

  if (start_with_object_furthest_from_first)
    do_start_with_object_furthest_from_first(first_selected);
  else if (start_with_object_furthest_from_everything)
    do_start_with_object_furthest_from_everything(first_selected);

  Spread_Object & fp0 = *(pool[first_selected]);
  fp0.set_selected();
  do_object_has_been_selected(first_selected);

  output << smiles_tag << fp0.smiles() << ">\n";
  output << identifier_tag << fp0.id() << ">\n";

  const Smiles_ID_Dist & sid = fp0.nsn();

  if (already_selected_molecules_present || previously_computed_nn_distance_tag.length())
  {
    output << smiles_tag << sid.smiles() << ">\n";
    output << identifier_tag << sid.id() << ">\n";
    output << distance_tag << static_cast<float>(sid.distance() * pool[first_selected]->scale()) << ">\n";
  }
  else
  {
    output << smiles_tag << "*>\n";
    output << identifier_tag << "*>\n";
    output << distance_tag << "1>\n";
  }
  output << "|\n";

  int number_selected = 1;

  if (! initialise_threads())
  {
    cerr << "Cannot initialise threads\n";
    return 0;
  }

  while (number_selected < number_to_select)
  {
    if (! post_to_array_of_semaphores(time_for_child_to_wake1, nthreads))
      return 0;

    wait_for_children();

//  cerr << "Back from max distance determination\n";

    similarity_type_t maxdist = tc[0].furthest_distance;
    int ichoose = tc[0].most_distant;

    for (int i = 1; i < nthreads; i++)
    {
//    cerr << "Thread " << i << " found " << tc[i].furthest_distance << endl;
      if (tc[i].furthest_distance > maxdist)
      {
        maxdist = tc[i].furthest_distance;
        ichoose = tc[i].most_distant;
      }
    }

    assert (ichoose >= 0);

    Spread_Object & fpsel = *(pool[ichoose]);

    fpsel.set_selected();

    for (int i = 0; i < nthreads; i++)  // tell each of the threads the identity of the selected fp
    {
      tc[i].sel = ichoose;
    }

    if (! post_to_array_of_semaphores(time_for_child_to_wake2, nthreads))
      return 0;

//  Do output while children computing

    output << smiles_tag << fpsel.smiles() << ">\n";
    output << identifier_tag << fpsel.id() << ">\n";
    const Smiles_ID_Dist & sid = fpsel.nsn();
    output << smiles_tag << sid.smiles() << ">\n";
    output << identifier_tag << sid.id() << ">\n";
    if (static_cast<float>(1.0) != fpsel.scale())
      output << "SCALE<" << fpsel.scale() << ">\n";
    output << distance_tag << fpsel.distance() << ">\n";    // the sid object does not know about any scaling of the distance
    output << "|\n";

    if (verbose > 1)
      cerr << "Selected " << number_selected << " '" << fpsel.id() << "' distance " << fpsel.distance() << " NSN '" << sid.id() << "'\n";

    output.write_if_buffer_holds_more_than(32768);

    nearest_selected_neighbour_distance.extra(fpsel.distance());

    number_selected++;

    wait_for_children();

    if (0 == number_selected % squeeze_frequency)
      squeeze_selected_items();
  }

// all the threads to destroy themselves

  for (int i = 0; i < nthreads; i++)
  {
    tc[i].thread_number = -1;
  }

  post_to_array_of_semaphores(time_for_child_to_wake1, nthreads);

  output.flush();

  if (verbose)
    cerr << "Returning with " << number_selected << " items selected\n";

  return number_selected;
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Usage <options> <input_file>\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -n <number>      specify how many items to select\n";
  cerr << " -A <file>        specify file of already selected items\n";
  cerr << " -N <tag>         gfp_nearneighbours has been run and initial distances are in <tag>\n";
  cerr << " -N col=<c>       initial distances are column <c> of the name\n";
  cerr << " -p COL=<col>     distance scaling factor is column <col> of name\n";
  cerr << " -p <tag>         specify distance scaling factor in <tag>\n";
  cerr << " -p FILE=<fname>  distance scaling factors in <fname>\n";
  cerr << " -p oknoscale     ok if not all items have a distance scaling factor - will default to 1.0\n";
  cerr << " -r <number>      report progress of initial distance assignments\n";
//cerr << " -E ...           options for choosing first selected molecule, enter '-E help'\n";
//cerr << " -C <number>      first selected is <number> times furthest from first fingerprint (DC)\n";
//cerr << " -R               first selected is random\n";
//cerr << " -w               first selected is one with highest scaling factor\n";
//cerr << " -f               first molecule is furthest away from everything else\n";
  cerr << " -j               recompute initial distances if no near neighbours\n";
  cerr << " -h <number>      number of threads to use\n";
  cerr << " -X ...           obscure options\n";
  cerr << " -i <tag>         specify identifier tag in pool\n";
  cerr << " -I <tag>         specify identifier tag in input file\n";
  cerr << " -F,-P,-W ...     gfp options, enter '-F help' for details\n";
//cerr << " -V <...>         Tversky conditions, enter '-V help' for details\n";
  cerr << " -v               verbose output\n";

  exit(rc);
}

/*static int
parse_filter_and_file (const const_IWSubstring & a,
                       const_IWSubstring & fname,
                       IW_TDT_Filter & filter)
{
  int i1 = a.index(",FILTER:");
  assert (i1 > 0);

  const_IWSubstring fstring = a;
  fstring.remove_leading_chars(i1 + 8);

  fname = a;
  fname.iwtruncate(i1);

  cerr << "filter '" << fstring << "', fname '" << fname << "'\n";

  if (! filter.build_from_string(fstring))
  {
    cerr << "Cannot initialise filter from '" << fstring << "'\n";
    return 0;
  }

  return 1;
}*/

static int
destroy_array_of_semaphores (sem_t * sem,
                             int nsem)
{
  int rc = 1;

  for (int i = 0; i < nsem; i++)
  {
    if (0 != sem_destroy(&(sem[i])))
    {
      perror("sem_destroy error");
      rc = 0;
    }
  }

  return rc;
}

static void
display_misc_options (std::ostream & os)
{
  os << " -X sq=<n>      squeeze out selected items every <n> items selected\n";
  os << " -X blurr=<d>   blurr distances to <d> resolution\n";

  exit(1);
}

static void
display_select_first_item_options (std::ostream & os,
                                   int rc = 0)
{
  os << " -C first        first item selected is first item in input\n";
  os << " -C rand         randomly choose first item to be selected\n";
  os << " -C hscale       first item selected is one with highest scale (-p)\n";
  os << " -C ff=<n>       first item selected is <n> times furthest from first\n";

  exit(rc);
}

static int
fpobj_spread (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vs:n:I:i:A:r:p:C:fRF:P:W:Q:jN:V:wh:X:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  if (need_to_call_initialise_fingerprints(cl))
  {
    if (! initialise_fingerprints(cl, verbose))
    {
      cerr << "Cannot initialise GFP options\n";
      usage(23);
    }
  }
  else if (! initialise_fingerprints(cl[0], verbose))
  {
    cerr << "Cannot initialise fingerprints from '" << cl[0] << "'\n";
    return 11;
  }

  if (cl.option_present('V'))
  {
    if (! tversky.parse_command_line(cl, 'V', verbose))
    {
      cerr << "Cannot initialise Tversky parameters\n";
      return 8;
    }
  }

  if (! cl.option_present('h'))
  {
    cerr << "Must specify number of threads via the -h option\n";
    usage(3);
  }

  if (cl.option_present('h'))
  {
    if (! cl.value('h', nthreads) || nthreads < 1)
    {
      cerr << "The number of threads to use (-h) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will use " << nthreads << " threads\n";

    time_for_child_to_wake1 = new sem_t[nthreads];
    time_for_child_to_wake2 = new sem_t[nthreads];

    if (! already_selected_molecules_present)
    {
      if (! initialise_array_of_semaphores(time_for_child_to_wake1, nthreads))
        return 3;

      if (0 != sem_init(&child_done, 0, 0))
      {
        perror("sem init error");
        return 6;
      }
    }

    if (! initialise_array_of_semaphores(time_for_child_to_wake2, nthreads))
      return 8;

    tc = new struct Thread_Communication[nthreads];
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 1)
    {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }

    if (! allocate_pool())
      return 83;
  }

// We need to be careful with the -i and -I options. Remember
// that the pool is built first

  if (cl.option_present('i'))
  {
    const_IWSubstring id;
    (void) cl.value('i', id);

    set_identifier_tag(id);

    if (verbose)
      cerr << "Identifiers in pool tagged as '" << id << "'\n";
  }

  int scaling_factor_specified = 0;

  if (cl.option_present('p'))
  {
    IWString tag, fname, col;

    int i = 0;
    const_IWSubstring p;

    while (cl.value('p', p, i++))
    {
      if (p.starts_with("FILE=") && 0 == fname.length())
      {
        p.remove_leading_chars(5);
        fname = p;
      }
      else if (p.starts_with("COL=") && 0 == col.length())
      {
        p.remove_leading_chars(4);
        col = p;
      }
      else if ("oknoscale" == p)
      {
        set_every_object_must_have_a_scale_factor(0);
      }
      else if (0 == tag.length())   // tag can only be specified once
        tag = p;
      else
      {
        cerr << "Unrecognised -p qualifier '" << p << "'\n";
        usage(4);
      }
    }

    if (0 == tag.length() && 0 == fname.length() && 0 == col.length())
    {
      cerr << "Must specify either tag, file name or column for the weighting factor\n";
      usage(3);
    }

    if (tag.length() && fname.length() && col.length())
    {
      cerr << "Must specify just one tag, just one file name or a column for the weighting factor\n";
      usage(3);
    }

    if (tag.length())
    {
      set_scale_tag(tag);

      if (verbose)
        cerr << "The scale factor will be the '" << tag << "' dataitem\n";
    }
    else if (fname.length())
    {
      if (! read_scaling_data(fname.null_terminated_chars(), verbose))
      {
        cerr << "Cannot read scaling data from '" << fname << "'\n";
        return 3;
      }
    }
    else if (col.length())
    {
      int c;
      if (! col.numeric_value(c) || c < 2)
      {
        cerr << "Invalid scaling factor column '" << col << "'\n";
        return 4;
      }

      set_scaling_factor_column(c);
    }

    scaling_factor_specified = 1;
  }

  int previous_distance_column = -1;

  if (cl.option_present('N'))
  {
    IWString fname;

    int attributes_specified = 0;

    int i = 0;
    const_IWSubstring n;

    while (cl.value('N', n, i++))
    {
      if (n.starts_with("FILE="))
      {
        n.remove_leading_chars(5);
        fname = n;
        attributes_specified++;
      }
      else if (n.starts_with("COL="))
      {
        n.remove_leading_chars(4);
        if (! n.numeric_value(previous_distance_column) || previous_distance_column < 1)
        {
          cerr << "The previous distance column must be a valid column number '" << n << "'\n";
          return 4;
        }
        attributes_specified++;
        previous_distance_column--;
      }
      else if (0 == previously_computed_nn_distance_tag.length())
      {
        previously_computed_nn_distance_tag = n;
        if (verbose)
          cerr << "Previously computed near neighbour distances in the '" << previously_computed_nn_distance_tag << "' tag\n";

        if (! previously_computed_nn_distance_tag.ends_with('<'))
          previously_computed_nn_distance_tag += '<';

        attributes_specified++;
      }
      else
      {
        cerr << "Unrecognised -N qualifier '" << n << "'\n";
        usage(3);
      }
    }

    if (1 != attributes_specified)
    {
      cerr << "Can specify just one of FILE=, COL= or tag for previously computed distances\n";
      usage(3);
    }

    if (fname.length())
    {
      if (! previously_computed_distances.read_data(fname))
      {
        cerr << "Cannot read previously computed distances from '" << fname << "'\n";
        return 4;
      }

      if (verbose)
        cerr << "Read " << previously_computed_distances.size() << " previously computed nn distances from '" << fname << "'\n";
    }
  }

// build the pool

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! build_pool(cl[i], previously_computed_nn_distance_tag))
    {
      cerr << "Yipes, cannot build pool from '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  pool_size = poolptr + 1;

  if (previous_distance_column >= 0)
  {
    for (int i = 0; i < pool_size; i++)
    {
      const_IWSubstring id = pool[i]->id().word(previous_distance_column);

      float d;
      if (! id.numeric_value(d) || d < static_cast<float>(0.0))
      {
        cerr << "Invalid distance '" << id << "' from '" << pool[i]->id() << "'\n";
        return 3;
      }

      pool[i]->set_distance(d);
    }
  }


  if (verbose && scaling_factor_specified)
  {
    const Accumulator<float> & sfs = scale_factor_statistics();
    cerr << sfs.n() << " of " << pool_size << " pool objects had demerit/scale factors\n";

    if (sfs.n() > 0)
    {
      cerr << "Scale factors between " << sfs.minval() << " and " << sfs.maxval();
      if (sfs.n() > 1)
        cerr << ", ave " << sfs.average();
      cerr << endl;
    }
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

  if (cl.option_present('r'))
  {
    if (! cl.value('r', report_establish_initial_distances) || report_establish_initial_distances < 1)
    {
      cerr << "The -r option must be followed by a whole positive number\n";
      usage(18);
    }

    if (verbose)
      cerr << "Will report initial neighbour assignments every " << report_establish_initial_distances << " fingerprints\n";
  }

  if (cl.option_present('j'))
  {
    retest_no_neighbours = 1;
    if (verbose)
      cerr << "Molecules with no initial nearest neighbour will be recomputed without window\n";
  }

// There can be any number of already present files.

  if (cl.option_present('A'))
  {
    set_scaling_factor_column(-1);   // turn off scaling factor stuff
    set_scale_tag("");
    set_every_object_must_have_a_scale_factor(0);

    const_IWSubstring fname;
    int i = 0;
    while (cl.value('A', fname, i++))
    {
      if (! establish_initial_distances(fname))
      {
        cerr << "Cannot establish initial distances from '" << fname << "'\n";
        return 54;
      }
    }

    already_selected_molecules_present = 1;

    for (int i = 0; i < nthreads; i++)
    { 
      pstc[i].thread_number = -1;
      pstc[i].previously_selected_fp = nullptr;
    }

    post_to_array_of_semaphores(time_for_child_to_wake1, nthreads);
    for (int i = 0; i < nthreads; i++)
    {
      sem_wait(&child_done);
    }

    delete [] pstc;
  }
  else
  {
    already_selected_molecules_present = 0;
  }

//There must be only one means of selecting the first molecule
  
  if (0 == cl.option_count('C'))
    ;
  else if (cl.option_count('C') > 1)
  {
    cerr << "Only one means of selecting the first item (-C) is alowed\n";
    usage(3);
  }
  else
  {
    if (already_selected_molecules_present)
    {
      cerr << "The -C option only works in the absence of the -A option\n";
      usage(18);
    }

    const_IWSubstring c = cl.string_value('C');
    if ("first" == c)
    {
      if (verbose)
        cerr << "First item in input will be selected first\n";
    }
    else if (c.starts_with("ff="))
    {
      c.remove_leading_chars(3);
      if (! c.numeric_value(start_with_object_furthest_from_first) || start_with_object_furthest_from_first < 1)
      {
        cerr << "Start furthest from first must be a whole +ve number\n";
        display_select_first_item_options(cerr, 1);
      }
    }
    else if (c.starts_with("rand"))
    {
      choose_first_item_randomly = 1;
      if (verbose)
        cerr << "Will choose first fingerprint at random\n";
    }
    else if (c.starts_with("hscale"))
    {
      if (! cl.option_present('p'))
      {
        cerr << "To start with the molecule with the highest scaling factor, must specify a scaling factor (-p)\n";
        usage(4);
      }

      first_item_is_one_with_highest_scale_factor = 1;
      if (verbose)
        cerr << "First selected will be molecule with highest scale factor\n";
    }
    else if ("help" == c)
    {
      display_select_first_item_options(cerr);
    }
    else
    {
      cerr << "Unrecognised -C qualifier '" << c << "'\n";
      display_select_first_item_options(cerr, 2);
    }
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
      cerr << "You asked for " << n << " molecules, but pool only contains " << pool_size << ". Shortened\n";
      n = pool_size;
    }

    number_to_select = n;
    if (verbose)
      cerr << number_to_select << " molecules will be selected\n";
  }
  else
    number_to_select = pool_size;

  if (cl.option_present('X'))
  {
    int i = 0;
    const_IWSubstring x;
    std::random_device rd;
    std::uniform_real_distribution<float> u(0.0, 1.0);
    while (cl.value('X', x, i++))
    {
      if (x.starts_with("blurr="))
      {
        x.remove_leading_chars(6);
        if (! x.numeric_value(blurr_distances) || blurr_distances < static_cast<float>(0.0))
        {
          cerr << "The blurr distances option must be a non negative number\n";
          usage(4);
        }

        if (verbose)
        {
          cerr << "Distance blurring factor set to " << blurr_distances << '\n';
          for (int i = 0; i < 5; i++)
          {
            similarity_type_t r = u(rd);
            cerr << "distance " << r << " becomes " << do_blurring(r, blurr_distances) << '\n';
          }
        }
      }
      else if (x.starts_with("sq="))
      {
        x.remove_leading_chars(3);
        if (! x.numeric_value(squeeze_frequency) || squeeze_frequency < 1)
        {
          cerr << "Squeeze frequency must be a whole +ve number\n";
          return 2;
        }
        if (verbose)
          cerr << "Will squeeze out selected items every " << squeeze_frequency << " items selected\n";
      }
      else if ("help" == x)
        display_misc_options(cerr);
      else
      {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        display_misc_options(cerr);
      }
    }
  }

  IWString_and_File_Descriptor output(1);

  (void) fpobj_spread(output);

  if (nthreads)
  {
    destroy_array_of_semaphores(time_for_child_to_wake1, nthreads);
    destroy_array_of_semaphores(time_for_child_to_wake2, nthreads);
    (void) sem_destroy(&child_done);
  }

  if (verbose)
  {
    cerr << "Nearest previously selected item distances between " << nearest_selected_neighbour_distance.minval() << " and " << nearest_selected_neighbour_distance.maxval();
    if (nearest_selected_neighbour_distance.n() > 1)
      cerr << " ave " << nearest_selected_neighbour_distance.average();
    cerr << endl;
  }

//delete [] pool;     leave this out for efficiency

  cerr << "Output can be processed with nplotnn\n";

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = fpobj_spread(argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status(stderr);
#endif

  return rc;
}
