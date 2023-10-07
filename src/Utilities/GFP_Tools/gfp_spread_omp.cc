/*
  Spread implementation
*/

#include <omp.h>
#include <stdlib.h>

#include <iostream>
#include <random>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/numeric_data_from_file.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Utilities/GFP_Tools/spread_v2.h"
#include "Utilities/GFP_Tools/tversky.h"

using std::cerr;

static int verbose = 0;

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

static Report_Progress report_establish_initial_distances;

/*
  Our pool is an array of FP objects
*/

static Spread_Object* pool = nullptr;

static int pool_size = 0;

static int poolptr = 0;  // the next item in the pool to be filled

static int number_to_select = 0;

static similarity_type_t stop_once_distance_drops_below = 0.0f;

/*
  We keep track of the distances of the nearest selected items
*/

static Accumulator<similarity_type_t> nearest_selected_neighbour_distance;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");

static IWString previously_computed_nn_distance_tag;

static int previous_computed_distance_column = -1;

static Numeric_Data_From_File<float> previously_computed_distances;

/*
  How often do we squeeze out selected items
*/

static int squeeze_pool = 0;

/*
  The normal output of the near neighbour programmes leaves $SMI and PCN tag in the
  output file. Those would confuse things, so they must be transformed to something
  else. We assume that has been done
*/

static IWString nn_smiles_tag("NNSMI<");
static IWString nn_id_tag("NNID<");

static float blurr_distances = static_cast<float>(0.0);

static similarity_type_t longest_distance_recognised =
    static_cast<similarity_type_t>(0.0);

static int
get_previously_computed_nearest_neighbour(
    Spread_Object& p, const IW_STL_Hash_Map_float& previously_computed_distances)
{
  IW_STL_Hash_Map_float::const_iterator f = previously_computed_distances.find(p.id());

  if (f ==
      previously_computed_distances.end()) {  // OK if no previously computed distance
    return 1;
  }

  p.set_nearest_previously_selected_neighbour("C", "UNK", (*f).second);

  return 1;
}

static int
get_previously_computed_nearest_neighbour(
    const IW_TDT& tdt, Spread_Object& p,
    const IWString& previously_computed_nn_distance_tag)
{
  similarity_type_t d;
  if (!tdt.dataitem_value(previously_computed_nn_distance_tag, d)) {
    return 0;
  }

  IWString nnsmiles;
  if (!tdt.dataitem_value(nn_smiles_tag, nnsmiles)) {
    return 0;
  }

  IWString nnid;
  if (!tdt.dataitem_value(nn_id_tag, nnid)) {
    return 0;
  }

  p.set_nearest_previously_selected_neighbour(nnsmiles, nnid, d);

  // cerr << "Set distance to previous " << p.distance() << '\n';

  return 1;
}

static int
build_pool(iwstring_data_source& input,
           const IWString& previously_computed_nn_distance_tag)
{
  assert(pool_size > 0);
  // cerr << "Pool ptr " << poolptr << ", pool size " << pool_size << '\n';
  assert(poolptr >= 0 && poolptr < pool_size);

  int items_with_previously_computed_distances = 0;

  IW_TDT tdt;
  while (tdt.next(input)) {
    int fatal;
    if (!pool[poolptr].construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      continue;
    }

    if (previously_computed_nn_distance_tag.length()) {
      get_previously_computed_nearest_neighbour(
          tdt, pool[poolptr],
          previously_computed_nn_distance_tag);  // should really check for fatal errors
    } else if (previously_computed_distances.size()) {
      if (!get_previously_computed_nearest_neighbour(pool[poolptr],
                                                     previously_computed_distances)) {
        return 0;
      }
    } else if (previous_computed_distance_column > 0) {
      if (!pool[poolptr].set_distance_to_previously_selected_from_column(
              previous_computed_distance_column)) {
        cerr << "Cannot set previously selected distance based on column '"
             << pool[poolptr].id() << "'\n";
        return 0;
      }
    }

    if (pool[poolptr].has_a_nearest_selected_neighbour()) {
      items_with_previously_computed_distances++;
    }

    poolptr++;

    if (poolptr >= pool_size) {
      if (verbose) {
        cerr << "Pool is full, max " << pool_size << '\n';
      }
      break;
    }
  }

  poolptr--;

  if (verbose) {
    cerr << "Pool now contains " << (poolptr + 1) << " objects\n";
    if (previously_computed_nn_distance_tag.length()) {
      cerr << items_with_previously_computed_distances
           << " items had previously computed distances\n";
    }
  }

  return 1;
}

static int
allocate_pool()
{
  assert(pool_size > 0);
  assert(NULL == pool);

  pool = new Spread_Object[pool_size];
  if (NULL == pool) {
    cerr << "Yipes, could not allocate pool of size " << pool_size << '\n';
    return 62;
  }

  if (verbose) {
    cerr << "system sized to " << pool_size << '\n';
  }

  return 1;
}

static int
build_pool(const char* fname, const IWString& previously_computed_nn_distance_tag)
{
  iwstring_data_source input(fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (0 == pool_size) {
    pool_size = input.count_records_starting_with("PCN<");

    if (0 == pool_size) {
      cerr << "Zero occurrences of '"
           << "PCN<"
           << "' in '" << fname << "'\n";
      return 0;
    }

    if (!allocate_pool()) {
      return 0;
    }
  }

  return build_pool(input, previously_computed_nn_distance_tag);
}

static int
do_squeeze_pool(Spread_Object* pool, int& pool_size)
{
  int iptr = 0;
  for (int i = 0; i < pool_size; i++) {
    if (pool[i].selected()) {
      continue;
    }

    if (iptr != i) {
      pool[iptr] = pool[i];
    }

    iptr++;
  }

  pool_size = iptr;

  return 1;
}

/*
  After establishing initial distances, some pool members don't have a nearest selected
  neighbour. re-scan the pool
*/

#ifdef NOT_USED_QWLEQWE
static int
rescan_for_no_neighbours(iwstring_data_source& input, IW_TDT_Filter& filter)
{
  resizable_array<Spread_Object*> to_scan;
  to_scan.resize(pool_size);

  for (int i = 0; i < pool_size; i++) {
    if (!pool[i].has_a_nearest_selected_neighbour()) {
      to_scan.add(&(pool[i]));
    }
  }

  const int nts = to_scan.number_elements();

  if (0 == nts) {
    return 1;
  }

  if (verbose) {
    cerr << "After reading previously selected, " << nts << " items with no neighbour\n";
  }

  if (!input.seekg(0)) {
    cerr << "rescan_for_no_neighbours: yipes, cannot seek back to beginning of file\n";
    return 0;
  }

  IW_TDT tdt;

  while (tdt.next(input)) {
    int fatal;
    Spread_Object fp;
    if (!fp.construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        cerr << "Cannot build fingerpint\n" << tdt << '\n';
        return 0;
      }

      continue;
    }

    for (int i = 0; i < to_scan.number_elements(); i++) {
      Spread_Object* p = to_scan[i];

      if (tversky.active()) {
        p->object_has_been_selected(fp, tversky);
      } else {
        p->object_has_been_selected(fp);
      }
    }
  }

  return 1;
}
#endif

static int
establish_initial_distances(iwstring_data_source& input)
{
  int ntdt = 0;

  IW_TDT tdt;
  while (tdt.next(input)) {
    int fatal;
    Spread_Object fp;
    if (!fp.construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        cerr << "Cannot build fingerpint\n" << tdt << '\n';
        return 0;
      }

      continue;
    }

    ntdt++;

    //  Run this fingerprint against everything in the pool

#pragma omp parallel
    {
#pragma omp for
      for (int i = 0; i < pool_size; i++) {
        if (tversky.active()) {
          pool[i].object_has_been_selected(fp, tversky);
        } else if (static_cast<similarity_type_t>(0.0) != longest_distance_recognised) {
          pool[i].object_has_been_selected_max_distance(fp, longest_distance_recognised);
        } else {
          pool[i].object_has_been_selected(fp);
        }
      }
    }

    if (report_establish_initial_distances()) {
      cerr << "Established initial distances " << ntdt << '\n';
    }
  }

  if (0 == ntdt) {
    cerr << "establish_initial_distances:: warning, no TDT's read\n";
  }

  if (verbose > 1) {
    cerr << "INitial distances established\n";
    for (int i = 0; i < pool_size; i++) {
      cerr << "Pool " << i << " is " << pool[i].id() << " dist " << pool[i].distance()
           << '\n';
    }
  }

  return 1;
}

static int
establish_initial_distances(const const_IWSubstring& fname)
{
  static int first_call = 1;

  if (first_call) {
    for (int i = 0; i < pool_size; i++) {
      pool[i].set_distance(2.0);
    }

    first_call = 0;
  }

  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open already selected file '" << fname << "'\n";
    return 0;
  }

  if (verbose) {
    cerr << "Establishing initial distances wrt '" << fname << "'\n";
  }

  return establish_initial_distances(input);
}

/*
  We use the atom count windows in an unusual manner.
  Usually atom count windows are used to NOT compare molecules
  that have very different atom counts.
  But when doing maximum dis-similarity selections, we actually
  want to find things at large distances. Therefore, invert
  the use of the atom count window. Initially, only consider
  molecules that have very different atom counts. When it
  seems no longer possible to do that, stop using the atom
  count windows
*/

static int use_atom_count_window = 0;

static void
do_object_has_been_selected_tversky(int isel)
{
  Spread_Object& fpsel = pool[isel];

  for (int i = 0; i < pool_size; i++) {
    if (pool[i].selected() || i == isel) {
      continue;
    }

    if (static_cast<float>(0.0) != blurr_distances) {
      pool[i].object_has_been_selected(fpsel, tversky, blurr_distances);
    } else {
      pool[i].object_has_been_selected(fpsel, tversky);
    }
  }

  return;
}

static void
do_object_has_been_selected_no_blurring(int isel)
{
  Spread_Object& fpsel = pool[isel];

  for (int i = 0; i < pool_size; i++) {
    if (pool[i].selected() || i == isel) {
      continue;
    }

    pool[i].object_has_been_selected(fpsel);
  }

  return;
}

static void
do_object_has_been_selected_with_blurring(int isel)
{
  Spread_Object& fpsel = pool[isel];

  for (int i = 0; i < pool_size; i++) {
    if (pool[i].selected() || i == isel) {
      continue;
    }

    pool[i].object_has_been_selected(fpsel, blurr_distances);
  }

  return;
}

static void
do_object_has_been_selected_with_distance_cutoff(int isel)
{
  Spread_Object& fpsel = pool[isel];
  for (int i = 0; i < pool_size; i++) {
    if (pool[i].selected() || i == isel) {
      continue;
    }

    pool[i].object_has_been_selected_max_distance(fpsel, longest_distance_recognised);
  }

  return;
}

static void
do_object_has_been_selected(int isel)
{
  if (tversky.active()) {
    do_object_has_been_selected_tversky(isel);
  } else if (static_cast<float>(0.0) != blurr_distances) {
    do_object_has_been_selected_with_blurring(isel);
  } else if (static_cast<float>(0.0) != longest_distance_recognised) {
    do_object_has_been_selected_with_distance_cutoff(isel);
  } else {
    do_object_has_been_selected_no_blurring(isel);
  }

  return;
}

static similarity_type_t
compute_the_distance(Spread_Object& fp1, Spread_Object& fp2)
{
  if (tversky.active()) {
    return static_cast<similarity_type_t>(1.0) -
           fp1.IW_General_Fingerprint::tversky(fp2, tversky);
  }

  return fp1.IW_General_Fingerprint::distance(fp2);
}

static int
choose_largest_previously_computed_distance()
{
  int rc = -1;
  similarity_type_t dmax = static_cast<similarity_type_t>(0.0);

  for (int i = 0; i < pool_size; i++) {
    similarity_type_t d = pool[i].distance();

    if (static_cast<similarity_type_t>(1.0) == d) {
      continue;
    }

    if (d > dmax) {
      rc = i;
      dmax = d;
    }
  }

  if (rc < 0) {
    cerr << "Warning, none of " << pool_size
         << " items have previously computed distances!\n";
    rc = 0;
  }

  return rc;
}

/*
  No checks as to whether things have scaling factors or not
*/

static int
item_with_highest_scale_factor()
{
  float highest_scale = pool[0].scale();
  int rc = 0;

  for (int i = 1; i < pool_size; i++) {
    if (pool[i].scale() > highest_scale) {
      highest_scale = pool[i].scale();
      rc = i;
    }
  }

  return rc;
}

/*
  Since Tversky is asymmetric, we should probably do both loops 1-pool_size
*/

static int
do_start_with_object_furthest_from_everything(int& istart)
{
  int id_of_further_distance_encountered = -1;
  similarity_type_t furthest_distance_encountered = 0.0;

  for (int i = 0; i < pool_size; i++) {
    Spread_Object& pi = pool[i];

    for (int j = i + 1; j < pool_size; j++) {
      similarity_type_t d = compute_the_distance(pi, pool[j]);
      if (d > furthest_distance_encountered) {
        furthest_distance_encountered = d;
        id_of_further_distance_encountered = i;
      }
    }
  }

  istart = id_of_further_distance_encountered;

  if (verbose) {
    cerr << "Starting with '" << pool[id_of_further_distance_encountered].id()
         << "' dist " << furthest_distance_encountered << '\n';
  }

  return 1;
}

static int
do_start_with_object_furthest_from_first(int& istart)
{
  resizable_array<int> already_done;
  already_done.resize(start_with_object_furthest_from_first);
  similarity_type_t furthest_distance_encountered = 0.0;
  int id_of_further_distance_encountered = 0;

  for (int i = 0; i < start_with_object_furthest_from_first; i++) {
    already_done.add(istart);

    Spread_Object& fp0 = pool[istart];

    int furthest_away = -1;

    similarity_type_t d0 = 0.0;
    for (int j = 1; j < pool_size; j++) {
      similarity_type_t d = compute_the_distance(pool[j], fp0);

      if (d <= d0) {
        continue;
      }

      if (j == istart || already_done.contains(j)) {
        continue;
      }

      d0 = d;
      furthest_away = j;
    }

    assert(furthest_away > 0);

    if (verbose) {
      cerr << "Furthest from first fingerprint is " << furthest_away << " '"
           << pool[furthest_away].id() << "', distance " << d0 << '\n';
    }

    istart = furthest_away;

    if (d0 > furthest_distance_encountered) {
      furthest_distance_encountered = d0;
      id_of_further_distance_encountered = istart;
    }
  }

  istart = id_of_further_distance_encountered;

  if (verbose) {
    cerr << "Starting with '" << pool[id_of_further_distance_encountered].id()
         << "' dist " << furthest_distance_encountered << '\n';
  }

  return 1;
}

static int
furthest_from_already_selected()
{
  int rc = 0;
  similarity_type_t maxd = pool[0].distance();

  for (int i = 1; i < pool_size; i++) {
    similarity_type_t d = pool[i].distance();
    //  cerr << "furthest_from_already_selected check " << d << " cmp " << maxd << '\n';

    if (d <= maxd) {
      continue;
    }

    maxd = d;
    rc = i;
  }

  // cerr << "Longest dfistance " << maxd << " item " << rc << '\n';

  return rc;
}

static int
fpobj_spread(IWString_and_File_Descriptor& output)
{
  int first_selected;

  if (choose_first_item_randomly) {
    std::random_device rd;
    std::uniform_int_distribution<int> u(0, pool_size - 1);
    first_selected = u(rd);
  } else if (already_selected_molecules_present) {
    first_selected = furthest_from_already_selected();
  } else if (previously_computed_nn_distance_tag.length()) {
    first_selected = choose_largest_previously_computed_distance();
  } else if (first_item_is_one_with_highest_scale_factor) {
    first_selected = item_with_highest_scale_factor();
  } else {
    first_selected = 0;
  }

  if (start_with_object_furthest_from_first) {
    do_start_with_object_furthest_from_first(first_selected);
  } else if (start_with_object_furthest_from_everything) {
    do_start_with_object_furthest_from_everything(first_selected);
  }

  Spread_Object& fp0 = pool[first_selected];
  fp0.set_selected();
  do_object_has_been_selected(first_selected);

  if (verbose > 1) {
    cerr << "First selected '" << fp0.id() << "'\n";
  }

  output << smiles_tag << fp0.smiles() << ">\n";
  output << identifier_tag << fp0.id() << ">\n";

  const Smiles_ID_Dist& sid = fp0.nsn();

  if (already_selected_molecules_present ||
      previously_computed_nn_distance_tag.length()) {
    output << smiles_tag << sid.smiles() << ">\n";
    output << identifier_tag << sid.id() << ">\n";
    output << distance_tag
           << static_cast<float>(sid.distance() * pool[first_selected].scale()) << ">\n";
  } else {
    output << smiles_tag << "*>\n";
    output << identifier_tag << "*>\n";
    output << distance_tag << "1>\n";
  }
  output << "|\n";

  int number_selected = 1;

  while (number_selected < number_to_select) {
    int max_pool = -1;
    float max_pool_value = -1.0f;  // just something less than 0
#pragma omp parallel
    {
      int max_pool_private = -1;
      float max_pool_value_private = -1.0f;  // again need to be less than 0
#pragma omp for nowait
      for (int i = 0; i < pool_size; i++) {  // find max value within thread
        if (pool[i].selected()) {
          continue;
        }

        if (pool[i].distance() > max_pool_value_private) {
          max_pool_private = i;
          max_pool_value_private = pool[i].distance();
        }
      }
#pragma omp critical
      {  // reduce across threads
        if (max_pool_value_private > max_pool_value) {
          max_pool = max_pool_private;
          max_pool_value = max_pool_value_private;
        }  // pragma omp critical
      }
    }  // pragma omp parallel

    const int ichoose = max_pool;
    if (ichoose < 0) {
      break;
    }

    //  cerr << "next " << ichoose << " '" << pool[ichoose].id() << "' sel " <<
    //  pool[ichoose].selected() << '\n';

    assert(ichoose >= 0);

    Spread_Object& fpsel = pool[ichoose];

    output << smiles_tag << fpsel.smiles() << ">\n";
    output << identifier_tag << fpsel.id() << ">\n";
    const Smiles_ID_Dist& sid = fpsel.nsn();
    output << smiles_tag << sid.smiles() << ">\n";
    output << identifier_tag << sid.id() << ">\n";
    if (static_cast<float>(1.0) != fpsel.scale()) {
      output << "SCALE<" << fpsel.scale() << ">\n";
    }
    output << distance_tag << fpsel.distance()
           << ">\n";  // the sid object does not know about any scaling of the distance
    output << "|\n";

    if (verbose > 1) {
      cerr << "Selected " << number_selected << " '" << fpsel.id() << "' (index "
           << ichoose << ") distance " << fpsel.distance() << " NSN '" << sid.id()
           << "'\n";
    }

    output.write_if_buffer_holds_more_than(32768);

    nearest_selected_neighbour_distance.extra(fpsel.distance());

    number_selected++;

    fpsel.set_selected();

#pragma omp parallel
    {
#pragma omp for
      for (int i = 0; i < pool_size; i++) {
        if (pool[i].selected()) {
          continue;
        }

        if (static_cast<float>(0.0) != blurr_distances) {
          pool[i].object_has_been_selected(fpsel, blurr_distances);
        } else if (static_cast<similarity_type_t>(0.0) != longest_distance_recognised) {
          pool[i].object_has_been_selected_max_distance(fpsel,
                                                        longest_distance_recognised);
        } else {
          pool[i].object_has_been_selected(fpsel);
        }
      }
    }  // pragma omp parallel

    if (squeeze_pool && 0 == number_selected % squeeze_pool) {
      do_squeeze_pool(pool, pool_size);
    }

    if (fpsel.distance() < stop_once_distance_drops_below) {
      break;
    }
  }

  if (verbose) {
    cerr << "Returning with " << number_selected << " items selected\n";
  }

  return number_selected;
}

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Usage <options> <input_file>\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -n <number>      specify how many items to select\n";
  cerr << " -t <dist>        stop selection once distance drops below <dist>\n";
  cerr << " -A <file>        specify file of already selected items\n";
  cerr << " -N <tag>         gfp_nearneighbours has been run and initial distances are in <tag>\n";
  cerr << " -p COL=<col>     distance scaling factor is column <col> of name\n";
  cerr << " -p <tag>         specify distance scaling factor in <tag>\n";
  cerr << " -p FILE=<fname>  distance scaling factors in <fname>\n";
  cerr << " -p oknoscale     ok if not all items have a distance scaling factor - will default to 1.0\n";
  cerr << " -r <number>      report progress of initial distance assignments\n";
  cerr << " -S ...           options for specifying first item selected, enter '-S help'\n";
  cerr << " -M ...           miscellaneous options, enter '-M help' for info\n";
  cerr << " -b <number>      \"blurr\" distances to <number> resolution per range\n";
//cerr << " -i <tag>         specify identifier tag in pool\n";
//cerr << " -I <tag>         specify identifier tag in input file\n";
  cerr << " -F,-P,...        gfp options, enter '-F help' for details (note -W not supported here)\n";
  cerr << " -V <...>         Tversky conditions, enter '-V help' for details\n";
  cerr << " -C <n>           number of OMP workers\n";
  cerr << " -v               verbose output\n";
  // clang-format on

  exit(rc);
}

static void
display_miscellaneous_options(std::ostream& os)
{
  os << " -M recomp       recompute distance if no nbrs found\n";
  os << " -M nscale       include scale factor of nbr with scale\n";
  os << " -M squeeze=nnn  squeeze out selected molecules evern <nnn> steps\n";
  os << " -M ldist=<dist> all distances truncated to <dist>\n";

  exit(3);
}

static void
display_first_item_selection_options(std::ostream& os)
{
  os << " -S rand         randomly choose first item\n";
  os << " -S hsf          start with item with highest scale factor\n";
  os << " -S furthest     start with item furthest from all other items\n";
  os << " -S fff=nnn      start with item nnn times furthest from first item\n";

  exit(1);
}

static int
fpobj_spread(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vs:n:I:i:A:r:p:F:P:Q:N:V:b:S:M:C:t:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('W')) {
    use_atom_count_window = 1;

    if (verbose) {
      cerr << "Atom count window specified, special processing for spread\n";
    }
  }

  if (cl.option_present('V')) {
    if (!tversky.parse_command_line(cl, 'V', verbose)) {
      cerr << "Cannot initialise Tversky parameters\n";
      return 8;
    }
  }

  if (cl.option_present('C')) {
    int nworkers;
    if (!cl.value('C', nworkers) || nworkers < 1) {
      cerr << "The number of OMP workers(-C) must be a whole +ve number\n";
      usage(1);
    }

    omp_set_num_threads(nworkers);
    //
    // IWString tmp;
    // tmp << nworkers;
    // __cilkrts_set_param("nworkers",tmp.null_terminated_chars());
  }

  if (cl.option_present('M')) {
    int i = 0;
    IWString m;
    while (cl.value('M', m, i++)) {
      m.to_lowercase();
      if ("nscale" == m) {
        if (!cl.option_present('p')) {
          cerr << "NO scaling data specified (-p) so cannot use 'nscale'\n";
          return 3;
        }

        set_include_scale_of_nsn_with_scale(1);
        if (verbose) {
          cerr << "Will include scaling factor of nearest selected nbr with scale\n";
        }
      } else if (m.starts_with("squeeze=")) {
        m.remove_leading_chars(8);
        if (!m.numeric_value(squeeze_pool) || squeeze_pool < 1) {
          cerr << "The squeeze= directive must be a +ve number\n";
          exit(3);
        }

        if (verbose) {
          cerr << "Already selected molecules will be squeezed out every " << squeeze_pool
               << " selections\n";
        }
      } else if (m.starts_with("ldist=")) {
        m.remove_leading_chars(6);
        if (!m.numeric_value(longest_distance_recognised) ||
            longest_distance_recognised <= 0.0f || longest_distance_recognised > 1.0f) {
          cerr << "The longest recognised distance directive 'ldist=' must be a valid "
                  "distance\n";
          exit(2);
        }

        if (verbose) {
          cerr << "Long distances truncated to " << longest_distance_recognised << '\n';
        }
      } else if ("help" == m) {
        display_miscellaneous_options(cerr);
      } else {
        cerr << "Unrecognised -M qualifier '" << m << "\n";
        display_miscellaneous_options(cerr);
      }
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  if (need_to_call_initialise_fingerprints(cl)) {
    if (!initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise GFP options\n";
      usage(23);
    }
  } else if (!initialise_fingerprints(cl[0], verbose)) {
    cerr << "Cannot initialise fingerprints from '" << cl[0] << "'\n";
    return 11;
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 1) {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }

    if (!allocate_pool()) {
      return 83;
    }
  }

  // We need to be careful with the -i and -I options. Remember
  // that the pool is built first

  if (cl.option_present('i')) {
    const_IWSubstring id;
    (void)cl.value('i', id);

    set_identifier_tag(id);

    if (verbose) {
      cerr << "Identifiers in pool tagged as '" << id << "'\n";
    }
  }

  int scaling_factor_specified = 0;

  if (cl.option_present('p')) {
    IWString tag, fname, col;

    int i = 0;
    const_IWSubstring p;

    while (cl.value('p', p, i++)) {
      if (p.starts_with("FILE=") && 0 == fname.length()) {
        p.remove_leading_chars(5);
        fname = p;
      } else if (p.starts_with("COL=") && 0 == col.length()) {
        p.remove_leading_chars(4);
        col = p;
      } else if ("oknoscale" == p) {
        set_every_object_must_have_a_scale_factor(0);
      } else if (0 == tag.length()) {  // tag can only be specified once
        tag = p;
      } else {
        cerr << "Unrecognised -p qualifier '" << p << "'\n";
        usage(4);
      }
    }

    if (0 == tag.length() && 0 == fname.length() && 0 == col.length()) {
      cerr << "Must specify either tag, file name or column for the weighting factor\n";
      usage(3);
    }

    if (tag.length() && fname.length() && col.length()) {
      cerr << "Must specify just one tag, just one file name or a column for the "
              "weighting factor\n";
      usage(3);
    }

    if (tag.length()) {
      set_scale_tag(tag);

      if (verbose) {
        cerr << "The scale factor will be the '" << tag << "' dataitem\n";
      }
    } else if (fname.length()) {
      if (!read_scaling_data(fname.null_terminated_chars(), verbose)) {
        cerr << "Cannot read scaling data from '" << fname << "'\n";
        return 3;
      }
    } else if (col.length()) {
      int c;
      if (!col.numeric_value(c) || c < 2) {
        cerr << "Invalid scaling factor column '" << col << "'\n";
        return 4;
      }

      set_scaling_factor_column(c);
    }

    scaling_factor_specified = 1;
  }

  if (cl.option_present('N')) {
    int nset = 0;

    IWString fname;

    int i = 0;
    const_IWSubstring n;

    while (cl.value('N', n, i++)) {
      if (n.starts_with("FILE=")) {
        n.remove_leading_chars(5);
        fname = n;
        nset++;
      } else if (n.starts_with("COL=")) {
        n.remove_leading_chars(4);
        if (!n.numeric_value(previous_computed_distance_column) ||
            previous_computed_distance_column < 1) {
          cerr << "The previously computed distance column must be a whole +ve number\n";
          return 3;
        }
        if (verbose) {
          cerr << "Previously computed near neighbour distances in column "
               << previous_computed_distance_column << '\n';
        }

        previous_computed_distance_column--;
        nset++;
      } else if (0 == previously_computed_nn_distance_tag.length()) {
        previously_computed_nn_distance_tag = n;
        if (verbose) {
          cerr << "Previously computed near neighbour distances in the '"
               << previously_computed_nn_distance_tag << "' tag\n";
        }

        if (!previously_computed_nn_distance_tag.ends_with('<')) {
          previously_computed_nn_distance_tag += '<';
        }

        nset++;
      } else {
        cerr << "Unrecognised -N qualifier '" << n << "'\n";
        usage(3);
      }
    }

    if (nset > 1) {
      cerr << "Can specify just one of FILE=, COL= or tag for previously computed "
              "distances\n";
      usage(3);
    }

    if (fname.length()) {
      if (!previously_computed_distances.read_data(fname)) {
        cerr << "Cannot read previously computed distances from '" << fname << "'\n";
        return 4;
      }

      if (verbose) {
        cerr << "Read " << previously_computed_distances.size()
             << " previously computed nn distances from '" << fname << "'\n";
      }
    }
  }

  // build the pool

  for (int i = 0; i < cl.number_elements(); i++) {
    if (!build_pool(cl[i], previously_computed_nn_distance_tag)) {
      cerr << "Yipes, cannot build pool from '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  pool_size = poolptr + 1;

  if (0 == squeeze_pool) {
    squeeze_pool = pool_size / 10 + 1;
  }

  if (verbose && scaling_factor_specified) {
    const Accumulator<float>& sfs = scale_factor_statistics();
    cerr << sfs.n() << " of " << pool_size << " pool objects had demerit/scale factors\n";

    if (sfs.n() > 0) {
      cerr << "Scale factors between " << sfs.minval() << " and " << sfs.maxval();
      if (sfs.n() > 1) {
        cerr << ", ave " << sfs.average();
      }
      cerr << '\n';
    }
  }

  // Now that the pool is built, we can switch identifiers if needed

  if (cl.option_present('I')) {
    const_IWSubstring id;
    cl.value('I', id);

    set_identifier_tag(id);

    if (verbose) {
      cerr << "Identifiers in input tagged as '" << id << "'\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_establish_initial_distances.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise reporter for progress on initial distances\n";
      usage(1);
    }
  }

  // There can be any number of already present files.
  // Each can have its own filter if the name starts with FILTER:

  if (cl.option_present('A')) {
    set_scaling_factor_column(-1);  // turn off scaling factor stuff
    set_scale_tag("");
    set_every_object_must_have_a_scale_factor(0);

    const_IWSubstring a;
    int i = 0;
    while (cl.value('A', a, i++)) {
      const_IWSubstring fname = a;  // by default

      if (!establish_initial_distances(fname)) {
        cerr << "Cannot establish initial distances from '" << fname << "'\n";
        return 54;
      }
    }

    already_selected_molecules_present = 1;
  } else {
    already_selected_molecules_present = 0;
  }

  if (cl.option_present('S')) {
    if (cl.option_count('S') > 1) {
      cerr << "There can be only one means of selecting the first item (-S)\n";
      usage(3);
    }

    if (already_selected_molecules_present) {
      cerr << "Cannot specify how to select first molecule if already selected molecules "
              "present (-A)\n";
      usage(3);
    }

    IWString s = cl.string_value('S');
    s.to_lowercase();

    if ("rand" == s) {
      choose_first_item_randomly = 1;
      if (verbose) {
        cerr << "Will choose the first item randomly";
      }
    } else if ("furthest" == s) {
      start_with_object_furthest_from_everything = 1;
      if (verbose) {
        cerr << "Will start with the item furthest from all other items\n";
      }
    } else if ("hsf" == s) {
      first_item_is_one_with_highest_scale_factor = 1;
      if (verbose) {
        cerr << "First selected will be molecule with highest scale factor\n";
      }
    } else if (s.starts_with("fff=")) {
      s.remove_leading_chars(4);
      if (!s.numeric_value(start_with_object_furthest_from_first) ||
          start_with_object_furthest_from_first < 1) {
        cerr << "The times furthest from first option 'fff=nnn' must be a +ve number\n";
        display_first_item_selection_options(cerr);
      }
      start_with_object_furthest_from_first = 1;

      if (verbose) {
        cerr << "Will start with the molecule furthest from the first molecule in the "
                "set\n";
      }
    } else if ("help" == s) {
      display_first_item_selection_options(cerr);
    } else {
      cerr << "Unrecognised first item selection directive '" << s << "'\n";
      display_first_item_selection_options(cerr);
    }
  }
  if (cl.option_present('b')) {
    if (!cl.value('b', blurr_distances) || blurr_distances < static_cast<float>(0.0)) {
      cerr << "The blurr distances option (-b) must be a non negative number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Distance blurring factor set to " << blurr_distances << '\n';
      std::random_device rd;
      std::uniform_real_distribution<similarity_type_t> u(0.0f, 1.0f);

      for (int i = 0; i < 5; i++) {
        similarity_type_t r = u(rd);
        cerr << "distance " << r << " becomes " << do_blurring(r, blurr_distances)
             << '\n';
      }
    }
  }

  if (cl.option_present('n')) {
    int n;
    if (!cl.value('n', n) || n < 1) {
      cerr << "the -n option must be followed by a whole positive number\n";
      usage(13);
    }

    if (n > pool_size) {
      cerr << "You asked for " << n << " molecules, but pool only contains " << pool_size
           << ". Shortened\n";
      n = pool_size;
    }

    number_to_select = n;
    if (verbose) {
      cerr << number_to_select << " molecules will be selected\n";
    }
  } else {
    number_to_select = pool_size;
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', stop_once_distance_drops_below) ||
        stop_once_distance_drops_below < 0.0f || stop_once_distance_drops_below >= 1.0f) {
      cerr << "The stop selection distance option (-t) must be a valid distance\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will stop selection once distance drops below "
           << stop_once_distance_drops_below << '\n';
    }
  }

  IWString_and_File_Descriptor output(1);

  (void)fpobj_spread(output);

  if (verbose) {
    cerr << "Nearest previously selected item distances between "
         << nearest_selected_neighbour_distance.minval() << " and "
         << nearest_selected_neighbour_distance.maxval();
    if (nearest_selected_neighbour_distance.n() > 1) {
      cerr << " ave " << nearest_selected_neighbour_distance.average();
    }
    cerr << '\n';
  }

  delete[] pool;  // leave this out for efficiency

  cerr << "Output can be processed with nplotnn\n";

  return 0;
}

int
main(int argc, char** argv)
{
  int rc = fpobj_spread(argc, argv);

  return rc;
}
