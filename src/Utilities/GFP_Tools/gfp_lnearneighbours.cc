/*
  Near neighbour programme where the pool is huge, and must be
  read sequentially.

  The intial implementation was with the Dow Agro database of over
  1 million molecules. We read our molecules of interest into memory,
  and then scan
*/

#include <iostream>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwqsort/iwqsort.h"

// #define FB_ENTROPY_WEIGHTED_FPS

#include "Utilities/GFP_Tools/nearneighbours.pb.h"
#include "gfp.h"
#include "neighbour_list.h"
#include "nndata.h"
#include "sparse_collection.h"
#include "tversky.h"

using std::cerr;
using std::endl;

static int fingerprints_processed = 0;

static Tversky tversky_parameters;

static IWString identifier_tag("PCN<");

static IWString distance_tag("DIST<");

static IWString smiles_tag("$SMI<");

static IWString number_neighbours_tag;

static bool max_distance_includes_max = false;

static int min_number_neighbours = 0;

static int do_not_compare_molecules_with_themselves = 0;

static int sort_neighbour_list_at_end = 0;

Set_of_Sparse_Fingerprint_Collection_Profile sfcp;

static Report_Progress report_progress;

static IWDigits iwdigits(200);

static int remove_leading_zeros_from_identifiers = 0;

static int chop_dos_carriage_returns = 0;

/*
  Distance thresholds
*/

static similarity_type_t upper_distance_threshold = 2.0;
static similarity_type_t lower_distance_threshold = -1.0;

/*
  We can also have a threshold associated with each molecule
*/

static IWString upper_distance_threshold_tag;

static int items_with_max_distance_threshold_specified = 0;

static int equal_weight_tanimoto_distance_measure = 0;

static int write_neighbours_as_protos = 0;

#ifdef PERHAPS_IMPLEMENT_SOMETIME
// This is the functionality of gfp_lnearneighbours_ndx which
// was deprecated.
// I think the behaviour can be done via nplotnn or nn2csv or...

// By default the PCN of each neighbour is its name.
// But we can also write it as its index in the haystack.
static int write_neighbours_as_indices = 0;
#endif

/*
  When searching against a sorted haystack, each fingerprint needs to
  have a status of not_yet_started, active or done
*/

#define NOT_YET_STARTED 0
#define SEARCHING_ACTIVE 1
#define SEARCHING_DONE 2

class NN_Object : public FP_and_Smiles,
                  public Neighbour_List<similarity_type_t, FP_and_Smiles, FP_and_Smiles> {
 private:
  similarity_type_t _upper_distance_threshold;

  int _search_status;

  // private functions

 public:
  NN_Object();

  void consider(FP_and_Smiles&);
  void consider(FP_and_Smiles& rhs, similarity_type_t must_be_larger_than);

  int construct_from_tdt(IW_TDT&, int&);

  int search_status() const {
    return _search_status;
  }

  void set_search_status(int s) {
    _search_status = s;
  }

  // Create a nnbr::NearNeighbours proto. Note that generation
  // is influenced by several global variables.
  nnbr::NearNeighbours MakeProto();

  // Write .nn to `output`.
  int WriteNNData(IWString_and_File_Descriptor& output);
};

NN_Object::NN_Object() {
  _upper_distance_threshold = upper_distance_threshold;

  _search_status = NOT_YET_STARTED;

  return;
}

static similarity_type_t zero_distance = static_cast<similarity_type_t>(0.0);

/*
  Nov 2001, I'm getting NaN's occasionally. Don't yet know why, so let's have the ability
  to ignore these things
*/

static int keep_going_after_fatal_error = 1;

static int fatal_errors_encountered = 0;

void
NN_Object::consider(FP_and_Smiles& rhs) {
  similarity_type_t d;

  // clang-format off
  //cerr << "equal_weight_tanimoto_distance_measure " << equal_weight_tanimoto_distance_measure << " values " << IW_General_Fingerprint::equal_weight_tanimoto(rhs) << ' ' << IW_General_Fingerprint::tanimoto(rhs) << " btw " << rhs.id() << " and " << id() << endl;
  // clang-format on

  if (equal_weight_tanimoto_distance_measure) {
    d = static_cast<similarity_type_t>(1.0) -
        IW_General_Fingerprint::equal_weight_tanimoto(rhs);
  } else if (tversky_parameters.active()) {
    d = static_cast<similarity_type_t>(1.0) -
        IW_General_Fingerprint::tversky(rhs, tversky_parameters);
  } else {
    d = static_cast<similarity_type_t>(1.0) - IW_General_Fingerprint::tanimoto(rhs);
  }

  if (d < static_cast<similarity_type_t>(0.0) ||
      d > static_cast<similarity_type_t>(1.0)) {
    cerr << "Bad news, invalid distance " << d << " between '" << _id << "' and '"
         << rhs.id() << "'\n";
    if (keep_going_after_fatal_error) {
      fatal_errors_encountered++;
      return;
    }

    abort();
  }

  bool inRange;
  if (max_distance_includes_max) {
    inRange = d > lower_distance_threshold && d <= _upper_distance_threshold;
  } else {
    inRange = d > lower_distance_threshold && d < _upper_distance_threshold;
  }

  if (inRange) {
    ;
  } else if (min_number_neighbours > 0) {
    if (number_neighbours() < min_number_neighbours) {
      ;
    } else if (d < distance_of_furthest_neighbour()) {
      ;
    } else {
      return;
    }
  } else {
    return;
  }

  if (do_not_compare_molecules_with_themselves && zero_distance == d && _id == rhs.id()) {
    return;
  }

  extra(rhs, d);

  // If we have a min number neighbours, we must guard against having the list
  // of neighbours grow without bounds

  // clang-format off
  //cerr << "nbrs " << number_neighbours() << " max d " << distance_of_furthest_neighbour() << " max " << _upper_distance_threshold << endl;
  // clang-format on

  if (min_number_neighbours > 0 && number_neighbours() > min_number_neighbours) {
    bool outOfRange;
    if (max_distance_includes_max) {
      outOfRange = distance_of_furthest_neighbour() > _upper_distance_threshold;
    } else {
      outOfRange = distance_of_furthest_neighbour() >= _upper_distance_threshold;
    }

    if (outOfRange) {
      Neighbour_List<similarity_type_t, FP_and_Smiles, FP_and_Smiles>::shrink(
          min_number_neighbours);
    }
  }
  return;
}

void
NN_Object::consider(FP_and_Smiles& rhs, similarity_type_t must_be_larger_than) {
  similarity_type_t d;

  if (!IW_General_Fingerprint::tanimoto(rhs, must_be_larger_than, d)) {
    return;
  }

  if (d < static_cast<similarity_type_t>(0.0) ||
      d > static_cast<similarity_type_t>(1.0)) {
    cerr << "Bad news, invalid distance " << d << " between '" << _id << "' and '"
         << rhs.id() << "'\n";
    if (keep_going_after_fatal_error) {
      fatal_errors_encountered++;
      return;
    }

    abort();
  }

  d = static_cast<similarity_type_t>(1.0) - d;

  bool inRange;
  if (max_distance_includes_max) {
    inRange = d > lower_distance_threshold && d <= _upper_distance_threshold;
  } else {
    inRange = d > lower_distance_threshold && d < _upper_distance_threshold;
  }

  if (inRange) {
    ;
  } else if (min_number_neighbours > 0) {
    if (number_neighbours() < min_number_neighbours) {
      ;
    } else if (d < distance_of_furthest_neighbour()) {
      ;
    } else {
      return;
    }
  } else {
    return;
  }

  // if (d <= lower_distance_threshold || d >= _upper_distance_threshold)
  //   return;

  if (do_not_compare_molecules_with_themselves && zero_distance == d && _id == rhs.id()) {
    return;
  }

  extra(rhs, d);

  return;
}

int
NN_Object::construct_from_tdt(IW_TDT& tdt, int& fatal) {
  if (!FP_and_Smiles::construct_from_tdt(tdt, fatal)) {
    return 0;
  }

  if (0 == upper_distance_threshold_tag.length()) {
    return 1;
  }

  const_IWSubstring u;
  if (!tdt.dataitem(upper_distance_threshold_tag, u)) {  // ok to be unspecified
    return 1;
  }

  u.remove_leading_chars(upper_distance_threshold_tag.length());
  if (u.starts_with('<')) {
    u.remove_leading_chars(1);
  }
  u.chop();

  if (!u.numeric_value(_upper_distance_threshold) || _upper_distance_threshold < 0.0 ||
      _upper_distance_threshold > 1.0) {
    cerr << "Invalid upper distance threshold value '" << u << "'\n";
    return 0;
  }

  items_with_max_distance_threshold_specified++;

  return 1;
}

nnbr::NearNeighbours
NN_Object::MakeProto() {
  nnbr::NearNeighbours result;

  result.set_name(_id.data(), _id.length());
  if (!smiles_tag.empty()) {
    result.set_smiles(_smiles.data(), _smiles.length());
  }

  if (sort_neighbour_list_at_end) {
    sort_neighbour_list();
  }

  if (min_number_neighbours > 0 && upper_distance_threshold <= 1.0) {
    remove_distant_neighbours(min_number_neighbours, upper_distance_threshold);
  }

  for (const Smiles_ID_Dist* sid : _neighbours) {
    nnbr::Nbr* nbr = result.add_nbr();

    const IWString& id = sid->id();
    nbr->set_id(id.data(), id.length());
    if (!smiles_tag.empty()) {
      const IWString& smi = sid->smiles();
      nbr->set_smi(smi.data(), smi.length());
    }
    nbr->set_dist(sid->distance());
  }

  return result;
}

int
NN_Object::WriteNNData(IWString_and_File_Descriptor& output) {
  if (smiles_tag.length()) {
    output << smiles_tag << _smiles << ">\n";
  }
  output << identifier_tag << _id;

  if (0 == number_neighbours_tag.length()) {
    output << ' ';
    iwdigits.append_number(output, number_neighbours());
    output << ">\n";
  } else {
    output << ">\n" << number_neighbours_tag;
    iwdigits.append_number(output, number_neighbours());
    output << ">\n";
  }

  if (sort_neighbour_list_at_end) {
    sort_neighbour_list();
  }

  if (min_number_neighbours > 0 && upper_distance_threshold <= 1.0) {
    remove_distant_neighbours(min_number_neighbours, upper_distance_threshold);
  }

  write(output);
  output << "|\n";

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

/*
  Our pool is an array of NN_Object objects
*/

static NN_Object* pool = nullptr;

static int pool_size = 0;

static int pool_is_sorted = 0;
static int number_active_fingerprints = 0;
static similarity_type_t must_be_larger_than = static_cast<float>(0.0);

static int similarity_computations_performed = 0;

const char* prog_name = nullptr;

static int verbose = 0;

/*
  During output we need to know which identifier is being used
*/

static const_IWSubstring output_identifier("PCN<");

/*
  Since we will want to go back and get members of the pool, we must make the
  file from which they come available throughout programme execution.
*/

static iwstring_data_source pool_data_source;

static int
build_pool(iwstring_data_source& input) {
  assert(pool_size > 0);

  int i = 0;

  IW_TDT tdt;
  while (tdt.next(input)) {
    int fatal;
    if (!pool[i].construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      continue;
    }

#ifdef FB_ENTROPY_WEIGHTED_FPS
    if (!pool[i].convert_01_fingerprint_to_integer_molecular_properties()) {
      return 0;
    }
#endif

    i++;

    if (i >= pool_size) {
      cerr << "Pool is full, max " << pool_size << endl;
      break;
    }
  }

  if (verbose) {
    cerr << i << " fingerprint objects added to pool\n";
  }

  pool_size = i;

  if (remove_leading_zeros_from_identifiers) {
    for (int i = 0; i < pool_size; i++) {
      NN_Object& pi = pool[i];

      if (!pi.id().starts_with('0')) {
        continue;
      }

      const_IWSubstring tmp(pi.id());
      tmp.remove_leading_chars('0');
      pi.set_id(tmp);
    }
  }

  return 1;
}

static int
build_pool(const const_IWSubstring& fname) {
  IWString tmp(fname);

  if (!pool_data_source.open(tmp))  // method is non-const on its argument!
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (chop_dos_carriage_returns) {
    pool_data_source.set_dos(1);
  }

  if (0 == pool_size)  // must grep the file to find out how many
  {
    assert(nullptr == pool);

    pool_size = pool_data_source.count_records_starting_with(identifier_tag);
    if (0 == pool_size) {
      cerr << "Yipes, cannot find any '" << identifier_tag << "' in the input\n";
      return 0;
    }

    if (verbose) {
      cerr << "Input contains " << pool_size << " fingerprints\n";
    }

    pool = new NN_Object[pool_size];

    if (nullptr == pool) {
      cerr << "Yipes, cannot allocate space for " << pool_size << " fingerprints\n";
      return 0;
    }

    number_active_fingerprints = pool_size;
  }

  return build_pool(pool_data_source);
}

template void resizable_array_base<Smiles_ID_Dist*>::iwqsort<SID_Comparator>(
    SID_Comparator&);
template void iwqsort<Smiles_ID_Dist*, SID_Comparator>(Smiles_ID_Dist**, int,
                                                       SID_Comparator&);
template void iwqsort<Smiles_ID_Dist*, SID_Comparator>(Smiles_ID_Dist**, int,
                                                       SID_Comparator&, void*);
template void compare_two_items<Smiles_ID_Dist*, SID_Comparator>(Smiles_ID_Dist**,
                                                                 SID_Comparator&, void*);
template void swap_elements<Smiles_ID_Dist*>(Smiles_ID_Dist*&, Smiles_ID_Dist*&, void*);
template void move_in_from_right<Smiles_ID_Dist*, SID_Comparator>(Smiles_ID_Dist**, int&,
                                                                  int&, SID_Comparator&);
template void move_in_from_left<Smiles_ID_Dist*, SID_Comparator>(Smiles_ID_Dist**, int&,
                                                                 int&, int,
                                                                 SID_Comparator&, void*);

// GH template class Neighbour_List<similarity_type_t, FP_and_Smiles, FP_and_Smiles>;

static int previous_natoms = 0;

static void
nearneighbours_sorted_pool(FP_and_Smiles& fp) {
  int n = fp.natoms();

  if (n == previous_natoms) {
    ;
  } else if (n > previous_natoms) {
    previous_natoms = n;
  } else {
    cerr << "Supposedly sorted haystack not sorted, previous " << previous_natoms
         << " now " << n << " atoms, " << fingerprints_processed
         << " fingerprints read, '" << fp.id() << "'\n";
    abort();
  }

  for (int i = 0; i < pool_size; i++) {
    NN_Object& pi = pool[i];

    if (SEARCHING_DONE == pi.search_status()) {
      continue;
    }

    if (can_be_compared(pi, fp)) {
      if (verbose > 1 && SEARCHING_ACTIVE != pi.search_status()) {
        cerr << "After " << fingerprints_processed << " fingerprints, pool item " << i
             << " (" << pi.natoms() << ") activated (" << fp.natoms() << ") "
             << number_active_fingerprints << " active fingerprints\n";
      }

      pi.set_search_status(SEARCHING_ACTIVE);

      pool[i].consider(fp, must_be_larger_than);

      similarity_computations_performed++;
    } else if (SEARCHING_ACTIVE == pi.search_status())  // was active, now turn off
    {
      pi.set_search_status(SEARCHING_DONE);
      number_active_fingerprints--;
      if (verbose > 1) {
        cerr << "After " << fingerprints_processed << " fingerprints, pool item " << i
             << " (" << pi.natoms() << ") turned off, (" << fp.natoms() << ") "
             << number_active_fingerprints << " remaining\n";
      }
      if (0 == number_active_fingerprints) {
        return;
      }
    }
  }

  return;
}

static void
nearneighbours(FP_and_Smiles& fp) {
  for (int i = 0; i < pool_size; i++) {
    if (!can_be_compared(pool[i], fp)) {
      continue;
    }

    pool[i].consider(fp);
  }

  return;
}

static int
compute_neighbours_stored(const NN_Object* pool, int pool_size) {
  int rc = 0;
  for (int i = 0; i < pool_size; i++) {
    rc += pool[i].number_neighbours();
  }

  return rc;
}

static void
do_remove_leading_zeros_from_identifiers(FP_and_Smiles& fp) {
  if (!fp.id().starts_with('0')) {
    return;
  }

  const_IWSubstring tmp(fp.id());
  tmp.remove_leading_chars('0');

  fp.set_id(tmp);

  return;
}

//  It goes somewhere in the list

static int
nearneighbours(iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  IW_TDT tdt;
  while (tdt.next(input) && number_active_fingerprints > 0) {
    FP_and_Smiles fp;
    int fatal;

    if (!fp.construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        return 0;
      }

      continue;
    }

#ifdef FB_ENTROPY_WEIGHTED_FPS
    if (!fp.convert_01_fingerprint_to_integer_molecular_properties()) {
      return 0;
    }
#endif

    if (remove_leading_zeros_from_identifiers) {
      do_remove_leading_zeros_from_identifiers(fp);
    }

    fingerprints_processed++;

    if (sfcp.active()) {
      fp.convert_to_non_sparse_forms(sfcp);
    }

    if (pool_is_sorted) {
      nearneighbours_sorted_pool(fp);
    } else {
      nearneighbours(fp);
    }

    if (report_progress()) {
      int n = compute_neighbours_stored(pool, pool_size);
      cerr << "Processed " << fingerprints_processed << " fingerprints, storing " << n
           << " neighbours (approx " << ((n * 100) / 1000000) << " MB)\n";
    }
  }

  return 1;
}

static int nearneighbours(const const_IWSubstring& fname,
                          IWString_and_File_Descriptor& output);

/*
  The input is a file containing a list of files to process
*/

static int
nearneighbours_list_of_files(iwstring_data_source& input,
                             IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    buffer.strip_leading_blanks();
    buffer.strip_trailing_blanks();

    if (0 == buffer.length()) {
      continue;
    }

    if (buffer.starts_with('#')) {
      continue;
    }

    if (!nearneighbours(buffer, output)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
nearneighbours(const const_IWSubstring& fname, IWString_and_File_Descriptor& output) {
  const_IWSubstring myfname = fname;

  if (fname.starts_with("F:")) {
    myfname.remove_leading_chars(2);
  }

  iwstring_data_source input(myfname);
  if (!input.ok()) {
    cerr << "Cannot open input '" << myfname << "'\n";
    return 0;
  }

  if (chop_dos_carriage_returns) {
    input.set_dos(1);
  }

  if (fname.starts_with("F:")) {
    return nearneighbours_list_of_files(input, output);
  }

  if (verbose) {
    cerr << "Processing '" << fname << "'\n";
  }

  return nearneighbours(input, output);
}

// clang-format off
static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << "Finds near neighbours\n";
  cerr << "Usage " << prog_name << " ... -p <needles> <haystack>\n";
  cerr << " -p <file>        specify file against which input is to be compared (needles)\n";
  cerr << " -s <number>      specify max pool size\n";
  cerr << " -n <number>      specify how many neighbours to find\n";
  cerr << " -T <distance>    specify upper distance threshold\n";
  cerr << " -m <number>      the minimum number of neighbours to find\n";
  cerr << " -T TAG=tag       upper distance threshold for each molecule in TAG<>\n";
  cerr << " -t <distance>    specify lower distance threshold\n";
  cerr << " -i <dataitem>    specify identifier dataitem in pool\n";
  cerr << " -I <dataitem>    specify identifier dataitem in input file\n";
  cerr << " -e <dataitem>    specify pool object dataitems to be echo'd (default $SMI and PCN)\n";
  cerr << " -r <number>      report progress every <number> fingerprints\n";
  cerr << " -D <tag>         tag for distances (default '" << distance_tag << "')\n";
  cerr << " -h               discard neighbours with zero distance and the same ID as the target\n";
  cerr << " -B <qualifier>   various other options, enter '-B help' for details\n";
  cerr << " -F ...           gfp options, enter '-F help' for details\n";
  cerr << " -V ...           Tversky specification, enter '-V help' for details\n";
  cerr << " -g <dist>        Abandon distance computation if any component > dist \n";
  cerr << " -k               generate nnbr::NearNeighbours textproto output\n";
  cerr << " -k               generate nnbr::NearNeighbours textproto output\n";
  cerr << " -v               verbose output\n";
// clang-format on

  exit (rc);
}

// clang-format on

static void
UpdateCounters(const NN_Object& needle, extending_resizable_array<int>& neighbours,
               Accumulator<float> closest_neighbour_distance) {
  const int nbrs = needle.number_neighbours();

  neighbours[nbrs]++;
  if (nbrs) {
    closest_neighbour_distance.extra(needle.distance_of_closest_neighbour());
  }
}

static int
nearneighbours(int argc, char** argv) {
  Command_Line cl(argc, argv, "vs:n:p:I:i:e:t:T:r:F:P:W:Q:V:hB:K:D:N:m:g:E:dk");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('d')) {
    chop_dos_carriage_returns = 1;

    if (verbose) {
      cerr << "Will chop DOS carriage returns\n";
    }
  }

  // initialise_string_distances();

  if (need_to_call_initialise_fingerprints(cl)) {
    if (!initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise GFP options\n";
      usage(23);
    }
  }

  if (cl.option_present('V')) {
    if (!tversky_parameters.parse_command_line(cl, 'V', verbose)) {
      cerr << "Cannot initialise Tversky conditions\n";
      usage(7);
    }
  }

  if (cl.option_present('t')) {
    similarity_type_t t;
    if (!cl.value('t', t) || t < 0.0 || t > 1.0) {
      cerr << "The -t option must be followed by a valid similarity value\n";
      usage(12);
    }

    lower_distance_threshold = t;

    if (verbose) {
      cerr << "Lower distance threshold set to " << lower_distance_threshold << endl;
    }
  }

  // And make sure we process the -T option before the pool is built in case
  // UPPER_DISTANCE_THRESHOLD_TAG is set. Also _upper_distance_threshold is set in the
  // constructor

  if (cl.option_present('T')) {
    int i = 0;
    const_IWSubstring t;
    while (cl.value('T', t, i++)) {
      if (t.starts_with("TAG=")) {
        upper_distance_threshold_tag = t;
        upper_distance_threshold_tag.remove_leading_chars(4);

        if (verbose) {
          cerr << "The upper distance threshold for each molecule in the '"
               << upper_distance_threshold_tag << "' tag\n";
        }

        continue;
      }

      if (!t.numeric_value(upper_distance_threshold) || upper_distance_threshold < 0.0 ||
          upper_distance_threshold > 1.0) {
        cerr << "The -T option must be followed by a valid similarity value\n";
        usage(12);
      }

      if (verbose) {
        cerr << "Upper distance threshold set to " << upper_distance_threshold << endl;
      }
    }
  }

  assert(lower_distance_threshold <= upper_distance_threshold);

#ifdef FB_ENTROPY_WEIGHTED_FPS
  if (!cl.option_present('E')) {
    cerr << "This version only works with entropy weighted fingerprints, must use -E "
            "option\n";
    usage(3);
  } else {
    if (!initialise_entropy_weighting_for_fingerprints(cl, 'E', verbose)) {
      cerr << "Cannot initialise entropy weighted fingerprint (-E)\n";
      return 0;
    }

    if (verbose) {
      cerr << "Entropy weighted fingerprints\n";
    }
  }

#endif

  if (cl.option_present('s')) {
    int nfp;
    if (!cl.value('s', nfp) || nfp < 1) {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }

    pool = new NN_Object[nfp];
    if (nullptr == pool) {
      cerr << "Yipes, could not allocate pool of size " << nfp << endl;
      return 62;
    }

    pool_size = nfp;
    number_active_fingerprints = pool_size;

    if (verbose) {
      cerr << "system sized to " << nfp << endl;
    }
  }

  // We need to be careful with the -i and -I options. Remember
  // that the pool is built first

  if (cl.option_present('i')) {
    (void)cl.value('i', identifier_tag);

    if (!identifier_tag.ends_with('<')) {
      identifier_tag += '<';
    }

    set_identifier_tag(identifier_tag);

    if (verbose) {
      cerr << "Identifiers in pool tagged as '" << identifier_tag << "'\n";
    }
  }

  if (cl.option_present('h')) {
    do_not_compare_molecules_with_themselves = 1;

    if (verbose) {
      cerr
          << "Will discard neighbours with zero distance and the same id as the target\n";
    }
  }

  if (cl.option_present('B')) {
    const_IWSubstring b;
    int i = 0;
    while (cl.value('B', b, i++)) {
      if ("nofatal" == b) {
        keep_going_after_fatal_error = 1;

        if (verbose) {
          cerr << "Will ignore bad distances!!!\n";
        }
      } else if ("nosmiles" == b) {
        smiles_tag.resize(0);

        if (verbose) {
          cerr << "Will discard all smiles\n";
        }
      } else if ("maxdistanceinclusive" == b) {
        max_distance_includes_max = true;

        if (verbose) {
          cerr << "Upper distance threshold is inclusive of the distance\n";
        }
      } else if (b.starts_with("nbrtag=")) {
        number_neighbours_tag = b;
        number_neighbours_tag.remove_leading_chars(7);

        if (verbose) {
          cerr << "Neighbour count written with tag '" << number_neighbours_tag << "'\n";
        }

        if (!number_neighbours_tag.ends_with('<')) {
          number_neighbours_tag += '<';
        }
      } else if ("bignn" == b) {
        sort_neighbour_list_at_end = 1;
        if (verbose) {
          cerr << "Algorithm for large neighbour lists\n";
        }
      } else if ("rmzero" == b) {
        remove_leading_zeros_from_identifiers = 1;
        if (verbose) {
          cerr << "Leading 0's removed from all identifiers\n";
        }
      } else if ("ewt" == b) {
        equal_weight_tanimoto_distance_measure = 1;
        if (verbose) {
          cerr << "Distance metric is equal weight Tanimoto\n";
        }
      } else if ("help" == b) {
        cerr << "The following -B qualifiers are recognised\n";

        cerr << " -B nofatal    ignore otherwise fatal errors\n";
        cerr << " -B nosmiles   discard neighbour smiles - things run faster and consume "
                "less memory\n";
        cerr << " -B nbrtag=<tag> write the number of neighbours for each target in "
                "<tag>\n";
        cerr << " -B bignn      use algorithm optimum for large neighour lists\n";
        cerr << " -B rmzero     remove leading zero's from identifiers\n";
        cerr << " -B ewt        distance metric is equal weight Tanimoto\n";
        cerr << " -B maxdistanceinclusive    distnaces up to and including the max "
                "distance are considered matches\n";
        return 0;
      } else {
        cerr << "Unrecognised -B qualifier '" << b << "'\n";
        usage(7);
      }
    }
  }

  if (cl.option_present('D')) {
    cl.value('D', distance_tag);

    if (verbose) {
      cerr << "Distances written with tag " << distance_tag << "'\n";
    }

    if (!distance_tag.ends_with('<')) {
      distance_tag.add('<');
    }
  }

  if (!cl.option_present('p')) {
    cerr << "Must specify a file of \"needles\" via the -p option\n";
    usage(5);
  }

  if (1 != cl.option_count('p')) {
    cerr << "Only one file of needles (-p) may be specified\n";
    usage(6);
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  if (cl.option_present('g')) {
    if (!cl.value('g', must_be_larger_than) || must_be_larger_than < 0.0 ||
        must_be_larger_than >= 1.0) {
      cerr << "The abandon distance computation option (-g) value must be a valid "
              "distance\n";
      usage(4);
    }
    pool_is_sorted = 1;

    if (verbose) {
      cerr << "Assuming pool sorted, each similarity component must be at least "
           << must_be_larger_than << endl;
    }

    must_be_larger_than = 1.0 - must_be_larger_than;  // convert to similarity
  }

  set_include_newlines_in_tdt(0);

  if (cl.option_present('K')) {
    if (!parse_sparse_to_dense_fingerprint_specifications(cl, 'K', verbose)) {
      cerr << "Cannot parse sparse to dense fingerprint specifications (-K)\n";
      usage(5);
    }
  }

  if (cl.option_present('k')) {
    write_neighbours_as_protos = 1;
    if (verbose) {
      cerr << "Will generate output in proto form\n";
    }
  }

  if (cl.option_present('p')) {
    const_IWSubstring fname = cl.string_value('p');

    if (!build_pool(fname)) {
      cerr << "Cannot build pool from '" << fname << "'\n";
      return 76;
    }

    if (0 == pool_size) {
      cerr << "Pool is empty, cannot continue\n";
      return 9;
    }

    if (verbose && upper_distance_threshold_tag.length()) {
      cerr << items_with_max_distance_threshold_specified << " of " << pool_size
           << " fingerprints had upper distance threshold individually specified\n";
    }
  }

  if (cl.option_present('K')) {
    for (int i = 0; i < pool_size; i++) {
      sfcp.build_profile(pool[i]);
    }

    cerr << "Profiling done\n";
    sfcp.report(cerr);
    sfcp.finished_profiling(verbose);

    for (int i = 0; i < pool_size; i++) {
      pool[i].convert_to_non_sparse_forms(sfcp);
    }
  }

  if (keep_going_after_fatal_error) {
    for (int i = 0; i < pool_size; i++) {
      pool[i].set_keep_going_after_fatal_error(1);
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

  // If they didn't ask for anything to be echo'd, just use defaults

  if (cl.option_present('e')) {
    cl.value('e', output_identifier);
    if (verbose) {
      cerr << "Identifier dataitem in output is '" << output_identifier << "'\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "The -r option must be followed by a positive whole number\n";
      usage(18);
    }
  }

  int neighbours_to_find = -1;
  if (cl.option_present('n')) {
    if (!cl.value('n', neighbours_to_find) || neighbours_to_find < 1) {
      cerr << "the -n option must be followed by a whole positive number\n";
      usage(13);
    }
  }

  if (neighbours_to_find < 0 &&
      !cl.option_present(
          'T'))  // otherwise it would keep the whole haystack as neighbours
  {
    neighbours_to_find = 1;

    if (verbose) {
      cerr << "By default, will find just one neighbour\n";
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', min_number_neighbours) || min_number_neighbours < 1) {
      cerr
          << "The minimum number of neighbours to find (-m) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose) {
      cerr << "Will find a minimum of " << min_number_neighbours
           << " neighbours for each needle\n";
    }

    if (neighbours_to_find > 0 && min_number_neighbours > neighbours_to_find) {
      cerr << "Inconsistent specification of neighbours_to_find " << neighbours_to_find
           << " and min_number_neighbours " << min_number_neighbours << endl;
      usage(4);
    }
  }

  if (neighbours_to_find > 0) {
    for (int i = 0; i < pool_size; i++) {
      pool[i].set_neighbours_to_find(neighbours_to_find);
    }

    if (verbose) {
      cerr << "A maximum of " << neighbours_to_find
           << " neighbours of each molecule will be found\n";
    }
  }

  if (sort_neighbour_list_at_end) {
    for (int i = 0; i < pool_size; i++) {
      pool[i].set_sort_neighbour_list_at_end(1);
    }
  }

  IWString_and_File_Descriptor output(1);

  for (int i = 0; i < cl.number_elements(); i++) {
    (void)nearneighbours(cl[i], output);
  }

  if (verbose) {
    cerr << fingerprints_processed << " fingerprints processed\n";
  }

  if (0 == fingerprints_processed) {
    cerr << "No fingerprints processed, no output\n";
    return 4;
  }

  if (fatal_errors_encountered) {
    cerr << "WARNING!, " << fatal_errors_encountered
         << " should have been fatal errors encountered\n";
  }

  extending_resizable_array<int> neighbours;
  Accumulator<similarity_type_t> closest_neighbour_distance;

  if (write_neighbours_as_protos) {
    for (int i = 0; i < pool_size; i++) {
      NN_Object& nni = pool[i];
      nnbr::NearNeighbours proto = nni.MakeProto();
      if (verbose) {
        UpdateCounters(nni, neighbours, closest_neighbour_distance);
      }
      gfp::WriteNNData(proto, output);
    }

  } else {
    set_default_iwstring_float_concatenation_precision(3);

    for (int i = 0; i < pool_size; i++) {
      NN_Object& nni = pool[i];

      nni.WriteNNData(output);
      if (verbose) {
        UpdateCounters(nni, neighbours, closest_neighbour_distance);
      }
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Closest neighbours between " << closest_neighbour_distance.minval()
         << " and " << closest_neighbour_distance.maxval();
    if (closest_neighbour_distance.n() > 1) {
      cerr << ", average " << closest_neighbour_distance.average();
    }
    cerr << '\n';

    for (int i = 0; i < neighbours.number_elements(); i++) {
      if (neighbours[i]) {
        cerr << neighbours[i] << " molecules had " << i << " neighbours\n";
      }
    }

    if (pool_is_sorted) {
      cerr << similarity_computations_performed
           << " similarity computations done on sorted pool\n";
    }
  }

  return 0;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = nearneighbours(argc, argv);

#ifdef USE_IWMALLOC
  terse_malloc_status(stderr);
#endif

  return rc;
}
