/*
  Near neighbour programme where the pool is huge, and must be
  read sequentially.
*/

#include <iostream>
#include <memory>

// GH Update code to fix the build under mac

#include "re2/re2.h"

#define USE_OMP
#ifdef USE_OMP
#include <omp.h>
#endif

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
// GH #define NEIGHBOUR_LIST_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Utilities/GFP_Tools/nearneighbours.pb.h"
#include "gfp_standard.h"
#include "nndata.h"

using std::cerr;
using std::endl;

// Must be defined before neighbour_list.h
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");
static IWString smiles_tag("$SMI<");

#define NEIGHBOUR_LIST_IMPLEMENTATION
#include "neighbour_list.h"
#include "tversky.h"

template class Neighbour_List<float, FP_and_Smiles, FP_and_Smiles>;

static int fingerprints_processed = 0;

static Tversky tversky_parameters;

static int min_number_neighbours = 0;

static int do_not_compare_molecules_with_themselves = 0;

static int sort_neighbour_list_at_end = 0;

static Report_Progress report_progress;

static IWDigits iwdigits(200);

static int strip_leading_zeros_from_identifiers = 0;

/*
  Distance thresholds
*/

static similarity_type_t upper_distance_threshold = 2.0;
static similarity_type_t lower_distance_threshold = -1.0;

static int equal_weight_tanimoto_distance_measure = 0;

/*
  Nov 2001, I'm getting NaN's occasionally. Don't yet know why, so let's have the ability
  to ignore these things
*/

static int keep_going_after_fatal_error = 0;

static int fatal_errors_encountered = 0;

static int write_neighbours_as_protos = 0;

class Needle : public Neighbour_List<similarity_type_t, Smiles_ID_Dist, GFP_Standard> {
 private:
  IWString _smiles;
  IWString _id;
  similarity_type_t _maxd;

  //  private functions

  void _consider(const similarity_type_t, const IWString& smiles, IWString& id);
  void _consider_min_nbrs(const similarity_type_t, const IWString& smiles, IWString& id);

 public:
  Needle();

  int build(const IW_TDT&);

  const IWString& smiles() const {
    return _smiles;
  }

  const IWString& id() const {
    return _id;
  }

  void consider(similarity_type_t, const IWString& smiles, IWString& id);

  // Create a nnbr::NearNeighbours proto from our data.
  nnbr::NearNeighbours MakeProto();

  // Write .nn data to `output`.
  int WriteNNData(IWString_and_File_Descriptor& output);
};

Needle::Needle() {
  _maxd = 1.0f;

  return;
}

int
Needle::build(const IW_TDT& tdt) {
  if (!tdt.dataitem_value(smiles_tag, _smiles) ||
      !tdt.dataitem_value(identifier_tag, _id)) {
    cerr << "Needle::build:cannot extract smiles and/or id from TDT\n";
    return 0;
  }

  if (strip_leading_zeros_from_identifiers) {
    _id.remove_leading_chars('0');
  }

  return 1;
}

void
Needle::consider(similarity_type_t d, const IWString& smiles,
                 IWString& id)  // not const
{
  if (static_cast<similarity_type_t>(0.0f) == d &&
      do_not_compare_molecules_with_themselves) {
    if (_id == id) {
      return;
    }

    if (strip_leading_zeros_from_identifiers && id.starts_with('0')) {
      id.remove_leading_chars('0');
      if (id == _id) {
        return;
      }
    }
  }

  if (min_number_neighbours > 0) {
    _consider_min_nbrs(d, smiles, id);
  } else {
    _consider(d, smiles, id);
  }

  return;
}

void
Needle::_consider(const similarity_type_t d, const IWString& smiles, IWString& id) {
  if (d < lower_distance_threshold) {
    return;
  }

  if (d > upper_distance_threshold) {
    return;
  }

  extra(smiles, id, d);
}

void
Needle::_consider_min_nbrs(const similarity_type_t d, const IWString& smiles,
                           IWString& id)  // not const
{
  assert(min_number_neighbours > 0);

  const int nn = number_neighbours();

  if (nn < min_number_neighbours)  // store them all until we get to the min
  {
    extra(smiles, id, d);
    return;
  }

  // Now we have the min number of nbrs, what do we do with this one?

  if (d < lower_distance_threshold) {
    return;
  }

  if (d > distance_of_furthest_neighbour()) {
    return;
  }

  if (0 == extra(smiles, id, d)) {
    return;
  }

  // If we have a min number neighbours, we must guard against having the list
  // of neighbours grow without bounds

  // cerr << "nbrs " << number_neighbours() << " max d " <<
  // distance_of_furthest_neighbour() << " max " << _upper_distance_threshold << endl;

  if (nn + 1 > min_number_neighbours &&
      distance_of_furthest_neighbour() > upper_distance_threshold) {
    this->shrink(min_number_neighbours);
  }

  return;
}

// Make a NearNeighbours proto based on the current state.
// Note that this interacts with several global variables, and
// the state may be changed as a result.
nnbr::NearNeighbours
Needle::MakeProto() {
  nnbr::NearNeighbours result;
  if (smiles_tag.length()) {
    result.set_smiles(_smiles.data(), _smiles.length());
  }
  result.set_name(_id.data(), _id.length());

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
Needle::WriteNNData(IWString_and_File_Descriptor& output) {
  if (smiles_tag.length()) {
    output << smiles_tag << smiles() << ">\n";
  }
  output << identifier_tag << id();

  output << ' ';
  iwdigits.append_number(output, number_neighbours());
  output << ">\n";

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
  Our pool is an array of fingerprints
  We maintain a separate array of neighbour list information
*/

static GFP_Standard* pool = nullptr;
static Needle* needle = nullptr;

static int pool_size = 0;

const char* prog_name = nullptr;

static int verbose = 0;

static const_IWSubstring output_identifier("PCN<");

static int
build_fingerprint(IW_TDT& tdt, GFP_Standard& fp, int check_tags = 0) {
  IW_General_Fingerprint gfp;

  int fatal;
  if (!gfp.construct_from_tdt(tdt, fatal)) {
    cerr << "Cannot read fingerprint\n";
    return 0;
  }

  if (check_tags) {
    if (!standard_fingerprints_present()) {
      return 0;
    }
  }

  fp.build_molecular_properties(gfp.molecular_properties_integer());
  fp.build_iw(gfp[0]);
  fp.build_mk(gfp[1]);
  fp.build_mk2(gfp[2]);

  return 1;
}

static int
build_pool(iwstring_data_source& input) {
  IW_TDT tdt;

  IWString tmp;

  int ndx = 0;

  for (; tdt.next(input), ndx < pool_size; ndx++) {
    if (!needle[ndx].build(tdt)) {
      return 0;
    }

    if (!build_fingerprint(tdt, pool[ndx], 0 == ndx)) {
      return 0;
    }
  }

  if (verbose) {
    cerr << ndx << " fingerprint objects added to pool\n";
  }

  pool_size = ndx;

  return 1;
}

static int
build_pool(const const_IWSubstring& fname) {
  IWString tmp(fname);

  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);

  if (0 == pool_size)  // must grep the file to find out how many
  {
    assert(nullptr == pool);

    pool_size = input.count_records_starting_with(identifier_tag);
    if (0 == pool_size) {
      cerr << "Yipes, cannot find any '" << identifier_tag << "' in the input\n";
      return 0;
    }

    if (verbose) {
      cerr << "Input contains " << pool_size << " fingerprints\n";
    }

    pool = new GFP_Standard[pool_size];
    needle = new Needle[pool_size];

    if (nullptr == pool) {
      cerr << "Yipes, cannot allocate space for " << pool_size << " fingerprints\n";
      return 0;
    }
  }

  return build_pool(input);
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

static void
nearneighbours(GFP_Standard& fp, similarity_type_t* d) {
#pragma omp parallel for schedule(dynamic, 256)
  for (int i = 0; i < pool_size; i++) {
    //  if (! can_be_compared(pool[i], fp))
    //    d[i] = -1.0f;
    //  else
    d[i] = static_cast<similarity_type_t>(1.0f) - pool[i].tanimoto(fp);
  }

  return;
}

static int
compute_neighbours_stored(const Needle* pool, int pool_size) {
  int rc = 0;
  for (int i = 0; i < pool_size; i++) {
    rc += needle[i].number_neighbours();
  }

  return rc;
}

//  It goes somewhere in the list

static int
nearneighbours(iwstring_data_source& input, IWString_and_File_Descriptor& output) {
  similarity_type_t* d = new similarity_type_t[pool_size];
  std::unique_ptr<similarity_type_t[]> free_d(d);

  IW_TDT tdt;
  while (tdt.next(input)) {
    IWString smiles, id;
    if (!tdt.dataitem_value(smiles_tag, smiles) ||
        !tdt.dataitem_value(identifier_tag, id)) {
      cerr << "Cannot extract smiles and/or id\n";
      return 0;
    }

    GFP_Standard fp;
    if (!build_fingerprint(tdt, fp)) {
      return 0;
    }

    fingerprints_processed++;

    nearneighbours(fp, d);

#pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < pool_size; ++i) {
      needle[i].consider(d[i], smiles, id);
    }

    if (report_progress()) {
      int n = compute_neighbours_stored(needle, pool_size);
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

  input.set_dos(1);

  if (fname.starts_with("F:")) {
    return nearneighbours_list_of_files(input, output);
  }

  if (verbose) {
    cerr << "Processing '" << fname << "'\n";
  }

  return nearneighbours(input, output);
}

static void
usage(int rc) {
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
  cerr << " -n <number>      specify how many neighbours to find\n";
  cerr << " -T <distance>    specify upper distance threshold\n";
  cerr << " -t <distance>    specify lower distance threshold\n";
  cerr << " -m <number>      the minimum number of neighbours to find\n";
  cerr << " -r <number>      report progress every <number> fingerprints\n";
  cerr << " -h               discard neighbours with zero distance and the same ID as the target\n";
  cerr << " -k               generate nnbr::NearNeighbours textproto output\n";
  // cerr << " -B <qualifier>   various other options, enter '-B help' for details\n";
  // cerr << " -V ...           Tversky specification, enter '-V help' for details\n";
  cerr << " -v               verbose output\n";
// clang-format on

  exit(rc);
}

static void
UpdateCounters(const Needle& needle, extending_resizable_array<int>& neighbours,
               Accumulator<float> closest_neighbour_distance) {
  const int nbrs = needle.number_neighbours();

  neighbours[nbrs]++;
  if (nbrs) {
    closest_neighbour_distance.extra(needle.distance_of_closest_neighbour());
  }
}

static int
nearneighbours(int argc, char** argv) {
  Command_Line cl(argc, argv, "vn:p:t:T:r:V:hB:N:m:zk");

  if (cl.unrecognised_options_encountered()) {
    if (cl.option_present('F') || cl.option_present('P')) {
      cerr << "NOte that traditional gfp options like -F and -P do not work with "
              "*standard* tools\n";
      return 1;
    }

    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  // initialise_string_distances();

  initialise_properties_ratios();

  set_default_iwstring_float_concatenation_precision(4);

  set_report_fingerprint_status(0);

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

  if (cl.option_present('T')) {
    int i = 0;
    const_IWSubstring t;
    while (cl.value('T', t, i++)) {
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

  if (cl.option_present('h')) {
    do_not_compare_molecules_with_themselves = 1;

    if (verbose) {
      cerr
          << "Will discard neighbours with zero distance and the same id as the target\n";
    }
  }

  if (cl.option_present('k')) {
    write_neighbours_as_protos = 1;
    if (verbose) {
      cerr << "Will generate output in proto form\n";
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
      } else if ("bignn" == b) {
        sort_neighbour_list_at_end = 1;
        if (verbose) {
          cerr << "Algorithm for large neighbour lists\n";
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
        cerr << " -B bignn      use algorithm optimum for large neighour lists\n";
        cerr << " -B ewt        distance metric is equal weight Tanimoto\n";
        return 0;
      } else {
        cerr << "Unrecognised -B qualifier '" << b << "'\n";
        usage(7);
      }
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

  set_include_newlines_in_tdt(0);

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
      needle[i].set_neighbours_to_find(neighbours_to_find);
    }

    if (verbose) {
      cerr << "A maximum of " << neighbours_to_find
           << " neighbours of each molecule will be found\n";
    }
  }

  if (sort_neighbour_list_at_end) {
    for (int i = 0; i < pool_size; i++) {
      needle[i].set_sort_neighbour_list_at_end(1);
    }
  }

  if (cl.option_present('z')) {
    strip_leading_zeros_from_identifiers = 1;

    if (verbose) {
      cerr << "Leading zero's stripped from identifiers\n";
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
      Needle& nni = needle[i];

      nnbr::NearNeighbours proto = nni.MakeProto();

      if (verbose) {
        UpdateCounters(nni, neighbours, closest_neighbour_distance);
      }
      gfp::WriteNNData(proto, output);
    }
  } else {
    for (int i = 0; i < pool_size; i++) {
      Needle& nni = needle[i];

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
