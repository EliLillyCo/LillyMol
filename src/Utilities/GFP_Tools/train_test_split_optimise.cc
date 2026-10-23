// Use the results of a nearest neighbour calculation to optimize
// the distances across a train/test split one or more train/test splits.

#include <stdlib.h>
#include <omp.h>

#include <algorithm>
#include <ctime>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "absl/container/flat_hash_map.h"

#define REPORT_PROGRESS_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"

#include "Utilities/GFP_Tools/train_test_split_activity.h"

#ifdef BUILD_BAZEL
#include "Utilities/GFP_Tools/nearneighbours.pb.h"
#else
#include "nearneighbours.pb.h"
#endif

namespace train_test_split {

using std::cerr;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << R"(Uses TFDataRecord nn data to generate optimized train/test splits
Randomly assigns train and test. Then iteratively swaps pairs if the sawp increases the inter-set distance.
Reads a set of precomputed nearest neighbours, computed to a fixed distance.
Make this as large as can reasonably be computed, but do not include unrealistic distances.
gfp_nearneighbours_single_file_tbb -T 0.35 -h 16 -S file.nn.tfdata file.gfp
train_test_split_opt -f 0.85 -n 1 -S OPT -o 2000000 -r 20000 -x 20000 -h 8 file.nn.tfdata.

 -f <train_fraction> fraction of data in the train split
 -n <nsplit>         number of splits needed.
 -S <stem>           write splits to <stem>, <stem>R<n> and <stem>E<n>.
 -c                  optimise for closness between train and test - the easiest possible split.
 -X                  also write cross split summary stats - expensive to compute.
 -C <fname>          weighted samples. A mapping from ID to number of instances.
 -A <fname>          activity file with header and id/activity columns. Repeat for multiple profiles.
 -A buckets=<n>      set default activity bucket count for subsequent -A files.
 -A <fname>:buckets=<n>[:quantile|:width] per-file activity bucketisation.
 -Y <weight>         weight for activity distribution optimisation; required with -A.
 -o <nopt>           number of optimisation steps to try per split.
 -s <sec>            run each split for <sec> seconds.
 -r <n>              report progress every <n> steps.
 -x <n>              abandon optimisation of <n> steps since last switch accepted.
 -h <n>              number of OMP threads to use during exact calculation.
 -Z ...              miscellaneous options, enter '-Z help' for info.
 -v                  verbose output.
)";
// clang-format on
  
  ::exit(rc);
}

void
DisplayMiscOptions(int rc) {
  cerr << R"( -Z seed=<n>        seed the random number generator for reproducible splits.
)";

  ::exit(rc);
}

// Neighbour information consists of an ID and an integer distance.
struct Nbr {
  // This will be the index of this neighbour into the overall array of needles,
  // which is external to this class.
  uint32_t _id;
  // This value could be narrower, even uint8_t would be fine.
  uint32_t _dist;
};

// Each needle knows its current state, and info about its neighbours.
// Deliberately minimal information
class Needle {
  private:
    // The needle does not need to know its id since that is the same as its position
    // in the overall array of needles.

    // Whether or not this item is in the training set.
    int _train;

    // the neighbours.
    uint32_t _number_nbrs;
    Nbr* _nbrs;

    uint32_t _count;

  public:
    Needle();
    ~Needle();

    int Build(const nnbr::NearNeighbours& proto,
              const std::unordered_map<std::string, uint32_t>& id_to_ndx,
              float upper_distance_limit);

    uint32_t number_neighbours() {
      return _number_nbrs;
    }

    int in_train() const {
      return _train;
    }
    void set_train(int s) {
      _train = s;
    }

    void invert_train() {
      _train = !_train;
    }

    uint32_t count() const {
      return _count;
    }
    void set_count(uint32_t s) {
      _count = s;
    }

    uint32_t Neighbour(uint32_t ndx) const {
      return _nbrs[ndx]._id;
    }

    uint32_t NeighbourDistance(uint32_t ndx) const {
      return _nbrs[ndx]._dist;
    }

    void SetNbrDistances(uint32_t* nbrdist, uint32_t max_distance) const;
    void SetNbrDistances(uint32_t* nbrdist) const;

    void WriteNbrs(std::ostream& output) const;

    // The distance of the nearest neighbour. If no neighbours, 101.0 is returned.
    // Note that the discretised distance is returned, the caller will need to convert
    // to a distance value in [0,1].
    float ClosestDistance() const;
};

Needle::Needle() {
  _train = 0;
  _number_nbrs = 0;
  _nbrs = nullptr;
  _count = 1;
}

Needle::~Needle() {
  if (_nbrs != nullptr) {
    delete [] _nbrs;
  }
}

uint32_t
DistanceToInt(float d) {
  if (d == 0.0f) {
    return 0;
  }

  return static_cast<uint32_t>(d * 100.0f) + 1;
}

// Since neighbours store their id's as indices, we need a mapping
// from names to indices `id_to_ndx`.
// `upper_distance_limit` does not work, not sure why, it breaks
// the optimisation.
int
Needle::Build(const nnbr::NearNeighbours& proto,
              const std::unordered_map<std::string, uint32_t>& id_to_ndx,
              float upper_distance_limit) {
  if (proto.nbr_size() == 0) {
    return 1;
  }

  int ok_nbrs = 0;
  if (upper_distance_limit == std::numeric_limits<float>::max()) {
    ok_nbrs = proto.nbr_size();
  } else {
    for (const nnbr::Nbr& nbr : proto.nbr()) {
      if (nbr.dist() <= upper_distance_limit) {
        ++ok_nbrs;
      }
    }

    if (ok_nbrs == 0) {
      return 1;
    }
  }

  _number_nbrs = ok_nbrs;
  _nbrs = new Nbr[_number_nbrs];

  int ndx = 0;
  for (const nnbr::Nbr& nbr : proto.nbr()) {
    const auto iter = id_to_ndx.find(nbr.id());
    if (iter == id_to_ndx.end()) {
      cerr << "Needle::Build:no ndx for '" << nbr.id() << "'\n";
      return 0;
    }

    _nbrs[ndx]._id = iter->second;
    _nbrs[ndx]._dist = DistanceToInt(nbr.dist());
    ++ndx;
    if (ndx >= ok_nbrs) {
      break;
    }
  }

  return 1;
}


// For each neighbour update the distance stored in `nbrdist`

void
Needle::SetNbrDistances(uint32_t* nbrdist) const {
  // For efficiency we duplicate the two loops so we do not
  // have an if statement deep inside the loop.
  if (_count == 1) [[ likely ]] {
    for (uint32_t i = 0; i < _number_nbrs; ++i) {
      const Nbr& nbr = _nbrs[i];
      uint32_t id = nbr._id;
      nbrdist[id] = nbr._dist;
    }
  } else {
    for (uint32_t i = 0; i < _number_nbrs; ++i) {
      const Nbr& nbr = _nbrs[i];
      uint32_t id = nbr._id;
      nbrdist[id] = nbr._dist * _count;
    }
  }
}


#ifdef NO_LONGER_USED
// For each neighbour update the distance stored in `nbrdist`
void
Needle::SetNbrDistances(uint32_t* nbrdist, uint32_t max_distance) const {
  for (uint32_t i = 0; i < _number_nbrs; ++i) {
    const Nbr& nbr = _nbrs[i];
    uint32_t id = nbr._id;
    if (nbrdist[id] == max_distance) {
      nbrdist[id] = nbr._dist;
    } else {
      nbrdist[id] = 0;
    }
  }
}
#endif

// Used for debugging.
void
Needle::WriteNbrs(std::ostream& output) const {
  output << "train " << _train << '\n';
  for (uint32_t i = 0; i < _number_nbrs; ++i) {
    output << "  " << i << " id " << _nbrs[i]._id << " dist " << _nbrs[i]._dist << '\n';
  }
}

float
Needle::ClosestDistance() const {
  if (_number_nbrs == 0) {
    return 101.0;
  }

  return _nbrs[0]._dist;
}

class Metropolis {
  private:
  public:
};

// Handle the overall optimisation of train/test splits.
class Optimise {
  private:
    int _verbose;

    // The number of splits to form, the -n option.
    int _nsplits;

    uint32_t _number_needles;
    Needle* _needle;

    // Does not work, do not used.
    float _upper_distance_threshold;

    // The -f option.
    float _fraction_train;

    absl::flat_hash_map<IWString, uint32_t> _id_to_count;

    std::mt19937 _rng;

    std::unique_ptr<std::uniform_int_distribution<uint32_t>> _uniform;

    // The sum of the distances between train and test. We want to maximise this.
    // If each distance is scaled to 100, then the largest job would
    // be N*(N-1)/2*100>.
    // Therefore use 64 bit int.
    uint64_t _current_score;

    // The number of optimisation steps per split, the -o option.
    uint32_t _nopt;

    // If running for a fixed time per split, the -t option.
    uint32_t _seconds;

    // Are we optimising for separation between train and test or closeness
    bool _optimise_for_separation;

    // For each pair of items where we know a distance, store that.
    // Key will be i * _number_needles + j where i > j.
    std::unordered_map<uint64_t, uint32_t> _distance;

    // One more than the max value of the values in _distance
    uint32_t _max_distance;

    int _nthreads;

    // Rather than storing smiles and id with the Needle, we deliberately
    // make a global array of both, hoping to keep the Needle class very
    // small.
    std::string* _smiles;
    std::string* _name;

    // If the -A option is specified.
    ActivityProfileSet _activity;

    // Global weight for mixing activity-bucket balance into the objective.
    uint32_t _activity_objective_weight;

    // Per-molecule multiplicity from the -C option. Usually all values are 1.
    std::vector<uint32_t> _sample_weight;

    // Keep track of how many times each needle appears in train.
    int* _times_in_train;

    // The file name stem for the output files, the -S option.
    IWString _stem;

    // This can be expensive so it is optional.
    int _write_cross_split_summary;

    // Activated by the -r option.
    Report_Progress _report_progress;

    // If nothing is accepted for this many steps, abandon optimisation.
    uint32_t _abandon_if;

  // Private functions.
    uint32_t GatherIdentifiers(const char* fname,
                            std::unordered_map<std::string, std::uint32_t>& id_to_ndx);
    uint32_t ReadNeighbours(const char* fname,
                            const std::unordered_map<std::string, uint32_t>& id_to_ndx);
    int ReadCounts(IWString& fname);
    int ReadCounts(iwstring_data_source& input);
    int ReadCount(const const_IWSubstring& buffer);

    uint32_t GetCount(const IWString& id) const;

    int RandomSplit();

    uint64_t FormKey(uint32_t i, uint32_t j) const;

    std::optional<uint32_t> Distance(uint32_t i, uint32_t j) const;
    uint32_t DistanceOrMax(uint32_t i, uint32_t j) const;

    uint64_t RecomputeCurrentScore();

    float AveDistanceAcrossSplit() const;
    std::tuple<float, uint32_t> AveDistNumberMax() const;

    std::tuple<uint32_t, uint32_t> ChooseTwo();

    int64_t CombinedObjectiveDelta(int64_t distance_delta, uint32_t out_of_train,
                                   uint32_t into_train) const;
    uint64_t ScoreAfterDelta(uint64_t score, int64_t delta) const;

    int MakeSplit(int split);

    int WriteSplit(int split) const;
    int WriteSmiles(int train, IWString& fname) const;
    int WriteSmiles(int train, IWString_and_File_Descriptor& output) const;
    int WriteSubset(int train, IWString& fname) const;
    int WriteSubset(int train, IWString_and_File_Descriptor& output) const;
    int WriteCrossSplitSummary(IWString& fname) const;
    int WriteCrossSplitSummary(IWString_and_File_Descriptor& output) const;
    int WriteRandomSplitStatus() const;
    void AccumulateStats();

    int MaybeReportActivityDistribution(std::ostream& output) const;
    std::vector<int> CurrentTrainMembership() const;

  public:
    Optimise();
    ~Optimise();

    int Initialise(const Command_Line& cl);

    int ReportSize(std::ostream& output) const;

    int ReadFingerprints(const char* fname);
    int BuildActivityProfiles();

    int Doit();

    int Report(std::ostream& output) const;
};

Optimise::Optimise() {
  _verbose = 0;
  _nsplits = 10;
  _fraction_train = 0.80;
  _nopt = 1000;

  _upper_distance_threshold = std::numeric_limits<float>::max();

  _max_distance = 0;

  _optimise_for_separation = true;

  _number_needles = 0;
  _needle = nullptr;

  _nthreads = 0;

  _activity_objective_weight = 0;

  _smiles = nullptr;
  _name = nullptr;
  _times_in_train = nullptr;

  _abandon_if = 0;

  _write_cross_split_summary = 0;

  std::random_device rd;
  _rng.seed(rd());
}

Optimise::~Optimise() {
  if (_needle != nullptr) {
    delete [] _needle;
  }
  if (_smiles != nullptr) {
    delete [] _smiles;
  }
  if (_name != nullptr) {
    delete [] _name;
  }

  delete [] _times_in_train;
}

int
Optimise::Initialise(const Command_Line& cl) {
  _verbose = cl.option_present('v');

  if(cl.option_present('n')) {
    if (! cl.value('n', _nsplits) || _nsplits < 1) {
      cerr << "Invalid nsplit (-n)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will create " << _nsplits << " splits\n";
    }
  }

  if (cl.option_present('f')) {
    if (! cl.value('f', _fraction_train) || _fraction_train <= 0.0f ||
         _fraction_train >= 1.0f) {
      cerr << "Optimise::Initialise:Invalid fraction train (-f)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will place " << _fraction_train << " of items in the training set\n";
    }
  }

  if (cl.option_present('c')) {
    _optimise_for_separation = false;
    if (_verbose) {
      cerr << "Will optimise for closeness between train and test\n";
    }
  }

  const int seconds_option_present = cl.option_present('s') || cl.option_present('t');
  if (cl.option_present('s') && cl.option_present('t')) {
    cerr << "The -s and -t options both specify seconds, use only one\n";
    return 0;
  }

  if (cl.option_present('o') && seconds_option_present) {
    cerr << "The -o (nopt) and -s (seconds) options are mutually incompatible\n";
    return 0;
  }

  if (cl.option_present('o')) {
    if (! cl.value('o', _nopt) || _nopt < 1) {
      cerr << "Optimise::Initialise:Invalid nopt (-o)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will perform " << _nopt << " optimisation attempts\n";
    }
  }

  if (cl.option_present('x')) {
    if (! cl.value('x', _abandon_if) || _abandon_if < 100) {  // 100 is arbitrary
      cerr << "The abandon if too long since last successful switch (-x) must be a number > 100\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will abandon optimistion of " << _abandon_if << " unsuccessful steps have been tried\n";
    }
  }

  if (seconds_option_present) {
    const char seconds_option = cl.option_present('s') ? 's' : 't';
    if (! cl.value(seconds_option, _seconds) || _seconds < 1) {
      cerr << "Optimise::Initialise:Invalid seconds (-" << seconds_option << ")\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will perform optimisation of each split for at least " << _seconds << " seconds\n";
    }

    _nopt = std::numeric_limits<uint32_t>::max();
  }

  if (cl.option_present('r')) {
    if (! _report_progress.initialise(cl, 'r', _verbose)) {
      cerr << "Cannot initialise progress reporting (-r)\n";
      return 0;
    }
  }

  if (cl.option_present('T')) {
    if (! cl.value('T', _upper_distance_threshold) || _upper_distance_threshold <= 0.0f ||
        _upper_distance_threshold > 1.0f) {
      cerr << "Invalid upper distance threshold (-T)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will drop nbrs at distances > " << _upper_distance_threshold << '\n';
    }
  }

  if (cl.option_present('C')) {
    IWString fname = cl.string_value('C');
    if (! ReadCounts(fname)) {
      cerr << "Cannot read count data '" << fname << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Read " << _id_to_count.size() << " count data records from '" << fname << "'\n";
    }
  }

  if (cl.option_present('S')) {
    cl.value('S', _stem);
    if (_verbose) {
      cerr << "Split files created with stem '" << _stem << "'\n";
    }
  }

  if (cl.option_present('h')) {
    if (! cl.value('h', _nthreads) || _nthreads < 1) {
      cerr << "The number of OMP threads (-h) must be a whole +ve number\n";
      return 0;
    }
    omp_set_num_threads(_nthreads);
    if (_verbose) {
      cerr << "Will use " << _nthreads << " OMP threads for exact calculation\n";
    }
  }

  if (cl.option_present('X')) {
    _write_cross_split_summary = 1;
    if (_verbose) {
      cerr << "Will write cross split summary files\n";
    }
  }

  if (cl.option_present('Z')) {
    const_IWSubstring z;
    for (int i = 0; cl.value('Z', z, i); ++i) {
      if (z == "help") {
        DisplayMiscOptions(0);
      }

      if (z.starts_with("seed=")) {
        z.remove_leading_chars(5);
        uint32_t seed;
        if (! z.numeric_value(seed)) {
          cerr << "Invalid random number seed (-Z seed=) '" << z << "'\n";
          return 0;
        }
        _rng.seed(seed);
        if (_verbose) {
          cerr << "Random number seed " << seed << '\n';
        }
      } else {
        cerr << "Unrecognised -Z qualifier '" << z << "'\n";
        DisplayMiscOptions(1);
      }
    }
  }

  if (cl.option_present('A')) {
    for (int i = 0; i < cl.option_count('A'); ++i) {
      const const_IWSubstring directive = cl.string_value('A', i);
      if (! _activity.AddDirective(directive)) {
        cerr << "Cannot process activity directive '" << directive << "'\n";
        return 0;
      }
    }
    if (_verbose) {
      cerr << "Read " << cl.option_count('A') << " activity directives\n";
    }
  }

  if (cl.option_present('Y')) {
    if (! cl.value('Y', _activity_objective_weight) || _activity_objective_weight == 0) {
      cerr << "The activity objective weight (-Y) must be a whole positive number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Activity distribution objective weight " << _activity_objective_weight << '\n';
    }
  }

  if (cl.option_present('A') && ! _activity.has_pending()) {
    cerr << "The activity option (-A) must specify at least one activity file\n";
    return 0;
  }

  if (_activity.has_pending() && _activity_objective_weight == 0) {
    cerr << "Activity files specified with -A, but no activity objective weight specified with -Y\n";
    return 0;
  }

  if (_activity_objective_weight > 0 && ! _activity.has_pending()) {
    cerr << "The activity objective weight (-Y) requires at least one activity file (-A)\n";
    return 0;
  }

  return 1;
}

int
Optimise::ReportSize(std::ostream& output) const {
  output << "Optimise:holds " << _number_needles << " needles\n";
  extending_resizable_array<uint32_t> nbr_count;
  for (uint32_t i = 0; i < _number_needles; ++i) {
    uint32_t n = _needle[i].number_neighbours();
    ++nbr_count[n];
  }

  for (uint32_t i = 0; i < nbr_count.size(); ++i) {
    if (nbr_count[i]) {
      output << nbr_count[i] << " items had " << i << " nbrs\n";
    }
  }

  output << "Distance matrix holds " << _distance.size() << " items, max dist " << _max_distance << '\n';

  return 1;
}

int
Optimise::ReadCounts(IWString& fname) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Optimise::ReadCounts:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadCounts(input);
}

int
Optimise::ReadCounts(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Optimise::ReadCounts:cannot read header record\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! ReadCount(buffer)) {
      cerr << "Optimise::ReadCounts:invalid data '" << buffer << "'\n";
      return 0;
    }
  }

  return _id_to_count.size();
}

int
Optimise::ReadCount(const const_IWSubstring& buffer) {
  static constexpr char kInputSeparator = ' ';

  int i = 0;
  IWString id;
  const_IWSubstring scount;
  if (! buffer.nextword(id, i, kInputSeparator) ||
      ! buffer.nextword(scount, i, kInputSeparator) ||
      id.empty() || scount.empty()) {
    cerr << "Optimise::ReadCount:cannot split\n";
    return 0;
  }

  uint32_t count;
  if (! scount.numeric_value(count) || count == 0) {
    cerr << "Optimise::ReadCount:invalid count " << count << '\n';
    return 0;
  }

  _id_to_count[id] = count;

  return 1;
}

std::optional<uint32_t>
Optimise::Distance(uint32_t i, uint32_t j) const {
  uint64_t key = FormKey(i, j);

  const auto iter = _distance.find(key);
  if (iter == _distance.end()) {
    return std::nullopt;
  }

  return iter->second;
}

uint32_t
Optimise::DistanceOrMax(uint32_t i, uint32_t j) const {
  const uint64_t key = FormKey(i, j);

  const auto iter = _distance.find(key);
  if (iter == _distance.end()) {
    // cerr << "No distance for " << i << ',' << j << '\n';
    return _max_distance;
  }

  return iter->second;
}

int
Optimise::RandomSplit() {
  std::bernoulli_distribution dist(_fraction_train);
  uint32_t in_train = 0;
  for (uint32_t i = 0; i < _number_needles; ++i) {
    if (dist(_rng)) {
      _needle[i].set_train(1);
      ++in_train;
    } else {
      _needle[i].set_train(0);
    }
  }

  if (_verbose) {
    cerr << "RandomSplit puts " << in_train << " in train\n";
  }

  // do some balancing.
  uint32_t expected_in_train = static_cast<uint32_t>(_number_needles * _fraction_train);
  if (in_train > expected_in_train) {
    uint32_t delta = in_train - expected_in_train;
    for (uint32_t i = 0; i < delta; ++i) {
      uint32_t x = (*_uniform)(_rng);
      if (_needle[x].in_train()) {
        _needle[x].set_train(0);
        --in_train;
      } else {
        --i;
      }
    }
  } else if (in_train < expected_in_train) {
    uint32_t delta = expected_in_train - in_train;
    for (uint32_t i = 0; i < delta; ++i) {
      uint32_t x = (*_uniform)(_rng);
      if (! _needle[x].in_train()) {
        _needle[x].set_train(1);
        ++in_train;
      } else {
        --i;
      }
    }
  }

  if (_verbose) {
    cerr << "RandomSplit after balancing " << in_train << " in train\n";
  }

#ifdef DEBUG_SWAP_ITEMS
  for (uint32_t i = 0; i < _number_needles; ++i) {
    cerr << " train " << i << ' ' << _needle[i].in_train() << '\n';
  }
#endif

  return 1;
}

uint64_t
Optimise::RecomputeCurrentScore() {
  _current_score = 0;

  uint64_t accumulate = 0;
  cerr << " Begin RecomputeCurrentScore\n";
#pragma omp parallel for reduction(+:accumulate)
  for (uint32_t i = 0; i < _number_needles; ++i) {
    const int itrain = _needle[i].in_train();
    for (uint32_t j = i + 1; j < _number_needles; ++j) {
      const int jtrain = _needle[j].in_train();
      if (itrain == jtrain) {
        continue;
      }

      const uint32_t d = DistanceOrMax(i, j);
      // cerr << " " << i << " and " << j << " opposite, dist " << d << '\n';

      accumulate += d;
    }
  }
  _current_score = accumulate;

  cerr << "Current score " << _current_score << '\n';

  return _current_score;
}

// Reading fingerprints is done twice. On the first pass, we determine
// the number of fingerprints, and establish the mapping from id to index.
int
Optimise::ReadFingerprints(const char* fname) {
  std::unordered_map<std::string, std::uint32_t> id_to_ndx;
  if (! GatherIdentifiers(fname, id_to_ndx)) {
    cerr << "Optimise::ReadFingerprints:cannot read '" << fname << "'\n";
    return 0;
  }

  return ReadNeighbours(fname, id_to_ndx);
}

uint32_t
Optimise::GatherIdentifiers(const char* fname,
                            std::unordered_map<std::string, std::uint32_t>& id_to_ndx) {
  iw_tf_data_record::TFDataReader reader;
  if (! reader.Open(fname)) {
    cerr << "Optimise::GatherIdentifiers:cannot open '" << fname << "'\n";
    return 0;
  }

  uint32_t ndx;
  for (ndx = 0; ; ++ndx) {
    std::optional<nnbr::NearNeighbours> needle = reader.ReadProto<nnbr::NearNeighbours>();
    if (! needle) {
      break;
    }
    // No checking for duplicates.
    id_to_ndx[needle->name()] = ndx;
  }

  if (id_to_ndx.size() != ndx) {
    cerr << "Optimise::GatherIdentifiers:size mismatch, read " << ndx << " records, but only " << id_to_ndx.size() << 
            " unique identifiers\n";
    return 0;
  }

  return ndx;
}

uint64_t
Optimise::FormKey(uint32_t i, uint32_t j) const {
  if (i < j) {
    return j * _number_needles + i;
  }

  return i * _number_needles + j;
}

// Choose 2 different members of the set that are across the train/split line.
std::tuple<uint32_t, uint32_t>
Optimise::ChooseTwo() {
  const uint32_t i1 = (*_uniform)(_rng);
  const int itrain  = _needle[i1].in_train();

  uint32_t i2 = (*_uniform)(_rng);
  int jtrain = _needle[i2].in_train();

  while (jtrain == itrain || i1 == i2) {
    i2 = (*_uniform)(_rng);
    jtrain = _needle[i2].in_train();
  }

  return std::make_tuple(i1, i2);
}

uint32_t
Optimise::GetCount(const IWString& id) const {
  if (_id_to_count.empty()) {
    return 1;
  }

  const auto iter = _id_to_count.find(id);
  if (iter == _id_to_count.end()) {
    cerr << "Optimise::GetCount:no count for '" << id << "'\n";
    return 1;
  }

  return iter->second;
}


uint32_t
Optimise::ReadNeighbours(const char* fname,
                         const std::unordered_map<std::string, uint32_t>& id_to_ndx) {
  iw_tf_data_record::TFDataReader reader;
  if (! reader.Open(fname)) {
    cerr << "Optimise::ReadNeighbours:cannot open '" << fname << "'\n";
    return 0;
  }

  _number_needles = id_to_ndx.size();
  _needle = new Needle[_number_needles];
  _smiles = new std::string[_number_needles];
  _name = new std::string[_number_needles];
  _times_in_train = new_int(_number_needles);
  _sample_weight.resize(_number_needles, 1);

  _uniform = std::make_unique<std::uniform_int_distribution<uint32_t>>(0, _number_needles - 1);

  uint32_t ndx;
  for (ndx = 0; ; ++ndx) {
    std::optional<nnbr::NearNeighbours> needle = reader.ReadProto<nnbr::NearNeighbours>();
    if (! needle) {
      break;
    }

    _smiles[ndx] = needle->smiles();
    // Do not store multi-token names, truncate to first token.
    IWString tmp = needle->name();
    tmp.truncate_at_first(' ');
    _name[ndx] = tmp.AsString();
    if (!_needle[ndx].Build(*needle, id_to_ndx, _upper_distance_threshold)) {
      cerr << "Optimise::ReadNeighbours:cannot process " << needle->ShortDebugString() << "\n";
      return 0;
    }

    const uint32_t count1 = GetCount(tmp);
    _needle[ndx].set_count(count1);
    _sample_weight[ndx] = count1;

    for (const nnbr::Nbr& nbr : needle->nbr()) {
      const auto iter = id_to_ndx.find(nbr.id());
      uint64_t key = FormKey(ndx, iter->second);
      const uint64_t count2 = GetCount(nbr.id());

      _distance[key] = DistanceToInt(nbr.dist()) * (count1 * count2);
    }
  }

  _max_distance = 0;
  for (const auto& [_, dist] : _distance) {
    if (dist > _max_distance) {
      _max_distance = dist;
    }
  }

  ++_max_distance;

#ifdef ECHO_DISTANCE
  for (const auto [k,v] : _distance) {
    cerr << " k " << k << " dist " << v << '\n';
  }
#endif

  if (ndx != _number_needles) {
    cerr << "Optimise::ReadNeighbours:size mismatch, expected " << _number_needles << " got " << ndx << '\n';
    return 0;
  }

#ifdef DEBUG_READ_NEIGHBOURS
  for (uint64_t q = 0; q < _number_needles; ++q) {
    cerr << "q " << q << ' ' << _name[q] << '\n';
    _needle[q].WriteNbrs(cerr);
  }
#endif

  return _number_needles;
}

int
Optimise::BuildActivityProfiles() {
  if (! _activity.has_pending()) {
    return 1;
  }

  std::vector<std::string> id;
  id.reserve(_number_needles);
  for (uint32_t i = 0; i < _number_needles; ++i) {
    id.push_back(_name[i]);
  }

  if (! _activity.Build(id, _sample_weight, _verbose)) {
    cerr << "Optimise::BuildActivityProfiles:cannot build activity profiles\n";
    return 0;
  }

  if (_verbose) {
    cerr << "Read " << _activity.number_profiles() << " activity profiles\n";
  }

  return 1;
}

std::vector<int>
Optimise::CurrentTrainMembership() const {
  std::vector<int> result(_number_needles);
  for (uint32_t i = 0; i < _number_needles; ++i) {
    result[i] = _needle[i].in_train();
  }

  return result;
}

// Return true if we are now > `seconds` seconds from `tzero`.
bool
BreakForTime(std::time_t tzero, uint32_t seconds) {
  return (std::time(nullptr) - tzero) >= seconds;
}

// Generate _nsplits splits.
int
Optimise::Doit() {
  for (int i = 0; i < _nsplits; ++i) {
    MakeSplit(i);
  }

  return 1;
}

// Return a signed diff between two unsigned numbers.
int64_t
Diff(const uint64_t v1, const uint64_t v2) {
  if (v1 > v2) {
    return v1 - v2;
  }
  return - static_cast<int64_t>(v2 - v1);
}

int64_t
Optimise::CombinedObjectiveDelta(int64_t distance_delta, uint32_t out_of_train,
                                 uint32_t into_train) const {
  if (_activity_objective_weight == 0) {
    return distance_delta;
  }

  const int64_t activity_delta = _activity.DeltaForSwap(out_of_train, into_train);
  return CombinedActivityDistanceDelta(distance_delta, activity_delta,
                                       _activity_objective_weight);
}

uint64_t
Optimise::ScoreAfterDelta(uint64_t score, int64_t delta) const {
  if (delta >= 0) {
    return score + delta;
  }

  return score - static_cast<uint64_t>(-delta);
}

// #define DEBUG_MAKE_SPLIT

// Starting with a random split, optimize it and write when done.
// This is still not correct.
// If run with verbose, after the optimisation the full score is recomputed
// and there are differences.
// Failures are rare, I cannot generate any failures with a random set of
// 10k molecules from Chembl.
// But 100k reliably generates failures.
// The resulting splits seem to be just fine.
int
Optimise::MakeSplit(int split) {
  std::time_t tzero = 0;
  if (_seconds) {
    tzero = std::time(nullptr);
  }

  RandomSplit();
  if (_activity.active()) {
    const std::vector<int> in_train = CurrentTrainMembership();
    uint64_t total_sample_weight = 0;
    uint64_t train_sample_weight = 0;
    for (uint32_t i = 0; i < _number_needles; ++i) {
      total_sample_weight += _sample_weight[i];
      if (in_train[i]) {
        train_sample_weight += _sample_weight[i];
      }
    }
    if (! _activity.InitialiseSplit(in_train, _sample_weight, total_sample_weight, train_sample_weight)) {
      return 0;
    }
  }
  uint64_t score = RecomputeCurrentScore();
  const uint64_t starting_score = score;
  if (_verbose) {
    const auto [across_split, number_at_max] = AveDistNumberMax();
    cerr << "Split " << split << " starting score " << score <<
            " ave dist " << across_split << " at max " << number_at_max << '\n';
  }

  // FIrst split, write the random split.
  if (split == 0 && _verbose) {
    cerr << "WRiting random split status\n";
    WriteRandomSplitStatus();
  }

  std::unique_ptr<uint32_t[]> id_dist = std::make_unique<uint32_t[]>(_number_needles);

  uint32_t steps_accepted = 0;
  uint32_t last_successful_switch = 0;
  cerr << "Begin " << _nopt << " optimisation steps\n";
#ifdef DEBUG_MAKE_SPLIT
  cerr << "Will perform " << _nopt << " optimisations\n";
  for (int i = 0; i < _number_needles; ++i) {
    cerr << i << " train " << _needle[i].in_train() << '\n';
  }
#endif

  for (uint32_t j = 0; j < _nopt; ++j) {
    auto [i1, i2] = ChooseTwo();
    const uint32_t out_of_train = _needle[i1].in_train() ? i1 : i2;
    const uint32_t into_train = _needle[i1].in_train() ? i2 : i1;
    _needle[i1].invert_train();
    _needle[i2].invert_train();

    const int train1 = _needle[i1].in_train();
    const int train2 = _needle[i2].in_train();
#ifdef DEBUG_MAKE_SPLIT
    cerr << "Selected " << i1 << " and " << i2 << '\n';
    for (uint32_t i = 0; i < _number_needles; ++i) {
      cerr << i << " train " << _needle[i].in_train() << '\n';
    }
#endif

    std::fill_n(id_dist.get(), _number_needles, _max_distance);
    _needle[i1].SetNbrDistances(id_dist.get());
    id_dist[i1] = 0;
    id_dist[i2] = 0;
#ifdef DEBUG_MAKE_SPLIT
    cerr << "After setting needle1 ";
    for (uint32_t q = 0; q < _number_needles; ++q) {
      cerr << ' ' << id_dist[q];
    }
    cerr << '\n';
#endif

    int64_t delta = 0;   // a signed quantity
#pragma omp parallel for reduction(+:delta)
    for (uint32_t i = 0; i < _number_needles; ++i) {
//    if (i == i1 || i == i2) [[ unlikely ]] {   not needed because the id_dist entries have been zero'd
//      continue;
//    }

      int32_t d = id_dist[i];

      // If both on the same side now, were previously on opposite sides.
      if (_needle[i].in_train() == train1) {
        delta += -d;
      } else  {  // different sides now, previously same.
        delta += d;
      }
      // cerr << "  delta updated to " << delta << '\n';
    }
#ifdef DEBUG_MAKE_SPLIT
    cerr << "At end of i1 delta " << delta << '\n';
#endif

    std::fill_n(id_dist.get(), _number_needles, _max_distance);
    _needle[i2].SetNbrDistances(id_dist.get());
    id_dist[i1] = 0;
    id_dist[i2] = 0;
#ifdef DEBUG_MAKE_SPLIT
    cerr << "After setting needle2 ";
    for (uint32_t q = 0; q < _number_needles; ++q) {
      cerr << ' ' << id_dist[q];
    }
    cerr << '\n';
#endif

#pragma omp parallel for reduction(+:delta)
    for (uint32_t i = 0; i < _number_needles; ++i) {
//    if (i == i1 || i == i2) [[ unlikely ]] {
//      continue;
//    }

      int32_t d = id_dist[i];
      // cerr <<  "  adj " << i2 << " to " << n << " dist " << d << " train2 " << train2 << " nbr " << _needle[n].in_train() << '\n';

      // If both on the same side now, were previously on opposite sides.
      if (_needle[i].in_train() == train2) {
        delta += -d;
      } else  {  // different sides now, previously same.
        delta += d;
      }
      // cerr << "   2nd loop, delta " << delta << '\n';
    }
#ifdef DEBUG_MAKE_SPLIT
    cerr << "At end of i2 delta " << delta << '\n';
#endif

    const int64_t combined_delta = CombinedObjectiveDelta(delta, out_of_train, into_train);
    const uint64_t new_score = ScoreAfterDelta(score, delta);
#ifdef DEBUG_SWAP_ITEMS
    cerr << "new_score " << new_score << " cmp " << score << " distance_delta " << delta <<
            " combined_delta " << combined_delta << '\n';
    for (uint32_t y = 0; y < _number_needles; ++y) {
      cerr << " needle " << y << " train " << _needle[y].in_train() << '\n';
    }
#endif

    if (_report_progress() && j > 0) {
      static uint32_t accepted_last_report = 0;
      uint32_t accepted_this_period = steps_accepted - accepted_last_report;
      accepted_last_report = steps_accepted;

      cerr << split << ' ' << j << " score " << new_score << " cmp " << starting_score <<
      " accepted " << steps_accepted << " delta " << accepted_this_period <<
      " distance_delta " << delta;
      if (_activity.active()) {
        cerr << " combined_delta " << combined_delta;
      }
      cerr << ' ' << iwmisc::Fraction<float>(steps_accepted, j) <<
      " last successful " << last_successful_switch << '\n';
      MaybeReportActivityDistribution(cerr);
    }

    bool accept_change = false;
    if (_optimise_for_separation && combined_delta > 0) {
      accept_change = true;
    } else if (! _optimise_for_separation && combined_delta < 0) {
      accept_change = true;
    }
    if (accept_change) {
      score = new_score;
      if (_activity.active()) {
        _activity.ApplySwap(out_of_train, into_train);
      }
      ++steps_accepted;
      last_successful_switch = j;
    } else {
      _needle[i1].invert_train();
      _needle[i2].invert_train();
    }

#ifdef DEBUG_SWAP_ITEMS
    auto r = RecomputeCurrentScore();
    cerr << j << " compute " << r << " new_score " << new_score << '\n';
    new_score = r;
#endif

    // Heuristics to stop this being checked too often.
    if (_seconds > 0 && j > 500 && j % 1000 == 0 && BreakForTime(tzero, _seconds)) [[ unlikely ]] {
      cerr << "Stopped for  timeout " << j << " steps completed\n";
      break;
    }

    if (_abandon_if > 0 && (j - last_successful_switch > _abandon_if)) [[ unlikely ]] {
      cerr << "Abandon optimisation, steps " << j << " last successful " << last_successful_switch << '\n';
      cerr << "Accepted " << steps_accepted << " switches\n";
      break;
    }
  }

  if (_verbose) {
    if (_activity.active()) {
      // If activity data is present, the score includes the influence of activity, so recomputation is not meaningful.
    } else {
      RecomputeCurrentScore();
      cerr << "Writing split " << split << " score " << score << " computed " << _current_score << " diff " << Diff(score, _current_score) << '\n';
    }

    const auto [across_split, number_at_max] = AveDistNumberMax();
    cerr << "starting_score " << starting_score << " score " << _current_score <<
            " improvement " << Diff(_current_score, starting_score) <<
            " across split " << AveDistanceAcrossSplit() << " at max " << number_at_max << '\n';
    cerr << "Split " << split << " accepted " << steps_accepted << " of " << _nopt << " steps\n";
  }

  AccumulateStats();

  return WriteSplit(split);
}

// Accumulate how many times each item appears in the training split.
void
Optimise::AccumulateStats() {
  for (uint32_t i = 0; i < _number_needles; ++i) {
    if (_needle[i].in_train()) {
      ++_times_in_train[i];
    }
  }
}

int
Optimise::MaybeReportActivityDistribution(std::ostream& output) const {
  if (! _activity.active()) {
    return 1;
  }

  output << "activity penalty " << _activity.current_penalty() << '\n';

  return output.good();
}

// Write the various split files, train and test identifiers, and smiles files.
// Distribution file as well.
int
Optimise::WriteSplit(int split) const {
  IWString fname;
  fname << _stem << 'R' << split;
  WriteSubset(1, fname);

  fname.resize_keep_storage(0);
  fname << _stem << 'E' << split;
  WriteSubset(0, fname);

  fname.resize_keep_storage(0);
  fname << _stem << 'R' << split << ".smi";
  WriteSmiles(1, fname);

  fname.resize_keep_storage(0);
  fname << _stem << 'E' << split << ".smi";
  WriteSmiles(0, fname);

  if (_write_cross_split_summary) {
    fname.resize_keep_storage(0);
    fname << _stem << "_stats" << split << ".txt";
    WriteCrossSplitSummary(fname);
  }

  return 1;
}

int
Optimise::WriteRandomSplitStatus() const {
  if (! _write_cross_split_summary) {
    return 1;
  }

  IWString fname;
  fname << _stem << "_stats_rand.txt";
  return WriteCrossSplitSummary(fname);
}

int
Optimise::WriteSubset(int train, IWString& fname) const {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "Optimise::WriteSubset:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteSubset(train, output);
}

// Writes out the _name[i] items that match `train`.
int
Optimise::WriteSubset(int train, IWString_and_File_Descriptor& output) const {
  for (uint32_t i = 0; i < _number_needles; ++i) {
    if (_needle[i].in_train() == train) {
      output << _name[i] << '\n';
      output.write_if_buffer_holds_more_than(8192);
    }
  }

  return 1;
}

int
Optimise::WriteSmiles(int train, IWString& fname) const {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "Optimise::WriteSubset:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteSmiles(train, output);
}

int
Optimise::WriteSmiles(int train, IWString_and_File_Descriptor& output) const {
  for (uint32_t i = 0; i < _number_needles; ++i) {
    if (_needle[i].in_train() == train) {
      output << _smiles[i] << ' ' << _name[i] << '\n';
      output.write_if_buffer_holds_more_than(8192);
    }
  }

  return 1;
}

int
Optimise::WriteCrossSplitSummary(IWString& fname) const {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "Optimise::WriteCrossSplitSummary:cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteCrossSplitSummary(output);
}

int
Optimise::WriteCrossSplitSummary(IWString_and_File_Descriptor& output) const {

  extending_resizable_array<int> count;

  count.reserve(_max_distance + 1);

  cerr << "WriteCrossSplitSummary\n";
  for (uint32_t i = 0; i < _number_needles; ++i) {
    int itrain = _needle[i].in_train();
    for (uint32_t j = i + 1; j < _number_needles; ++j) {
      if (_needle[j].in_train() == itrain) {
        continue;
      }
      uint32_t d = DistanceOrMax(i, j);
      ++count[d];
    }
  }
  cerr << "WriteCrossSplitSummary finished\n";

  constexpr char kSep = ' ';

  output << "Dist" << kSep << "Count" << kSep << "Cumulative" << '\n';

  uint32_t sum = 0;
  for (int i = 0; i < count.number_elements(); ++i) {
    if (count[i]) {
      sum += count[i];
      output << iwmisc::Fraction<float>(i, 100) << kSep << count[i] << kSep << sum << '\n';
    }
  }

  return 1;
}

float
Optimise::AveDistanceAcrossSplit() const {
  uint64_t sum = 0;
  uint64_t n = 1;

#pragma omp parallel for reduction(+:sum, n)
  for (uint64_t i = 0; i < _number_needles; ++i) {
    const int itrain = _needle[i].in_train();
    for (uint64_t j = i + 1; j <_number_needles; ++j) {
      if (_needle[j].in_train() == itrain) {
        continue;
      }
      uint32_t d = DistanceOrMax(i, j);
      sum += d;
      ++n;
    }
  }

  // Do the calc as a double, will be cast to float on return.
  return iwmisc::Fraction<double>(sum, n);
}

// Return the average distance of cross-split pairs that are not at the
// maximum distance, and the number that are at the maximum distance.
std::tuple<float, uint32_t>
Optimise::AveDistNumberMax() const {
  uint64_t number_max = 0;
  uint64_t sum = 0;
  uint64_t n = 0;

#pragma omp parallel for reduction(+:sum, number_max)
  for (uint64_t i = 0; i < _number_needles; ++i) {
    const int itrain = _needle[i].in_train();
    for (uint64_t j = i + 1; j <_number_needles; ++j) {
      if (_needle[j].in_train() == itrain) {
        continue;
      }
      uint32_t d = DistanceOrMax(i, j);
      if (d == _max_distance) {
        ++number_max;
      } else {
        sum += d;
        ++n;
      }
    }
  }

  float ave = iwmisc::Fraction<double>(sum, n);

  return std::tuple(ave, number_max);
}

int
Optimise::Report(std::ostream& output) const {
  static constexpr char kSep = ' ';

  output << "generated " << _nsplits << " splits\n";
  output << "ID" << kSep << "TimesInTrain" << kSep << "NNDist" << '\n';

  for (uint32_t i = 0; i < _number_needles; ++i) {
    output << _name[i] << kSep << _times_in_train[i] << kSep << (_needle[i].ClosestDistance() / 100.0f) << '\n';
  }

  return output.good();
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vf:S:n:o:r:T:s:t:x:C:h:XA:Y:Z:c");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_present('v');

  if (! cl.option_present('S')) {
    cerr << "Must specify output file name stem via the -S option\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient aguments\n";
    Usage(1);
  }

  Optimise optimise;
  if (! optimise.Initialise(cl)) {
    cerr << "Cannot initialise calculation\n";
    return 1;
  }

  if (! optimise.ReadFingerprints(cl[0])) {
    cerr << "Cannot read TFDataRecord nearneighbour proto data '" << cl[0] << "'\n";
    return 1;
  }

  if (! optimise.BuildActivityProfiles()) {
    cerr << "Cannot build activity profiles\n";
    return 1;
  }

  if (verbose) {
    optimise.ReportSize(cerr);
  }

  optimise.Doit();

  if (verbose) {
    optimise.Report(cerr);
  }

  return 0;
}

}  // namespace train_test_split

int
main(int argc, char** argv) {
  int rc = train_test_split::Main(argc, argv);
  return rc;
}
