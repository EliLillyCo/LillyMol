// Use the results of a nearest neighbour calculation to optimize
// one or more train/test splits.
// In its current state this is not working. Needs debugging.
//  TODO:ianwatson fix this.

#include <stdlib.h>

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

#define REPORT_PROGRESS_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Utilities/GFP_Tools/nearneighbours.pb.h"

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
  cerr << "Uses TFDataRecord nn data to generate optimized train/test splits\n";
  cerr << " -f <train_fraction> fraction of data in the train split\n";
  cerr << " -n <nsplit>         number of splits needed\n";
  cerr << " -S <stem>           write splits to <stem>, <stem>R and <stem>E\n";
  //cerr << " -s ...              stratified sampling...\n";
  cerr << " -o <nopt>           number of optimisation steps to try per split\n";
  cerr << " -t <sec>            run each split for <sec> seconds\n";
  cerr << " -r <n>              report progress every <n> steps\n";
  cerr << " -T <dist>           discard neighbours <dist>\n";
  cerr << " -v                  verbose output\n";
// clang-format on
  
  ::exit(rc);
}

// Neighbour information consists of an ID and an integer distance.
struct Nbr {
  // This will be the index of this neighbour into the overall array of needles.
  uint32_t _id;
  // This value could be narrower, even uint8_t would be fine.
  uint32_t _dist;
};

// Each needle knows its current state, and info about its neighbours.
// Deliberately minimal information
class Needle {
  private:
    // The needle does not need to know its id since that is the same as its position
    // in the overall raay of needles.

    // Whether or not this item is in the training set.
    int _train;

    // the neighbours.
    uint32_t _number_nbrs;
    Nbr* _nbrs;

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

    uint32_t Neighbour(uint32_t ndx) const {
      return _nbrs[ndx]._id;
    }
};

Needle::Needle() {
  _train = 0;
  _number_nbrs = 0;
  _nbrs = nullptr;
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
      ok_nbrs = 1;
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

class Optimise {
  private:
    int _verbose;

    int _nsplits;

    uint32_t _number_needles;
    Needle* _needle;

    float _upper_distance_threshold;

    float _fraction_train;

    std::mt19937 _rng;

    std::unique_ptr<std::uniform_int_distribution<uint32_t>> _uniform;

    uint32_t _current_score;

    uint32_t _nopt;

    uint32_t _seconds;

    // For each pair of items where we know a distance, store that.
    // Key will be i * _number_needles + j where i > j.
    std::unordered_map<uint64_t, uint32_t> _distance;

    // One more than the max value of the values in _distance
    uint32_t _max_distance;

    // Rather than storing smiles and id with the Needle, we deliberately
    // make a global array of both, hoping to keep the Needle class very
    // small.
    std::string* _smiles;
    std::string* _name;

    IWString _stem;

    Report_Progress _report_progress;

  // Private functions.
    uint32_t GatherIdentifiers(const char* fname,
                            std::unordered_map<std::string, std::uint32_t>& id_to_ndx);
    uint32_t ReadNeighbours(const char* fname,
                            const std::unordered_map<std::string, uint32_t>& id_to_ndx);
    int RandomSplit();

    uint64_t FormKey(uint32_t i, uint32_t j) const;

    std::optional<uint32_t> Distance(uint32_t i, uint32_t j) const;
    uint32_t DistanceOrMax(uint32_t i, uint32_t j) const;

    uint32_t RecomputeCurrentScore();

    std::tuple<uint32_t, uint32_t> ChooseTwo();

    int MakeSplit(int split);

    int WriteSplit(int split) const;
    int WriteSmiles(int train, IWString& fname) const;
    int WriteSmiles(int train, IWString_and_File_Descriptor& output) const;
    int WriteSubset(int train, IWString& fname) const;
    int WriteSubset(int train, IWString_and_File_Descriptor& output) const;

  public:
    Optimise();
    ~Optimise();

    int Initialise(const Command_Line& cl);

    int ReportSize(std::ostream& output) const;

    int ReadFingerprints(const char* fname);

    int Doit();
};

Optimise::Optimise() {
  _verbose = 0;
  _nsplits = 10;
  _fraction_train = 0.80;
  _nopt = 1000;

  _upper_distance_threshold = std::numeric_limits<float>::max();

  _max_distance = 0;

  _number_needles = 0;
  _needle = nullptr;

  _smiles = nullptr;
  _name = nullptr;

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

  if (cl.option_present('o') && cl.option_present('s')) {
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

  if (cl.option_present('s')) {
    if (! cl.value('s', _seconds) || _seconds < 1) {
      cerr << "Optimise::Initialise:Invalid seconds (-s)\n";
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

  if (cl.option_present('S')) {
    cl.value('S', _stem);
    if (_verbose) {
      cerr << "Split files created with stem '" << _stem << "'\n";
    }
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

  output << "Distance matrix holds " << _distance.size() << " items, max " << _max_distance << '\n';

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

uint32_t
Optimise::RecomputeCurrentScore() {
  _current_score = 0;
  for (uint32_t i = 0; i < _number_needles; ++i) {
    const int itrain = _needle[i].in_train();
    for (uint32_t j = i + 1; j < _number_needles; ++j) {
      const int jtrain = _needle[j].in_train();
      if (itrain == jtrain) {
        continue;
      }

      const uint32_t d = DistanceOrMax(i, j);
      //cerr << " " << i << " and " << j << " opposite, dist " << d << '\n';

      _current_score += d;
    }
  }

  //cerr << "Current score " << _current_score << '\n';

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

  _uniform = std::make_unique<std::uniform_int_distribution<uint32_t>>(0, _number_needles - 1);

  uint32_t ndx;
  for (ndx = 0; ; ++ndx) {
    std::optional<nnbr::NearNeighbours> needle = reader.ReadProto<nnbr::NearNeighbours>();
    if (! needle) {
      break;
    }

    _smiles[ndx] = needle->smiles();
    _name[ndx] = needle->name();
    if (!_needle[ndx].Build(*needle, id_to_ndx, _upper_distance_threshold)) {
      cerr << "Optimise::ReadNeighbours:cannot process " << needle->ShortDebugString() << "\n";
      return 0;
    }

    for (const nnbr::Nbr& nbr : needle->nbr()) {
      const auto iter = id_to_ndx.find(nbr.id());
      uint64_t key = FormKey(ndx, iter->second);
      _distance[key] = DistanceToInt(nbr.dist());
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

  return _number_needles;
}

// Return true if we are now > `seconds` seconds from `tzero`.
bool
BreakForTime(std::time_t tzero, uint32_t seconds) {
  return (std::time(nullptr) - tzero) >= seconds;
}

int
Optimise::Doit() {
  for (int i = 0; i < _nsplits; ++i) {
    MakeSplit(i);
  }

  return 1;
}

// Starting with a random split, optimize it and write when done.
int
Optimise::MakeSplit(int split) {
  std::time_t tzero = 0;
  if (_seconds) {
    tzero = std::time(nullptr);
  }

  RandomSplit();
  uint64_t score = RecomputeCurrentScore();
  const uint64_t starting_score = score;
  if (_verbose) {
    cerr << "Split " << split << " starting score " << score << '\n';
  }

  uint32_t steps_accepted = 0;
  cerr << "Will perform " << _nopt << " optimisations\n";
  for (uint32_t j = 0; j < _nopt; ++j) {
    auto [i1, i2] = ChooseTwo();
    _needle[i1].invert_train();
    _needle[i2].invert_train();

    const int train1 = _needle[i1].in_train();
    const int train2 = _needle[i2].in_train();
    //cerr << "Selected " << i1 << " and " << i2 << " has " << _needle[i1].number_neighbours() << " nbrs\n";

    int32_t delta = 0;   // a signed quantity.
    uint32_t number_nbrs = _needle[i1].number_neighbours();
    for (uint32_t x = 0; x < number_nbrs; ++x) {
      uint32_t n = _needle[i1].Neighbour(x);
      if (n == i2) {
        continue;
      }
      uint32_t d = DistanceOrMax(i2, n);
      // cerr <<  "  adj " << i1 << " to " << n << " dist " << d << " train1 " << train1 << " nbr " << _needle[n].in_train() << '\n';

      // If both on the same side now, were previously on opposite sides.
      if (_needle[n].in_train() == train2) {
        delta -= d;
      } else  {  // different sides now, previously same.
        delta += d;
      }
      // cerr << "  delta updated to " << delta << '\n';
    }

    number_nbrs = _needle[i2].number_neighbours();
    for (uint32_t x = 0; x < number_nbrs; ++x) {
      uint32_t n = _needle[i2].Neighbour(x);
      if (n == i1) {
        continue;
      }
      uint32_t d = DistanceOrMax(i1, n);
      // cerr <<  "  adj " << i2 << " to " << n << " dist " << d << " train2 " << train2 << " nbr " << _needle[n].in_train() << '\n';

      // If both on the same side now, were previously on opposite sides.
      if (_needle[n].in_train() == train1) {
        delta -= d;
      } else  {  // different sides now, previously same.
        delta += d;
      }
      // cerr << "   2nd loop, delta " << delta << '\n';
    }
    // cerr << "Score " << score << " delta " << delta << '\n';
    uint64_t new_score = score + delta;
#ifdef DEBUG_SWAP_ITEMS
    cerr << "new_score " << new_score << " cmp " << score << '\n';
    for (uint32_t y = 0; y < _number_needles; ++y) {
      cerr << " needle " << y << " train " << _needle[y].in_train() << '\n';
    }
#endif

    if (_report_progress() && j > 0) {
      cerr << split << ' ' << j << " score " << new_score << " cmp " << starting_score << " accepted " << steps_accepted << " " << iwmisc::Fraction<float>(steps_accepted, j) << '\n';
    }

    // Heuristics to stop this being checked too often.
    if (_seconds > 0 && j > 500 && j % 1000 == 0 && BreakForTime(tzero, _seconds)) {
      break;
    }

    if (new_score == score) {
      continue;
    }
    // If better, sets more separated, always accept.
    if (new_score > score) {
      score = new_score;
      ++steps_accepted;
      continue;
    }
    // Worse, revert. Generate random number...
    if (new_score < score) {
      _needle[i1].invert_train();
      _needle[i2].invert_train();
      continue;
    }
  }

  cerr << "Writing split " << split << '\n';
  WriteSplit(split);

  cerr << "Split " << split << " accepted " << steps_accepted << " of " << _nopt << " steps\n";

  return 1;
}

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

  return 1;
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
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vf:S:n:o:r:T:s:t:");
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

  if (verbose) {
    optimise.ReportSize(cerr);
  }

  optimise.Doit();

  return 0;
}

}  // namespace train_test_split

int
main(int argc, char** argv) {
  int rc = train_test_split::Main(argc, argv);
  return rc;
}
