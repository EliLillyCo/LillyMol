#include <cstdint>
#include <cstdlib>

#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.h"

namespace gfp_leader_buckets {

using std::cerr;

constexpr int kDefaultMinScanForThreading = 10000;

struct DistanceHit {
  int pool_index;
  similarity_type_t distance;
};

class Candidate : public IW_General_Fingerprint {
 private:
  enum class State : uint8_t {
    kAvailable,
    kLeader,
    kClusterMember,
    kPreviouslySelected,
  };

  State _state = State::kAvailable;
  IWString _smiles;
  IWString _bucket_label;
  int _bucket = -1;
  similarity_type_t _distance = 0.0f;

 public:
  int ConstructFromTdt(IW_TDT& tdt, const IWString& smiles_tag, const IWString& id_tag,
                       int& fatal);

  int selected() const {
    return _state != State::kAvailable;
  }

  int available() const {
    return _state == State::kAvailable;
  }

  void MarkLeader() {
    _state = State::kLeader;
    _distance = 0.0f;
  }

  void MarkClusterMember(similarity_type_t d) {
    _state = State::kClusterMember;
    _distance = d;
  }

  void MarkPreviouslySelected(similarity_type_t d) {
    _state = State::kPreviouslySelected;
    _distance = d;
  }

  const IWString& smiles() const {
    return _smiles;
  }

  const IWString& bucket_label() const {
    return _bucket_label;
  }

  int bucket() const {
    return _bucket;
  }

  void set_bucket(int b, const IWString& label) {
    _bucket = b;
    _bucket_label = label;
  }

  similarity_type_t distance() const {
    return _distance;
  }
};

int
Candidate::ConstructFromTdt(IW_TDT& tdt, const IWString& smiles_tag,
                            const IWString& id_tag, int& fatal) {
  if (! IW_General_Fingerprint::construct_from_tdt(tdt, fatal)) {
    return 0;
  }

  const_IWSubstring value;
  if (! tdt.dataitem_value(smiles_tag, value)) {
    cerr << "Candidate::ConstructFromTdt:missing smiles tag '" << smiles_tag << "'\n";
    fatal = 1;
    return 0;
  }
  _smiles = value;

  if (id().empty()) {
    cerr << "Candidate::ConstructFromTdt:missing identifier tag '" << id_tag << "'\n";
    fatal = 1;
    return 0;
  }

  return 1;
}

class Bucket {
 private:
  IWString _label;
  int _total = 0;
  int _leaders = 0;
  int _cluster_members = 0;
  int _previously_selected = 0;
  int _available = 0;
  std::vector<int> _members;
  uint32_t _cursor = 0;

 public:
  explicit Bucket(const IWString& label) : _label(label) {}

  const IWString& label() const {
    return _label;
  }

  int total() const {
    return _total;
  }

  int leaders() const {
    return _leaders;
  }

  int cluster_members() const {
    return _cluster_members;
  }

  int previously_selected() const {
    return _previously_selected;
  }

  int available() const {
    return _available;
  }

  void AddMember(int ndx) {
    _members.push_back(ndx);
    ++_total;
    ++_available;
  }

  void DecrementAvailable() {
    --_available;
  }

  void LeaderSelected() {
    ++_leaders;
  }

  void ClusterMemberSelected() {
    ++_cluster_members;
  }

  void PreviouslySelected() {
    ++_previously_selected;
  }

  int NextAvailable(const resizable_array_p<Candidate>& pool);
};

int
Bucket::NextAvailable(const resizable_array_p<Candidate>& pool) {
  while (_cursor < _members.size()) {
    const int ndx = _members[_cursor];
    if (pool[ndx]->available()) {
      return ndx;
    }
    ++_cursor;
  }

  return -1;
}

class BucketLeaderOptions {
 private:
  int _verbose = 0;
  IWString _smiles_tag = "$SMI<";
  IWString _id_tag = "PCN<";
  IWString _distance_tag = "DIST<";
  IWString _bucket_output_tag = "BCKT<";

  similarity_type_t _threshold = -1.0f;
  similarity_type_t _previously_selected_threshold = -1.0f;
  int _previously_selected_threshold_specified = 0;
  int _max_leaders = std::numeric_limits<int>::max();
  int _nthreads = 1;
  int _report_progress = 0;
  int _min_scan_for_threading = kDefaultMinScanForThreading;
  IWString _bucket_summary_fname;

  int _bucket_column = -1;
  IWString _bucket_file;
  std::unordered_map<std::string, IWString> _bucket_from_id;

 public:
  int Initialise(Command_Line& cl);

  int verbose() const {
    return _verbose;
  }

  const IWString& smiles_tag() const {
    return _smiles_tag;
  }

  const IWString& id_tag() const {
    return _id_tag;
  }

  const IWString& distance_tag() const {
    return _distance_tag;
  }

  const IWString& bucket_output_tag() const {
    return _bucket_output_tag;
  }

  similarity_type_t threshold() const {
    return _threshold;
  }

  similarity_type_t previously_selected_threshold() const {
    return _previously_selected_threshold_specified ? _previously_selected_threshold : _threshold;
  }

  int max_leaders() const {
    return _max_leaders;
  }

  int nthreads() const {
    return _nthreads;
  }

  int report_progress() const {
    return _report_progress;
  }

  int min_scan_for_threading() const {
    return _min_scan_for_threading;
  }

  const IWString& bucket_summary_fname() const {
    return _bucket_summary_fname;
  }

  int AssignBucket(Candidate& candidate) const;

 private:
  int ReadBucketFile(const IWString& fname);
  int InitialiseMiscellaneousOptions(Command_Line& cl);
};

void
DisplayMiscellaneousOptions(std::ostream& output) {
  output << " -D minscan=<n>  run serial scans unless at least <n> candidates remain\n";
  output << "                 default " << kDefaultMinScanForThreading << '\n';
}

int
BucketLeaderOptions::ReadBucketFile(const IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open bucket file '" << fname << "'\n";
    return 0;
  }

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    buffer.strip_leading_blanks();
    buffer.strip_trailing_blanks();
    if (buffer.empty() || buffer.starts_with('#')) {
      continue;
    }

    if (buffer.nwords() != 2) {
      cerr << "Invalid bucket file record line " << input.lines_read() << " '" << buffer << "'\n";
      return 0;
    }

    const_IWSubstring id, bucket;
    buffer.word(0, id);
    buffer.word(1, bucket);
    _bucket_from_id.emplace(std::string(id.data(), id.length()), IWString(bucket));
  }

  if (_bucket_from_id.empty()) {
    cerr << "No bucket assignments read from '" << fname << "'\n";
    return 0;
  }

  if (_verbose) {
    cerr << "Read " << _bucket_from_id.size() << " bucket assignments from '" << fname << "'\n";
  }

  return 1;
}

int
BucketLeaderOptions::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('I')) {
    cl.value('I', _id_tag);
    if (! _id_tag.ends_with('<')) {
      _id_tag << '<';
    }
    set_identifier_tag(_id_tag);
  }

  if (cl.option_present('B')) {
    if (cl.option_count('B') > 1) {
      cerr << "Only one -B bucket source is supported\n";
      return 0;
    }

    const_IWSubstring b = cl.string_value('B');
    if (b.starts_with("col=")) {
      b.remove_leading_chars(4);
      if (! b.numeric_value(_bucket_column) || _bucket_column < 1) {
        cerr << "Invalid bucket column '" << b << "'\n";
        return 0;
      }
      --_bucket_column;
      if (_verbose) {
        cerr << "Bucket labels in name column " << (_bucket_column + 1) << '\n';
      }
    } else if (b.starts_with("file=")) {
      _bucket_file = b;
      _bucket_file.remove_leading_chars(5);
      if (! ReadBucketFile(_bucket_file)) {
        return 0;
      }
    } else {
      _bucket_file = b;
      if (! ReadBucketFile(_bucket_file)) {
        return 0;
      }
    }
  } else {
    cerr << "Must specify a bucket source via -B col=<n> or -B <file>\n";
    return 0;
  }

  if (cl.option_present('t')) {
    if (! cl.value('t', _threshold) || _threshold < 0.0f || _threshold > 1.0f) {
      cerr << "The sphere exclusion threshold (-t) must be in [0,1]\n";
      return 0;
    }
  } else {
    cerr << "Must specify sphere exclusion threshold via -t\n";
    return 0;
  }

  if (cl.option_present('a')) {
    if (! cl.value('a', _previously_selected_threshold) ||
        _previously_selected_threshold < 0.0f || _previously_selected_threshold > 1.0f) {
      cerr << "The previously selected threshold (-a) must be in [0,1]\n";
      return 0;
    }
    _previously_selected_threshold_specified = 1;
  }

  if (cl.option_present('n')) {
    if (! cl.value('n', _max_leaders) || _max_leaders < 1) {
      cerr << "The max leaders option (-n) must be a positive integer\n";
      return 0;
    }
  }

  if (cl.option_present('j')) {
    if (! cl.value('j', _nthreads) || _nthreads < 1) {
      cerr << "The thread count option (-j) must be a positive integer\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will use " << _nthreads << " threads for distance scans\n";
    }
  }

  if (cl.option_present('r')) {
    if (! cl.value('r', _report_progress) || _report_progress < 1) {
      cerr << "The report progress option (-r) must be a positive integer\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will report progress every " << _report_progress << " leaders\n";
    }
  }

  if (cl.option_present('Y')) {
    _bucket_summary_fname = cl.string_value('Y');
    if (_verbose) {
      cerr << "Will write bucket summary to '" << _bucket_summary_fname << "'\n";
    }
  }

  if (cl.option_present('D') && ! InitialiseMiscellaneousOptions(cl)) {
    return 0;
  }

  return 1;
}

int
BucketLeaderOptions::InitialiseMiscellaneousOptions(Command_Line& cl) {
  const_IWSubstring d;
  for (int i = 0; cl.value('D', d, i); ++i) {
    if (d == "help") {
      DisplayMiscellaneousOptions(cerr);
      ::exit(0);
    }

    if (d.starts_with("minscan=")) {
      d.remove_leading_chars(8);
      if (! d.numeric_value(_min_scan_for_threading) || _min_scan_for_threading < 0) {
        cerr << "Invalid -D minscan value '" << d << "'\n";
        return 0;
      }
      if (_verbose) {
        cerr << "Will use serial distance scans below " << _min_scan_for_threading
             << " candidates\n";
      }
    } else {
      cerr << "Unrecognised -D qualifier '" << d << "'\n";
      DisplayMiscellaneousOptions(cerr);
      return 0;
    }
  }

  return 1;
}

int
BucketLeaderOptions::AssignBucket(Candidate& candidate) const {
  IWString label;

  if (_bucket_column >= 0) {
    const IWString& id = candidate.id();
    const_IWSubstring token;
    if (! id.word(_bucket_column, token)) {
      cerr << "Cannot extract bucket column " << (_bucket_column + 1) << " from '"
           << id << "'\n";
      return 0;
    }
    label = token;
  } else {
    auto iter = _bucket_from_id.find(std::string(candidate.id().data(), candidate.id().length()));
    if (iter == _bucket_from_id.end()) {
      cerr << "No bucket assignment for '" << candidate.id() << "'\n";
      return 0;
    }
    label = iter->second;
  }

  // bucket number assigned later, once labels have been compacted.
  candidate.set_bucket(-1, label);
  return 1;
}

class BucketLeader {
 private:
  BucketLeaderOptions _options;
  resizable_array_p<Candidate> _pool;
  resizable_array_p<Bucket> _buckets;
  std::unordered_map<std::string, int> _bucket_number;
  int _first_unselected = 0;
  int _leaders_found = 0;
  int _discarded_by_previously_selected = 0;

 public:
  int Initialise(Command_Line& cl) {
    return _options.Initialise(cl);
  }

  int BuildPool(const char* fname);
  int ApplyPreviouslySelected(const char* fname);
  int Run(IWString_and_File_Descriptor& output);
  int WriteBucketSummary() const;
  int Report(std::ostream& output) const;

 private:
  int AvailableCount() const;
  void MaybeReportProgress() const;
  int BuildPool(iwstring_data_source& input);
  int AddCandidate(Candidate* candidate);
  int BucketNumber(const IWString& label);
  int ChooseNextLeader() const;
  int FormCluster(int leader, std::vector<int>& cluster_members);
  int MarkLeader(int ndx);
  int MarkClusterMember(int ndx, similarity_type_t d);
  int MarkPreviouslySelected(int ndx, similarity_type_t d);
  void AdvanceFirstUnselected();
  int FindWithinThreshold(IW_General_Fingerprint& reference, similarity_type_t threshold,
                          std::vector<DistanceHit>& hits) const;
  void FindWithinThresholdRange(IW_General_Fingerprint& reference,
                                similarity_type_t threshold, int begin, int end,
                                std::vector<DistanceHit>& hits) const;
  int WriteCluster(int leader, const std::vector<int>& cluster_members,
                   IWString_and_File_Descriptor& output) const;
};

int
BucketLeader::BucketNumber(const IWString& label) {
  const std::string key(label.data(), label.length());
  auto iter = _bucket_number.find(key);
  if (iter != _bucket_number.end()) {
    return iter->second;
  }

  const int result = _buckets.number_elements();
  _bucket_number.emplace(key, result);
  _buckets.add(new Bucket(label));
  return result;
}

int
BucketLeader::AddCandidate(Candidate* candidate) {
  if (! _options.AssignBucket(*candidate)) {
    return 0;
  }

  const int b = BucketNumber(candidate->bucket_label());
  candidate->set_bucket(b, candidate->bucket_label());

  const int ndx = _pool.number_elements();
  _pool.add(candidate);
  _buckets[b]->AddMember(ndx);
  return 1;
}

int
BucketLeader::BuildPool(iwstring_data_source& input) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    Candidate* candidate = new Candidate;
    int fatal = 0;
    if (! candidate->ConstructFromTdt(tdt, _options.smiles_tag(), _options.id_tag(), fatal)) {
      delete candidate;
      if (fatal) {
        cerr << "Cannot parse TDT\n" << tdt;
        return 0;
      }
      continue;
    }

    if (! AddCandidate(candidate)) {
      delete candidate;
      cerr << "Cannot assign bucket to TDT\n" << tdt;
      return 0;
    }
  }

  if (_pool.empty()) {
    cerr << "No fingerprints read\n";
    return 0;
  }

  if (_options.verbose()) {
    cerr << "Read " << _pool.number_elements() << " fingerprints in " << _buckets.number_elements()
         << " buckets\n";
  }

  return 1;
}

int
BucketLeader::BuildPool(const char* fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return BuildPool(input);
}

int
BucketLeader::MarkLeader(int ndx) {
  Candidate& candidate = *_pool[ndx];
  if (! candidate.available()) {
    return 0;
  }

  candidate.MarkLeader();
  Bucket& bucket = *_buckets[candidate.bucket()];
  bucket.DecrementAvailable();
  bucket.LeaderSelected();
  AdvanceFirstUnselected();
  return 1;
}

int
BucketLeader::MarkClusterMember(int ndx, similarity_type_t d) {
  Candidate& candidate = *_pool[ndx];
  if (! candidate.available()) {
    return 0;
  }

  candidate.MarkClusterMember(d);
  Bucket& bucket = *_buckets[candidate.bucket()];
  bucket.DecrementAvailable();
  bucket.ClusterMemberSelected();
  return 1;
}

int
BucketLeader::MarkPreviouslySelected(int ndx, similarity_type_t d) {
  Candidate& candidate = *_pool[ndx];
  if (! candidate.available()) {
    return 0;
  }

  candidate.MarkPreviouslySelected(d);
  Bucket& bucket = *_buckets[candidate.bucket()];
  bucket.DecrementAvailable();
  bucket.PreviouslySelected();
  ++_discarded_by_previously_selected;
  return 1;
}

void
BucketLeader::AdvanceFirstUnselected() {
  while (_first_unselected < _pool.number_elements() &&
         _pool[_first_unselected]->selected()) {
    ++_first_unselected;
  }
}

int
BucketLeader::ChooseNextLeader() const {
  int best_bucket = -1;
  for (int i = 0; i < _buckets.number_elements(); ++i) {
    const Bucket& bucket = *_buckets[i];
    if (bucket.available() == 0) {
      continue;
    }

    if (best_bucket < 0) {
      best_bucket = i;
      continue;
    }

    const Bucket& best = *_buckets[best_bucket];
    if (bucket.leaders() < best.leaders() ||
        (bucket.leaders() == best.leaders() && bucket.available() < best.available())) {
      best_bucket = i;
    }
  }

  if (best_bucket < 0) {
    return -1;
  }

  return _buckets[best_bucket]->NextAvailable(_pool);
}

void
BucketLeader::FindWithinThresholdRange(IW_General_Fingerprint& reference,
                                       similarity_type_t threshold, int begin, int end,
                                       std::vector<DistanceHit>& hits) const {
  hits.clear();
  for (int i = begin; i < end; ++i) {
    if (_pool[i]->selected()) {
      continue;
    }
    if (! can_be_compared(reference, *_pool[i])) {
      continue;
    }

    const similarity_type_t d = reference.IW_General_Fingerprint::distance(*_pool[i]);
    if (d <= threshold) {
      hits.push_back(DistanceHit{i, d});
    }
  }
}

int
BucketLeader::FindWithinThreshold(IW_General_Fingerprint& reference,
                                  similarity_type_t threshold,
                                  std::vector<DistanceHit>& hits) const {
  hits.clear();

  const int pool_size = _pool.number_elements();
  if (_first_unselected >= pool_size) {
    return 1;
  }

  const int items_to_scan = pool_size - _first_unselected;
  const int nthreads = std::min(_options.nthreads(), items_to_scan);
  if (nthreads <= 1 || items_to_scan < _options.min_scan_for_threading()) {
    FindWithinThresholdRange(reference, threshold, _first_unselected, pool_size, hits);
    return 1;
  }

  std::vector<std::vector<DistanceHit>> per_thread_hits(nthreads);
  std::vector<std::thread> threads;
  threads.reserve(nthreads - 1);

  for (int t = 1; t < nthreads; ++t) {
    const int begin = _first_unselected + (items_to_scan * t) / nthreads;
    const int end = _first_unselected + (items_to_scan * (t + 1)) / nthreads;

    threads.emplace_back([this, &reference, threshold, begin, end, &per_thread_hits, t]() {
      // Each worker gets a private fingerprint copy and a private hit vector.
      // Candidate and bucket state are read-only until all threads have joined.
      IW_General_Fingerprint local_reference(reference);
      FindWithinThresholdRange(local_reference, threshold, begin, end, per_thread_hits[t]);
    });
  }

  // The caller handles one slice. This avoids creating a worker just so the main
  // thread can immediately wait, which is especially costly with -j 2.
  const int begin = _first_unselected;
  const int end = _first_unselected + items_to_scan / nthreads;
  IW_General_Fingerprint local_reference(reference);
  FindWithinThresholdRange(local_reference, threshold, begin, end, per_thread_hits[0]);

  for (std::thread& thread : threads) {
    thread.join();
  }

  for (std::vector<DistanceHit>& thread_hits : per_thread_hits) {
    hits.insert(hits.end(), thread_hits.begin(), thread_hits.end());
  }

  return 1;
}

int
BucketLeader::FormCluster(int leader, std::vector<int>& cluster_members) {
  cluster_members.clear();

  std::vector<DistanceHit> hits;
  if (! FindWithinThreshold(*_pool[leader], _options.threshold(), hits)) {
    return 0;
  }

  for (const DistanceHit& hit : hits) {
    if (hit.pool_index == leader) {
      continue;
    }
    if (MarkClusterMember(hit.pool_index, hit.distance)) {
      cluster_members.push_back(hit.pool_index);
    }
  }

  AdvanceFirstUnselected();
  return 1;
}

int
BucketLeader::WriteCluster(int leader, const std::vector<int>& cluster_members,
                           IWString_and_File_Descriptor& output) const {
  const Candidate& centre = *_pool[leader];
  output << _options.smiles_tag() << centre.smiles() << ">\n";
  output << _options.id_tag() << centre.id() << ">\n";
  output << "CLUSTER<" << _leaders_found << ">\n";
  output << "CSIZE<" << (cluster_members.size() + 1) << ">\n";
  output << _options.bucket_output_tag() << centre.bucket_label() << ">\n";

  for (int ndx : cluster_members) {
    const Candidate& member = *_pool[ndx];
    output << _options.smiles_tag() << member.smiles() << ">\n";
    output << _options.id_tag() << member.id() << ">\n";
    output << _options.distance_tag() << member.distance() << ">\n";
    output << _options.bucket_output_tag() << member.bucket_label() << ">\n";
  }

  output << "|\n";
  output.write_if_buffer_holds_more_than(32768);
  return output.good();
}

int
BucketLeader::ApplyPreviouslySelected(const char* fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open previously selected file '" << fname << "'\n";
    return 0;
  }

  const similarity_type_t threshold = _options.previously_selected_threshold();
  IW_TDT tdt;
  while (tdt.next(input)) {
    IW_General_Fingerprint fp;
    int fatal = 0;
    if (! fp.construct_from_tdt(tdt, fatal)) {
      if (fatal) {
        cerr << "Cannot parse previously selected TDT\n" << tdt;
        return 0;
      }
      continue;
    }

    std::vector<DistanceHit> hits;
    if (! FindWithinThreshold(fp, threshold, hits)) {
      return 0;
    }
    for (const DistanceHit& hit : hits) {
      MarkPreviouslySelected(hit.pool_index, hit.distance);
    }
    AdvanceFirstUnselected();
  }

  if (_options.verbose()) {
    cerr << "Discarded " << _discarded_by_previously_selected
         << " fingerprints by previously selected file(s)\n";
  }

  return _first_unselected < _pool.number_elements();
}

int
BucketLeader::AvailableCount() const {
  int result = 0;
  for (const Bucket* bucket : _buckets) {
    result += bucket->available();
  }

  return result;
}

void
BucketLeader::MaybeReportProgress() const {
  const int report_progress = _options.report_progress();
  if (report_progress == 0 || (_leaders_found % report_progress) != 0) {
    return;
  }

  cerr << "gfp_leader_buckets:formed " << _leaders_found << " leaders, "
       << AvailableCount() << " candidates remain\n";
}

int
BucketLeader::Run(IWString_and_File_Descriptor& output) {
  std::vector<int> cluster_members;

  while (_leaders_found < _options.max_leaders()) {
    const int leader = ChooseNextLeader();
    if (leader < 0) {
      break;
    }

    if (! MarkLeader(leader)) {
      cerr << "Internal error, cannot mark leader " << leader << "\n";
      return 0;
    }

    if (! FormCluster(leader, cluster_members)) {
      return 0;
    }

    if (! WriteCluster(leader, cluster_members, output)) {
      return 0;
    }

    ++_leaders_found;
    MaybeReportProgress();
  }

  return 1;
}


int
BucketLeader::WriteBucketSummary() const {
  const IWString& fname = _options.bucket_summary_fname();
  if (fname.empty()) {
    return 1;
  }

  IWString tmp(fname);
  IWString_and_File_Descriptor output;
  if (! output.open(tmp.null_terminated_chars())) {
    cerr << "Cannot open bucket summary file '" << fname << "'\n";
    return 0;
  }

  const char sep = fname.ends_with(".csv") ? ',' : ' ';
  output << "bucket" << sep << "total" << sep << "leaders" << sep
         << "cluster_members" << sep << "previously_selected" << sep
         << "available\n";

  for (const Bucket* bucket : _buckets) {
    output << bucket->label() << sep << bucket->total() << sep << bucket->leaders()
           << sep << bucket->cluster_members() << sep << bucket->previously_selected()
           << sep << bucket->available() << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }

  output.flush();
  return output.good();
}

int
BucketLeader::Report(std::ostream& output) const {
  output << "Found " << _leaders_found << " leaders\n";
  if (_discarded_by_previously_selected) {
    output << "Discarded " << _discarded_by_previously_selected
           << " fingerprints by previously selected file(s)\n";
  }

  for (const Bucket* bucket : _buckets) {
    output << "Bucket '" << bucket->label() << "' total " << bucket->total()
           << " leaders " << bucket->leaders()
           << " cluster_members " << bucket->cluster_members()
           << " previously_selected " << bucket->previously_selected()
           << " available " << bucket->available() << '\n';
  }

  return output.good();
}

void
Usage(int rc) {
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif

  cerr << R"(Performs bucket-balanced leader sphere exclusion clustering on fingerprints.

Usage: gfp_leader_buckets -t <dist> -B col=<n>|-B <file> [options] input.gfp

 -t <dist>      sphere exclusion distance threshold.
 -B col=<n>     bucket label is token <n> in the identifier/name field.
 -B <file>      two column id bucket mapping file; file=<fname> also accepted.
 -n <n>         maximum number of leaders to write.
 -j <n>         number of threads for distance scans, default 1.
 -r <n>         report progress every <n> leaders formed.
 -Y <fname>     write bucket summary; .csv suffix gives comma-separated output.
 -A <file>      previously selected fingerprints; candidates within threshold are discarded.
 -a <dist>      distance threshold for -A files; default is -t value.
 -I <tag>       identifier tag, default PCN<.
 -D ...         miscellaneous options, enter '-D help' for details.
 -F -P -W       standard GFP options, enter '-F help' for details.
 -v             verbose output.
)";

  ::exit(rc);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vF:P:W:t:B:n:j:r:Y:A:a:I:D:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (cl.number_elements() != 1) {
    cerr << "gfp_leader_buckets expects exactly one input file\n";
    Usage(1);
  }

  if (need_to_call_initialise_fingerprints(cl)) {
    if (! initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise fingerprint options\n";
      Usage(1);
    }
  }

  BucketLeader bucket_leader;
  if (! bucket_leader.Initialise(cl)) {
    cerr << "Cannot initialise bucket leader options\n";
    Usage(1);
  }

  if (! bucket_leader.BuildPool(cl[0])) {
    cerr << "Cannot build pool from '" << cl[0] << "'\n";
    return 1;
  }

  IWString fname;
  for (int i = 0; cl.value('A', fname, i); ++i) {
    if (! bucket_leader.ApplyPreviouslySelected(fname.null_terminated_chars())) {
      cerr << "Cannot process previously selected file '" << fname << "'\n";
      return 1;
    }
  }

  IWString_and_File_Descriptor output(1);
  if (! bucket_leader.Run(output)) {
    cerr << "Fatal error during clustering\n";
    return 1;
  }

  output.flush();

  if (! bucket_leader.WriteBucketSummary()) {
    return 1;
  }

  if (verbose) {
    bucket_leader.Report(cerr);
  }

  return 0;
}

}  // namespace gfp_leader_buckets

int
main(int argc, char** argv) {
  return gfp_leader_buckets::Main(argc, argv);
}
