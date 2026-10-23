#include "Utilities/GFP_Tools/train_test_split_activity.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <utility>

#include "absl/container/flat_hash_map.h"

#include "Foundational/iwmisc/iw_tabular_data.h"
#include "Foundational/iwstring/absl_hash.h"

namespace train_test_split {

using std::cerr;

namespace {

char
SeparatorForFileName(const IWString& fname) {
  if (fname.ends_with(".csv")) {
    return ',';
  }

  return ' ';
}

int
ParseBucketisation(const const_IWSubstring& value, ActivityBucketisation& bucketisation) {
  if (value == "quantile" || value == "quantiles" || value == "q") {
    bucketisation = ActivityBucketisation::kQuantile;
    return 1;
  }
  if (value == "width" || value == "equal_width" || value == "range" || value == "w") {
    bucketisation = ActivityBucketisation::kEqualWidth;
    return 1;
  }

  cerr << "Unrecognised activity bucketisation '" << value << "'\n";
  return 0;
}

int
ApplyDirective(const const_IWSubstring& token, ActivityProfileOptions& options) {
  if (token.starts_with("buckets=")) {
    const_IWSubstring value(token);
    value.remove_leading_chars(8);
    if (! value.numeric_value(options.buckets) || options.buckets == 0 ||
        options.buckets > std::numeric_limits<uint8_t>::max()) {
      cerr << "Invalid activity bucket count '" << token << "'\n";
      return 0;
    }
    return 1;
  }

  if (token.starts_with("weight=") || token.starts_with("objective_weight=")) {
    const_IWSubstring value(token);
    if (token.starts_with("weight=")) {
      value.remove_leading_chars(7);
    } else {
      value.remove_leading_chars(17);
    }
    if (! value.numeric_value(options.objective_weight) || options.objective_weight == 0) {
      cerr << "Invalid activity profile objective weight '" << token << "'\n";
      return 0;
    }
    return 1;
  }
  if (token.starts_with("tolerance=")) {
    const_IWSubstring value(token);
    value.remove_leading_chars(10);
    if (! value.numeric_value(options.tolerance) || options.tolerance > 100) {
      cerr << "Invalid activity tolerance '" << token << "'\n";
      return 0;
    }
    return 1;
  }

  if (token.starts_with("tol=")) {
    const_IWSubstring value(token);
    value.remove_leading_chars(4);
    if (! value.numeric_value(options.tolerance) || options.tolerance > 100) {
      cerr << "Invalid activity tolerance '" << token << "'\n";
      return 0;
    }
    return 1;
  }

  if (token.starts_with("bucket=")) {
    const_IWSubstring value(token);
    value.remove_leading_chars(7);
    return ParseBucketisation(value, options.bucketisation);
  }

  if (token.starts_with("bucketisation=")) {
    const_IWSubstring value(token);
    value.remove_leading_chars(14);
    return ParseBucketisation(value, options.bucketisation);
  }

  return ParseBucketisation(token, options.bucketisation);
}

uint64_t
AbsDiff(uint64_t v1, uint64_t v2) {
  return v1 > v2 ? v1 - v2 : v2 - v1;
}

}  // namespace

int64_t
CombinedActivityDistanceDelta(int64_t distance_delta, int64_t activity_delta,
                              uint32_t activity_objective_weight) {
  return distance_delta + static_cast<int64_t>(activity_objective_weight) * activity_delta;
}

int
ActivityProfileSet::AddDirective(const const_IWSubstring& directive) {
  resizable_array_p<const_IWSubstring> tokens;
  if (directive.split(tokens, ':') == 0) {
    cerr << "ActivityProfileSet::AddDirective:empty directive\n";
    return 0;
  }

  ActivityProfileOptions options = _defaults;
  const_IWSubstring first(*tokens[0]);

  if (first.contains('=')) {
    for (int i = 0; i < tokens.number_elements(); ++i) {
      if (! ApplyDirective(*tokens[i], _defaults)) {
        return 0;
      }
    }
    return 1;
  }

  options.fname = first;
  for (int i = 1; i < tokens.number_elements(); ++i) {
    if (! ApplyDirective(*tokens[i], options)) {
      return 0;
    }
  }

  if (options.fname.empty()) {
    cerr << "ActivityProfileSet::AddDirective:no file name in '" << directive << "'\n";
    return 0;
  }

  _pending.push_back(std::move(options));
  return 1;
}

int
ActivityProfileSet::Build(const std::vector<std::string>& id,
                          const std::vector<uint32_t>& sample_weight,
                          int verbose) {
  _profile.reserve(_pending.size());
  for (const ActivityProfileOptions& options : _pending) {
    ActivityProfile profile;
    if (! profile.Build(options, id, sample_weight, verbose)) {
      return 0;
    }
    _profile.push_back(std::move(profile));
  }

  return 1;
}

int
ActivityProfile::Build(const ActivityProfileOptions& options,
                       const std::vector<std::string>& id,
                       const std::vector<uint32_t>& sample_weight,
                       int verbose) {
  if (id.size() != sample_weight.size()) {
    cerr << "ActivityProfile::Build:size mismatch " << id.size() << " vs " << sample_weight.size() << '\n';
    return 0;
  }

  _options = options;
  _sample_weight = &sample_weight;

  IW_Tabular_Data<double> reader;
  reader.set_has_header(1);
  reader.set_first_column_is_identifier(1);

  const char sep = SeparatorForFileName(_options.fname);
  if (! reader.build(_options.fname.null_terminated_chars(), sep)) {
    cerr << "ActivityProfile::Build:cannot read activity data from '" << _options.fname << "'\n";
    return 0;
  }
  if (reader.ncols() == 0) {
    cerr << "ActivityProfile::Build:no activity columns in '" << _options.fname << "'\n";
    return 0;
  }

  absl::flat_hash_map<IWString, uint32_t> id_to_ndx;
  id_to_ndx.reserve(id.size());
  for (uint32_t i = 0; i < id.size(); ++i) {
    id_to_ndx[IWString(id[i])] = i;
  }
  if (id_to_ndx.size() != id.size()) {
    cerr << "ActivityProfile::Build:duplicate molecule identifiers supplied\n";
    return 0;
  }

  std::vector<double> activity(id.size(), 0.0);
  std::vector<uint8_t> seen(id.size(), 0);
  absl::flat_hash_map<IWString, int> activity_ids_seen;
  activity_ids_seen.reserve(reader.nrows());
  uint32_t matched = 0;
  _activity_records_not_used = 0;

  for (int i = 0; i < reader.nrows(); ++i) {
    const IWString& activity_id = reader.ids()[i];
    if (activity_ids_seen.contains(activity_id)) {
      cerr << "ActivityProfile::Build:duplicate activity for '" << activity_id << "'\n";
      return 0;
    }
    activity_ids_seen[activity_id] = 1;

    const auto iter = id_to_ndx.find(activity_id);
    if (iter == id_to_ndx.end()) {
      ++_activity_records_not_used;
      continue;
    }

    const uint32_t ndx = iter->second;
    seen[ndx] = 1;
    activity[ndx] = reader.value(i, 0);
    ++matched;
  }

  _observed = matched;
  _missing = id.size() - matched;
  if (_observed == 0) {
    cerr << "ActivityProfile::Build:no activity values in '" << _options.fname <<
            "' matched the molecule set\n";
    return 0;
  }

  _scaled.assign(id.size(), kMissingScaled);
  double minval = std::numeric_limits<double>::max();
  double maxval = -std::numeric_limits<double>::max();
  for (uint32_t i = 0; i < activity.size(); ++i) {
    if (! seen[i]) {
      continue;
    }
    minval = std::min(minval, activity[i]);
    maxval = std::max(maxval, activity[i]);
  }

  if (minval == maxval) {
    cerr << "ActivityProfile::Build:all activity values in '" << _options.fname <<
            "' are " << minval << ", no bucketisation possible\n";
    return 0;
  }

  const double range = maxval - minval;
  for (uint32_t i = 0; i < activity.size(); ++i) {
    if (! seen[i]) {
      continue;
    }
    int scaled = static_cast<int>(std::lround(100.0 * (activity[i] - minval) / range));
    scaled = std::clamp(scaled, 0, 100);
    _scaled[i] = static_cast<uint8_t>(scaled);
  }

  if (verbose) {
    cerr << "ActivityProfile::Build '" << _options.fname << "' observed " << _observed <<
            " missing " << _missing << " unused " << _activity_records_not_used <<
            " min " << minval << " max " << maxval << " buckets " << _options.buckets <<
            " weight " << _options.objective_weight << " tolerance " << _options.tolerance << "%\n";
  }

  return AssignBuckets(activity, seen, verbose);
}

int
ActivityProfile::AssignBuckets(const std::vector<double>& activity,
                               const std::vector<uint8_t>& seen,
                               int verbose) {
  _bucket.assign(activity.size(), kMissingBucket);
  _total_bucket_count.assign(_options.buckets, 0);
  _train_bucket_count.assign(_options.buckets, 0);
  _ideal_bucket_count.assign(_options.buckets, 0);
  _tolerance.assign(_options.buckets, 0);

  if (_options.bucketisation == ActivityBucketisation::kEqualWidth) {
    if (! AssignEqualWidthBuckets(seen)) {
      return 0;
    }
  } else if (! AssignQuantileBuckets(activity, seen)) {
    return 0;
  }

  return CheckBucketPopulation(verbose);
}

int
ActivityProfile::AssignEqualWidthBuckets(const std::vector<uint8_t>& seen) {
  for (uint32_t i = 0; i < _scaled.size(); ++i) {
    if (! seen[i]) {
      continue;
    }
    const uint32_t bucket = std::min<uint32_t>(_options.buckets - 1,
                            (static_cast<uint32_t>(_scaled[i]) * _options.buckets) / 101);
    _bucket[i] = bucket;
    ++_total_bucket_count[bucket];
  }

  return 1;
}

int
ActivityProfile::AssignQuantileBuckets(const std::vector<double>& activity,
                                        const std::vector<uint8_t>& seen) {
  std::vector<std::pair<double, uint32_t>> sorted;
  sorted.reserve(_observed);
  for (uint32_t i = 0; i < activity.size(); ++i) {
    if (seen[i]) {
      sorted.emplace_back(activity[i], i);
    }
  }

  std::sort(sorted.begin(), sorted.end());

  uint32_t begin = 0;
  while (begin < sorted.size()) {
    uint32_t end = begin + 1;
    while (end < sorted.size() && sorted[end].first == sorted[begin].first) {
      ++end;
    }

    const uint32_t bucket = std::min<uint32_t>(_options.buckets - 1,
                                (begin * _options.buckets) / sorted.size());
    for (uint32_t i = begin; i < end; ++i) {
      _bucket[sorted[i].second] = bucket;
      ++_total_bucket_count[bucket];
    }
    begin = end;
  }

  return 1;
}

int
ActivityProfile::CheckBucketPopulation(int verbose) const {
  uint32_t populated = 0;
  uint32_t empty = 0;
  for (uint32_t i = 0; i < _options.buckets; ++i) {
    if (_total_bucket_count[i] == 0) {
      ++empty;
    } else {
      ++populated;
    }
  }

  if (verbose) {
    cerr << "ActivityProfile::Build bucket counts for '" << _options.fname << "'";
    for (uint32_t i = 0; i < _options.buckets; ++i) {
      cerr << ' ' << i << ':' << _total_bucket_count[i];
    }
    cerr << '\n';
  }

  if (empty > 0) {
    cerr << "ActivityProfile::Build:warning '" << _options.fname << "' has " << empty <<
            " empty buckets\n";
  }

  if (populated <= 1) {
    cerr << "ActivityProfile::Build:only " << populated << " populated bucket in '" <<
            _options.fname << "'\n";
    return 0;
  }

  if (_observed < 2 * _options.buckets) {
    cerr << "ActivityProfile::Build:warning '" << _options.fname << "' has only " <<
            _observed << " observations for " << _options.buckets << " buckets\n";
  }

  return 1;
}

void
ActivityProfile::PrecomputeBucketTargets() {
  for (uint32_t i = 0; i < _options.buckets; ++i) {
    _ideal_bucket_count[i] = _total_bucket_count[i] * _target_train_sample_weight;
    _tolerance[i] = (_options.tolerance == 0)
                    ? 0
                    : (_ideal_bucket_count[i] * _options.tolerance) / 100;
  }
}

uint64_t
ActivityProfile::PenaltyForBucket(uint32_t bucket, uint64_t train_bucket_count) const {
  if (_total_sample_weight == 0) {
    return 0;
  }

  const uint64_t actual = train_bucket_count * _total_sample_weight;
  const uint64_t diff = AbsDiff(actual, _ideal_bucket_count[bucket]);

  if (diff <= _tolerance[bucket]) {
    return 0;
  }

  return diff - _tolerance[bucket];
}

uint64_t
ActivityProfile::Penalty(const std::vector<uint64_t>& train_bucket_count) const {
  uint64_t penalty = 0;
  for (uint32_t i = 0; i < _options.buckets; ++i) {
    penalty += PenaltyForBucket(i, train_bucket_count[i]);
  }

  return penalty;
}

int
ActivityProfile::InitialiseSplit(const std::vector<int>& in_train,
                                 const std::vector<uint32_t>& sample_weight,
                                 uint64_t total_sample_weight,
                                 uint64_t train_sample_weight) {
  if (in_train.size() != _bucket.size() || sample_weight.size() != _bucket.size()) {
    cerr << "ActivityProfile::InitialiseSplit:size mismatch\n";
    return 0;
  }

  std::fill(_train_bucket_count.begin(), _train_bucket_count.end(), 0);
  _total_sample_weight = 0;
  uint64_t observed_train_sample_weight = 0;

  for (uint32_t i = 0; i < _bucket.size(); ++i) {
    const uint32_t b = _bucket[i];
    if (b == kMissingBucket) {
      continue;
    }

    _total_sample_weight += sample_weight[i];
    if (in_train[i]) {
      _train_bucket_count[b] += sample_weight[i];
      observed_train_sample_weight += sample_weight[i];
    }
  }

  if (total_sample_weight > 0 && train_sample_weight > 0) {
    _target_train_sample_weight = (_total_sample_weight * train_sample_weight) / total_sample_weight;
  } else {
    _target_train_sample_weight = observed_train_sample_weight;
  }

  PrecomputeBucketTargets();
  _penalty = Penalty(_train_bucket_count);
  return 1;
}

uint64_t
ActivityProfile::RecomputePenalty(const std::vector<int>& in_train,
                                  const std::vector<uint32_t>& sample_weight,
                                  uint64_t total_sample_weight,
                                  uint64_t train_sample_weight) {
  InitialiseSplit(in_train, sample_weight, total_sample_weight, train_sample_weight);
  return _penalty;
}

int64_t
ActivityProfile::DeltaForSwap(uint32_t out_of_train, uint32_t into_train) const {
  const uint32_t b1 = _bucket[out_of_train];
  const uint32_t b2 = _bucket[into_train];
  if (b1 == kMissingBucket && b2 == kMissingBucket) {
    return 0;
  }

  uint64_t old_penalty = 0;
  uint64_t new_penalty = 0;

  if (b1 != kMissingBucket) {
    old_penalty += PenaltyForBucket(b1, _train_bucket_count[b1]);
    const uint64_t new_count = _train_bucket_count[b1] - (*_sample_weight)[out_of_train];
    new_penalty += PenaltyForBucket(b1, new_count);
  }

  if (b2 != kMissingBucket && b2 != b1) {
    old_penalty += PenaltyForBucket(b2, _train_bucket_count[b2]);
    const uint64_t new_count = _train_bucket_count[b2] + (*_sample_weight)[into_train];
    new_penalty += PenaltyForBucket(b2, new_count);
  } else if (b2 != kMissingBucket) {
    const uint64_t new_count = _train_bucket_count[b2] - (*_sample_weight)[out_of_train] + (*_sample_weight)[into_train];
    new_penalty = PenaltyForBucket(b2, new_count);
  }

  return static_cast<int64_t>(old_penalty) - static_cast<int64_t>(new_penalty);
}

void
ActivityProfile::ApplySwap(uint32_t out_of_train, uint32_t into_train) {
  const uint32_t b1 = _bucket[out_of_train];
  if (b1 != kMissingBucket) {
    _train_bucket_count[b1] -= (*_sample_weight)[out_of_train];
  }

  const uint32_t b2 = _bucket[into_train];
  if (b2 != kMissingBucket) {
    _train_bucket_count[b2] += (*_sample_weight)[into_train];
  }

  _penalty = Penalty(_train_bucket_count);
}

int
ActivityProfileSet::InitialiseSplit(const std::vector<int>& in_train,
                                    const std::vector<uint32_t>& sample_weight,
                                    uint64_t total_sample_weight,
                                    uint64_t train_sample_weight) {
  for (ActivityProfile& profile : _profile) {
    if (! profile.InitialiseSplit(in_train, sample_weight, total_sample_weight, train_sample_weight)) {
      return 0;
    }
  }

  return 1;
}

uint64_t
ActivityProfileSet::RecomputePenalty(const std::vector<int>& in_train,
                                     const std::vector<uint32_t>& sample_weight,
                                     uint64_t total_sample_weight,
                                     uint64_t train_sample_weight) {
  uint64_t result = 0;
  for (ActivityProfile& profile : _profile) {
    result += profile.objective_weight() * profile.RecomputePenalty(in_train, sample_weight, total_sample_weight, train_sample_weight);
  }

  return result;
}

int64_t
ActivityProfileSet::DeltaForSwap(uint32_t out_of_train, uint32_t into_train) const {
  int64_t result = 0;
  for (const ActivityProfile& profile : _profile) {
    result += static_cast<int64_t>(profile.objective_weight()) * profile.DeltaForSwap(out_of_train, into_train);
  }

  return result;
}

void
ActivityProfileSet::ApplySwap(uint32_t out_of_train, uint32_t into_train) {
  for (ActivityProfile& profile : _profile) {
    profile.ApplySwap(out_of_train, into_train);
  }
}

uint64_t
ActivityProfileSet::current_penalty() const {
  uint64_t result = 0;
  for (const ActivityProfile& profile : _profile) {
    result += profile.objective_weight() * profile.current_penalty();
  }

  return result;
}

}  // namespace train_test_split
