#ifndef UTILITIES_GFP_TOOLS_TRAIN_TEST_SPLIT_ACTIVITY_H_
#define UTILITIES_GFP_TOOLS_TRAIN_TEST_SPLIT_ACTIVITY_H_

#include <cstdint>
#include <iosfwd>
#include <limits>
#include <string>
#include <vector>

#include "Foundational/iwstring/iwstring.h"

namespace train_test_split {

enum class ActivityBucketisation {
  kQuantile,
  kEqualWidth,
};

struct ActivityProfileOptions {
  IWString fname;
  uint32_t buckets = 10;
  uint32_t objective_weight = 1;
  // Percent tolerance around each bucket's ideal train count. Zero means exact.
  uint32_t tolerance = 0;
  ActivityBucketisation bucketisation = ActivityBucketisation::kQuantile;
};

int64_t CombinedActivityDistanceDelta(int64_t distance_delta, int64_t activity_delta,
                                      uint32_t activity_objective_weight);

class ActivityProfile {
  private:
    static constexpr uint32_t kMissingBucket = std::numeric_limits<uint32_t>::max();
    static constexpr uint8_t kMissingScaled = std::numeric_limits<uint8_t>::max();

    ActivityProfileOptions _options;
    std::vector<uint8_t> _scaled;
    std::vector<uint32_t> _bucket;
    std::vector<uint64_t> _total_bucket_count;
    std::vector<uint64_t> _train_bucket_count;
    std::vector<uint64_t> _ideal_bucket_count;
    std::vector<uint64_t> _tolerance;
    const std::vector<uint32_t>* _sample_weight = nullptr;
    uint64_t _total_sample_weight = 0;
    uint64_t _target_train_sample_weight = 0;
    uint64_t _observed = 0;
    uint64_t _missing = 0;
    uint64_t _activity_records_not_used = 0;
    uint64_t _penalty = 0;

    int AssignBuckets(const std::vector<double>& activity,
                      const std::vector<uint8_t>& seen, int verbose);
    int AssignQuantileBuckets(const std::vector<double>& activity,
                              const std::vector<uint8_t>& seen);
    int AssignEqualWidthBuckets(const std::vector<uint8_t>& seen);
    int CheckBucketPopulation(int verbose) const;

    void PrecomputeBucketTargets();
    uint64_t PenaltyForBucket(uint32_t bucket, uint64_t train_bucket_count) const;
    uint64_t Penalty(const std::vector<uint64_t>& train_bucket_count) const;

  public:
    int Build(const ActivityProfileOptions& options,
              const std::vector<std::string>& id,
              const std::vector<uint32_t>& sample_weight,
              int verbose = 0);

    int InitialiseSplit(const std::vector<int>& in_train,
                        const std::vector<uint32_t>& sample_weight,
                        uint64_t total_sample_weight = 0,
                        uint64_t train_sample_weight = 0);
    uint64_t RecomputePenalty(const std::vector<int>& in_train,
                              const std::vector<uint32_t>& sample_weight,
                              uint64_t total_sample_weight = 0,
                              uint64_t train_sample_weight = 0);
    int64_t DeltaForSwap(uint32_t out_of_train, uint32_t into_train) const;
    void ApplySwap(uint32_t out_of_train, uint32_t into_train);

    bool missing(uint32_t i) const {
      return _bucket[i] == kMissingBucket;
    }
    uint32_t bucket(uint32_t i) const {
      return _bucket[i];
    }
    uint8_t scaled(uint32_t i) const {
      return _scaled[i];
    }
    uint64_t current_penalty() const {
      return _penalty;
    }
    uint32_t nbuckets() const {
      return _options.buckets;
    }
    uint32_t objective_weight() const {
      return _options.objective_weight;
    }
    uint64_t observed() const {
      return _observed;
    }
    uint64_t missing_count() const {
      return _missing;
    }
    uint64_t activity_records_not_used() const {
      return _activity_records_not_used;
    }
    const IWString& fname() const {
      return _options.fname;
    }
};

class ActivityProfileSet {
  private:
    std::vector<ActivityProfileOptions> _pending;
    std::vector<ActivityProfile> _profile;
    ActivityProfileOptions _defaults;

  public:
    int AddDirective(const const_IWSubstring& directive);
    int Build(const std::vector<std::string>& id,
              const std::vector<uint32_t>& sample_weight,
              int verbose = 0);

    int InitialiseSplit(const std::vector<int>& in_train,
                        const std::vector<uint32_t>& sample_weight,
                        uint64_t total_sample_weight = 0,
                        uint64_t train_sample_weight = 0);
    uint64_t RecomputePenalty(const std::vector<int>& in_train,
                              const std::vector<uint32_t>& sample_weight,
                              uint64_t total_sample_weight = 0,
                              uint64_t train_sample_weight = 0);
    int64_t DeltaForSwap(uint32_t out_of_train, uint32_t into_train) const;
    void ApplySwap(uint32_t out_of_train, uint32_t into_train);

    bool active() const {
      return !_profile.empty();
    }
    bool has_pending() const {
      return !_pending.empty();
    }
    int number_profiles() const {
      return _profile.size();
    }
    uint64_t current_penalty() const;
};

}  // namespace train_test_split

#endif  // UTILITIES_GFP_TOOLS_TRAIN_TEST_SPLIT_ACTIVITY_H_
