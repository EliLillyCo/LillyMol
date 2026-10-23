#include "Utilities/GFP_Tools/train_test_split_activity.h"

#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {

using testing::ElementsAre;

std::string
TempFileName(const char* stem) {
  const char* tmpdir = std::getenv("TEST_TMPDIR");
  if (tmpdir == nullptr) {
    tmpdir = "/tmp";
  }

  std::string result(tmpdir);
  result.append("/");
  result.append(stem);
  return result;
}

void
WriteFile(const std::string& fname, const char* contents) {
  std::ofstream output(fname);
  output << contents;
}

TEST(TestActivityProfile, EqualWidthWhitespace) {
  const std::string fname = TempFileName("activity_whitespace.txt");
  WriteFile(fname, "id activity\nA 0\nB 10\nC 20\nD 30\n");

  train_test_split::ActivityProfileOptions options;
  options.fname = fname.c_str();
  options.buckets = 2;
  options.bucketisation = train_test_split::ActivityBucketisation::kEqualWidth;

  train_test_split::ActivityProfile profile;
  const std::vector<std::string> id = {"A", "B", "C", "D"};
  const std::vector<uint32_t> sample_weight = {1, 1, 1, 1};
  ASSERT_TRUE(profile.Build(options, id, sample_weight));

  EXPECT_THAT(std::vector<uint32_t>({profile.bucket(0), profile.bucket(1),
                                    profile.bucket(2), profile.bucket(3)}),
              ElementsAre(0, 0, 1, 1));
  EXPECT_THAT(std::vector<uint32_t>({profile.scaled(0), profile.scaled(1),
                                    profile.scaled(2), profile.scaled(3)}),
              ElementsAre(0, 33, 67, 100));
}

TEST(TestActivityProfile, QuantileCsvKeepsTiesTogether) {
  const std::string fname = TempFileName("activity_quantile.csv");
  WriteFile(fname, "id,activity\nA,1\nB,1\nC,3\nD,4\nE,5\n");

  train_test_split::ActivityProfileOptions options;
  options.fname = fname.c_str();
  options.buckets = 3;

  train_test_split::ActivityProfile profile;
  const std::vector<std::string> id = {"A", "B", "C", "D", "E"};
  const std::vector<uint32_t> sample_weight = {1, 1, 1, 1, 1};
  ASSERT_TRUE(profile.Build(options, id, sample_weight));

  EXPECT_EQ(profile.bucket(0), profile.bucket(1));
  EXPECT_THAT(std::vector<uint32_t>({profile.bucket(0), profile.bucket(1),
                                    profile.bucket(2), profile.bucket(3),
                                    profile.bucket(4)}),
              ElementsAre(0, 0, 1, 1, 2));
}

TEST(TestActivityProfile, WeightedSwapDeltaMatchesRecompute) {
  const std::string fname = TempFileName("activity_weighted.txt");
  WriteFile(fname, "id activity\nA 0\nB 10\nC 20\nD 30\n");

  train_test_split::ActivityProfileOptions options;
  options.fname = fname.c_str();
  options.buckets = 2;
  options.bucketisation = train_test_split::ActivityBucketisation::kEqualWidth;

  train_test_split::ActivityProfile profile;
  const std::vector<std::string> id = {"A", "B", "C", "D"};
  const std::vector<uint32_t> sample_weight = {3, 1, 2, 1};
  ASSERT_TRUE(profile.Build(options, id, sample_weight));

  std::vector<int> in_train = {1, 0, 1, 0};
  ASSERT_TRUE(profile.InitialiseSplit(in_train, sample_weight, 7, 5));
  const uint64_t initial_penalty = profile.current_penalty();
  const int64_t delta = profile.DeltaForSwap(0, 1);

  profile.ApplySwap(0, 1);
  in_train[0] = 0;
  in_train[1] = 1;
  const uint64_t recomputed = profile.RecomputePenalty(in_train, sample_weight, 7, 5);

  EXPECT_EQ(profile.current_penalty(), recomputed);
  EXPECT_EQ(static_cast<int64_t>(initial_penalty) -
            static_cast<int64_t>(profile.current_penalty()), delta);
}

TEST(TestActivityProfileSet, StickyAndPerFileDirectives) {
  const std::string fname1 = TempFileName("activity_sticky1.txt");
  const std::string fname2 = TempFileName("activity_sticky2.txt");
  WriteFile(fname1, "id activity\nA 0\nB 1\nC 2\nD 3\n");
  WriteFile(fname2, "id activity\nA 0\nB 1\nC 2\nD 3\n");

  train_test_split::ActivityProfileSet profiles;
  ASSERT_TRUE(profiles.AddDirective("buckets=4"));
  std::string directive1(fname1);
  directive1.append(":buckets=2:width");
  ASSERT_TRUE(profiles.AddDirective(directive1.c_str()));
  ASSERT_TRUE(profiles.AddDirective(fname2.c_str()));

  const std::vector<std::string> id = {"A", "B", "C", "D"};
  const std::vector<uint32_t> sample_weight = {1, 1, 1, 1};
  ASSERT_TRUE(profiles.Build(id, sample_weight));
  EXPECT_TRUE(profiles.active());
  EXPECT_EQ(profiles.number_profiles(), 2);
}


TEST(TestActivityProfile, MissingValuesAreIgnored) {
  const std::string fname = TempFileName("activity_missing.txt");
  WriteFile(fname, "id activity\nA 0\nC 20\nD 30\nX 99\n");

  train_test_split::ActivityProfileOptions options;
  options.fname = fname.c_str();
  options.buckets = 2;
  options.bucketisation = train_test_split::ActivityBucketisation::kEqualWidth;

  train_test_split::ActivityProfile profile;
  const std::vector<std::string> id = {"A", "B", "C", "D"};
  const std::vector<uint32_t> sample_weight = {1, 1, 1, 1};
  ASSERT_TRUE(profile.Build(options, id, sample_weight));

  EXPECT_FALSE(profile.missing(0));
  EXPECT_TRUE(profile.missing(1));
  EXPECT_FALSE(profile.missing(2));
  EXPECT_FALSE(profile.missing(3));
  EXPECT_EQ(profile.observed(), 3u);
  EXPECT_EQ(profile.missing_count(), 1u);
  EXPECT_EQ(profile.activity_records_not_used(), 1u);

  std::vector<int> in_train = {1, 1, 0, 0};
  ASSERT_TRUE(profile.InitialiseSplit(in_train, sample_weight, 4, 2));
  const uint64_t initial_penalty = profile.current_penalty();

  // Moving missing B out of train and observed C into train only changes C's bucket.
  const int64_t delta = profile.DeltaForSwap(1, 2);
  profile.ApplySwap(1, 2);
  in_train[1] = 0;
  in_train[2] = 1;
  EXPECT_EQ(profile.current_penalty(), profile.RecomputePenalty(in_train, sample_weight, 4, 2));
  EXPECT_EQ(static_cast<int64_t>(initial_penalty) -
            static_cast<int64_t>(profile.current_penalty()), delta);
}

TEST(TestActivityProfile, DuplicateActivityIdsAreFatal) {
  const std::string fname = TempFileName("activity_duplicate.txt");
  WriteFile(fname, "id activity\nA 0\nA 1\nB 2\nC 3\n");

  train_test_split::ActivityProfileOptions options;
  options.fname = fname.c_str();
  options.buckets = 2;

  train_test_split::ActivityProfile profile;
  const std::vector<std::string> id = {"A", "B", "C"};
  const std::vector<uint32_t> sample_weight = {1, 1, 1};
  EXPECT_FALSE(profile.Build(options, id, sample_weight));
}

TEST(TestActivityProfile, OnePopulatedBucketIsFatal) {
  const std::string fname = TempFileName("activity_one_bucket.txt");
  WriteFile(fname, "id activity\nA 1\nB 1\nC 1\nD 1\n");

  train_test_split::ActivityProfileOptions options;
  options.fname = fname.c_str();
  options.buckets = 4;

  train_test_split::ActivityProfile profile;
  const std::vector<std::string> id = {"A", "B", "C", "D"};
  const std::vector<uint32_t> sample_weight = {1, 1, 1, 1};
  EXPECT_FALSE(profile.Build(options, id, sample_weight));
}

TEST(TestActivityProfile, ToleranceCreatesDeadband) {
  const std::string fname = TempFileName("activity_tolerance.txt");
  WriteFile(fname, "id activity\nA 0\nB 10\nC 20\nD 30\n");

  train_test_split::ActivityProfileOptions options;
  options.fname = fname.c_str();
  options.buckets = 2;
  options.tolerance = 100;
  options.bucketisation = train_test_split::ActivityBucketisation::kEqualWidth;

  train_test_split::ActivityProfile profile;
  const std::vector<std::string> id = {"A", "B", "C", "D"};
  const std::vector<uint32_t> sample_weight = {1, 1, 1, 1};
  ASSERT_TRUE(profile.Build(options, id, sample_weight));

  std::vector<int> in_train = {1, 1, 0, 0};
  ASSERT_TRUE(profile.InitialiseSplit(in_train, sample_weight, 4, 2));
  EXPECT_EQ(profile.current_penalty(), 0u);
}


TEST(TestActivityProfileSet, ObjectiveWeightScalesDelta) {
  const std::string fname = TempFileName("activity_objective_weight.txt");
  WriteFile(fname, "id activity\nA 0\nB 10\nC 20\nD 30\n");

  train_test_split::ActivityProfile profile;
  train_test_split::ActivityProfileOptions options;
  options.fname = fname.c_str();
  options.buckets = 2;
  options.bucketisation = train_test_split::ActivityBucketisation::kEqualWidth;

  const std::vector<std::string> id = {"A", "B", "C", "D"};
  const std::vector<uint32_t> sample_weight = {1, 1, 1, 1};
  ASSERT_TRUE(profile.Build(options, id, sample_weight));

  std::vector<int> in_train = {1, 1, 0, 0};
  ASSERT_TRUE(profile.InitialiseSplit(in_train, sample_weight, 4, 2));
  const int64_t unweighted_delta = profile.DeltaForSwap(0, 2);

  train_test_split::ActivityProfileSet profiles;
  std::string directive(fname);
  directive.append(":buckets=2:width:objective_weight=3");
  ASSERT_TRUE(profiles.AddDirective(directive.c_str()));
  ASSERT_TRUE(profiles.Build(id, sample_weight));
  ASSERT_TRUE(profiles.InitialiseSplit(in_train, sample_weight, 4, 2));

  EXPECT_EQ(profiles.DeltaForSwap(0, 2), 3 * unweighted_delta);
}


TEST(TestCombinedActivityDistanceDelta, CombinesSignedDeltas) {
  EXPECT_EQ(train_test_split::CombinedActivityDistanceDelta(10, 0, 100), 10);
  EXPECT_EQ(train_test_split::CombinedActivityDistanceDelta(10, 2, 3), 16);
  EXPECT_EQ(train_test_split::CombinedActivityDistanceDelta(10, -2, 3), 4);
  EXPECT_EQ(train_test_split::CombinedActivityDistanceDelta(-10, 4, 3), 2);
  EXPECT_EQ(train_test_split::CombinedActivityDistanceDelta(10, 4, 0), 10);
}

}  // namespace
