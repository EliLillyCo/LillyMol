// Test suite for the feture scaler class

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"
#include "google/protobuf/text_format.h"

#define FEATURE_SCALER_IMPLEMENTATION
#include "scaler.h"

namespace {
TEST(TestFeatureScaler, TestInRange1) {
  const std::string proto_text = R"pb(
    min: -1
    max: 1
    )pb";
  FeatureScaling::FeatureScaling proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(proto_text, &proto));

  feature_scaler::FeatureScaler<float> scaler = feature_scaler::FeatureScaler<float>::Build(proto);

  EXPECT_FLOAT_EQ(scaler.Scale(-1.0), 0.0);
  EXPECT_FLOAT_EQ(scaler.Scale(-0.8), 0.1);
  EXPECT_FLOAT_EQ(scaler.Scale(-0.5), 0.25);
  EXPECT_FLOAT_EQ(scaler.Scale(0.0), 0.5);
  EXPECT_FLOAT_EQ(scaler.Scale(1.0), 1.0);
}

TEST(TestFeatureScaler, TestExtrapolate) {
  const std::string proto_text = R"pb(
    min: -3
    max: -1
    )pb";
  FeatureScaling::FeatureScaling proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(proto_text, &proto));

  feature_scaler::FeatureScaler<float> scaler = feature_scaler::FeatureScaler<float>::Build(proto);

  EXPECT_FLOAT_EQ(scaler.Scale(-3.0), 0.0);
  EXPECT_FLOAT_EQ(scaler.Scale(-2.0), 0.5);
  EXPECT_FLOAT_EQ(scaler.Scale(-1.0), 1.0);

  EXPECT_FLOAT_EQ(scaler.Scale(-4.0), -0.5);
  EXPECT_FLOAT_EQ(scaler.Scale(-5.0), -1.0);

  EXPECT_FLOAT_EQ(scaler.Scale(0.0), 1.5);
  EXPECT_FLOAT_EQ(scaler.Scale(1.0), 2.0);
}

TEST(TestFeatureScaler, TestTruncate) {
  const std::string proto_text = R"pb(
    min: -3
    max: 4
    )pb";
  FeatureScaling::FeatureScaling proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(proto_text, &proto));

  feature_scaler::FeatureScaler<float> scaler = feature_scaler::FeatureScaler<float>::Build(proto);
  scaler.set_truncate_out_of_range(true);
  EXPECT_FLOAT_EQ(scaler.Scale(-4.0), 0.0);
  EXPECT_FLOAT_EQ(scaler.Scale(-3.0), 0.0);
  EXPECT_FLOAT_EQ(scaler.Scale(4.0), 1.0);
  EXPECT_FLOAT_EQ(scaler.Scale(5.0), 1.0);
}

TEST(TestFeatureScaler, TestScaleBack1) {
  const std::string proto_text = R"pb(
    min: 3
    max: 8
    )pb";

  FeatureScaling::FeatureScaling proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(proto_text, &proto));

  feature_scaler::FeatureScaler<float> scaler = feature_scaler::FeatureScaler<float>::Build(proto);

  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(0.0), 3.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(0.5), 5.5);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(1.0), 8.0);
}

TEST(TestFeatureScaler, TestScaleBack2) {
  const std::string proto_text = R"pb(
    min: -5
    max: 5
    )pb";

  FeatureScaling::FeatureScaling proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(proto_text, &proto));

  feature_scaler::FeatureScaler<float> scaler = feature_scaler::FeatureScaler<float>::Build(proto);

  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(0.0), -5.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(0.5), 0.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(1.0), 5.0);

  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(-0.1), -6.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(1.1), 6.0);
}

TEST(TestFeatureScaler, TestScaleBackTruncate) {
  const std::string proto_text = R"pb(
    min: 5
    max: 9
    )pb";

  FeatureScaling::FeatureScaling proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(proto_text, &proto));

  feature_scaler::FeatureScaler<float> scaler = feature_scaler::FeatureScaler<float>::Build(proto);
  scaler.set_truncate_out_of_range(true);

  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(0.0), 5.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(0.5), 7.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(1.0), 9.0);

  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(-0.1), 5.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(1.1), 9.0);
}

TEST(FeatureScaler, Test11ScalingTruncated) {
  const std::string proto_text = R"pb(
    min: 5
    max: 9
    range_type: R11
  )pb";

  FeatureScaling::FeatureScaling proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(proto_text, &proto));

  feature_scaler::FeatureScaler<float> scaler = feature_scaler::FeatureScaler<float>::Build(proto);
  scaler.set_truncate_out_of_range(true);

  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(5.0)), 5.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(7.0)), 7.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(9.0)), 9.0);

  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(4.0)), 5.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(10.0)), 9.0);
}

TEST(FeatureScaler, Test11ScalingNotTruncated) {
  const std::string proto_text = R"pb(
    min: 5
    max: 9
    range_type: R11
  )pb";

  FeatureScaling::FeatureScaling proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(proto_text, &proto));

  feature_scaler::FeatureScaler<float> scaler = feature_scaler::FeatureScaler<float>::Build(proto);
  scaler.set_truncate_out_of_range(false);

  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(5.0)), 5.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(7.0)), 7.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(9.0)), 9.0);

  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(4.0)), 4.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(10.0)), 10.0);
}

TEST(FeatureScaler, Test11ScalingNotTruncated11SetViaAPI) {
  const std::string proto_text = R"pb(
    min: 5
    max: 9
  )pb";

  FeatureScaling::FeatureScaling proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(proto_text, &proto));

  feature_scaler::FeatureScaler<float> scaler = feature_scaler::FeatureScaler<float>::Build(proto);
  scaler.set_truncate_out_of_range(false);
  scaler.Initialise11Scaling();

  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(5.0)), 5.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(7.0)), 7.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(9.0)), 9.0);

  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(4.0)), 4.0);
  EXPECT_FLOAT_EQ(scaler.ScaleBackToOrignalRange(scaler.Scale(10.0)), 10.0);
}

}  // namespace

