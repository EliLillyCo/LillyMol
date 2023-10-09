// Tests for ring_extraction

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "ring_ext_rep.h"

namespace {

using ring_replacement::RingHash;

TEST(TestRingSizeHash, Test4A) {
  extending_resizable_array<int> aliph, arom;
  aliph[4] = 1;
  EXPECT_EQ(RingHash(aliph, arom), 40);
}

TEST(TestRingSizeHash, Test5A) {
  extending_resizable_array<int> aliph, arom;
  aliph[5] = 1;
  EXPECT_EQ(RingHash(aliph, arom), 50);
}

TEST(TestRingSizeHash, Test6A) {
  extending_resizable_array<int> aliph, arom;
  aliph[6] = 1;
  EXPECT_EQ(RingHash(aliph, arom), 60);
}

TEST(TestRingSizeHash, Test7A) {
  extending_resizable_array<int> aliph, arom;
  aliph[7] = 1;
  EXPECT_EQ(RingHash(aliph, arom), 70);
}

TEST(TestRingSizeHash, Test4a) {
  extending_resizable_array<int> aliph, arom;
  arom[4] = 1;
  EXPECT_EQ(RingHash(aliph, arom), 41);
}

TEST(TestRingSizeHash, Test5a) {
  extending_resizable_array<int> aliph, arom;
  arom[5] = 1;
  EXPECT_EQ(RingHash(aliph, arom), 51);
}

TEST(TestRingSizeHash, Test4A4A) {
  extending_resizable_array<int> aliph, arom;
  aliph[4] = 2;
  EXPECT_EQ(RingHash(aliph, arom), 4040);
}

TEST(TestRingSizeHash, Test4A4a) {
  extending_resizable_array<int> aliph, arom;
  aliph[4] = 1;
  arom[4] = 1;
  EXPECT_EQ(RingHash(aliph, arom), 4041);
}

TEST(TestRingSizeHash, Test4A5a) {
  extending_resizable_array<int> aliph, arom;
  aliph[4] = 1;
  arom[5] = 1;
  EXPECT_EQ(RingHash(aliph, arom), 4051);
}

TEST(TestRingSizeHash, Test4A5a6A) {
  extending_resizable_array<int> aliph, arom;
  aliph[4] = 1;
  aliph[6] = 1;
  arom[5] = 1;
  EXPECT_EQ(RingHash(aliph, arom), 405160);
}

}  // namespace
