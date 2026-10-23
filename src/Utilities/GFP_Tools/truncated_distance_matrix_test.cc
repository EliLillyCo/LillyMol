#include "Utilities/GFP_Tools/truncated_distance_matrix.h"

#include <filesystem>
#include <string>

#include <unistd.h>

#include "gtest/gtest.h"

#include "Foundational/data_source/tfdatarecord.h"
#include "Utilities/GFP_Tools/nearneighbours.pb.h"

namespace {

using truncated_distance_matrix::ProtoType;
using truncated_distance_matrix::Storage;
using truncated_distance_matrix::TruncatedDistanceMatrix;

std::string
TmpFileName(const char* stem) {
  const auto dir = std::filesystem::temp_directory_path();
  return (dir / (std::string(stem) + std::to_string(::getpid()) + ".tfrecord")).string();
}

void
AddNbr(nnbr::NearNeighbours& proto, const std::string& id, float dist) {
  nnbr::Nbr* nbr = proto.add_nbr();
  nbr->set_id(id);
  nbr->set_dist(dist);
}

void
AddNbr(nnbr::NearNeighboursIndices& proto, uint32_t id, float dist) {
  nnbr::NbrNdx* nbr = proto.add_nbr();
  nbr->set_id(id);
  nbr->set_dist(dist);
}

std::string
WriteNameProtoFile() {
  const std::string fname = TmpFileName("tdm_name");
  iw_tf_data_record::TFDataWriter writer;
  EXPECT_TRUE(writer.Open(fname.c_str()));

  nnbr::NearNeighbours a;
  a.set_name("A");
  a.set_smiles("C");
  AddNbr(a, "B", 0.25f);
  AddNbr(a, "C", 0.50f);
  EXPECT_TRUE(writer.WriteSerializedProto(a));

  nnbr::NearNeighbours b;
  b.set_name("B");
  b.set_smiles("N");
  AddNbr(b, "A", 0.25f);
  EXPECT_TRUE(writer.WriteSerializedProto(b));

  nnbr::NearNeighbours c;
  c.set_name("C");
  c.set_smiles("O");
  AddNbr(c, "A", 0.50f);
  EXPECT_TRUE(writer.WriteSerializedProto(c));

  nnbr::NearNeighbours d;
  d.set_name("D");
  d.set_smiles("F");
  EXPECT_TRUE(writer.WriteSerializedProto(d));

  writer.Close();
  return fname;
}

std::string
WriteIndexedProtoFile() {
  const std::string fname = TmpFileName("tdm_index");
  iw_tf_data_record::TFDataWriter writer;
  EXPECT_TRUE(writer.Open(fname.c_str()));

  nnbr::NearNeighboursIndices a;
  a.set_name("A");
  AddNbr(a, 1, 0.25f);
  AddNbr(a, 2, 0.50f);
  EXPECT_TRUE(writer.WriteSerializedProto(a));

  nnbr::NearNeighboursIndices b;
  b.set_name("B");
  AddNbr(b, 0, 0.25f);
  EXPECT_TRUE(writer.WriteSerializedProto(b));

  nnbr::NearNeighboursIndices c;
  c.set_name("C");
  AddNbr(c, 0, 0.50f);
  EXPECT_TRUE(writer.WriteSerializedProto(c));

  nnbr::NearNeighboursIndices d;
  d.set_name("D");
  EXPECT_TRUE(writer.WriteSerializedProto(d));

  writer.Close();
  return fname;
}

void
TestNameProto(Storage storage) {
  const std::string fname = WriteNameProtoFile();
  TruncatedDistanceMatrix dm;
  ASSERT_TRUE(dm.Build(fname, storage, ProtoType::kNearNeighbours));

  EXPECT_EQ(dm.size(), 4u);
  EXPECT_EQ(dm.number_distances(), 2u);
  ASSERT_TRUE(dm.Index("A"));
  EXPECT_EQ(*dm.Index("A"), 0u);
  EXPECT_EQ(dm.Name(2), "C");

  ASSERT_TRUE(dm.Distance(0, 1));
  EXPECT_NEAR(*dm.Distance(0, 1), 0.25, 1.0 / 255.0);
  ASSERT_TRUE(dm.Distance(1, 0));
  EXPECT_NEAR(*dm.Distance(1, 0), 0.25, 1.0 / 255.0);
  ASSERT_TRUE(dm.Distance(0, 2));
  EXPECT_NEAR(*dm.Distance(0, 2), 0.50, 1.0 / 255.0);

  EXPECT_FALSE(dm.Distance(1, 2));
  EXPECT_NEAR(dm.DistanceOrDefault(1, 2), dm.MaxStoredDistance(), 1.0e-6);
  EXPECT_NEAR(dm.DistanceOrDefault(3, 1), dm.MaxStoredDistance(), 1.0e-6);
  EXPECT_NEAR(dm.DistanceOrDefault(3, 3), 0.0, 1.0e-6);

  EXPECT_FALSE(dm.SetDefaultDistance(0.25f));
  EXPECT_TRUE(dm.SetDefaultDistance(1.0f));
  EXPECT_NEAR(dm.DistanceOrDefault(1, 2), 1.0, 1.0e-6);

  const std::vector<float> batch = dm.DistancesOrDefault({0, 1, 3}, {1, 2, 3});
  ASSERT_EQ(batch.size(), 3u);
  EXPECT_NEAR(batch[0], 0.25, 1.0 / 255.0);
  EXPECT_NEAR(batch[1], 1.0, 1.0e-6);
  EXPECT_NEAR(batch[2], 0.0, 1.0e-6);

  std::filesystem::remove(fname);
}

}  // namespace

TEST(TruncatedDistanceMatrix, RowSparseNameProto) {
  TestNameProto(Storage::kRowSparse);
}

TEST(TruncatedDistanceMatrix, RowHashNameProto) {
  TestNameProto(Storage::kRowHash);
}

TEST(TruncatedDistanceMatrix, IndexedProto) {
  const std::string fname = WriteIndexedProtoFile();
  TruncatedDistanceMatrix dm;
  ASSERT_TRUE(dm.Build(fname, Storage::kRowSparse, ProtoType::kNearNeighboursIndices));

  EXPECT_EQ(dm.size(), 4u);
  EXPECT_EQ(dm.number_distances(), 2u);
  EXPECT_EQ(dm.Name(1), "B");
  ASSERT_TRUE(dm.Distance(2, 0));
  EXPECT_NEAR(*dm.Distance(2, 0), 0.50, 1.0 / 255.0);

  std::filesystem::remove(fname);
}


TEST(TruncatedDistanceMatrix, OneByteDuplicateDifferenceStoresSmallerValue) {
  const std::string fname = TmpFileName("tdm_one_byte");
  iw_tf_data_record::TFDataWriter writer;
  ASSERT_TRUE(writer.Open(fname.c_str()));

  nnbr::NearNeighbours a;
  a.set_name("A");
  AddNbr(a, "B", 10.0f / 255.0f);
  ASSERT_TRUE(writer.WriteSerializedProto(a));

  nnbr::NearNeighbours b;
  b.set_name("B");
  AddNbr(b, "A", 11.0f / 255.0f);
  ASSERT_TRUE(writer.WriteSerializedProto(b));

  writer.Close();

  TruncatedDistanceMatrix dm;
  ASSERT_TRUE(dm.Build(fname, Storage::kRowSparse, ProtoType::kNearNeighbours));
  EXPECT_EQ(dm.number_distances(), 1u);
  EXPECT_EQ(dm.duplicate_distances_differing_by_one(), 1u);
  ASSERT_TRUE(dm.Distance(0, 1));
  EXPECT_NEAR(*dm.Distance(0, 1), 10.0f / 255.0f, 1.0e-6);

  std::filesystem::remove(fname);
}

TEST(TruncatedDistanceMatrix, ConflictingDuplicateFails) {
  const std::string fname = TmpFileName("tdm_conflict");
  iw_tf_data_record::TFDataWriter writer;
  ASSERT_TRUE(writer.Open(fname.c_str()));

  nnbr::NearNeighbours a;
  a.set_name("A");
  AddNbr(a, "B", 0.25f);
  ASSERT_TRUE(writer.WriteSerializedProto(a));

  nnbr::NearNeighbours b;
  b.set_name("B");
  AddNbr(b, "A", 0.75f);
  ASSERT_TRUE(writer.WriteSerializedProto(b));

  writer.Close();

  TruncatedDistanceMatrix dm;
  EXPECT_FALSE(dm.Build(fname, Storage::kRowSparse, ProtoType::kNearNeighbours));

  std::filesystem::remove(fname);
}
