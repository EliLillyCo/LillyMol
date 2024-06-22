// Tests for substructure_atom

#include <string>

#include "aromatic.h"
#include "substructure.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {

using testing::UnorderedElementsAre;

TEST(TestSubstructureAtom, TestOr) {
  const std::string string_proto = R"pb(
query {
  query_atom {
    id: 0
    atom_properties {
      atomic_number: 6
    }
    atom_properties {
      atomic_number: 7
      logical_operator: SS_OR
    }
  }
}
  )pb";

  SubstructureSearch::SubstructureQuery proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  Substructure_Query query;
  ASSERT_TRUE(query.ConstructFromProto(proto));

  Molecule m1;
  ASSERT_TRUE(m1.build_from_smiles("C"));
  EXPECT_EQ(query.substructure_search(&m1), 1);

  Molecule m2;
  ASSERT_TRUE(m2.build_from_smiles("N"));
  EXPECT_EQ(query.substructure_search(&m2), 1);
}

TEST(TestSubstructureAtom, TestAnd) {
  const std::string string_proto = R"pb(
query {
  query_atom {
    id: 0
    atom_properties {
      atomic_number: 6
    }
    atom_properties {
      aromatic: true
      logical_operator: SS_AND
    }
  }
}
  )pb";

  SubstructureSearch::SubstructureQuery proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  Substructure_Query query;
  ASSERT_TRUE(query.ConstructFromProto(proto));

  Molecule m1;
  ASSERT_TRUE(m1.build_from_smiles("C"));
  EXPECT_EQ(query.substructure_search(&m1), 0);

  Molecule m2;
  ASSERT_TRUE(m2.build_from_smiles("c1ccccc1"));
  EXPECT_EQ(query.substructure_search(&m2), 6);
}

TEST(TestSubstructureAtom, TestImpossibleAnd) {
  const std::string string_proto = R"pb(
query {
  query_atom {
    id: 0
    atom_properties {
      atomic_number: 6
    }
    atom_properties {
      atomic_number: 7
      logical_operator: SS_AND
    }
  }
}
  )pb";

  SubstructureSearch::SubstructureQuery proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  Substructure_Query query;
  ASSERT_TRUE(query.ConstructFromProto(proto));

  Molecule m1;
  ASSERT_TRUE(m1.build_from_smiles("C"));
  EXPECT_EQ(query.substructure_search(&m1), 0);

  Molecule m2;
  ASSERT_TRUE(m2.build_from_smiles("CN"));
  EXPECT_EQ(query.substructure_search(&m2), 0);

  Molecule m3;
  ASSERT_TRUE(m3.build_from_smiles("c1ccccc1"));
  EXPECT_EQ(query.substructure_search(&m3), 0);

  Molecule m4;
  ASSERT_TRUE(m4.build_from_smiles("c1ncccc1"));
  EXPECT_EQ(query.substructure_search(&m4), 0);
}

TEST(TestSubstructureAtom, TestXor) {
  const std::string string_proto = R"pb(
query {
  query_atom {
    id: 0
    atom_properties {
      atomic_number: 6
    }
    atom_properties {
      aromatic: true
      logical_operator: SS_XOR
    }
  }
}
  )pb";

  SubstructureSearch::SubstructureQuery proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  Substructure_Query query;
  ASSERT_TRUE(query.ConstructFromProto(proto));

  Molecule m1;
  ASSERT_TRUE(m1.build_from_smiles("C"));
  EXPECT_EQ(query.substructure_search(&m1), 1);

  Molecule m2;
  ASSERT_TRUE(m2.build_from_smiles("c1ccccc1"));
  EXPECT_EQ(query.substructure_search(&m2), 0);
}

}  // namespace
