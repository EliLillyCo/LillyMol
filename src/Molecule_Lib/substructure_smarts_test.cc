#include <iostream>
#include <string>

#include "aromatic.h"
#include "substructure.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {

using std::cerr;
using std::endl;

using testing::UnorderedElementsAre;
using testing::Contains;
using testing::ElementsAre;
using testing::WhenSorted;
using testing::Eq;

class TestSubstructureSmarts : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    std::string _string_proto;

    IWString _smiles;

    SubstructureSearch::SubstructureQuery _proto;

    Substructure_Query _query;

    Substructure_Results _sresults;

    Molecule _m;

    static constexpr int _ntest = 10;

  protected:
    void _WriteQuery(const char * fname) {
      IWString tmp(fname);
      _query.write_msi(tmp);
    }

    const Set_of_Atoms _FirstAtomEachEmbedding()
    {
      Set_of_Atoms to_be_returned;
      for (const auto * e : _sresults.embeddings()) {
        to_be_returned.add(e->item(0));
      }

      return to_be_returned;
    }

    bool _DoPerumationsTests(const int expected);
};

void
TestSubstructureSmarts::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

bool TestSubstructureSmarts::_DoPerumationsTests(const int expected)
{
  const IWString initial_smiles = _m.smiles();

  for (int i = 0; i < _ntest; ++i) {
    const IWString smiles = _m.random_smiles();
    if (! _m.build_from_smiles(smiles)) 
    {
      cerr << "_DoPerumationsTests:invalid smiles '" << smiles << "'\n";
      return false;
    }

    const int hits = _query.substructure_search(_m, _sresults);

    if (expected == hits)
      continue;

    cerr << "_DoPerumationsTests:permutation test failure, expected " << expected << " got " << hits << " hits\n";
    cerr << "Initial smiles " << initial_smiles << endl;
    cerr << "Random  smiles " << smiles << endl;

    return false;
  }

  return true;
}

TEST_F(TestSubstructureSmarts, TestD0)
{
  _string_proto = R"(
    query {
      smarts: "[CD0]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C.CC.CC(C)C.CC(C)(C)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_matches = 1;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_matches);

  const auto matched = _FirstAtomEachEmbedding();

  EXPECT_THAT(matched, ElementsAre(0));

  EXPECT_TRUE(_DoPerumationsTests(expected_matches));
}

TEST_F(TestSubstructureSmarts, TestD1)
{
  _string_proto = R"(
    query {
      smarts: "[CD1]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C.CC.CC(C)C.CC(C)(C)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_matches = 9;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_matches);

  const Set_of_Atoms matched = _FirstAtomEachEmbedding();

  EXPECT_THAT(matched, WhenSorted(ElementsAre(1, 2, 3, 5, 6, 7, 9, 10, 11)));

  EXPECT_TRUE(_DoPerumationsTests(expected_matches));
}

TEST_F(TestSubstructureSmarts, TestD2)
{
  _string_proto = R"(
    query {
      smarts: "[CD2]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C.CCC.CC(C)CC.CC(C)(C)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_matches = 2;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_matches);

  const Set_of_Atoms matched = _FirstAtomEachEmbedding();

  EXPECT_THAT(matched, UnorderedElementsAre(2, 7));

  EXPECT_TRUE(_DoPerumationsTests(expected_matches));
}

TEST_F(TestSubstructureSmarts, TestD3)
{
  _string_proto = R"(
    query {
      smarts: "[CD3]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C.CCC.CC(C)CC.CC(C)(C)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_matches = 1;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_matches);

  const Set_of_Atoms matched = _FirstAtomEachEmbedding();

  EXPECT_THAT(matched, UnorderedElementsAre(5));

  EXPECT_TRUE(_DoPerumationsTests(expected_matches));
}

TEST_F(TestSubstructureSmarts, TestD4)
{
  _string_proto = R"(
    query {
      smarts: "[CD4]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C.CCC.CC(C)CC.CC(C)(C)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_matches = 1;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_matches);

  const Set_of_Atoms matched = _FirstAtomEachEmbedding();

  EXPECT_THAT(matched, UnorderedElementsAre(10));

  EXPECT_TRUE(_DoPerumationsTests(expected_matches));
}

TEST_F(TestSubstructureSmarts, TestNotD2)
{
  _string_proto = R"(
    query {
      smarts: "[!CD2]"   # interpreted as (Not Carbon) And (degree 2)
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC.CNC.CN.CC(C)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_matches = 1;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_matches);

  const Set_of_Atoms matched = _FirstAtomEachEmbedding();

  EXPECT_THAT(matched[0], Eq(3));

  EXPECT_TRUE(_DoPerumationsTests(expected_matches));
}

TEST_F(TestSubstructureSmarts, TestNotD2Recursive)
{
  _string_proto = R"(
    query {
      smarts: "[!$([CD2])]"   # interpreted as Not (Carbon and degree 2)
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC.CNC.CCC.CN.CC(C)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_matches = 13;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_matches);

  const Set_of_Atoms matched = _FirstAtomEachEmbedding();

  EXPECT_THAT(matched, UnorderedElementsAre(0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13));

  EXPECT_TRUE(_DoPerumationsTests(expected_matches));
}

TEST_F(TestSubstructureSmarts, TestCarbonNotDegree2)
{
  _string_proto = R"(
    query {
      smarts: "[C!D2]"   # interpreted as (Carbon) and (not degree 2)
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC.CNC.CCC.CN.CC(C)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_matches = 11;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_matches);

  const Set_of_Atoms matched = _FirstAtomEachEmbedding();

  EXPECT_THAT(matched, UnorderedElementsAre(0, 1, 2, 4, 5, 7, 8, 10, 11, 12, 13));

  EXPECT_TRUE(_DoPerumationsTests(expected_matches));
}

TEST_F(TestSubstructureSmarts, TestCarbonOrDegree2)
{
  _string_proto = R"(
    query {
      smarts: "[C,D2]"   # interpreted as (Carbon) or (degree 2)
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC.CNC.CCO.CN.SC(F)Cl";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_matches = 9;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_matches);

  const Set_of_Atoms matched = _FirstAtomEachEmbedding();

  EXPECT_THAT(matched, UnorderedElementsAre(0, 1, 2, 3, 4, 5, 6, 8, 11));

  EXPECT_TRUE(_DoPerumationsTests(expected_matches));
}

TEST_F(TestSubstructureSmarts, TestCarbonXOrDegree2)
{
  _string_proto = R"(
    query {
      smarts: "[C^D2]"   # interpreted as (Carbon) xor (degree 2)
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C.CCC.CNC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_matches = 6;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_matches);

  const Set_of_Atoms matched = _FirstAtomEachEmbedding();

  EXPECT_THAT(matched, UnorderedElementsAre(0, 1, 3, 4, 5, 6));

  EXPECT_TRUE(_DoPerumationsTests(expected_matches));
}

// It is unclear what Cr5 should be interpreted as.
// The SMARTS specification says nothing about the position in
// the smarts of the isotope. But by convention, everyone
// writes the isotope first. For now, we allow the isotope
// to be anywhere, but may revisit this.
TEST_F(TestSubstructureSmarts, TestCr5)
{
  _string_proto = R"(
    query {
      smarts: "[Cr5]"   # Isotope 5 of Chromium
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CCCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  _smiles = "C1CCCC1[Cr]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  _smiles = "C1CCCC1[5Cr]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  _smiles = "C1CCCC1[5Cr]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  _string_proto = R"(

    query {
      smarts: "[C&r5]"     # Carbon in a 5 membered ring
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  Substructure_Query query2;

  EXPECT_TRUE(query2.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  ASSERT_EQ(query2.substructure_search(_m, _sresults), 5);
}

}  // namespace
