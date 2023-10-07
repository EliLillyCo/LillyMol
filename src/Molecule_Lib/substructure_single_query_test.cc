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

class TestSubstructure : public testing::Test
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
TestSubstructure::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

bool TestSubstructure::_DoPerumationsTests(const int expected)
{
  const IWString initial_smiles = _m.smiles();

  for (int i = 0; i < _ntest; ++i) {
    const IWString smiles = _m.random_smiles();
    if (! _m.build_from_smiles(smiles)) {
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

TEST_F(TestSubstructure, TestOneEmbeddingPerStartAtom)
{
  _string_proto = R"(
    query {
      one_embedding_per_start_atom: true
      smarts: "C(F)(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F.FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestUniqueEmbeddingsOnly)
{
  _string_proto = R"(
    query {
      unique_embeddings_only: true
      smarts: "FC(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestMaxMatchesToFine)
{
  _string_proto = R"(
    query {
      max_matches_to_find: 3
      smarts: "FC(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPerumationsTests(3));
}

TEST_F(TestSubstructure, TestSubtractFromRc)
{
  _string_proto = R"(
    query {
      subtract_from_rc: 2
      smarts: "FC(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 22);

  EXPECT_TRUE(_DoPerumationsTests(22));
}

TEST_F(TestSubstructure, TestHitsNeeded)
{
  _string_proto = R"(
    query {
      min_hits_needed: 15
      smarts: "FC(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 24);

  EXPECT_TRUE(_DoPerumationsTests(24));
}

TEST_F(TestSubstructure, TestNoSymmetricHits)
{
  _string_proto = R"(
    query {
      perceive_symmetric_equivalents: false
      smarts: "FC(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestEmbeddingsDoNotOverlap1)
{
  _string_proto = R"(
    query {
      embeddings_do_not_overlap: true
      smarts: "FC(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();


  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

// Hmmm, was going to use benzene as a target molecule, with the query cc, but then
// I realised this is undefined. There might be either 2 or 3 possible embeddings that
// do not violate the overlapping embeddings constraint. Doing it optimally would be
// hard, so I just do not run that test.
TEST_F(TestSubstructure, TestEmbeddingsDoNotOverlapAmbiguous)
{
  _string_proto = R"(
    query {
      embeddings_do_not_overlap: true
      smarts: "cc"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  // currently the software returns 3, but there is no guarantee of that.
//EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);
}

TEST_F(TestSubstructure, TestEmbeddingsDoNotOverlapEasy)
{
  _string_proto = R"(
    query {
      embeddings_do_not_overlap: true
      smarts: "O=C-N"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();


  _smiles = "CC(=O)NC(=O)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestAllHitsInSameFragment)
{
  _string_proto = R"(
    query {
      all_hits_in_same_fragment: true
      smarts: "O=C-N"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();


  _smiles = "CC(=O)NC(=O)C.CC(=O)NC(=O)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestOnlyMatchLargestFragment)
{
  _string_proto = R"(
    query {
      only_match_largest_fragment: true
      smarts: "O=C-N"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC(=O)NC(=O)C.CC(=O)NC(=O)CC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  Set_of_Atoms matched;
  for (const auto* e : _sresults.embeddings()) {
    matched += *e;
  }

  EXPECT_THAT(matched, UnorderedElementsAre(8, 9, 10, 10, 11, 12));

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestOnlyMatchLargestFragmentNoMatchesInFirst)
{
  _string_proto = R"(
    query {
      only_match_largest_fragment: true
      smarts: "O=C-N"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC(=O)NC(=O)CC.CCCCCCCCCCCCCCCCCCCCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  _DoPerumationsTests(0);
}

TEST_F(TestSubstructure, TestImplicitRingCondition1)
{
  _string_proto = R"(
    query {
      implicit_ring_condition: 1
      smarts: "C-C-C-C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 8);

  EXPECT_TRUE(_DoPerumationsTests(8));
}

TEST_F(TestSubstructure, TestImplicitRingCondition0)
{
  _string_proto = R"(
    query {
      implicit_ring_condition: 0
      smarts: "C-C-C-C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestDistanceBetweenHitsMax)
{
  _string_proto = R"(
    query {
      max_distance_between_hits: 3
      smarts: "F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Fc1c(F)cc(F)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms matches = _FirstAtomEachEmbedding();

  EXPECT_THAT(matches, UnorderedElementsAre(0, 3));

//EXPECT_TRUE(_DoPerumationsTests(2)); // permutation tests not done because of ordering effects
}

// Here again results are ambiguous and not tested.
// The algorithm works by removing hits that are too close
// to a previous match. but since the order in which matches
// are discovered is undefined, there is no real way of 
// telling how the winnowing will work.
TEST_F(TestSubstructure, TestDistanceBetweenHitsMin)
{
  _string_proto = R"(
    query {
      min_distance_between_hits: 4
      smarts: "F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Fc1c(F)cc(F)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms matches = _FirstAtomEachEmbedding();

  EXPECT_THAT(matches, Contains(6));

  EXPECT_TRUE(_DoPerumationsTests(2));
}

// This might be flaky. I think this should work with distance_between_hits:5
// but it did not. Should be investigated.
TEST_F(TestSubstructure, TestDistanceBetweenHitsMatchedAtom)
{
  _string_proto = R"(
    query {
      distance_between_hits_ncheck: 2
      max_distance_between_hits: 5
      smarts: "OC(=O)c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms matches = _FirstAtomEachEmbedding();

  EXPECT_THAT(matches, UnorderedElementsAre(0, 7));

//EXPECT_TRUE(_DoPerumationsTests(2));  permutation tests not done
}

TEST_F(TestSubstructure, TestAttachedHeteroatomCount)
{
  _string_proto = R"(
    query {
      attached_heteroatom_count: 0
      smarts: "OC(=O)c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms matches = _FirstAtomEachEmbedding();

  EXPECT_THAT(matches, UnorderedElementsAre(0, 7));

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestRingAtomsMatched)
{
  _string_proto = R"(
    query {
      ring_atoms_matched: 1
      smarts: "OC(=O)c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPerumationsTests(3));
}

TEST_F(TestSubstructure, TestHeteroatomsMatched)
{
  _string_proto = R"(
    query {
      heteroatoms_matched: 2
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);

  EXPECT_TRUE(_DoPerumationsTests(4));
}

TEST_F(TestSubstructure, TestHeteroatomsInMoleculeNoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        heteroatoms_in_molecule: 3
      }
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestHeteroatomsInMoleculeMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        heteroatoms_in_molecule: 9
      }
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);

  EXPECT_TRUE(_DoPerumationsTests(4));
}

TEST_F(TestSubstructure, TestNatomsNo)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        natoms: 9
      }
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNatomsMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        natoms: 19
      }
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);

  EXPECT_TRUE(_DoPerumationsTests(4));
}

TEST_F(TestSubstructure, TestNringsNoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        nrings: 2
      }
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNringsMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        nrings: 1
      }
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);

  EXPECT_TRUE(_DoPerumationsTests(4));
}

TEST_F(TestSubstructure, TestNcon)
{
  _string_proto = R"(
    query {
      ncon: 2
      smarts: "FC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FCC.FC(C)C.FC(C)(C)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms matched = _FirstAtomEachEmbedding();

  EXPECT_EQ(matched[0], 3);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNcon2)
{
  _string_proto = R"(
    query {
      ncon: 2
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NC.CNC.CNCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(6, 7));

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestHeteroatomsInMoleculeNo)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        heteroatoms_in_molecule: 2
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestHeteroatomsInMoleculeMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        heteroatoms_in_molecule: 2
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCN";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestNringsNo)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        nrings: 1
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCN";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNrings0OK)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        nrings: 0
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCN";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestNrings1)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        nrings: 1
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestFusedRingsNone)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        fused_rings: 1
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestFusedRingsZeroOk)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        fused_rings: 0
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestFusedRings2)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        fused_rings: 2
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CNC1CC2C1CC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestStronglyFusedRingsNo)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        min_strongly_fused_rings: 2
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1CCC12C3C4C1C5C2C3C45";  // Cubane has no rings that share more than 1 bond.

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestStronglyFusedRingsOK)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        min_strongly_fused_rings: 2
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1C1C3CC2CC(CC1C2)C3";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}


TEST_F(TestSubstructure, TestIsolatedRingsNotIsolated)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        isolated_rings: 1
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c12cccc1c(CN)ccn2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestIsolatedRingsNoRing)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        isolated_rings: 1
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestIsolatedRingsOK)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        isolated_rings: 1
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCc1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestIsolatedRingObjectsNoRings)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        ring_systems: 1
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestIsolatedRingObjectsNotMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        ring_systems: 2
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCc1cccc2ccncc12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestIsolatedRingObjectsMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        ring_systems: 1
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCc1cccc2ccncc12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestAromaticAtomsNone)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        aromatic_atoms: 0
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestAromaticAtoms0NoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        aromatic_atoms: 0
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCCc1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestAromaticAtomsOK)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        aromatic_atoms: 6
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCCc1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestHeteroatomsMatches)
{
  _string_proto = R"(
    query {
      heteroatoms: [9, 17]
      attached_heteroatom_count: 1
      heteroatoms_matched: 1
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCF";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestHeteroatomsNoMatches)
{
  _string_proto = R"(
    query {
      heteroatoms: [9, 17]
      attached_heteroatom_count: 1
      heteroatoms_matched: 1
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestElementHitsNeededMatches)
{
  _string_proto = R"(
    query {
      element_hits_needed {
        atomic_number: 7
        hits_needed: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestElementHitsNeededNoMatches)
{
  _string_proto = R"(
    query {
      element_hits_needed {
        atomic_number: 7
        hits_needed: 2
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestElementsNeededNone)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        elements_needed {
          atomic_number: 7
          hits_needed: 1
        }
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestElementsNeededNumberNoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        elements_needed {
          atomic_number: 7
          hits_needed: 2
        }
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestElementsNeededNumberMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        elements_needed {
          atomic_number: 7
          hits_needed: 1
        }
        elements_needed {
          atomic_number: 6
          hits_needed: 3
        }
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNonAromaticRingsNone)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        non_aromatic_rings: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNonAromaticRingsNoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        non_aromatic_rings: 2
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNonAromaticRingsMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        non_aromatic_rings: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNumberIsotopicAtomsNone)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        number_isotopic_atoms: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNumberIsotopicAtomsNoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        number_isotopic_atoms: 2
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1[3CH2]C1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNumberIsotopicAtomsMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        number_isotopic_atoms: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1[3CH2]C1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNumberFragmentsNoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        number_fragments: 2
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNumberFragmentsMatch1)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        number_fragments: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNumberFragmentsMatch2)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        number_fragments: 2
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1CC1.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestAtomsInSpinachAllSpinach)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        atoms_in_spinach: 7
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCCC.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestAtomsInSpinachMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        atoms_in_spinach: 3
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestInterRingAtomsNoneMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        inter_ring_atoms: 0
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestInterRingAtomsZero)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        inter_ring_atoms: 0
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1C1CC1.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestInterRingAtomsOne)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        inter_ring_atoms: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1CC1CC1.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestUnmatchedAtomsNoMatch)
{
  _string_proto = R"(
    query {
      unmatched_atoms: 1
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1C.[Br-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestUnmatchedAtomsMatch)
{
  _string_proto = R"(
    query {
      unmatched_atoms: 5
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1C.[Br-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNetFormalChargeZeroMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        net_formal_charge: 0
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNetFormalChargeZeroNoMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        net_formal_charge: 0
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNetFormalChargePositive)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        net_formal_charge: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "[NH3+]CC1CCC1C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNetFormalChargeNegative)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        net_formal_charge: -1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1C.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestMinFractionAtomsMatchedFailOne)
{
  _string_proto = R"(
    query {
      min_fraction_atoms_matched: 0.34
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1CO.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestMinFractionAtomsMatchedMatchOne)
{
  _string_proto = R"(
    query {
      min_fraction_atoms_matched: 0.330
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1CO.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestMinFractionAtomsMatchedMatchTwo)
{
  _string_proto = R"(
    query {
      min_fraction_atoms_matched: 0.660
      smarts: "N-CC.CCO";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1CO.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestMinFractionAtomsMatchedMaxViolated)
{
  _string_proto = R"(
    query {
      max_fraction_atoms_matched: 0.330
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1CO.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestRespectInitialAlllNumbering)
{
  _string_proto = R"(
    query {
      respect_initial_atom_numbering: true
      query_atom {
        id: 5
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 6
        }
        single_bond: 5
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1CO.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms * e = _sresults.embedding(0);

  EXPECT_THAT(*e, testing::ElementsAre(-1, -1, -1, 1, -1, 0));

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNoMatchedAtomsBetweenMatches)
{
  _string_proto = R"(
    query {
      unique_embeddings_only: true
      no_matched_atoms_between {
        a1: 1
        a2: 4
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 8
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 8
        }
        double_bond: 1
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 8
        }
      }
      query_atom {
        id: 4
        atom_properties {
          atomic_number: 6
        }
        single_bond: 3
      }
      query_atom {
        id: 5
        atom_properties {
          atomic_number: 8
        }
        double_bond: 4
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "COC(=O)CCCCC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNoMatchedAtomsBetweenNoMatch)
{
  _string_proto = R"(
    query {
      unique_embeddings_only: true
      no_matched_atoms_between {
        a1: 1
        a2: 4
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 8
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 8
        }
        double_bond: 1
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 8
        }
      }
      query_atom {
        id: 4
        atom_properties {
          atomic_number: 6
        }
        single_bond: 3
      }
      query_atom {
        id: 5
        atom_properties {
          atomic_number: 8
        }
        double_bond: 4
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC(=O)OCCCCC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestLinkAtomsMatches)
{
  _string_proto = R"(
    query {
      link_atoms {
        a1: 1
        a2: 2
        distance: 3
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 9
        }
        single_bond: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCCCF";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestLinkAtomsNoMatches)
{
  _string_proto = R"(
    query {
      link_atoms {
        a1: 1
        a2: 2
        min_distance: 4
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 9
        }
        single_bond: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCCCF";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

}  // namespace
