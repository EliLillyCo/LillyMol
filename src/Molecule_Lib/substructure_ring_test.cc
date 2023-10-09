#include <iostream>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "aromatic.h"
#include "substructure.h"

namespace {

using std::cerr;
using std::endl;

using testing::UnorderedElementsAre;
using testing::Contains;

class TestSubstructureRing : public testing::Test
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

    bool _DoPermutationTests(const int expected);
};

void
TestSubstructureRing::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

bool TestSubstructureRing::_DoPermutationTests(const int expected)
{
  const IWString initial_smiles = _m.smiles();

  for (int i = 0; i < _ntest; ++i) {
    const IWString smiles = _m.random_smiles();
    if (! _m.build_from_smiles(smiles)) 
    {
      cerr << "_DoPermutationTests:invalid smiles '" << smiles << "'\n";
      return false;
    }

    const int hits = _query.substructure_search(_m, _sresults);

    if (expected == hits)
      continue;

    cerr << "_DoPermutationTests:permutation test failure, expected " << expected << " got " << hits << " hits\n";
    cerr << "Initial smiles " << initial_smiles << endl;
    cerr << "Random  smiles " << smiles << endl;

    return false;
  }

  return true;
}

TEST_F(TestSubstructureRing, TestRingHitsNeededMatches)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          hits_needed: 2
        }
        ring_size: 3
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1CC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 7);

  EXPECT_TRUE(_DoPermutationTests(7));
}

TEST_F(TestSubstructureRing, TestRingHitsNeededNoMatches)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          min_hits_needed: 3
        }
        ring_size: 3
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1CC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRing, TestAttachedHeteroatomCountMatches)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          attached_heteroatom_count: 1
        }
        ring_size: 5
        aromatic: true
      }
      smarts: "c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ncnc1N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPermutationTests(3));
}

TEST_F(TestSubstructureRing, TestAttachedHeteroatomCountNoMatches)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          attached_heteroatom_count: 0
        }
        ring_size: 5
        aromatic: true
      }
      smarts: "c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ncnc1N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRing, TestHeteroatomCountMatches)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          heteroatom_count: 2
        }
        ring_size: 5
        aromatic: true
      }
      smarts: "c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ncnc1N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPermutationTests(3));
}

TEST_F(TestSubstructureRing, TestHeteroatomCountNoMatches)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          heteroatom_count: 1
        }
        ring_size: 5
        aromatic: true
      }
      smarts: "c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ncnc1N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRing, TestNconMatchesIsolated)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          ncon: 1
        }
        ring_size: 5
        aromatic: true
      }
      smarts: "c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ncnc1N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPermutationTests(3));
}

TEST_F(TestSubstructureRing, TestNconMatchesFused)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          ncon: 2
        }
        ring_size: 5
        aromatic: true
      }
      smarts: "c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  // Ring fusion bonds count as exocyclic.
  _smiles = "c12cnnc2ccc(N)c1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 7);

  EXPECT_TRUE(_DoPermutationTests(7));
}

TEST_F(TestSubstructureRing, TestallHitsInSameFragmentMatches)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          all_hits_in_same_fragment: true
          hits_needed: 2
        }
        ring_size: 5
        aromatic: true
      }
      smarts: "c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  // Ring fusion bonds count as exocyclic.
  _smiles = "n1occc1Cc1ncoc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 6);

  EXPECT_TRUE(_DoPermutationTests(6));
}

TEST_F(TestSubstructureRing, TestallHitsInSameFragmentNoMatches)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          all_hits_in_same_fragment: true
          hits_needed: 2
        }
        ring_size: 5
        aromatic: true
      }
      smarts: "c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  // Ring fusion bonds count as exocyclic.
  _smiles = "n1occc1.c1ncoc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRing, TestWithinRingUnsaturation)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          within_ring_unsaturation: 2
          hits_needed: 2
        }
        ring_size: 5
        aromatic: true
      }
      smarts: "c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "n1occc1.c1ncoc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 6);

  EXPECT_TRUE(_DoPermutationTests(6));
}

TEST_F(TestSubstructureRing, TestLargestNumberOfBondsSharedWithAnotherRingIsolated)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          largest_number_of_bonds_shared_with_another_ring: 0
        }
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC1CC1.C1CC12CC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 9);

  EXPECT_TRUE(_DoPermutationTests(9));
}

TEST_F(TestSubstructureRing, TestLargestNumberOfBondsSharedWithAnotherRing1)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          largest_number_of_bonds_shared_with_another_ring: 1
        }
      }
      smarts: "c:c(:c):c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "n1ccc2ccccc12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 6);

  EXPECT_TRUE(_DoPermutationTests(6));
}

TEST_F(TestSubstructureRing, TestLargestNumberOfBondsSharedWithAnotherRing2)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          largest_number_of_bonds_shared_with_another_ring: 2
        }
      }
      smarts: "C[CD3](C)C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C3CC2CC(CC1C2)C3";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 24);

  EXPECT_TRUE(_DoPermutationTests(24));
}

TEST_F(TestSubstructureRing, TestLargestNumberOfBondsSharedWithAnotherRingNoMatch)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          largest_number_of_bonds_shared_with_another_ring: 1
        }
      }
      smarts: "C[CD3](C)C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C3CC2CC(CC1C2)C3";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRing, TestPiElectrionsZero)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          atoms_with_pi_electrons: 0
        }
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1.c1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPermutationTests(3));
}

TEST_F(TestSubstructureRing, TestPiElectrionsBenzene)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          atoms_with_pi_electrons: 6
        }
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1.c1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPermutationTests(3));
}

TEST_F(TestSubstructureRing, TestStronglyFusedRingNeighbours)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          strongly_fused_ring_neighbours: 2
        }
      }
      smarts: "[CD3]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1.c1ccccc1.C1C3CC2CC(CC1C2)C3";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);

  EXPECT_TRUE(_DoPermutationTests(4));
}

TEST_F(TestSubstructureRing, TestEnvironmentZeroMatchLess)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          environment: "<2c-F"
        }
      }
      smarts: "[cD2]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1.c1ccccc1.C1C3CC2CC(CC1C2)C3";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 6);

  EXPECT_TRUE(_DoPermutationTests(6));
}

TEST_F(TestSubstructureRing, TestEnvironmentOkMatchLess)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          environment: "<2c-F"
        }
      }
      smarts: "[cD3]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1.c1cc(F)ccc1.C1C3CC2CC(CC1C2)C3";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPermutationTests(1));
}

TEST_F(TestSubstructureRing, TestEnvironmentAnd)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          environment: "<2c-F&&1c-N"
        }
      }
      smarts: "[CD3]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1.c1cc(F)ccc1N.C1C3CC2CC(CC1C2)C3";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);

  EXPECT_TRUE(_DoPermutationTests(4));
}

TEST_F(TestSubstructureRing, TestEnvironmentOr)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          environment: "<2c-F||1c-N"
        }
      }
      smarts: "[CD3]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1.c1ccccc1N.C1C3CC2CC(CC1C2)C3";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);

  EXPECT_TRUE(_DoPermutationTests(4));
}

TEST_F(TestSubstructureRing, TestEnvironmentXor1)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          environment: "<2c-F^^1c-N"
        }
      }
      smarts: "[CD2]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  // Note that the 3 membered ring satisfies the environment specification.
  // Shows how this can be tricky.
  _smiles = "c1ccc(F)cc1N.C1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPermutationTests(3));
}

TEST_F(TestSubstructureRing, TestEnvironmentXor2)
{
  _string_proto = R"(
    query {
      ring_specifier {
        aromatic: true
        base {
          environment: "<2c-F^^1c-N"
        }
      }
      smarts: "[CD2]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  // Because the ring was specified as aromatic, the 3 membered ring can no longer match.
  _smiles = "c1ccc(F)cc1N.C1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRing, TestEnvironmentXor3)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          environment: "<2c-F^^1c-N"
        }
      }
      smarts: "[CD2]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ccc(F)cc1N.CCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRing, TestEnvironmentEnvironmentCanMatchInRingAtoms)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          environment_can_match_in_ring_atoms: true
          environment: "1ncN"
        }
      }
      smarts: "[ND1]c(:c):n"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Nc1ncccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPermutationTests(1));
}

TEST_F(TestSubstructureRing, TestEnvironmentEnvironmentCanMatchInRingAtomsX)
{
  _string_proto = R"(
    query {
      ring_specifier {
        base {
          environment_can_match_in_ring_atoms: true
          environment: "1ncN"
        }
      }
      smarts: "Nc"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  // Remember the ring matching is separate from the atom matching.
  _smiles = "Nc1ncccc1.Nc1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPermutationTests(2));
}

TEST_F(TestSubstructureRing, TestFusedNoRing)
{
  _string_proto = R"(
    query {
      ring_specifier {
        fused: 0
      }
      smarts: "*"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRing, TestFused0)
{
  _string_proto = R"(
    query {
      ring_specifier {
        fused: 0
      }
      smarts: "*"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPermutationTests(3));
}

TEST_F(TestSubstructureRing, TestFused2)
{
  _string_proto = R"(
    query {
      ring_specifier {
        fused: 2
        base {
          hits_needed: 1
        }
      }
      smarts: "*"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2CC12.C1C3CC2CC(CC1C2)C3";
  _smiles = "C1C2C3CC321";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 5);

  EXPECT_TRUE(_DoPermutationTests(5));
}

TEST_F(TestSubstructureRing, TestFusedAdamantane2)
{
  _string_proto = R"(
    query {
      ring_specifier {
        fused: 2
        base {
          hits_needed: 3
        }
      }
      smarts: "*"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C3CC2CC(CC1C2)C3";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 10);

  EXPECT_TRUE(_DoPermutationTests(10));
}

TEST_F(TestSubstructureRing, TestFused3)
{
  _string_proto = R"(
    query {
      ring_specifier {
        fused: 3
        ring_size: 4
        base {
          hits_needed: 1
        }
      }
      smarts: "[Cx4]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC2C34CCC3(CC4)C21";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPermutationTests(2));
}

TEST_F(TestSubstructureRing, TestFusedAromaticNeighborsNotFusedAtAll)
{
  _string_proto = R"(
    query {
      ring_specifier {
        fused_aromatic_neighbours: 1
      }
      smarts: "*"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRing, TestFusedAromaticNeighborsZero)
{
  _string_proto = R"(
    query {
      ring_specifier {
        fused_aromatic_neighbours: 0
      }
      smarts: "*"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPermutationTests(3));
}

TEST_F(TestSubstructureRing, TestFusedAromaticNeighborsMatches2)
{
  _string_proto = R"(
    query {
      ring_specifier {
        fused_aromatic_neighbours: 2
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1.c1ccc-2c(c1)-c3c2cccc3";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPermutationTests(3));
}

TEST_F(TestSubstructureRing, TestFusedAromaticNeighborsMatches1)
{
  _string_proto = R"(
    query {
      ring_specifier {
        fused_aromatic_neighbours: 1
      }
      smarts: "[Cx2H2D2T0G0AX4+0]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1.o1ccc2ccccc12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPermutationTests(3));
}

TEST_F(TestSubstructureRing, TestFusedNonAromaticNoRingsMin)
{
  _string_proto = R"(
    query {
      ring_specifier {
        min_fused_non_aromatic_neighbours: 1
      }
      smarts: "[*]";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRing, TestFusedNonAromaticNoRingsNeverMatches)
{
  _string_proto = R"(
    query {
      ring_specifier {
        fused_non_aromatic_neighbours: 0
      }
      smarts: "[*]";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRing, TestFusedNonAromaticRingsNeedToMatch)
{
  _string_proto = R"(
    query {
      ring_specifier {
        max_fused_non_aromatic_neighbours: 1
      }
      smarts: "[*]";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPermutationTests(3));
}

TEST_F(TestSubstructureRing, TestFusedNonAromaticNeighborsNoSpiro)
{
  _string_proto = R"(
    query {
      ring_specifier {
        fused_non_aromatic_neighbours: 1
      }
      smarts: "[aT2]";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N12N=CC=C1NCC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPermutationTests(1));
}

} // namespace
