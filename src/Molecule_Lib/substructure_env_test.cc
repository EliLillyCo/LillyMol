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

class TestSubstructureEnv : public testing::Test
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
TestSubstructureEnv::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

bool TestSubstructureEnv::_DoPerumationsTests(const int expected)
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

TEST_F(TestSubstructureEnv, TestEnvAttachmentPoint)
{
  _string_proto = R"(
    query {
      smarts: "C"
      environment {
        smarts: "F"
        attachment {
          attachment_point: 0
          btype: SS_SINGLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructureEnv, TestEnvAttachmentPointMultipleBonds)
{
  _string_proto = R"(
    query {
      smarts: "C"
      environment {
        smarts: "N"
        attachment {
          attachment_point: 0
          btype: SS_SINGLE_BOND
          btype: SS_DOUBLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCC=N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms matches = _FirstAtomEachEmbedding();

  EXPECT_THAT(matches, UnorderedElementsAre(1, 4));

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructureEnv, TestEnvAttachmentPointSubstructureBond)
{
  _string_proto = R"(
    query {
      smarts: "C"
      environment {
        smarts: "N"
        attachment {
          substructure_bond: "-,= 0"
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCC=N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms matches = _FirstAtomEachEmbedding();

  EXPECT_THAT(matches, UnorderedElementsAre(1, 4));

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructureEnv, TestEnvironmentMustMatchUnmatchedAtomsEasy)
{
  _string_proto = R"(
    query {
      smarts: "C"
      environment_must_match_unmatched_atoms: true
      environment {
        smarts: "N"
        attachment {
          substructure_bond: "-,= 0"
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCC=N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms matches = _FirstAtomEachEmbedding();

  EXPECT_THAT(matches, UnorderedElementsAre(1, 4));

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructureEnv, TestEnvironmentMustMatchUnmatchedAtomsNoMatch)
{
  _string_proto = R"(
    query {
      smarts: "C"
      environment_must_match_unmatched_atoms: true
      environment {
        smarts: "NC"
        attachment {
          substructure_bond: "-,= 0"
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCC=N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructureEnv, TestEnvironmentMustMatchUnmatchedAtomsMatchesAlreadyMatched)
{
  _string_proto = R"(
    query {
      smarts: "C"
      environment_must_match_unmatched_atoms: false
      environment {
        smarts: "N~C"
        attachment {
          substructure_bond: "-,= 0"
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCC=N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructureEnv, TestEnvMatchesShareAttachmentPoints)
{
  _string_proto = R"(
    query {
      smarts: "C"
      env_matches_can_share_attachment_points: false
      environment {
        smarts: "N"
        attachment {
          substructure_bond: "- 0"
        }
      }
      environment {
        smarts: "N"
        attachment {
          substructure_bond: "= 0"
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NC(=N)C(C#N)CC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructureEnv, TestEnvMatchesShareAttachmentPointsTrue)
{
  _string_proto = R"(
    query {
      smarts: "C"
      env_matches_can_share_attachment_points: true
      environment {
        smarts: "N"
        attachment {
          substructure_bond: "- 0"
        }
      }
      environment {
        smarts: "N"
        attachment {
          substructure_bond: "= 0"
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NC(=N)C(C#N)CC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructureEnv, TestHitsNeeded)
{
  _string_proto = R"(
    query {
      smarts: "C"
      env_matches_can_share_attachment_points: true
      environment {
        smarts: "N"
        hits_needed: 2
        attachment {
          substructure_bond: "-,= 0"
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NC(=N)C(C#N)CC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructureEnv, TestHitsNeededMultipleSites)
{
  _string_proto = R"(
    query {
      smarts: "c"
      env_matches_can_share_attachment_points: true
      environment {
        smarts: "F"
        hits_needed: 2
        attachment {
          attachment_point: 0
          btype: SS_SINGLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Fc1cc(F)ccNc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructureEnv, TestHitsNeededMultipleSitesCorrect)
{
  _string_proto = R"(
    query {
      smarts: "c1ccccc1"
      env_matches_can_share_attachment_points: true
      environment {
        smarts: "F"
        hits_needed: 2
        attachment {
          attachment_point: [0, 2]
          btype: SS_SINGLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Fc1cc(F)cc(N)c1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

// This time env_matches_can_share_attachment_points is in the environment.
TEST_F(TestSubstructureEnv, TestHitsNeededMultipleSitesCorrectEnv)
{
  _string_proto = R"(
    query {
      smarts: "c1ccccc1"
      environment {
        env_matches_can_share_attachment_points: true
        smarts: "F"
        hits_needed: 2
        attachment {
          attachment_point: [0, 2]
          btype: SS_SINGLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Fc1cc(F)cc(N)c1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructureEnv, TestNoOtherSubstituentsAllowed)
{
  _string_proto = R"(
    query {
      smarts: "C1CCOCC1"
      environment {
        no_other_substituents_allowed: true
        smarts: "N"
        hits_needed: 1
        attachment {
          attachment_point: 0
          btype: SS_SINGLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NC(OC)1CCOCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructureEnv, TestNoOtherSubstituentsAllowedEnvMatchesOK)
{
  _string_proto = R"(
    query {
      smarts: "C1CCOCC1"
      environment {
        no_other_substituents_allowed: true
        smarts: "N"
        attachment {
          attachment_point: 0
          btype: SS_SINGLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  // The extra substituent also matches the environment.
  _smiles = "NC(NC)1CCOCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructureEnv, TestHydrogenOnlyH)
{
  _string_proto = R"(
    query {
      smarts: "C1CCOCC1"
      environment {
        hydrogen_ok: true
        smarts: "N"
        hits_needed: 1
        attachment {
          attachment_point: 0
          btype: SS_SINGLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CCOCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructureEnv, TestHydrogenNotHydrogen)
{
  _string_proto = R"(
    query {
      smarts: "C1CCOCC1"
      environment {
        hydrogen_ok: true
        smarts: "N"
        hits_needed: 1
        attachment {
          attachment_point: 0
          btype: SS_SINGLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  // Rather than an environment match, there is a non-H that matches
  _smiles = "CC1CCOCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructureEnv, Testmax_env_matches_per_anchor)
{
  _string_proto = R"(
    query {
      smarts: "CCC"
      environment {
        smarts: "F"
        max_env_matches_per_anchor: 1
        hits_needed: 2
        attachment {
          attachment_point: [0, 2]
          btype: SS_SINGLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)CCF";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructureEnv, TestMultipleAnd)
{
  _string_proto = R"(
    query {
      smarts: "CCC"
      environment {
        and_id: 1
        smarts: "F"
        attachment {
          attachment_point: [0, 2]
          btype: SS_SINGLE_BOND
        }
      }
      environment {
        and_id: 1
        attachment {
          attachment_point: 1
          btype: SS_SINGLE_BOND
        }
        query_atom {
          id: 3
          atom_properties {
            atomic_number : 6
            ncon: 1
          }
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)C(C)CF";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructureEnv, TestMultipleOr)
{
  _string_proto = R"(
    query {
      smarts: "c1cnccc1"
      environment {
        or_id: 1
        smarts: "F"
        attachment {
          attachment_point: [0, 1, 2, 3]
          btype: SS_SINGLE_BOND
        }
      }
      environment {
        or_id: 1
        attachment {
          attachment_point: [0, 1, 2, 3, 4, 5]
          btype: SS_SINGLE_BOND
        }
        query_atom {
          id: 6
          atom_properties {
            atomic_number : 8
            ncon: 1
          }
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Fc1cnccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPerumationsTests(1));

  _smiles = "Oc1cnccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructureEnv, TestRejection)
{
  _string_proto = R"(
    query {
      smarts: "c1cnccc1"
      environment_no_match {
        smarts: "F"
        attachment {
          attachment_point: [0, 1, 2, 3, 4, 5]
          btype: SS_SINGLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Fc1cnccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));

  _smiles = "Oc1cnccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructureEnv, TestNumericEnv1) {
  _string_proto = R"(
    query {
      smarts: "c1cnccc1"
      environment {
        smarts: ">1F"
        attachment {
          attachment_point: [0, 1, 2, 3, 4, 5]
          btype: SS_SINGLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Fc1cnccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));

  _smiles = "Fc1cnc(F)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructureEnv, TestNumericEnv2) {
  _string_proto = R"(
    query {
      smarts: "c1cnccc1"
      environment {
        smarts: "{-2}F"
        attachment {
          attachment_point: [0, 1, 2, 3, 4, 5]
          btype: SS_SINGLE_BOND
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Fc1cnccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));

  _smiles = "Fc1cnc(F)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_TRUE(_DoPerumationsTests(2));

  _smiles = "Fc1cnc(F)c(F)c1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPerumationsTests(0));

}

}  // namespace
