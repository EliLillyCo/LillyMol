#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "aromatic.h"
#include "substructure.h"

namespace {

//using google::protobuf::TextFormat::ParseFromString;

using testing::UnorderedElementsAre;
using testing::ElementsAre;

class TestSubstructure : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    std::string _string_proto;

    IWString _smiles;

    Substructure_Query _query;

    Substructure_Results _sresults;

    Molecule _m;

  protected:
    void _WriteQuery(const char * fname) {
      IWString tmp(fname);
      _query.write_msi(tmp);
    }
};

void
TestSubstructure::SetUp()
{
  set_global_aromaticity_type(Daylight);
  reset_aromatic_file_scope_variables();
}

TEST_F(TestSubstructure, SingleAtom)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  EXPECT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  ASSERT_TRUE(_m.build_from_smiles("C"));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

}

TEST_F(TestSubstructure, SingleAtomMultipleMatches)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  EXPECT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  for (int i = 0; i < 20; ++i) {
    _smiles << "C";
    ASSERT_TRUE(_m.build_from_smiles(_smiles));
    EXPECT_EQ(_query.substructure_search(_m, _sresults), i + 1);
  }
}

TEST_F(TestSubstructure, SingleAtomNot)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        match_as_match: false
        atom_properties {
          atomic_number : 6
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  EXPECT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestComposite)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
    query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  EXPECT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestCompositeAnd)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
    logexp: SS_AND
    query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  EXPECT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
  _smiles = "NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestCompositeLowPriorityAnd)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
    logexp: SS_LP_AND
    query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  EXPECT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
  _smiles = "NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestCompositeXor)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
    logexp: SS_XOR
    query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "CN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestCompositeEachComponentMatch)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_symbol : "C"
        }
      }
    }
    match_each_component: 1
    query {
      query_atom {
        id: 0
        atom_properties {
          atomic_symbol : "N"
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "CN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
}

TEST_F(TestSubstructure, TestRespectInitialAtomNumbering)
{
  _string_proto = R"(query {
      respect_initial_atom_numbering: true
      query_atom {
        id: 0
        initial_atom_number: 2
        atom_properties {
          atomic_symbol : "Na"
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "[Na]";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 3);

  EXPECT_EQ(INVALID_ATOM_NUMBER, e->item(0));
  EXPECT_EQ(INVALID_ATOM_NUMBER, e->item(1));
  EXPECT_EQ(0, e->item(2));
}

TEST_F(TestSubstructure, TestMatchAsMatchNoMatch)
{
  _string_proto = R"(query {
      query_atom {
        match_as_match: false
        id: 0
        atom_properties {
          atomic_symbol : "I"
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "[I]";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestMatchAsMatchMatch)
{
  _string_proto = R"(query {
      query_atom {
        match_as_match: false
        id: 0
        atom_properties {
          atomic_symbol : "I"
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "[I]C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 1);
  EXPECT_EQ(e->item(0), 1);

  _smiles = "C[I]";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 1);
  EXPECT_EQ(e->item(0), 0);
}

TEST_F(TestSubstructure, TestSymmetricMatches)
{
  _string_proto = R"(query {
      perceive_symmetric_equivalents: false
      query_atom {
        id: 0
        atom_properties {
          atomic_symbol : "C"
        }
      }
      query_atom {
        id: 3
        single_bond: 0
        atom_properties {
          atomic_number: 6
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  _query.set_perceive_symmetry_equivalent_matches(true);
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
}

TEST_F(TestSubstructure, TestUniqueEmbeddingsOnly)
{
  _string_proto = R"(query {
      unique_embeddings_only: true
      smarts: "C(F)(F)(F)F"
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FC(F)(F)F";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 5);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
}

TEST_F(TestSubstructure, TestIncludeAtomInEmbeddingSmarts)
{
  // When build from smarts, respect_initial_atom_numbering is in effect.
  _string_proto = R"(query {
      unique_embeddings_only: true
      smarts: "[/IWxC](F)(F)(F)F"
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FC(F)(F)F";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 5);
  EXPECT_EQ(INVALID_ATOM_NUMBER, e->item(0));
  EXPECT_EQ(_m.atomic_number(e->item(1)), 9);
  _query.set_respect_initial_atom_numbering(0);
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 4);
}

TEST_F(TestSubstructure, TestAllVariantsMatch)
{
  // When build from smarts, respect_initial_atom_numbering is in effect.
  _string_proto = R"(query {
      smarts: "*(F)(F)(F)F"
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FC(F)(F)F";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 24);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 5);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
}

TEST_F(TestSubstructure, TestOneEmebeddingPerStartAtom)
{
  // When build from smarts, respect_initial_atom_numbering is in effect.
  _string_proto = R"(query {
      unique_embeddings_only: true
      smarts: "c1ccccc1"
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 6);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
}

TEST_F(TestSubstructure, TestIncludeAtomInEmbeddingProto)
{
  _string_proto = R"(query {
      unique_embeddings_only: true
      query_atom {
        id: 0
        include_in_embedding: false
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 9
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 9
        }
        single_bond: 0
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 9
        }
        single_bond: 0
      }
      query_atom {
        id: 4
        atom_properties {
          atomic_number: 9
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FC(F)(F)F";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 4);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 9);
}

#ifdef OR_ID_NOT_IMPLEMENTED
TEST_F(TestSubstructure, TestAtomOrId)
{
  _string_proto = R"(query {
      unique_embeddings_only: true
      query_atom {
        id: 0
        or_id: 75
        atom_properties {
          atomic_symbol : "F"
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FC(F)(F)F";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
  for (int i = 0; i < 4; ++i) {
    const Set_of_Atoms * e = _sresults.embedding(0);
    EXPECT_EQ(e->number_elements(), 1);
  }
}
#endif

TEST_F(TestSubstructure, TestInitialAtomNumberSparse)
{
  _string_proto = R"(query {
      respect_initial_atom_numbering: true
      query_atom {
        id: 0
        initial_atom_number: 8
        atom_properties {
          atomic_symbol : "O"
        }
      }
      query_atom {
        id: 1
        initial_atom_number: 2
        atom_properties {
          atomic_symbol : "C"
        }
        double_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "NC(=O)C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 9);
  EXPECT_EQ(e->item(0), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(1), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(2), 1);
  EXPECT_EQ(e->item(3), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(4), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(5), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(6), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(7), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(8), 2);
}

TEST_F(TestSubstructure, TestRingIdMatches)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        ring_id: 8
        atom_properties {
          atomic_symbol: "C"
          aromatic: true
        }
      }
      query_atom {
        id: 1
        ring_id: 8
        atom_properties {
          atomic_symbol : "N"
          aromatic: true
        }
        aromatic_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1ncccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 2);
}

TEST_F(TestSubstructure, TestRingIdNotMatch)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        ring_id: 2
        atom_properties {
          atomic_symbol: "C"
          aromatic: true
        }
      }
      query_atom {
        id: 1
        ring_id: 8
        atom_properties {
          atomic_symbol : "N"
          aromatic: true
        }
        aromatic_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1ncccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestFusedSystemIdMatch)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        fused_system_id: 3
        atom_properties {
          atomic_symbol: "C"
          aromatic: true
        }
      }
      query_atom {
        id: 1
        fused_system_id: 3
        atom_properties {
          atomic_symbol : "N"
          aromatic: true
        }
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_symbol : "C"
          aromatic: false
        }
        single_bond: 1
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1nnc2cn(C)ccc12";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 6);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 3);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
}

TEST_F(TestSubstructure, TestFusedSystemIdNoMatch)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        fused_system_id: 3
        atom_properties {
          atomic_symbol: "C"
          aromatic: true
        }
      }
      query_atom {
        id: 1
        fused_system_id: 9
        atom_properties {
          atomic_symbol : "N"
          aromatic: true
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1ccccc1c1ncccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 6);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
  EXPECT_EQ(_m.atomic_number(e->item(1)), 7);
}

TEST_F(TestSubstructure, TestFragmentIdMatchesSameFrag)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        fragment_id: 3
        atom_properties {
          atomic_number: 7
          aromatic: false
        }
      }
      query_atom {
        id: 1
        fragment_id: 3
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestFragmentIdNoMatches)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        fragment_id: 3
        atom_properties {
          atomic_number: 7
          aromatic: false
        }
      }
      query_atom {
        id: 1
        fragment_id: 9
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestFragmentIdMatchesDifferentFrags)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        fragment_id: 3
        atom_properties {
          atomic_number: 7
          aromatic: false
        }
      }
      query_atom {
        id: 1
        fragment_id: 9
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C.N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestBondSmartsSingle)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
          aromatic: false
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
        bond_smarts: "- 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestBondSmartsMultipleAttach)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
          aromatic: false
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
        bond_smarts: "- 0"
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
        bond_smarts: "- 0 1"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
}

TEST_F(TestSubstructure, TestBondSmartsDouble)
{
  _string_proto = R"(query {
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
        bond_smarts: "= 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C=NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestBondSmartsTriple)
{
  _string_proto = R"(query {
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
        bond_smarts: "# 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C#NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestBondSmartsArom)
{
  _string_proto = R"(query {
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
        bond_smarts: ": 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1nccnc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
}

TEST_F(TestSubstructure, TestBondSmartsAny)
{
  _string_proto = R"(query {
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
        bond_smarts: "~ 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C-N.C=N.C#N.c1ncccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 5);
}

TEST_F(TestSubstructure, TestBondSmartsRing)
{
  _string_proto = R"(query {
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
        bond_smarts: "@ 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C-N.C1NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
}


TEST_F(TestSubstructure, TestBondSmartsDoubleAndRing)
{
  _string_proto = R"(query {
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
        bond_smarts: "@= 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C=N.C1C=NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestBondSmartsDoubleNotRing)
{
  _string_proto = R"(query {
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
        bond_smarts: "!@= 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C=N=C.C1C=NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
}

TEST_F(TestSubstructure, TestBondSmartsSingleOrDouble)
{
  _string_proto = R"(query {
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
        bond_smarts: "-,= 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C=N-C.C1C=NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
}

TEST_F(TestSubstructure, TestBondSmartsNonRing)
{
  _string_proto = R"(query {
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
        bond_smarts: "!@ 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();


  _smiles = "C-N.C1NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestBondSmartsXor)
{
  _string_proto = R"(query {
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
        bond_smarts: "@^= 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1NC1.C1C=NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);
}

TEST_F(TestSubstructure, TestBondSmartsSingleDoubleAndRing)
{
  _string_proto = R"(query {
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
        bond_smarts: "-,=;@ 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1NC1.C1C=NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
}

TEST_F(TestSubstructure, TestQueryBondSimple)
{
  _string_proto = R"(query {
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
        query_bond {
          btype: SS_SINGLE_BOND
          other_end: 0
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1NC1.C1C=NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);
}

TEST_F(TestSubstructure, TestQueryBondComplex)
{
  _string_proto = R"(query {
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
        query_bond {
          btype: SS_SINGLE_BOND
          btype: SS_TRIPLE_BOND
          other_end: 0
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "NC#N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
}

TEST_F(TestSubstructure, TestPreferenceDefault)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
        preference {
          ncon: 2
          preference_value: 5
        }
        preference {
          ncon: 1
          preference_value: 1
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  // By default, low preference hit are removed.

  _smiles = "NC.NCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(_m.ncon(e->item(0)), 2);

  _smiles = "NCC.NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  e = _sresults.embedding(0);
  EXPECT_EQ(_m.ncon(e->item(0)), 2);
}

TEST_F(TestSubstructure, TestPreferenceSortPreferences)
{
  _string_proto = R"(query {
      sort_by_preference_value: true
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
        preference {
          ncon: 2
          preference_value: 5
        }
        preference {
          ncon: 1
          preference_value: 1
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  // By default, low preference hit are removed.

  _smiles = "NC.NCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(_m.ncon(e->item(0)), 2);

  _smiles = "NCC.NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
  e = _sresults.embedding(0);
  EXPECT_EQ(_m.ncon(e->item(0)), 2);
}

TEST_F(TestSubstructure, TestPreferenceSumPreferences)
{
  _string_proto = R"(query {
      sort_by_preference_value: true
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
        preference {
          ncon: 2
          preference_value: 1
        }
        preference {
          ncon: 1
          preference_value: 5
        }
        preference {
          ring_bond_count: 2
          preference_value: 6
        }
        sum_all_preference_hits: true
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "NC.NCC.N1CC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(_m.ncon(e->item(0)), 2);
  EXPECT_EQ(_m.ring_bond_count(e->item(0)), 2);

  _smiles = "N1CC1.NCC.NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
  e = _sresults.embedding(0);
  EXPECT_EQ(_m.ncon(e->item(0)), 2);
  EXPECT_EQ(_m.ring_bond_count(e->item(0)), 2);
}

TEST_F(TestSubstructure, TestElementsNeededMatches)
{
  _string_proto = R"(query {
      required_molecular_properties {
        elements_needed {
          atomic_number: 9
          hits_needed: 1
        }
        elements_needed {
          atomic_number: 6
          min_hits_needed: 1
        }
        elements_needed {
          atomic_number: 7
          max_hits_needed: 1
        }
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FCN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestElementsNeededNoMatches)
{
  _string_proto = R"(query {
      required_molecular_properties {
        elements_needed {
          atomic_number: 9
          hits_needed: 2
        }
        elements_needed {
          atomic_number: 6
          min_hits_needed: 1
        }
        elements_needed {
          atomic_number: 7
          max_hits_needed: 1
        }
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FCN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestElementHitsNeededMatches)
{
  _string_proto = R"(query {
      required_molecular_properties {
        elements_needed {
          atomic_number: 6
          min_hits_needed: 1
        }
        elements_needed {
          atomic_number: 7
          max_hits_needed: 1
        }
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
      element_hits_needed {
        atomic_number: 9
        hits_needed: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FCN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestElementHitsNeededNoMatches)
{
  _string_proto = R"(query {
      required_molecular_properties {
        elements_needed {
          atomic_number: 6
          min_hits_needed: 1
        }
        elements_needed {
          atomic_number: 7
          max_hits_needed: 1
        }
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
      element_hits_needed {
        atomic_number: 9
        hits_needed: 1
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FCN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestChiral4Matches)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
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
          atomic_number: 7
        }
        single_bond: 1
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 9
        }
        single_bond: 1
      }
      query_atom {
        id: 4
        atom_properties {
          atomic_number: 6
        }
        single_bond: 1
      }
      query_atom {
        id: 5
        atom_properties {
          atomic_number: 6
        }
        single_bond: 4
      }
      chiral_centre {
        center: 1
        top_front {
          atom_number: 0
        }
        top_back {
          atom_number: 2
        }
        left_down {
          atom_number: 3
        }
        right_down {
          atom_number: 4
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C[C@@](F)(N)CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "C[C@](F)(N)CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestChiral3Matches)
{
  _string_proto = R"(query {
      smarts: "C[CH](N)CC";
      chiral_centre {
        center: 1
        top_front {
          atom_number: 0
        }
        top_back {
          atom_number: 2
        }
        left_down {
          h_or_lp: HYDROGEN
        }
        right_down {
          atom_number: 3
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C[C@@H](N)CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "C[C@H](N)CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}


TEST_F(TestSubstructure, TestAnyLength)
{
  _string_proto = R"(query {
      smarts: "C#[#{hello}D2]=[#{world}D>1]";
    }
  )";

  set_auto_create_new_elements(1);
  set_atomic_symbols_can_have_arbitrary_length(1);

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "[Cr][world]=[hello]#CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  set_auto_create_new_elements(0);
  set_atomic_symbols_can_have_arbitrary_length(0);
}

TEST_F(TestSubstructure, TestLowerCaseUppercaseElement) {
  set_auto_create_new_elements(1);
  set_atomic_symbols_can_have_arbitrary_length(1);

  ASSERT_TRUE(_query.create_from_smarts("[#{qQ}]"));
  ASSERT_TRUE(_m.build_from_smiles("[Qq][qQ][qQ]"));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  set_auto_create_new_elements(0);
  set_atomic_symbols_can_have_arbitrary_length(0);
}

TEST_F(TestSubstructure, TestCompsite)
{
  _string_proto = R"(query {
      smarts: "CC";
    }
    query {
      smarts: "NN";
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1cccnc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  _smiles = "NN.CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
  EXPECT_EQ(_sresults.number_embeddings(), 2);

  for (const auto* e : _sresults.embeddings())
  {
    EXPECT_EQ(e->number_elements(), 2);
    EXPECT_THAT(*e, UnorderedElementsAre(2, 3));
  }

  _smiles = "CC.c1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
  EXPECT_EQ(_sresults.number_embeddings(), 2);

  for (const auto* e : _sresults.embeddings())
  {
    EXPECT_EQ(e->number_elements(), 2);
    EXPECT_THAT(*e, UnorderedElementsAre(0, 1));
  }
}

TEST_F(TestSubstructure, TestCompsiteXor)
{
  _string_proto = R"(query {
      smarts: "CC";
    }
    logexp: SS_XOR
    query {
      smarts: "NN";
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1cccnc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  _smiles = "NN.CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  _smiles = "NN.C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
  EXPECT_EQ(_sresults.number_embeddings(), 2);

  for (const auto* e : _sresults.embeddings())
  {
    EXPECT_EQ(e->number_elements(), 2);
    EXPECT_THAT(*e, UnorderedElementsAre(0, 1));
  }

  _smiles = "CC.NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
  // A surprising result, no embeddings. This is correct because
  // the last query done, looking for NN, failed. But at that point
  // the result of the XOR was available. But at that time,
  // _sresults was holding a non match.
  EXPECT_EQ(_sresults.number_embeddings(), 0);
}

TEST_F(TestSubstructure, TestManyMatches)
{
  _string_proto = R"(query {
      smarts: "CC(C)(C)c1ccc(C(C)(C)C)cc1";
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CC(C)(C)c1ccc(C(C)(C)C)cc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 144);
}

TEST_F(TestSubstructure, TestSameEmbedding)
{
  _string_proto = R"(query {
      smarts: "[R]";
      compress_embeddings: true
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CC(C)(C)c1ccc(C(C)(C)C)cc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(4, 5, 6, 7, 12, 13));
}

TEST_F(TestSubstructure, TestUnmatchedAtomsAttachedZero)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
          ncon: 2
          aromatic: true
        }
      }
      query_atom {
        id: 1
          atom_properties {
          atomic_number: 6
          ncon: 3
          aromatic: true
        }
        aromatic_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 6
          ncon: 2
          aromatic: true
        }
        aromatic_bond: 1
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
        unmatched_atoms_attached: 0
        single_bond: 1
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1ccccc1C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0, 4, 5, 6));
}

TEST_F(TestSubstructure, TestUnmatchedAtomsAttachedOne)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
          ncon: 2
          aromatic: true
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
          ncon: 3
          aromatic: true
        }
        aromatic_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 6
          ncon: 2
          aromatic: true
        }
        aromatic_bond: 1
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
        unmatched_atoms_attached: 1
        single_bond: 1
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1ccccc1CC";  // Has an unmatched atom attached.
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0, 4, 5, 6));
}

TEST_F(TestSubstructure, TestUnmatchedAtomsAttachedMinOne)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
          ncon: 2
          aromatic: true
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
          ncon: 3
          aromatic: true
        }
        aromatic_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 6
          ncon: 2
          aromatic: true
        }
        aromatic_bond: 1
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
        min_unmatched_atoms_attached: 1
        single_bond: 1
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1ccccc1CC";  // Has an unmatched atom attached.
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0, 4, 5, 6));
}

TEST_F(TestSubstructure, TestUnmatchedAtomsAttachedMaxOne)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
          ncon: 2
          aromatic: true
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
          ncon: 3
          aromatic: true
        }
        aromatic_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 6
          ncon: 2
          aromatic: true
        }
        aromatic_bond: 1
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
        max_unmatched_atoms_attached: 1
        single_bond: 1
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1ccccc1CC";  // Has an unmatched atom attached.
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0, 4, 5, 6));
}

TEST_F(TestSubstructure, TestUnmatchedAtomsAttachedRing)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 8
          ncon: 2
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
          ncon: 2
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 7
          ncon: 2
        }
        unmatched_atoms_attached: 5
        single_bond: 1
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "O1CNCC(C)C1C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0, 1, 2));
}

TEST_F(TestSubstructure, TestUnmatchedAtomsAttachedLogicalExpOr)
{
  // C-N-
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
          ncon: 1
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
          ncon: 2
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 7
        }
        atom_properties {
          atomic_number: 6
          logical_operator: SS_OR
        }
        min_unmatched_atoms_attached: 2
        single_bond: 1
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CNC";   // 0 unmatched atoms attached.
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  _smiles = "CNCC";   // 1 unmatched atoms attached.
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
  return;

  _smiles = "CNCCO";   // 2 unmatched atoms attached, C
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0, 1, 2));

  _smiles = "CNNCO";   // 2 unmatched atoms attached, N
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0, 1, 2));
}

TEST_F(TestSubstructure, TestRecursiveMatchesNotMatched)
{
  // C-N-
  _string_proto = R"(query {
      smarts: "C[$(C(=O)[!C])]"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CCC(=O)N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

// By default, atom environments can match already matched atoms.
// Test turning that feature off.
TEST_F(TestSubstructure, TestRecursiveMatchesAlreadyMatched)
{
  // C-N-
  _string_proto = R"(query {
      smarts: "C[$(C(=O)C)]"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CCC(=O)N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  set_atom_environment_only_matches_unmatched_atoms(1);
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestRecursiveBenzene)
{
  _string_proto = R"(query {
    smarts: "[$([$([acr6]%64[$(c%25[cR1][0#6H1r6R1D2X3x2+0a][$([cH1a]%51[0#6+0a][0#6H1r6R1D2X3x2+0a]c[#6aH][$(cc)]%51)][#6aH][cX3H]%25)][x2Hc][c^n][ch][0#6H1r6R1D2X3x2+0a]%64)]8[$([x2Hc]%82[#6aH][$([0#6H1r6R1D2X3x2+0a]%65c[$(c9[ac][cR1][x2Hc][$(cc)][acr6]9)][$([$(cc)]%97[x2Hc]c[c^n][cX3H][cR1]%97)][ch][$(cc)]%65)][$(c%97[acr6][$([$(cc)]%70[0#6H1r6R1D2X3x2+0a][ac][x2Hc][acr6]c%70)]c[0#6H1r6R1D2X3x2+0a][$([acr6]%77[0#6+0a][ch][$(cc)][cH1a][cX3H]%77)]%97)][$([cH1a]%41[0#6+0a][ac][$([cR1]%54[$(cc)][$([ch]%87[$([#6aH]%86c[c^n][cH1a][acr6][$([$([cR1]%90[0#6+0a][c^n][ac]c[x2Hc]%90)]%45[ac][cR1][$([ch]%41[0#6+0a][cH1a][cR1]c[#6aH]%41)][$([cX3H]%72[cH1a][x2Hc][$(c%26[$(cc)][ch][cH1a][$([ac]8[#6aH][x2Hc][c^n][cH1a]c8)][0#6+0a]%26)][$([x2Hc]%86[cX3H][$(cc)]c[c^n][0#6H1r6R1D2X3x2+0a]%86)][$(cc)]%72)][#6aH]%45)]%86)][cX3H][c^n][#6aH][$([ch]%62[$(cc)]c[cX3H][ac][0#6+0a]%62)]%87)][cH1a][acr6][ch]%54)][#6aH][$([$([x2Hc]%39[cH1a][ch][0#6H1r6R1D2X3x2+0a][ac][$([ac]%31[c^n][0#6+0a][#6aH][acr6][ch]%31)]%39)]%61[$([0#6H1r6R1D2X3x2+0a]%74[cX3H][acr6][c^n][0#6+0a][ch]%74)][$([$(cc)]%29[cR1][cH1a]c[0#6+0a][0#6H1r6R1D2X3x2+0a]%29)]c[ch]c%61)]%41)][$([$(cc)]4[0#6+0a][$([$([$(c%48[ch][ac][cX3H][$(cc)][0#6+0a]%48)]%29[cH1a]c[c^n][x2Hc][0#6+0a]%29)]9c[acr6][ac][ch][$([#6aH]%12[cX3H]cc[cH1a][ac]%12)]9)][acr6][x2Hc][ac]4)]%82)][$([$([$([ch]%83[$([$([0#6+0a]%24[ac][cH1a]c[#6aH][c^n]%24)]%67[x2Hc]c[0#6+0a][$(cc)][0#6H1r6R1D2X3x2+0a]%67)][$([cR1]%45[$(cc)]c[cH1a][ch][c^n]%45)]cc[0#6+0a]%83)]%23[0#6+0a][$([0#6H1r6R1D2X3x2+0a]%70[cH1a][x2Hc][ch][0#6+0a][$([ch]%18[cH1a][cX3H][c^n][ac]c%18)]%70)][$([ac]%15[acr6][$([cX3H]%40c[c^n][$(c%55[cX3H][cH1a]c[c^n][ch]%55)][ch][0#6+0a]%40)][$([cR1]%28c[ch][c^n][0#6+0a][0#6H1r6R1D2X3x2+0a]%28)][ch][cR1]%15)][$(cc)][$([c^n]%43[#6aH][0#6+0a][cX3H][cH1a][$([acr6]%87c[#6aH][cH1a]c[0#6H1r6R1D2X3x2+0a]%87)]%43)]%23)]7[$([$(cc)]%91[$([0#6+0a]%86[c^n][$([c^n]%90[cX3H][$([0#6+0a]5[$(cc)][acr6][ac][cH1a]c5)][#6aH][x2Hc][0#6+0a]%90)][cH1a][cX3H][x2Hc]%86)]c[acr6][0#6H1r6R1D2X3x2+0a][$([c^n]9[0#6+0a][ch][ac]c[x2Hc]9)]%91)][cR1][$(cc)][$(c%97[$([acr6]%96[ac][0#6H1r6R1D2X3x2+0a]c[cX3H][#6aH]%96)][$([$(cc)]%63[ac]c[#6aH][x2Hc][cR1]%63)][$([0#6H1r6R1D2X3x2+0a]%57[x2Hc][ch]c[c^n][ac]%57)][acr6][$([x2Hc]%76[acr6]c[$([cH1a]%59[x2Hc][$(cc)][0#6H1r6R1D2X3x2+0a][cR1]c%59)][ch][cR1]%76)]%97)][$([0#6H1r6R1D2X3x2+0a]%14[x2Hc][ac][cX3H][c^n][cH1a]%14)]7)][$([$([$(c%73[cX3H][$([$([acr6]%68c[x2Hc][#6aH][ac][ch]%68)]%64[cX3H][cH1a][#6aH][$(cc)][0#6+0a]%64)][acr6][ac][$([x2Hc]%70[$(cc)]c[ac][cX3H][cH1a]%70)]%73)]%48[$([c^n]%79[$([$([ac]9[$(cc)][cX3H][0#6+0a]cc9)]%10c[0#6+0a][0#6H1r6R1D2X3x2+0a][x2Hc][#6aH]%10)][ch][cH1a][ac][#6aH]%79)][$(c%31[0#6+0a][c^n][$([ch]%71[$(cc)][#6aH]c[cR1][cH1a]%71)]c[$(cc)]%31)][$([#6aH]%56[cH1a]c[0#6+0a][cR1]c%56)][c^n][cX3H]%48)]%32[$([$(cc)]%95[$([$(cc)]%93[#6aH]c[$([cH1a]%74[acr6][$(cc)][#6aH]c[cR1]%74)][$([0#6H1r6R1D2X3x2+0a]%84[#6aH][0#6+0a][$(cc)][$([0#6+0a]%58[ac][c^n][cR1][#6aH]c%58)][ac]%84)][c^n]%93)][c^n][ac][$([cH1a]%78[$([acr6]%11[#6aH][0#6H1r6R1D2X3x2+0a][c^n][cH1a][ch]%11)][x2Hc][ac]c[$(cc)]%78)][cR1]%95)][$([0#6H1r6R1D2X3x2+0a]%54[x2Hc][ac][$([x2Hc]%78[acr6][$([c^n]%30[x2Hc][$([x2Hc]%50[#6aH][$(cc)][$([c^n]%48[#6aH][ac][x2Hc][$([ch]%82[cH1a][cR1][#6aH]c[acr6]%82)][ch]%48)][cH1a][ac]%50)][$(cc)][ac][#6aH]%30)][$([x2Hc]%45[$([ac]%71[#6aH][cH1a]c[ch][acr6]%71)][ch][cH1a][cR1][cX3H]%45)][$([0#6+0a]%51[cH1a]c[x2Hc][cX3H][0#6H1r6R1D2X3x2+0a]%51)][0#6H1r6R1D2X3x2+0a]%78)][$([cH1a]%42c[c^n][cX3H][$(cc)][x2Hc]%42)][$([$([$([$([$([$(c%38[cH1a][#6aH][cX3H][ac][$(cc)]%38)]%17[x2Hc][cX3H]c[$([#6aH]%92[0#6H1r6R1D2X3x2+0a]c[x2Hc][cX3H][cR1]%92)][#6aH]%17)]%23[$([#6aH]%28[cH1a]c[$([c^n]%59[ac][cR1][x2Hc][#6aH][acr6]%59)]c[c^n]%28)][ac][cX3H][cR1]c%23)]%28[cR1][$([ch]%51[#6aH]c[$(cc)][cR1][cH1a]%51)][$([$(cc)]%10[cH1a][cX3H][ch]c[$([cR1]%63[#6aH][cH1a][$(cc)][x2Hc][c^n]%63)]%10)][$([$([0#6+0a]%67[$(cc)][$([cH1a]%76[0#6H1r6R1D2X3x2+0a][cX3H][$(cc)][0#6+0a][ch]%76)][cH1a][$([0#6+0a]%90[acr6][cR1][ch][#6aH][0#6H1r6R1D2X3x2+0a]%90)]c%67)]%89c[ac][cH1a][$(cc)][ch]%89)][0#6H1r6R1D2X3x2+0a]%28)]%28[ch][#6aH][cH1a][$([cX3H]%58[$([cH1a]%12[ac][0#6H1r6R1D2X3x2+0a][cX3H]c[#6aH]%12)][0#6H1r6R1D2X3x2+0a][$([$(c%36[ch][$(cc)][cH1a][0#6H1r6R1D2X3x2+0a][acr6]%36)]%89[$(cc)][acr6][0#6+0a][x2Hc][cR1]%89)][#6aH][cR1]%58)][ac]%28)]%94c[$(c%12[ch][$([$(cc)]%15[$(c%97[0#6+0a][x2Hc][cR1][#6aH][$([0#6+0a]%94[c^n][#6aH][cH1a][acr6][ac]%94)]%97)][cX3H][$([$(cc)]%15c[$([0#6+0a]%93[cR1][x2Hc][acr6][ac][cH1a]%93)][#6aH][$([0#6H1r6R1D2X3x2+0a]%23[c^n][$([cX3H]%36[ch][cH1a][c^n][cR1]c%36)]c[#6aH][cH1a]%23)][$([x2Hc]%49[$([c^n]%23[acr6][ch]cc[cX3H]%23)][cR1][0#6+0a][cX3H][$(cc)]%49)]%15)][c^n][ch]%15)][cX3H][$([cR1]%15[ch][$([ac]%48[x2Hc][0#6H1r6R1D2X3x2+0a][$(cc)][c^n][cR1]%48)][acr6][cH1a][ac]%15)][$(c%50[ac][acr6][cH1a][0#6H1r6R1D2X3x2+0a][cR1]%50)]%12)][$([$([x2Hc]%15c[cH1a][$([$(cc)]%22[x2Hc]c[ac][cR1][cX3H]%22)][0#6+0a][cR1]%15)]%60c[0#6H1r6R1D2X3x2+0a]c[$(cc)][$([$([$(cc)]%98[$([ch]%24[0#6+0a][$([acr6]%24[$(cc)][cH1a][0#6+0a]c[#6aH]%24)][0#6H1r6R1D2X3x2+0a][cR1][$([$([cH1a]%23[ch][cR1][$([c^n]%25[x2Hc][$(cc)][ch][ac][cX3H]%25)][$([x2Hc]%35[ac][0#6+0a][acr6][cH1a][c^n]%35)]c%23)]%70[0#6H1r6R1D2X3x2+0a][ac][ch][$(cc)][acr6]%70)]%24)][cX3H][$([0#6H1r6R1D2X3x2+0a]3cc[ch][x2Hc][$([#6aH]%73[$([ac]5[c^n][cH1a][x2Hc][acr6][0#6H1r6R1D2X3x2+0a]5)][cH1a][0#6+0a][ac][cX3H]%73)]3)][#6aH][$([0#6H1r6R1D2X3x2+0a]%29c[cX3H][ch][$(cc)][$([0#6+0a]%15[#6aH][cX3H][0#6H1r6R1D2X3x2+0a][x2Hc][cH1a]%15)]%29)]%98)]%46[cH1a][$([0#6H1r6R1D2X3x2+0a]%56[$([0#6+0a]%51[cR1][$([cX3H]%82[ch][0#6H1r6R1D2X3x2+0a][cR1]c[c^n]%82)][cH1a][#6aH][$([$(c%56[$([cX3H]%61[$(cc)][$(c%75[0#6+0a][x2Hc][ch][$(cc)][ac]%75)]c[cH1a][cR1]%61)]c[$([cH1a]%87[0#6+0a][c^n][x2Hc][cX3H]c%87)][c^n][cX3H]%56)]%11c[cR1][#6aH][0#6+0a][ac]%11)]%51)][$([ch]%59[x2Hc][ac][0#6H1r6R1D2X3x2+0a][acr6][c^n]%59)][cX3H][ac][$([acr6]%91[c^n]c[ch][cX3H][x2Hc]%91)]%56)][0#6H1r6R1D2X3x2+0a][$([0#6+0a]%19[c^n][cX3H][$([ac]%46[acr6][0#6H1r6R1D2X3x2+0a][x2Hc][c^n]c%46)][acr6][x2Hc]%19)][cX3H]%46)]%60)][acr6][ch]%94)]%54)][acr6][ac][x2Hc]%32)][$([cX3H]%87[cH1a][0#6+0a][x2Hc][$(c%31[$([$([#6aH]%93[ac][0#6H1r6R1D2X3x2+0a][0#6+0a][$([$([ch]%73[#6aH][$([$(c%97[c^n][x2Hc][0#6H1r6R1D2X3x2+0a][0#6+0a][acr6]%97)]%72[0#6+0a][$([0#6+0a]%12[0#6H1r6R1D2X3x2+0a][$(cc)][ac][cH1a]c%12)][ch][0#6H1r6R1D2X3x2+0a][acr6]%72)][x2Hc][$([ac]%23[0#6+0a][$([cR1]%61[cX3H]c[ch][#6aH][cH1a]%61)][cH1a][0#6H1r6R1D2X3x2+0a][cX3H]%23)][$(cc)]%73)]%92[$([$([#6aH]%85[cR1][0#6H1r6R1D2X3x2+0a][cH1a][0#6+0a][cX3H]%85)]%16[x2Hc]c[acr6]c[#6aH]%16)][ch][cH1a][$([cH1a]%31[$(cc)][#6aH][0#6+0a][$([0#6H1r6R1D2X3x2+0a]%21[$(cc)]c[x2Hc][c^n][acr6]%21)][x2Hc]%31)][$(c%19[c^n][$([cX3H]7[x2Hc][#6aH][acr6][ch][cR1]7)][cH1a][ch][cR1]%19)]%92)][$([0#6H1r6R1D2X3x2+0a]%34[acr6][x2Hc]c[0#6+0a]c%34)]%93)]4[acr6][$([cX3H]%89[0#6+0a][c^n][$(cc)][#6aH][$([#6aH]%26[0#6+0a][cH1a][acr6][0#6H1r6R1D2X3x2+0a][x2Hc]%26)]%89)][ch][$(cc)][$([x2Hc]%13[ch][acr6][$([$(cc)]%26[0#6H1r6R1D2X3x2+0a][0#6+0a]c[x2Hc][cH1a]%26)][0#6H1r6R1D2X3x2+0a][ac]%13)]4)][$([x2Hc]%97[c^n]c[cX3H][$(cc)][cR1]%97)][#6aH][$([$([0#6H1r6R1D2X3x2+0a]%66c[ac][x2Hc][0#6+0a][cR1]%66)]%38c[x2Hc]c[$([$([cX3H]%98[c^n][$([ch]%78[#6aH][0#6+0a]c[0#6H1r6R1D2X3x2+0a][$(cc)]%78)][$([#6aH]%30[x2Hc]c[cH1a][acr6][ac]%30)][cH1a][#6aH]%98)]%46[$(cc)][0#6+0a][ch][$([cR1]%19[acr6][#6aH][0#6+0a]c[ac]%19)][0#6H1r6R1D2X3x2+0a]%46)][$([$(c%18[$([ac]%62c[$(cc)][$(c%18[cH1a][0#6+0a][0#6H1r6R1D2X3x2+0a][ch][cX3H]%18)][x2Hc]c%62)][acr6][ch][cH1a][0#6+0a]%18)]%11c[$(c%24[x2Hc][#6aH][c^n][acr6][cX3H]%24)][0#6H1r6R1D2X3x2+0a][ac][cH1a]%11)]%38)][$([0#6H1r6R1D2X3x2+0a]%74[cH1a][acr6][x2Hc][$([$([cR1]%74[#6aH][x2Hc][cX3H][$(c%52[cR1]c[$(cc)][cH1a][0#6+0a]%52)][$(cc)]%74)]%48c[cH1a][$([0#6+0a]%27[cH1a][c^n][acr6]c[#6aH]%27)][c^n][$(cc)]%48)][c^n]%74)]%31)][$(cc)]%87)]c8)]"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  set_atom_environment_only_matches_unmatched_atoms(0);

  _smiles = "c1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 6);
}

TEST_F(TestSubstructure, TestNoSaveMatchedAtoms)
{
  _string_proto = R"(query {
    save_matched_atoms: false
    smarts: "C!@C"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
  EXPECT_EQ(_sresults.number_embeddings(), 0);
}

TEST_F(TestSubstructure, TestOneEmbeddingPerStartAtom)
{
  _string_proto = R"(query {
    one_embedding_per_start_atom: true
    smarts: "[CD4]C"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CC(C)(C)C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  EXPECT_EQ(_sresults.number_embeddings(), 1);
}

TEST_F(TestSubstructure, TestMaxMatchesToFind)
{
  _string_proto = R"(query {
    max_matches_to_find: 3
    smarts: "[CD4]C"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CC(C)(C)C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);
  EXPECT_EQ(_sresults.number_embeddings(), 3);
  EXPECT_EQ(_sresults.embedding(0)->item(0), 1);
  EXPECT_EQ(_sresults.embedding(1)->item(0), 1);
  EXPECT_EQ(_sresults.embedding(2)->item(0), 1);
}

TEST_F(TestSubstructure, TestSubtractFromRc)
{
  _string_proto = R"(query {
    subtract_from_rc: 3
    smarts: "[$(cF)]"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "Fc1c(F)c(F)c(F)c(F)c1F";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);
  EXPECT_EQ(_sresults.number_embeddings(), 6);
  EXPECT_EQ(_sresults.embedding(0)->item(0), 1);
  EXPECT_EQ(_sresults.embedding(1)->item(0), 2);
  EXPECT_EQ(_sresults.embedding(2)->item(0), 4);
}

TEST_F(TestSubstructure, TestImplicitRingConditionChainOK)
{
  _string_proto = R"(query {
    implicit_ring_condition: 0
    smarts: "CCC"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
  EXPECT_EQ(_sresults.number_embeddings(), 2);
  EXPECT_THAT(*_sresults.embedding(0), ElementsAre(0, 1, 2));
}

TEST_F(TestSubstructure, TestImplicitRingConditionChainNotMet)
{
  _string_proto = R"(query {
    implicit_ring_condition: 1
    smarts: "CCC"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

#ifdef UNSURE_ABOUT_THIS_TEST_NEEDS_INVESTIGATION
TEST_F(TestSubstructure, TestImplicitRingConditionRingMet)
{
  _string_proto = R"(query {
    implicit_ring_condition: 1
    smarts: "CCC"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CC1CC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
  _query.set_do_not_perceive_symmetry_equivalent_matches(1);
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}
#endif

TEST_F(TestSubstructure, TestTwoPiElectronsAromatic) {
  _string_proto = R"(query {
    smarts: "[ND1H2]-c(:c):c=O"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  set_allow_two_electron_systems_to_be_aromatic(1);
//_smiles = "CN(C)C1=C(N)C(=O)C1=O";
  _smiles = "CN(C)c1c(N)c(=O)c1=O";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestRidSame1) {
  _string_proto = R"(query {
    smarts: "[/IWrid1].[/IWrid1]"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1CC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 6);
}

TEST_F(TestSubstructure, TestRidSame2) {
  _string_proto = R"(query {
    smarts: "[/IWrid1].[/IWrid1]"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1CC1CC1CC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 12);
}

TEST_F(TestSubstructure, TestRidDifferentNumberOneRing) {
  _string_proto = R"(query {
    smarts: "[/IWrid1].[/IWrid2]"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1CC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestRidDifferentNumberTwoRing) {
  _string_proto = R"(query {
    smarts: "[/IWrid1R].[/IWrid2R]"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1CC1CC1CC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 18);
}

// This leaves out the R specification in the smarts. In that case
// the non-ring Carbon atom also counts as not in the same ring.
TEST_F(TestSubstructure, TestRidChainAtomIsNotInTheSameRing) {
  _string_proto = R"(query {
    smarts: "[/IWrid1].[/IWrid2]"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1CC1CC1CC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 30);
}

TEST_F(TestSubstructure, TestGlobalIDRing) {
  _string_proto = R"(query {
    ring_specifier {
      base {
        heteroatom_count: 1
        set_global_id: 2
      }
      aromatic: true
    }
    smarts: "[F,Br]-[/IWgid2a]"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "Fc1ccc(cc1)CCc1cnc(Br)cc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->size(), 2);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 35);
}

TEST_F(TestSubstructure, TestGlobalIDRingSysMatches) {
  _string_proto = R"(query {
    ring_system_specifier {
      base {
        set_global_id: 3
      }
      ring_size: [4, 6]
      aromatic_ring_count: 1
    }
    smarts: "[F,Br]-[/IWgid3a]"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "Fc1ccc2c(c1)C(C2)c1ccc2ccc(Cl)cc2c1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->size(), 2);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 9);
}

TEST_F(TestSubstructure, TestGlobalIDRingSysNoMatch) {
  _string_proto = R"(query {
    ring_system_specifier {
      base {
        set_global_id: 8
      }
      ring_size: [4, 6]
      aromatic_ring_count: 1
    }
    smarts: "[F,Br]-[/IWgid3a]"
  }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "Fc1ccc2c(c1)C(C2)c1ccc2ccc(Cl)cc2c1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

struct ProtoMolMatches {
  std::string proto;
  IWString smiles;
  int expected;
};

std::ostream&
operator<<(std::ostream& output, const ProtoMolMatches& pmm) {
  output << pmm.proto << ' ' << " smiles " << pmm.smiles << " expected " << pmm.expected;

  return output;
}

TEST_F(TestSubstructure, TestAndWithTwoZeros) {
  ASSERT_TRUE(_query.create_from_smarts("0[<23C]&&0[>23C]"));

  _smiles = "[22C]";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
  _smiles = "[24C]";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
  _smiles = "[23C]";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

const std::string aminunch = R"pb(
query {
  query_atom {
    id: 0
    atom_properties {
      atomic_number: 7
      formal_charge: 0
      unsaturation: 0
      aromatic: false
      max_aryl: 1
    }
  }
  environment_no_match {
    attachment {
      attachment_point: 0
      btype: SS_SINGLE_BOND
    }
    query_atom {
      id: 2
      atom_properties {
        atomic_number: [6, 7, 16]
      }
    }
    query_atom {
      id: 3
      atom_properties {
        atomic_number: [8, 16]
      }
      bond_smarts: "= 2"
    }
  }
  environment_no_match {
    attachment {
      attachment_point: 0
      btype: SS_SINGLE_BOND
    }
    query_atom {
      id: 4
      atom_properties {
        atomic_number: [6, 7, 8]
      }
    }
    query_atom {
      id: 5
      atom_properties {
        atomic_number: 7
        formal_charge: 1
      }
      query_bond {
        other_end: 4
        btype: [SS_DOUBLE_BOND, SS_AROMATIC_BOND]
      }
    }
  }
}
)pb";

const std::string ncon1 = R"pb(
query {
  smarts: "c1ccccc1"
  ncon: 1
}
)pb";

const std::string substituent_min_natoms2 = R"pb(
query {
  ring_specifier {
    base {
      substituent {
        min_natoms: 2
        set_global_id: 3
      }
    }
  }
  smarts: "[/IWgid3]"
}
)pb";
const std::string substituent_min_natoms4 = R"pb(
query {
  ring_specifier {
    base {
      substituent {
        min_natoms: 4
        set_global_id: 3
      }
    }
  }
  smarts: "[/IWgid3]"
}
)pb";
const std::string substituent_min_natoms7 = R"pb(
query {
  ring_specifier {
    base {
      substituent {
        min_natoms: 7
        set_global_id: 3
      }
    }
  }
  smarts: "[/IWgid3]"
}
)pb";
const std::string substituent_nrings0 = R"pb(
query {
  ring_specifier {
    base {
      substituent {
        nrings: 0
        set_global_id: 3
      }
    }
  }
  smarts: "[/IWgid3]"
}
)pb";
const std::string substituent_min_nrings1 = R"pb(
query {
  ring_specifier {
    base {
      substituent {
        min_nrings: 1
        set_global_id: 3
      }
    }
  }
  smarts: "[/IWgid3]"
}
)pb";
const std::string substituent_min_nrings1b = R"pb(
query {
  ring_specifier {
    base {
      substituent {
        min_nrings: 1
        set_global_id: 3
      }
      environment: "O"
    }
  }
  smarts: "[/IWgid3]"
}
)pb";
const std::string substituent_length2 = R"pb(
query {
  ring_specifier {
    base {
      substituent {
        length: 2
        set_global_id: 3
      }
    }
  }
  smarts: "[/IWgid3]"
}
)pb";
const std::string substituent_isotope1 = R"pb(
query {
  ring_specifier {
    base {
      substituent {
        set_global_id: 3
        required_smarts: "[1]"
      }
    }
  }
  smarts: "[/IWgid3]"
}
)pb";
const std::string substituent_different_global_ids_23 = R"pb(
query {
  ring_specifier {
    base {
      set_global_id: 2
      substituent {
        set_global_id: 3
        natoms: 2
      }
    }
  }
  smarts: "[/IWgid3]"
}
)pb";
const std::string substituent_different_global_ids_33 = R"pb(
query {
  ring_specifier {
    base {
      set_global_id: 2
      substituent {
        set_global_id: 3
        natoms: 2
      }
    }
  }
  smarts: "[/IWgid2]"
}
)pb";

class TestSubstructureP: public testing::TestWithParam<ProtoMolMatches> {
  protected:
    SubstructureSearch::SubstructureQuery _proto;
    Substructure_Query _query;
    Molecule _mol;
    Substructure_Results _sresults;
};

TEST_P(TestSubstructureP, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  // cerr << "Testing " << params.smiles << " expecting " << params.expected << '\n';
  EXPECT_EQ(_query.substructure_search(_mol, _sresults), params.expected) << params;
}
INSTANTIATE_TEST_SUITE_P(TestSubstructureP, TestSubstructureP, testing::Values(
  ProtoMolMatches{aminunch, "CN", 1},
  ProtoMolMatches{aminunch, "C[N+H3]", 0},
//ProtoMolMatches{aminunch, "ON", 0},
  ProtoMolMatches{aminunch, "C=N", 0},
  ProtoMolMatches{aminunch, "C#N", 0},
  ProtoMolMatches{aminunch, "Nc1ccccc1", 1},
  ProtoMolMatches{aminunch, "c1ccccc1Nc1ccccc1", 0},
  ProtoMolMatches{aminunch, "c1cnccc1", 0},
  // First environment_no_match
  ProtoMolMatches{aminunch, "CNC=O", 0},
  ProtoMolMatches{aminunch, "CNC=S", 0},
  ProtoMolMatches{aminunch, "CNN=O", 0},
  ProtoMolMatches{aminunch, "CNN=S", 0},
  ProtoMolMatches{aminunch, "CNS=O", 0},

  ProtoMolMatches{ncon1, "c1ccccc1", 0},
  ProtoMolMatches{ncon1, "c1ccccc1C", 12},
  ProtoMolMatches{ncon1, "Cc1ccccc1C", 0},

  ProtoMolMatches{substituent_min_natoms4, "CCCCCCC1CC1", 6},
  ProtoMolMatches{substituent_min_natoms7, "CCCCCCC1CC1", 0},
  ProtoMolMatches{substituent_min_natoms2, "CCC1CC1CC", 4},
  ProtoMolMatches{substituent_nrings0, "CCCCCCC1CC1", 6},
  ProtoMolMatches{substituent_min_nrings1,  "CCCCCCC1CC1", 0},
  ProtoMolMatches{substituent_min_nrings1b,  "C1OC1CCCC1CC1", 6},
  ProtoMolMatches{substituent_length2,  "CCCCCCC1CC1C", 0},
  ProtoMolMatches{substituent_length2,  "CCC1CC1C", 2},
  ProtoMolMatches{substituent_length2,  "CC(C)(C)C1CC1C", 4},
  ProtoMolMatches{substituent_length2,  "CCC(C)(C)C1CC1C", 0},
  ProtoMolMatches{substituent_isotope1,  "CCC(C)(C)C1CC1C", 0},
  ProtoMolMatches{substituent_isotope1,  "C[1CH2](C)(C)C1CC1CN", 4},
  ProtoMolMatches{substituent_isotope1,  "C[1CH2](C)(C)C1CC1C[1NH2]", 6},
  ProtoMolMatches{substituent_different_global_ids_23,  "CCC1CC1C", 2},
  ProtoMolMatches{substituent_different_global_ids_33,  "CCC1CC1C", 3},
  
  ProtoMolMatches{R"pb(
query {
  smarts: "NC=O"
  substituent {
    heteroatom_count: 0
  }
}
)pb", "NC(=O)C", 1},
  
  ProtoMolMatches{R"pb(
query {
  smarts: "NC=O"
  substituent {
    heteroatom_count: 0
  }
}
)pb", "NC(=O)CCCO", 0},
  
  ProtoMolMatches{R"pb(
query {
  smarts: "NC=O"
  substituent {
    heteroatom_count: 1
  }
}
)pb", "NC(=O)CCCO", 1},
 
  ProtoMolMatches{R"pb(
query {
  smarts: "NC=O"
  substituent {
    heteroatom_count: 1
    length: 4
    natoms: 4
    nrings: 0
    unsaturation_count: 0
    hits_needed: 1
  }
}
)pb", "NC(=O)CCCO", 1},
 
  ProtoMolMatches{R"pb(
query {
  smarts: "CCCC"
  unique_embeddings_only: true
  substituent {
    heteroatom_count: 1
    length: 1
    natoms: 1
    nrings: 0
    unsaturation_count: 0
    hits_needed: 4
  }
}
)pb", "FC(F)(F)CCCF", 1}
));

struct SmilesSmartsMatches {
  IWString smiles;
  IWString smarts;
  int expected;
};
class TestSubstructureSmartsP: public testing::TestWithParam<SmilesSmartsMatches>  {
  protected:
    Substructure_Query _query;
    Molecule _mol;
    Substructure_Results _sresults;
};

TEST_P(TestSubstructureSmartsP, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  //std::cerr << "Testing " << params.smiles << " smarts " << params.smarts << " expecting " << params.expected << '\n';
  EXPECT_EQ(_query.substructure_search(_mol, _sresults), params.expected);
}
INSTANTIATE_TEST_SUITE_P(TestSubstructureSmartsP, TestSubstructureSmartsP, testing::Values(
  SmilesSmartsMatches{"C", "1C", 1},
  SmilesSmartsMatches{"C", "0C", 0},
  SmilesSmartsMatches{"C", "2C", 0},
  SmilesSmartsMatches{"CC", "1C", 0},
  SmilesSmartsMatches{"CC", "2C", 2},
  SmilesSmartsMatches{"CCC", "2C", 0},
  SmilesSmartsMatches{"CCC", ">1C", 3},
  SmilesSmartsMatches{"CCC", ">2C", 3},
  SmilesSmartsMatches{"CCC", "3C", 3},
  SmilesSmartsMatches{"CCC", ">3C", 0},
  SmilesSmartsMatches{"CCC", "<3C", 0},
  SmilesSmartsMatches{"CCC", "<4C", 3},
  SmilesSmartsMatches{"CCC", "1-2C", 0},
  SmilesSmartsMatches{"CCC", "1-3C", 3},
  SmilesSmartsMatches{"CCC", "1-4C", 3},
  SmilesSmartsMatches{"CCC", "2-4C", 3},
  SmilesSmartsMatches{"CCC", "3-4C", 3},
  SmilesSmartsMatches{"CCC", "3-3C", 3},
  SmilesSmartsMatches{"CCC", "4-28C", 0},

  // Aromatic sulphur tests
  SmilesSmartsMatches{"S1=CC=CC=C1", "s", 1},
  SmilesSmartsMatches{"[SH]1=CC=CC=C1", "s", 1},
  SmilesSmartsMatches{"CS1=CC=CC=C1", "s", 1}
));

const std::string inter_ring_region1 = R"pb(
query {
  inter_ring_region {
  }
}
)pb";

class TestSubstructureInterRingRegion: public testing::TestWithParam<ProtoMolMatches>  {
  protected:
    SubstructureSearch::SubstructureQuery _proto;
    Substructure_Query _query;
    Molecule _mol;
    Substructure_Results _sresults;
};

TEST_P(TestSubstructureInterRingRegion, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  //std::cerr << "Testing " << params.smiles << " smarts " << params.smarts << " expecting " << params.expected << '\n';
  EXPECT_EQ(_query.substructure_search(_mol, _sresults), params.expected);
}
INSTANTIATE_TEST_SUITE_P(TestSubstructureInterRingRegion, TestSubstructureInterRingRegion, testing::Values(
  ProtoMolMatches{inter_ring_region1, "C1CC1C1CC1", 0}
));

const std::string fraction_matched = R"pb(
query {
  min_all_matches_fraction_atoms_matched: 0
  max_all_matches_fraction_atoms_matched: 0
}
)pb";


struct SmiFractions {
  IWString smiles;
  IWString smarts;
  float min_fraction = {};
  float max_fraction = {};
  int expected = {};
};

class TestAllEmbeddingsFractionMatched: public testing::TestWithParam<SmiFractions>  {
  protected:
    SubstructureSearch::SubstructureQuery _proto;
    Substructure_Query _query;
    Molecule _mol;
    Substructure_Results _sresults;
};

TEST_P(TestAllEmbeddingsFractionMatched, Tests) {
  const auto params = GetParam();

  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));

  SubstructureSearch::SingleSubstructureQuery* q = _proto.add_query();
  q->set_smarts(params.smarts.data(), params.smarts.length());
  q->set_min_all_matches_fraction_atoms_matched(params.min_fraction);
  q->set_max_all_matches_fraction_atoms_matched(params.max_fraction);

  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  //std::cerr << "Testing " << params.smiles << " smarts " << params.smarts << " min " << params.min_fraction << " max " << params.max_fraction << " expecting " << params.expected << '\n';
  EXPECT_EQ(_query.substructure_search(_mol, _sresults), params.expected);
}
INSTANTIATE_TEST_SUITE_P(TestAllEmbeddingsFractionMatched, TestAllEmbeddingsFractionMatched, testing::Values(
  SmiFractions{"CC", "C", 0, 0, 2},
  SmiFractions{"CCC", "C", 0, 0, 3},
  SmiFractions{"CCC", "C", 1, 0, 3},
  SmiFractions{"CCC", "C", 0.5, 0, 3},
  SmiFractions{"CCC", "C", 0, 1, 3},
  SmiFractions{"CCC", "C", 0.1, 0.0, 3},
  SmiFractions{"CCC", "C", 0.2, 0.0, 3},
  // unspecified max gets set to 1.
  SmiFractions{"CCC", "C", 0.9, 0.0, 3},
  SmiFractions{"CNN", "C", 0.0, 0.0, 1},
  SmiFractions{"CNN", "C", 0.33, 0.0, 1},
  SmiFractions{"CNN", "C", 0.34, 0.0, 0},
  SmiFractions{"CNN", "C", 0.0, 0.33, 0},
  SmiFractions{"CNN", "C", 0.0, 0.34, 1},
  SmiFractions{"CNN", "C", 0.33, 0.34, 1}
));


class TestComponentGrouping: public testing::TestWithParam<SmilesSmartsMatches>  {
  protected:
    Substructure_Query _query;
    Molecule _mol;
    Substructure_Results _sresults;
};

TEST_P(TestComponentGrouping, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  //std::cerr << "Testing " << params.smiles << " smarts " << params.smarts << " expecting " << params.expected << '\n';
  EXPECT_EQ(_query.substructure_search(_mol, _sresults), params.expected) <<
     "mol " << params.smiles << " smt " << params.smarts << " expected " << params.expected;
}
INSTANTIATE_TEST_SUITE_P(TestComponentGrouping, TestComponentGrouping, testing::Values(
  SmilesSmartsMatches{"CCCC", "C.C", 12},
  SmilesSmartsMatches{"CCCC", "(C.C)", 12},
  SmilesSmartsMatches{"CCCC", "(C).(C)", 0},
  SmilesSmartsMatches{"CCCC.CCCC", "(C).(C)", 32},
  SmilesSmartsMatches{"CCCC", "(C).C", 12},
  SmilesSmartsMatches{"CCCC.CCCC", "(C).(C).C", 192}
));

}  // namespace
