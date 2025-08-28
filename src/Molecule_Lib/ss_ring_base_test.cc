// Tests for functions that happen to reside in ss_ring_base.
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "molecule.h"
#include "substructure.h"
#include "target.h"

namespace {

using testing::UnorderedElementsAre;
using testing::UnorderedElementsAreArray;

class TestSubstituent : public testing::Test
{
  protected:

  protected:
    std::string _string_proto;
    SubstructureSearch::SubstructureQuery _proto;
    IWString _smiles;
    Molecule _mol;
    Substructure_Query _query;
    Substructure_Results _sresults;
};

TEST_F(TestSubstituent, TestAtomCountMatches) {
  _string_proto = R"pb(
query {
  ring_specifier {
    aromatic: false
    base {
      substituent {
        natoms: 1;
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_EQ(_sresults.embedding(0)->item(0), 3);
}

TEST_F(TestSubstituent, TestAtomCountNoMatch) {
  _string_proto = R"pb(
query {
  ring_specifier {
    aromatic: false
    base {
      substituent {
        natoms: 2;
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0);
}

TEST_F(TestSubstituent, TestAtomCountMultipleMatches) {
  _string_proto = R"pb(
query {
  ring_specifier {
    aromatic: false
    base {
      substituent {
        natoms: 1;
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "CC1CC1C";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 2);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0));
  EXPECT_THAT(*_sresults.embedding(1), UnorderedElementsAre(4));
}

TEST_F(TestSubstituent, TestSSMatches) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_specifier {
    aromatic: false
    base {
      substituent {
        required_smarts: "[CD1]";
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "CCC1CC1C";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({0, 1, 5}));
}

TEST_F(TestSubstituent, TestSSMatchesNumericQualifier) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_specifier {
    aromatic: false
    base {
      substituent {
        required_smarts: "2[CD1]";
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "CC(C)C1CC1C";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({0, 1, 2}));
}

TEST_F(TestSubstituent, TestSSDisqualifying) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_specifier {
    aromatic: false
    base {
      substituent {
        disqualifying_smarts: "O";
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "CC(C)C1CC1C";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({0, 1, 2, 6}));
}

TEST_F(TestSubstituent, TestNrings0) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 0
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c12ccccc1cc(CC)cc2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({8, 9}));
}

TEST_F(TestSubstituent, TestNrings0NoMatch) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 0
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c12ccccc1cc(CCC1CC1)cc2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0);
}


TEST_F(TestSubstituent, TestNrings1Match) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 1
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c12ccccc1cc(CCC1CC1)cc2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({8, 9, 10, 11, 12}));
}

TEST_F(TestSubstituent, TestNrings1NoMatch) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 2
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c12ccccc1cc(CCC1CC1)cc2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0);
}

TEST_F(TestSubstituent, TestNrings1NoMatchTooManyRings) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 1
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c12ccccc1cc(CCC1CC1C1CC1)cc2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0);
}

TEST_F(TestSubstituent, TestNrings1NoSubstituents) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 1
      }
    }
  }
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "S1N=NC2=C1SC=C2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0);
}

TEST_F(TestSubstituent, TestNringsNoMatchNoRing) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 1
      }
    }
  }
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "N1=C2C=NNC2=C(C)N1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0);
}

TEST_F(TestSubstituent, TestNringsSeparateRings) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 1
    rings_in_system: 1
    base {
      substituent {
        nrings: 2
        required_smarts: "[/IWfsid1].[/IWfsid2]"
        set_global_id: 2
      }
    }
  }
  smarts: "[/IWgid2]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c1cncc1C(C1CC1)(C1CC1)";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({5, 6, 7, 8, 9, 10, 11}));
}
TEST_F(TestSubstituent, TestNringsFusedRingsShouldNotMatch) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 1
    rings_in_system: 1
    base {
      substituent {
        nrings: 2
        required_smarts: "[/IWfsid1R].[/IWfsid2R]"
        set_global_id: 2
      }
    }
  }
  smarts: "[/IWgid2]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c1cncc1CC12CC1C2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0);
}

TEST_F(TestSubstituent, InterRingRegion1) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    natoms: 1
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({3}));
}

TEST_F(TestSubstituent, InterRingRegion1NoMatch) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    natoms: 1
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1CCC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0);
}

TEST_F(TestSubstituent, InterRingRegionOkLength2) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    length: 2
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({3}));
}

TEST_F(TestSubstituent, InterRingRegionOkLength3) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    length: 3
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({3, 7}));
}

TEST_F(TestSubstituent, InterRingRegionOkRequiredLength3) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    required_length: [2, 3, 3]
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({3, 7}));
}

TEST_F(TestSubstituent, InterRingRegionRequiredLengthNotMatched) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    required_length: [2, 3, 4]
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0);
}

TEST_F(TestSubstituent, InterRingRegionRingConnections) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    ring_connections: 2
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)CC1CC1CCCC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({11, 12, 13}));
}

TEST_F(TestSubstituent, InterRingRegionRingHitsNeededOk) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    ring_connections: 2
    hits_needed: 2
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)C1CC1CC1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({10, 14}));
}

TEST_F(TestSubstituent, InterRingRegionRingHitsNeededTooSmall) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    ring_connections: 2
    hits_needed: 1
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)C1CC1CC1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0);
}

TEST_F(TestSubstituent, InterRingRegionRingHitsNeededTooLarge) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    ring_connections: 2
    hits_needed: 3
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)C1CC1CC1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0);
}

TEST_F(TestSubstituent, InterRingRegionRingRequiredSmarts) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    ring_connections: 2
    hits_needed: 2
    required_smarts: "[R]!@[AD2]!@[R]"
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)C1CC1CC1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({10, 14}));
}

TEST_F(TestSubstituent, InterRingRegionRingDisqualifyingSmarts) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    disqualifying_smarts: "[R]!@[AD2]!@[R]"
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)C1CC1CC1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({3}));
}

struct ProtoSmilesExpected {
  std::string proto;
  IWString smiles;
  int expected;
};

std::ostream&
operator<<(std::ostream& output, const ProtoSmilesExpected& pse) {
  output << pse.proto << ' ' << pse.smiles << " expect " << pse.expected;

  return output;
}

class TestSubstructureP: public testing::TestWithParam<ProtoSmilesExpected> {
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
  ProtoSmilesExpected{R"pb(
query {
  ring_specifier {
    base {
      substituent {
        natoms: 1;
        heteroatom_count: 0
        set_global_id: 1
        hits_needed: 2
      }
    }
  }
  smarts: "[/IWgid1C]"
}
)pb", "C1CC1(C)C", 2},

  ProtoSmilesExpected{R"pb(
query {
  substituent {
    max_natoms: 3
    heteroatom_count: 0
    no_other_substituents_allowed: true
  }
  smarts: "[CD3T2](=O)-[ND2]"
}
)pb", "NCCC(=O)NCCC", 0}

));


}  // namespace
