#include <iostream>
#include <string>
#include <vector>

// Tests for Substructre_Atom_Specifier

#include "aromatic.h"
#include "substructure.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {

using std::cerr;
using std::endl;

using testing::UnorderedElementsAre;

//using google::protobuf::TextFormat::ParseFromString;

class TestSubstructureSpec : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    std::string _string_proto;

    IWString _smiles;

    Substructure_Query _query;

    Substructure_Results _sresults;

    Molecule _m;

    SubstructureSearch::SubstructureQuery _proto;

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
};

void
TestSubstructureSpec::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

TEST_F(TestSubstructureSpec, AtomicNumber)
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

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  ASSERT_TRUE(_m.build_from_smiles("C"));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructureSpec, AtomicSymbol)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_symbol: "C"
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CN";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms * e = _sresults.embedding(0);

  EXPECT_EQ(e->number_elements(), 1);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
}

TEST_F(TestSubstructureSpec, Ncon1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          ncon: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C.NCN.COC.FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  for (int i = 0; i < 2; ++i) {
    const Set_of_Atoms * e = _sresults.embedding(i);
    EXPECT_EQ(e->number_elements(), 1);
    EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
    EXPECT_EQ(_m.ncon(e->item(0)), 1);
    EXPECT_EQ(_m.attached_heteroatom_count(e->item(0)), 1);
  }

  const std::vector<int> matched_atoms = {_sresults.embedding(0)->item(0),
                                          _sresults.embedding(1)->item(0)};
  EXPECT_THAT(matched_atoms, UnorderedElementsAre(4, 6));
}

TEST_F(TestSubstructureSpec, MinNcon1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          min_ncon: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C.NCN.COC.NC(N)N.C(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  const std::vector<int> matched_atoms = {_sresults.embedding(0)->item(0),
                                          _sresults.embedding(1)->item(0)};
  EXPECT_THAT(matched_atoms, UnorderedElementsAre(8, 11));
}

TEST_F(TestSubstructureSpec, Ncon2)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          ncon2: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC(C(F)(F)F)C(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms* e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 0);
}

TEST_F(TestSubstructureSpec, Nbonds)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          nbonds: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC=N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms* e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 1);
}

TEST_F(TestSubstructureSpec, FormalChargePos)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          formal_charge: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "[N+]CC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms* e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 0);
}

TEST_F(TestSubstructureSpec, FormalChargeNeg)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 8
          formal_charge: -1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "[O-]C(=O)c1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms* e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 0);
}

TEST_F(TestSubstructureSpec, FormalChargeZero)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          formal_charge: 0
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "[N+]CN";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms* e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 2);
}

TEST_F(TestSubstructureSpec, NringsPositive)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 8
          nrings: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC1OC1.O1C2CCC1CC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms* e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 2);
}

TEST_F(TestSubstructureSpec, RingBondCount)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 8
          ring_bond_count: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC1OC1.O1C2CCC1CC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);
  Set_of_Atoms matched;
  for (const auto* e : _sresults.embeddings()) {
    matched.add(e->item(0));
  }

  EXPECT_THAT(matched, UnorderedElementsAre(2, 4));
}

TEST_F(TestSubstructureSpec, RingSize)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          ring_size: 5
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CC1.N1CCC1.[1N]1CCCC1.N1CCCCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 7);
  EXPECT_EQ(_m.isotope(e->item(0)), 1);
}

TEST_F(TestSubstructureSpec, Hcount)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          hcount: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 1);
}

TEST_F(TestSubstructureSpec, Aromatic)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          aromatic: true
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ncccc1C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 5);
}

TEST_F(TestSubstructureSpec, Chirality)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          chirality: true
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C[C@H](N)CC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->item(0), 1);
}

TEST_F(TestSubstructureSpec, AromaticRingSize)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          aromatic_ring_size: 6
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "n1ccc2cnccc12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 5);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();

  EXPECT_THAT(e, UnorderedElementsAre(3, 4, 6, 7, 8));
}

TEST_F(TestSubstructureSpec, AliphaticRingSize)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          aliphatic_ring_size: 4
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1OC1.C1CCC1.C1CCCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 4);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();

  EXPECT_THAT(e, UnorderedElementsAre(3, 4, 5, 6));
}

TEST_F(TestSubstructureSpec, AttachedHeteroatomCount)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          attached_heteroatom_count: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C.CC.C(C)C.C(C)(C)C.N[1CH](N)N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();

  EXPECT_EQ(_m.isotope(e[0]), 1);
}

TEST_F(TestSubstructureSpec, LonePairCount1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          lone_pair_count: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NC.CNC.CN(C)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 3, 6));
}

TEST_F(TestSubstructureSpec, LonePairCount2)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 8
          lone_pair_count: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "O.OC.COC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 1, 4));
}

TEST_F(TestSubstructureSpec, Unsaturation1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 8
          unsaturation: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "O.OC.COC.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_EQ(e[0], 8);
}

TEST_F(TestSubstructureSpec, Unsaturation2)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          unsaturation: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CCN.C=C.C#N.C=C=N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(5, 8));
}

TEST_F(TestSubstructureSpec, DaylightX1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          daylight_x: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "[CH].CN.[C]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_EQ(e[0], 0);
}

TEST_F(TestSubstructureSpec, DaylightX2)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          daylight_x: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C[C]C.[CH].CN.[C]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_EQ(e[0], 1);
}

TEST_F(TestSubstructureSpec, DaylightX3)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          daylight_x: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N.CN.[N+]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 2));
}

TEST_F(TestSubstructureSpec, Isotope)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 92
          isotope: [3, 4]
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "[3U].[4U].[U]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 1));
}

TEST_F(TestSubstructureSpec, Aryl)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          aromatic: false
          aryl: [1, 2]
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCNc1ccccc1.c1cccnc1Nc1ncccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(2, 15));
}


TEST_F(TestSubstructureSpec, Vinyl)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          vinyl: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CNC.COC.CC(=O)C.Cc1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(6, 9));
}

TEST_F(TestSubstructureSpec, FusedSystemSize)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          fused_system_size: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1.N1C2CCC1CC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_EQ(e[0], 4);
}

TEST_F(TestSubstructureSpec, HeteroatomsInRing)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          heteroatoms_in_ring: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1.N1C2CCC1CC2.C1NCNC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(12, 14));
}

TEST_F(TestSubstructureSpec, SpinachOnlyNoRings1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          match_spinach_only: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N.CN";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructureSpec, SpinachOnlyNoRings0)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          match_spinach_only: 0
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N.CN.C1CC1NC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);
}

TEST_F(TestSubstructureSpec, SpinachOnly)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          match_spinach_only: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NC1CCN1.N1CCC1.N1C2CCC1CC2.C1NCNC1.C1CC1NC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_EQ(e[0], 0);
}

TEST_F(TestSubstructureSpec, ScaffoldBondsAttachedToRing1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          scaffold_bonds_attached_to_ring: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1CC1CC1.C1CC1CCC1CC1.C1CC1C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 12);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 14));
}

TEST_F(TestSubstructureSpec, ScaffoldBondsAttachedToRing2)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          scaffold_bonds_attached_to_ring: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1CC1CC1CCC1CC1F.C1CC1CC1CC1.C1CC1CCC1CC1.C1CC1C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(4, 5, 6));
}

TEST_F(TestSubstructureSpec, SymmetryDegree1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 9
          symmetry_degree: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FCC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_EQ(e[0], 0);
}

TEST_F(TestSubstructureSpec, SymmetryDegree3)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 9
          symmetry_degree: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FCC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(3,4,5));
}

TEST_F(TestSubstructureSpec, SymmetryDegree6)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 9
          symmetry_degree: 6
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)c1ccc(C(F)(F)F)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 6);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 2, 3, 9, 10, 11));
}

TEST_F(TestSubstructureSpec, SymmetryDegree3b)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 9
          symmetry_degree: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)c1ccc(C(F)(F)F)cc1F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 6);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 2, 3, 9, 10, 11));
}

TEST_F(TestSubstructureSpec, SymmetryGroup1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 9
          symmetry_degree: 3
          symmetry_group: 75
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 9
          symmetry_degree: 3
          symmetry_group: 75
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 6);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(2, 2, 3, 3, 4, 4));
}

TEST_F(TestSubstructureSpec, TestXor)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
        atom_properties {
          ncon: 2
          nrings: 0
          logical_operator: SS_XOR
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
          aromatic: true
        }
        single_bond: 0
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ccc(N)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(3, 4));

  _smiles = "c1ccc(CC)cc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(3, 4));

  // Both parts of the XOR are true.
  _smiles = "c1ccc(NC)cc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructureSpec, TestAtomType)
{
  _string_proto = R"(query {
      atom_type: "UST:Y"
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
          atom_type: 6007
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
          aromatic: true
          atom_type: 3001
        }
        single_bond: 0
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ccc(N)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(3, 4));
}

TEST_F(TestSubstructureSpec, TestAtomTypeGroupMatches)
{
  _string_proto = R"(query {
      atom_type: "UST:C"
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
          aryl: 1
        }
        atom_type_group: 3
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
          aryl: 1
        }
        atom_type_group: 3
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Cc1ccc(N)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0, 5));
}

TEST_F(TestSubstructureSpec, TestAtomTypeGroupNoMatch)
{
  _string_proto = R"(query {
      atom_type: "UST:C"
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
          aromatic: false
          aryl: 1
        }
        atom_type_group: 3
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
          aromatic: false
          aryl: 1
        }
        atom_type_group: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Cc1ccc(N)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructureSpec, TestAtomTypeGroupOrProblem)
{
  _string_proto = R"(query {
      atom_type: "UST:Y"
      query_atom {
        id: 0
        atom_properties {
          atomic_number: [6, 7]
          aromatic: false
        }
        atom_type_group: 3
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : [6, 7]
          aromatic: false
          aryl: 1
        }
        atom_type_group: 2
        single_bond: 0
      }
    }
  )";

//cerr << "Building query from proto\n";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CNc1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  if (0 == _query.substructure_search(_m, _sresults))
    cerr << "No matches, hit " << _query.max_query_atoms_matched_in_search() << " atoms\n";
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0, 1));

  _smiles = "NNc1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  _smiles = "CCc1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST(TestSmartsNumericQualifier, NoClosingBrace) {
  const char * s("{foobar");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_FALSE(substructure_spec::SmartsNumericQualifier(s, nchars, result));
}

TEST(TestSmartsNumericQualifier, HasMin) {
  const char * s(">2");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), nchars);
  Min_Max_Specifier<int> expected;
  expected.set_min(3);
  EXPECT_EQ(result, expected);
}

TEST(TestSmartsNumericQualifier, HasMax) {
  const char * s("<4");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), nchars);
  Min_Max_Specifier<int> expected;
  expected.set_max(3);
  EXPECT_EQ(result, expected);
}

TEST(TestSmartsNumericQualifier, HasValue) {
  const char * s("4");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), nchars);
  Min_Max_Specifier<int> expected;
  expected.add(4);
  EXPECT_EQ(result, expected);
}

TEST(TestSmartsNumericQualifier, HasValue2) {
  const char * s("46");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), nchars);
  Min_Max_Specifier<int> expected;
  expected.add(46);
  EXPECT_EQ(result, expected);
}

TEST(TestSmartsNumericQualifier, HasValue3) {
  const char * s("317");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), nchars);
  Min_Max_Specifier<int> expected;
  expected.add(317);
  EXPECT_EQ(result, expected);
}

TEST(TestSmartsNumericQualifier, NotANumber) {
  const char * s("V");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), 0);
}

TEST(TestSmartsNumericQualifier, RangeNoEndEmpty) {
  const char * s("{");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), 0);
}

TEST(TestSmartsNumericQualifier, RangeNoEndNotEmpty) {
  const char * s("{4");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), 0);
}

TEST(TestSmartsNumericQualifier, RangeMinOnly) {
  const char * s("{4-}");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), nchars);
  Min_Max_Specifier<int> expected;
  expected.set_min(4);
  EXPECT_EQ(result, expected);
}

TEST(TestSmartsNumericQualifier, RangemaxOnly) {
  const char * s("{-10}");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), nchars);
  Min_Max_Specifier<int> expected;
  expected.set_max(10);
  EXPECT_EQ(result, expected);
}

TEST(TestSmartsNumericQualifier, InvalidRange) {
  const char * s("{3-2}");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), 0);
}

TEST(TestSmartsNumericQualifier, ValidRange) {
  const char * s("{3-10}");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), nchars);
  Min_Max_Specifier<int> expected;
  expected.set_min(3);
  expected.set_max(10);
  EXPECT_EQ(result, expected);
}

TEST(TestSmartsNumericQualifier, JustANumber) {
  const char * s("46");
  const int nchars = strlen(s);
  Min_Max_Specifier<int> result;
  EXPECT_EQ(substructure_spec::SmartsNumericQualifier(s, nchars, result), nchars);
  Min_Max_Specifier<int> expected;
  expected.add(46);
  EXPECT_EQ(result, expected);
}

struct SmilesSmartsNhits {
  IWString smiles;
  IWString smarts;
  int nhits;
};

class TestRanges : public testing::TestWithParam<SmilesSmartsNhits> {
  protected:
    Substructure_Query _query;
    Molecule _m;
};

TEST_P(TestRanges, TestH) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  //cerr << "TestingH '" << params.smiles << "' smarts '" << params.smarts << " xpt " << params.nhits << '\n';
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestRanges, TestRanges, testing::Values(
  SmilesSmartsNhits{"C", "[CH]", 0},
  SmilesSmartsNhits{"C", "[CH0]", 0},
  SmilesSmartsNhits{"C", "[CH1]", 0},
  SmilesSmartsNhits{"C", "[CH2]", 0},
  SmilesSmartsNhits{"C", "[CH3]", 0},
  SmilesSmartsNhits{"C", "[CH{0-}]", 1},
  SmilesSmartsNhits{"C", "[CH{1-}]", 1},
  SmilesSmartsNhits{"C", "[CH{2-}]", 1},
  SmilesSmartsNhits{"C", "[CH{3-}]", 1},
  SmilesSmartsNhits{"C", "[CH{4-}]", 1},
  SmilesSmartsNhits{"C", "[CH{-0}]", 0},
  SmilesSmartsNhits{"C", "[CH{-1}]", 0},
  SmilesSmartsNhits{"C", "[CH{-2}]", 0},
  SmilesSmartsNhits{"C", "[CH{-3}]", 0},
  SmilesSmartsNhits{"C", "[CH{-4}]", 1},
  SmilesSmartsNhits{"C", "[CH{0-4}]", 1},
  SmilesSmartsNhits{"C", "[CH{1-4}]", 1},
  SmilesSmartsNhits{"C", "[CH{2-4}]", 1},
  SmilesSmartsNhits{"C", "[CH{3-4}]", 1},
  SmilesSmartsNhits{"C", "[CH{4-4}]", 1},
  SmilesSmartsNhits{"C", "[CH{4-5}]", 1},
  SmilesSmartsNhits{"C", "[CH{5-5}]", 0},
  SmilesSmartsNhits{"C", "[CH<1]", 0},
  SmilesSmartsNhits{"C", "[CH<2]", 0},
  SmilesSmartsNhits{"C", "[CH<3]", 0},
  SmilesSmartsNhits{"C", "[CH<4]", 0},
  SmilesSmartsNhits{"C", "[CH<5]", 1},
  SmilesSmartsNhits{"C", "[CH>0]", 1},
  SmilesSmartsNhits{"C", "[CH>1]", 1},
  SmilesSmartsNhits{"C", "[CH>2]", 1},
  SmilesSmartsNhits{"C", "[CH>3]", 1},
  SmilesSmartsNhits{"C", "[CH>4]", 0}
));

TEST_P(TestRanges, TestD) {
  const auto params = GetParam();
//cerr << "TestingD '" << params.smiles << "' smarts '" << params.smarts << " xpt " << params.nhits << '\n';
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestRangesD, TestRanges, testing::Values(
  SmilesSmartsNhits{"C", "[CD0]", 1},
  SmilesSmartsNhits{"C", "[CD1]", 0},
  SmilesSmartsNhits{"C", "[CD2]", 0},
  SmilesSmartsNhits{"C", "[CD3]", 0},
  SmilesSmartsNhits{"C", "[CD>0]", 0},
  SmilesSmartsNhits{"C", "[CD>1]", 0},
  SmilesSmartsNhits{"C", "[CD<1]", 1},
  SmilesSmartsNhits{"C", "[CD<2]", 1},
  SmilesSmartsNhits{"C", "[CD{-0}]", 1},
  SmilesSmartsNhits{"C", "[CD{-1}]", 1},
  SmilesSmartsNhits{"C", "[CD{-2}]", 1},
  SmilesSmartsNhits{"C", "[CD{0-2}]", 1},
  SmilesSmartsNhits{"C", "[CD{1-2}]", 0},
  SmilesSmartsNhits{"C", "[CD{2-2}]", 0},
  SmilesSmartsNhits{"C-C", "[CD1]", 2},
  SmilesSmartsNhits{"C-C", "[CD>0]", 2},
  SmilesSmartsNhits{"C-C", "[CD>1]", 0},
  SmilesSmartsNhits{"C-C", "[CD<1]", 0},
  SmilesSmartsNhits{"C-C", "[CD<2]", 2},
  SmilesSmartsNhits{"C-C", "[CD{-1}]", 2},
  SmilesSmartsNhits{"C-C", "[CD{-2}]", 2},
  SmilesSmartsNhits{"C-C", "[CD{1-}]", 2},
  SmilesSmartsNhits{"C-C", "[CD{2-}]", 0},
  SmilesSmartsNhits{"C-C", "[CD{0-1}]", 2},
  SmilesSmartsNhits{"C-C", "[CD{0-2}]", 2},
  SmilesSmartsNhits{"C-C", "[CD{1-3}]", 2},
  SmilesSmartsNhits{"C-C", "[CD{2-3}]", 0},
  SmilesSmartsNhits{"C-C-C", "[CD0]", 0},
  SmilesSmartsNhits{"C-C-C", "[CD{0-}]", 3},
  SmilesSmartsNhits{"C-C-C", "[CD{0-}]", 3},
  SmilesSmartsNhits{"C-C-C", "[CD{1-}]", 3},
  SmilesSmartsNhits{"C-C-C", "[CD{2-}]", 1},
  SmilesSmartsNhits{"C-C-C", "[CD{1-2}]", 3},
  SmilesSmartsNhits{"C-C-C", "[CD{1-1}]", 2},
  SmilesSmartsNhits{"C-C-C", "[CD{2-2}]", 1},
  SmilesSmartsNhits{"C-C-C", "[CD<2]", 2},
  SmilesSmartsNhits{"C-C-C", "[CD>1]", 1},
  SmilesSmartsNhits{"C-C-C", "[CD1,CD2]", 3},
  SmilesSmartsNhits{"C-C-C", "[C;D1,D2]", 3},
  SmilesSmartsNhits{"C-C-C", "[CD{1-2},CD{1-2}]", 3}
));

TEST_P(TestRanges, TestR) {
  const auto params = GetParam();
//cerr << "TestingR '" << params.smiles << "' smarts '" << params.smarts << " xpt " << params.nhits << '\n';
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestRangesR, TestRanges, testing::Values(
  SmilesSmartsNhits{"C", "[CR]", 0},
  SmilesSmartsNhits{"C", "[CR0]", 1},
  SmilesSmartsNhits{"C", "[CR>0]", 0},
  SmilesSmartsNhits{"C", "[CR<1]", 1},
  SmilesSmartsNhits{"C", "[CR{-1}]", 1},
  SmilesSmartsNhits{"C", "[CR{1-}]", 0},
  SmilesSmartsNhits{"C", "[CR{0-0}]", 1},
  SmilesSmartsNhits{"C1CC1", "[CR{0-0}]", 0},
  SmilesSmartsNhits{"C1CC1", "[CR{1-1}]", 3},
  SmilesSmartsNhits{"C1CC1", "[CR{1-}]", 3},
  SmilesSmartsNhits{"C1CC1", "[CR0]", 0},
  SmilesSmartsNhits{"C1CC1", "[CR]", 3},
  SmilesSmartsNhits{"C1CC1", "[CR1]", 3},
  SmilesSmartsNhits{"C1CC1", "[CR>0]", 3},
  SmilesSmartsNhits{"C1CC1", "[CR>1]", 0},
  SmilesSmartsNhits{"C1CC1", "[CR<1]", 0},
  SmilesSmartsNhits{"C1CC1", "[CR<2]", 3}
));

TEST_P(TestRanges, Testr) {
  const auto params = GetParam();
//cerr << "Testingr '" << params.smiles << "' smarts '" << params.smarts << " xpt " << params.nhits << '\n';
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestRangesr, TestRanges, testing::Values(
  SmilesSmartsNhits{"C", "[r]", 0},
  SmilesSmartsNhits{"C", "[r>2]", 0},
  SmilesSmartsNhits{"C", "[r{3-}]", 0},
  SmilesSmartsNhits{"C1CC1", "[r]", 3},
  SmilesSmartsNhits{"C1CC1", "[r3]", 3},
  SmilesSmartsNhits{"C1CC1", "[r>2]", 3},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r3]", 3},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r4]", 4},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r5]", 5},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r>3]", 7},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r>4]", 5},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r<5]", 5},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r{3-3}]", 3},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r{4-4}]", 4},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r{5-5}]", 5},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r{3-}]", 8},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r{4-}]", 7},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r{5-}]", 5},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r{3-4}]", 5},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r{-4}]", 5},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r{-5}]", 8},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[r{-6}]", 8}
));

TEST_P(TestRanges, Testx) {
  const auto params = GetParam();
//cerr << "Testingx '" << params.smiles << "' smarts '" << params.smarts << " xpt " << params.nhits << '\n';
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestRangesx, TestRanges, testing::Values(
  SmilesSmartsNhits{"C", "[x]", 0},
  SmilesSmartsNhits{"C", "[x0]", 1},
  SmilesSmartsNhits{"C", "[x{0-}]", 1},
  SmilesSmartsNhits{"C1CC1", "[x2]", 3},
  SmilesSmartsNhits{"C1CC1", "[x3]", 0},
  SmilesSmartsNhits{"C1CC1", "[x>2]", 0},
  SmilesSmartsNhits{"C1CC1", "[x<3]", 3},
  SmilesSmartsNhits{"C1CC1", "[r{3-}]", 3},
  SmilesSmartsNhits{"C1CC1", "[r{-3}]", 3},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[x3]", 4},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[x2]", 4},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[x>1]", 8},
  SmilesSmartsNhits{"C12CC2C3CCCC31", "[x{2-}]", 8}
));

TEST_P(TestRanges, Testv) {
  const auto params = GetParam();
//cerr << "Testingv '" << params.smiles << "' smarts '" << params.smarts << " xpt " << params.nhits << '\n';
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestRangesv, TestRanges, testing::Values(
  SmilesSmartsNhits{"C", "[v]", 0},
  SmilesSmartsNhits{"CC", "[v]", 0},
  SmilesSmartsNhits{"CC", "[v4]", 2},
  SmilesSmartsNhits{"C=C", "[v2]", 0},
  SmilesSmartsNhits{"C=C", "[v3]", 0},
  SmilesSmartsNhits{"C#C", "[v3]", 0},
  SmilesSmartsNhits{"C#C", "[v4]", 2},
  SmilesSmartsNhits{"C", "[v0]", 0},
  SmilesSmartsNhits{"C", "[v{0-}]", 1},
  SmilesSmartsNhits{"CC", "[v{0-}]", 2},
  SmilesSmartsNhits{"CC", "[v{1-}]", 2},
  SmilesSmartsNhits{"C=C", "[v{1-}]", 2},
  SmilesSmartsNhits{"C=C", "[v{2-}]", 2},
  SmilesSmartsNhits{"C=C", "[v{3-}]", 2},
  SmilesSmartsNhits{"C#C", "[v{3-}]", 2}
));

TEST_P(TestRanges, TestG) {
  const auto params = GetParam();
//cerr << "TestingG '" << params.smiles << "' smarts '" << params.smarts << " xpt " << params.nhits << '\n';
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestRangesG, TestRanges, testing::Values(
  SmilesSmartsNhits{"C", "[G0]", 1},
  SmilesSmartsNhits{"CC", "[G0]", 2},
  SmilesSmartsNhits{"C=C", "[G1]", 2},
  SmilesSmartsNhits{"C#C", "[G2]", 2},
  SmilesSmartsNhits{"C", "[G{0-}]", 1},
  SmilesSmartsNhits{"C", "[G{1-}]", 0},
  SmilesSmartsNhits{"C=C", "[G{1-}]", 2},
  SmilesSmartsNhits{"C=C", "[G{2-}]", 0},
  SmilesSmartsNhits{"C#C", "[G{2-}]", 2}
));

// Got through this a fair way and realised that these instantiations
// seemed to be instantiating many instances of the original class.
// Seems I should have made a separate class for each test??
// Not worried about it now, but next time...
class TestRangesT : public testing::TestWithParam<SmilesSmartsNhits> {
  protected:
    Substructure_Query _query;
    Molecule _m;
};

TEST_P(TestRangesT, TestT) {
  const auto params = GetParam();
//cerr << "TestingT '" << params.smiles << "' smarts '" << params.smarts << " xpt " << params.nhits << '\n';
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestRangesT, TestRangesT, testing::Values(
  SmilesSmartsNhits{"C", "[T0]", 1},
  SmilesSmartsNhits{"CC", "[T0]", 2},
  SmilesSmartsNhits{"CCC", "[T0]", 3},
  SmilesSmartsNhits{"CCO", "[T0]", 2},
  SmilesSmartsNhits{"CCO", "[T1]", 1},
  SmilesSmartsNhits{"OO", "[T0]", 0},
  SmilesSmartsNhits{"OO", "[T1]", 2},
  SmilesSmartsNhits{"OC(N)O", "[T0]", 3},
  SmilesSmartsNhits{"OC(N)O", "[T1]", 0},
  SmilesSmartsNhits{"OC(N)O", "[T2]", 0},
  SmilesSmartsNhits{"OC(N)O", "[T3]", 1}
));

class TestRangesSpiro : public testing::TestWithParam<SmilesSmartsNhits> {
  protected:
    Substructure_Query _query;
    Molecule _m;
};

TEST_P(TestRangesSpiro, TestSpiro) {
  const auto params = GetParam();
  // std::cerr << "TestingSpiro '" << params.smiles << "' smarts '" << params.smarts << " xpt " << params.nhits << '\n';
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestRangesSpiro, TestRangesSpiro, testing::Values(
  SmilesSmartsNhits{"C", "[/IWspiro]", 0},
  SmilesSmartsNhits{"C1CC1C1CC1", "[/IWspiroC]", 0},
  SmilesSmartsNhits{"C1CC12CC2", "[/IWspiroC]", 1},
  SmilesSmartsNhits{"C12(N3CC4(C)C(=O)C(C)(C3)CN1C4)C1=C(C=CC(=C1)Br)NC2=O", "[/IWspiroC]", 1}
));

class TestInvalidSmarts : public testing::TestWithParam<IWString> {
  protected:
    Substructure_Query _query;
};

TEST_P(TestInvalidSmarts, TestInvalidSmarts) {
  const auto params = GetParam();
//cerr << "Testing bad smarts " << params << "\n";
  EXPECT_FALSE(_query.create_from_smarts(params));
}
INSTANTIATE_TEST_SUITE_P(TestInvalidSmarts, TestInvalidSmarts, testing::Values(
  IWString{"[CD]"},
  IWString{"[r1]"},
  IWString{"[r2]"},
  IWString{"[r<2]"},
  IWString{"[r<3]"},
  IWString{"[r0]"},
  IWString{"[r{0-}]"},
  IWString{"[r{1-}]"},
  IWString{"[r{2-}]"},
  IWString{"[r{0-0}]"},
  IWString{"[D{}]"},
  IWString{"[D{ }]"},
  IWString{"[D{3}]"},
  IWString{"[D{q-}]"},
  IWString{"[D{-q}]"},
  IWString{"[D{-3 }]"},
  IWString{"[D{ 3-}]"},
  IWString{"[G]"}
));

class TestCipStereo : public testing::TestWithParam<SmilesSmartsNhits> {
  protected:
    Substructure_Query _query;
    Molecule _m;
};

TEST_P(TestCipStereo, TestCipStereo) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_query.create_from_smarts(params.smarts));
  EXPECT_EQ(_query.substructure_search(&_m), params.nhits);
}
INSTANTIATE_TEST_SUITE_P(TestCipStereo, TestCipStereo, testing::Values(
  SmilesSmartsNhits{"I[C@H](Br)F",  "[/IWcipS]", 1},
  SmilesSmartsNhits{"Br[C@H](F)I",  "[/IWcipS]", 1},
  SmilesSmartsNhits{"I[C@H](Br)F",  "[/IWcipS]", 1},
  SmilesSmartsNhits{"Br[C@@H](I)F", "[/IWcipS]", 1},
  SmilesSmartsNhits{"F[C@H](I)Br",  "[/IWcipS]", 1},
  SmilesSmartsNhits{"I[C@@H](F)Br", "[/IWcipS]", 1},
  SmilesSmartsNhits{"F[C@@H](Br)I", "[/IWcipS]", 1},
  SmilesSmartsNhits{"I[C@H](Br)F",  "[/IWcipR]", 0},
  SmilesSmartsNhits{"Br[C@H](F)I",  "[/IWcipR]", 0},
  SmilesSmartsNhits{"I[C@H](Br)F",  "[/IWcipR]", 0},
  SmilesSmartsNhits{"Br[C@@H](I)F", "[/IWcipR]", 0},
  SmilesSmartsNhits{"F[C@H](I)Br",  "[/IWcipR]", 0},
  SmilesSmartsNhits{"I[C@@H](F)Br", "[/IWcipR]", 0},
  SmilesSmartsNhits{"F[C@@H](Br)I", "[/IWcipR]", 0},

  SmilesSmartsNhits{"F[C@@H](Br)I", "[/IWcipSC]", 1},
  SmilesSmartsNhits{"F[C@@H](Br)I", "[/IWcipSCD3]", 1},
  SmilesSmartsNhits{"F[C@@H](Br)I", "Br[/IWcipSCD3](F)I", 1},

  SmilesSmartsNhits{"F[C@@H](I)Br", "[/IWcipR]", 1},
  SmilesSmartsNhits{"Br[C@H](I)F",  "[/IWcipR]", 1},
  SmilesSmartsNhits{"Br[C@@H](F)I", "[/IWcipR]", 1},
  SmilesSmartsNhits{"F[C@H](Br)I",  "[/IWcipR]", 1},
  SmilesSmartsNhits{"I[C@@H](Br)F", "[/IWcipR]", 1},
  SmilesSmartsNhits{"I[C@H](F)Br",  "[/IWcipR]", 1},

  SmilesSmartsNhits{"F[C@@H](I)Br", "[/IWcipS]", 0},
  SmilesSmartsNhits{"Br[C@H](I)F",  "[/IWcipS]", 0},
  SmilesSmartsNhits{"Br[C@@H](F)I", "[/IWcipS]", 0},
  SmilesSmartsNhits{"F[C@H](Br)I",  "[/IWcipS]", 0},
  SmilesSmartsNhits{"I[C@@H](Br)F", "[/IWcipS]", 0},
  SmilesSmartsNhits{"I[C@H](F)Br",  "[/IWcipS]", 0}
));

}  // namespace
