// Tester for the No_Matched_Atoms
#include <iostream>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "substructure.h"

namespace {

using std::cerr;
using std::endl;

using testing::Eq;
using testing::ElementsAre;
using testing::IsEmpty;
using testing::Property;
using testing::ResultOf;

class TestNMABToken : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    const_IWSubstring s;

    resizable_array_p<NMAB_Token> tokens;
};

void
TestNMABToken::SetUp() 
{
}

TEST_F(TestNMABToken, Empty) {
  s = "{}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, Incomplete1) {
  s = "{";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, Incomplete2) {
  s = "}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, Incomplete3) {
  s = "a}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, Incomplete4) {
  s = "{bb";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, SingleNumber) {
  s = "{3}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);

  EXPECT_EQ(tokens[0]->numbers()[0], 3);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, ElementsAre(3)));
}

TEST_F(TestNMABToken, InvalidRange1) {
  s = "{3-}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, InvalidRange2) {
  s = "{3-2}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, InvalidRange3) {
  s = "{-2}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, Range1) {
  s = "{3-3}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, ElementsAre(3)));
}

TEST_F(TestNMABToken, Range35) {
  s = "{3-5}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, ElementsAre(3, 4, 5)));
}

TEST_F(TestNMABToken, LessThan) {
  s = "{<3}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::relational, -1));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, ElementsAre(3)));
}

TEST_F(TestNMABToken, GreaterThan) {
  s = "{<0}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::relational, -1));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, ElementsAre(0)));
//EXPECT_THAT(*tokens[0], ResultOf([](const NMAB_Token& nmbt) { return nmbt.number1();}, Eq(0)));
}

TEST_F(TestNMABToken, JustOperator) {
  s = "{&}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);

  s = "{,}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);

  s = "{^}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);

  s = "{;}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, JustRelational) {
  s = "{<}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);

  s = "{>}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, JustOpAndRelational) {
  s = "{&<}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);

  s = "{;>}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, EmptySmarts) {
  s = "{[]}";
  cerr << "first\n";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
  cerr << "second\n";
  s = "{1[]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
  cerr << "third\n";
  s = "{>1[]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
  s = "{,>1[]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
  s = "{,4-5[]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, ValidSmarts1) {
  s = "{[c]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::op, Eq(0)));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::relational, 0));
// Not sure why these fail. I want to test negative, not a specific -ve value.
//EXPECT_THAT(*tokens[0], Property(&NMAB_Token::number1, testing::Lt(0)));
//EXPECT_THAT(*tokens[0], Property(&NMAB_Token::number2, testing::Lt(0)));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, IsEmpty()));

  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::smarts, Eq("c")));
}

TEST_F(TestNMABToken, ValidSmarts2) {
  s = "{[cc]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::op, Eq(0)));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::relational, 0));
// Not sure why these fail. I want to test negative, not a specific -ve value.
//EXPECT_THAT(*tokens[0], Property(&NMAB_Token::number1, testing::Lt(0)));
//EXPECT_THAT(*tokens[0], Property(&NMAB_Token::number2, testing::Lt(0)));

  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::smarts, Eq("cc")));
}

TEST_F(TestNMABToken, NoFirstOp) {
  s = "{;3,4}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, NoOperator) {
  s = "{[c]3}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, Combinations1) {
  s = "{<2[cc],3-5[n]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 2);

  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::op, Eq(0)));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::relational, -1));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, ElementsAre(2)));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::smarts, Eq("cc")));

  EXPECT_THAT(*tokens[1], Property(&NMAB_Token::op, Eq(2)));
  EXPECT_THAT(*tokens[1], Property(&NMAB_Token::relational, 0));
  EXPECT_THAT(*tokens[1], Property(&NMAB_Token::numbers, ElementsAre(3, 4, 5)));
  EXPECT_THAT(*tokens[1], Property(&NMAB_Token::smarts, Eq("n")));
}

class TestNMAB : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    Single_Substructure_Query _query;
    Molecule _m;
    Substructure_Results _sresults;
};

void
TestNMAB::SetUp()
{
}

TEST_F(TestNMAB, TestBadSmarts1) {
  EXPECT_EQ(_query.create_from_smarts("C..."), false);
}

TEST_F(TestNMAB, TestBadSmarts2) {
  EXPECT_EQ(_query.create_from_smarts("..."), false);
}

TEST_F(TestNMAB, TestBadSmarts3) {
  EXPECT_EQ(_query.create_from_smarts("...C"), false);
}

TEST_F(TestNMAB, EstersMatch) {
  ASSERT_TRUE(_query.create_from_smarts("COC(=O)...C(=O)OC"));
  ASSERT_TRUE(_m.build_from_smiles("COC(=O)CCCCCC(=O)OC"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, EstersNoMatch) {
  ASSERT_TRUE(_query.create_from_smarts("COC(=O)...C(=O)OC"));
  ASSERT_TRUE(_m.build_from_smiles("CC(=O)OCCCCCC(=O)OC"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 0);
}

TEST_F(TestNMAB, Match1) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1Distance) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{1}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1NoMatchDistance) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{2}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 0);
}

TEST_F(TestNMAB, Match1NoMatchDistanceLongPath) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{3}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 0);
}

TEST_F(TestNMAB, Match1MatchSmartsAll) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsCount) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{1[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsGeater) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{>0[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsLess) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{<3[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsNotAllSame1) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{>0[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cnc(O)cc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsNotAllSame2) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{>0[n]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cnc(O)cc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
  ASSERT_TRUE(_m.build_from_smiles("Oc1ccc(O)nc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsNotAllSame3) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{>0[n],>0[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cnc(O)cc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
  ASSERT_TRUE(_m.build_from_smiles("Oc1ccc(O)nc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsExactCount) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{1[n];1[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cnc(O)cc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
  ASSERT_TRUE(_m.build_from_smiles("Oc1ccc(O)nc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsOr) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{[n,c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cnc(O)cc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
  ASSERT_TRUE(_m.build_from_smiles("Oc1ccc(O)nc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, AtomsOrSeparation) {
  ASSERT_TRUE(_query.create_from_smarts("O...{>5,>2[C]}O"));
  ASSERT_TRUE(_m.build_from_smiles("ONNNNNNO"));
  _query.set_find_unique_embeddings_only(1);
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 1);
  ASSERT_TRUE(_m.build_from_smiles("OCNCCO"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 1);
}

TEST_F(TestNMAB, AllRingAtoms) {
  ASSERT_TRUE(_query.create_from_smarts("O...{[R]}O"));
  ASSERT_TRUE(_m.build_from_smiles("ONNNNNNO"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 0);
  ASSERT_TRUE(_m.build_from_smiles("OC1C2C(O)C12"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, AllRingAtomsButNotBonds) {
  ASSERT_TRUE(_query.create_from_smarts("O...{[R]}O"));
  ASSERT_TRUE(_m.build_from_smiles("OC1CC1C2CC2O"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, NoRingAtomsBetween) {
  ASSERT_TRUE(_query.create_from_smarts("[R]-!@*...{1;0[R]}*-!@[R]"));
  ASSERT_TRUE(_m.build_from_smiles("C1CC1CCC1CC1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 0);
}

TEST_F(TestNMAB, ZeroLengthRegion) {
  ASSERT_TRUE(_query.create_from_smarts("[>0]...{[R]}[>0]"));
  ASSERT_TRUE(_m.build_from_smiles("[1CH3][1CH2]C1CC(F)(F)C1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 0);

}

struct ProtoMolMatches {
  std::string proto;
  IWString smiles;
  int expected;
};

class TestRegionsP: public testing::TestWithParam<ProtoMolMatches> {
  protected:
    SubstructureSearch::SubstructureQuery _proto;
    Substructure_Query _query;
    Molecule _mol;
    Substructure_Results _sresults;
};

TEST_P (TestRegionsP, Tests) {
  const auto params = GetParam();
  cerr << "Calling TextFormat\n";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  cerr << "Next ConstructFromProto\n";
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  // cerr << "Testing " << params.smiles << " expecting " << params.expected << '\n';
  EXPECT_EQ(_query.substructure_search(_mol, _sresults), params.expected) << 
      "mismatch " << params.smiles << " region " << _proto.query(0).region(0).ShortDebugString();
}
INSTANTIATE_TEST_SUITE_P(TestRegionsP, TestRegionsP, testing::Values(
  ProtoMolMatches{R"pb(
query {
  smarts: "[CD1].[CD1]"
  unique_embeddings_only: true
  region {
    atom: [0, 1]
    natoms: 2
  }
}
)pb", "CCCC", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "[N].[N]"
  unique_embeddings_only: true
  region {
    atom: [0, 1]
    natoms: 2
  }
}
)pb", "CNCCNC", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "[N].[N]"
  unique_embeddings_only: true
  region {
    atom: [0, 1]
    natoms: 3
  }
}
)pb", "CCNCCCNCC", 1},

  // A failing case because the atoms are in a ring.
  ProtoMolMatches{R"pb(
query {
  smarts: "[N].[N]"
  unique_embeddings_only: true
  region {
    atom: [0, 1]
    min_natoms: 1
  }
}
)pb", "N1CCNCC1", 0},

  ProtoMolMatches{R"pb(
query {
  smarts: "[N].[N]"
  unique_embeddings_only: true
  region {
    atom: [0, 1]
    natoms: 7
  }
}
)pb", "CCNCC1C(CC1)CCNCC", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "[N].[N]"
  unique_embeddings_only: true
  region {
    atom: [0, 1]
    nrings: 1
  }
}
)pb", "CCNCC1C(C1)CNCC", 1}

));

class TestNearbyAtomsP: public testing::TestWithParam<ProtoMolMatches> {
  protected:
    SubstructureSearch::SubstructureQuery _proto;
    Substructure_Query _query;
    Molecule _mol;
    Substructure_Results _sresults;
};

TEST_P (TestNearbyAtomsP, Tests) {
  const auto params = GetParam();
  cerr << "Calling TextFormat\n";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  cerr << "Next ConstructFromProto\n";
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  // cerr << "Testing " << params.smiles << " expecting " << params.expected << '\n';
  EXPECT_EQ(_query.substructure_search(_mol, _sresults), params.expected) << 
      "mismatch " << params.smiles << " from " << _proto.ShortDebugString() << '\n';
}
INSTANTIATE_TEST_SUITE_P(TestNearbyAtomsP, TestNearbyAtomsP, testing::Values(

  ProtoMolMatches{R"pb(
query {
  smarts: "O"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    min_bonds_between: 2
    max_bonds_between: 4
  }
}
)pb", "OCCN", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "O"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    min_bonds_between: 3
    max_bonds_between: 4
  }
}
)pb", "OCCN", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "O"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    min_bonds_between: 4
    max_bonds_between: 5
  }
}
)pb", "OCCN", 0},

  ProtoMolMatches{R"pb(
query {
  smarts: "n1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    min_bonds_between: 2
    max_bonds_between: 3
    matched_atom: 0
  }
}
)pb", "n1c(CN)cccc1", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "n1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    bonds_between: 2
    matched_atom: 0
  }
}
)pb", "n1c(CN)cccc1", 0},

  ProtoMolMatches{R"pb(
query {
  smarts: "n1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    bonds_between: 3
    matched_atom: 0
  }
}
)pb", "n1c(CN)cccc1", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "n1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    bonds_between: 3
    matched_atom: 0
    hits_needed: 1
  }
}
)pb", "n1c(CN)cccc1", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "n1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    bonds_between: 3
    matched_atom: 0
    hits_needed: 2
  }
}
)pb", "n1c(CN)cccc1", 0},

  ProtoMolMatches{R"pb(
query {
  smarts: "n1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    bonds_between: 3
    matched_atom: 0
    min_hits_needed: 1
  }
}
)pb", "n1c(CN)cccc1", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "n1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    bonds_between: 3
    matched_atom: 0
    min_hits_needed: 1
    rejection: true
  }
}
)pb", "n1c(CN)cccc1", 0},

  ProtoMolMatches{R"pb(
query {
  smarts: "Fc1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "F"
    bonds_between: 1
    can_overlap_matched_atoms: false
  }
}
)pb", "Fc1ccccc1", 0},

  ProtoMolMatches{R"pb(
query {
  smarts: "Fc1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "F"
    bonds_between: 1
    can_overlap_matched_atoms: true
  }
}
)pb", "Fc1ccccc1", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "Fc1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    max_bonds_between: 2
    matched_atom: [3, 4, 5]
  }
}
)pb", "Fc1ccccc1", 0},

  ProtoMolMatches{R"pb(
query {
  smarts: "Fc1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    max_bonds_between: 2
    matched_atom: [3, 4, 5]
  }
}
)pb", "Fc1cc(N)ccc1", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "Fc1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    max_bonds_between: 2
    matched_atom: [3, 4, 5]
  }
}
)pb", "Fc1cc(CN)ccc1", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "Fc1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    max_bonds_between: 2
    matched_atom: [3, 4, 5]
  }
}
)pb", "Fc1cc(CCN)ccc1", 0},

  ProtoMolMatches{R"pb(
query {
  smarts: "Fc1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    max_bonds_between: 2
    matched_atom: [3, 4, 5]
    min_hits_needed: 2
  }
}
)pb", "Fc1cc(CN)cc(CN)c1", 1},

  ProtoMolMatches{R"pb(
query {
  smarts: "Fc1ccccc1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "N"
    max_bonds_between: 2
    matched_atom: [3, 4, 5]
    min_hits_needed: 2
    max_hits_needed: 2
  }
}
)pb", "Fc1cc(CN)c(N)c(CN)c1", 0},

  // Wow, look at what this can do.
  // Can specify that some number of the matched atoms must be a given type.
  // We always had element_hits_needed, but this seems more precise.
  ProtoMolMatches{R"pb(
query {
  smarts: "Fc1caaac1"
  unique_embeddings_only: true
  nearby_atoms {
    smarts: "n"
    max_bonds_between: 0
    matched_atom: [3, 4, 5]
    hits_needed: 2
    can_overlap_matched_atoms: true
  }
}
)pb", "Fc1cncnc1", 1}

));

}  // namespace
