// Tester for minor_changes

#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Tools/minor_changes.h"

namespace {

using namespace google::protobuf;
using testing::UnorderedElementsAreArray;

class TestMinorChanges : public testing::Test {
  protected:
    minor_changes::Options _options;

    std::string _string_proto;

    Molecule _m;

    minor_changes_data::MinorChangesData _proto;

    resizable_array_p<Molecule> _results;
};

TEST_F(TestMinorChanges, TestNotActivated) {
  ASSERT_TRUE(_m.build_from_smiles("CCC"));

  EXPECT_GT(_options.Process(_m, _results), 0);
}

TEST_F(TestMinorChanges, TestSingleToDouble_1) {
  ASSERT_TRUE(_m.build_from_smiles("CCC"));

  _string_proto = R"pb(
    single_to_double_bond: true
)pb";

  ASSERT_TRUE(TextFormat::ParseFromString(_string_proto, &_proto));

  _options.SetConfig(_proto);

  EXPECT_EQ(_options.Process(_m, _results), 1);

  EXPECT_EQ(_results.size(), 1);

  EXPECT_EQ(_results[0]->smiles(), "C=CC");
}

TEST_F(TestMinorChanges, TestDoubleTosingle_1) {
  ASSERT_TRUE(_m.build_from_smiles("CC=C"));

  _string_proto = R"pb(
    double_to_single_bond: true
)pb";

  ASSERT_TRUE(TextFormat::ParseFromString(_string_proto, &_proto));

  _options.SetConfig(_proto);

  EXPECT_EQ(_options.Process(_m, _results), 1);

  ASSERT_EQ(_results.size(), 1);

  EXPECT_EQ(_results[0]->smiles(), "CCC");
}

TEST_F(TestMinorChanges, TestNoProcessAmide) {
  ASSERT_TRUE(_m.build_from_smiles("CC(=O)N"));

  _string_proto = R"pb(
    double_to_single_bond: true
)pb";

  ASSERT_TRUE(TextFormat::ParseFromString(_string_proto, &_proto));

  _options.SetConfig(_proto);

  EXPECT_EQ(_options.Process(_m, _results), 0);

  ASSERT_EQ(_results.size(), 0);
}

struct InputData {
  IWString smiles;
  std::string proto;
  std::vector<IWString> expected;
};

class TestMinorChangesP: public testing::TestWithParam<InputData> {
  protected:
    minor_changes::Options _options;

    std::string _string_proto;

    Molecule _mol;

    minor_changes_data::MinorChangesData _proto;

    resizable_array_p<Molecule> _results;
};
TEST_P(TestMinorChangesP, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  ASSERT_TRUE(TextFormat::ParseFromString(params.proto, &_proto));
  std::cerr << "Starting with " << params.smiles << '\n';
  _options.SetConfig(_proto);
  EXPECT_EQ(_options.Process(_mol, _results), params.expected.size());
  if (params.expected.empty()) {
    return;
  }

  std::vector<IWString> smiles;
  for (Molecule* m : _results) {
    smiles.push_back(m->unique_smiles());
    std::cerr << "Added " << m->unique_smiles() << '\n';
  }
  EXPECT_THAT(smiles, UnorderedElementsAreArray(params.expected));
}
INSTANTIATE_TEST_SUITE_P(TestMinorChangesP, TestMinorChangesP, testing::Values(
  // Do not lower a double bond in an amide. #0
  InputData{"CC(=O)N", R"pb(
    double_to_single_bond: true
)pb",
    {}},

  // Do not lower a double bond in a sulfonamide. #1
  InputData{"CS(=O)(=O)N", R"pb(
    double_to_single_bond: true
)pb",
    {}},

  // Do not lower a double bond in a nitro. #2
  InputData{"CN(=O)(=O)C", R"pb(
    double_to_single_bond: true
)pb",
    {}},

  // Do not lower a double bond in a sulphone. #3
  InputData{"CS(=O)(=O)C", R"pb(
    double_to_single_bond: true
)pb",
    {}},

  // Something with no double bonds should be unchanged. #4
  InputData{"CCC", R"pb(
    double_to_single_bond: true
)pb",
    {}},

  // Triple bonds always get lowered. #5
  InputData{"CC#C", R"pb(
    double_to_single_bond: true
)pb",
    {"C=CC"}},

  // Triple bonds always get lowered. #6
  InputData{"CC#N", R"pb(
    double_to_single_bond: true
)pb",
    {"N=CC"}},

  // Do not lower an amidine. #7
  InputData{"CC(=N)N", R"pb(
    double_to_single_bond: true
)pb",
    {}},

  // unspiro. #8
  InputData{"C1CC12CC2", R"pb(
    unspiro: true
)pb",
    {"C1CC1C1CC1"}},

  // make three ring OK. #9
  InputData{"CCCC", R"pb(
    make_three_membered_rings: true
)pb",
    {"CC1CC1"}},

  // No three membered rings with heteroatoms. #10
  InputData{"COCC", R"pb(
    make_three_membered_rings: true
)pb",
    {}},

  // Change carbon to Nitrogen, heteratom prevents. #11
  InputData{"COCC", R"pb(
    change_carbon_to_nitrogen: true
)pb",
    {}},

  // Change carbon to Nitrogen. #12
  InputData{"CCCC", R"pb(
    change_carbon_to_nitrogen: true
)pb",
    {"NCCC", "CCNC"}},

  // Change carbon to Nitrogen, not in a ring. #13
  InputData{"C1CC1", R"pb(
    change_carbon_to_nitrogen: true
)pb",
    {}},

  // Change carbon to Nitrogen, aromatic ring. #14
  InputData{"C1=CC=CC=C1", R"pb(
    change_carbon_to_nitrogen: true
)pb",
    {"[n]1ccccc1"}},

  // Change carbon to Nitrogen, aromatic ring, but not 3 connected. #15
  InputData{"C1=CC=CC=C1C", R"pb(
    change_carbon_to_nitrogen: true
)pb",
    {"Cc1[n]cccc1", "Cc1ccc[n]c1", "Cc1cc[n]cc1", "Nc1ccccc1"}},

  // Change carbon to Oxygen. #16
  InputData{"CC", R"pb(
    change_carbon_to_oxygen: true
)pb",
    {"OC"}},

  // Change carbon to Oxygen, heteroatom inhibits. #17
  InputData{"CNC", R"pb(
    change_carbon_to_oxygen: true
)pb",
    {}},

  // Change carbon to Oxygen, beta heteroatom inhibits. #18
  InputData{"CCN", R"pb(
    change_carbon_to_oxygen: true
)pb",
    {}},

  // Change carbon to Oxygen, two connected only. #19
  InputData{"CC(C)C", R"pb(
    change_carbon_to_oxygen: true
)pb",
    {"OC(C)C"}},

  // Change Nitrogen to Carbon. #20
  InputData{"NC", R"pb(
    change_nitrogen_to_carbon: true
)pb",
    {"CC"}},

  // Change Nitrogen to Carbon. Quats neutralised #21
  InputData{"C[N+](C)(C)C", R"pb(
    change_nitrogen_to_carbon: true
)pb",
    {"CC(C)(C)C"}},

  // Change Nitrogen to Carbon. Avoid exocyclic to arom #22
  InputData{"N=C1NC=CC=C1", R"pb(
    change_nitrogen_to_carbon: true
)pb",
    {}},

  // Insert CH2. #23
  InputData{"CC", R"pb(
    insert_ch2: true
)pb",
    {"CCC"}},

  // Remove CH2. #24
  InputData{"CCC", R"pb(
    remove_ch2: true
)pb",
    {"CC"}},

  // Remove CH2. Heteratoms inhibit. #25
  InputData{"NCN", R"pb(
    remove_ch2: true
)pb",
    {}},

  // Remove CH2. Ring membership inhibits. #26
  InputData{"C1CCC1", R"pb(
    remove_ch2: true
)pb",
    {}},

  // Destroy aromatic rings. not aromatic #27
  InputData{"C1CCC1", R"pb(
    destroy_aromatic_rings: true
)pb",
    {}},

  // Destroy aromatic rings. works #28
  InputData{"c1ncccc1", R"pb(
    destroy_aromatic_rings: true
)pb",
    {"N1CCCCC1"}},

  // Destroy aromatic rings. each ring separately #29
  InputData{"c1cccc2ncnc12", R"pb(
    destroy_aromatic_rings: true
)pb",
    {"C1Nc2ccccc2N1", "N1=CNC2C1CCCC2"}},

  // Destroy aromatic ring systems. isolated ring #30
  InputData{"c1ncccc1", R"pb(
    destroy_aromatic_ring_systems: true
)pb",
    {"N1CCCCC1"}},

  // Destroy aromatic ring systems. benzimidazole #31
  InputData{"c1cccc2ncnc12", R"pb(
    destroy_aromatic_ring_systems: true
)pb",
    {"C1NC2CCCCC2N1"}},

  // Swap adjacent atoms. parent not regenerated #32
  InputData{"CC", R"pb(
    swap_adjacent_atoms: true
)pb",
    {}},

  // Swap adjacent atoms. same not regenerated #33
  InputData{"CCC", R"pb(
    swap_adjacent_atoms: true
)pb",
    {}},

  // Swap adjacent atoms. works terminal atom #34
  InputData{"COC", R"pb(
    swap_adjacent_atoms: true
)pb",
    {"OCC"}},

  // Swap adjacent atoms. works internal atom, no adj heteroatoms #35
  InputData{"COCN", R"pb(
    swap_adjacent_atoms: true
)pb",
    {"OCCN"}},

  // single_to_double_bond. do not change 5 valent N #36
  InputData{"N1(=C(C)CCC1(C)C)=O CHEMBL325242", R"pb(
    double_to_single_bond: true
)pb",
    {}},

  // destroy_aromatic_ring. neutralize charges #37
  InputData{"[N+]1(=C(C=C(SC)SC)C=CC=C1C)C", R"pb(
    destroy_aromatic_rings: true
)pb",
    {"CC1N(C)C(C=C(SC)SC)CCC1"}},

  // destroy_aromatic_ring_systems. neutralize charges #38
  InputData{"[N+]1(=C(C=C(SC)SC)C=CC=C1C)C", R"pb(
    destroy_aromatic_ring_systems: true
)pb",
    {"CC1N(C)C(C=C(SC)SC)CCC1"}},

  // destroy_aromatic_ring_systems. no process n=O groups #39
  InputData{"C1=CC2=N(=O)N=C(N(=C2C=C1)=O)NCC CHEMBL214041", R"pb(
    destroy_aromatic_rings: true
)pb",
    {"CCNc1[n][n](=O)c2c([n]1=O)CCCC2"}},

  // destroy_aromatic_ring_systems. keep charge in unchanged ring #40
  InputData{"C1=C[N+]2=C(N=C1)N(CCCC)[C@](O)(C1=CC=CC=C1)C2 CHEMBL1851918", R"pb(
    destroy_aromatic_rings: true
)pb",
    {"CCCCN1C2=[N+](CCCN2)C[C@@]1(O)c1ccccc1", "CCCCN1c2[n]ccc[n+]2C[C@@]1(O)C1CCCCC1"}},

  // destroy_aromatic_ring_systems. no change five valent N #41
  InputData{"N1=C(C2=CC=CC=C2)N=C2C(=N(=O)C3=C(N2C(C)C)C=CC=C3)C1=O", R"pb(
    destroy_aromatic_rings: true
)pb",
    {"O=c1[n]c([n]c2[n](C(C)C)c3c([n](=O)c12)cccc3)C1CCCCC1", "O=c1[n]c([n]c2[n](C(C)C)c3c([n](=O)c12)CCCC3)c1ccccc1", "O=C1NC(NC2N(C(C)C)c3c(N(=O)=C12)cccc3)c1ccccc1"}},

  // make_three_membered_rings. complex existing 3 membered ring #42
  InputData{"C12C3C(CC(C1)CC3)CCC2", R"pb(
    make_three_membered_rings: true
)pb",
    {"C1CC23C4CCC(C2)CC34C1", "C1CC2CC34CCC2C3(C1)C4", "C1CC2C3C4CC(CC34C1)C2",
     "C1C2CC34CC3CC1C4CC2", "C1CC2C3C4CCC23C(C4)C1", "C1CC2C34CC3C(C2)CC4C1",
     "C1C2CC3C4(C1)C(CC3)C4C2", "C1CC2C3CCC4C2C4C3C1", "C1CC2CC3C4C3CC2C4C1",
     "C1C2C3C4C(CC2)C4C1CC3", "C1CC2C3C4CC4(C2)CC3C1", "C1C2C3C4CCC(C3)CC4C12"}},

  // add_fragments. simple case #43
  InputData{"CC", R"pb(
    add_fragments: true
    fragment: "S"
)pb",
    {"SCC"}},

  // add_fragments. with atom typing matches #44
  InputData{"CC", R"pb(
    add_fragments: true
    fragment: "[3001SH2]"
    atype: "UST:AY"
)pb",
    {"[3001SH]CC"}},

  // insert_fragments. fails because not enough H on S #44
  InputData{"CC", R"pb(
    insert_fragments: true
    bivalent_fragment: "[3001SH2]C"
    atype: "UST:AY"
)pb",
    {}},

  // insert_fragments. works #45
  InputData{"CC", R"pb(
    insert_fragments: true
    bivalent_fragment: "[3001CH3]S"
    atype: "UST:AY"
)pb",
    {"S[3001CH](C)C"}},

  // insert_fragments. works 1 atom inserted #46
  InputData{"CCC", R"pb(
    replace_inner_fragments: true
    bivalent_fragment: "[3001CH3]S"
    atype: "UST:AY"
)pb",
    {"S[3001CH](C)C"}},

  // insert_fragments. works 2 atoms inserted #47
  InputData{"CCCCO", R"pb(
    replace_inner_fragments: true
    bivalent_fragment: "[3001CH3](S)[3001CH3]N"
    atype: "UST:AY"
)pb",
    {"S[3001CH](C)[3001CH](N)CO", "S[3001CH](CO)[3001CH](N)C"}},

  // replace_terminal_fragments. #48
  InputData{"CCCCO", R"pb(
    replace_terminal_fragments: true
    fragment: "CCN"
    atype: "UST:AY"
)pb",
    {"OCCCCCN", "OCCCCN", "NCCCCC", "NCCCCCC"}},

  // remove_fragment. #49
  InputData{"CCCCO", R"pb(
    remove_fragment: [1, 2]
)pb",
    {"OCCC", "CCCC"}}
));

}  // namespace
