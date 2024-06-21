// Tests for dicer_api

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "dicer_api.h"

namespace dicer_api {
namespace {

using dicer_api::Dicer;
using testing::UnorderedElementsAre;
using testing::Pair;

TEST(TestEthane, TestEthaneNoBreakCarbonCarbon) {
  Dicer dicer;
  Molecule ethane;
  ASSERT_TRUE(ethane.build_from_smiles("CC ethane"));
  std::unordered_map<std::string, uint32_t> result;
  EXPECT_EQ(dicer.Dice(ethane, result), 0);
}

TEST(TestEthane, TestEthaneokBreakCarbonCarbon) {
  Dicer dicer;
  dicer.set_break_cc_bonds(1);

  Molecule ethane;
  ASSERT_TRUE(ethane.build_from_smiles("CC ethane"));

  std::unordered_map<std::string, uint32_t> result;
  ASSERT_EQ(dicer.Dice(ethane, result), 1);
  EXPECT_THAT(result, UnorderedElementsAre(Pair("C", 2)));
}

TEST(TestMethanol, TestDefault) {
  Dicer dicer;

  Molecule methanol;
  ASSERT_TRUE(methanol.build_from_smiles("CO methanol"));

  std::unordered_map<std::string, uint32_t> result;
  ASSERT_EQ(dicer.Dice(methanol, result), 2);
  EXPECT_THAT(result, UnorderedElementsAre(Pair("C", 1), Pair("O", 1)));
}

TEST(TestMethanol, TestIsotope) {
  Dicer dicer;
  constexpr isotope_t kIsotope = 4;

  dicer.set_label_join_points(kIsotope);

  Molecule methanol;
  ASSERT_TRUE(methanol.build_from_smiles("CO methanol"));

  std::unordered_map<std::string, uint32_t> result;
  ASSERT_EQ(dicer.Dice(methanol, result), 2);
  EXPECT_THAT(result, UnorderedElementsAre(Pair("[4CH4]", 1), Pair("[4OH2]", 1)));
}

class TestDicerApi : public testing::Test {
  protected:
    Dicer _dicer;
    Molecule _m;
    std::unordered_map<std::string, uint32_t> _result;
};

TEST_F(TestDicerApi, TestEthanolNoMinFragSize) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);

  ASSERT_TRUE(_m.build_from_smiles("CCO methanol"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 2);
  EXPECT_THAT(_result, UnorderedElementsAre(Pair("[8CH3]C", 1), Pair("[8OH2]", 1)));
}

TEST_F(TestDicerApi, TestEthanolMinFragSize) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_min_fragment_size(2);

  ASSERT_TRUE(_m.build_from_smiles("CCO methanol"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 1);
  EXPECT_THAT(_result, UnorderedElementsAre(Pair("[8CH3]C", 1)));
}

TEST_F(TestDicerApi, TestEthanolmaxFragSize) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_max_fragment_size(1);

  ASSERT_TRUE(_m.build_from_smiles("CCO methanol"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 1);
  EXPECT_THAT(_result, UnorderedElementsAre(Pair("[8OH2]", 1)));
}

TEST_F(TestDicerApi, TestEthanolBreakCCBondsDepth1) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_break_cc_bonds(1);

  ASSERT_TRUE(_m.build_from_smiles("CCO methanol"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 4);
  EXPECT_THAT(_result, UnorderedElementsAre(Pair("[8CH4]", 1), Pair("[8CH3]C", 1), Pair("O[8CH3]", 1),
                        Pair("[8OH2]", 1)));
}

TEST_F(TestDicerApi, TestEthanolBreakCCBondsDepth2) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_break_cc_bonds(1);
  _dicer.set_max_bonds_to_break(2);

  ASSERT_TRUE(_m.build_from_smiles("CCO methanol"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 4);
  EXPECT_THAT(_result, UnorderedElementsAre(Pair("[8CH4]", 4), Pair("[8CH3]C", 1), Pair("O[8CH3]", 1),
                        Pair("[8OH2]", 2)));
}

TEST_F(TestDicerApi, TestCF3) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_break_cc_bonds(1);

  ASSERT_TRUE(_m.build_from_smiles("CC(F)(F)F x"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 2);
  EXPECT_THAT(_result, UnorderedElementsAre(Pair("[8CH4]", 1), Pair("F[8CH](F)F", 1)));
}

TEST_F(TestDicerApi, TestAmide) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_break_cc_bonds(1);

  ASSERT_TRUE(_m.build_from_smiles("C1=CC=CC=C1C(=O)NC x"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 4);
  EXPECT_THAT(_result, UnorderedElementsAre(Pair("[8cH]1ccccc1", 1), Pair("O=[8CH]NC", 1),
                Pair("[8CH4]", 1), Pair("O=C([8NH2])c1ccccc1", 1)));
}

TEST_F(TestDicerApi, TestSulfonAmide) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_break_cc_bonds(1);

  ASSERT_TRUE(_m.build_from_smiles("C1=CC=CC=C1S(=O)(=O)NC x"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 4);
  EXPECT_THAT(_result, UnorderedElementsAre(Pair("[8cH]1ccccc1", 1), Pair("CN[8SH](=O)=O", 1),
                Pair("[8CH4]", 1), Pair("[8NH2]S(=O)(=O)c1ccccc1", 1)));
}

TEST_F(TestDicerApi, TestGuanidine) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_break_cc_bonds(1);

  ASSERT_TRUE(_m.build_from_smiles("CC(=N)(N)NC x"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 3);
  EXPECT_THAT(_result, UnorderedElementsAre(Pair("[8NH2]C(=N)(N)C", 1), Pair("N[8C](=N)NC", 1),
                Pair("[8CH4]", 2)));
}

TEST_F(TestDicerApi, Test1AminoEthanol) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_break_cc_bonds(1);

  ASSERT_TRUE(_m.build_from_smiles("C(C)(O)N x"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 6);
  EXPECT_THAT(_result, UnorderedElementsAre(Pair("O[8CH2]C", 1),
                        Pair("[8OH2]", 1),
                        Pair("N[8CH2]C", 1),
                        Pair("[8CH4]", 1),
                        Pair("[8NH3]", 1),
                        Pair("O[8CH2]N", 1)));
}

TEST_F(TestDicerApi, Test1AminoEthanol2Bonds) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_break_cc_bonds(1);
  _dicer.set_max_bonds_to_break(2);

  ASSERT_TRUE(_m.build_from_smiles("C(C)(O)N x"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 9);
  EXPECT_THAT(_result, UnorderedElementsAre(
      Pair("O[8CH2]C", 1), 
      Pair("N[8CH2]C", 1), 
      Pair("[8CH4]", 3), 
      Pair("O[8CH3]", 2), 
      Pair("[8OH2]", 3), 
      Pair("[8CH3]C", 2), 
      Pair("N[8CH3]", 2), 
      Pair("[8NH3]", 3), 
      Pair("O[8CH2]N", 1)));
}

TEST_F(TestDicerApi, Test1AminoEthanol3Bonds) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_break_cc_bonds(1);
  _dicer.set_max_bonds_to_break(3);

  ASSERT_TRUE(_m.build_from_smiles("C(C)(O)N x"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 9);
  EXPECT_THAT(_result, UnorderedElementsAre(
      Pair("O[8CH2]C", 1), 
      Pair("N[8CH2]C", 1), 
      Pair("[8CH4]", 11), 
      Pair("O[8CH3]", 2), 
      Pair("[8OH2]", 5), 
      Pair("[8CH3]C", 2), 
      Pair("N[8CH3]", 2), 
      Pair("[8NH3]", 5), 
      Pair("O[8CH2]N", 1)));
}

TEST_F(TestDicerApi, TestSmarts1) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_break_cc_bonds(1);
  _dicer.set_max_bonds_to_break(3);
  ASSERT_TRUE(_dicer.AddBreakBondSmarts("C-F"));

  ASSERT_TRUE(_m.build_from_smiles("CC(F)(F)F x"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 4);
  EXPECT_THAT(_result, UnorderedElementsAre(
    Pair("[8FH]", 15), 
    Pair("[8CH3]C", 6),
    Pair("F[8CH2]C", 6),
    Pair("C[8CH](F)F", 3)));
}

TEST_F(TestDicerApi, TestSmarts1NoSymm) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  // Note that this does not do anything. Once queries are
  // specified, the hard coded rules are not consulted.
  //_dicer.set_break_cc_bonds(1);

  _dicer.set_max_bonds_to_break(3);
  ASSERT_TRUE(_dicer.AddBreakBondSmarts("C-F"));
  ASSERT_TRUE(_dicer.AddBreakBondSmarts("C-C"));
  _dicer.set_perceive_symmetry_equivalent_matches(0);

  ASSERT_TRUE(_m.build_from_smiles("CC(F)(F)F x"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 5);
  EXPECT_THAT(_result, UnorderedElementsAre(
        Pair("C[8CH](F)F", 1),
        Pair("[8FH]", 2), 
        Pair("F[8CH2]F", 2), 
        Pair("F[8CH](F)F", 1), 
        Pair("[8CH4]", 2)));
}

TEST_F(TestDicerApi, TestAspirin1) {
  constexpr isotope_t kIsotope = 8;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_break_cc_bonds(1);

  ASSERT_TRUE(_m.build_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O aspirin"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 6);
  EXPECT_THAT(_result, UnorderedElementsAre(
        Pair("O[8CH]=O", 1), 
        Pair("O=C(Oc1[8cH]cccc1)C", 1), 
        Pair("OC(=O)c1[8cH]cccc1", 1), 
        Pair("OC(=O)c1c(O[8CH]=O)cccc1", 1), 
        Pair("[8OH]C(=O)C", 1), 
        Pair("[8CH4]", 1)));
}

TEST_F(TestDicerApi, TestRecap1) {
  constexpr isotope_t kIsotope = 1;

  _dicer.set_label_join_points(kIsotope);
  _dicer.set_work_like_recap(1);
  ASSERT_TRUE(_m.build_from_smiles("O=C(N(C)C1=C(N=C2N1C=C(C=C2)C(=O)NCCOC1=CC=C(OC)C=C1)CC)CC1=CC=CC=C1 CHEMBL2354634"));

  ASSERT_TRUE(_dicer.AddBreakBondSmarts("[R]!@*"));

  ASSERT_EQ(_dicer.Dice(_m, _result), 7);
  EXPECT_THAT(_result, UnorderedElementsAre(
        Pair("[1CH3]C", 1),
        Pair("[1OH]CCN[1CH]=O", 1),
        Pair("[1cH]1cc[1cH]cc1", 1),
        Pair("[n]1c2[n](c[1cH]cc2)[1cH][1cH]1", 1),
        Pair("O=C([1NH]C)[1CH3]", 1),
        Pair("[1cH]1ccccc1", 1),
        Pair("[1OH]C", 1)));
}

}  // namespace
}  // namespace dicer_api
