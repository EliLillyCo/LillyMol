// Tests for aromaticity settings

#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "aromatic.h"
#include "molecule.h"
#include "smiles.h"

namespace {

using std::cerr;

struct SmilesSettingsExpected {
  IWString smiles;
  int aromaticity;
  int min_aromatic_ring_size;
  int max_aromatic_ring_size;
  int two_electrons_aromatic;

  IWString expected;
};

class TestAromaticity: public testing::TestWithParam<SmilesSettingsExpected> {
  protected:
    Molecule _m;

    void SetUp();
};

void
TestAromaticity::SetUp() {
  reset_aromatic_file_scope_variables();
}

TEST_P(TestAromaticity, Various) {
  const auto params = GetParam();
  set_global_aromaticity_type(params.aromaticity);
  set_default_unique_smiles_aromaticity(params.aromaticity);
  set_min_aromatic_ring_size(params.min_aromatic_ring_size);
  set_max_aromatic_ring_size(params.max_aromatic_ring_size);
  set_allow_two_electron_systems_to_be_aromatic(params.two_electrons_aromatic);
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(_m.unique_smiles(), params.expected);
}
INSTANTIATE_TEST_SUITE_P(TestSettings, TestAromaticity, testing::Values(
  SmilesSettingsExpected{"O=C1C=CC1=O", Daylight, 4, 6, 0, "O=C1C(=O)C=C1"},
  SmilesSettingsExpected{"O=C1C=CC1=O", Daylight, 4, 6, 1, "O=c1c(=O)cc1"},
  SmilesSettingsExpected{"O=C1C=CC1=O", Daylight, 5, 6, 1, "O=C1C(=O)C=C1"},
  SmilesSettingsExpected{"N1=CC=CC=CC=C1", Daylight, 5, 6, 1, "N1=CC=CC=CC=C1"},
  SmilesSettingsExpected{"N1=CC=CC=CC=C1", Daylight, 5, 8, 1, "N1=CC=CC=CC=C1"},
  SmilesSettingsExpected{"N1=CC=CC=CC=C1", ANY_EVEN_NUMBER_OF_PI_ELECTRONS, 5, 9, 1, "[n]1ccccccc1"},
  SmilesSettingsExpected{"N1=CC=CC=CC=C1", EVERYTHING_HAS_A_PI_ELECTRON, 5, 9, 1, "[n]1ccccccc1"}
));

class TestRoundTrip : public testing::TestWithParam<IWString> {
  protected:
    Molecule _m1;
    Molecule _m2;
};

//  Interesting smiles that have been found various places.
TEST_P(TestRoundTrip, Various) {
  const auto params = GetParam();
  ASSERT_TRUE(_m1.build_from_smiles(params)) << "Cannot parse " << params;
  EXPECT_TRUE(_m2.build_from_smiles(_m1.unique_smiles()));
  EXPECT_EQ(_m1.unique_smiles(), _m2.unique_smiles()) << " not equal " <<
                _m1.unique_smiles() << '\n';
                _m2.unique_smiles();
}
INSTANTIATE_TEST_SUITE_P(TestRoundTrip, TestRoundTrip, testing::Values(
  // https://issueexplorer.com/issue/rdkit/rdkit/4701
  IWString{"[C:12]12([CH:62]([CH3:65])[c:61]3[cH:64][cH:67][cH:68][cH:66][cH:63]3)[CH:20]4[c:30]5[c:40]6[c:49]7[c:57]8[c:60]([c:59]9[c:55]([c:47]([c:44]([c:52]9[c:51]([c:43]%10[c:35]%11[c:25]%12[c:19]%13%14)[c:53]8[c:45]%11[c:39]6[c:29]4%13)[c:34]([c:24]%15[c:15]%16[c:7]%17[c:3]%18%19)[c:33]%10[c:23]%16[c:16]%12[c:8]%18[c:11]%14[c:5]1%20)[c:37]([c:36]%21[c:26]%22[c:18]%23[c:10]%24[c:13]%25[c:6]%26%27)[c:27]%15[c:17]%22[c:9]%17[c:4]%24[c:1]%19[c:2]%20%26)[c:54]([c:46]%21[c:38]%28[c:28]%23[c:21]%25%29)[c:56]%30[c:48]%28[c:41]%31[c:31]%29[c:22]%32[c:14]2%27)[c:58]%30[c:50]7[c:42]%31[c:32]5%32"},
  IWString{"[CH3:36][c:30]1[cH:35][cH:41][cH:37][c:31]([c:27]1[N+:26]2=[C:29]([CH-:33][C:40](=[N+:32]([Al+:28]2)[c:39]3[c:42]([cH:45][cH:49][cH:47][c:43]3[CH3:48])[CH3:46])[CH3:44])[CH3:34])[CH3:38]"},
  IWString{"C1Cn2cn[nH]c2=N1"},  // From GBD
  IWString{"N1=c2[nH][n]c[n]2CC1"}, // From GBD, input is unique form.
  IWString{"C1CN2C=NNC2=N1"},  // From GBD, input is kekule form.
  IWString{"C=c1cc(C)sc1=NC(=O)C"},  // exocyclic amide test
  IWString{"Nc1cnc2cn(C)nc2n1"}  // initially finding a bad Kekule form
));

// There are smiles that are really hard to read, but which fail unique smiles tests.
// Just test to see if they can be interpreted, hopefully later fix problems and
// move them to TestRoundTrip
class TestReadOnly : public testing::TestWithParam<IWString> {
  protected:
    Molecule _m;
};

TEST_P(TestReadOnly, Various) {
  const auto& params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params));
}
INSTANTIATE_TEST_SUITE_P(TestReadOnly, TestReadOnly, testing::Values(
  IWString{"Clc1c2c(c(Cl)c(c1Cl)Cl)c1[n]3[Cu][n]4c(=NC5=NC(=Nc23)c2c(Cl)c(c(c(c25)Cl)Cl)Cl)c2c(c4=NC3=NC(=N1)c1c(Cl)c(c(c(c13)Cl)Cl)Cl)c(Cl)c(c(c2Cl)Cl)Cl"}
));

#ifdef FIGURE_THESE_OUT
these are currently problematic, not sure what to do yet...
O=[3SH]1=c2c(=NC(=N1)N1CCN(C(=O)c3ccccc3)CC1)cc[3cH]c2
O=C(c1ccccc1)N1CCN(C2=N[3SH](=O)=c3c(=N2)cc[3cH]c3)CC1
#endif

TEST(TestReadAromaticSulphur, TestReadAromaticSulphur) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("[sH]1ccccc1"));
}

}  //namespace
