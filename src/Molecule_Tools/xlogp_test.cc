// Tester for xlogp
#include <iostream>
#include <memory>
#include <span>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Tools/xlogp.h"

namespace {

using testing::ElementsAreArray;

#ifdef NOT_NEEDED__
std::ostream&
operator<< (std::ostream& os, const std::span<int>& values) {
  for (int v : values) {
    os << ' ' << v;
  }

  return os;
}
#endif

IWString
ToString(const std::span<int>& values) {
  IWString rc;
  for (int v : values) {
    rc << ' ' << v;
  }

  return rc;
}

struct SmilesExpected {
  IWString smiles;
  double expected;
  std::vector<int> atypes;
};
class TestXlogpP: public testing::TestWithParam<SmilesExpected> {
  protected:
    Molecule _m;
    std::unique_ptr<int[]> _atype;
};

TEST_P(TestXlogpP, WithoutCorrections) {
  const auto& params = GetParam();
  xlogp::ForTestingSetApplyCorrections(0);

  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  _atype.reset(new int[_m.natoms()]);
  std::optional<double> x = xlogp::XLogP(_m, _atype.get());
  ASSERT_NE(x, std::nullopt);
  EXPECT_NEAR(*x, params.expected, 0.001) << _m.smiles() << 
        " got " << *x << " expect " << params.expected;
  std::span<int> values(_atype.get(), _m.natoms());
  EXPECT_THAT(values, ElementsAreArray(params.atypes)) << _m.smiles() << ToString(values);
}
INSTANTIATE_TEST_SUITE_P(TestXlogpPNoCorrections, TestXlogpP, testing::Values(
  SmilesExpected{"CC", 1.244, {1,1}},
  SmilesExpected{"c1ccccc1", 1.962, {26,26,26,26,26,26}},
  SmilesExpected{"c1ncccc1", 0.653, {27,57,27,26,26,26}},
  SmilesExpected{"CC(F)(F)F", 1.402, {1,16,72,72,72}},
  SmilesExpected{"Fc1ccccc1", 2.064, {72, 30, 26, 26, 26, 26, 26}},
  SmilesExpected{"Clc1ccccc1", 2.581, {73, 30, 26, 26, 26, 26, 26}},
  SmilesExpected{"Brc1ccccc1", 2.758, {74, 30, 26, 26, 26, 26, 26}},
  SmilesExpected{"Ic1ccccc1", 3.06, {75, 30, 26, 26, 26, 26, 26}},
  // BioByte 1.885 Marvin 1.91
  SmilesExpected{"O=N(=O)c1ccccc1", 1.835, {80, 80, 80, 30, 26, 26, 26, 26, 26}},
  // BioByte 3.305 Marvin 3.00
  SmilesExpected{"S=C=Nc1ccccc1", 3.197, {78, 78, 78, 30, 26, 26, 26, 26, 26}},
  // BioByte 1.575 Marvin 1.83
  SmilesExpected{"N#Cc1ccccc1", 1.681, {77, 77, 29, 26, 26, 26, 26, 26}},
  // BioByte 1.885 Marvin 1.91
  SmilesExpected{"O=N(=O)c1ccccc1", 1.835, {80, 80, 80, 30, 26, 26, 26, 26, 26}},
  // Biobyte 2.025, Marvin 2.05
  SmilesExpected{"O=Nc1ccccc1", 1.648, {79, 79, 30, 26, 26, 26, 26, 26}},
  //Biobyte 3.305 Marvin 3.00
  SmilesExpected{"S=C=Nc1ccccc1", 3.197, {78, 78, 78, 30, 26, 26, 26, 26, 26}},
  // Biobyte -0.65 Marvin -0.81  POOR MATCH
  SmilesExpected{"N1NC(=O)C=C1", 0.055, {60, 64, 24, 44, 19, 20}},
  // Biobyte 2.641 Marvin 2.49
  SmilesExpected{"Cc1ccccc1", 2.243, {2, 29, 26, 26, 26, 26, 26}},
  // Biobyte -0.664 Marvin -0.63
  SmilesExpected{"CN", -0.533, {3, 46}},
  // Biobyte 0.915 Marvin 1.14
  SmilesExpected{"Nc1ccccc1", 1.214, {47, 30, 26, 26, 26, 26, 26}},
  // Biobyte -0.208 Marvin 0.11
  SmilesExpected{"O=C(C)C", 0.192, {44, 24, 2, 2}},
  // Biobyte 0.018 Marvin 0.19
  SmilesExpected{"CN(C)C", 0.314, {3, 51, 3, 3}},
  // Biobyte -1.350 Marvin -0.064
  SmilesExpected{"CN(C)NC", -0.581, {3, 52, 3, 50, 3}},
  // Biobyte -1.114 Marvin -1.03
  SmilesExpected{"CC(=O)N", -0.813, {2, 24, 44, 63}},
  // Biobyte -1.078 Marvin -0.81
  SmilesExpected{"CC(=O)NC", -0.323, {2, 24, 44, 64, 3}},
  // Biobyte -0.801 Marvin -0.58 POOR MATCH
  SmilesExpected{"CC(=O)N(C)C", -0.184, {2, 24, 44, 65, 3, 3}},
  // Biobyte 0.646 Marvin 0.84
  SmilesExpected{"CS", 0.755, {3, 66}},
  // Biobyte 1.788 Marvin 1.75
  SmilesExpected{"CSC", 0.985, {3, 67, 3}},
  // Biobyte 0.322 Marvin 1.00
  SmilesExpected{"S=C(C)C", 0.432, {69, 24, 2, 2}},
  // Biobyte -1.378 Marvin -1.41
  SmilesExpected{"CS(=O)C", -1.082, {3, 70, 45, 3}},
  // Biobyte -1.498 Marvin -1.30
  SmilesExpected{"CS(=O)(=O)C", -0.428, {3, 71, 45, 45, 3}},
  // Biobyte -0.440 Marvin -0.37
  SmilesExpected{"CC=N", 0.072, {2, 21, 53}},
  // Biobyte 0.432 Marvin -0.15
  SmilesExpected{"CC=NC", -0.017, {2, 21, 53, 3}},
  // Biobyte 0.032 Marvin -0.44
  SmilesExpected{"CC=NNC", -0.016, {2, 21, 54, 50, 3}},
  //  Biobyte 0.632 Marvin 0.20
  SmilesExpected{"CN=NC", 0.556, {3, 55, 55, 3}},
  // Biobyte -0.048 Marvin 0.17 POOR MATCH
  SmilesExpected{"CN=NNC", -1.021, {3, 55, 56, 50, 3}},
  // Biobyte -0.518 MATCH -0.19
  SmilesExpected{"CNC", 0.000, {3, 49, 3}},
  // Biobyte -1.315 Marvin -0.09
  SmilesExpected{"OC(=O)NC(=O)O", -0.698, {39, 25, 44, 64, 25, 44, 39}},
  // Biobyte -0.235 Marvin -0.16
  SmilesExpected{"CCO", 0.017, {1, 6, 38}},
  // Biobyte -0.048 MATCH 0.12
  SmilesExpected{"COC", 0.311, {3, 41, 3}},
  // Biobyte -0.312 Marvin -0.82
  SmilesExpected{"CC(=O)N(C)C(=O)C", -0.255, {2, 24, 44, 65, 3, 24, 44, 2}},
  // Biobyte -0.312 Marvin 0.05
  SmilesExpected{"c1ncncc1", 0.056, {27, 57, 28, 57, 27, 26}},
  // Biobyte 2.68 Marvin 2.08
  SmilesExpected{"CC(C)C", 1.963, {1, 8, 1, 1}},
  // Biobyte 2.281 Marvin 1.80 
  SmilesExpected{"CCC", 1.695, {1, 4, 1}},
  // Biobyte 3.699 Marvin 3.38
  SmilesExpected{"CCCc1ccccc1", 3.111, {1, 4, 5, 29, 26, 26, 26, 26, 26}},
  // Biobyte -0.135 Marvin -0.27
  SmilesExpected{"NCC", -0.12, {46, 6, 1}},
  // Biobyte -1.565 Marvin -1.31
  SmilesExpected{"NCN", -1.327, {46, 7, 46}},
  // Biobyte 3.569 Marvin 3.22
  SmilesExpected{"CC(C)c1ccccc1", 3.089, {1, 9, 1, 29, 26, 26, 26, 26, 26}},
  // Biobyte 0.174 Marvin 0.15
  SmilesExpected{"CC(C)N", 0.383, {1, 10, 1, 46}},
  // Biobyte -1.256 Marvin -1.08
  SmilesExpected{"CC(N)N", -0.766, {1, 11, 46, 46}},
  // Biobyte -2.687 Marvin -1.51
  SmilesExpected{"NC(N)N", -1.878, {46, 11, 46, 46}},
  // Biobyte 3.079 Marvin 2.38
  SmilesExpected{"CC(C)(C)C", 2.11, {1, 12, 1, 1, 1}},
  // Biobyte 3.968 Marvin 3.52
  SmilesExpected{"CC(C)(C)N", 0.778, {1, 14, 1, 1, 46}},
  // Biobyte -0.947 Marvin -0.83   POOR MATCH
  SmilesExpected{"CC(C)(N)N", -0.132, {1, 15, 1, 46, 46}},
  //Biobyte -2.378 Marvin -1.26
  SmilesExpected{"CC(N)(N)N", -1.547, {1, 16, 46, 46, 46}},
  // Biobyte -4.794 Marvin -1.69
  SmilesExpected{"NC(N)(N)N", -2.322, {46, 17, 46, 46, 46}},
  // Biobyte -0.194 Marvin -0.22
  SmilesExpected{"CC(=O)O", -0.097, {2, 24, 44, 39}},
  // Biobyte -4.794 Marvin -1.69
  SmilesExpected{"NC(N)(N)N", -2.322, {46, 17, 46, 46, 46}},
  // Biobyte 1.427 Marvin 0.90
  SmilesExpected{"CC(F)(F)C", 1.834, {1, 15, 72, 72, 1}},
  // Biobyte 1.116 Marvin 1.46
  SmilesExpected{"CC(F)(F)F", 1.402, {1, 16, 72, 72, 72}},
  // Biobyte 1.025 Marvin 0.73
  SmilesExpected{"CCF", 0.863, {1, 6, 72}},
  // Figure 3 from paper
  // Biobyte 1.277 Marvin 1.17
  SmilesExpected{"NC(=O)c1ccccc1O", 0.444, {63, 24, 44, 29, 26, 26, 26, 26, 30, 39}}
));

// There should be a way of doing this without creating another class...
class TestXlogpPWithCorrections: public testing::TestWithParam<SmilesExpected> {
  protected:
    Molecule _m;
    std::unique_ptr<int[]> _atype;
};

TEST_P(TestXlogpPWithCorrections, WithCorrections) {
  const auto& params = GetParam();

  xlogp::ForTestingSetApplyCorrections(1);

  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  _atype.reset(new int[_m.natoms()]);
  std::optional<double> x = xlogp::XLogP(_m, _atype.get());
  ASSERT_NE(x, std::nullopt);
  EXPECT_NEAR(*x, params.expected, 0.001) << _m.smiles() << 
        " got " << *x << " expect " << params.expected << " C " << _m.name();
  std::span<int> values(_atype.get(), _m.natoms());
  EXPECT_THAT(values, ElementsAreArray(params.atypes)) << _m.smiles() << ToString(values) << ' ' << _m.name();
}
INSTANTIATE_TEST_SUITE_P(TestXlogpPWithCorrections, TestXlogpPWithCorrections,  testing::Values(
  // Biobyte 1.752 Marvin 1.35
  SmilesExpected{"CC", 1.624, {1,1}},
  // Biobyte 2.142 Marvin 1.97
  SmilesExpected{"c1ccccc1", 1.962, {26,26,26,26,26,26}},
  SmilesExpected{"c1ncccc1", 0.653, {27,57,27,26,26,26}},
  // Biobyte 1.116 Marvin 1.46
  SmilesExpected{"Fc1ccccc1", 2.064, {72, 30, 26, 26, 26, 26, 26}},
  SmilesExpected{"Clc1ccccc1", 2.581, {73, 30, 26, 26, 26, 26, 26}},
  SmilesExpected{"Brc1ccccc1", 2.758, {74, 30, 26, 26, 26, 26, 26}},
  SmilesExpected{"Ic1ccccc1", 3.06, {75, 30, 26, 26, 26, 26, 26}},
  // BioByte 1.885 Marvin 1.91
  SmilesExpected{"O=N(=O)c1ccccc1", 1.925, {80, 80, 80, 30, 26, 26, 26, 26, 26}},
  // BioByte 3.305 Marvin 3.00
  SmilesExpected{"S=C=Nc1ccccc1", 3.197, {78, 78, 78, 30, 26, 26, 26, 26, 26}},
  // BioByte 1.575 Marvin 1.83
  SmilesExpected{"N#Cc1ccccc1", 1.681, {77, 77, 29, 26, 26, 26, 26, 26}},
  // BioByte 1.885 Marvin 1.91
  SmilesExpected{"O=N(=O)c1ccccc1", 1.925, {80, 80, 80, 30, 26, 26, 26, 26, 26}},
  // Biobyte 2.025, Marvin 2.05
  SmilesExpected{"O=Nc1ccccc1", 1.648, {79, 79, 30, 26, 26, 26, 26, 26}},
  //Biobyte 3.305 Marvin 3.00
  SmilesExpected{"S=C=Nc1ccccc1", 3.197, {78, 78, 78, 30, 26, 26, 26, 26, 26}},
  // Biobyte -0.65 Marvin -0.81  POOR MATCH
  SmilesExpected{"N1NC(=O)C=C1", 0.055, {60, 64, 24, 44, 19, 20}},
  // Biobyte 2.641 Marvin 2.49
  SmilesExpected{"Cc1ccccc1", 2.433, {2, 29, 26, 26, 26, 26, 26}},
  // Biobyte -0.664 Marvin -0.63
  SmilesExpected{"CN methylamine", -0.533, {3, 46}},
  // Biobyte 0.915 Marvin 1.14
  SmilesExpected{"Nc1ccccc1", 1.214, {47, 30, 26, 26, 26, 26, 26}},
  // Biobyte -0.208 Marvin 0.11
  SmilesExpected{"O=C(C)C", 0.192, {44, 24, 2, 2}},
  // Biobyte 0.018 Marvin 0.19
  SmilesExpected{"CN(C)C", 0.314, {3, 51, 3, 3}},
  // Biobyte -1.350 Marvin -0.064
  SmilesExpected{"CN(C)NC", -0.581, {3, 52, 3, 50, 3}},
  // Biobyte -1.114 Marvin -1.03
  SmilesExpected{"CC(=O)N", -0.813, {2, 24, 44, 63}},
  // Biobyte -1.078 Marvin -0.81
  SmilesExpected{"CC(=O)NC", -0.323, {2, 24, 44, 64, 3}},
  // Biobyte -0.801 Marvin -0.58 POOR MATCH
  SmilesExpected{"CC(=O)N(C)C", -0.184, {2, 24, 44, 65, 3, 3}},
  // Biobyte 0.646 Marvin 0.84
  SmilesExpected{"CS", 0.755, {3, 66}},
  // Biobyte 1.788 Marvin 1.75
  SmilesExpected{"CSC", 0.985, {3, 67, 3}},
  // Biobyte 0.322 Marvin 1.00
  SmilesExpected{"S=C(C)C", 0.432, {69, 24, 2, 2}},
  // Biobyte -1.378 Marvin -1.41
  SmilesExpected{"CS(=O)C", -1.082, {3, 70, 45, 3}},
  // Biobyte -1.498 Marvin -1.30
  SmilesExpected{"CS(=O)(=O)C", -0.428, {3, 71, 45, 45, 3}},
  // Biobyte -0.440 Marvin -0.37
  SmilesExpected{"CC=N", 0.072, {2, 21, 53}},
  // Biobyte 0.432 Marvin -0.15
  SmilesExpected{"CC=NC", -0.017, {2, 21, 53, 3}},
  // Biobyte 0.032 Marvin -0.44
  SmilesExpected{"CC=NNC", -0.016, {2, 21, 54, 50, 3}},
  //  Biobyte 0.632 Marvin 0.20
  SmilesExpected{"CN=NC", 0.556, {3, 55, 55, 3}},
  // Biobyte -0.048 Marvin 0.17 POOR MATCH
  SmilesExpected{"CN=NNC", -1.021, {3, 55, 56, 50, 3}},
  // Biobyte -0.518 MATCH -0.19
  SmilesExpected{"CNC", 0.000, {3, 49, 3}},
  // Biobyte -1.315 Marvin -0.09
  SmilesExpected{"OC(=O)NC(=O)O", -0.698, {39, 25, 44, 64, 25, 44, 39}},
  // Biobyte -0.235 Marvin -0.16
  SmilesExpected{"CCO", 0.017, {1, 6, 38}},
  // Biobyte -0.048 MATCH 0.12
  SmilesExpected{"COC", 0.311, {3, 41, 3}},
  // Biobyte -0.312 Marvin -0.82
  SmilesExpected{"CC(=O)N(C)C(=O)C", -0.255, {2, 24, 44, 65, 3, 24, 44, 2}},
  // Biobyte -0.312 Marvin 0.05
  SmilesExpected{"c1ncncc1", 0.056, {27, 57, 28, 57, 27, 26}},
  // Biobyte 2.68 Marvin 2.08
  SmilesExpected{"CC(C)C", 2.723, {1, 8, 1, 1}},
  // Biobyte 2.281 Marvin 1.80 
  SmilesExpected{"CCC", 2.265, {1, 4, 1}},
  // Biobyte 3.699 Marvin 3.38
  SmilesExpected{"CCCc1ccccc1", 3.681, {1, 4, 5, 29, 26, 26, 26, 26, 26}},
  // Biobyte -0.135 Marvin -0.27
  SmilesExpected{"NCC", -0.120, {46, 6, 1}},
  // Biobyte -1.565 Marvin -1.31
  SmilesExpected{"NCN", -1.327, {46, 7, 46}},
  // Biobyte 3.569 Marvin 3.22
  SmilesExpected{"CC(C)c1ccccc1", 3.659, {1, 9, 1, 29, 26, 26, 26, 26, 26}},
  // Biobyte 0.174 Marvin 0.15
  SmilesExpected{"CC(C)N", 0.383, {1, 10, 1, 46}},
  // Biobyte -1.256 Marvin -1.08
  SmilesExpected{"CC(N)N", -0.766, {1, 11, 46, 46}},
  // Biobyte -2.687 Marvin -1.51
  SmilesExpected{"NC(N)N", -1.878, {46, 11, 46, 46}},
  // Biobyte 3.079 Marvin 2.38
  SmilesExpected{"CC(C)(C)C", 3.060, {1, 12, 1, 1, 1}},
  // Biobyte 3.968 Marvin 3.52  POOR MATCH
  SmilesExpected{"CC(C)(C)N", 0.778, {1, 14, 1, 1, 46}},
  // Biobyte -0.947 Marvin -0.83   POOR MATCH
  SmilesExpected{"CC(C)(N)N", -0.132, {1, 15, 1, 46, 46}},
  //Biobyte -2.378 Marvin -1.26
  SmilesExpected{"CC(N)(N)N", -1.547, {1, 16, 46, 46, 46}},
  // Biobyte -4.794 Marvin -1.69
  SmilesExpected{"NC(N)(N)N", -2.322, {46, 17, 46, 46, 46}},
  // Biobyte -0.194 Marvin -0.22
  SmilesExpected{"CC(=O)O acetic acid", -0.097, {2, 24, 44, 39}},
  // Biobyte -4.794 MATCH -1.69
  SmilesExpected{"NC(N)(N)N", -2.322, {46, 17, 46, 46, 46}},
  // Biobyte 1.427 Marvin 0.90  POOR MATCH
  SmilesExpected{"CC(F)(F)C", 1.914, {1, 15, 72, 72, 1}},
  // Biobyte 1.116 Marvin 1.46
  SmilesExpected{"CC(F)(F)F", 1.482, {1, 16, 72, 72, 72}},
  // Biobyte 1.025 Marvin 0.73
  SmilesExpected{"CCF", 0.863, {1, 6, 72}},
  // Figure 3 from paper
  // Biobyte 1.277 Marvin 1.17
  SmilesExpected{"NC(=O)c1ccccc1O", 1.044, {63, 24, 44, 29, 26, 26, 26, 26, 30, 39}},
  // Biobyte 1.268 Marvin 1.11
  SmilesExpected{"C=C", 1.354, {18, 18}},
  // Biobyte 0.923 Marvin 0.99
  SmilesExpected{"C#CC propyne", 0.943, {35, 36, 2}},
  // Biobyte1.7971.373 Marvin 1.49
  SmilesExpected{"C=CC propylene", 1.644, {18, 19, 2}},
  // Biobyte 2.236 Marvin 1.78
  SmilesExpected{"C1CC1 cyclopropane", 1.923, {4, 4, 4}},
  // Biobyte 2.281 Marvin 1.80
  SmilesExpected{"CCC propane", 2.265, {1, 4, 1}},
  // Biobyte 1.902 Marvin 1.63
  SmilesExpected{"C=CC=C 1,3-butadiene", 2.296, {18, 19, 19, 18}},
  // Biobyte 4.780 Marvin 4.09
  SmilesExpected{"C1=CC=C(C=C1)C#CC2=CC=CC=C2 diphenylacetylene", 4.568, {26, 26, 26, 29, 26, 26, 36, 36, 29, 26, 26, 26, 26, 26}},
  // Biobyte 2.412 Marvin 2.12
  SmilesExpected{"C#CC1=CC=CC=C1 ethynylbenzene", 2.384, {35, 36, 29, 26, 26, 26, 26, 26,}},
  // Biobyte 2.851 Marvin 2.70
  SmilesExpected{"C1C=CC2=CC=CC=C21 indene", 3.145, {5, 19, 19, 29, 26, 26, 26, 26, 29}},
  // Biobyte -0.764 Marvin -0.52
  SmilesExpected{"CO methanol", -0.396, {3, 38}},
  // Biobyte 1.974 Marvin 2.18
  SmilesExpected{"CC1=CC=C(C=C1)O p-cresol", 1.869, {2, 29, 26, 26, 30, 26, 26, 39}},
  // Biobyte 2.702 Marvin 2.13 Paper 2.03
  SmilesExpected{"C1=CC=C2C(=C1)C=CO2 benzofuran", 2.034, {26, 26, 26, 30, 29, 26, 19, 20, 43}},
  // Biobyte 3.172 Marvin 2.85
  SmilesExpected{"C1=CC=C2C(=C1)C=CS2 benzothiophene 2.67", 2.671, {26, 26, 26, 30, 29, 26, 19, 20, 68}},
  // Biobyte -0.005 Marvin 0.10
  SmilesExpected{"C(=O)(C(F)(F)F)N trifluoroacetanide -0.12", -0.259, {24, 44, 16, 72, 72, 72, 63}},

  // Biobyte 1.299 Marvin 1.57
  SmilesExpected{"CCC=O propionaldehyde", 0.438, {1, 5, 21, 44}},
  // Biobyte 0.655 Marvin 0.82
  SmilesExpected{"C1=CC=C(C=C1)C(=O)N benzamide", 0.818, {26, 26, 26, 29, 26, 26, 24, 44, 63}},
  // Biobyte 1.299 Marvin 1.57
  SmilesExpected{"CN1C2CCC1CC(C2)OC(=O)C(CO)C3=CC=CC=C3 atropine", 2.308, {3, 51, 10, 4, 4, 10, 4, 10, 4, 41, 24, 44, 9, 6, 38, 29, 26, 26, 26, 26, 26}},
  // Biobyte 1.486 Marvin 1.76
  SmilesExpected{"CC(C)NCC(COC1=CC=C(C=C1)CCOC)O metoprolol", 1.877, {1, 10, 1, 49, 6, 10, 6, 41, 30, 26, 26, 29, 26, 26, 5, 6, 41, 3, 38}},
  // Biobyte 2.652 Marvin 2.69
  SmilesExpected{"CC(C)NCC(COC1=CC=CC=C1CC=C)O alprenolol", 2.543, {1, 10, 1, 49, 6, 10, 6, 41, 30, 26, 26, 26, 26, 29, 5, 19, 18, 38}},
  // Biobyte 2.067 Marvin 1.81
  // Paper reports 2.61
  SmilesExpected{"C1=CC=C(C=C1)N2C(=CC=N2)NS(=O)(=O)C3=CC=C(C=C3)N sulfaphenazole", 2.783, {26, 26, 26, 30, 26, 26, 62, 23, 19, 21, 54, 64, 71, 45, 45, 30, 26, 26, 30, 26, 26, 47}},
  // Biobyte 2.305 Marvin 2.92
  SmilesExpected{"C1=CC=C(C=C1)C2=NC(C(=O)NC3=C2C=C(C=C3)Cl)O oxezepam", 2.870, {26, 26, 26, 29, 26, 26, 24, 53, 11, 24, 44, 64, 30, 29, 26, 30, 26, 26, 73, 38}},
  // Biobyte -0.109 Marvin 0.43
  SmilesExpected{"CC(C)NCC(COC1=CC=C(C=C1)CC(=O)N)O atenolol", 0.656, {1, 10, 1, 49, 6, 10, 6, 41, 30, 26, 26, 29, 26, 26, 5, 24, 44, 63, 38}},
  // Biobyte 4.466 Marvin 5.04
  // Paper reports 3.79. This seems correct. TODO:ianwatson
  SmilesExpected{"CC(C)C(CCCN(C)CCC1=CC(=C(C=C1)OC)OC)(C#N)C2=CC(=C(C=C2)OC)OC verapamil 3.79", 5.933, {1, 8, 1, 13, 4, 4, 6, 51, 3, 6, 5, 29, 26, 30, 30, 26, 26, 41, 3, 41, 3, 77, 77, 29, 26, 30, 30, 26, 26, 41, 3, 41, 3}},
  // Biobyte 0.981 Marvin 0.981
  // Paper reports 0.72. TODO:ianwatson
  SmilesExpected{"COC1=CC(=CC(=C1OC)OC)CC2=CN=C(N=C2N)N trimethoprim", 1.364, {3, 41, 30, 26, 29, 26, 30, 30, 41, 3, 41, 3, 5, 29, 27, 57, 33, 57, 32, 47, 47}},
  // Biobyte 3.371 Marvin 2.44
  SmilesExpected{"CCCN(CCC)S(=O)(=O)C1=CC=C(C=C1)C(=O)O probenecid", 2.459, {1, 4, 6, 65, 6, 4, 1, 71, 45, 45, 30, 26, 26, 29, 26, 26, 24, 44, 39}},
  // Biobyte 2.538 Marvin 1.88
  SmilesExpected{"CCN(CC)CCOC(=O)C1=CC=C(C=C1)N procaine", 1.845, {1, 6, 51, 6, 1, 6, 6, 41, 24, 44, 29, 26, 26, 30, 26, 26, 47}},
  // Biobyte 1.423 Marvin 0.95 Paper 1.27
  SmilesExpected{"CCN(CC)CCNC(=O)C1=CC=C(C=C1)N procainamide", 1.282, {1, 6, 51, 6, 1, 6, 6, 64, 24, 44, 29, 26, 26, 30, 26, 26, 47}},
  // Biobyte  2.861 Marvin 2.68 Paper 2.75
  SmilesExpected{"C1=CC=C2C(=C1)C=C3C=CC=C(C3=N2)N 4-aminoacridine", 2.755, {26, 26, 26, 34, 34, 26, 26, 34, 26, 26, 26, 30, 34, 57, 47}},
  // Biobyte 1.565 Marvin 1.26 Paper 1.38
  SmilesExpected{"C1=CC=C2C(=C1)NC=N2 benzimidazole 1.38", 1.380, {26, 26, 26, 30, 30, 26, 60, 21, 53}},
  // Biobyte 1.175 Marvin 1.22
  // Problem is aromaticity definitions. TODO:ianwatson
  SmilesExpected{"C1=CC2=C(NC=C2)N=C1 7-azaindole 1.28", 1.281, {26, 26, 29, 32, 60, 20, 19, 57, 27}},
  // Biobyte 0.288 Marvin 0.43
  SmilesExpected{"c1nccc1 pyrrole 0.80", 0.802, {20, 60, 20, 19, 19}},
  // Biobyte 1.318 Marvin 1.11
  SmilesExpected{"c1occc1 furan 0.65", 0.649, {20, 43, 20, 19, 19}},
  // Biobyte  4.95 Marvin 4.28 Paper 5.16
  SmilesExpected{"C1=CC2=C3C(=C1)C=CC4=CC=CC(=C43)C=C2 pyrene 5.16", 5.166, {26, 26, 34, 34, 34, 26, 26, 26, 34, 26, 26, 26, 34, 34, 26, 26}},
  // Biobyte -0.208 Marvin 0.11
  SmilesExpected{"CC(=O)C acetone 0.19", 0.192, {2, 24, 44, 2}},
  // Biobyte 0.185 Marvin 0.54
  SmilesExpected{"CCC#N propionitrile 0.47", 0.467, {1, 5, 77, 77}},
  // Biobyte 1.788 Marvin 1.75
  SmilesExpected{"c1sccc1 thiophene 1.29", 1.286, {20, 68, 20, 19, 19}},
  // Biobyte 1.283 Marvin 0.88 Paper 1.46
  // The assigments seem to be good, not sure why the discrepancy. TODO:ianwatson
  SmilesExpected{"C1=CC(=CC=C1[C@@H](O)[C@@H](CO)NC(=O)C(Cl)Cl)N(=O)=O chloramphenicol 1.46", 0.966, {26, 26, 30, 26, 26, 29, 10, 38, 10, 6, 38, 64, 24, 44, 11, 73, 73, 80, 80, 80}},
  // Biobyte -1.00 Marvin -0.44
  SmilesExpected{"C1=C2C(=CC(=C1Cl)S(=O)(=O)N)S(=O)(=O)N=CN2 chlorothiazide -0.58", 0.246, {26, 30, 30, 26, 30, 30, 73, 71, 45, 45, 63, 71, 45, 45, 54, 21, 58}},
  // Biobyte 5.50 Marvin 4.54
  SmilesExpected{"CN(C)CCCN1C2=CC=CC=C2SC2=C1C=C(Cl)C=C2 chlorpromazine 4.91", 5.022, {3, 51, 3, 6, 4, 6, 62, 30, 26, 26, 26, 26, 30, 67, 30, 30, 26, 30, 73, 26, 26}},
  // Biobyte 0.190 Marvin -0.11 Paper 0.20
  // Not sure what the discrepancy is. TODO:ianwatson
  SmilesExpected{"CC1=C(N=CN1)CSCCNC(=NC)NC#N cimetidine 0.20", 0.279, {2, 23, 23, 53, 21, 60, 6, 67, 6, 6, 49, 25, 53, 3, 58, 77, 77}},
  // Biobyte 2.961 Marvin 3.08 Paper 2.99
  // Likely problems with the N atom classifications TODO::ianwatson
  SmilesExpected{"CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C1=CC=CC=C1 diazepam 2.99", 2.988, {3, 65, 24, 44, 6, 53, 24, 29, 30, 26, 26, 30, 26, 73, 29, 26, 26, 26, 26, 26}},
  // Biobyte 3.647 Marvin 2.73 Paper 3.14
  SmilesExpected{"CC(=O)OC1C(SC2=CC=CC=C2N(C1=O)CCN(C)C)C1=CC=C(C=C1)OC diltiazem 3.14", 3.141, {2, 24, 44, 41, 10, 10, 67, 30, 26, 26, 26, 26, 30, 65, 24, 44, 6, 6, 51, 3, 3, 29, 26, 26, 30, 26, 26, 41, 3}},
  //Biobyte 2.89 Marvin 3.27  Paper 3.74
  // The discrepancy is hard to understand. The atom types seem quite straightforward
  // but there is a significant difference. TODO:ianwatson
  SmilesExpected{"N(C)CCOC(C1=CC=CC=C1)C1=CC=CC=C1 diphenhydramine 3.74", 3.439, {49, 3, 6, 6, 41, 10, 29, 26, 26, 26, 26, 26, 29, 26, 26, 26, 26, 26}},
  // Biobyte 2.875 Marvin 3.00
  // This is problematic. Without the halogen correction the sum is 3.678 and
  // if we subtract 3*-0.26 we get 2.89. But other times the halogen correction
  // does not get applied that way.
  SmilesExpected{"ClC(Cl)(Cl)Cl CCl4 2.89", 3.418, {73, 17, 73, 73, 73}},
  // Biobyte 1.952 Marvin 1.83
  // We get close to 1.83 if we take the raw sum 2.622- 3*-0.26.
  SmilesExpected{"C(Cl)(Cl)Cl Chloroform 1.83", 2.362, {11, 73, 73, 73}},
  //  0.368 Marvin 0.42
  SmilesExpected{"FCF difluoromethane 0.71", 0.719, {72, 7, 72}},
  // Biobyte  1.249 1.29
  SmilesExpected{"ClCCl dichloromethane 1.413", 1.413, {73, 7, 73}},
  // Biobyte 0.496 Marvin 0.37
  SmilesExpected{"CF methylfluoride 0.45", 0.450, {3, 72}},
  // Biobyte 3.48 Marvin 2.62
  SmilesExpected{"ClC(Cl)=C(Cl)Cl tetrachloroethylene 3.30", 3.304, {73, 23, 73, 23, 73, 73}},
  // Biobyte 2.85 Marvin 2.50
  // This is very mysterious. This molecule is very simple, but I cannot
  // discern how the paper gets 3.04. Without the halogen correction the sum is 2.594.
  // That means the sum increases by 0.44 as a result of the halogen correction.
  // Cannot figure out what they did. TODO:ianwatson
  SmilesExpected{"C(C(F)(F)Cl)(F)(F)Cl 1,2-dichlorotctrafluoroethane 3.04", 2.754, {16, 16, 72, 72, 73, 72, 72, 73}},
  // Biobyte 1.735 Marvin 2.08
  // This is hard to understand. Without the halogen correction the sum is 1.56 and so
  // again the halogen correction is 0.44. This time there are just F atoms, before there
  // were F and Cl. Very hard to understand.
  SmilesExpected{"FC(F)(F)C(F)(F)F hexafluoroethane 2.01", 1.72, {72, 16, 72, 72, 16, 72, 72, 72}},
  // Biobyte 0.898 Marvin 0.65
  SmilesExpected{"CC(F)F 1,1-difluoroethane 1.27", 1.28, {1, 11, 72, 72}},
  // Biobyte 1.778 Marvin 1.52
  // Hard to understand how they get 2.16. The sum before correctios is 2.234
  // and if we subtract 0.26 get the answer below. TODO:ianwatson
  SmilesExpected{"CC(Cl)Cl 1,1-diChloroethane 2.16", 1.974, {1, 11, 73, 73}},
  // Biobyte 3.568 Marvin 3.18
  SmilesExpected{"Clc1ccc(Cl)cc1 p-dichlorobenzene 3.20", 3.20, {73, 30, 26, 26, 30, 73, 26, 26}},
  // Biobyte 3.849 Marvin 3.66
  // Very significant differences here. Very hard to understand
  // why. TODO:ianwatson
  SmilesExpected{"C1CN(CCC1(O)C1=CC=C(Cl)C=C1)CCCC(=O)C1=CC=C(F)C=C1 haloperidol 3.57", 4.365, {4, 6, 51, 6, 4, 14, 38, 29, 26, 26, 30, 73, 26, 26, 6, 4, 5, 24, 44, 29, 26, 26, 30, 72, 26, 26}}
));


}  // namespace
