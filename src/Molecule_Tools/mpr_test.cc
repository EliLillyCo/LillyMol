// Tests for mpr

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Foundational/cmdline/cmdline.h"
#include "Molecule_Tools/mpr.h"

namespace {

using testing::ElementsAreArray;

struct PropertyData {
  std::vector<const char*> argv;
  const char* options;
  IWString smiles;
  int nfeatures;
  std::vector<int> result;
};

class TestMPR: public testing::TestWithParam<PropertyData> {
  protected:
    Molecule _m;
    mpr::MolecularPropertiesGenerator _mpg;
};

TEST_P(TestMPR, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles)) << "invalid smiles " << params.smiles;

  const int argc = params.argv.size();
  ASSERT_GT(argc, 0);
  char ** argv = new char*[argc];
  std::unique_ptr<char *[]> free_argv(argv);
  int ndx = 0;
  for (const auto& c : params.argv) {
    argv[ndx] = const_cast<char*>(c);
    ++ndx;
  }
  Command_Line cl(argc, argv, params.options);

  ASSERT_TRUE(_mpg.Initialise(cl, 'O'));
  int nfeatures = _mpg.number_features();
  ASSERT_EQ(nfeatures, params.nfeatures);

  std::unique_ptr<int[]> res = std::make_unique<int[]>(nfeatures);
  ASSERT_TRUE(_mpg.GenerateMolecularProperties(_m, res.get()));

  // Copy _res to a vector so we can use ElementsAreArray.
  std::vector<int> tmp(res.get(), res.get() + nfeatures);
  EXPECT_THAT(tmp, ElementsAreArray(params.result)) << "Mismatch " << _m.smiles();
}
// We deliberately do not test '-O none' by itself. Because _mpg will
// use a default set of features if none are specified and a calculation is requested.
INSTANTIATE_TEST_SUITE_P(TestMPR, TestMPR, testing::Values(
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "natoms"}, "O:", "C", 1, {1} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "lrsz"}, "O:", "C1CC1", 1, {3} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "nrings"}, "O:", "C1CC1C", 1, {1} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "ringatom"}, "O:", "C1CCC1", 1, {4} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "aroma"}, "O:", "C1CCC1", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "aroma"}, "O:", "c1ccccc1O", 1, {6} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "frac"}, "O:", "C1CC1", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "frac"}, "O:", "C1CC12CC2", 1, {1} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "frac"}, "O:", "C12CCC1CC2", 1, {2} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "htroatom"}, "O:", "C", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "htroatom"}, "O:", "CN", 1, {1} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "unsatura"}, "O:", "CN", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "unsatura"}, "O:", "C=C", 1, {2} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "unsatura"}, "O:", "c1ccccc1", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "mxdst"}, "O:", "C", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "mxdst"}, "O:", "c1ccccc1", 1, {3} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "scaffold"}, "O:", "C", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "scaffold"}, "O:", "c1ccccc1", 1, {6} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "scaffold"}, "O:", "Cc1cc(C)ccc1", 1, {6} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "lrsysz"}, "O:", "C", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "lrsysz"}, "O:", "C1CC1", 1, {1} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "lrsysz"}, "O:", "C12CCC2CC1", 1, {2} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "linkatom"}, "O:", "C", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "linkatom"}, "O:", "CC", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "linkatom"}, "O:", "C1CC1C1CC1", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "linkatom"}, "O:", "C1CC1CC1CC1", 1, {1} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "rotbond"}, "O:", "CC", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "rotbond"}, "O:", "CCC", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "rotbond"}, "O:", "CCC", 1, {0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "rotbond"}, "O:", "CCCC", 1, {1} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "rotbond"}, "O:", "CCCC(F)(F)F", 1, {1} },

  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "natoms", "-O", "nrings"},
                "O:", "C", 2, {1, 0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "all"},
                "O:", "C", 13, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} },
  PropertyData{std::vector<const char*>{"-", "-O", "none", "-O", "all"},
                "O:", "CC(CNC(=O)C1=CC=C(Cl)C=C1I)N1CCCC1 CHEMBL3469286", 13,
                { 19, 6, 2, 11, 6, 0, 5, 2, 11, 16, 1, 5, 4 } }
));

}  // namespace
