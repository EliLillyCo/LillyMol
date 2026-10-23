// Tests for nvrtspsa


#include <optional>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Tools/nvrtspsa.h"

namespace {

struct Data {
  IWString smiles;
  double expected;
};

class NvrtsPSATest : public testing::TestWithParam<Data> {
  protected:
    Molecule _m;

};

TEST_P(NvrtsPSATest, Tests) {
  const auto& params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  nvrtspsa::NovartisPolarSurfaceArea tpsa;
  tpsa.set_convert_to_charge_separated(1);
  tpsa.set_zero_for_all_phosphorus_atoms(1);
  tpsa.set_zero_for_all_sulphur_atoms(1);

  std::optional<double> result = tpsa.PolarSurfaceArea(_m);
  ASSERT_TRUE(result) << _m.smiles();
  EXPECT_NEAR(*result, params.expected, 0.001) << _m.smiles();
}
INSTANTIATE_TEST_SUITE_P(NvrtsPSATest, NvrtsPSATest, testing::Values(
  // This value is quite different from RdKit due to aromaticity.
  Data{"O=C1N(C2=C3C(=CC=CC3=CC=C2)N1)CC CHEMBL1397386", 37.79},
  Data{"C1(=CC=CN1)N(=O)=O CHEMBL1651441", 58.93},
  // This is different from RDKit, not sure why
  Data{"C12=C(C3=NCCN3C(=O)N1)N(C)C=N2 CHEMBL326789", 67.97},
  // Different from RDKit, maybe aromaticity.
  Data{"C12=NCCN1C1=C(C3=CC=CC=C23)C=CC=C1 CHEMBL1181670", 17.29}
));

}  // namespace
