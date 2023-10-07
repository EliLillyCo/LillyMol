// Tests for Cahn Ingold Prelog chirality

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "molecule.h"

namespace {
struct RS {
  IWString smiles;
  // Do we expect a valid result or not.
  // For example, if we test an atom that does not
  // have a chiral centre, we do not get a result.
  bool expect_result;
  atom_number_t zatom;
  // The expected value
  CahnIngoldPrelog rs;
};

class TestCip: public testing::TestWithParam<RS> {
  protected:
    Molecule _m;
};

TEST_P(TestCip, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles)) << "Bad smiles " << params.smiles;
  std::optional<CahnIngoldPrelog> cip = _m.CahnIngoldPrelogValue(params.zatom);
  if (! params.expect_result) {
    EXPECT_EQ(cip, std::nullopt);
    return;
  }

  EXPECT_THAT(_m.CahnIngoldPrelogValue(params.zatom), testing::Optional(params.rs)) <<
    "Mismatch in " << params.smiles << " atom " << params.zatom << " expected " << params.rs;
}
INSTANTIATE_TEST_SUITE_P(TestCip, TestCip, testing::Values(
  RS{"C", false, 0, CahnIngoldPrelog::kUnspecified},

  RS{"I[C@H](Br)F", true, 1, CahnIngoldPrelog::S},
  RS{"Br[C@H](F)I", true, 1, CahnIngoldPrelog::S},
  RS{"I[C@H](Br)F", true, 1, CahnIngoldPrelog::S},
  RS{"Br[C@@H](I)F", true, 1, CahnIngoldPrelog::S},
  RS{"F[C@H](I)Br", true, 1, CahnIngoldPrelog::S},
  RS{"I[C@@H](F)Br", true, 1, CahnIngoldPrelog::S},
  RS{"F[C@@H](Br)I", true, 1, CahnIngoldPrelog::S},

  RS{"F[C@@H](I)Br", true, 1, CahnIngoldPrelog::R},
  RS{"Br[C@H](I)F", true, 1, CahnIngoldPrelog::R},
  RS{"Br[C@@H](F)I", true, 1, CahnIngoldPrelog::R},
  RS{"F[C@H](Br)I", true, 1, CahnIngoldPrelog::R},
  RS{"I[C@@H](Br)F", true, 1, CahnIngoldPrelog::R},
  RS{"I[C@H](F)Br", true, 1, CahnIngoldPrelog::R}
));

}  // namespace
