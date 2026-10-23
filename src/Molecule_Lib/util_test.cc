#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/set_of_atoms.h"
#include "Molecule_Lib/util.h"

namespace {

using testing::ElementsAre;

Molecule
MolFromSmiles(const char* smiles) {
  Molecule result;
  EXPECT_TRUE(result.build_from_smiles(smiles)) << smiles;
  return result;
}

TEST(SetAtomsWithinRadius, RadiusZeroSetsOnlySeeds) {
  Molecule m = MolFromSmiles("CCCC");
  Set_of_Atoms seeds;
  seeds << 1;

  std::vector<int> destination(m.natoms(), 0);
  EXPECT_EQ(lillymol::SetAtomsWithinRadius(m, seeds, 0, 7, destination.data()), 1);

  EXPECT_THAT(destination, ElementsAre(0, 7, 0, 0));
}

TEST(SetAtomsWithinRadius, RadiusOneFromSingleSeed) {
  Molecule m = MolFromSmiles("CCCC");
  Set_of_Atoms seeds;
  seeds << 1;

  std::vector<int> destination(m.natoms(), 0);
  EXPECT_EQ(lillymol::SetAtomsWithinRadius(m, seeds, 1, 3, destination.data()), 3);

  EXPECT_THAT(destination, ElementsAre(3, 3, 3, 0));
}

TEST(SetAtomsWithinRadius, RadiusTwoFromMultipleSeeds) {
  Molecule m = MolFromSmiles("CCCCCC");
  Set_of_Atoms seeds;
  seeds << 1;
  seeds << 5;

  std::vector<int> destination(m.natoms(), 0);
  EXPECT_EQ(lillymol::SetAtomsWithinRadius(m, seeds, 2, 1, destination.data()), 6);

  EXPECT_THAT(destination, ElementsAre(1, 1, 1, 1, 1, 1));
}

TEST(SetAtomsWithinRadius, BranchesAreReachedByBondDistance) {
  Molecule m = MolFromSmiles("CC(C)CC");
  Set_of_Atoms seeds;
  seeds << 1;

  std::vector<int> destination(m.natoms(), 0);
  EXPECT_EQ(lillymol::SetAtomsWithinRadius(m, seeds, 1, 1, destination.data()), 4);

  EXPECT_THAT(destination, ElementsAre(1, 1, 1, 1, 0));
}

TEST(SetAtomsWithinRadius, CountsOnlyChangedEntries) {
  Molecule m = MolFromSmiles("CCCC");
  Set_of_Atoms seeds;
  seeds << 1;

  std::vector<int> destination = {0, 5, 0, 0};
  EXPECT_EQ(lillymol::SetAtomsWithinRadius(m, seeds, 1, 5, destination.data()), 2);

  EXPECT_THAT(destination, ElementsAre(5, 5, 5, 0));
}

TEST(SetAtomsWithinRadius, NegativeRadiusIsNoOp) {
  Molecule m = MolFromSmiles("CCCC");
  Set_of_Atoms seeds;
  seeds << 1;

  std::vector<int> destination(m.natoms(), 0);
  EXPECT_EQ(lillymol::SetAtomsWithinRadius(m, seeds, -1, 1, destination.data()), 0);

  EXPECT_THAT(destination, ElementsAre(0, 0, 0, 0));
}

}  // namespace
