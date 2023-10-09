// Tests for the dihedral angle functions

#include <array>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "molecule.h"

namespace {

using testing::FloatEq;

struct CoordsAngle {
  // Coordinates for 4 atoms.
  std::array<float, 3> coords[4];
  // The expected dihedral angle
  float expected_angle;
  // A rotation axis. The dihedral should be invariant to
  // rotations about this axis.
  std::array<float, 3> rot;
};

class TestDihedral : public testing::TestWithParam<CoordsAngle> {
  protected:
    Molecule mol;
};

TEST_P(TestDihedral, TestAngles) {
  mol.build_from_smiles("CCCC");
  const auto params = GetParam();

  constexpr double kAbsDiff = 0.0001;

  for (int i = 0; i < 4; ++i) {
    mol.setxyz(i, params.coords[i][0], params.coords[i][1], params.coords[i][2]);
  }

  const float angle = mol.dihedral_angle(0, 1, 2, 3) * RAD2DEG;
  // std::cerr << " Initial config angle " << angle << '\n';
  EXPECT_NEAR(angle, params.expected_angle, kAbsDiff);

  Coordinates rot(params.rot[0], params.rot[1], params.rot[2]);
  rot.normalise();

  // We will do a series of rotations by 30 degrees.
  constexpr double k30 = 30.0 * RAD2DEG;
  for (int i = 1; i < 12; ++i) {
    // Reset coordinates each time to avoid float inaccuracies.
    for (int j = 0; j < 4; ++j) {
      mol.setxyz(j, params.coords[j][0], params.coords[j][1], params.coords[j][2]);
    }
    mol.rotate_atoms(rot, i * k30);
    const float angle = mol.dihedral_angle(0, 1, 2, 3) * RAD2DEG;
    // std::cerr << " i " << i << " angle " << angle << '\n';
    EXPECT_NEAR(angle, params.expected_angle, kAbsDiff);
  }
}
INSTANTIATE_TEST_SUITE_P(TestDihedral, TestDihedral, testing::Values(
  CoordsAngle{{{-1.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, 0.0, {1,2,3}},
  CoordsAngle{{{-1.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, 180.0, {-1,2,3}},
  CoordsAngle{{{-1.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 1.0}}, 90.0, {0,2,-1}},
  CoordsAngle{{{-1.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, -1.0}}, 90.0, {-3,-2,-1}}
));

class TestSignedDihedral : public testing::TestWithParam<CoordsAngle> {
  protected:
    Molecule mol;
};

TEST_P(TestSignedDihedral, TestAngles) {
  mol.build_from_smiles("CCCC");
  const auto params = GetParam();

  constexpr double kAbsDiff = 0.0001;

  for (int i = 0; i < 4; ++i) {
    mol.setxyz(i, params.coords[i][0], params.coords[i][1], params.coords[i][2]);
  }

  const float angle = mol.signed_dihedral_angle(0, 1, 2, 3) * RAD2DEG;
  // std::cerr << " Signed: Initial config angle " << angle << '\n';
  EXPECT_NEAR(angle, params.expected_angle, kAbsDiff);

  Coordinates rot(params.rot[0], params.rot[1], params.rot[2]);
  rot.normalise();

  // We will do a series of rotations by 30 degrees.
  constexpr double k30 = 30.0 * RAD2DEG;
  for (int i = 1; i < 12; ++i) {
    // Reset coordinates each time to avoid float inaccuracies.
    for (int j = 0; j < 4; ++j) {
      mol.setxyz(j, params.coords[j][0], params.coords[j][1], params.coords[j][2]);
    }
    mol.rotate_atoms(rot, i * k30);
    const float angle = mol.signed_dihedral_angle(0, 1, 2, 3) * RAD2DEG;
    // std::cerr << " i " << i << " angle " << angle << '\n';
    EXPECT_NEAR(angle, params.expected_angle, kAbsDiff);
  }
}
INSTANTIATE_TEST_SUITE_P(TestSignedDihedral, TestSignedDihedral, testing::Values(
  CoordsAngle{{{-1.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, 0.0, {1,2,3}},
  CoordsAngle{{{-1.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, 180.0, {-1,2,3}},
  CoordsAngle{{{-1.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 1.0}}, 90.0, {0,2,-1}},
  CoordsAngle{{{-1.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, -1.0}}, -90.0, {-3,-2,-1}},

  CoordsAngle{{{-1.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.01}}, -0.57293868, {1,2,3}},
  CoordsAngle{{{-1.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, -0.01}}, 0.57293868, {1,2,3}}
));
}
