// Tests for rotate_atoms

#include <memory>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "molecule.h"

namespace {

TEST(TestRotateAtoms, Test1) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("N{{-0.0178,1.4648,0.0101}}C{{0.0021,-0.0041,0.002}}C{{0.8136,-0.5019,1.1702}}(=O{{1.3252,0.2894,1.9346}})S{{1.0136,-2.1839,1.4188}}"));

  std::unique_ptr<int[]> fragment_membership = m.fragment_membership();

  int moving_frag = fragment_membership[0];

  std::vector<float> initial_distances;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      initial_distances.push_back(m.distance_between_atoms(i, j));
    }
  }

  std::vector<float> initial_angles;
  for (int i = 0; i < matoms; ++i) {
    for (int j = 0; j < matoms; ++j) {
      if (i == j) {
        continue;
      }
      for (int k = 0; k < matoms; ++k) {
        if (k == i || k == j) {
          continue;
        }
        initial_angles.push_back(m.bond_angle(i, j, k));
      }
    }
  }

  m.translate_atoms(-m[0], fragment_membership.get(), moving_frag);
  EXPECT_FLOAT_EQ(m.x(0), 0.0f);
  EXPECT_FLOAT_EQ(m.y(0), 0.0f);
  EXPECT_FLOAT_EQ(m.z(0), 0.0f);

  const float theta = -2.42896;
  Space_Vector<float> axis(-0.576296,-0.0579064,-0.492563);
  axis.normalise();
  m.rotate_atoms(axis, theta, fragment_membership.get(), moving_frag);

  int ndx = 0;
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      float d = m.distance_between_atoms(i, j);
      EXPECT_FLOAT_EQ(d, initial_distances[ndx]) << i << ' ' << j;
      ++ndx;
    }
  }

  ndx = 0;
  for (int i = 0; i < matoms; ++i) {
    for (int j = 0; j < matoms; ++j) {
      if (i == j) {
        continue;
      }
      for (int k = 0; k < matoms; ++k) {
        if (k == i || k == j) {
          continue;
        }
        float theta = m.bond_angle(i, j, k);
        EXPECT_NEAR(theta, initial_angles[ndx], 0.0001);
        ++ndx;
      }
    }
  }
}

}  // namespace
