// Tests for ring finding
#include <algorithm>
#include <array>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "molecule.h"
#include "path.h"

namespace {

TEST(TestRing, Unfused) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1CC1"));
  EXPECT_EQ(m.nrings(), 1);
  const Ring* r = m.ringi(0);
  EXPECT_EQ(r->size(), 3);
  EXPECT_FALSE(r->is_fused());
  // Ring numbers are not guaranteed.
  EXPECT_EQ(r->ring_number(), 0);
  EXPECT_EQ(r->fragment_membership(), 0);
  EXPECT_EQ(r->fused_ring_neighbours(), 0);
}

TEST(TestRing, Fused1) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1C2CC12"));
  EXPECT_EQ(m.nrings(), 2);
  for (const Ring* r : m.sssr_rings()) {
    EXPECT_EQ(r->size(), 3);
    EXPECT_TRUE(r->is_fused());
    EXPECT_EQ(r->fused_ring_neighbours(), 1);
  }

  std::array<int, 4> visited;
  std::fill_n(visited.begin(), 4, 0);
  for (const Ring* r : m.sssr_rings()) {
    int nbrs = r->fused_ring_neighbours();
    for (int i = 0; i < nbrs; ++i) {
      const Ring* nbr = r->fused_neighbour(i);
      nbr->increment_vector(visited.data(), 1);
    }
  }
  EXPECT_THAT(visited, testing::ElementsAre(1, 2, 1, 2));

  // Make sure that a molecule derived from a copy constructor also
  // works.
  Molecule m2(m);
  EXPECT_EQ(m2.nrings(), m.nrings());
  std::fill_n(visited.begin(), 4, 0);
  for (const Ring* r : m2.sssr_rings()) {
    int nbrs = r->fused_ring_neighbours();
    for (int i = 0; i < nbrs; ++i) {
      const Ring* nbr = r->fused_neighbour(i);
      nbr->increment_vector(visited.data(), 1);
    }
  }
  EXPECT_THAT(visited, testing::ElementsAre(1, 2, 1, 2));

  // Make sure that copied Molecule's also work.
  Molecule m3;
  m3 = m;
  std::fill_n(visited.begin(), 4, 0);
  for (const Ring* r : m3.sssr_rings()) {
    int nbrs = r->fused_ring_neighbours();
    for (int i = 0; i < nbrs; ++i) {
      const Ring* nbr = r->fused_neighbour(i);
      nbr->increment_vector(visited.data(), 1);
    }
  }
  EXPECT_THAT(visited, testing::ElementsAre(1, 2, 1, 2));
}

}  // namespace
