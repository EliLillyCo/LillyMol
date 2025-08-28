// Tests for ring finding.

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "molecule.h"
#include "path.h"

namespace {

// This is a surprising molecule. Even though it is a cage-like structure
// none of the component rings are strongly fused - all rings share just one
// bond with another ring.
TEST(TestFused, TestFused) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1C2C3C1N1C2C31"));
  EXPECT_EQ(m.nrings(), 4);

  for (int i = 0; i < 4; ++i) {
    const Ring* ri = m.ringi(i);
    EXPECT_TRUE(ri->is_fused());
  }

  EXPECT_EQ(m.non_sssr_rings(), 0);

  for (int i = 0; i < 4; ++i) {
    const Ring* ri = m.ringi(i);
    EXPECT_EQ(ri->largest_number_of_bonds_shared_with_another_ring(), 1);
    EXPECT_EQ(ri->strongly_fused_ring_neighbours(), 0);
  }

}

}  // namespace
