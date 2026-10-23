// Tests for Molecule::discern_chirality_from_wedge_bonds() (mdl.cc).
//
// _discern_chirality_from_wedge_bond_4() used to reduce a stereocentre's
// neighbours to a single 2D "rotation" sign without ever consulting which
// neighbour actually carried the wedge. That is only valid when the
// neighbours span more than 180 degrees around the centre; it silently gives
// the same (sometimes wrong) answer whenever two of them are drawn
// (anti)collinear through the centre - a routine occurrence at ring
// fusions/bridgeheads and in ordinary zig-zag chains, not a rare pathology.

#include <fstream>
#include <string>

#include "gtest/gtest.h"

#include "Foundational/data_source/iwstring_data_source.h"

#include "molecule.h"

namespace {

Molecule
ReadSdf(const char* sdf, const char* stem) {
  const std::string fname = testing::TempDir() + stem;
  std::ofstream output(fname);
  output << sdf;
  output.close();

  iwstring_data_source input(fname.c_str());
  EXPECT_TRUE(input.good());

  Molecule m;
  EXPECT_TRUE(m.read_molecule_ds(input, FILE_TYPE_SDF));
  return m;
}

// The minimal reproducer from the original bug report: a single sp3 carbon
// with three explicit substituents fanned out below it, differing only in
// which bond carries the hash. The two molfiles must be enantiomers.
TEST(WedgeBondChirality, WedgeAtomIdentityDeterminesHandedness) {
  static constexpr char kHashToF[] = R"(minimal-hash-to-F
  hand written                    2D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.5000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990   -0.7500    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  6  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END
$$$$
)";

  static constexpr char kHashToCl[] = R"(minimal-hash-to-Cl
  hand written                    2D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.5000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990   -0.7500    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  6  0  0  0
  1  4  1  0  0  0  0
M  END
$$$$
)";

  Molecule m1 = ReadSdf(kHashToF, "/hash_to_f.sdf");
  Molecule m2 = ReadSdf(kHashToCl, "/hash_to_cl.sdf");

  EXPECT_NE(m1.unique_smiles(), m2.unique_smiles())
      << "moving the hash to a different bond must change the stereoisomer";
}

// Two of a 3-connected stereocentre's substituents drawn directly opposite
// each other through the centre - here because they are the xy-projection
// of a genuine tetrahedral arrangement, the configuration an ordinary
// zig-zag chain regularly produces. Flipping the wedge direction while
// holding every coordinate fixed must flip the enantiomer.
TEST(WedgeBondChirality, CollinearInPlaneNeighboursWedgeUp) {
  static constexpr char kWedgeUp[] = R"(degenerate-wedge-up
  hand written                    2D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5130    1.4095    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END
$$$$
)";

  static constexpr char kWedgeDown[] = R"(degenerate-wedge-down
  hand written                    2D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5130    1.4095    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  6  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END
$$$$
)";

  Molecule up = ReadSdf(kWedgeUp, "/degenerate_up.sdf");
  Molecule down = ReadSdf(kWedgeDown, "/degenerate_down.sdf");

  EXPECT_NE(up.unique_smiles(), down.unique_smiles())
      << "Cl and Br are drawn collinear through the stereocentre; only the "
         "wedge on F distinguishes the two enantiomers";

  // Cross check against ground truth that does not go anywhere near 2D
  // wedge interpretation: the same molecule, genuinely 3D, with F lifted
  // above the page and an explicit H completing a real tetrahedron, read
  // with discern_chirality_from_3d_structure(). F's wedge-up 2D depiction
  // above is exactly the xy-projection of this 3D structure.
  static constexpr char k3D[] = R"(tetrahedral-3d
  hand written                    3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5130    1.4095    1.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000   -1.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000   -1.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.5130   -1.4095    1.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
)";
  Molecule m3d = ReadSdf(k3D, "/degenerate_3d.sdf");
  ASSERT_TRUE(m3d.discern_chirality_from_3d_structure());
  m3d.remove_explicit_hydrogens();

  EXPECT_EQ(up.unique_smiles(), m3d.unique_smiles())
      << "the 2D wedge-up depiction is a projection of this 3D structure, so "
         "they must resolve to the same stereoisomer";
}

}  // namespace
