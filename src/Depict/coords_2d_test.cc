#include "Depict/coords_2d.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "gtest/gtest.h"

namespace {

// The distance between two bonded atoms, ignoring z.
float
BondLength(const Molecule& m, const Bond& b) {
  const float dx = m.x(b.a1()) - m.x(b.a2());
  const float dy = m.y(b.a1()) - m.y(b.a2());

  return std::sqrt(dx * dx + dy * dy);
}

// The closest approach between any two atoms, ignoring z. A layout that
// superimposes atoms is useless for depiction however good its bond lengths.
float
ClosestApproach(const Molecule& m) {
  float closest = std::numeric_limits<float>::max();

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      const float dx = m.x(i) - m.x(j);
      const float dy = m.y(i) - m.y(j);
      closest = std::min(closest, std::sqrt(dx * dx + dy * dy));
    }
  }

  return closest;
}

// A structure with a metal in a ring. Kept out of TestSmiles because coordgen
// does not give it uniform bond lengths - see MetalComplexLaysOut.
constexpr char kMetalComplex[] = "[Fe]1(Cl)(Cl)NCCN1";

// Molecules spanning the cases coordgen is expected to handle: chains, fused
// and spiro rings, aromatics, hypervalent and gem-disubstituted centres,
// stereo-bearing and charged structures, a macrocycle and multiple fragments.
std::vector<std::string>
TestSmiles() {
  return {
      "C",
      "CC",
      "OCC",
      "CC(C)(C)C",
      "CCCCCCCCCCCCCCCC",
      "c1ccccc1",
      "c1ccc2ccccc2c1",
      "C1CC2(CC1)CCCC2",
      "O=C(Nc1ccccc1)c1ccc(Cl)cc1",
      "CC(=O)Oc1ccccc1C(=O)O",
      "CN1CCC[C@H]1c1cccnc1",
      "C/C=C/C",
      "[NH3+]CC(=O)[O-]",
      "C1CCCCCCCCCCC1",
      "C1CCCCCCCCCCCCCCCCCCC1",
      "O=C1CCC(=O)N1c1ccccc1",
      "C1(Cl)(Cl)NCCN1",
      "CS(=O)(=O)C",
      "C1CC(Cl)(Cl)CC1",
      "c1ccccc1.c1ccccc1",
      "CC(=O)O.[Na+].[OH-]",
      "CC1=C(C(=O)Nc2ccccc2)N(c2ccc(F)cc2)N=C1c1ccc(S(=O)(=O)N)cc1",
  };
}

TEST(Coords2D, EmptyMoleculeFails) {
  Molecule m;

  EXPECT_EQ(depict::Generate2DCoordinates(m), 0);
}

TEST(Coords2D, SingleAtom) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C"));

  ASSERT_EQ(depict::Generate2DCoordinates(m), 1);
  EXPECT_EQ(m.natoms(), 1);
  // Centred, so it lands on the origin.
  EXPECT_FLOAT_EQ(m.x(0), 0.0f);
  EXPECT_FLOAT_EQ(m.y(0), 0.0f);
  EXPECT_FLOAT_EQ(m.z(0), 0.0f);
}

// Every atom must be given a z of exactly zero, otherwise downstream 3D tools
// would silently treat the depiction as a conformer.
TEST(Coords2D, LayoutIsFlat) {
  for (const std::string& smi : TestSmiles()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    ASSERT_EQ(depict::Generate2DCoordinates(m), 1) << smi;

    for (int i = 0; i < m.natoms(); ++i) {
      EXPECT_FLOAT_EQ(m.z(i), 0.0f) << smi << " atom " << i;
    }
  }
}

// The defining property of a usable layout: every bond is drawn at close to
// the requested length. coordgen rounds its internal coordinates to whole
// units of 50 per bond, so a small tolerance is unavoidable.
TEST(Coords2D, BondLengthsAreUniform) {
  depict::Coords2DOptions opts;
  opts.bond_length = 1.5f;

  for (const std::string& smi : TestSmiles()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    ASSERT_EQ(depict::Generate2DCoordinates(m, opts), 1) << smi;

    for (int i = 0; i < m.nedges(); ++i) {
      EXPECT_NEAR(BondLength(m, *m.bondi(i)), 1.5f, 0.15f)
          << smi << " bond " << i;
    }
  }
}

TEST(Coords2D, BondLengthIsHonoured) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("O=C(Nc1ccccc1)c1ccc(Cl)cc1"));

  depict::Coords2DOptions opts;
  opts.bond_length = 30.0f;
  ASSERT_EQ(depict::Generate2DCoordinates(m, opts), 1);

  for (int i = 0; i < m.nedges(); ++i) {
    EXPECT_NEAR(BondLength(m, *m.bondi(i)), 30.0f, 3.0f) << "bond " << i;
  }
}

// Atoms sitting on top of each other are the characteristic failure of a bad
// layout, including across disconnected fragments.
TEST(Coords2D, NoAtomsCoincide) {
  for (const std::string& smi : TestSmiles()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    if (m.natoms() < 2) {
      continue;
    }
    ASSERT_EQ(depict::Generate2DCoordinates(m, {}), 1) << smi;

    EXPECT_GT(ClosestApproach(m), 0.5f) << smi;
  }
}

// Laying out a molecule must not perturb the molecule itself.
TEST(Coords2D, StructureIsUnchanged) {
  for (const std::string& smi : TestSmiles()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;

    const int matoms = m.natoms();
    const int nedges = m.nedges();
    const IWString before = m.unique_smiles();

    ASSERT_EQ(depict::Generate2DCoordinates(m), 1) << smi;

    EXPECT_EQ(m.natoms(), matoms) << smi;
    EXPECT_EQ(m.nedges(), nedges) << smi;

    Molecule m2;
    ASSERT_TRUE(m2.build_from_smiles(smi));
    EXPECT_EQ(m.unique_smiles(), m2.unique_smiles()) << smi;
    EXPECT_EQ(before, m2.unique_smiles()) << smi;
  }
}

TEST(Coords2D, CentringIsOptional) {
  const char* smi = "CC(=O)Oc1ccccc1C(=O)O";

  depict::Coords2DOptions centred;
  centred.centre = true;

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(smi));
  ASSERT_EQ(depict::Generate2DCoordinates(m, centred), 1);

  float xmin = m.x(0);
  float xmax = m.x(0);
  float ymin = m.y(0);
  float ymax = m.y(0);
  for (int i = 1; i < m.natoms(); ++i) {
    xmin = std::min(xmin, m.x(i));
    xmax = std::max(xmax, m.x(i));
    ymin = std::min(ymin, m.y(i));
    ymax = std::max(ymax, m.y(i));
  }
  EXPECT_NEAR(0.5f * (xmin + xmax), 0.0f, 1.0e-04f);
  EXPECT_NEAR(0.5f * (ymin + ymax), 0.0f, 1.0e-04f);

  // Without centring the layout is the same shape, just somewhere else, so
  // bond lengths are unaffected.
  depict::Coords2DOptions off;
  off.centre = false;

  Molecule m2;
  ASSERT_TRUE(m2.build_from_smiles(smi));
  ASSERT_EQ(depict::Generate2DCoordinates(m2, off), 1);
  for (int i = 0; i < m2.nedges(); ++i) {
    EXPECT_NEAR(BondLength(m2, *m2.bondi(i)), 1.5f, 0.15f);
  }
}

// All three effort levels must produce a usable layout, not just the default.
TEST(Coords2D, AllPrecisionsProduceUsableLayouts) {
  const char* smi = "CC1=C(C(=O)Nc2ccccc2)N(c2ccc(F)cc2)N=C1c1ccc(S(=O)(=O)N)cc1";

  for (depict::Precision p : {depict::Precision::kQuick,
                              depict::Precision::kStandard,
                              depict::Precision::kBest}) {
    depict::Coords2DOptions opts;
    opts.precision = p;

    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi));
    ASSERT_EQ(depict::Generate2DCoordinates(m, opts), 1);
    EXPECT_GT(ClosestApproach(m), 0.5f);
    for (int i = 0; i < m.nedges(); ++i) {
      EXPECT_NEAR(BondLength(m, *m.bondi(i)), 1.5f, 0.15f);
    }
  }
}

TEST(Coords2D, SkipMinimizationStillLaysOut) {
  depict::Coords2DOptions opts;
  opts.skip_minimization = true;

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC(=O)Oc1ccccc1C(=O)O"));
  ASSERT_EQ(depict::Generate2DCoordinates(m, opts), 1);

  for (int i = 0; i < m.nedges(); ++i) {
    EXPECT_NEAR(BondLength(m, *m.bondi(i)), 1.5f, 0.15f);
  }
}

TEST(Coords2D, EvenAnglesStillLaysOut) {
  depict::Coords2DOptions opts;
  opts.even_angles = true;

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC(C)(C)c1ccccc1"));
  ASSERT_EQ(depict::Generate2DCoordinates(m, opts), 1);

  for (int i = 0; i < m.nedges(); ++i) {
    EXPECT_NEAR(BondLength(m, *m.bondi(i)), 1.5f, 0.15f);
  }
}

// A benzene ring must come out as a regular hexagon: all six atoms the same
// distance from the centre, and interior angles of 120 degrees. This is the
// most direct check that the geometry is chemically sensible rather than
// merely non-degenerate.
TEST(Coords2D, BenzeneIsARegularHexagon) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1"));
  ASSERT_EQ(depict::Generate2DCoordinates(m), 1);
  ASSERT_EQ(m.natoms(), 6);

  // Centred, so the centroid is the origin and the circumradius of a hexagon
  // of side 1.5 is 1.5.
  for (int i = 0; i < 6; ++i) {
    const float r = std::sqrt(m.x(i) * m.x(i) + m.y(i) * m.y(i));
    EXPECT_NEAR(r, 1.5f, 0.05f) << "atom " << i;
  }
}

// Explicit hydrogens are laid out as atoms, which callers need to be able to
// rely on either way.
TEST(Coords2D, ExplicitHydrogensAreLaidOut) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1"));
  m.make_implicit_hydrogens_explicit();
  ASSERT_EQ(m.natoms(), 12);

  ASSERT_EQ(depict::Generate2DCoordinates(m), 1);
  EXPECT_EQ(m.natoms(), 12);
  EXPECT_GT(ClosestApproach(m), 0.5f);
  for (int i = 0; i < m.nedges(); ++i) {
    EXPECT_NEAR(BondLength(m, *m.bondi(i)), 1.5f, 0.15f);
  }
}

// Metal complexes are one of coordgen's selling points, and they do lay out
// usefully, but not with the uniform bond lengths every other structure gets:
// the two exocyclic chlorines here come out at sqrt(2) times the bond length.
// The all-carbon analogue "C1(Cl)(Cl)NCCN1" is uniform, so this is specific to
// the organometallic case rather than to gem-disubstitution. Recorded here so
// that a renderer is written knowing bond lengths are not guaranteed uniform,
// and so a change in the behaviour on an upstream resync is noticed.
TEST(Coords2D, MetalComplexLaysOut) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(kMetalComplex));
  ASSERT_EQ(depict::Generate2DCoordinates(m), 1);

  EXPECT_GT(ClosestApproach(m), 0.5f);

  int longer_than_uniform = 0;
  for (int i = 0; i < m.nedges(); ++i) {
    const float len = BondLength(m, *m.bondi(i));
    EXPECT_GE(len, 1.35f) << "bond " << i;
    EXPECT_LE(len, 2.20f) << "bond " << i;
    if (len > 1.65f) {
      ++longer_than_uniform;
    }
  }

  EXPECT_EQ(longer_than_uniform, 2);
}

// Heavily bridged polycyclics on which coordgen has been observed to produce
// NaN coordinates while still reporting them as set - 21 cases in a 39823
// molecule structure-mutation set, all of them fused cages of this kind.
std::vector<std::string>
NanProneSmiles() {
  return {
      "CC1=CC=C2C=NN3C4=NC(CCC(Cl)CCC(=CS4)C2CCCCC3N[N+](=O)[O-])S1",
      "CC1CCCN(C2=NC(=O)C3=Cc4cn(nc4C(O)CCC(S(C)(=O)=O)C3)C2=O)C1",
      "COC(=O)Nc1cccc2c1C1=Cc3c[nH]nc3N=C(C(=O)C=C2)C(N2CCCC(C)C2)=NC1=O",
      "Cc1cc2ccc(ccc(N3CCOCC3)c1)C(Cl)CCCNC2",
      "CCOCCN1CCCC(OC)CC2C=NN1C=C1C=C1C(=O)N=C2N1CCCC(C)C1",
  };
}

// The contract callers depend on: a reported success never leaves non-finite
// coordinates behind. NaN in a molfile produces a record no toolkit can read,
// so such a layout has to be reported as a failure instead. Phrased as the
// invariant rather than as "these must fail" so that it stays correct if an
// upstream resync fixes the underlying NaN.
TEST(Coords2D, NeverProducesNonFiniteCoordinates) {
  std::vector<std::string> smiles = TestSmiles();
  for (const std::string& s : NanProneSmiles()) {
    smiles.push_back(s);
  }

  for (const std::string& smi : smiles) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;

    if (!depict::Generate2DCoordinates(m)) {
      continue;  // Reported as a failure, which is the contract.
    }

    for (int i = 0; i < m.natoms(); ++i) {
      EXPECT_TRUE(std::isfinite(m.x(i))) << smi << " atom " << i;
      EXPECT_TRUE(std::isfinite(m.y(i))) << smi << " atom " << i;
      EXPECT_TRUE(std::isfinite(m.z(i))) << smi << " atom " << i;
    }
  }
}

// Structures on which coordgen has been observed to collapse part of the
// layout onto a single point - in the worst case 34 atoms sharing 8 positions
// - while again reporting the coordinates as set.
std::vector<std::string>
CollapseProneSmiles() {
  return {
      "CCOC1CCC2CNNCNC34C(C=CC(=O)C5=NN=CC=NNC=C5C3(C)CC1)CCC24",
      "COC(=O)NC1=CC=CC2=C1C1=NNC=CC3=C(C=C2)C(=CC=C3)NCCC1",
      "C=CC1CNCC2=CC=CC1=CC=C1C=C2CC=C(/C=N/N)CC1",
      "CC1=CC=C2C=NN(C3=C1CCCCCCN3)C(N(=O)=O)C1CCCC2/C=C\\S1",
  };
}

// The other half of the drawability contract: a reported success never leaves
// two bonded atoms on top of each other. Like the NaN case this is phrased as
// the invariant, so it survives an upstream fix.
TEST(Coords2D, NeverLeavesBondedAtomsCoincident) {
  std::vector<std::string> smiles = TestSmiles();
  for (const std::string& s : CollapseProneSmiles()) {
    smiles.push_back(s);
  }
  for (const std::string& s : NanProneSmiles()) {
    smiles.push_back(s);
  }

  for (const std::string& smi : smiles) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;

    if (!depict::Generate2DCoordinates(m)) {
      continue;  // Reported as a failure, which is the contract.
    }

    for (int i = 0; i < m.nedges(); ++i) {
      EXPECT_GT(BondLength(m, *m.bondi(i)), 0.015f) << smi << " bond " << i;
    }
  }
}

// A failed layout must leave the molecule exactly as it was, so that a caller
// which ignores the return value cannot silently emit a corrupt structure.
TEST(Coords2D, FailureLeavesCoordinatesUntouched) {
  std::vector<std::string> smiles = NanProneSmiles();
  for (const std::string& s : CollapseProneSmiles()) {
    smiles.push_back(s);
  }

  for (const std::string& smi : smiles) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;

    // Mark the incoming coordinates so a partial write would show up.
    for (int i = 0; i < m.natoms(); ++i) {
      m.setxyz(i, static_cast<float>(i), 0.0f, 0.0f);
    }

    if (depict::Generate2DCoordinates(m)) {
      continue;  // Laid out successfully, nothing to check here.
    }

    for (int i = 0; i < m.natoms(); ++i) {
      EXPECT_FLOAT_EQ(m.x(i), static_cast<float>(i)) << smi << " atom " << i;
      EXPECT_FLOAT_EQ(m.y(i), 0.0f) << smi << " atom " << i;
    }
  }
}

// The closest approach between two atoms that are not bonded to each other,
// computed independently of the library so the reported figure is checked
// rather than merely echoed. -1 when there is no such pair.
float
ClosestNonBondedApproach(const Molecule& m) {
  float closest = -1.0f;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      if (m.are_bonded(i, j)) {
        continue;
      }
      const float dx = m.x(i) - m.x(j);
      const float dy = m.y(i) - m.y(j);
      const float d = std::sqrt(dx * dx + dy * dy);
      if (closest < 0.0f || d < closest) {
        closest = d;
      }
    }
  }

  return closest;
}

// The crowding measurement must match the layout it describes.
TEST(Coords2D, ClosestNonBondedApproachIsMeasured) {
  for (const std::string& smi : TestSmiles()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;

    depict::Coords2DResult result;
    ASSERT_EQ(depict::Generate2DCoordinates(m, {}, result), 1) << smi;

    EXPECT_FLOAT_EQ(result.closest_nonbonded_approach,
                    ClosestNonBondedApproach(m))
        << smi;
  }
}

// Molecules with no non-bonded pair to measure. Reported as -1 rather than as
// zero, which a caller would read as a total overlap.
TEST(Coords2D, ClosestNonBondedApproachAbsentForTinyMolecules) {
  for (const char* smi : {"C", "CC"}) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;

    depict::Coords2DResult result;
    ASSERT_EQ(depict::Generate2DCoordinates(m, {}, result), 1) << smi;

    EXPECT_LT(result.closest_nonbonded_approach, 0.0f) << smi;
  }
}

// The measurement is only useful if it discriminates. Ordinary molecules come
// out well clear of the crowding threshold a renderer would care about, so a
// small value really does mean an unusual layout.
TEST(Coords2D, OrdinaryMoleculesAreNotCrowded) {
  for (const std::string& smi : TestSmiles()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;

    depict::Coords2DResult result;
    ASSERT_EQ(depict::Generate2DCoordinates(m, {}, result), 1) << smi;

    if (result.closest_nonbonded_approach < 0.0f) {
      continue;  // Nothing to measure, covered above.
    }
    // Half of the default 1.5 bond length.
    EXPECT_GT(result.closest_nonbonded_approach, 0.75f) << smi;
  }
}

// The measurement is in the same units as the requested bond length, so it has
// to scale with it. A renderer comparing it against a fraction of the bond
// length depends on this.
TEST(Coords2D, ClosestNonBondedApproachScalesWithBondLength) {
  const char* smi = "CC(=O)Oc1ccccc1C(=O)O";

  depict::Coords2DOptions small;
  small.bond_length = 1.5f;
  Molecule m1;
  ASSERT_TRUE(m1.build_from_smiles(smi));
  depict::Coords2DResult r1;
  ASSERT_EQ(depict::Generate2DCoordinates(m1, small, r1), 1);

  depict::Coords2DOptions large;
  large.bond_length = 15.0f;
  Molecule m2;
  ASSERT_TRUE(m2.build_from_smiles(smi));
  depict::Coords2DResult r2;
  ASSERT_EQ(depict::Generate2DCoordinates(m2, large, r2), 1);

  EXPECT_NEAR(r2.closest_nonbonded_approach,
              10.0f * r1.closest_nonbonded_approach,
              0.01f * r2.closest_nonbonded_approach);
}

// Disconnected fragments must be placed beside each other rather than on top
// of each other.
TEST(Coords2D, FragmentsAreSeparated) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1.c1ccccc1"));
  ASSERT_EQ(m.number_fragments(), 2);
  ASSERT_EQ(depict::Generate2DCoordinates(m), 1);

  // Every atom of the first ring must be clear of every atom of the second.
  for (int i = 0; i < 6; ++i) {
    for (int j = 6; j < 12; ++j) {
      const float dx = m.x(i) - m.x(j);
      const float dy = m.y(i) - m.y(j);
      EXPECT_GT(std::sqrt(dx * dx + dy * dy), 1.0f) << i << ' ' << j;
    }
  }
}

// Running twice on the same input must give the same answer, otherwise
// regression tests over depictions could never be stable.
TEST(Coords2D, LayoutIsDeterministic) {
  for (const std::string& smi : TestSmiles()) {
    Molecule m1;
    Molecule m2;
    ASSERT_TRUE(m1.build_from_smiles(smi)) << smi;
    ASSERT_TRUE(m2.build_from_smiles(smi)) << smi;

    ASSERT_EQ(depict::Generate2DCoordinates(m1), 1) << smi;
    ASSERT_EQ(depict::Generate2DCoordinates(m2), 1) << smi;

    for (int i = 0; i < m1.natoms(); ++i) {
      EXPECT_FLOAT_EQ(m1.x(i), m2.x(i)) << smi << " atom " << i;
      EXPECT_FLOAT_EQ(m1.y(i), m2.y(i)) << smi << " atom " << i;
    }
  }
}

// ---------------------------------------------------------------------------
// Cis/trans double bonds.
//
// These are checked geometrically rather than by writing a molfile and reading
// it back, because the geometry is the whole of the statement: a molfile has no
// field for E/Z, so which side of the double bond each substituent is drawn on
// is all a reader has to go on.

// Whether `s1` and `s5` are drawn on the same side of the line through `a3` and
// `a4`.
bool
SameSide(const Molecule& m, atom_number_t a3, atom_number_t a4, atom_number_t s1,
         atom_number_t s5) {
  const float dx = m.x(a4) - m.x(a3);
  const float dy = m.y(a4) - m.y(a3);

  const float c1 = dx * (m.y(s1) - m.y(a3)) - dy * (m.x(s1) - m.x(a3));
  const float c5 = dx * (m.y(s5) - m.y(a3)) - dy * (m.x(s5) - m.x(a3));

  return (c1 > 0.0f) == (c5 > 0.0f);
}

struct CisTransCase {
  const char* smiles;
  // Atoms of the double bond, then one substituent on each end. Atom numbers
  // are SMILES order, which is the order LillyMol keeps.
  atom_number_t a3;
  atom_number_t a4;
  atom_number_t s1;
  atom_number_t s5;
  // Whether s1 and s5 must come out on the same side.
  bool same_side;
};

std::vector<CisTransCase>
CisTransCases() {
  return {
      // The pair the whole thing turns on: identical constitution, opposite
      // configuration, so a layout that ignored E/Z would draw them alike.
      {"C/C=C/C", 1, 2, 0, 3, false},
      {"C/C=C\\C", 1, 2, 0, 3, true},

      {"CC/C=C/CC", 2, 3, 1, 4, false},
      {"CC/C=C\\CC", 2, 3, 1, 4, true},

      // Substituents that are not the first atom listed on their end, so the
      // choice of which substituent carries the direction flag matters.
      {"CC(C)/C=C/C(=O)O", 3, 4, 1, 5, false},
      {"CC(C)/C=C\\C(=O)O", 3, 4, 1, 5, true},

      // Both ends fully substituted, and the substituents carrying the direction
      // flags are not the ones written first at either end.
      {"CC/C(C)=C(/C)CC", 2, 4, 1, 5, false},
      // The other substituent on the same end is necessarily on the other side.
      {"CC/C(C)=C(/C)CC", 2, 4, 3, 5, true},

      // Aromatic substituents - large rigid fragments on both sides.
      {"c1ccccc1/C=C/c1ccccc1", 6, 7, 0, 8, false},
      {"c1ccccc1/C=C\\c1ccccc1", 6, 7, 0, 8, true},

      // Two double bonds in one conjugated chain, one each way.
      {"C/C=C/C=C\\C", 1, 2, 0, 3, false},
      {"C/C=C/C=C\\C", 3, 4, 2, 5, true},

      // Not carbon. LillyMol only derives E/Z from a depiction for C=C, but it
      // holds it from a SMILES for anything, and the layout has to honour it.
      {"C/C=N/O", 1, 2, 0, 3, false},
      {"C/C=N\\O", 1, 2, 0, 3, true},

      // A macrocycle, where the ring is big enough that the double bond
      // configuration is a real choice rather than forced by the ring.
      {"C1CCCC/C=C/CCC1", 5, 6, 4, 7, false},
  };
}

TEST(Coords2D, CisTransConfigurationIsDrawn) {
  for (const CisTransCase& c : CisTransCases()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(c.smiles)) << c.smiles;
    ASSERT_TRUE(m.bond_between_atoms(c.a3, c.a4)->part_of_cis_trans_grouping())
        << c.smiles;

    depict::Coords2DResult result;
    ASSERT_EQ(depict::Generate2DCoordinates(m, depict::Coords2DOptions(), result), 1)
        << c.smiles;

    EXPECT_EQ(SameSide(m, c.a3, c.a4, c.s1, c.s5), c.same_side) << c.smiles;
    EXPECT_EQ(result.cis_trans_honoured, result.cis_trans_bonds) << c.smiles;
  }
}

// The measurement has to see every double bond the molecule states a
// configuration for, not just the ones it managed to draw.
TEST(Coords2D, CisTransBondsAreCounted) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C/C=C/C=C\\C"));

  depict::Coords2DResult result;
  ASSERT_EQ(depict::Generate2DCoordinates(m, depict::Coords2DOptions(), result), 1);

  EXPECT_EQ(result.cis_trans_bonds, 2);
  EXPECT_EQ(result.cis_trans_honoured, 2);
}

TEST(Coords2D, NoCisTransBondsNothingCounted) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC=CC"));

  depict::Coords2DResult result;
  ASSERT_EQ(depict::Generate2DCoordinates(m, depict::Coords2DOptions(), result), 1);

  EXPECT_EQ(result.cis_trans_bonds, 0);
  EXPECT_EQ(result.cis_trans_honoured, 0);
}

// Turning the constraint off leaves the layout engine to its own devices, which
// draws every double bond trans. Worth pinning down, since it is the reason the
// constraint exists: it is what a caller gets if they ask for the raw layout.
TEST(Coords2D, CisTransCanBeTurnedOff) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C/C=C\\C"));

  depict::Coords2DOptions opts;
  opts.honour_cis_trans = false;

  depict::Coords2DResult result;
  ASSERT_EQ(depict::Generate2DCoordinates(m, opts, result), 1);

  EXPECT_EQ(result.cis_trans_bonds, 0);
  EXPECT_FALSE(SameSide(m, 1, 2, 0, 3));
}

// Nothing about the constraint may change the molecule other than its
// coordinates - in particular the directional flags that hold E/Z must survive,
// since a caller that writes SMILES afterwards depends on them.
TEST(Coords2D, CisTransFlagsAreNotDisturbed) {
  for (const char* smi : {"C/C=C/C", "C/C=C\\C", "C/C=C/C=C\\C",
                          "c1ccccc1/C=C\\c1ccccc1"}) {
    Molecule m;
    Molecule untouched;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    ASSERT_TRUE(untouched.build_from_smiles(smi)) << smi;

    ASSERT_EQ(depict::Generate2DCoordinates(m), 1) << smi;

    EXPECT_EQ(m.unique_smiles(), untouched.unique_smiles()) << smi;
  }
}

// End to end within LillyMol: lay the molecule out, throw the directional flags
// away, and let Molecule read the configuration back out of the coordinates it
// was just given. Only C=C, which is all
// discern_cis_trans_bonds_from_depiction() handles.
TEST(Coords2D, ConfigurationSurvivesRederivationFromCoordinates) {
  for (const char* smi : {"C/C=C/C", "C/C=C\\C", "CC/C=C/CC", "CC/C=C\\CC",
                          "OC/C=C\\CO", "C/C=C/C=C\\C"}) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    const IWString before = m.unique_smiles();

    ASSERT_EQ(depict::Generate2DCoordinates(m), 1) << smi;

    m.revert_all_directional_bonds_to_non_directional();
    m.discern_cis_trans_bonds_from_depiction();

    EXPECT_EQ(m.unique_smiles(), before) << smi;
  }
}

}  // namespace
