#include "Depict/wedge_2d.h"

#include <string>
#include <vector>

#include "Molecule_Lib/chiral_centre.h"

#include "Depict/coords_2d.h"

#include "gtest/gtest.h"

namespace {

// A wedge as it will be written to a molfile: the bond in its stored order,
// which is the order the bond block uses, and which way it points.
struct Wedge {
  atom_number_t a1;
  atom_number_t a2;
  int up;  // 1 for a wedge, 0 for a hash
};

std::vector<Wedge>
WedgesOf(const Molecule& m) {
  std::vector<Wedge> result;

  for (int i = 0; i < m.nedges(); ++i) {
    const Bond* b = m.bondi(i);
    if (b->is_wedge_up()) {
      result.push_back({b->a1(), b->a2(), 1});
    } else if (b->is_wedge_down()) {
      result.push_back({b->a1(), b->a2(), 0});
    }
  }

  return result;
}

// Coordinates are a precondition for wedge assignment, so the tests generate a
// layout first, exactly as make_2d_coordinates does.
void
LayOut(Molecule& m) {
  ASSERT_EQ(depict::Generate2DCoordinates(m), 1) << m.smiles();
}

// Give `to` the layout of `from`. Used to compare two stereoisomers on
// identical geometry, so that any difference in the wedges assigned is
// attributable to the chirality rather than to the two being laid out
// differently.
void
CopyCoordinates(const Molecule& from, Molecule& to) {
  ASSERT_EQ(from.natoms(), to.natoms());
  for (int i = 0; i < from.natoms(); ++i) {
    to.setxyz(i, from.x(i), from.y(i), from.z(i));
  }
}

// Molecules whose every chiral centre really is stereogenic, so all of them
// must end up wedged. Acyclic and ring centres, one and several centres per
// molecule, a quaternary centre, and a bridged system.
std::vector<std::string>
StereogenicSmiles() {
  return {
      "C[C@H](N)C(=O)O",                            // L-alanine
      "C[C@@H](N)C(=O)O",                           // D-alanine
      "C[C@H](O)CC",                                // butan-2-ol
      "N[C@@H](Cc1ccccc1)C(=O)O",                   // phenylalanine
      "C[C@H](N)[C@@H](C)O",                        // two centres
      "O[C@@H]1CCCC[C@H]1N",                        // aminocyclohexanol
      "O[C@H]1CCCC[C@H]1N",                         // its diastereomer
      "C[C@](N)(O)CC",                              // quaternary, acyclic
      "C[C@H]1C[C@@H]2CC[C@H]1C2",                  // bridged bicyclic
      "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",  // glucose
  };
}

// Marked chiral in the SMILES but not stereogenic: both ring branches from each
// centre are identical, so there is no stereochemistry to draw.
std::vector<std::string>
NotStereogenicSmiles() {
  return {
      "C[C@H]1CC[C@@H](C)CC1",  // cis-1,4-dimethylcyclohexane
      "O[C@H]1CC[C@H](O)CC1",   // trans-cyclohexane-1,4-diol
      "C[C@]1(O)CCCCC1",        // 1-methylcyclohexan-1-ol
  };
}

TEST(Wedge2D, EmptyMolecule) {
  Molecule m;

  depict::WedgeResult result;
  EXPECT_EQ(depict::AssignWedgeBonds(m, result), 0);
  EXPECT_EQ(result.chiral_centres, 0);
}

TEST(Wedge2D, NoChiralCentresNoWedges) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC(=O)Oc1ccccc1C(=O)O"));
  LayOut(m);

  depict::WedgeResult result;
  EXPECT_EQ(depict::AssignWedgeBonds(m, result), 0);
  EXPECT_EQ(result.chiral_centres, 0);
  EXPECT_TRUE(WedgesOf(m).empty());
}

TEST(Wedge2D, EveryStereocentreIsWedged) {
  for (const std::string& smi : StereogenicSmiles()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    LayOut(m);
    const int ncentres = m.chiral_centres();
    ASSERT_GT(ncentres, 0) << smi;

    depict::WedgeResult result;
    EXPECT_EQ(depict::AssignWedgeBonds(m, result), ncentres) << smi;
    EXPECT_EQ(result.chiral_centres, ncentres) << smi;
    EXPECT_EQ(result.wedged, ncentres) << smi;
    EXPECT_EQ(result.unresolved, 0) << smi;
    EXPECT_EQ(result.not_stereogenic, 0) << smi;
    EXPECT_EQ(static_cast<int>(WedgesOf(m).size()), ncentres) << smi;
  }
}

// The narrow end of a wedge sits on the stereocentre - that is what makes the
// wedge a statement about that atom - and a molfile puts the bond's first atom
// first. Readers that take the first atom as the stereocentre, as RDKit does,
// discard a wedge written the other way round, so the stereochemistry silently
// fails to travel. Checked as its own property because LillyMol's own reader
// compensates for the reversed form and so cannot detect it.
TEST(Wedge2D, NarrowEndIsAtTheStereocentre) {
  for (const std::string& smi : StereogenicSmiles()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    LayOut(m);
    ASSERT_GT(depict::AssignWedgeBonds(m), 0) << smi;

    for (const Wedge& w : WedgesOf(m)) {
      EXPECT_NE(m.chiral_centre_at_atom(w.a1), nullptr)
          << smi << " wedge " << w.a1 << " to " << w.a2;
    }
  }
}

// A wedge that is shared between two centres, or a second wedge at one centre,
// makes the drawing ambiguous.
TEST(Wedge2D, OneWedgePerCentreAndPerBond) {
  for (const std::string& smi : StereogenicSmiles()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    LayOut(m);
    ASSERT_GT(depict::AssignWedgeBonds(m), 0) << smi;

    std::vector<int> wedges_at_atom(m.natoms(), 0);
    for (const Wedge& w : WedgesOf(m)) {
      ++wedges_at_atom[w.a1];
      ++wedges_at_atom[w.a2];
    }
    for (int i = 0; i < m.natoms(); ++i) {
      EXPECT_LE(wedges_at_atom[i], 1) << smi << " atom " << i;
    }
  }
}

// "Either" is the molfile's way of saying the stereochemistry is unknown, which
// is not what a molecule with known chirality should be written as.
TEST(Wedge2D, WedgesAreUpOrDownNeverEither) {
  for (const std::string& smi : StereogenicSmiles()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    LayOut(m);
    ASSERT_GT(depict::AssignWedgeBonds(m), 0) << smi;

    for (int i = 0; i < m.nedges(); ++i) {
      EXPECT_FALSE(m.bondi(i)->is_wedge_either()) << smi << " bond " << i;
    }
  }
}

// Assigning wedges is a statement about how the molecule is drawn, not a change
// to the molecule. In particular it must not perturb the chirality it is
// describing, which the unique smiles would show.
TEST(Wedge2D, StructureAndChiralityUnchanged) {
  std::vector<std::string> smiles = StereogenicSmiles();
  for (const std::string& s : NotStereogenicSmiles()) {
    smiles.push_back(s);
  }

  for (const std::string& smi : smiles) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    LayOut(m);

    const int matoms = m.natoms();
    const int nedges = m.nedges();
    const IWString before = m.unique_smiles();

    depict::AssignWedgeBonds(m);

    EXPECT_EQ(m.natoms(), matoms) << smi;
    EXPECT_EQ(m.nedges(), nedges) << smi;
    EXPECT_EQ(m.unique_smiles(), before) << smi;
  }
}

TEST(Wedge2D, Idempotent) {
  for (const std::string& smi : StereogenicSmiles()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    LayOut(m);

    ASSERT_GT(depict::AssignWedgeBonds(m), 0) << smi;
    const std::vector<Wedge> first = WedgesOf(m);

    depict::WedgeResult again;
    ASSERT_GT(depict::AssignWedgeBonds(m, again), 0) << smi;
    const std::vector<Wedge> second = WedgesOf(m);

    ASSERT_EQ(first.size(), second.size()) << smi;
    for (uint32_t i = 0; i < first.size(); ++i) {
      EXPECT_EQ(first[i].a1, second[i].a1) << smi;
      EXPECT_EQ(first[i].a2, second[i].a2) << smi;
      EXPECT_EQ(first[i].up, second[i].up) << smi;
    }
  }
}

// A wedge that came in with the molecule was drawn against whatever
// coordinates it arrived with, so after a fresh layout it is no longer a true
// statement. All of them go, including ones on bonds this code would not have
// chosen.
TEST(Wedge2D, WedgesFromTheInputAreDiscarded) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C[C@H](N)C(=O)O"));
  LayOut(m);

  // Not a bond at the stereocentre, so nothing here would clear it otherwise.
  const atom_number_t carboxyl_c = 3;
  const atom_number_t hydroxyl_o = 5;
  ASSERT_TRUE(m.are_bonded(carboxyl_c, hydroxyl_o));
  ASSERT_TRUE(m.set_wedge_bond_between_atoms(carboxyl_c, hydroxyl_o, 1));
  ASSERT_EQ(WedgesOf(m).size(), 1u);

  depict::WedgeResult result;
  ASSERT_EQ(depict::AssignWedgeBonds(m, result), 1);

  const std::vector<Wedge> wedges = WedgesOf(m);
  ASSERT_EQ(wedges.size(), 1u);
  EXPECT_NE(m.chiral_centre_at_atom(wedges[0].a1), nullptr);
  EXPECT_FALSE(m.bond_between_atoms(carboxyl_c, hydroxyl_o)->is_wedge_any());
}

// The sharpest check that the direction is derived rather than guessed: two
// enantiomers, drawn identically, must differ by exactly the sense of the
// wedge. A rule that got the geometry wrong in a way that cancelled out would
// give both the same direction and be caught here.
TEST(Wedge2D, EnantiomersGetOppositeWedges) {
  Molecule l_ala;
  Molecule d_ala;
  ASSERT_TRUE(l_ala.build_from_smiles("C[C@H](N)C(=O)O"));
  ASSERT_TRUE(d_ala.build_from_smiles("C[C@@H](N)C(=O)O"));
  ASSERT_NE(l_ala.unique_smiles(), d_ala.unique_smiles());

  LayOut(l_ala);
  CopyCoordinates(l_ala, d_ala);

  ASSERT_EQ(depict::AssignWedgeBonds(l_ala), 1);
  ASSERT_EQ(depict::AssignWedgeBonds(d_ala), 1);

  const std::vector<Wedge> lw = WedgesOf(l_ala);
  const std::vector<Wedge> dw = WedgesOf(d_ala);
  ASSERT_EQ(lw.size(), 1u);
  ASSERT_EQ(dw.size(), 1u);

  // Same drawing, so the same bond carries the wedge.
  EXPECT_EQ(lw[0].a1, dw[0].a1);
  EXPECT_EQ(lw[0].a2, dw[0].a2);
  EXPECT_NE(lw[0].up, dw[0].up);
}

// Diastereomers differ at one centre only, so drawn identically their wedges
// must differ at that centre and agree at the other.
TEST(Wedge2D, DiastereomersDifferAtOneCentre) {
  Molecule a;
  Molecule b;
  ASSERT_TRUE(a.build_from_smiles("O[C@@H]1CCCC[C@H]1N"));
  ASSERT_TRUE(b.build_from_smiles("O[C@H]1CCCC[C@H]1N"));
  ASSERT_NE(a.unique_smiles(), b.unique_smiles());

  LayOut(a);
  CopyCoordinates(a, b);

  ASSERT_EQ(depict::AssignWedgeBonds(a), 2);
  ASSERT_EQ(depict::AssignWedgeBonds(b), 2);

  const std::vector<Wedge> aw = WedgesOf(a);
  const std::vector<Wedge> bw = WedgesOf(b);
  ASSERT_EQ(aw.size(), 2u);
  ASSERT_EQ(bw.size(), 2u);

  int differ = 0;
  for (uint32_t i = 0; i < aw.size(); ++i) {
    ASSERT_EQ(aw[i].a1, bw[i].a1);
    ASSERT_EQ(aw[i].a2, bw[i].a2);
    if (aw[i].up != bw[i].up) {
      ++differ;
    }
  }
  EXPECT_EQ(differ, 1);
}

// An atom can be marked chiral in a SMILES without being stereogenic. Drawing
// no wedge is the right answer, and must be reported as such rather than as a
// failure to express something.
TEST(Wedge2D, NotStereogenicIsNotAFailure) {
  for (const std::string& smi : NotStereogenicSmiles()) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    LayOut(m);
    const int ncentres = m.chiral_centres();
    ASSERT_GT(ncentres, 0) << smi;

    depict::WedgeResult result;
    EXPECT_EQ(depict::AssignWedgeBonds(m, result), 0) << smi;
    EXPECT_EQ(result.not_stereogenic, ncentres) << smi;
    EXPECT_EQ(result.unresolved, 0) << smi;
    EXPECT_TRUE(WedgesOf(m).empty()) << smi;
  }
}

// Every centre is accounted for in exactly one of the three outcomes, so a
// caller reporting these numbers cannot lose one.
TEST(Wedge2D, EveryCentreIsAccountedFor) {
  std::vector<std::string> smiles = StereogenicSmiles();
  for (const std::string& s : NotStereogenicSmiles()) {
    smiles.push_back(s);
  }

  for (const std::string& smi : smiles) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smi)) << smi;
    LayOut(m);

    depict::WedgeResult result;
    depict::AssignWedgeBonds(m, result);

    EXPECT_EQ(result.chiral_centres,
              result.wedged + result.not_stereogenic + result.unresolved)
        << smi;
    EXPECT_EQ(result.wedged, static_cast<int>(WedgesOf(m).size())) << smi;
  }
}

// With no layout every atom is at the origin, so there is no direction to
// reason about. Nothing may be wedged - a wedge drawn on collapsed coordinates
// would be a statement with no meaning - and it must not be treated as though
// there were nothing to express.
TEST(Wedge2D, NoCoordinatesNothingWedged) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C[C@H](N)C(=O)O"));
  ASSERT_EQ(m.chiral_centres(), 1);

  depict::WedgeResult result;
  EXPECT_EQ(depict::AssignWedgeBonds(m, result), 0);
  EXPECT_EQ(result.chiral_centres, 1);
  EXPECT_EQ(result.wedged, 0);
  EXPECT_EQ(result.unresolved, 1);
  EXPECT_TRUE(WedgesOf(m).empty());
}

// Cis/trans is held on the single bonds beside a double bond, using the same
// field a wedge would occupy. Taking one of those for a wedge would throw the
// E/Z away, so those bonds are off limits.
TEST(Wedge2D, CisTransBondsAreNotUsed) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C/C=C/[C@H](N)O"));
  LayOut(m);
  const IWString before = m.unique_smiles();

  depict::WedgeResult result;
  ASSERT_EQ(depict::AssignWedgeBonds(m, result), 1);

  // The E/Z is still there.
  EXPECT_EQ(m.unique_smiles(), before);

  for (const Wedge& w : WedgesOf(m)) {
    const Bond* b = m.bond_between_atoms(w.a1, w.a2);
    EXPECT_FALSE(b->part_of_cis_trans_grouping())
        << "wedge " << w.a1 << " to " << w.a2;
  }
}

}  // namespace
