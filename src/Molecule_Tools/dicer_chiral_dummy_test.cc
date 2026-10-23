#include <memory>
#include <vector>

#include "Molecule_Lib/chiral_centre.h"

#include "Molecule_Tools/dicer_lib.h"
#include "gtest/gtest.h"

namespace dicer_lib {
namespace {

TEST(ChiralDummy, RoundTripBothEnantiomersAndCoordinationTypes) {
  for (const char* smiles :
       {"Cl[C@H](F)Br", "Cl[C@@H](F)Br", "Cl[C@](F)(Br)I", "Cl[C@@](F)(Br)I"}) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smiles));
    const IWString original = m.unique_smiles();
    const int natoms = m.natoms();
    ASSERT_TRUE(BreakBondWithChiralDummies(m, 1, 0));
    ASSERT_EQ(m.natoms(), natoms + 1);
    EXPECT_EQ(m.chiral_centres(), 1);
    EXPECT_EQ(m.atomic_number(natoms), 0);
    EXPECT_EQ(m.hcount(natoms), 0);
    EXPECT_FALSE(InitialAtomNumber(m.atomi(natoms)));
    ASSERT_TRUE(m.stereo_preserving_substitute(1, natoms, 0));
    m.remove_atom(natoms);
    EXPECT_EQ(m.unique_smiles(), original);
  }
}

TEST(ChiralDummy, BothEndpointsKeepTheirNeighbourSlots) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("Cl[C@H](F)[C@@H](Br)I"));
  const Chiral_Centre c1(*m.chiral_centre_at_atom(1));
  const Chiral_Centre c2(*m.chiral_centre_at_atom(3));
  const int n = m.natoms();
  ASSERT_TRUE(BreakBondWithChiralDummies(m, 1, 3));
  ASSERT_EQ(m.chiral_centres(), 2);
  EXPECT_EQ(m.natoms(), n + 2);
  EXPECT_TRUE(m.are_bonded(1, n));
  EXPECT_TRUE(m.are_bonded(3, n + 1));
  EXPECT_FALSE(m.are_bonded(1, 3));
  Chiral_Centre restored1(*m.chiral_centre_at_atom(1));
  Chiral_Centre restored2(*m.chiral_centre_at_atom(3));
  ASSERT_TRUE(restored1.change_atom_number(n, 3));
  ASSERT_TRUE(restored2.change_atom_number(n + 1, 1));
  EXPECT_EQ(restored1.top_front(), c1.top_front());
  EXPECT_EQ(restored1.top_back(), c1.top_back());
  EXPECT_EQ(restored1.left_down(), c1.left_down());
  EXPECT_EQ(restored1.right_down(), c1.right_down());
  EXPECT_EQ(restored2.top_front(), c2.top_front());
  EXPECT_EQ(restored2.top_back(), c2.top_back());
  EXPECT_EQ(restored2.left_down(), c2.left_down());
  EXPECT_EQ(restored2.right_down(), c2.right_down());
}

TEST(ChiralDummy, SecondCapRejectsWholeCutWithoutMutation) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("Cl[C@H](F)[C@@H](Br)I"));
  ASSERT_TRUE(BreakBondWithChiralDummies(m, 3, 5));
  const IWString original = m.unique_smiles();
  const int natoms = m.natoms();
  EXPECT_FALSE(BreakBondWithChiralDummies(m, 1, 3));
  EXPECT_EQ(m.natoms(), natoms);
  EXPECT_EQ(m.unique_smiles(), original);
}

TEST(ChiralDummy, OriginalWildcardIsNotAnAddedCap) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("*[C@H](F)Cl"));
  std::vector<USPVPTR> provenance(m.natoms());
  initialise_atom_pointers(m, nullptr, provenance.data());
  ASSERT_TRUE(BreakBondWithChiralDummies(m, 1, 3));
  EXPECT_EQ(m.natoms(), 5);
}

TEST(ChiralDummy, AchiralCutsAddNothingIncludingMultipleBonds) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC=O"));
  EXPECT_TRUE(BreakBondWithChiralDummies(m, 1, 2));
  EXPECT_TRUE(BreakBondWithChiralDummies(m, 0, 1));
  EXPECT_EQ(m.natoms(), 3);
}

TEST(ChiralDummy, ReconstructedSubsetsMatchCappedComponents) {
  Molecule parent;
  ASSERT_TRUE(parent.build_from_smiles("Cl[C@H](F)Br"));
  const IWString original = parent.unique_smiles();
  const int subset[] = {0, 1, 1, 1};
  int xref[4];
  Molecule fragment, complement;
  ASSERT_TRUE(CreateChiralDummySubsets(parent, subset, fragment, complement, xref));
  EXPECT_EQ(fragment.natoms(), 4);
  EXPECT_EQ(fragment.chiral_centres(), 1);
  EXPECT_EQ(complement.natoms(), 1);
  EXPECT_EQ(parent.unique_smiles(), original);
  Molecule broken(parent);
  ASSERT_TRUE(BreakBondWithChiralDummies(broken, 1, 0));
  resizable_array_p<Molecule> components;
  broken.create_components(components);
  bool found = false;
  for (Molecule* component : components) {
    if (component->unique_smiles() == fragment.unique_smiles()) {
      found = true;
    }
  }
  EXPECT_TRUE(found);
  const int double_cut[] = {0, 1, 0, 1};
  EXPECT_FALSE(CreateChiralDummySubsets(parent, double_cut, fragment, complement, xref));
}
}  // namespace
}  // namespace dicer_lib
