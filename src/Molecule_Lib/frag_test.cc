// Tester for fragment related things.
#include <iostream>
#include <memory>

#include "molecule.h"

//#include "googlemock/include/gmock/gmock.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {

using testing::UnorderedElementsAre;
using testing::UnorderedElementsAreArray;

TEST(TestFrags, NoAtoms) {
  Molecule m;
  EXPECT_EQ(m.number_fragments(), 0);
}

TEST(TestFrags, SingleFrags) {
}

struct SmilesNfrag {
  IWString smiles;
  int number_fragments;
};

class TestNumberFragments : public testing::TestWithParam<SmilesNfrag> {
};

TEST_P(TestNumberFragments, TestCounts) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.smiles));
  EXPECT_EQ(m.number_fragments(), params.number_fragments);
}

INSTANTIATE_TEST_SUITE_P(TestCounts, TestNumberFragments, testing::Values(
   SmilesNfrag{"C", 1},
   SmilesNfrag{"CC", 1},
   SmilesNfrag{"CCC", 1},
   SmilesNfrag{"CC(C)(C)C", 1},
   SmilesNfrag{"C1CC1", 1},
   SmilesNfrag{"C1CC12CC2", 1},
   SmilesNfrag{"C1CC1C1CC1", 1},
   SmilesNfrag{"C1C2CC12", 1},
   SmilesNfrag{"C12C3C4C1C5C2C3C45", 1},

   SmilesNfrag{"C.C", 2},
   SmilesNfrag{"C.C.C", 3},
   SmilesNfrag{"C.C.C.C", 4},
   SmilesNfrag{"C.C.C.C.C", 5},
   SmilesNfrag{"C1CC1.C1CC1", 2}
));

struct SmilesSameFrag {
  IWString smiles;
  atom_number_t a1;
  atom_number_t a2;

  bool same;
};

class TestFragmentMembership : public testing::TestWithParam<SmilesSameFrag> {
};

TEST_P(TestFragmentMembership, TestFragMembership) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.smiles));
  if (params.same) {
    EXPECT_EQ(m.fragment_membership(params.a1), m.fragment_membership(params.a2));
  } else {
    EXPECT_NE(m.fragment_membership(params.a1), m.fragment_membership(params.a2));
  }
}
INSTANTIATE_TEST_SUITE_P(TestFragMembership, TestFragmentMembership, testing::Values(
   SmilesSameFrag{"CC", 0, 1, true},
   SmilesSameFrag{"C.C", 0, 1, false},
   SmilesSameFrag{"C.C.C", 0, 2, false},
   SmilesSameFrag{"C1CC1.CC", 3, 4, true}
));

TEST(TestCreateSubset, TestCreateSubset) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("FCNO"));

  for (int i = 0; i < m.natoms(); ++i) {
    Set_of_Atoms to_keep;
    to_keep << i;
    Molecule subset = m.create_subset(to_keep);
    EXPECT_EQ(subset.natoms(), 1);
    EXPECT_EQ(subset.atomic_number(0), m.atomic_number(i));
  }

  for (int i = 0; i < m.natoms(); ++i) {
    for (int j = i + 1; j < m.natoms(); ++j)  {
      Set_of_Atoms to_keep;
      to_keep << i << j;
      Molecule subset = m.create_subset(to_keep);
      EXPECT_EQ(subset.natoms(), 2);
      std::vector<atomic_number_t> in_parent {m.atomic_number(i), m.atomic_number(j)};
      EXPECT_THAT(in_parent,
        UnorderedElementsAre(subset.atomic_number(0), subset.atomic_number(1)));
    }
  }

}

TEST(TestAtomsInFragment, FirstCall) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C.CC"));
  EXPECT_EQ(m.atoms_in_fragment(0), 1);
  EXPECT_EQ(m.atoms_in_fragment(1), 2);
}

TEST(TestRingsInFragment, TestRingsInFragment) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCC.C1CC1.C12CC1C2"));
  EXPECT_EQ(m.rings_in_fragment(0), 0);
  EXPECT_EQ(m.rings_in_fragment(1), 1);
  EXPECT_EQ(m.rings_in_fragment(2), 2);
}

struct ForCreateSubsets {
  IWString smiles;
  std::vector<int> subset;
  std::vector<IWString> expected;
  int f0_chiral_centres = 0;
  int f1_chiral_centres = 0;
  int f0_implicit_hydrogen_count = 0;
  int f1_implicit_hydrogen_count = 0;
};
class TestCreateSubsets : public testing::TestWithParam<ForCreateSubsets> {
};

TEST_P(TestCreateSubsets, TestCreateSubsets) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.smiles));
  const int matoms = m.natoms();
  ASSERT_EQ(static_cast<int>(params.subset.size()), m.natoms());

  std::unique_ptr<int[]> subset = std::make_unique<int[]>(matoms);
  for (int i = 0; i < matoms; ++i) {
    subset[i] = params.subset[i];
  }

  std::unique_ptr<int[]> xref = std::make_unique<int[]>(matoms);

  Molecule f0, f1;
  EXPECT_TRUE(m.create_subsets(subset.get(), f0, f1, xref.get()));

  std::vector<IWString> fragment_smiles;
  fragment_smiles.push_back(f0.smiles());
  fragment_smiles.push_back(f1.smiles());

  EXPECT_THAT(fragment_smiles, UnorderedElementsAreArray(params.expected)) << f0.smiles() << ' ' << f1.smiles();
  EXPECT_EQ(f0.chiral_centres(), params.f0_chiral_centres);
  EXPECT_EQ(f1.chiral_centres(), params.f1_chiral_centres);
  if (params.f0_chiral_centres > 0) {
    ASSERT_NE(f0.chiral_centre_in_molecule_not_indexed_by_atom_number(0), nullptr);
    EXPECT_EQ(f0.chiral_centre_in_molecule_not_indexed_by_atom_number(0)->implicit_hydrogen_count(),
              params.f0_implicit_hydrogen_count);
  }
  if (params.f1_chiral_centres > 0) {
    ASSERT_NE(f1.chiral_centre_in_molecule_not_indexed_by_atom_number(0), nullptr);
    EXPECT_EQ(f1.chiral_centre_in_molecule_not_indexed_by_atom_number(0)->implicit_hydrogen_count(),
              params.f1_implicit_hydrogen_count);
  }
}
INSTANTIATE_TEST_SUITE_P(TestCreateSubsets, TestCreateSubsets, testing::Values(
   ForCreateSubsets{"CC", {0, 1}, {"C", "C"}},
   ForCreateSubsets{"CCC", {0, 0, 1}, {"CC", "C"}},
   ForCreateSubsets{"CCCO", {0, 1, 1, 0}, {"C.O", "CC"}},
   ForCreateSubsets{"F[C@](Cl)(Br)I", {0, 0, 0, 0, 0}, {"F[C@](Cl)(Br)I", "."}, 1, 0},
   ForCreateSubsets{"F[C@](Cl)(Br)I", {0, 0, 0, 0, 1}, {"F[C@H](Cl)Br", "I"}, 1, 0, 1, 0},
   ForCreateSubsets{"F[C@](Cl)(Br)I", {1, 1, 1, 1, 0}, {"I", "F[C@H](Cl)Br"}, 0, 1, 0, 1},
   ForCreateSubsets{"F[C@](Cl)(Br)I", {0, 1, 1, 1, 1}, {"F", "Cl[C@@H](Br)I"}, 0, 1, 0, 1},
   ForCreateSubsets{"CCCC[C@](F)(Cl)Br", {1, 1, 1, 1, 0, 0, 0, 0}, {"F[C@@H](Cl)Br", "CCCC"}, 1, 0, 1, 0}
));

}  // namespace

