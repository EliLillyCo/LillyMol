// Tester for fragment related things.
#include <iostream>

#include "molecule.h"

//#include "googlemock/include/gmock/gmock.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {

using std::cerr;
using std::endl;

using testing::UnorderedElementsAre;

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
  cerr << "Building from " << params.smiles << '\n';
  ASSERT_TRUE(m.build_from_smiles(params.smiles));
  cerr << "Built molecule from " << params.smiles << '\n';
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
  cerr << "Building from " << params.smiles << '\n';
  ASSERT_TRUE(m.build_from_smiles(params.smiles));
  cerr << "Built molecule from " << params.smiles << '\n';
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
}  // namespace

