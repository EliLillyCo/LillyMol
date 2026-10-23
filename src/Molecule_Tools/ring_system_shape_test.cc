#include "Molecule_Tools/ring_system_shape.h"

#include <algorithm>
#include <memory>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {

using ring_system_shape::AnalyseRingSystemShape;
using ring_system_shape::NonRingBranchPointCount;
using ring_system_shape::RingAtomBranchPointCount;
using ring_system_shape::RingSystemShape;
using ring_system_shape::RingSystemShapeClass;

RingSystemShape
AnalyseFirstRingSystem(Molecule& m, isotope_t only_process_isotope = 0) {
  std::unique_ptr<int[]> ring_system_membership = std::make_unique<int[]>(m.natoms());
  const int number_ring_systems =
      m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership.get());
  EXPECT_GT(number_ring_systems, 0);

  RingSystemShape result;
  EXPECT_TRUE(AnalyseRingSystemShape(m, ring_system_membership.get(), 1,
                                     only_process_isotope, result));
  return result;
}

int
NonRingBranchPointCount(Molecule& m) {
  std::unique_ptr<int[]> ring_system_membership = std::make_unique<int[]>(m.natoms());
  m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership.get());

  return ring_system_shape::NonRingBranchPointCount(m, ring_system_membership.get());
}

int
RingAtomBranchPointCount(Molecule& m, isotope_t only_process_isotope = 0) {
  std::unique_ptr<int[]> ring_system_membership = std::make_unique<int[]>(m.natoms());
  m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership.get());

  return ring_system_shape::RingAtomBranchPointCount(m, ring_system_membership.get(),
                                                        only_process_isotope);
}

const ring_system_shape::RingSystemSpan*
FindRingSystemSpan(const RingSystemShape& shape, atom_number_t from) {
  for (const ring_system_shape::RingSystemSpan& ring_system_span :
       shape.ring_system_spans) {
    if (ring_system_span.from == from) {
      return &ring_system_span;
    }
  }

  return nullptr;
}

void
LabelRodLikeExitPair(Molecule& m, isotope_t isotope) {
  RingSystemShape shape = AnalyseFirstRingSystem(m);

  for (const ring_system_shape::RingSystemSpan& span1 : shape.ring_system_spans) {
    for (const atom_number_t atom2 : span1.farthest_atoms) {
      const ring_system_shape::RingSystemSpan* span2 = FindRingSystemSpan(shape, atom2);
      if (span2 == nullptr) {
        continue;
      }
      if (std::find(span2->farthest_atoms.begin(), span2->farthest_atoms.end(),
                    span1.from) == span2->farthest_atoms.end()) {
        continue;
      }

      m.set_isotope(span1.from, isotope);
      m.set_isotope(atom2, isotope);
      return;
    }
  }

  ADD_FAILURE() << "No rod-like exit pair found";
}

atom_number_t
LabelExitAtomWithAttachmentCount(Molecule& m, int attachment_count, isotope_t isotope) {
  RingSystemShape shape = AnalyseFirstRingSystem(m);
  std::vector<int> attachment_count_by_atom(m.natoms(), 0);
  for (const ring_system_shape::RingSystemAttachment& attachment : shape.attachments) {
    ++attachment_count_by_atom[attachment.ring_atom];
  }

  for (int i = 0; i < m.natoms(); ++i) {
    if (attachment_count_by_atom[i] == attachment_count) {
      m.set_isotope(i, isotope);
      return i;
    }
  }

  ADD_FAILURE() << "No ring atom with requested attachment count found";
  return kInvalidAtomNumber;
}

TEST(RingSystemShape, TerminalBenzeneIsNotApplicable) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("Cc1ccccc1"));

  const RingSystemShape shape = AnalyseFirstRingSystem(m);

  EXPECT_EQ(shape.shape_class, RingSystemShapeClass::kNotApplicable);
  EXPECT_EQ(shape.attachments.size(), 1u);
  EXPECT_EQ(shape.ring_system_spans.size(), 1u);
  EXPECT_EQ(shape.ring_system_spans[0].max_separation, 3);
  EXPECT_EQ(shape.ring_system_spans[0].farthest_atoms.size(), 1u);
}

TEST(RingSystemShape, OrthoBenzeneIsNotRodLike) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("Cc1ccccc1C"));

  const RingSystemShape shape = AnalyseFirstRingSystem(m);

  EXPECT_EQ(shape.shape_class, RingSystemShapeClass::kNotRodLike);
  EXPECT_EQ(shape.attachments.size(), 2u);
  EXPECT_EQ(shape.ring_system_spans.size(), 2u);
  EXPECT_EQ(shape.observed_separation, 1);
  EXPECT_EQ(shape.rod_deficit, 2);
  for (const ring_system_shape::RingSystemSpan& ring_system_span :
       shape.ring_system_spans) {
    EXPECT_EQ(ring_system_span.max_separation, 3);
  }
}

TEST(RingSystemShape, MetaBenzeneIsNotRodLike) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("Cc1cccc(C)c1"));

  const RingSystemShape shape = AnalyseFirstRingSystem(m);

  EXPECT_EQ(shape.shape_class, RingSystemShapeClass::kNotRodLike);
  EXPECT_EQ(shape.attachments.size(), 2u);
  EXPECT_EQ(shape.ring_system_spans.size(), 2u);
  EXPECT_EQ(shape.observed_separation, 2);
  EXPECT_EQ(shape.rod_deficit, 1);
}

TEST(RingSystemShape, ParaBenzeneIsRodLike) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("Cc1ccc(C)cc1"));

  const RingSystemShape shape = AnalyseFirstRingSystem(m);

  EXPECT_EQ(shape.shape_class, RingSystemShapeClass::kRodLike);
  EXPECT_EQ(shape.observed_separation, 3);
  EXPECT_EQ(shape.rod_deficit, 0);
  ASSERT_EQ(shape.attachments.size(), 2u);
  ASSERT_EQ(shape.ring_system_spans.size(), 2u);

  const atom_number_t a1 = shape.attachments[0].ring_atom;
  const atom_number_t a2 = shape.attachments[1].ring_atom;

  const ring_system_shape::RingSystemSpan* span1 = FindRingSystemSpan(shape, a1);
  ASSERT_NE(span1, nullptr);
  EXPECT_EQ(span1->max_separation, 3);
  EXPECT_THAT(span1->farthest_atoms, testing::ElementsAre(a2));

  const ring_system_shape::RingSystemSpan* span2 = FindRingSystemSpan(shape, a2);
  ASSERT_NE(span2, nullptr);
  EXPECT_EQ(span2->max_separation, 3);
  EXPECT_THAT(span2->farthest_atoms, testing::ElementsAre(a1));
}

TEST(RingSystemShape, TrisubstitutedBenzeneIsMultiSubstituted) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("Cc1cc(C)cc(C)c1"));

  const RingSystemShape shape = AnalyseFirstRingSystem(m);

  EXPECT_EQ(shape.shape_class, RingSystemShapeClass::kMultiSubstituted);
  EXPECT_EQ(shape.attachments.size(), 3u);
  EXPECT_EQ(shape.ring_system_spans.size(), 3u);
}

TEST(RingSystemShape, IsotopeFilterSelectsLabelledRodLikeExitPair) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("Cc1ccc(C)c(C)c1"));

  const RingSystemShape unfiltered = AnalyseFirstRingSystem(m);
  ASSERT_EQ(unfiltered.shape_class, RingSystemShapeClass::kMultiSubstituted);
  ASSERT_EQ(unfiltered.attachments.size(), 3u);

  LabelRodLikeExitPair(m, 7);
  const RingSystemShape filtered = AnalyseFirstRingSystem(m, 7);

  EXPECT_EQ(filtered.shape_class, RingSystemShapeClass::kRodLike);
  EXPECT_EQ(filtered.attachments.size(), 2u);
  EXPECT_EQ(filtered.ring_system_spans.size(), 2u);
  EXPECT_EQ(filtered.rod_deficit, 0);
}

TEST(RingAtomBranchPointCount, IsotopeFilterOnlyCountsLabelledRingAtoms) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1(CC)(CCC)CCC(C)CC1"));
  EXPECT_EQ(RingAtomBranchPointCount(m), 1);

  const atom_number_t terminal_attachment = LabelExitAtomWithAttachmentCount(m, 1, 3);
  ASSERT_NE(terminal_attachment, kInvalidAtomNumber);
  EXPECT_EQ(RingAtomBranchPointCount(m, 3), 0);

  m.set_isotope(terminal_attachment, 0);
  const atom_number_t branched_attachment = LabelExitAtomWithAttachmentCount(m, 2, 3);
  ASSERT_NE(branched_attachment, kInvalidAtomNumber);
  EXPECT_EQ(RingAtomBranchPointCount(m, 3), 1);
}

TEST(RingSystemShape, MultipleSubstituentsOnOneRingAtomCountAsOneExitPoint) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1(C)(F)CCCCC1"));

  const RingSystemShape shape = AnalyseFirstRingSystem(m);

  EXPECT_EQ(shape.shape_class, RingSystemShapeClass::kNotApplicable);
  EXPECT_EQ(shape.attachments.size(), 2u);
  EXPECT_EQ(shape.ring_system_spans.size(), 1u);
  EXPECT_EQ(shape.ring_system_spans[0].max_separation, 3);
}

TEST(RingSystemShape, TerminalCarbonylOxygenIsNotAnExitPoint) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("O=C1CCCCC1"));

  const RingSystemShape shape = AnalyseFirstRingSystem(m);

  EXPECT_EQ(shape.shape_class, RingSystemShapeClass::kNotApplicable);
  EXPECT_EQ(shape.attachments.size(), 0u);
  EXPECT_EQ(shape.ring_system_spans.size(), 0u);
}

TEST(RingSystemShape, TerminalCarbonylDoesNotDistortRodDeficit) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC1CCC(C)CC1=O"));

  const RingSystemShape shape = AnalyseFirstRingSystem(m);

  EXPECT_EQ(shape.shape_class, RingSystemShapeClass::kRodLike);
  EXPECT_EQ(shape.attachments.size(), 2u);
  EXPECT_EQ(shape.ring_system_spans.size(), 2u);
  EXPECT_EQ(shape.observed_separation, 3);
  EXPECT_EQ(shape.rod_deficit, 0);
}

TEST(RingSystemShape, BiphenylRingAttachmentIsExternal) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1-c1ccccc1"));

  std::unique_ptr<int[]> ring_system_membership = std::make_unique<int[]>(m.natoms());
  const int number_ring_systems =
      m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership.get());
  ASSERT_EQ(number_ring_systems, 2);

  RingSystemShape shape1;
  ASSERT_TRUE(AnalyseRingSystemShape(m, ring_system_membership.get(), 1, shape1));
  EXPECT_EQ(shape1.shape_class, RingSystemShapeClass::kNotApplicable);
  EXPECT_EQ(shape1.attachments.size(), 1u);
  EXPECT_EQ(shape1.ring_system_spans.size(), 1u);
  EXPECT_EQ(ring_system_membership[shape1.attachments[0].outside_atom], 2);

  RingSystemShape shape2;
  ASSERT_TRUE(AnalyseRingSystemShape(m, ring_system_membership.get(), 2, shape2));
  EXPECT_EQ(shape2.shape_class, RingSystemShapeClass::kNotApplicable);
  EXPECT_EQ(shape2.attachments.size(), 1u);
  EXPECT_EQ(shape2.ring_system_spans.size(), 1u);
  EXPECT_EQ(ring_system_membership[shape2.attachments[0].outside_atom], 1);
}

TEST(RingSystemShape, InvalidRingSystemIdFails) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("Cc1ccccc1"));
  std::unique_ptr<int[]> ring_system_membership = std::make_unique<int[]>(m.natoms());
  m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership.get());

  RingSystemShape shape;
  EXPECT_FALSE(AnalyseRingSystemShape(m, ring_system_membership.get(), 2, shape));
  EXPECT_EQ(shape.shape_class, RingSystemShapeClass::kInvalid);
}

TEST(RingSystemShape, Adamantane3) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC1C3CC2CC(C)(CC1C2)C3"));
  std::unique_ptr<int[]> ring_system_membership = std::make_unique<int[]>(m.natoms());
  int number_ring_systems =
      m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership.get());
  ASSERT_EQ(number_ring_systems, 1);

  RingSystemShape shape1;
  ASSERT_TRUE(AnalyseRingSystemShape(m, ring_system_membership.get(), 1, shape1));
  EXPECT_EQ(shape1.shape_class, RingSystemShapeClass::kNotRodLike);
  EXPECT_EQ(shape1.attachments.size(), 2u);
  EXPECT_EQ(shape1.ring_system_spans.size(), 2u);
  EXPECT_EQ(ring_system_membership[shape1.attachments[0].ring_atom], 1);
  EXPECT_EQ(ring_system_membership[shape1.attachments[0].outside_atom], 0);
  EXPECT_EQ(ring_system_membership[shape1.attachments[1].ring_atom], 1);
  EXPECT_EQ(ring_system_membership[shape1.attachments[1].outside_atom], 0);
  EXPECT_EQ(shape1.attachments[0].ring_atom, 1);
  EXPECT_EQ(shape1.attachments[0].outside_atom, 0);
  EXPECT_EQ(shape1.attachments[1].ring_atom, 6);
  EXPECT_EQ(shape1.attachments[1].outside_atom, 7);

  EXPECT_EQ(shape1.ring_system_spans.size(), 2u);
  EXPECT_EQ(shape1.ring_system_spans[0].from, 1);
  EXPECT_EQ(shape1.ring_system_spans[0].max_separation, 4);
}

TEST(RingSystemShape, Adamantane4) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC1C3CC2C(C)C(CC1C2)C3"));
  std::unique_ptr<int[]> ring_system_membership = std::make_unique<int[]>(m.natoms());
  int number_ring_systems =
      m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership.get());
  ASSERT_EQ(number_ring_systems, 1);

  RingSystemShape shape1;
  ASSERT_TRUE(AnalyseRingSystemShape(m, ring_system_membership.get(), 1, shape1));
  EXPECT_EQ(shape1.shape_class, RingSystemShapeClass::kRodLike);
  EXPECT_EQ(shape1.attachments.size(), 2u);
  EXPECT_EQ(shape1.ring_system_spans.size(), 2u);
  EXPECT_EQ(ring_system_membership[shape1.attachments[0].ring_atom], 1);
  EXPECT_EQ(ring_system_membership[shape1.attachments[0].outside_atom], 0);
  EXPECT_EQ(ring_system_membership[shape1.attachments[1].ring_atom], 1);
  EXPECT_EQ(ring_system_membership[shape1.attachments[1].outside_atom], 0);

  EXPECT_EQ(shape1.attachments[0].ring_atom, 1);
  EXPECT_EQ(shape1.attachments[0].outside_atom, 0);
  EXPECT_EQ(shape1.attachments[1].ring_atom, 5);
  EXPECT_EQ(shape1.attachments[1].outside_atom, 6);

  EXPECT_EQ(shape1.ring_system_spans.size(), 2u);
  EXPECT_EQ(shape1.ring_system_spans[0].from, 1);
  EXPECT_EQ(shape1.ring_system_spans[0].max_separation, 4);
}

TEST(NonRingBranchPointCount, StraightAlkaneHasNoBranches) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCCCCC"));

  EXPECT_EQ(NonRingBranchPointCount(m), 0);
}

TEST(NonRingBranchPointCount, TertButylIsNotSubstantialBranching) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC(C)(C)C"));

  EXPECT_EQ(NonRingBranchPointCount(m), 0);
}

TEST(NonRingBranchPointCount, CarbonylCarbonIsNotSubstantialBranching) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC(=O)NC"));

  EXPECT_EQ(NonRingBranchPointCount(m), 0);
}

TEST(NonRingBranchPointCount, SulfoneSulfurIsNotSubstantialBranching) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CS(=O)(=O)C"));

  EXPECT_EQ(NonRingBranchPointCount(m), 0);
}

TEST(NonRingBranchPointCount, TrifluoromethylIsNotSubstantialBranching) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC(F)(F)F"));

  EXPECT_EQ(NonRingBranchPointCount(m), 0);
}

TEST(NonRingBranchPointCount, ThreeLargeSubstituentsCountAsBranching) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCC(CC)CC"));

  EXPECT_EQ(NonRingBranchPointCount(m), 1);
}

TEST(NonRingBranchPointCount, RingAtomsAreSkipped) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1(C)(C)CCCCC1"));

  EXPECT_EQ(NonRingBranchPointCount(m), 0);
}

TEST(RingAtomBranchPointCount, SingleSubstantialRingSubstituentIsNotBranching) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCC1CCCCC1"));

  EXPECT_EQ(RingAtomBranchPointCount(m), 0);
}

TEST(RingAtomBranchPointCount, TwoTerminalSingleAtomSubstituentsAreNotBranching) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1(C)(F)CCCCC1"));

  EXPECT_EQ(RingAtomBranchPointCount(m), 0);
}

TEST(RingAtomBranchPointCount, EthylAndPropylOnOneRingAtomIsBranching) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1(CC)(CCC)CCCCC1"));

  EXPECT_EQ(RingAtomBranchPointCount(m), 1);

  const RingSystemShape shape = AnalyseFirstRingSystem(m);
  EXPECT_EQ(shape.shape_class, RingSystemShapeClass::kNotApplicable);
}

TEST(RingAtomBranchPointCount, BranchedRingAtomSuppressesRodLikeMoleculePolicy) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1(CC)(CCC)CCC(C)CC1"));

  EXPECT_EQ(RingAtomBranchPointCount(m), 1);

  const RingSystemShape shape = AnalyseFirstRingSystem(m);
  EXPECT_EQ(shape.shape_class, RingSystemShapeClass::kRodLike);
}

}  // namespace
