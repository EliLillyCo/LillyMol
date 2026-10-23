#include <memory>

#include "gtest/gtest.h"

#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Tools/iwecfp_lib.h"

namespace iwecfp {
namespace {

TEST(IwecfpSubset, NullSubsetUsesExistingPath) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));

  atype_t atype[3] = {6, 6, 8};

  Iwecfp fp1;
  Sparse_Fingerprint_Creator sfc1;
  EXPECT_EQ(fp1.Fingerprint(m, atype, &sfc1), FingerprintResult::kOk);

  Iwecfp fp2;
  Sparse_Fingerprint_Creator sfc2;
  EXPECT_EQ(fp2.Fingerprint(m, atype, nullptr, &sfc2), FingerprintResult::kOk);

  EXPECT_EQ(sfc1, sfc2);
}

TEST(IwecfpSubset, EmptySubsetGeneratesNoFingerprint) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));

  atype_t atype[3] = {6, 6, 8};
  int include_atom[3] = {0, 0, 0};

  Iwecfp fp;
  Sparse_Fingerprint_Creator sfc;
  EXPECT_EQ(fp.Fingerprint(m, atype, include_atom, &sfc), FingerprintResult::kNoStartAtoms);
  EXPECT_EQ(sfc.nbits(), 0u);
}

TEST(IwecfpSubset, AllAtomsIncludedMatchesFullFingerprint) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));

  atype_t atype[3] = {6, 6, 8};
  int include_atom[3] = {1, 1, 1};

  Iwecfp fp1;
  Sparse_Fingerprint_Creator sfc1;
  EXPECT_EQ(fp1.Fingerprint(m, atype, &sfc1), FingerprintResult::kOk);

  Iwecfp fp2;
  Sparse_Fingerprint_Creator sfc2;
  EXPECT_EQ(fp2.Fingerprint(m, atype, include_atom, &sfc2), FingerprintResult::kOk);

  EXPECT_EQ(sfc1, sfc2);
}

TEST(IwecfpSubset, ExpansionDoesNotCrossExcludedAtoms) {
  Molecule parent;
  ASSERT_TRUE(parent.build_from_smiles("CCC"));

  // Only the terminal atoms are included. The middle atom must be invisible.
  atype_t parent_atype[3] = {6, 6, 8};
  int include_atom[3] = {1, 0, 1};

  Iwecfp subset_fp;
  Sparse_Fingerprint_Creator subset_sfc;
  EXPECT_EQ(subset_fp.Fingerprint(parent, parent_atype, include_atom, &subset_sfc),
            FingerprintResult::kOk);

  Molecule disconnected;
  ASSERT_TRUE(disconnected.build_from_smiles("C.C"));

  atype_t disconnected_atype[2] = {6, 8};

  Iwecfp disconnected_fp;
  Sparse_Fingerprint_Creator disconnected_sfc;
  EXPECT_EQ(disconnected_fp.Fingerprint(disconnected, disconnected_atype,
                                        &disconnected_sfc), FingerprintResult::kOk);

  EXPECT_EQ(subset_sfc, disconnected_sfc);
}

}  // namespace
}  // namespace iwecfp
