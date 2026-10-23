#include "Molecule_Lib/iwmfingerprint.h"

#include "Molecule_Lib/molecule.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {

TEST(TestIWMFingerprint, ExplicitOptionsControlBitCount) {
  const int original_nbits = iwmfingerprint_nbits();

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("CCO ethanol"));

  IWMFingerprintOptions options;
  options.bits_per_fingerprint = 512;
  IWMFingerprint fp512(options);
  ASSERT_TRUE(fp512.construct_fingerprint(mol));
  EXPECT_EQ(fp512.nbits(), 512);

  set_iwmfingerprint_nbits(1024);
  IWMFingerprint fp_default;
  ASSERT_TRUE(fp_default.construct_fingerprint(mol));
  EXPECT_EQ(fp_default.nbits(), 1024);

  IWMFingerprint fp512_again(options);
  ASSERT_TRUE(fp512_again.construct_fingerprint(mol));
  EXPECT_EQ(fp512_again.nbits(), 512);

  set_iwmfingerprint_nbits(original_nbits);
}

TEST(TestIWMFingerprint, WidthConstructorControlsBitVector) {
  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("CCO ethanol"));

  IWMFingerprint fp512(512);
  ASSERT_TRUE(fp512.construct_fingerprint(mol));
  EXPECT_EQ(fp512.nbits(), 512);
  ASSERT_NE(fp512.vector(), nullptr);

  int sum = 0;
  for (int i = 0; i < fp512.nbits(); ++i) {
    sum += fp512.vector()[i];
  }
  EXPECT_GT(sum, 0);
}

TEST(TestIWMFingerprint, LegacySetterControlsDefaultOptions) {
  const int original_nbits = iwmfingerprint_nbits();

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("CCO ethanol"));

  set_iwmfingerprint_nbits(256);
  IWMFingerprint fp;
  ASSERT_TRUE(fp.construct_fingerprint(mol));
  EXPECT_EQ(fp.nbits(), 256);

  set_iwmfingerprint_nbits(original_nbits);
}

}  // namespace
