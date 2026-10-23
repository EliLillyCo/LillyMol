#include "Utilities/GFP_Tools/gfp_context.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "Foundational/iwbits/fixed_bit_vector.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/molecule.h"

#include "Utilities/GFP_Tools/fixed_size_counted_fingerprint.h"
#include "Utilities/GFP_Tools/sparsefp.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {

using gfp_context::GFPList;
using testing::ElementsAre;
using testing::FloatNear;

IWString
FingerprintRecord(const IWString& tag, const std::vector<int>& bits) {
  fixed_bit_vector::FixedBitVector fp(128);
  for (int bit : bits) {
    fp.set_bit(bit);
  }

  IWString result;
  result << tag << fp.DaylightAsciiRepresentationIncludingNsetInfo() << ">\n";
  return result;
}

IWString
SparseRecord() {
  Sparse_Fingerprint_Creator sfc;
  sfc.hit_bit(100, 2);
  sfc.hit_bit(100000, 3);

  IWString result;
  sfc.daylight_ascii_form_with_counts_encoded("NCA<", result);
  result << '\n';
  return result;
}

IWString
MolecularPropertiesRecord() {
  return FingerprintRecord("MPR<", {0, 3, 5, 11, 17, 31, 47, 63});
}

IWString
FixedCountedRecord() {
  const unsigned char counts[] = {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0};
  Fixed_Size_Counted_Fingerprint_uchar fp;
  fp.construct_from_array(counts, 16);

  std::ostringstream output;
  fp.write_daylight_ascii_representation(output, "FCA<", 1);

  return IWString(output.str().c_str());
}

std::string
TempFileName() {
  const char* tmpdir = std::getenv("TEST_TMPDIR");
  if (tmpdir == nullptr) {
    tmpdir = "/tmp";
  }

  std::string result(tmpdir);
  result.append("/gfp_context_test.gfp");
  return result;
}

std::string
CarbonTempFileName() {
  const char* tmpdir = std::getenv("TEST_TMPDIR");
  if (tmpdir == nullptr) {
    tmpdir = "/tmp";
  }

  std::string result(tmpdir);
  result.append("/gfp_context_carbon_test.gfp");
  return result;
}

std::string
TestDataFile(const char* fname) {
  const char* test_srcdir = std::getenv("TEST_SRCDIR");
  if (test_srcdir != nullptr) {
    for (const char* path : {"/_main/Utilities/GFP_Tools/testdata/", "/Utilities/GFP_Tools/testdata/"}) {
      std::string candidate(test_srcdir);
      candidate.append(path);
      candidate.append(fname);
      if (std::filesystem::exists(candidate)) {
        return candidate;
      }
    }
  }

  std::string candidate("Utilities/GFP_Tools/testdata/");
  candidate.append(fname);
  return candidate;
}

std::string
WriteCarbonGfpFile() {
  const std::string fname = CarbonTempFileName();

  std::ofstream output(fname);
  output << "$SMI<C>\n";
  output << "PCN<Methane>\n";
  output << "FCTS<.E..........2;1;1;1;1>\n";
  output << "|\n";
  output << "$SMI<CC>\n";
  output << "PCN<Ethane>\n";
  output << "FCTS<.U..........2;1;1;1;1>\n";
  output << "|\n";

  return fname;
}

std::string
WriteShuffledStandardTagFile() {
  const char* tmpdir = std::getenv("TEST_TMPDIR");
  if (tmpdir == nullptr) {
    tmpdir = "/tmp";
  }

  std::string fname(tmpdir);
  fname.append("/gfp_context_shuffled_standard_test.gfp");

  IWString contents;
  contents << "$SMI<C>\n";
  contents << "PCN<methane>\n";
  contents << MolecularPropertiesRecord();
  contents << FingerprintRecord("FPMK2<", {1, 2});
  contents << FingerprintRecord("FPIW<", {3, 4});
  contents << FingerprintRecord("FPMK<", {5, 6});
  contents << "|\n";

  std::ofstream output(fname);
  output << contents;

  return fname;
}

std::string
WriteTestGfpFile() {
  const std::string fname = TempFileName();

  IWString contents;
  contents << "$SMI<C>\n";
  contents << "PCN<methane>\n";
  contents << FingerprintRecord("FPA<", {0, 1, 2});
  contents << FingerprintRecord("FPB<", {10, 11});
  contents << SparseRecord();
  contents << MolecularPropertiesRecord();
  contents << FixedCountedRecord();
  contents << "|\n";
  contents << "$SMI<CC>\n";
  contents << "PCN<ethane>\n";
  contents << FingerprintRecord("FPA<", {0, 1});
  contents << FingerprintRecord("FPB<", {10, 11});
  contents << SparseRecord();
  contents << MolecularPropertiesRecord();
  contents << FixedCountedRecord();
  contents << "|\n";
  contents << "$SMI<O>\n";
  contents << "PCN<water>\n";
  contents << FingerprintRecord("FPA<", {20, 21});
  contents << FingerprintRecord("FPB<", {30, 31});
  contents << SparseRecord();
  contents << MolecularPropertiesRecord();
  contents << FixedCountedRecord();
  contents << "|\n";

  std::ofstream output(fname);
  output << contents;

  return fname;
}

TEST(TestGFPContext, ReadFileDiscoversContextAndComputesDistances) {
  const std::string fname = WriteTestGfpFile();

  GFPList gfp;
  ASSERT_TRUE(gfp.ReadFile(fname.c_str()));

  EXPECT_EQ(gfp.size(), 3);
  EXPECT_EQ(gfp.smiles(0), "C");
  EXPECT_EQ(gfp.id(1), "ethane");
  EXPECT_THAT(gfp.context().Tags(), ElementsAre("FCA<", "FPA<", "FPB<", "MPR<", "NCA<"));

  EXPECT_THAT(gfp.Distance(0, 0), FloatNear(0.0f, 1.0e-6f));
  EXPECT_THAT(gfp.Distance(0, 1), FloatNear(1.0f / 15.0f, 1.0e-6f));
  EXPECT_THAT(gfp.Distance(1, 0), FloatNear(gfp.Distance(0, 1), 1.0e-6f));

  const uint64_t hash = gfp.context().context_hash();
  ASSERT_TRUE(gfp.mutable_context().SetWeight("FPB<", 0.0f));
  EXPECT_EQ(gfp.context().context_hash(), hash);
  EXPECT_THAT(gfp.Distance(0, 1), FloatNear(1.0f / 12.0f, 1.0e-6f));

  std::vector<IWString> only_fpa;
  only_fpa.emplace_back("FPA<");
  ASSERT_TRUE(gfp.mutable_context().UseOnly(only_fpa));
  EXPECT_THAT(gfp.Distance(0, 1), FloatNear(1.0f / 3.0f, 1.0e-6f));

  gfp.mutable_context().UseAll();
  EXPECT_THAT(gfp.Distance(0, 1), FloatNear(1.0f / 12.0f, 1.0e-6f));
}

TEST(TestGFPContext, NearestNeighbours) {
  const std::string fname = WriteTestGfpFile();

  GFPList gfp;
  ASSERT_TRUE(gfp.ReadFile(fname.c_str()));

  const auto nbrs = gfp.NearestNeighbours(0, 2);
  ASSERT_EQ(nbrs.size(), 2u);
  EXPECT_EQ(nbrs[0].index, 1);
  EXPECT_THAT(nbrs[0].distance, FloatNear(1.0f / 15.0f, 1.0e-6f));
  EXPECT_EQ(nbrs[1].index, 2);
  EXPECT_THAT(nbrs[1].distance, FloatNear(0.4f, 1.0e-6f));
}

TEST(TestGFPContext, FixedCountedCarbonDistance) {
  const std::string fname = WriteCarbonGfpFile();

  GFPList gfp;
  ASSERT_TRUE(gfp.ReadFile(fname.c_str()));

  EXPECT_EQ(gfp.size(), 2);
  EXPECT_THAT(gfp.context().Tags(), ElementsAre("FCTS<"));
  EXPECT_THAT(gfp.Distance(0, 0), FloatNear(0.0f, 1.0e-6f));
  EXPECT_THAT(gfp.Distance(0, 1), FloatNear(0.5f, 1.0e-6f));
  EXPECT_THAT(gfp.Distance(1, 0), FloatNear(0.5f, 1.0e-6f));
}

TEST(TestGFPContext, StandardGfpGoldenDistances) {
  GFPList gfp;
  ASSERT_TRUE(gfp.ReadFile(TestDataFile("rand10.standard.gfp").c_str()));

  ASSERT_GT(gfp.size(), 3);
  ASSERT_EQ(gfp.id(0), "CHEMBL3460651");
  ASSERT_EQ(gfp.id(1), "CHEMBL3460651.a");
  ASSERT_EQ(gfp.id(3), "CHEMBL1417367");

  EXPECT_THAT(gfp.Distance(0, 1), FloatNear(0.0421f, 0.0001f));
  EXPECT_THAT(gfp.Distance(1, 0), FloatNear(0.0421f, 0.0001f));
  EXPECT_THAT(gfp.Distance(3, 0), FloatNear(0.499f, 0.001f));
  EXPECT_THAT(gfp.Distance(0, 3), FloatNear(0.499f, 0.001f));
}

TEST(TestGFPContext, NearestNeighboursWithinDistance) {
  const std::string fname = WriteTestGfpFile();

  GFPList gfp;
  ASSERT_TRUE(gfp.ReadFile(fname.c_str()));

  const auto close = gfp.NearestNeighboursWithinDistance(0, 0.1f);
  ASSERT_EQ(close.size(), 1u);
  EXPECT_EQ(close[0].index, 1);
  EXPECT_THAT(close[0].distance, FloatNear(1.0f / 15.0f, 1.0e-6f));

  const auto none = gfp.NearestNeighboursWithinDistance(0, 0.01f);
  EXPECT_TRUE(none.empty());

  const auto all = gfp.NearestNeighboursWithinDistance(0, 1.0f);
  ASSERT_EQ(all.size(), 2u);
  EXPECT_EQ(all[0].index, 1);
  EXPECT_EQ(all[1].index, 2);
}

TEST(TestGFPContext, StandardContextCanonicalOrderMatchesShuffledFile) {
  GFPList from_file;
  ASSERT_TRUE(from_file.ReadFile(WriteShuffledStandardTagFile().c_str()));

  auto standard = std::make_shared<gfp_context::GFPContext>();
  ASSERT_TRUE(standard->BuildStandard());

  EXPECT_EQ(from_file.context().context_hash(), standard->context_hash());
  EXPECT_THAT(from_file.context().Tags(),
              ElementsAre("FPIW<", "FPMK<", "FPMK2<", "MPR<"));
}

TEST(TestGFPContext, BuildFromSpecsMatchesStandard) {
  gfp_context::GFPContext standard;
  ASSERT_TRUE(standard.BuildStandard());

  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::IWMFingerprint(),
      gfp_context::GFPGeneratorSpec::MACCSKeys(true),
      gfp_context::GFPGeneratorSpec::MolecularProperties()};
  gfp_context::GFPContext from_specs;
  ASSERT_TRUE(from_specs.BuildFromSpecs(specs));

  EXPECT_EQ(from_specs.context_hash(), standard.context_hash());
  EXPECT_THAT(from_specs.Tags(), ElementsAre("FPIW<", "FPMK<", "FPMK2<", "MPR<"));

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("CCO ethanol"));
  gfp_context::GFPFingerprint fp1;
  ASSERT_TRUE(standard.Fingerprint(mol, fp1));
  gfp_context::GFPFingerprint fp2;
  ASSERT_TRUE(from_specs.Fingerprint(mol, fp2));
  EXPECT_THAT(standard.Distance(fp1, fp2), FloatNear(0.0f, 1.0e-6f));
}

TEST(TestGFPContext, BuildFromSpecsIsOrderIndependent) {
  std::vector<gfp_context::GFPGeneratorSpec> specs1 = {
      gfp_context::GFPGeneratorSpec::IWMFingerprint(),
      gfp_context::GFPGeneratorSpec::MACCSKeys(true),
      gfp_context::GFPGeneratorSpec::MolecularProperties()};
  std::vector<gfp_context::GFPGeneratorSpec> specs2 = {
      gfp_context::GFPGeneratorSpec::MolecularProperties(),
      gfp_context::GFPGeneratorSpec::MACCSKeys(true),
      gfp_context::GFPGeneratorSpec::IWMFingerprint()};

  gfp_context::GFPContext context1;
  ASSERT_TRUE(context1.BuildFromSpecs(specs1));
  gfp_context::GFPContext context2;
  ASSERT_TRUE(context2.BuildFromSpecs(specs2));

  EXPECT_EQ(context1.context_hash(), context2.context_hash());
  EXPECT_THAT(context1.Tags(), ElementsAre("FPIW<", "FPMK<", "FPMK2<", "MPR<"));
  EXPECT_THAT(context2.Tags(), ElementsAre("FPIW<", "FPMK<", "FPMK2<", "MPR<"));
}

TEST(TestGFPContext, BuildFromSpecsRejectsDuplicateTags) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::IWMFingerprint(),
      gfp_context::GFPGeneratorSpec::IWMFingerprint()};
  gfp_context::GFPContext context;
  EXPECT_FALSE(context.BuildFromSpecs(specs));
}

TEST(TestGFPContext, XLogPReplicates) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::XLogP(4)};
  gfp_context::GFPContext context;
  ASSERT_TRUE(context.BuildFromSpecs(specs));
  EXPECT_THAT(context.Tags(), ElementsAre("NCXLOGP4<"));

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("CCO ethanol"));
  gfp_context::GFPFingerprint fp;
  ASSERT_TRUE(context.Fingerprint(mol, fp));

  const Sparse_Fingerprint& sparse = fp.sparse(0);
  ASSERT_EQ(sparse.nbits(), 4);
  const int count = sparse.count_for_bit(0);
  EXPECT_GT(count, 0);
  EXPECT_EQ(sparse.count_for_bit(1), count);
  EXPECT_EQ(sparse.count_for_bit(2), count);
  EXPECT_EQ(sparse.count_for_bit(3), count);
  EXPECT_EQ(sparse.count_for_bit(4), 0);
  EXPECT_THAT(context.Distance(fp, fp), FloatNear(0.0f, 1.0e-6f));
}

TEST(TestGFPContext, XLogPRejectsInvalidReplicates) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::XLogP(0)};
  gfp_context::GFPContext context;
  EXPECT_FALSE(context.BuildFromSpecs(specs));
}

TEST(TestGFPContext, TPSAReplicates) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::TPSA(4)};
  gfp_context::GFPContext context;
  ASSERT_TRUE(context.BuildFromSpecs(specs));
  EXPECT_THAT(context.Tags(), ElementsAre("NCTPSA4<"));

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("CCO ethanol"));
  gfp_context::GFPFingerprint fp;
  ASSERT_TRUE(context.Fingerprint(mol, fp));

  const Sparse_Fingerprint& sparse = fp.sparse(0);
  ASSERT_EQ(sparse.nbits(), 4);
  const int count = sparse.count_for_bit(0);
  EXPECT_GT(count, 0);
  EXPECT_EQ(sparse.count_for_bit(1), count);
  EXPECT_EQ(sparse.count_for_bit(2), count);
  EXPECT_EQ(sparse.count_for_bit(3), count);
  EXPECT_EQ(sparse.count_for_bit(4), 0);
  EXPECT_THAT(context.Distance(fp, fp), FloatNear(0.0f, 1.0e-6f));
}

TEST(TestGFPContext, TPSARejectsInvalidReplicates) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::TPSA(0)};
  gfp_context::GFPContext context;
  EXPECT_FALSE(context.BuildFromSpecs(specs));
}

TEST(TestGFPContext, FormulaFingerprint) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::FormulaFingerprint()};
  gfp_context::GFPContext context;
  ASSERT_TRUE(context.BuildFromSpecs(specs));
  EXPECT_THAT(context.Tags(), ElementsAre("FCFML<"));

  Molecule ethane;
  ASSERT_TRUE(ethane.build_from_smiles("CC ethane"));
  gfp_context::GFPFingerprint fp1;
  ASSERT_TRUE(context.Fingerprint(ethane, fp1));

  Molecule propane;
  ASSERT_TRUE(propane.build_from_smiles("CCC propane"));
  gfp_context::GFPFingerprint fp2;
  ASSERT_TRUE(context.Fingerprint(propane, fp2));

  EXPECT_THAT(context.Distance(fp1, fp1), FloatNear(0.0f, 1.0e-6f));
  EXPECT_GT(context.Distance(fp1, fp2), 0.0f);
  EXPECT_LT(context.Distance(fp1, fp2), 1.0f);
}

TEST(TestGFPContext, CATSFingerprint) {
  if (std::getenv("LILLYMOL_HOME") == nullptr) {
    GTEST_SKIP() << "LILLYMOL_HOME not set";
  }

  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::CATS(4, true)};
  gfp_context::GFPContext context;
  ASSERT_TRUE(context.BuildFromSpecs(specs));
  EXPECT_THAT(context.Tags(), ElementsAre("NCCATS4<"));

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("CCO ethanol"));
  gfp_context::GFPFingerprint fp;
  ASSERT_TRUE(context.Fingerprint(mol, fp));
  EXPECT_THAT(context.Distance(fp, fp), FloatNear(0.0f, 1.0e-6f));
}

TEST(TestGFPContext, CATSSuppressedHydrophobicPairsTag) {
  if (std::getenv("LILLYMOL_HOME") == nullptr) {
    GTEST_SKIP() << "LILLYMOL_HOME not set";
  }

  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::CATS(4, false)};
  gfp_context::GFPContext context;
  ASSERT_TRUE(context.BuildFromSpecs(specs));
  EXPECT_THAT(context.Tags(), ElementsAre("NCCATSP4<"));
}

TEST(TestGFPContext, CATSRejectsInvalidMaxPathLength) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::CATS(0, true)};
  gfp_context::GFPContext context;
  EXPECT_FALSE(context.BuildFromSpecs(specs));
}

TEST(TestGFPContext, AtomPairFingerprint) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::AtomPair(1, 3, IWString("UST:Y"), false)};
  gfp_context::GFPContext context;
  ASSERT_TRUE(context.BuildFromSpecs(specs));
  EXPECT_THAT(context.Tags(), ElementsAre("NCAP1M3USTY<"));

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("CCO ethanol"));
  gfp_context::GFPFingerprint fp;
  ASSERT_TRUE(context.Fingerprint(mol, fp));
  EXPECT_GT(fp.sparse(0).nbits(), 0);
  EXPECT_THAT(context.Distance(fp, fp), FloatNear(0.0f, 1.0e-6f));
}

TEST(TestGFPContext, AtomPairOutOfRangeTag) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::AtomPair(0, 2, IWString("UST:Y"), true)};
  gfp_context::GFPContext context;
  ASSERT_TRUE(context.BuildFromSpecs(specs));
  EXPECT_THAT(context.Tags(), ElementsAre("NCAPT0M2USTY<"));
}

TEST(TestGFPContext, AtomPairRejectsInvalidSpecs) {
  gfp_context::GFPContext context;
  std::vector<gfp_context::GFPGeneratorSpec> negative_min = {
      gfp_context::GFPGeneratorSpec::AtomPair(-1, 3, IWString("UST:Y"), false)};
  EXPECT_FALSE(context.BuildFromSpecs(negative_min));

  std::vector<gfp_context::GFPGeneratorSpec> max_less_than_min = {
      gfp_context::GFPGeneratorSpec::AtomPair(4, 3, IWString("UST:Y"), false)};
  EXPECT_FALSE(context.BuildFromSpecs(max_less_than_min));

  std::vector<gfp_context::GFPGeneratorSpec> bad_atom_type = {
      gfp_context::GFPGeneratorSpec::AtomPair(1, 3, IWString("UST:Y-"), false)};
  EXPECT_FALSE(context.BuildFromSpecs(bad_atom_type));
}

TEST(TestGFPContext, RingSubstitutionFingerprint) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::RingSubstitution()};
  gfp_context::GFPContext context;
  ASSERT_TRUE(context.BuildFromSpecs(specs));
  EXPECT_THAT(context.Tags(), ElementsAre("NCRS<"));

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("c1ccccc1C toluene"));
  gfp_context::GFPFingerprint fp;
  ASSERT_TRUE(context.Fingerprint(mol, fp));
  EXPECT_GT(fp.sparse(0).nbits(), 0);
  EXPECT_THAT(context.Distance(fp, fp), FloatNear(0.0f, 1.0e-6f));
}

TEST(TestGFPContext, ECFingerprintDefaultAtomType) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::ECFingerprint(3, IWString("UST:Z"))};
  gfp_context::GFPContext context;
  ASSERT_TRUE(context.BuildFromSpecs(specs));
  EXPECT_THAT(context.Tags(), ElementsAre("NCEC3USTZ<"));

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("CCO ethanol"));
  gfp_context::GFPFingerprint fp;
  ASSERT_TRUE(context.Fingerprint(mol, fp));
  EXPECT_GT(fp.sparse(0).nbits(), 0);
  EXPECT_THAT(context.Distance(fp, fp), FloatNear(0.0f, 1.0e-6f));
}

TEST(TestGFPContext, ECFingerprintAllowsDifferentAtomTypes) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::ECFingerprint(3, IWString("UST:AY")),
      gfp_context::GFPGeneratorSpec::ECFingerprint(3, IWString("C"))};
  gfp_context::GFPContext context;
  ASSERT_TRUE(context.BuildFromSpecs(specs));
  EXPECT_THAT(context.Tags(), ElementsAre("NCEC3C<", "NCEC3USTAY<"));

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("c1ccccc1 benzene"));
  gfp_context::GFPFingerprint fp;
  ASSERT_TRUE(context.Fingerprint(mol, fp));
  EXPECT_GT(fp.sparse(0).nbits(), 0);
  EXPECT_GT(fp.sparse(1).nbits(), 0);
}

TEST(TestGFPContext, ECFingerprintRejectsInvalidSpecs) {
  gfp_context::GFPContext context;
  std::vector<gfp_context::GFPGeneratorSpec> negative_radius = {
      gfp_context::GFPGeneratorSpec::ECFingerprint(-1, IWString("UST:Z"))};
  EXPECT_FALSE(context.BuildFromSpecs(negative_radius));

  std::vector<gfp_context::GFPGeneratorSpec> bad_atom_type_chars = {
      gfp_context::GFPGeneratorSpec::ECFingerprint(3, IWString("UST:A-Y"))};
  EXPECT_FALSE(context.BuildFromSpecs(bad_atom_type_chars));
}

TEST(TestGFPContext, ALogPReplicates) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::ALogP(4)};
  gfp_context::GFPContext context;
  ASSERT_TRUE(context.BuildFromSpecs(specs));
  EXPECT_THAT(context.Tags(), ElementsAre("NCALOGP4<"));

  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("CCO ethanol"));
  gfp_context::GFPFingerprint fp;
  ASSERT_TRUE(context.Fingerprint(mol, fp));

  const Sparse_Fingerprint& sparse = fp.sparse(0);
  ASSERT_EQ(sparse.nbits(), 4);
  const int count = sparse.count_for_bit(0);
  EXPECT_GT(count, 0);
  EXPECT_EQ(sparse.count_for_bit(1), count);
  EXPECT_EQ(sparse.count_for_bit(2), count);
  EXPECT_EQ(sparse.count_for_bit(3), count);
  EXPECT_EQ(sparse.count_for_bit(4), 0);
  EXPECT_THAT(context.Distance(fp, fp), FloatNear(0.0f, 1.0e-6f));
}

TEST(TestGFPContext, ALogPRejectsInvalidReplicates) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::ALogP(0)};
  gfp_context::GFPContext context;
  EXPECT_FALSE(context.BuildFromSpecs(specs));
}

TEST(TestGFPContext, MaccsCanOmitLevel2) {
  std::vector<gfp_context::GFPGeneratorSpec> specs = {
      gfp_context::GFPGeneratorSpec::MACCSKeys(false)};
  gfp_context::GFPContext context;
  ASSERT_TRUE(context.BuildFromSpecs(specs));

  EXPECT_THAT(context.Tags(), ElementsAre("FPMK<"));
  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("CCO ethanol"));
  gfp_context::GFPFingerprint fp;
  ASSERT_TRUE(context.Fingerprint(mol, fp));
  EXPECT_THAT(context.Distance(fp, fp), FloatNear(0.0f, 1.0e-6f));
}

TEST(TestGFPContext, StandardFingerprintGeneration) {
  auto context = std::make_shared<gfp_context::GFPContext>();
  ASSERT_TRUE(context->BuildStandard());
  EXPECT_THAT(context->Tags(), ElementsAre("FPIW<", "FPMK<", "FPMK2<", "MPR<"));
  EXPECT_TRUE(context->can_generate_fingerprints());

  Molecule ethane;
  ASSERT_TRUE(ethane.build_from_smiles("CC ethane"));
  Molecule propane;
  ASSERT_TRUE(propane.build_from_smiles("CCC propane"));

  gfp_context::GFPFingerprint fp1;
  ASSERT_TRUE(context->Fingerprint(ethane, fp1));
  gfp_context::GFPFingerprint fp2;
  ASSERT_TRUE(context->Fingerprint(propane, fp2));

  EXPECT_THAT(context->Distance(fp1, fp1), FloatNear(0.0f, 1.0e-6f));
  const float d12 = context->Distance(fp1, fp2);
  EXPECT_GE(d12, 0.0f);
  EXPECT_LE(d12, 1.0f);
  EXPECT_THAT(context->Distance(fp2, fp1), FloatNear(d12, 1.0e-6f));
}

TEST(TestGFPContext, StandardListAddAndQueryFingerprint) {
  auto gfp = GFPList::Standard();
  ASSERT_NE(gfp, nullptr);

  Molecule ethane;
  ASSERT_TRUE(ethane.build_from_smiles("CC ethane"));
  Molecule propane;
  ASSERT_TRUE(propane.build_from_smiles("CCC propane"));
  Molecule butane;
  ASSERT_TRUE(butane.build_from_smiles("CCCC butane"));

  ASSERT_TRUE(gfp->Add(ethane));
  ASSERT_TRUE(gfp->Add(propane));
  ASSERT_TRUE(gfp->Add(butane));

  EXPECT_EQ(gfp->size(), 3);
  EXPECT_EQ(gfp->id(0), "ethane");
  EXPECT_THAT(gfp->Distance(0, 0), FloatNear(0.0f, 1.0e-6f));

  gfp_context::GFPFingerprint query;
  ASSERT_TRUE(gfp->mutable_context().Fingerprint(propane, query));
  EXPECT_THAT(gfp->Distance(query, 1), FloatNear(0.0f, 1.0e-6f));

  const auto hits = gfp->NearestNeighbours(query, 2);
  ASSERT_EQ(hits.size(), 2u);
  EXPECT_EQ(hits[0].index, 1);
  EXPECT_THAT(hits[0].distance, FloatNear(0.0f, 1.0e-6f));
}

TEST(TestGFPContext, StandardListAddMoleculesWithoutMetadata) {
  auto gfp = GFPList::Standard();
  ASSERT_NE(gfp, nullptr);

  Molecule ethane;
  ASSERT_TRUE(ethane.build_from_smiles("CC ethane"));
  Molecule propane;
  ASSERT_TRUE(propane.build_from_smiles("CCC propane"));
  Molecule butane;
  ASSERT_TRUE(butane.build_from_smiles("CCCC butane"));

  std::vector<Molecule*> molecules = {&ethane, &propane, &butane};
  ASSERT_TRUE(gfp->AddMolecules(molecules));

  EXPECT_EQ(gfp->size(), 3);
  EXPECT_FALSE(gfp->metadata_stored());
  EXPECT_THAT(gfp->Distance(1, 1), FloatNear(0.0f, 1.0e-6f));

  Molecule pentane;
  ASSERT_TRUE(pentane.build_from_smiles("CCCCC pentane"));
  EXPECT_FALSE(gfp->Add(pentane));
}

TEST(TestGFPContext, StandardListFromMoleculesCanStoreMetadata) {
  Molecule ethane;
  ASSERT_TRUE(ethane.build_from_smiles("CC ethane"));
  Molecule propane;
  ASSERT_TRUE(propane.build_from_smiles("CCC propane"));

  std::vector<Molecule*> molecules = {&ethane, &propane};
  auto gfp = GFPList::StandardFromMolecules(molecules, true, true);
  ASSERT_NE(gfp, nullptr);

  EXPECT_EQ(gfp->size(), 2);
  EXPECT_TRUE(gfp->metadata_stored());
  EXPECT_EQ(gfp->id(0), "ethane");
  EXPECT_EQ(gfp->smiles(1), "CCC");
}

}  // namespace
