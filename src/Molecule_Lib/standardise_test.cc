#include <algorithm>

#include "googletest/include/gtest/gtest.h"

#include "aromatic.h"
#include "molecule.h"
#include "smiles.h"
#include "standardise.h"

namespace {


class TestStandardisation : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    Chemical_Standardisation _chemical_standardisation;

    IWString _smiles;
    Molecule _m1;
    Molecule _m2;
};

void
TestStandardisation::SetUp()
{
  set_global_aromaticity_type(Daylight);
  set_unique_smiles_legacy_atom_ordering(true);
}

TEST_F(TestStandardisation, EmptyMolecule)
{
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestAcidYes)
{
  _smiles = "CC(=O)[O-]";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  EXPECT_EQ(_m1.molecular_formula(), "C2O2H3");
  EXPECT_EQ(_m1.smiles(), "CC(=O)[O-]");

  _chemical_standardisation.Activate(CS_ACID, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.molecular_formula(), "C2O2H4");
  EXPECT_EQ(_m1.smiles(), "CC(=O)O");
  EXPECT_EQ(_m1.unique_smiles(), "OC(=O)C");
}

TEST_F(TestStandardisation, TestChargedImidazole)
{
  _smiles = "CN1C=C[N+](CC)=C1";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _m2 = _m1;
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  EXPECT_EQ(_m1.molecular_formula(), "C6N2H11");
  EXPECT_EQ(_m1.smiles(), "CN1C=C[N+](=C1)CC");

  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  _smiles = "CCN1C=C[N+](C)=C1";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_m1.molecular_formula(), "C6N2H11");
  EXPECT_EQ(_m1.unique_smiles(), "C[n+]1c[n](CC)cc1");

  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "C[n]1c[n+](CC)cc1");

  _m1.invalidate_smiles();

  // The transformed molecule should not change

  constexpr int replicates = 10;
  for (int i = 0; i < replicates; ++i) {
    const IWString & smiles = _m1.random_smiles();
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smiles));
    _chemical_standardisation.process(m);
    EXPECT_EQ(m.unique_smiles(), "C[n]1c[n+](CC)cc1");
  }

  // Random variants should all get transformed to the same form.

  for (int i = 0; i < replicates; ++i) {
    const IWString & smiles = _m2.random_smiles();
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smiles));
    _chemical_standardisation.process(m);
    EXPECT_EQ(m.unique_smiles(), "C[n]1c[n+](CC)cc1");
  }
}

TEST_F(TestStandardisation, TestMisdrawnSulfonamideNoChange)
{
  _smiles = "CCS(=O)(=O)NC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_MSDSA, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestMisdrawnSulfonamideChanges)
{
  _smiles = "CCS(=O)(O)=NC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_MSDSA, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "CNS(=O)(=O)CC");
}

TEST_F(TestStandardisation, TestEnoltoKetoYes)
{
  _smiles = "CC(O)=C";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_KETO_ENOL, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "O=C(C)C");
}

TEST_F(TestStandardisation, TestEnoltKetoNoBecauseRing)
{
  _smiles = "C1C(O)=CC1";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_KETO_ENOL, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestEnoltKetoNoBecauseAdjacentKeto)
{
  _smiles = "CC(O)=CC(=O)C";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_KETO_ENOL, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestEnoltKetoNoBecauseAdjacentUnsaturation)
{
  _smiles = "C(=O)(O)C(=CC(=O)C1=CC=CC=C1OC)O";   //  CHEMBL4171775
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestEnoltKetoNoBecauseComplexInterdependency)
{
  _smiles = "C1=CC(=CC(=C1O)C(=O)CC(=O)C=CC1=CC=C(O)C=C1)Cl"; // CHEMBL4208282
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestEnoltKetoNoBecauseAdjacentUnsaturationDoubleBond)
{
  _smiles = "S(C1=NC2=CC(=CC=C2N1)C)CC(O)=C(C(=N)C)C#N"; // CHEMBL3197234
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestChargedPyrazole)
{
  _smiles = "[N+]1(=C(C)C=CN1CC1OC(=O)C(C1)(C1=CC=CC=C1)C1=CC=CC=C1)CC";  // CHEMBL140300
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _m2 = _m1;
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);

  EXPECT_EQ(_m1.unique_smiles(), "O=C1OC(CC1(c1ccccc1)c1ccccc1)C[n+]1[n](c(cc1)C)CC");

  _m1.invalidate_smiles();

  // The transformed molecule should not change

  constexpr int replicates = 10;
  for (int i = 0; i < replicates; ++i) {
    const IWString & smiles = _m1.random_smiles();
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smiles));
    _chemical_standardisation.process(m);
    EXPECT_EQ(m.unique_smiles(), "O=C1OC(CC1(c1ccccc1)c1ccccc1)C[n+]1[n](c(cc1)C)CC");
  }
  // Random variants of the starting molecule should all end up the same.
  for (int i = 0; i < replicates; ++i) {
    const IWString & smiles = _m2.random_smiles();
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smiles));
    _chemical_standardisation.process(m);
    EXPECT_EQ(m.unique_smiles(), "O=C1OC(CC1(c1ccccc1)c1ccccc1)C[n+]1[n](c(cc1)C)CC");
  }
}

// Since the pyrazole algorithm depends on the atom ordering, do the same test with
// different atom orderings.
TEST_F(TestStandardisation, TestChargedPyrazoleIncreasing) {
  _smiles = "ON1C=CC=[N+]1N";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);

  EXPECT_EQ(_m1.unique_smiles(), "O[n+]1[n](N)ccc1");
}

TEST_F(TestStandardisation, TestChargedPyrazoledecreasing) {
  _smiles = "N[N+]1=CC=CN1O";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);

  EXPECT_EQ(_m1.unique_smiles(), "O[n+]1[n](N)ccc1");
}

TEST_F(TestStandardisation, TestChargedImidazole3ConnectedNplus) {
  _smiles = "[O-][n]1cco[n+]1=C";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestHRemoval) {
  _smiles = "CC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _m1.make_implicit_hydrogens_explicit();
  _chemical_standardisation.Activate(CS_XH, /*verbose*/ false);
  EXPECT_GT(_chemical_standardisation.process(_m1), 0);
}

#ifdef REDUNDANT_TEST
// Because of other transformations, this guard
// is never made.
TEST_F(TestStandardisation, NoChargedPyrazolones) {
  _smiles = "[O-]C1=C2C3=[N+](N1)C23";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}
#endif

TEST_F(TestStandardisation, PyrazoloneNoProcess) {
  _smiles = "C1=N(=O)NC(=C1)O";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

// First do a test without canonicalization, and then
// repeat with canonicalisation.
TEST_F(TestStandardisation, TestNoUsmi) {
  _smiles = "N1=C(N)C=CN=C1S";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _smiles = "N1=C(S)N=CC=C1N";
  ASSERT_TRUE(_m2.build_from_smiles(_smiles));

  _chemical_standardisation.activate_all();
  EXPECT_GT(_chemical_standardisation.process(_m1), 0);
  EXPECT_GT(_chemical_standardisation.process(_m2), 0);
  EXPECT_NE(_m1.unique_smiles(), _m2.unique_smiles());
}

// First do a test without canonicalization, and then
// repeat with canonicalisation.
TEST_F(TestStandardisation, TestUsmi) {
  _smiles = "N1=C(N)C=CN=C1S";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _smiles = "N1=C(S)N=CC=C1N";
  ASSERT_TRUE(_m2.build_from_smiles(_smiles));

  _chemical_standardisation.activate_all();
  _chemical_standardisation.set_convert_to_canonical_order(standardise::Canonicalise::kReinterpretSmiles);
  EXPECT_GT(_chemical_standardisation.process(_m1), 0);
  EXPECT_GT(_chemical_standardisation.process(_m2), 0);
  EXPECT_EQ(_m1.unique_smiles(), _m2.unique_smiles());
}

TEST_F(TestStandardisation, TestCminus) {
  _smiles = "O=C1[C-](C(=O)C2=CC=CC=C12)C(=O)C1=CC=CC=C1";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "O=C(C1C(=O)c2c(C1=O)cccc2)c1ccccc1");
}

TEST_F(TestStandardisation, TestSulfonylUrea) {
  _smiles = "CNC(=N)S";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "S=C(NC)N") << 
        "sulfonyl urea not converted " << _m1.unique_smiles() << '\n';
}

TEST_F(TestStandardisation, TestSulfonylUreaRing) {
  _smiles = "N1=C(C(C)CN1C(S)=N)C1=CC=CC(=C1)Cl CHEMBL90428";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "Clc1cc(C2=NN(C(=S)N)CC2C)ccc1") << 
        "sulfonyl urea not converted " << _m1.unique_smiles() << '\n';
}

TEST_F(TestStandardisation, Test124Triazole) {
  _smiles = "SC1=NC(=C(C)N=N1)S CHEMBL3272977";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "Cc1[n][nH]c(=S)[nH]c1=S") << 
        "124 triazole not converted " << _m1.unique_smiles() << '\n';
}

}  // namespace
