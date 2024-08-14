#include <algorithm>
#include <vector>

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

#ifdef THIS_CASE_DOES_NOT_EXIST
Jan 2024.
Looking at Chembl, this case seems not to exist. Turn off for now.
The existing standardisation does not change [N+]#[C-].
TEST_F(TestStandardisation, ReverseReversedCyano) {
  _smiles = "CCN#C";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  constexpr int kVerbose = 0;
  _chemical_standardisation.Activate(CS_REVERSE_NV5, kVerbose);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "CC[N+]#[C-]") << _m1.unique_smiles() << " not match";
}
#endif

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

TEST_F(TestStandardisation, TestExternalNoSmilesNoSmarts) {
  const_IWSubstring buffer = R"(
)";
  const_IWSubstring fname_not_used;

  EXPECT_FALSE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
}

TEST_F(TestStandardisation, TestExternalNoSmarts) {
  const_IWSubstring buffer = R"(
  smiles: "O=[N+]-[O-]"
)";
  const_IWSubstring fname_not_used;

  EXPECT_FALSE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
}

TEST_F(TestStandardisation, TestExternalNoSmiles) {
  const_IWSubstring buffer = R"(
  smarts: "O=N=O"
)";
  const_IWSubstring fname_not_used;

  EXPECT_FALSE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
}

TEST_F(TestStandardisation, TestExternalNitroToChargeSeparated) {
  const_IWSubstring buffer = R"(
  smiles: "O=[N+]-[O-]"
  smarts: "O=N=O"
)";
  const_IWSubstring fname_not_used;

  EXPECT_TRUE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
  ASSERT_TRUE(_chemical_standardisation.active());

  _smiles = "CN(=O)=O";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_TRUE(_chemical_standardisation.process(_m1));

  EXPECT_EQ(_m1.unique_smiles(), "O=[N+]([O-])C") << "GOt smiles " << _m1.smiles();
}

TEST_F(TestStandardisation, TestExternalChargeSeparatedToNitro) {
  const_IWSubstring buffer = R"(
  smarts: "O=[N+]-[O-]"
  smiles: "O=N=O"
)";
  const_IWSubstring fname_not_used;

  EXPECT_TRUE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
  ASSERT_TRUE(_chemical_standardisation.active());

  _smiles = "C[N+](=O)[O-]";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_TRUE(_chemical_standardisation.process(_m1));

  EXPECT_EQ(_m1.unique_smiles(), "CN(=O)=O") << "GOt smiles " << _m1.smiles();
}

TEST_F(TestStandardisation, TestExternalNeutraliseAcid) {
  const_IWSubstring buffer = R"(
  smarts: "[O-][C,S]=O"
  smiles: "O-*=O"
  name: "acid"
)";
  const_IWSubstring fname_not_used;

  EXPECT_TRUE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
  ASSERT_TRUE(_chemical_standardisation.active());
  _chemical_standardisation.set_append_string_depending_on_what_changed(1);

  _smiles = "[O-]C=O";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _m1.set_name("foo");
  EXPECT_TRUE(_chemical_standardisation.process(_m1));

  EXPECT_EQ(_m1.unique_smiles(), "OC=O") << "GOt smiles " << _m1.smiles();
  EXPECT_EQ(_m1.name(), "foo STD:acid") << "Name mismatch got " << _m1.name();
}

TEST_F(TestStandardisation, TestExternalChargedAcid) {
  const_IWSubstring buffer = R"(
  smarts: "[OH][C,S]=O"
  smiles: "[O-]-*=O"
  name: "acid"
)";
  const_IWSubstring fname_not_used;

  EXPECT_TRUE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
  ASSERT_TRUE(_chemical_standardisation.active());
  _chemical_standardisation.set_append_string_depending_on_what_changed(1);

  _smiles = "OC=O";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _m1.set_name("foo");
  EXPECT_TRUE(_chemical_standardisation.process(_m1));

  EXPECT_EQ(_m1.unique_smiles(), "O=C[O-]") << "GOt smiles " << _m1.smiles();
  EXPECT_EQ(_m1.name(), "foo STD:acid") << "Name mismatch got " << _m1.name();
}

struct ForStd {
  std::vector<IWString> directives;
  IWString smiles;
  IWString expected;
};

class TestStandardisationP: public testing::TestWithParam<ForStd> {
  protected:
    Chemical_Standardisation _chemical_standardisation;
    Molecule _m;
};

TEST_P(TestStandardisationP, Tests) {
  const auto params = GetParam();

  static constexpr int kVerbose = 0;

  for (const IWString& directive : params.directives) {
    ASSERT_TRUE(_chemical_standardisation.Activate(directive, kVerbose));
  }
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_chemical_standardisation.process(_m));
  EXPECT_EQ(_m.unique_smiles(), params.expected) << "got " << 
            _m.unique_smiles() << " expected " << params.expected;
}
INSTANTIATE_TEST_SUITE_P(TestStandardisationP, TestStandardisationP, testing::Values(
  ForStd{{"rvnv5"}, "N1(=NC(=N(=O)C2=CC(=CC=C12)OCCCN1CCOCC1)CC)=O CHEMBL553213", 
         "CCc1[n][n+]([O-])c2c([n+]1[O-])cc(OCCCN1CCOCC1)cc2"},
  ForStd{{"isotope"}, "[2H]-C", "C[H]"},
  ForStd{{"isotope", "all"}, "[2H]-C", "C"}
));

}  // namespace
