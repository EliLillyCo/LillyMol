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

TEST_F(TestStandardisation, TestKetoToEnolNoTerminalDoubleBonds) {
  _smiles = "CC(=O)C";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.set_keto_enol_preference(standardisation::KetoEnol::kToEnol);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  constexpr int kVerbose = false;

  _chemical_standardisation.Activate(CS_KETO_ENOL, kVerbose);

  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestKetoToEnolSimple) {
  _smiles = "CC(=O)CC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.set_keto_enol_preference(standardisation::KetoEnol::kToEnol);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  constexpr int kVerbose = false;

  _chemical_standardisation.Activate(CS_KETO_ENOL, kVerbose);

  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "OC(=CC)C") << "from " << _smiles << " got "
                << _m1.unique_smiles();
}

TEST_F(TestStandardisation, TestKetoToEnolAdjacentNotProcessed) {
  _smiles = "OC(=O)CC(=O)CC(=O)C=CC(O)=O CHEMBL1743220";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.set_keto_enol_preference(standardisation::KetoEnol::kToEnol);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  constexpr int kVerbose = false;

  _chemical_standardisation.Activate(CS_KETO_ENOL, kVerbose);

  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestKetoToEnolNotResolved) {
  _smiles = "CCC(=O)CC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.set_keto_enol_preference(standardisation::KetoEnol::kToEnol);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  constexpr int kVerbose = false;

  _chemical_standardisation.Activate(CS_KETO_ENOL, kVerbose);

  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestKetoToEnolHigherConnectivity) {
  _smiles = "CC(C)C(=O)CC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.set_keto_enol_preference(standardisation::KetoEnol::kToEnol);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  constexpr int kVerbose = false;

  _chemical_standardisation.Activate(CS_KETO_ENOL, kVerbose);

  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "OC(=C(C)C)CC") << "from " << _smiles << " get " << 
            _m1.unique_smiles();
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
  _chemical_standardisation.set_convert_to_canonical_order(standardisation::Canonicalise::kReinterpretSmiles);
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

struct ExternalDirective {
  IWString directive;
  IWString smiles;
  int expected_rc;
  IWString expected_result;
};

class TestExternalDirective : public testing::TestWithParam<ExternalDirective> {
  protected:
    Chemical_Standardisation _chemical_standardisation;
    Molecule _m;
};

TEST_P(TestExternalDirective, Tests) {
  const auto params = GetParam();

  ASSERT_TRUE(_m.build_from_smiles(params.smiles));

  const_IWSubstring fname_not_used;
  EXPECT_TRUE(_chemical_standardisation.ReadExternalProto(params.directive, fname_not_used));
  EXPECT_EQ(_chemical_standardisation.process(_m), params.expected_rc);
  if (params.expected_rc > 0) {
    EXPECT_EQ(_m.unique_smiles(), params.expected_result) << " got " <<
                _m.unique_smiles() << " expected " << params.expected_result;
  }
}
INSTANTIATE_TEST_SUITE_P(TestExternalDirective, TestExternalDirective, testing::Values(
  ExternalDirective({R"(smarts: "[ND1H1]=[CD3]-[OD1H]" smiles: "N-C=O" name: "amide")",
                    "C1C2CCC1C(C2C(=N)O)C(=N)O CHEMBL1717199", 1,
                    "O=C(N)C1C(C(=O)N)C2CC1CC2"}),
  ExternalDirective({R"(smarts: "[NR0]=c1:[cD2]:[cD2]:[nH]:c:c1" smiles: "N-C=CC=N" name: "para-amino")",
                    "N(C(=NC(C)C)O)S(=O)(=O)C1=CNC=CC1=NC1=CC(C)=CC=C1 CHEMBL1148", 1,
                    "OC(=NC(C)C)NS(=O)(=O)c1c(Nc2cc(C)ccc2)cc[n]c1"})
));


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
  ASSERT_TRUE(_chemical_standardisation.process(_m)) << "doing " << params.directives[0] <<
              ' ' << params.smiles;
  EXPECT_EQ(_m.unique_smiles(), params.expected) << "got " << 
            _m.unique_smiles() << " expected " << params.expected;

  for (int i = 0; i < 10; ++i) {
    _m.build_from_smiles(params.smiles);
    const IWString s = _m.random_smiles();
    ASSERT_TRUE(_m.build_from_smiles(s));
    ASSERT_TRUE(_chemical_standardisation.process(_m)) << "doing " << params.directives[0] <<
                ' ' << params.smiles << " rsmi " << s;
    EXPECT_EQ(_m.unique_smiles(), params.expected) << "got " << 
              _m.unique_smiles() << " expected " << params.expected << " rsmi " << s;
  }
}
INSTANTIATE_TEST_SUITE_P(TestStandardisationP, TestStandardisationP, testing::Values(
  ForStd{{"rvnv5"}, "N1(=NC(=N(=O)C2=CC(=CC=C12)OCCCN1CCOCC1)CC)=O CHEMBL553213", 
         "CCc1[n][n+]([O-])c2c([n+]1[O-])cc(OCCCN1CCOCC1)cc2"},
  ForStd{{"isotope"}, "[2H]-C", "C[H]"},
  ForStd{{"isotope", "all"}, "[2H]-C", "C"},

  ForStd{{"to2ap"}, "O=C(NN)c1cc[nH]c(=N)c1 CHEMBL3091876", "O=C(NN)c1cc([n]cc1)N"},
  ForStd{{"to2ap"}, "N=C(N)N=c1sc([n][nH]1)c1c(C)cccc1 CHEMBL1188079", "NC(=N)Nc1sc(c2c(C)cccc2)[n][n]1"},

  ForStd{{"isoxazole"}, "C12=C(SN=C1-O)CCCC2NC CHEMBL150489", "O=c1[nH]sc2c1C(NC)CCC2"},

  ForStd{{"oxopyrimidine"}, "C1(=C(NC(=NC1=O)C)NCC1=CN=CC=C1)C#N CHEMBL17125", "O=c1[nH]c(C)[n]c(NCc2c[n]ccc2)c1C#N"},
  ForStd{{"oxopyrimidine"}, "C1CCCC2=C1NC(=NC2=O)N=C(NC1=CC=C(C=C1)OC)N CHEMBL4932203",
                            "O=c1[nH]c([n]c2c1CCCC2)N=C(Nc1ccc(OC)cc1)N"},

//We are no longer standardising fused pyrazoles.
//ForStd{{"pyrazole"}, "C1=C2C(=NC(=NC2=NN1)C)O CHEMBL154781",
//              "Oc1c2c[n][nH]c2[n]c(C)[n]1"}

  ForStd{{"indoleh"}, "BrC1=CNC2=NC=NC2=C1 CHEMBL1562708",
                            "Brc1c[n]c2[n]c[nH]c2c1"},
  ForStd{{"indoleh"}, "CC1=NC(=C2C=CN=C2N1)NC(=O)[C@@H]1CCCOC1 CHEMBL5442493",
                      "O=C(Nc1c2c([nH]cc2)[n]c(C)[n]1)[C@H]1COCCC1"},
// Turned off because of the change to connectivity testing.
// Re-enable once that is resolved.
//ForStd{{"guan"}, "C(=NC1=CC=CN=C1)(NC#N)N1CCCCC1 CHEMBL86956",
//                          "N#CN=C(N1CCCCC1)Nc1c[n]ccc1"},

  ForStd{{"rvnv5o"}, "C1(=N(=O)C=CC=C1)O CHEMBL3898657",
                            "Oc1[n+]([O-])cccc1"},
  ForStd{{"rvnv5o"}, "N1(=C(C)CCC1C)=O CHEMBL332704",
                            "[O-][N+]1=C(CCC1C)C"}
));

struct ForStdKetoEnol {
  std::vector<IWString> directives;
  standardisation::KetoEnol keto_enol;
  IWString smiles;
  IWString expected;
};

class TestStandardisationKetoEnol : public testing::TestWithParam<ForStdKetoEnol> {
  protected:
    Chemical_Standardisation _chemical_standardisation;
    Molecule _m;
};

TEST_P(TestStandardisationKetoEnol, Tests) {
  const auto& params = GetParam();

  static constexpr int kVerbose = 1;

  for (const IWString& directive : params.directives) {
    ASSERT_TRUE(_chemical_standardisation.Activate(directive, kVerbose));
  }
  _chemical_standardisation.set_keto_enol_preference(params.keto_enol);

  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_TRUE(_chemical_standardisation.process(_m));

  EXPECT_EQ(_m.unique_smiles(), params.expected) << "got " << 
            _m.unique_smiles() << " expected " << params.expected << ' ' <<
            _m.name();

  ASSERT_TRUE(_m.valence_ok()) << "invalid valence " << _m.smiles() << ' ' <<
              _m.name();

  for (int i = 0; i < 10; ++i) {
    _m.build_from_smiles(params.smiles);
    const IWString s = _m.random_smiles();
    _m.build_from_smiles(s);
    ASSERT_TRUE(_chemical_standardisation.process(_m));
    EXPECT_EQ(_m.unique_smiles(), params.expected) << " from " << s << " got "
                << _m.unique_smiles() << " expected " << params.expected << '\n';
  }
}
INSTANTIATE_TEST_SUITE_P(TestStandardisationKetoEnol, TestStandardisationKetoEnol, testing::Values(
  ForStdKetoEnol{{"isoxazole"}, standardisation::KetoEnol::kToEnol, "C12=C(SNC1=O)CCCC2NC CHEMBL150489",
                 "Oc1[n]sc2c1C(NC)CCC2"},
  ForStdKetoEnol{{"isoxazole"}, standardisation::KetoEnol::kToKeto, "Oc1[n]sc2c1C(NC)CCC2 CHEMBL150489",
                 "O=c1[nH]sc2c1C(NC)CCC2"},

  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToEnol, "C1=CNNC1=O",
                 "Oc1[n][nH]cc1"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToKeto, "Oc1[n][nH]cc1",
                 "O=c1[nH][nH]cc1"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToEnol, "C1(=C(O)C(=O)NN1CCC)O CHEMBL366101.toEnol",
                "CCC[n]1[n]c(c(c1O)O)O"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToKeto, "C1=CC(=CC=C1C1=NN(C)C(=C1)O)Cl CHEMBL4743178",
                "Clc1ccc(c2[nH][n](c(=O)c2)C)cc1"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToKeto, "C(=O)(O)C1=CC=NC(=C1)N1N=CC=C1O CHEMBL4287317",
                "OC(=O)c1cc([n]2[nH]ccc2=O)[n]cc1"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToKeto, "N1N=C(O)C=C1 CHEMBL4227850.toKeto",
                "O=c1[nH][nH]cc1"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToEnol, "O=c1[nH][nH]cc1 CHEMBL4227850.ToEnol",
                "Oc1[n][nH]cc1"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToKeto, "OC1=C(C)C(=NN1)C CHEMBL1476106.toKeto",
                "O=c1[nH][nH]c(c1C)C"},
  // Note that the round trip does not return to the starting molecule. We would then
  // need to run the molecule through pyrazole standardisation.
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToEnol, "O=c1[nH][nH]c(c1C)C CHEMBL1476106.toEnol",
                "Oc1[n][nH]c(c1C)C"},

  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToKeto, "C1(=CC(=NN1)O)N toKeto",
                "O=c1cc(N)[nH][nH]1"},
  // Another transformation changes the exocyclic bond
  ForStdKetoEnol{{"pyrazolone", "enol-fused"}, standardisation::KetoEnol::kToKeto, "OC1=NN(C)C2=NC(C)=CC(C)=C12 toKeto",
                "O=c1[nH][n](c2[n]c(cc(c12)C)C)C"}


));

struct NoChange {
  std::vector<IWString> directives;
  IWString smiles;
};

class TestStandardisationNoChange : public testing::TestWithParam<NoChange> {
  protected:
    Chemical_Standardisation _chemical_standardisation;
    Molecule _m;
};

TEST_P(TestStandardisationNoChange, Tests) {
  const auto params = GetParam();

  static constexpr int kVerbose = 0;

  for (const IWString& directive : params.directives) {
    ASSERT_TRUE(_chemical_standardisation.Activate(directive, kVerbose));
  }
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_FALSE(_chemical_standardisation.process(_m)) << params.directives[0] << " changed " <<
                params.smiles;
}
INSTANTIATE_TEST_SUITE_P(TestStandardisationNoChange, TestStandardisationNoChange, testing::Values(
  NoChange{{"to2ap"}, "C12=NC(=N)C=CN1[C@H]1O[C@H]([C@H](O)[C@H]1O2)CO CHEMBL4303543"},
  NoChange{{"to2ap"}, "CC=CC1=C(C(=CC(=C1)CC1=CNC(=N)NC1=N)CCC)OC CHEMBL528943"},
  NoChange{{"indoleh"}, "CC1=CC(=O)N2N=CC=C2N1 CHEMBL4525057"},
  NoChange{{"guan"}, "C(=N)(NC1=C2C(=CC=C1)C=CC=C2)N1CCCCC1 CHEMBL91145"}
));

}  // namespace
