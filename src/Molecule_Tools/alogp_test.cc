// Tests for alogp

#include <filesystem>
#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/alogp.h"

namespace {

#ifdef OLD_VERSION_NOT_USED__
class TestALogP : public testing::Test {
  protected:
    IWString _smiles;

    Molecule _m;

    // We can compare type atom types.
    IWString _expected_smiles;

  // protected functions.
    void SetUp() override;

    Charge_Assigner _charge_assigner;

    alogp::ALogP _alogp;

  public:
};

void
TestALogP::SetUp() {
  const char* test_srcdir = getenv("TEST_SRCDIR");
  // Some diagnostic stuff to see the relationship between TEST_SRCDIR
  // and where our files actually show up.
#ifdef DEBUG_FILE_PATHS
  std::string qq(test_srcdir);
  qq += "/../alogp_test.runfiles/charge_assigner/";
  for (auto const& dir_entry : std::filesystem::directory_iterator{qq}) {
    std::cerr << dir_entry << '\n';
  }
#endif

  IWString queries_file(test_srcdir);
  queries_file << "/../alogp_test.runfiles/charge_assigner/queries";

  IWString cmd;
  cmd << "F:" << queries_file;

  if (! _charge_assigner.build(cmd)) {
    std::cerr << "TestALogP::SetUp:cannot initialise charge assigner '" << cmd << "'\n";
  } else {
    // std::cerr << "Charge assigner initialised '" << cmd << "'\n";
    _charge_assigner.set_apply_charges_to_molecule(1);
  }

  _alogp.set_label_with_atom_type(1);
}

TEST_F(TestALogP, TestMethane) {
  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  _charge_assigner.process(_m);
  std::optional<double> a = _alogp.LogP(_m);
  ASSERT_NE(a, std::nullopt);
  EXPECT_NEAR(*a, 0.636, 0.001);
}
#endif

struct SmilesExpected {
  IWString smiles;
  float alogp;
  IWString labelled_smiles;
};

class TestAlogpP: public testing::TestWithParam<SmilesExpected> {
  protected:
    Molecule _mol;
    alogp::ALogP _alogp;

    Charge_Assigner _charge_assigner;

  protected:
    void SetUp();
};

void
TestAlogpP::SetUp() {
  const char* test_srcdir = getenv("TEST_SRCDIR");
  // Some diagnostic stuff to see the relationship between TEST_SRCDIR
  // and where our files actually show up.
#ifdef DEBUG_FILE_PATHS
  std::string qq(test_srcdir);
  qq += "/../alogp_test.runfiles/charge_assigner/";
  for (auto const& dir_entry : std::filesystem::directory_iterator{qq}) {
    std::cerr << dir_entry << '\n';
  }
#endif

  IWString queries_file(test_srcdir);
  queries_file << "/../alogp_test.runfiles/charge_assigner/queries";

  IWString cmd;
  cmd << "F:" << queries_file;

  if (! _charge_assigner.build(cmd)) {
    std::cerr << "TestALogP::SetUp:cannot initialise charge assigner '" << cmd << "'\n";
  } else {
    std::cerr << "Charge assigner initialised '" << cmd << "'\n";
    _charge_assigner.set_apply_charges_to_molecule(1);
  }

  _alogp.set_label_with_atom_type(1);
  _alogp.set_use_alcohol_for_acid(1);
  _alogp.set_rdkit_charged_nitrogen(1);
  _alogp.set_rdkit_phoshoric_acid_hydrogen(1);
}

TEST_P(TestAlogpP, TestMolecules) {
  const auto& params = GetParam();

  ASSERT_TRUE(_mol.build_from_smiles(params.smiles)) << "bad smiles " << params.smiles;

  // _charge_assigner.process(_mol);

  std::optional<double> a = _alogp.LogP(_mol);
  ASSERT_NE(a, std::nullopt);
  EXPECT_NEAR(*a, params.alogp, 0.001) << params.smiles <<
        " expect " << static_cast<float>(params.alogp) <<
        "\ncomp " << *a << 
        ' ' << _mol.aromatic_smiles();
  if (! params.labelled_smiles.empty()) {
    EXPECT_EQ(_mol.aromatic_smiles(), params.labelled_smiles) << params.labelled_smiles <<
        " got\n" << _mol.aromatic_smiles() << ' ' << params.alogp << ' ' <<
        _mol.name();
  }
}
INSTANTIATE_TEST_SUITE_P(TestAlogpP, TestAlogpP, testing::Values(
  SmilesExpected{"CC ethane", 1.026, ""},
  SmilesExpected{"CCC propane", 1.416, ""},
  SmilesExpected{"CCCC butane", 1.806, ""},
  SmilesExpected{"C1CC1 cyclopropane", 1.170, ""},
  SmilesExpected{"CC(C)C isobutane", 1.662, ""},
  SmilesExpected{"CC(C)(C)C neopentane", 2.052, ""},
  SmilesExpected{"c1ccccc1 benzene", 1.687, ""},
  SmilesExpected{"CO methanol", -0.392, "[3CH3][50OH]"},
  SmilesExpected{"c1ccccc1O phenol", 1.392, "[18cH]1[18cH][18cH][18cH][18cH][23c]1[50OH]"},
  SmilesExpected{"c1ccccc1S benzenethiol", 1.975, "[18cH]1[18cH][18cH][18cH][18cH][24c]1[68SH]"},
  SmilesExpected{"c1ccccc1F fluoro-benzene", 1.826, "[18cH]1[18cH][18cH][18cH][18cH][14c]1[62F]"},
  SmilesExpected{"c1ccccc1Cl chloro-benzene", 2.340, "[18cH]1[18cH][18cH][18cH][18cH][15c]1[63Cl]"},
  SmilesExpected{"c1ccccc1Br bromo-benzene", 2.449, "[18cH]1[18cH][18cH][18cH][18cH][16c]1[64Br]"},
  SmilesExpected{"c1ccccc1I iodo-benzene", 2.291, "[18cH]1[18cH][18cH][18cH][18cH][17c]1[65I]"},
  SmilesExpected{"CC(F)(F)F 1,1,1-trofluoroethane", 1.569, "[1CH3][4C]([62F])([62F])[62F]"},
  SmilesExpected{"COC dimethyl ether", 0.263, "[3CH3][51O][3CH3]"},
  SmilesExpected{"C12=C(C=CC=C1)C=CN2 indole", 2.168, "[19c]12[19c]([18cH][18cH][18cH][18cH]1)[18cH][18cH][44nH]2"},
  SmilesExpected{"c1ccccc1c1ccccc1 biphenyl", 3.354, "[18cH]1[18cH][18cH][18cH][18cH][20c]1[20c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"OC(C)C isopropyl alcohol", 0.387, "[50OH][4CH]([1CH3])[1CH3]"},
  SmilesExpected{"o1cccc1 furan", 1.280, "[49o]1[18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"s1cccc1 thiophene", 1.748, "[70s]1[18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"CC(=O)O acetic acid", 0.091, "[1CH3][5C](=[57O])[50OH]"},
  SmilesExpected{"CC=O acetaldehyde", 0.205, "[1CH3][5CH]=[57O]"},
  SmilesExpected{"CC#N acetonitrile", 0.530, "[1CH3][7C]#[42N]"},
  SmilesExpected{"C=C methene", 0.802, "[6CH2]=[6CH2]"},
  SmilesExpected{"CC=N ethanimine", 0.656, "[1CH3][5CH]=[40NH]"},
  SmilesExpected{"CCc1ccccc1 ethylbenzene", 2.249, "[1CH3][10CH2][21c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"CC(C)c1ccccc1 isopropylbenezene", 2.810, "[1CH3][11CH]([1CH3])[21c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"CC(C)(C)c1ccccc1 t-butyl benzene", 2.984, "[1CH3][12C]([1CH3])([1CH3])[21c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"C=Cc1ccccc1 styrene", 2.330, "[6CH2]=[26CH][21c]1[18cH][18cH][18cH][18cH][18cH]1"},

  // We only get concordance with RDKit if we assume 2 Hydrogens on the N+.
  SmilesExpected{"CN methylamine", -0.425, "[3CH3][34NH2]"},
  SmilesExpected{"CNC dimethyl methylamine", -0.164, "[3CH3][35NH][3CH3]"},
  SmilesExpected{"n1ccccc1 pyridine", 1.082, "[44n]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"[nH]1cccc1 pyrole", 1.015, "[44nH]1[18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"Nc1ccccc1 aniline", 1.269, "[36NH2][22c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"CNc1ccccc1 N-methylaniline", 1.728, "[3CH3][37NH][22c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"CC=NC N-methylethanimine", 0.707, "[1CH3][5CH]=[39N][3CH3]"},
  SmilesExpected{"CN(C)C trimethylamine", 0.178, "[3CH3][40N]([3CH3])[3CH3]"},
  SmilesExpected{"CN(C)c1ccccc1 N,N-dimethylaniline", 1.753, "[3CH3][41N]([3CH3])[22c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"CNC(=O)C N-METHYLACETAMIDE", -0.248, "[3CH3][35NH][5C](=[57O])[1CH3]"},
  SmilesExpected{"CNC(=O)NC 1,3-DIMETHYLUREA", -0.455, "[3CH3][35NH][5C](=[59O])[35NH][3CH3]"},
  SmilesExpected{"O=C1NN=CN1 CHEMBL1865594", -0.902, "[56O]=[25c]1[44nH][44n][18cH][44nH]1"},
  SmilesExpected{"C1(=O)NC=CC=C1 CHEMBL662", 0.375, "[25c]1(=[56O])[44nH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"O=S(=O)NCC Ethanesulfonamide", -0.878, "[54O]=[69SH](=[54O])[35NH][3CH2][1CH3]"},
  SmilesExpected{"C12=C(SNC1=O)CCNC2 CHEMBL171241", 0.082, "[21c]12[21c]([70s][44nH][25c]1=[56O])[10CH2][3CH2][35NH][10CH2]2"},
  SmilesExpected{"C1=CC=C2C(=C1)NC=N2 benzimidazole", 1.563, "[18cH]1[18cH][18cH][19c]2[19c]([18cH]1)[44nH][18cH][44n]2"},
  SmilesExpected{"S1C(=C(C=C1)NC(N)=N)C(=O)OC CHEMBL4299981", 0.840, "[70s]1[21c]([22c]([18cH][18cH]1)[37NH][5C]([34NH2])=[40NH])[5C](=[58O])[51O][3CH3]"},
  SmilesExpected{"C(=O)(C1=CC=CC(=C1)Br)NCCN CHEMBL128615", 1.138, "[5C](=[58O])([21c]1[18cH][18cH][18cH][16c]([18cH]1)[64Br])[35NH][3CH2][3CH2][34NH2]"},
  SmilesExpected{"O=N(=O)C1=CC(=CC(=C1)C(=O)N)N(=O)=O CHEMBL1437065", 0.602, "[53O]=[46N](=[53O])[22c]1[18cH][22c]([18cH][21c]([18cH]1)[5C](=[58O])[34NH2])[46N](=[53O])=[53O]"},
  SmilesExpected{"C(=O)(C1=CC=C(C=C1)CC(C)C)NO CHEMBL439659", 2.004, "[5C](=[58O])([21c]1[18cH][18cH][21c]([18cH][18cH]1)[10CH2][2CH]([1CH3])[1CH3])[35NH][50OH]"},
  SmilesExpected{"C1(=C(N=C(C)C2=C1N=CN(C2=C)CC)OC)C#N CHEMBL1836266", 2.236, "[21c]1([23c]([44n][21c]([8CH3])[21c]2[22c]1[39N]=[5CH][40N]([26C]2=[6CH2])[3CH2][1CH3])[52O][3CH3])[7C]#[42N]"},
  SmilesExpected{"C12=CC(=NN1C(C)CN(C2=O)C1=CC=CC(=C1)OC)COC1=NC=C(Cl)C=C1 CHEMBL3617639", 3.741, "[21c]12[18cH][21c]([44n][44n]1[11CH]([1CH3])[3CH2][41N]([5C]2=[58O])[22c]1[18cH][18cH][18cH][23c]([18cH]1)[52O][3CH3])[10CH2][52O][23c]1[44n][18cH][15c]([63Cl])[18cH][18cH]1"},
  SmilesExpected{"S(=O)(=O)(N(C)CC(=O)NCC1=CC=CC=N1)C1=CC(=CC=C1OC)C CHEMBL1736135", 1.336, "[69S](=[54O])(=[54O])([40N]([3CH3])[3CH2][5C](=[57O])[35NH][10CH2][21c]1[18cH][18cH][18cH][18cH][44n]1)[24c]1[18cH][21c]([18cH][18cH][23c]1[52O][3CH3])[8CH3]"},
  SmilesExpected{"C12(C)C3(C(CC1C1CCC4=CC(=O)C=CC4(C)C1(F)C(O)C2)CN(C)O3)C(=O)COC(=O)C CHEMBL441963", 2.331, "[2C]12([1CH3])[4C]3([2CH]([1CH2][2CH]1[2CH]1[1CH2][1CH2][6C]4=[6CH][5C](=[57O])[6CH]=[6CH][2C]4([1CH3])[4C]1([62F])[4CH]([50OH])[1CH2]2)[3CH2][40N]([3CH3])[51O]3)[5C](=[57O])[3CH2][51O][5C](=[57O])[1CH3]"},
  SmilesExpected{"C(=O)(N(NC(=O)C(C)NC(=O)C(C)NC(=O)N1CCNCC1)CC(=O)N)C1C(O1)C(=O)N(CC1=CC=CC=C1)CC1=CC=CC=C1 CHEMBL584157", -1.164, "[5C](=[57O])([40N]([35NH][5C](=[57O])[4CH]([1CH3])[35NH][5C](=[57O])[4CH]([1CH3])[35NH][5C](=[59O])[40N]1[3CH2][3CH2][35NH][3CH2][3CH2]1)[3CH2][5C](=[57O])[34NH2])[4CH]1[4CH]([51O]1)[5C](=[57O])[40N]([10CH2][21c]1[18cH][18cH][18cH][18cH][18cH]1)[10CH2][21c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"C1=C(C=CC(=C1)CN(SN1CCCCC1)N(C(=O)C1=CC(=CC(=C1)C)C)C(C)(C)C)CC CHEMBL2228842", 6.573, "[18cH]1[21c]([18cH][18cH][21c]([18cH]1)[10CH2][40N]([68S][40N]1[3CH2][1CH2][1CH2][1CH2][3CH2]1)[40N]([5C](=[58O])[21c]1[18cH][21c]([18cH][21c]([18cH]1)[8CH3])[8CH3])[4C]([1CH3])([1CH3])[1CH3])[10CH2][1CH3]"},
  SmilesExpected{"C1(=O)C2=C(C(=O)C3=C1C=C(C=C3OCCC(C)C)NC(=O)CN=N#N)C(=CC=C2)OCCC(C)C CHEMBL4541126", 5.561, "[5C]1(=[58O])[21c]2[21c]([5C](=[58O])[21c]3[21c]1[18cH][22c]([18cH][23c]3[52O][3CH2][1CH2][2CH]([1CH3])[1CH3])[37NH][5C](=[57O])[3CH2][39N]=[47N]#[47N])[23c]([18cH][18cH][18cH]2)[52O][3CH2][1CH2][2CH]([1CH3])[1CH3]"},
  // Substantial difference from RDKit due to aromaticity differences. I think we are more correct.
  // SmilesExpected{"S(C1=NC2=C(C3=NC(C(=O)N13)CCC(=O)NCC1=CC=CO1)C=CC=C2)CC(=O)NC1=CC=C(OCC)C=C1 CHEMBL1707125", 5.561, "[5C]1(=[58O])[21c]2[21c]([5C](=[58O])[21c]3[21c]1[18cH][22c]([18cH][23c]3[52O][3CH2][1CH2][2CH]([1CH3])[1CH3])[37NH][5C](=[57O])[3CH2][39N]=[47N]#[47N])[23c]([18cH][18cH][18cH]2)[52O][3CH2][1CH2][2CH]([1CH3])[1CH3]"},
  // SmilesExpected{"[N+]12(CCCCC1)CN1C(=CC(=C(C1=N2)C(O)=O)C)C CHEMBL4764125", 5.561, "[5C]1(=[58O])[21c]2[21c]([5C](=[58O])[21c]3[21c]1[18cH][22c]([18cH][23c]3[52O][3CH2][1CH2][2CH]([1CH3])[1CH3])[37NH][5C](=[57O])[3CH2][39N]=[47N]#[47N])[23c]([18cH][18cH][18cH]2)[52O][3CH2][1CH2][2CH]([1CH3])[1CH3]"},
  SmilesExpected{"N1(=C2C(=N(=O)C=C1)C=CC=C2)=O CHEMBL2104626", 0.107, "[45n]1([19c]2[19c]([45n](=[53O])[18cH][18cH]1)[18cH][18cH][18cH][18cH]2)=[53O]"},
  // RDKit does not classify the Hydrogens on the P-[OH] groups as acidic. That
  // seems wrong. Using our own value here.
  SmilesExpected{"P(=O)(O)(O)C1=CC=CC(=C1)CP(=O)(O)O CHEMBL149333", 0.167, "[67P](=[61O])([50OH])([50OH])[13c]1[18cH][18cH][18cH][21c]([18cH]1)[10CH2][67P](=[61O])([50OH])[50OH]"},
  SmilesExpected{"S(=O)(=O)(N=S(C)CCCC)C1=CC=C(C)C=C1 CHEMBL1533581", 2.916, "[69S](=[54O])(=[54O])([39N]=[69S]([3CH3])[3CH2][1CH2][1CH2][1CH3])[24c]1[18cH][18cH][21c]([8CH3])[18cH][18cH]1"},
  // RDKit gets atom 0 wrong, it assigns C20 when it should be C19. Probably aromaticity.
  SmilesExpected{"C12=C(N=C(NCC3=CC=CC=C3)NC1=NC(=N2)C(C)C)NC1=CC=C(Cl)C=C1 CHEMBL318880", 5.485, "[19c]12[22c]([44n][22c]([37NH][10CH2][21c]3[18cH][18cH][18cH][18cH][18cH]3)[44nH][19c]1[44n][21c]([44n]2)[11CH]([1CH3])[1CH3])[37NH][22c]1[18cH][18cH][15c]([63Cl])[18cH][18cH]1"},
  SmilesExpected{"O=S(=O)(N)C1=CC=C(C(=O)C2=C(N)N=C(S2)NC2=CC=C(C=C2)S(=O)(=O)N)C=C1 CHEMBL2377821", 0.995, "[54O]=[69S](=[54O])([34NH2])[24c]1[18cH][18cH][21c]([5C](=[58O])[21c]2[22c]([36NH2])[44n][22c]([70s]2)[37NH][22c]2[18cH][18cH][24c]([18cH][18cH]2)[69S](=[54O])(=[54O])[34NH2])[18cH][18cH]1"},
  SmilesExpected{"C1(=C(C(F)(F)F)C=CC(=N1)OC1=CC(=CC=C1)C(=N)N)OC1=CC=CC(=C1)C(=N)N CHEMBL50714", 4.253, "[23c]1([21c]([12C]([62F])([62F])[62F])[18cH][18cH][23c]([44n]1)[52O][23c]1[18cH][21c]([18cH][18cH][18cH]1)[5C](=[40NH])[34NH2])[52O][23c]1[18cH][18cH][18cH][21c]([18cH]1)[5C](=[40NH])[34NH2]"},
  // Aromaticity differences prevent this matching RDKit (0.977)
  SmilesExpected{"C1(C)(C)N(=C2C=CC(=CC2=N1=O)COC1=CC=C(C=C1)C=NNC(=S)NCC=C)=O CHEMBL3347314", 1.916, "[4C]1([1CH3])([1CH3])[48N](=[5C]2[6CH]=[6CH][6C](=[6CH][5C]2=[48N]1=[53O])[3CH2][52O][23c]1[18cH][18cH][21c]([18cH][18cH]1)[5CH]=[39N][35NH][5C](=[68S])[35NH][3CH2][6CH]=[6CH2])=[53O]"},
  // Aromaticity differences prevent this matching RDKit (0.689)
  SmilesExpected{"C(CCNCCCN)N1C(=O)C2=CC=C3C4=C(C=CC(=C24)C1=O)C(=O)N(CCCNCCCN)C3=O CHEMBL2079466", -0.2774, "[10CH2]([1CH2][3CH2][35NH][3CH2][1CH2][3CH2][34NH2])[44n]1[25c](=[56O])[19c]2[18cH][18cH][19c]3[19c]4[19c]([18cH][18cH][19c]([19c]24)[25c]1=[56O])[25c](=[56O])[44n]([10CH2][1CH2][3CH2][35NH][3CH2][1CH2][3CH2][34NH2])[25c]3=[56O]"},
  // Aromaticity differences prevent this matching RDKit (0.220)
  SmilesExpected{"C1=CC=CC(=C1)NC(=O)C(=O)NCC1=CC=C(O1)C=C1C(=O)N2C(=NCC2)S1 CHEMBL221991", 1.820, "[18cH]1[18cH][18cH][18cH][22c]([18cH]1)[37NH][5C](=[57O])[5C](=[57O])[35NH][10CH2][21c]1[18cH][18cH][21c]([49o]1)[26CH]=[6C]1[5C](=[57O])[40N]2[5C](=[39N][3CH2][3CH2]2)[68S]1"},
  SmilesExpected{"C1(=S)C=CSS1 CHEMBL368700", 2.539, "[28c]1(=[68S])[18cH][18cH][70s][70s]1"}
));

}  // namespace
