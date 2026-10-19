// Tests for Molecule

#include <fstream>
#include <memory>
#include <span>
#include <string>
#include <vector>
#include <unordered_map>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"


#include "Foundational/data_source/iwstring_data_source.h"

#include "molecule.h"
#include "moleculeio.h"
#include "path.h"

namespace {

using testing::UnorderedElementsAreArray;
using testing::FloatNear;
using testing::DoubleNear;
using testing::ElementsAre;
using testing::UnorderedElementsAre;

class ExtraTextInfoScope {
  private:
    int _read_extra_text_info;
    int _write_extra_text_info;

  public:
    ExtraTextInfoScope() :
        _read_extra_text_info(moleculeio::read_extra_text_info()),
        _write_extra_text_info(moleculeio::write_extra_text_info()) {
    }

    ~ExtraTextInfoScope() {
      moleculeio::set_read_extra_text_info(_read_extra_text_info);
      moleculeio::set_write_extra_text_info(_write_extra_text_info);
    }
};

std::vector<std::string>
ExtraTextInfo(const Molecule& m) {
  std::vector<std::string> result;
  result.reserve(m.number_records_text_info());
  for (int i = 0; i < m.number_records_text_info(); ++i) {
    result.push_back(m.text_info(i).AsString());
  }

  return result;
}

static std::string
WriteSdfFixture(const char* stem) {
  static constexpr char kSdf[] = R"(methane
  LillyMol  0808262D

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0
M  END
>  <ID>
CHEMBL123

$$$$
)";

  const std::string fname = testing::TempDir() + stem;
  std::ofstream output(fname);
  EXPECT_TRUE(output.good());
  output << kSdf;
  return fname;
}

TEST(MoleculeExtraTextInfo, ReaderLocalSdfTextInfo) {
  ExtraTextInfoScope scope;
  moleculeio::set_read_extra_text_info(0);

  const std::string keep_fname = WriteSdfFixture("/keep_extra_text_info.sdf");
  iwstring_data_source keep_input(keep_fname.c_str());
  ASSERT_TRUE(keep_input.good());

  Molecule_Read_Options keep_options;
  keep_options.mdl_file_supporting_material().set_read_extra_text_info(1);

  Molecule keep;
  ASSERT_TRUE(keep.read_molecule_ds(keep_input, FILE_TYPE_SDF, keep_options));
  EXPECT_THAT(ExtraTextInfo(keep), ElementsAre(">  <ID>", "CHEMBL123", ""));

  const std::string discard_fname = WriteSdfFixture("/discard_extra_text_info.sdf");
  iwstring_data_source discard_input(discard_fname.c_str());
  ASSERT_TRUE(discard_input.good());

  Molecule_Read_Options discard_options;
  discard_options.mdl_file_supporting_material().set_read_extra_text_info(0);

  Molecule discard;
  ASSERT_TRUE(discard.read_molecule_ds(discard_input, FILE_TYPE_SDF, discard_options));
  EXPECT_EQ(discard.number_records_text_info(), 0);
}

TEST(MoleculeExtraTextInfo, SdfTextInfoRoundTripAndManipulation) {
  static constexpr char kSdf[] = R"(methane
  LillyMol  0808262D

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0
M  END
>  <ID>
CHEMBL123

>  <MULTI LINE>
first line
second line

$$$$
)";

  ExtraTextInfoScope scope;
  moleculeio::set_read_extra_text_info(1);

  const std::string fname = testing::TempDir() + "/extra_text_info.sdf";
  {
    std::ofstream output(fname);
    ASSERT_TRUE(output.good());
    output << kSdf;
  }

  iwstring_data_source input(fname.c_str());
  ASSERT_TRUE(input.good());

  Molecule m;
  ASSERT_TRUE(m.read_molecule_ds(input, FILE_TYPE_SDF));
  EXPECT_THAT(ExtraTextInfo(m), ElementsAre(
      ">  <ID>",
      "CHEMBL123",
      "",
      ">  <MULTI LINE>",
      "first line",
      "second line",
      ""));

  IWString written;
  EXPECT_GT(m.write_extra_text_info(written), 0);
  EXPECT_EQ(written.AsString(),
            ">  <ID>\nCHEMBL123\n\n>  <MULTI LINE>\nfirst line\nsecond line\n\n");

  Molecule copied(m);
  EXPECT_THAT(ExtraTextInfo(copied), ElementsAre(
      ">  <ID>",
      "CHEMBL123",
      "",
      ">  <MULTI LINE>",
      "first line",
      "second line",
      ""));

  Molecule assigned;
  assigned = m;
  EXPECT_THAT(ExtraTextInfo(assigned), ElementsAre(
      ">  <ID>",
      "CHEMBL123",
      "",
      ">  <MULTI LINE>",
      "first line",
      "second line",
      ""));

  Molecule moved;
  moved = std::move(assigned);
  EXPECT_THAT(ExtraTextInfo(moved), ElementsAre(
      ">  <ID>",
      "CHEMBL123",
      "",
      ">  <MULTI LINE>",
      "first line",
      "second line",
      ""));

  moved.discard_extra_text_info();
  EXPECT_EQ(moved.number_records_text_info(), 0);

  EXPECT_EQ(moved.add_extra_text_info("added text"), 1);
  EXPECT_THAT(ExtraTextInfo(moved), ElementsAre("added text"));
}

class TestSubstructure : public testing::Test
{
  protected:
    IWString _smiles;
    Molecule _m;
};

TEST_F(TestSubstructure, TestSpaceAtEndOfSmiles) {
  _smiles = "C ";
  EXPECT_TRUE(_m.build_from_smiles(_smiles));
}

TEST_F(TestSubstructure, TestSwapAtoms1) {
  _smiles = "CC=N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  const IWString initial_usmi = _m.unique_smiles();
  _m.swap_atoms(0, 2);
  EXPECT_EQ(initial_usmi, _m.unique_smiles());
}


TEST_F(TestSubstructure, TestSwapAtomsChiral) {
  _smiles = "F[C@H](C)N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  const IWString initial_usmi = _m.unique_smiles();
  _m.swap_atoms(0, 2);
  EXPECT_EQ(initial_usmi, _m.unique_smiles());
  _m.swap_atoms(0, 3);
  EXPECT_EQ(initial_usmi, _m.unique_smiles());
  _m.swap_atoms(1, 3);
  EXPECT_EQ(initial_usmi, _m.unique_smiles());
  _m.swap_atoms(2, 3);
  EXPECT_EQ(initial_usmi, _m.unique_smiles());
}

TEST_F(TestSubstructure, TestRenumberAtoms) {
  _smiles = "CNO";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  // new_number[i] is the new atom number for atom i in the starting molecule.
  static constexpr int kNewNumber[] = {2, 0, 1};
  EXPECT_EQ(_m.renumber_atoms(kNewNumber), 1);
  EXPECT_EQ(_m.atomic_number(0), 7);
  EXPECT_EQ(_m.atomic_number(1), 8);
  EXPECT_EQ(_m.atomic_number(2), 6);

  static constexpr int kDuplicate[] = {0, 0, 1};
  EXPECT_EQ(_m.renumber_atoms(kDuplicate), 0);

  static constexpr int kOutOfRange[] = {0, 1, 3};
  EXPECT_EQ(_m.renumber_atoms(kOutOfRange), 0);
}

// Tests the "ok_atom_number" method which tells if a given atom index is valid for a given molecule

TEST_F(TestSubstructure, TestAtomNumberIndex1) {
  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_TRUE(_m.ok_atom_number(0));
  EXPECT_FALSE(_m.ok_atom_number(-1));
  EXPECT_FALSE(_m.ok_atom_number(1));
  EXPECT_FALSE(_m.ok_atom_number(2));
  EXPECT_FALSE(_m.ok_atom_number(3));
}


TEST_F(TestSubstructure, TestAtomNumberIndex2) {
  _smiles = "CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_TRUE(_m.ok_atom_number(0));
  EXPECT_TRUE(_m.ok_atom_number(1));
  EXPECT_FALSE(_m.ok_atom_number(-1));
  EXPECT_FALSE(_m.ok_atom_number(2));
  EXPECT_FALSE(_m.ok_atom_number(3));
  EXPECT_FALSE(_m.ok_atom_number(4));
}


TEST_F(TestSubstructure, TestAtomNumberIndex3) {
  _smiles = "CCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_TRUE(_m.ok_atom_number(0));
  EXPECT_TRUE(_m.ok_atom_number(1));
  EXPECT_TRUE(_m.ok_atom_number(2));
  EXPECT_FALSE(_m.ok_atom_number(-1));
  EXPECT_FALSE(_m.ok_atom_number(3));
  EXPECT_FALSE(_m.ok_atom_number(4));
}


TEST_F(TestSubstructure, TestAtomNumberIndex4) {
  _smiles = "N#Cc1ccc(-c2coc3cc(O)ccc3c2=O)cc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_TRUE(_m.ok_atom_number(0));
  EXPECT_TRUE(_m.ok_atom_number(1));
  EXPECT_TRUE(_m.ok_atom_number(2));
  EXPECT_TRUE(_m.ok_atom_number(3));
  EXPECT_TRUE(_m.ok_atom_number(4));
  EXPECT_TRUE(_m.ok_atom_number(5));
  EXPECT_TRUE(_m.ok_atom_number(6));
  EXPECT_TRUE(_m.ok_atom_number(7));
  EXPECT_TRUE(_m.ok_atom_number(8));
  EXPECT_TRUE(_m.ok_atom_number(9));
  EXPECT_TRUE(_m.ok_atom_number(10));
  EXPECT_TRUE(_m.ok_atom_number(11));
  EXPECT_TRUE(_m.ok_atom_number(12));
  EXPECT_TRUE(_m.ok_atom_number(13));
  EXPECT_TRUE(_m.ok_atom_number(14));
  EXPECT_TRUE(_m.ok_atom_number(15));
  EXPECT_TRUE(_m.ok_atom_number(16));
  EXPECT_TRUE(_m.ok_atom_number(17));
  EXPECT_TRUE(_m.ok_atom_number(18));
  EXPECT_TRUE(_m.ok_atom_number(19));
  EXPECT_FALSE(_m.ok_atom_number(-1));
  EXPECT_FALSE(_m.ok_atom_number(20));
  EXPECT_FALSE(_m.ok_atom_number(21));
}

// Tests the "ok_2_atoms" method which tells if a given set of two atom indices is valid for a given molecule

TEST_F(TestSubstructure, TestTwoAtomNumberIndex1) {
  _smiles = "CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_TRUE(_m.ok_2_atoms(0,1));
  EXPECT_FALSE(_m.ok_2_atoms(0,0));
  EXPECT_FALSE(_m.ok_2_atoms(0,2));
  EXPECT_FALSE(_m.ok_2_atoms(1,2));
  EXPECT_FALSE(_m.ok_2_atoms(1,1));
}

// Tests the "ok_3_atoms" method which tells if a given set of three atom indicies is valid for a given molecule

TEST_F(TestSubstructure, TestThreeAtomNumberIndex1) {
  _smiles = "CCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_TRUE(_m.ok_3_atoms(0,1,2));
  EXPECT_FALSE(_m.ok_3_atoms(0,0,0));
  EXPECT_FALSE(_m.ok_3_atoms(1,1,1));
  EXPECT_FALSE(_m.ok_3_atoms(2,2,2));
  EXPECT_FALSE(_m.ok_3_atoms(3,3,3));
  EXPECT_FALSE(_m.ok_3_atoms(0,1,4));
}

// Tests the "ok_4_atoms" method which tells if a given set of four atom indicies is valid for a given molecule

TEST_F(TestSubstructure, TestFourAtomNumberIndex1) {
  _smiles = "COCCl";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_TRUE(_m.ok_4_atoms(0,1,2,3));
  EXPECT_FALSE(_m.ok_4_atoms(0,0,0,0));
  EXPECT_FALSE(_m.ok_4_atoms(1,1,1,1));
  EXPECT_FALSE(_m.ok_4_atoms(2,2,2,2));
  EXPECT_FALSE(_m.ok_4_atoms(3,3,3,3));
  EXPECT_FALSE(_m.ok_4_atoms(4,4,4,4));
  EXPECT_FALSE(_m.ok_4_atoms(1,2,3,5));
}

// Parameterized Test for the "nrings" method which determines the
// number of rings in a given molecule.

struct SmilesNrings {
  IWString smiles; // input smiles
  int nrings; // input expected number of rings
};

class TestRings : public testing::TestWithParam<SmilesNrings> {
  protected:
    Molecule _m;
};

TEST_P(TestRings, TestRingCount1) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(_m.nrings(), params.nrings);
}
INSTANTIATE_TEST_SUITE_P(TestRings, TestRings, testing::Values(
  SmilesNrings{"CC=N", 0},
  SmilesNrings{"O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl", 2},
  SmilesNrings{"C1CC1", 1},
  SmilesNrings{"C1CCC1", 1},
  SmilesNrings{"C12CC1C2", 2},
  SmilesNrings{"C123CC1C2C3", 3},
  SmilesNrings{"C1C2CC3CC1CC(C2)C3", 3},
  SmilesNrings{"C1CC1.C1CC1", 2}
));

// Parameterized Test for the "natoms" method that determines the
// number of atoms in a given molecule

struct SmilesNatoms {
  IWString smiles; // input smiles
  int natoms; // input number of expected atoms
};

class TestAtomCount : public testing::TestWithParam<SmilesNatoms> {
  protected:
    Molecule _m;
};

TEST_P(TestAtomCount, TestAtomCount1) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(_m.natoms(), params.natoms);
}
INSTANTIATE_TEST_SUITE_P(TestAtomCount, TestAtomCount, testing::Values(
  SmilesNatoms{"CC=N", 3},
  SmilesNatoms{"C", 1},
  SmilesNatoms{"CC", 2},
  SmilesNatoms{"CCC", 3},
  SmilesNatoms{"COCCl", 4},
  SmilesNatoms{"N#Cc1ccc(-c2coc3cc(O)ccc3c2=O)cc1", 20},
  SmilesNatoms{"C[N+](C)(CCOc1ccccc1)Cc1ccccc1.O=C([O-])c1cc2ccccc2cc1O", 33},
  SmilesNatoms{"CC(C)=CCC/C(C)=C/CC/C(C)=C/CSC[C@H](NC(=O)CCCCCN1CCN(C)CC1)C(=O)NC1CCCC1", 41},
  SmilesNatoms{"CCOc1ccc(Cc2cc([C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c3c(c2Cl)CC(C)(C)O3)cc1", 33},
  SmilesNatoms{"Cc1ccccc1N1CCN(CN2C(=O)CC(Cc3ccccc3)C2=O)CC1", 28}
));

// Parameterized Test for the "molecular_formula" and "isis_like_molecular_formula_dot_between_fragments" methonds that determine the correct molecular formulae for a given molecule

struct SmilesFormulae {
  IWString smiles; // input smiles
  IWString molform; // input expected MF
  IWString isis_molform; // input expected ISISMF
  float mol_weight; // input expected molecular weight
  double ex_mass; // input expected exact mass
};

class TestMolecularFormula : public testing::TestWithParam<SmilesFormulae> {
  protected:
    Molecule _m;
};

TEST_P(TestMolecularFormula, TestMolecularFormula1) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(_m.molecular_formula(), params.molform);
  IWString isis_molform;
  _m.isis_like_molecular_formula_dot_between_fragments(isis_molform);
  EXPECT_EQ(isis_molform, params.isis_molform);
  EXPECT_THAT(_m.molecular_weight(), FloatNear(params.mol_weight, 0.001));
  EXPECT_THAT(_m.exact_mass(), DoubleNear(params.ex_mass, 0.000001));
}
INSTANTIATE_TEST_SUITE_P(TestMolecularFormula, TestMolecularFormula, testing::Values(
  SmilesFormulae{"C.C", "C2H8", "CH4.CH4", 32.08492, 32.06260028},
  SmilesFormulae{"C=O.O", "CO2H4", "CH2O.H2O", 48.04126, 48.0211294},
  SmilesFormulae{"c1ccccc1", "C6H6", "C6H6", 78.11184, 78.04695021},
  SmilesFormulae{"[NH4+].OCC#N", "C2N2OH7", "H4N.C2H3NO", 75.08986, 75.05583788},
  SmilesFormulae{"O.O", "O2H4", "H2O.H2O", 36.03056, 36.0211294},
  SmilesFormulae{"C12C3C4C1C5C2C3C45", "C8H8", "C8H8", 104.1491, 104.0626003},
  SmilesFormulae{"C1C2CC3CC1CC(C2)C3", "C10H16", "C10H16", 136.234, 136.1252006},
  SmilesFormulae{"C1CC12CC2", "C5H8", "C5H8", 68.11702, 68.06260028},
  SmilesFormulae{"S(=O)(=O)(F)C", "CO2SFH3", "CH3FO2S", 98.09673, 97.98377828},
  SmilesFormulae{"O1CNCC1", "C3NOH7", "C3H7NO", 73.09382, 73.05276388},
  SmilesFormulae{"O(C)CCN", "C3NOH9", "C3H9NO", 75.1097, 75.06841395},
  SmilesFormulae{"C(C)(C)[N+]#[C-]", "C4NH7", "C4H7N", 69.10512, 69.05784925},
  SmilesFormulae{"C(=C)(C)C=O", "C4OH6", "C4H6O", 70.08984, 70.04186484},
  SmilesFormulae{"C(=O)(O)C=O", "C2O3H2", "C2H2O3", 74.03548, 74.00039396},
  SmilesFormulae{"NCC=CCl", "C3NClH6", "C3H6ClN", 91.53918, 91.01887693},
  SmilesFormulae{"NCC#CBr", "C3NBrH4", "C3H4BrN", 133.9746, 132.9527102},
  SmilesFormulae{"CS(O)(=O)=O", "CO3SH4", "CH4O3S", 96.10567, 95.98811473},
  SmilesFormulae{"NCC=C=C", "C4NH7", "C4H7N", 69.10512, 69.05784925}
));

// Parameterized Test for the "add_bond" and "are_bonded" methods that
// determine if a bond can be added between two given atoms in a given
// molecule.

struct SmilesAtomNumbers {
  IWString smiles; // input smiles
  int atom_1; // input for a given atom
  int atom_2; // input for a given atom
};

class TestAddBond : public testing::TestWithParam<SmilesAtomNumbers> {
  protected:
    Molecule _m;
};

TEST_P(TestAddBond, TestAddBond1) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_TRUE(_m.add_bond(params.atom_1, params.atom_2, SINGLE_BOND));
  EXPECT_TRUE(_m.are_bonded(params.atom_1, params.atom_2));
}
INSTANTIATE_TEST_SUITE_P(TestAddBond, TestAddBond, testing::Values(
  SmilesAtomNumbers{"CCC", 0, 2},
  SmilesAtomNumbers{"O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl", 0, 17},
  SmilesAtomNumbers{"COCCl", 1, 3},
  SmilesAtomNumbers{"N#Cc1ccc(-c2coc3cc(O)ccc3c2=O)cc1", 8, 15},
  SmilesAtomNumbers{"C[N+](C)(CCOc1ccccc1)Cc1ccccc1.O=C([O-])c1cc2ccccc2cc1O", 12, 31},
  SmilesAtomNumbers{"CC(C)=CCC/C(C)=C/CC/C(C)=C/CSC[C@H](NC(=O)CCCCCN1CCN(C)CC1)C(=O)NC1CCCC1", 0, 34},
  SmilesAtomNumbers{"CCOc1ccc(Cc2cc([C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c3c(c2Cl)CC(C)(C)O3)cc1", 0, 27},
  SmilesAtomNumbers{"Cc1ccccc1N1CCN(CN2C(=O)CC(Cc3ccccc3)C2=O)CC1", 4, 17}
));

// Parameterized Test for the "connections" methods that determines
// the atoms that any given atom is connected to in a given molecule
// onetest if for the connections method that returns an array of the
// connected atoms, and one that takes an array and fills it with the
// connected atoms

struct SmilesAtomNeighbors {
  IWString smiles; // input smiles
  int atom_1; // input for a given atom
  int nconnections; // input for the number of expected connections for atom_1
  Set_of_Atoms conn; //array that will hold the set of connected atoms
};

class TestAtomNeighbors : public testing::TestWithParam<SmilesAtomNeighbors> {
  protected:
    Molecule _m;
};

TEST_P(TestAtomNeighbors, TestAtomNeighbors1) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  Set_of_Atoms conn;
  EXPECT_EQ(_m.connections(params.atom_1, conn), params.nconnections);
  EXPECT_THAT(conn, UnorderedElementsAreArray(params.conn));
}
INSTANTIATE_TEST_SUITE_P(TestAtomNeighbors, TestAtomNeighbors, testing::Values(
  SmilesAtomNeighbors{"CCC", 0, 1, {1}},
  SmilesAtomNeighbors{"C1CC1", 0, 2, {1,2}},
  SmilesAtomNeighbors{"C1CC1", 1, 2, {0,2}},
  SmilesAtomNeighbors{"C1CC1", 2, 2, {0,1}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 0, 3, {1,3,5}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 1, 3, {6,2,0}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 2, 3, {7,1,3}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 3, 3, {2,0,4}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 4, 3, {3,5,7}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 5, 3, {0,6,4}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 6, 3, {1,5,7}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 7, 3, {6,2,4}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 0, 2, {1,5}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 1, 3, {0,2,8}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 2, 2, {1,3}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 3, 3, {2,4,9}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 4, 2, {3,5}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 5, 3, {0,4,6}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 6, 2, {7,5}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 7, 3, {6,8,9}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 8, 2, {1,7}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 9, 2, {3,7}},
  SmilesAtomNeighbors{"C1CC12CC2", 0, 2, {1,2}},
  SmilesAtomNeighbors{"C1CC12CC2", 1, 2, {0,2}},
  SmilesAtomNeighbors{"C1CC12CC2", 2, 4, {1,0,4,3}},
  SmilesAtomNeighbors{"C1CC12CC2", 3, 2, {4,2}},
  SmilesAtomNeighbors{"C1CC12CC2", 4, 2, {3,2}}
));

TEST_P(TestAtomNeighbors, TestAtomNeighbors2) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  Set_of_Atoms nbrs;
  EXPECT_EQ(_m.connections(params.atom_1, nbrs), params.nconnections);
  EXPECT_THAT(params.conn, UnorderedElementsAreArray(nbrs));
}
INSTANTIATE_TEST_SUITE_P(TestAtomNeighbors3, TestAtomNeighbors, testing::Values(
  SmilesAtomNeighbors{"CCC", 0, 1, {1}},
  SmilesAtomNeighbors{"C1CC1", 0, 2, {1,2}},
  SmilesAtomNeighbors{"C1CC1", 1, 2, {0,2}},
  SmilesAtomNeighbors{"C1CC1", 2, 2, {0,1}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 0, 3, {1,3,5}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 1, 3, {6,2,0}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 2, 3, {7,1,3}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 3, 3, {2,0,4}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 4, 3, {3,5,7}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 5, 3, {0,6,4}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 6, 3, {1,5,7}},
  SmilesAtomNeighbors{"C12C3C4C1C5C2C3C45", 7, 3, {6,2,4}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 0, 2, {1,5}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 1, 3, {0,2,8}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 2, 2, {1,3}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 3, 3, {2,4,9}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 4, 2, {3,5}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 5, 3, {0,4,6}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 6, 2, {7,5}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 7, 3, {6,8,9}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 8, 2, {1,7}},
  SmilesAtomNeighbors{"C1C2CC3CC1CC(C2)C3", 9, 2, {3,7}},
  SmilesAtomNeighbors{"C1CC12CC2", 0, 2, {1,2}},
  SmilesAtomNeighbors{"C1CC12CC2", 1, 2, {0,2}},
  SmilesAtomNeighbors{"C1CC12CC2", 2, 4, {1,0,4,3}},
  SmilesAtomNeighbors{"C1CC12CC2", 3, 2, {4,2}},
  SmilesAtomNeighbors{"C1CC12CC2", 4, 2, {3,2}}
));

// Parameterized test for the "bonds_between" method that determines
// the number of bonds between two given atoms

struct SmilesAtomsBonds {
  IWString smiles; // input smiles
  int atom_1; // input for a given atom
  int atom_2; // input for a given atom
  int nbonds; // input for the number of expected bonds between atom 1 and atom 2
};

class TestBondsBetween : public testing::TestWithParam<SmilesAtomsBonds> {
  protected:
    Molecule _m;
};

TEST_P(TestBondsBetween, TestBondsBetween1) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(_m.bonds_between(params.atom_1, params.atom_2), params.nbonds);
}
INSTANTIATE_TEST_SUITE_P(TestBondsBetween, TestBondsBetween, testing::Values(
  SmilesAtomsBonds{"CC", 0, 1, 1},
  SmilesAtomsBonds{"CCC", 0, 1, 1},
  SmilesAtomsBonds{"CCC", 0, 2, 2}
));

TEST(TestAllAtomsBetween, DifferentiatesMultipleShortestPaths) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1CCC1"));

  Set_of_Atoms one_path;
  EXPECT_EQ(m.atoms_between(0, 2, one_path), 2);
  ASSERT_EQ(one_path.number_elements(), 1);

  Set_of_Atoms all_paths;
  EXPECT_EQ(m.AllAtomsBetween(0, 2, all_paths), 2);
  EXPECT_THAT(all_paths, UnorderedElementsAre(1, 3));
}

TEST(TestAllAtomsBetween, BreadthFirstOrder) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1CCCCC1"));

  Set_of_Atoms all_paths;
  EXPECT_EQ(m.AllAtomsBetween(0, 3, all_paths), 4);
  EXPECT_THAT(all_paths, ElementsAre(1, 5, 2, 4));
}

TEST(TestAllAtomsBetween, AdjacentAtomsReturnEmptySet) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCC"));

  Set_of_Atoms all_paths;
  EXPECT_EQ(m.AllAtomsBetween(0, 1, all_paths), 0);
  EXPECT_TRUE(all_paths.empty());
}

TEST(TestFsid, TestFsid) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C12CC2C1C1C2CC12.C12CC2C1C1C2CC12.C1CC1"));
  std::unordered_map<int, int> fsid;
  for (const Ring* r : m.sssr_rings()) {
    const int f = r->fused_system_identifier();
    fsid[f] += 1;
  }

  EXPECT_EQ(fsid.size(), 5u);
}

TEST(TestNameConcat, TestNameConcat1) {
  Molecule m;
  m << "foo";
  EXPECT_EQ(m.name(), "foo");
  m << "bar";
  EXPECT_EQ(m.name(), "foobar");
}

TEST(TestGetSetCoordinates, TestGetSetCoordinates1) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C{{-0.0178,1.4608,0.0101}}N{{0.0021,-0.0041,0.002}}1C{{0.0222,-0.6766,1.1648}}=C{{0.037,-2.0379,1.1987}}C{{0.0412,-2.7759,0.0216}}=C{{0.0254,-2.0993,-1.2002}}2C{{0.0056,-0.6775,-1.1921}}1=C{{-0.0067,-0.2866,-2.5234}}C{{0.0055,-1.4481,-3.2818}}=N{{0.0243,-2.4979,-2.4818}}2"));
  std::unique_ptr<float[]> initial_coords = m.GetCoordinates();
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    m.setxyz(i, 0.0, 0.1, 0.2);
  }

  m.SetCoordinates(initial_coords.get());
  for (int i = 0; i < matoms; ++i) {
    EXPECT_EQ(m.x(i), initial_coords[3 * i]);
    EXPECT_EQ(m.y(i), initial_coords[3 * i + 1]);
    EXPECT_EQ(m.z(i), initial_coords[3 * i + 2]);
  }
}

struct SmilesNsys {
  IWString smiles;
  int number_ring_systems;
};

class TestNumberRingSystems : public testing::TestWithParam<SmilesNsys> {
  protected:
    Molecule _m;
};

TEST_P(TestNumberRingSystems, TestNumberRingSystems) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(_m.number_ring_systems(), params.number_ring_systems) <<
    params.smiles << " expected " << params.number_ring_systems <<
    " got " << _m.number_ring_systems();
}
INSTANTIATE_TEST_SUITE_P(TestNumberRingSystems, TestNumberRingSystems, testing::Values(
  SmilesNsys{"CC", 0},
  SmilesNsys{"C1CC1", 1},
  SmilesNsys{"C1CC1C1CC1", 2},
  SmilesNsys{"C12CC1C2", 1},
  SmilesNsys{"C1CC1C2CC2", 2}
));

struct SmilesAtomValue {
  IWString smiles;
  std::vector<atom_number_t> atom;
  std::vector<int> expected;
};

class TestSaturation : public testing::TestWithParam<SmilesAtomValue> {
  protected:
    Molecule _m;
};
TEST_P(TestSaturation, TestSaturation) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  for (uint32_t i = 0; i < params.atom.size(); ++i) {
    EXPECT_EQ(_m.saturated(params.atom[i]), params.expected[i]);
  }
}

INSTANTIATE_TEST_SUITE_P(TestSaturation, TestSaturation, testing::Values(
  SmilesAtomValue{"CCC", {0, 1, 2}, {1, 1, 1}},
  SmilesAtomValue{"CC=C", {0, 1, 2}, {1, 0, 0}},
  SmilesAtomValue{"CNC(=O)C", {0, 1, 2, 3, 4}, {1, 1, 0, 0, 1}},
  SmilesAtomValue{"CNS(=O)(=O)C", {0, 1, 2, 3, 4, 5}, {1, 1, 0, 0, 0, 1}}
));

TEST(TestAtomOrderInSmiles, Test1) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("CN"));
  const int matoms = m.natoms();

  EXPECT_EQ(m.smiles(), "CN");

  const resizable_array<atom_number_t>& aoit1 = m.atom_order_in_smiles();
  for (int i = 0; i < matoms; ++i) {
    EXPECT_EQ(aoit1[i], i);
  }

  EXPECT_EQ(m.smiles_starting_with_atom(1), "NC");
  const resizable_array<atom_number_t>& aoit2 = m.atom_order_in_smiles();
  EXPECT_EQ(aoit2[0], 1);
  EXPECT_EQ(aoit2[1], 0);
}

struct SmilesAndSpinach {
  IWString smiles;
  isotope_t isotope;
  std::vector<int> expected_spinach;
  // If `isotope` is specified.
  IWString expected_smiles;
};

class TestSpinach : public testing::TestWithParam<SmilesAndSpinach> {
  protected:
    Molecule _m;
};
TEST_P(TestSpinach, TestSpinach) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));

  const int matoms = _m.natoms();

  std::unique_ptr<int[]> spinach = std::make_unique<int[]>(matoms);
  _m.IdentifySpinachLabel(spinach.get(), params.isotope);
  const std::span<int> gotspan(spinach.get(), spinach.get() + matoms);
  EXPECT_THAT(gotspan, testing::ElementsAreArray(params.expected_spinach));
  EXPECT_EQ(_m.smiles(), params.expected_smiles) << " got " << _m.smiles();
}
INSTANTIATE_TEST_SUITE_P(TestSpinach, TestSpinach, testing::Values(
  SmilesAndSpinach{"CC", 1, {1, 1}, "CC"},
  SmilesAndSpinach{"C1CC1", 1, {0, 0, 0}, "C1CC1"},
  SmilesAndSpinach{"C1CC1=O", 1, {0, 0, 0, 0}, "C1CC1=O"},
  SmilesAndSpinach{"CC1CC1=O", 1, {1, 0, 0, 0, 0}, "[1CH3][1CH]1CC1=O"} 
));

TEST(BondTypeChange, RemoveObsoleteAlkeneDirectionFlags) {
  for (const char* smiles : {"C/C=C/C", "C/C=C\\C"}) {
    Molecule mol;
    ASSERT_TRUE(mol.build_from_smiles(smiles));
    ASSERT_TRUE(mol.set_bond_type_between_atoms(1, 2, SINGLE_BOND));
    Molecule expected;
    ASSERT_TRUE(expected.build_from_smiles("CCCC"));
    EXPECT_EQ(mol.unique_smiles(), expected.unique_smiles());
    EXPECT_FALSE(mol.bond_between_atoms(0, 1)->is_directional());
    EXPECT_FALSE(mol.bond_between_atoms(2, 3)->is_directional());
  }
}

TEST(BondTypeChange, PreserveNeighbouringAlkeneGeometry) {
  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("C/C=C/C=C/C"));
  ASSERT_TRUE(mol.set_bond_type_between_atoms(1, 2, SINGLE_BOND));
  Molecule expected;
  ASSERT_TRUE(expected.build_from_smiles("CCC/C=C/C"));
  EXPECT_EQ(mol.unique_smiles(), expected.unique_smiles());
  EXPECT_FALSE(mol.bond_between_atoms(0, 1)->is_directional());
  EXPECT_TRUE(mol.bond_between_atoms(2, 3)->is_directional());
  EXPECT_TRUE(mol.bond_between_atoms(3, 4)->part_of_cis_trans_grouping());
}

TEST(DeleteAllFragmentsExcept, KeepsSelectedFragmentAndReturnsRemovedAtomCount) {
  for (int frag = 0; frag < 3; ++frag) {
    Molecule mol;
    ASSERT_TRUE(mol.build_from_smiles("CC.O.N mixture"));
    const char* expected[] = {"CC", "O", "N"};
    EXPECT_EQ(mol.delete_all_fragments_except(frag), frag == 0 ? 2 : 3);
    EXPECT_EQ(mol.unique_smiles(), expected[frag]);
    EXPECT_EQ(mol.number_fragments(), 1);
    EXPECT_EQ(mol.name(), "mixture");
  }
}

TEST(DeleteAllFragmentsExcept, RejectsInvalidInputWithoutChangingMolecule) {
  for (const char* smiles : {"", "CC", "CC.O.N"}) {
    Molecule mol;
    if (*smiles) ASSERT_TRUE(mol.build_from_smiles(smiles));
    mol.set_name("invalid_input");
    const IWString before = mol.unique_smiles();
    const int natoms = mol.natoms();
    const int nf = mol.number_fragments();
    for (int frag : {-1, nf}) {
      EXPECT_EQ(mol.delete_all_fragments_except(frag), 0);
      EXPECT_EQ(mol.unique_smiles(), before);
      EXPECT_EQ(mol.natoms(), natoms);
    }
    if (nf <= 1) {
      EXPECT_EQ(mol.delete_all_fragments_except(0), 0);
      EXPECT_EQ(mol.unique_smiles(), before);
    }
  }
}

}  // namespace
