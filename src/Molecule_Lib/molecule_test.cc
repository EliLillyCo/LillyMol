// Tests for Molecule

#include <unordered_map>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"


#include "molecule.h"
#include "path.h"

namespace {

using testing::UnorderedElementsAreArray;
using testing::FloatNear;
using testing::DoubleNear;

class TestSubstructure : public testing::Test
{
  protected:
    IWString _smiles;
    Molecule _m;
};

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

TEST(TestFsid, TestFsid) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C12CC2C1C1C2CC12.C12CC2C1C1C2CC12.C1CC1"));
  std::unordered_map<int, int> fsid;
  for (const Ring* r : m.sssr_rings()) {
    const int f = r->fused_system_identifier();
    fsid[f] += 1;
  }

  EXPECT_EQ(fsid.size(), 5);
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

}  // namespace
