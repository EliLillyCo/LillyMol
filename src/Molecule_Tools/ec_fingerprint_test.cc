#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "gtest/gtest.h"

#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"

#include "ec_fingerprint.h"

namespace ec_fingerprint {
namespace {

using std::cerr;

class TestECFingerPrint : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    IWString _smiles;
    Molecule _m;
    int _matoms;
    Atom_Typing_Specification _atom_typing;
    ProduceFingerprint _fp;
    ECFingerprint _ecfp;
    IWString _fingerprint;
};

void
TestECFingerPrint::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

TEST_F(TestECFingerPrint, EmptyMolecule)
{
  atom_type_t * atype = new atom_type_t[1]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  EXPECT_EQ(_ecfp.Fingerprint(_m, nullptr, atype, _fp), 0);
  const Sparse_Fingerprint_Creator& sfc = _fp.sfc();
  EXPECT_EQ(sfc.nbits(), 0l);
}

TEST_F(TestECFingerPrint, TestSingleAtom)
{
  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  ASSERT_EQ(_matoms, 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);

  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_ecfp.Fingerprint(_m, nullptr, atype, _fp));

  const Sparse_Fingerprint_Creator& sfc = _fp.sfc();
  EXPECT_EQ(sfc.nbits(), 1l);

  sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  cerr << _fingerprint << "\n";
//EXPECT_EQ(_fingerprint, "..+u..2.....2");
  EXPECT_EQ(_fingerprint, "...xPE2.....2");
}

TEST_F(TestECFingerPrint, TestOneAtomExcluded)
{
  _smiles = "CN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  const int matoms = _m.natoms();

  int * include_atom = new int[matoms]; std::unique_ptr<int[]> free_include_atom(include_atom);
  include_atom[0] = 1;
  include_atom[1] = 0;

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);

  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_ecfp.Fingerprint(_m, include_atom, atype, _fp));

  const Sparse_Fingerprint_Creator& sfc = _fp.sfc();

  EXPECT_EQ(sfc.nbits(), 1l);

  sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  cerr << "fp is " << _fingerprint << '\n';
//EXPECT_EQ(_fingerprint, "..+u..2.....2");
  EXPECT_EQ(_fingerprint, "...xPE2.....2");
}

TEST_F(TestECFingerPrint, TestSingleAtomCType) 
{
  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 1);

  ASSERT_TRUE(_atom_typing.build("C"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);

  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_ecfp.Fingerprint(_m, nullptr, atype, _fp));

  const Sparse_Fingerprint_Creator& sfc = _fp.sfc();

  EXPECT_EQ(sfc.nbits(), 1l);

  sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  cerr << _fingerprint << '\n';
  EXPECT_EQ(_fingerprint, "...nF.2.....2");
}

TEST_F(TestECFingerPrint, TestFragments) 
{
  _smiles = "C.C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 2);
  EXPECT_EQ(_m.number_fragments(), 2);

  ASSERT_TRUE(_atom_typing.build("C"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_ecfp.Fingerprint(_m, nullptr, atype, _fp));
  const Sparse_Fingerprint_Creator& sfc = _fp.sfc();

  EXPECT_EQ(sfc.nbits(), 1l);

  sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  cerr << "fp " << _fingerprint << '\n';
  EXPECT_EQ(_fingerprint, "...nF.6.....2");
}

TEST_F(TestECFingerPrint, TestBenzene) 
{
  _smiles = "c1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 6);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("C"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_ecfp.Fingerprint(_m, nullptr, atype, _fp));

  const Sparse_Fingerprint_Creator& sfc = _fp.sfc();

  EXPECT_EQ(sfc.nbits(), 4l);

  sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "...mF0LDCGtn7azuyCRLR.M4+UM.2");
}

TEST_F(TestECFingerPrint, TestButane) 
{
  _smiles = "CCCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 4);
  EXPECT_EQ(_m.number_fragments(), 1);
  _m.compute_aromaticity();

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_ecfp.Fingerprint(_m, nullptr, atype, _fp));

  const Sparse_Fingerprint_Creator& sfc = _fp.sfc();
  EXPECT_EQ(sfc.nbits(), 6l);

  sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "...xPJ9gSrtKWUubfFDTsEE0.U9KajCmrbDpT.60....2");
}

TEST_F(TestECFingerPrint, TestMinRadiusTooLong)
{
  _smiles = "CCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 3);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _ecfp.set_min_radius(4);
  ASSERT_EQ(_ecfp.Fingerprint(_m, nullptr, atype, _fp), 1);  // even though no bits found

  const Sparse_Fingerprint_Creator& sfc = _fp.sfc();

  EXPECT_EQ(sfc.nbits(), 0l);
}

TEST_F(TestECFingerPrint, TestMinMaxSame)
{
  _smiles = "CCCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 4);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _ecfp.set_min_radius(1);
  _ecfp.set_max_radius(1);
  ASSERT_TRUE(_ecfp.Fingerprint(_m, nullptr, atype, _fp));

  const Sparse_Fingerprint_Creator& sfc = _fp.sfc();

  EXPECT_EQ(sfc.nbits(), 2l);

  sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "JccCduoHry20.U..3");
}

TEST_F(TestECFingerPrint, TestMissingCentreAtom)
{
  _smiles = "CCC(CC)CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  int * include_atom = new int[_matoms]; std::unique_ptr<int[]> free_include_atom(include_atom);

  std::fill_n(include_atom, _matoms, 1);
  include_atom[2] = 0;

  EXPECT_EQ(_matoms, 7);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_ecfp.Fingerprint(_m, include_atom, atype, _fp));

  const Sparse_Fingerprint_Creator& sfc = _fp.sfc();

  EXPECT_EQ(sfc.nbits(), 2l);

  sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "...xPJO81eQ4+U..3");
}

TEST_F(TestECFingerPrint, TestCoverage)
{
  _smiles = "C1(=CN(C)C(=C1)C)C1=CSC(=N1)NC(N)=N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 16);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:ACHY"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _ecfp.set_max_radius(2);
  AtomMapCoverage coverage;
  ASSERT_TRUE(_ecfp.Fingerprint(_m, nullptr, atype, coverage));

  coverage.FingerprintingComplete(_m);

  EXPECT_EQ(_m.smiles(), "[C:6]1(=[CH:4][N:6]([CH3:2])[C:6](=[CH:4]1)[CH3:2])[C:6]1=[CH:4][S:4][C:6](=[N:4]1)[NH:4][C:6]([NH2:2])=[NH:2]");
}

TEST_F(TestECFingerPrint, TestCoverage2)
{
  _smiles = "CN1C=CN=C1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 6);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:CYP"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _ecfp.set_max_radius(2);
  AtomMapCoverage coverage;
  ASSERT_TRUE(_ecfp.Fingerprint(_m, nullptr, atype, coverage));

  coverage.FingerprintingComplete(_m);

  EXPECT_EQ(_m.smiles(), "[CH3:2][N:6]1[CH:4]=[CH:4][N:4]=[CH:4]1");
}

const std::vector<std::string> aspirin = {
  "N1=CN(C)C2=C1N(C)C(=O)N(C)C2=O",
  "N1(C2=C(N(C(=O)N(C)C2=O)C)N=C1)C",
  "N1=CN(C)C2=C1N(C)C(=O)N(C2=O)C",
  "C1(=O)N(C2=C(N(C=N2)C)C(=O)N1C)C",
  "O=C1N(C)C2=C(N(C=N2)C)C(=O)N1C",
  "C12=C(C(=O)N(C)C(=O)N1C)N(C=N2)C",
  "C12=C(N=CN1C)N(C)C(=O)N(C)C2=O",
  "N1(C)C=NC2=C1C(=O)N(C)C(=O)N2C",
  "C12=C(N(C=N1)C)C(=O)N(C)C(=O)N2C",
  "C1=NC2=C(N1C)C(=O)N(C(=O)N2C)C",
  "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
  "C1=NC2=C(C(=O)N(C(=O)N2C)C)N1C",
  "N1(C)C2=C(N(C=N2)C)C(=O)N(C)C1=O",
  "O=C1N(C(=O)C2=C(N1C)N=CN2C)C",
  "C1(=O)N(C)C2=C(C(=O)N1C)N(C=N2)C",
  "C1=NC2=C(C(=O)N(C)C(=O)N2C)N1C",
  "N1(C)C(=O)N(C2=C(N(C=N2)C)C1=O)C",
  "C1(=O)N(C)C(=O)C2=C(N=CN2C)N1C",
  "C1(=O)N(C)C(=O)N(C)C2=C1N(C=N2)C",
  "C12=C(N=CN1C)N(C)C(=O)N(C2=O)C",
  "N1(C)C(=O)N(C)C2=C(N(C=N2)C)C1=O",
  "N1(C)C(=O)N(C2=C(C1=O)N(C=N2)C)C",
  "O=C1C2=C(N=CN2C)N(C(=O)N1C)C",
  "CN1C2=C(N(C=N2)C)C(=O)N(C)C1=O",
  "CN1C(=O)N(C)C2=C(C1=O)N(C)C=N2",
  "N1(C2=C(N=C1)N(C)C(=O)N(C2=O)C)C",
  "CN1C(=O)N(C)C2=C(N(C=N2)C)C1=O",
  "CN1C(=O)N(C)C(=O)C2=C1N=CN2C",
  "CN1C(=O)N(C2=C(N(C)C=N2)C1=O)C",
  "C1(=O)N(C(=O)N(C2=C1N(C=N2)C)C)C",
  "N1(C)C(=O)C2=C(N(C)C1=O)N=CN2C",
  "CN1C2=C(N(C(=O)N(C)C2=O)C)N=C1",
  "O=C1N(C)C(=O)N(C2=C1N(C=N2)C)C",
  "N1(C(=O)N(C)C(=O)C2=C1N=CN2C)C",
  "N1(C(=O)C2=C(N=CN2C)N(C)C1=O)C",
  "N1(C)C2=C(N(C)C=N2)C(=O)N(C1=O)C",
  "N1(C2=C(N=C1)N(C)C(=O)N(C)C2=O)C",
  "N1(C)C(=O)C2=C(N=CN2C)N(C1=O)C",
  "C1=NC2=C(N1C)C(=O)N(C)C(=O)N2C",
  "C1(=O)N(C)C(=O)N(C2=C1N(C)C=N2)C",
  "O=C1N(C)C2=C(C(=O)N1C)N(C=N2)C",
  "C1(=O)N(C)C(=O)N(C)C2=C1N(C)C=N2",
  "N1(C)C(=O)N(C)C2=C(N(C)C=N2)C1=O",
  "C1(=O)N(C2=C(C(=O)N1C)N(C=N2)C)C",
  "N1(C(=O)C2=C(N(C1=O)C)N=CN2C)C",
  "C1(=O)N(C(=O)N(C)C2=C1N(C)C=N2)C",
  "O=C1N(C(=O)C2=C(N=CN2C)N1C)C",
  "N1(C(=O)C2=C(N=CN2C)N(C1=O)C)C",
  "N1(C(=O)N(C)C2=C(C1=O)N(C)C=N2)C",
  "C12=C(N(C)C=N1)C(=O)N(C(=O)N2C)C",
  "N1=CN(C2=C1N(C(=O)N(C2=O)C)C)C",
  "C1(=O)N(C)C(=O)C2=C(N1C)N=CN2C",
  "N1(C(=O)N(C2=C(C1=O)N(C=N2)C)C)C",
  "O=C1N(C)C(=O)N(C)C2=C1N(C)C=N2",
  "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
  "C12=C(N=CN1C)N(C(=O)N(C2=O)C)C",
  "C12=C(C(=O)N(C(=O)N1C)C)N(C)C=N2",
  "N1(C)C=NC2=C1C(=O)N(C(=O)N2C)C",
  "N1(C)C(=O)N(C(=O)C2=C1N=CN2C)C",
  "O=C1N(C(=O)N(C)C2=C1N(C=N2)C)C",
  "CN1C(=O)C2=C(N=CN2C)N(C1=O)C",
  "CN1C(=O)C2=C(N(C1=O)C)N=CN2C",
  "C1(=O)C2=C(N(C)C(=O)N1C)N=CN2C",
  "O=C1N(C2=C(C(=O)N1C)N(C)C=N2)C",
  "O=C1C2=C(N(C)C(=O)N1C)N=CN2C",
  "C12=C(N(C)C(=O)N(C)C1=O)N=CN2C",
  "O=C1C2=C(N=CN2C)N(C)C(=O)N1C",
  "O=C1C2=C(N(C(=O)N1C)C)N=CN2C",
  "N1(C(=O)C2=C(N(C)C1=O)N=CN2C)C",
  "O=C1N(C)C(=O)C2=C(N1C)N=CN2C",
  "C12=C(N(C(=O)N(C1=O)C)C)N=CN2C",
  "C1(=O)N(C(=O)C2=C(N1C)N=CN2C)C",
  "N1(C2=C(N(C)C=N2)C(=O)N(C1=O)C)C",
  "CN1C2=C(N=C1)N(C)C(=O)N(C)C2=O",
  "N1(C)C(=O)N(C)C2=C(C1=O)N(C=N2)C",
  "N1(C=NC2=C1C(=O)N(C)C(=O)N2C)C",
  "N1(C=NC2=C1C(=O)N(C(=O)N2C)C)C",
  "CN1C2=C(N=C1)N(C(=O)N(C2=O)C)C",
  "N1=CN(C)C2=C1N(C(=O)N(C)C2=O)C",
  "CN1C(=O)N(C(=O)C2=C1N=CN2C)C",
  "C1(=O)N(C)C2=C(N(C=N2)C)C(=O)N1C",
  "N1(C)C(=O)C2=C(N=CN2C)N(C)C1=O",
  "CN1C2=C(N=C1)N(C)C(=O)N(C2=O)C",
  "N1(C2=C(N(C)C(=O)N(C)C2=O)N=C1)C",
  "C1(=O)N(C(=O)C2=C(N=CN2C)N1C)C",
  "N1(C(=O)N(C)C2=C(N(C)C=N2)C1=O)C",
  "N1(C)C(=O)N(C)C(=O)C2=C1N=CN2C",
  "CN1C2=C(N(C)C=N2)C(=O)N(C1=O)C",
  "CN1C2=C(C(=O)N(C1=O)C)N(C)C=N2",
  "C12=C(N(C(=O)N(C)C1=O)C)N=CN2C",
  "C1(=O)N(C2=C(C(=O)N1C)N(C)C=N2)C",
  "C1(=O)C2=C(N=CN2C)N(C(=O)N1C)C",
  "O=C1N(C(=O)N(C2=C1N(C)C=N2)C)C",
  "N1(C)C2=C(N(C)C(=O)N(C2=O)C)N=C1",
  "N1=CN(C2=C1N(C)C(=O)N(C2=O)C)C",
  "CN1C(=O)N(C2=C(C1=O)N(C=N2)C)C",
  "C12=C(N=CN1C)N(C(=O)N(C)C2=O)C",
  "CN1C2=C(C(=O)N(C)C1=O)N(C=N2)C",
  "N1(C2=C(N(C=N2)C)C(=O)N(C)C1=O)C",
  "N1=CN(C)C2=C1N(C(=O)N(C2=O)C)C",
  "C1(=O)N(C(=O)N(C2=C1N(C)C=N2)C)C",
  "C1(=O)C2=C(N(C(=O)N1C)C)N=CN2C",
  "O=C1N(C(=O)N(C)C2=C1N(C)C=N2)C",
  "N1(C2=C(N(C)C(=O)N(C2=O)C)N=C1)C",
  "CN1C2=C(N(C(=O)N(C2=O)C)C)N=C1",
  "O=C1N(C)C2=C(C(=O)N1C)N(C)C=N2",
  "O=C1N(C)C(=O)C2=C(N=CN2C)N1C",
  "N1(C2=C(C(=O)N(C1=O)C)N(C)C=N2)C",
  "O=C1N(C)C(=O)N(C2=C1N(C)C=N2)C",
  "N1(C)C2=C(C(=O)N(C1=O)C)N(C)C=N2",
  "C1(=O)N(C)C2=C(N(C)C=N2)C(=O)N1C",
  "CN1C(=O)C2=C(N=CN2C)N(C)C1=O",
  "C1(=O)N(C)C(=O)N(C2=C1N(C=N2)C)C",
  "C1(=O)C2=C(N=CN2C)N(C)C(=O)N1C",
  "C1(=O)N(C2=C(N(C)C=N2)C(=O)N1C)C",
  "N1(C(=O)N(C(=O)C2=C1N=CN2C)C)C",
  "CN1C2=C(N(C)C(=O)N(C)C2=O)N=C1",
  "N1=CN(C2=C1N(C)C(=O)N(C)C2=O)C",
  "N1(C(=O)N(C)C2=C(N(C=N2)C)C1=O)C",
  "N1(C)C2=C(C(=O)N(C)C1=O)N(C=N2)C",
  "CN1C(=O)C2=C(N(C)C1=O)N=CN2C",
  "N1(C)C2=C(N=C1)N(C(=O)N(C)C2=O)C",
  "O=C1N(C)C(=O)N(C)C2=C1N(C=N2)C",
  "C1(=O)N(C(=O)N(C)C2=C1N(C=N2)C)C",
  "N1(C)C2=C(N=C1)N(C(=O)N(C2=O)C)C",
  "O=C1N(C2=C(N(C)C=N2)C(=O)N1C)C",
  "N1(C2=C(N=C1)N(C(=O)N(C)C2=O)C)C",
  "C12=C(N(C)C(=O)N(C1=O)C)N=CN2C",
  "N1(C)C(=O)C2=C(N(C1=O)C)N=CN2C",
  "N1(C2=C(N=C1)N(C(=O)N(C2=O)C)C)C",
  "O=C1N(C2=C(N(C=N2)C)C(=O)N1C)C",
  "CN1C(=O)N(C2=C(C1=O)N(C)C=N2)C",
  "CN1C(=O)N(C)C2=C(C1=O)N(C=N2)C",
  "N1(C2=C(N(C(=O)N(C2=O)C)C)N=C1)C",
  "O=C1N(C2=C(C(=O)N1C)N(C=N2)C)C",
  "N1=CN(C2=C1N(C(=O)N(C)C2=O)C)C",
  "C1(=O)N(C)C2=C(C(=O)N1C)N(C)C=N2",
  "CN1C2=C(N=C1)N(C(=O)N(C)C2=O)C",
  "N1(C2=C(C(=O)N(C)C1=O)N(C=N2)C)C",
  "N1(C)C(=O)N(C2=C(N(C)C=N2)C1=O)C",
  "CN1C2=C(N(C)C(=O)N(C2=O)C)N=C1",
  "N1(C)C(=O)N(C2=C(C1=O)N(C)C=N2)C",
  "N1(C(=O)N(C2=C(N(C)C=N2)C1=O)C)C",
  "N1(C(=O)N(C2=C(C1=O)N(C)C=N2)C)C",
  "O=C1N(C)C2=C(N(C)C=N2)C(=O)N1C",
  "O=C1N(C(=O)N(C2=C1N(C=N2)C)C)C",
  "N1(C(=O)N(C)C2=C(C1=O)N(C=N2)C)C",
  "N1(C(=O)N(C2=C(N(C=N2)C)C1=O)C)C",
  "N1(C)C2=C(N=C1)N(C)C(=O)N(C)C2=O",
  "N1(C)C(=O)N(C)C2=C(C1=O)N(C)C=N2",
  "N1(C)C2=C(N=C1)N(C)C(=O)N(C2=O)C",
  "CN1C(=O)N(C2=C(N(C=N2)C)C1=O)C",
  "CN1C(=O)N(C)C2=C(N(C)C=N2)C1=O",
  "N1(C)C2=C(N(C)C(=O)N(C)C2=O)N=C1",
  "N1(C)C2=C(N(C(=O)N(C2=O)C)C)N=C1",
  "N1(C)C2=C(N(C(=O)N(C)C2=O)C)N=C1"
};

bool
SameAtomTypes(const std::unordered_map<atom_type_t, int>& h1,
              const std::unordered_map<atom_type_t, int>& h2) {
  if (h1.size() != h2.size()) {
    return false;
  }

  for (const auto& [b1, c1] : h1) {
    const auto iter = h2.find(b1);
    if (iter == h2.end()) {
      std::cerr << "bit " << b1 << " not in h2\n";
      return false;
    }

    if (c1 != iter->second) {
      std::cerr << "Count mismatch on bit " << b1 << " counts " << c1 << " and " << iter->second << '\n';
      return false;
    }
  }

  return true;
}

bool
SameBitsAndCounts(const Sparse_Fingerprint_Creator& fp1,
                  const Sparse_Fingerprint_Creator& fp2) {
  const auto& bits1 = fp1.bits_found();
  const auto& bits2 = fp2.bits_found();

  if (bits1.size() != bits2.size()) {
    std::cerr << "Bit count mismatch " << bits1.size() << " and " << bits2.size() << '\n';
    return false;
  }

  for (const auto& [b, c] : bits1) {
    const auto iter = bits2.find(b);
    if (iter == bits2.end()) {
      cerr << "Bit " << b << " not found in b2\n";
      return false;
    }

    if (c != iter->second) {
      std::cerr << "Count mismatch on bit " << b << " counts " << c << " and " << iter->second << '\n';
      return false;
    }
  }

  return true;
}

TEST_F(TestECFingerPrint, TestRandomOrder) {
  const int nsmiles = aspirin.size();

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(aspirin[0]));

  Atom_Typing_Specification atom_typing;
  ASSERT_TRUE(atom_typing.build("UST:ACY"));

  const int matoms = m.natoms();

  std::unique_ptr<atom_type_t[]> atype = std::make_unique<atom_type_t[]>(matoms);
  ASSERT_TRUE(atom_typing.assign_atom_types(m, atype.get()));

  std::unordered_map<atom_type_t, int> atype_hash;
  for (int i = 0; i < matoms; ++i) {
    ++atype_hash[atype[i]];
  }

  ECFingerprint fp_creator;

  ProduceFingerprint bits;
  fp_creator.Fingerprint(m, nullptr, atype.get(), bits);

  for (int i = 1; i < nsmiles; ++i) {
    Molecule variant;
    EXPECT_TRUE(variant.build_from_smiles(aspirin[i])) << aspirin[i];
    ASSERT_EQ(matoms, variant.natoms());
    ASSERT_EQ(m.unique_smiles(), variant.unique_smiles());

    std::unique_ptr<atom_type_t[]> atype = std::make_unique<atom_type_t[]>(variant.natoms());
    ASSERT_TRUE(atom_typing.assign_atom_types(variant, atype.get()));

    std::unordered_map<atom_type_t, int> variant_atype_hash;
    for (int j = 0; j < matoms; ++j) {
      ++variant_atype_hash[atype[j]];
    }
    EXPECT_TRUE(SameAtomTypes(atype_hash, variant_atype_hash));

    ProduceFingerprint variant_bits;
    fp_creator.Fingerprint(variant, nullptr, atype.get(), variant_bits);

    EXPECT_TRUE(SameBitsAndCounts(bits.sfc(), variant_bits.sfc()));
  }
}

}  // namespace
}  // namespace ec_fingerprint
