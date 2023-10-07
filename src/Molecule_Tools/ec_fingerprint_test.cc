#include <algorithm>
#include <iostream>

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

}  // namespace
}  // namespace ec_fingerprint
