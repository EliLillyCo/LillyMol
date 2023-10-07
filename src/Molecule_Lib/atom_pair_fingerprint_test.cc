#include <algorithm>

#include "googletest/include/gtest/gtest.h"

#include "atom_pair_fingerprint.h"

#include "aromatic.h"
#include "atom_typing.h"

#include "Foundational/iwmisc/sparse_fp_creator.h"

namespace atom_pair_fingerprint {
namespace {


class TestAPFingerprint : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    IWString _smiles;
    Molecule _m;
    int _matoms;
    Atom_Typing_Specification _atom_typing;
    AtomPairFingerprint _apfp;
    Sparse_Fingerprint_Creator _sfp;
    IWString _fingerprint;
};

void
TestAPFingerprint::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

TEST_F(TestAPFingerprint, EmptyMolecule)
{
  atom_type_t * atype = new atom_type_t[1]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  EXPECT_EQ(_apfp.Fingerprint(_m, nullptr, atype, _sfp), 0);
  EXPECT_EQ(_sfp.nbits(), 0l);
}

TEST_F(TestAPFingerprint, TestSingleAtom)
{
  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  ASSERT_EQ(_matoms, 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);

  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_EQ(_apfp.Fingerprint(_m, nullptr, atype, _sfp), 0);

  EXPECT_EQ(_sfp.nbits(), 0l);
}

TEST_F(TestAPFingerprint, TestOneAtomExcluded)
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

  ASSERT_EQ(_apfp.Fingerprint(_m, include_atom, atype, _sfp), 0);

  EXPECT_EQ(_sfp.nbits(), 0l);
}

TEST_F(TestAPFingerprint, TestSingleAtomCType) 
{
  _smiles = "CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 2);

  ASSERT_TRUE(_atom_typing.build("C"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);

  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_apfp.Fingerprint(_m, nullptr, atype, _sfp));

  EXPECT_EQ(_sfp.nbits(), 1l);

  _sfp.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "1NS7n.2.....2");
}

TEST_F(TestAPFingerprint, TestFragments) 
{
  _smiles = "C.C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 2);
  EXPECT_EQ(_m.number_fragments(), 2);

  ASSERT_TRUE(_atom_typing.build("C"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _apfp.set_min_separation(1);
  ASSERT_EQ(_apfp.Fingerprint(_m, nullptr, atype, _sfp), 0);

  EXPECT_EQ(_sfp.nbits(), 0l);
}

TEST_F(TestAPFingerprint, TestBenzene) 
{
  _smiles = "c1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 6);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("C"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_apfp.Fingerprint(_m, nullptr, atype, _sfp));

  EXPECT_EQ(_sfp.nbits(), 3l);

  _sfp.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "+cK1v.OCsaU4a21Y+UM1....1");
}

TEST_F(TestAPFingerprint, TestButane) 
{
  _smiles = "CCCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 4);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_apfp.Fingerprint(_m, nullptr, atype, _sfp));

  EXPECT_EQ(_sfp.nbits(), 3l);

  _sfp.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "VmxGH6Qsiz45EWKK.k6+....1");
}

TEST_F(TestAPFingerprint, TestMinSeparationTooLong)
{
  _smiles = "CCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 3);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _apfp.set_min_separation(4);
  ASSERT_EQ(_apfp.Fingerprint(_m, nullptr, atype, _sfp), 0);

  EXPECT_EQ(_sfp.nbits(), 0l);
}

TEST_F(TestAPFingerprint, TestMinMaxSame)
{
  _smiles = "CCCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 4);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _apfp.set_min_separation(3);
  _apfp.set_max_separation(3);
  ASSERT_TRUE(_apfp.Fingerprint(_m, nullptr, atype, _sfp));

  EXPECT_EQ(_sfp.nbits(), 1l);

  _sfp.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "Vo6ZZU2.....2");
}

TEST_F(TestAPFingerprint, TestMissingCentreAtom)
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

  ASSERT_TRUE(_apfp.Fingerprint(_m, include_atom, atype, _sfp));

  EXPECT_EQ(_sfp.nbits(), 4l);

  _sfp.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "VmxGH6Qsiz45EWKKVoiDCkA1+UA.2");
}

TEST_F(TestAPFingerprint, TestOutOfRangeOff)
{
  _smiles = "CCCCCCCCCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 10);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:ACHY"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _apfp.set_max_separation(5);
  ASSERT_TRUE(_apfp.Fingerprint(_m, nullptr, atype, _sfp));

  EXPECT_EQ(_sfp.nbits(), 10l);

  _sfp.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "mx2tnAjFCib9qeZZmxefbkQ0+U99t+XymyEQJQjhW7T9vMo9+E62.gjqy119xjr+.k6.....1");
}

TEST_F(TestAPFingerprint, TestOutOfRangeOn)
{
  _smiles = "CCCCCCCCCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 10);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:ACHY"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _apfp.set_max_separation(5);
  _apfp.set_include_out_of_range_separations(true);
  ASSERT_TRUE(_apfp.Fingerprint(_m, nullptr, atype, _sfp));

  EXPECT_EQ(_sfp.nbits(), 13l);

  _sfp.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "mx2tnAjFCib9qeZZmxefbkQ0+U99t+XymyEQJQjhW7T9vMo9+E62.gjqy119xjr+n.+bmQk.PbQ1.UA4qBMXCU2.....2");
}

TEST_F(TestAPFingerprint, TestCoverage)
{
  _smiles = "C1(=CN(C)C(=C1)C)C1=CSC(=N1)NC(N)=N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 16);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:ACHY"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _apfp.set_examine_coverage(true);
  _apfp.set_max_separation(8);
  ASSERT_TRUE(_apfp.Fingerprint(_m, nullptr, atype, _sfp));

  EXPECT_EQ(_sfp.nbits(), 87l);

  _sfp.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "2BXLNF1NFkIEuvoL2Cj5zk2+.E2EuxDz2FXhZl2VHhsFB1xY.E2+.F2oJqEkOcO6A4eATX+yJbM+.E2+A8ElKH0YakBH35BqVyTKuE2+.E45txPiVyTNuMTlAYi5yfxK.E2+.MTvPkm6+.LzW.EFXsUBPxY+.E2+W+Pf0sUUEsq667KzW0bYFk6+.E87YZYgWN7Q26aGM.W7YaHb.EE+.caGNCe7YaSdWN8TSMaPkkM+.E20WNj6ncaPo9u7axdwWNjOUUI+.E47b2yUWOIgs6aZELG7dJ.O.E60.MaZK2W7dJVLWOJVG6aZzwQ+.E2+WOuWGcaigWe7fgKaWOvEfUA0.E47i.0IWPUvGcasGCm7i3Tg.E2+.MatM+K7kNCKWQ5+Dcb+kJk+.U2+WQ5HDcbIgS87pAbWeiKnVU2+.E4etSsJeiwLOefj5peevmYM.E20.Ofj8Fuey6EcejWE28fsbeo+.U60ejWawug032uf.ju8ekhlT.2+.U8f0tS6ekiLcOg9dcWf5cUg.E6+.egSbGnMdjGCq9+e7hWtoFU0.E2+qAn8rxXLby1MrwZL.U60....1");

  EXPECT_EQ(_m.smiles(), "[C:15]1(=[CH:15][N:15]([CH3:13])[C:15](=[CH:15]1)[CH3:13])[C:15]1=[CH:15][S:15][C:15](=[N:15]1)[NH:15][C:15]([NH2:13])=[NH:13]");
}

}  // namespace
}  // namespace atom_pair_fingerprint
