#include <algorithm>

#include "gtest/gtest.h"

#include "linear_fingerprint.h"

#include "aromatic.h"
#include "atom_typing.h"

#include "Foundational/iwmisc/sparse_fp_creator.h"

namespace linear_fingerprint {
namespace {


class TestLinearFingerprint : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    IWString _smiles;
    Molecule _m;
    int _matoms;
    Atom_Typing_Specification _atom_typing;
    Sparse_Fingerprint_Creator _sfc;
    LinearFingerprintGenerator _lfp;
    IWString _fingerprint;
};

void
TestLinearFingerprint::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

TEST_F(TestLinearFingerprint, EmptyMolecule)
{
  atom_type_t * atype = new atom_type_t[1]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  EXPECT_EQ(_lfp.Fingerprint(_m, nullptr, atype, _sfc), 0);
  EXPECT_EQ(_sfc.nbits(), 0l);
}

TEST_F(TestLinearFingerprint, TestSingleAtom)
{
  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  ASSERT_EQ(_matoms, 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);

  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_lfp.Fingerprint(_m, nullptr, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 1l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "TCF1F.2.....2");
}

TEST_F(TestLinearFingerprint, TestOneAtomExcluded)
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

  ASSERT_TRUE(_lfp.Fingerprint(_m, include_atom, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 1l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "TCF1F.2.....2");
}

TEST_F(TestLinearFingerprint, TestSingleAtomCType) 
{
  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 1);

  ASSERT_TRUE(_atom_typing.build("C"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);

  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_lfp.Fingerprint(_m, nullptr, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 1l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "28KAE.2.....2");
}

TEST_F(TestLinearFingerprint, TestFragments) 
{
  _smiles = "C.C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 2);
  EXPECT_EQ(_m.number_fragments(), 2);

  ASSERT_TRUE(_atom_typing.build("C"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_lfp.Fingerprint(_m, nullptr, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 1l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "28KAE.6.....2");
}

TEST_F(TestLinearFingerprint, TestBenzene) 
{
  _smiles = "c1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 6);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("C"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_lfp.Fingerprint(_m, nullptr, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 6l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "+TscE.saIWEKHbqE5bOeV.M4+UMabhY.9gQ7+.M4....2");
}

TEST_F(TestLinearFingerprint, TestButane) 
{
  _smiles = "CCCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 4);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_lfp.Fingerprint(_m, nullptr, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 4l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, ".2qmqbnYEoFz8bJezUREC.22.UA.2");
}

TEST_F(TestLinearFingerprint, TestMinPathTooLong)
{
  _smiles = "CCCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 4);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _lfp.set_min_length(4);
  ASSERT_TRUE(!_lfp.Fingerprint(_m, nullptr, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 0l);
}

TEST_F(TestLinearFingerprint, TestMinMaxSame)
{
  _smiles = "CCCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 4);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _lfp.set_min_length(1);
  _lfp.set_max_length(1);
  ASSERT_TRUE(_lfp.Fingerprint(_m, nullptr, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 1l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "zUREC.A.....2");
}

TEST_F(TestLinearFingerprint, TestPathSameAsMol)
{
  _smiles = "CCC(CC)CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 7);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _lfp.set_min_length(4);
  _lfp.set_max_length(4);
  ASSERT_TRUE(_lfp.Fingerprint(_m, nullptr, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 1l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "UL26W.A.....2");
}

TEST_F(TestLinearFingerprint, TestMissingCentreAtom)
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

  ASSERT_TRUE(_lfp.Fingerprint(_m, include_atom, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 2l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, "TCF1FDs5I1U4.k..3");
}

TEST_F(TestLinearFingerprint, TestRingBitsWithNoRings)
{
  _smiles = "CCC(CC)(CC)CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 9);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _lfp.set_fingerprint_ring_presence(true);
  ASSERT_TRUE(_lfp.Fingerprint(_m, nullptr, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 5l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, ".2qmqbnYEoFz8bJeUL26W.k70UPy+p.s0.......1");
}

TEST_F(TestLinearFingerprint, TestRingCaseWithoutRingsRequested)
{
  _smiles = "CCC1(CC)CC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 7);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  ASSERT_TRUE(_lfp.Fingerprint(_m, nullptr, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 5l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, ".2qmqbnYEoFz8bJeUL26W.c50ULy+p.s+k......1");
}

TEST_F(TestLinearFingerprint, TestRingCaseWithRingsRequested)
{
  _smiles = "CCC1(CC)CC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 7);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _lfp.set_fingerprint_ring_presence(true);
  ASSERT_TRUE(_lfp.Fingerprint(_m, nullptr, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 6l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, ".2qmqYdIndJ8XwIpTCF1F.20.URz8bJezUREC.20....2");
}

TEST_F(TestLinearFingerprint, TestPathsCanCross)
{
  _smiles = "CCC1(CC)CC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  _matoms = _m.natoms();

  EXPECT_EQ(_matoms, 7);
  EXPECT_EQ(_m.number_fragments(), 1);

  ASSERT_TRUE(_atom_typing.build("UST:Y"));

  atom_type_t * atype = new atom_type_t[_matoms]; std::unique_ptr<atom_type_t[]> free_atype(atype);
  ASSERT_TRUE(_atom_typing.assign_atom_types(_m, atype));

  _lfp.set_paths_can_cross(true);
  ASSERT_TRUE(_lfp.Fingerprint(_m, nullptr, atype, _sfc));

  EXPECT_EQ(_sfc.nbits(), 8l);

  _sfc.daylight_ascii_form_with_counts_encoded(_fingerprint);
  EXPECT_EQ(_fingerprint, ".2qmqU8IRbE2qtg4TCF1F+.4.URz8bJeUL26W6Crz7vy+p.s0UY2+k..1");
}

}  // namespace
}  // namespace linear_fingerprint
