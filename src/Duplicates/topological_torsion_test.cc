
// Tests for topological_torsion

#include <iostream>

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

#include "Molecule_Lib/atom_typing.h"
#include "topological_torsion.h"

namespace {
using topological_torsion::TorsionOptions;

using std::cerr;
using std::endl;

class TestTopologicalTorsion : public testing::Test {
  protected:
    Molecule _m;

    Atom_Typing_Specification _atom_typing;

    std::unique_ptr<atom_type_t[]> _atom_type;

    std::unique_ptr<int[]> _include_atom;

    int AssignAtomTypes(const char* ust);

    TorsionOptions _options;
};

int
TestTopologicalTorsion::AssignAtomTypes(const char * ust) {
  if (! _atom_typing.build(ust)) {
    cerr << "TestTopologicalTorsion:Cannot assign type '" << ust << "'\n";
    return 0;
  }

  const int natoms = _m.natoms();

  if (natoms > 0) {
    _atom_type.reset(new atom_type_t[natoms]);
    _include_atom.reset(new int[natoms]);
  }

  if (! _atom_typing.assign_atom_types(_m, _atom_type.get())) {
    cerr << "AssignAtomTypes:cannot assign '" << ust << "'\n";
    return 0;
  }

  return 1;
}

TEST_F(TestTopologicalTorsion, EmptyMolecule) {
  Sparse_Fingerprint_Creator sfc = TopologicalTorsion(_m, nullptr, nullptr, _options);
  EXPECT_EQ(sfc.nbits(), 0);
}

TEST_F(TestTopologicalTorsion, OneAtom) {
  ASSERT_TRUE(_m.build_from_smiles("C"));
  ASSERT_TRUE(AssignAtomTypes("UST:Y"));

  Sparse_Fingerprint_Creator sfc = TopologicalTorsion(_m, _atom_type.get(), nullptr, _options);
  EXPECT_EQ(sfc.nbits(), 0);
}

TEST_F(TestTopologicalTorsion, TwoAtom) {
  ASSERT_TRUE(_m.build_from_smiles("CC"));
  ASSERT_TRUE(AssignAtomTypes("UST:Y"));

  Sparse_Fingerprint_Creator sfc = TopologicalTorsion(_m, _atom_type.get(), nullptr, _options);
  EXPECT_EQ(sfc.nbits(), 0);
}

TEST_F(TestTopologicalTorsion, ThreeAtom) {
  ASSERT_TRUE(_m.build_from_smiles("CCC"));
  ASSERT_TRUE(AssignAtomTypes("UST:Y"));

  Sparse_Fingerprint_Creator sfc = TopologicalTorsion(_m, _atom_type.get(), nullptr, _options);
  EXPECT_EQ(sfc.nbits(), 0);
}

TEST_F(TestTopologicalTorsion, FourAtom) {
  ASSERT_TRUE(_m.build_from_smiles("CCCC"));
  ASSERT_TRUE(AssignAtomTypes("UST:Y"));

  Sparse_Fingerprint_Creator sfc = TopologicalTorsion(_m, _atom_type.get(), nullptr, _options);
  EXPECT_EQ(sfc.nbits(), 1);
}

TEST_F(TestTopologicalTorsion, OnePath) {
  ASSERT_TRUE(_m.build_from_smiles("ONCF"));
  ASSERT_TRUE(AssignAtomTypes("UST:Y"));

  Sparse_Fingerprint_Creator sfc = TopologicalTorsion(_m, _atom_type.get(), nullptr, _options);
  EXPECT_EQ(sfc.nbits(), 1);

  Molecule m2;
  ASSERT_TRUE(m2.build_from_smiles("FCNO"));
  ASSERT_TRUE(_atom_typing.assign_atom_types(m2, _atom_type.get()));
  Sparse_Fingerprint_Creator sfc2 = TopologicalTorsion(m2, _atom_type.get(), nullptr, _options);
  EXPECT_EQ(sfc2.nbits(), 1);

  const auto& bits1 = sfc.bits_found();
  const auto& bits2 = sfc2.bits_found();
  const auto f1 = bits1.cbegin();
  const auto f2 = bits2.cbegin();
  EXPECT_EQ(f1->first, f2->first);
}

TEST_F(TestTopologicalTorsion, DifferentBonds1) {
  ASSERT_TRUE(_m.build_from_smiles("C=CCC"));
  ASSERT_TRUE(AssignAtomTypes("UST:Y"));

  Sparse_Fingerprint_Creator sfc = TopologicalTorsion(_m, _atom_type.get(), nullptr, _options);
  EXPECT_EQ(sfc.nbits(), 1);

  Molecule m2;
  ASSERT_TRUE(m2.build_from_smiles("CCC=C"));
  ASSERT_TRUE(_atom_typing.assign_atom_types(m2, _atom_type.get()));
  Sparse_Fingerprint_Creator sfc2 = TopologicalTorsion(m2, _atom_type.get(), nullptr, _options);
  EXPECT_EQ(sfc2.nbits(), 1);

  const auto& bits1 = sfc.bits_found();
  const auto& bits2 = sfc2.bits_found();
  const auto f1 = bits1.cbegin();
  const auto f2 = bits2.cbegin();
  EXPECT_EQ(f1->first, f2->first);
}

TEST_F(TestTopologicalTorsion, DifferentBonds2) {
  ASSERT_TRUE(_m.build_from_smiles("C#CCC"));
  ASSERT_TRUE(AssignAtomTypes("UST:Y"));

  Sparse_Fingerprint_Creator sfc = TopologicalTorsion(_m, _atom_type.get(), nullptr, _options);
  EXPECT_EQ(sfc.nbits(), 1);

  Molecule m2;
  ASSERT_TRUE(m2.build_from_smiles("CCC#C"));
  ASSERT_TRUE(_atom_typing.assign_atom_types(m2, _atom_type.get()));
  Sparse_Fingerprint_Creator sfc2 = TopologicalTorsion(m2, _atom_type.get(), nullptr, _options);
  EXPECT_EQ(sfc2.nbits(), 1);

  const auto& bits1 = sfc.bits_found();
  const auto& bits2 = sfc2.bits_found();
  const auto f1 = bits1.cbegin();
  const auto f2 = bits2.cbegin();
  EXPECT_EQ(f1->first, f2->first);
}

TEST_F(TestTopologicalTorsion, Aromatic) {
  ASSERT_TRUE(_m.build_from_smiles("C1=CC=CC=C1"));
  ASSERT_TRUE(AssignAtomTypes("UST:Y"));

  Sparse_Fingerprint_Creator sfc = TopologicalTorsion(_m, _atom_type.get(), nullptr, _options);
  EXPECT_EQ(sfc.nbits(), 1);
}

TEST_F(TestTopologicalTorsion, Ring) {
  ASSERT_TRUE(_m.build_from_smiles("C1CC1"));
  ASSERT_TRUE(AssignAtomTypes("UST:Y"));

  _options.fingerprint_3_membered_rings = true;
  Sparse_Fingerprint_Creator sfc = TopologicalTorsion(_m, _atom_type.get(), nullptr, _options);
  EXPECT_EQ(sfc.nbits(), 1);
}

class TestMultipleSmiles : public testing::Test {
  protected:
    std::vector<Molecule> _mols;

    Atom_Typing_Specification _atom_typing;

    std::unique_ptr<atom_type_t[]> _atom_type;

    std::unique_ptr<int[]> _include_atom;

    int InitialiseAtomTyping(const char* ust, int natoms);

    TorsionOptions _options;
};

int
TestMultipleSmiles::InitialiseAtomTyping(const char * ust, int natoms) {
  if (! _atom_typing.build(ust)) {
    cerr << "TestTopologicalTorsion:Cannot assign type '" << ust << "'\n";
    return 0;
  }

  if (natoms > 0) {
    _atom_type.reset(new atom_type_t[natoms]);
    _include_atom.reset(new int[natoms]);
  }

  return 1;
}

TEST_F(TestMultipleSmiles, TestRingAllSingleBonds) {
  _mols.reserve(6);
  std::vector<IWString> smiles {"C1ON1", "C1NO1", "O1NC1", "O1CN1", "N1OC1", "N1CO1"};
  for (const auto& smi : smiles) {
    Molecule m;
    m.build_from_smiles(smi);
    _mols.emplace_back(std::move(m));
  }

  InitialiseAtomTyping("TT", _mols[0].natoms());

  _options.fingerprint_3_membered_rings = true;

  std::vector<Sparse_Fingerprint_Creator> sfcs;
  sfcs.reserve(_mols.size());
  for (Molecule& m : _mols) {
    ASSERT_TRUE(_atom_typing.assign_atom_types(m, _atom_type.get()));
    sfcs.emplace_back(std::move(TopologicalTorsion(m, _atom_type.get(), nullptr, _options)));
  }

  uint32_t bit = 0;
  for (const Sparse_Fingerprint_Creator& sfc : sfcs) {
    EXPECT_EQ(sfc.nbits(), 1);
    const auto bits = sfc.bits_found();
    const auto f = bits.cbegin();
    if (bit == 0) {
      bit = f->first;
    } else {
      EXPECT_EQ(f->first, bit);
    }
  }
}

TEST_F(TestMultipleSmiles, TestRingWithDoubleBond) {
  _mols.reserve(3);
  std::vector<IWString> smiles {"C1=CN1", "C=1NC=1", "N1C=C1"};
  for (const auto& smi : smiles) {
    Molecule m;
    m.build_from_smiles(smi);
    _mols.emplace_back(std::move(m));
  }
  InitialiseAtomTyping("TT", _mols[0].natoms());

  _options.fingerprint_3_membered_rings = true;

  std::vector<Sparse_Fingerprint_Creator> sfcs;
  sfcs.reserve(_mols.size());
  for (Molecule& m : _mols) {
    ASSERT_TRUE(_atom_typing.assign_atom_types(m, _atom_type.get()));
    sfcs.emplace_back(std::move(TopologicalTorsion(m, _atom_type.get(), nullptr, _options)));
  }

  uint32_t bit = 0;
  for (const Sparse_Fingerprint_Creator& sfc : sfcs) {
    EXPECT_EQ(sfc.nbits(), 1);
    const auto bits = sfc.bits_found();
    const auto f = bits.cbegin();
    if (bit == 0) {
      bit = f->first;
    } else {
      EXPECT_EQ(f->first, bit);
    }
  }
}

// Return the keys from a map.
std::vector<uint32_t>
GetKeys(const IW_Hash_Map<unsigned int, int>& hash) {
  std::vector<uint32_t> result;
  result.reserve(hash.size());
  for (auto [key, value] : hash) {
    result.push_back(key);
  }
  return result;
}

TEST_F(TestMultipleSmiles, TestFourMemberedRingSingle) {
  std::vector<IWString> smiles {"C1ONS1", "C1SNO1", "O1NSC1", "O1CSN1", "N1OCS1", "N1SCO1", "S1CON1", "S1NOC1"};
  _mols.reserve(smiles.size());
  for (const auto& smi : smiles) {
    Molecule m;
    m.build_from_smiles(smi);
    _mols.emplace_back(std::move(m));
  }

  for (size_t i = 1; i < _mols.size(); ++i) {
    EXPECT_EQ(_mols[i-1].unique_smiles(), _mols[i].unique_smiles());
  }
  InitialiseAtomTyping("TT", _mols[0].natoms());

  std::vector<Sparse_Fingerprint_Creator> sfcs;
  sfcs.reserve(_mols.size());
  for (Molecule& m : _mols) {
    ASSERT_TRUE(_atom_typing.assign_atom_types(m, _atom_type.get()));
    sfcs.emplace_back(std::move(TopologicalTorsion(m, _atom_type.get(), nullptr, _options)));
  }

  std::vector<uint32_t> expected;
  for (const Sparse_Fingerprint_Creator& sfc : sfcs) {
    EXPECT_EQ(sfc.nbits(), 4);
    const auto bits = sfc.bits_found();
    if (expected.empty()) {
      expected = GetKeys(bits);
    } else {
      const std::vector<uint32_t> got = GetKeys(bits);
      EXPECT_THAT(expected, testing::UnorderedElementsAreArray(got));
    }
  }
}

}  // namespace
