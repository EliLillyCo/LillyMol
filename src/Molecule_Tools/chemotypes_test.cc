#include "Molecule_Tools/chemotypes.h"

#include <algorithm>
#include <memory>
#include <vector>

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/substructure.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

namespace {

std::unique_ptr<int[]>
RingSystemMembership(Molecule& m, int& number_ring_systems) {
  std::unique_ptr<int[]> ring_system = std::make_unique<int[]>(m.natoms());
  number_ring_systems =
      m.label_atoms_by_ring_system_including_spiro_fused(ring_system.get());
  return ring_system;
}

Substructure_Query*
AddSmarts(resizable_array_p<Substructure_Query>& queries, const char* smarts) {
  Substructure_Query* query = new Substructure_Query;
  EXPECT_TRUE(query->create_from_smarts(smarts));
  queries.add(query);
  return query;
}

void
ExpectSameEmbedding(const Set_of_Atoms& lhs, const Set_of_Atoms& rhs) {
  EXPECT_THAT(lhs, testing::ElementsAreArray(rhs));
}

chemotypes::ChemotypeQueryMatch
MatchWithSeedAtom(Molecule& m, atom_number_t seed_atom) {
  chemotypes::ChemotypeQueryMatch result;
  result.query_index = 0;
  result.hits = 1;
  result.embedding = Set_of_Atoms({seed_atom});
  result.seed_atom = seed_atom;
  result.ring_system.resize(m.natoms(), 0);
  m.label_atoms_by_ring_system_including_spiro_fused(result.ring_system.data());
  result.seed_ring_system = result.ring_system[seed_atom];
  return result;
}

int
CountMaskedAtoms(const std::vector<int>& mask) {
  return std::count_if(mask.begin(), mask.end(), [](int value) { return value != 0; });
}



TEST(Chemotypes, FirstChemotypeQueryMatchReportsNoQueryMatch) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1"));

  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[N]");

  chemotypes::ChemotypeQueryMatch result;
  EXPECT_EQ(chemotypes::FirstChemotypeQueryMatch(m, queries, result),
            chemotypes::ChemotypeQueryMatchStatus::kNoQueryMatch);
  EXPECT_EQ(result.query_index, -1);
}

TEST(Chemotypes, FirstChemotypeQueryMatchRequiresRingAtomInEmbedding) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CC"));

  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[C]");

  chemotypes::ChemotypeQueryMatch result;
  EXPECT_EQ(chemotypes::FirstChemotypeQueryMatch(m, queries, result),
            chemotypes::ChemotypeQueryMatchStatus::kMatchedQueryNoRingAtom);
  EXPECT_EQ(result.query_index, 0);
  EXPECT_EQ(result.hits, 2);
  EXPECT_EQ(result.seed_atom, kInvalidAtomNumber);
  EXPECT_EQ(result.seed_ring_system, 0);
}

TEST(Chemotypes, FirstChemotypeQueryMatchUsesFirstMatchingQuery) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("Cc1ccccc1"));

  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[N]");
  AddSmarts(queries, "[c]");
  AddSmarts(queries, "[C]");

  chemotypes::ChemotypeQueryMatch result;
  ASSERT_EQ(chemotypes::FirstChemotypeQueryMatch(m, queries, result),
            chemotypes::ChemotypeQueryMatchStatus::kMatched);
  EXPECT_EQ(result.query_index, 1);
  EXPECT_GT(result.hits, 0);
  ASSERT_NE(result.seed_atom, kInvalidAtomNumber);
  EXPECT_EQ(result.seed_ring_system, result.ring_system[result.seed_atom]);
}

TEST(Chemotypes, FirstChemotypeQueryMatchUsesFirstEmbedding) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("Cc1ccccc1"));

  resizable_array_p<Substructure_Query> queries;
  Substructure_Query* query = AddSmarts(queries, "[c]");

  Substructure_Results expected;
  ASSERT_GT(query->substructure_search(m, expected), 0u);

  chemotypes::ChemotypeQueryMatch result;
  ASSERT_EQ(chemotypes::FirstChemotypeQueryMatch(m, queries, result),
            chemotypes::ChemotypeQueryMatchStatus::kMatched);
  ExpectSameEmbedding(result.embedding, *expected.embedding(0));
  EXPECT_EQ(result.seed_atom, expected.embedding(0)->item(0));
}

TEST(Chemotypes, FirstChemotypeQueryMatchAllowsMultipleEmbeddingsInOneRingSystem) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1"));

  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[c]");

  chemotypes::ChemotypeQueryMatch result;
  ASSERT_EQ(chemotypes::FirstChemotypeQueryMatch(m, queries, result),
            chemotypes::ChemotypeQueryMatchStatus::kMatched);
  EXPECT_EQ(result.seed_ring_system, 1);
}

TEST(Chemotypes, FirstChemotypeQueryMatchRejectsMultipleRingSystems) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1-c1ccccc1"));

  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[c]");

  chemotypes::ChemotypeQueryMatch result;
  EXPECT_EQ(chemotypes::FirstChemotypeQueryMatch(m, queries, result),
            chemotypes::ChemotypeQueryMatchStatus::kAmbiguousQueryMatches);
}

TEST(Chemotypes, FirstChemotypeQueryMatchCanUseFirstEmbeddingAcrossRingSystems) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1-c1ccccc1"));

  resizable_array_p<Substructure_Query> queries;
  Substructure_Query* query = AddSmarts(queries, "[c]");

  Substructure_Results expected;
  ASSERT_GT(query->substructure_search(m, expected), 0u);

  chemotypes::ChemotypeQueryMatch result;
  constexpr bool kChooseFirstEmbedding = true;
  ASSERT_EQ(chemotypes::FirstChemotypeQueryMatch(m, queries, result,
                                                kChooseFirstEmbedding),
            chemotypes::ChemotypeQueryMatchStatus::kMatched);
  ExpectSameEmbedding(result.embedding, *expected.embedding(0));
}


TEST(Chemotypes, ChemotypeAtomMaskSeedRingOnly) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1-c1ccccc1"));

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 0);
  const std::vector<int> mask = chemotypes::ChemotypeAtomMask(m, match, 0);

  EXPECT_EQ(CountMaskedAtoms(mask), 6);
  for (atom_number_t atom = 0; atom < m.natoms(); ++atom) {
    EXPECT_EQ(mask[atom], match.ring_system[atom] == match.seed_ring_system ?
                              chemotypes::kChemotypeCoreAtom :
                              chemotypes::kChemotypeNotKept);
  }
}

TEST(Chemotypes, ChemotypeAtomMaskIncludesRequestedAdjacentRingSystem) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1-c1ccccc1"));

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 0);
  const std::vector<int> mask = chemotypes::ChemotypeAtomMask(m, match, 1);

  EXPECT_EQ(CountMaskedAtoms(mask), 12);
  EXPECT_THAT(mask, testing::Each(chemotypes::kChemotypeCoreAtom));
}

TEST(Chemotypes, ChemotypeAtomMaskDoesNotIncludeRingBeyondInterveningRing) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1CCc1ccccc1CCc1ccccc1"));

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 0);
  const std::vector<int> mask = chemotypes::ChemotypeAtomMask(m, match, 2);

  EXPECT_EQ(CountMaskedAtoms(mask), 14);
  for (atom_number_t atom = 0; atom < m.natoms(); ++atom) {
    if (match.ring_system[atom] == match.ring_system[m.natoms() - 1]) {
      EXPECT_EQ(mask[atom], 0);
    }
  }
}

TEST(Chemotypes, ChemotypeAtomMaskHonoursMaxDistance) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1CCc1ccccc1"));

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 0);
  EXPECT_EQ(CountMaskedAtoms(chemotypes::ChemotypeAtomMask(m, match, 1, 2)), 6);
  EXPECT_EQ(CountMaskedAtoms(chemotypes::ChemotypeAtomMask(m, match, 1, 3)), 14);
}

TEST(Chemotypes, ChemotypeAtomMaskIncludesMultipleAdjacentRingSystems) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1cc(-c2ccccc2)ccc1-c1ccccc1"));

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 0);
  EXPECT_EQ(CountMaskedAtoms(chemotypes::ChemotypeAtomMask(m, match, 1)), 12);
  EXPECT_EQ(CountMaskedAtoms(chemotypes::ChemotypeAtomMask(m, match, 2)), 18);
}

TEST(Chemotypes, ChemotypeAtomMaskCanIncludeTiesAtAdjacentRingCutoff) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1cc(-c2ccccc2)ccc1-c1ccccc1"));

  chemotypes::ChemotypeOptions options;
  options.adjacent_ring_systems_to_include = 1;
  options.include_tied_adjacent_ring_systems = 1;

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 0);
  EXPECT_EQ(CountMaskedAtoms(chemotypes::ChemotypeAtomMask(m, match, options)), 18);
}

TEST(Chemotypes, ChemotypeAtomMaskInvalidMatchReturnsAllZero) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1"));

  chemotypes::ChemotypeQueryMatch match;
  match.seed_ring_system = 1;
  const std::vector<int> mask = chemotypes::ChemotypeAtomMask(m, match, 1);

  EXPECT_EQ(mask.size(), static_cast<size_t>(m.natoms()));
  EXPECT_THAT(mask, testing::Each(0));
}



TEST(Chemotypes, ChemotypeAtomMaskCanReuseScratchStorage) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1CCc1ccccc1"));

  chemotypes::ChemotypeOptions options;
  options.adjacent_ring_systems_to_include = 1;
  chemotypes::ChemotypeScratch scratch;

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 0);
  EXPECT_EQ(CountMaskedAtoms(chemotypes::ChemotypeAtomMask(m, match, options, scratch)), 14);
  EXPECT_EQ(CountMaskedAtoms(chemotypes::ChemotypeAtomMask(m, match, options, scratch)), 14);
}

TEST(Chemotypes, ChemotypeAtomMaskCanSuppressLinkerAtoms) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1CCc1ccccc1"));

  chemotypes::ChemotypeOptions options;
  options.adjacent_ring_systems_to_include = 1;
  options.include_linker_atoms = 0;

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 0);
  EXPECT_EQ(CountMaskedAtoms(chemotypes::ChemotypeAtomMask(m, match, options)), 12);
}

TEST(Chemotypes, ChemotypeAtomMaskIncludesExocyclicDoubleBondedAtoms) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("O=C1CCCCC1"));

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 1);
  const std::vector<int> mask = chemotypes::ChemotypeAtomMask(m, match, 0);

  EXPECT_EQ(CountMaskedAtoms(mask), 7);
  EXPECT_EQ(mask[0], chemotypes::kChemotypeCoreAtom);
}


TEST(Chemotypes, ChemotypeAtomMaskIncludesTerminalDoubleBondOnLinkerAtom) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1C(=O)c1ccccc1"));

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 0);
  const std::vector<int> mask = chemotypes::ChemotypeAtomMask(m, match, 1);

  EXPECT_EQ(CountMaskedAtoms(mask), 14);
  EXPECT_EQ(mask[7], chemotypes::kChemotypeCoreAtom);
}

TEST(Chemotypes, ChemotypeAtomMaskCanSuppressExocyclicDoubleBondedAtoms) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("O=C1CCCCC1"));

  chemotypes::ChemotypeOptions options;
  options.include_exocyclic_double_bonded_atoms = 0;

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 1);
  const std::vector<int> mask = chemotypes::ChemotypeAtomMask(m, match, options);

  EXPECT_EQ(CountMaskedAtoms(mask), 6);
  EXPECT_EQ(mask[0], 0);
}

TEST(Chemotypes, ChemotypeAtomMaskIncludesOptionalAttachedAtoms) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CN1CCCCC1"));

  chemotypes::ChemotypeOptions options;
  options.include_attached_atoms = 1;

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 2);
  const std::vector<int> mask = chemotypes::ChemotypeAtomMask(m, match, options);

  EXPECT_EQ(CountMaskedAtoms(mask), 7);
  EXPECT_EQ(mask[0], chemotypes::kChemotypeAttachedAtom);
}

TEST(Chemotypes, ChemotypeAtomMaskCanIgnoreSinglyConnectedAttachedAtoms) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CN1CCCCC1"));

  chemotypes::ChemotypeOptions options;
  options.include_attached_atoms = 1;
  options.ignore_singly_connected_attached_atoms = 1;

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 2);
  const std::vector<int> mask = chemotypes::ChemotypeAtomMask(m, match, options);

  EXPECT_EQ(CountMaskedAtoms(mask), 6);
  EXPECT_EQ(mask[0], 0);
}



TEST(Chemotypes, ChemotypeAtomMaskAttachedAtomsDoNotPullInUnselectedRing) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1-c1ccccc1"));

  chemotypes::ChemotypeOptions options;
  options.include_attached_atoms = 1;

  const chemotypes::ChemotypeQueryMatch match = MatchWithSeedAtom(m, 0);
  const std::vector<int> mask = chemotypes::ChemotypeAtomMask(m, match, options);

  EXPECT_EQ(CountMaskedAtoms(mask), 6);
  for (atom_number_t atom = 0; atom < m.natoms(); ++atom) {
    if (match.ring_system[atom] != match.seed_ring_system) {
      EXPECT_EQ(mask[atom], 0);
    }
  }
}

TEST(Chemotypes, ReduceToChemotypeRemovesAtomsNotInMask) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("COc1ccccc1CCc1ccccc1F"));

  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[c]");

  chemotypes::ChemotypeOptions options;
  options.adjacent_ring_systems_to_include = 1;
  options.choose_first_embedding = true;
  chemotypes::ChemotypeScratch scratch;
  chemotypes::ChemotypeQueryMatch match;

  ASSERT_EQ(chemotypes::ReduceToChemotype(m, queries, options, scratch, match),
            chemotypes::ChemotypeQueryMatchStatus::kMatched);
  EXPECT_EQ(m.natoms(), 14);
}

TEST(Chemotypes, ReduceToChemotypeReportsNoQueryMatch) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1"));

  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[N]");

  chemotypes::ChemotypeOptions options;
  chemotypes::ChemotypeScratch scratch;
  chemotypes::ChemotypeQueryMatch match;

  EXPECT_EQ(chemotypes::ReduceToChemotype(m, queries, options, scratch, match),
            chemotypes::ChemotypeQueryMatchStatus::kNoQueryMatch);
  EXPECT_EQ(m.natoms(), 6);
}


TEST(Chemotypes, ReduceToChemotypeLabelsRingExitPointWithFixedIsotope) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CN1CCCCC1"));

  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[N]");

  constexpr isotope_t kExitPointIsotope = 99;
  chemotypes::ChemotypeOptions options;
  options.isotope_for_exit_points = kExitPointIsotope;
  chemotypes::ChemotypeScratch scratch;
  chemotypes::ChemotypeQueryMatch match;

  ASSERT_EQ(chemotypes::ReduceToChemotype(m, queries, options, scratch, match),
            chemotypes::ChemotypeQueryMatchStatus::kMatched);

  ASSERT_EQ(m.natoms(), 6);
  int labelled_atoms = 0;
  for (atom_number_t atom = 0; atom < m.natoms(); ++atom) {
    if (m.isotope(atom) == kExitPointIsotope) {
      ++labelled_atoms;
      EXPECT_EQ(m.atomic_number(atom), static_cast<atomic_number_t>(7));
    } else {
      EXPECT_EQ(m.isotope(atom), static_cast<isotope_t>(0));
    }
  }
  EXPECT_EQ(labelled_atoms, 1);
}

TEST(Chemotypes, ReduceToChemotypeLabelsNonTerminalAttachmentAtomsByAtomType) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCN1CCCCC1"));

  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[N]");

  chemotypes::ChemotypeOptions options;
  chemotypes::ChemotypeScratch scratch;
  chemotypes::ChemotypeQueryMatch match;
  Atom_Typing_Specification atom_typing;
  atom_typing.set_atom_type(IWATTYPE_Z);

  ASSERT_EQ(chemotypes::ReduceToChemotype(m, queries, options, scratch, match,
                                          &atom_typing),
            chemotypes::ChemotypeQueryMatchStatus::kMatched);

  ASSERT_EQ(m.natoms(), 7);
  int labelled_atoms = 0;
  for (atom_number_t atom = 0; atom < m.natoms(); ++atom) {
    if (m.isotope(atom)) {
      ++labelled_atoms;
      // The labelled atom was non-terminal in the parent molecule, but becomes
      // terminal after the rest of the ethyl group is removed.
      EXPECT_EQ(m.ncon(atom), 1);
      EXPECT_EQ(m.atomic_number(atom), static_cast<atomic_number_t>(6));
      EXPECT_EQ(m.isotope(atom), static_cast<isotope_t>(6));
    }
  }
  EXPECT_EQ(labelled_atoms, 1);
}

TEST(Chemotypes, ReduceToChemotypeDoesNotLabelSinglyConnectedAttachmentAtoms) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CN1CCCCC1"));

  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[N]");

  chemotypes::ChemotypeOptions options;
  chemotypes::ChemotypeScratch scratch;
  chemotypes::ChemotypeQueryMatch match;
  Atom_Typing_Specification atom_typing;
  atom_typing.set_atom_type(IWATTYPE_Z);

  ASSERT_EQ(chemotypes::ReduceToChemotype(m, queries, options, scratch, match,
                                          &atom_typing),
            chemotypes::ChemotypeQueryMatchStatus::kMatched);

  ASSERT_EQ(m.natoms(), 7);
  for (atom_number_t atom = 0; atom < m.natoms(); ++atom) {
    EXPECT_EQ(m.isotope(atom), static_cast<isotope_t>(0)) << " atom " << atom;
  }
}

TEST(Chemotypes, ReduceToChemotypeDoesNotLabelTerminalDoubleBondedAtoms) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("O=C1CCCCC1"));

  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[R]");

  chemotypes::ChemotypeOptions options;
  chemotypes::ChemotypeScratch scratch;
  chemotypes::ChemotypeQueryMatch match;
  Atom_Typing_Specification atom_typing;
  atom_typing.set_atom_type(IWATTYPE_Z);

  ASSERT_EQ(chemotypes::ReduceToChemotype(m, queries, options, scratch, match,
                                          &atom_typing),
            chemotypes::ChemotypeQueryMatchStatus::kMatched);

  ASSERT_EQ(m.natoms(), 7);
  for (atom_number_t atom = 0; atom < m.natoms(); ++atom) {
    EXPECT_EQ(m.isotope(atom), static_cast<isotope_t>(0)) << " atom " << atom;
  }
}

TEST(Chemotypes, BiphenylRingsAreDirectlyAdjacent) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1-c1ccccc1"));

  int number_ring_systems;
  std::unique_ptr<int[]> ring_system = RingSystemMembership(m, number_ring_systems);
  ASSERT_EQ(number_ring_systems, 2);

  const int seed_ring_system = ring_system[0];
  Set_of_Atoms embedding({0});
  const std::vector<chemotypes::AdjacentRingSystem> adjacent =
      chemotypes::DirectlyAdjacentRingSystems(m, ring_system.get(), seed_ring_system,
                                              embedding);

  ASSERT_EQ(adjacent.size(), 1u);
  EXPECT_NE(adjacent[0].ring_system, seed_ring_system);
  EXPECT_EQ(adjacent[0].distance, 1);
}

TEST(Chemotypes, MaxDistanceLimitsAdjacentRingSystems) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1CCc1ccccc1"));

  int number_ring_systems;
  std::unique_ptr<int[]> ring_system = RingSystemMembership(m, number_ring_systems);
  ASSERT_EQ(number_ring_systems, 2);

  const int seed_ring_system = ring_system[0];
  Set_of_Atoms embedding({0});
  EXPECT_TRUE(chemotypes::DirectlyAdjacentRingSystems(m, ring_system.get(),
                                                      seed_ring_system,
                                                      embedding, 2).empty());

  const std::vector<chemotypes::AdjacentRingSystem> adjacent =
      chemotypes::DirectlyAdjacentRingSystems(m, ring_system.get(), seed_ring_system,
                                              embedding, 3);
  ASSERT_EQ(adjacent.size(), 1u);
  EXPECT_EQ(adjacent[0].distance, 3);
}

TEST(Chemotypes, RingBeyondInterveningRingIsNotDirectlyAdjacent) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1CCc1ccccc1CCc1ccccc1"));

  int number_ring_systems;
  std::unique_ptr<int[]> ring_system = RingSystemMembership(m, number_ring_systems);
  ASSERT_EQ(number_ring_systems, 3);

  const int seed_ring_system = ring_system[0];
  const int terminal_ring_system = ring_system[m.natoms() - 1];
  ASSERT_NE(terminal_ring_system, seed_ring_system);

  Set_of_Atoms embedding({0});
  const std::vector<chemotypes::AdjacentRingSystem> adjacent =
      chemotypes::DirectlyAdjacentRingSystems(m, ring_system.get(), seed_ring_system,
                                              embedding);

  ASSERT_EQ(adjacent.size(), 1u);
  EXPECT_NE(adjacent[0].ring_system, terminal_ring_system);
}


TEST(Chemotypes, LargerAdjacentRingPreferenceUsesAtomsDownTheCandidateBond) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1C(CC1CC1)CC1CC1C"));

  int number_ring_systems;
  std::unique_ptr<int[]> ring_system = RingSystemMembership(m, number_ring_systems);
  ASSERT_EQ(number_ring_systems, 3);

  const int seed_ring_system = ring_system[0];
  Set_of_Atoms embedding({0});
  const std::vector<chemotypes::AdjacentRingSystem> adjacent =
      chemotypes::DirectlyAdjacentRingSystems(m, ring_system.get(), seed_ring_system,
                                              embedding, 0,
                                              true);  // prefer larger adjacent ring system.

  ASSERT_EQ(adjacent.size(), 2u);
  EXPECT_EQ(adjacent[0].distance, adjacent[1].distance);
  EXPECT_GT(adjacent[0].atoms_down_the_bond, adjacent[1].atoms_down_the_bond);
}

TEST(Chemotypes, LargerAdjacentRingPreferenceIsIndependentOfInputOrder) {
  Molecule e3b;
  ASSERT_TRUE(e3b.build_from_smiles("c1ccccc1C(CC1CC1)CC1CC1C"));
  Molecule e3c;
  ASSERT_TRUE(e3c.build_from_smiles("c1ccccc1C(CC1CC1C)CC1CC1"));

  int number_ring_systems1;
  std::unique_ptr<int[]> ring_system1 = RingSystemMembership(e3b, number_ring_systems1);
  ASSERT_EQ(number_ring_systems1, 3);
  int number_ring_systems2;
  std::unique_ptr<int[]> ring_system2 = RingSystemMembership(e3c, number_ring_systems2);
  ASSERT_EQ(number_ring_systems2, 3);

  Set_of_Atoms embedding({0});
  const std::vector<chemotypes::AdjacentRingSystem> adjacent1 =
      chemotypes::DirectlyAdjacentRingSystems(e3b, ring_system1.get(), ring_system1[0],
                                              embedding, 0, true);
  const std::vector<chemotypes::AdjacentRingSystem> adjacent2 =
      chemotypes::DirectlyAdjacentRingSystems(e3c, ring_system2.get(), ring_system2[0],
                                              embedding, 0, true);

  ASSERT_EQ(adjacent1.size(), 2u);
  ASSERT_EQ(adjacent2.size(), 2u);
  EXPECT_EQ(adjacent1[0].atoms_down_the_bond, adjacent2[0].atoms_down_the_bond);
  EXPECT_GT(adjacent1[0].atoms_down_the_bond, adjacent1[1].atoms_down_the_bond);
  EXPECT_GT(adjacent2[0].atoms_down_the_bond, adjacent2[1].atoms_down_the_bond);
}

TEST(Chemotypes, QueryAtomOrderBreaksEquivalentAdjacentRingChoices) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1cc(-c2ccccc2)ccc1-c1ccccc1"));

  int number_ring_systems;
  std::unique_ptr<int[]> ring_system = RingSystemMembership(m, number_ring_systems);
  ASSERT_EQ(number_ring_systems, 3);

  const int seed_ring_system = ring_system[0];
  Set_of_Atoms no_query_atoms;
  const std::vector<chemotypes::AdjacentRingSystem> adjacent =
      chemotypes::DirectlyAdjacentRingSystems(m, ring_system.get(), seed_ring_system,
                                              no_query_atoms);
  ASSERT_EQ(adjacent.size(), 2u);
  ASSERT_NE(adjacent[0].seed_atom, adjacent[1].seed_atom);

  Set_of_Atoms embedding({adjacent[1].seed_atom, adjacent[0].seed_atom});
  const std::vector<chemotypes::AdjacentRingSystem> reordered =
      chemotypes::DirectlyAdjacentRingSystems(m, ring_system.get(), seed_ring_system,
                                              embedding);

  ASSERT_EQ(reordered.size(), 2u);
  EXPECT_EQ(reordered[0].seed_atom, adjacent[1].seed_atom);
  EXPECT_EQ(reordered[0].candidate_atom, adjacent[1].candidate_atom);
}

TEST(Chemotypes, InvalidSeedRingSystemReturnsEmptyResult) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("c1ccccc1-c1ccccc1"));

  int number_ring_systems;
  std::unique_ptr<int[]> ring_system = RingSystemMembership(m, number_ring_systems);
  ASSERT_EQ(number_ring_systems, 2);

  Set_of_Atoms embedding({0});
  EXPECT_TRUE(chemotypes::DirectlyAdjacentRingSystems(m, ring_system.get(), 0,
                                                      embedding).empty());
}

}  // namespace

TEST(Chemotypes, ComplementLabelsBothEndpointsBeyondRetainedAttachment) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCN1CCCCC1"));
  m.set_name("parent");
  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[N]");
  chemotypes::ChemotypeOptions options;
  options.include_attached_atoms = 1;
  options.isotope_for_exit_points = 99;
  chemotypes::ChemotypeScratch scratch;
  chemotypes::ChemotypeQueryMatch match;
  Molecule complement;
  ASSERT_EQ(chemotypes::ReduceToChemotype(m, queries, options, scratch, match,
                                        nullptr, &complement),
            chemotypes::ChemotypeQueryMatchStatus::kMatched);
  ASSERT_EQ(m.natoms(), 7);
  ASSERT_EQ(complement.natoms(), 1);
  EXPECT_EQ(complement.name(), "parent");
  EXPECT_EQ(complement.isotope(0), static_cast<isotope_t>(99));
  int labelled = 0;
  for (int i = 0; i < m.natoms(); ++i) {
    if (m.isotope(i) == 99) {
      ++labelled;
      EXPECT_EQ(m.atomic_number(i), 6);
      EXPECT_FALSE(m.is_ring_atom(i));
    }
  }
  EXPECT_EQ(labelled, 1);
}

TEST(Chemotypes, ComplementDisconnectedPiecesAndEmptyComplement) {
  for (const char* smiles : {"CN1CCC(C)CC1", "N1CCCCC1"}) {
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smiles));
    const int original_atoms = m.natoms();
    resizable_array_p<Substructure_Query> queries;
    AddSmarts(queries, "[N]");
    chemotypes::ChemotypeOptions options;
    options.isotope_for_exit_points = 99;
    chemotypes::ChemotypeScratch scratch;
    chemotypes::ChemotypeQueryMatch match;
    Molecule complement;
    ASSERT_EQ(chemotypes::ReduceToChemotype(m, queries, options, scratch, match,
                                          nullptr, &complement),
              chemotypes::ChemotypeQueryMatchStatus::kMatched);
    EXPECT_EQ(m.natoms() + complement.natoms(), original_atoms);
    int labelled = 0;
    for (int i = 0; i < m.natoms(); ++i) {
      labelled += m.isotope(i) == 99;
    }
    EXPECT_EQ(labelled, original_atoms - 6);
    for (int i = 0; i < complement.natoms(); ++i) {
      EXPECT_EQ(complement.isotope(i), static_cast<isotope_t>(99));
    }
    if (! complement.empty()) {
      EXPECT_EQ(complement.number_fragments(), 2);
    }
  }
}

TEST(Chemotypes, FixedExitIsotopeWithoutComplementOutput) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCN1CCCCC1"));
  resizable_array_p<Substructure_Query> queries;
  AddSmarts(queries, "[N]");
  chemotypes::ChemotypeOptions options;
  options.include_attached_atoms = 1;
  options.isotope_for_exit_points = 99;
  chemotypes::ChemotypeScratch scratch;
  chemotypes::ChemotypeQueryMatch match;
  ASSERT_EQ(chemotypes::ReduceToChemotype(m, queries, options, scratch, match),
            chemotypes::ChemotypeQueryMatchStatus::kMatched);
  int labelled = 0;
  for (int i = 0; i < m.natoms(); ++i) {
    labelled += m.isotope(i) == 99;
  }
  EXPECT_EQ(labelled, 1);
}
