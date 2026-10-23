#include "Molecule_Tools/chemotypes.h"

#include <algorithm>
#include <limits>
#include <optional>
#include <queue>
#include <utility>
#include <vector>

#include "Molecule_Lib/atom.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/bond.h"
#include "Molecule_Lib/substructure.h"

namespace chemotypes {
namespace {

int
SeedSystemDistance(Molecule& m, const int* ring_system, int seed_ring_system,
                   atom_number_t from, atom_number_t to) {
  if (from == to) {
    return 0;
  }

  const int matoms = m.natoms();
  std::vector<int> distance(matoms, -1);
  std::queue<atom_number_t> to_process;

  distance[from] = 0;
  to_process.push(from);

  while (! to_process.empty()) {
    const atom_number_t atom = to_process.front();
    to_process.pop();

    const int next_distance = distance[atom] + 1;
    for (const Bond* bond : m[atom]) {
      const atom_number_t other = bond->other(atom);
      if (ring_system[other] != seed_ring_system) {
        continue;
      }
      if (distance[other] >= 0) {
        continue;
      }
      if (other == to) {
        return next_distance;
      }

      distance[other] = next_distance;
      to_process.push(other);
    }
  }

  return std::numeric_limits<int>::max();
}

std::vector<int>
EmbeddingDistances(Molecule& m, const int* ring_system, int seed_ring_system,
                   const Set_of_Atoms& embedding, atom_number_t seed_atom) {
  std::vector<int> result;
  for (atom_number_t atom : embedding) {
    if (ring_system[atom] != seed_ring_system) {
      continue;
    }

    result.push_back(SeedSystemDistance(m, ring_system, seed_ring_system, atom,
                                        seed_atom));
  }

  return result;
}

std::pair<int, int>
FirstMatchedAtomPriority(const std::vector<int>& distances) {
  if (distances.empty()) {
    return {std::numeric_limits<int>::max(), std::numeric_limits<int>::max()};
  }

  return {0, distances[0]};
}

bool
Precedes(const AdjacentRingSystem& lhs, const AdjacentRingSystem& rhs,
         bool prefer_larger_adjacent_ring_system) {
  if (lhs.distance != rhs.distance) {
    return lhs.distance < rhs.distance;
  }
  if (prefer_larger_adjacent_ring_system &&
      lhs.atoms_down_the_bond != rhs.atoms_down_the_bond) {
    return lhs.atoms_down_the_bond > rhs.atoms_down_the_bond;
  }
  if (lhs.matched_atom_rank != rhs.matched_atom_rank) {
    return lhs.matched_atom_rank < rhs.matched_atom_rank;
  }
  if (lhs.distances_from_matched_atoms != rhs.distances_from_matched_atoms) {
    return lhs.distances_from_matched_atoms < rhs.distances_from_matched_atoms;
  }
  if (lhs.distance_from_matched_atom != rhs.distance_from_matched_atom) {
    return lhs.distance_from_matched_atom < rhs.distance_from_matched_atom;
  }
  if (lhs.seed_atom != rhs.seed_atom) {
    return lhs.seed_atom < rhs.seed_atom;
  }
  if (lhs.candidate_atom != rhs.candidate_atom) {
    return lhs.candidate_atom < rhs.candidate_atom;
  }

  return lhs.ring_system < rhs.ring_system;
}

void
ResizeAdjacentRingSystems(std::vector<AdjacentRingSystem>& adjacent, int requested,
                          int include_ties) {
  if (requested <= 0) {
    adjacent.clear();
    return;
  }

  if (requested >= static_cast<int>(adjacent.size())) {
    return;
  }

  if (! include_ties) {
    adjacent.resize(requested);
    return;
  }

  const int cutoff_distance = adjacent[requested - 1].distance;
  int keep = requested;
  while (keep < static_cast<int>(adjacent.size()) &&
         adjacent[keep].distance == cutoff_distance) {
    ++keep;
  }

  adjacent.resize(keep);
}

void
MaybeUpdate(std::vector<AdjacentRingSystem>& result, const AdjacentRingSystem& candidate,
            bool prefer_larger_adjacent_ring_system) {
  for (AdjacentRingSystem& existing : result) {
    if (existing.ring_system != candidate.ring_system) {
      continue;
    }

    if (Precedes(candidate, existing, prefer_larger_adjacent_ring_system)) {
      existing = candidate;
    }
    return;
  }

  result.push_back(candidate);
}

void
RecordCandidate(Molecule& m, const int* ring_system, int seed_ring_system,
                const Set_of_Atoms& embedding, atom_number_t seed_atom,
                atom_number_t candidate_atom, atom_number_t first_outside_atom,
                atom_number_t atom_before_candidate, int distance,
                bool prefer_larger_adjacent_ring_system,
                std::vector<AdjacentRingSystem>& result) {
  std::vector<int> distances = EmbeddingDistances(m, ring_system, seed_ring_system,
                                                  embedding, seed_atom);
  const auto [matched_atom_rank, distance_from_matched_atom] =
      FirstMatchedAtomPriority(distances);

  AdjacentRingSystem candidate;
  candidate.ring_system = ring_system[candidate_atom];
  candidate.distance = distance;
  candidate.seed_atom = seed_atom;
  candidate.candidate_atom = candidate_atom;
  candidate.first_outside_atom = first_outside_atom;
  candidate.atom_before_candidate = atom_before_candidate;
  if (atom_before_candidate != kInvalidAtomNumber) {
    std::vector<int> down_the_bond(m.natoms(), 0);
    const std::optional<int> atoms_down_the_bond =
        m.DownTheBond(atom_before_candidate, candidate_atom, down_the_bond.data());
    if (atoms_down_the_bond) {
      candidate.atoms_down_the_bond = *atoms_down_the_bond + distance - 1;
    }
  }
  candidate.matched_atom_rank = matched_atom_rank;
  candidate.distance_from_matched_atom = distance_from_matched_atom;
  candidate.distances_from_matched_atoms = std::move(distances);

  MaybeUpdate(result, candidate, prefer_larger_adjacent_ring_system);
}

void
FindAdjacentFromSeedAtom(Molecule& m, const int* ring_system, int seed_ring_system,
                         const Set_of_Atoms& embedding, atom_number_t seed_atom,
                         int max_distance, bool prefer_larger_adjacent_ring_system,
                         std::vector<AdjacentRingSystem>& result) {
  const int matoms = m.natoms();
  std::vector<int> distance(matoms, -1);
  std::vector<atom_number_t> first_outside_atom(matoms, kInvalidAtomNumber);
  std::queue<atom_number_t> to_process;

  for (const Bond* bond : m[seed_atom]) {
    const atom_number_t other = bond->other(seed_atom);
    if (ring_system[other] == seed_ring_system) {
      continue;
    }

    if (ring_system[other] > 0) {
      RecordCandidate(m, ring_system, seed_ring_system, embedding, seed_atom,
                      other, other, seed_atom, 1,
                      prefer_larger_adjacent_ring_system, result);
      continue;
    }

    distance[other] = 1;
    first_outside_atom[other] = other;
    to_process.push(other);
  }

  while (! to_process.empty()) {
    const atom_number_t atom = to_process.front();
    to_process.pop();

    if (max_distance > 0 && distance[atom] >= max_distance) {
      continue;
    }

    const int next_distance = distance[atom] + 1;
    for (const Bond* bond : m[atom]) {
      const atom_number_t other = bond->other(atom);
      if (ring_system[other] == seed_ring_system) {
        continue;
      }

      if (ring_system[other] > 0) {
        RecordCandidate(m, ring_system, seed_ring_system, embedding, seed_atom,
                        other, first_outside_atom[atom], atom, next_distance,
                        prefer_larger_adjacent_ring_system, result);
        continue;
      }

      if (distance[other] >= 0) {
        continue;
      }

      distance[other] = next_distance;
      first_outside_atom[other] = first_outside_atom[atom];
      to_process.push(other);
    }
  }
}


bool
IsTerminalExocyclicDoubleBondedAtom(Molecule& m, const Bond& bond,
                                    atom_number_t atom) {
  if (! bond.is_double_bond()) {
    return false;
  }
  if (m.ncon(atom) != 1) {
    return false;
  }

  const atomic_number_t z = m.atomic_number(atom);
  return z == 6 || z == 7 || z == 8;
}

int
SelectedRingSystemCount(const std::vector<int>& selected_ring_system) {
  return std::count(selected_ring_system.begin(), selected_ring_system.end(), 1);
}

void
AddLinkerComponents(Molecule& m, const std::vector<int>& ring_system,
                    const std::vector<int>& selected_ring_system,
                    std::vector<int>& mask, ChemotypeScratch& scratch) {
  const int matoms = m.natoms();
  scratch.PrepareForNewMatch(ring_system, selected_ring_system.size());

  for (atom_number_t atom = 0; atom < matoms; ++atom) {
    if (scratch.visited[atom]) {
      continue;
    }

    scratch.PrepareForNewLinkerComponent(selected_ring_system.size());
    scratch.to_process.push(atom);

    while (! scratch.to_process.empty()) {
      const atom_number_t current = scratch.to_process.front();
      scratch.to_process.pop();
      if (scratch.visited[current]) {
        continue;
      }
      scratch.visited[current] = 1;
      scratch.component.push_back(current);

      for (const Bond* bond : m[current]) {
        const atom_number_t other = bond->other(current);
        if (ring_system[other] > 0) {
          if (ring_system[other] < static_cast<int>(scratch.adjacent_ring_system.size()) &&
              selected_ring_system[ring_system[other]]) {
            scratch.adjacent_ring_system[ring_system[other]] = 1;
          }
          continue;
        }

        if (! scratch.visited[other]) {
          scratch.to_process.push(other);
        }
      }
    }

    if (SelectedRingSystemCount(scratch.adjacent_ring_system) < 2) {
      continue;
    }

    for (atom_number_t a : scratch.component) {
      mask[a] = kChemotypeCoreAtom;
    }
  }
}

void
AddExocyclicDoubleBondedAtoms(Molecule& m, const std::vector<int>& ring_system,
                              std::vector<int>& mask) {
  const int matoms = m.natoms();
  for (atom_number_t atom = 0; atom < matoms; ++atom) {
    if (mask[atom] || m.ncon(atom) != 1) {
      continue;
    }

    const Atom& a = m[atom];
    const Bond* bond = a[0];
    if (! IsTerminalExocyclicDoubleBondedAtom(m, *bond, atom)) {
      continue;
    }

    const atom_number_t other = bond->other(atom);
    if (mask[other]) {
      mask[atom] = kChemotypeCoreAtom;
    }
  }
}

void
AddAttachedAtoms(Molecule& m, const std::vector<int>& ring_system,
                 const ChemotypeOptions& options, std::vector<int>& mask) {
  const int matoms = m.natoms();

  for (atom_number_t atom = 0; atom < matoms; ++atom) {
    if (! mask[atom] || ring_system[atom] == 0) {
      continue;
    }

    for (const Bond* bond : m[atom]) {
      const atom_number_t other = bond->other(atom);
      if (mask[other]) {
        continue;
      }
      if (ring_system[other] > 0) {
        continue;
      }
      if (IsTerminalExocyclicDoubleBondedAtom(m, *bond, other)) {
        continue;
      }
      if (options.ignore_singly_connected_attached_atoms && m.ncon(other) == 1) {
        continue;
      }

      mask[other] = kChemotypeAttachedAtom;
    }
  }
}


void
ApplyExitPointIsotopes(Molecule& m, const std::vector<int>& keep,
                       isotope_t isotope) {
  if (isotope == 0) {
    return;
  }

  // Label both original endpoints before extracting either subset. The final
  // mask includes optional attachment atoms, so cuts need not occur at rings.
  const int matoms = m.natoms();
  for (atom_number_t atom = 0; atom < matoms; ++atom) {
    if (! keep[atom]) {
      continue;
    }
    for (const Bond* bond : m[atom]) {
      const atom_number_t other = bond->other(atom);
      if (! keep[other]) {
        m.set_isotope(atom, isotope);
        m.set_isotope(other, isotope);
      }
    }
  }
}

int
ApplyAtomTypeIsotopes(Molecule& m, Atom_Typing_Specification& atom_typing) {
  const int matoms = m.natoms();
  std::vector<uint32_t> atom_type(matoms, 0);
  if (! atom_typing.assign_atom_types(m, atom_type.data())) {
    return 0;
  }

  for (atom_number_t atom = 0; atom < matoms; ++atom) {
    m.set_isotope(atom, static_cast<isotope_t>(atom_type[atom]));
  }

  return 1;
}

// Atom type labelling is via attached atoms. Core chemotype atoms lose
// their isotopes, while one-hop attachment atoms retain atom-type isotopes.
void
UnsetIsotopesAddAdjacent(Molecule& m, std::vector<int>& keep) {
  const int matoms = m.natoms();
  for (atom_number_t atom = 0; atom < matoms; ++atom) {
    // Singly connected atoms never have isotopic labels.
    if (m.ncon(atom) == 1) {
      m.set_isotope(atom, 0);
      continue;
    }

    if (keep[atom] != kChemotypeCoreAtom) {
      continue;
    }

    m.set_isotope(atom, 0);

    for (const Bond* b : m[atom]) {
      const atom_number_t other = b->other(atom);
      if (keep[other] == kChemotypeNotKept) {
        keep[other] = kChemotypeAttachedAtom;
      }
    }
  }
}

}  // namespace



void
ChemotypeScratch::PrepareForNewMatch(const std::vector<int>& ring_system,
                                     int ring_system_count) {
  visited.resize(ring_system.size());
  for (size_t i = 0; i < ring_system.size(); ++i) {
    visited[i] = ring_system[i] > 0;
  }

  adjacent_ring_system.resize(ring_system_count);
  component.clear();
  while (! to_process.empty()) {
    to_process.pop();
  }
}

void
ChemotypeScratch::PrepareForNewLinkerComponent(int ring_system_count) {
  component.clear();
  adjacent_ring_system.assign(ring_system_count, 0);
  while (! to_process.empty()) {
    to_process.pop();
  }
}

int
SetMatchFromEmbedding(int query_index, int hits, const Set_of_Atoms& embedding,
                      const std::vector<int>& ring_system, int required_ring_system,
                      ChemotypeQueryMatch& result) {
  for (atom_number_t atom : embedding) {
    if (ring_system[atom] == 0) {
      continue;
    }
    if (required_ring_system > 0 && ring_system[atom] != required_ring_system) {
      continue;
    }

    result.query_index = query_index;
    result.hits = hits;
    result.embedding = embedding;
    result.ring_system = ring_system;
    result.seed_atom = atom;
    result.seed_ring_system = ring_system[atom];
    return 1;
  }

  return 0;
}

ChemotypeQueryMatchStatus
FirstChemotypeQueryMatch(Molecule& m, resizable_array_p<Substructure_Query>& queries,
                         ChemotypeQueryMatch& result, bool choose_first_embedding) {
  result = ChemotypeQueryMatch();

  for (int i = 0; i < queries.number_elements(); ++i) {
    Substructure_Query* query = queries[i];
    if (query == nullptr) {
      continue;
    }

    Substructure_Results sresults;
    const int hits = query->substructure_search(m, sresults);
    if (hits == 0) {
      continue;
    }

    std::vector<int> ring_system(m.natoms(), 0);
    m.label_atoms_by_ring_system_including_spiro_fused(ring_system.data());

    result.query_index = i;
    result.hits = hits;
    result.embedding = *sresults.embedding(0);
    result.ring_system = ring_system;

    if (choose_first_embedding) {
      if (SetMatchFromEmbedding(i, hits, *sresults.embedding(0), ring_system, 0, result)) {
        return ChemotypeQueryMatchStatus::kMatched;
      }
      return ChemotypeQueryMatchStatus::kMatchedQueryNoRingAtom;
    }

    std::vector<int> ring_system_seen(m.natoms() + 1, 0);
    int unique_ring_systems = 0;
    int only_ring_system = 0;
    for (const Set_of_Atoms* embedding : sresults.embeddings()) {
      for (atom_number_t atom : *embedding) {
        const int rs = ring_system[atom];
        if (rs == 0 || ring_system_seen[rs]) {
          continue;
        }
        ring_system_seen[rs] = 1;
        only_ring_system = rs;
        ++unique_ring_systems;
        if (unique_ring_systems > 1) {
          return ChemotypeQueryMatchStatus::kAmbiguousQueryMatches;
        }
      }
    }

    if (unique_ring_systems == 0) {
      return ChemotypeQueryMatchStatus::kMatchedQueryNoRingAtom;
    }

    for (const Set_of_Atoms* embedding : sresults.embeddings()) {
      if (SetMatchFromEmbedding(i, hits, *embedding, ring_system, only_ring_system, result)) {
        return ChemotypeQueryMatchStatus::kMatched;
      }
    }

    return ChemotypeQueryMatchStatus::kMatchedQueryNoRingAtom;
  }

  return ChemotypeQueryMatchStatus::kNoQueryMatch;
}

std::vector<AdjacentRingSystem>
DirectlyAdjacentRingSystems(Molecule& m, const int* ring_system, int seed_ring_system,
                            const Set_of_Atoms& embedding, int max_distance,
                            bool prefer_larger_adjacent_ring_system) {
  std::vector<AdjacentRingSystem> result;
  if (seed_ring_system <= 0) {
    return result;
  }

  const int matoms = m.natoms();
  for (atom_number_t atom = 0; atom < matoms; ++atom) {
    if (ring_system[atom] != seed_ring_system) {
      continue;
    }
    FindAdjacentFromSeedAtom(m, ring_system, seed_ring_system, embedding, atom,
                             max_distance, prefer_larger_adjacent_ring_system, result);
  }

  std::sort(result.begin(), result.end(),
            [prefer_larger_adjacent_ring_system](const AdjacentRingSystem& lhs,
                                                 const AdjacentRingSystem& rhs) {
              return Precedes(lhs, rhs, prefer_larger_adjacent_ring_system);
            });
  return result;
}


std::vector<int>
ChemotypeAtomMask(Molecule& m, const ChemotypeQueryMatch& match,
                  const ChemotypeOptions& options, ChemotypeScratch& scratch) {
  const int matoms = m.natoms();
  std::vector<int> result(matoms, kChemotypeNotKept);

  if (match.seed_ring_system <= 0 ||
      match.ring_system.size() != static_cast<size_t>(matoms)) {
    return result;
  }

  int max_ring_system = match.seed_ring_system;
  for (int r : match.ring_system) {
    if (r > max_ring_system) {
      max_ring_system = r;
    }
  }

  std::vector<int> selected_ring_system(max_ring_system + 1, 0);
  selected_ring_system[match.seed_ring_system] = 1;

  std::vector<AdjacentRingSystem> adjacent = DirectlyAdjacentRingSystems(
      m, match.ring_system.data(), match.seed_ring_system, match.embedding,
      options.max_distance, options.prefer_larger_adjacent_ring_system);

  const int adjacent_ring_systems_to_include =
      std::max(0, options.adjacent_ring_systems_to_include);
  ResizeAdjacentRingSystems(adjacent, adjacent_ring_systems_to_include,
                            options.include_tied_adjacent_ring_systems);

  for (const AdjacentRingSystem& ring_system : adjacent) {
    if (ring_system.ring_system >= 0 &&
        ring_system.ring_system < static_cast<int>(selected_ring_system.size())) {
      selected_ring_system[ring_system.ring_system] = 1;
    }
  }

  for (atom_number_t atom = 0; atom < matoms; ++atom) {
    const int r = match.ring_system[atom];
    if (r > 0 && r < static_cast<int>(selected_ring_system.size()) &&
        selected_ring_system[r]) {
      result[atom] = kChemotypeCoreAtom;
    }
  }

  if (options.include_linker_atoms) {
    AddLinkerComponents(m, match.ring_system, selected_ring_system, result, scratch);
  }

  if (options.include_exocyclic_double_bonded_atoms) {
    AddExocyclicDoubleBondedAtoms(m, match.ring_system, result);
  }

  if (options.include_attached_atoms) {
    AddAttachedAtoms(m, match.ring_system, options, result);
  }

  return result;
}


std::vector<int>
ChemotypeAtomMask(Molecule& m, const ChemotypeQueryMatch& match,
                  const ChemotypeOptions& options) {
  ChemotypeScratch scratch;
  return ChemotypeAtomMask(m, match, options, scratch);
}

std::vector<int>
ChemotypeAtomMask(Molecule& m, const ChemotypeQueryMatch& match,
                  int adjacent_ring_systems_to_include, int max_distance) {
  ChemotypeOptions options;
  options.adjacent_ring_systems_to_include = adjacent_ring_systems_to_include;
  options.max_distance = max_distance;
  return ChemotypeAtomMask(m, match, options);
}


ChemotypeQueryMatchStatus
ReduceToChemotype(Molecule& m, resizable_array_p<Substructure_Query>& queries,
                  const ChemotypeOptions& options, ChemotypeScratch& scratch,
                  ChemotypeQueryMatch& match,
                  Atom_Typing_Specification* atom_typing, Molecule* complement) {
  if (complement != nullptr) {
    *complement = Molecule();
  }
  if (atom_typing != nullptr && ! ApplyAtomTypeIsotopes(m, *atom_typing)) {
    return ChemotypeQueryMatchStatus::kAtomTypingFailed;
  }

  if (const ChemotypeQueryMatchStatus status = FirstChemotypeQueryMatch(
          m, queries, match, options.choose_first_embedding);
      status != ChemotypeQueryMatchStatus::kMatched) {
    return status;
  }

  std::vector<int> keep = ChemotypeAtomMask(m, match, options, scratch);
  if (atom_typing != nullptr) {
    UnsetIsotopesAddAdjacent(m, keep);
  } else if (options.isotope_for_exit_points != 0) {
    ApplyExitPointIsotopes(m, keep, options.isotope_for_exit_points);
  }

  if (complement != nullptr) {
    *complement = m;
    std::vector<int> outside(keep.size());
    std::transform(keep.begin(), keep.end(), outside.begin(),
                   [](int value) { return value == kChemotypeNotKept; });
    complement->remove_atoms(outside.data(), 0);
  }
  m.remove_atoms(keep.data(), 0);
  return ChemotypeQueryMatchStatus::kMatched;
}

}  // namespace chemotypes
