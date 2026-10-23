#ifndef MOLECULE_TOOLS_CHEMOTYPES_H_
#define MOLECULE_TOOLS_CHEMOTYPES_H_

#include <queue>
#include <vector>

#include "Foundational/iwaray/iwaray.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/set_of_atoms.h"

class Atom_Typing_Specification;
class Substructure_Query;

namespace chemotypes {

// Values in the atom mask returned by ChemotypeAtomMask. Nonzero values are
// retained by Molecule::remove_atoms(mask.data(), 0).
static constexpr int kChemotypeNotKept = 0;
static constexpr int kChemotypeCoreAtom = 1;
static constexpr int kChemotypeAttachedAtom = 2;

enum class ChemotypeQueryMatchStatus {
  kMatched,
  kNoQueryMatch,
  kMatchedQueryNoRingAtom,
  kAmbiguousQueryMatches,
  kAtomTypingFailed,
};

struct ChemotypeOptions {
  int adjacent_ring_systems_to_include = 0;
  int max_distance = 0;
  int include_linker_atoms = 1;
  int include_exocyclic_double_bonded_atoms = 1;
  int include_attached_atoms = 0;
  int ignore_singly_connected_attached_atoms = 0;
  int include_tied_adjacent_ring_systems = 0;
  bool prefer_larger_adjacent_ring_system = false;
  bool choose_first_embedding = false;
  isotope_t isotope_for_exit_points = 0;
};


struct ChemotypeScratch {
  std::vector<int> visited;
  std::vector<atom_number_t> component;
  std::vector<int> adjacent_ring_system;
  std::queue<atom_number_t> to_process;

  void PrepareForNewMatch(const std::vector<int>& ring_system,
                          int ring_system_count);
  void PrepareForNewLinkerComponent(int ring_system_count);
};

struct ChemotypeQueryMatch {
  int query_index = -1;
  int hits = 0;
  Set_of_Atoms embedding;
  atom_number_t seed_atom = kInvalidAtomNumber;
  int seed_ring_system = 0;
  std::vector<int> ring_system;
};

// Scan queries in order and return a chemotype-defining match. By default, the
// first query that matches must identify exactly one ring system across all its
// embeddings. With choose_first_embedding, the first embedding of the first
// matching query is used even if the query matches multiple ring systems.
ChemotypeQueryMatchStatus FirstChemotypeQueryMatch(
    Molecule& m, resizable_array_p<Substructure_Query>& queries,
    ChemotypeQueryMatch& result, bool choose_first_embedding = false);

struct AdjacentRingSystem {
  int ring_system = 0;
  int distance = 0;
  atom_number_t seed_atom = kInvalidAtomNumber;
  atom_number_t candidate_atom = kInvalidAtomNumber;
  int matched_atom_rank = 0;
  int distance_from_matched_atom = 0;
  atom_number_t first_outside_atom = kInvalidAtomNumber;
  atom_number_t atom_before_candidate = kInvalidAtomNumber;
  int atoms_down_the_bond = 0;

  // Distances from seed_atom to matched atoms in embedding order. This lets
  // callers make deterministic choices in symmetric cases by query atom order.
  std::vector<int> distances_from_matched_atoms;
};

// Return ring systems directly adjacent to `seed_ring_system`. Paths may start
// from any atom in `seed_ring_system`, but internal atoms must be non-ring
// atoms. The traversal stops when it reaches a different ring system, so ring
// systems beyond an intervening ring system are not adjacent.
//
// The returned systems are sorted by chemical closeness, then by the order of
// atoms in `embedding` so query atom order controls ambiguous choices. If
// `prefer_larger_adjacent_ring_system` is true, same-distance candidates are
// ordered by the number of atoms on the candidate side, including linker atoms
// from the seed ring system to the candidate and atoms down the bond entering
// the candidate ring system, before deterministic query/atom order tie-breakers.
std::vector<AdjacentRingSystem> DirectlyAdjacentRingSystems(
    Molecule& m, const int* ring_system, int seed_ring_system,
    const Set_of_Atoms& embedding, int max_distance = 0,
    bool prefer_larger_adjacent_ring_system = false);


// Return an atom mask for the chemotype core. The mask includes all atoms in
// the seed ring system and, if requested, the first N directly adjacent ring
// systems according to DirectlyAdjacentRingSystems. Linker atoms, exocyclic
// double-bonded atoms and one-hop attachment atoms are controlled by options.
std::vector<int> ChemotypeAtomMask(Molecule& m, const ChemotypeQueryMatch& match,
                                   const ChemotypeOptions& options,
                                   ChemotypeScratch& scratch);

std::vector<int> ChemotypeAtomMask(Molecule& m, const ChemotypeQueryMatch& match,
                                   const ChemotypeOptions& options);

std::vector<int> ChemotypeAtomMask(Molecule& m, const ChemotypeQueryMatch& match,
                                   int adjacent_ring_systems_to_include,
                                   int max_distance = 0);


// Find the first matching query, build the chemotype atom mask, optionally
// label terminal attachment atoms by atom type, and remove atoms outside the
// mask from `m`. Atom typing, when specified, is applied before any atoms are
// removed so full-molecule context is available. If supplied, complement receives
// the excluded atoms with the parent name; it is empty on a failed match or when
// all atoms are retained. Fixed exit isotopes label both endpoints of each cut.
ChemotypeQueryMatchStatus ReduceToChemotype(
    Molecule& m, resizable_array_p<Substructure_Query>& queries,
    const ChemotypeOptions& options, ChemotypeScratch& scratch,
    ChemotypeQueryMatch& match, Atom_Typing_Specification* atom_typing = nullptr,
    Molecule* complement = nullptr);

}  // namespace chemotypes

#endif  // MOLECULE_TOOLS_CHEMOTYPES_H_
