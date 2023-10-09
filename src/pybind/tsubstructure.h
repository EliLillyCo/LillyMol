#ifndef PYBIND_TSUBSTRUCTURE_H_
#define PYBIND_TSUBSTRUCTURE_H_

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"

namespace pybind_substructure {

struct TSubstructure {
  bool reduce_to_largest_fragment = false;

  bool unique_embeddings_only = false;
  bool find_one_embedding_per_root_atom = false;
  bool perceive_symmetry_equivalent_matches = true;
  uint32_t max_matches_to_find = 0;
  bool make_implicit_hydrogens_explicit = false;

  isotope_t isotope = 0;

  std::vector<Substructure_Query*> query;

  bool must_match_all_queries = false;

  // Private functions
  void ApplyQueryConstraints();

  int Preprocess(Molecule& m);

  bool SubstructureSearchAllMatch(Molecule& m);
  bool SubstructureSearchAnyMatch(Molecule& m);

 public:
  TSubstructure();
  ~TSubstructure();

  // Get queries, same as the -q option to tsubstructure.
  bool ReadQueries(std::string& directive);

  bool AddQueryFromSmarts(const std::string& smarts);

  // Return true if there is a match.
  // If `must_match_all_queries` is true, then all queries must
  // be matched. Otherwise return true if any query matches.
  bool SubstructureSearch(const std::string& smi);
  bool SubstructureSearch(Molecule& mol);

  // Return a vector containing the number of times each query
  // matches the input.
  std::vector<int> NumofMatches(const std::string& smiles);
  std::vector<int> NumofMatches(Molecule& m);

  // Return an isotopically labelled smiles where each matched
  // atom, across all queries, will be isotopically labelled `isotope`.
  // Note that if no queries match, an empty string is returned.
  std::string LabelMatchedAtoms(const std::string& smi);
  std::string LabelMatchedAtoms(Molecule& m);

  // For each query, lists of the matched atoms.
  // If a particular query does not match the input, an empty vector is
  // returned for that query.
  // for example if the queries are "C" and "N", passing "CC" will
  // result in [ [0, 1], [] ] as the result.
  // Note that what is returned is very similar to what is labelled
  // by LabelMatchedAtoms.
  // Note that the same atom may appear multiple times in a given
  // match.
  std::vector<std::vector<uint32_t>> MatchedAtoms(const std::string& smi);
  std::vector<std::vector<uint32_t>> MatchedAtoms(Molecule& m);
};

};  // namespace pybind_substructure

#endif  // PYBIND_TSUBSTRUCTURE_H_
