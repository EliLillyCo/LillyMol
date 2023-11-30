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

  // For each query the number of times it has matched a molecule.
  std::vector<uint32_t> query_matched;

  bool must_match_all_queries = false;

  // When LabelMatchedAtoms is called with a string, what kind of smiles
  // string is returned.
  bool labeled_smiles_are_unique = false;

  // When LabelMatchedAtoms is called, we can number the atoms
  // with the (query_number + 1) that hits - must be offset since
  // isotope 0 is invisible.
  bool label_by_query_number = false;

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
  bool AddQueriesFromSmarts(const std::vector<std::string>& smarts);

  void set_unique_embeddings_only(bool s) {
    unique_embeddings_only = s;
    ApplyQueryConstraints();
  }

  void set_find_one_embedding_per_root_atom(bool s) {
    find_one_embedding_per_root_atom = s;
    ApplyQueryConstraints();
  }

  void set_perceive_symmetry_equivalent_matches(bool s) {
    perceive_symmetry_equivalent_matches = s;
    ApplyQueryConstraints();
  }

  void set_max_matches_to_find(int s) {
    max_matches_to_find = s;
    ApplyQueryConstraints();
  }

  void set_reduce_to_largest_fragment(bool s) {
    reduce_to_largest_fragment = s;
  }

  void set_make_implicit_hydrogens_explicit(bool s) {
    make_implicit_hydrogens_explicit = s;
  }

  void set_labeled_smiles_are_unique(bool s) {
    labeled_smiles_are_unique = s;
  }

  void set_label_by_query_number(bool s) {
    label_by_query_number = s;
  }

  const std::vector<uint32_t> QueryMatched() const {
    return query_matched;
  }

  uint32_t number_queries() const {
    return query.size();
  }

  std::vector<std::string> query_names() const;

  // Return true if there is a match.
  // If `must_match_all_queries` is true, then all queries must
  // be matched. Otherwise return true if any query matches.
  bool SubstructureSearch(const std::string& smi);
  bool SubstructureSearch(Molecule& mol);

  // Return a vector containing the number of times each query
  // matches the input.
  std::vector<int> NumofMatches(const std::string& smiles);
  std::vector<int> NumofMatches(Molecule& m);

  // Perform substructure queries on `mol` and apply isotopic label
  // `isotope` to the matched atoms. 
  std::string LabelMatchedAtoms(const std::string& smi);
  // The matched atoms are marked with isotopic labels, and the number
  // of queries that match is returned.
  int LabelMatchedAtoms(Molecule& m);

  // For each molecule return the number of matches in each query.
  std::vector<std::vector<uint32_t>> NumberMatches(const std::vector<std::string>& smiles);
  // The argument is non-const because preprocessing may change the molecule - strip
  // to largest fragment, add implicit hydrogens ...
  // and substructure searching a molecule is non-const.
  std::vector<std::vector<uint32_t>> NumberMatches(std::vector<Molecule>& mols);

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
