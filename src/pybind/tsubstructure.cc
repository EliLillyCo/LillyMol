// Substructure Query object for pybond

#include "tsubstructure.h"

#include <iostream>

#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/target.h"

namespace pybind_substructure {

using std::cerr;

TSubstructure::TSubstructure() {
}

TSubstructure::~TSubstructure() {
  for (Substructure_Query* q : query) {
    delete q;
  }
}

bool
TSubstructure::ReadQueries(std::string& directive) {
  static constexpr int kVerbose = 1;

  const const_IWSubstring tmp(directive);
  resizable_array_p<Substructure_Query> local_queries;
  if (!process_cmdline_token('*', tmp, local_queries, kVerbose)) {
    cerr << "Cannot process '" << directive << "'\n";
    return false;
  }

  for (Substructure_Query* q : local_queries) {
    query.push_back(q);
  }

  local_queries.resize_no_delete(0);

  ApplyQueryConstraints();

  query_matched.resize(query.size(), 0);

  return true;
}

// For all queries in `query` apply all known query constraints.
void
TSubstructure::ApplyQueryConstraints() {
  for (Substructure_Query* q : query) {
    if (unique_embeddings_only) {
      q->set_find_unique_embeddings_only(1);
    }
    if (find_one_embedding_per_root_atom) {
      q->set_find_one_embedding_per_atom(1);
    }
    if (!perceive_symmetry_equivalent_matches) {
      q->set_perceive_symmetry_equivalent_matches(0);
    }
    if (max_matches_to_find > 0) {
      q->set_max_matches_to_find(max_matches_to_find);
    }
  }
}

bool
TSubstructure::AddQueryFromSmarts(const std::string& smarts) {
  std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
  if (!q->CreateFromSmarts(smarts)) {
    cerr << "TSubstructure::AddQueryFromSmarts:invalid smarts '" << smarts << "'\n";
    return false;
  }

  query.push_back(q.release());

  query_matched.resize(query.size(), 0);

  // Only the last query read needs to be set.
  ApplyQueryConstraints();

  return true;
}

bool
TSubstructure::AddQueriesFromSmarts(const std::vector<std::string>& smarts) {
  for (const std::string& smt : smarts) {
    if (! AddQueryFromSmarts(smt)) {
      cerr << "TSubstructure::AddQueriesFromSmarts:invalid smarts '" << smt << "'\n";
      return false;
    }
  }

  return true;
}

int
TSubstructure::Preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }
  if (make_implicit_hydrogens_explicit) {
    m.make_implicit_hydrogens_explicit();
  }

  return 1;
}

bool
TSubstructure::SubstructureSearch(Molecule& m) {
  Preprocess(m);

  Molecule_to_Match target(&m);

  if (query.empty()) {
    cerr << "No query has been loaded\n";
    return false;
  }

  if (must_match_all_queries) {
    return SubstructureSearchAllMatch(m);
  } else {
    return SubstructureSearchAnyMatch(m);
  }
}

bool
TSubstructure::SubstructureSearchAnyMatch(Molecule& m) {
  Preprocess(m);

  Molecule_to_Match target(&m);

  for (uint32_t i = 0; i < query.size(); ++i) {
    if (query[i]->substructure_search(target)) {
      ++query_matched[i];
      return true;
    }
  }
  return false;
}

bool
TSubstructure::SubstructureSearchAllMatch(Molecule& m) {
  Preprocess(m);
  Molecule_to_Match target(&m);

  for (uint32_t i = 0; i < query.size(); ++i) {
    if (query[i]->substructure_search(target)) {
      ++query_matched[i];
    } else {
      return false;
    }
  }

  return true;
}

bool
TSubstructure::SubstructureSearch(const std::string& smiles) {
  Molecule m;
  if (!m.build_from_smiles(smiles)) {
    std::cerr << "Invalid smiles '" << smiles << "'\n";
    return 0;
  }

  return SubstructureSearch(m);
}

std::vector<int>
TSubstructure::NumofMatches(const std::string& smiles) {
  Molecule m;
  if (!m.build_from_smiles(smiles)) {
    std::cerr << "Invalid smiles '" << smiles << "'\n";
  }

  return NumofMatches(m);
}

std::vector<int>
TSubstructure::NumofMatches(Molecule& m) {
  Preprocess(m);

  Molecule_to_Match target(&m);
  std::vector<int> res;

  const uint32_t nq = query.size();
  res.reserve(nq);
  for (uint32_t i = 0; i < nq; ++i) {
    int nhits = query[i]->substructure_search(target);
    res.push_back(nhits);
    if (nhits) {
      ++query_matched[i];
    }
  }

  return res;
}

std::string
TSubstructure::LabelMatchedAtoms(const std::string& smi) {
  if (isotope == 0) {
    cerr << "TSubstructure::LabelMatchedAtoms:isotope not specified, set to 1\n";
    isotope = 1;
  }

  Molecule m;
  if (! m.build_from_smiles(smi)) {
    cerr << "TSubstructure::LabelMatchedAtoms:invalid smiles '" << smi << "'\n";
    return std::string("");
  }

  LabelMatchedAtoms(m);

  return m.smiles().AsString();
}

int
TSubstructure::LabelMatchedAtoms(Molecule& m) {
  Preprocess(m);

  int got_result = 0;

  Molecule_to_Match target(&m);
  for (uint32_t i = 0; i < query.size(); ++i) {
    Substructure_Results sresults;
    if (query[i]->substructure_search(target, sresults) == 0) {
      continue;
    }

    ++query_matched[i];

    ++got_result;

    // Loop through all the embeddings, and through each atom
    // in the embedding and set the isotope.
    for (const Set_of_Atoms* e : sresults.embeddings()) {
      for (atom_number_t a : *e) {
        if (label_by_query_number) {
          m.set_isotope(a, i + 1);
        } else {
          m.set_isotope(a, isotope);
        }
      }
    }
  }

  return got_result;
}

std::vector<std::vector<uint32_t>>
TSubstructure::MatchedAtoms(const std::string& smi) {
  Molecule m;
  if (! m.build_from_smiles(smi)) {
    cerr << "TSubstructure::MatchedAtoms:invalid smiles '" << smi << "'\n";
    return std::vector<std::vector<uint32_t>>(query.size());
  }

  return MatchedAtoms(m);
}

std::vector<std::vector<uint32_t>>
TSubstructure::MatchedAtoms(Molecule& m) {
  std::vector<std::vector<uint32_t>> result(query.size());

  Preprocess(m);

  Molecule_to_Match target(&m);

  const int nqueries = query.size();
  for (int i = 0; i < nqueries; ++i) {
    Substructure_Results sresults;
    if (query[i]->substructure_search(target, sresults) == 0) {
      continue;
    }

    ++query_matched[i];

    // Loop through all the embeddings, and through each atom
    // in the embedding and set the isotope.
    for (const Set_of_Atoms* e : sresults.embeddings()) {
      for (atom_number_t a : *e) {
        result[i].push_back(a);
      }
    }
  }

  return result;
}

std::vector<std::vector<uint32_t>>
TSubstructure::NumberMatches(const std::vector<std::string>& smiles) {
  const uint32_t nmols = smiles.size();
  std::vector<Molecule> mols(nmols);
  for (uint32_t i = 0; i < nmols; ++i) {
    if (! mols[i].build_from_smiles(smiles[i])) {
      cerr << "TSubstructure::NumberMatches:invalid smiles '" << smiles[i] << "'\n";
      mols[i].resize(0);
    }
  }

  return NumberMatches(mols);
}

std::vector<std::vector<uint32_t>>
TSubstructure::NumberMatches(std::vector<Molecule>& mols) {
  const uint32_t nmols = mols.size();

  std::vector<std::vector<uint32_t>> result(nmols);

  const int nq = query.size();

  for (uint32_t i = 0; i < nmols; ++i) {
    Preprocess(mols[i]);

    Molecule_to_Match target(&(mols[i]));

    result[i].reserve(nq);
    for (int j = 0; j < nq; ++j) {
      const int nhits = query[j]->substructure_search(target);
      result[i].push_back(nhits);
      if (nhits) {
        ++query_matched[j];
      }
    }
  }

  return result;
}

std::vector<std::string>
TSubstructure::query_names() const {
  std::vector<std::string> result(query.size());
  for (const Substructure_Query* q : query) {
    result.push_back(q->comment().AsString());
  }

  return result;
}


}  // namespace pybind_substructure
