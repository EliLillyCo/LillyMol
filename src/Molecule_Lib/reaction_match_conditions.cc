#include <stdlib.h>

#include <iostream>
#include <limits>

#include "reaction_match_conditions.h"

using std::cerr;

// Don't forget to update the assignment operator if new values added.
// Also IsDefault will need to be updated.

static constexpr int default_find_unique_embeddings = 1;
static constexpr int default_ignore_not_reacting = 0;
static constexpr int default_process_hit_number = -1;
static constexpr int default_max_matches_to_find = 0;
static constexpr int default_one_embedding_per_start_atom = 0;
static constexpr int default_ignore_symmetry_related_matches = 0;
static constexpr int default_embeddings_can_overlap = 1;
static constexpr int default_suppress_if_more_than_this_many_substructure_search_hits = std::numeric_limits<int>::max();
static constexpr int default_issue_sidechain_no_match_warnings = 1;

Match_Conditions::Match_Conditions() {
  _verbose = 0;

  _find_unique_embeddings = default_find_unique_embeddings;

  _ignore_not_reacting = default_ignore_not_reacting;

  _process_hit_number = default_process_hit_number;

  _max_matches_to_find = default_max_matches_to_find;

  _one_embedding_per_start_atom = default_one_embedding_per_start_atom;

  _ignore_symmetry_related_matches = default_ignore_symmetry_related_matches;

  _embeddings_can_overlap = default_embeddings_can_overlap;

  _suppress_if_more_than_this_many_substructure_search_hits =
    default_suppress_if_more_than_this_many_substructure_search_hits;

  _issue_sidechain_no_match_warnings = default_issue_sidechain_no_match_warnings;

  return;
}

int
Match_Conditions::IsDefault() const {
  if (_find_unique_embeddings != default_find_unique_embeddings) {
    return 0;
  }
  if (_ignore_not_reacting != default_ignore_not_reacting) {
    return 0;
  }
  if (_process_hit_number != default_process_hit_number) {
    return 0;
  }
  if (_max_matches_to_find != default_max_matches_to_find) {
    return 0;
  }
  if (_one_embedding_per_start_atom != default_one_embedding_per_start_atom) {
    return 0;
  }
  if (_ignore_symmetry_related_matches != default_ignore_symmetry_related_matches) {
    return 0;
  }
  if (_embeddings_can_overlap != default_embeddings_can_overlap) {
    return 0;
  }
  if (_suppress_if_more_than_this_many_substructure_search_hits != 
                default_suppress_if_more_than_this_many_substructure_search_hits) {
    return 0;
  }
  if (_issue_sidechain_no_match_warnings != default_issue_sidechain_no_match_warnings) {
    return 0;
  }

  return 1;
}

Match_Conditions&
Match_Conditions::operator=(const Match_Conditions& rhs) {
  _verbose = rhs._verbose;

  _find_unique_embeddings = rhs._find_unique_embeddings;

  _ignore_not_reacting = rhs._ignore_not_reacting;

  _process_hit_number = rhs._process_hit_number;

  _one_embedding_per_start_atom = rhs._one_embedding_per_start_atom;

  _ignore_symmetry_related_matches = rhs._ignore_symmetry_related_matches;

  _embeddings_can_overlap = rhs._embeddings_can_overlap;

  _multiple_match_string = rhs._multiple_match_string;

  return *this;
}

static constexpr int default_make_new_reagent_for_each_hit = 0;
static constexpr int default_strip_reagents_to_largest_fragment = 0;

Sidechain_Match_Conditions::Sidechain_Match_Conditions() {
  _make_new_reagent_for_each_hit = default_make_new_reagent_for_each_hit;

  _strip_reagents_to_largest_fragment = default_strip_reagents_to_largest_fragment;

  _active = 0;

  return;
}

int
Sidechain_Match_Conditions::IsDefault() const {
  if (! Match_Conditions::IsDefault()) {
    return 0;
  }

  if (_make_new_reagent_for_each_hit != default_make_new_reagent_for_each_hit) {
    return 0;
  }
  if (_strip_reagents_to_largest_fragment != default_strip_reagents_to_largest_fragment) {
    return 0;
  }

  return 1;
}

Sidechain_Match_Conditions&
Sidechain_Match_Conditions::operator=(const Sidechain_Match_Conditions& rhs) {
  Match_Conditions::operator=(rhs);

  _make_new_reagent_for_each_hit = rhs._make_new_reagent_for_each_hit;

  _max_matches_to_find = rhs._max_matches_to_find;

  _strip_reagents_to_largest_fragment = rhs._strip_reagents_to_largest_fragment;

  return *this;
}

static constexpr int default_enumerate_scaffold_hits_individually = 0;
// Zero has special meaning.
static constexpr int default_combinatorial_expansion_of_scaffold_hits = -1;

Scaffold_Match_Conditions::Scaffold_Match_Conditions() {
  _enumerate_scaffold_hits_individually = default_enumerate_scaffold_hits_individually;
  _combinatorial_expansion_of_scaffold_hits = default_combinatorial_expansion_of_scaffold_hits;

  return;
}

int
Scaffold_Match_Conditions::IsDefault() const {
  if (! Match_Conditions::IsDefault()) {
    return 0;
  }
  if (_enumerate_scaffold_hits_individually != default_enumerate_scaffold_hits_individually) {
    return 0;
  }
  if (_combinatorial_expansion_of_scaffold_hits != default_combinatorial_expansion_of_scaffold_hits) {
    return 0;
  }

  return 1;
}
