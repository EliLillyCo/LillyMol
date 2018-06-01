#include <stdlib.h>
#include <limits>

#include "reaction_match_conditions.h"

// Don't forget to update the assignment operator. 

Match_Conditions::Match_Conditions ()
{
  _verbose = 0;

  _find_unique_embeddings = 1;

  _ignore_not_reacting = 0;

  _process_hit_number = -1;

  _one_embedding_per_start_atom = 0;

  _ignore_symmetry_related_matches = 0;

  _suppress_if_more_than_this_many_substructure_search_hits = std::numeric_limits<int>::max();

  return;
}

Match_Conditions &
Match_Conditions::operator= (const Match_Conditions & rhs)
{
  _verbose = rhs._verbose;

  _find_unique_embeddings = rhs._find_unique_embeddings;

  _ignore_not_reacting = rhs._ignore_not_reacting;

  _process_hit_number = rhs._process_hit_number;

  _one_embedding_per_start_atom = rhs._one_embedding_per_start_atom;

  _ignore_symmetry_related_matches = rhs._ignore_symmetry_related_matches;

  _multiple_match_string = rhs._multiple_match_string;

  return *this;
}

Sidechain_Match_Conditions::Sidechain_Match_Conditions ()
{
  _make_new_reagent_for_each_hit = 0;

  _max_matches_to_find = 0;

  _strip_reagents_to_largest_fragment = 0;

  return;
}

Sidechain_Match_Conditions &
Sidechain_Match_Conditions::operator= (const Sidechain_Match_Conditions & rhs)
{
  Match_Conditions::operator= (rhs);

  _make_new_reagent_for_each_hit = rhs._make_new_reagent_for_each_hit;

  _max_matches_to_find = rhs._max_matches_to_find;

  _strip_reagents_to_largest_fragment = rhs._strip_reagents_to_largest_fragment;

  return *this;
}

Scaffold_Match_Conditions::Scaffold_Match_Conditions ()
{
  _enumerate_scaffold_hits_individually = 0;
  _combinatorial_expansion_of_scaffold_hits = -1;    // zero has special meaning

  return;
}
