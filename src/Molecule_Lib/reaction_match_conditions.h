#ifndef MATCH_CONDITIONS_H
#define MATCH_CONDITIONS_H

#include "iwstring.h"

/*
  We need a means of passing around what to do with the different
  possible match results
*/

class Match_Conditions
{
  protected:
    int _verbose;

    int _ignore_not_reacting;

    int _process_hit_number;

//  When we use just one of several possible sidechain matches, we can have
//  some text appended to the name to indicate that

    IWString _multiple_match_string;

    int _find_unique_embeddings;

    int _one_embedding_per_start_atom;
    
    int _ignore_symmetry_related_matches;

    int _ignore_multiple_matches_involving_atoms_not_changing;

//  Mar 2015. Suppress any reagent that returns multiple substructure hits

    int _suppress_if_more_than_this_many_substructure_search_hits;

  public:
    Match_Conditions ();

    Match_Conditions & operator = (const Match_Conditions &);

    int verbose () const { return _verbose;}
    void set_verbose (int v) { _verbose = v;}

    int ignore_not_reacting () const { return _ignore_not_reacting;}
    void set_ignore_not_reacting (int i) { _ignore_not_reacting = i;}

    int find_unique_embeddings_only () const { return _find_unique_embeddings; }
    void set_find_unique_embeddings_only (int i) { _find_unique_embeddings = i;}

    int process_hit_number () const { return _process_hit_number;}
    void set_process_hit_number (int p) { _process_hit_number = p;}

    int one_embedding_per_start_atom () const { return _one_embedding_per_start_atom;}
    void set_one_embedding_per_start_atom (int p) { _one_embedding_per_start_atom = p;}

    int ignore_symmetry_related_matches () const { return _ignore_symmetry_related_matches;}
    void set_ignore_symmetry_related_matches (int p) { _ignore_symmetry_related_matches = p;}

    const IWString & multiple_match_string () const { return _multiple_match_string;}
    void set_multiple_match_string (const IWString & s) { _multiple_match_string = s;}

    void set_ignore_multiple_substucture_matches (int s) { _suppress_if_more_than_this_many_substructure_search_hits = s;}
    int suppress_if_more_than_this_many_substructure_search_hits () const { return _suppress_if_more_than_this_many_substructure_search_hits;}
};

class Scaffold_Match_Conditions : public Match_Conditions
{
  private:

// If there are multiple matches in the scaffold, we may want to make a separate product from each

    int _enumerate_scaffold_hits_individually;
    int _combinatorial_expansion_of_scaffold_hits;

  public:
    Scaffold_Match_Conditions ();

    int enumerate_scaffold_hits_individually () const { return _enumerate_scaffold_hits_individually;}
    void set_enumerate_scaffold_hits_individually (int e) { _enumerate_scaffold_hits_individually = e;}

    void set_combinatorial_expansion_of_scaffold_hits(const int s) { _combinatorial_expansion_of_scaffold_hits = s;}
    int combinatorial_expansion_of_scaffold_hits() const { return _combinatorial_expansion_of_scaffold_hits;}
};

class Sidechain_Match_Conditions : public Match_Conditions
{
  private:

//  In the sidechain, if there are multiple matches, we can make regio-isomers

    int _make_new_reagent_for_each_hit;

    int _max_matches_to_find;

    int _strip_reagents_to_largest_fragment;

  public:
    Sidechain_Match_Conditions ();

    Sidechain_Match_Conditions & operator= (const Sidechain_Match_Conditions &);

    int make_new_reagent_for_each_hit () const { return _make_new_reagent_for_each_hit;}
    void set_make_new_reagent_for_each_hit (int m) { _make_new_reagent_for_each_hit = m;}

    int max_matches_to_find () const { return _max_matches_to_find;}
    void set_max_matches_to_find (int p) { _max_matches_to_find = p;}

    int strip_reagents_to_largest_fragment () const { return _strip_reagents_to_largest_fragment;}
    void set_strip_reagents_to_largest_fragment (int p) { _strip_reagents_to_largest_fragment = p;}
};


#endif
