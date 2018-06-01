#include <stdlib.h>

#include "substructure.h"
#include "msi_object.h"
#include "target.h"

void
Substructure_Query::_default_values ()
{
  _each_component_search = 0;

  return;
}

Substructure_Query::Substructure_Query ()
{
  _default_values ();
}

Substructure_Query::Substructure_Query (const char * cc)
{
  _default_values ();

  set_comment(cc);
}

Substructure_Query::Substructure_Query (const const_IWSubstring & cc)
{
  _default_values();

  set_comment(cc);
}

Substructure_Query::Substructure_Query (const IWString & cc)
{
  _default_values();

  set_comment(cc);
}

Substructure_Query::~Substructure_Query ()
{

  return;
}

//#define SHOW_SSQ_OK

int
Substructure_Query::ok () const
{
#ifdef SHOW_SSQ_OK
  cerr << "Checking OK for Substructure Query with " << _number_elements << " queries\n";
#endif

  if (! _operator.ok())
    return 0;

  for (int i = 0; i < _number_elements; i++)
    if (! _things[i]->ok())
      return 0;

  return 1;
}

int
Substructure_Query::debug_print (std::ostream & os) const
{
  os << "Substructure Query ";
  if (1 == _number_elements)
    os << "single query only\n";
  else
  {
    os << _number_elements << " components";
    _operator.debug_print(os);
  }

  IWString indentation = "   ";
  for (int i = 0; i < _number_elements; i++)
  {
    os << " Query " << i << endl;
    _things[i]->debug_print(os, indentation);
  }

  return os.good();
}

int
Substructure_Query::terse_details (std::ostream & os) const
{
  os << "Substructure Query ";
  if (1 == _number_elements)
    os << "single query only\n";
  else
  {
    os << _number_elements << " components\n";
    _operator.debug_print(os);
  }

  IWString indentation = "   ";
  for (int i = 0; i < _number_elements; i++)
  {
    os << " Query " << i << endl;
    _things[i]->terse_details(os, indentation);
  }

  return os.good();
}

int
Substructure_Query::numeric_value (double & d, int ndx) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->numeric_value(d, ndx))
      return 1;
  }

  return 0;
}

void
Substructure_Query::set_numeric_value (double d, int ndx) 
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_numeric_value(d, ndx);
  }

  return;
}

void
Substructure_Query::add_numeric_value (double d)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->add_numeric_value(d);
  }

  return;
}

void
Substructure_Query::discard_all_numeric_values ()
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->discard_all_numeric_values();
  }

  return;
}

int
Substructure_Query::unique_numbers_from_initial_atom_numbers ()
{
  for (int i = 0; i < _number_elements; i++)
  {
    if (! _things[i]->unique_numbers_from_initial_atom_numbers())
    {
      return 0;
    }
  }

  return 1;
}

//#define DEBUG_SUBSTRUCTURE_SEARCH

int
Substructure_Query::_substructure_search (Molecule_to_Match & target,
                                          Substructure_Results & results)
{
  if (_each_component_search)
    return substructure_search_do_each_component(target, results);

  _operator.reset();

#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  cerr << "Query evaluating " << _number_elements << " components\n";
#endif

  for (int i = 0; i < _number_elements; i++)
  {
    if (i > 0 && ! _operator.result_needed(i))    // perhaps as a result of an OR operator
      continue;

    int tmp = _things[i]->substructure_search(target, results);

#ifdef DEBUG_SUBSTRUCTURE_SEARCH
    cerr << "Return code from i = " << i << " is " << tmp << endl;
#endif

    _operator.set_result(i, tmp);

    int result;
    if (! _operator.evaluate(result))     // need to evaluate more components
      continue;

//  We have enough results for a final answer

#ifdef DEBUG_SUBSTRUCTURE_SEARCH
    cerr << "Operator successfully evaluated, result " << result << endl;
#endif

    if (0 == result)       // not a match
      return 0;

    return tmp;     // the whole function rc is nhits for last component to match
  }

#ifdef DEBUG_SUBSTRUCTURE_SEARCH
  cerr << "No match to composite query\n";
#endif

  return 0;
}

int
Substructure_Query::substructure_search (Molecule_to_Match & target,
                                         Substructure_Results & results)
{
  assert (ok());
  assert (target.ok());

  return _substructure_search(target, results);
}

int
Substructure_Query::substructure_search (Molecule & m,
                                         Substructure_Results & results)
{
  Molecule_to_Match target(&m);

  return substructure_search(target, results);
}

int
Substructure_Query::substructure_search (Molecule * m,
                                         Substructure_Results & results)
{
  return substructure_search(*m, results);
}
int
Substructure_Query::substructure_search (Molecule_to_Match & target)
{
  Substructure_Results results;

  return substructure_search(target, results);
}

int
Substructure_Query::substructure_search (Molecule * m)
{
  assert (ok());
  assert (m->ok());

  Molecule_to_Match target(m);
  Substructure_Results results;

  return substructure_search(target, results);
}

int
Substructure_Query::substructure_search_do_each_component (Molecule_to_Match & target,
                                                           Substructure_Results & sresults)
{
  sresults.initialise(target.natoms());

  for (int i = 0; i < _number_elements; i++)
  {
    Substructure_Results tmp;

    if (0 == _things[i]->substructure_search(target, tmp))
      continue;

    sresults.add_embeddings(tmp);
  }

  return sresults.number_embeddings();
}

int
Substructure_Query::set_find_one_embedding_per_atom (int s)
{
  assert (ok());

  if (0 == _number_elements)
  {
    cerr << "Substructure_Query::set_find_one_embedding_per_atom: no components\n";
    return 0;
  }

//cerr << "Substructure_Query::set_find_one_embedding_per_atom:has " << _number_elements << " components\n";

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_find_one_embedding_per_atom(s);
  }

  return 1;
}

int
Substructure_Query::set_max_matches_to_find (int m)
{
  assert (ok());

  if (0 == _number_elements)
  {
    cerr << "Substructure_Query::set_max_matches_to_find: no components\n";
    return 0;
  }

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_max_matches_to_find(m);
  }

  return 1;
}

int
Substructure_Query::set_min_matches_to_find (int m)
{
  assert (ok());

  if (0 == _number_elements)
  {
    cerr << "Substructure_Query::set_min_matches_to_find: no components\n";
    return 0;
  }

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_min_matches_to_find(m);
  }

  return 1;
}

int
Substructure_Query::set_find_unique_embeddings_only (int s)
{
  assert (ok());

  if (0 == _number_elements)
  {
    cerr << "Substructure_Query::set_find_unique_embeddings_only: no queries\n";
    return 0;
  }

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_find_unique_embeddings_only(s);
  }

  return 1;
}

int
Substructure_Query::set_do_not_perceive_symmetry_equivalent_matches (int s)
{
  assert (ok());

  if (0 == _number_elements)
  {
    cerr << "Substructure_Query::set_do_not_perceive_symmetry_equivalent_matches: no queries\n";
    return 0;
  }

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_do_not_perceive_symmetry_equivalent_matches(s);
  }

  return 1;
}

int
Substructure_Query::max_query_atoms_matched_in_search () const
{
  assert (ok());

  iwmax<int> rc(0);
  for (int i = 0; i < _number_elements; i++)
  {
    rc.extra(_things[i]->max_query_atoms_matched_in_search());
  }

  return rc.maxval();
}

int
Substructure_Query::max_atoms_in_query ()
{
  assert (ok());

  iwmax<int> rc(0);
  for (int i = 0; i < _number_elements; i++)
  {
    rc.extra(_things[i]->max_atoms_in_query());
  }

  return rc.maxval();
}

int
Substructure_Query::print_environment_matches (std::ostream & os) const
{
  int environment_present = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->environment_present())
    {
      environment_present = 1;
      break;
    }
  }

  if (! environment_present)
    return os.good();

  os << "Environment matches for query '" << _comment << "'\n";

  for (int i = 0; i < _number_elements; i++)
  {
    os << " Query component " << i << endl;
    _things[i]->print_environment_matches(os);
  }

  return os.good();
}

int
Substructure_Query::add (Single_Substructure_Query * q,
                         int logop)
{
  resizable_array_p<Single_Substructure_Query>::add(q);

  if (_number_elements > 1)
    _operator.add_operator(logop);

  return 1;
}       
Substructure_Atom *
Substructure_Query::query_atom_with_initial_atom_number (int a) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    Substructure_Atom * rc = _things[i]->query_atom_with_initial_atom_number(a);
    if (NULL != rc)
      return rc;
  }

  return NULL;
}
Substructure_Atom *
Substructure_Query::query_atom_with_atom_map_number (int a) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    Substructure_Atom * rc = _things[i]->query_atom_with_atom_map_number(a);
    if (NULL != rc)
      return rc;
  }

  return NULL;
}


int
Substructure_Query::highest_initial_atom_number () const
{
  int rc = -1;

  for (int i = 0; i < _number_elements; i++)
  {
    int h = _things[i]->highest_initial_atom_number();

    if (h > rc)
      rc = h;
  }

  return rc;
}

int
Substructure_Query::highest_atom_map_number () const
{
  int rc = -1;

  for (int i = 0; i < _number_elements; i++)
  {
    int h = _things[i]->highest_atom_map_number();

    if (h > rc)
      rc = h;
  }

  return rc;
}

void
Substructure_Query::set_respect_initial_atom_numbering (int s)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_respect_initial_atom_numbering(s);
  }

  return;
}

void
Substructure_Query::set_ncon (const Min_Max_Specifier<int> & n)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_ncon(n);
  }

  return;
}

void
Substructure_Query::set_distance_between_hits (const Min_Max_Specifier<int> & n)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_distance_between_hits(n);
  }

  return;
}

void
Substructure_Query::set_only_keep_matches_in_largest_fragment (int s)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_only_keep_matches_in_largest_fragment(s);
  }

  return;
}

/*
  We leave out the possibility of query components beyond 0 being
  set differently, why????
*/

int
Substructure_Query::only_keep_matches_in_largest_fragment () const
{
  if (0 == _number_elements)
    return 0;

  return _things[0]->only_keep_matches_in_largest_fragment();
}

void
Substructure_Query::set_embeddings_do_not_overlap (int s)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_embeddings_do_not_overlap(s);
  }

  return;
}

int
Substructure_Query::set_save_matched_atoms(int s)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_save_matched_atoms(s);
  }

  return _number_elements;
}

int
Substructure_Query::set_min_atoms_to_match (int s)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_min_atoms_to_match(s);
  }

  return _number_elements;
}

int
Substructure_Query::set_max_atoms_to_match (int s)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->set_max_atoms_to_match(s);
  }

  return _number_elements;
}

void
Substructure_Query::identify_atom_numbers (extending_resizable_array<int> & a) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->identify_atom_numbers(a);
  }

  return;
}

void
Substructure_Query::identify_atom_map_numbers (extending_resizable_array<int> & a) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->identify_atom_map_numbers(a);
  }

  return;
}

const Substructure_Bond *
Substructure_Query::bond_between_atoms (int a1, int a2) const
{
  const Substructure_Bond * rc;

  for (int i = 0; i < _number_elements; i++)
  {
    rc = _things[i]->bond_between_atoms(a1, a2);

    if (NULL != rc)
      return rc;
  }

  return NULL;
}

const Substructure_Bond *
Substructure_Query::bond_between_atom_map_numbers (int a1, int a2) const
{
  const Substructure_Bond * rc;

  for (int i = 0; i < _number_elements; i++)
  {
    rc = _things[i]->bond_between_atom_map_numbers(a1, a2);

    if (NULL != rc)
      return rc;
  }

  return NULL;
}

void
Substructure_Query::assign_unique_id_from_atom_number_if_set (extending_resizable_array<int> & numbers_in_use)
{
  for (auto i = 0; i < _number_elements; ++i)
  {
    _things[i]->assign_unique_id_from_atom_number_if_set(numbers_in_use);
  }

  return;
}
int
Substructure_Query::assign_atom_map_numbers(int & amap)
{
  int rc = 0;
  for (auto i = 0; i < _number_elements; ++i)
  {
    rc += _things[i]->assign_atom_map_numbers(amap);
  }

  return rc;
}

int
Substructure_Query::print_connectivity_graph(std::ostream & output) const
{
  output << "Substructure_Query::print_connectivity_graph:has " << _number_elements << " components\n";

  for (int i = 0; i < _number_elements; ++i)
  {
    _things[i]->print_connectivity_graph(output);
  }

  return 1;
}
