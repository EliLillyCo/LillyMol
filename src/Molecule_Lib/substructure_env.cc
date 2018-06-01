/*
  All the functions for matching a query environment
*/

#include <stdlib.h>
#include <limits>
#include <memory>

#include "misc.h"
#include "substructure.h"
#include "target.h"

Substructure_Environment::Substructure_Environment ()
{
  _unique_id = 0;

  _or_id = 0;
  _and_id = 0;

  _no_other_substituents_allowed = 0;

  _max_environment_matches_per_attachment_point = std::numeric_limits<int>::max();

  _environments_can_share_attachment_points = 1;

  _max_matches_to_find = std::numeric_limits<int>::max();

  _hydrogen_ok_as_environment_match = 0;
 
  return;
}

int
Substructure_Environment::ok () const
{
  if (! resizable_array_p<Substructure_Atom>::ok())
    return 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (! _things[i]->ok())
      return 0;
  }

  return 1;
}

/*
  debug_print and terse_details need to share some code
*/

int
Substructure_Environment::_print_common_info (std::ostream & os,
                                const IWString & indentation) const
{
  int np = _possible_parents.number_elements();

  os << indentation << "Environment atom with " << np << " possible parent";
  if (np > 1)
    os << 's';
  os << ':';

  for (int i = 0; i < np; i++)
  {
    const Substructure_Atom * p = _possible_parents[i];
    os << ' ' << p->unique_id();
  }

  os << " and " << _number_elements << " components\n";

  if (_query_environment_match_as_rejection)
    os << indentation << "Match as rejection is " << _query_environment_match_as_rejection << endl;

  if (_and_id)
    os << indentation << " and_id " << _and_id << endl;

  if (_or_id)
    os << indentation << " or_id " << _or_id << endl;

  if (_no_other_substituents_allowed)
    os << indentation << " no other substituents allowed\n";

  return os.good();
}

int
Substructure_Environment::terse_details (std::ostream & os,
                                         const IWString & indentation) const
{
  _print_common_info(os, indentation);

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->terse_details(os, indentation + "  ");
  }

  return os.good();
}

int
Substructure_Environment::debug_print (std::ostream & os,
                                       const IWString & indentation) const
{
  _print_common_info(os, indentation);

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->recursive_debug_print(os, indentation + "  ");
  }

  return os.good();
}

void
Substructure_Environment::assign_unique_atom_numbers (int & id)
{
  _unique_id = id++;
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->assign_unique_atom_numbers(id);
  }

  return;
}

int
Substructure_Environment::attributes_specified ()
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
    rc += _things[i]->attributes_specified();

  return rc;
}

int
Substructure_Environment::involves_aromatic_bond_specifications (int & r) const
{
  for (int i = 0; i < _number_elements; i++)
    if (_things[i]->involves_aromatic_bond_specifications(r))
      return 1;

  return 0;
}

Substructure_Environment_Match::Substructure_Environment_Match ()
{
  _query_environment_match_as_rejection = 0;
}

Substructure_Environment_Rejection::Substructure_Environment_Rejection()
{
  _query_environment_match_as_rejection = 1;
}

//#define DEBUG_ENVIRONMENT_SEARCH

/*
  Initially I set up the environment search to only allow one match per
  attachment point, but in Apr 98, I needed to change that. Hopefully
  this still works...
*/

int
Substructure_Atom::environment_search (Query_Atoms_Matched & matched_query_atoms,
                                       int * previously_matched_atoms)
{
  assert (0 == matched_query_atoms.number_elements());

  int nc = _children.number_elements();

#ifdef DEBUG_ENVIRONMENT_SEARCH
  cerr << "Atom " << _unique_id << " beginning environment_search, " << nc << " children\n";
#endif

  int rc = 0;       // the number of hits

  while (move_to_next_match_from_current_anchor(previously_matched_atoms, matched_query_atoms))
  {
#ifdef DEBUG_ENVIRONMENT_SEARCH
  cerr << "Env base atom matches atom " << _current_hold_atom->atom_number() << endl;
#endif

    if (0 == nc)    // no children, have a match right here
    {
      rc++;
      continue;
    }

    matched_query_atoms.resize_keep_storage(0);
    add_your_children(matched_query_atoms);
    int atom_to_process = 0;

    while (atom_to_process >= 0)
    {
      Substructure_Atom * a = (Substructure_Atom *) matched_query_atoms[atom_to_process];

#ifdef DEBUG_ENVIRONMENT_SEARCH
      cerr << "Atom to process is " << atom_to_process << " which is " << a->unique_id() << endl;
#endif

      if (! a->move_to_next_match_from_current_anchor(previously_matched_atoms, matched_query_atoms))
      {
#ifdef DEBUG_ENVIRONMENT_SEARCH
        cerr << "Move to next failed for query environment atom " << a->unique_id() << endl;
#endif

        a->remove_your_children(matched_query_atoms, previously_matched_atoms);

        if (a->or_id() && atom_to_process < matched_query_atoms.number_elements() - 1 &&
            a->or_id() == matched_query_atoms[atom_to_process + 1]->or_id())
          matched_query_atoms.remove_item(atom_to_process);
        else
          atom_to_process--;
      }
      else
      {
#ifdef DEBUG_ENVIRONMENT_SEARCH
        cerr << "Move to next match succeeded " << 
                a->unique_id() << "(" << a->atom_number_matched() << "), or = " << a->or_id() <<
                " atom to process = " << atom_to_process << " matched = " << matched_query_atoms.number_elements() << endl;
#endif
        int orid = a->or_id();
        if (orid)
          remove_atoms_with_same_or(matched_query_atoms, atom_to_process + 1, orid);
  
        a->add_your_children(matched_query_atoms);   // does nothing if already added

        atom_to_process++;

        if (atom_to_process >= matched_query_atoms.number_elements())      // the == condition is the only one which should ever happen
        {
          rc++;
          break;        // break from while (atom_to_process >= 0) loop, we are only interested in one embedding per start atom
        }
      }
    }
  }

#ifdef DEBUG_ENVIRONMENT_SEARCH
  cerr << "environment search returning " << rc << endl;
#endif

  return rc;
}

//#define DEBUG_SS_ENV_MATCHES

/*
  Matches a query environment.
  This is complicated by the presence of the presence of the _no_other_substituents_allowed
  attribute. In that case, we must check all possible attachment points, in
  order to make sure that there are no non-matching groups substituted
  at those points
*/

int
Substructure_Environment::matches (int * previously_matched_atoms,
                                   int * anchor_atom_available)
{
  Query_Atoms_Matched qam;
  qam.resize(20);    // 20 seems pretty large

  int nhits = 0;
  int np = _possible_parents.number_elements();

  int parents_with_unmatched_connections = 0;

#ifdef DEBUG_SS_ENV_MATCHES
  cerr << "Environment object " << _unique_id << " with " << np << " attach points and " << _number_elements << " components matching\n";
#endif

  for (int i = 0; i < np; i++)    // loop over possible parents
  {
    Substructure_Atom * p = _possible_parents[i];
    Target_Atom * a = p->current_hold_atom();
    if (NULL == a)
      continue;

#ifdef DEBUG_SS_ENV_MATCHES
    cerr << "Trying possible parent " << i << " atom " << a->atom_number() << " uc = " << p->unmatched_connections(previously_matched_atoms) << " anchor avail " << anchor_atom_available[a->atom_number()] << endl;
    const Molecule * m = a->m();
    for (int x = 0; x < m->natoms(); ++x)
    {
      cerr << " previously_matched_atoms[" << x << "] = " << previously_matched_atoms[x] << endl;
    }
#endif

    if (! anchor_atom_available[a->atom_number()])
      continue;

    const int parent_unmatched_connections = p->unmatched_connections(previously_matched_atoms);

    if (0 == parent_unmatched_connections)
    {
#ifdef DEBUG_SS_ENV_MATCHES
      if (_hydrogen_ok_as_environment_match && a->hcount())
        cerr << "HYDROGEN matches as environment match\n";
#endif
      if (_hydrogen_ok_as_environment_match && a->hcount())
        nhits++;
      continue;
    }

    parents_with_unmatched_connections++;

    int found_match_this_parent = 0;

    for (int j = 0; j < _number_elements; j++)   // now loop over all members of the environment
    {
      Substructure_Atom * aj = _things[j];

#ifdef DEBUG_SS_ENV_MATCHES
      cerr << "Preparing component " << j << endl;
#endif

      aj->prepare_for_matching(a);
      aj->set_parent(p, &_bond);

      qam.resize_keep_storage(0);

      int esearch = aj->environment_search(qam, previously_matched_atoms);

      aj->recursive_release_hold();
      aj->no_parent();

#ifdef DEBUG_SS_ENV_MATCHES
      cerr << "Environment component " << j << " matches " << esearch << " times\n";
#endif

      if (esearch)
      {
        if (esearch > _max_environment_matches_per_attachment_point)
          esearch = _max_environment_matches_per_attachment_point;

        _matches[j] += esearch;

        nhits += esearch;
        anchor_atom_available[a->atom_number()] = _environments_can_share_attachment_points;

        if (! _hits_needed.is_set() && 0 == _no_other_substituents_allowed)
          return 1;

        if (nhits >= _max_matches_to_find)
          break;

//      Found a match at this attachment point. Break to move to the next attachment point
//      unless we might be making sure there is nothing else at this location

        found_match_this_parent += esearch;
        if (0 == _no_other_substituents_allowed)    // no need to keep looking
          break;
      }
    }

#ifdef DEBUG_SS_ENV_MATCHES
    cerr << "Found " << found_match_this_parent << " matches from parent\n";
#endif

    if (_no_other_substituents_allowed)
    {
      if (found_match_this_parent != parent_unmatched_connections)   // means something else at this parent
        return 0;
    }
  }

#ifdef DEBUG_SS_ENV_MATCHES
  cerr << "After checking components, nhits = " << nhits << ", pwuc = " << parents_with_unmatched_connections << endl;
  cerr << "_hits_needed.matches? " << _hits_needed.matches(nhits) << endl;
#endif

  if (nhits > 0)
    return _hits_needed.matches(nhits);

  if (_hits_needed.is_set() && _hits_needed.matches(nhits))
    return 1;

// If we got no hits, but all the possible parents were fully matched,
// that counts as an OK match

  if (0 == parents_with_unmatched_connections && _no_other_substituents_allowed)
    return 1;

  return 0;
}

/*
  Try to match a group (consisting of possibly only one) environment_atoms.
  If any of them match, return 1
*/

int
Single_Substructure_Query::_query_environment_or_group_matched (int * previously_matched_atoms,
                         const resizable_array<Substructure_Environment *> & or_group,
                         int * anchor_atom_available)
{
  const int no = or_group.number_elements();

//cerr << "Single_Substructure_Query::_query_environment_or_group_matched:checking or group with " << no << " components\n";

  for (int i = 0; i < no; i++)
  {
    Substructure_Environment * e = or_group[i];
    if (e->matches(previously_matched_atoms, anchor_atom_available))
      return 1;
  }

  return 0;    // none of them matched.
}

//#define DEBUG_AND_GROUP

/*
  Try to match an AND group (consisting of possibly only one) environment_atoms.
  If all match, return 1
*/

int
Single_Substructure_Query::_query_environment_and_group_matched (int * previously_matched_atoms,
                         const resizable_array<Substructure_Environment *> & zgroup,
                         int * anchor_atom_available)
{
  int no = zgroup.number_elements();

#ifdef DEBUG_AND_GROUP
  cerr << "Checking and group of " << no << " environments\n";
#endif

  for (int i = 0; i < no; i++)
  {
    Substructure_Environment * e = zgroup[i];
    if (! e->matches(previously_matched_atoms, anchor_atom_available))
    {
#ifdef DEBUG_AND_GROUP
      cerr << "Member " << i << " of and group size " << no << " did not match\n";
#endif

      return 0;
    }
#ifdef DEBUG_AND_GROUP
    cerr << "Member " << i << " of " << no << " and group found match\n";
#endif
  }

#ifdef DEBUG_AND_GROUP
  cerr << "And grouping of " << no << " members all matched\n";
#endif

  return 1;    // they all matched.
}

/*
  Try to match a group (consisting of possibly only one) environment_atoms.
  If any of them match, return 1

  Not used right now. Hard to know what to do with this. A more flexible
  approach would be to impose a requirement on the number of matches,
  rather than what we have now.
*/

int
Single_Substructure_Query::_query_environment_xor_group_matched (int * previously_matched_atoms,
                         const resizable_array<Substructure_Environment *> & zgroup,
                         int * anchor_atom_available)
{
  int nhits = 0;

  int no = zgroup.number_elements();
  for (int i = 0; i < no; i++)
  {
    Substructure_Environment * e = zgroup[i];
    if (e->matches(previously_matched_atoms, anchor_atom_available))
    {
      nhits++;
      if (nhits > 1)
        return 0;
    }
  }

  return nhits;
}

/*
  The environment rejections are processed AND group at a time.
  Every AND group must match.
*/

//#define DEBUG_ENVIRONMENT_REJECTIONS

int
Single_Substructure_Query::_environment_rejections_matched (const int atoms_in_target_molecule, 
                                       int * previously_matched_atoms,
                                       int * env_already_done)
{
  int nr = _environment_rejections.number_elements();

#ifdef DEBUG_ENVIRONMENT_REJECTIONS
  cerr << "Checking " << nr << " rejections\n";
#endif

  int * anchor_atom_available = new_int(atoms_in_target_molecule, 1); std::unique_ptr<int[]> free_anchor_atom_available(anchor_atom_available);

  for (int i = 0; i < nr; i++)
  {
    if (env_already_done[i])
      continue;

    resizable_array<Substructure_Environment *> z_group;

    Substructure_Environment * e = _environment_rejections[i];
    z_group.add(e);

    if (e->and_id())
    {
      for (int j = i + 1; j < nr; j++)
      {
        Substructure_Environment * ej = _environment_rejections[j];
        if (ej->and_id() == e->and_id())
        {
          z_group.add(ej);
          env_already_done[j] = 1;
        }
      }
    }

#ifdef DEBUG_ENVIRONMENT_REJECTIONS
    cerr << "Grouping contains " << z_group.number_elements() << " components\n";
#endif

    if (_query_environment_and_group_matched(previously_matched_atoms, z_group, anchor_atom_available))
      return 1;
  }

  return 0;   // none of the group matched
}

//#define DEBUG_QUERY_ENVIRONMENT_ALSO_MATCHED

/*
  By default, every component of the environment must match.
  The environment is processed OR group at a time.
  Every OR group must match.
*/

int
Single_Substructure_Query::_query_environment_also_matched (const int atoms_in_target_molecule,
                                           int * previously_matched_atoms,
                                           int * env_already_done)
{
  int ne = _environment.number_elements();

#ifdef DEBUG_QUERY_ENVIRONMENT_ALSO_MATCHED
  cerr << "Checking " << ne << " enviromment groups\n";
#endif

  int * anchor_atom_available = new_int(atoms_in_target_molecule, 1); std::unique_ptr<int[]> free_anchor_atom_available(anchor_atom_available);

  for (int i = 0; i < ne; i++)
  {
    if (env_already_done[i])
      continue;

    resizable_array<Substructure_Environment *> z_group;

    Substructure_Environment * e = _environment[i];
    z_group.add(e);

    if (e->or_id())
    {
      for (int j = i + 1; j < ne; j++)
      {
        if (env_already_done[j])
          continue;

        Substructure_Environment * ej = _environment[j];
        if (ej->or_id() == e->or_id())
        {
          z_group.add(ej);
          env_already_done[j] = 1;
        }
      }
    }

    if (! _query_environment_or_group_matched(previously_matched_atoms, z_group, anchor_atom_available))
      return 0;
  }

  return 1;
}

/*
  There are lots of interesting aspects to evaluating the query environment.
*/

int
Single_Substructure_Query::_query_environment_also_matched(Query_Atoms_Matched & matched_query_atoms, 
                                      int atoms_in_target_molecule)
{
  int ne = _environment.number_elements();

  int nr = _environment_rejections.number_elements();

  if (0 == ne && 0 == nr)
    return 1;

#ifdef DEBUG_QUERY_ENVIRONMENT_ALSO_MATCHED
  cerr << "Single_Substructure_Query::_query_environment_also_matched:checking " << ne << " match and " << nr << " environment rejections\n";
#endif

  int * previously_matched = new_int(atoms_in_target_molecule); std::unique_ptr<int[]> free_previously_matched(previously_matched);

//cerr << "Single_Substructure_Query::_query_environment_also_matched: _environment_must_match_unmatched_atoms " << _environment_must_match_unmatched_atoms << endl;

  if (_environment_must_match_unmatched_atoms)
  {
    int na = matched_query_atoms.number_elements();
    for (int i = 0; i < na; i++)
    {
      const Substructure_Atom * a = matched_query_atoms[i];
      const int j = a->atom_number_matched();

#ifdef DEBUG_QUERY_ENVIRONMENT_ALSO_MATCHED
      cerr << "Query atom " << a->unique_id() << " matched with atom " << j;
      if (0 != previously_matched[j])
        cerr << " YIPES, that atom not marked as matched";
      cerr << endl;
#endif

      assert (j >= 0 && 0 == previously_matched[j]);
      previously_matched[j] = 1;
    }
  }

// Really just needs to be sized the max of (ne and nr)

  int * env_already_done = new_int(ne + nr); std::unique_ptr<int[]> free_env_already_done(env_already_done);

#ifdef DEBUG_QUERY_ENVIRONMENT_ALSO_MATCHED
  cerr << " ne " << ne << " and nr " << nr << endl;
#endif

  if (ne)
  {
    if (! _query_environment_also_matched(atoms_in_target_molecule, previously_matched, env_already_done))
    {
      _no_match_to_environment++;
      return 0;
    }
  }

  if (nr)   // either no environment matches, or we matched the environment
  {
    if (ne)                                 // if already used, reset to 0
      set_vector(env_already_done, nr, 0);

    if (_environment_rejections_matched(atoms_in_target_molecule, previously_matched, env_already_done))
    {
      _match_to_environemt_rejection++;
      return 0;
    }
  }

  return 1;
}

int
Substructure_Environment::print_environment_matches (std::ostream & os) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    os << "    component " << i << " matched ";
    if (i < _matches.number_elements())
      os << _matches.item(i);
    else
      os << '0';
      
    os << " times\n";
  }

  return os.good();
}
