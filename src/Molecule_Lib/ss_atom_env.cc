#include <stdlib.h>
#include <memory>

/*
  Implementation of the Substructure_Atom_Environment class
*/

#include "substructure.h"
#include "tokenise_atomic_smarts.h"
#include "target.h"
#include "misc.h"

static int environment_only_matches_unmatched_atoms = 0;

void
set_atom_environment_only_matches_unmatched_atoms (int s)
{
  environment_only_matches_unmatched_atoms = s;

  return;
}

Substructure_Atom_Environment::Substructure_Atom_Environment ()
{
  _all_components_single_atoms = 1;
}

int
Substructure_Atom_Environment::ok () const
{
  if (1 == _number_elements && 0 == _operator.number_operators())
    ;
  else if (_number_elements == _operator.number_operators())    // true during building
    ;
  else if (_number_elements != _operator.number_results())
    return 0;

  return _operator.ok();
}

int
Substructure_Atom_Environment::debug_print (std::ostream & os) const
{
  os << "Details on Substructure_Atom_Environment object with " << _number_elements << " components\n";

  _operator.debug_print(os);

  if (! ok())
    os << "Warning, NOT OK\n";

  return os.good();
}

int
Substructure_Atom_Environment::involves_aromatic_bond_specifications (int & _need_to_compute_ring_membership) const
{ 
  for (int i = 0; i < _number_elements; i++)
  {
    if (_things[i]->involves_aromatic_bond_specifications(_need_to_compute_ring_membership))
      return 1;
  }

  return 0;
}

//#define DEBUG_SS_ATOM_ENV_CREATE_FROM_SMARTS

int
Substructure_Atom_Environment::create_from_smarts (const Atomic_Smarts_Component & env)
{
  assert (env.starts_with("$(") && env.ends_with(')'));

  const_IWSubstring mysmarts = env;      // so we can make changes

  mysmarts.remove_leading_chars(2);
  mysmarts.chop();

#ifdef DEBUG_SS_ATOM_ENV_CREATE_FROM_SMARTS
  cerr << "Substructure_Atom_Environment::create_from_smarts: parsing " << env << endl;
#endif

  Substructure_Atom * a = new Substructure_Atom;

  if (! a->parse_smarts_specifier(mysmarts))
  {
    cerr << "Substructure_Atom_Environment::create_from_smarts: cannot parse '" << env << "'\n";
    delete a;
    return 0;
  }

  add (a);

  if (a->number_children())
    _all_components_single_atoms = 0;

// Unfortunately complex initialising the operator. Remember that an 
// environment component comes here with the operator following it

  _operator.set_unary_operator(_operator.number_results() - 1, env.unary_operator());

  if (IW_LOGEXP_UNDEFINED != env.op())
    _operator.add_operator(env.op());

#ifdef DEBUG_SS_ATOM_ENV_CREATE_FROM_SMARTS
  cerr << "After parsing smarts\n";
  debug_print(cerr);
#endif

  assert(ok());

  return 1;
}

//#define DEBUG_SS_ATOM_ENV_MATCHES

int
Substructure_Atom_Environment::matches (Target_Atom & target_atom,
                                        int atoms_in_target,
                                        const int * already_matched)
{
#ifdef DEBUG_SS_ATOM_ENV_MATCHES
  cerr << "Substructure_Atom_Environment::matches: " << _number_elements << " components, trying to match atom " << target_atom.atom_number() << endl;
#endif

  _operator.reset();

  int * tmp;
/*
  Oct 99, this turns out to create all kinds of problems.  Just make a
  vector for all searches

  if (_all_components_single_atoms)
    tmp = NULL;
  else */

  tmp = new int[atoms_in_target]; std::unique_ptr<int[]> free_tmp(tmp);

  int rc = _matches(target_atom, atoms_in_target, already_matched, tmp);

#ifdef DEBUG_SS_ATOM_ENV_MATCHES
  cerr << "Substructure_Atom_Environment::matches: returning " << rc << endl;
#endif

  return rc;
}

int
Substructure_Atom_Environment::_matches (Target_Atom & target_atom,
                                         int atoms_in_target,
                                         const int * already_matched_by_query,
                                         int * already_matched)
{
  assert (0 == already_matched_by_query[target_atom.atom_number()]);

  for (int i = 0; i < _number_elements; i++)
  {
    if (! _operator.result_needed(i))
      continue;

#ifdef DEBUG_SS_ATOM_ENV_MATCHES
    cerr << "Substructure_Atom_Environment::_matches: trying to match component " << i << endl;
#endif

    set_vector(already_matched, atoms_in_target, 0);

    if ( _match_component(target_atom, i, already_matched_by_query, already_matched))
    {
#ifdef DEBUG_SS_ATOM_ENV_MATCHES
      cerr << "Component " << i << " matches\n";
#endif
      _operator.set_result(i, 1);
    }
    else
    {
#ifdef DEBUG_SS_ATOM_ENV_MATCHES
      cerr << "Component " << i << " does not match\n";
#endif
      _operator.set_result(i, 0);
    }

#ifdef DEBUG_SS_ATOM_ENV_MATCHES
    cerr << "Result for component " << i << " is " << _operator.result(i) << endl;
    _operator.debug_print(cerr);
#endif

    _things[i]->recursive_release_hold();

    int rc;
    if (_operator.evaluate(rc))
      return rc;
  }

  return 0;
}

int
Substructure_Atom_Environment::_match_component (Target_Atom & target_atom,
                                                 int which_component,
                                                 const int * already_matched_by_query,
                                                 int * already_matched)
{
#ifdef DEBUG_SS_ATOM_ENV_MATCHES
  cerr << "Substructure_Atom_Environment::_match_component: component " << which_component << " trying to match atom " << target_atom.atom_number() << endl;
#endif

  Substructure_Atom * root_atom = _things[which_component];

  assert (0 == already_matched[target_atom.atom_number()]);

  if (! root_atom->matches(target_atom, already_matched_by_query))
    return 0;

  if (0 == root_atom->number_children())
    return 1;

  if (environment_only_matches_unmatched_atoms)
  {
    int na = target_atom.atoms_in_target_molecule();

    for (int j = 0; j < na; j++)
    {
      already_matched[j] = already_matched_by_query[j];
#ifdef DEBUG_SS_ATOM_ENV_MATCHES
      cerr << "Atom " << j << " already_matched " << already_matched[j] << endl;
#endif
    }
  }

  root_atom->set_hold(&target_atom, already_matched);

  Query_Atoms_Matched matched_atoms;
  matched_atoms.resize(23);         // just an arbitrary sizing

  root_atom->add_your_children(matched_atoms);

  int atom_to_process = 0;
  while (atom_to_process >= 0)
  {
    Substructure_Atom * a = (Substructure_Atom *) matched_atoms[atom_to_process];

    if (! a->move_to_next_match_from_current_anchor(already_matched, matched_atoms))
    {
#ifdef DEBUG_SS_ATOM_ENV_MATCHES
      cerr << "Move to next failed for atom " << a->unique_id() << endl;
#endif

      a->remove_your_children(matched_atoms, already_matched);
      if (a->or_id() && atom_to_process < matched_atoms.number_elements() - 1 &&
          a->or_id() == matched_atoms[atom_to_process + 1]->or_id())
        matched_atoms.remove_item(atom_to_process);
      else
        atom_to_process--;
    }
    else
    {
#ifdef DEBUG_SS_ATOM_ENV_MATCHES
      cerr << "Move to next match succeeded " << 
              a->unique_id() << "(" << a->atom_number_matched() <<
              "), or = " << a->or_id() <<
              " atom to process = " << atom_to_process << " matched = " << matched_atoms.number_elements() << endl;
#endif
      int orid = a->or_id();
      if (orid)
        remove_atoms_with_same_or(matched_atoms, atom_to_process + 1, orid);

      a->add_your_children(matched_atoms);   // does nothing if already added

      if (atom_to_process == matched_atoms.number_elements() - 1)    // all atoms matched
        return 1;

      atom_to_process++;
    }
  }

  return 0;
}
