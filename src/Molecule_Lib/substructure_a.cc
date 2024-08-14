#include <ctype.h>
#include <assert.h>
#include <iostream>
#include <memory>

#define IWARAY_EACH_IMPLEMENTATION

#include "Foundational/iwmisc/misc.h"
#include "substructure.h"
#include "target.h"
#include "misc2.h"
#include "smiles.h"

using std::cerr;

void
Substructure_Atom::_default_values()
{
  _unique_id = -1;
  _initial_atom_number = -1;
  _current_hold_atom = nullptr;

  _atom_map_number = -1;

  _global_match_id = -1;

  _or_id = 0;

  _ring_id = 0;
  _fused_system_id = 0;

  _include_in_embedding = 1;

  _bond_to_parent     = nullptr;
  _parent             = nullptr;

  _con = -1;
  _anchor = nullptr;

  _sum_all_preference_hits = 0;

  _match_as_match_or_rejection = 1;

  _fragment_id = 0;

  _atom_type = 0;

  _atom_type_group = 0;

  _on_queue_for_matching = 0;

  return;
}

Substructure_Atom::Substructure_Atom()
{
  _default_values();

  return;
}

Substructure_Atom::Substructure_Atom(const Element * e) :
                      Substructure_Atom_Specifier(e)
{
  _default_values();

  return;
}

Substructure_Atom::~Substructure_Atom()
{
//cerr << "Deleting Substructure atom " << _unique_id << " at " << this << '\n';

  if (-5 == _unique_id)
    cerr << "Deleting already deleted Substructure_Atom\n";

  assert(ok());

  _unique_id = -5;

  return;
}

/*
  The first time we assign the unique_id, retain that as _INITIAL_ATOM_NUMBER -
  unless it has already been set previously (by the smarts parser for example)
*/

void
Substructure_Atom::set_unique_id(int i)
{
  assert (i >= 0);

  if (_unique_id < 0 && _initial_atom_number < 0)
    _initial_atom_number = i;

  _unique_id = i;

  return;
}

int 
Substructure_Atom::unique_id_from_initial_atom_number()
{
  if (_initial_atom_number < 0)
  {
    cerr << "Substructure_Atom::unique_id_from_initial_atom_number: initial atom number not set\n";
    return 0;
  }

  _unique_id = _initial_atom_number;

  for (Substructure_Atom* c : _children)
  {
    c->unique_id_from_initial_atom_number();
  }

  return 1;
}

int
Substructure_Atom::collect_all_atoms(extending_resizable_array<Substructure_Atom *> & completed)
{
  completed[_unique_id] = this;

  int rc = 1;
  for (Substructure_Atom* c : _children)
  {
    c->collect_all_atoms(completed);
  }

  return rc;
}

int
Substructure_Atom::highest_initial_atom_number() const
{
  int rc = _initial_atom_number;

  for (Substructure_Atom* c : _children)
  {
    int tmp = c->highest_initial_atom_number();
    if (tmp > rc)
      rc = tmp;
  }

  return rc;
}

int
Substructure_Atom::highest_atom_map_number() const
{
  int rc = _atom_map_number;

  for (const Substructure_Atom* c: _children)
  {
    const int tmp = c->highest_atom_map_number();
    if (tmp > rc)
      rc = tmp;
  }

  return rc;
}

int
Substructure_Atom::number_descendants() const
{
  int rc = 1;

  for (const Substructure_Atom* c : _children)
  {
    rc += c->number_descendants();
  }

  return rc;
}

int
Substructure_Atom::preferences_present() const
{
  if (_preference_value)
    return 1;

  if (_preferences.number_elements())
    return 1;

  for (Substructure_Atom_Specifier* c : _components)
  {
    if (c->preference_value())
      return 1;
  }

  for (const Substructure_Atom* c : _children)
  {
    if (c->preferences_present())
      return 1;
  }

  return 0;
}

int
Substructure_Atom::set_preference_value(int s)
{
  if (_components.empty())
  {
    cerr << "Substructure_Atom::set_preference_value:no components\n";
    return 0;
  }

  _components[0]->set_preference_value(s);

  return 1;
}

int
Substructure_Atom::ring_ids_present() const
{
  if (_ring_id > 0)
    return 1;

  for (const Substructure_Atom* c : _children)
  {
    if (c->ring_ids_present())
      return 1;
  }

  return 0;
}

int
Substructure_Atom::fused_system_ids_present() const
{
  if (_fused_system_id > 0)
    return 1;

  for (const Substructure_Atom* c : _children)
  {
    if (c->fused_system_ids_present())
      return 1;
  }

  return 0;
}

int
Substructure_Atom::fragment_ids_present() const
{
  if (_fragment_id > 0)
    return 1;

  for (const Substructure_Atom* c : _children)
  {
    if (c->fragment_ids_present())
      return 1;
  }

  return 0;
}

/*
  during detection of implicit rings, we need to know whether or not one
  Substructure_Atom is in the ring closure bond list of another atom
*/

int
Substructure_Atom::has_ring_closure_bond_to(const Substructure_Atom * a) const
{
  for (const Substructure_Bond* b : _bonds)
  {
    if (a == b->a())
      return 1;
  }

  return 0;
}

int
Substructure_Atom::count_attributes_specified() 
{
  int rc = Substructure_Atom_Specifier::count_attributes_specified();

  for (auto* c : _components)
  {
    rc += c->count_attributes_specified();
  }

  for (auto* c : _children)
  {
    rc += c->count_attributes_specified();
  }

  for (auto* e : _environment)
  {
    rc += e->count_attributes_specified();
  }

  return rc;
}

int
Substructure_Atom::spinach_match_specified() const
{
  if (Substructure_Atom_Specifier::spinach_match_specified())
    return 1;

  for (const Substructure_Atom* c : _children)
  {
    if (c->spinach_match_specified())
      return 1;
  }

  return 0;
}

int
Substructure_Atom::fragment_id (int & result) const
{
  result = 0;   // make sure initialised

  if (0 != _fragment_id)
  {
    result = _fragment_id;
    return 1;
  }

  for (const Substructure_Atom* c : _children)
  {
    if (c->fragment_id(result))
      return 1;
  }

  return 0;      // no fragment_id with this atom or its descendants
}

/*
  Is a given query atom connected to another?
*/

/*int
Substructure_Atom::bonded_to (int j) const
{
  assert (ok());

  int nb = _bonds.number_elements();
  for (int i = 0; i < nb; i++)
  {
    Substructure_Atom * s = other(i);
    if (s->unique_id() == j)
      return 1;
  }

  return 0;
}*/

/*
  Is this query atom bonded to S, and if so, return the bond which makes
  that connection.
*/

/*Substructure_Bond *
Substructure_Atom::bond_to (const Substructure_Atom * s) const
{
  int nb = _bonds.number_elements();
  for (int i = 0; i < nb; i++)
  {
    Substructure_Bond * b = _bonds[i];
    if (s == b->other(this))
      return b;
  }

  return nullptr;
}*/

//#define SHOW_SUBSTRUCTURE_ATOM_OK

int
Substructure_Atom::ok () const
{
#ifdef SHOW_SUBSTRUCTURE_ATOM_OK
  if (nullptr != _parent && nullptr == _bond_to_parent)
  {
    cerr << "Substructure_Atom::ok: parent is not null, but bond to parent is, atom " << _unique_id << '\n';
  }
#endif

  if (nullptr != _parent)
  {
    if (nullptr == _bond_to_parent)
      return 0;
  }

#ifdef SHOW_SUBSTRUCTURE_ATOM_OK
  if (nullptr == _parent && nullptr != _bond_to_parent)
  {
    cerr << "Yikes, has no parent but _bond_to_parent exists!!\n";
    _bond_to_parent->debug_print(cerr, ' ');
  }
#endif

  if (nullptr == _parent && nullptr != _bond_to_parent)
    return 0;

#ifdef SHOW_SUBSTRUCTURE_ATOM_OK
  if (_current_hold_atom)
    cerr << "Atom " << _unique_id << " checking hold atom " << _current_hold_atom->ok() << '\n';
#endif

  if (_current_hold_atom && ! _current_hold_atom->ok())
    return 0;

#ifdef SHOW_SUBSTRUCTURE_ATOM_OK
  if (! _operator.ok())
    cerr << "Operator not OK\n";
#endif

  if (! _operator.ok())
    return 0;

  if (_components.empty() && 1 == _operator.number_results())
    ;
  else if (_components.number_elements() == _operator.number_results() - 1)    // true only during building
    ;
  else if (_components.number_elements() != _operator.number_results())
  {
#ifdef SHOW_SUBSTRUCTURE_ATOM_OK
    cerr << "Components " << _components.number_elements() << " operator mismatch " << _operator.number_results() << '\n';
#endif
    return 0;
  }

#ifdef SHOW_SUBSTRUCTURE_ATOM_OK
  cerr << "Checking underlying specifier\n";
#endif

  return Substructure_Atom_Specifier::ok();
}

int
Substructure_Atom::ok_recursive() const
{
  if (! ok())
    return 0;

  for (const Substructure_Atom* c : _children)
  {
    if (! c->ok_recursive())
      return 0;
  }

  return 1;
}

int
Substructure_Atom::debug_print(std::ostream & os, const IWString & indentation) const
{
  assert (os.good());

  os << indentation << "Substructure Atom " << _unique_id << " (initial " << _initial_atom_number << ")";
  if (_or_id)
    os << " (or id " << _or_id << ")";

  os << " with " << _children.number_elements() << " children, " << 
        _components.number_elements() << " components\n";

  os << indentation << "Operators ";
  _operator.debug_print(os);
  os << '\n';

  if (! ok())
    os << indentation << "Warning, OK fails\n";

  Substructure_Atom_Specifier::debug_print(os, indentation);

  if (nullptr != _anchor)
  {
    os << indentation << " Currently anchored at atom " << _anchor->atom_number() << ", _con = " << _con << '\n';
  }
  else
  {
    os << indentation << " Anchor not set\n";
  }

  if (nullptr != _current_hold_atom)
    os << indentation << " Currently matched with atom " << _current_hold_atom;
  else
    os << indentation << " Not currently matched";

  if (_parent)
    os << ", child of " << _parent->unique_id();
  else
    os << ", no parent";

  IWString ind;
  ind << indentation << "  ";

  int nb = _bonds.number_elements();
  if (0 == nb)
    os << ", no bonds\n";
  else
  {
    os << ", " << nb << " bonds\n";

    for (int i = 0; i < nb; i++)
    {
      Substructure_Bond * b = _bonds[i];
      os << indentation << "    " << i << " ";
      b->debug_print(os, ind);

      Substructure_Atom * other_query_atom = b->a();
      if (nullptr == other_query_atom)
        os << indentation << "Very strange, atom at other end of bond is NULL, perhaps I am an environment atom\n";
      else
      {
        const Target_Atom * other_atom = other_query_atom->current_hold_atom();

        if (nullptr == other_atom)
          os << indentation << "    Substructure atom at other end (" << other_query_atom->unique_id() << ") not matched\n";
        else
          os << indentation << "    Substructure atom at other end (" << other_query_atom->unique_id() << ") is matched with atom " <<
                other_query_atom->current_hold_atom()->atom_number() << '\n';
      }
    }
  }

  int nc = _components.number_elements();
  for (int i = 0; i < nc; i++)
  {
    const Substructure_Atom_Specifier * c = _components[i];
    os << indentation << "  Component " << i << '\n';
    c->debug_print(os, ind);
  }

  if (0 == _include_in_embedding)
    os << indentation << "  Not included in embedding\n";

  return os.good();
}

int
Substructure_Atom::recursive_debug_print (std::ostream & os, 
                              const IWString & indentation) const
{
  debug_print(os, indentation);

  int nc = _children.number_elements();
  for (int i = 0; i < nc; i++)
  {
    os << indentation << "Child " << i << " of " << _unique_id << '\n';
    _children[i]->recursive_debug_print(os, indentation + "  ");
  }

  return os.good();
}

int
Substructure_Atom::terse_details (std::ostream & os,
                      const IWString & indentation) const
{
  assert (os.good());

  os << indentation << "Substructure Atom " << _unique_id << " (initial " << _initial_atom_number << ")";
  if (_or_id)
    os << " (or id " << _or_id << ")";

  os << " with " << _children.number_elements() << " children";

  if (_components.number_elements())
    os << " and " << _components.number_elements() << " components";

  os << '\n';

  if (! ok())
    os << "Warning, OK fails\n";

  os << indentation << "Substructure atom specifications\n";
  if (_element.number_elements())
  {
    os << indentation << " atomic number:";
    for (int i = 0; i < _element.number_elements(); i++)
    {
      os << ' ' << _element[i]->atomic_number();
    }
    os << '\n';
  }
  if (_ncon.is_set())
    os << indentation << " ncon " << _ncon << '\n';
  if (_nbonds.is_set())
    os << indentation << " nbonds " << _nbonds << '\n';
  if (_ncon2.is_set())
    os << indentation << " ncon2 " << _ncon2 << '\n';
  if (_formal_charge.is_set())
    os << indentation << " Formal charge " << _formal_charge << '\n';
  if (_nrings.is_set())
    os << indentation << " Nrings " << _nrings << '\n';
  if (_hcount.is_set())
    os << indentation << " Hcount " << _hcount << '\n';
  if (SUBSTRUCTURE_NOT_SPECIFIED != _aromaticity)
    os << indentation << " Aromaticity = " << _aromaticity << '\n';

  if (nullptr != _current_hold_atom)
    os << indentation << "Currently matched with atom " << _current_hold_atom << '\n';

  if (_bonds.number_elements())
    os << _bonds.number_elements() << " ring closure bonds";

  if (nullptr != _anchor)
  {
    os << indentation << "Currently anchored at atom " << _anchor->atom_number() << ", _con = " << _con << '\n';
  }

  int nc = _components.number_elements();
  for (int i = 0; i < nc; i++)
  {
    const Substructure_Atom_Specifier * c = _components[i];
    os << indentation << "  Component " << i << '\n';
    c->terse_details(os, indentation + "  ");
  }

  return os.good();
}

int
Substructure_Atom::recursive_terse_details(std::ostream & os,
                      const IWString & indentation) const
{
  assert (os.good());

  terse_details(os, indentation);

  int nc = _children.number_elements();
  for (int i = 0; i < nc; i++)
  {
    os << indentation << "Child " << i << " of " << _unique_id << '\n';
    _children[i]->recursive_terse_details(os, indentation + "  ");
  }

  return os.good();
}

int 
Substructure_Atom::print_connectivity_graph(std::ostream & output) const
{
  const int uid = 10000 * _unique_id + 1000 * _initial_atom_number + 100 * _atom_map_number;
  output << "atom " << _unique_id << " initial " << _initial_atom_number << " map " << _atom_map_number << ", UID " << uid << ", " << _children.number_elements() << " children. ";
  if (nullptr == _bond_to_parent)
    output << " Root\n";
  else if (_bonds.number_elements() > 1)
  {
    output << _bonds.number_elements() << " bonds\n";
    for (int i = 0; i < _bonds.number_elements(); ++i)
    {
      output << " bond to " << _bonds[i]->a()->unique_id() << '\n';
    }
  }
    
  for (int i = 0; i < _children.number_elements(); ++i)
  {
    output << " child " << i << " of " << uid << '\n';
    _children[i]->print_connectivity_graph(output);
  }

  return 1;
}

int
Substructure_Atom::print_current_hold (std::ostream & os) const
{
  assert (ok());
  
  os << "Query atom " << _unique_id;
  if (nullptr == _current_hold_atom)
  {
    os << " not bound\n";
    return 0;
  }

  os << " hold is " << _current_hold_atom->atom_number() << '\n';
  return 1;
}

int
Substructure_Atom::involves_rings() const
{
  assert (nullptr == "This is not working");

  if (Substructure_Atom_Specifier::involves_rings())
    return 1;

  for (const Substructure_Atom_Specifier* c : _components)
  {
    if (c->involves_rings())
      return 1;
  }

  return 0;
}

//#define DEBUG_HOLD

/*
  This Substructure_Atom is matched with atom A
*/

int
Substructure_Atom::set_hold (Target_Atom * a, int * already_matched)
{
//assert (ok());

  assert (nullptr == _current_hold_atom);

  int anum = a->atom_number();
  assert (0 == already_matched[anum]);

#ifdef DEBUG_HOLD
  cerr << "Query atom " << _unique_id << " set hold to " << anum << '\n';
#endif

  _current_hold_atom = a;

  already_matched[anum] = 1;

  return 1;
}

int
Substructure_Atom::release_hold (int * already_matched)
{
//assert (ok());

  assert (nullptr != _current_hold_atom);

  assert (_current_hold_atom->ok());

#ifdef DEBUG_HOLD
  cerr << "Atom " << _unique_id << " being released from atom " << _current_hold_atom->atom_number() << '\n';
#endif

  int anum = _current_hold_atom->atom_number();
  assert (already_matched[anum]);

  already_matched[anum] = 0;

  _current_hold_atom = nullptr;

  return 1;
}

int
Substructure_Atom::recursive_release_hold()
{
  if (nullptr == _current_hold_atom)
    return 1;

  _current_hold_atom = nullptr;

//#define IW_HAVE_LAMBDA
#ifdef IW_HAVE_LAMBDA
  _children.each([] (Substructure_Atom & c) { c.recursive_release_hold();});
  _environment.each([] (Substructure_Atom & e) { e.recursive_release_hold();});
#else
  for (int i = _children.number_elements() - 1; i >= 0; --i)
  {
    _children[i]->recursive_release_hold();
  }

  for (int i = _environment.number_elements() - 1; i >= 0; --i)
  {
    _environment[i]->recursive_release_hold();
  }
#endif

  return 1;
}

atom_number_t
Substructure_Atom::atom_number_matched () const
{
  assert (_current_hold_atom);

  return _current_hold_atom->atom_number();
}

/*
  A molecule has been made from the query. Each component of the query
  needs to ensure that it does not violate any aspect of the molecule.

  Note that the Molecule may have 0 == nrings for any atom, but the
  query itself may specify a positive ring value, so all we really
  can do is to check that any max_nrings specified is at least as
  large as nr
*/

int
Substructure_Atom::_adjust_nrings (int nr)
{
//  There is nothing we can do in this case. Even if the atom in the query
//  does not appear in a ring, the query can ultimately be embedded in a ring.

  if (0 == nr && _nrings.is_set())
    return 1;
  
  int mxnr;
  if (! _nrings.max(mxnr))
    return 1;

  if (mxnr < nr)
  {
    cerr << "Substructure_Atom::_adjust_nrings: nrings violation\n";
    cerr << "Perceived " << nr << " rings in query, max is " << mxnr << '\n';
    if (! _nrings.adjust_to_accommodate(nr))
    {
      cerr << "Could not adjust nr " << _nrings << '\n';
    }
  }


  return 1;
}

/*
  The query has been found to be cyclic. Each Substructure_Atom must ensure
  that its own _ring_sizes specifications are consistent with the query.
*/

int
Substructure_Atom::_adjust_ring_sizes (const List_of_Ring_Sizes & ring_sizes_perceived)
{
  if (ring_sizes_perceived.empty())
    return 1;

  if (! _ring_size.is_set())
  {
    for (int r : ring_sizes_perceived)
    {
      _ring_size.add(r);
    }

    return 1;
  }

  for (const auto r : ring_sizes_perceived) {
    _ring_size.add_if_not_already_present(r);
  }

  return 1;
}

/*
  A ring has been detected in the query. Each substructure atom must check to 
  ensure its own consistency with the ring system in the query.

  There are two parts of ring specification
    _nrings and _ring_sizes
*/

int
Substructure_Atom::adjust_ring_specifications (int nr, 
                                  const List_of_Ring_Sizes & ring_sizes)
{
  assert (nr >= 0);

  if (! _adjust_nrings(nr))
    return 0;

  if (! _adjust_ring_sizes(ring_sizes))
    return 0;

  return 1;
}

/*
  Note that this is ambiguous. We are really returning the
  first specification of an atomic number encountered in the query atom.
  So a query like [CD3,ND2] would be poorly served by this.
  Be watching for problems...
*/

int
Substructure_Atom::atomic_number (atomic_number_t & z) const
{
  if (Substructure_Atom_Specifier::atomic_number(z))
    return 1;

  for (Substructure_Atom_Specifier* c : _components)
  {
    if (c->atomic_number(z))
      return 1;
  }

  return 0;
}

/*
  Substructure_Atom::min_ncon is only used by ::create_molecule. It returns an
  estimate of the minimum connectivity associated with an atom.

  A floor on the return code is the actual number of connections present in the
  query.

  If a minimum has been specified, return that value.
  Otherwise, check to see if there is an explicit connectivity available.
*/

int
Substructure_Atom::min_ncon () const
{
  int minimum_possible = _children.number_elements();

  int tmp = Substructure_Atom_Specifier::min_ncon();
  if (tmp > minimum_possible)
    minimum_possible = tmp;

  for (const Substructure_Atom_Specifier* c : _components)
  {
    int tmp = c->min_ncon();
    if (tmp > minimum_possible)
      minimum_possible = tmp;
  }

  return minimum_possible;
}


int
Substructure_Atom::min_nbonds () const
{
  int minimum_possible = min_ncon();

  int tmp = Substructure_Atom_Specifier::min_nbonds();
  if (tmp > minimum_possible)
    minimum_possible = tmp;

  for (const Substructure_Atom_Specifier* c : _components)
  {
    int tmp = c->min_nbonds();
    if (tmp > minimum_possible)
      minimum_possible = tmp;
  }

  return minimum_possible;
}

/*
  In order to match a query atom, we first check any properties of the
  Substructure_Atom itself, followed by the individual
  Substructure_Atom_Specifiers
*/

// #define DEBUG_ATOM_MATCHES

int
Substructure_Atom::_matches(Target_Atom & target, const int * already_matched)
{
  assert (nullptr == _current_hold_atom);
  assert (0 == already_matched[target.atom_number()]);

  if (_attributes_specified > 0) {
    const int m = Substructure_Atom_Specifier::matches(target);

#ifdef DEBUG_ATOM_MATCHES
    cerr << "Matching atom with " << _components.number_elements() << " components, underlying specifier match " << m << " match or rej " << _match_as_match_or_rejection << '\n';
#endif
  
    if (m && 0 == _match_as_match_or_rejection) {    // matches, but we are a rejection criterion
      return 0;
    }

    if (0 == m && 0 != _match_as_match_or_rejection) {    // no match, but we must match
      return 0;
    }

    if (0 == m) {
      return ! _match_as_match_or_rejection;
    }
  }

#ifdef DEBUG_ATOM_MATCHES
  cerr << "Query atom " << _unique_id << " and " << target.atom_number() <<
          " not rejected by global conditions\n";
#endif       

  _preference_value_current_match = _preference_value;

// Process all the components. Note that with more complex sets of operators
// than just OR, it is not entirely clear what to do with the preference
// value. For now, we sum it for all components which match - until the
// expression is evaluated that is. Note also the potential for problems
// with components which match as exclusions. Beware....

  int nc = _components.number_elements();

  if (nc)
  {
    int result = -1;

    _operator.reset();
#ifdef DEBUG_ATOM_MATCHES
    cerr << "Starting evaluation of " << nc << " components\n";
    _operator.debug_print(cerr);
#endif

    for (int i = 0; i < nc; i++)
    {
//    cerr << "Do we need a result from component " << i << "? " << _operator.result_needed(i) << '\n';

      if (! _operator.result_needed(i))
        continue;
  
      int tmp = _components[i]->matches(target);
      _operator.set_result(i, tmp);

#ifdef DEBUG_ATOM_MATCHES
      cerr << "Component " << i << " result is " << tmp << '\n';
#endif

      if (tmp)
        _preference_value_current_match += _components[i]->preference_value();

      if (_operator.evaluate(result))
      {
        if (0 == result)       // the expression is false
          return ! _match_as_match_or_rejection;

        break;
      }
#ifdef DEBUG_ATOM_MATCHES
      cerr << "Could not evaluate operator ";
      _operator.debug_print(cerr);
#endif
    }
  }

  if (_environment.number_elements())
  {
#ifdef DEBUG_ATOM_MATCHES
    cerr << "Checking " << _environment.number_elements() << " environment components\n";
#endif
    if (! _environment.matches(target, target.atoms_in_target_molecule(), already_matched))
      return ! _match_as_match_or_rejection;

#ifdef DEBUG_ATOM_MATCHES
    cerr << "Environment matches, continuing\n";
#endif
  }

#ifdef DEBUG_ATOM_MATCHES
  cerr << "Query atom " << _unique_id << " trying to match " << target.atom_number() << " environment OK\n";
  cerr << "_global_match_id " << _global_match_id << " cmp target " << target.global_id() << '\n';
#endif

  if (_global_match_id > 0 && _global_match_id != target.global_id()) {
    return ! _match_as_match_or_rejection;
  }

// sum any numeric preferences.
// Any zero preference value means reject the match!

  if (_preferences.empty()) {
    return 1;
  }

//#define DEBUG_PREFERENCE_STUFF
#ifdef DEBUG_PREFERENCE_STUFF
  cerr << "Checking " << np << " preferences, sum " << _sum_all_preference_hits << 
          " matched atom " << target.atom_number() << " type " << target.atomic_number() << '\n';
#endif

  for (Substructure_Atom_Specifier* a : _preferences)
  {
    if (a->matches(target))
    {
#ifdef DEBUG_PREFERENCE_STUFF
      cerr << "Preference component " << i << " matches, value " << a->preference_value() << '\n';
#endif
      if (0 == a->preference_value())
        return ! _match_as_match_or_rejection;

      _preference_value_current_match += a->preference_value();

//    Do we sum all preference hits, or just take the first one

      if (0 == _sum_all_preference_hits)
        break;
    }
  }

#ifdef DEBUG_PREFERENCE_STUFF
  cerr << "Preference sum " << _preference_value_current_match << '\n';
#endif

  return 1;
}

/*
  Public interface for _matches
*/

int
Substructure_Atom::matches(Target_Atom & target, const int * already_matched)
{
#ifdef DEBUG_ATOM_MATCHES
  cerr << "Query atom " << _unique_id;
  if (_parent)
    cerr << " (anchor matched with " << _parent->current_hold_atom()->atom_number() << ")";
  cerr << " trying to match atom " << target.atom_number() << ' ' << target.m()->atomic_symbol(target.atom_number()) << '\n';
#endif

  int rc = _matches(target, already_matched);

#ifdef DEBUG_ATOM_MATCHES
  if (rc)
  {
    cerr << "Returning match: query atom " << _unique_id << " with " <<
             target.atom_number() << " OK\n";
  }
#endif

  return rc;
}

/*
  We are thinking of matching THIS atom with target_atom A.
  We must first decide whether or not all the ring closure bonds
  associated with THIS can be satisfied.
*/

int
Substructure_Atom::_ring_closure_bonds_match(Target_Atom * a) const
{
  int nb = _bonds.number_elements();
  if (1 == nb)     // just the bond to parent, no ring closure bonds
    return 1;

  for (int i = 1; i < nb; i++)    // start with i = 1
  {
    Substructure_Bond * b = _bonds[i];
    const Substructure_Atom * s = b->a();
    if (! s->is_matched())
    {
      cerr << "Substructure_Atom::_ring_closure_bonds_match: bond ordering error\n";
      cerr << "Query atom " << _unique_id << " being placed, but query atom " <<
               s->unique_id() << " not yet matched\n";
      assert (nullptr == "Change bond ordering in query file\n");
      return 0;
    }

    Bond_and_Target_Atom * zbond;
    if (nullptr == (zbond = a->bond_to_atom(s->atom_number_matched())))    // not connected
      return 0;

    if (! b->matches(*zbond))
      return 0;
  }

  return 1;    // none of the bonds did not match
}

int
Substructure_Atom::set_anchor (Target_Atom * a)
{
//assert (a->ok());

  _con = 0;

  _anchor = a;

  _anchor_ncon = _anchor->ncon();

  return 1;
}

/*
  This is called whenever an atom is placed on the pending list
  All conditions are initialised to 'nothing'
*/

int
Substructure_Atom::prepare_for_matching (Target_Atom * new_anchor)
{
  _release_hold();

  set_anchor(new_anchor);

  return 1;
}

void
Substructure_Atom::_release_hold ()
{
  _current_hold_atom = nullptr;

  _preference_value_current_match = 0;

  return;
}

/*
  set_parent () is used by Substructure_Environment. In that case, we don't
  add the bond to our array - it is owned by someone else.
  If the ADD_BOND_TO_BONDS_ARRAY is set, then we assume ownership of the bond
*/

void
Substructure_Atom::set_parent (Substructure_Atom * a, 
                               Substructure_Bond * b,
                               int add_bond_to_bonds_array)
{
  _parent = a;
  _bond_to_parent = b;

  if (nullptr == b->a())
    b->set_atom(a);

  if (add_bond_to_bonds_array)
    _bonds.add(b);

  return;
}

/*
  We are proposing to match THIS with TARGET. 
  Before we do that, we must make sure that no _ring_id's would be
  violated.
  Examine all the already placed atoms and check their _ring_id values.
  When set, check for consistency with THIS
*/

int
Substructure_Atom::_match_all_ring_ids (Target_Atom * target,
                      const Query_Atoms_Matched & matched_atoms)
{
#ifdef DEBUG_MATCH_ALL_RING_IDS
  cerr << "Checking " << na << " matched atons for ring id's\n";
#endif

  for (const Substructure_Atom* a : matched_atoms)
  {
    if (0 == a->ring_id())
      continue;

    Target_Atom * am = a->current_hold_atom();
    if (nullptr == am)     // we have covered all the matched atoms
      return 1;

    if (_ring_id == a->ring_id())
    {
      if (! target->in_same_rings(am))
        return 0;
    }
    else if (target->in_same_rings(am))
      return 0;
  }

  return 1;
}


/*
  A Substructure_Atom which may or may not already be matched, must move
  on to the next atom.

  Its parent must be matched for this to work!!

  When we have multi-root queries, we must allow the possibility for
  there to be no anchor, in which case, this always fails
*/

//#define DEBUG_MOVE_TO_NEXT_FROM_ANCHOR

int
Substructure_Atom::move_to_next_match_from_current_anchor (int * already_matched,
                             const Query_Atoms_Matched & matched_atoms)
{
//assert (ok());

#ifdef DEBUG_MOVE_TO_NEXT_FROM_ANCHOR
  cerr << "Query atom " << _unique_id << " moving.";
  if (_current_hold_atom)
    cerr << " current hold atom " << _current_hold_atom->atom_number() << '\n';
  else
    cerr << " no current hold atom\n";
#endif

  if (_current_hold_atom)
    release_hold(already_matched);

  if (nullptr == _anchor)
    return 0;

#ifdef DEBUG_MOVE_TO_NEXT_FROM_ANCHOR
  cerr << "move_to_next_match_from_current_anchor:: query atom " << _unique_id << " anchored at " << _parent->current_hold_atom()->atom_number() << '\n';
  cerr << "_con = " << _con << " and _anchor_ncon = " << _anchor_ncon << '\n';
#endif

// Loop through all remaining connections to anchor

  for ( ; _con < _anchor_ncon; _con++)
  {
    Bond_and_Target_Atom & bata = _anchor->other(_con);

    Target_Atom * a = bata.other();

#ifdef DEBUG_MOVE_TO_NEXT_FROM_ANCHOR
    cerr << "Query atom " << _unique_id << ", _con " << _con << ", may go to atom " << a->atom_number() << " atn " << a->atomic_number() << '\n';
    if (already_matched[a->atom_number()])
      cerr << ". Nope, that one's matched\n";
    else
      cerr << '\n';
#endif

    if (already_matched[a->atom_number()])
      continue;

#ifdef DEBUG_MOVE_TO_NEXT_FROM_ANCHOR
    if (! _bond_to_parent->matches(bata))
    {
      cerr << "Bond type mismatch\n";
      _bond_to_parent->debug_print(cerr, "MISMATCH");
    }
#endif

    if (! _bond_to_parent->matches(bata))
      continue;

#ifdef DEBUG_MOVE_TO_NEXT_FROM_ANCHOR
    if (! matches(*a, already_matched))
      cerr << "Cannot match that atom\n";
#endif

    if (! matches(*a, already_matched))
      continue;

    // If only the bond to the parent, no need to check ring closure bonds.
    if (_bonds.number_elements() == 1) {
    } else if (! _ring_closure_bonds_match(a)) {
      continue;
    }

//  Looking good. Check any _ring_id values

    if (_ring_id && ! _match_all_ring_ids(a, matched_atoms))
      continue;

//  We have a match!

    set_hold(a, already_matched);

    _con++;

#ifdef DEBUG_MOVE_TO_NEXT_FROM_ANCHOR
    cerr << "move_to_next_match: query atom " << _unique_id << ", _con " << _con << " matched with atom " << 
            a->atom_number() << '\n';
#endif

    return 1;
  }

// One might think that once we had exhausted all the connections to a
// given anchor, we should invalidate that anchor, but instead, we need
// to reset in case someone is trying a different arrangement of the
// connected atoms about that anchor

  _con = 0;      // get ready for another attempt with same anchor

#ifdef DEBUG_MOVE_TO_NEXT_FROM_ANCHOR
  cerr << "Move_to_next_match_from_current_anchor: returning 0, query atom " << _unique_id << '\n';
#endif

  return 0;
}

/*
  This function was written for the fingerprint routines.
  They need to know whether or not this query atom is specified 
  as being in a ring or not
*/

int
Substructure_Atom::determine_ring_or_non_ring (int & result) const
{
  if (! _nrings.is_set())
    return 0;

  const int nelements = _nrings.number_elements();

  if (nelements == 1) {
    result = _nrings[0];
    return 1;
  }

  int tmp;
  if (_nrings.min(tmp) && tmp > 0) {
    result = tmp;
    return 1;
  }

  // must be just a max value specified.
  if (nelements == 0) {
    return 0;
  }

// There are multiple values for nrings. If they are all > 0, then ok.
// Fingerprint just needs to know ring or non ring.

  for (int i = 0; i < nelements; ++i) {
    if (_nrings[i] == 0) {
      return 0;
    }
  }

// All _nrings values were > 0

  result = _nrings[0];
  return 1;
}

/*
  This function was written for the fingerprint routines.
  Does an explicit ncon value exist for this atom.
*/

int
Substructure_Atom::ncon_specification (int & result) const
{
  if (! _ncon.is_set())
    return 0;

  if (_ncon.number_elements() > 1)
    return 0;

  if (1 == _ncon.number_elements())
  {
    result = _ncon[0];
    return 1;
  }

  return 0;
}

/*
  When making a molecule from a query, we need to know whether or not there
  is a unique hcount specification for an atom
*/

int
Substructure_Atom::hcount_specification (int & result) const
{
  if (! _hcount.is_set())
    return 0;

  if (_hcount.number_elements() > 1)
    return 0;

  if (1 == _hcount.number_elements())
  {
    result = _hcount[0];
    return 1;
  }

  return 0;
}

int
Substructure_Atom::formal_charge_specification (formal_charge_t & result) const
{
  if (! _formal_charge.is_set())
    return 0;

  if (_formal_charge.number_elements() > 1)
    return 0;

  if (1 == _formal_charge.number_elements())
  {
    result = _formal_charge[0];
    return 1;
  }

  return 0;
}

// why do both of these methods exist?

int
Substructure_Atom::formal_charge (formal_charge_t & fc) const
{
  assert (ok());

  if (! _formal_charge.is_set())
    return 0;

  if (_formal_charge.number_elements())
  {
    fc = _formal_charge[0];
    return 1;
  }

//  Should probably do a check for min and max. Have never used them, so...

  return 0;
}

/*
  A Substructure_Query object can do a consistency check upon itself.
  All components must support a check_internal_consistency function.
*/

int
Substructure_Atom::check_internal_consistency(int connections) const
{
  (void) connections;
  return 1;
}

int
Substructure_Atom::ring_sizes_specified(resizable_array<int> & ring_sizes) const
{
  int n = _aromatic_ring_size.number_elements();
  for (int i = 0; i < n; ++i) {
    ring_sizes.add_if_not_already_present(_aromatic_ring_size[i]);
  }

  n = _aliphatic_ring_size.number_elements();
  for (int i = 0; i < n; ++i) {
    ring_sizes.add_if_not_already_present(_aliphatic_ring_size[i]);
  }

  return ring_sizes.number_elements();
}

/*
  Nov 97: I had previously ignored OR conditions between the children,
  but ran into a problem with that. Therefore we take special care
  when we have an OR condition.

  Note that we don't include the environment.
*/

int
Substructure_Atom::min_atoms_for_match() const
{
//assert (ok ());

  int nc = _children.number_elements();

  if (0 == nc)
    return 1;

  int * already_done = new_int(nc); std::unique_ptr<int[]> free_already_done(already_done);
  
  int rc = 1;     // this atom
  for (int i = 0; i < nc; i++)
  {
    if (already_done[i])
      continue;

    int tmp = _children[i]->min_atoms_for_match();
    if (0 == _children[i]->or_id())
    {
      rc += tmp;
      continue;
    }

//  Now things get complicated. This child is in an OR relationship

    for (int j = i + 1; j < nc; j++)
    {
      if (already_done[j])
        continue;

      if (_children[j]->or_id() != _children[i]->or_id())
        continue;

      already_done[j] = 1;

      int tmpj = _children[j]->min_atoms_for_match();
      if (tmpj < tmp)
        tmp = tmpj;
    }

    rc += tmp;
  }

  return rc;
}

int
Substructure_Atom::_set_implicit_hydrogens (Molecule & m,
                                            atom_number_t a) const
{
  int ih;
  if (! hcount_specification(ih))
  {
//  cerr << "Setting implicit hcount for atom " << i << " to zero\n";
    m.set_implicit_hydrogens(a, 0);
  }
  else
    m.set_implicit_hydrogens(a, ih);

  return 1;
}

/*
  Called during create_molecule_from_query
*/

int
Substructure_Atom::_fill_min_ncon (Molecule & m,
                                   atom_number_t a) const
{
  int connection_shortage = min_ncon() - m.ncon(a);
  if (0 == connection_shortage)
    return 1;

// A negative connection shortage means that this query atom does not know
// about its full connectivity. The query should probably be changed to
// reflect this.

  if (connection_shortage < 0)
  {
    cerr << "Substructure_Atom::_fill_min_ncon: Query atom " << _unique_id << 
            " connections not specified fully\n";
    cerr << "Min known is " << min_ncon() << " but molecule has " << m.ncon(a) << '\n';
    cerr << "Suggest you change the query to reflect that\n";
  }

  int bond_shortage = min_nbonds() - m.nbonds(a);
  if (bond_shortage <= 0)
    bond_shortage = 0;

// cerr << "Atom " << a << " cs = " << connection_shortage << " bs = " << bond_shortage << '\n';
  for (int i = 0; i < connection_shortage; i++)
  {
    Atom * h = new Atom(0);
    m.add(h);
    bond_type_t bt = SINGLE_BOND;
    if (bond_shortage > connection_shortage)
    {
      bt = DOUBLE_BOND;
      bond_shortage--;
    }
    m.add_bond(a, m.natoms() - 1, bt);
  }

  return 1;
}

//#define DEBUG_ADD_BOND_TO_MOLECULE

int
Substructure_Atom::_add_bond_to_molecule (Molecule & m,
                                          const Substructure_Atom * other_atom,
                                          const Substructure_Bond * b) const
{
  bond_type_t types_matched = b->types_matched();

  bond_type_t bt;

  if (SINGLE_BOND & types_matched)
    bt = SINGLE_BOND;
  else if (DOUBLE_BOND & types_matched)
    bt = DOUBLE_BOND;
  else if (TRIPLE_BOND & types_matched)
    bt = TRIPLE_BOND;
  else
    bt = SINGLE_BOND;

  assert (other_atom->ok());

#ifdef DEBUG_ADD_BOND_TO_MOLECULE
  cerr << "_add_bond_to_molecule: molecule has " << m.natoms() << " atoms\n";
  cerr << "Adding bond between " << _unique_id << " and " << other_atom->unique_id() << '\n';
#endif

  if (_unique_id == other_atom->unique_id())
  {
    cerr << "Substructure_Atom::_add_bond_to_molecule: yipes, atoms the same " << _unique_id << '\n';
    return 0;
  }

  int matoms = m.natoms();

  if (_unique_id >= matoms || other_atom->unique_id() >= matoms)
  {
    cerr << "Substructure_Atom::_add_bond_to_molecule: yipes, atom numbers out of range\n";
    cerr << "molecule contains " << matoms << " atoms, id's are " << _unique_id << " and " << other_atom->unique_id() << '\n';
    return 0;
  }

  if (m.are_bonded(_unique_id, other_atom->unique_id()))
  {
    cerr << "Substructure_Atom::_add_bond_to_molecule: yipes, atoms " << _unique_id << " and " << other_atom->unique_id() << " already bonded\n";
    return 0;
  }

  return m.add_bond(_unique_id, other_atom->unique_id(), bt);
}

int
Substructure_Atom::_add_bond_to_molecule (Molecule & m,
                                          const Substructure_Bond * b) const
{
  const Substructure_Atom * other_atom = b->a();

  return _add_bond_to_molecule(m, other_atom, b);
}

Atom *
Substructure_Atom::create_atom () const
{
  if (_element.number_elements())
    return new Atom(_element[0]);

  return new Atom(0);
}

int
Substructure_Atom::create_molecule(Molecule & m, 
                                   int fill_min_ncon, 
                                   int set_implicit_hydrogens) const
{
  assert(ok());

  Atom * a = create_atom();

  int my_atom_number = m.natoms();

  formal_charge_t fc;
  if (formal_charge_specification(fc))
    a->set_formal_charge(fc);

  m.add(a);

// Note that we don't do anything about OR values in the children

  for (const Substructure_Atom* c : _children)
  {
    if (! c->create_molecule(m, fill_min_ncon, set_implicit_hydrogens))
    {
      cerr << "Substructure_Atom::create_molecule: atom " << _unique_id << " atom number " << my_atom_number << " failed\n";
      return 0;
    }
  }

  for (const Substructure_Bond * b : _bonds)
  {
    if (! _add_bond_to_molecule(m, b))
    {
      cerr << "Substructure_Atom::create_molecule: atom " << _unique_id << " atom number " << my_atom_number << " adding " << _bonds.number_elements() << " bonds\n";
      return 0;
    }
  }

  if (fill_min_ncon)
    _fill_min_ncon(m, my_atom_number);

  if (set_implicit_hydrogens)
    _set_implicit_hydrogens(m, my_atom_number);

  return 1;
}

void
Substructure_Atom::assign_unique_atom_numbers (int & id)
{
  _unique_id = id++;

  for (Substructure_Atom* c : _children)
  {
    c->assign_unique_atom_numbers(id);
  }

  return;
}

void
Substructure_Atom::assign_unique_id_from_atom_number_if_set (extending_resizable_array<int> & numbers_in_use)
{
//cerr << "Substructure_Atom::assign_unique_id_from_atom_number_if_set:_initial_atom_number " << _initial_atom_number << " _unique_id " << _unique_id << '\n';

  if (_initial_atom_number >= 0) {     // use it
    _unique_id = _initial_atom_number;
    ++numbers_in_use[_unique_id];
  } else     // find an unused number to assign
  {
    int next_to_assign = -1;
    for (unsigned int i = 0; i < numbers_in_use.size(); ++i)    // find current max index
    {
      if (numbers_in_use[i] > 0)     // number I is already in use
        continue;

      next_to_assign = i;
      break;
    }

    if (next_to_assign > 0)
      next_to_assign++;
    else
      next_to_assign = numbers_in_use.size();

    _unique_id = next_to_assign;
    numbers_in_use[next_to_assign]++;
  }

  // THe environment is separate from the main query.
  for (Substructure_Atom* e : _environment) {
    extending_resizable_array<int> for_env;
    e->assign_unique_id_from_atom_number_if_set(for_env);
  }

  for (unsigned int i = 0; i < _children.size(); ++i)
  {
    _children[i]->assign_unique_id_from_atom_number_if_set(numbers_in_use);
  }

  return;
}

/*
  Make sure we return the number of children added
*/

int
Substructure_Atom::add_your_children(Query_Atoms_Matched & alist)
{
  for (Substructure_Atom* c : _children)
  {
    c->prepare_for_matching(_current_hold_atom);   // tell children their new anchor
    alist.add_if_not_already_present(c);
  }

  return _children.number_elements();
}

/*
  A particular matching does not work, and this atom is being taken
  off the list of matched atoms.
  First, we remove any children we had placed on the pending list.
  Then, we release any hold we may have.
  Finally, we invalidate any info about a previous anchor
*/
int
Substructure_Atom::remove_your_children(Query_Atoms_Matched & atoms,
                                        int * already_matched)
{
  if (_current_hold_atom) {
    release_hold(already_matched);
  }

  const int nc = _children.number_elements();
  if (0 == nc) {
    return 1;
  }

  if (1 == nc)
    atoms.remove_first(_children[0]);
  else if (2 == nc)
    atoms.remove_two_items(_children[0], _children[1]);
  else
  {
    for (Substructure_Atom* c : _children)
    {
      atoms.remove_first(c);
    }
  }


  return 1;
}

/*
  A substrucure query object needs to know whether or not any of its
  components depend on knowing the aromaticity of a bond.
  We return 1 if we do. If one of our bonds depends on ring membership,
  we set that
*/

int
Substructure_Atom::involves_aromatic_bond_specifications(int & need_rings) const
{
  if (_aromatic_ring_size.is_set())
    return 1;

  if (_aliphatic_ring_size.is_set())
    return 1;

  if (_aromaticity != SUBSTRUCTURE_NOT_SPECIFIED)
    return 1;

  for (const auto * a : _components) {
    if (a->aromaticity() != SUBSTRUCTURE_NOT_SPECIFIED)
      return 1;
  }

  for (Substructure_Bond* b : _bonds)
  {
    if (b->involves_aromatic_bond_specification(need_rings))
      return 1;
  }

  for (const Substructure_Atom* c : _children) {
    if (c->involves_aromatic_bond_specifications(need_rings))
      return 1;
  }

  return _environment.involves_aromatic_bond_specifications(need_rings);
}

int
Substructure_Atom::unmatched_atoms_attached_specified() const
{
  int rc = 0;
  if (_unmatched_atoms_attached.is_set())
    rc = 1;

  for (const Substructure_Atom * c : _children)
  {
    rc += c->unmatched_atoms_attached_specified();
  }

  return rc;
}

void
Substructure_Atom::notify_extra_child(Substructure_Atom * a)
{
  assert (ok());
  assert (a && a->ok());

  _children.add(a);

  return;
}

int 
Substructure_Atom::smarts(IWString & s) const
{
  assert (ok());

  s = "not implemented";

  return 1;
}

/*static int
transfer_atoms_with_or_value (Target_Atom * target, 
                              resizable_array<Substructure_Atom *> & env,
                              Query_Atoms_Matched & matched_atoms,
                              int orid)
{
  int ne = env.number_elements ();
  int rc = 0;

  for (int i = 0; i < ne; i++)
  {
    Substructure_Atom * a = env[i];
    if (orid != a->or_id ())
      continue;

//  a->recursive_release_hold ();
    a->prepare_for_matching (target);

    env.remove_item (i);
    i--;
    ne--;
    rc++;
    matched_atoms.add (a);
  }

  return rc;
}*/

int
Substructure_Atom::numeric_value (double & result) const
{
  return _numeric_value.value(result);
}

/*int 
Substructure_Atom::smarts(IWString & s) const
{
  assert (ok ());

  s = '[';       // for simplicity we always put brackets

  int na = _atomic_number.number_elements();
  for (int i = 0; i < _atomic_number.number_elements(); i++)
  {
    if (i)
      s += ',';
    s.append_number(_atomic_number[i]);
  }

  int top_level_attributes = count_attributes_specified();
  int nc = _components.number_elements();

  s += ']';

  return 1;
}*/

/*
*/

int
Substructure_Atom::unmatched_connections(const int * already_matched) const
{
  assert (_current_hold_atom);

  const Atom * a = _current_hold_atom->atom ();

  int acon = a->number_elements ();

  int rc = 0;
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(_current_hold_atom->atom_number(), i);
    if (0 == already_matched[j])
      rc++;
  }

  return rc;
}

//#define DEBUG_DETERMINE_START_STOP

/*
  Just before doing a substructure search, we can restrict the
  range for searches by querying the target for atom types

  Note that if all atomic number info is in components, they
  will be ignored and the whole molecule will be searched
*/

int
Substructure_Atom::determine_start_stop(const Molecule_to_Match & target,
                       int & istart, int & istop) const
{
  int matoms = target.natoms();

  istart = 0;
  istop = matoms;

  if (!_match_as_match_or_rejection)  // Too complex otherwise.
    return 1;

// If we have no info about atomic numbers, we must search the whole target
  if (_element.empty()) {
    return 1;
  }

  const int na = _element.number_elements();

  assert (_element_unique_id.number_elements() == na);

#ifdef DEBUG_DETERMINE_START_STOP
  cerr << "Determining start/stop points, matoms = " << matoms << " na = " << na << '\n';
  for (int i = 0; i < na; i++)
  {
    cerr << " z = " << _element[i]->atomic_number();
  }
  if (na > 0)
    cerr << '\n';
#endif

// If we have no info about atomic numbers, we must search the whole target

  if (0 == na)
    return 1;

  atomic_number_t z = _element[0]->atomic_number();

  if (z < 0)    // we have one or more non periodic table items, too hard
    return 1;

  istart = target.first(z);

#ifdef DEBUG_DETERMINE_START_STOP
  cerr << "First occurrence of atomic number " << z << " at " << istart << '\n';
#endif

  if (istart >= 0)
  {
    istop = target.last(z) + 1;
    if (1 == na)
      return 1;
  }

// Now the more complex case where we have multiple atomic numbers

  for (int i = 1; i < na; i++)   // we already did item[0] above
  {
    z = _element[i]->atomic_number();

    int j = target.first(z);
    if (j < 0)
      continue;

    if (istart < 0 || j < istart)
      istart = j;

    j = target.last(z) + 1;
    if (j > istop)
      istop = j;
  }

  return (istart >= 0);
}

#ifdef VERSION_USING_ELEMENT

int
Substructure_Atom::determine_start_stop(const Molecule_to_Match & target,
                       int & istart, int & istop) const
{
  const int matoms = target.natoms();

  int na = _element.number_elements();

#ifdef DEBUG_DETERMINE_START_STOP
  cerr << "Determining start/stop points, matoms = " << matoms << " na = " << na << '\n';
  for (int i = 0; i < na; i++)
  {
    cerr << " z = " << _element[i]->atomic_number();
  }
  if (na > 0)
    cerr << '\n';
#endif

// If we have no info about atomic numbers, we must search the whole target

  if (0 == na)
  {
    istart = 0;
    istop = matoms;
    return 1;
  }

  atomic_number_t z = _element[0]->atomic_number();

  istart = target.first(z);

#ifdef DEBUG_DETERMINE_START_STOP
  cerr << "First occurrence of atomic number " << z << " at " << istart << '\n';
#endif

  if (istart >= 0)
  {
    istop = target.last(z) + 1;
    if (1 == na)
      return 1;
  }
  else
    istop = matoms;

// Now the more complex case where we have multiple atomic numbers

  for (int i = 1; i < na; i++)   // we already did item[0] above
  {
    z = _element[i]->atomic_number();

    int j = target.first(z);
    if (j < 0)
      continue;

    if (istart < 0 || j < istart)
      istart = j;

    j = target.last(z) + 1;
    if (j > istop)
      istop = j;
  }

  return (istart >= 0);
}
#endif

Substructure_Atom *
Substructure_Atom::query_atom_with_initial_atom_number(atom_number_t a)
{
  if (_initial_atom_number == a)
    return this;

  for (Substructure_Atom* c : _children)
  {
    Substructure_Atom * rc = c->query_atom_with_initial_atom_number(a);
    if (nullptr != rc)
      return rc;
  }

  return nullptr;
}

Substructure_Atom *
Substructure_Atom::query_atom_with_atom_map_number(atom_number_t a)
{
  if (_atom_map_number == a)
    return this;

  for (Substructure_Atom* c : _children)
  {
    Substructure_Atom * rc = c->query_atom_with_atom_map_number(a);
    if (nullptr != rc)
      return rc;
  }

  return nullptr;
}

int
Substructure_Atom::is_bonded_to(atom_number_t o) const
{
  if (nullptr == _parent)
    ;
  else if (0 == _parent->unique_id())
    return 1;

  for (const Substructure_Bond* b : _bonds)
  {
    if (o == b->a()->unique_id())
      return 1;
  }

  return 0;
}
int
Substructure_Atom::add_ncon_preference_object(int n, int p)
{
  Substructure_Atom_Specifier * s = new Substructure_Atom_Specifier;
  s->set_ncon(n);
  s->set_preference_value(p);

  _preferences.add(s);

  return 1;
}

void
Substructure_Atom::identify_atom_numbers(extending_resizable_array<int> & a) const
{
//cerr << "Substructure_Atom::identify_atom_numbers:my number " << _initial_atom_number << '\n';

  if (_initial_atom_number >= 0)
    a[_initial_atom_number]++;

  for (Substructure_Atom* c : _children)
  {
    c->identify_atom_numbers(a);
  }

  return;
}

void
Substructure_Atom::identify_atom_map_numbers(extending_resizable_array<int> & a) const
{
//cerr << "Substructure_Atom::identify_atom_numbers:my number " << _initial_atom_number << '\n';

  if (_atom_map_number >= 0)
    a[_atom_map_number]++;

  for(Substructure_Atom* c : _children)
  {
    c->identify_atom_map_numbers(a);
  }

  return;
}

/*
  This is made complicated because of the Substructure_Atom datastructure.
  When there is a ring closure, only the second Substructure_Atom knows about that
  bond. So in all cases, if we do nt find a match at a given atom, the 
  children must be searched, since somewhere down the 'tree' there may
  be a ring closure bond pointing back to this atom
*/

//#define DEBUG_BOND_BETWEEN_ATOM_MAP_NUMBERS

const Substructure_Bond *
Substructure_Atom::bond_between_atom_map_numbers(int a1, int a2) const
{
#ifdef DEBUG_BOND_BETWEEN_ATOM_MAP_NUMBERS
  cerr << "Substructure_Atom::bond_between_atom_map_numbers:looking for " << a1 << " and " << a2 << " I am " << _atom_map_number << '\n';
#endif

  const Substructure_Bond * rc = _bond_between_atom_map_numbers(a1, a2);

  if (nullptr != rc)
    return rc;

  for (const Substructure_Atom* c : _children)  // go looking in the children.
  {
    rc = c->bond_between_atom_map_numbers(a1, a2);
    if (nullptr != rc)
      return rc;
  }

  return nullptr;
}

const Substructure_Bond *
Substructure_Atom::_bond_between_atom_map_numbers(int a1, int a2) const
{
  if (a1 == _atom_map_number) 
    ;
  else if (a2 == _atom_map_number)    // make sure a1 is associated with THIS
    iwswap(a1, a2);
  else    // we are neither a1 or a2, caller will check children
    return nullptr;

  assert (a1 == _atom_map_number);

  for (const Substructure_Atom* c : _children)
  {
#ifdef DEBUG_BOND_BETWEEN_ATOM_MAP_NUMBERS
    cerr << "Substructure_Atom::bond_between_atom_map_numbers:from " << _atom_map_number << " child is " << c->atom_map_number() << " CHILD\n";
#endif

    if (a2 != c->atom_map_number())
      continue;

    return c->bond_to_parent();
  }

// what about ring closure bonds

#ifdef DEBUG_BOND_BETWEEN_ATOM_MAP_NUMBERS
  cerr << "N = " << _bonds.number_elements() << '\n';
#endif

  for (const Substructure_Bond* b : _bonds)
  {
    const Substructure_Atom * a = b->a();

#ifdef DEBUG_BOND_BETWEEN_ATOM_MAP_NUMBERS
    cerr << "Substructure_Atom::bond_between_atom_map_numbers:from " << _atom_map_number << " to " << a->atom_map_number() << " BOND\n";
#endif

    if (a2 == a->atom_map_number())
      return b;
  }

  return nullptr;
}

/*
  same as with bond_between_atom_map_numbers, we need to be very
  careful about ring closure bonds
*/

const Substructure_Bond *
Substructure_Atom::bond_between_atoms (int a1, int a2) const    
{
  const Substructure_Bond * rc = _bond_between_atoms(a1, a2);
  if (nullptr != rc)
    return rc;

  for (const Substructure_Atom* c : _children)
  {
    rc = c->bond_between_atoms(a1, a2);
    if (nullptr != rc)
      return rc;
  }

  return nullptr;
}

const Substructure_Bond *
Substructure_Atom::_bond_between_atoms(int a1, int a2) const
{
  if (a1 == _initial_atom_number)
    ;
  else if (a2 == _initial_atom_number)
    iwswap(a1, a2);
  else
  {
    return nullptr;
  }
 
  for (const Substructure_Atom* c : _children)
  {
    if (a2 != c->initial_atom_number())
      continue;

    return c->bond_to_parent();
  }

  const int n = _bonds.number_elements();

  for (int i = 1; i < n; ++i)
  {
    const Substructure_Bond * b = _bonds[i];

    const Substructure_Atom * a = b->a();

    if (a2 == a->initial_atom_number())
      return b;
  }

  return nullptr;
}

int
Substructure_Atom::query_atom_with_isotope (int iso) const
{
  if (_isotope.contains(iso))
    return _initial_atom_number;

  for (auto i = 0; i < _components.number_elements(); ++i)
  {
    auto c = _components[i]->isotope();

    if (c.contains(iso))
      return _initial_atom_number;
  }

  for (auto i = 0; i < _children.number_elements(); ++i)
  {
    auto rc = _children[i]->query_atom_with_isotope(iso);

    if (rc >= 0)
      return rc;
  }

  return -1;
}

void
Substructure_Atom::adjust_initial_atom_numbers(const int * xref)
{
  _initial_atom_number = xref[_initial_atom_number];

  for (Substructure_Atom* c : _children)
  {
    c->adjust_initial_atom_numbers(xref);
  }

  return;
}

int
Substructure_Atom::symmetry_groups_present() const
{
//cerr << "Substructure_Atom::symmetry_groups_present:checking " << _symmetry_group << ' ' << this << '\n';

  if (_symmetry_group > 0)
    return 1;

  for (const auto* c : _components)
  {
    if (c->symmetry_group() > 0)
      return 1;
  }

  for (const auto* c : _children)
  {
    if (c->symmetry_groups_present())
      return 1;
  }

  return 0;
}

int
Substructure_Atom::atom_type_groups_present() const
{
  int rc = 0;

  if (_atom_type_group > 0)
    rc = 1;

  for (const auto* c : _children)
  {
    rc += c->atom_type_groups_present();
  }

  return rc;
}

int
Substructure_Atom::first_symmetry_group() const
{
  if (_symmetry_group > 0)
    return _symmetry_group;

  for (const Substructure_Atom_Specifier* c: _components)
  {
    if (c->symmetry_group() > 0)
      return c->symmetry_group();
  }

  return 0;
}

int
Substructure_Atom::assign_atom_map_numbers(int & amap)
{
  int rc = 0;

//cerr << "Substructure_Atom::assign_atom_map_numbers:atom " << _unique_id << " map " << _atom_map_number << " do I need " << amap << '\n';

  if (_atom_map_number < 0)
  {
    amap++;
    _atom_map_number = amap;
    rc++;
  }

  for (Substructure_Atom* c : _children)
  {
    rc += c->assign_atom_map_numbers(amap);
  }

  return rc;
}

int
UnmatchedAtomsAttached(const Molecule& m,
                       int * already_matched,
                       const atom_number_t zatom)
{
  int rc = 0;

  for (const Bond* b : *m.atomi(zatom))
  {
    const atom_number_t i = b->other(zatom);
    if (already_matched[i])
      continue;
    already_matched[i] = 1;
    rc++;
    rc += UnmatchedAtomsAttached(m, already_matched, i);
  }

  return rc;
}

int
Substructure_Atom::unmatched_atoms_attached_matches(const Query_Atoms_Matched& matched_query_atoms,
                                         const int * already_matched,
                                         const Molecule_to_Match& target)
{
  if (! _unmatched_atoms_attached.is_set())
    return 1;

  if (nullptr == _current_hold_atom)  // Not sure if this can happen.
    return 1;

  const atom_number_t matched = _current_hold_atom->atom_number();

  const Molecule& m = *target.molecule();

  int * tmp = new int[m.natoms()]; std::unique_ptr<int[]> free_tmp(tmp);
  std::copy_n(already_matched, m.natoms(), tmp);

  const int u = UnmatchedAtomsAttached(m, tmp, matched);

  return _unmatched_atoms_attached.matches(u);
}

int
Substructure_Atom::add_component(Substructure_Atom_Specifier* spec) {
  _components.add(spec);

  if (_components.size() > 1) {
    _operator.add_operator(IW_LOGEXP_AND);
  }

  return 1;
}
