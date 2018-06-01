#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

/*
  Special processing for the no matched atoms between directive
*/

#include "misc.h"

#include "mdl_molecule.h"
#include "substructure.h"
#include "target.h"

static void
fill_matched_atoms_array(const Query_Atoms_Matched & matched_atoms,
                          int * matched)
{
  int qam = matched_atoms.number_elements();

  for (int i = 0; i < qam; i++)
  {
    const Substructure_Atom * a = matched_atoms[i];
    if (! a->include_in_embedding())
      continue;

    matched[a->atom_number_matched()] = 1;
  }

  return;
}

//#define DEBUG_NO_MATCHED_ATOM_PATH

/*
  Can we find a path with no matched atoms between two atoms
*/

static int
no_matched_atom_path(const Atom ** atoms, atom_number_t destination,
                      int * already_covered_this_path,
                      const int * matched,
                      atom_number_t my_atom)
{
//assert (0 == matched[my_atom]);   NOT TRUE, the first atom will be on the matched atom list

  assert (0 == already_covered_this_path[my_atom]);

  already_covered_this_path[my_atom] = 1;    // put ourselves on the path

  const Atom * a = atoms[my_atom];

#ifdef DEBUG_NO_MATCHED_ATOM_PATH
  cerr << "no_matched_atom_path continues with atom " << my_atom << " type " << a->atomic_symbol() << " ncon " << a->ncon() << " destination " << destination << endl;
#endif

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(my_atom, i);

    if (destination == j)     // bingo. Do this test first as destination will be set in matched[]
      return 1;

    if (already_covered_this_path[j])
      continue;

    if (matched[j])
      continue;

    if (1 == atoms[j]->ncon())     // terminal atom, no sense following it
      continue;

    if (no_matched_atom_path(atoms, destination, already_covered_this_path, matched, j))
      return 1;
  }

  already_covered_this_path[my_atom] = 0;    // take ourselves off the path

  return 0;
}

//#define DEBUG_NO_MATCHED_ATOMS_BETWEEN_SATISFIED

int
Single_Substructure_Query::_no_matched_atoms_between_satisfied(Query_Atoms_Matched & matched_atoms) const
{
  int nnmab = _no_matched_atoms_between.number_elements();

#ifdef DEBUG_NO_MATCHED_ATOMS_BETWEEN_SATISFIED
  cerr << "Testing " << nnmab << " matched atoms between constraints\n";
#endif

  Molecule * m = matched_atoms[0]->current_hold_atom()->m();

  int matoms = m->natoms();

  int * matched = new_int(matoms); std::unique_ptr<int[]> free_matched(matched);

  fill_matched_atoms_array(matched_atoms, matched);

  const Atom ** atoms = new const Atom *[matoms]; std::unique_ptr<const Atom *[]> free_atoms(atoms);
  m->atoms(atoms);

  int * tmp = new_int(matoms);    // used in constructing the path
  std::unique_ptr<int[]> free_tmp(tmp);

  for (int i = 0; i < nnmab; i++)
  {
    const Bond * b = _no_matched_atoms_between[i];

//  the numbers in the Bond object refer to query atoms. Need to convert these to atom
//  numbers in the molecule

//  const Substructure_Atom * qa1 = matched_atoms[b->a1()];
//  const Substructure_Atom * qa2 = matched_atoms[b->a2()];

    assert (matched_atoms[b->a1()]->is_matched());
    assert (matched_atoms[b->a2()]->is_matched());

    atom_number_t a1 = matched_atoms[b->a1()]->current_hold_atom()->atom_number();
    atom_number_t a2 = matched_atoms[b->a2()]->current_hold_atom()->atom_number();

#ifdef DEBUG_NO_MATCHED_ATOMS_BETWEEN_SATISFIED
    cerr << "Testing constraint between matched atoms " << b->a1() << " and " << b->a2() << endl;
    cerr << "Corresponding to atoms " << a1 << " and " << a2 << " in the molecule\n";
#endif

    if (m->fragment_membership(a1) != m->fragment_membership(a2))
      return 1;

    if (! no_matched_atom_path(atoms, a1, tmp, matched, a2))
      return 0;
  }

  return 1;    // all constraints satisfied
}

int
Single_Substructure_Query::_distance_between_root_atoms_satisfied(Query_Atoms_Matched & matched_atoms) const
{
  assert (_distance_between_root_atoms.is_set());

  int nr = _root_atoms.number_elements();

  assert (nr > 1);     // doesn't make sense otherwise

  Molecule * m = matched_atoms[0]->current_hold_atom()->m();

  for (int i = 0; i < nr; i++)
  {
    const Substructure_Atom * ri = _root_atoms[i];
    assert (ri->is_matched());

    atom_number_t ai = ri->current_hold_atom()->atom_number();

    for (int j = i + 1; j < nr; j++)
    {
      const Substructure_Atom * rj = _root_atoms[j];
      assert (rj->is_matched());

      atom_number_t aj = rj->current_hold_atom()->atom_number();

      int d = m->bonds_between(ai, aj);
      if (! _distance_between_root_atoms.matches(d))
        return 0;
    }
  }

  return 1;
}

Link_Atom::Link_Atom()
{
  _a1 = INVALID_ATOM_NUMBER;
  _bt = INVALID_BOND_TYPE;
  _a2 = INVALID_ATOM_NUMBER;

  _bond_topology = -1;

  _mdl_atom_data = NULL;

  return;
}

Link_Atom::Link_Atom(const Link_Atom & rhs)
{
  _a1 = rhs._a1;
  _bt = rhs._bt;
  _a2 = rhs._a2;

  _bond_topology = rhs._bond_topology;

  _symbol = rhs._symbol;

  _e = NULL;

  if (NULL != _mdl_atom_data)
    delete _mdl_atom_data;

  _mdl_atom_data = rhs._mdl_atom_data;

  _d = rhs._d;

  assert(_d.ok());

  return;
}

int
Link_Atom::debug_print(std::ostream & os) const
{
  if (_a1 < 0)
  {
    os << "Link_Atom::debug_print:not initialised\n";
    return os.good();
  }

  os << "Link_Atom::debug_print:between matched atoms " << _a1 << " and " << _a2;
  if (_bond_topology >= 0)
    os << " bond topology " << _bond_topology;

  if (_symbol.length())
    os << ", symbol '" << _symbol << "'";

  os << '\n';
  
  return _d.debug_print(os);
}

Link_Atom_Current_State::Link_Atom_Current_State()
{
  _lhs = INVALID_ATOM_NUMBER;
  _rhs = INVALID_ATOM_NUMBER;
  _bt = INVALID_BOND_TYPE;

  return;
}

int
Link_Atom_Current_State::initialise(Link_Atom const* l)
{
  _lhs = l->a1();
  _rhs = l->a2();

  const bond_type_t b = l->btype();

  if (IS_SINGLE_BOND(b))
    _bt = SINGLE_BOND;
  else if (IS_DOUBLE_BOND(b))
    _bt = DOUBLE_BOND;
  else if (IS_TRIPLE_BOND(b))
    _bt = TRIPLE_BOND;
  else
    _bt = SINGLE_BOND;

  return 1;
}

//#define DEBUG_LINK_ATOMS_SATISFIED

int
Single_Substructure_Query::_link_atoms_satisfied(Query_Atoms_Matched & matched_atoms) const
{
#ifdef DEBUG_LINK_ATOMS_SATISFIED
  cerr << "Testing " << _link_atom.number_elements() << " link atoms\n";
#endif

  for (int i = 0; i < _link_atom.number_elements(); i++)
  {
    const Link_Atom * l = _link_atom[i];

    if (! _link_atom_satisfied(*l, matched_atoms))
      return 0;
  }

  return 1;
}

int
Single_Substructure_Query::_link_atom_satisfied(const Link_Atom & l,
                                                Query_Atoms_Matched & matched_atoms) const
{
#ifdef DEBUG_LINK_ATOMS_SATISFIED
  cerr << "Matched atoms";
  for (int i = 0; i < matched_atoms.number_elements(); i++)
  {
    const Substructure_Atom * a = matched_atoms[i];

    cerr << ' ' << a->current_hold_atom()->atom_number();
  }
  cerr << endl;
#endif

  Molecule * m = matched_atoms[0]->current_hold_atom()->m();

  const atom_number_t a1 = matched_atoms[l.a1()]->current_hold_atom()->atom_number();
  const atom_number_t a2 = matched_atoms[l.a2()]->current_hold_atom()->atom_number();

  if (m->fragment_membership(a1) != m->fragment_membership(a2))
  {
//  cerr << "Single_Substructure_Query::_link_atom_satisfied:atoms in different fragments\n";
    return 0;
  }

  int d = m->bonds_between(a1, a2) - 1;   // Apr 2004, change from bonds between to atoms between

#ifdef DEBUG_LINK_ATOMS_SATISFIED
  cerr << "Link atoms " << l.a1() << " and " << l.a2() << endl;

  cerr << "matched atoms " << a1 << " and " << a2 << " are " << d << " bonds apart, satisfies separation " << l.satisfies_separation(d) << endl;
#endif

  return l.satisfies_separation(d);
}

int
Link_Atom::swap_atoms(atom_number_t n1, atom_number_t n2)
{
  if (_a1 == n1)
  {
    _a1 = n2;
    return 1;
  }

  if (_a2 == n1)
  {
    _a2 = n2;
    return 1;
  }

  return 1;
}

int
Link_Atom::atom_is_being_removed(atom_number_t zatom)
{
  assert (zatom != _a1);
  assert (zatom != _a2);

  if (_a1 > zatom)
    _a1--;

  if (_a2 > zatom)
    _a2--;

  return 1;
}

int
Link_Atom::initialise_from_mdl_record(const const_IWSubstring & buffer,
                                       int matoms,
                                       atom_number_t & a)
{
  if (buffer.nwords() < 6)
  {
    cerr << "Link_Atom::initialise_from_mdl_record:must have at least 6 words\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);    // M
  buffer.nextword(token, i);    // LIN
  buffer.nextword(token, i);    // sqeuential number, not of interest

  buffer.nextword(token, i);   // atom number

  if (! token.numeric_value(a) || a < 1 || a > matoms)
  {
    cerr << "MDL_Molecule::_parse_link_record:invalid central atom\n";
    return 0;
  }

  a--;    // convert to our numbering

  buffer.nextword(token, i);  //   range

  int range;
  if (token.numeric_value(range) && range > 0)
  {
    _d.set_max(range);
  }
  else
  {
    _d.reset();

    if (_d.initialise(token))
    {
      for (int j = 0; j < _d.number_elements(); j++)    // convert from atom count to bond count
      {
        _d[j] = _d[j] + 1;
      }
    }
    else
    {
      cerr << "Link_Atom::initialise_from_mdl_record:invalid range '" << token << "'\n";
      return 0;
    }
  }

  assert (_d.ok());

  buffer.nextword(token, i);
  if (! token.numeric_value(_a1) || _a1 < 1)
  {
    cerr << "Link_Atom::initialise_from_mdl_record:invalid a1 '" << token << "'\n";
    return 0;
  }

  _a1--;

  buffer.nextword(token, i);
  if (! token.numeric_value(_a2) || _a2 < 1)
  {
    cerr << "Link_Atom::initialise_from_mdl_record:invalid a2 '" << token << "'\n";
    return 0;
  }

  _a2--;

//cerr << "Link_Atom:initialised, between " << _a1 << " and " << _a2 << endl;

  return 1;
}

int
Link_Atom::set_symbol(const const_IWSubstring & s)
{
  _symbol = s;

  _e = get_element_from_symbol_no_case_conversion(s);   // may return NULL, A, Q, [CD4]

  return 1;
}

int
Single_Substructure_Query::add_no_matched_atoms_between_initial_atom_numbers(int n1,
                                                                             int n2)
{
  Substructure_Atom * a1 = query_atom_with_initial_atom_number(n1);
  Substructure_Atom * a2 = query_atom_with_initial_atom_number(n2);

  if (NULL == a1 || NULL == a2)
  {
    cerr << "Single_Substructure_Query::add_no_matched_atoms_between_initial_atom_numbers:no query atoms found\n";
    cerr << "a1 = " << a1 << " or a2 = " << a2 << endl;
    return 0;
  }

  Bond * b = new Bond(a1->unique_id(), a2->unique_id(), SINGLE_BOND);

  _no_matched_atoms_between.add(b);

  return _no_matched_atoms_between.number_elements();
}

int
Link_Atom::create_next_variant(MDL_Molecule & m, Link_Atom_Current_State & lacs) const
{
  int maxdist;
  (void) _d.max(maxdist);

  if (lacs.number_placed() == maxdist)
    return 0;

  int natoms = m.natoms();

// If this is a clean element, just use it, otherwise transfer query specifications

  if (NULL != _e)              // atom was C, O, etc...
    m.add(_e);            
  else if (_symbol.length())              // atom was *, A, Q, [SD2], etc...
    m.add_atom_based_on_symbol(_symbol);
  else
    m.add_atom_based_on_symbol("*");    // how could this happen?

#ifdef DEBUG_CREATE_NEXT_VARIANT
  cerr << "Start with '" << m.smiles() << "'\n";
  cerr << " adding bond " << lacs.lhs() << " and " << natoms << endl;
#endif

// When nothing has been placed, there will be a gap between lhs and rhs.
// Otherwise we need to remove the existing bond to rhs and insert our new atom

  if (0 == lacs.number_placed())
    m.add_bond(lacs.lhs(), natoms, lacs.btype_for_molecule(), _bt, _bond_topology);
  else
  {
    m.remove_bond_between_atoms(lacs.last_atom_placed(), lacs.rhs());
    m.add_bond(lacs.last_atom_placed(), natoms, lacs.btype_for_molecule(), _bt, _bond_topology);
  }

  m.add_bond(natoms, lacs.rhs(), lacs.btype_for_molecule(), _bt, _bond_topology);

#ifdef DEBUG_CREATE_NEXT_VARIANT
  cerr << "create_next_variant created '" << m.smiles() << endl;
#endif

  lacs.add(natoms);

  if (! m.arrays_allocated())
    m.build(m);

  if (NULL != _mdl_atom_data)
    m.set_mdl_atom_data(natoms, _mdl_atom_data);

  return 1;
}

/*
  We are deleting the atom which was the link atom. We need to transfer
  any query specifications that may have been associated with it.
*/

int
Link_Atom::set_mdl_atom_data(const MDL_Atom_Data * a)
{
  if (NULL == _mdl_atom_data)
    _mdl_atom_data = new MDL_Atom_Data(*a);
  else
    (*_mdl_atom_data) = *a;

  return 1;
}
