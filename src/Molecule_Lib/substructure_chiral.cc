#include <stdlib.h>

#include "substructure.h"
#include "target.h"


Substructure_Chiral_Centre::Substructure_Chiral_Centre()
{
  _numeric = NULL;

  _centre = NULL;
  _top_front = NULL;
  _top_back = NULL;
  _left_down = NULL;
  _right_down = NULL;

  _number_explicit_connections = 0;

  return;
}

Substructure_Chiral_Centre::~Substructure_Chiral_Centre()
{
  if (NULL != _numeric)
    delete _numeric;

  return;
}
int
Substructure_Chiral_Centre::debug_print (std::ostream & os) const
{
  os << "Substructure_Chiral_Centre::debug_print:has " << _number_explicit_connections << " explicit connections\n";
  if (NULL != _numeric)
  {
    os << "Numeric is object allocated\n";
    _numeric->debug_print (cerr);
  }

  if (NULL != _centre)
    os << " centre " << _centre->unique_id() << endl;
  if (NULL != _top_front)
    os << " top_front " << _top_front->unique_id() << endl;
  if (NULL != _top_back)
    os << " top_back " << _top_back->unique_id() << endl;
  if (NULL != _left_down)
    os << " left_down " << _left_down->unique_id() << endl;
  if (NULL != _right_down)
    os << " right_down " << _right_down->unique_id() << endl;

  return os.good ();
}

void
Substructure_Chiral_Centre::set_top_front(const Substructure_Atom * s)
{
  _top_front = s;
  if (NULL != _top_front)
    _number_explicit_connections++;

  return;
}

void
Substructure_Chiral_Centre::set_top_back(const Substructure_Atom * s)
{
  _top_back = s;
  if (NULL != _top_back)
    _number_explicit_connections++;

  return;
}

void
Substructure_Chiral_Centre::set_left_down(const Substructure_Atom * s)
{
  _left_down = s;
  if (NULL != _left_down)
    _number_explicit_connections++;

  return;
}

void
Substructure_Chiral_Centre::set_right_down(const Substructure_Atom * s)
{
  _right_down = s;
  if (NULL != _right_down)
    _number_explicit_connections++;

  return;
}

int
Single_Substructure_Query::_build_chirality_specification_from_msi_attribute (const IWString & s)
{
  if (5 != s.nwords())
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:must be 5 tokens '" << s << "'\n";
    return 0;
  }

//cerr << "Building from '" << s << "'\n";

  Substructure_Chiral_Centre * c = new Substructure_Chiral_Centre;

  resizable_array<int> numbers_encountered;   // make sure no duplicates

  int i = 0;
  const_IWSubstring token;

  s.nextword(token, i);

  int uid;
  if (! token.numeric_value(uid) || uid < 0)
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:invalid centre '" << s << "'\n";
    return 0;
  }

  const Substructure_Atom * a = query_atom_with_initial_atom_number(uid);
  if (NULL == a)
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:no centre atom '" << s << "'\n";
    return 0;
  }

  c->set_centre(a);

  numbers_encountered.add(uid);

  s.nextword(token, i);

  if ("H" == token)
    numbers_encountered.add_if_not_already_present(-9);
  else if (! token.numeric_value(uid) || uid < 0)
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:invalid TF '" << s << "'\n";
    return 0;
  }
  else if (NULL == (a = query_atom_with_initial_atom_number(uid)))
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:no TF atom '" << s << "'\n";
    return 0;
  }
  else
  {
    c->set_top_front(a);
    numbers_encountered.add_if_not_already_present(uid);
  }

  s.nextword(token, i);

  if ("H" == token)
    numbers_encountered.add_if_not_already_present(-9);
  else if (! token.numeric_value(uid) || uid < 0)
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:invalid TB '" << s << "'\n";
    return 0;
  }
  else if (NULL == (a = query_atom_with_initial_atom_number(uid)))
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:no TB atom '" << s << "'\n";
    return 0;
  }
  else
  {
    c->set_top_back(a);
    numbers_encountered.add_if_not_already_present(uid);
  }

  s.nextword(token, i);

  if ("H" == token)
    numbers_encountered.add_if_not_already_present(-9);
  else if (! token.numeric_value(uid) || uid < 0)
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:invalid LD '" << s << "'\n";
    return 0;
  }
  else if (NULL == (a = query_atom_with_initial_atom_number(uid)))
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:no LD atom '" << s << "'\n";
    return 0;
  }
  else
  {
    c->set_left_down(a);
    numbers_encountered.add(uid);
  }

  s.nextword(token, i);

  if ("H" == token)
    numbers_encountered.add_if_not_already_present(-9);
  else if (! token.numeric_value(uid) || uid < 0)
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:invalid RD '" << s << "'\n";
    return 0;
  }
  else if (NULL == (a = query_atom_with_initial_atom_number(uid)))
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:no RD atom '" << s << "'\n";
    return 0;
  }
  else
  {
    c->set_right_down(a);
    numbers_encountered.add_if_not_already_present(uid);
  }

  if (5 != numbers_encountered.number_elements())
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:duplicate atom numbers '" << s << "'\n";
    return 0;
  }

  _chirality.add(c);

  return 1;
}

static int
get_matched_atom_number (const Substructure_Atom * a,
                         atom_number_t & matched_atom,
                         atom_number_t what_to_add_if_not_matched)
{
  if (NULL == a)
  {
    matched_atom = what_to_add_if_not_matched;
    return 0;
  }

  matched_atom = a->current_hold_atom()->atom_number();
  return 1;
}

//#define DEBUG_SUBSTURE_CHIRAL_CENTRE_IS_MATCHED

/*
  We need to figure out of two sets of 3 numbers are in the same direction or not.
  Some of the matches may be negative, which corresponds to implicit hydrogens
  We may need to rotate them.
*/

static int
check_match(atom_number_t north1,
            atom_number_t sw1,
            atom_number_t se1,
            atom_number_t north2,
            atom_number_t sw2,
            atom_number_t se2)
{
#ifdef DEBUG_SUBSTURE_CHIRAL_CENTRE_IS_MATCHED
  cerr << "check_match:checking: north1 " << north1 << " sw1 " << sw1 << " se1 " << se1 << endl;
  cerr << "        compare with: north2 " << north2 << " sw2 " << sw2 << " se2 " << se2 << endl;
#endif

  if (north1 == north2)
      return sw1 == sw2;

  if (north1 == se2)
    return sw1 == north2;

  if (north1 == sw2)
    return se1 == north2;

#ifdef DEBUG_SUBSTURE_CHIRAL_CENTRE_IS_MATCHED
  cerr << "check_match:no match: north1 " << north1 << " sw1 " << sw1 << " se1 " << se1 << endl;
  cerr << "        compare with: north2 " << north2 << " sw2 " << sw2 << " se2 " << se2 << endl;
#endif

  return 0;
}

int
Single_Substructure_Query::_chiral_atoms_matched (Query_Atoms_Matched & matched_atoms,
                                    Molecule_to_Match & target_molecule) const
{
  const Molecule * m = target_molecule.molecule();

  if (NULL == m)
  {
//  cerr << "Single_Substructure_Query::_chiral_atoms_matched:ignoring chirality\n";
    return 1;
  }

//cerr << "Checking " << _chirality.number_elements() << " chirality specifications\n";

  for (int i = 0; i < _chirality.number_elements(); i++)
  {
    if (! _chirality[i]->is_matched (target_molecule.molecule()))
      return 0;
  }

  return 1;         // all chirality specifications OK
}


//#define DEBUG_SUBSTURE_CHIRAL_CENTRE_IS_MATCHED
#ifdef DEBUG_SUBSTURE_CHIRAL_CENTRE_IS_MATCHED
static int
write_null_or_current_hold_atom (const Substructure_Atom * a)
{
  if (NULL == a)
    return -1;
  else
    return a->current_hold_atom()->atom_number();

  return -2;
}

#endif


int
Substructure_Chiral_Centre::is_matched (const Molecule * m) const
{
  assert (NULL != _centre);

  if (NULL == _centre->current_hold_atom())    // hard to imagine
    return 0;

  atom_number_t centre_atom = _centre->current_hold_atom()->atom_number();

#ifdef DEBUG_SUBSTURE_CHIRAL_CENTRE_IS_MATCHED
  cerr << "Centre matched with " << centre_atom << endl;
#endif

  const Chiral_Centre * c = m->chiral_centre_at_atom(centre_atom);
  if (NULL == c)
    return 0;

#ifdef DEBUG_SUBSTURE_CHIRAL_CENTRE_IS_MATCHED
  cerr << "Got chiral centre on matched atom\n";
#endif

  int matched_connections = 0;

#ifdef DEBUG_SUBSTURE_CHIRAL_CENTRE_IS_MATCHED
  cerr << " tf " << write_null_or_current_hold_atom(_top_front) << " tb " << write_null_or_current_hold_atom(_top_back) << " ld " << write_null_or_current_hold_atom(_left_down) << " rd " << write_null_or_current_hold_atom(_right_down) << endl;
#endif

  atom_number_t what_to_add_if_not_matched;

  if (4 == m->ncon(centre_atom))
    what_to_add_if_not_matched = INVALID_ATOM_NUMBER;
  else if (const_cast<Molecule *>(m)->hcount(centre_atom))
    what_to_add_if_not_matched = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else
    what_to_add_if_not_matched = CHIRAL_CONNECTION_IS_LONE_PAIR;

//cerr << "Extra is " << what_to_add_if_not_matched << endl;

  atom_number_t tf;
  if (get_matched_atom_number(_top_front, tf, what_to_add_if_not_matched))
    matched_connections++;

  atom_number_t tb;
  if (get_matched_atom_number(_top_back, tb, what_to_add_if_not_matched))
    matched_connections++;

  atom_number_t ld;
  if (get_matched_atom_number(_left_down, ld, what_to_add_if_not_matched))
    matched_connections++;

  atom_number_t rd;
  if (get_matched_atom_number(_right_down, rd, what_to_add_if_not_matched))
    matched_connections++;

#ifdef DEBUG_SUBSTURE_CHIRAL_CENTRE_IS_MATCHED
  cerr << "Found " << matched_connections << " matched connections, compare " << _number_explicit_connections << endl;
  cerr << "tf " << tf << " tb " << tb << " ld " << ld << " rd " << rd << endl;
  c->debug_print(cerr);
  cerr << _number_explicit_connections << " explicit connections, target " << c->number_connections_specified() << endl;
#endif

  if (matched_connections != _number_explicit_connections)
    return 0;

// If we only have 3 atoms matched, the first thing is to identify the
// unmatched atom

  if (3 == _number_explicit_connections && 4 == c->number_atoms_specified())
  {
    if (c->top_front() == tf || c->top_front() == tb || c->top_front() == ld || c->top_front() == rd)
      ;
    else
    {
      if (tf < 0)
        return check_match(tb, ld, rd, c->top_back(), c->left_down(), c->right_down());
      else if (tb < 0)
        return check_match(tf, rd, ld, c->top_back(), c->left_down(), c->right_down());
      else if (ld < 0)
        return check_match(tf, tb, rd, c->top_back(), c->left_down(), c->right_down());
      else if (rd < 0)
        return check_match(tf, ld, tb, c->top_back(), c->left_down(), c->right_down());
    }

    if (c->top_back() == tf || c->top_back() == tb || c->top_back() == ld || c->top_back() == rd)
      ;
    else
    {
      if (tf < 0)
        return check_match(tb, ld, rd, c->top_front(), c->right_down(), c->left_down());
      else if (tb < 0)
        return check_match(tf, rd, ld, c->top_front(), c->right_down(), c->left_down());
      else if (ld < 0)
        return check_match(tf, tb, rd, c->top_front(), c->right_down(), c->left_down());
      else if (rd < 0)
        return check_match(tf, ld, tb, c->top_front(), c->right_down(), c->left_down());
    }

    if (c->left_down() == tf || c->left_down() == tb || c->left_down() == ld || c->left_down() == rd)
      ;
    else
    {
      if (tf < 0)
        return check_match(tb, ld, rd, c->top_front(), c->top_back(), c->right_down());
      else if (tb < 0)
        return check_match(tf, rd, ld, c->top_front(), c->top_back(), c->right_down());
      else if (ld < 0)
        return check_match(tf, tb, rd, c->top_front(), c->top_back(), c->right_down());
      else if (rd < 0)
        return check_match(tf, ld, tb, c->top_front(), c->top_back(), c->right_down());
    }

    if (c->right_down() == tf || c->right_down() == tb || c->right_down() == ld || c->right_down() == rd)
      ;
    else
    {
      if (tf < 0)
        return check_match(tb, ld, rd, c->top_front(), c->left_down(), c->top_back());
      else if (tb < 0)
        return check_match(tf, rd, ld, c->top_front(), c->left_down(), c->top_back());
      else if (ld < 0)
        return check_match(tf, tb, rd, c->top_front(), c->left_down(), c->top_back());
      else if (rd < 0)
        return check_match(tf, ld, tb, c->top_front(), c->left_down(), c->top_back());
    }

    cerr << "Substructure_Chiral_Centre:should not come to here\n";
    return 0;
  }

#ifdef DEBUG_SUBSTURE_CHIRAL_CENTRE_IS_MATCHED
  cerr << "LINE " << __LINE__ << " tb " << tb << " tf " << tf << endl;
#endif

// At this stage, we have assembled all the matched atoms. Now just do
// the comparison. Make sure we pass the arguments in anti-clockwise order

  if (tb < 0)
  {
    if (tf == c->top_front())
      return check_match(tb, ld, rd, c->top_back(), c->left_down(), c->right_down());
    else if (tf == c->top_back())
      return check_match(tb, ld, rd, c->top_front(), c->right_down(), c->left_down());
    else if (tf == c->left_down())
      return check_match(tb, ld, rd, c->top_front(), c->top_back(), c->right_down());
    else if (tf == c->right_down())
      return check_match(tb, ld, rd, c->top_front(), c->left_down(), c->top_back());
  }
  else if (tf >= 0)
  {
    if (tf == c->top_front())
      return check_match(tb, ld, rd, c->top_back(), c->left_down(), c->right_down());
    else if (tf == c->top_back())
      return check_match(tb, ld, rd, c->top_front(), c->right_down(), c->left_down());
    else if (tf == c->left_down())
      return check_match(tb, ld, rd, c->top_front(), c->top_back(), c->right_down());
    else if (tf == c->right_down())
      return check_match(tb, ld, rd, c->top_front(), c->left_down(), c->top_back());
  }
  else if (tb >= 0)
  {
    if (tb == c->top_front())
      return check_match(tf, rd, ld, c->top_back(), c->left_down(), c->right_down());
    else if (tb == c->top_back())
      return check_match(tf, rd, ld, c->top_front(), c->right_down(), c->left_down());
    else if (tb == c->left_down())
      return check_match(tf, rd, ld, c->top_front(), c->top_back(), c->right_down());
    else if (tb == c->right_down())
      return check_match(tf, rd, ld, c->top_front(), c->left_down(), c->top_back());
  }
  else
  {
    cerr << "Substructure_Chiral_Centre::is_matched:huh, what's matched\n";
    return 0;
  }

  cerr << "Should not come to here\n";

  return 0;
}


/*int
Substructure_Chiral_Centre::construct_from_msi_attribute (const msi_attribute & att)
{
  assert (NULL == _numeric);

  const IWString & s = att.stringval();

  int nw = s.nwords();

  if (nw < 4)
  {
    cerr << "Substructure_Chiral_Centre::construct_from_msi_attribute:invalid chirality '" << s << "'\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  s.nextword(token, i);

  atom_number_t zatom;
  if (! token.numeric_value(zatom) || zatom < 0)
  {
    cerr << "Substructure_Chiral_Centre::construct_from_msi_attribute:invalid centre '" << s << "'\n";
    return 0;
  }

  _numeric = new Chiral_Centre (zatom);

  _number_explicit_connections = 0;

  int nc = 0;    // don't set _number_explicit_connections here, it will be set when pointers are assigned

  s.nextword(token, i);

  if ('H' == token)
    _numeric->set_top_front(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  else if (! token.numeric_value(zatom) || zatom < 0)
  {
    cerr << "Substructure_Chiral_Centre::construct_from_msi_attribute:invalid TF '" << s << "'\n";
    return 0;
  }
  else
  {
    _numeric->set_top_front(zatom);
    nc++;
  }

  s.nextword(token, i);

  if ('H' == token)
    _numeric->set_top_back(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  else if (! token.numeric_value(zatom) || zatom < 0)
  {
    cerr << "Substructure_Chiral_Centre::construct_from_msi_attribute:invalid TB '" << s << "'\n";
    return 0;
  }
  else
  {
    _numeric->set_top_back(zatom);
    nc++;
  }

  s.nextword(token, i);

  if ('H' == token)
    _numeric->set_left_down(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  else if (! token.numeric_value(zatom) || zatom < 0)
  {
    cerr << "Substructure_Chiral_Centre::construct_from_msi_attribute:invalid LD '" << s << "'\n";
    return 0;
  }
  else
  {
    _numeric->set_left_down(zatom);
    nc++;
  }

  s.nextword(token, i);

  if ('H' == token)
    _numeric->set_right_down(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  else if (! token.numeric_value(zatom) || zatom < 0)
  {
    cerr << "Substructure_Chiral_Centre::construct_from_msi_attribute:invalid RD '" << s << "'\n";
    return 0;
  }
  else
  {
    _numeric->set_right_down(zatom);
    nc++;
  }


  if (nc < 3)
  {
    cerr << "Substructure_Chiral_Centre::construct_from_msi_attribute:not enough connections '" << s << "'\n";
    return 0;
  }

  return 1;
}*/

int
Substructure_Chiral_Centre::write_msi (std::ostream & os, 
                                     const IWString & ind,
                                     const char * attribute_name) const
{
  if (NULL == _centre)
  {
    cerr << "Substructure_Chiral_Centre::write_msi:centre not initialised\n";
    return 0;
  }

  os << ind << "  (A C " << attribute_name << " \"";

  os << _centre->unique_id () << ' ';

  if (NULL == _top_front)
    os << "H ";
  else
    os << _top_front->unique_id() << ' ';

  if (NULL == _top_back)
    os << "H ";
  else
    os << _top_back->unique_id() << ' ';

  if (NULL == _left_down)
    os << "H ";
  else
    os << _left_down->unique_id() << ' ';

  if (NULL == _right_down)
    os << "H";
  else
    os << _right_down->unique_id();

  os << "\")\n";

  return os.good();
}

int
Single_Substructure_Query::_add_chiral_centre (const Molecule & m,
                                               const Chiral_Centre & c)
{
  Substructure_Chiral_Centre * tmp = new Substructure_Chiral_Centre();

  atom_number_t zatom = c.a();

  const Substructure_Atom * a = query_atom_with_initial_atom_number(zatom);
  if (NULL == a)
  {
    cerr << "Single_Substructure_Query::_add_chiral_centre:no centre atom\n";
    return 0;
  }

  tmp->set_centre(a);

  zatom = c.top_front();

  if (zatom < 0)
    ;
  else if (NULL == (a = query_atom_with_initial_atom_number(zatom)))
  {
    cerr << "Single_Substructure_Query::_add_chiral_centre:no TF atom " << zatom << endl;
    return 0;
  }
  else
    tmp->set_top_front(a);

  zatom = c.top_back();

  if (zatom < 0)
    ;
  else if (NULL == (a = query_atom_with_initial_atom_number(zatom)))
  {
    cerr << "Single_Substructure_Query::_add_chiral_centre:no TB atom " << zatom << endl;
    return 0;
  }
  else
    tmp->set_top_back(a);

  zatom = c.left_down();

  if (zatom < 0)
    ;
  else if (NULL == (a = query_atom_with_initial_atom_number(zatom)))
  {
    cerr << "Single_Substructure_Query::_add_chiral_centre:no LD atom " << zatom << endl;
    return 0;
  }
  else
    tmp->set_left_down(a);

  zatom = c.right_down();

  if (zatom < 0)
    ;
  else if (NULL == (a = query_atom_with_initial_atom_number(zatom)))
  {
    cerr << "Single_Substructure_Query::_add_chiral_centre:no RD atom " << zatom << endl;
    return 0;
  }
  else
    tmp->set_right_down(a);

  _chirality.add(tmp);

  return 1;
}

/*
  Note that we are not checking for one of the connections being a lone
  pair. Too rare...
*/

/*int
Substructure_Chiral_Centre::extract_connections(const Chiral_Centre & c)
{
  assert (NULL == _numeric);

  _numeric = new Chiral_Centre(c.a());

  _numeric->set_top_front(c.top_front());
  _numeric->set_top_back(c.top_back());
  _numeric->set_left_down(c.left_down());
  _numeric->set_right_down(c.right_down());

  return 1;
}*/

int
Substructure_Atom::first_chirality_value_in_any_component() const
{
  if (_chirality > 0)
    return _chirality;

  for (int i = 0; i < _components.number_elements(); i++)
  {
    if (_components[i]->chirality() > 0)
      return _components[i]->chirality();
  }

  return -1;
}

int
Substructure_Chiral_Centre::invert()
{
  assert (NULL != _centre);

  const Substructure_Atom * tmp = _top_front;
  _top_front = _top_back;
  _top_back = tmp;

  return 1;
}

int
Single_Substructure_Query::_build_chirality_specifications_from_atom_chiral_info()
{
  for (int i = 0; i < _root_atoms.number_elements(); i++)
  {
    if (! _build_chirality_specifications_from_atom_chiral_info(_root_atoms[i]))
    {
      cerr << "Single_Substructure_Query::_build_chirality_specifications_from_atom_chiral_info:invalid chirality\n";
      return 0;
    }
  }

  return 1;
}

int
Single_Substructure_Query::_build_chirality_specifications_from_atom_chiral_info(const Substructure_Atom * r)
{
//cerr << "_build_chirality_specifications_from_atom_chiral_info:chirality " << r->first_chirality_value_in_any_component() << endl;
  if (r->first_chirality_value_in_any_component() > 0)
  {
    if (! _build_chirality_specification_from_atom_chiral_info(r))
    {
      cerr << "Single_Substructure_Query::_build_chirality_specifications_from_atom_chiral_info:invalid chirality, atom " << r->unique_id() << "\n";
      return 0;
    }
  }

  for (int i = 0; i < r->number_children(); i++)
  {
    Substructure_Atom * ci = r->child(i);
    if (! _build_chirality_specifications_from_atom_chiral_info(ci))
      return 0;
  }

  return 1;
}

int
Single_Substructure_Query::_build_chirality_specification_from_atom_chiral_info(const Substructure_Atom * a)
{
  assert (a->first_chirality_value_in_any_component() > 0);

//cerr << "Building chirality information for atom " << a->initial_atom_number() << endl;

  Substructure_Chiral_Centre * tmp = new Substructure_Chiral_Centre;

  tmp->set_centre(a);

  const Substructure_Atom * p = a->parent();

  if (NULL != p)
    tmp->set_top_front(p);

  int nc = a->number_children() + a->nbonds() - 1;   // children + ring closure bonds

  if (nc < 1 || nc > 4)
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_atom_chiral_info:not enough connections, nc = " << nc << endl;
    cerr << "Atom " << a->initial_atom_number() << " children " << a->number_children() << " bonds " << a->nbonds() << " parent? " << (NULL != p) << " ring closure bonds " << a->number_ring_closure_bonds() << endl;
    return 0;
  }

  int ndx = 0;

  for (int i = 1; i < a->number_ring_closure_bonds(); i++)   // first one is bond to parent
  {
    const Substructure_Bond * b = a->ring_closure_bond(i);

    const Substructure_Atom * a = b->a();

    if (0 == ndx)
      tmp->set_top_back(a);
    else if (1 == ndx)
      tmp->set_left_down(a);
    else if (2 == ndx)
      tmp->set_right_down(a);

    ndx++;

#ifdef DEBUG_BUILD_CHIRALITY_SPECIFICATION
    cerr << "processing ring closures, ndx " << ndx << endl;
    tmp->debug_print(cerr);
#endif
  }

  for (int i = 0; i < a->number_children(); i++)
  {
    const Substructure_Atom * c = a->child(i);
    if (0 == ndx)
      tmp->set_top_back(c);
    else if (1 == ndx)
      tmp->set_left_down(c);
    else if (2 == ndx)
      tmp->set_right_down(c);

    ndx++;

#ifdef DEBUG_BUILD_CHIRALITY_SPECIFICATION
    cerr << "processing children, atom " << c->initial_atom_number() << ", ndx " << ndx << endl;
    tmp->debug_print(cerr);
#endif
  }

#ifdef DEBUG_BUILD_CHIRALITY_SPECIFICATION
  cerr << "At end of processing ndx " << ndx << endl;
  tmp->debug_print (cerr);
#endif

// We need to check to see if any other atoms have a ring closure bond to this atom, thereby
// completing the chiral centre

  int ring_closure_bonds_present = 0;

  if (3 == ndx)   // done
    ;
  else
  {
    for (int i = 0; i < a->number_children(); i++)
    {
      const Substructure_Atom * c = a->child(i);

      if (0 == c->complete_chiral_centre_via_ring_closure(*tmp, ndx))
        continue;

      ring_closure_bonds_present++;

      if (3 == ndx)
        break;
    }
  }

#ifdef DEBUG_BUILD_CHIRALITY_SPECIFICATION
  cerr << "After looking for ring closures in children\n";
  tmp->debug_print(cerr);

  cerr << "Check inversion " << a->first_chirality_value_in_any_component() << endl;
#endif

/*
  No, this is not right. We do not deal with the case of two ring closures coming
  into the atom. Too messy! This does not get used much...
*/

  if (ring_closure_bonds_present)
  {
    if (1 == a->first_chirality_value_in_any_component())
      tmp->invert();
  }
  else if (2 == a->first_chirality_value_in_any_component())    // atom was @@
    tmp->invert();

  _chirality.add(tmp);

//cerr << "Now have " << _chirality.number_elements() << " chirality specifications\n";

  return 1;
}

int
Substructure_Atom::complete_chiral_centre_via_ring_closure (Substructure_Chiral_Centre & c,
                                           int & ndx) const
{
  int rc = 0;

  for (int i = 1; i < _bonds.number_elements(); i++)    // omit bond to parent
  {
    const Substructure_Bond * b = _bonds[i];

    const Substructure_Atom * o = b->a();

    if (o->initial_atom_number() != c.centre()->initial_atom_number())
      continue;

    if (1 == ndx)
      c.set_left_down(this);
    else if (2 == ndx)
      c.set_right_down(this);

    rc++;

    ndx++;
    if (ndx > 2)
      return 1;
  }

  for (int i = 0; i < _children.number_elements(); i++)
  {
    if (! _children[i]->complete_chiral_centre_via_ring_closure(c, ndx))
      continue;

    rc++;
    if (ndx > 2)
      return rc;
  }

  return rc;
}

/*
  The chirality objects just have their atom numbers set, no pointers yet
*/

/*int 
Single_Substructure_Query::_initialise_chirality_info()
{
//cerr << "Initialising " << _chirality.number_elements() << " chirality specifications\n";

  for (int i = 0; i < _chirality.number_elements(); i++)
  {
    if (! _initialise_chirality_info(_chirality[i]))
    {
      cerr << "Single_Substructure_Query::initialise_chirality_info:bad chirality info\n";
      return 0;
    }
  }

  return 1;
}*/

/*
  The Substructure_Chiral_Centre object knows the query atom numbers in
  its _numeric attribute. We need to convert those to pointers to
  Substructure_Atoms
*/

/*int
Single_Substructure_Query::_initialise_chirality_info(Substructure_Chiral_Centre * c)
{
  const Chiral_Centre * numeric = c->numeric_chiral_centre();

  atom_number_t zatom = numeric->a();

  const Substructure_Atom * centre_query_atom = query_atom_with_initial_atom_number (zatom);
  if (NULL == centre_query_atom)
  {
    cerr << "Single_Substructure_Query::_initialise_chirality_info:no atom number " << zatom << "'\n";
    return 0;
  }

  c->set_centre(centre_query_atom);

#define DEBUG_INITIALISE_CHIRALITY_INFO
#ifdef DEBUG_INITIALISE_CHIRALITY_INFO
  cerr << "Initialising from ";
  numeric->debug_print(cerr);
#endif

  const Substructure_Atom * a;    // scope here for efficiency

  zatom = numeric->top_front();
  if (zatom >= 0)
  {
    a = query_atom_with_initial_atom_number (zatom);
    if (NULL == a)
    {
      cerr << "Single_Substructure_Query::_initialise_chirality_info:no atom number " << zatom << "'\n";
      return 0;
    }
    c->set_top_front(a);
  }

  zatom = numeric->top_back();
  if (zatom >= 0)
  {
    a = query_atom_with_initial_atom_number (zatom);
    if (NULL == a)
    {
      cerr << "Single_Substructure_Query::_initialise_chirality_info:no atom number " << zatom << "'\n";
      return 0;
    }
    c->set_top_back(a);
  }

  zatom = numeric->left_down();
  if (zatom >= 0)
  {
    a = query_atom_with_initial_atom_number (zatom);
    if (NULL == a)
    {
      cerr << "Single_Substructure_Query::_initialise_chirality_info:no atom number " << zatom << "'\n";
      return 0;
    }
    c->set_left_down(a);
  }

  zatom = numeric->right_down();
  if (zatom >= 0)
  {
    a = query_atom_with_initial_atom_number (zatom);
    if (NULL == a)
    {
      cerr << "Single_Substructure_Query::_initialise_chirality_info:no atom number " << zatom << "'\n";
      return 0;
    }
    c->set_right_down(a);
  }

  cerr << "Chirality value for centre atom " << centre_query_atom->first_chirality_value_in_any_component() << endl;
  if (1 == centre_query_atom->first_chirality_value_in_any_component())
    c->invert();

  return 1;
}*/
