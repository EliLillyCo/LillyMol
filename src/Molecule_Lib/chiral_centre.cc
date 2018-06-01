#include <stdlib.h>

#include <assert.h>
#include <limits>

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "molecule.h"
#include "mdl.h"
#include "chiral_centre.h"
#include "misc2.h"
#include "iwminmax.h"

static int _automatically_add_implicit_hydrogen_to_incomplete_chiral_centre = 0;

void
set_automatically_add_implicit_hydrogen_to_incomplete_chiral_centre (int s)
{
  _automatically_add_implicit_hydrogen_to_incomplete_chiral_centre = s;

  return;
}

void
Chiral_Centre::_default_values()
{
  _a = INVALID_ATOM_NUMBER;

  _top_front = _top_back = _left_down = _right_down = INVALID_ATOM_NUMBER;

  return;
}

Chiral_Centre::Chiral_Centre (atom_number_t a)
{
  _default_values();

  _a = a;

  _chirality_known = 0;

  return;
}

Chiral_Centre::Chiral_Centre(const Chiral_Centre & rhs)
{
  _a = rhs._a;
  _top_front = rhs._top_front;
  _top_back = rhs._top_back;
  _left_down = rhs._left_down;
  _right_down = rhs._right_down;

  return;
}

Chiral_Centre::~Chiral_Centre()
{
  assert(ok());

  _default_values();
}

static void
print_atom (std::ostream & os,
            const char * s,
            atom_number_t a)
{
  os << s;
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == a)
  {
    os << "implicit hydrogen";
    return;
  }

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == a)
  {
    os << "lone pair";
    return;
  }

  os << a;

  return;
}

int
Chiral_Centre::debug_print (std::ostream & os) const
{
  os << "Chiral centre at atom " << _a << endl;
  if (! ok())
    os << "Warning, OK function fails\n";

  print_atom(os, "    Top Front = ", _top_front);
  os << endl;

  print_atom(os, "        Top Back  = ", _top_back);
  os << endl;

  print_atom(os, "Left Bottom = ", _left_down);
  os << "    ";
  print_atom(os, "Right Bottom = ", _right_down);
  os << endl;

  return os.good();
}

void
Molecule::_print_atom_and_type (std::ostream & os,
                              const char * s,
                              atom_number_t centre,
                              atom_number_t a) const
{
  print_atom(os, s, a);

  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == a)
    ;
  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == a)
    ;
  else if (INVALID_ATOM_NUMBER == a)
    ;
  else if (a >= 0 && a < _number_elements)
  {
    os << " (" << atomic_symbol(a) << ")";
    if (INVALID_ATOM_NUMBER != centre && ! are_bonded(centre, a))
      os << " NOT BONDED!!";
  }
  else
    os << " Invalid atom number (" << a << ')';

  return;
}

int
Molecule::print_chiral_centre_details (const Chiral_Centre * c, std::ostream & os) const
{
  atom_number_t centre = c->a();

  _print_atom_and_type(os, "Chiral centre at atom ", INVALID_ATOM_NUMBER, centre);
  os << ' ' << ncon(c->a()) << " connections";
  if (c->chirality_known())
    os << " chirality known\n";
  else
    os << " chirality unknown\n";

  if (! ok())
    os << "Warning, OK function fails\n";

  _print_atom_and_type(os, "Top Front = ", centre, c->top_front());
  os << endl;

  _print_atom_and_type(os, "Top Back = ", centre, c->top_back());
  os << endl;

  _print_atom_and_type(os, "Left Down = ", centre, c->left_down());
  os << " ";
  _print_atom_and_type(os, "Right Down = ", centre, c->right_down());
  os << endl;

  return os.good();
}

/*
  Since these things get built gradually, just about the only thing
  we can say for sure is that the centre atom must be defined.

  We also add a test for at most one of the connected atoms being an
  implicit hydrogen
*/

int
Chiral_Centre::ok() const
{
  if (INVALID_ATOM_NUMBER == _a)
    return 0;

// Top_front must be set before Top_back

  if (INVALID_ATOM_NUMBER != _top_back &&
      INVALID_ATOM_NUMBER == _top_front)
    return 0;

  if (implicit_hydrogen_count() > 1)
    return 0;

  if (lone_pair_count() > 1)
    return 0;

  return 1;
}

int
Molecule::valid_chiral_centre (const Chiral_Centre * c) const
{
  if (! ok())
    return 0;

  atom_number_t centre = c->a();

  if (centre < 0 || centre >= _number_elements)
    return 0;

  const Atom * a = _things[centre];

  int acon = a->ncon();

  if (4 == acon)
  {
    if (c->implicit_hydrogen_count() || c->lone_pair_count())
    return 0;
  }
  else if (3 == acon)
  {
//  cerr << "IH " << c->implicit_hydrogen_count() << " LP " << c->lone_pair_count() << endl;
    if (1 != (c->implicit_hydrogen_count() + c->lone_pair_count()))
      return 0;
  }
  else       // something with fewer than 3 connections cannot be a chiral centre - well, it could have a lone pair and an implicit Hydrogen...
    return 0;

  atom_number_t top_front = c->top_front();

  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == top_front)
    ;
  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == top_front)
    ;
  else if (INVALID_ATOM_NUMBER == top_front)
    ;
  else if (! a->is_bonded_to (top_front))
    return 0;

  atom_number_t top_back = c->top_back();

  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == top_back)
    ;
  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == top_back)
    ;
  else if (INVALID_ATOM_NUMBER == top_back)
    ;
  else if (! a->is_bonded_to(top_back))
    return 0;

  atom_number_t left_down = c->left_down();

  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == left_down)
    ;
  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == left_down)
    ;
  else if (INVALID_ATOM_NUMBER == left_down)
    ;
  else if (! a->is_bonded_to(left_down))
    return 0;

  atom_number_t right_down = c->right_down();
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == right_down)
    ;
  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == right_down)
    ;
  else if (INVALID_ATOM_NUMBER == right_down)
    ;
  else if (! a->is_bonded_to(right_down))
    return 0;

  return 1;
}

/*
  Is a member of a chiral_centre object unset?
*/

static int
atom_not_set (int zatom)
{
  if (zatom >= 0)    // a valid atom number
    return 0;

  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == zatom)
    return 0;

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == zatom)
    return 0;

// Anything else is unspecified

  if (INVALID_ATOM_NUMBER == zatom)
    return 1;

  return 1;    // anything else < 0
}

/*
  Note that the assert(ok()) is at the end of the function.
*/

int
Chiral_Centre::complete() const
{
//debug_print(cerr);
  int number_not_set = 0;

  number_not_set += atom_not_set(_top_front);
  number_not_set += atom_not_set(_top_back);
  number_not_set += atom_not_set(_left_down);
  number_not_set += atom_not_set(_right_down);

//cerr << number_not_set << " not set\n";
  if (number_not_set)
    return 0;

  assert(ok());

  return 1;
}

int
Chiral_Centre::number_connections_specified () const
{
  int rc = ! atom_not_set(_top_front) +
           ! atom_not_set(_top_back)  +
           ! atom_not_set(_left_down) +
           ! atom_not_set(_right_down);

  return rc;
}

int
Chiral_Centre::number_atoms_specified () const
{
  return (_top_front >= 0) +
           (_top_back >= 0) +
           (_left_down >= 0) +
           (_right_down >= 0);
}

int
Chiral_Centre::involves (atom_number_t at) const
{
  assert(INVALID_ATOM_NUMBER != at);
  assert(ok());

  if (at == _a)
    return 1;

  if (_top_front == at)
    return 1;
  if (_top_back == at)
    return 1;
  if (_left_down == at)
    return 1;
  if (_right_down == at)
    return 1;

  return 0;
}

/*
  Does a chiral centre involve two different atoms?
*/

int
Chiral_Centre::involves (atom_number_t a1, atom_number_t a2) const
{
  assert(INVALID_ATOM_NUMBER != a1);
  assert(INVALID_ATOM_NUMBER != a2);
  assert(a1 != a2);

  assert(ok());

  if (_a == a1)
    ;
  else if (_top_front == a1)
    ;
  else if (_top_back == a1)
    ;
  else if (_left_down == a1)
    ;
  else if (_right_down == a1)
    ;
  else
    return 0;

  if (_a == a2)
    return 1;
  if (_top_front == a2)
    return 1;
  if (_top_back == a2)
    return 1;
  if (_left_down == a2)
    return 1;
  if (_right_down == a2)
    return 1;

  return 0;
}

int
Chiral_Centre::make_copy (const Chiral_Centre & rhs, const int * xref)
{
  assert(rhs.complete());

  _a = xref[rhs._a];

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == rhs._top_front)
    _top_front = CHIRAL_CONNECTION_IS_LONE_PAIR;
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == rhs._top_front)
    _top_front = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else
    _top_front = xref[rhs._top_front];

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == rhs._top_back)
    _top_back = CHIRAL_CONNECTION_IS_LONE_PAIR;
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == rhs._top_back)
    _top_back = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else
    _top_back = xref[rhs._top_back];

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == rhs._left_down)
    _left_down = CHIRAL_CONNECTION_IS_LONE_PAIR;
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == rhs._left_down)
    _left_down = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else
    _left_down = xref[rhs._left_down];

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == rhs._right_down)
    _right_down = CHIRAL_CONNECTION_IS_LONE_PAIR;
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == rhs._right_down)
    _right_down = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else
    _right_down = xref[rhs._right_down];

  _chirality_known = rhs._chirality_known;

  return ok();
}

int
Chiral_Centre::invert()
{
  assert(complete());

  atom_number_t tmp = _left_down;

  _left_down = _right_down;
  _right_down = tmp;

  return 1;
}

/*
  Most likely someone is deleting all the H atoms.
*/

int
Chiral_Centre::convert_to_implicit_hydrogen (atom_number_t h)
{
//cerr << "Chiral_Centre::convert_to_implicit_hydrogen:atom " << h << " was an explicit hydrogen\n";
  if (h == _top_front)
    _top_front = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else if (h == _top_back)
    _top_back = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else if (h == _left_down)
    _left_down = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else if (h == _right_down)
    _right_down = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else
  {
    cerr << "Chiral_Centre::convert_to_implicit_hydrogen: does not involve " << h << endl;
    debug_print (cerr);
    assert(NULL == "This should not happen");
  }

  return 1;
}

//#define DEBUG_ADJUST_CHIRAL_CENTRES_FOR_LOSS_OF_ATOM

int
Molecule::_adjust_chiral_centres_for_loss_of_atom (atom_number_t a,
                                                   int was_hydrogen)
{
#ifdef DEBUG_ADJUST_CHIRAL_CENTRES_FOR_LOSS_OF_ATOM
  cerr << "Removing atom " << a << " Hydrogen(?) = " << was_hydrogen << " need to check " << _chiral_centres.number_elements() << " existing chiral centres\n";
#endif

  int nc = _chiral_centres.number_elements();
  int rc = 0;
  for (int i = 0; i < nc; i++)
  {
    Chiral_Centre * cc = _chiral_centres[i];

    if (! cc->involves(a))
    {
      cc->adjust_for_loss_of_atom(a);
      continue;
    }

#ifdef DEBUG_ADJUST_CHIRAL_CENTRES_FOR_LOSS_OF_ATOM
    cerr << "Chiral centre " << i << " involves atom " << a << endl;
    if (was_hydrogen && 0 == cc->implicit_hydrogen_count())
      cerr << "Will merely convert to implicit hydrogen\n";
#endif

    if (was_hydrogen && 0 == cc->implicit_hydrogen_count())
    {
      cc->convert_to_implicit_hydrogen(a);
      cc->adjust_for_loss_of_atom(a);
    }
    else
    {
      _things[cc->a()]->set_implicit_hydrogens_known(0);   // hopefully true in most cases

      _chiral_centres.remove_item(i);
      i--;
      nc--;
      rc++;
    }
  }

  return rc;
}

int
Molecule::_stereo_centre_hydrogens_become_implicit (Chiral_Centre * c)
{
  if (1 == atomic_number(c->top_front()))
    c->set_top_front(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  else if (1 == atomic_number(c->top_back()))
    c->set_top_back(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  else if (1 == atomic_number(c->left_down()))
    c->set_left_down(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  else if (1 == atomic_number(c->right_down()))
    c->set_right_down(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  else
  {
    c->debug_print(cerr);
    assert(NULL == "How could this happen?");
  }

  return 1;
}

/*
  Frequently when removing explicit hydrogens, we nevertheless want
  to preserve stereo info. Therefore before deleting the hydrogen,
  we can ask the molecule to convert all explicit H's to implicit H's
*/

int
Molecule::stereo_centre_hydrogens_become_implicit()
{
  assert(ok());

  int rc = 0;
  int nc = _chiral_centres.number_elements();
  for (int i = 0; i < nc; i++)
  {
    Chiral_Centre * c = _chiral_centres[i];
    if (1 == c->implicit_hydrogen_count())   // can have only 1 implicit H
      continue;

    rc += _stereo_centre_hydrogens_become_implicit(c);
  }

  return rc;
}

int
Chiral_Centre::set_top_front (atom_number_t ntf)
{
  assert(ok());

  if (ntf == _a ||
      ntf == _top_back ||
      ntf == _left_down ||
      ntf == _right_down)
  {
    cerr << "Chiral_Centre::set_top_front: atom " << ntf << " not allowed\n";
    debug_print(cerr);
    return 0;
  }

  _top_front = ntf;

  return 1;
}

int
Chiral_Centre::set_top_back (atom_number_t ntb)
{
  assert(ok());

  if (ntb == _a ||
      ntb == _top_front ||
      ntb == _left_down ||
      ntb == _right_down)
  {
    cerr << "Chiral_Centre::set_top_back: atom " << ntb << " not allowed\n";
    debug_print(cerr);
    return 0;
  }

  _top_back = ntb;

  return 1;
}

int
Chiral_Centre::set_left_down (atom_number_t nld)
{
  assert(ok());
  if (nld == _a ||
      nld == _top_front ||
      nld == _top_back ||
      nld == _right_down)
  {
    cerr << "Chiral_Centre::set_left_down: atom " << nld << " not allowed\n";
    debug_print(cerr);
    return 0;
  }

  _left_down = nld;

  return 1;
}

int
Chiral_Centre::set_right_down (atom_number_t nrd)
{
  assert(ok());
  if (nrd == _a ||
      nrd == _top_front ||
      nrd == _top_back ||
      nrd == _left_down)
  {
    cerr << "Chiral_Centre::set_right_down: atom " << nrd << " not allowed\n";
    debug_print(cerr);
    return 0;
  }

  _right_down = nrd;

  return 1;
}

void
Chiral_Centre::adjust_for_loss_of_atom (atom_number_t l)
{
  if (_a > l)
    _a--;

  if (_top_front > l)
    _top_front--;
  if (_top_back > l)
    _top_back--;
  if (_left_down > l)
    _left_down--;
  if (_right_down > l)
    _right_down--;

  return;
}

int
Chiral_Centre::adjust_atom_numbers (const int * xref)
{
  if (_a >= 0)
    _a = xref[_a];

  if (_top_front >= 0)
    _top_front = xref[_top_front];
  if (_top_back >= 0)
    _top_back = xref[_top_back];
  if (_left_down >= 0)
    _left_down = xref[_left_down];
  if (_right_down >= 0)
    _right_down = xref[_right_down];

  return ok();
}

/*
  When creating a molecule, there are cases where one of the "connections"
  to a chiral centre is a lone pair. This function is called whenever
  that situation has been found to be true.

  Top front stays the same in all cases.

  Either _left_down or _right_down must be unknown.

  Consider the two cases

  _left_down:

          tb                                   Lp
          |                                     |
          |               becomes               |
         / \                                   / \
        /   \                                 /   \
      Unk    rd                              rd    tb

  _right_down

          tb                                   Lp
          |                                     |
          |               becomes               |
         / \                                   / \
        /   \                                 /   \
       ld   Unk                              tb    ld

*/

int
Chiral_Centre::centre_atom_has_a_lone_pair()
{
  if (INVALID_ATOM_NUMBER == _left_down)
  {
    _left_down = _right_down;
    _right_down = _top_back;
    _top_back = CHIRAL_CONNECTION_IS_LONE_PAIR;
  }
  else if (INVALID_ATOM_NUMBER == _right_down)
  {
    _right_down = _left_down;
    _left_down = _top_back;
    _top_back = CHIRAL_CONNECTION_IS_LONE_PAIR;
  }
  else
  {
    cerr << "Chiral_Centre::centre_atom_has_a_lone_pair: Both left_down and right_down set\n";
    debug_print(cerr);
//  abort();
    cerr << "IGNORED - check your molecule!\n";
    return 0;
  }

  return 1;
}

//#define DEBUG_SMI_PROCESS_NEW_CHIRAL_CENTRE

/*
  In processing a smiles, we have encountered a chiral centre.
  At this stage, that atom must be singly connected.
*/

int
Molecule::_smi_process_new_chiral_centre (Chiral_Centre * c,
                                          int hcount) const
{
  atom_number_t a = c->a();
  assert(a == _number_elements - 1);
  assert(hcount >= 0 && hcount <= 1);

  int acon = _things[_number_elements - 1]->ncon();

#ifdef DEBUG_SMI_PROCESS_NEW_CHIRAL_CENTRE
  cerr << "Processing new chiral centre, atom " << a << ", hcount = " << hcount << endl;
  cerr << "Atom has " << acon << " connections\n";
  if (1 == acon)
    cerr << "Single connection to atom " << other(a, 0) << endl;
#endif

  if (0 == acon)
  {
    if (hcount)
      c->set_top_front(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
/*  else
    {
      cerr << "Molecule::_smi_process_new_chiral_centre: zero hcount to acon = 0\n";
      cerr << "Smiles starts with zero connected atom with no implicit H atoms\n";
      cerr << "Don't know how to parse these...\n";
      return 0;
    }*/
  }
  else if (1 == acon)
  {
    atom_number_t j = other(a, 0);
    c->set_top_front(j);
    if (hcount)
      c->set_top_back(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  }
  else
  {
    cerr << "Molecule::_smi_process_new_chiral_centre: atom " << a << " ncon = " << acon << endl;
    assert(NULL == "How could this happen?");
  }

  c->set_chirality_known(1);    // chiral centres derived from smiles are of known chirality

#ifdef DEBUG_SMI_PROCESS_NEW_CHIRAL_CENTRE
  cerr << "New chiral centre constructed\n";
  c->debug_print(cerr);
#endif

  return 1;
}

Chiral_Centre *
Molecule::chiral_centre_at_atom (atom_number_t a) const
{
  int nc = _chiral_centres.number_elements();
  for (int i = 0; i < nc; i++)
  {
    Chiral_Centre * cc = _chiral_centres[i];
    if (a == cc->a())
      return cc;
  }

  return NULL;
}

Chiral_Centre *
Molecule::remove_no_delete_chiral_centre_at_atom (atom_number_t zatom)
{
  int nc = _chiral_centres.number_elements();

  for (int i = 0; i < nc; i++)
  {
    Chiral_Centre * c = _chiral_centres[i];

    if (zatom != c->a())
      continue;

    _chiral_centres.remove_no_delete(i);

    _set_modified();

    return c;
  }

  return NULL;
}

Chiral_Centre *
Molecule::chiral_centre_in_molecule_not_indexed_by_atom_number (int i) const
{
  assert(ok());
  assert(_chiral_centres.ok_index(i));

  return _chiral_centres[i];
}

//#define DEBUG_SMI_LAST_ATOM_IS_PART_OF_CHIRAL_CENTRE

/*
  In parsing a smiles, we have found an atom bonded to a chiral centre.
*/

int
Molecule::_smi_atom_bonded_to_chiral_centre (atom_number_t previous_atom,
                                             int previous_atom_chiral_count,
                                             atom_number_t atom_bonded_to_chiral_centre)
{
// First fetch the chiral centre anchored on PREVIOUS_ATOM

  Chiral_Centre * c = chiral_centre_at_atom(previous_atom);
  if (NULL == c)
  {
    cerr << "Molecule::_smi_last_atom_is_part_of_chiral_centre: no chiral atom found\n";
    cerr << "Atom " << previous_atom << endl;
    return 0;
  }

#ifdef DEBUG_SMI_LAST_ATOM_IS_PART_OF_CHIRAL_CENTRE
  cerr << "Got another connection to chiral atom " << previous_atom << endl;
  c->debug_print(cerr);
#endif

// In the case where the first atom in the molecule has chirality by no Hydrogen, top_front will be unset

  if (INVALID_ATOM_NUMBER == c->top_front())
  {
    c->set_top_front(atom_bonded_to_chiral_centre);
    return 1;
  }

// The second atom of the specification is always top_back

  if (INVALID_ATOM_NUMBER == c->top_back())
  {
    c->set_top_back(atom_bonded_to_chiral_centre);
    return 1;
  }

  if (c->complete())
  {
    if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == c->top_front())
      c->set_top_front(atom_bonded_to_chiral_centre);
    else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == c->top_back())
      c->set_top_back(atom_bonded_to_chiral_centre);
    else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == c->left_down())
      c->set_left_down(atom_bonded_to_chiral_centre);
    else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == c->right_down())
      c->set_right_down(atom_bonded_to_chiral_centre);
    else
    {
      cerr << "Molecule::_smi_last_atom_is_part_of_chiral_centre: chiral specifier is complete\n";
      cerr << "Atom " << atom_bonded_to_chiral_centre << endl;
      c->debug_print(cerr);
      return 0;
    }

    atom_number_t a = c->a();

    _things[a]->set_implicit_hydrogens(0, 1);

    cerr << "Molecule::_smi_last_atom_is_part_of_chiral_centre:replaced implicit Hydrogen with connection\n";
    cerr << "Atom " << atom_bonded_to_chiral_centre << endl;
    return 1;
  }

// The atom will be either left_down or right_down

  if (1 == previous_atom_chiral_count)     // type @, counterclockwise
  {
    if (INVALID_ATOM_NUMBER == c->left_down())    // this is the 3rd atom
    {
      if (INVALID_ATOM_NUMBER != c->right_down())
      {
        cerr << "Molecule::_smi_last_atom_is_part_of_chiral_centre:already full!\n";
        c->debug_print(cerr);
        return 0;
      }

      c->set_left_down(atom_bonded_to_chiral_centre);
    }
    else     // this is the 4th atom
    {
      if (INVALID_ATOM_NUMBER != c->right_down())
      {
        cerr << "Molecule::_smi_last_atom_is_part_of_chiral_centre:already full!\n";
        c->debug_print(cerr);
        return 0;
      }

      c->set_right_down(atom_bonded_to_chiral_centre);
    }
  }
  else if (2 == previous_atom_chiral_count)    // type @@, going clockwise
  {
    if (INVALID_ATOM_NUMBER == c->right_down())    // 3rd atom
    {
      if (INVALID_ATOM_NUMBER != c->left_down())
      {
        cerr << "Molecule::_smi_last_atom_is_part_of_chiral_centre:already full\n";
        c->debug_print(cerr);
        return 0;
      }

      c->set_right_down(atom_bonded_to_chiral_centre);
    }
    else    // this is the 4th atom of the chiral centre
    {
      if (INVALID_ATOM_NUMBER != c->left_down())
      {
        cerr << "Molecule::_smi_last_atom_is_part_of_chiral_centre:already full\n";
        c->debug_print(cerr);
        return 0;
      }

      c->set_left_down(atom_bonded_to_chiral_centre);
    }
  }
  else
  {
    cerr << "Bad chiral count " << previous_atom_chiral_count << endl;
    iwabort();
  }

#ifdef DEBUG_SMI_LAST_ATOM_IS_PART_OF_CHIRAL_CENTRE
  cerr << "After adding atom ";
  c->debug_print(cerr);
#endif

  return 1;
}

//#define DEBUG_CHECK_FOR_INCOMPLETE_CHIRAL_SPECIFICATIONS

/*
  Checking for incomplete chiral specifications became more difficult
  once we allow lone pairs.
  When only H atoms fill out the valence, the connectivity is known
  when the SMILES is read in. But when a lone pair fills out a chiral
  centre, that cannot be discerned until after the molecule is 
  constructed.
*/

int
Molecule::_check_for_incomplete_chiral_specifications (Chiral_Centre * c)
{
  atom_number_t zatom = c->a();

#ifdef DEBUG_CHECK_FOR_INCOMPLETE_CHIRAL_SPECIFICATIONS
  cerr << "Checking for incomplete chiral specification\n";
  c->debug_print(cerr);
#endif

  if (atom_not_set(c->top_front()) ||
      atom_not_set(c->top_back())  ||
      atom_not_set(c->left_down()) ||
      atom_not_set(c->right_down()))
  {
    int notset = atom_not_set(c->top_front()) + atom_not_set(c->top_back()) + atom_not_set(c->left_down()) + atom_not_set(c->right_down());
    if (1 == notset && _automatically_add_implicit_hydrogen_to_incomplete_chiral_centre && 0 == c->implicit_hydrogen_count())
    {
      cerr << "Molecule::_check_for_incomplete_chiral_specifications:filling missing implicit hydrogen on atom " << c->a() <<endl;
      if (atom_not_set(c->top_front()))
        c->set_top_front(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
      else if (atom_not_set(c->top_back()))
        c->set_top_back(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
      else if (atom_not_set(c->left_down()))
        c->set_left_down(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
      else if (atom_not_set(c->right_down()))
        c->set_right_down(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
      else
      {
        cerr << "Molecule::_check_for_incomplete_chiral_specifications:nset inconsistency\n";
        return 0;
      }

      _things[zatom]->set_implicit_hydrogens(1, 1);   // 2nd arg means over-ride known value

      return 1;
    }

    cerr << "Molecule::_check_for_incomplete_chiral_specifications: incomplete chiral centre, type " << atomic_symbol(c->a()) << endl;
    c->debug_print(cerr);
    return 0;
  }

  int acon = _things[zatom]->ncon();

  if (acon < 3 || acon > 4)
  {
    cerr << "Molecule::_check_for_incomplete_chiral_specifications: Atom " << zatom <<
            " has " << acon << " connections\n";
    return 0;
  }

  if (4 == acon)
  {
    cerr << "Molecule::_check_for_incomplete_chiral_specifications: Atom " << zatom <<
            " has 4 connections, but incomplete!!!\n";
    print_chiral_centre_details(c, cerr);
    return 0;
  }

  assert(3 == acon);

  if (1 != implicit_hydrogens(c->a()))
  {
    cerr << "Molecule::_check_for_incomplete_chiral_specifications: Atom " << zatom <<
            " has " << implicit_hydrogens(c->a()) << " implicit hydrogens\n";
    return 0;
  }

  if (INVALID_ATOM_NUMBER == c->left_down())
    c->set_left_down(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  else if (INVALID_ATOM_NUMBER == c->right_down())
    c->set_right_down(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  else
  {
    assert(NULL == "this should not happen");
  }

  return 1;
}

int
Molecule::_check_for_incomplete_chiral_specifications()
{
  int rc = 1;    // assume OK unless we run into a problem

//cerr << "Molecule::_check_for_incomplete_chiral_specifications:have " << _chiral_centres.number_elements() << " chiral centres to check\n";

  for (int i = _chiral_centres.number_elements() - 1; i >= 0; i--)
  {
    Chiral_Centre * c = _chiral_centres[i];

    int lp;
    if (lone_pair_count(c->a(), lp) && 1 == lp)
    {
      const Atom * a = _things[c->a()];
      if (16 == a->atomic_number() && 4 == a->ncon())   // not sure what chirality this is
        ;
      else
        c->centre_atom_has_a_lone_pair();
    }

    if (c->complete())
      continue;

    if (_check_for_incomplete_chiral_specifications(c))
      continue;

    atomic_number_t cz = _things[c->a()]->atomic_number();

    if (7 == cz || 14 == cz)    // not sure what do do with Chiral Silicon or Nitrogen
    {
      cerr << "Chiral centre at atomic number " << cz << " ignored\n";
      _chiral_centres.remove_item(i);
      continue;
    }

    _chiral_centres.remove_item(i);

    if (! ignore_incorrect_chiral_input())
      rc = 0;
  }

//cerr << "Molecule::_check_for_incomplete_chiral_specifications:at end of processing have " << _chiral_centres.number_elements() << " chiral centres\n";

  return rc;

}

/*
  When forming the smiles for a chiral atom, we need to order
  the connections. If a connection is an implicit hydrogen, it
  is considered the "first" neighbour.
  If A is an implicit H, we return a number which is (hopefully) less
  than zorder of all other connections.
  If the connection to A is a ring opening, give it an order corresponding
  to which ring opening it is. The rule is that the first ring opening
  gets the low number
*/

static int
assign_zorder (atom_number_t a,
               const int * zorder,
               const resizable_array<atom_number_t> & ring_openings,
               const resizable_array<atom_number_t> & ring_closures)
{
  int nro = ring_openings.number_elements();
  int nrc = ring_closures.number_elements();

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == a)
    return -(nro + nrc) - 5;   // lower than anything else

  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == a)
    return -(nro + nrc) - 1;  // lower than anything else returned here

// When constructing a smiles, ring closures are done before openings

  int i = ring_closures.index(a);
  if (i >= 0)
    return -(nro + nrc) + i;

  i = ring_openings.index(a);
  if (i >= 0)
    return -nro + i;

// At this stage, atom A is neither a ring closure, nor a ring opening

  return zorder[a];
}

//#define DEBUG_APPEND_SMILES_CHIRALITY_SYMBOL

/*
  In building a smiles, the path has identified three atoms,
  SOUTH_WEST, NORTH and SOUTH_EAST as being the next three
  connections to _A.
  Do these three atoms actually look like SW, N, SE

  We must be particularly careful about implicit hydrogens, because
  the standard is that the implicit H is the FIRST connection.
*/

int
Chiral_Centre::_discern_clockwise (IWString & smiles,
                                   const int * zorder,
                                   atom_number_t south_west,
                                   atom_number_t north,
                                   atom_number_t south_east,
                                   const resizable_array<atom_number_t> & ring_openings,
                                   const resizable_array<atom_number_t> & ring_closures) const
{
  const int zsw    = assign_zorder(south_west, zorder, ring_openings, ring_closures);
  const int znorth = assign_zorder(north, zorder, ring_openings, ring_closures);
  const int zse    = assign_zorder(south_east, zorder, ring_openings, ring_closures);

#ifdef DEBUG_APPEND_SMILES_CHIRALITY_SYMBOL
  cerr << "South West = " << south_west << " zsw = " << zsw << endl;
  cerr << "North      = " << north      << " zn  = " << znorth << endl;
  cerr << "South East = " << south_east << " zse = " << zse << endl;
#endif

  if (zsw < znorth && znorth < zse)
    smiles += '@';
  else if (znorth < zse && zse < zsw)
    smiles += '@';
  else if (zse < zsw && zsw < znorth)
    smiles += '@';

  return 1;
}

/*
  A smiles is being constructed
  We have passed from PREVIOUS_ATOM to _a. The ordering is as given
  by ZORDER. If atom _A has any ring openings, the neighbouring atoms 
  at the end of the ring openings are in RING_OPENINGS
*/

int
Chiral_Centre::append_smiles_chirality_symbol (IWString & smiles,
                         const int * zorder,
                         atom_number_t previous_atom,
                         const resizable_array<const Bond *> & ring_opening_bonds,
                         const resizable_array<atom_number_t> & ring_closures) const
{
  resizable_array<atom_number_t> ring_openings;

  const int nrob = ring_opening_bonds.number_elements();
  ring_openings.resize(nrob);
  for (int i = 0; i < nrob; i++)
  {
    const Bond * b = ring_opening_bonds[i];
    ring_openings.add(b->other(_a));
  }

  assert(_a != previous_atom);

  smiles += '@';    // there will be at least one of these, from now on,
                    // we just need to identify the clockwise (@@) centres

// Looking down the PREVIOUS_ATOM - _A bond, the other bonds look like:

/*          North
//            |
//            |
//           / \
//          /   \
//         SW    SE
*/

// Pass to discern_clockwise SOUTH_WEST, NORTH, SOUTH_EAST

#ifdef DEBUG_APPEND_SMILES_CHIRALITY_SYMBOL
  cerr << "Append_smiles_chirality_symbol from " << previous_atom << " (order " << zorder[previous_atom] << ") to atom " << _a << " (order " << zorder[_a] << ")";
  if (ring_openings.number_elements())
  {
    cerr << " Ring openings to:";
    for (int i = 0; i < ring_openings.number_elements(); i++)
      cerr << " " << ring_openings[i];
  }
  cerr << endl;
  if (ring_closures.number_elements())
  {
    cerr << " Ring closures to:";
    for (int i = 0; i < ring_closures.number_elements(); i++)
      cerr << " " << ring_closures[i];
  }
  cerr << endl;
#endif

  if (previous_atom == top_front())
    _discern_clockwise(smiles, zorder, _left_down, _top_back, _right_down, ring_openings, ring_closures);
  else if (previous_atom == top_back())
    _discern_clockwise(smiles, zorder, _right_down, _top_front, _left_down, ring_openings, ring_closures);
  else if (previous_atom == left_down())
    _discern_clockwise(smiles, zorder, _top_front, _right_down, _top_back, ring_openings, ring_closures);
  else if (previous_atom == right_down())
    _discern_clockwise(smiles, zorder, _top_back, _left_down, _top_front, ring_openings, ring_closures);

  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_front)
    _discern_clockwise(smiles, zorder, _left_down, _top_back, _right_down, ring_openings, ring_closures);
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_back)
    _discern_clockwise(smiles, zorder, _right_down, _top_front, _left_down, ring_openings, ring_closures);
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _left_down)
    _discern_clockwise(smiles, zorder, _top_front, _right_down, _top_back, ring_openings, ring_closures);
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _right_down)
    _discern_clockwise(smiles, zorder, _top_back, _left_down, _top_front, ring_openings, ring_closures);

  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_front)
    _discern_clockwise(smiles, zorder, _left_down, _top_back, _right_down, ring_openings, ring_closures);
  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_back)
    _discern_clockwise(smiles, zorder, _right_down, _top_front, _left_down, ring_openings, ring_closures);
  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == _left_down)
    _discern_clockwise(smiles, zorder, _top_front, _right_down, _top_back, ring_openings, ring_closures);
  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == _right_down)
    _discern_clockwise(smiles, zorder, _top_back, _left_down, _top_front, ring_openings, ring_closures);
  else
  {
    cerr << "Chiral_Centre::append_smiles_chirality_symbol: huh, previous atom is " << previous_atom << endl;
    cerr << "Smiles is '" << smiles << "'\n";
    debug_print(cerr);
//  assert(NULL == "this should not happen");
  }

  return 1;
}

static int
do_orientation (unsigned int r1, unsigned int r2, unsigned int r3)
{
//cerr << "Doing ranking of " << r1 << " " << r2 << " " << r3 << endl;

  if (r1 > r2 && r2 > r3)
    return 1;

  if (r2 > r3 && r3 > r1)
    return 1;

  if (r3 > r1 && r1 > r2)
    return 1;

  if (r1 == r2 || r2 == r3 || r1 == r3)
    return 0;

  return -1;
}

static int
do_orientation (const unsigned int * rank,
                atom_number_t a1,
                atom_number_t a2)
{
  int r1 = rank[a1];
  int r2 = rank[a2];

  if (r1 == r2)
    return 0;

  if (r1 < r2)
    return -1;

  return 1;
}

/*
  The chiral centre involves a lone pair. Examine the remaining 3 connections.
  Lone pairs must be treated specially
*/

static int
do_orientation (const unsigned int * rank,
                atom_number_t north,
                atom_number_t sw,
                atom_number_t se)
{
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == north)
    return do_orientation(rank, sw, se);
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == sw)
    return do_orientation(rank, se, north);
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == se)
    return do_orientation(rank, north, sw);

// No lone pairs, do normal processing

  return do_orientation(rank[north], rank[sw], rank[se]);
}

//#define DEBUG_CC_ORIENTATION
#ifdef DEBUG_CC_ORIENTATION

static void
write_rank (const char * which_rank,
            int zatom,
            const unsigned int * rank)
{
  cerr << which_rank;
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == zatom || CHIRAL_CONNECTION_IS_LONE_PAIR == zatom)
    cerr << "*";
  else
    cerr << rank[zatom];

  return;
}

#endif

/*
  When unique smiles determinations are being done with chirality, we need
  a means getting an indication of the orientation at a centre
*/

int
Chiral_Centre::orientation (const unsigned int * rank) const
{
#ifdef DEBUG_CC_ORIENTATION
  cerr << "Orientation";
  write_rank(" tf ", _top_front, rank);
  write_rank(" tb ", _top_back, rank);
  write_rank(" ld ", _left_down, rank);
  write_rank(" rd ", _right_down, rank);
  cerr << endl;
#endif

// Deal with implicit hydrogens first

  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_front)
    return do_orientation(rank, _top_back, _left_down, _right_down);
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_back)
    return do_orientation(rank, _top_front, _right_down, _left_down);
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _right_down)
    return do_orientation(rank, _left_down, _top_back, _top_front);
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _left_down)
    return do_orientation(rank, _top_front, _top_back, _right_down);

// then lone pairs

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_front)
    return do_orientation(rank[_top_back], rank[_left_down], rank[_right_down]);

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_back)
    return do_orientation(rank[_top_front], rank[_right_down], rank[_left_down]);

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _left_down)
    return do_orientation(rank[_top_front], rank[_top_back], rank[_right_down]);

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _right_down)
    return do_orientation(rank[_top_back], rank[_top_front], rank[_left_down]);

// Looks like all connections are actual atoms, fetch their ranks

  int rtf = rank[_top_front];
  int rtb = rank[_top_back];
  int rld = rank[_left_down];
  int rrd = rank[_right_down];

// If any are duplicates, we cannot resolve the centre

  if (rtf == rtb || rtf == rld || rtf == rrd || rtb == rld || rtb == rrd || rld == rrd)
    return 0;

// All ranks are different, now this gets ugly

  if (rtf > rtb && rtf > rld && rtf > rrd)
    return do_orientation(rtb, rld, rrd);
  if (rtb > rtf && rtb > rld && rtb > rrd)
    return do_orientation(rtf, rrd, rld);
  if (rld > rtf && rld > rtb && rld > rrd)
    return do_orientation(rtf, rtb, rrd);
  if (rrd > rtf && rrd > rtb && rrd > rld)
    return do_orientation(rtb, rtf, rld);

  assert(NULL == "Should not come to here");

  return 0;
}

/*static int
do_chiral_influence (const unsigned int * rank,
                     atom_number_t left,
                     atom_number_t right, 
                     atom_number_t a1,
                     atom_number_t a2)
{
  unsigned int rl = rank[left];
  unsigned int rr = rank[right];

  if (rl == rr)
    return 0;

  if (rl < rr)
    return -1;
  else
    return 1;
}*/

/*
  The chiral centre had just a single lone pair or just a single implicit hydrogen
*/

/*static int
do_chiral_influence (const unsigned int r1, const unsigned int r2)
{
  if (r1 == r2)
    return 0;

  if (r1 < r2)
    return -1;
  else
    return 1;
}*/

/*
  The chiral centre had both a lone pair and an implicit hydrogen. Those
  two atoms have defined an orientation
*/

/*static int
do_chiral_influence (atom_number_t left,
                     atom_number_t right,
                     atom_number_t a1,
                     atom_number_t a2)
{
  if (a1 == left && a2 == right)
   return 1;

  if (a2 == left && a1 == right)
   return -1;

  abort();

  return 0;
}*/

/*
  One connection has been fixed
*/

/*static int
do_chiral_influence (const unsigned int * rank,
                     atom_number_t north,
                     atom_number_t sw,
                     atom_number_t se,
                     atom_number_t a1,
                     atom_number_t a2)
{

// All three of north, se and sw are actual atoms. One of them isn't a1 or a2

  if (a1 != north && a2 != north)
    return do_chiral_influence(se, sw, a1, a2);
  if (a1 != se && a2 != se)
    return do_chiral_influence(sw, north, a1, a2);
  if (a1 != sw && a2 != sw)
    return do_chiral_influence(north, se, a1, a2);

  cerr << "do_chiral_influence: yipes, no match found\n";
  abort();
  return 0;
}*/

/*
  We come here when the chiral centre has an implicit hydrogen.
  Check to see if it also has a lone pair
*/

/*static int
do_chiral_influence_check_lp (const unsigned int * rank,
                              atom_number_t north,
                              atom_number_t sw,
                              atom_number_t se,
                              atom_number_t a1,
                              atom_number_t a2)
{
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == north)
    return do_chiral_influence(rank, se, sw, a1, a2);
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == sw)
    return do_chiral_influence(rank, north, se, a1, a2);
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == se)
    return do_chiral_influence(rank, sw, north, a1, a2);

// no special treatment needed

  return do_chiral_influence (rank, north, sw, se, a1, a2);
}*/

/*
  Used during unique determinations
*/

/*int
Chiral_Centre::influence (const unsigned int * rank,
                          atom_number_t a1, atom_number_t a2) const
{
  assert(rank[a1] == rank[a2]);    // these are the atoms we are trying to resolve

  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_front)
    return do_chiral_influence_check_lp (rank, _top_back, _left_down, _right_down, a1, a2);
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_back)
    return do_chiral_influence_check_lp (rank, _top_front, _right_down, _left_down, a1, a2);
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _right_down)
    return do_chiral_influence_check_lp (rank, _top_front, _top_back, _left_down, a1, a2);
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _left_down)
    return do_chiral_influence_check_lp (rank, _top_front, _top_back, _right_down, a1, a2);

// then lone pairs

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_front)
    return do_chiral_influence (rank, _top_back, _left_down, _right_down, a1, a2);

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_back)
    return do_chiral_influence (rank, _top_front, _right_down, _left_down, a1, a2);

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _right_down)
    return do_chiral_influence (rank, _top_back, _top_front, _left_down, a1, a2);

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _left_down)
    return do_chiral_influence (rank, _top_front, _top_back, _right_down, a1, a2);

// All the atoms are actual atoms. Match the atoms up. Put the A1 match
// along the Z axis and look for anticlockwise from the A2 atom

  if (_top_front == a1)
  {
    if (_top_back == a2)
      return do_chiral_influence (rank[_left_down], rank[_right_down]);
    if (_left_down == a2)
      return do_chiral_influence (rank[_right_down], rank[_top_back]);
    if (_right_down == a2)
      return do_chiral_influence (rank[_top_back], rank[_left_down]);
  }

  if (_top_back == a1)
  {
    if (_top_front == a2)
      return do_chiral_influence (rank[_right_down], rank[_left_down]);
    if (_left_down == a2)
      return do_chiral_influence (rank[_top_front], rank[_right_down]);
    if (_right_down == a2)
      return do_chiral_influence (rank[_left_down], rank[_top_front]);
  }

  if (_left_down == a1)
  {
    if (_top_front == a2)
      return do_chiral_influence (rank[_top_back], rank[_right_down]);
    if (_top_back == a2)
      return do_chiral_influence (rank[_right_down], rank[_top_front]);
    if (_right_down == a2)
      return do_chiral_influence (rank[_top_front], rank[_top_back]);
  }

  if (_right_down == a1)
  {
    if (_top_front == a2)
      return do_chiral_influence (rank[_left_down], rank[_top_back]);
    if (_top_back == a2)
      return do_chiral_influence (rank[_top_front], rank[_left_down]);
    if (_left_down == a2)
      return do_chiral_influence (rank[_top_front], rank[_top_back]);
  }

  abort();

  return 0;
}*/

static unsigned int
get_rank_for_connection (const unsigned int * rank, atom_number_t a)
{
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == a)
    return std::numeric_limits<unsigned int>::max();

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == a)
    return std::numeric_limits<unsigned int>::max() - 1;

  return  rank[a];
}

static int
do_influence (const unsigned int * rank,
              atom_number_t north,
              atom_number_t sw,
              atom_number_t se)
{
  unsigned int rank_north = get_rank_for_connection(rank, north);
  unsigned int rank_sw    = get_rank_for_connection(rank, sw);
  unsigned int rank_se    = get_rank_for_connection(rank, se);

  if (rank_north > rank_sw && rank_sw > rank_se)
    return 1;

  if (rank_sw > rank_se && rank_se > rank_north)
    return 1;

  if (rank_se > rank_north && rank_north > rank_sw)
    return 1;

  if (rank_north < rank_sw && rank_sw < rank_se)
    return -1;

  if (rank_sw < rank_se && rank_se < rank_north)
    return -1;

  if (rank_se < rank_north && rank_north < rank_sw)
    return -1;

  return 0;
}

/*
  Atom ZATOM is part of the centre.
*/

int
Chiral_Centre::influence (const unsigned int * rank,
                          atom_number_t zatom) const
{
  if (_top_front == zatom)
    return do_influence(rank, _top_back, _left_down, _right_down);
  if (_top_back == zatom)
    return do_influence(rank, _top_front, _right_down, _left_down);
  if (_left_down == zatom)
    return do_influence(rank, _top_front, _top_back, _right_down);
  if (_right_down == zatom)
    return do_influence(rank, _top_front, _left_down, _top_back);

  return 0;
}

/*
  What is the highest atom number less than BELOW.
  Note that this works only because both CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN
  and INVALID_ATOM_NUMBER are negative quantities
*/

/*int
Chiral_Centre::_highest_valid_atom_number (atom_number_t below) const
{
  iwmax<atom_number_t> amax (-1);

  if (_top_front < below)
    amax.extra (_top_front);
  if (_top_back < below)
    amax.extra (_top_back);
  if (_left_down < below)
    amax.extra (_left_down);
  if (_right_down < below)
    amax.extra (_right_down);

  return amax.maxval();
}*/

int
Chiral_Centre::mdl_stereo_centre_value (atom_number_t sw,
                                        atom_number_t north,
                                        atom_number_t se) const
{
  if (sw < north && north < se)
    return 1;
  if (north < se && se < sw)
    return 1;
  if (se < sw && sw < north)
    return 1;

  return 2;
}

/*
  the caller is writing a MOLFILE, and needs to know whether to put
  a 1, 2 or 3
*/

int
Chiral_Centre::mdl_stereo_centre_value() const
{
  if (0 == _chirality_known)
    return 3;

/*
// MDL wants the highest atom number to be pointing into the page.
// Therefore this is the opposite of what we had for smiles
// Pass the arguments SW N SE
//
//             N
//             |
//            / \
//          SW   SE

*/

  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_front)
    return mdl_stereo_centre_value(_right_down, _top_back, _left_down);
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_back)
    return mdl_stereo_centre_value(_left_down, _top_front, _right_down);
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _left_down)
    return mdl_stereo_centre_value(_top_front, _top_back, _right_down);
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _right_down)
    return mdl_stereo_centre_value(_top_back, _top_front, _left_down);

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_front)
    return mdl_stereo_centre_value(_right_down, _top_back, _left_down);
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_back)
    return mdl_stereo_centre_value(_left_down, _top_front, _right_down);
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _left_down)
    return mdl_stereo_centre_value(_top_front, _top_back, _right_down);
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _right_down)
    return mdl_stereo_centre_value(_top_back, _top_front, _left_down);

// No implicit H

  if (_top_front > _top_back && _top_front > _left_down && _top_front > _right_down)
    return mdl_stereo_centre_value(_right_down, _top_back, _left_down);
  if (_top_back > _top_front && _top_back > _left_down && _top_back > _right_down)
    return mdl_stereo_centre_value(_left_down, _top_front, _right_down);
  if (_left_down > _top_front && _left_down > _top_back && _left_down > _right_down)
    return mdl_stereo_centre_value(_top_front, _top_back, _right_down);
  if (_right_down > _top_front && _right_down > _top_back && _right_down > _left_down)
    return mdl_stereo_centre_value(_top_back, _top_front, _left_down);

  cerr << "Chiral_Centre::mdl_stereo_centre_value: Huh??\n";
  debug_print(cerr);
  iwabort();

  return 0;
}

/*
  Note that implicit_hydrogen_count and lone_pair_count do not call ok
  as ok() calls them!
*/

int
Chiral_Centre::implicit_hydrogen_count() const
{
  int rc = 0;

  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_front)
    rc++;
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_back)
    rc++;
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _left_down)
    rc++;
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _right_down)
    rc++;

  return rc;
}

int
Chiral_Centre::lone_pair_count() const
{
  int rc = 0;

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_front)
    rc++;
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_back)
    rc++;
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _left_down)
    rc++;
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _right_down)
    rc++;

  return rc;
}

/*
  An external process has determined that atom A can be a chiral
  centre. Note that we create it with unspecified chirality

  Note that when calling this from mdl.cc, the atom will have 0
  connections, so we have a flag which suppresses errors when
  0 == ncon

*/

Chiral_Centre *
Molecule::create_chiral_centre (atom_number_t zatom,
                                int zero_connections_ok)
{
  assert(NULL == chiral_centre_at_atom(zatom));    // cannot already be one at atom A

  Atom * a = _things[zatom];

  int lp = 0;

  int acon = a->ncon();
  if (0 == acon && zero_connections_ok)
    ;
  else if (4 == acon)
    ;
  else if (4 == acon + a->implicit_hydrogens())
    ;
  else if (a->lone_pair_count(lp) && 1 == lp)
    ;
  else
  {
    cerr << "Molecule::create_chiral_centre: atom " << zatom << " (" << a->atomic_symbol() << ") has " << acon <<
            " connections\n";
    debug_print(cerr);
    return NULL;
  }

  Chiral_Centre * c = new Chiral_Centre(zatom);
  if (acon)
  {
    c->set_top_front(a->other(zatom, 0));
    c->set_top_back (a->other(zatom, 1));
    c->set_left_down(a->other(zatom, 2));
    if (4 == acon)
      c->set_right_down(a->other(zatom, 3));
    else if (lp)
      c->set_right_down(CHIRAL_CONNECTION_IS_LONE_PAIR);
    else
      c->set_right_down(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  }

  _chiral_centres.add(c);

  _set_modified(zatom);

  return c;
}

//#define DEBUG_COMPLETE_CHIRAL_CENTRE

/*
  In setting the chirality from a molfile, we must follow the
  directions in the MOLFILE manual
  In their parlance, atom 4 is always top_back
  and atom 3 is always top_front

  Note that we treat an implicit hydrogen and a lone pair the same

  I checked with ISIS draw and if you have a Deuterium and an implicit Hydrogen,
  the implicit Hydrogen is still considered the highest numbered atom

  Sulphur is particularly problematic. -S(=O)- seems SP3 hybridised.
  Directional bonds are used near S(=O)(=O) groups to indicate the bond
  angle through the Sulphur. Easiest to just allow Sulphur.
  Phosphorus seems to be the same way...
*/

int
Molecule::_complete_chiral_centre_from_mdl_files (Chiral_Centre * c,
                                                  const MDL_File_Supporting_Material & mdlfos)
{
#ifdef DEBUG_COMPLETE_CHIRAL_CENTRE
  cerr << "Molecule::_complete_chiral_centre_from_mdl_files:completing\n";
  c->debug_print(cerr);
#endif

  atom_number_t a = c->a();

  const Atom * atom_a = _things[a];

  int acon = atom_a->ncon();

  if (acon == atom_a->nbonds())   // good, fully saturated
    ;
  else if (6 == atom_a->atomic_number())  // the most common case
    ;
  else if (16 == atom_a->atomic_number())  // many different possibilities
    ;
  else if (15 == atom_a->atomic_number())   // may indeed be correct
    ;
  else
  {
    cerr << "Molecule::_complete_chiral_centre_from_mdl_files:unsaturated chiral atom, type " << atom_a->atomic_symbol() << "\n";
    c->debug_print(cerr);
    return 0;
  }

  Set_of_Atoms con;

  int explicit_hydrogen_present = 0;
  int hydrogen_isotope_present = 0;    // I suppose they could have D and T and an implicit Hydrogen

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = atom_a->other(a, i);
    if (1 == _things[j]->atomic_number() && mdlfos.mdl_read_h_correct_chiral_centres())
    {
      if (_things[j]->isotope())
      {
        hydrogen_isotope_present++;
        con.add(j);
      }
      else
      {
        explicit_hydrogen_present++;
        con.add(_number_elements + j);    // special atom number that will be larger than any other atom number
      }
    }
    else
      con.add(j);
  }

  int ih = _things[a]->implicit_hydrogens();

#ifdef DEBUG_COMPLETE_CHIRAL_CENTRE
  cerr << "Completing chiral centre at atom " << a << " ncon = " << acon << ", ih = " << ih << endl;
#endif

// An explicit Hydrogen and an implicit Hydrogen is an error - but, what about a partially Deuterated atom!!

  if (ih + explicit_hydrogen_present > 1)
  {
    cerr << "Molecule::_complete_chiral_centre_from_mdl_files: atom " << a << " has " << ih << " implicit hydrogens and " << explicit_hydrogen_present << " explicit\n";
    cerr << "type " << _things[a]->atomic_symbol() << ", " << _things[a]->ncon() << " connections\n";
    return 0;
  }

  int tcon = acon + ih;

// Three connections and a lone pair is fine.

  int lone_pairs = 0;
  if (3 == tcon)
  {
    if ( ! lone_pair_count(a, lone_pairs))
      return 0;
      
    if (1 != lone_pairs)
      return 0;
  }
  else if (acon < 3 || tcon < 3 || tcon > 4)
  {
    cerr << "Molecule::_complete_chiral_centre_from_mdl_files: atom " << a << 
            " has " << acon << " connections\n";
    cerr << "tcon = " << tcon << " type " << _things[a]->atomic_symbol() << endl;

    return 0;
  }

  con.sort(int_comparitor_larger);

  if (explicit_hydrogen_present)
  {
    for (int i = 0; i < acon; i++)
    {
      atom_number_t j = con[i];
      if (j >= _number_elements)
      {
        con[i] = j - _number_elements;
        break;
      }
    }
  }

#ifdef DEBUG_COMPLETE_CHIRAL_CENTRE
  cerr << "The following atom numbers are bonded:";
  for (int i = 0; i < con.number_elements(); i++)
  {
    int j = con[i];
    cerr << ' ' << j << " (" << _things[j]->atomic_symbol() << ')';
  }
  cerr << endl;
  cerr << " ih " << ih << " lone_pairs " << lone_pairs << endl;
#endif

  if (2 == acon && 1 == ih && 1 == lone_pairs)
  {
    c->set_top_front(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
    c->set_top_back(CHIRAL_CONNECTION_IS_LONE_PAIR);
  }
  else
  {
    c->set_top_front(con[2]);

    if (4 == acon)
      c->set_top_back(con[3]);
    else if (ih)
      c->set_top_back(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
    else if (lone_pairs)
      c->set_top_back(CHIRAL_CONNECTION_IS_LONE_PAIR);
    else
    {
      cerr << "Molecule::_complete_chiral_centre_from_mdl_files:three connections by no H or LP\n";
      return 0;
    }
  }

#ifdef DEBUG_COMPLETE_CHIRAL_CENTRE
  cerr << "Input chirality = " << c->chirality_known() << " top done\n";
  c->debug_print(cerr);
  cerr << "Other atoms are 0: " << other(a, 0) << " and 1 " << other(a, 1) << endl;
#endif

// As coded, I believe this to be incorrect, but Concord seems to interpret
// this this way.

  if (1 == c->chirality_known())     // clockwise
  {
    c->set_left_down (con[1]);
    c->set_right_down(con[0]);
  }
  else if (2 == c->chirality_known())  // anti-clockwise
  {
    c->set_left_down (con[0]);
    c->set_right_down(con[1]);
  }
  else
  {
    c->debug_print(cerr);
    assert(NULL == "Unknown chirality_known value");
  }

  return 1;
}

/*
  We have just finished parsing a MOLFILE, and there were
  stereo centres in there. All mdl.cc does is identify the
  stereo centres, we need to check them and fill them out
*/


int
Molecule::_complete_chiral_centres_from_mdl_files(const MDL_File_Supporting_Material & mdlfos)
{
  int rc = 1;    // assume OK until we find a problem

//cerr << "Molecule has " << nc << " chiral centres\n";

  for (int i = _chiral_centres.number_elements() - 1; i >= 0; i--)
  {
    Chiral_Centre * c = _chiral_centres[i];

    if (c->complete())
      continue;

    if (_complete_chiral_centre_from_mdl_files(c, mdlfos))
      continue;

    if (ignore_incorrect_chiral_input())
    {
      cerr << "Discarding invalid chiral centre on atom " << c->a() << " '" << smarts_equivalent_for_atom(c->a()) << "'\n";

      _chiral_centres.remove_item(i);
    }
    else
      rc = 0;
  }

  return rc;
}

int
Molecule::invert_chirality_on_atom (atom_number_t a)
{
  Chiral_Centre * c = chiral_centre_at_atom(a);
  assert(c);

  c->invert();

  invalidate_smiles();

  _symmetry_class_and_canonical_rank.invalidate();

  return 1;
}

int
Molecule::_check_chiral_centres() const
{
  int nc = _chiral_centres.number_elements();
  for (int i = 0; i < nc; i++)
  {
    if (! valid_chiral_centre(_chiral_centres[i]))
    {
      cerr << "Yipes, chiral centre " << i << " is bad, " << smarts_equivalent_for_atom(_chiral_centres[i]->a()) << endl;
      print_chiral_centre_details(_chiral_centres[i], cerr);
      return 0;
    }
  }

  return 1;
}

//#define DEBUG_GOT_RING_OPENING_BOND

/*
  We need a way of representing the fact that one or more of the atoms
  of a chiral_centre are tied to ring closures, and so their atom
  numbers are not yet known
*/

#define CHIRAL_CENTRE_PENDING_RING_CLOSURE(r) (-1000000 - (r))

/*
  In parsing a smiles, a ring opening is found at a chiral atom.
  If 1 == chiral_count, it is of type '@',
  If 2 == chiral_count, it is of type '@@',
*/

int
Chiral_Centre::got_ring_opening_bond (int ring_number, int chiral_count)
{
#ifdef DEBUG_GOT_RING_OPENING_BOND
  cerr << "Ring " << ring_number << " opens from chiral centre ";
  debug_print(cerr);
#endif

  if (INVALID_ATOM_NUMBER == _top_front)
  {
    _top_front = CHIRAL_CENTRE_PENDING_RING_CLOSURE(ring_number);
    return 1;
  }

  if (INVALID_ATOM_NUMBER == _top_back)
  {
    _top_back = CHIRAL_CENTRE_PENDING_RING_CLOSURE(ring_number);
    return 1;
  }

// Jan 2003. Jibo had a molecule with [C@H]23 where 2 was a ring closure and 3 was a ring
// opening. In that case, there will be just one slot open

  if (1 == chiral_count && INVALID_ATOM_NUMBER == _left_down)
  {
    _left_down = CHIRAL_CENTRE_PENDING_RING_CLOSURE(ring_number);
    return 1;
  }

  if (2 == chiral_count && INVALID_ATOM_NUMBER == _right_down)
  {
    _right_down = CHIRAL_CENTRE_PENDING_RING_CLOSURE(ring_number);
    return 1;
  }

  assert(1 == ((INVALID_ATOM_NUMBER == _left_down) + (INVALID_ATOM_NUMBER == _right_down)));    // just one of these should be unset

  if (INVALID_ATOM_NUMBER == _left_down)
  {
    _left_down = CHIRAL_CENTRE_PENDING_RING_CLOSURE(ring_number);
    return 1;
  }

  if (INVALID_ATOM_NUMBER == _right_down)
  {
    _right_down = CHIRAL_CENTRE_PENDING_RING_CLOSURE(ring_number);
    return 1;
  }

  cerr << "Chiral_Centre::got_ring_opening_bond: chiral count = " << chiral_count << endl;
  debug_print(cerr);

  return 0;
}

//#define DEBUG_GOT_RING_CLOSURE_BOND

/*
  In parsing a smiles, a ring closing bond has been encountered
  somewhere in the smiles. The smiles parser detects that atom _A
  is a chiral centre at the other end of that bond.

  This function provides the info that the atom at the other end
  of ring number RING_NUMBER is atom A
*/

int
Chiral_Centre::got_ring_closure_bond (int ring_number, atom_number_t a)
{
  assert(INVALID_ATOM_NUMBER != a);

#ifdef DEBUG_GOT_RING_CLOSURE_BOND
  cerr << "Atom " << a << " closes ring " << ring_number << endl;
  debug_print(cerr);
#endif

  if (_top_front == CHIRAL_CENTRE_PENDING_RING_CLOSURE(ring_number))
  {
    _top_front = a;
    return 1;
  }

  if (_top_back == CHIRAL_CENTRE_PENDING_RING_CLOSURE(ring_number))
  {
    _top_back = a;
    return 1;
  }

  if (_left_down == CHIRAL_CENTRE_PENDING_RING_CLOSURE(ring_number))
  {
    _left_down = a;
    return 1;
  }

  if (_right_down == CHIRAL_CENTRE_PENDING_RING_CLOSURE(ring_number))
  {
    _right_down = a;
    return 1;
  }

  cerr << "Chiral_Centre::got_ring_closure_bond: not pending ring " << ring_number << endl;
  debug_print(cerr);
  iwabort();
  return 0;
}

static int
do_copy (atom_number_t & a,
         atom_number_t initial,
         int diff)
{
  if (initial >= 0)     // valid atom number, not a lone pair or implicit Hydrogen
  {
    a = initial - diff;
    assert(a >= 0);
  }
  else
    a = initial;

  return 1;
}

int
Chiral_Centre::make_copy (const Chiral_Centre * c2)
{
  assert(c2->ok());

  _chirality_known = c2->_chirality_known;

  int diff = c2->_a - _a;

  do_copy(_top_back, c2->_top_back, diff);
  do_copy(_top_front, c2->_top_front, diff);
  do_copy(_left_down, c2->_left_down, diff);
  do_copy(_right_down, c2->_right_down, diff);

  return 1;
}

/*
  Someone is making implicit hydrogens explicit, and our implicit
  hydrogen now has an atom number
*/

int
Chiral_Centre::implicit_hydrogen_is_now_atom_number (atom_number_t a)
{
  if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_back)
    _top_back = a;
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_front)
    _top_front = a;
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _left_down)
    _left_down = a;
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _right_down)
    _right_down = a;
  else
  {
    cerr << "Chiral_Centre::implicit_hydrogen_is_now_atom_number: huh, no implicit hydrogen\n";
    debug_print(cerr);
    cerr << "Proposed atom was " << a << endl;

    return 0;
  }

  return 1;
}

/*
  Similarly a lone pair has been substituted
*/

int
Chiral_Centre::lone_pair_is_now_atom_number (atom_number_t a)
{
  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_back)
    _top_back = a;
  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_front)
    _top_front = a;
  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == _left_down)
    _left_down = a;
  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == _right_down)
    _right_down = a;
  else
  {
    cerr << "Chiral_Centre::implicit_hydrogen_is_now_atom_number: huh, no implicit hydrogen\n";
    debug_print(cerr);
    cerr << "Proposed atom was " << a << endl;

    return 0;
  }

  return 1;
}

/*
  sometimes we want to make a given atom appear in a given position
  within the chiral centre.
*/

int
Chiral_Centre::make_top_front (atom_number_t ntf)
{
  if (ntf == _top_front)
    return 1;

  if (ntf == _top_back)    // lave right_down alone
  {
    _top_back = _left_down;
    _left_down = _top_front;

    _top_front = ntf;

    return 1;
  }

  if (ntf == _left_down)     // leave top back alone
  {
    _left_down = _right_down;
    _right_down = _top_front;

    _top_front = ntf;

    return 1;
  }

  if (ntf == _right_down)    // leave top back alone
  {
    _right_down = _left_down;
    _left_down = _top_front;

    _top_front = ntf;

    return 1;
  }

  cerr << "Chiral_Centre::make_top_front: I don't use atom " << ntf << endl;
  debug_print(cerr);
  iwabort();
  
  return 0;
}

int
Chiral_Centre::atom_is_now_implicit_hydrogen (atom_number_t a)
{
  if (_top_front == a)
    _top_front = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else if (_top_back == a)
    _top_back = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else if (_left_down == a)
    _left_down = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else if (_right_down == a)
    _right_down = CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN;
  else
  {
    cerr << "Chiral_Centre::atom_is_now_implicit_hydrogen: I don't include atom " << a << endl;
    return 0;
  }

  return 1;
}

int
Chiral_Centre::atom_is_now_lone_pair (atom_number_t a)
{
  if (_top_front == a)
    _top_front = CHIRAL_CONNECTION_IS_LONE_PAIR;
  else if (_top_back == a)
    _top_back = CHIRAL_CONNECTION_IS_LONE_PAIR;
  else if (_left_down == a)
    _left_down = CHIRAL_CONNECTION_IS_LONE_PAIR;
  else if (_right_down == a)
    _right_down = CHIRAL_CONNECTION_IS_LONE_PAIR;
  else
  {
    cerr << "Chiral_Centre::atom_is_now_lone_pair: I don't include atom " << a << endl;
    return 0;
  }

  return 1;
}

/*
  this function is mainly used by Molecule::create_subset.

  Are all the atoms which form this chiral centre in the subset

  Note the somewhat unusual convention that subset atoms will
  have a subset[] value >= 0
*/

int
Chiral_Centre::all_atoms_in_subset (const int * subset, int id) const
{
  assert(complete());

  if (subset[_a] < 0)
    return 0;

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_front)
    ;
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_front)
    ;
  else if (id != subset[_top_front])
    return 0;

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _top_back)
    ;
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _top_back)
    ;
  else if (id != subset[_top_back])
    return 0;

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _left_down)
    ;
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _left_down)
    ;
  else if (id != subset[_left_down])
    return 0;

  if (CHIRAL_CONNECTION_IS_LONE_PAIR == _right_down)
    ;
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == _right_down)
    ;
  else if (id != subset[_right_down])
    return 0;

  return 1;    // yep, all our atoms are in the subset
}

int
Molecule::remove_chiral_centre_at_atom (atom_number_t a)
{
  assert(ok());
  assert(a >= 0 && a < _number_elements);

  for (int i = 0; i < _chiral_centres.number_elements(); i++)
  {
    if (a == _chiral_centres[i]->a())
    {
      _chiral_centres.remove_item(i);
      _things[a]->set_implicit_hydrogens_known(0);

      _set_modified();

      return 1;
    }
  }

  cerr << "Molecule::remove_chiral_centre: no chiral centre on atom " << a << endl;

  return 0;
}

int
Molecule::remove_all_chiral_centres()
{
  int rc = _chiral_centres.number_elements();
  if (0 == rc)
    return 0;

  for (int i = 0; i < rc; i++)
  {
    const Chiral_Centre * ci = _chiral_centres[i];
    atom_number_t a = ci->a();
    _things[a]->set_implicit_hydrogens_known(0);
  }

  _chiral_centres.resize(0);

  _set_modified();

  return rc;
}

int
Molecule::_add_chiral_centre_checking_for_duplicate (Chiral_Centre * c)
{
  atom_number_t a = c->a();

  int nc = _chiral_centres.number_elements();
  for (int i = 0; i < nc; i++)
  {
    if (a != _chiral_centres[i]->a())
      continue;

    _chiral_centres.remove_item(i);
    break;
  }

  _chiral_centres.add(c);

  return 1;
}

/*
  Should really check to make sure that there isn't already a chiral centre
  centered on the main atom of this new chiral centre....
*/

int
Molecule::add_chiral_centre (Chiral_Centre * c,
                             int check_for_existing)
{
  assert(valid_chiral_centre(c));

  int rc;
  if (check_for_existing && _chiral_centres.number_elements())
    rc = _add_chiral_centre_checking_for_duplicate(c);
  else
  {
    _chiral_centres.add(c);
    rc = 1;
  }

  invalidate_smiles();

  _symmetry_class_and_canonical_rank.invalidate();

  return rc;
}

int
Molecule::involved_in_chiral_centre (int * s, int v) const
{
  int n = _chiral_centres.number_elements();

  for (int i = 0; i < n; i++)
  {
    _chiral_centres[i]->set_vector(s, v);
  }

  return n;
}

int
Molecule::at_centre_of_chiral_centre (int * s, int v) const
{
  int n = _chiral_centres.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = _chiral_centres[i]->a();

    assert(j >= 0 && j < _number_elements);

    s[j] = v;
  }

  return n;
}

int
Chiral_Centre::set_vector (int * s, int v) const
{
  int rc = 0;

  if (_top_front >= 0)
  {
    s[_top_front] = v;
    rc++;
  }

  if (_top_back >= 0)
  {
    s[_top_back] = v;
    rc++;
  }

  if (_left_down >= 0)
  {
    s[_left_down] = v;
    rc++;
  }

  if (_right_down >= 0)
  {
    s[_right_down] = v;
    rc++;
  }

  return rc;
}

static int
atom_matches (atom_number_t & target, atom_number_t zold, atom_number_t znew)
{
  if (target == zold)
  {
    target = znew;
    return 1;
  }

  return 0;
}

int
Chiral_Centre::change_atom_number (atom_number_t zold, atom_number_t znew)
{
  assert(ok());

  if (atom_matches(_a, zold, znew))
    return 1;
  if (atom_matches(_top_back, zold, znew))
    return 1;
  if (atom_matches(_top_front, zold, znew))
    return 1;
  if (atom_matches(_left_down, zold, znew))
    return 1;
  if (atom_matches(_right_down,  zold, znew))
    return 1;

  cerr << "Chiral_Centre::change_atom_number: no involvement with atom " << zold << endl;
  debug_print(cerr);

  return 0;
}

/*
  Kind of strange function. Atom number ZATOM is moving to be the last atom in the
  molecule. Therefore we need to change the atom numbers of everything larger than
  ZATOM, and change any reference to ZATOM
*/

static void
do_move_atom_to_end_of_atom_list (atom_number_t & target, atom_number_t zatom, int atoms_in_molecule)
{
  if (zatom == target)      // goes to the end of the list
    target = atoms_in_molecule - 1;
  else if (CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN == target)
    ;
  else if (CHIRAL_CONNECTION_IS_LONE_PAIR == target)
    ;
  else if (target > zatom)      // atom ZATOM has been moved out, shift TARGET down
    target--;

  return;
}

int
Chiral_Centre::move_atom_to_end_of_atom_list (atom_number_t zatom, int atoms_in_molecule)
{
  do_move_atom_to_end_of_atom_list(_a,          zatom, atoms_in_molecule);
  do_move_atom_to_end_of_atom_list(_top_front,  zatom, atoms_in_molecule);
  do_move_atom_to_end_of_atom_list(_top_back,   zatom, atoms_in_molecule);
  do_move_atom_to_end_of_atom_list(_left_down,  zatom, atoms_in_molecule);
  do_move_atom_to_end_of_atom_list(_right_down, zatom, atoms_in_molecule);

  return 1;
}

static void
check_atom_numbers_swapped (atom_number_t & target,
                            atom_number_t a1,
                            int & done1,
                            atom_number_t a2,
                            int & done2)
{
  if (! done1 && target == a1)
  {
    target = a2;
    done1 = 1;
  }
  else if (! done2 && target == a2)
  {
    target = a1;
    done2 = 1;
  }

  return;
}

int
Chiral_Centre::atom_numbers_are_swapped (atom_number_t a1, atom_number_t a2)
{
  assert(a1 != a2);

  int done1 = 0;
  int done2 = 0;

#ifdef DEBUG_ATOM_NUMBERS_ARE_SWAPPED
  cerr << "Before swapping " << a1 << " and " << a2 << endl;
  debug_print(cerr);
#endif

  check_atom_numbers_swapped(_a, a1, done1, a2, done2);
  check_atom_numbers_swapped(_top_front, a1, done1, a2, done2);
  check_atom_numbers_swapped(_top_back, a1, done1, a2, done2);
  check_atom_numbers_swapped(_left_down, a1, done1, a2, done2);
  check_atom_numbers_swapped(_right_down, a1, done1, a2, done2);

#ifdef DEBUG_ATOM_NUMBERS_ARE_SWAPPED
  cerr << "After swapping " << a1 << " and " << a2 << endl;
  debug_print(cerr);
#endif

  return 1;
}

atom_number_t
Chiral_Centre::next_atom (int & i) const
{
  if (0 == i)
  {
    i++;
    if (_top_front >= 0)
      return _top_front;
  }

  if (1 == i)
  {
    i++;
    if (_top_back >= 0)
      return _top_back;
  }

  if (2 == i)
  {
    i++;
    if (_left_down >= 0)
      return _left_down;
  }

  if (3 == i)
  {
    i++;
    if (_right_down >= 0)
      return _right_down;
  }

  return INVALID_ATOM_NUMBER;
}

/*
  Initial implementation from Jibo. Use most of the concepts from his code,
  but just a minimal implementation here. If full functionality needed,
  use his programme.
*/

int
Molecule::_discern_chirality_from_3d_structure(atom_number_t zatom)
{
  const Atom * a = _things[zatom];

  Set_of_Atoms nbr;

  a->connections(zatom, nbr);

  int ncon = nbr.number_elements();

  if (ncon < 3)   // maybe the case of a lone pair and an implicit Hydrogen!!?
  {
    cerr << "Molecule::_discern_chirality_from_3d_structure:cannot process atom with 2 connections, '" << name() << "'\n";
    return 0;
  }

  Chiral_Centre * c = new Chiral_Centre(zatom);

  Coordinates coord[5];   // coords of centre atom and attachments

  coord[0] = *a;

  c->set_top_front(nbr[0]);
  coord[1] = *(_things[nbr[0]]);

  c->set_top_back(nbr[1]);
  coord[2] = *(_things[nbr[1]]);

  c->set_left_down(nbr[2]);
  coord[3] = *(_things[nbr[2]]);

  if (4 == ncon)
  {
    c->set_right_down(nbr[3]);
    coord[4] = *(_things[nbr[3]]);   // never used
  }
  else if (implicit_hydrogens(zatom))
    c->set_right_down(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
  else
    c->set_right_down(CHIRAL_CONNECTION_IS_LONE_PAIR);

// First translate all atoms to the centre

  for (int i = 1; i < 4; i++)
  {
    coord[i] -= coord[0];
  }

  if (coord[1].length() < 1.0e-03 ||
      coord[2].length() < 1.0e-03 ||
      coord[3].length() < 1.0e-03)   // how could this happen?
  {
    cerr << "Zero length bond in chirality determination, '" << name() << "'\n";
    return 0;
  }

  coord[1].cross_product(coord[2]);

  angle_t theta = coord[1].angle_between(coord[3]);

  if (theta > M_PI * 0.5)
    c->invert();

  add_chiral_centre(c);

//cerr << "Angle is " << theta << endl;

  return 1;
}

/*static int
do_print_canonical_ranking(Molecule & m, 
                           std::ostream & os)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    os << i << " " << a->atomic_symbol() << " " << a->ncon() << " connnections, rank " << m.canonical_rank(i) << endl;
  }

  return 1;
}*/

int
Molecule::discern_chirality_from_3d_structure()
{
  if (highest_coordinate_dimensionality() < 3)
  {
    cerr << "discern_chirality_from_3d_structure:cannot process non 3D molecule '" << name() << "'\n";
    return 0;
  }

  remove_all_chiral_centres();

  int rc = 1;

  for (int i = 0; i < _number_elements; i++)
  {
    const Atom * ai = _things[i];

    if (ai->ncon() < 3)
      continue;

    if (6 != ai->atomic_number())
      continue;

     if (ai->ncon() < ai->nbonds())
       continue;

    if (! is_actually_chiral(*this, i))
      continue;
    
    if (! _discern_chirality_from_3d_structure(i))
      rc = 0;
  }

  return rc;
}

void
reset_chiral_centre_file_scope_variables()
{
  _automatically_add_implicit_hydrogen_to_incomplete_chiral_centre = 0;
}
void
Chiral_Centre::new_atom_numbers (const int * xref)
{
  if (_a >= 0 && xref[_a] != _a)
    _a = xref[_a];
  if (_top_back >= 0 && xref[_top_back] != _top_back)
    _top_back = xref[_top_back];
  if (_top_front >= 0 && xref[_top_front] != _top_front)
    _top_front = xref[_top_front];
  if (_left_down >= 0 && xref[_left_down] != _left_down)
    _left_down = xref[_left_down];
  if (_right_down >= 0 && xref[_right_down] != _right_down)
    _right_down = xref[_right_down];

  return;
}
