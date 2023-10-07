#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <iostream>
#include <memory>

#define COMPILING_CTB

#include "Foundational/iwmisc/misc.h"

#include "misc2.h"
#include "molecule.h"
#include "path_scoring.h"

using std::cerr;
using std::endl;

Atomic_Numbers_Encounterd::Atomic_Numbers_Encounterd()
{
  initialise();

  return;
}

static atomic_number_t
atomic_number_from_element(const Element * e)
{
  const IWString s = e->symbol();

  if (! isupper(s[0]))
  {
//  cerr << "Strange element symbol '" << s << "'\n";
    return 26 + 26 * 26 + 1;
  }

  atomic_number_t z = s[0] - 'A';
  if (1 == s.length())
    return z;

  if (! islower(s[1]))
  {
//  cerr << "Strange element symbol '" << s << "'\n";
    return 26 + 26 * 26 + 1;
  }

  z = 26 + 26 * z + s[1] - 'a';

  return z;
}

std::ostream &
operator << (std::ostream & os, const Atomic_Numbers_Encounterd & ane)
{
  os << "Atomic numbers";

  for (int i = 0; i < SIZE_OF_ATOMIC_NUMBERS_ENCOUNTERED_ARRAY; i += 11)
  {
    for (int j = 0; j < 11; j++)
    {
      if (0 == ane._found[i + j])
        continue;

      os << ' ' << (i / 11);

      if (0 == j)
        os << " single-1";
      if (1 == j)
        os << " single-2";
      else if (2 == j)
        os << " single-3";
      else if (3 == j)
        os << " single-4";
      else if (4 == j)
        os << " aromatic-2";
      else if (5 == j)
        os << " aromatic-3";
      else if (6 == j)
        os << " double-1";
      else if (7 == j)
        os << " double-2";
      else if (8 == j)
        os << " double-3";
      else if (9 == j)
        os << " triple-1";
      else if (10 == j)
        os << " triple-2";

//    if (ane._found[i + j] > 1)
        os << '(' << ane._found[i + j] << ')';
    }
  }

  return os;
}

std::ostream &
operator << (std::ostream & os, const Path_Scoring & ps)
{
  ps.debug_print(os);

  return os;
}

void
Atomic_Numbers_Encounterd::initialise()
{
  set_vector(_found, SIZE_OF_ATOMIC_NUMBERS_ENCOUNTERED_ARRAY, 0);

  return;
}

/*int
Atomic_Numbers_Encounterd::extra(const Bond * b, atomic_number_t z)
{
  assert (z >= 0 && z <= HIGHEST_ATOMIC_NUMBER);

  int i;
  if (b->is_aromatic())
    i = 1;
  else if (b->is_single_bond())
    i = 0;
  else if (b->is_double_bond())
    i = 2;
  else if (b->is_triple_bond())
    i = 3;
  else       // Huh?
    i = 0;

  _found[z * 11 + i]++;

  return 1;
}*/

int
Atomic_Numbers_Encounterd::extra(const Element * e)
{
  const atomic_number_t z = atomic_number_from_element(e);

  _found[z * 11]++;

  return 1;
}

/*
  See the header files for the layout within an atomic number grouping
*/

int
Atomic_Numbers_Encounterd::extra(const Bond * b,
                                 const Element * e,
                                 int ncon)
{
  assert (ncon > 0);

  int i;

  if (b->is_aromatic())
  {
    assert (ncon >= 2);

    i = 4 + ncon - 2;
  }
  else if (b->is_single_bond())
  {
    i = ncon - 1;
  }
  else if (b->is_double_bond())
  {
    i = 6 + ncon - 1;
  }
  else if (b->is_triple_bond())
  {
    i = 9 + ncon - 1;
  }
  else    // huh!
  {
    i = ncon - 1;      // pretend it is a single bond
  }

  const atomic_number_t z = atomic_number_from_element(e);

  _found[z * 11 + i]++;

  return 1;
}

int
Atomic_Numbers_Encounterd::operator < (const Atomic_Numbers_Encounterd & rhs) const
{
  for (int i = SIZE_OF_ATOMIC_NUMBERS_ENCOUNTERED_ARRAY - 1; i >= 0; i--)
  {
    if (_found[i] == rhs._found[i])
      continue;
    else if (_found[i] < rhs._found[i])
      return 1;
    else
      return 0;
  }

  return 0;         // they must be the same
}

int
Atomic_Numbers_Encounterd::operator != (const Atomic_Numbers_Encounterd & rhs) const
{
  for (int i = SIZE_OF_ATOMIC_NUMBERS_ENCOUNTERED_ARRAY - 1; i >= 0; i--)
  {
    if (_found[i] != rhs._found[i])
      return 1;       // these two are different
  }

  return 0;      // they must be the same
}

int
Atomic_Numbers_Encounterd::operator == (const Atomic_Numbers_Encounterd & rhs) const
{
  for (int i = SIZE_OF_ATOMIC_NUMBERS_ENCOUNTERED_ARRAY - 1; i >= 0; i--)
  {
    if (_found[i] != rhs._found[i])
      return 0;        // these two are different
  }

  return 1;      // these two are the same
}

/*
  Return -1, 0 or 1 depending in a comparison between the two objects
*/

int
Atomic_Numbers_Encounterd::compare (const Atomic_Numbers_Encounterd & rhs) const
{
  for (int i = SIZE_OF_ATOMIC_NUMBERS_ENCOUNTERED_ARRAY - 1; i >= 0; i--)
  {
//  if (_found[i] != rhs._found[i])
//    cerr << " i = " << i << " lhs " << _found[i] << " and rhs " << rhs._found[i] << endl;

    if (_found[i] < rhs._found[i])
      return 1;
    else if (_found[i] > rhs._found[i])
      return -1;
  }

  return 0;     // they must be identical
}

int
Path_Scoring::compare (const Path_Scoring & rhs) const
{
  if (_number_elements < rhs._number_elements)
    return -1;
  else if (_number_elements > rhs._number_elements)
    return 1;

  if (_first_bond < rhs._first_bond)
    return -1;
  else if (_first_bond > rhs._first_bond)
    return 1;

  for (int i = 0; i < _number_elements; i++)
  {
    int tmp = _things[i]->compare(*(rhs._things[i]));

    if (0 != tmp)
      return tmp;
  }

  return 0;    // must be the same
}

Path_Scoring::Path_Scoring()
{
  _nsteps = 0;

  _first_bond = 0;

  _assigned_rank = -1;

  _edge_atom.resize(20);       // hopefully large enough for most molecules

  return;
}

int
Path_Scoring::debug_print (std::ostream & os) const
{
  os << "Path_Scoring object with " << _number_elements << " steps\n";

  for (int i = 0; i < _number_elements; i++)
  {
    const Atomic_Numbers_Encounterd & ane = *(_things[i]);

    os << ' ' << i << "   " << ane << endl;
  }

  return os.good();
}

/*
  Initialise a path scoring object with atom number A - which should
  be the same as atom AT.
  Note that we don't make any allowance for multiple bonds along this
  first path. Maybe we should - here would be the place to do it
*/

int
Path_Scoring::initialise (atom_number_t a, const Atom * at)
{
  resize_keep_storage(0);

  _start_atom = a;

  if (0 == _elements_allocated)
    resize(10);

// We can only advance beyond atom AT if it has more than 1 connection

  if (at->ncon() > 1)
  {
    if (0 == _edge_atom.elements_allocated())
      _edge_atom.resize(10);

    _edge_atom.add(a);
  }
  else     // just one connection, we cannot advance
    _edge_atom.resize_keep_storage(0);

  Atomic_Numbers_Encounterd * ane = new Atomic_Numbers_Encounterd;
  ane->extra(at->element());

  add(ane);
  
  return 1;
}

int
Path_Scoring::operator == (const Path_Scoring & rhs) const
{
  if (_nsteps != rhs._nsteps)
    return 0;

  if (_number_elements != rhs._number_elements)
    return 0;

  for (int i = 0; i < _number_elements; i++)
  {
    if (*(_things[i]) == *(rhs._things[i]))
      continue;

    return 0;
  }
  
  return 1;
}

int
Path_Scoring::operator != (const Path_Scoring & rhs) const
{
  if (_number_elements != rhs._number_elements)
    return 1;              // these are in fact not equal

  for (int i = 0; i < _number_elements; i++)
  {
    if (*(_things[i]) != *(rhs._things[i]))
      return 1;
  }

  return 0;     // These must be the same
}

//#define DEBUG_ADVANCE

/*
  If we are on the first step, initialise will have put an Atomic_Numbers_Encounterd
  object in our array. Otherwise we need to allocate a new one for this step
*/

int
Path_Scoring::advance (Atom * const * atom,
                       const int * claimed)
{
  int ne = _edge_atom.number_elements();
  if (0 == ne)
    return 1;

  Atomic_Numbers_Encounterd * ane = new Atomic_Numbers_Encounterd;

  add(ane);

  Set_of_Atoms new_edge_atoms;    // easier than doing an in-place replacement into _edge_atom

  if (new_edge_atoms.elements_allocated() < ne * 2 + 5)
    new_edge_atoms.resize(ne * 2 + 5);

#ifdef DEBUG_ADVANCE
  cerr << "Advancing " << ne << " edge atoms, nsteps = " << _nsteps << endl;
#endif

  for (int i = 0; i < ne; i++)
  {
    int j = _edge_atom[i];

#ifdef DEBUG_ADVANCE
    cerr << "Atom " << j << " is an edge atom\n";
#endif

    assert (1 == claimed[j]);

    const Atom * a = atom[j];

    int acon = a->ncon();

    for (int k = 0; k < acon; k++)
    {
      const Bond * b = a->item(k);

      atom_number_t l = b->other(j);

#ifdef DEBUG_ADVANCE
      cerr << "Atom " << l << " is attached, claimed = " << claimed[l] << endl;
#endif

      if (1 == claimed[l])      // someone else already got this one
        continue;

//    Add the extra atom here. We may have two paths converging on the same atom,
//    but one doing it via a single bond, the other via a double bond. 

      ane->extra(b, atom[l]->element(), atom[l]->ncon());

      if (new_edge_atoms.contains(l))     // we already got this off some other atom
        continue;

      if (atom[l]->ncon() > 1)    // if atom L is terminal, we cannot advance past it
        new_edge_atoms.add(l);

#ifdef DEBUG_ADVANCE
      if (_nsteps > 0)
        cerr << "Added atomic number " << atom[l]->atomic_number() << endl;
#endif
    }
  }

  _edge_atom.resize_keep_storage(0);

  if (new_edge_atoms.empty())
    return 0;

  _nsteps++;

  _edge_atom = new_edge_atoms;

  return _edge_atom.number_elements();
}

int
Path_Scoring::update_claimed(int * claimed) const
{
  _edge_atom.set_vector(claimed, 1);

  return 1;
}

/*
  This is called when there is a tolerance on E/Z determinations.
  The high priority atoms have been identified and we have the path

  a       d            a                     a
   \     /              \                     \
    b = c       OR       b = c         OR      b = c - d
                              \
                               d

  We check to see if any of the angles are within TOLERANCE.
*/ 

static int
straight_bond (const Coordinates & ab,
               const Coordinates & bc,
               const Coordinates & cd,
               const angle_t tolerance)
{
  angle_t theta = ab.angle_between_unit_vectors(bc);

//cerr << "Theta1 " << theta << endl;

  if (fabs(theta) <= tolerance)
    return 1;

  theta = bc.angle_between_unit_vectors(cd);

//cerr << "Theta2 " << theta << endl;

  if (fabs(theta) < tolerance)
    return 1;

  return 0;
}

/*
  Look for singly bonded atoms to atom A.
*/

static int
identify_attached_atoms (const Atom * a,
                         atom_number_t zatom,
                         atom_number_t & a1,
                         atom_number_t & a2)
{
  a1 = INVALID_ATOM_NUMBER;
  a2 = INVALID_ATOM_NUMBER;

  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_single_bond())
      continue;

    atom_number_t j = b->other(zatom);

    if (INVALID_ATOM_NUMBER == a1)
      a1 = j;
    else
    {
      a2 = j;
      return 1;
    }
  }

  return (INVALID_ATOM_NUMBER != a1);
}

/*
  We have identified a bond as possibly being part of a cis-trans grouping.
  We need to perceive the directionality associated with that bond.
*/

/*static void
set_directionality (const Bond & b,
                    atom_number_t a,
                    int & direction)
{
  if (! b.is_directional())
  {
    direction = 0;
    return;
  }

  if (a == b.a1())
  {
    if (b.is_directional_up())
      direction = 1;
    else if (b.is_directional_down())
      direction = -1;
  }
  else if (a == b.a2())
  {
    if (b.is_directional_up())
      direction = -1;
    else if (b.is_directional_down())
      direction = 1;
  }

  return;
}*/

/*
  Atom A is at one end of a cis-trans bond. Identify the single bonds
  that are attached
*/

static int
identify_attached_bonds (const Atom * a,
                         Bond * & b12,
                         Bond * & b13)
{
  b12 = nullptr;
  b13 = nullptr;     // may not exist

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_double_bond())
      continue;

    if (nullptr == b12)
      b12 = const_cast<Bond *>(b);
    else if (nullptr == b13)
    {
      b13 = const_cast<Bond *>(b);
      return 2;
    }
  }

  return (nullptr != b12);
}

/*static int
identify_attached_atoms (const Atom * a,
                         atom_number_t zatom,
                         atom_number_t & a1,
                         int & a1_direction,
                         atom_number_t & a2,
                         int & a2_direction)
{
  a1 = INVALID_ATOM_NUMBER;
  a2 = INVALID_ATOM_NUMBER;

  a2_direction = 0;    // a1_direction will always be set

  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_single_bond())
      continue;

    atom_number_t j = b->other(zatom);

    if (INVALID_ATOM_NUMBER == a1)
    {
      a1 = j;
      set_directionality(*b, zatom, a1_direction);
    }
    else
    {
      a2 = j;
      set_directionality(*b, zatom, a2_direction);
      return 1;
    }
  }

  return (INVALID_ATOM_NUMBER != a1);
}*/

/*
  Someone is interested in an E/Z classification.
  We must first identify the atoms involved

  a1        a5
    \      /
     a3==a4
    /      \
   a2       a6

  A3 and A4 define the double bond. ALHS and ARHS will be set to
  the two high priority atoms connected to each end of the double
  bond.

  Note that some of these may be missing. We return 0 if we either
  cannot identify or cannot resolve ALHS and ARHS
*/

int
Molecule::identify_ez_atoms (atom_number_t a3, atom_number_t a4,
                             atom_number_t & alhs, atom_number_t & arhs) const
{
  assert (ok_2_atoms(a3, a4));

  alhs = INVALID_ATOM_NUMBER;
  arhs = INVALID_ATOM_NUMBER;

  const Atom * aa3 = _things[a3];
  const Atom * aa4 = _things[a4];

  int a3con = aa3->ncon();
  int a4con = aa4->ncon();

  if (a3con < 2 || a3con > 3 || a4con < 2 || a4con > 3)
  {
//  cerr << "Molecule::score_ez_bond: atom " << a3 << " (" << a3con << " connections) or atom " << a4 << " (" << a4con << " connections) cannot form a E/Z bond\n";
    return 0;
  }

  bond_type_t bt;
  if (! aa3->is_bonded_to(a4, bt) || ! IS_DOUBLE_BOND(bt))
  {
    cerr << "Molecule::score_ez_bond: atoms " << a3 << " and " << a4 << " not bonded via a double bond\n";
    return 0;
  }

// Identify the atoms attached

  atom_number_t a1, a2;
  if (! identify_attached_atoms(aa3, a3, a1, a2))
    return 0;

  assert (INVALID_ATOM_NUMBER != a1);

  atom_number_t a5, a6;
  if (! identify_attached_atoms(aa4, a4, a5, a6))
    return 0;

  assert (INVALID_ATOM_NUMBER != a5);

  arhs = alhs = INVALID_ATOM_NUMBER;

// Handle the cases of just one bond off the end

  if (INVALID_ATOM_NUMBER == a2)
    alhs = a1;
  if (INVALID_ATOM_NUMBER == a6)
    arhs = a5;

// If both atoms are already known, we must have a double bond that has
// just one substituent at each end *-*=*-*

  if (INVALID_ATOM_NUMBER != alhs && INVALID_ATOM_NUMBER != arhs)
    return 1;

  int * already_done = new_int(_number_elements); std::unique_ptr<int[]> free_already_done(already_done);

  already_done[a1] = 1;
  if (INVALID_ATOM_NUMBER != a2)
    already_done[a2] = 1;

  already_done[a3] = 1;
  already_done[a4] = 1;

  already_done[a5] = 1;
  if (INVALID_ATOM_NUMBER != a6)
    already_done[a6] = 1;

//#define DEBUG_SCORE_EZ

  int rc = 1;

  if (INVALID_ATOM_NUMBER == alhs)     // 2 substituents, need to resolve
  {
    Path_Scoring ps[2];
    ps[0].initialise(a1, _things[a1]);
    ps[1].initialise(a2, _things[a2]);

    int lhs = _score_ez_bond(already_done, ps[0], ps[1]);

#ifdef DEBUG_SCORE_EZ
  cerr << "After scoring atom " << a3 << " (atoms " << a1 << " and " << a2 << ") lhs " << lhs << endl;
#endif

    if (0 == lhs)      // could not be resolved, double bond cannot be E/Z
      rc = 0;
    else if (lhs < 0)
      alhs = a1;
    else
      alhs = a2;
  }

  if (0 == rc)    // lhs not resolved, no point scoring the other end
    ;
  else if (INVALID_ATOM_NUMBER == arhs)
  {
    Path_Scoring ps[2];

    ps[0].initialise(a5, _things[a5]);
    ps[1].initialise(a6, _things[a6]);

    int rhs = _score_ez_bond(already_done, ps[0], ps[1]);

#ifdef DEBUG_SCORE_EZ
  cerr << "After scoring atom " << a4 << " (atoms " << a5 << " and " << a6 << ") rhs " << rhs << endl;
#endif

    if (0 == rhs)     // rhs could not be resolved
      rc = 0;
    else if (rhs < 0)
      arhs = a5;
    else
      arhs = a6;
  }

  return rc;
}

//#define DEBUG_EZ_GEOMETRY

/*
  What is the geometric arrangement of A1 and A4 across the line
  defined by A2-A3
*/

int
Molecule::ez_by_geometry (atom_number_t a1,
                          atom_number_t a2,
                          atom_number_t a3,
                          atom_number_t a4,
                          angle_t tolerance) const
{
#ifdef DEBUG_EZ_GEOMETRY
  cerr << "Checking geometry for atoms " << a1 << ", " << a2 << ", " << a3 << " and " << a4 << endl;
#endif

  const Atom * aa1 = _things[a1];
  const Atom * aa2 = _things[a2];
  const Atom * aa3 = _things[a3];
  const Atom * aa4 = _things[a4];

//assert (0.0 == aa1->z() && 0.0 == aa2->z() && 0.0 == aa3->z() && 0.0 == aa4->z());

  Coordinates ab = *aa1 - *aa2;
  Coordinates bc = *aa2 - *aa3;
  Coordinates cd = *aa3 - *aa4;

  if (static_cast<coord_t>(0.0) == ab.x() && static_cast<coord_t>(0.0) == ab.y() && static_cast<coord_t>(0.0) == ab.z())
  {
    cerr << "Molecule::ez_by_geometry:zero distance between atoms, '" << _molecule_name << "'\n";
    abort();
    return 0;
  }

#ifdef DEBUG_EZ_GEOMETRY
  cerr << "AB " << ab << endl;
  cerr << "BC " << bc << endl;
  cerr << "CD " << cd << endl;
#endif

  ab.normalise();
  bc.normalise();
  cd.normalise();

  if (tolerance > static_cast<angle_t>(0.0))
  {
    if (straight_bond(ab, bc, cd, tolerance))
      return 0;
  }

  Coordinates x1 = bc;
  x1.cross_product(ab);

  Coordinates x2 = bc;
  x2.cross_product(cd);

  angle_t a = x1.angle_between(x2);

#ifdef DEBUG_EZ_GEOMETRY
  cerr << "Checking result " << x1 << " and " << x2 << " angle " << (a * RAD2DEG) << endl;
#endif

  if (fabs(a) <= (M_PI * 0.5))    // pointing in roughly the same direction
    return -1;
  else
    return 1;

#ifdef OLD_TWO_DIMENSIONAL_CODE
  if (x1.z() < 0.0 && x2.z() < 0.0)
  {
    return -1;
  }
  else if (x1.z() > 0.0 && x2.z() > 0.0)
  {
    return -1;
  }
  else if (x1.z() > 0.0 && x2.z() < 0.0)
  {
    return 1;
  }
  else if (x1.z() < 0.0 && x2.z() > 0.0)
  {
    return 1;
  }
  else
  {
    cerr << "Molecule::::_ez_by_geometry: atoms " << a1 << ", " << a2 << ", " << a3 << ", " << a4 << " possibly collinear\n";
    return 0;
  }
#endif
}

/*
  Score TWO path scoring objects

  Return -1 if the first one is dominant,
  return  1 if the 2nd   one is dominant
  return  0 if they are equivalent
*/

int
Molecule::_score_ez_bond (int * already_done,
                          Path_Scoring & ps1, Path_Scoring & ps2) const
{
  while (1)
  {
    int compare = ps1.compare(ps2);

//#define DEBUG_SCORE_EZ_BOND
#ifdef DEBUG_SCORE_EZ_BOND
    cerr << "Compare = " << compare << endl;
#endif

    if (compare < 0)
      return -1;
    else if (compare > 0)
      return 1;

#ifdef DEBUG_SCORE_EZ_BOND
    cerr << "Active " << ps1.active() << " and " << ps2.active() << endl;
#endif

    if (ps1.active() && ps2.active())    // not yet resolved, need to keep going
      ;
    else if (! ps1.active() && ! ps2.active())    // both done, not resolved
      return 0;
    else if (! ps1.active())
      return 1;
    else if (! ps2.active())
      return -1;

    ps1.advance(_things, already_done);
    ps2.advance(_things, already_done);

    ps1.update_claimed(already_done);
    ps2.update_claimed(already_done);
  }

  return 0;     // they cannot be resolved
}


/*
  The primary comparison is on the number of steps taken
*/

static int
path_scoring_comparitor (Path_Scoring * const * pps1, Path_Scoring * const * pps2)
{
  Path_Scoring * ps1 = *pps1;
  Path_Scoring * ps2 = *pps2;

  int n1 = ps1->nsteps();
  int n2 = ps2->nsteps();

  if (n1 < n2)        // ps1 has taken fewer steps
    return 1;
  else if (n1 > n2)   // ps2 has taken fewer steps
    return -1;
  else                // they have both taken the same number of steps
    return ps1->compare(*ps2);
}

//#define DEBUG_RESOLVED

/*
  For a set to be resolved, all items must be different.
  Sort, then look for identical neighbours

  If we come to a pair of objects which are equal, but both inactive,
  then we set the STOPPED variable
*/

int
resolved(resizable_array_p<Path_Scoring> & ps,
         int & stopped)
{
  stopped = 0;

  const int np = ps.number_elements();

#ifdef DEBUG_RESOLVED
  cerr << "Are " << np << " path scoring objects resolved\n";
#endif

  if (1 == np)    // hard to imagine this happening
    return 1;

  if (2 == np)
  {
    if (*(ps[0]) != *(ps[1]))    // are different, resolved
      return 1;

    if (! ps[0]->active() || ! ps[1]->active())
      stopped = 1;

    return 0;
  }

  ps.sort(path_scoring_comparitor);

  for (int i = 1; i < np; i++)
  {
    const Path_Scoring & prev = *(ps[i - 1]);
    const Path_Scoring & curr = *(ps[i]);

#ifdef DEBUG_RESOLVED
    cerr << "Comparing " << prev << endl;
    cerr << "     with " << curr << endl;
#endif

    if (curr != prev)     // these two are resolved
      continue;

#ifdef DEBUG_RESOLVED
    cerr << "They are equal\n";
#endif

//  Two adjacent items are the same. If they are both inactive, we are stopped

    if (! curr.active() && ! prev.active())
      stopped = 1;
  
    return 0;     // items not resolved
  }

  return 1;    // yes, we are resolved
}

/*
  We have two atoms at the end of a bond. Set a directionality indicator if needed
*/

static int
discern_directionality (atom_number_t a1, 
                        atom_number_t a2,
                        const Bond & b)
{
  if (INVALID_ATOM_NUMBER == a2)
    return 0;

  if (! b.is_directional())
    return 0;

  if (b.is_directional_up())
  {
    if (b.a1() == a1)
      return 1;
    else
      return -1;
  }
  else
  {
    if (b.a1() == a1)
      return -1;
    else
      return 1;
  }
}

//#define DEBUG_DISCERN_CIS_TRANS_BOND_FROM_DEPICTION

/*
  This was written for use when reading mdl files

  a1        a5
    \      /
     a3==a4
    /      \
  a2        a6

  We only return 0 if there is an error
*/

int
Molecule::_discern_cis_trans_bond_from_depiction (Bond * b)
{
#ifdef  DEBUG_DISCERN_CIS_TRANS_BOND_FROM_DEPICTION
  cerr << "Looking for cis-trans bond involving " << (*b) << endl;
#endif

  if (b->part_of_cis_trans_grouping())   // already done
    return 0;

  atom_number_t a3 = b->a1();
  atom_number_t a4 = b->a2();

  if (in_same_ring(a3, a4))
    return 1;

  atom_number_t alhs, arhs;    // not used

  if (! identify_ez_atoms(a3, a4, alhs, arhs))      // redundancy between this call and the two calls to identify_attached_atoms
    return 1;

  Bond * b31;
  Bond * b32;

  if (! identify_attached_bonds(_things[a3], b31, b32))
    return 1;

  Bond * b45;
  Bond * b46;

  if (! identify_attached_bonds(_things[a4], b45, b46))
    return 0;

  atom_number_t a1 = b31->other(a3);
  int a1_direction = discern_directionality(a3, a1, *b31);

  atom_number_t a2;
  int a2_direction;
  if (nullptr == b32)
  {
    a2 = INVALID_ATOM_NUMBER;
    a2_direction = 0;
  }
  else 
  {
    a2 = b32->other(a3);
    a2_direction = discern_directionality(a3, a2, *b32);
  }

  atom_number_t a5 = b45->other(a4);
  int a5_direction = discern_directionality(a4, a5, *b45);

  atom_number_t a6;
  int a6_direction;
  if (nullptr == b46)
  {
    a6 = INVALID_ATOM_NUMBER;
    a6_direction = 0;
  }
  else
  {
    a6 = b46->other(a4);
    a6_direction = discern_directionality(a4, a6, *b46);
  }

// Look for any inconsistent directionality specification

  if (INVALID_ATOM_NUMBER == a2)     // just one connection, nothing to check
    ;
  else if (0 == a1_direction && 0 == a2_direction)   // neither bond directional, great
    ;
  else if (a1_direction == a2_direction)
  {
    cerr << "Molecule::_discern_cis_trans_bond_from_depiction:inconsistent directionality\n";
    cerr << "Atoms " << a1 << " (" << a1_direction << ") - " << a3 << " - " << a2 << " (" <<  a2_direction << ")\n";
    return 0;
  }

  if (INVALID_ATOM_NUMBER == a6)            // just one connection, nothing to check
    ;
  else if (0 == a5_direction && 0 == a6_direction)   // neither bond directional, great
    ;
  else if (a5_direction == a6_direction)
  {
    cerr << "Molecule::_discern_cis_trans_bond_from_depiction:inconsistent directionality\n";
    cerr << "Atoms " << a5 << " (" << a5_direction << ") - " << a4 << " - " << a6 << " (" << a6_direction << ")\n";
    return 0;
  }

// Now we need to discern the relative orientation around the bond

  int ez = ez_by_geometry(a1, a3, a4, a5, static_cast<angle_t>(10.0 * DEG2RAD));

#ifdef DEBUG_DISCERN_CIS_TRANS_BOND_FROM_DEPICTION
  cerr << "Atoms " << a1 << "-" << a3 << "=" << a4 << "-" << a5 << " ez = " << ez << endl;
#endif

  if (0 == ez)
    return 1;

// Remember that ez_by_geometry returns 1 for Z (same side) and -1 for E (opposite side)

  if (ez < 0)      // E configuration, A1 and A5 on opposite sides. Switch A5 and A6
  {
    std::swap(a5, a6);
    std::swap(a5_direction, a6_direction);
    Bond * btmp = b46;
    b46 = b45;
    b45 = btmp;
  }

#ifdef DEBUG_DISCERN_CIS_TRANS_BOND_FROM_DEPICTION
  cerr << "Directions:\n";
  cerr << "a1 " << a1 << ' ' << a1_direction << " a2 " << a2 << ' ' << a2_direction << " a5 " << a5 << ' ' << a5_direction << " a6 " << a6 << ' ' << a6_direction << endl;
#endif

// Check the various pairs of atoms across the double bond (1,5), (2,5), (2,6) and (1,6)

  if (0 == a1_direction || 0 == a5_direction)
    ;
  else if (a1_direction != a5_direction)
  {
    cerr << "Molecule::_discern_cis_trans_bond_from_depiction:inconsistent 1,5 directionality across bond\n";
    cerr << "Atoms " << a1 << " (" << a1_direction << ") - " << a3 << " = " << a4 << " - " << a5 << " (" << a5_direction << ") molecule '" << _molecule_name << "'\n";
    return 0;
  }

  if (INVALID_ATOM_NUMBER == a2 || 0 == a5_direction || 0 == a2_direction)
    ;
  else if (a2_direction == a5_direction)
  {
    cerr << "Molecule::_discern_cis_trans_bond_from_depiction:inconsistent 2,5 directionality across bond\n";
    cerr << "Atoms " << a2 << " (" << a2_direction << ") - " << a3 << " = " << a4 << " - " << a5 << " (" << a5_direction << ") molecule '" << _molecule_name << "'\n";
    return 0;
  }

  if (INVALID_ATOM_NUMBER == a2 || INVALID_ATOM_NUMBER == a6 || 0 == a2_direction || 0 == a6_direction)
    ;
  else if (a2_direction != a6_direction)
  {
    cerr << "Molecule::_discern_cis_trans_bond_from_depiction:inconsistent directionality across bond\n";
    cerr << "Atoms " << a2 << " - " << a3 << " = " << a4 << " - " << a6 << " molecule '" << _molecule_name << "'\n";
    return 0;
  }

  if (INVALID_ATOM_NUMBER == a6 || 0 == a1_direction || 0 == a6_direction)
    ;
  else if (a6_direction == a1_direction)
  {
    cerr << "Molecule::_discern_cis_trans_bond_from_depiction:inconsistent directionality across bond\n";
    cerr << "Atoms " << a1 << " - " << a3 << " = " << a4 << " - " << a6 << " molecule '" << _molecule_name << "'\n";
    return 0;
  }

// Now our checking is done, start setting things

  if (0 == a1_direction && 0 == a2_direction && 0 == a5_direction && 0 == a6_direction)   // hopefully the most common case
  {
    b31->set_directional_up(a3, a1);
    if (nullptr != b32)
      b32->set_directional_down(a3, a2);
    if (nullptr != b45)
      b45->set_directional_up(a4, a5);
    if (nullptr != b46)
      b46->set_directional_down(a4, a6);
  }
  else     // we assume that just one of the directions is set
  {

#define D56_FIVE_ON_TOP 1
#define D56_SIX_ON_TOP -1
#define D12_ONE_ON_TOP 1
#define D12_TWO_ON_TOP -1

    int d12 = 0;
    int d56 = 0;

    if (a1_direction < 0)
    {
      d12 = D12_TWO_ON_TOP;
      d56 = D56_SIX_ON_TOP;
    }
    else if (a1_direction > 0)
    {
      d12 = D12_ONE_ON_TOP;
      d56 = D56_FIVE_ON_TOP;
    }
    else if (a2_direction < 0)
    {
      d12 = D12_ONE_ON_TOP;
      d56 = D56_FIVE_ON_TOP;
    }
    else if (a2_direction > 0)
    {
      d12 = D12_TWO_ON_TOP;
      d56 = D56_SIX_ON_TOP;
    }
    else if (a5_direction < 0)
    {
      d12 = D12_TWO_ON_TOP;
      d56 = D56_SIX_ON_TOP;
    }
    else if (a5_direction > 0)
    {
      d12 = D12_ONE_ON_TOP;
      d56 = D56_FIVE_ON_TOP;
    }
    else if (a6_direction < 0)
    {
      d12 = D12_ONE_ON_TOP;
      d56 = D56_FIVE_ON_TOP;
    }
    else if (a6_direction > 0)
    {
      d12 = D12_TWO_ON_TOP;
      d56 = D56_SIX_ON_TOP;
    }

    if (D12_ONE_ON_TOP == d12)
    {
      b31->set_directional_up(a3, a1);
      if (nullptr != b32)
        b32->set_directional_down(a3, a2);
    }
    else
    {
      b31->set_directional_down(a3, a1);
      if (nullptr != b32)
        b32->set_directional_up(a3, a2);
    }

    if (D56_FIVE_ON_TOP == d56)
    {
      if (nullptr != b45)
        b45->set_directional_up(a4, a5);
      if (nullptr != b46)
        b46->set_directional_down(a4, a6);
    }
    else
    {
      if (nullptr != b45)
        b45->set_directional_down(a4, a5);
      if (nullptr != b46)
        b46->set_directional_up(a4, a6);
    }
  }

/*
  WRONG. Did not account for any pre-existing direction if it showed up in a5 or a6

  if (0 == a1_direction && 0 == a2_direction)
  {
    b31->set_directional_up (a3, a1);
    if (nullptr != b32)
      b32->set_directional_down (a3, a2);
  }
  else if (0 != a2_direction)    
  {
    if (a2_direction < 0)
      b31->set_directional_up (a3, a1);
    else
      b31->set_directional_down (a3, a1);
  }
  else if (nullptr != b32)
  {
    if (a1_direction < 0)
      b32->set_directional_up (a3, a2);
    else
      b32->set_directional_down (a3, a2);
  }

  if (0 == a5_direction && 0 == a6_direction)
  {
    if (nullptr != b45)
      b45->set_directional_up (a4, a5);
    if (nullptr != b46)
      b46->set_directional_down (a4, a6);
  }
  else if (0 != a6_direction)    
  {
    if (nullptr == b45)
      ;
    else if (a6_direction < 0)
      b45->set_directional_up (a4, a5);
    else
      b45->set_directional_down (a4, a5);
  }
  else if (nullptr != b46)
  {
    if (a5_direction < 0)
      b46->set_directional_up (a4, a6);
    else
      b46->set_directional_down (a4, a6);
  }*/

  b->set_part_of_cis_trans_grouping(1);

// Are there any double bonds joined to any of the atoms attached

  if (! _extend_cis_trans_system(a1))
    return 0;

  if (INVALID_ATOM_NUMBER != a2)
  {
    if (! _extend_cis_trans_system(a2))
      return 0;
  }

  if (INVALID_ATOM_NUMBER != a5)
  {
    if (! _extend_cis_trans_system(a5))
      return 0;
  }

  if (INVALID_ATOM_NUMBER != a6)
  {
    if (! _extend_cis_trans_system(a6))
      return 0;
  }

  return 1;
}

/*
  Atom ZATOM is attached to (but not part of) a double bond that has cis-trans stereochemistry.
  Is it part of another cis-trans grouping?
  We only return 0 if there is an error
*/

int
Molecule::_extend_cis_trans_system (atom_number_t zatom)
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  if (acon == a->nbonds())
    return 1;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    if (b->part_of_cis_trans_grouping())   // already done
      continue;

    if (b->is_cis_trans_either_double_bond())
      continue;

    return _discern_cis_trans_bond_from_depiction(const_cast<Bond *>(b));
  }

  return 1;
}

int
Molecule::_discern_cis_trans_bond_from_depiction (atom_number_t zatom)
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  if (acon == a->nbonds())
    return 1;

  int rc = 1;

  for (int i = 0; i < acon; i++)
  {
    Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    if (b->is_cis_trans_either_double_bond())
      continue;

    if (! _discern_cis_trans_bond_from_depiction(b))
      rc = 0;
  }

  return rc;
}

int
assign_ranks (resizable_array_p<Path_Scoring> & ps)
{
  ps.sort(path_scoring_comparitor);

  int np = ps.number_elements();

  int rank_to_assign = 1;

  ps[0]->set_assigned_rank(rank_to_assign);

  for (int i = 1; i < np; i++)
  {
    const Path_Scoring * prev = ps[i - 1];
    Path_Scoring * pi = ps[i];

    int comparison = pi->compare(*prev);

    if (0 != comparison)
      rank_to_assign++;

    pi->set_assigned_rank(rank_to_assign);
  }

  return 1;
}

int
Molecule::remove_invalid_directional_bonds()
{
  int rc = 0;    // number of invalid cis-trans bonds we remove

  int * already_done = nullptr;

  int nb = _bond_list.number_elements();

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];

    if (! b->is_double_bond())
      continue;

    if (! b->part_of_cis_trans_grouping())
      continue;

    atom_number_t a3 = b->a1();
    atom_number_t a4 = b->a2();

    int remove_it = 0;

    if (1 == _things[a3]->ncon())    // definitely bad
      remove_it = 1;
    if (2 == _things[a3]->ncon())
      ;
    else
    {
      if (nullptr == already_done)
        already_done = new_int(_number_elements);

      atom_number_t a1, a2;
      identify_attached_atoms(_things[a3], a3, a1, a2);

      Path_Scoring ps[2];
      ps[0].initialise(a1, _things[a1]);
      ps[1].initialise(a2, _things[a2]);
      already_done[a1] = already_done[a2] = 1;

      if (0 == _score_ez_bond(already_done, ps[0], ps[1]))
        remove_it = 1;
    }

    if (remove_it)
      ;
    else if (1 == _things[a4]->ncon())
      remove_it = 1;
    else if (2 == _things[a4]->ncon())
      ;
    else
    {
      atom_number_t a5, a6;
      identify_attached_atoms(_things[a4], a4, a5, a6);

      Path_Scoring ps[2];
      ps[0].initialise(a5, _things[a5]);
      ps[1].initialise(a6, _things[a6]);

      if (nullptr == already_done)
        already_done = new_int(_number_elements);

      already_done[a5] = already_done[a6] = 1;

      if (0 == _score_ez_bond(already_done, ps[0], ps[1]))
        remove_it = 1;
    }

//  cerr << "Bond between " << a3 << " and " << a4 << " remove_it " << remove_it << endl;

    if (! remove_it)
      continue;

    const_cast<Bond *>(b)->set_part_of_cis_trans_grouping(0);

//  Any single bonds attached that are not part of another cis-trans bond must be notified

    _cis_trans_bond_has_been_invalidated(a3);
    _cis_trans_bond_has_been_invalidated(a4);
    rc++;
  }

  if (nullptr == already_done)
    delete [] already_done;

  return rc;
}

/*
  A cis-trans bond involving atom zatom has been invalidated
  Notify any atoms singly bonded to ZATOM
*/

int
Molecule::_cis_trans_bond_has_been_invalidated (atom_number_t zatom)
{
  int rc = 0;

  const Atom * a = _things[zatom];

  for (int i = 0; i < a->ncon(); i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_single_bond())
      continue;

    atom_number_t j = b->other(zatom);

    const Atom * aj = _things[j];

    int jcon = aj->ncon();

    int found_another_cis_trans_bond = 0;

    for (int k = 0; k < jcon; k++)
    {
      const Bond * b = aj->item(k);

      if (! b->is_double_bond())
        continue;

      if (b->part_of_cis_trans_grouping())
      {
        found_another_cis_trans_bond = 1;
        break;
      }
    }

    if (! found_another_cis_trans_bond)
    {
      const_cast<Bond *>(b)->set_not_directional();
      rc = 1;
    }
  }

  return rc;
}
