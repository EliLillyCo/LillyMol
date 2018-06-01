#include <stdlib.h>
#include <assert.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "misc2.h"

#include "molecule.h"

/*
  NOTE THIS IS OBSOLETE, AND NO LONGER USED

  Scoring a centre for chirality is somewhat involved.  A (possibly)
  chiral centre has a number of paths from it which must be explored. 
  Each path must know its "value" by whatever means that is computed. 

  It must also know which atoms it has covered, array _counted.  It
  also keeps track of the atoms it is processing in the _atom_list
  array.  Atoms are added to the end of the _atom_list as they are
  found.  Contrast this with the _counted array in which the entries
  are boolean, ordered by atom number.

  When I figure it out, we will have a pointer function which returns
  a value for an element.
*/

class Path_Scoring_Object
{
  friend
    operator << (ostream &, const Path_Scoring_Object &);

  private:
    int * _atom_list;
    int * _counted;
    int _ncounted;
    int _start_of_last_generation;
    int _score;
    int _bonds_counted;
    int _stopped;

//  This function returns a value for _score, based upon an element type

    int (* _evaluate_score) (const Element *);

  public:
    Path_Scoring_Object (Molecule *, atom_number_t, atom_number_t);
    ~Path_Scoring_Object ();

    int ok () const;
    int debug_print (ostream &) const;

    int score () const { return _score;}
    int stopped () const { return _stopped;}

    int advance_to_next_generation (Molecule *, const int *);
    int update_global_counted_array (int, int *);
};

// Root atom is the atom being evaluated

Path_Scoring_Object::Path_Scoring_Object (Molecule * m, atom_number_t root_atom,
                                          atom_number_t start_atom)
{
  assert (m->are_bonded (root_atom, start_atom));

  _atom_list = new_int (m->natoms ());
  _counted   = new_int (m->natoms ());

  _atom_list[0] = root_atom;
  _atom_list[1] = start_atom;
  _ncounted = 2;
  _counted[root_atom] = 1;
  _counted[start_atom] = 1;

  _start_of_last_generation = 1;

  _bonds_counted = 1;

  _score = m->atomic_number (start_atom);

  _stopped = 0;

//_evaluate_score = evaluate_score;
//_score = (m->elementi (start_atom))->_evaluate_score ();

  return;
}

Path_Scoring_Object::~Path_Scoring_Object ()
{
  assert (ok ());

  if (-98 == _start_of_last_generation)
    cerr << "Freeing already deleted Path_Scoring_Object\n";

  _start_of_last_generation = -98;
  _ncounted = -99;

  delete _atom_list;
  delete _counted;

  _score = -9;

  _stopped = 1;

 return;
}

int
Path_Scoring_Object::ok () const
{
  if (NULL == _atom_list)
    return 0;

  if (NULL == _counted)
    return 0;

// When the object is initialised, we always set two members as known.

  if (_ncounted < 2)
    return 0;

  if (_start_of_last_generation < 0)
    return 0;

  if (_start_of_last_generation >= _ncounted && ! _stopped)
    return 0;

  return 1;
}

int
Path_Scoring_Object::debug_print (ostream & os) const
{
  assert (os.good ());

  if (! ok ())
    os << "Path_Scoring_Object:: warning, ok fails\n";

  os << "Path_Scoring_Object:: included " << _ncounted << " atoms, last generation starts at " << 
        _start_of_last_generation << ", score = " << _score << endl;
  os << _bonds_counted << " bonds counted\n";

  os << "Atoms";
  for (int i = 0; i < _ncounted; i++)
    os << " " << _atom_list[i];

  os << endl;

  return 1;
}

/*
  Advance to the next generation.
  Array COUNTED will come from the Atom_Scoring_Object invoking us. We ignore any atom
  for which _counted is set.

  Note that we do not update the counted array. Function update_global_counted_array is
  called to do that. Consider the case of a five membered ring.
    1    2
     ____
 0  /    |
    \____|
    3    4

  and starting with atom 0. the paths would advance to atoms 1 and 3. Everything is fine.
  Then, the upper path advances to atom 2. IF we updated the counted array, the lower
  path could not then advance, and would stop at atom 3, and record the paths as
  different (which they are not). 

  The array indexing here is a little tricky.
*/

int
Path_Scoring_Object::advance_to_next_generation (Molecule * m, const int * counted)
{
  assert (ok ());
  assert (NULL != counted);

//cerr << "Path scoring object advancing, start of last = " << _start_of_last_generation << " counted = " << _ncounted << endl;

  if (_stopped)
    return 0;

  int istop = _ncounted;
  int rc = 0;
  for (int i = _start_of_last_generation; i < istop; i++)
  {
    atom_number_t a = _atom_list[i];
    _score += m->implicit_hydrogens (a);    // don't forget our implicit hudrogens.

    int acon = m->ncon (a);
    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = m->other (a, j);
//    cerr << "Examining atom " << k << ", counted = " << counted[k] << endl;
      if (counted[k] || _counted[k])     // the Atom_Scoring_Object knows about this atom already
        continue;

      rc++;
      _score += m->atomic_number (k);

      _atom_list[_ncounted] = k;
      _ncounted++;
      _counted[k] = 1;
    }
  }

  _start_of_last_generation = istop;

  if (0 == rc)
    _stopped = 1;
  else
    _bonds_counted++;

  return rc;
}

/*
  We update the global array of counted atoms. Note that our array is a list of
  atom numbers, whereas the global array is a boolean.
*/

int
Path_Scoring_Object::update_global_counted_array (int matoms, int * counted)
{
  assert (ok ());
  assert (NULL != counted);

  return 1;
  for (int i = 0; i < matoms; i++)
    counted[i] = _counted[i];

  return 1;
}

/*
  An atom keeps track of all Path_Scoring_Object's emanating from its start atom.
*/

template class resizable_array_p<Path_Scoring_Object>;
template class resizable_array_base<Path_Scoring_Object *>;

class Atom_Scoring_Object
{
  friend
    operator << (ostream &, const Atom_Scoring_Object &);

  private:
    atom_number_t _root_atom;
    int * _counted;
    resizable_array_p<Path_Scoring_Object> _paths;
    int _score;
    int _bonds_counted;

//  private functions
  
    void _default_values (Molecule *, atom_number_t);
    int  _update_global_counted_array (int);

  public:
    Atom_Scoring_Object (Molecule *, int);
    ~Atom_Scoring_Object ();

    int ok () const;
    int debug_print (ostream & = cerr) const;

    int bonds_counted () const { return _bonds_counted;}
    int score () const { return _score;}

    int advance_to_next_bonds (Molecule *);
    int can_advance () const;
    int all_scores_different () const;
};

void
Atom_Scoring_Object::_default_values (Molecule * m, atom_number_t a)
{
  assert (OK_ATOM_NUMBER (m, a));

  _root_atom = a;

  int matoms = m->natoms ();

  _score = 0;

  _counted = new_int (matoms);
  _counted[a] = 1;
  _score = m->atomic_number (a);

  int acon = m->ncon (a);
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = m->other (a, i);
    _counted[j] = 1;
    Path_Scoring_Object * p = new Path_Scoring_Object (m, a, j);
    _paths.add (p);
    _score += p->score ();
  }

  _bonds_counted = 1;

  _update_global_counted_array (matoms);

  return;
}

Atom_Scoring_Object::Atom_Scoring_Object (Molecule * m, atom_number_t a)
{
  _default_values (m, a);

  return;
}

Atom_Scoring_Object::~Atom_Scoring_Object ()
{
  assert (ok ());

  if (-654 == _bonds_counted)
    cerr << "Atom_Scoring_Object:: deleting an already deleted object\n";

  delete _counted;
  _counted = NULL;

  _score = -76;

  _root_atom = INVALID_ATOM_NUMBER;

  _bonds_counted = -654;

  return;
}

int
Atom_Scoring_Object::ok () const
{
  if (NULL == _counted)
    return 0;

  if (_score < 0)
    return 0;

  if (_bonds_counted < 0)
    return 0;

  if (_root_atom < 0)
    return 0;

  return 1;
}

int
Atom_Scoring_Object::debug_print (ostream & os) const
{
  assert (os.good ());

  if (! ok ())
    os << "Atom_Scoring_Object:: warning, ok fails\n";

  os << "Atom_Scoring_Object:: for atom " << _root_atom << ", " << _bonds_counted <<
        " bonds counted, score = " << _score << endl;

  int np = _paths.number_elements ();
  os << np << " bonds to start atom\n";
  for (int i = 0; i < np; i++)
    _paths[i]->debug_print (os);

  return 1;
}

int
Atom_Scoring_Object::_update_global_counted_array (int matoms)
{
  int np = _paths.number_elements ();
  for (int i = 0; i < np; i++)
  {
    Path_Scoring_Object * p = _paths[i];
    p->update_global_counted_array (matoms, _counted);
  }

  return 1;
}

int
Atom_Scoring_Object::advance_to_next_bonds (Molecule * m)
{
  assert (ok ());
  assert (OK_MOLECULE (m));

  int np = _paths.number_elements ();
  _score = 0;
//cerr << "Atom scoring object advancing, np = " << np << endl;
  int rc = 0;
  for (int j = 0; j < np; j++)
  {
    Path_Scoring_Object * p = _paths[j];
    int tmp = p->advance_to_next_generation (m, _counted);
    if (tmp)
    {
      rc++;
      _score += p->score ();
    }
  }
  _bonds_counted++;
  _update_global_counted_array (m->natoms ());

  return rc;
}

/*int Atom_Scoring_Object::score () const
{
  assert (ok ());
  assert (bonds_counted <= _bonds_counted);

  int rc = 0;
  int np = _paths.number_elements ();
  for (int i = 0; i < np; i++)
  {
    const Path_Scoring_Object * p = _paths[i];
    rc += p->score ();
  }

  return rc;
}*/

int
Atom_Scoring_Object::all_scores_different () const
{
  int np = _paths.number_elements ();

  for (int i = 0; i < np; i++)
  {
    const Path_Scoring_Object * pi = _paths[i];
    for (int j = i + 1; j < np; j++)
    {
      const Path_Scoring_Object * pj = _paths[j];
      if (pi->score () == pj->score ())
        return 0;
    }
  }

  return 1;
}

/*
  We need to know whether or not scoring can usefully progress to the next bond.

  If there are stopped scores which are identical, then no use continuing.
  If all paths are stopped, no use continuing.
*/

int
Atom_Scoring_Object::can_advance () const
{
  assert (ok ());

  int np = _paths.number_elements ();

  int stopped_paths = 0;
  for (int i = 0; i < np; i++)
  {
    const Path_Scoring_Object * pi = _paths[i];
    if (pi->stopped ())
      stopped_paths++;

    for (int j = i + 1; j < np; j++)
    {
      const Path_Scoring_Object * pj = _paths[j];
      if (pi->score () == pj->score () && pi->stopped () && pj->stopped ())
        return 0;
    }
  }

  if (np == stopped_paths)
    return 0;

  return 1;
}

int
Molecule::is_possible_chiral_centre (atom_number_t a)
{
  assert (OK_ATOM_NUMBER (this, a));

  int ih = implicit_hydrogens (a);
  if (ih > 1)
    return 0;

  int acon = ncon (a);
  if (acon + ih < 4)
    return 0;

  if (acon + ih > 4)                // don't know how to handle these yet - hopefully never needed.
    return 0; 

// If there are two hydrogens present, or an implicit hydrogen and an explicit one, then no
// chiral centre.

  int eh = explicit_hydrogens (a);
  if (eh + ih > 1)
    return 0;

// Now, we have a centre with more than 4 connections, and one or fewer implicit hydrogens.
// We must score it.

  Atom_Scoring_Object score (this, a);
//cerr << "After creation ";
//score.debug_print (cerr);

  while (1)
  {
    if (! score.advance_to_next_bonds (this))
      return 0;

//  cerr << "After advance to " << score.bonds_counted () << endl;
//  score.debug_print (cerr);

    if (score.all_scores_different ())
      return 1;

    if (! score.can_advance ())
      return 0;
  }
}

/*
  This relies on the equivalence of chiral_type_t and int types
*/

int
Molecule::_allocate_chirality_storage ()
{
  assert (NULL == _chirality);

  _chirality = new_int (_number_elements, CHIRALITY_NOT_DETERMINED);

  return 1;
}

chiral_type_t
Molecule::chirality (atom_number_t a)
{
  assert (OK_ATOM_NUMBER (this, a));

  if (NULL != _chirality && CHIRALITY_NOT_DETERMINED != _chirality[a])
    return _chirality[a];

  if (NULL == _chirality)
    _allocate_chirality_storage ();

  if (is_possible_chiral_centre (a))
    return _chirality[a] = CHIRAL_BUT_UNSPECIFIED;
  else
    return _chirality[a] = NON_CHIRAL;
}

int
Molecule::perceive_chirality ()
{
  assert (ok ());


  int rc = 0;    // we return the number of chiral centres found

  int matoms = natoms ();

  if (NULL == _chirality)
    _allocate_chirality_storage ();

  for (int i = 0; i < matoms; i++)
  {
    if (CHIRALITY_NOT_DETERMINED == _chirality[i])
    {
      if (is_possible_chiral_centre (i))
      {
        _chirality[i] = CHIRAL_BUT_UNSPECIFIED;
        rc++;
      }
      else
        _chirality[i] = NON_CHIRAL;
    }
    else if (NON_CHIRAL != _chirality[i])
      rc++;
  }

  return rc;
}

int
Molecule::set_chirality (atom_number_t a, chiral_type_t c)
{
  assert (OK_ATOM_NUMBER (this, a));

  if (NULL == _chirality)
    _allocate_chirality_storage ();

  invalidate_smiles ();    // perhaps we should do _set_modified??

  _chirality[a] = c;

  return 1;
}

int
Molecule::unspecified_chiral_centres ()
{
  assert (ok ());

  perceive_chirality ();

  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
    if (CHIRAL_BUT_UNSPECIFIED == _chirality[i])
      rc++;

  return rc;
}
