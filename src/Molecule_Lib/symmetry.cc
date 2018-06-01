#include <iostream>
#include <memory>
using std::cerr;
using std::endl;


//#define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#include "misc.h"

#include "molecule.h"
#include "path.h"
#include "misc2.h"

#include "symmetry.h"

void
Symmetry_Info::_default_values ()
{
  _matoms = 0;
  _symmetry_class = NULL;
  _degree = NULL;

  return;
}

Symmetry_Info::Symmetry_Info ()
{
  _default_values();
}

Symmetry_Info::~Symmetry_Info ()
{
  assert (ok());

  DELETE_IF_NOT_NULL(_symmetry_class);

  DELETE_IF_NOT_NULL(_degree);

  return;
}

int
Symmetry_Info::ok () const
{
  if (_number_elements && NULL == _symmetry_class)
    return 0;

  if (_symmetry_class && 0 == _number_elements)
    return 0;

  if (NULL == _symmetry_class && 0 == _number_elements && 0 == _matoms)
    return 1;

  if (_number_elements > _matoms)
    return 0;

  return 1;
}

int
Symmetry_Info::debug_print (std::ostream & os) const
{
  os << "Info on symmetry info for molecule with " << _matoms << " atoms and " <<
        _number_elements << " symmetry classes\n";

  for (int i = 0; i < _number_elements; i++)
  {
    const Set_of_Atoms * a = _things[i];

    os << "Symmetry class " << i << " atoms " << (*a) << endl;
  }

  int ng = _symmetry_groupings.number_elements();
  if (ng)
  {
    os << ng << " symmetric groupings\n";
    for (int i = 0; i < ng; i++)
    {
      const Symmetric_Atoms * s = _symmetry_groupings[i];
      os << " grouping " << i << " contains " << s->number_elements();
      if (s->number_elements() == _matoms)
        os << " (all)";
      os << " atoms, degree " << s->degree() << endl;

      os << " Atoms";
      for (int j = 0; j < s->number_elements(); j++)
        os << ' ' << s->item(j);

      os << endl;
    }
  }

  return os.good();
}

/*
  Atom J was initially identified as non symmetric. Is it bonded to
  symmetric atoms?

  We return the lowest number of connections with the same symmetry number

  For example, imagine we have something 5 connected. Three atoms
  are in the same symmetry class, and the other two in a different
  class. We need to return 2 in that case, as any symmetry grouping
  containing this atom will have, at most, twofold symmetry.
*/

/*int
Symmetry_Info::_determine_degree (Molecule & m, atom_number_t zatom)
{
  const Atom * a = m.atomi (zatom);
  int acon = a->ncon ();
  if (1 == acon)     // cannot be a centre
    return 0;

  resizable_array<int> symmetry_classes_found;
  symmetry_classes_found.resize (acon);

// For efficiency, we check the case of all values the same as we go along

  int all_the_same = 1;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other (zatom, i);
    if (_symmetry_class[j])
    {
      symmetry_classes_found.add (_symmetry_class[j]);
      if (all_the_same)
      {
        int nf = symmetry_classes_found.number_elements ();
        if (nf > 1 && symmetry_classes_found[nf - 2] != symmetry_classes_found[nf - 1])
          all_the_same = 0;
      }
    }
  }

  int nscf = symmetry_classes_found.number_elements ();

// If no symmetric atoms are bonded, this is not symmetric

  if (0 == nscf)
    return 0;

// It is impossible to have just one symmetric atom bonded to a so-far non-symmetric one

  assert (nscf > 1);

// The case of two symmetric attachments will be common

  if (2 == nscf)
    return 2;

// If all neighbours have the same symmetry class, that's easy...

  if (all_the_same)
    return nscf;

// Now things get ugly. We have 3 or more symmetric connections, not all the same.
// consider nscf == 3. We could have three connections the same, or 2 + 1
// As 3 is a very common case, do it without a sort

  if (3 == nscf)
  {
    if (symmetry_classes_found[0] != symmetry_classes_found[1])
      return 2;
    if (symmetry_classes_found[0] != symmetry_classes_found[2])
      return 2;

    return 3;     // they are all the same
  }

// Work out the smallest number of symmetry classes the same. Note that
// we skip over case of just one member of that class

  symmetry_classes_found.sort (int_comparitor_larger);

  int prev = -1;
  int min_number_found = nscf + 1;
  int number_found = 0;
  for (int i = 0; i < nscf; i++)
  {
    int isym = symmetry_classes_found[i];
    if (isym != prev)
    {
      if (number_found > 1 && number_found < min_number_found)
        min_number_found = number_found;

      prev = isym;
      number_found = 1;
    }
    else
      number_found++;
  }

  if (number_found > 1 && number_found < min_number_found)
    min_number_found = number_found;

// We must have found at least two the same

  assert (min_number_found < nscf);

  return min_number_found;
}*/

//#define DEBUG_SYMINFO_INITIALISE

int
Symmetry_Info::initialise (Molecule & m)
{
  _matoms = m.natoms ();

  _symmetry_class = new int[_matoms];

   copy_vector (_symmetry_class, m.symmetry_classes (), _matoms);

  _degree = new_int (_matoms);

#ifdef DEBUG_SYMINFO_INITIALISE
  cerr << "Symmetry classes\n";
  for (int i = 0; i < _matoms; i++)
  {
    cerr << "Atom " << i << " in class " << _symmetry_class[i] << endl;
  }
#endif

  int atoms_assigned = 0;

  for (int i = 1; i <= _matoms; i++)
  {
    int jstart = locate_item_in_array(i, _matoms, _symmetry_class);
    if (jstart < 0)
      continue;

    Set_of_Atoms * a = new Set_of_Atoms;
    a->resize(3);       // just a guess, yes, this is wasteful...
    for (int j = jstart; j < _matoms; j++)
    {
      if (_symmetry_class[j] == i)
        a->add(j);
    }

//  and update the _degree array while we are at it.

    int na = a->number_elements();

    atoms_assigned += na;

    if (na > 1)
      a->set_vector(_degree, na);

    assert(a->number_elements());

    add(a);

//  cerr << "Found group " << (*a) << endl;

    if (atoms_assigned == _matoms)
      break;
  }

  assert (atoms_assigned <= _matoms);

  if (! _compute_symmetry_groupings(m))
    return 0;

  assert (ok());

  return _number_elements;
}

int
Symmetry_Info::_compute_symmetry_groupings (Molecule & m,
                               const int * symmetric_neighbours,
                               int * already_done)
{
  for (int i = 0; i < _number_elements; i++)
  {
    const Set_of_Atoms * s = _things[i];
    int ns = s->number_elements();
    if (1 == ns)     // a possibly asymmetric atom
      continue;

    if (already_done[s->item(0)])
      continue;

    Symmetric_Atoms * sa = new Symmetric_Atoms;

    sa->build(m, s->item(0), ns, symmetric_neighbours, _degree, already_done);

    _symmetry_groupings.add(sa);
  }

  return _symmetry_groupings.number_elements();
}

/*
  Private function called only from initialise
*/

int
Symmetry_Info::_compute_symmetry_groupings (Molecule & m)
{
  if (0 == _number_elements)    // if there is no symmetry, we are done
    return 0;

  int * symmetric_neighbours = new_int(_matoms); std::unique_ptr<int[]> free_symmetric_neighbours(symmetric_neighbours);

#ifdef SLOWVERSIONWITHATOMS
  for (int i = 0; i < _matoms; i++)
  {
    const Atom * a = m.atomi(i);
    int acon = a->ncon();
    int sn = 0;
    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

      if (_degree[k])
        sn++;
    }

    symmetric_neighbours[i] = sn;
  }
#endif

  for (int i = m.nedges() - 1; i >= 0; i--)
  {
    const Bond * b = m.bondi(i);

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (_degree[a1])
      symmetric_neighbours[a2]++;
    if (_degree[a2])
      symmetric_neighbours[a1]++;
  }

  int * already_done = new_int(_matoms); std::unique_ptr<int[]> free_already_done(already_done);

  return _compute_symmetry_groupings(m, symmetric_neighbours, already_done);
}

/*
  We usually aren't interested in CF3 and T-Butyl, and other
  symmetries where 3 atoms are all joined to a common base.
*/

int
Symmetry_Info::remove_trivial_cf3 (const Molecule & m)
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    const Set_of_Atoms * s = _things[i];

    int ns = s->number_elements();
    if (3 != ns)
      continue;

    const Atom * a = m.atomi(s->item(0));
    if (1 != a->ncon())
      continue;

    atom_number_t anchor = a->other(s->item(0), 0);

    int remove_it = 1;
    for (int j = 1; j < ns; j++)
    {
      atom_number_t k = s->item(j);

      if (anchor != m.other(k, 0))
      {
        remove_it = 0;
        break;
      }
    }

    if (remove_it)
    {
      remove_item(i);
      _remove_symmetry_grouping(*s);
      i--;
      rc++;
    }
  }

  return rc;
}

static int
is_part_of_benzene_ring (const Set_of_Atoms & s,
                         const resizable_array<const Ring *> & benzene_rings)
{
  for (int i = 0; i < benzene_rings.number_elements(); i++)
  {
    const Ring * ri = benzene_rings[i];
    if (ri->contains(s[0]) && ri->contains(s[1]))    // only really need to check the first one
      return 1;
  }

  return 0;
}

/*
  We are looking for a benzene ring with just 1 connection outside the ring
*/

static int
is_benzene_ring (const Molecule & m,
                 const Ring & ri)
{
  assert (6 == ri.number_elements());

  int three_connected_ring_atoms = 0;

  for (int j = 0; j < 6; j++)
  {
    atom_number_t k = ri[j];
    if (6 != m.atomic_number(k))
      return 0;

    if (3 == m.ncon(k))
    {
      three_connected_ring_atoms++;
      if (three_connected_ring_atoms > 1)
        return 0;
    }
  }

  return (1 == three_connected_ring_atoms);
}

int
Symmetry_Info::remove_benzene (Molecule & m)
{
  if (m.natoms() <= 6)
    return 0;

  int nr = m.nrings();
  if (0 == nr)
    return 0;

  m.compute_aromaticity_if_needed();

// Look for benzene rings

  resizable_array<const Ring *> benzene_rings;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (6 != ri->number_elements())
      continue;

    if (ri->is_fused())
      continue;

    if (! ri->is_aromatic())
      continue;

    if (is_benzene_ring(m, *ri))
      benzene_rings.add(ri);
  }

  if (0 == benzene_rings.number_elements())
    return 0;

#ifdef DEBUG_REMOVE_BENZENE
  cerr << "Found " << benzene_rings.number_elements() << " benzene rings\n";
#endif

// Now remove all the 2 member symmetry groupings involving these rings

  int rc = 0;
  for (int i = _number_elements - 1; i >= 0; i--)
  {
    const Set_of_Atoms * s = _things[i];

#ifdef DEBUG_REMOVE_BENZENE
//  cerr << "Examining " << (*s) << endl;
#endif
    
    int ns = s->number_elements();

    if (2 != ns)
      continue;

    if (! m.is_aromatic(s->item(0)))
      continue;

    if (! is_part_of_benzene_ring(*s, benzene_rings))
      continue;

#ifdef DEBUG_REMOVE_BENZENE
    cerr << "Yes, part of benzene ring\n";
#endif

    remove_item(i);
    rc++;
  }

  if (rc)
    _remove_symmetry_groupings_that_are_just_benzene(benzene_rings, rc);

  return rc;
}

Symmetric_Atoms::Symmetric_Atoms ()
{
  _degree = 0;
}

//#define DEBUG_BUILD_SYMMETRIC_ATOMS

/*
  Build large scale collections of symmetry related atoms

  SYMMETRY_CLASS is what comes back from the molecule object.

  Particularly interesting are atoms which are not symmetry
  related to any other atom (0 == symmetry_class[i]), but
  which are connected to symmetric atoms. On entry SYMMETRY_COUNT
  must contain the number of symmetric atoms attached to each
  atom in the molecule

  This is buggy. The problem is that the real symmetry of the molecule
  isn't really known until this procedure completes, but this is the
  procedure that builds the real symmetry. 

  This shows the problem

  N(C1=CC=C(F)C=C1)(S(=O)(=O)C1=CC(Cl)=CC(=C1)Cl)CC(OCC)OCC PBCHM2800580

  Ideally we would include the central N atom with the symmetry grouping
  that includes the Sulphone and the benzene ring, but we haven't perceived
  that yet.
*/

int
Symmetric_Atoms::_build (Molecule & m,
                         atom_number_t zatom,
                         const int * symmetric_neighbours,
                         const int * degree,
                         int myflag,
                         int * already_done)
{
  assert (0 == already_done[zatom]);

  already_done[zatom] = myflag;
  add(zatom);

  if (0 == degree[zatom])   // must be two-fold symmetric (not symmetric by itself)
    ;
  else if (degree[zatom] < _degree)
    _degree = degree[zatom];    // the lowest degree of any atom component

  const Atom * a = m.atomi(zatom);
  int acon = a->ncon();

// First compute the number of unprocessed neighbours with no symmetry

  int non_symmetric_neighbours = 0;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (myflag == already_done[j])
      continue;

    if (0 == degree[j])
      non_symmetric_neighbours++;
  }

#ifdef DEBUG_BUILD_SYMMETRIC_ATOMS
  cerr << "Symmetric_Atoms::_build: continuing with atom " << zatom << " degree " << degree[zatom] << ", " << non_symmetric_neighbours << " non symmetric neighbours\n";
#endif

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);

    if (myflag == already_done[j])
      continue;

    if (degree[j] && 0 == already_done[j])    // all symmetric attachments to a symmetric atom are included
    {
#ifdef DEBUG_BUILD_SYMMETRIC_ATOMS
      cerr << "Continuing with neighbour " << j << " degree " << degree[j] << endl;
#endif

      _build(m, j, symmetric_neighbours, degree, myflag, already_done);
      continue;
    }

    int jcon = m.ncon(j);

//  We have a singly connected atom hanging off

    if (1 == jcon && 1 == non_symmetric_neighbours)
    {
#ifdef DEBUG_BUILD_SYMMETRIC_ATOMS
    if (1 == jcon)
      cerr << "neighbour " << j << " is singly connected\n";
#endif

      add(j);
      continue;
    }

// Now the more difficult case where atom J is not symmetric by itself
// - nor is atom I. Two connected neighbours are always added

//  Continue down a long chain

    if (1 == non_symmetric_neighbours && 2 == jcon && 0 == b->nrings() && 0 == already_done[j])
    {
#ifdef DEBUG_BUILD_SYMMETRIC_ATOMS
      cerr << " neighbour " << j << " is doubly connected\n";
#endif

      _degree = 2;    // down to two-fold symmetry
      _build(m, j, symmetric_neighbours, degree, myflag, already_done);
      continue;
    }

//  If atom J is bonded to some (other) symmetric atoms, we continue with it

    if (symmetric_neighbours[j] && 0 == already_done[j])
    {
#ifdef DEBUG_BUILD_SYMMETRIC_ATOMS
    if (symmetric_neighbours[j])
      cerr << "Neighbour " << j << " has " << symmetric_neighbours[j] << " symmetric neighbours\n";
#endif

      _build(m, j, symmetric_neighbours, degree, myflag, already_done);
      continue;
    }

//  If we have come to a ring, that isn't symmetric, just add the atom

    if (0 == b->nrings() && m.nrings(j))
    {
#ifdef DEBUG_BUILD_SYMMETRIC_ATOMS
      cerr << "neighbour " << j << " is at a ring\n";
#endif

      add(j);
      continue;
    }

//  We may be at a stopping point

    if (1 == non_symmetric_neighbours && 0 == b->nrings())
    {
      _degree = 2;
      add(j);
      continue;
    }
  }

  return _number_elements;
}

int
Symmetric_Atoms::build (Molecule & m,
                        atom_number_t zatom,
                        int initial_degree,
                        const int * symmetric_neighbours,
                        const int * degree,
                        int * already_done)
{
  int matoms = m.natoms();

  resize(matoms);

  _degree = initial_degree;

  return _build(m, zatom, symmetric_neighbours, degree, zatom + 1, already_done);
}

/*
  We are trying to figure out if Symmetric_Atoms grouping R2 is a benzene
  ring that covers R1
*/

static int
no_more_than_N_atoms_mismatched (const Set_of_Atoms & r1,
                                 const Set_of_Atoms & r2,
                                 int max_mismatches)
{
  int nr1 = r1.number_elements();

  assert (nr1 < r2.number_elements());

  int mismatches = 0;

  for (int i = 0; i < nr1; i++)
  {
    if (r2.contains(r1[i]))
      continue;

    mismatches++;
    if (mismatches > max_mismatches)
      return 0;
  }

  return 1;
}


static int
is_benzene_ring (const resizable_array<const Ring *> & benzene_rings,
                 const Set_of_Atoms & r)
{
  for (int i = 0; i < benzene_rings.number_elements(); i++)
  {
    const Ring * bri = benzene_rings[i];

    if (no_more_than_N_atoms_mismatched(*bri, r, 1))
      return 1;
  }

  return 0;
}

int
Symmetry_Info::_remove_symmetry_groupings_that_are_just_benzene (const resizable_array<const Ring *> & benzene_rings,
                                       int nb)
{
  int rc = 0;

  for (int i = _symmetry_groupings.number_elements() - 1; i >= 0; i--)
  {
    const Symmetric_Atoms * s = _symmetry_groupings[i];

    if (7 != s->number_elements())
      continue;

    if (is_benzene_ring(benzene_rings, *s))
    {
      _symmetry_groupings.remove_item(i);
      rc++;
      if (rc == nb)
        break;
    }
  }

  return rc;
}

/*
  We are discarding symmetry info about something like a CF3 or T butyl.
  Find the corresponding Symmetric_Atoms and remove it
*/

int
Symmetry_Info::_remove_symmetry_grouping (const Set_of_Atoms & s)
{
  assert (s.number_elements());

  for (int i = 0; i < _symmetry_groupings.number_elements(); i++)
  {
    const Symmetric_Atoms * smgpi = _symmetry_groupings[i];

    if (no_more_than_N_atoms_mismatched(s, *smgpi, 1))
    {
      _symmetry_groupings.remove_item(i);
      return 1;
    }
  }

  cerr << "Symmetry_Info::_remove_symmetry_grouping: could not remove group for " << s << endl;
  cerr << s.number_elements() << " atoms\n";

  return 0;
}

static int
all_atoms_separated_enough (Molecule & m,
                            const Set_of_Atoms & s,
                            int nb)
{
  int n = s.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = s[i];

    for (int k = i + 1; k < n; k++)
    {
      if (m.bonds_between(j, s[k]) < nb)
        return 0;
    }
  }

  return 1;   // all separations OK
}

/*int
Symmetry_Info::remove_symmetry_closer_than (Molecule & m,
                                            int nb)
{
  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
  {
    const Set_of_Atoms * si = _things[i];

    int ns = s->number_elements ();

    if (all_atoms_separated_enough (m, *s, nb))
      continue;
  }

}*/
