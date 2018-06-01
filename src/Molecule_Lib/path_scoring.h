#ifndef PATH_SCORING_H
#define PATH_SCORING_H

class Atom;
#include "iwmtypes.h"
#include "molecule.h"

/*
  For R and S determinations and E/Z determinations, we need a means
  of scanning through the molecule.

  One of the underlying objects is a list of atomic numbers encountered.

  Sept 99. This is incorrect. Ran into a case (LLY 484066) where just
  keeping track of atomic numbers is not good enough, as that doesn't
  distinguish between a singly bonded and multiply bonded atom.

  For this reason, we reserve 4 slots for each atomic number.
    4 * z = single bond
    4 * z + 1 = aromatic bond
    4 * z + 2 = double bond
    4 * z + 3 = triple bond

  Also, I find that we must implement a set of Atomic_Numbers_Encounterd
  for each "layer" in the R/S determination. Otherwise we cannot distinguish
  between C-N and N-C

  Mar 2000. Still not right.
  We must account for the connectivity of the atoms encountered. John Tomer
  sent me examples of asymmetric norboranes that failed.

  For example

     C1CC2C=C(C1C2)C(O)=O PBCHM281894

  the path scoring runs around the ring and decided it is done. But, it
  hasn't realised that one of the bonds at the end of the double bond is
  3 connected, and the other is 2 connected.

  Therefore we need to take the connectivity into account. We allow
  up to 4 connections. Therefore, each atom needs 11 slots

    11 * z      = single_bond, 1 connection
    11 * z + 1  = single_bond, 2 connections
    11 * z + 2  = single_bond, 3 connections
    11 * z + 3  = single_bond, 4 connections
    11 * z + 4  = aromatic_bond, 2 connections
    11 * z + 5  = aromatic_bond, 3 connections
    11 * z + 6  = double_bond, 1 connections
    11 * z + 7  = double_bond, 2 connections
    11 * z + 8  = double_bond, 3 connections
    11 * z + 9  = triple_bond, 1 connections
    11 * z + 10 = triple_bond, 2 connections

    Aug 02. What if we try to do a chirality determination on a molecule
    that contains non-periodic table elements? Since these can be encountered
    by any kind of bond, we use a pseudo atomic number that is computed
    from the atomic symbol.

    Allow one number for each single letter element, one for each two letter element,
    and one for everything else
*/

#define NUMBER_OF_POSSIBLE_ELEMENTS (26 + 26*26 + 1 + 1)   // needs to be 1 larger than the highest possible value

#define SIZE_OF_ATOMIC_NUMBERS_ENCOUNTERED_ARRAY (11 * NUMBER_OF_POSSIBLE_ELEMENTS)

class Atomic_Numbers_Encounterd
{
  friend
    std::ostream & operator << (std::ostream &, const Atomic_Numbers_Encounterd &);

  private:
    atomic_number_t _found[SIZE_OF_ATOMIC_NUMBERS_ENCOUNTERED_ARRAY];

  public:
    Atomic_Numbers_Encounterd ();

    int  extra (const Element *);
    int  extra (const Bond * b, atomic_number_t z);
    int  extra (const Bond * b, const Element * e, int ncon);

    void initialise ();

    int compare (const Atomic_Numbers_Encounterd &) const;

    int  operator <  (const Atomic_Numbers_Encounterd &) const;
    int  operator >  (const Atomic_Numbers_Encounterd &) const;
    int  operator != (const Atomic_Numbers_Encounterd &) const;
    int  operator == (const Atomic_Numbers_Encounterd &) const;
};

/*
  The molecule will allocate one of these objects for each atom to be
  differentiated.
  First call initialise - to give it an atom with which to start

  while (unresolved)
  {
    call advance.
    test for differentiated
    call update_claimed
  }

  We keep a Atomic_Numbers_Encounterd object for each step we take.

  During a scan through a molecule, different paths may become
  blocked at different times. We need to record the number of
  steps each path scoring object has taken

  Feb 2015. In some cases, the first bond matters.
*/

class Path_Scoring: public resizable_array_p<Atomic_Numbers_Encounterd>
{
  friend
    std::ostream & operator << (std::ostream &, const Path_Scoring &);

  private:
    atom_number_t _start_atom;

    int _first_bond;

    int _assigned_rank;

//  As we advance through the molecule, we need to keep track of the atoms
//  which are at the edge of our advance.

    Set_of_Atoms _edge_atom;

    int _nsteps;

  public:
    Path_Scoring ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int operator == (const Path_Scoring &) const;
    int operator != (const Path_Scoring &) const;

    int compare (const Path_Scoring &) const;

    int active () const { return _edge_atom.number_elements ();}

    int nsteps () const { return _nsteps;}

    int initialise (atom_number_t, const Atom *);
    void set_first_bond (int s) { _first_bond = s;}

    atom_number_t start_atom () const { return _start_atom;}

    int assigned_rank () const { return _assigned_rank;}
    void set_assigned_rank (int r) { _assigned_rank = r;}

    int advance (Atom * const *, const int *);

//  After all Path_Scoring objects have advanced, they need to update the
//  global array of atoms claimed

    int update_claimed (int *) const;
};

/*
  Resolving a set of these involves sorting the array and then looking
  through the array to see if all neighbours are equal
*/

extern int resolved (resizable_array_p<Path_Scoring> &, int &);

extern int assign_ranks (resizable_array_p<Path_Scoring> & ps);

#endif
