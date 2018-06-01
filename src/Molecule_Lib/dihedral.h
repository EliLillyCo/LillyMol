#ifndef DIHEDRAL_H
#define DIHEDRAL_H

/*
  We often need a specifier for dihedral atoms, and the atom
  number at which the moving atoms defined by the dihedral angle
  begin within a molecule.
*/

#include <iostream>

#include "iwmtypes.h"

#define DIHEDRAL_MAGIC -973151

class Dihedral_Atoms
{
  protected:
    atom_number_t  _a1, _a2, _a3, _a4;
    magic_number_t _magic;

  public:
    Dihedral_Atoms ();
    Dihedral_Atoms (atom_number_t, atom_number_t, atom_number_t, atom_number_t);
    ~Dihedral_Atoms ();

    void invalidate ();
    int ok () const;
    int debug_print (ostream &) const;

    void set (atom_number_t, atom_number_t, atom_number_t, atom_number_t);

    atom_number_t a1 () const;
    atom_number_t a2 () const;
    atom_number_t a3 () const;
    atom_number_t a4 () const;
};

#define OK_MOL_DIHEDRAL(m, d) \
          ( (d)->ok () && \
	    OK_4_ATOMS ((m), (d)->a1 (), (d)->a2 (), \
	                     (d)->a3 (), (d)->a4 ()) )

/*
  We often need not only the specifier or a dihedral angle, but also
  info on which is the first atom in the molecule to move.
*/

#define TWIST_MAGIC -117722

#define OK_TWIST(t) (NULL != (t) && (t)->ok ())

#define OK_MOL_TWIST(m, t) \
          ( OK_MOL_DIHEDRAL((m), (t)) && \
	    (t).ok ());
	    

class Twist_Atoms : public Dihedral_Atoms
{
  private:
    magic_number_t _magic;
    atom_number_t  _frag_start;

  public:
    Twist_Atoms (atom_number_t, atom_number_t, atom_number_t, atom_number_t,
                 atom_number_t);
    ~Twist_Atoms ();

    boolean ok () const;
    int     debug_print (ostream &) const;

    atom_number_t frag_start () const;
};

extern ostream &
operator << (ostream &, const Twist_Atoms &);

#endif
