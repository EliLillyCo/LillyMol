#ifndef IW_SYMMETRY_H
#define IW_SYMMETRY_H

/*
  When perceiving whole molecule symmetry, and large scale symmetry
  within the molecule, we need a structure to hold that.
*/

class Symmetric_Atoms : public Set_of_Atoms
{
  private:
    int _degree;

// private functions

    int _build (Molecule &, atom_number_t, const int *, const int *, int, int *);

  public:
    Symmetric_Atoms ();

    int degree () const { return _degree;}

    int build (Molecule &, atom_number_t, int, const int *, const int *, int *);
};

/*
  When studying symmetry classes in a molecule, we need a helper object
  to hold more info about the molecule

  Each symmetry class gets its own set of atoms
*/

class Symmetry_Info : public resizable_array_p<Set_of_Atoms>
{
  private:
    int _matoms;

    int * _symmetry_class;    // from the molecule

//  Degree will be the number of atoms in the same symmetry class (3 for each
//  of the F atoms in CF3, or 6 for benzene carbon's)

    int * _degree;

    resizable_array_p<Symmetric_Atoms> _symmetry_groupings;

//  private functions

    void _default_values ();
    int  _determine_degree (Molecule & m, atom_number_t zatom);

    int _compute_symmetry_groupings (Molecule & m, const int * symmetric_neighbours, int * already_done);
    int _compute_symmetry_groupings (Molecule & m);

    int _remove_symmetry_grouping (const Set_of_Atoms &);
    int _remove_symmetry_groupings_that_are_just_benzene (const resizable_array<const Ring *> & benzene_rings, int nb);

  public:
    Symmetry_Info ();
    ~Symmetry_Info ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int initialise (Molecule &);

    int remove_trivial_cf3 (const Molecule &);
    int remove_benzene (Molecule &);

    int symmetry_groupings () const { return _symmetry_groupings.number_elements ();}
    const Symmetric_Atoms * symmetry_grouping (int i) const { return _symmetry_groupings[i];}

    void remove_symmetry_grouping (int i) { _symmetry_groupings.remove_item (i);}

    const int * symmetry_class () const { return _symmetry_class;}

    const int * degree () const { return _degree;}
};

#endif
