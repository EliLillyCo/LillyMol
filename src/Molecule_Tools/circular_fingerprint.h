#ifndef CIRCULAR_FINGERPRINT_GENERATOR_H
#define CIRCULAR_FINGERPRINT_GENERATOR_H

#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Molecule_Lib/molecule.h"

/*
  Class for generating circular fingerprints for a molecule.
  Note that much of the storage is specific to the molecule
  being processed, so no thread safety here
*/

class Circular_Fingerprint_Generator 
{
  private:
    int _additive;
    int _matoms;
    int _allocated;
    int _min_radius;
    int _max_radius;
    const Atom ** _atom;
    int * _atype;
    int * _processing_status;

//  private functions

    void _free_dynamic_arrays ();
    int _initialise_molecule (const Molecule & m, const int * atype, const int * include_these_atoms);
    int _generate_shells (int radius, unsigned int sum_so_far, Sparse_Fingerprint_Creator & sfc);
    void _increment (unsigned int & sum_so_far, int bc, int atom_constant) const;

  public:
    Circular_Fingerprint_Generator();
    ~Circular_Fingerprint_Generator();

    void set_min_radius(int s) { _min_radius = s;}
    void set_max_radius(int s) { _max_radius = s;}

    void set_additive (int s) { _additive = s;}

    int generate_fingerprint (const Molecule & m, const int * atype, const int * include_these_atoms, Sparse_Fingerprint_Creator &);
};

void set_default_circular_fingerprint_additive (int);

#endif
