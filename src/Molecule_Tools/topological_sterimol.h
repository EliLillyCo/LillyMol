#ifndef TOPOSTERIMOL_H
#define TOPOSTERIMOL_H

#include <iostream>

#include "Molecule_Lib/molecule.h"

/*
  We use the Sterimol object to hold the results of the 3D computations.
  We need another object to hold the topological results
*/

#define TOPOLOGICAL_STERIMOL_NATOMS 0
#define TOPOLOGICAL_STERIMOL_HETERO 1
#define TOPOLOGICAL_STERIMOL_DONOR 2
#define TOPOLOGICAL_STERIMOL_ACCEPT 3
#define TOPOLOGICAL_STERIMOL_NRINGS 4
#define TOPOLOGICAL_STERIMOL_AROM 5
#define TOPOLOGICAL_STERIMOL_MXDIST 6
#define TOPOLOGICAL_STERIMOL_MXACCDIST 7
#define TOPOLOGICAL_STERIMOL_MXDONDIST 8
#define TOPOLOGICAL_STERIMOL_POS 9
#define TOPOLOGICAL_STERIMOL_NEG 10
#define TOPOLOGICAL_STERIMOL_UNSAT 11
#define TOPOLOGICAL_STERIMOL_CIGAR 12
#define TOPOLOGICAL_STERIMOL_AVE_DIST_10 13
#define TOPOLOGICAL_STERIMOL_ATOMS_AT_1 14
#define TOPOLOGICAL_STERIMOL_ATOMS_AT_2 15
#define TOPOLOGICAL_STERIMOL_ATOMS_AT_3 16
#define TOPOLOGICAL_STERIMOL_MIN_AROM 17
#define TOPOLOGICAL_STERIMOL_MAX_AROM 18
#define TOPOLOGICAL_STERIMOL_MIN_HTRO 19
#define TOPOLOGICAL_STERIMOL_MAX_HTRO 20

#define NTOPOLOGICAL_STERIMOL 21

#define TOPOLOGICAL_STERIMOL_Q0 0
#define TOPOLOGICAL_STERIMOL_Q1 1
#define TOPOLOGICAL_STERIMOL_MAXQ 2
#define TOPOLOGICAL_STERIMOL_MINQ 3
#define TOPOLOGICAL_STERIMOL_QEXT 4
#define TOPOLOGICAL_STERIMOL_TABS 5

#define NTOPO_STERIMOL_CHARGE 6

class Topological_Sterimol
{
 private:
  int _descriptor[NTOPOLOGICAL_STERIMOL];

  charge_t _partial_charge_descriptor[NTOPO_STERIMOL_CHARGE];

  int _do_partial_charge_descriptors;

 public:
  Topological_Sterimol();

  void
  set_do_partial_charge_descriptors(int s)
  {
    _do_partial_charge_descriptors = s;
  }

  int write_descriptors(std::ostream&) const;

  int compute_descriptors(Molecule&, const int* don_acc, int*, atom_number_t,
                          atom_number_t);
};

#endif
