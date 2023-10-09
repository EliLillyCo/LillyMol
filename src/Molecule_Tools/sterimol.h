#ifndef MOLECULE_TOOLS_STERIMOL_H
#define MOLECULE_TOOLS_STERIMOL_H

#include "Foundational/accumulator/accumulator.h"

#include "Molecule_Lib/iwmtypes.h"

class Molecule;

typedef float sterimol_t;
class Atom;

#define STERIMOL_LONGEST_DISTANCE 0
#define STERIMOL_LONGEST_PERP_DIST 1
#define STERIMOL_YCROSS 2
#define STERIMOL_ZCROSS 3
#define STERIMOL_AVE_Y 4
#define STERIMOL_AVE_Z 5
#define STERIMOL_YZPROD 6
#define STERIMOL_Y_ANGLE 7
#define STERIMOL_Z_ANGLE 8
#define STERIMOL_ANGLE_RATIO 9
#define STERIMOL_LARGE_CONE 10
#define STERIMOL_OVERALL_CONE 11

#define NSTERIMOL 12

#define STERMIOL_SPATIAL_AVERAGE 0

#define NSTERIMOL_PARTIAL_CHARGE 1

class Sterimol
{
 private:
  sterimol_t _b[NSTERIMOL];

  charge_t _partial_charge_descriptor[NSTERIMOL_PARTIAL_CHARGE];

  // If we are looking for distances from a given atom, we encode those separately

  Accumulator<sterimol_t> _from_base_atom;

  int _do_partial_charge_descriptors;

  int _atoms_in_fragment;

  //  private functions

  int _compute_distances_from_a1(int matoms, const int* process, atom_number_t a1,
                                 const Atom** atom);

 public:
  Sterimol();

  void
  set_do_partial_charge_descriptors(int s)
  {
    _do_partial_charge_descriptors = s;
  }

  int
  do_partial_charge_descriptors() const
  {
    return _do_partial_charge_descriptors;
  }

  void reset();

  int copy_if_longer(const Sterimol&);

  sterimol_t operator[](int) const;
  sterimol_t& operator[](int);

  int write_descriptors(std::ostream&) const;

  int do_computation(Molecule&);
  int do_computation(Molecule&, const int* don_acc, const int*,
                     atom_number_t = INVALID_ATOM_NUMBER);
};

extern int sterimol(Molecule&, Sterimol&);

extern int write_sterimol_descriptor_headers(std::ostream&, const IWString&, int = 0);

#endif  // MOLECULE_TOOLS_STERIMOL_H
