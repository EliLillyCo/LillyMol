#ifndef MOLECULE_TOOLS_DICER_LIB_H_
#define MOLECULE_TOOLS_DICER_LIB_H_

#include <optional>

#include "Molecule_Lib/atom.h"
#include "Molecule_Lib/iwmtypes.h"
#include "Molecule_Lib/molecule.h"

namespace dicer_lib {

/*
  We can store various things with the user specified void pointer of the parent.
  This gets passed down to the fragments
*/

struct USPVPTR
{
  public:
    atom_number_t _atom_number_in_parent;
    uint32_t _atom_type_in_parent;

  public:
    USPVPTR();

    void set_atom_number_in_parent(const atom_number_t s) { _atom_number_in_parent = s;}
    atom_number_t atom_number_in_parent() const { return _atom_number_in_parent;}

    void set_atom_type_in_parent(const uint32_t s) { _atom_type_in_parent = s;}
    uint32_t atom_type_in_parent() const { return _atom_type_in_parent;}

};  // USPVPTR

// Given that `a` has a user_specified_void_ptr of USPVPTR, return the atom
// number in the parent molecule.
std::optional<atom_number_t> InitialAtomNumber(const Atom* a);

// Given a set of atom types, assign the user_specified_void_ptr attribute for
// each atom in `m`.
int
initialise_atom_pointers(Molecule & m,
                         const uint32_t * atype,
                         USPVPTR * uspvptr);

// Cut a bond, replacing a departing neighbour at each specified chiral
// endpoint with a terminal *. Caps have no parent identity. Returns false,
// without changing m, if a chiral endpoint already has a cap or needs a non-single cap bond.
bool BreakBondWithChiralDummies(Molecule& m, atom_number_t a1, atom_number_t a2);

// Rebuild a parent partition with the same capping semantics as recursive cuts.
// xref has space for the original parent atoms only. False discards the pair.
bool CreateChiralDummySubsets(Molecule& parent, const int* subset,
                             Molecule& fragment, Molecule& complement, int* xref);

// Return true if `b` is the bond leading to a CF3 or t-butyl group.
int is_cf3_like(const Molecule & m, const Bond & b);
int is_amide(const Molecule & m, atom_number_t a1, atom_number_t a2);

} // namespace dicer_lib

#endif  // MOLECULE_TOOLS_DICER_LIB_H_
