#ifndef DEPICT_COORDS_2D_H
#define DEPICT_COORDS_2D_H

// Generation of 2D coordinates suitable for depiction.
//
// The layout engine is Schrodinger's coordgen (BSD-3-Clause), an optional
// dependency built by src/build_linux.sh when BUILD_SCHRODINGER_2D is set. None
// of it appears in this header: coordgen's headers #define unguarded macros such
// as BONDLENGTH and MACROCYCLE, so they are confined to coords_2d.cc and callers
// need neither the include path nor the macros.

#include "Molecule_Lib/molecule.h"

namespace depict {

// How hard the layout engine works to resolve clashes. Higher settings cost
// more time and only differ on molecules that are hard to lay out; for most
// drug-like molecules all three give the same answer.
enum class Precision {
  kQuick,
  kStandard,
  kBest
};

class Coords2DOptions {
  public:
    // Length assigned to a bond in the generated layout. 1.5 is the length
    // used by MDL molfiles and is what most viewers expect.
    float bond_length = 1.5f;

    Precision precision = Precision::kStandard;

    // Skip the force field refinement that follows the initial construction.
    // Faster, but leaves clashes that would otherwise be resolved.
    bool skip_minimization = false;

    // Spread substituents evenly around an atom rather than preferring the
    // idealised angles of the parent ring or chain.
    bool even_angles = false;

    // Translate the result so its centre of geometry is the origin.
    bool centre = true;

    // Lay double bonds out so that the cis/trans configuration the molecule
    // already holds is what the drawing says. A molfile has no field for E/Z -
    // it lives only in the coordinates - so without this the layout engine
    // draws every double bond trans, silently turning Z into E. Turn it off
    // only to reproduce the raw engine layout.
    bool honour_cis_trans = true;
};

// Measurements of the layout that was produced, for callers wanting to know
// how crowded it is - a renderer might use this to flag a depiction, or to
// decide how much room atom labels have. Only filled in when generation
// succeeded.
//
// These are measured here rather than taken from the layout engine.
// sketcherMinimizer::runGenerateCoordinates() does return a bool that looks
// like it reports a clean pose, but it is the verdict of an intermediate stage,
// taken before the remedial clash-avoidance and minimization passes run. It is
// therefore false for aspirin and for hexadecane, whose final layouts are
// perfectly clean, and cannot be used as a quality signal.
class Coords2DResult {
  public:
    // Distance between the closest pair of atoms that are not bonded to each
    // other, in the same units as Coords2DOptions::bond_length. Small values
    // mean a crowded drawing: atoms overlapping their neighbours' labels.
    // Comfortable layouts are around a bond length or more; below about a third
    // of one the drawing is hard to read. Negative when there is no such pair,
    // i.e. fewer than three atoms.
    float closest_nonbonded_approach = -1.0f;

    // Double bonds whose cis/trans configuration the molecule stated, and which
    // therefore constrained the layout. Zero when
    // Coords2DOptions::honour_cis_trans is false, or when the molecule holds no
    // E/Z - note that LillyMol only holds E/Z from a SMILES, or from a molfile
    // read with -i dctb.
    int cis_trans_bonds = 0;

    // Of those, the ones the finished layout actually draws the right way round,
    // measured from the coordinates. The rest have been drawn as the other
    // isomer and will be read back as such; see Generate2DCoordinates for when
    // that happens.
    int cis_trans_honoured = 0;
};

// Replace the coordinates of `m` with a generated 2D layout. Every atom is
// given z == 0. Connectivity, atom order and atom count are unchanged.
//
// Returns 1 on success and 0 if no usable layout was produced, in which case
// the coordinates of `m` are left untouched. A success is guaranteed to have
// finite coordinates and no two bonded atoms at the same position, so the
// result can always be written out and read back. Failure is rare, and means
// one of:
//   - an empty molecule;
//   - a single connected fragment containing more than 40 rings, a coordgen
//     limit;
//   - coordgen produced NaN or collapsed part of the structure onto one point.
//     Only seen on heavily bridged polycyclics, and it reports the coordinates
//     as set in both cases, so these are checked for explicitly. Over a 39823
//     molecule structure-mutation set 47 molecules failed this way; over 1893
//     ChEMBL molecules, none.
//
// Note that a success is not a promise of a good layout, only of a drawable
// one. Bond lengths are uniform for ordinary molecules but cannot be for
// bridged ring systems, where 2D drawings are distorted by nature. Use the
// overload taking a Coords2DResult to find out how crowded the result is.
//
// Explicit hydrogens are laid out as atoms in their own right, which is
// usually not what a depiction wants. Callers that want implicit hydrogens
// should remove them before calling.
//
// Aromatic bonds are passed to the layout engine in their Kekule form. A
// molecule for which no Kekule form could be found still lays out, with those
// bonds treated as single, which affects the layout only marginally.
//
// Cis/trans double bonds are laid out to match the configuration `m` holds,
// which is not something the layout engine does on its own - left to itself it
// draws every double bond trans. A bond it cannot honour is left as it laid it
// out rather than being reported as a failure, since the rest of the layout is
// still correct; the count is in Coords2DResult. That happens for a double bond
// inside a ring of eight atoms or fewer, where the ring geometry decides the
// configuration and there is nothing to choose, and for a double bond whose two
// ends carry substituents the engine considers equivalent.
//
// Over the 1020 molecules of test/gfp_naive_bayesian/case_1/in/train.smi that
// state a configuration, laying out and reading back with RDKit gives the stated
// configuration for all of them; with honour_cis_trans off, 370 of them come
// back as the other isomer. A further 16 molecules have a double bond whose
// configuration they do not state and which the layout necessarily gives one to
// - unavoidable, since a drawing has to put the substituents somewhere.
//
// Note that LillyMol only holds E/Z if it was told: from a SMILES always, but
// from a 2D molfile only when it was read with `-i dctb`, which is off by
// default. A molecule read from a molfile without it has no configuration to
// honour and nothing here will invent one.
int Generate2DCoordinates(Molecule& m, const Coords2DOptions& opts,
                          Coords2DResult& result);

// Overloads for callers that do not want the measurements. They are cheaper:
// the closest approach is quadratic in the atom count, so it is only computed
// when asked for.
int Generate2DCoordinates(Molecule& m, const Coords2DOptions& opts);
int Generate2DCoordinates(Molecule& m);

}  // namespace depict

#endif  // DEPICT_COORDS_2D_H
