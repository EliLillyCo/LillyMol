#ifndef DEPICT_WEDGE_2D_H
#define DEPICT_WEDGE_2D_H

// Assignment of wedge and hash bonds to a 2D depiction, so that the chirality
// LillyMol holds internally survives being written to a molfile.
//
// This is needed because LillyMol writes chirality as the MDL atom parity
// field, which most other software ignores - RDKit reports such an atom as
// having no assigned stereochemistry. What toolkits do read is the wedge/hash
// flag in the bond block, interpreted against the 2D coordinates. Without
// wedges a depiction is stereochemically silent no matter how correct the
// molecule is inside LillyMol.
//
// Nothing here needs the layout engine, only coordinates that are already
// present, so it does not depend on coordgen.

#include "Molecule_Lib/molecule.h"

namespace depict {

// What happened, for callers that want to report or test it.
class WedgeResult {
  public:
    // Chiral centres found in the molecule.
    int chiral_centres = 0;

    // Centres given a wedge or hash bond that reproduces their chirality.
    int wedged = 0;

    // Centres where chirality perception saw no stereocentre at all, whichever
    // way a wedge was tried. These are normal and not a problem: a SMILES can
    // mark an atom as chiral when it is not stereogenic, as in
    // 1,4-dimethylcyclohexane, whose two ring branches are identical. There is
    // no stereochemistry to draw, so none is drawn.
    int not_stereogenic = 0;

    // Centres that are stereogenic but could not be given a wedge that
    // reproduces their chirality. Their stereochemistry will not reach software
    // that reads wedges. See AssignWedgeBonds for when this happens.
    int unresolved = 0;
};

// Give every chiral centre in `m` a wedge or hash bond, chosen so that reading
// the 2D coordinates together with the wedge recovers the chirality `m` already
// holds. Every wedge and hash bond already on the molecule is cleared first, so
// this is idempotent, and a molecule read from a molfile and then re-laid out
// does not keep wedges that were drawn against the old coordinates. Cis/trans
// directionality on single bonds is left alone, since that is how E/Z is held.
//
// `m` must already have 2D coordinates; with all atoms at the origin there is
// no geometry to reason about and every centre is left unresolved.
//
// Returns the number of centres wedged. A centre is left unresolved, rather
// than being given a wedge that might state the wrong stereochemistry, when:
//   - it has no single bond to a neighbour that can carry a wedge;
//   - its 2D geometry is degenerate, e.g. two neighbours superimposed;
//   - the verification pass below could not confirm the assignment.
//
// The direction of each wedge is not derived from a hand-written parity rule.
// A wedge says "this neighbour is toward the viewer", so the candidate is
// tested by lifting that neighbour out of the plane, asking Molecule to derive
// chirality from the resulting three dimensional geometry, and comparing the
// handedness with what `m` started with; the direction that agrees is the one
// used. That leans on Molecule::discern_chirality_from_3d_structure(), which is
// well exercised by every 3D file LillyMol reads, and deliberately avoids
// Molecule::_discern_chirality_from_wedge_bond_4(), which its own source
// comment describes as buggy. Every assignment is then verified the same way
// before being kept, so a wrong wedge is reported as unresolved rather than
// written out.
int AssignWedgeBonds(Molecule& m, WedgeResult& result);

int AssignWedgeBonds(Molecule& m);

}  // namespace depict

#endif  // DEPICT_WEDGE_2D_H
