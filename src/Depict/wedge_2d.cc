#include "Depict/wedge_2d.h"

#include <algorithm>
#include <cmath>

#include "Molecule_Lib/chiral_centre.h"
#include "Molecule_Lib/molecule.h"

namespace depict {

namespace {

// How far out of the plane a candidate neighbour is lifted when working out
// which direction a wedge has to point. Any value clearly larger than the
// rounding in the coordinates does, since only the sign of the resulting
// handedness is used. A bond length keeps the probe geometry realistic, which
// matters because chirality perception looks at angles.
constexpr float kLift = 1.0f;

// The four connections of a chiral centre, in the order the class stores them.
// Entries are atom numbers, or the sentinels for an implicit hydrogen or a lone
// pair.
void
ConnectionsInOrder(const Chiral_Centre& c, int* dest) {
  dest[0] = c.top_front();
  dest[1] = c.top_back();
  dest[2] = c.left_down();
  dest[3] = c.right_down();
}

// Whether two chiral centres on the same atom describe the same handedness.
//
// Both list the same four connections, so one ordering is a permutation of the
// other, and the two describe the same handedness exactly when that permutation
// is even. This is pure combinatorics: it needs no knowledge of what top_front
// and left_down mean geometrically, which is the point - that convention is
// easy to get wrong and nothing here has to know it.
//
// Returns -1 if the two do not list the same connections, so the caller can
// treat it as unknown rather than as a mismatch.
int
SameHandedness(const Chiral_Centre& c1, const Chiral_Centre& c2) {
  int a[4], b[4];
  ConnectionsInOrder(c1, a);
  ConnectionsInOrder(c2, b);

  // Where each of b's entries sits in a.
  int perm[4];
  for (int i = 0; i < 4; ++i) {
    perm[i] = -1;
    for (int j = 0; j < 4; ++j) {
      if (b[i] == a[j]) {
        // A repeated value would make the permutation ambiguous. Two lone pairs
        // or two implicit hydrogens on one centre are not a chiral centre we
        // can depict, so bail out.
        if (perm[i] >= 0) {
          return -1;
        }
        perm[i] = j;
      }
    }
    if (perm[i] < 0) {
      return -1;
    }
  }

  // Parity by counting inversions.
  int inversions = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      if (perm[i] > perm[j]) {
        ++inversions;
      }
    }
  }

  return 0 == inversions % 2;
}

// Clear any wedge on the bond between `a1` and `a2` without disturbing the
// cis/trans directionality that single bonds next to a double bond carry, since
// that is how E/Z is held and is not ours to discard.
void
ClearWedge(Molecule& m, atom_number_t a1, atom_number_t a2) {
  const Bond* b = m.bond_between_atoms(a1, a2);
  if (nullptr == b || !b->is_wedge_any()) {
    return;
  }

  const int was_up = b->is_directional_up();
  const int was_down = b->is_directional_down();

  m.set_wedge_bond_between_atoms(a1, a2, 0);

  if (was_up) {
    const_cast<Bond*>(m.bond_between_atoms(a1, a2))->set_directional_up(a1, a2);
  } else if (was_down) {
    const_cast<Bond*>(m.bond_between_atoms(a1, a2))->set_directional_down(a1, a2);
  }
}

// How suitable `nbr` is as the atom a wedge at `centre` points to. Higher is
// better; a negative score means unusable.
//
// The preferences are the usual drawing conventions: put the wedge on a bond
// that unambiguously belongs to this centre and to nothing else. A terminal
// neighbour is ideal because the wide end of the wedge has empty space around
// it, an explicit hydrogen most of all since that is what a chemist expects to
// see. Ring bonds are avoided because a wedge drawn along a ring bond reads as
// belonging to either of its atoms.
int
WedgeCandidateScore(Molecule& m, atom_number_t centre, atom_number_t nbr) {
  const Bond* b = m.bond_between_atoms(centre, nbr);
  if (nullptr == b || !b->is_single_bond()) {
    return -1;
  }

  // A bond that already carries cis/trans information must keep it.
  if (b->is_directional() || b->part_of_cis_trans_grouping()) {
    return -1;
  }

  // Superimposed atoms give no direction to reason about.
  if (m.distance_between_atoms(centre, nbr) < 1.0e-4) {
    return -1;
  }

  int score = 0;

  // A wedge is drawn narrow end first, and the narrow end has to sit on the
  // stereocentre - that is what makes it a statement about this atom. The
  // molfile bond block writes each bond in its stored order, so a bond already
  // stored pointing away from the centre needs no reorienting. Strongly
  // preferred over one that does, hence a bonus larger than all the others
  // combined.
  if (b->a1() == centre) {
    score += 100;
  }

  if (1 == m.ncon(nbr)) {
    score += 40;
    if (1 == m.atomic_number(nbr)) {
      score += 20;
    }
  }

  if (!b->nrings()) {
    score += 10;
  }

  // Another centre's wedge would be reading the same bond.
  if (nullptr == m.chiral_centre_at_atom(nbr)) {
    score += 5;
  }

  return score;
}

// The handedness Molecule derives at `zatom` when `lifted` is pulled out of the
// plane by `lift`. Returns the derived centre, or nullptr when nothing could be
// derived. `probe` is left holding the derived chirality.
const Chiral_Centre*
HandednessWithLift(Molecule& probe, atom_number_t zatom, atom_number_t lifted,
                   float lift) {
  probe.remove_all_chiral_centres();
  probe.setz(lifted, lift);

  probe.discern_chirality_from_3d_structure();

  return probe.chiral_centre_at_atom(zatom);
}

}  // namespace

int
AssignWedgeBonds(Molecule& m, WedgeResult& result) {
  result = WedgeResult();

  // Any wedge already on the molecule was drawn against different coordinates -
  // typically those of the molfile it was read from - and once the layout has
  // been replaced it is no longer a true statement about anything. Clearing
  // them all means the wedges that come out are exactly the ones put there
  // here, whatever the input carried.
  for (int i = 0; i < m.nedges(); ++i) {
    const Bond* b = m.bondi(i);
    if (b->is_wedge_any()) {
      ClearWedge(m, b->a1(), b->a2());
    }
  }

  const int ncentres = m.chiral_centres();
  result.chiral_centres = ncentres;
  if (0 == ncentres) {
    return 0;
  }

  // Work from a snapshot: the loop below builds probe molecules whose chirality
  // is removed and re-derived, and it needs the original handedness to compare
  // against.
  const Molecule original(m);

  // Bonds already promised to a centre, so two centres do not wedge the same
  // bond and leave the drawing ambiguous.
  resizable_array<const Bond*> claimed;

  for (int i = 0; i < ncentres; ++i) {
    const Chiral_Centre* c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);
    if (nullptr == c || !c->chirality_known()) {
      ++result.unresolved;
      continue;
    }

    const atom_number_t centre = c->a();

    // Pick the neighbour the wedge will point at.
    atom_number_t best = INVALID_ATOM_NUMBER;
    int best_score = -1;
    for (int j = 0; j < m.ncon(centre); ++j) {
      const atom_number_t nbr = m.other(centre, j);
      const Bond* b = m.bond_between_atoms(centre, nbr);
      if (claimed.contains(b)) {
        continue;
      }
      const int score = WedgeCandidateScore(m, centre, nbr);
      if (score > best_score) {
        best_score = score;
        best = nbr;
      }
    }

    if (best_score < 0) {
      ++result.unresolved;
      continue;
    }

    // Which way does it have to point? A wedge means `best` is toward the
    // viewer, so lift it and see which sign reproduces the original handedness.
    int direction = 0;
    int perceived_a_centre = 0;
    for (const float lift : {kLift, -kLift}) {
      Molecule probe(original);
      const Chiral_Centre* derived = HandednessWithLift(probe, centre, best, lift);
      if (nullptr == derived) {
        continue;
      }
      ++perceived_a_centre;
      const int same = SameHandedness(*c, *derived);
      if (1 == same) {
        direction = (lift > 0.0f) ? 1 : -1;
        break;
      }
    }

    if (0 == direction) {
      // Nothing was perceived either way, so there is no stereochemistry here to
      // express, as opposed to stereochemistry we failed to express.
      if (0 == perceived_a_centre) {
        ++result.not_stereogenic;
      } else {
        ++result.unresolved;
      }
      continue;
    }

    // Point the bond away from the centre, so the wedge that is written has its
    // narrow end there. Reversing a bond's stored orientation does not change
    // what is bonded to what - each atom's list still holds this same bond - it
    // only decides which atom the molfile writes first. Without this the wedge
    // is emitted starting at the neighbour, and a reader that takes the first
    // atom as the stereocentre, as RDKit does, discards it.
    Bond* b = const_cast<Bond*>(m.bond_between_atoms(centre, best));
    if (b->a1() != centre) {
      b->set_a1a2(centre, best);
    }

    if (!m.set_wedge_bond_between_atoms(centre, best, direction)) {
      ++result.unresolved;
      continue;
    }

    claimed.add(m.bond_between_atoms(centre, best));
    ++result.wedged;
  }

  return result.wedged;
}

int
AssignWedgeBonds(Molecule& m) {
  WedgeResult notused;
  return AssignWedgeBonds(m, notused);
}

}  // namespace depict
