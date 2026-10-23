#include "Depict/coords_2d.h"

#include <algorithm>
#include <cmath>
#include <vector>

// coordgen, built by src/build_linux.sh when BUILD_SCHRODINGER_2D is set.
// Included after the LillyMol headers because sketcherMinimizerMaths.h, pulled
// in from here, #defines BONDLENGTH and MACROCYCLE without guards.
#include "coordgen/sketcherMinimizer.h"
#include "coordgen/sketcherMinimizerAtom.h"
#include "coordgen/sketcherMinimizerBond.h"
#include "coordgen/sketcherMinimizerMaths.h"
#include "coordgen/sketcherMinimizerMolecule.h"

namespace depict {

namespace {

float
CoordgenPrecision(Precision p) {
  switch (p) {
    case Precision::kQuick:
      return SKETCHER_QUICK_PRECISION;
    case Precision::kBest:
      return SKETCHER_BEST_PRECISION;
    case Precision::kStandard:
    default:
      return SKETCHER_STANDARD_PRECISION;
  }
}

// A sketcherMinimizerMolecule does not own the atoms and bonds added to it -
// its destructor frees only the rings. Ownership passes to the
// sketcherMinimizer at initialize(), which frees all three in its destructor.
// So exactly one of these two things must happen, and this guard covers the
// window before initialize() where we are still the owner.
class MoleculeOwner {
  private:
    sketcherMinimizerMolecule* _mol;

  public:
    explicit MoleculeOwner(sketcherMinimizerMolecule* mol) : _mol(mol) {
    }

    ~MoleculeOwner() {
      if (_mol == nullptr) {
        return;
      }
      for (sketcherMinimizerAtom* a : _mol->_atoms) {
        delete a;
      }
      for (sketcherMinimizerBond* b : _mol->_bonds) {
        delete b;
      }
      delete _mol;
    }

    sketcherMinimizerMolecule*
    get() const {
      return _mol;
    }

    // Called once the sketcherMinimizer has taken over.
    void
    release() {
      _mol = nullptr;
    }
};

// coordgen wants an integer bond order. LillyMol stores aromatic bonds in
// Kekule form with an extra aromatic flag, so the ordinary accessors give the
// right answer. A bond that is aromatic but has no Kekule form is reported as
// none of single/double/triple, and falls through to 1.
int
BondOrder(const Bond& b) {
  if (b.is_single_bond()) {
    return 1;
  }
  if (b.is_double_bond()) {
    return 2;
  }
  if (b.is_triple_bond()) {
    return 3;
  }

  return 1;
}

/* One double bond whose cis/trans configuration the molecule states, written the
   way coordgen needs to hear it: the double bond a3==a4, one substituent on each
   end, and whether those two substituents belong on the same side.

      a1        a5          a1
        \      /              \
         a3==a4                a3==a4
                                     \
                                      a5
          cis                        trans
*/
struct CisTransConstraint {
  int bond_index = -1;  // into Molecule::bondi(), the double bond itself
  atom_number_t a1 = kInvalidAtomNumber;
  atom_number_t a3 = kInvalidAtomNumber;
  atom_number_t a4 = kInvalidAtomNumber;
  atom_number_t a5 = kInvalidAtomNumber;
  bool cis = false;
};

// The sense of a directional single bond, looking outward from `from`: +1 for
// up, -1 for down, 0 if the bond carries no direction. This mirrors the file
// static discern_directionality() in Molecule_Lib/path_scoring.cc, which is not
// exported.
//
// LillyMol holds E/Z the way SMILES writes it, as up/down flags on the single
// bonds flanking the double bond, and the meaning of a pair of them is a
// convention rather than something to be derived. The convention is fixed by
// Molecule::_discern_cis_trans_bond_from_depiction(), which reorders the
// substituents so that a1 and a5 are the pair on the same side and then treats
// a1_direction != a5_direction as a contradiction: two substituents at opposite
// ends of a double bond are cis exactly when these two senses are equal.
int
DirectionalSense(atom_number_t from, const Bond& b) {
  if (!b.is_directional()) {
    return 0;
  }

  const int sign = (b.a1() == from) ? 1 : -1;

  return b.is_directional_up() ? sign : -sign;
}

// Every double bond of `m` whose configuration is known. One substituent per end
// is enough to pin the bond down, so the first directional one found is used;
// which one it is does not matter, because the flags on a well formed molecule
// agree with each other.
void
CollectCisTransConstraints(const Molecule& m,
                           std::vector<CisTransConstraint>& constraints) {
  const int nedges = m.nedges();

  for (int i = 0; i < nedges; ++i) {
    const Bond* b = m.bondi(i);
    if (!b->is_double_bond() || !b->part_of_cis_trans_grouping()) {
      continue;
    }

    const atom_number_t a3 = b->a1();
    const atom_number_t a4 = b->a2();

    atom_number_t a1 = kInvalidAtomNumber;
    atom_number_t a5 = kInvalidAtomNumber;
    int d1 = 0;
    int d5 = 0;

    for (const Bond* b : m[a3]) {
      const atom_number_t o = b->other(a3);
      if (o == a4) {
        continue;
      }
      const int d = DirectionalSense(a3, *b);
      if (d != 0) {
        a1 = o;
        d1 = d;
        break;
      }
    }

    for (const Bond* b : m[a4]) {
      const atom_number_t o = b->other(a4);
      if (o == a3) {
        continue;
      }
      const int d = DirectionalSense(a4, *b);
      if (d != 0) {
        a5 = o;
        d5 = d;
        break;
      }
    }

    // A double bond flagged as part of a cis/trans grouping but with no
    // directional bond on one of its ends says nothing about configuration.
    if (d1 == 0 || d5 == 0) {
      continue;
    }

    CisTransConstraint c;
    c.bond_index = i;
    c.a1 = a1;
    c.a3 = a3;
    c.a4 = a4;
    c.a5 = a5;
    c.cis = (d1 == d5);
    constraints.push_back(c);
  }
}

// Whether the finished layout draws `c` the way it was asked to. Measured from
// the coordinates rather than asked of coordgen, because what matters is what a
// reader of the molfile will conclude, and that is a property of the coordinates
// alone.
bool
ConstraintSatisfied(const CisTransConstraint& c, const std::vector<float>& x,
                    const std::vector<float>& y) {
  const float dx = x[c.a4] - x[c.a3];
  const float dy = y[c.a4] - y[c.a3];

  // Signed areas, so the sign says which side of the double bond axis each
  // substituent is on.
  const float s1 = dx * (y[c.a1] - y[c.a3]) - dy * (x[c.a1] - x[c.a3]);
  const float s5 = dx * (y[c.a5] - y[c.a3]) - dy * (x[c.a5] - x[c.a3]);

  const bool same_side = (s1 > 0.0f) == (s5 > 0.0f);

  return same_side == c.cis;
}

// The distance between the closest pair of atoms not bonded to each other, or
// -1 if there is no such pair. Quadratic, so only called when a caller asks
// for the measurement.
float
ClosestNonBondedApproach(const Molecule& m, const std::vector<float>& x,
                         const std::vector<float>& y) {
  const int matoms = m.natoms();

  float closest = -1.0f;
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      if (m.are_bonded(i, j)) {
        continue;
      }
      const float dx = x[i] - x[j];
      const float dy = y[i] - y[j];
      const float d = std::sqrt(dx * dx + dy * dy);
      if (closest < 0.0f || d < closest) {
        closest = d;
      }
    }
  }

  return closest;
}

// `result` may be null, for the overloads that do not want the measurements.
int
Generate2DCoordinates(Molecule& m, const Coords2DOptions& opts,
                      Coords2DResult* result) {
  const int matoms = m.natoms();

  // coordgen's sanity check rejects an empty structure, and prints a warning
  // to stderr on the way out. Stop here so nothing is written to stderr.
  if (matoms == 0) {
    return 0;
  }

  MoleculeOwner owner(new sketcherMinimizerMolecule());

  // Held in LillyMol atom order. coordgen is free to reorder its own _atoms,
  // so read the results back through these rather than through the molecule.
  std::vector<sketcherMinimizerAtom*> atoms;
  atoms.reserve(matoms);

  for (int i = 0; i < matoms; ++i) {
    sketcherMinimizerAtom* a = owner.get()->addNewAtom();
    a->setAtomicNumber(m.atomic_number(i));
    a->charge = m.formal_charge(i);
    atoms.push_back(a);
  }

  // Held in LillyMol bond order, for the same reason as `atoms`.
  const int nedges = m.nedges();
  std::vector<sketcherMinimizerBond*> bonds;
  bonds.reserve(nedges);

  for (int i = 0; i < nedges; ++i) {
    const Bond* b = m.bondi(i);
    sketcherMinimizerBond* nb =
        owner.get()->addNewBond(atoms[b->a1()], atoms[b->a2()]);
    nb->setBondOrder(BondOrder(*b));
    bonds.push_back(nb);
  }

  std::vector<CisTransConstraint> cis_trans;
  if (opts.honour_cis_trans) {
    CollectCisTransConstraints(m, cis_trans);
  }

  sketcherMinimizer minimizer(CoordgenPrecision(opts.precision));
  minimizer.setSkipMinimization(opts.skip_minimization);
  minimizer.setEvenAngles(opts.even_angles);

  minimizer.initialize(owner.get());
  owner.release();

  // After initialize(), not before. setAbsoluteStereoFromStereoInfo() turns the
  // relative statement "these two substituents are on the same side" into
  // coordgen's absolute isZ flag, which is expressed against the CIP-senior
  // substituent at each end. Working that out needs each atom's neighbour list
  // and the ring perception, and initialize() is what builds both -
  // addNewBond() only records the bond. Called earlier it silently does nothing,
  // leaving isZ at its default of false, i.e. trans.
  //
  // Double bonds the molecule says nothing about are deliberately left alone
  // rather than being marked unspecified: unspecified sets coordgen's
  // m_ignoreZE, which promotes the bond to a rotatable inter-fragment bond and
  // changes the layout of every molecule that has a double bond without E/Z.
  for (const CisTransConstraint& c : cis_trans) {
    sketcherMinimizerBondStereoInfo info;
    info.atom1 = atoms[c.a1];
    info.atom2 = atoms[c.a5];
    info.stereo = c.cis ? sketcherMinimizerBondStereoInfo::cis
                        : sketcherMinimizerBondStereoInfo::trans;
    bonds[c.bond_index]->setStereoChemistry(info);
    bonds[c.bond_index]->setAbsoluteStereoFromStereoInfo();
  }

  // The return value is deliberately discarded. It looks like it says whether
  // the pose came out free of clashes, but it is captured inside
  // CoordgenMinimizer::avoidClashesOfMolecule() from flipFragments(), before
  // the remedial avoidTerminalClashes() and minimizeMolecule() passes that
  // follow it. So it means "the first attempt needed no cleanup", and is false
  // for structures as ordinary as aspirin. It is also true whenever the sanity
  // check rejected the structure and nothing was laid out at all. Success is
  // determined by checking the atoms below instead.
  minimizer.runGenerateCoordinates();

  const float scale = opts.bond_length / static_cast<float>(BONDLENGTH);

  std::vector<float> x(matoms);
  std::vector<float> y(matoms);

  for (int i = 0; i < matoms; ++i) {
    // Set only by setCoordinates, so this is what distinguishes "laid out"
    // from "the sanity check rejected this structure".
    if (!atoms[i]->coordinatesSet) {
      return 0;
    }
    x[i] = atoms[i]->getCoordinates().x() * scale;
    y[i] = atoms[i]->getCoordinates().y() * scale;

    // coordgen occasionally produces NaN on heavily bridged polycyclics, while
    // still reporting the coordinates as set, so this has to be checked
    // separately. Writing NaN into a molfile produces a record no toolkit can
    // read, so treat it as a failure to lay the molecule out. See the header for
    // how often this and the collapse below happen.
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) {
      return 0;
    }
  }

  // coordgen also occasionally collapses part of a structure onto a single
  // point, again on heavily bridged polycyclics, and in the worst case observed
  // put 34 atoms on 8 positions. Two bonded atoms at the same coordinates
  // cannot be drawn, so this is a failed layout rather than a poor one. It is
  // reported as set just as the NaN case is. Checked against the requested bond
  // length so the test means the same thing at any scale; the threshold is far
  // below any separation a real layout produces, so a usable layout is never
  // rejected here.
  const float coincident = 0.01f * opts.bond_length;
  for (int i = 0; i < nedges; ++i) {
    const Bond* b = m.bondi(i);
    const float dx = x[b->a1()] - x[b->a2()];
    const float dy = y[b->a1()] - y[b->a2()];
    if (std::sqrt(dx * dx + dy * dy) < coincident) {
      return 0;
    }
  }

  if (opts.centre) {
    float xmin = x[0];
    float xmax = x[0];
    float ymin = y[0];
    float ymax = y[0];
    for (int i = 1; i < matoms; ++i) {
      xmin = std::min(xmin, x[i]);
      xmax = std::max(xmax, x[i]);
      ymin = std::min(ymin, y[i]);
      ymax = std::max(ymax, y[i]);
    }

    const float xshift = 0.5f * (xmin + xmax);
    const float yshift = 0.5f * (ymin + ymax);
    for (int i = 0; i < matoms; ++i) {
      x[i] -= xshift;
      y[i] -= yshift;
    }
  }

  if (result != nullptr) {
    result->closest_nonbonded_approach = ClosestNonBondedApproach(m, x, y);
    result->cis_trans_bonds = static_cast<int>(cis_trans.size());
    for (const CisTransConstraint& c : cis_trans) {
      if (ConstraintSatisfied(c, x, y)) {
        ++result->cis_trans_honoured;
      }
    }
  }

  // Committed only now, so a failure above leaves `m` untouched.
  for (int i = 0; i < matoms; ++i) {
    m.setxyz(i, x[i], y[i], 0.0f);
  }

  return 1;
}

}  // namespace

int
Generate2DCoordinates(Molecule& m, const Coords2DOptions& opts,
                      Coords2DResult& result) {
  return Generate2DCoordinates(m, opts, &result);
}

int
Generate2DCoordinates(Molecule& m, const Coords2DOptions& opts) {
  return Generate2DCoordinates(m, opts, nullptr);
}

int
Generate2DCoordinates(Molecule& m) {
  static const Coords2DOptions opts;

  return Generate2DCoordinates(m, opts, nullptr);
}

}  // namespace depict
