#ifndef MOLECULE_LIB_GEOMETRIC_CONSTRAINTS_H
#define MOLECULE_LIB_GEOMETRIC_CONSTRAINTS_H
// Geometry matching conditions.

#include "molecule.h"
#include "space_vector.h"
#include "Molecule_Lib/geometric_constraints.pb.h"

namespace geometric_constraints {

class AllowedRange {
  private:
    float _min_value;
    float _max_value;
  public:
    AllowedRange();

    int BuildFromProto(const GeometricConstraints::Range& proto);

    // Is this initialised?
    int Active() const {
      return _min_value <= _max_value;
    }

    int set_range(float d1, float d2) {
      _min_value = d1;
      _max_value = d2;
      return 1;
    }

    int Matches(float x) const {
      return x >= _min_value && x <= _max_value;
    }

    friend std::ostream& operator<< (std::ostream& output, const AllowedRange& arange);
};

// All geometric constraints consist of a set of atom numbers and an allowed range
// for the quantity defined by those atoms.
class ConstraintBaseClass {
  protected:
    // All constraints have 2 or more atoms.
    resizable_array<int> _atoms;
    // All constraints have a range of allowed values.
    AllowedRange _allowed_range;
  public:
    int SetRange(float d1, float d2) {
      return _allowed_range.set_range(d1, d2);
    }

    // Add our atom numbers to `atoms_present` if not already there.
    int AtomNumbersPresent(resizable_array<int>& atoms_present) const;

    friend std::ostream& operator<< (std::ostream& output, const ConstraintBaseClass& constraint);

    // The range must be set and the atoms unique. Expensive.
    int IsValid() const;

    // Inexpensive check that does not ensure internal validity.
    int Active() const { return _allowed_range.Active();}

    int SetAtoms(int s1, int s2);
    int SetAtoms(int s1, int s2, int s3);
    int SetAtoms(int s1, int s2, int s3, int s4);

    // Mostly for testing and debugging.
    const resizable_array<int> Atoms() const {
      return _atoms;
    }
};

// this is instantiated from the corresponding Message in the proto.
// Note, not currently implemented.
struct LineOfSight {
  vdw::VdwType vdw_type = vdw::VdwType::kUndefined;
  float delta = 0.0f;

  public:
    int BuildFromProto(const GeometricConstraints::LineOfSight& proto);

    int Active() const {
      return vdw_type != vdw::VdwType::kUndefined;
    }
};

// Classes for matching geometric constraints. They are all very similar.
// All contain a Matches method that has two signatures.
// The two Matches function differ in how the atom numbers (_a1, _a2...) are interpreted.
// Without an embedding, the atom numbers are just atom numbers in the molecule.
// With an embedding, they are matched atom numbers from the embedding.
// Note too that the Active methods just check the first of the atom numbers for
// validity.
// The IsValid method provides more expensive checking.

// Class specifying a geometric distance between atoms constraint.
class DistanceConstraint : public ConstraintBaseClass {
  private:
    float _tube_radius;

    int MatchesWithEmptyTube(Molecule& m, atom_number_t a1, atom_number_t a2) const;

  public:
    DistanceConstraint();

    int BuildFromProto(const GeometricConstraints::Distance& proto);

    int Matches(Molecule& m);
    int Matches(Molecule& m, const Set_of_Atoms& embedding);

    friend std::ostream& operator<< (std::ostream& output, const DistanceConstraint& constraint);
};

// Class specifying constaints on a bond angle.
class BondAngleConstraint : public ConstraintBaseClass{
  private:

  public:
    int BuildFromProto(const GeometricConstraints::BondAngle& proto);

    int set_range_in_degrees(float mn, float mx) {
      return _allowed_range.set_range(iwmisc::Deg2Rad(mn), iwmisc::Deg2Rad(mx));
    }
    int set_range_in_radians(float mn, float mx) {
      return _allowed_range.set_range(mn, mx);
    }

    int Matches(const Molecule& m);
    int Matches(const Molecule& m, const Set_of_Atoms& embedding);

    friend std::ostream& operator<< (std::ostream& output, const BondAngleConstraint& constraint);
};

// Class specifying constaints on a torsion angle.
class TorsionAngleConstraint : public ConstraintBaseClass {
  private:

  public:
    int BuildFromProto(const GeometricConstraints::TorsionAngle& proto);

    int set_range_in_degrees(float mn, float mx) {
      return _allowed_range.set_range(iwmisc::Deg2Rad(mn), iwmisc::Deg2Rad(mx));
    }
    int set_range_in_radians(float mn, float mx) {
      return _allowed_range.set_range(mn, mx);
    }

    int Matches(const Molecule& m);
    int Matches(const Molecule& m, const Set_of_Atoms& embedding);

    friend std::ostream& operator<< (std::ostream& output, const TorsionAngleConstraint& constraint);
};

class PositionalConstraint : public ConstraintBaseClass {
  private:
    GeometricConstraints::SpatialConstraintMatchType _match_type;

    Space_Vector<double> _coords;

  // private functions
    template <typename P> int BuildFromProto(const P& proto);
  
    int MatchesAverage(const Molecule& m, const Set_of_Atoms& embedding) const;
    int MatchesAll(const Molecule& m, const Set_of_Atoms& embedding) const;
    int MatchesAny(const Molecule& m, const Set_of_Atoms& embedding) const;
  public:
    PositionalConstraint();

    int BuildFromProto(const GeometricConstraints::PositiveSpatialConstraint& proto);
    int BuildFromProto(const GeometricConstraints::NegativeSpatialConstraint& proto);

    int Matches(const Molecule& m) const;
    int Matches(const Molecule& m, const Set_of_Atoms& embedding) const;

    friend std::ostream& operator<< (std::ostream& output, const PositionalConstraint& constraint);
};

// A common situation is to have an arbitrary number of constraints of different kinds.
// The _number_to_match parameter governs matching.
// it specifies the number of queries must match - set to the sum of
// distance + angle + torsion conditions to force all to match.
// Set to 1 for an 'or' type match. But note that there is no ability to 
// force matching at least one item from each type of constraint.
class SetOfGeometricConstraints {
  private:
    // Whether or not any data is in here.
    int _active;

    resizable_array_p<DistanceConstraint> _distances;
    resizable_array_p<BondAngleConstraint> _bond_angles;
    resizable_array_p<TorsionAngleConstraint> _torsion_angles;
    resizable_array_p<PositionalConstraint> _positive_spatial_constraints;
    resizable_array_p<PositionalConstraint> _negative_spatial_constraints;

    int _number_to_match;

    // Private functions.

    int _number_constraints() const;

  public:
    SetOfGeometricConstraints();

    int BuildFromProto(const GeometricConstraints::SetOfConstraints& proto);
    int BuildProto(GeometricConstraints::SetOfConstraints& proto) const;

    void DebugPrint(std::ostream& output) const;

    // Return a unique list of all the atom numbers mentioned in any component.
    resizable_array<int> AtomNumbersPresent() const;

    int IsValid() const;

    int Active() const { return _active ;}

    // The only reason this is non-const is that the distance constraint
    // may choose to rotate the molecule to make computing the tube exclusion
    // zone easier.
    int Matches(Molecule& m) const;
    int Matches(Molecule& m, const Set_of_Atoms& embedding) const;
};

}  // namespace geometric_constraints
#endif // MOLECULE_LIB_GEOMETRIC_CONSTRAINTS_H
