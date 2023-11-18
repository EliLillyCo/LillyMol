#include <memory>

#include "geometric_constraints.h"

namespace geometric_constraints {

using std::cerr;

// Initialize with arbitrary values that will prevent matching anything.
AllowedRange::AllowedRange() {
  _min_value = 1.0;
  _max_value = 0.0;
}

std::ostream&
operator<< (std::ostream& output, const AllowedRange& arange) {
  output << "AllowedRange [" << arange._min_value << ',' << arange._max_value << ']';
  return output;
}

int
AllowedRange::BuildFromProto(const GeometricConstraints::Range& proto) {
  _min_value = proto.min();
  _max_value = proto.max();

  // If an upper bound is not specified, assume infinite.
  if (_min_value > 0.0 && _max_value == 0.0) {
    _max_value = std::numeric_limits<float>::max();
  }

  if (_min_value > _max_value) {
    cerr << "AllowedRange::BuildFromProto:invalid range " << *this << '\n';
    return 0;
  }

  return 1;
}

int
AllowedRange::BuildProto(GeometricConstraints::Range& proto) const {
  proto.set_min(_min_value);
  proto.set_max(_max_value);

  return 1;
}


int
ConstraintBaseClass::IsValid() const {
  // Unset is a valid state.
  if (!_allowed_range.Active() && _atoms.empty()) {
    return 1;
  }

  const int natoms = _atoms.number_elements();

  if (natoms < 2)  // The DistanceConstraint has 2 atoms.
    return 0;
  if (natoms > 4)  // Torsion has 4 atoms.
    return 0;

  // We have atoms specified, but no _allowed_range.
  if (!_allowed_range.Active()) {
    return 0;
  }

  // Are all the atom numbers unique?
  for (int i = 0; i < natoms; ++i) {
    int a1 = _atoms[i];
    for (int j = i + 1; j < natoms; ++j) {
      if (_atoms[j] == a1)
        return 0;
    }
  }

  return 1;
}

int
ConstraintBaseClass::SetAtoms(int a1, int a2) {
  return (_atoms.add_if_not_already_present(a1)> 0) +
         (_atoms.add_if_not_already_present(a2)> 0);
}

int
ConstraintBaseClass::SetAtoms(int a1, int a2, int a3) {
  return (_atoms.add_if_not_already_present(a1)> 0) +
         (_atoms.add_if_not_already_present(a2)> 0) +
         (_atoms.add_if_not_already_present(a3)> 0);
}

int
ConstraintBaseClass::SetAtoms(int a1, int a2, int a3, int a4) {
  return (_atoms.add_if_not_already_present(a1)> 0) +
         (_atoms.add_if_not_already_present(a2)> 0) +
         (_atoms.add_if_not_already_present(a3)> 0) +
         (_atoms.add_if_not_already_present(a4)> 0);
}

int
ConstraintBaseClass::AtomNumbersPresent(resizable_array<int>& atom_numbers) const {
  int rc = 0;
  for (int a : _atoms) {
    rc += atom_numbers.add_if_not_already_present(a);
  }
  return rc;
}

std::ostream& operator<<(std::ostream& output, const ConstraintBaseClass& constraint) {
  for (int i : constraint._atoms) {
    output << ' ' << i;
  }
  output << ' ' << constraint._allowed_range;

  return output;
}

DistanceConstraint::DistanceConstraint() {
  _tube_radius = 0.0f;
}

int
DistanceConstraint::BuildFromProto(const GeometricConstraints::Distance& proto) {
  if (! _allowed_range.BuildFromProto(proto.range())) {
    cerr << "DistanceConstraint::BuildFromProto:invalid range " << proto.ShortDebugString() << '\n';
    return 0;
  }
  if (SetAtoms(proto.a1(), proto.a2()) != 2) {
    cerr << "DistanceConstraint::BuildFromProto:invalid atom numbers " << proto.ShortDebugString() << '\n';
    return 0;
  }

  _tube_radius = proto.tube_radius();

  if (! IsValid()) {
    cerr << "DistanceConstraint::BuildFromProto:invalid " << *this << '\n';
    return 0;
  }

  return 1;
}

std::ostream& operator<< (std::ostream& output, const DistanceConstraint& constraint) {
  const ConstraintBaseClass& me = constraint;
  output << "DistanceConstraint:atoms " << me;
  return output;
}

int
DistanceConstraint::Matches(Molecule& m) {
  if (! _allowed_range.Matches(m.distance_between_atoms(_atoms[0], _atoms[1]))) {
    return 0;
  }

  if (_tube_radius <= 0.0f) {
    return 1;
  }

  // First save the coordinates.
  std::unique_ptr<float[]> coords = m.GetCoordinates();
  const Coordinates c1 = m.get_coords(_atoms[0]);
  m.translate_atoms(-c1.x(), -c1.y(), -c1.z());

  return 0;
}

int
DistanceConstraint::Matches(Molecule& m, const Set_of_Atoms& embedding) {
  if (!_allowed_range.Matches(m.distance_between_atoms(embedding[_atoms[0]], embedding[_atoms[1]]))) {
    return 0;
  }

  if (_tube_radius <= 0.0f) {
    return 1;
  }

  // First save the coordinates.
  std::unique_ptr<float[]> coords = m.GetCoordinates();
  int rc = MatchesWithEmptyTube(m, embedding[_atoms[0]], embedding[_atoms[1]]);
  m.SetCoordinates(coords.get());
  return rc;
}

int
DistanceConstraint::MatchesWithEmptyTube(Molecule& m, atom_number_t a1, atom_number_t a2) const {
  // Shift `m` so that `a1` is at the origin.
  const Coordinates c1 = m.get_coords(a1);
  m.translate_atoms(-c1.x(), -c1.y(), -c1.z());

  // And rotate so `a2` is along the X axis.
  Coordinates x(static_cast<coord_t>(1.0), static_cast<coord_t>(0.0), static_cast<coord_t>(0.0));
  const Atom* atom2 = m.atomi(a2);

  // concepts lifted from tshadow...
  Coordinates r(atom2->x(), atom2->y(), atom2->z());

  r.normalise();
  angle_t theta = x.angle_between_unit_vectors(r);

  r.cross_product(x);
  r.normalise();

  m.rotate_atoms(r, theta);

  const float tube_rhs = m.x(a2);
  // Take the square to avoid calculating square roots.
  const float r2 = _tube_radius * _tube_radius;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (i == a1 || i == a2) {
      continue;
    }
    const Atom* a = m.atomi(i);
    if (a->x() < 0.0f || a->x() > tube_rhs) {
      continue;
    }
    const float r = a->y() * a->y() + a->z() * a->z();
    if (r < r2) {
      return 0;
    }
  }

  return 1;
}

int
BondAngleConstraint::BuildFromProto(const GeometricConstraints::BondAngle& proto) {
  if (! _allowed_range.BuildFromProto(proto.range())) {
    cerr << "DistanceConstraint::BuildFromProto:invalid range " << proto.ShortDebugString() << '\n';
    return 0;
  }
  _allowed_range *= DEG2RAD;

  if (SetAtoms(proto.a1(), proto.a2(), proto.a3()) != 3) {
    cerr << "BondAngle::BuildFromProto:invalid atom numbers " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! IsValid()) {
    cerr << "BondAngleConstraint::BuildFromProto:invalid " << *this << '\n';
    return 0;
  }

  return 1;
}

std::ostream&
operator<< (std::ostream& output, const BondAngleConstraint& constraint) {
  const ConstraintBaseClass& me = constraint;
  output << "BondAngleConstraint:atoms " <<  me;
  return output;
}

int
BondAngleConstraint::Matches(const Molecule& m) {
//std::cerr << "Angle is " << m.bond_angle(_a1, _a2, _a3) << " radians " << _allowed_range << '\n';
  return _allowed_range.Matches(m.bond_angle(_atoms[0], _atoms[1], _atoms[2]));
}

int
BondAngleConstraint::Matches(const Molecule& m, const Set_of_Atoms& embedding) {
  return _allowed_range.Matches(m.bond_angle(embedding[_atoms[0]], embedding[_atoms[1]], embedding[_atoms[2]]));
}

int
TorsionAngleConstraint::BuildFromProto(const GeometricConstraints::TorsionAngle& proto) {
  if (! _allowed_range.BuildFromProto(proto.range())) {
    cerr << "DistanceConstraint::BuildFromProto:invalid range " << proto.ShortDebugString() << '\n';
    return 0;
  }

  _allowed_range *= DEG2RAD;

  if (SetAtoms(proto.a1(), proto.a2(), proto.a3(), proto.a4()) != 4) {
    cerr << "TorsionAngle::BuildFromProto:invalid atom numbers " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! IsValid()) {
    cerr << "TorsionAngle::BuildFromProto:invalid " << *this << '\n';
    return 0;
  }

  return 1;
}

std::ostream&
operator<< (std::ostream& output, const TorsionAngleConstraint& constraint) {
  const ConstraintBaseClass& me = constraint;
  output << "TorsionAngleConstraint:atoms " << me;
  return output;
}

int
TorsionAngleConstraint::Matches(const Molecule& m) {
  return _allowed_range.Matches(m.dihedral_angle(_atoms[0], _atoms[1], _atoms[2], _atoms[3]));
}

int
TorsionAngleConstraint::Matches(const Molecule& m, const Set_of_Atoms& embedding) {
  return _allowed_range.Matches(m.dihedral_angle(embedding[_atoms[0]], embedding[_atoms[1]], embedding[_atoms[2]], embedding[_atoms[3]]));
}

SetOfGeometricConstraints::SetOfGeometricConstraints() {
  _active = 0;
  _number_to_match = 0;
}

// The number of constraints across all different kinds.
int
SetOfGeometricConstraints::_number_constraints() const {
  return _distances.number_elements() +
         _bond_angles.number_elements() +
         _torsion_angles.number_elements() +
         _positive_spatial_constraints.number_elements() +
         _negative_spatial_constraints.number_elements() +
         _angle_between_vectors.number_elements();
}

int
SetOfGeometricConstraints::IsValid() const {
  int nset = _number_constraints();
  if (_number_to_match > nset) {
    return 0;
  }

  if (nset == 0 &&_number_to_match > 0) {
    return 0;
  }

  return 1;
}

int
SetOfGeometricConstraints::BuildFromProto(const GeometricConstraints::SetOfConstraints& proto) {
  if (proto.number_to_match() > 0) {
    _number_to_match = proto.number_to_match();
  }
  for (const auto& dist : proto.distance()) {
    std::unique_ptr<DistanceConstraint> c = std::make_unique<DistanceConstraint>();
    if (! c->BuildFromProto(dist)) {
      cerr << "SetOfConstraints::BuildFromProto:bad distance " << dist.ShortDebugString() << '\n';
      return 0;
    }
    _distances << c.release();
  }

  for (const auto& angle : proto.bond_angle()) {
    std::unique_ptr<BondAngleConstraint> c = std::make_unique<BondAngleConstraint>();
    if (! c->BuildFromProto(angle)) {
      cerr << "SetOfConstraints::BuildFromProto:bad angle " << angle.ShortDebugString() << '\n';
      return 0;
    }
    _bond_angles << c.release();
  }

  for (const auto& torsion : proto.torsion_angle()) {
    std::unique_ptr<TorsionAngleConstraint> c = std::make_unique<TorsionAngleConstraint>();
    if (! c->BuildFromProto(torsion)) {
      cerr << "SetOfConstraints::BuildFromProto:bad torsion " << torsion.ShortDebugString() << '\n';
      return 0;
    }

    _torsion_angles << c.release();
  }

  for (const auto& spatial: proto.positive_spatial_constraint()) {
    std::unique_ptr<PositionalConstraint> p = std::make_unique<PositionalConstraint>();
    if (! p->BuildFromProto(spatial)) {
      cerr << "SetOfConstraints::BuildFromProto:bad spatial " << spatial.ShortDebugString() << '\n';
      return 0;
    }

    _positive_spatial_constraints << p.release();
  }

  for (const auto& spatial: proto.negative_spatial_constraint()) {
    std::unique_ptr<PositionalConstraint> p = std::make_unique<PositionalConstraint>();
    if (! p->BuildFromProto(spatial)) {
      cerr << "SetOfConstraints::BuildFromProto:bad spatial " << spatial.ShortDebugString() << '\n';
      return 0;
    }

    _negative_spatial_constraints << p.release();
  }

  for (const auto& vector_angle : proto.angle_between_vectors()) {
    std::unique_ptr<AngleBetweenVectors> s = std::make_unique<AngleBetweenVectors>();
    if (! s->BuildFromProto(vector_angle)) {
      cerr << "SetOfConstraints::BuildFromProto:bad vector angle " << vector_angle.ShortDebugString() << '\n';
      return 0;
    }

    _angle_between_vectors << s.release();
  }

  _active = 1;
  if (_number_to_match == 0) {
    _number_to_match = _number_constraints();
  }

  return IsValid();
}


PositionalConstraint::PositionalConstraint() {
}

template <typename P>
int
PositionalConstraint::BuildFromProto(const P& proto) {
  if (! _allowed_range.BuildFromProto(proto.range())) {
    cerr << "PositionalConstraint::BuildFromProto:invalid range " << proto.ShortDebugString() << '\n';
    return 0;
  }

  for (auto atom : proto.atom()) {
    _atoms.add_if_not_already_present(atom);
  }

  if (_atoms.empty()) {
    cerr << "PositionalConstraint::BuildFromProto:no atoms\n";
    return 0;
  }

  _match_type = proto.match_type();

  _coords.set_x(proto.position().x());
  _coords.set_y(proto.position().y());
  _coords.set_z(proto.position().z());

  return 1;
}

int
PositionalConstraint::BuildFromProto(const GeometricConstraints::PositiveSpatialConstraint& proto) {
  return BuildFromProto<GeometricConstraints::PositiveSpatialConstraint>(proto);
}
int
PositionalConstraint::BuildFromProto(const GeometricConstraints::NegativeSpatialConstraint& proto) {
  return BuildFromProto<GeometricConstraints::NegativeSpatialConstraint>(proto);
}

int
PositionalConstraint::Matches(const Molecule& m,
                              const Set_of_Atoms& embedding) const {
  switch (_match_type) {
    case GeometricConstraints::AVERAGE:
      return MatchesAverage(m, embedding);
    case GeometricConstraints::ANY:
      return MatchesAny(m, embedding);
    case GeometricConstraints::ALL:
      return MatchesAll(m, embedding);
    default:
      cerr << "PositionalConstraint::Matches:no match type\n";
      return 0;
  }
}

int
PositionalConstraint::MatchesAverage(const Molecule& m,
                              const Set_of_Atoms& embedding) const {
  Space_Vector<double> acc;
  for (int a : embedding) {
    const Atom& atom = m.atom(a);
    const Space_Vector<double> tmp(atom.x(), atom.y(), atom.z());
    acc += tmp;
  }

  acc /= static_cast<double>(embedding.size());

  double dist = _coords.distance(acc);

  return _allowed_range.Matches(dist);
}

int
PositionalConstraint::MatchesAny(const Molecule& m,
                              const Set_of_Atoms& embedding) const {
  for (int a : embedding) {
    const Atom& atom = m.atom(a);
    const Space_Vector<double> tmp(atom.x(), atom.y(), atom.z());
    double d = _coords.distance(tmp);
    if (_allowed_range.Matches(d)) {
      return 1;
    }
  }

  return 0;
}

int
PositionalConstraint::MatchesAll(const Molecule& m,
                              const Set_of_Atoms& embedding) const {
  for (int a : embedding) {
    const Atom& atom = m.atom(a);
    const Space_Vector<double> tmp(atom.x(), atom.y(), atom.z());
    double d = _coords.distance(tmp);
    if (! _allowed_range.Matches(d)) {
      return 0;
    }
  }

  return 1;
}

AngleBetweenVectors::AngleBetweenVectors() {
  _base1 = -1;
  _end1 = -1;

  _base2 = -1;
  _end2 = -1;
}

int
AngleBetweenVectors::BuildFromProto(const GeometricConstraints::AngleBetweenVectors& proto) {
  if (! proto.has_v1() || ! proto.has_v2() || ! proto.has_angle()) {
    cerr << "AngleBetweenVectors::BuildFromProto:must specify all components\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  _base1 = proto.v1().base();
  _end1 = proto.v1().end();

  _base2 = proto.v2().base();
  _end2 = proto.v2().end();

  if (_base1 == _end1 || _base2 == _end2) {
    cerr << "AngleBetweenVectors::BuildFromProto:invalid matched atoms\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! _allowed_range.BuildFromProto(proto.angle())) {
    cerr << "AngleBetweenVectors::BuildFromProto:invalid angle\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  // the human specifies things in degrees.
  _allowed_range *= DEG2RAD;
  // cerr << "after building " << _allowed_range << '\n';

  return 1;
}

int
AngleBetweenVectors::BuildProto(GeometricConstraints::AngleBetweenVectors& proto) const {
  proto.mutable_v1()->set_base(_base1);
  proto.mutable_v1()->set_end(_end1);

  proto.mutable_v2()->set_base(_base2);
  proto.mutable_v2()->set_end(_end2);

  return _allowed_range.BuildProto(*proto.mutable_angle());  // does not convert back to degrees...TODO:ianwatson
}

int
AngleBetweenVectors::Matches(const Molecule& m, const Set_of_Atoms& embedding) const {
  atom_number_t base1 = embedding[_base1];
  atom_number_t end1 = embedding[_end1];
  Space_Vector<float> v1 = m[end1] - m[base1];

  atom_number_t base2 = embedding[_base2];
  atom_number_t end2 = embedding[_end2];
  Space_Vector<float> v2 = m[end2] - m[base2];

  // cerr << "angle is " << (v1.angle_between(v2) * RAD2DEG) << " matches " << _allowed_range.Matches(v1.angle_between(v2)) << '\n';
  return _allowed_range.Matches(v1.angle_between(v2));
}

int
SetOfGeometricConstraints::Matches(Molecule& m) const {
  int matched_here = 0;
  for (auto c : _distances) {
    if (c->Matches(m)) {
      matched_here++;
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    } else {
      return 0;
    }
  }
  for (auto c : _bond_angles) {
    if (c->Matches(m)) {
      matched_here++;
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    } else {
      return 0;
    }
  }
  for (auto c : _torsion_angles) {
    if (c->Matches(m)) {
      matched_here++;
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    } else {
      return 0;
    }
  }

  return 0;
}

int
SetOfGeometricConstraints::Matches(Molecule& m, const Set_of_Atoms& embedding) const {
  // cerr << "SetOfGeometricConstraints::Matches " << m.smiles() << " need " << _number_to_match << '\n';
  int matched_here = 0;
#ifdef DEBUG_CONSTRAINTS_MATCHES
  std::cerr << "Checking " << _distances.size() << " distance_constraints\n";
#endif
  for (auto c : _distances) {
    if (c->Matches(m, embedding)) {
      matched_here++;
#ifdef DEBUG_CONSTRAINTS_MATCHES
      cerr << "distance constaint matched, now " << matched_here << " matches\n";
#endif
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    }
#ifdef DEBUG_CONSTRAINTS_MATCHES
    else {
      cerr << "Constraint " << *c << " did not match\n";
      const auto atoms = c->Atoms();
      std::cerr << "dist " << m.distance_between_atoms(embedding[atoms[0]], embedding[atoms[1]]) << '\n';
    }
#endif

  }
  for (auto c : _bond_angles) {
    if (c->Matches(m, embedding)) {
      matched_here++;
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    }
  }
  for (auto c : _torsion_angles) {
    if (c->Matches(m, embedding)) {
      matched_here++;
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    }
  }
  for (auto &c : _positive_spatial_constraints) {
    if (c->Matches(m, embedding)) {
      ++matched_here;
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    }
  }
  for (auto &c : _negative_spatial_constraints) {
    if (c->Matches(m, embedding)) {
      return 0;
    }
    ++matched_here;
  }

  for (auto& v : _angle_between_vectors) {
    if (!v->Matches(m, embedding)) {
      return 0;
    }
    ++matched_here;
  }

  if (matched_here >= _number_to_match) {
    return matched_here;
  }

#ifdef DEBUG_CONSTRAINTS_MATCHES
  std::cerr << "only matched " << matched_here << " items\n";
#endif

  return 0;
}

template <typename T>
void
write_constraints(const resizable_array_p<T>& constraints, std::ostream& output) {
  for (const T * constraint : constraints) {
    output << ' ' << *constraint;
  }
  output << '\n';
}
void
SetOfGeometricConstraints::DebugPrint(std::ostream& output) const {
  output << "SetOfGeometricConstraints\n";
  if (_distances.size()) {
    output << " distances";
    write_constraints(_distances, output);
  }
  if (_bond_angles.size()) {
    output << " bond angles";
    write_constraints(_bond_angles, output);
  }
  if (_torsion_angles.size()) {
    output << " torsion angles";
    write_constraints(_torsion_angles, output);
  }
  output << "must match " << _number_to_match << '\n';
}

resizable_array<int>
SetOfGeometricConstraints::AtomNumbersPresent() const {
  resizable_array<int> result;
  for (const auto * c : _distances) {
    c->AtomNumbersPresent(result);
  }
  for (const auto * c : _bond_angles) {
    c->AtomNumbersPresent(result);
  }
  for (const auto * c : _torsion_angles) {
    c->AtomNumbersPresent(result);
  }

  return result;
}

int
SetOfGeometricConstraints::BuildProto(GeometricConstraints::SetOfConstraints&) const {
  // TODO:implement this
  cerr << "SetOfConstraints::BuildProto:implement\n";
  return 1;
}

}  // namespace geometric_constraints
