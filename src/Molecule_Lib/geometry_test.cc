// Tests for various geometry related functions.

#include <random>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "aromatic.h"
#include "substructure.h"

namespace {

struct GeometryTestInput {
  IWString smiles;
  atom_number_t atom;
  float dist;
  Coordinates coords;
};

class GeometryTestBed: public testing::TestWithParam<GeometryTestInput>  {
  protected:
    Molecule _m;
    Coordinates _result;
    Coordinates _start_atom;
};

TEST_P(GeometryTestBed, AlongX1a) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles)) << "invalid smiles " << params.smiles;
  ASSERT_TRUE(_m.LocationOfSubstituent(params.atom, params.dist, _result));
  if (params.coords.x() != 0.0 && params.coords.y() != 0.0 && params.coords.z() != 0.0) {
    EXPECT_FLOAT_EQ(params.coords.x(), _result.x()) << " X mismatch, expect " << params.coords.x() << " got " << _result.x();
    EXPECT_FLOAT_EQ(params.coords.y(), _result.y()) << " Y mismatch, expect " << params.coords.y() << " got " << _result.y();
    EXPECT_FLOAT_EQ(params.coords.z(), _result.z()) << " Z mismatch, expect " << params.coords.z() << " got " << _result.z();
  }

  // `result` must be at the right distance from params.zatom
  _m.get_coordinates(params.atom, _start_atom);
  float d = _start_atom.distance(_result);
  EXPECT_NEAR(d, params.dist, 0.001);

  // All bonded atoms must be further away than params.dist

  Space_Vector<float> v12 = _result - _m[params.atom];

  const Atom& a = _m[params.atom];
  for (const Bond* b : a) {
    atom_number_t o = b->other(params.atom);
    const Atom& ao = _m[o];

    // A not very stringent test.
    EXPECT_GT(_result.distance(ao), params.dist);

    // Make sure we are increasing the distance to the other atoms by a reasonable amount.

    /*  result
    //  |
    //  |  v12
    //  |
    //  params.atom
    //   \
    //     \   v23
    //      \
    //        o                 */
    const Space_Vector<float> v23 = ao - _m[params.atom];
    // In realistic situtions it might be the tetrahedral bond 108 degrees.
    EXPECT_GT(v12.angle_between(v23), DEG2RAD * 100.0) << " bad angle " <<
                (v12.angle_between(v23) * RAD2DEG);
  }

  // Never a good idea to use random numbers in a unit test, I know...
  // And yes, we know this is not uniformly sampling over rotation space.
  std::random_device rd;
  std::uniform_real_distribution<float> ubox(-1.0, 1.0);
  std::uniform_real_distribution<float> upi(-M_PI, M_PI);

  for (int i = 0; i < 10; ++i) {
    Coordinates rot(ubox(rd), ubox(rd), ubox(rd));
    rot.normalise();
    const angle_t angle = upi(rd);
    _m.translate_atoms(ubox(rd), ubox(rd), ubox(rd));
    _m.rotate_atoms(rot, angle);

    ASSERT_TRUE(_m.LocationOfSubstituent(params.atom, params.dist, _result));

    v12 = _result - _m[params.atom];

    _m.get_coordinates(params.atom, _start_atom);
    float d = _start_atom.distance(_result);
    EXPECT_NEAR(d, params.dist, 0.001) << "distance not right, got " << d << 
        "expected " << params.dist;

    for (const Bond* b : _m[params.atom]) {
      atom_number_t o = b->other(params.atom);
      const Space_Vector<float> v23 = _m[o] - _m[params.atom];
      EXPECT_GT(v12.angle_between(v23), DEG2RAD * 100.0) << " bad angle " <<
                (v12.angle_between(v23) * RAD2DEG);
    }
  }
}
INSTANTIATE_TEST_SUITE_P(GeometryTestBed, GeometryTestBed, testing::Values(
  GeometryTestInput{"C{{0,0,0}}C{{1,0,0}}", 1, 1, {2.0, 0.0, 0.0}},
  GeometryTestInput{"C{{0,0,0}}C{{-1,0,0}}", 1, 1, {-2.0, 0.0, 0.0}},
  GeometryTestInput{"C{{0,0,0}}C{{0,1,0}}", 1, 1, {0.0, 2.0, 0.0}},
  GeometryTestInput{"C{{0,0,0}}C{{0,-1,0}}", 1, 1, {0.0, -2.0, 0.0}},
  GeometryTestInput{"C{{0,0,0}}C{{0,0,1}}", 1, 1, {0.0, 0.0, 2.0}},
  GeometryTestInput{"C{{0,0,0}}C{{0,0,-1}}", 1, 1, {0.0, 0.0, -2.0}},
  GeometryTestInput{"C{{0,0,0}}C{{1,1,1}}", 1, 1, {1.57735, 1.57735, 1.57735}},
  GeometryTestInput{"C{{1,1,1}}C{{2,2,2}}", 1, 1, {2.57735, 2.57735, 2.57735}},
  GeometryTestInput{"C{{-0.0187,1.5258,0.0104}}C{{0.0021,-0.0041,0.002}}C{{0.7182,-0.4975,-1.2568}}", 1, 1, {0,0,0}},
  GeometryTestInput{"C{{-0.0187,1.5258,0.0104}}C{{0.0021,-0.0041,0.002}}(C{{0.7182,-0.4975,-1.2568}})C{{0.7421,-0.5109,1.2415}}", 1, 1, {0,0,0}}
));

}  // namespace
