// Tester for Space_Vector class

#include <memory>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "space_vector.h"

namespace {

using testing::ElementsAre;

template <typename T>
class TestSpaceVector : public testing::Test {
};

using MyTypes = ::testing::Types<float, double>;
TYPED_TEST_SUITE(TestSpaceVector, MyTypes);

TYPED_TEST(TestSpaceVector, TestSetGet) {
  Space_Vector<TypeParam> v;
  v.setxyz(1.0, 2.0, 0.0);
  EXPECT_EQ(v.x(), 1.0);
  EXPECT_EQ(v.y(), 2.0);
  EXPECT_EQ(v.z(), 0.0);
}
TYPED_TEST(TestSpaceVector, TestGetArray) {
  Space_Vector<TypeParam> v;
  v.setxyz(1.0, 2.0, 0.0);
  std::unique_ptr<double[]> values = std::make_unique<double[]>(3);
  v.getxyz(values.get());
  EXPECT_EQ(values.get()[0], 1.0);
  EXPECT_EQ(values.get()[1], 2.0);
  EXPECT_EQ(values.get()[2], 0.0);
}

TYPED_TEST(TestSpaceVector, TestReadOK) {
  const const_IWSubstring input("0,1.1,2.2");
  Space_Vector<TypeParam> v;
  ASSERT_TRUE(v.read(input, ','));
  EXPECT_FLOAT_EQ(v.x(), 0.0);
  EXPECT_FLOAT_EQ(v.y(), 1.1);
  EXPECT_FLOAT_EQ(v.z(), 2.2);
}
TYPED_TEST(TestSpaceVector, TestAdd) {
  Space_Vector<TypeParam> v(1.0, 2.2, 3.3);
  v.add(-1.0, -2.2, 3.3);
  EXPECT_FLOAT_EQ(v.x(), 0.0);
  EXPECT_FLOAT_EQ(v.y(), 0.0);
  EXPECT_FLOAT_EQ(v.z(), 6.6);
}
TYPED_TEST(TestSpaceVector, TestTranslate) {
  Space_Vector<TypeParam> v(1.0, 2.2, 3.3);
  v.translate(1.0, 2.2, 3.3);
  EXPECT_FLOAT_EQ(v.x(), 2.0);
  EXPECT_FLOAT_EQ(v.y(), 4.4);
  EXPECT_FLOAT_EQ(v.z(), 6.6);
}
TYPED_TEST(TestSpaceVector, TestNormalise) {
  TypeParam x = 2.0;
  TypeParam y = 1.0;
  TypeParam z = 3.0;
  Space_Vector<TypeParam> v(x, y, z);
  TypeParam norm = v.norm();
  EXPECT_FLOAT_EQ(norm, sqrt(4 + 1 + 9));
  v.normalise();
  EXPECT_FLOAT_EQ(v.x(), x/norm);
  EXPECT_FLOAT_EQ(v.y(), y/norm);
  EXPECT_FLOAT_EQ(v.z(), z/norm);
}
TYPED_TEST(TestSpaceVector, TestOpEqual) {
  TypeParam x = 2.0;
  TypeParam y = 1.0;
  TypeParam z = 3.0;
  const Space_Vector<TypeParam> v1(x, y, z);
  const Space_Vector<TypeParam> v2 = v1;
  EXPECT_EQ(v1, v2);
}
TYPED_TEST(TestSpaceVector, TestOpPlusScalar) {
  TypeParam x = 2.0;
  TypeParam y = 1.0;
  TypeParam z = 3.0;
  Space_Vector<TypeParam> v(x, y, z);
  TypeParam add = static_cast<TypeParam>(1.5);
  v += add;
  EXPECT_FLOAT_EQ(v.x(), x + add);
  EXPECT_FLOAT_EQ(v.y(), y + add);
  EXPECT_FLOAT_EQ(v.z(), z + add);
}
TYPED_TEST(TestSpaceVector, TestOpPlusVector) {
  TypeParam x = 2.0;
  TypeParam y = -1.0;
  TypeParam z = 3.0;
  Space_Vector<TypeParam> v1(x, y, z);
  Space_Vector<TypeParam> v2 = v1 + 1.0;
  Space_Vector<TypeParam> v3 = v1 + v2;
  EXPECT_FLOAT_EQ(v3.x(), 2 * x + 1);
  EXPECT_FLOAT_EQ(v3.y(), 2 * y + 1);
  EXPECT_FLOAT_EQ(v3.z(), 2 * z + 1);
}
TYPED_TEST(TestSpaceVector, TestOpSubScalar) {
  TypeParam x = 2.0;
  TypeParam y = -1.0;
  TypeParam z = 3.0;
  Space_Vector<TypeParam> v1(x, y, z);
  Space_Vector<TypeParam> v2 = v1 - 1.0;
  EXPECT_FLOAT_EQ(v2.x(), x - 1);
  EXPECT_FLOAT_EQ(v2.y(), y - 1);
  EXPECT_FLOAT_EQ(v2.z(), z - 1);
}
TYPED_TEST(TestSpaceVector, TestOpSubVector) {
  TypeParam x = 2.0;
  TypeParam y = -1.0;
  TypeParam z = 3.0;
  Space_Vector<TypeParam> v1(x, y, z);
  Space_Vector<TypeParam> v2 = v1 + 1.0;
  Space_Vector<TypeParam> v3 = v1 - v2;
  EXPECT_FLOAT_EQ(v3.x(), -1);
  EXPECT_FLOAT_EQ(v3.y(), -1);
  EXPECT_FLOAT_EQ(v3.z(), -1);
}
TYPED_TEST(TestSpaceVector, TestOpMultiplyZero) {
  TypeParam x = 2.0;
  TypeParam y = -1.0;
  TypeParam z = 3.0;
  Space_Vector<TypeParam> v1(x, y, z);
  Space_Vector<TypeParam> v2 = v1 * static_cast<TypeParam>(0);
  EXPECT_FLOAT_EQ(v2.x(), 0);
  EXPECT_FLOAT_EQ(v2.y(), 0);
  EXPECT_FLOAT_EQ(v2.z(), 0);
}
TYPED_TEST(TestSpaceVector, TestOpMultiplyNonZero) {
  TypeParam x = 2.0;
  TypeParam y = -1.0;
  TypeParam z = 3.0;
  Space_Vector<TypeParam> v1(x, y, z);
  TypeParam mult = -4.2;
  Space_Vector<TypeParam> v2 = v1 / static_cast<TypeParam>(mult);
  EXPECT_FLOAT_EQ(v2.x(), x / mult);
  EXPECT_FLOAT_EQ(v2.y(), y / mult);
  EXPECT_FLOAT_EQ(v2.z(), z / mult);
}
TYPED_TEST(TestSpaceVector, TestOpMinus) {
  TypeParam x = 2.0;
  TypeParam y = -1.0;
  TypeParam z = 3.0;
  Space_Vector<TypeParam> v1(x, y, z);
  Space_Vector<TypeParam> v2 = - v1;
  EXPECT_FLOAT_EQ(v2.x(), -x);
  EXPECT_FLOAT_EQ(v2.y(), -y);
  EXPECT_FLOAT_EQ(v2.z(), -z);
}
TYPED_TEST(TestSpaceVector, TestOpMultiplyEqualsSame) {
  TypeParam x = 2.0;
  TypeParam y = -1.0;
  TypeParam z = 3.0;
  const Space_Vector<TypeParam> v1(x, y, z);
  const Space_Vector<TypeParam> v2(v1);
  Space_Vector<TypeParam> v3 = v1 * v2;
  EXPECT_FLOAT_EQ(v3.x(), 0);
  EXPECT_FLOAT_EQ(v3.y(), 0);
  EXPECT_FLOAT_EQ(v3.z(), 0);
}
TYPED_TEST(TestSpaceVector, TestOpMultiplyEqualsOpposite) {
  TypeParam x = 2.0;
  TypeParam y = -1.0;
  TypeParam z = 3.0;
  Space_Vector<TypeParam> v1(x, y, z);
  Space_Vector<TypeParam> v2 = -v1;
  Space_Vector<TypeParam> v3 = v1 * v2;
  EXPECT_FLOAT_EQ(v3.x(), 0);
  EXPECT_FLOAT_EQ(v3.y(), 0);
  EXPECT_FLOAT_EQ(v3.z(), 0);
  v2.normalise();
  v3 = v1 * v2;
  EXPECT_NEAR(v3.norm(), 0, 1.0e-07);
}
TYPED_TEST(TestSpaceVector, TestOpMultiplyccw) {
  Space_Vector<TypeParam> v1(1, 0, 0);
  Space_Vector<TypeParam> v2(0, 1, 0);
  Space_Vector<TypeParam> v3 = v1 * v2;
  Space_Vector<TypeParam> expected(0, 0, 1);
  EXPECT_EQ(v3, expected);
  EXPECT_FLOAT_EQ(v1.angle_between(v2), 0.5 * M_PI);
  EXPECT_FLOAT_EQ(v2.angle_between(v1), 0.5 * M_PI);
  EXPECT_FLOAT_EQ(expected.angle_between(v1), 0.5 * M_PI);
  EXPECT_FLOAT_EQ(expected.angle_between(v2), 0.5 * M_PI);
}
TYPED_TEST(TestSpaceVector, TestOpMultiplymcw) {
  Space_Vector<TypeParam> v1(1, 0, 0);
  Space_Vector<TypeParam> v2(0, 1, 0);
  Space_Vector<TypeParam> v3 = v2 * v1;
  Space_Vector<TypeParam> expected(0, 0, -1);
  EXPECT_EQ(v3, expected);
}
TYPED_TEST(TestSpaceVector, TestCrossProductccw) {
  Space_Vector<TypeParam> v1(1, 0, 0);
  Space_Vector<TypeParam> v2(0, 1, 0);
  v1.cross_product(v2);
  Space_Vector<TypeParam> expected(0, 0, 1);
  EXPECT_EQ(v1, expected);
}
TYPED_TEST(TestSpaceVector, TestCrossProductcw) {
  Space_Vector<TypeParam> v1(1, 0, 0);
  Space_Vector<TypeParam> v2(0, -1, 0);
  v1.cross_product(v2);
  Space_Vector<TypeParam> expected(0, 0, -1);
  EXPECT_EQ(v1, expected);
}
TYPED_TEST(TestSpaceVector, TestToResizableArray) {
  const TypeParam x = 0;
  const TypeParam y = 1;
  const TypeParam z = 2;
  Space_Vector<TypeParam> v(x, y, z);
  EXPECT_THAT(v.ToResizableArray(), ElementsAre(x, y, z));
}
TYPED_TEST(TestSpaceVector, TestFormUnitVector) {
  Space_Vector<TypeParam> v1(0.6, -1.7, 2.32);
  Space_Vector<TypeParam> v2(1.6, 1.2, -3.2);
  const Space_Vector<TypeParam> norm = v1.form_unit_vector(v2);
  EXPECT_FLOAT_EQ(norm.norm(), 1);
}
TYPED_TEST(TestSpaceVector, TestCloserThan) {
  const TypeParam x = 0;
  const TypeParam y = 1;
  const TypeParam z = 2;
  Space_Vector<TypeParam> v1(x, y, z);
  EXPECT_TRUE(v1.closer_than(v1, 0.01));
  Space_Vector<TypeParam> v2(v1);
  v2.z() += 1.0;
  EXPECT_TRUE(v1.closer_than(v2, 1.01));
  EXPECT_FALSE(v1.closer_than(v2, 0.99));
}
TYPED_TEST(TestSpaceVector, TestAngleBetween45) {
  Space_Vector<TypeParam> v1(1, 0, 0);
  EXPECT_FLOAT_EQ(v1.angle_between(v1), 0.0);
  Space_Vector<TypeParam> v2(1, 1, 0);
  EXPECT_FLOAT_EQ(v1.angle_between(v2), M_PI * 0.25);
  v2.setxyz(0, 1, 0);
  EXPECT_FLOAT_EQ(v1.angle_between(v2), M_PI * 0.5);
  v2.setxyz(-1, 1, 0);
  EXPECT_FLOAT_EQ(v1.angle_between(v2), M_PI * 0.75);
  v2.setxyz(-1, 0, 0);
  EXPECT_FLOAT_EQ(v1.angle_between(v2), M_PI);
  v2.setxyz(-1, -1, 0);
  EXPECT_FLOAT_EQ(v1.angle_between(v2), M_PI * 0.75);
}


}  // namespace
