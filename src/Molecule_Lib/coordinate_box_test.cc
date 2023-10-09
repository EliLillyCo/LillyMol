// Test the coordinate box idea.

#include <random>

#include "coordinate_box.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {

using coordinate_box::CoordinateBox;
using testing::ElementsAre;
using testing::FloatEq;
using testing::FloatNear;

TEST(TestCoordinateBox, BadBoxSpec1) {
  const const_IWSubstring s;
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecTooManyTokens) {
  const const_IWSubstring s = "1,2,3,4";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecBadCellSize1) {
  const const_IWSubstring s = "0.0,12.0,3.0";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecBadCellSize2) {
  const const_IWSubstring s = "-0.2,12.0,3.0";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecBadX) {
  const const_IWSubstring s = "0.1,0.0,3.0";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecBadY) {
  const const_IWSubstring s = "0.1,1.0,0.0";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecNotAMultipleX) {
  const const_IWSubstring s = "0.1,1.05,1.0";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecNotAMultipleY) {
  const const_IWSubstring s = "0.1,1.00,0.93";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, OriginIsZero) {
  const const_IWSubstring s = "0.1,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  EXPECT_EQ(box.CellNumber(0.0, 0.0, 0.0), 0);
}
#ifdef CELL_NUMBER_CHECK_IN_BOUNDS
TEST(TestCoordinateBox, OutOfRangeXneg) {
  const const_IWSubstring s = "0.1,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  EXPECT_LT(box.CellNumber(-0.1, 0.0, 0.0), 0);
}
TEST(TestCoordinateBox, OutOfRangeX) {
  const const_IWSubstring s = "0.1,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  EXPECT_LT(box.CellNumber(1.1, 0.0, 0.0), 0);
}
TEST(TestCoordinateBox, OutOfRangeYneg) {
  const const_IWSubstring s = "0.1,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  EXPECT_LT(box.CellNumber(0.1, -0.1, 0.0), 0);
}
TEST(TestCoordinateBox, OutOfRangeY) {
  const const_IWSubstring s = "0.1,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  EXPECT_LT(box.CellNumber(0.5, 1.01, 0.0), 0);
}
TEST(TestCoordinateBox, OutOfRangeZneg) {
  const const_IWSubstring s = "0.1,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  EXPECT_LT(box.CellNumber(0.1, 0.7, -1.001), 0);
}
#endif
TEST(TestCoordinateBox, AllCloseX) {
  const const_IWSubstring s = "0.001,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  for (int i = 0; i < 1000; ++i) {
    const float x = i * 0.001;
    const float y = 0.0;
    const float z = 0.0;
    Space_Vector<float> initial(x, y, z);
    const int cell_number = box.CellNumber(x, y, z);
    ASSERT_GE(cell_number, 0);
    const Space_Vector<float> back = box.CoordinatesAsVector<float>(cell_number);
    EXPECT_LT(back.distance(initial), 0.001);
  }
}
TEST(TestCoordinateBox, AllCloseY) {
  const const_IWSubstring s = "0.001,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  for (int i = 0; i < 1000; ++i) {
    const float x = 0.0;
    const float y = i * 0.001;
    const float z = 0.0;
    Space_Vector<float> initial(x, y, z);
    const int cell_number = box.CellNumber(x, y, z);
    ASSERT_GE(cell_number, 0);
    const Space_Vector<float> back = box.CoordinatesAsVector<float>(cell_number);
    EXPECT_LT(back.distance(initial), 0.001);
  }
}
TEST(TestCoordinateBox, AllCloseZ) {
  const const_IWSubstring s = "0.001,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  for (int i = 0; i < 1000; ++i) {
    const float x = 0.0;
    const float y = 0.0;
    const float z = i * 0.001;
    Space_Vector<float> initial(x, y, z);
    const int cell_number = box.CellNumber(x, y, z);
    ASSERT_GE(cell_number, 0);
    const Space_Vector<float> back = box.CoordinatesAsVector<float>(cell_number);
    EXPECT_LT(back.distance(initial), 0.001);
  }
}
TEST(TestCoordinateBox, ArbitraryPositions) {
  const const_IWSubstring s = "0.001,1.00,2.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  for (int i = 0; i < 37; ++i) {
    const float x = static_cast<double>(i) / 37.0;
    for (int j = 0; j < 57; ++j) {
      const float y = static_cast<double>(j) / 57.0;
      for (int k = 0; k < 87; ++k) {
        const float z = static_cast<double>(k) / 87.0;
        Space_Vector<float> initial(x, y, z);
        const int cell_number = box.CellNumber(x, y, z);
        ASSERT_GE(cell_number, 0);
        const Space_Vector<float> back = box.CoordinatesAsVector<float>(cell_number);
        EXPECT_LT(back.distance(initial), 0.001);
      }
    }
  }
}

TEST(TestCellToLayer, TestCellToLayer) {
  using coordinate_box::CellToLayer;

  EXPECT_EQ(CellToLayer(0), 0);
  for (int i = 1; i < 27; ++i) {
    EXPECT_EQ(CellToLayer(i), 1);
  }

  for (int i = 27; i < 124; ++i) {
    EXPECT_EQ(CellToLayer(i), 2);
  }

  for (int i = 125; i < 342; ++i) {
    EXPECT_EQ(CellToLayer(i), 3);
  }

  for (int i = 343; i < 728; ++i) {
    EXPECT_EQ(CellToLayer(i), 4);
  }

  for (int layer = 5; layer < 100; ++layer) {
    const int n = 2 * layer + 1;
    const int previous_layer = layer - 1;
    const int previous_n = 2 * previous_layer + 1;
    int cells_in_box = n * n * n;
    int cells_in_previous_box = previous_n * previous_n * previous_n;
    for (int i = cells_in_previous_box; i < cells_in_box; ++i) {
      EXPECT_EQ(CellToLayer(i), layer);
    }
  }
}

TEST(TestConcentricBox, TestCellsLayer1) {
  constexpr double dx = 1.0;
  coordinate_box::ConcentricBox box(dx);

  Space_Vector<float> coords(0.0, 0.0, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 0);

  coords.setxyz(-0.51, -0.51, -0.51);
  EXPECT_EQ(box.CellNumber(coords), 1);

  coords.setxyz(-0.51, -0.51, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 2);

  coords.setxyz(-0.51, -0.51, 0.51);
  EXPECT_EQ(box.CellNumber(coords), 3);

  coords.setxyz(-0.51, 0.0, -0.51);
  EXPECT_EQ(box.CellNumber(coords), 4);

  coords.setxyz(-0.51, 0.0, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 5);

  coords.setxyz(-0.51, 0.0, 0.51);
  EXPECT_EQ(box.CellNumber(coords), 6);

  coords.setxyz(-0.51, 0.51, -0.51);
  EXPECT_EQ(box.CellNumber(coords), 7);

  coords.setxyz(-0.51, 0.51, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 8);

  coords.setxyz(-0.51, 0.51, 0.51);
  EXPECT_EQ(box.CellNumber(coords), 9);

  // RHS
  coords.setxyz(0.51, -0.51, -0.51);
  EXPECT_EQ(box.CellNumber(coords), 10);

  coords.setxyz(0.51, -0.51, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 11);

  coords.setxyz(0.51, -0.51, 0.51);
  EXPECT_EQ(box.CellNumber(coords), 12);

  coords.setxyz(0.51, 0.0, -0.51);
  EXPECT_EQ(box.CellNumber(coords), 13);

  coords.setxyz(0.51, 0.0, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 14);

  coords.setxyz(0.51, 0.0, 0.51);
  EXPECT_EQ(box.CellNumber(coords), 15);

  coords.setxyz(0.51, 0.51, -0.51);
  EXPECT_EQ(box.CellNumber(coords), 16);

  coords.setxyz(0.51, 0.51, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 17);

  coords.setxyz(0.51, 0.51, 0.51);
  EXPECT_EQ(box.CellNumber(coords), 18);

  // middle layer, bottom.

  coords.setxyz(0.0, -0.51, -0.51);
  EXPECT_EQ(box.CellNumber(coords), 19);

  coords.setxyz(0.0, -0.51, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 20);

  coords.setxyz(0.0, -0.51, 0.51);
  EXPECT_EQ(box.CellNumber(coords), 21);

  // Middle layer, top.

  coords.setxyz(0.0, 0.51, -0.51);
  EXPECT_EQ(box.CellNumber(coords), 23);

  coords.setxyz(0.0, 0.51, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 24);

  coords.setxyz(0.0, 0.51, 0.51);
  EXPECT_EQ(box.CellNumber(coords), 25);

  // Back and front.

  coords.setxyz(0.0, 0.0, 0.51);
  EXPECT_EQ(box.CellNumber(coords), 22);

  coords.setxyz(0.0, 0.0, -0.51);
  EXPECT_EQ(box.CellNumber(coords), 26);
}

TEST(TestConcentricBox, TestCellsLayer2) {
  constexpr double dx = 1.0;
  coordinate_box::ConcentricBox box(dx);

  Space_Vector<float> coords(-1.8, -1.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 27);

  // Cells on the LHS face, x = -1.8

  coords.setxyz(-1.8, -1.8, -0.8);
  EXPECT_EQ(box.CellNumber(coords), 28);

  coords.setxyz(-1.8, -1.8, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 29);

  coords.setxyz(-1.8, -1.8, 0.8);
  EXPECT_EQ(box.CellNumber(coords), 30);

  coords.setxyz(-1.8, -1.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 31);

  coords.setxyz(-1.8, -0.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 32);

  coords.setxyz(-1.8, -0.8, -0.8);
  EXPECT_EQ(box.CellNumber(coords), 33);

  coords.setxyz(-1.8, -0.8, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 34);

  coords.setxyz(-1.8, -0.8, 0.8);
  EXPECT_EQ(box.CellNumber(coords), 35);

  coords.setxyz(-1.8, -0.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 36);

  coords.setxyz(-1.8, 0.0, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 37);

  coords.setxyz(-1.8, 0.0, -0.8);
  EXPECT_EQ(box.CellNumber(coords), 38);

  coords.setxyz(-1.8, 0.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 42);

  coords.setxyz(-1.8, 1.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 47);

  coords.setxyz(-1.8, 1.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 51);

  // RHS face, x = 1.8

  coords.setxyz(1.8, -1.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 52);

  coords.setxyz(1.8, -1.8, -0.8);
  EXPECT_EQ(box.CellNumber(coords), 53);

  coords.setxyz(1.8, -1.8, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 54);

  coords.setxyz(1.8, -1.8, 0.8);
  EXPECT_EQ(box.CellNumber(coords), 55);

  coords.setxyz(1.8, -1.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 56);

  coords.setxyz(1.8, -0.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 57);

  coords.setxyz(1.8, -0.8, -0.8);
  EXPECT_EQ(box.CellNumber(coords), 58);

  coords.setxyz(1.8, 0.0, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 62);

  coords.setxyz(1.8, 0.0, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 64);

  coords.setxyz(1.8, 1.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 76);

  // Move one layer right, x = -0.8

  coords.setxyz(-0.8, -1.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 77);

  coords.setxyz(-0.8, -1.8, -0.8);
  EXPECT_EQ(box.CellNumber(coords), 78);

  coords.setxyz(-0.8, -1.8, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 79);

  coords.setxyz(-0.8, -1.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 81);

  coords.setxyz(-0.8, -0.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 82);

  coords.setxyz(-0.8, 0.0, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 83);

  coords.setxyz(-0.8, 0.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 84);

  // Top face

  coords.setxyz(-0.8, 1.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 85);

  coords.setxyz(-0.8, 1.8, -0.8);
  EXPECT_EQ(box.CellNumber(coords), 86);

  coords.setxyz(-0.8, 1.8, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 87);

  coords.setxyz(-0.8, 1.8, 0.8);
  EXPECT_EQ(box.CellNumber(coords), 88);

  coords.setxyz(-0.8, 1.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 89);

  // Down the back

  coords.setxyz(-0.8, -0.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 90);

  coords.setxyz(-0.8, 0.0, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 91);

  coords.setxyz(-0.8, 0.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 92);

  // Next layer, x = 0.0

  coords.setxyz(0.0, -1.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 93);

  coords.setxyz(0.0, -1.8, -0.8);
  EXPECT_EQ(box.CellNumber(coords), 94);

  coords.setxyz(0.0, -1.8, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 95);

  coords.setxyz(0.0, -1.8, 0.8);
  EXPECT_EQ(box.CellNumber(coords), 96);

  coords.setxyz(0.0, -1.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 97);

  // up the front.

  coords.setxyz(0.0, -0.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 98);

  coords.setxyz(0.0, 0.0, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 99);

  coords.setxyz(0.0, 0.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 100);

  // Across the top

  coords.setxyz(0.0, 1.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 101);

  coords.setxyz(0.0, 1.8, -0.8);
  EXPECT_EQ(box.CellNumber(coords), 102);

  coords.setxyz(0.0, 1.8, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 103);

  coords.setxyz(0.0, 1.8, 0.8);
  EXPECT_EQ(box.CellNumber(coords), 104);

  coords.setxyz(0.0, 1.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 105);

  // Down the back

  coords.setxyz(0.0, -0.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 106);

  coords.setxyz(0.0, 0.0, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 107);

  coords.setxyz(0.0, 0.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 108);

  // Next X layer, 0.8

  coords.setxyz(0.8, -1.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 109);

  coords.setxyz(0.8, -1.8, -0.8);
  EXPECT_EQ(box.CellNumber(coords), 110);

  coords.setxyz(0.8, -1.8, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 111);

  coords.setxyz(0.8, -1.8, 0.8);
  EXPECT_EQ(box.CellNumber(coords), 112);

  coords.setxyz(0.8, -1.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 113);

  // UP the front

  coords.setxyz(0.8, -0.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 114);

  coords.setxyz(0.8, 0.0, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 115);

  coords.setxyz(0.8, 0.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 116);

  // Across the top

  coords.setxyz(0.8, 1.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 117);

  coords.setxyz(0.8, 1.8, -0.8);
  EXPECT_EQ(box.CellNumber(coords), 118);

  coords.setxyz(0.8, 1.8, 0.0);
  EXPECT_EQ(box.CellNumber(coords), 119);

  coords.setxyz(0.8, 1.8, 0.8);
  EXPECT_EQ(box.CellNumber(coords), 120);

  coords.setxyz(0.8, 1.8, 1.8);
  EXPECT_EQ(box.CellNumber(coords), 121);

  // Down the back

  coords.setxyz(0.8, -0.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 122);

  coords.setxyz(0.8, 0.0, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 123);

  coords.setxyz(0.8, 0.8, -1.8);
  EXPECT_EQ(box.CellNumber(coords), 124);


  std::vector<float> xyz{-1.8, -0.8, 0.0, 0.8, 1.8};
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      for (int k = 0; k < 5; ++k) {
        coords.setxyz(xyz[i], xyz[j], xyz[k]);
      }
    }
  }
}

// Given a coordinate 'v' that gets passed through a ConcentricBox with
// resolution 1.0, what is the expected result when converted back to
// coordinates.
float
Expected(float v) {
  if (v < 0.0f) {
    return -static_cast<float>(static_cast<int>(-v + 0.499999));
  } else {
    return static_cast<float>(static_cast<int>(v + 0.499999));
  }
}

TEST(TestConcentricBox, TestLayer1) {
  constexpr double dx = 1.0;
  coordinate_box::ConcentricBox box(dx);

  Space_Vector<float> coords(0.0, 0.0, 0.0);
  uint64_t cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 0);

  coords.setxyz(0.499, 0.4999, 0.4999);
  EXPECT_EQ(box.CellNumber(coords), 0);
  coords.setxyz(-0.499, -0.4999, -0.4999);
  EXPECT_EQ(box.CellNumber(coords), 0);

  // Individual tests that walk around layer 1.

  coords.setxyz(-0.51, -0.51, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 2);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -1.0, 0.0));

  coords.setxyz(-0.51, -0.51, -0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 1);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -1.0, -1.0));

  coords.setxyz(-0.51, 0.0, -0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 4);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 0.0, -1.0));

  coords.setxyz(-0.51, 0.0, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 5);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 0.0, 0.0));

  coords.setxyz(-0.51, 0.0, 0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 6);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 0.0, 1.0));

  coords.setxyz(-0.51, -0.51, 0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 3);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -1.0, 1.0));

  coords.setxyz(-0.51, 0.51, 0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 9);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 1.0, 1.0));

  coords.setxyz(-0.51, 0.51, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 8);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 1.0, 0.0));

  coords.setxyz(-0.51, 0.51, -0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 7);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 1.0, -1.0));

  // Now the right face.

  coords.setxyz(0.51, -0.51, -0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 10);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, -1.0, -1.0));

  coords.setxyz(0.51, -0.51, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 11);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, -1.0, 0.0));

  coords.setxyz(0.51, -0.51, 0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 12);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, -1.0, 1.0));

  coords.setxyz(0.51, 0.0, -0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 13);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, 0.0, -1.0));

  coords.setxyz(0.51, 0.0, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 14);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, 0.0, 0.0));

  coords.setxyz(0.51, 0.0, 0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 15);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, 0.0, 1.0));

  coords.setxyz(0.51, 0.51, -0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 16);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, 1.0, -1.0));

  coords.setxyz(0.51, 0.51, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 17);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, 1.0, 0.0));

  coords.setxyz(0.51, 0.51, 0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 18);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, 1.0, 1.0));

  // Now start moving along the x axis.

  coords.setxyz(0.0, -0.51, -0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 19);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, -1.0, -1.0));

  coords.setxyz(0.0, -0.51, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 20);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, -1.0, 0.0));

  coords.setxyz(0.0, -0.51, 0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 21);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, -1.0, 1.0));

  coords.setxyz(0.0, 0.51, -0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 23);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, 1.0, -1.0));

  coords.setxyz(0.0, 0.51, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 24);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, 1.0, 0.0));

  coords.setxyz(0.0, 0.51, 0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 25);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, 1.0, 1.0));

  // Now things in the middle.

  coords.setxyz(0.0, 0.0, -0.51);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 26);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, 0.0, -1.0));
}

TEST(TestConcentricBox, TestRandom1) {
  constexpr double dx = 1.0;
  coordinate_box::ConcentricBox box(dx);

  std::default_random_engine rng;
  std::uniform_real_distribution<float> distribution(0.0, 10.0);

  constexpr int kNtest = 100;
  for (int i = 0; i < kNtest; ++i) {
    const float x = distribution(rng);
    const float y = distribution(rng);
    const float z = distribution(rng);
    Space_Vector<float> coords(x, y, z);

    const uint64_t cell = box.CellNumber(coords);

    const Space_Vector<float> returned = box.CellNumberToCoordinates<float>(cell);

    const float x_expected = Expected(x);
    const float y_expected = Expected(y);
    const float z_expected = Expected(z);

    EXPECT_THAT(returned.ToResizableArray(), ElementsAre(x_expected, y_expected, z_expected));
  }
}

TEST(TestConcentricBox, TestP1) {
  coordinate_box::ConcentricBox box(0.01);
  const Space_Vector<float> coords(0.0, 2.6, 15.1);
  const uint64_t cell = box.CellNumber(coords);
  const Space_Vector<float> returned = box.CellNumberToCoordinates<float>(cell);
  EXPECT_THAT(returned.ToResizableArray(), 
              ElementsAre(FloatNear(coords.x(), 1.0e-05),
              FloatNear(coords.y(), 0.5),
              FloatNear(coords.z(), 0.5)));
}

TEST(TestConcentricBox, TestRandom01) {
  coordinate_box::ConcentricBox box(0.01);

  std::default_random_engine rng;
  std::uniform_real_distribution<float> distribution(0.0, 20.0);

  constexpr int kNtest = 100;
  for (int i = 0; i < kNtest; ++i) {
    const float x = distribution(rng);
    const float y = distribution(rng);
    const float z = distribution(rng);
    Space_Vector<float> coords(x, y, z);

    const coordinate_box::LayerPosition lpil = box.Position(coords);

    const Space_Vector<float> returned = box.CellToCoordinates<float>(lpil);

    EXPECT_THAT(returned.ToResizableArray(),
                ElementsAre(FloatNear(x, 0.1), FloatNear(y, 0.1), FloatNear(z, 0.1)));
  }
}

TEST(TestConcentricBox, TestRandomDefResolution) {
  coordinate_box::ConcentricBox box;

  std::default_random_engine rng;
  std::uniform_real_distribution<float> distribution(0.0, 90.0);

  constexpr int kNtest = 1000;
  constexpr float abs_diff = 1.0e-03;
  for (int i = 0; i < kNtest; ++i) {
    const float x = distribution(rng);
    const float y = distribution(rng);
    const float z = distribution(rng);
    Space_Vector<float> coords(x, y, z);

    const coordinate_box::LayerPosition lpil = box.Position(coords);

    const Space_Vector<float> returned = box.CellToCoordinates<float>(lpil);

    EXPECT_THAT(returned.ToResizableArray(), 
                ElementsAre(FloatNear(x, abs_diff), FloatNear(y, abs_diff), FloatNear(z, abs_diff)));
    EXPECT_LT(returned.distance(coords), 1.0e-03);
  }
}

TEST(TestConcentricBox, TestLayer2) {
  constexpr double dx = 1.0;
  coordinate_box::ConcentricBox box(dx);

  Space_Vector<float> coords;
  coords.setxyz(-1.6, -1.6, -1.6);
  uint64_t cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 27);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -2.0, -2.0));

  coords.setxyz(-1.6, -1.6, -0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 28);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -2.0, -1.0));

  coords.setxyz(-1.6, -1.6, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 29);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -2.0, 0.0));

  coords.setxyz(-1.6, -1.6, 0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 30);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -2.0, 1.0));

  coords.setxyz(-1.6, -1.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 31);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -2.0, 2.0));

  coords.setxyz(-1.6, -0.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 32);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -1.0, -2.0));

  coords.setxyz(-1.6, -0.6, -0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 33);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -1.0, -1.0));

  coords.setxyz(-1.6, -0.6, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 34);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -1.0, 0.0));

  coords.setxyz(-1.6, -0.6, 0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 35);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -1.0, 1.0));

  coords.setxyz(-1.6, -0.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 36);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -1.0, 2.0));

  coords.setxyz(-1.6, 0.0, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 37);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 0.0, -2.0));

  coords.setxyz(-1.6, 0.0, -0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 38);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 0.0, -1.0));

  coords.setxyz(-1.6, 0.0, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 39);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 0.0, 0.0));

  coords.setxyz(-1.6, 0.0, 0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 40);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 0.0, 1.0));

  coords.setxyz(-1.6, 0.0, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 41);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 0.0, 2.0));

  coords.setxyz(-1.6, 0.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 42);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 1.0, -2.0));

  coords.setxyz(-1.6, 0.6, -0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 43);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 1.0, -1.0));

  coords.setxyz(-1.6, 0.6, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 44);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 1.0, 0.0));

  coords.setxyz(-1.6, 0.6, 0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 45);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 1.0, 1.0));

  coords.setxyz(-1.6, 0.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 46);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 1.0, 2.0));

  coords.setxyz(-1.6, 1.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 47);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 2.0, -2.0));

  coords.setxyz(-1.6, 1.6, -0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 48);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 2.0, -1.0));

  coords.setxyz(-1.6, 1.6, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 49);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 2.0, 0.0));

  coords.setxyz(-1.6, 1.6, 0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 50);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 2.0, 1.0));

  coords.setxyz(-1.6, 1.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 51);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 2.0, 2.0));

  // Right face

  coords.setxyz(1.6, -1.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 52);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, -2.0, -2.0));

  coords.setxyz(1.6, -1.6, -0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 53);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, -2.0, -1.0));

  coords.setxyz(1.6, -1.6, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 54);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, -2.0, 0.0));

  coords.setxyz(1.6, -1.6, 0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 55);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, -2.0, 1.0));

  coords.setxyz(1.6, -1.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 56);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, -2.0, 2.0));

  coords.setxyz(1.6, -0.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 57);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, -1.0, -2.0));

  coords.setxyz(1.6, -0.6, -0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 58);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, -1.0, -1.0));

  coords.setxyz(1.6, -0.6, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 59);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, -1.0, 0.0));

  coords.setxyz(1.6, -0.6, 0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 60);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, -1.0, 1.0));

  coords.setxyz(1.6, -0.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 61);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, -1.0, 2.0));

  coords.setxyz(1.6, 0.0, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 62);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 0.0, -2.0));

  coords.setxyz(1.6, 0.0, -0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 63);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 0.0, -1.0));

  coords.setxyz(1.6, 0.0, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 64);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 0.0, 0.0));

  coords.setxyz(1.6, 0.0, 0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 65);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 0.0, 1.0));

  coords.setxyz(1.6, 0.0, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 66);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 0.0, 2.0));

  coords.setxyz(1.6, 0.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 67);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 1.0, -2.0));

  coords.setxyz(1.6, 0.6, -0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 68);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 1.0, -1.0));

  coords.setxyz(1.6, 0.6, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 69);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 1.0, 0.0));

  coords.setxyz(1.6, 0.6, 0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 70);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 1.0, 1.0));

  coords.setxyz(1.6, 0.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 71);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 1.0, 2.0));

  coords.setxyz(1.6, 1.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 72);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 2.0, -2.0));

  coords.setxyz(1.6, 1.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 76);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 2.0, 2.0));

  // IN the middle

  coords.setxyz(-0.6, -1.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 77);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -2.0, -2.0));

  coords.setxyz(-0.6, -1.6, -0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 78);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -2.0, -1.0));

  coords.setxyz(-0.6, -1.6, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 79);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -2.0, 0.0));

  coords.setxyz(-0.6, -1.6, 0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 80);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -2.0, 1.0));

  coords.setxyz(-0.6, -1.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 81);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -2.0, 2.0));

  coords.setxyz(-0.6, -0.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 82);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -1.0, 2.0));

  coords.setxyz(-0.6, 0.0, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 83);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 0.0, 2.0));

  coords.setxyz(-0.6, 0.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 84);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 1.0, 2.0));

  coords.setxyz(-0.6, 1.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 85);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 2.0, -2.0));

  coords.setxyz(-0.6, 1.6, -0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 86);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 2.0, -1.0));

  coords.setxyz(-0.6, 1.6, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 87);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 2.0, 0.0));

  coords.setxyz(-0.6, 1.6, 0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 88);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 2.0, 1.0));

  coords.setxyz(-0.6, 1.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 89);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 2.0, 2.0));

  coords.setxyz(-0.6, -0.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 90);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -1.0, -2.0));

  coords.setxyz(-0.6, 0.0, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 91);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 0.0, -2.0));

  coords.setxyz(-0.6, 0.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 92);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 1.0, -2.0));

  // Next slice

  coords.setxyz(0.0, -1.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 93);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, -2.0, -2.0));

  coords.setxyz(0.0, -1.6, -0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 94);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, -2.0, -1.0));

  coords.setxyz(0.0, -1.6, 0.0);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 95);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, -2.0, 0.0));

  coords.setxyz(0.0, -1.6, 0.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 96);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, -2.0, 1.0));

  coords.setxyz(0.0, -0.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 98);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, -1.0, 2.0));

  coords.setxyz(0.0, 0.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 100);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, 1.0, 2.0));

  coords.setxyz(0.0, 1.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 101);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, 2.0, -2.0));

  coords.setxyz(0.0, 1.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 105);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, 2.0, 2.0));

  coords.setxyz(0.0, -0.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 106);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, -1.0, -2.0));

  coords.setxyz(0.0, 0.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 108);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(0.0, 1.0, -2.0));

  // Next slice

  coords.setxyz(0.6, -1.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 109);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, -2.0, -2.0));

  coords.setxyz(0.6, -1.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 113);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, -2.0, 2.0));

  coords.setxyz(0.6, -0.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 114);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, -1.0, 2.0));

  coords.setxyz(0.6, 0.0, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 115);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, 0.0, 2.0));

  coords.setxyz(0.6, 0.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 116);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, 1.0, 2.0));

  coords.setxyz(0.6, 1.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 117);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, 2.0, -2.0));

  coords.setxyz(0.6, 1.6, 1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 121);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, 2.0, 2.0));

  coords.setxyz(0.6, -0.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 122);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, -1.0, -2.0));

  coords.setxyz(0.6, 0.6, -1.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 124);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(1.0, 1.0, -2.0));
}

TEST(TestConcentricBox, TestLayer3) {
  constexpr double dx = 1.0;
  coordinate_box::ConcentricBox box(dx);

  Space_Vector<float> coords;
  coords.setxyz(-2.6, -2.6, -2.6);
  uint64_t cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 125);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-3.0, -3.0, -3.0));

  coords.setxyz(-2.6, -2.6, 2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 131);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-3.0, -3.0, 3.0));

  coords.setxyz(-2.6, 2.6, 2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 173);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-3.0, 3.0, 3.0));

  coords.setxyz(2.6, -2.6, -2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 174);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(3.0, -3.0, -3.0));

  coords.setxyz(2.6, 2.6, 2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 222);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(3.0, 3.0, 3.0));

  coords.setxyz(-1.6, -2.6, -2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 223);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -3.0, -3.0));

  coords.setxyz(-1.6, -2.6, 2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 229);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -3.0, 3.0));

  coords.setxyz(-1.6, -1.6, 2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 230);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -2.0, 3.0));

  coords.setxyz(-1.6, 1.6, 2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 234);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 2.0, 3.0));

  coords.setxyz(-1.6, 2.6, -2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 235);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 3.0, -3.0));

  coords.setxyz(-1.6, 2.6, 2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 241);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 3.0, 3.0));

  coords.setxyz(-1.6, -1.6, -2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 242);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, -2.0, -3.0));

  coords.setxyz(-1.6, 1.6, -2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 246);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-2.0, 2.0, -3.0));

  // Next X slice

  coords.setxyz(-0.6, -2.6, -2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 247);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -3.0, -3.0));

  coords.setxyz(-0.6, -2.6, 2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 253);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -3.0, 3.0));

  coords.setxyz(-0.6, -1.6, 2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 254);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, -2.0, 3.0));

  coords.setxyz(-0.6, 0.0, 2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 256);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 0.0, 3.0));

  coords.setxyz(-0.6, 1.6, 2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 258);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-1.0, 2.0, 3.0));

  coords.setxyz(1.6, 1.6, -2.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 342);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(2.0, 2.0, -3.0));
}

TEST(TestConcentricBox, TestLayer4) {
  constexpr double dx = 1.0;
  coordinate_box::ConcentricBox box(dx);

  Space_Vector<float> coords;
  coords.setxyz(-3.6, -3.6, -3.6);
  uint64_t cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 343);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-4.0, -4.0, -4.0));

  coords.setxyz(3.6, 3.6, 3.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 504);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(4.0, 4.0, 4.0));

  coords.setxyz(2.6, 2.6, 3.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 712);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(3.0, 3.0, 4.0));

  coords.setxyz(2.6, -2.6, -3.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 722);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(3.0, -3.0, -4.0));

  coords.setxyz(2.6, 2.6, -3.6);
  cell = box.CellNumber(coords);
  EXPECT_EQ(cell, 728);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(3.0, 3.0, -4.0));
}

TEST(TestConcentricBox, TestLayer5) {
  constexpr double dx = 1.0;
  coordinate_box::ConcentricBox box(dx);

  Space_Vector<float> coords;
  coords.setxyz(-4.6, -4.6, -4.6);
  uint64_t cell = box.CellNumber(coords);
  EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(), ElementsAre(-5.0, -5.0, -5.0));

  for (float x = -4.6f; x < 5.0f; x += 1.0f) {
    for (float y = -4.6f; y < 5.0f; y += 1.0f) {
      for (float z = -4.6f; z < 5.0f; z += 1.0f) {
        const Space_Vector<float> coords(x, y, z);
        float x_expected = Expected(x);
        float y_expected = Expected(y);
        float z_expected = Expected(z);
        const uint64_t cell = box.CellNumber(coords);
        EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(),
                    ElementsAre(x_expected, y_expected, z_expected));
      }
    }
  }
}

TEST(TestConcentricBox, TestLayer9) {
  constexpr double dx = 1.0;
  coordinate_box::ConcentricBox box(dx);

  for (float x = -9.6f; x < 9.0f; x += 1.0f) {
    for (float y = -9.6f; y < 9.0f; y += 1.0f) {
      for (float z = -9.6f; z < 9.0f; z += 1.0f) {
        const Space_Vector<float> coords(x, y, z);
        float x_expected = Expected(x);
        float y_expected = Expected(y);
        float z_expected = Expected(z);
        const uint64_t cell = box.CellNumber(coords);
        EXPECT_THAT(box.CellNumberToCoordinates<float>(cell).ToResizableArray(),
                    ElementsAre(x_expected, y_expected, z_expected));
      }
    }
  }
}

TEST(TestFromString, Empty) {
  const const_IWSubstring s;
  coordinate_box::LayerPosition layer_position;
  EXPECT_EQ(coordinate_box::FromString(s, layer_position), 0);
}

TEST(TestFromString, JustSeparator) {
  const const_IWSubstring s(":");
  coordinate_box::LayerPosition layer_position;
  EXPECT_EQ(coordinate_box::FromString(s, layer_position), 0);
}

TEST(TestFromString, EmptyLayer) {
  const const_IWSubstring s(":0");
  coordinate_box::LayerPosition layer_position;
  EXPECT_EQ(coordinate_box::FromString(s, layer_position), 0);
}

TEST(TestFromString, EmptyPosition) {
  const const_IWSubstring s("0:");
  coordinate_box::LayerPosition layer_position;
  EXPECT_EQ(coordinate_box::FromString(s, layer_position), 0);
}

TEST(TestFromString, BothZero) {
  const const_IWSubstring s("0:0");
  coordinate_box::LayerPosition layer_position;
  ASSERT_EQ(coordinate_box::FromString(s, layer_position), s.length());
  EXPECT_EQ(layer_position.layer, 0);
  EXPECT_EQ(layer_position.position_in_layer, 0);
}

TEST(TestFromString, BothZeroWithTrailing) {
  const const_IWSubstring s("0:0:::");
  coordinate_box::LayerPosition layer_position;
  EXPECT_EQ(coordinate_box::FromString(s, layer_position), 3);
  EXPECT_EQ(layer_position.layer, 0);
  EXPECT_EQ(layer_position.position_in_layer, 0);
}

TEST(TestFromString, ArbitraryNumbers) {
  for (uint32_t layer = 17; layer < 8000; layer += 103) {
    uint64_t within_layer = 123456;
    for (int i = 0; i < 100; ++i) {
      within_layer *= 987654;
      IWString s;
      s << layer << ':' << within_layer;
      coordinate_box::LayerPosition layer_position;
      EXPECT_EQ(coordinate_box::FromString(s, layer_position), s.length());
      EXPECT_EQ(layer_position.layer, layer);
      EXPECT_EQ(layer_position.position_in_layer, within_layer);
    }
  }
}

TEST(TestFromString, RoundTrip) {
  for (uint32_t layer = 17; layer < 8000; layer += 10399) {
    uint64_t within_layer = 1234567;
    for (int i = 0; i < 100; ++i) {
      within_layer *= 987654;
      coordinate_box::LayerPosition layer_position(layer, within_layer);
      IWString s;
      s << layer_position;
      coordinate_box::LayerPosition returned;
      EXPECT_EQ(coordinate_box::FromString(s, returned), s.length());
      EXPECT_EQ(returned.layer, layer);
      EXPECT_EQ(returned.position_in_layer, within_layer);
    }
  }
}

}  // namespace
