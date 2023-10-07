// Tests for FixedbitVector

#include <memory>

#include "fixed_bit_vector.h"

//#include "googletest/include/gtest/gtest.h"
#include <gtest/gtest.h>
#include "google/protobuf/text_format.h"

namespace {
using fixed_bit_vector::FixedBitVector;

// Really just a test that nothing crashes.
TEST(TestFixedBitVector, Empty) {
  FixedBitVector foo;
  EXPECT_EQ(1, 1);
}

TEST(TestFixedBitVector, ResizeFromNothing1) {
  FixedBitVector foo;
  foo.resize(3);
  EXPECT_EQ(foo.nbits(), 64);
  for (int b = 0; b < 64; ++b) {
    EXPECT_FALSE(foo.is_set(b));
  }
}

TEST(TestFixedBitVector, ResizeFromNothing2) {
  FixedBitVector foo;
  foo.resize(125);
  EXPECT_EQ(foo.nbits(), 64 * 2);
  for (int b = 0; b < 64 * 2; ++b) {
    EXPECT_FALSE(foo.is_set(b));
  }
}

TEST(TestFixedBitVector, ResizeFromNothingExact) {
  FixedBitVector foo;
  foo.resize(128);
  EXPECT_EQ(foo.nbits(), 64 * 2);
  for (int b = 0; b < 64 * 2; ++b) {
    EXPECT_FALSE(foo.is_set(b));
  }
}

TEST(TestFixedBitVector, BitsSpecifiedInConstructor) {
  FixedBitVector foo(231);
  EXPECT_EQ(foo.nbits(), 64 * 4);
  for (int b = 0; b < 64 * 4; ++b) {
    EXPECT_FALSE(foo.is_set(b));
  }
}

TEST(TestFixedBitVector, BitsSpecifiedInConstructorExact) {
  FixedBitVector foo(256);
  EXPECT_EQ(foo.nbits(), 64 * 4);
  for (int b = 0; b < 64 * 4; ++b) {
    EXPECT_FALSE(foo.is_set(b));
  }
}

TEST(TestFixedBitVector, TestSet) {
  constexpr int nbits = 256;
  FixedBitVector foo(nbits);
  for (int i = 0; i < nbits; ++i) {
    foo.set_bit(i);
    for (int j = 0; j < nbits; ++j) {
      if (i == j)
        EXPECT_TRUE(foo.is_set(j));
      else
        EXPECT_FALSE(foo.is_set(j));
    }
    foo.unset_bit(i);
  }
}

TEST(TestFixedBitVector, TestBitsInCommon) {
  constexpr int nbits = 128;
  FixedBitVector f1(nbits), f2(nbits);
  for (int i = 0; i < nbits; ++i) {
    f1.set_bit(i);
    for (int j = 0; j < nbits; ++j) {
      f2.set_bit(j);
      if (i == j) {
        EXPECT_EQ(f1.BitsInCommon(f2), 1);
      } else {
        EXPECT_EQ(f1.BitsInCommon(f2), 0);
      }
      f2.unset_bit(j);
    }
    f1.unset_bit(i);
  }
}

TEST(TestFixedBitVector, TestNsetSingleBot) {
  constexpr int nbits = 512;
  FixedBitVector foo(nbits);
  for (int i = 0; i < nbits; ++i) {
    foo.set_bit(i);
    EXPECT_EQ(foo.nset(), 1);
    foo.unset_bit(i);
  }
}

TEST(TestFixedBitVector, TestNsetAll) {
  constexpr int nbits = 2048;
  FixedBitVector foo(nbits);
  for (int i = 0; i < nbits; ++i) {
    foo.set_bit(i);
    EXPECT_EQ(foo.nset(), i + 1);
  }
}

TEST(TestLog2, TestSingleBitSet) {
  constexpr uint64_t one = 1;
  for (int i = 1; i < 64; ++i) {
    uint64_t v = one << i;
    EXPECT_EQ(fixed_bit_vector::first_bit_set(v), i);
  }
}

TEST(TestLog2, TestOtherValues) {
  constexpr uint64_t one = 1;
  for (int i = 1; i < 63; ++i) {
    uint64_t v = one << i;
    EXPECT_EQ(fixed_bit_vector::first_bit_set(v), i);
    v |= one << (i+1);
    EXPECT_EQ(fixed_bit_vector::first_bit_set(v), i);
    v |= one << (i-1);
    EXPECT_EQ(fixed_bit_vector::first_bit_set(v), i - 1);
  }
}

TEST(TestFixedBitVector, TestFirstBitSetSingleValue) {
  constexpr int nbits = 2048;
  FixedBitVector foo(nbits);
  for (int i = 0; i < nbits; ++i) {
    foo.set_bit(i);
    EXPECT_EQ(foo.FirstBitSet(), i);
    foo.reset();
  }
}

TEST(TestLog2, TestUnsetSingleValues) {
  constexpr uint64_t one = 1;
  EXPECT_EQ(fixed_bit_vector::first_unset_bit(0), 0);
  EXPECT_LT(fixed_bit_vector::first_unset_bit(std::numeric_limits<uint64_t>::max()), 0);

  for (int i = 1; i < 64; ++i) {
    uint64_t v = one << i;
    EXPECT_EQ(fixed_bit_vector::first_unset_bit(v), 0);
  }
}

TEST(TestLog2, TestUnsetSingleUnset) {
  constexpr uint64_t one = 1;
  for (int i = 0; i < 64; ++i) {
    uint64_t v = std::numeric_limits<uint64_t>::max();
    v ^= (one << i);
    EXPECT_EQ(fixed_bit_vector::first_unset_bit(v), i);
  }
}

TEST(TestLog2, TestUnsetManyUnset) {
  constexpr uint64_t one = 1;
  for (int i = 0; i < 64; ++i) {
    uint64_t v = std::numeric_limits<uint64_t>::max();
    v ^= (one << i);
    for (int j = i + 1; j < 64; j += 3) {
      v ^= (one << j);
    }
    EXPECT_EQ(fixed_bit_vector::first_unset_bit(v), i);
  }
}

TEST(TestFixedBitVector, TestFirstBitMultipleValues) {
  constexpr int nbits = 1024;
  for (int i = 0; i < nbits - 64; ++i) {
    FixedBitVector foo(nbits);
    foo.set_bit(i);
    foo.set_bit(i + 32);
    foo.set_bit(i + 33);
    EXPECT_EQ(foo.FirstBitSet(), i);
  }
}

TEST(TestFixedBitVector, TestOperatorLtLtSingleBit) {
  constexpr int nbits = 64;
  constexpr uint64_t one = 1;
  for (int i = 0; i < nbits; ++i) {
    FixedBitVector foo(nbits);
    foo.set_bit(i);
    std::stringstream ss;
    ss << foo;
    std::stringstream expected;
    expected << "FixedBitVector 64 " << std::hex << (one << i);
    EXPECT_EQ(ss.str(), expected.str());
  }
}

TEST(TestFixedBitVector, TestBuildFromHexEmpty) {
  FixedBitVector foo(64);
  foo.set_bit(34);
  IWString empty_hex;
  EXPECT_TRUE(foo.ConstructFromHex(empty_hex));
  EXPECT_EQ(foo.nbits(), 0);
}

TEST(TestFixedBitVector, TestBuildFromHexOneByte) {
  FixedBitVector foo(64);
  foo.set_bit(34);
  IWString hex("00");
  EXPECT_TRUE(foo.ConstructFromHex(hex));
  EXPECT_EQ(foo.nwords(), 1);
  EXPECT_EQ(foo.nset(), 0);
}

char
HexChar(int i) {
  if (i < 10) {
    return '0' + i;
  } else {
    return 'a' + i - 10;
  }
}

TEST(TestFixedBitVector, TestAllHex) {
  for (int i = 0; i < 16; ++i) {
    for (int j = 0; j < 16; ++j) {
      IWString zhex;
      zhex << HexChar(i);
      zhex << HexChar(j);
      FixedBitVector foo;
      ASSERT_TRUE(foo.ConstructFromHex(zhex));
      EXPECT_EQ(foo.nwords(), 1);
      if (j % 2 == 1) {
        EXPECT_TRUE(foo.is_set(0));
      } else {
        EXPECT_FALSE(foo.is_set(0));
      }
      EXPECT_EQ(16 * i + j, foo.bits()[0]);
    }
  }
}

TEST(TestFixedBitVector, TwoWordsFromHex) {
  IWString zhex("0000000000000000");
  FixedBitVector foo;
  ASSERT_TRUE(foo.ConstructFromHex(zhex));
  EXPECT_EQ(foo.nwords(), 1);
  zhex << "00";
  ASSERT_TRUE(foo.ConstructFromHex(zhex));
  EXPECT_EQ(foo.nwords(), 2);
  EXPECT_EQ(foo.nset(), 0);
}

TEST(TestFixedBitVector, TestDaylight) {
  resizable_array<int> bits_set;
  FixedBitVector foo(257);
  bits_set << 7 << 31 << 65 << 67 << 68 << 120 << 151 << 191 << 208 << 251 << 256;
  for (auto bit : bits_set) {
    foo.set_bit(bit);
  }

  const IWString daylight = foo.DaylightAsciiRepresentation();
  FixedBitVector round_trip;

  ASSERT_TRUE(round_trip.ConstructFromDaylightAsciiBitRep(daylight));
  EXPECT_EQ(foo, round_trip);
}

TEST(TestFixedBitVector, TestFold) {
  resizable_array<int> bits_set;
  FixedBitVector foo(257);
  bits_set << 7 << 31 << 65 << 67 << 68 << 120 << 151 << 191 << 208 << 251 << 256;
  foo.set_bits(bits_set);

  EXPECT_EQ(foo.nbits(), 320);
  foo.Fold(1);

  EXPECT_EQ(foo.nbits(), 160);

  for (auto bit : bits_set) {
    EXPECT_TRUE(foo.is_set(bit) || foo.is_set(bit / 2));
  }
}

TEST(TestFixedBitVector, TestConcatenate) {
  FixedBitVector foo1(128);
  FixedBitVector foo2(256);

  resizable_array<int> bits1, bits2;
  bits1 << 9 << 127 << 44 << 4 << 119 << 118 << 121 << 64 << 65 << 3 << 33;
  bits2 << 46 << 93 << 27 << 73 << 12 << 43 << 89 << 5 << 37 << 85 << 25;

  foo1.set_bits(bits1);
  foo2.set_bits(bits2);

  EXPECT_EQ(foo1.nbits(), 128);
  EXPECT_EQ(foo2.nbits(), 256);

  foo1 += foo2;
  EXPECT_EQ(foo1.nbits(), 128 + 256);
  for (auto b : bits1) {
    EXPECT_TRUE(foo1.is_set(b));
  }
  for (auto b : bits2) {
    EXPECT_TRUE(foo1.is_set(128 + b));
  }
}

TEST(TestFixedBitVector, TestAscii01) {
  IWString ascii;
  ascii.extend(128, '0');
  resizable_array<int> bits;
  bits << 0 << 1 << 2 << 33 << 46 << 93 << 27 << 117 << 116 << 127 << 3 << 99 << 101;
  for (int b : bits) {
    ascii[b] = '1';
  }
  FixedBitVector foo;
  ASSERT_TRUE(foo.ConstructFromAscii01Representation(ascii));
  EXPECT_EQ(foo.nbits(), 128);
  EXPECT_EQ(foo.nset(), bits.number_elements());

  for (int b : bits) {
    EXPECT_TRUE(foo.is_set(b));
  }
}

TEST(TestFixedBitVector, TestFlip) {
  resizable_array<int> bits;
  bits << 0 << 1 << 2 << 33 << 46 << 93 << 27 << 117 << 116 << 127 << 3 << 99 << 101;
  FixedBitVector foo(128);
  for (int b : bits) {
    foo.set_bit(b);
  }
  foo.Flip();
  for (int i = 0; i < 128; ++i) {
    if (bits.contains(i)) {
      EXPECT_FALSE(foo.is_set(i));
    } else {
      EXPECT_TRUE(foo.is_set(i));
    }
  }
}
TEST(TestFixedBitVector, TestOr) {
  resizable_array<int> common_bits;
  common_bits << 0 << 1 << 2 << 33 << 46 << 93 << 27 << 117 << 116;
  resizable_array<int> extra1, extra2;
  extra1 << 101 << 102 << 103 << 127;
  extra2 << 104 << 105 << 106 << 121;

  FixedBitVector foo1(128);
  FixedBitVector foo2(128);

  foo1.set_bits(common_bits);
  foo2.set_bits(common_bits);
  foo1.set_bits(extra1);
  foo2.set_bits(extra2);

  foo1.iwor(foo2);
  EXPECT_EQ(foo1.nset(), common_bits.number_elements() + extra1.number_elements() + extra2.number_elements());

  for (int b : common_bits) {
    EXPECT_TRUE(foo1.is_set(b));
  }
  for (int b : extra1) {
    EXPECT_TRUE(foo1.is_set(b));
  }
  for (int b : extra2) {
    EXPECT_TRUE(foo1.is_set(b));
  }
}
TEST(TestFixedBitVector, TestOrChanged) {
  resizable_array<int> common_bits;
  common_bits << 0 << 1 << 2 << 33 << 46 << 93 << 27 << 117 << 116;
  resizable_array<int> extra1, extra2;
  extra1 << 101 << 102 << 103 << 127;
  extra2 << 104 << 105 << 106 << 121;

  FixedBitVector foo1(128);
  FixedBitVector foo2(128);

  foo1.set_bits(common_bits);
  foo2.set_bits(common_bits);
  foo1.set_bits(extra1);
  foo2.set_bits(extra2);

  int changed;
  foo1.iwor(foo2, changed);
  EXPECT_TRUE(changed);

  FixedBitVector foo3(foo1);
  foo1.iwor(foo3, changed);
  EXPECT_FALSE(changed);
}

TEST(TestFixedBitVector, TestAnd) {
  resizable_array<int> common_bits;
  common_bits << 0 << 1 << 2 << 33 << 46 << 93 << 27 << 117 << 116;
  resizable_array<int> extra1, extra2;
  extra1 << 101 << 102 << 103 << 127;
  extra2 << 104 << 105 << 106 << 121;

  FixedBitVector foo1(128);
  FixedBitVector foo2(128);

  foo1.set_bits(common_bits);
  foo2.set_bits(common_bits);
  foo1.set_bits(extra1);
  foo2.set_bits(extra2);

  foo1.iwand(foo2);
  EXPECT_EQ(foo1.nset(), common_bits.number_elements());
  for (int b : common_bits) {
    EXPECT_TRUE(foo1.is_set(b));
  }
}

TEST(TestFixedBitVector, TestAndChanged) {
  resizable_array<int> common_bits;
  common_bits << 0 << 1 << 2 << 33 << 46 << 93 << 27 << 117 << 116;
  resizable_array<int> extra1, extra2;
  extra1 << 101 << 102 << 103 << 127;
  extra2 << 104 << 105 << 106 << 121;

  FixedBitVector foo1(128);
  FixedBitVector foo2(128);

  foo1.set_bits(common_bits);
  foo2.set_bits(common_bits);
  foo1.set_bits(extra1);
  foo2.set_bits(extra2);

  int changed = 0;
  foo1.iwand(foo2, changed);
  EXPECT_TRUE(changed);

  FixedBitVector foo3(foo1);
  foo1.iwand(foo3, changed);
  EXPECT_FALSE(changed);
}

TEST(TestFixedBitVector, TestXor) {
  resizable_array<int> common_bits;
  common_bits << 0 << 1 << 2 << 33 << 46 << 93 << 27 << 117 << 116;
  resizable_array<int> extra1, extra2;
  extra1 << 101 << 102 << 103 << 127;
  extra2 << 104 << 105 << 106 << 121;

  FixedBitVector foo1(128);
  FixedBitVector foo2(128);

  foo1.set_bits(common_bits);
  foo2.set_bits(common_bits);
  foo1.set_bits(extra1);
  foo2.set_bits(extra2);

  foo1.iwxor(foo2);
  EXPECT_EQ(foo1.nset(), extra1.number_elements() + extra2.number_elements());
  for (int b : common_bits) {
    EXPECT_FALSE(foo1.is_set(b));
  }
  for (int b : extra1) {
    EXPECT_TRUE(foo1.is_set(b));
  }
  for (int b : extra2) {
    EXPECT_TRUE(foo1.is_set(b));
  }
}

TEST(TestFixedBitVector, TestConstructFromTdtRecordNset) {
  resizable_array<int> common_bits;
  common_bits << 0 << 1 << 2 << 33 << 46 << 93 << 27 << 117 << 116;
  constexpr int nbits = 256;
  FixedBitVector foo1(nbits);
  foo1.set_bits(common_bits);
  IWString tdt_record;
  tdt_record << "FPX<" << foo1.DaylightAsciiRepresentationIncludingNsetInfo() << '>';

  FixedBitVector foo2;
  int nset2;
  ASSERT_TRUE(foo2.ConstructFromTdtRecordNset(tdt_record, nset2));
  EXPECT_EQ(foo2.nset(), common_bits.number_elements());
  EXPECT_EQ(foo1, foo2);
}

TEST(TestFixedBitVector, TestConstructFromAscii01RepresentationWithSpaces) {
  resizable_array<int> common_bits;
  common_bits << 0 << 1 << 33 << 46 << 93 << 27 << 117 << 116;

  IWString ascii;
  ascii.extend(256, ' ');
  for (int i = 0; i < 256; i += 2) {
    ascii[i] = '0';
  }
  for (int b : common_bits) {
    int position =  b * 2;
    ascii[position] = '1';
  }

  FixedBitVector foo;
  ASSERT_TRUE(foo.ConstructFromAscii01RepresentationWithSpaces(ascii));
  EXPECT_EQ(foo.nset(), common_bits.number_elements());
  for (int b : common_bits) {
    EXPECT_TRUE(foo.is_set(b));
  }
}

TEST(TestFixedBitVector, TestMake01NoSpace) {
  resizable_array<int> common_bits;
  common_bits << 0 << 1 << 33 << 46 << 93 << 27 << 117 << 116;

  FixedBitVector foo1(128);
  foo1.set_bits(common_bits);

  const IWString b01 = foo1.To01tf('1', '0');
  FixedBitVector foo2;
  ASSERT_TRUE(foo2.ConstructFromAscii01Representation(b01));
  EXPECT_EQ(foo1, foo2);
}

TEST(TestFixedBitVector, TestOperatorEquals) {
  resizable_array<int> common_bits;
  common_bits << 0 << 1 << 33 << 46 << 93 << 27 << 117 << 116;
  FixedBitVector foo1(256);
  for (int b : common_bits) {
    foo1.set_bit(b);
  }

  FixedBitVector foo2;
  foo2 = foo1;
  EXPECT_TRUE(foo1 == foo2);
  EXPECT_EQ(foo1, foo2);
}

TEST(TestFixedBitVector, TestNextOnBitEmpty) {
  FixedBitVector foo(1024);

  int number_retrieved = 0;
  int i = 0;
  int bit;
  while ((bit = foo.NextOnBit(i)) >= 0) {
    number_retrieved++;
  }

  EXPECT_EQ(number_retrieved, 0);
}

TEST(TestFixedBitVector, TestNextOnBitAllSet) {
  constexpr int nbits = 1024;
  FixedBitVector foo(nbits);
  for (int i = 0; i < nbits; ++i) {
    foo.set_bit(i);
  }

  int number_retrieved = 0;
  int i = 0;
  int bit;
  while ((bit = foo.NextOnBit(i)) >= 0) {
    number_retrieved++;
  }

  EXPECT_EQ(number_retrieved, nbits);
}

TEST(TestFixedBitVector, TestNextOnBit) {
  resizable_array<int> common_bits;
  common_bits << 0 << 1 << 27 << 33 << 46 << 93 << 117 << 116 << 127 << 200 << 255;
  FixedBitVector foo1(256);
  for (int b : common_bits) {
    foo1.set_bit(b);
  }
  EXPECT_EQ(foo1.nset(), common_bits.number_elements());

  int number_retrieved = 0;
  int i = 0;
  int bit;
  while ((bit = foo1.NextOnBit(i)) >= 0) {
    EXPECT_TRUE(common_bits.contains(bit));
    number_retrieved++;
  }

  EXPECT_EQ(number_retrieved, common_bits.number_elements());
}

TEST(TestFixedBitVector, TestFromCongiguousStorage) {
  resizable_array<int> common_bits;
  common_bits << 0 << 1 << 27 << 33 << 46 << 93 << 117 << 116 << 127 << 200 << 255;
  FixedBitVector foo1(256);
  foo1.set_bits(common_bits);
  const int nw = 256 / 64;

  std::unique_ptr<uint64_t[]> storage(new uint64_t[nw]);
  void * ptr = foo1.CopyToContiguousStorage(storage.get());
  EXPECT_EQ(ptr, storage.get() +  nw);

  FixedBitVector foo2;
  ptr = foo2.BuildFromContiguousStorage(storage.get(), nw);
  EXPECT_EQ(ptr, storage.get() + nw);
  EXPECT_EQ(foo1, foo2);
}

}  // namespace
