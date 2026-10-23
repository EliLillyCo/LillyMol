#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

#include "gtest/gtest.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Tools/iwdescr_lib.h"

namespace {

int
DescriptorIndex(const IWDescr& iwdescr, const char* name) {
  for (int i = 0; i < iwdescr.number_descriptors(); ++i) {
    if (iwdescr.descriptor_name(i) == name) {
      return i;
    }
  }

  return -1;
}

float
DescriptorValue(const IWDescr& iwdescr, const std::vector<float>& values, const char* name) {
  const int index = DescriptorIndex(iwdescr, name);
  EXPECT_GE(index, 0) << name;
  if (index < 0) {
    return std::numeric_limits<float>::quiet_NaN();
  }

  return values[index];
}

TEST(IWDescr, InitialiseAllAndProcess) {
  if (getenv("LILLYMOL_HOME") == nullptr) {
    GTEST_SKIP() << "LILLYMOL_HOME is not defined";
  }

  IWDescr iwdescr;
  ASSERT_TRUE(iwdescr.InitialiseAll());
  ASSERT_GT(iwdescr.number_descriptors(), 0);

  std::unordered_set<std::string> names;
  for (int i = 0; i < iwdescr.number_descriptors(); ++i) {
    const std::string name = iwdescr.descriptor_name(i).AsString();
    EXPECT_FALSE(name.empty()) << i;
    EXPECT_TRUE(names.insert(name).second) << name;
  }

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));
  std::vector<float> values(iwdescr.number_descriptors());
  ASSERT_TRUE(iwdescr.Process(m, values.data()));

  int finite_values = 0;
  for (float value : values) {
    finite_values += std::isfinite(value);
  }
  EXPECT_GT(finite_values, 0);
}

TEST(IWDescr, SpinachRotatableBondDescriptors) {
  if (getenv("LILLYMOL_HOME") == nullptr) {
    GTEST_SKIP() << "LILLYMOL_HOME is not defined";
  }

  IWDescr iwdescr;
  ASSERT_TRUE(iwdescr.InitialiseAll());

  Molecule biphenyl;
  ASSERT_TRUE(biphenyl.build_from_smiles("c1ccccc1-c1ccccc1"));
  std::vector<float> biphenyl_values(iwdescr.number_descriptors());
  ASSERT_TRUE(iwdescr.Process(biphenyl, biphenyl_values.data()));
  EXPECT_EQ(DescriptorValue(iwdescr, biphenyl_values, "rotbond"), 1.0f);
  EXPECT_EQ(DescriptorValue(iwdescr, biphenyl_values, "scafrotb"), 1.0f);
  EXPECT_EQ(DescriptorValue(iwdescr, biphenyl_values, "spchrotb"), 0.0f);
  EXPECT_EQ(DescriptorValue(iwdescr, biphenyl_values, "maxrotbgrp"), 1.0f);

  Molecule hexane;
  ASSERT_TRUE(hexane.build_from_smiles("CCCCCC"));
  std::vector<float> hexane_values(iwdescr.number_descriptors());
  ASSERT_TRUE(iwdescr.Process(hexane, hexane_values.data()));
  EXPECT_EQ(DescriptorValue(iwdescr, hexane_values, "rotbond"), 3.0f);
  EXPECT_EQ(DescriptorValue(iwdescr, hexane_values, "scafrotb"), 0.0f);
  EXPECT_EQ(DescriptorValue(iwdescr, hexane_values, "spchrotb"), 3.0f);
  EXPECT_EQ(DescriptorValue(iwdescr, hexane_values, "maxrotbgrp"), 3.0f);
}

}  // namespace
