// Tests for known_fragment_data

#include <filesystem>
#include <iostream>

#include "Foundational/iwstring/iwstring.h"
#include "Molecule_Tools/known_fragment_data.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {

using std::cerr;

// Write `contents` to `fname`.
int
Write(const IWString& contents,
      const char* fname) {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    std::cerr << "Write:cannot open '" << fname << "'\n";
    return 0;
  }
  
  output << contents;
  if (!contents.ends_with('\n')) {
    output << '\n';
  }

  return 1;
}

struct SmilesFragSpec {
  IWString smiles;

  // If we need data in a temp file.
  std::string fname;
  IWString file_contents;

  IWString result;
};

class TestFragments : public testing::TestWithParam<SmilesFragSpec> {
  protected:
    Molecule _m;
    Known_Fragment_Data _known_fragment_data;
};

TEST_P(TestFragments, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  auto path = std::filesystem::temp_directory_path();
  path.append(params.fname);

  ASSERT_TRUE(Write(params.file_contents, path.c_str()));
  ASSERT_TRUE(_known_fragment_data.ReadKnownSaltsSmartsFile(IWString(path)));
  _known_fragment_data.process(_m);
  EXPECT_EQ(_m.smiles(), params.result) << "Got " << _m.smiles() << 
        " expected " << params.result;
}
INSTANTIATE_TEST_SUITE_P(TestFragments, TestFragments, testing::Values(
  SmilesFragSpec{"C.[Na].[Cu]", "t1", "[Na,Cu]", "C"},
  SmilesFragSpec{"CC.C", "t2", "[CD1]C", "C"},
  SmilesFragSpec{"CC.C.CCC", "t3", "[CD1]C", "C.CCC"}
));


class TestKnownParents : public testing::TestWithParam<SmilesFragSpec> {
  protected:
    Molecule _m;
    Known_Fragment_Data _known_fragment_data;
};

TEST_P(TestKnownParents, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  auto path = std::filesystem::temp_directory_path();
  path.append(params.fname);

  ASSERT_TRUE(Write(params.file_contents, path.c_str()));
  ASSERT_TRUE(_known_fragment_data.read_known_parents(IWString(path)));
  _known_fragment_data.process(_m);
  EXPECT_EQ(_m.smiles(), params.result) << "Got " << _m.smiles() << 
        " expected " << params.result;
}
INSTANTIATE_TEST_SUITE_P(TestKnownParents, TestKnownParents, testing::Values(
  SmilesFragSpec{"C.CC", "t1", "O", "C.CC"},
  SmilesFragSpec{"C.CC", "t2", "C", "C"},
  SmilesFragSpec{"C.C.CC", "t3", "C", "C.C"}
));

class TestKnownSalts : public testing::TestWithParam<SmilesFragSpec> {
  protected:
    Molecule _m;
    Known_Fragment_Data _known_fragment_data;
};

TEST_P(TestKnownSalts, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  auto path = std::filesystem::temp_directory_path();
  path.append(params.fname);

  ASSERT_TRUE(Write(params.file_contents, path.c_str()));
  ASSERT_TRUE(_known_fragment_data.read_known_salts(IWString(path)));
  _known_fragment_data.process(_m);
  EXPECT_EQ(_m.smiles(), params.result) << "Got " << _m.smiles() << 
        " expected " << params.result;
}
INSTANTIATE_TEST_SUITE_P(TestKnownSalts, TestKnownSalts, testing::Values(
  SmilesFragSpec{"C.CC", "t1", "O", "C.CC"},
  SmilesFragSpec{"C.CC", "t2", "C", "CC"},
  SmilesFragSpec{"C.CC", "t3", "CC", "C"},
  SmilesFragSpec{"CC.C.CC", "t4", "CC", "C"}
));


}  // namespace
