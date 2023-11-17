// Tester for fileconv

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Lib/istream_and_type.h"

#include "fileconv_opts.h"

namespace {

using std::vector;

class FileconvTest : public testing::Test {
  protected:
    Molecule _m;
    fileconv::FileconvConfig _config;
    void SetUp();
    void TearDown();
};

void
FileconvTest::SetUp() {
}

void
FileconvTest::TearDown() {
}

TEST_F(FileconvTest, TestLongestPathPasses) {
  int argc = 3;
  const char * argv[] = {"fileconv", "-Y", "maxplen=4"};
  Command_Line cl(argc, const_cast<char**>(argv), "Y:");
  ASSERT_TRUE(_config.Build(cl));

  ASSERT_TRUE(_m.build_from_smiles("CCCCC"));
  fileconv::FileconvResult result = _config.Process(_m);
  EXPECT_FALSE(result.rejected);
}

TEST_F(FileconvTest, TestLongestPathRejected) {
  int argc = 3;
  const char * argv[] = {"fileconv", "-Y", "maxplen=4"};
  Command_Line cl(argc, const_cast<char**>(argv), "Y:");
  ASSERT_TRUE(_config.Build(cl));

  ASSERT_TRUE(_m.build_from_smiles("CCCCCC"));
  fileconv::FileconvResult result = _config.Process(_m);
  EXPECT_TRUE(result.rejected);
}

// Pattern is established, now write parameterized tests.

struct OptionsSmilesResult {
  // First inputs.
  std::vector<const char*> argv;
  const char * options;
  const char * smiles;
  const char * name;

  // And then the results.
  int expected_rejected;
  // For actions that change the name of the molecule.
  const char * expected_smiles;
  const char * expected_name;
};

class FileConvTestP : public testing::TestWithParam<OptionsSmilesResult> {
  protected:
    Molecule _m;
    fileconv::FileconvConfig _config;
};

TEST_P(FileConvTestP, FileConvTestP) {
  const auto& params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  _m.set_name(params.name);
  const int argc = params.argv.size();
  ASSERT_GT(argc, 0);
  char ** argv = new char*[argc];
  std::unique_ptr<char *[]> free_argv(argv);
  int ndx = 0;
  for (const auto& c : params.argv) {
    argv[ndx] = const_cast<char*>(c);
    ++ndx;
  }
  Command_Line cl(argc, argv, params.options);
  // std::cerr << "Testing initial smiles " << _m.smiles() << ' ' << params.name << '\n';
  ASSERT_TRUE(_config.Build(cl));
  fileconv::FileconvResult result = _config.Process(_m);
  EXPECT_EQ(result.rejected, params.expected_rejected) << " rej " << result.rejected << " expected " << params.expected_rejected;
  EXPECT_EQ(_m.smiles(), params.expected_smiles) << " smiles mismatch " << _m.smiles() << ' ' << params.expected_smiles;
  EXPECT_EQ(_m.name(), params.expected_name) << " name mismatch " << _m.name() << ' ' << params.expected_name;
}
INSTANTIATE_TEST_SUITE_P(FileConvTestP, FileConvTestP, testing::Values(
  OptionsSmilesResult{vector<const char*>{"_", "-c", "4"}, "c:", 
                      "CC", "name", 1, "CC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-c", "4"}, "c:",
                      "CCCC", "name", 0, "CCCC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-C", "4"}, "C:",
                      "CCCCCCC", "name", 1, "CCCCCCC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-C", "4"}, "C:",
                      "CCCC", "name", 0, "CCCC", "name"},

  OptionsSmilesResult{vector<const char*>{"_", "-p", "AMW"}, "p:", 
                      "OC(F)(Cl)C(P)(S)NI", "name", 0, "OC(F)(Cl)C(P)(S)NI", "name AMW = 303.4616"},

  OptionsSmilesResult{vector<const char*>{"_", "-Y", "zpad=4"}, "Y:",
                      "C", "name", 0, "C", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-Y", "zpad=9"}, "Y:",
                      "C", "name", 0, "C", "00000name"},
  OptionsSmilesResult{vector<const char*>{"_", "-Y", "zpad=-2"}, "Y:",
                      "C", "name", 0, "C", "me"},

  OptionsSmilesResult{vector<const char*>{"_", "-Y", "FHrmsqb"}, "Y:",
                      "[1C][1C][1C](N)N[1C](F)(F)F name", "name", 0, "[1CH3][1CH2][1CH](N)N[1C](F)(F)F", "name"},

  OptionsSmilesResult{vector<const char*>{"_", "-r", "1"}, "r:",
                      "C", "name", 1, "C", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-r", "1"}, "r:",
                      "C1CC1", "name", 0, "C1CC1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-r", "S:1"}, "r:",
                      "C1CC1", "name", 0, "C1CC1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-r", "S:2"}, "r:",
                      "C1CC1", "name", 1, "C1CC1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-r", "S:1"}, "r:",
                      "C1CC1.C1CC1", "name", 0, "C1CC1.C1CC1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-R", "S:1"}, "R:",
                      "C1CC1.C1CC1", "name", 1, "C1CC1.C1CC1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-R", "S:1"}, "R:",
                      "C1CC12CC2", "name", 1, "C1CC12CC2", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-R", "S:1", "-R", "spiro"}, "R:",
                      "C1CC12CC2", "name", 0, "C1CC12CC2", "name"},

  OptionsSmilesResult{vector<const char*>{"_", "-r", "aliph:1" }, "r:",
                      "c1ccccc1", "name", 1, "C1=CC=CC=C1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-r", "aliph:1" }, "r:",
                      "C1CC1", "name", 0, "C1CC1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-r", "aliph:2" }, "r:",
                      "C1CC1", "name", 1, "C1CC1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-r", "aliph:2" }, "r:",
                      "C1CC1C1CC1", "name", 0, "C1CC1C1CC1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-r", "arom:1" }, "r:",
                      "C1CC1", "name", 1, "C1CC1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-r", "arom:1" }, "r:",
                      "c1ccccc1", "name", 0, "C1=CC=CC=C1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-r", "arom:2" }, "r:",
                      "c1ccccc1", "name", 1, "C1=CC=CC=C1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-r", "arom:2" }, "r:",
                      "c1ccccc1.c1ccccc1", "name", 0, "C1=CC=CC=C1.C1=CC=CC=C1", "name"},

  OptionsSmilesResult{vector<const char*>{"_", "-R", "arom:0" }, "R:",
                      "C1CC1", "name", 0, "C1CC1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-R", "arom:0" }, "R:",
                      "c1ccccc1C", "name", 1, "C1=CC=CC=C1C", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-R", "arom:1" }, "R:",
                      "c1ccccc1C", "name", 0, "C1=CC=CC=C1C", "name"},

  OptionsSmilesResult{vector<const char*>{"_", "-R", "aliph:0" }, "R:",
                      "CCC", "name", 0, "CCC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-R", "aliph:0" }, "R:",
                      "C1CC1", "name", 1, "C1CC1", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-R", "aliph:1" }, "R:",
                      "C1CC1", "name", 0, "C1CC1", "name"},

  // Some fragment related tests. More are needed.
  OptionsSmilesResult{vector<const char*>{"_", "-f", "rmlarge"}, "f:",
                      "C.CC", "name", 0, "C", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "rmlarge"}, "f:",
                      "CCC.C.CC", "name", 0, "C.CC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "rmlarge=2"}, "f:",
                      "CCC.C.CC", "name", 0, "C", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "rmsmall"}, "f:",
                      "CCC.C.CC", "name", 0, "CCC.CC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "rmsmall=2"}, "f:",
                      "CCC.C.CC", "name", 0, "CCC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "keepsmall"}, "f:",
                      "CCC.C.CC", "name", 0, "C", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "keepsmall=2"}, "f:",
                      "CCC.C.CC", "name", 0, "C.CC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "smarts:[CD0]"}, "f:",
                      "CCC.C.CC", "name", 0, "C", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "smarts:[1C]"}, "f:",
                      "CCC.C.[1CH3]C", "name", 0, "[1CH3]C", "name"},

  OptionsSmilesResult{vector<const char*>{"_", "-f", "RMDUP"}, "f:",
                      "C.C", "name", 0, "C", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "RMDUP"}, "f:",
                      "C.C.C", "name", 0, "C", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "RMDUP"}, "f:",
                      "CC.CCC.C", "name", 0, "CC.CCC.C", "name"},

  OptionsSmilesResult{vector<const char*>{"_", "-f", "WINDOW=,3"}, "f:",
                      "CCC.C.CC.CCCC", "name", 0, "CCC.C.CC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "WINDOW=1,3"}, "f:",
                      "CCC.C.CC.CCCC", "name", 0, "CCC.C.CC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "WINDOW=2,4"}, "f:",
                      "CCC.C.CC.CCCC", "name", 0, "CCC.CC.CCCC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-f", "WINDOW=2,2"}, "f:",
                      "CCC.C.CC.CCCC.CC", "name", 0, "CC.CC", "name"},

  OptionsSmilesResult{vector<const char*>{"_", "-O", "def"}, "O:",
                      "CNOFSClBrI[Li][Na][K][Mg][Ca]", "name", 0, "CNOFSClBrI[Li][Na][K][Mg][Ca]", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-O", "def"}, "O:",
                      "B", "name", 1, "B", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-O", "okOBO"}, "O:",
                      "B(O)(O)B(O)O", "name", 0, "B(O)(O)B(O)O", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-O", "none", "-f", "lod"}, "f:O:",
                      "[Cu].CC", "name", 0, "CC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-O", "def", "-f", "lod"}, "f:O:",
                      "[Li].CC", "name", 0, "CC", "name"},

  OptionsSmilesResult{vector<const char*>{"_", "-h", "SMARTS:[CD1]"}, "h:",
                      "CC1CC1", "name", 0, "C([H])([H])([H])C1CC1", "name"},

  OptionsSmilesResult{vector<const char*>{"_", "-I", "add:3-CC"}, "I:",
                      "[3CH3]C", "name", 0, "[3CH2](C)CC", "name"}
));

}  // namespace
