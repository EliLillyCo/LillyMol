#include <unistd.h>

#include <initializer_list>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "argv_permutation.h"
#include "cmdline.h"

namespace {
void Check(std::initializer_list<const char*> input,
           std::initializer_list<const char*> expected, const char* spec = "vs:c:") {
  std::vector<std::string> storage(input.begin(), input.end());
  std::vector<char*> argv;
  for (auto& token : storage) {
    argv.push_back(token.data());
  }
  argv.push_back(nullptr);
  const int parse_argc = cmdline_internal::PermuteArguments(storage.size(), argv.data(), spec);
  EXPECT_LE(parse_argc, static_cast<int>(storage.size()));
  int i = 0;
  for (const char* token : expected) {
    EXPECT_STREQ(argv[i++], token);
  }
  EXPECT_EQ(argv[storage.size()], nullptr);
}

TEST(ArgvPermutation, InterleavedOptionsAndFiles) {
  Check({"tool", "a.smi", "-s", "c", "b.smi", "-v"},
        {"tool", "-s", "c", "-v", "a.smi", "b.smi"});
}
TEST(ArgvPermutation, ClusterAttachedValueAndNegativeValue) {
  Check({"tool", "a.smi", "-vsc", "-c", "-2", "-"},
        {"tool", "-vsc", "-c", "-2", "a.smi", "-"});
}
TEST(ArgvPermutation, EndOfOptions) {
  Check({"tool", "a.smi", "-v", "--", "-s", "c"},
        {"tool", "-v", "--", "a.smi", "-s", "c"});
}
TEST(ArgvPermutation, DelimiterAsRequiredValue) {
  Check({"tool", "a.smi", "-s", "--", "-v"},
        {"tool", "-s", "--", "-v", "a.smi"});
}
TEST(ArgvPermutation, MissingValueRemainsMissing) {
  Check({"tool", "a.smi", "-s"}, {"tool", "-s", "a.smi"});
}
TEST(ArgvPermutation, OptionalValueIsAttachedOnly) {
  Check({"tool", "a.smi", "-s", "b.smi", "-sc"},
        {"tool", "-s", "-sc", "a.smi", "b.smi"}, "s::");
}
TEST(ArgvPermutation, ExplicitOrdering) {
  Check({"tool", "a.smi", "-v"}, {"tool", "a.smi", "-v"}, "+v");
}
TEST(CommandLine, RepeatedInterleavedParsing) {
  for (int i = 0; i < 3; ++i) {
    char program[] = "tool", file[] = "file.smi", option[] = "-s", value[] = "c";
    char* argv[] = {program, file, option, value, nullptr};
    Command_Line cl(4, argv, "s:");
    ASSERT_EQ(cl.unrecognised_options_encountered(), 0);
    ASSERT_EQ(cl.option_count('s'), 1);
    EXPECT_STREQ(cl.option_value('s'), "c");
    ASSERT_EQ(cl.number_elements(), 1);
    EXPECT_STREQ(cl[0], "file.smi");
  }
}
}  // namespace

TEST(ArgvPermutation, MissingValueParsingBoundary) {
  char program[] = "tool", file[] = "a.smi", option[] = "-s";
  char* argv[] = {program, file, option, nullptr};
  EXPECT_EQ(cmdline_internal::PermuteArguments(3, argv, "s:"), 2);
  EXPECT_STREQ(argv[2], "a.smi");
}

TEST(ArgvPermutation, StopAtFirstOperandParser) {
  // GNU getopt with '+' reproduces BSD's stop-at-first-operand behavior.
  for (int i = 0; i < 3; ++i) {
    char program[] = "tool", file[] = "a.smi", option[] = "-s", value[] = "c";
    char* argv[] = {program, file, option, value, nullptr};
    const int argc = cmdline_internal::PermuteArguments(4, argv, "s:");
#if defined(__APPLE__)
    optreset = 1;
    optind = 1;
#else
    optind = 0;
#endif
    EXPECT_EQ(getopt(argc, argv, "+s:"), 's');
    EXPECT_STREQ(optarg, "c");
    EXPECT_EQ(getopt(argc, argv, "+s:"), -1);
    EXPECT_EQ(optind, 3);
    EXPECT_STREQ(argv[optind], "a.smi");
  }
}
TEST(ArgvPermutation, MissingValueIsReported) {
  char program[] = "tool", file[] = "a.smi", option[] = "-s";
  char* argv[] = {program, file, option, nullptr};
  const int argc = cmdline_internal::PermuteArguments(3, argv, "s:");
#if defined(__APPLE__)
  optreset = 1;
  optind = 1;
#else
  optind = 0;
#endif
  opterr = 0;
  EXPECT_EQ(getopt(argc, argv, ":s:"), ':');
  EXPECT_EQ(getopt(argc, argv, ":s:"), -1);
#if defined(__APPLE__)
  EXPECT_EQ(optind, 3);
#else
  EXPECT_EQ(optind, 2);
#endif
  // The operand is outside the reduced parsing range. BSD getopt can advance
  // optind beyond that range when reporting a missing value.
  EXPECT_STREQ(argv[argc], "a.smi");
}
