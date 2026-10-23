#include "cmdline.h"

#include <initializer_list>
#include <string>
#include <vector>

#include "gtest/gtest.h"

namespace {

// Own writable, null-terminated argv storage for as long as Command_Line uses
// its strings. Each initializer element is one argument, including spaces.
class Arguments {
 private:
  std::vector<std::string> _strings;
  std::vector<char*> _argv;

 public:
  Arguments(std::initializer_list<const char*> args)
      : _strings(args.begin(), args.end()) {
    for (std::string& arg : _strings) {
      _argv.push_back(arg.data());
    }
    _argv.push_back(nullptr);
  }

  int
  argc() const {
    return static_cast<int>(_strings.size());
  }

  char**
  argv() {
    return _argv.data();
  }
};

// Keep options before operands so these tests exercise the shared API without
// depending on getopt permutation or POSIXLY_CORRECT on the host platform.
TEST(CommandLine, NoArguments) {
  Arguments args{"tool"};
  Command_Line cl(args.argc(), args.argv(), "vn:");
  EXPECT_TRUE(cl.empty());
  EXPECT_EQ(cl.unrecognised_options_encountered(), 0);
  EXPECT_EQ(cl.option_present('v'), 0);
  EXPECT_EQ(cl.option_count('n'), 0);
  EXPECT_EQ(cl.option_value('n'), nullptr);
  EXPECT_TRUE(cl.string_value('n').empty());
}

TEST(CommandLine, PositionalArgumentsWithoutOptions) {
  Arguments args{"tool", "first.smi", "file with spaces.smi", ""};
  Command_Line cl(args.argc(), args.argv(), "");
  ASSERT_EQ(cl.number_elements(), 3);
  EXPECT_STREQ(cl[0], "first.smi");
  EXPECT_STREQ(cl[1], "file with spaces.smi");
  EXPECT_STREQ(cl[2], "");
  EXPECT_EQ(cl.unrecognised_options_encountered(), 0);
}

TEST(CommandLine, ClusteredAndRepeatedFlags) {
  Arguments args{"tool", "-vvq", "-v", "input.smi"};
  Command_Line cl(args.argc(), args.argv(), "vq");
  EXPECT_NE(cl.option_present('v'), 0);
  EXPECT_NE(cl.option_present('q'), 0);
  EXPECT_EQ(cl.option_count('v'), 3);
  EXPECT_EQ(cl.option_count('q'), 1);
  EXPECT_EQ(cl.option_value('v'), nullptr);
  EXPECT_TRUE(cl.string_value('v').empty());
  int value = 17;
  EXPECT_EQ(cl.value('v', value), 0);
  EXPECT_EQ(value, 17);
  ASSERT_EQ(cl.number_elements(), 1);
  EXPECT_STREQ(cl[0], "input.smi");
}

TEST(CommandLine, SeparateAndAttachedValuesInOccurrenceOrder) {
  Arguments args{"tool", "-s", "first", "-ssecond", "-vs", "third"};
  Command_Line cl(args.argc(), args.argv(), "vs:");
  ASSERT_EQ(cl.option_count('s'), 3);
  EXPECT_STREQ(cl.option_value('s'), "first");
  EXPECT_STREQ(cl.option_value('s', 1), "second");
  EXPECT_STREQ(cl.option_value('s', 2), "third");
  EXPECT_EQ(cl.option_value('s', 3), nullptr);
  EXPECT_EQ(cl.option_value('z'), nullptr);
  EXPECT_NE(cl.option_present('v'), 0);
  EXPECT_TRUE(cl.empty());
}

TEST(CommandLine, ValuesAreNotRetokenized) {
  Arguments args{"tool", "-s", "a b 'c' > output", "-s", ""};
  Command_Line cl(args.argc(), args.argv(), "s:");
  ASSERT_EQ(cl.option_count('s'), 2);
  EXPECT_STREQ(cl.option_value('s'), "a b 'c' > output");
  ASSERT_NE(cl.option_value('s', 1), nullptr);
  EXPECT_STREQ(cl.option_value('s', 1), "");
}

TEST(CommandLine, DashPrefixedRequiredValues) {
  Arguments args{"tool", "-s", "-v", "-n", "-12", "input.smi"};
  Command_Line cl(args.argc(), args.argv(), "vs:n:");
  EXPECT_EQ(cl.option_present('v'), 0);
  EXPECT_STREQ(cl.option_value('s'), "-v");
  int n = 0;
  ASSERT_NE(cl.value('n', n), 0);
  EXPECT_EQ(n, -12);
  EXPECT_EQ(cl.unrecognised_options_encountered(), 0);
}

TEST(CommandLine, EndOfOptionsAndStdin) {
  Arguments args{"tool", "-v", "--", "-s", "-", "input.smi"};
  Command_Line cl(args.argc(), args.argv(), "vs:");
  EXPECT_EQ(cl.option_count('v'), 1);
  EXPECT_EQ(cl.option_present('s'), 0);
  ASSERT_EQ(cl.number_elements(), 3);
  EXPECT_STREQ(cl[0], "-s");
  EXPECT_STREQ(cl[1], "-");
  EXPECT_STREQ(cl[2], "input.smi");
}

TEST(CommandLine, LoneDashIsAnOperand) {
  Arguments args{"tool", "-v", "-"};
  Command_Line cl(args.argc(), args.argv(), "v");
  ASSERT_EQ(cl.number_elements(), 1);
  EXPECT_STREQ(cl[0], "-");
  EXPECT_EQ(cl.unrecognised_options_encountered(), 0);
}

TEST(CommandLine, UnknownOptionsDoNotHideValidOnes) {
  Arguments args{"tool", "-z", "-v", "-s", "known", "input.smi"};
  Command_Line cl(args.argc(), args.argv(), "vs:");
  EXPECT_EQ(cl.unrecognised_options_encountered(), 1);
  EXPECT_EQ(cl.option_present('z'), 0);
  EXPECT_EQ(cl.option_count('v'), 1);
  EXPECT_STREQ(cl.option_value('s'), "known");
  ASSERT_EQ(cl.number_elements(), 1);
  EXPECT_STREQ(cl[0], "input.smi");
}

TEST(CommandLine, MissingRequiredValue) {
  Arguments args{"tool", "-v", "-s"};
  Command_Line cl(args.argc(), args.argv(), "vs:");
  EXPECT_EQ(cl.unrecognised_options_encountered(), 1);
  EXPECT_EQ(cl.option_count('v'), 1);
  EXPECT_EQ(cl.option_count('s'), 0);
  EXPECT_TRUE(cl.empty());
}

TEST(CommandLine, StringAccessors) {
  Arguments args{"tool", "-s", "first value", "-s", "second"};
  Command_Line cl(args.argc(), args.argv(), "s:");
  IWString iw;
  ASSERT_NE(cl.value('s', iw), 0);
  EXPECT_EQ(iw, "first value");
  const_IWSubstring substring;
  ASSERT_NE(cl.value('s', substring, 1), 0);
  EXPECT_EQ(substring, "second");
  char buffer[32];
  ASSERT_NE(cl.value('s', buffer, 1), 0);
  EXPECT_STREQ(buffer, "second");
  EXPECT_EQ(cl.string_value('s', 1), "second");
  EXPECT_TRUE(cl.string_value('s', 2).empty());
  std::string standard;
  ASSERT_NE(cl.value_as_std_string('s', standard), 0);
  EXPECT_EQ(standard, "first value");
  EXPECT_EQ(cl.std_string_value('s', 1), "second");
  EXPECT_TRUE(cl.std_string_value('s', 2).empty());
  ASSERT_NE(cl.value<std::string>('s', standard, 1), 0);
  EXPECT_EQ(standard, "second");
}

TEST(CommandLine, AbsentValuesLeaveDestinationsUnchanged) {
  Arguments args{"tool", "-s", "present"};
  Command_Line cl(args.argc(), args.argv(), "s:n:");
  int number = 42;
  EXPECT_EQ(cl.value('n', number), 0);
  EXPECT_EQ(number, 42);
  IWString iw("unchanged");
  EXPECT_EQ(cl.value('s', iw, 1), 0);
  EXPECT_EQ(iw, "unchanged");
  std::string standard("unchanged");
  EXPECT_EQ(cl.value_as_std_string('s', standard, 1), 0);
  EXPECT_EQ(standard, "unchanged");
}

template <typename T>
class NumericValues : public testing::Test {};

using NumericTypes = testing::Types<int, unsigned int, long, unsigned long, long long,
                                    unsigned long long, float, double>;
TYPED_TEST_SUITE(NumericValues, NumericTypes);

TYPED_TEST(NumericValues, ConvertsRepeatedValuesAndRejectsInvalidText) {
  Arguments args{"tool", "-n", "12", "-n", "34", "-n", "invalid", "-n", "12junk"};
  Command_Line cl(args.argc(), args.argv(), "n:");
  TypeParam value = 0;
  ASSERT_NE(cl.value('n', value), 0);
  EXPECT_EQ(value, static_cast<TypeParam>(12));
  ASSERT_NE(cl.value('n', value, 1), 0);
  EXPECT_EQ(value, static_cast<TypeParam>(34));
  EXPECT_EQ(cl.value('n', value, 2), 0);
  EXPECT_EQ(cl.value('n', value, 3), 0);
  EXPECT_EQ(cl.value('n', value, 4), 0);
  // Conversion failures are distinct from syntactically invalid options.
  EXPECT_EQ(cl.unrecognised_options_encountered(), 0);
}

TYPED_TEST(NumericValues, RejectsEmptyValue) {
  Arguments args{"tool", "-n", ""};
  Command_Line cl(args.argc(), args.argv(), "n:");
  TypeParam value = 17;
  EXPECT_EQ(cl.value('n', value), 0);
  EXPECT_EQ(value, static_cast<TypeParam>(17));
  EXPECT_EQ(cl.unrecognised_options_encountered(), 0);
}

TEST(CommandLine, FloatingPointValues) {
  Arguments args{"tool", "-n", "-1.25", "-n", "2.5e2"};
  Command_Line cl(args.argc(), args.argv(), "n:");
  double d = 0;
  ASSERT_NE(cl.value('n', d), 0);
  EXPECT_DOUBLE_EQ(d, -1.25);
  float f = 0;
  ASSERT_NE(cl.value('n', f, 1), 0);
  EXPECT_FLOAT_EQ(f, 250.0f);
}

TEST(CommandLine, AllValuesAppendInOccurrenceOrder) {
  Arguments args{"tool", "-s", "a", "-s", "b c"};
  Command_Line cl(args.argc(), args.argv(), "s:");
  resizable_array<const char*> values;
  values.add("existing");
  EXPECT_EQ(cl.all_values('s', values), 2);
  ASSERT_EQ(values.number_elements(), 3);
  EXPECT_STREQ(values[0], "existing");
  EXPECT_STREQ(values[1], "a");
  EXPECT_STREQ(values[2], "b c");
  EXPECT_EQ(cl.all_values('z', values), 0);
  EXPECT_EQ(values.number_elements(), 3);
}

TEST(CommandLine, AllStringValuesWithAndWithoutSplitting) {
  Arguments args{"tool", "-s", "a", "-s", "b  c", "-s", ""};
  Command_Line cl(args.argc(), args.argv(), "s:");
  resizable_array_p<IWString> whole;
  EXPECT_EQ(cl.all_values('s', whole), 3);
  ASSERT_EQ(whole.number_elements(), 3);
  EXPECT_EQ(*whole[0], "a");
  EXPECT_EQ(*whole[1], "b  c");
  EXPECT_TRUE(whole[2]->empty());
  resizable_array_p<IWString> split;
  EXPECT_EQ(cl.all_values('s', split, 1), 4);
  ASSERT_EQ(split.number_elements(), 4);
  EXPECT_EQ(*split[0], "a");
  EXPECT_EQ(*split[1], "b");
  EXPECT_EQ(*split[2], "c");
  EXPECT_TRUE(split[3]->empty());
}

TEST(CommandLine, SuccessiveParsersDoNotShareOptionsOrErrors) {
  for (int iteration = 0; iteration < 3; ++iteration) {
    {
      Arguments args{"tool", "-z", "-v", "-s", "first"};
      Command_Line cl(args.argc(), args.argv(), "vs:");
      EXPECT_EQ(cl.unrecognised_options_encountered(), 1);
      EXPECT_STREQ(cl.option_value('s'), "first");
    }
    {
      Arguments args{"other", "-n", "42", "second.smi"};
      Command_Line cl(args.argc(), args.argv(), "n:");
      EXPECT_EQ(cl.unrecognised_options_encountered(), 0);
      EXPECT_EQ(cl.option_present('v'), 0);
      int value = 0;
      ASSERT_NE(cl.value('n', value), 0);
      EXPECT_EQ(value, 42);
      ASSERT_EQ(cl.number_elements(), 1);
      EXPECT_STREQ(cl[0], "second.smi");
    }
  }
}

}  // namespace
