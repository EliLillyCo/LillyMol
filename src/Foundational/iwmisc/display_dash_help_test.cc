#include "Foundational/iwmisc/misc.h"

#include <iostream>

#include "gtest/gtest.h"

namespace {

Command_Line
MakeCommandLine(std::initializer_list<const char*> args) {
  std::vector<char*> argv;
  argv.reserve(args.size());
  for (const char* arg : args) {
    argv.push_back(const_cast<char*>(arg));
  }

  return Command_Line(static_cast<int>(argv.size()), argv.data(), "q:");
}

TEST(DisplayDashHelpIfRequested, IgnoresMissingOption) {
  Command_Line cl = MakeCommandLine({"prog"});

  int calls = 0;
  DisplayDashHelpIfRequested(cl, 'q', [&calls]() { ++calls; });

  EXPECT_EQ(calls, 0);
}

TEST(DisplayDashHelpIfRequested, IgnoresOtherValues) {
  Command_Line cl = MakeCommandLine({"prog", "-q", "SMARTS:c", "-q", "F:query.qry"});

  int calls = 0;
  DisplayDashHelpIfRequested(cl, 'q', [&calls]() { ++calls; });

  EXPECT_EQ(calls, 0);
}

TEST(DisplayDashHelpIfRequested, CallsFunctionAndExitsOnHelp) {
  EXPECT_EXIT(
      {
        Command_Line cl = MakeCommandLine({"prog", "-q", "SMARTS:c", "-q", "help"});
        DisplayDashHelpIfRequested(cl, 'q', []() { std::cerr << "query help\n"; });
        std::cerr << "not reached\n";
      },
      testing::ExitedWithCode(0), "query help");
}

TEST(DisplayDashHelpIfRequested, PassesZeroToFunctionTakingInt) {
  EXPECT_EXIT(
      {
        Command_Line cl = MakeCommandLine({"prog", "-q", "help"});
        DisplayDashHelpIfRequested(cl, 'q', [](int rc) { std::cerr << "rc " << rc << "\n"; });
        std::cerr << "not reached\n";
      },
      testing::ExitedWithCode(0), "rc 0");
}

}  // namespace
