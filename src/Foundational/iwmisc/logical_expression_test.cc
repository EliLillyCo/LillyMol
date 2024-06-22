// Tester for logical expression

#include <vector>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/iwstring.h"

#include "logical_expression.h"

namespace {
struct LogExpData {
  IWString representation;
  std::vector<int> result;
  int done_at;
  int expected_result;
};

class TestLogicalExpression: public testing::TestWithParam<LogExpData> {
  protected:
    IW_Logical_Expression _logexp;
    int _result;
};

IWString
AsString(const std::vector<int>& v) {
  IWString result;
  result << '[';
  for (uint32_t i = 0; i < v.size(); ++i) {
    if (i > 0) {
      result << ',';
    }
    result << v[i];
  }
  result << ']';

  return result;
}

TEST_P(TestLogicalExpression, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(_logexp.BuildFromString(params.representation));
  const int n = params.result.size();
  for (int i = 0; i < n; ++i) {
    _logexp.set_result(i, params.result[i]);
    if (i == params.done_at) {
      EXPECT_TRUE(_logexp.evaluate(_result));
      EXPECT_EQ(_result, params.expected_result) << params.representation << ' ' <<
                AsString(params.result);
      return;
    } else {
      EXPECT_FALSE(_logexp.evaluate(_result)) << "Got result at " << i << " expected " <<
                params.done_at << ' ' << params.representation << ' ' << AsString(params.result);
    }
  }
}
INSTANTIATE_TEST_SUITE_P(TestLogicalExpression, TestLogicalExpression, testing::Values(
  LogExpData{".", {1}, 0, 1},
  LogExpData{".", {0}, 0, 0},
  LogExpData{"!.", {1}, 0, 0},
  LogExpData{"!.", {0}, 0, 1},

  LogExpData{".&.", {0, 0}, 1, 0},
  LogExpData{".&.", {1, 0}, 1, 0},
  LogExpData{".&.", {0, 1}, 1, 0},
  LogExpData{".&.", {1, 1}, 1, 1},

  LogExpData{".,.", {0, 0}, 1, 0},
  LogExpData{".,.", {1, 0}, 0, 1},
  LogExpData{".,.", {0, 1}, 1, 1},
  LogExpData{".,.", {1, 1}, 0, 1},

  LogExpData{".^.", {0, 0}, 1, 0},
  LogExpData{".^.", {1, 0}, 1, 1},
  LogExpData{".^.", {0, 1}, 1, 1},
  LogExpData{".^.", {1, 1}, 1, 0},

  LogExpData{".;.", {0, 0}, 1, 0},
  LogExpData{".;.", {1, 0}, 1, 0},
  LogExpData{".;.", {0, 1}, 1, 0},
  LogExpData{".;.", {1, 1}, 1, 1},

  LogExpData{".,.;.", {0, 0, 0}, 1, 0},
  LogExpData{".,.;.", {1, 0, 0}, 2, 0},
  LogExpData{".,.;.", {1, 0, 1}, 2, 1},
  LogExpData{".,.;.", {0, 1, 0}, 2, 0},
  LogExpData{".,.;.", {0, 1, 1}, 2, 1},

  LogExpData{".;.,.", {0, 0, 0}, 0, 0},
  LogExpData{".;.,.", {0, 1, 0}, 0, 0},
  LogExpData{".;.,.", {0, 1, 1}, 0, 0},
  LogExpData{".;.,.", {1, 0, 0}, 2, 0},
  LogExpData{".;.,.", {1, 1, 0}, 1, 1},
  LogExpData{".;.,.", {1, 0, 1}, 2, 1},
  LogExpData{".;.,.", {1, 1, 1}, 1, 1},

  LogExpData{".;.;.", {0, 0, 0}, 0, 0},
  LogExpData{".;.;.", {1, 0, 0}, 1, 0},
  LogExpData{".;.;.", {1, 1, 0}, 2, 0},
  LogExpData{".;.;.", {1, 1, 1}, 2, 1},

  LogExpData{".&.;.", {0, 0, 0}, 0, 0},
  LogExpData{".&.;.", {1, 0, 0}, 1, 0},
  LogExpData{".&.;.", {1, 1, 0}, 2, 0},
  LogExpData{".&.;.", {1, 1, 1}, 2, 1},

  LogExpData{".;.&.", {0, 0, 0}, 0, 0},
  LogExpData{".;.&.", {1, 0, 0}, 1, 0},
  LogExpData{".;.&.", {1, 1, 0}, 2, 0},
  LogExpData{".;.&.", {1, 1, 1}, 2, 1}
));

}  // namespace
