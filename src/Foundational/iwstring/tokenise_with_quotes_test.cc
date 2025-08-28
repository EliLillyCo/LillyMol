
#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

#include "iwstring.h"

namespace {

using iwstring::TokeniseWithQuotes;

struct Data {
  IWString buffer;
  char sep;
  int ntokens;
  std::vector<const char*> expected;
};

class TestTokenise: public testing::TestWithParam<Data> {
  protected:
    resizable_array<int> _tstart;
    resizable_array<int> _tstop;
};

TEST_P(TestTokenise, TestTokenise) {
  const auto params = GetParam();
  EXPECT_EQ(TokeniseWithQuotes(params.buffer, params.sep, _tstart, _tstop), params.ntokens) <<
                params.buffer;

  // Expected failure encountered, cannot extract matching tokens.
  if (params.ntokens < 0) {
    return;
  }

  for (int i = 0; i < params.ntokens; ++i) {
    int b = _tstart[i];
    int e = _tstop[i];
    // std::cerr << "b " << b << " e " << e << '\n';
    const_IWSubstring token(params.buffer.rawdata() + b, e - b);
    EXPECT_EQ(params.expected[i], token) << i << " mismatch '" << params.expected[i] <<
                        "' got '" << token << "' in " << params.buffer;
  }
}
INSTANTIATE_TEST_SUITE_P(TestTokenise, TestTokenise, testing::Values(
  Data{"a,b", ',', 2, {"a", "b"}},
  Data{"aa,b", ',', 2, {"aa", "b"}},
  Data{"aa,bb", ',', 2, {"aa", "bb"}},
  Data{"aaa,bb", ',', 2, {"aaa", "bb"}},
  Data{R"("a","b")", ',', 2, {"a", "b"}},
  Data{R"("a a","b")", ',', 2, {"a a", "b"}},
  Data{R"("a a","b b")", ',', 2, {"a a", "b b"}},
  Data{R"(a,"b b")", ',', 2, {"a", "b b"}},
  Data{R"("a a",b)", ',', 2, {"a a", "b"}},
  Data{R"("a,a",b)", ',', 2, {"a,a", "b"}},
  Data{R"("a,a",,b)", ',', 3, {"a,a", "", "b"}},
  Data{R"("a,a",,b,)", ',', 4, {"a,a", "", "b", ""}},
  Data{R"(,"a,a",,b,)", ',', -1, {"", "a,a", "", "b", ""}}
));


}  // namespace
