// Tester for the string classes
#include <stdlib.h>

#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

#include "iwstring.h"

namespace {

TEST(TestConstructors, TestAll) {
  IWString foo;
  
  EXPECT_EQ(foo, "");

  EXPECT_EQ("", foo);

  EXPECT_FALSE(foo != "");

  EXPECT_TRUE(" " != foo);

  IWString bar("hello world");
  EXPECT_EQ(bar, "hello world");
  EXPECT_EQ("hello world", bar);

  IWString bg = 'h';
  EXPECT_EQ(bg, "h");
  EXPECT_EQ("h", bg);

  const_IWSubstring s1;
  s1 = "hello world";
  EXPECT_EQ(s1, "hello world");
  EXPECT_EQ("hello world", s1);

  s1 = bar;
  EXPECT_EQ(s1, "hello world");
  EXPECT_EQ("hello world", s1);

  const char * hello = "hello";
  IWString h = "hello";

  EXPECT_EQ(hello, h);
  EXPECT_EQ(h, hello);

  const_IWSubstring hs = h;
  EXPECT_EQ(hello, hs);

  foo = "abc def ghi";

  EXPECT_EQ(foo, "abc def ghi");

  EXPECT_EQ("abc def ghi", foo);

  EXPECT_EQ("def", foo.word(1));

  EXPECT_EQ("ghi", foo.word(2));

  hs = "blargldy glerf";

  EXPECT_EQ("blargldy glerf", hs);
}


TEST(TestIWString, TestAsString) {
  IWString s("hello");
  const std::string as_string = s.AsString();
  EXPECT_EQ(as_string, "hello");
}
TEST(TestConstIWSubstring, TestAsString) {
  const_IWSubstring s("hello");
  const std::string as_string = s.AsString();
  EXPECT_EQ(as_string, "hello");
}
TEST(TestAppendOperator, TestStd) {
  IWString s("hello");
  std::string_view space(" ");
  std::string world("world");
  s << space << world;
  EXPECT_EQ(s, "hello world");
}
TEST(TestEqualsOperator, TestStd) {
  IWString hello = "hello";
  std::string world;
  EXPECT_FALSE(iwstring::Equals(hello, world));
}

TEST(TestIWString, TestExpandEnvironmentVariablesNothing) {
  const IWString hello("hello world");
  EXPECT_EQ(hello.ExpandEnvironmentVariables(), hello);
}
TEST(TestIWString, TestExpandEnvironmentVariablesTooShort) {
  const IWString hello("${}");
  EXPECT_EQ(hello.ExpandEnvironmentVariables(), hello);
}
TEST(TestIWString, TestExpandEnvironmentVariablesNotSet) {
  const IWString hello("${NOTSET}");
  EXPECT_FALSE(hello.ExpandEnvironmentVariables());
}
TEST(TestIWString, TestExpandEnvironmentVariablesEmpty) {
  const IWString hello("xx${}foobar");
  EXPECT_FALSE(hello.ExpandEnvironmentVariables());
}
TEST(TestIWString, TestExpandEnvironmentVariablesNotClosed) {
  const IWString hello("xx${xyfoobar");
  EXPECT_FALSE(hello.ExpandEnvironmentVariables());
}

struct EnvData {
  // Pairs of SHELL_VAR=>value
  std::unordered_map<std::string, std::string> to_set;
  // something like "hello ${world}"
  IWString input;
  // something like "hello world"
  IWString expected;
};

class TestExpandEnv: public testing::TestWithParam<EnvData> {
  protected:
    // Set of shell variables to be unset
    std::unordered_set<std::string> _to_clear;

    void TearDown();
};

void
TestExpandEnv::TearDown() {
  for (const auto& vname : _to_clear) {
    ::unsetenv(vname.c_str());
  }
}

TEST_P(TestExpandEnv, TestExpandEnv) {
  const auto params = GetParam();
  for (const auto& [vname, value] : params.to_set) {
    _to_clear.emplace(vname);
    ::setenv(vname.c_str(), value.c_str(), 1);
  }

  std::optional<IWString> expanded = params.input.ExpandEnvironmentVariables();
  EXPECT_EQ(expanded, params.expected) << " expected " << params.expected <<
        " got '" << expanded->AsString();
}
INSTANTIATE_TEST_SUITE_P(TestExpandEnv, TestExpandEnv, testing::Values(
  EnvData{{}, "hello", "hello"},
  EnvData{{{"world", "world"}}, "hello ${world}", "hello world"},
  EnvData{{{"hello", "hello"}}, "${hello} world", "hello world"},
  EnvData{{{"hello", "welcome"}}, "${hello} world", "welcome world"},
  EnvData{{{"a", "abcdefghi"}}, "${a} world", "abcdefghi world"},
  EnvData{{{"a", "abcdefghi"}}, "xxx ${a} world", "xxx abcdefghi world"},
  EnvData{{{"abcdefghi", "a"}}, "${abcdefghi} world", "a world"},
  EnvData{{{"abcdefghi", "a"}}, "xxx ${abcdefghi} y", "xxx a y"},
  EnvData{{{"mm93", "marc marquez"}, 
           {"vr46", "valentino rossi"}},
           "motogp ${mm93} and ${vr46} greats",
           "motogp marc marquez and valentino rossi greats"},
  EnvData{{{"mm93", "marc marquez"}}, "hello $mm93 motogp", "hello $mm93 motogp"},
  EnvData{{{"mm93", "marquez marquez"}}, "hello $mm93}", "hello $mm93}"}
));

TEST(TestIWString, TestEnsureEndsWithEmpty) {
  IWString s;
  EXPECT_TRUE(s.EnsureEndsWith('a'));
  EXPECT_EQ(s, 'a');
}

template <typename T>
class EnsureEndsWithTest : public testing::Test {
 protected:
  void TestEmpty();
  void TestAlreadyEndsWith();
  void TestMustBeAdded();
};

template <typename T>
void
EnsureEndsWithTest<T>::TestEmpty() {
  IWString s;
  T extra("a");
  EXPECT_EQ(s.EnsureEndsWith(extra), 1);
  EXPECT_EQ(s, "a");
}

template <typename T>
void
EnsureEndsWithTest<T>::TestAlreadyEndsWith() {
  IWString s('a');
  T extra("a");
  EXPECT_EQ(s.EnsureEndsWith(extra), 0);
  EXPECT_EQ(s, 'a');
}

template <typename T>
void
EnsureEndsWithTest<T>::TestMustBeAdded() {
  IWString s('a');
  T extra("b");
  EXPECT_EQ(s.EnsureEndsWith(extra), 1);
  EXPECT_EQ(s, "ab");
}

using MyTypes = ::testing::Types<const char*, const IWString&, const const_IWSubstring&>;
TYPED_TEST_SUITE_P(EnsureEndsWithTest);

TYPED_TEST_P(EnsureEndsWithTest, StartsEmpty) {
  // Inside a test, refer to TypeParam to get the type parameter.
  // TypeParam n = 0;
  // Maybe something could be done to combine the const char* type into the template?
  // std::is_pointer...
  // if (std::is_integral<TypeParam>::value) {
  // }

  // You will need to use `this` explicitly to refer to fixture members.
  this->TestEmpty();
}

TYPED_TEST_P(EnsureEndsWithTest, AlreadyEndsWith) {
  // TypeParam n = 0;

  this->TestAlreadyEndsWith();
}

TYPED_TEST_P(EnsureEndsWithTest, MustBeAdded) { 
  // TypeParam n = 0;

  this->TestMustBeAdded();
}

REGISTER_TYPED_TEST_SUITE_P(EnsureEndsWithTest,
                            StartsEmpty, AlreadyEndsWith, MustBeAdded);

INSTANTIATE_TYPED_TEST_SUITE_P(My, EnsureEndsWithTest, MyTypes);

// typed test suites are too complicated...

TEST(TestIWString, TestStartsWithChar) {
  IWString foo('c');
  EXPECT_TRUE(foo.starts_with('c'));
  EXPECT_FALSE(foo.starts_with('b'));
  EXPECT_FALSE(foo.starts_with("cc"));
}

TEST(TestConst_IWSubstring, TestStartsWithChar) {
  const char* s = "c";
  const_IWSubstring foo(s);
  EXPECT_TRUE(foo.starts_with('c'));
  EXPECT_FALSE(foo.starts_with('b'));
  EXPECT_FALSE(foo.starts_with("cc"));
}

TEST(TestIWString, TestStartsWithCharStar) {
  IWString foo('c');
  EXPECT_TRUE(foo.starts_with("c"));
  EXPECT_TRUE(foo.starts_with("c", 1));
  EXPECT_TRUE(foo.starts_with("cx", 1));
  EXPECT_TRUE(foo.starts_with("cx", 1));
  EXPECT_FALSE(foo.starts_with("cx"));

  foo = "abc";
  EXPECT_TRUE(foo.starts_with("a"));
  EXPECT_TRUE(foo.starts_with("ab"));
  EXPECT_TRUE(foo.starts_with("abc"));
  EXPECT_TRUE(foo.starts_with("abc", 1));
  EXPECT_TRUE(foo.starts_with("abc", 2));
  EXPECT_TRUE(foo.starts_with("abc", 3));
  EXPECT_FALSE(foo.starts_with("abcd", 4));
  EXPECT_FALSE(foo.starts_with("x"));
}

TEST(TestConst_IWSubstring, TestStartsWithCharStar) {
  IWString foo('c');
  EXPECT_TRUE(foo.starts_with("c"));
  EXPECT_TRUE(foo.starts_with("c", 1));
  EXPECT_TRUE(foo.starts_with("cx", 1));
  EXPECT_TRUE(foo.starts_with("cx", 1));
  EXPECT_FALSE(foo.starts_with("x"));

  foo = "abc";
  EXPECT_TRUE(foo.starts_with("a"));
  EXPECT_TRUE(foo.starts_with("ab"));
  EXPECT_TRUE(foo.starts_with("abc"));
  EXPECT_FALSE(foo.starts_with("abcd"));
  EXPECT_FALSE(foo.starts_with("x"));
}

TEST(TestIWString, TestStartsWithIWString) {
  IWString foo('c');
  EXPECT_TRUE(foo.starts_with(foo));
  IWString bar("cx");
  EXPECT_TRUE(bar.starts_with(bar));
  EXPECT_FALSE(foo.starts_with(bar));
  EXPECT_TRUE(bar.starts_with(foo));

  foo = "abc";
  EXPECT_TRUE(foo.starts_with(foo));
  EXPECT_FALSE(foo.starts_with(bar));
  EXPECT_FALSE(bar.starts_with(foo));
  bar = "a";
  EXPECT_TRUE(foo.starts_with(bar));
  EXPECT_FALSE(bar.starts_with(foo));
  bar = "ab";
  EXPECT_TRUE(foo.starts_with(bar));
  EXPECT_FALSE(bar.starts_with(foo));
  bar = "abc";
  EXPECT_TRUE(foo.starts_with(bar));
  EXPECT_TRUE(bar.starts_with(foo));
  bar = "abcd";
  EXPECT_FALSE(foo.starts_with(bar));
  EXPECT_TRUE(bar.starts_with(foo));
}

TEST(TestConst_IWSubstring, TestStartsWithIWString) {
  const_IWSubstring foo("c");
  EXPECT_TRUE(foo.starts_with(foo));
  const_IWSubstring bar("cx");
  EXPECT_TRUE(bar.starts_with(bar));
  EXPECT_FALSE(foo.starts_with(bar));
  EXPECT_TRUE(bar.starts_with(foo));

  foo = "abc";
  EXPECT_TRUE(foo.starts_with(foo));
  EXPECT_FALSE(foo.starts_with(bar));
  EXPECT_FALSE(bar.starts_with(foo));
  bar = "a";
  EXPECT_TRUE(foo.starts_with(bar));
  EXPECT_FALSE(bar.starts_with(foo));
  bar = "ab";
  EXPECT_TRUE(foo.starts_with(bar));
  EXPECT_FALSE(bar.starts_with(foo));
  bar = "abc";
  EXPECT_TRUE(foo.starts_with(bar));
  EXPECT_TRUE(bar.starts_with(foo));
  bar = "abcd";
  EXPECT_FALSE(foo.starts_with(bar));
  EXPECT_TRUE(bar.starts_with(foo));
}

TEST(TextNextWord, TestSpaces) {
  const_IWSubstring hello_world("hello  world foo");
  int i = 0;
  const_IWSubstring token;
  ASSERT_TRUE(hello_world.nextword(token, i));
  EXPECT_EQ(token, "hello");
  ASSERT_TRUE(hello_world.nextword(token, i));
  EXPECT_EQ(token, "world");
  ASSERT_TRUE(hello_world.nextword(token, i));
  EXPECT_EQ(token, "foo");
}

TEST(TestStrcmp, TestAll) {
  IWString foo("hello");
  IWString bar("hello");
  
  EXPECT_EQ(foo.strcmp(bar), 0) << "hello.strcmp(hello)";

  EXPECT_EQ(bar.strncmp(foo, 4), 0) << "hello.strncmp(hello, 4";

  bar = "HELLO";

  EXPECT_EQ(foo.strcasecmp(bar), 0) << "HELLO.strcasecmp(hello)";

  EXPECT_EQ(bar.strcasecmp(foo), 0) << "HELLO.strcasecmp(hello)";
}

TEST(TestStrncat, TestAll) {
  IWString foo("hello");

  foo.strncat(" world", 6);

  EXPECT_EQ(foo, "hello world") << "hello.strncat( world, 6)";
}

TEST(TestTestSplitIntoDirectiveAndValue, TestAll) {
  IWString foo("mm=93");

  const_IWSubstring directive;
  int value;
  ASSERT_TRUE(foo.split_into_directive_and_value(directive, '=', value));

  EXPECT_EQ(directive, "mm");
  EXPECT_EQ(value, 93);
}

TEST(TestAppendNumber, TestAll) {
  IWString foo;
  foo.append_number(0);
  EXPECT_EQ(foo, "0") << "append number 0";

  foo.append_number(-8);

  EXPECT_EQ(foo, "0-8") << "append number -8";

  foo.append_number(123456789);

  EXPECT_EQ(foo, "0-8123456789") << "append number 123456789";

  foo = "";

  foo.append_number(-9);
  EXPECT_EQ(foo, "-9") << "append -9";

  foo.append_number(1000);

  EXPECT_EQ(foo, "-91000") << "append 1000";
}


TEST(TestIndex, TestAll) {
  IWString foo("hello");

  auto i = foo.rindex(' ');
  EXPECT_EQ(i, -1) << foo << " rindex space";

  foo = " hello";

  i = foo.rindex(' ');
  EXPECT_EQ(i, 0) << foo << " rindex";

  const_IWSubstring xx(foo);
  i = xx.rindex(' ');

  EXPECT_EQ(i, 0) << xx << " rindex";

  foo = "hello";
  i = foo.rindex('l');
  EXPECT_EQ(i, 3) << foo << " rindex l";

  i = foo.index('l');

  EXPECT_EQ(i, 2) << foo << " index l";

  i = foo.index('h');
  EXPECT_EQ(i, 0) << foo << " index h";
}


TEST(TestStartsWith, TestAll) {
  IWString foo;
  foo = "abcdef";

  EXPECT_TRUE(foo.starts_with ("abc")) << foo << " starts with abc";

  IWString bar("abcd");
  EXPECT_TRUE(foo.starts_with(bar)) << bar << " starts with abcd";

  const_IWSubstring bb = substr(bar, 0, 2);
  EXPECT_TRUE(foo.starts_with(bb));
}

TEST(TestCharacterEquality, TestAll) {
  IWString foo;

  foo = "f";

  EXPECT_EQ(foo, 'f');
  EXPECT_EQ('f', foo);

  EXPECT_TRUE(foo == 'f');
  EXPECT_TRUE('f' == foo);
  EXPECT_TRUE(foo == foo);

  EXPECT_FALSE(foo != 'f');
  EXPECT_FALSE('f' != foo);
  
  const_IWSubstring bar = "b";
  EXPECT_TRUE(bar == 'b') << bar << " == b";
  EXPECT_TRUE('b' == bar) << " b == " << bar;

  EXPECT_FALSE(bar != 'b');
  EXPECT_FALSE('b' != bar);
}

TEST(TestChop, TestAll) {
  IWString foo ("a;sldkfjas;dlfkj");

  foo.chop();
  EXPECT_EQ(foo, "a;sldkfjas;dlfk");

  foo.chop(2);
  EXPECT_EQ(foo, "a;sldkfjas;dl");
}

TEST(TestPlus, TestAll) {
  IWString foo("abc");
  foo += "def";

  EXPECT_EQ(foo, "abcdef");

  IWString bar = foo + "1";

  EXPECT_EQ(bar, "abcdef1");

  foo = "abc";
  bar = "111";

  IWString s = foo + "3" + bar;

  EXPECT_EQ(s, "abc3111");

  foo='3';
  s = "12" + foo;

  EXPECT_EQ(s, "123");
}

TEST(TestPlusPlus, TestAll) {
  const IWString foo = "hello world";

  const_IWSubstring bar;

  foo.from_to (0, foo.length() - 1, bar);

  EXPECT_EQ(bar, "hello world");
  IWString result;
  char c;
  while ((c = ++bar)) {
    result += c;
  }

  EXPECT_EQ(result, "hello world");

  bar = foo.from_to (0, foo.length() - 1);

  EXPECT_EQ(bar, "hello world");

  result.iwtruncate (0);
  while ((c = bar++)) {
    result += c;
  }

  EXPECT_EQ(result, "ello world");
}

TEST(TestAppend, TestAll) {
  IWString foo;

  foo += "hello";
  EXPECT_EQ(foo, "hello");

  foo += ' ';

  EXPECT_EQ(foo, "hello ");

  foo += "world";
  EXPECT_EQ(foo, "hello world");

  IWString bar = "hello " + foo;
  EXPECT_EQ(bar, "hello hello world");

  foo = " hello";
  const_IWSubstring world (" world");

  foo += world;

  EXPECT_EQ(foo, " hello world");

  foo = "";
  foo.append_with_spacer ("hello");

  EXPECT_EQ(foo, "hello");

  foo.append_with_spacer ("world");

  EXPECT_EQ(foo, "hello world");
}

TEST(TestRemoveToFirst, TestAll) {
  IWString foo("hello world");

  int chars_removed = foo.remove_up_to_first('e');
  EXPECT_EQ(chars_removed, 2);

  EXPECT_EQ(foo, "llo world");

  const_IWSubstring bar(foo);

  chars_removed = bar.remove_up_to_first('z');
  EXPECT_EQ(chars_removed, 0);

  chars_removed = bar.remove_up_to_first(' ');
  EXPECT_EQ(chars_removed, 4);
  EXPECT_EQ(bar, "world");
}

TEST(TestChange, TestAll) {
  IWString foo("hello world");

  foo.change(0, 4, "hello");

  EXPECT_EQ(foo, "hello world");

  foo.change(6, 10, "world");

  EXPECT_EQ(foo, "hello world");

  foo.change(3, 3, "q");

  EXPECT_EQ(foo, "helqo world");

  foo.change(2, 3, "ll");

  EXPECT_EQ(foo, "hello world");

  foo.change(6, 10, "potato");

  EXPECT_EQ(foo, "hello potato");

  foo.change(6, 11, "");

  EXPECT_EQ(foo, "hello ");

  foo.change(5, 5, " world");

  EXPECT_EQ(foo, "hello world");

  foo.change(9, 10, "ld torana");

  EXPECT_EQ(foo, "hello world torana");

  foo = "abcdefghijklmnopqrstuvwxyz0123456789";

  foo.change(5, 10, "x");

  EXPECT_EQ(foo, "abcdexlmnopqrstuvwxyz0123456789");
}

TEST(TestRemoveSuffix, TestAll) {
  IWString foo("hello world");

  ASSERT_FALSE(foo.remove_suffix()) << "Should not be able to remove suffix " << foo;

  foo = "hello.world";

  ASSERT_TRUE(foo.remove_suffix());

  EXPECT_EQ(foo, "hello");

  foo = ".smi";

  EXPECT_TRUE(foo.remove_suffix());

  EXPECT_EQ(foo, "");
}

TEST(TestCompareWighoutCase, TestAll) {
  IWString foo("ABCDEF");

  EXPECT_TRUE(foo.matches_ignore_case('A')) << "No match A";

  EXPECT_TRUE(foo.matches_ignore_case('a'));

  EXPECT_TRUE(foo.matches_ignore_case("ABCDEF"));

  EXPECT_TRUE(foo.matches_ignore_case("Abcdef"));
}

TEST(TestUnHtml, TestAll) {
  IWString foo = "abcdef";
  int t = foo.unhtml();
  EXPECT_EQ(t, 0);

  foo = "&abcdef";
  t = foo.unhtml();
  EXPECT_EQ(t, 0);

  foo = "  &foo;bar";
  t = foo.unhtml();
  EXPECT_EQ(t, 0);

  foo = "&nbsp;hello";
  t = foo.unhtml();
  EXPECT_EQ(foo, " hello");

  foo = "hello&nbsp;world";
  t = foo.unhtml();
  EXPECT_EQ(foo, "hello world");

  foo = "hello world&nbsp;";
  foo.unhtml();
  EXPECT_EQ(foo, "hello world ");

  foo = "&nbsp;&lt;&gt;&quot;&amp;&apos;";
  foo.unhtml();
  EXPECT_EQ(foo, " <>\"&'");
}

TEST(TestStrips, TestAll) {
  IWString foo ("hello world");

  foo.strip_leading_blanks();

  EXPECT_EQ(foo, "hello world");

  foo.strip_trailing_blanks();

  EXPECT_EQ(foo, "hello world");

  foo = " ";
  foo.strip_leading_blanks();
  EXPECT_EQ(foo, "");

  foo = " ";
  foo.strip_trailing_blanks();
  EXPECT_EQ(foo, "");

  foo = "          hello world  ";
  foo.strip_leading_blanks();
  EXPECT_EQ(foo, "hello world  ");

  foo = "   hello   world        ";
  foo.strip_trailing_blanks();
  EXPECT_EQ(foo, "   hello   world");

  const_IWSubstring bar = "hello world";
  bar.strip_leading_blanks();
  EXPECT_EQ(bar, "hello world");

  bar.strip_trailing_blanks();
  EXPECT_EQ(bar, "hello world");

  bar = "   ";
  bar.strip_leading_blanks();
  EXPECT_EQ(bar, "");

  bar = "       ";
  bar.strip_trailing_blanks();
  EXPECT_EQ(bar, "");

  bar = "    hello world   ";
  bar.strip_leading_blanks();
  EXPECT_EQ(bar, "hello world   ");

  bar.strip_trailing_blanks();
  EXPECT_EQ(bar, "hello world");
}


TEST(TestAppendHex, TestAll) {
  unsigned char y[256];
  for (int i = 0; i < 256; ++i) {
    y[i] = i;
  }

  IWString s;
  s.append_hex(y, 256);

  EXPECT_EQ(s, "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f202122232425262728292a2b2c2d2e2f303132333435363738393a3b3c3d3e3f404142434445464748494a4b4c4d4e4f505152535455565758595a5b5c5d5e5f606162636465666768696a6b6c6d6e6f707172737475767778797a7b7c7d7e7f808182838485868788898a8b8c8d8e8f909192939495969798999a9b9c9d9e9fa0a1a2a3a4a5a6a7a8a9aaabacadaeafb0b1b2b3b4b5b6b7b8b9babbbcbdbebfc0c1c2c3c4c5c6c7c8c9cacbcccdcecfd0d1d2d3d4d5d6d7d8d9dadbdcdddedfe0e1e2e3e4e5e6e7e8e9eaebecedeeeff0f1f2f3f4f5f6f7f8f9fafbfcfdfeff");
}

TEST(TestTruncate, TestAll) {
  IWString foo = "Hello world";

  foo.iwtruncate(foo.length());
  EXPECT_EQ(foo, "Hello world");

  foo.iwtruncate(5);
  EXPECT_EQ(foo, "Hello");

  foo << " world   ";

  foo.truncate_at_first(' ');
  EXPECT_EQ(foo, "Hello");

  foo = "abc.def.ghi";

  foo.truncate_at_last('.');
  EXPECT_EQ(foo, "abc.def");
}

TEST(TestInsert, TestAll) {
  IWString foo("helloworld");

  foo.insert(' ', 5);
  EXPECT_EQ(foo, "hello world");

  foo.insert(" everyone", 5);
  EXPECT_EQ(foo, "hello everyone world");

  foo = "world";
  foo.insert("ello ", 0);
  EXPECT_EQ(foo, "ello world");

  foo.insert('h', 0);
  EXPECT_EQ(foo, "hello world");

  IWString bar = "hello";
  foo = " world";

  foo.insert(bar, 0);
  EXPECT_EQ(foo, "hello world");

  bar = " the";
  foo.insert(bar, 5);
  EXPECT_EQ(foo, "hello the world");

  foo.insert(substr(bar, 1, 1), 3);
  EXPECT_EQ(foo, "heltlo the world");
}

TEST(TestErase, TestAll) {
  IWString foo = "hello world";

  foo.erase(0, 1);
  EXPECT_EQ(foo, "llo world");

  foo.erase(0, 0);
  EXPECT_EQ(foo, "lo world");

  foo.erase(1, 1);
  EXPECT_EQ(foo, "l world");

  foo.erase(1, 4);
  EXPECT_EQ(foo, "lld");

  foo.erase(2, 2);
  EXPECT_EQ(foo, "ll");

  foo = "abcdef";

  foo.erase(0);
  EXPECT_EQ(foo, "bcdef");

  foo.erase(1);
  EXPECT_EQ(foo, "bdef");
}


TEST(TestRemoveChars, TestAll) {
  IWString foo = "hello world";

  foo.remove_chars(0, 1);
  EXPECT_EQ(foo, "ello world");

  foo.remove_chars(1, 2);
  EXPECT_EQ(foo, "eo world");

  foo.remove_chars(foo.length() - 1, 1);
  EXPECT_EQ(foo, "eo worl");

  foo = "hello world";
  foo.remove_chars(0, foo.length());
  EXPECT_EQ(foo, "");

  foo = "hello world";

  foo.remove_from_to(0, 10);
  EXPECT_EQ(foo, "");

  foo = "hello world";
  foo.remove_from_to(1, 1);
  EXPECT_EQ(foo, "hllo world");

  foo.remove_from_to(5, 9);
  EXPECT_EQ(foo, "hllo ");
}

TEST(TestMatchesAtPosition, TestAll) {
  IWString foo("hello world");

  EXPECT_TRUE(foo.matches_at_position(0, "hello"));

  EXPECT_TRUE(foo.matches_at_position(1, "e"));

  EXPECT_TRUE(foo.matches_at_position(6, "world"));
}

TEST(TestRemoveLeadingChars, TestAll) {
  IWString foo("abcdef");

  foo.remove_leading_chars(1);

  EXPECT_EQ(foo, "bcdef");

  foo.remove_leading_chars(2);

  EXPECT_EQ(foo, "def");

  foo = "hello     world   to you and you";

  const_IWSubstring bar = foo;

  bar.remove_leading_words(1);

  EXPECT_EQ(bar, "world   to you and you");

  bar.remove_leading_words(2);

  EXPECT_EQ(bar, "you and you");

  foo = "abcdef";

  foo.remove_leading_words(1);

  EXPECT_EQ(foo, "");

  foo = "abcdef";

  foo.shift(3, ' ');

  EXPECT_EQ(foo, "   abcdef") << "shift(3)";

  foo.remove_leading_chars(4, 'q');

  EXPECT_EQ(foo, "bcdefqqqq");

  foo = "a b c d e";

  foo.remove_word(0);

  EXPECT_EQ(foo, "b c d e");

  foo.remove_word(3);

  EXPECT_EQ(foo, "b c d");

  foo.remove_word(1);

  EXPECT_EQ(foo, "b d");
}

TEST(TestSubstring, TestAll) {
  IWString foo = "abcdef";
  EXPECT_EQ(foo.substr (1, 3), "bcd");

  EXPECT_EQ(substr (foo, 0, 1), "a");

  EXPECT_EQ(foo.substr (0, 1), "a");

  EXPECT_EQ(substr (foo, 0, 2), "ab");

  EXPECT_EQ(foo.substr (0, 2), "ab");

  EXPECT_EQ(substr (foo, 2, 2), "cd");

  EXPECT_EQ(foo.substr (2, 2), "cd");

  EXPECT_EQ(foo.substr (1), "bcdef");

  EXPECT_EQ(foo.from_to (2, 3), "cd");

  EXPECT_EQ(substr (foo, 3), "def");

  EXPECT_EQ(substr (substr (foo, 2), 2), "ef");
}

TEST(TestFind, TestAll) {
  IWString foo = "hello world";
  int i = foo.find ("hello");
  EXPECT_EQ(i, 0);

  i = foo.find ("world");
  EXPECT_EQ(i, 6);

  IWString bar = "hello";
  i = foo.find (bar);
  EXPECT_EQ(i, 0);

  bar = "world";
  i = foo.find (bar);
  EXPECT_EQ(i, 6);

  const_IWSubstring s = foo.from_to (0, 4);
  i = foo.find (s);
  EXPECT_EQ(i, 0);

  s = foo.from_to (6, 10);
  i = foo.find (s);
  EXPECT_EQ(i, 6);
}

TEST(TestStrncpy, TestAll) {
  IWString foo;

  foo.strncpy("hello world", 11);

  EXPECT_EQ(foo, "hello world");

  foo.strncpy("hello world", 5);

  EXPECT_EQ(foo, "hello");

  const_IWSubstring bar;
  
  bar = foo.from_to (0, 4);
  EXPECT_EQ(bar, "hello");

  char buffer[9];

  bar.copy_to_char_array (buffer);

#ifdef NONONO
// this is invalid! bar points to a subset of foo, so we cannot assign it this way!!
   unfortunately, there is no autmatic guard against this
  foo = bar;

  EXPECT_EQ(foo, "hello");
#endif
}

TEST(TestStrncmp, TestAll) {
  IWString foo ("hello world");

  int tmp = foo.strncmp("hello world", ::strlen("hello world"));

  EXPECT_EQ(tmp, 0);

  const_IWSubstring bar = foo;

  tmp = bar.strncmp("hello world", ::strlen("hello world"));

  EXPECT_EQ(tmp, 0);
}

TEST(TestLooksLike, TestAll) {
  IWString foo ("he");

  EXPECT_TRUE(foo.looks_like ("hello world", 2));

  EXPECT_FALSE(foo.looks_like ("h", 1));

  EXPECT_TRUE(foo.looks_like ("he", 2));
}


TEST(TestFromTo, TestAll) {
  IWString foo = "0123456789";

  EXPECT_EQ(from_to (foo, 0, 0), "0");

  EXPECT_EQ(foo.from_to (0, 0), "0");

  EXPECT_EQ(from_to (foo, 0, 1), "01");

  EXPECT_EQ(foo.from_to (0, 1), "01");

  EXPECT_EQ(from_to (foo, 2, 8), "2345678");

  EXPECT_EQ(foo.from_to (2, 8), "2345678");

  EXPECT_EQ(foo.from_to (2, "6"), "23456");

  EXPECT_EQ(foo.before ('0'), "");
  EXPECT_EQ(foo.before ('1'), "0");
  EXPECT_EQ(foo.after  ('9'), "");
  EXPECT_EQ(foo.after  ('3'), "456789");

  const_IWSubstring bar;

  foo.from_to (1, 4, bar);
  EXPECT_EQ(bar, "1234");

  bar = "hello world";

  EXPECT_EQ(bar.before (' '), "hello");
  EXPECT_EQ(bar.after (' '),  "world");
}

TEST(TestCompressBlanks, TestAll) {
  IWString foo ("abcdef");
  foo.compress_blanks();
  EXPECT_EQ(foo, "abcdef");

  foo = "abc d e ";
  foo.compress_blanks();
  EXPECT_EQ(foo, "abc d e ");

  foo = " a b c";
  foo.compress_blanks();
  EXPECT_EQ(foo, " a b c");

  foo = "      ";
  foo.compress_blanks();
  EXPECT_EQ(foo, " ");

  foo = " a   b";
  foo.compress_blanks();
  EXPECT_EQ(foo, " a b");

  foo = "       a            b     zyz      ";
  foo.compress_blanks();
  EXPECT_EQ(foo, " a b zyz ");
}



TEST(TestEndsWith, TestAll) {
  IWString foo ("barbleqwerf");

  EXPECT_TRUE(foo.ends_with('f'));

  EXPECT_TRUE(foo.ends_with("f"));

  EXPECT_TRUE(foo.ends_with("qwerf"));

  EXPECT_TRUE(foo.ends_with("barbleqwerf"));

  EXPECT_TRUE(foo.ends_with(foo));

  IWString bar ("werf");
  EXPECT_TRUE(foo.ends_with(bar));
}

TEST(TestPrevword, TestAll) {
  IWString foo = "hello world from foo";
  int i = foo.length() - 1;
  IWString bar;
  IWString result;
  while (foo.prevword (bar, i))
  {
    if (result.length())
      result += ' ';
    result += bar;
  }

  EXPECT_EQ(result, "foo from world hello");

  foo = "  hello     world";
  i = foo.length() - 1;
  result = "";
  while (foo.prevword (bar, i))
  {
    if (result.length())
      result += ' ';
    result += bar;
  }

  EXPECT_EQ(result, "world hello");

  foo = "  hello world   ";
  i = foo.length() - 1;
  result = "";
  while (foo.prevword (bar, i))
  {
    if (result.length())
      result += ' ';
    result += bar;
  }

  EXPECT_EQ(result, "world hello");
}

TEST(TestNextword, TestAll) {
  IWString foo = "hello world from foo";
  int i = 0;
  IWString bar;
  IWString result;
  while (foo.nextword (bar, i))
  {
    if (result.length())
      result += ' ';
    result += bar;
  }

  EXPECT_EQ(foo, result.null_terminated_chars());

  const_IWSubstring sbar;
  result.iwtruncate(0);

  i = 0;
  while (foo.nextword (sbar, i))
  {
    if (result.length())
      result += ' ';
    result += sbar;
  }

  EXPECT_EQ(result, foo.null_terminated_chars());

  foo = "  hello     world";
  i = 0;
  result = "";
  while (foo.nextword (bar, i))
  {
    if (result.length())
      result += ' ';
    result += bar;
  }

  EXPECT_EQ(result, "hello world");

  foo = "  hello world   ";
  i = 0;
  result = "";
  while (foo.nextword (bar, i))
  {
    if (result.length())
      result += ' ';
    result += bar;
  }

  EXPECT_EQ(result, "hello world");

  foo = "d1";

  IWString d1;
  i = 0;
  foo.nextword (d1, i);

  EXPECT_EQ(d1, "d1");

  sbar = "d1";
  i = 0;
  sbar.nextword (d1, i);

  EXPECT_EQ(d1, "d1");

  foo = "hello world ";
  i = 0;
  foo.nextword (d1, i);
  foo.nextword (d1, i);

  EXPECT_EQ(d1, "world");
}

TEST(TestNextwordSingleDelimiter, TestAll) {
  IWString foo = "hello world from foo";
  int i = 0;
  const_IWSubstring token;
  IWString tmp;

  while (foo.nextword_single_delimiter(token, i, ' '))
  {
    tmp.append_with_spacer(token);
  }

  EXPECT_EQ(tmp, foo);

  foo = "hello  world";

  i = 0;
  foo.nextword_single_delimiter(token, i, ' ');

  EXPECT_EQ(token, "hello");

  foo.nextword_single_delimiter(token, i, ' ');

  EXPECT_EQ(token, "");

  foo.nextword_single_delimiter(token, i, ' ');

  EXPECT_EQ(token, "world");

  const_IWSubstring glerf ("  abc   def  ");
  i = 0;

  glerf.nextword_single_delimiter(token, i, ' ');
  EXPECT_EQ(token, "");

  glerf.nextword_single_delimiter(token, i, ' ');
  EXPECT_EQ(token, "");

  glerf.nextword_single_delimiter(token, i, ' ');
  EXPECT_EQ(token, "abc");

  glerf.nextword_single_delimiter(token, i, ' ');
  EXPECT_EQ(token, "");

  glerf.nextword_single_delimiter(token, i, ' ');
  EXPECT_EQ(token, "");

  glerf.nextword_single_delimiter(token, i, ' ');
  EXPECT_EQ(token, "def");

  glerf.nextword_single_delimiter(token, i, ' ');
  EXPECT_EQ(token, "");

  glerf.nextword_single_delimiter(token, i, ' ');
  EXPECT_EQ(token, "");

  glerf = "a,b,c,d,,,e,,g";

  i = 0;
  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "a");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "b");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "c");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "d");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "e");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "g");

  glerf = "a,b,,,,";
  i = 0;

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "a");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "b");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "");

  glerf.nextword_single_delimiter(token, i, ',');
  EXPECT_EQ(token, "");
}

TEST(TestNext, TestAll) {
  IWString foo = "01234012340123401234";

  const std::vector<int> expected{0, 5, 10, 15};

  std::vector<int> found;

  int j = 0;
  while (foo.next('0', j)) {
    found.push_back(j);
    ++j;
  }

  EXPECT_EQ(found, expected);
}

TEST(TestBaseName, TestAll) {
  IWString foo("foo");

  const_IWSubstring b;

  foo.iwbasename(b);

  EXPECT_EQ(b, "foo");

  foo = "foo/";

  foo.iwbasename(b);

  EXPECT_EQ(b, "foo");

  foo = "glerf//////";

  foo.iwbasename(b);

  EXPECT_EQ(b, "glerf");

  foo = "   ";

  foo.iwbasename(b);

  EXPECT_EQ(b, "   ");

  foo = "abc/def";

  foo.iwbasename(b);

  EXPECT_EQ(b, "def");

  const_IWSubstring bar("a/b/c/d/e/f");

  bar.iwbasename(b);

  EXPECT_EQ(b, "f");
}

TEST(TestStdString, TestAll) {
  std::string foo = "hello world";

  IWString s1(foo);

  EXPECT_EQ(s1, "hello world");

  const_IWSubstring s2(foo);

  EXPECT_EQ(s2, "hello world");

  foo = " hunga dunga";

  s1 = foo;

  EXPECT_EQ(s1, " hunga dunga");

  s2 = foo;

  EXPECT_EQ(s2, " hunga dunga");

  foo = " extra";

  s1 += foo;

  EXPECT_EQ(s1, " hunga dunga extra");

  s1 = "hello world";

  foo = s1;

  EXPECT_EQ(foo, "hello world");

  s2 = " &#^@ _";
}

TEST(Test_GsubString, TestAll) {
  IWString foo ("hello world");

  int i = foo.gsub (" ", "_");
  EXPECT_EQ(i, 1);

  foo.gsub ('_', ' ');

  i = foo.gsub("h", "H");
  EXPECT_EQ(i, 1);

  i = foo.gsub ("d", "D");
  EXPECT_EQ(i, 1);

  i = foo.gsub (" ", "");
  EXPECT_EQ(i, 1);

  foo = "hello world";

  i = foo.gsub ("o w", "O W");
  EXPECT_EQ(i, 1);

  foo = "hello world";

  i = foo.gsub ("hello", "HELLO");
  EXPECT_EQ(i, 1);

  i = foo.gsub ("world", "WORLD");
  EXPECT_EQ(i, 1);

  foo = "hello world";

  i = foo.gsub ("hello", "gday");
  EXPECT_EQ(i, 1);

  i = foo.gsub ("world", "mate");
  EXPECT_EQ(i, 1);

  i = foo.gsub ("day mat", "e");
  EXPECT_EQ(i, 1);

  foo = "hello world";
  i = foo.gsub ("hello world", "gday mate");
  EXPECT_EQ(i, 1);

  foo = "gday mate";

  i = foo.gsub ("gday", "hello");
  EXPECT_EQ(i, 1);

  i = foo.gsub ("mate", "world");
  EXPECT_EQ(i, 1);

  foo = "hello hello";

  i = foo.gsub ("hello", "gdbye");
  EXPECT_EQ(i, 2);

  foo = "hello hello hello hello";

  i = foo.gsub ("hello", "xx");
  EXPECT_EQ(i, 4);

  i = foo.gsub ("xx", "yyy");
  EXPECT_EQ(i, 4);

  foo = "";

  int n = 100;

  for (int i = 0; i < n; i++) {
    foo << " x";
  }

  i = foo.gsub ("x", "yz");
  EXPECT_EQ(foo.length(), 3*n);

  for (int i = 0; i < 3 * n; i++) {
    int j = i % 3;
    if (0 == j && ' ' == foo[i])
      ;
    else if (1 == j && 'y' == foo[i])
      ;
    else if (2 == j && 'z' == foo[i])
      ;
    else {
      EXPECT_TRUE(false) << "Improper character at i = " << i << " got '" << foo[i] << "'\n";
    }
  }
}

TEST(TestGsub, TestAll) {
  IWString foo = "hello world";
  int i = foo.gsub (' ', '_');
  EXPECT_EQ(i, 1);

  i = foo.gsub ('o', ' ', 1);
  EXPECT_EQ(1, i);

  i = foo.gsub ('o', ' ');
  EXPECT_EQ(1, i);

  i = foo.gsub (' ', 'o');
  EXPECT_EQ(2, i);
}

TEST(TestIsInt, TestAll) {
  IWString foo ("err");
  
  int tmp;
  EXPECT_EQ(foo.is_int (tmp), 0);

  foo = "";
  EXPECT_EQ(foo.is_int (tmp), 0);

  foo = "      ";
  EXPECT_EQ(foo.is_int (tmp), 0);

  foo = "+1";
  EXPECT_EQ(foo.is_int (tmp), 1);
  EXPECT_EQ(tmp, 1);

  foo = "234";
  EXPECT_EQ(foo.is_int (tmp), 1);
  EXPECT_EQ(tmp, 234);

  foo = "-009";
  EXPECT_EQ(foo.is_int (tmp), 1);
  EXPECT_EQ(tmp, -9);

  foo = "a1998";
  EXPECT_EQ(is_int (substr (foo, 1), tmp), 1);
  EXPECT_EQ(tmp, 1998);

  EXPECT_EQ(foo.substr (1).is_int (tmp), 1);
  EXPECT_EQ(tmp, 1998);

  foo = "00987";
  EXPECT_EQ(foo.is_int (tmp), 1);
  EXPECT_EQ(tmp, 987);

  foo = "0x1";
  unsigned int utmp;
  EXPECT_EQ(foo.is_hex (utmp), 1);
  EXPECT_EQ(int (utmp), 1);

  foo = "0x1234123";
  int expected_result = 0x1234123;
  EXPECT_EQ(foo.is_hex (utmp), 1);
  EXPECT_EQ(int (utmp), expected_result);
}


static IWString
produce_a_string (const char * s1,
                  const char * s2)
{
  IWString rc(s1);
  rc << s2;

  return rc;
}

TEST(TestReturnedString, TestAll) {
  IWString foo = produce_a_string("abc", "def");

  EXPECT_EQ(foo, "abcdef");

  foo << "1";

  EXPECT_EQ(foo, "abcdef1");
}

TEST(TestPerlSplit, TestAll) {
  IWString foo = "abc d e fghijkl m n OPQRST UVW XYZ";
  resizable_array_p<const_IWSubstring> tokens;

  int ntokens = foo.split(tokens);

  assert (ntokens == tokens.number_elements());

  EXPECT_EQ(*(tokens[0]), "abc");
  EXPECT_EQ(*(tokens[1]), "d");
  EXPECT_EQ(*(tokens[ntokens - 1]), "XYZ");

  foo = "abc/def/ ";
  ntokens = foo.split(tokens, '/');

  ASSERT_EQ(ntokens, tokens.number_elements());
  EXPECT_EQ(ntokens, 3);

  EXPECT_EQ(*(tokens[0]), "abc");
  EXPECT_EQ(*(tokens[1]), "def");
  EXPECT_EQ(*(tokens[2]), " ");

// Since the iwaray<> template cannot be resized, we do those tests within braces

  if (1)
  {
    iwaray<const_IWSubstring> t2;

    foo = "hello world   testing";
  
    ntokens = foo.split(t2, ' ');

    EXPECT_EQ(ntokens, 3);

    EXPECT_EQ(t2[0], "hello");
    EXPECT_EQ(t2[1], "world");
    EXPECT_EQ(t2[2], "testing");
  }

  if (1)
  {
    iwaray<const_IWSubstring> t2;

    foo = "w0;1;hello world;3;444;glerf;last one";

    ntokens = foo.split(t2, ';');

    EXPECT_EQ(ntokens, 7);

    EXPECT_EQ(t2[0], "w0");
    EXPECT_EQ(t2[1], "1");
    EXPECT_EQ(t2[2], "hello world");
    EXPECT_EQ(t2[3], "3");
    EXPECT_EQ(t2[4], "444");
    EXPECT_EQ(t2[5], "glerf");
    EXPECT_EQ(t2[6], "last one");
  }

  if (1)
  {
    foo = " hello     world     ";
    iwaray<IWString> tokens;

    int ntokens = foo.split(tokens, ' ');

    EXPECT_EQ(ntokens, 2);
    EXPECT_EQ(tokens[0], "hello");
    EXPECT_EQ(tokens[1], "world");
  }
}

TEST(TestOperators, TestAll) {
  IWString foo;

  foo << "Hello " << 3 << '5';

  EXPECT_EQ(foo, "Hello 35");

  foo = "Hello ";
  foo << 8 << " world";

  EXPECT_EQ(foo, "Hello 8 world");

  IWString bar = foo;

  EXPECT_EQ(foo, bar);

  EXPECT_FALSE(foo != bar);

  const_IWSubstring g = "hello world";
  const_IWSubstring h = g;

  EXPECT_EQ(g, h);

  EXPECT_FALSE(g != h);

#ifdef DOES_NOT_WORK_MEMORY_PROBLEM
  g = 'q';

  if (g == 'q')
    ;
  else
  {
    cerr << "const_IWSubstring::operator != (char) failed, line " << __LINE__ << "\n";
    return 0;
  }
#endif

  const char * rhs = nullptr;
  IWString nn = rhs;

  EXPECT_TRUE(nn.empty()) << " Null assignment";
}

TEST(TestNWords, TestAll) {
  IWString foo = "";
  int i = foo.nwords();
  EXPECT_EQ(i, 0);

  foo = " ";
  i = foo.nwords();
  EXPECT_EQ(i, 0);

  i = foo.nwords ('"');
  EXPECT_EQ(i, 1);

  foo = "ab";
  i = foo.nwords();
  EXPECT_EQ(i, 1);

  foo = " ab";
  i = foo.nwords();
  EXPECT_EQ(i, 1);

  foo = "ab ";
  i = foo.nwords();
  EXPECT_EQ(i, 1);

  foo = " a  b  c ";
  i = foo.nwords();
  EXPECT_EQ(i, 3);
}

TEST(TestOperatorAngleAngle, TestAll) {
  IWString foo ("hello");
  IWString bar ("world");

  foo << " " << bar;

  EXPECT_EQ(foo, "hello world");

  const_IWSubstring orl = bar.from_to (1, 3);

  foo << " " << orl;

  EXPECT_EQ(foo, "hello world orl");
}
TEST(TestAppendNumberWidth, TestAll) {
  IWString foo;
  append_number(foo, 1, 1);

  EXPECT_EQ(foo, "1");

  append_number(foo, 1, 2);

  EXPECT_EQ(foo, "1 1");

  append_number(foo, 876, 1);

  EXPECT_EQ(foo, "1 1876");

  append_number(foo, -312, 1);

  EXPECT_EQ(foo, "1 1876-312");

  foo.resize(0);

  append_number(foo, 1, 4);

  EXPECT_EQ(foo, "   1");

  append_number(foo, -14, 4);

  EXPECT_EQ(foo, "   1 -14");
}

TEST(TestSplit, TestAll) {
  IWString foo = "hello/world";
  const_IWSubstring before, after;

  int rc = foo.split(before, ' ', after);
  EXPECT_EQ(rc, 0);

  rc = foo.split(before, '/', after);
  EXPECT_EQ(rc, 2);

  EXPECT_EQ(before, "hello");
  EXPECT_EQ(after, "world");

  foo = "/blerfl";

  rc = foo.split(before, '/', after);
  EXPECT_EQ(rc, 1);
  EXPECT_EQ(before, "");
  EXPECT_EQ(after, "blerfl");

  foo = "gbwerdy?";
  rc = foo.split(before, '?', after);
  EXPECT_EQ(rc, 1);
  EXPECT_EQ(before, "gbwerdy");
  EXPECT_EQ(after, "");
}

TEST(TestWord, TestAll) {
  IWString foo = "abc def";
  EXPECT_EQ (foo.word (0), "abc");

  EXPECT_EQ (foo.word (1), "def");

  EXPECT_EQ (foo.word (2), "");

  IWString bar;

  foo = "     abc     def    ";

  foo.word (0, bar);
  EXPECT_EQ (bar, "abc");

  foo.word (1, bar);
  EXPECT_EQ (bar, "def");

  foo.word (2, bar);
  EXPECT_EQ (bar, "");

  foo.word(-1, bar);
  EXPECT_EQ(bar, "def");

  foo.word(-2, bar);
  EXPECT_EQ(bar, "abc");

  bar = foo.word(-1);
  EXPECT_EQ(bar, "def");
  bar = foo.word(-2);
  EXPECT_EQ(bar, "abc");

  bar = "hello ";

  bar.word (0, foo);
  EXPECT_EQ (foo, "hello");

  const_IWSubstring glerf = "	yabba	dabba do";
  IWString w;
  glerf.whitespace_delimited_word (0, w);

  EXPECT_EQ (w, "yabba");

  glerf.whitespace_delimited_word (1, w);
  EXPECT_EQ (w, "dabba");

  glerf.whitespace_delimited_word (2, w);
  EXPECT_EQ (w, "do");
}

TEST(TestKMG, TestAll) {
  IWString foo = "42";
  uint64_t i;
  ASSERT_TRUE(foo.NumericValueKMG(i));
  EXPECT_EQ(i, 42);

  foo = "42k";
  ASSERT_TRUE(foo.NumericValueKMG(i));
  EXPECT_EQ(i, 42000);

  foo = "42M";
  ASSERT_TRUE(foo.NumericValueKMG(i));
  EXPECT_EQ(i, 42000000);

  foo = "42G";
  ASSERT_TRUE(foo.NumericValueKMG(i));
  EXPECT_EQ(i, 42000000000);
}


}  // namespace
