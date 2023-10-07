#include <ostream>

#include "fetch_via_regexp.h"
#include "iwre2.h"

Fetch_via_Regexp::Fetch_via_Regexp()
{
  _n = 0;
  _rx = nullptr;
  _matches_found = nullptr;

  return;
}

Fetch_via_Regexp::~Fetch_via_Regexp()
{
  if (nullptr != _rx)
    delete [] _rx;

  if (nullptr != _matches_found)
    delete [] _matches_found;

  return;
}

int
Fetch_via_Regexp::build(Command_Line & cl, const char flag)
{
  assert (0 == _n);

  _n = cl.option_count(flag);

  if (0 == _n)
  {
    std::cerr << "Fetch_via_Regexp::build:no option '" << flag << "' specified\n";
    return 0;
  }

  _rx = new std::unique_ptr<re2::RE2>[_n];

  for (int i = 0; i < _n; ++i)
  {
    const_IWSubstring s;
    cl.value(flag, s);

    if (! iwre2::RE2Reset(_rx[i], s)) {
      std::cerr << "Fetch_via_Regexp::build:invalid regular expression '" << s << "'\n";
      return 0;
    }
  }

  _matches_found = new_int(_n);

  return _n;
}

int
Fetch_via_Regexp::matches(const IWString & s)
{
  const re2::StringPiece tmp(s.data(), s.length());
  for (int i = 0; i < _n; ++i)
  {
    if (iwre2::RE2PartialMatch(s, *_rx[i])) {
      _matches_found[i]++;
      return 1;
    }
  }

  return 0;
}

int
Fetch_via_Regexp::report(std::ostream & output) const
{
  output << "Fetch_via_Regexp::report\n";
  for (int i = 0; i < _n; ++i)
  {
    output << " rx   " << _rx[i]->pattern() << ' ' << _matches_found[i] << '\n';
  }

  return 1;
}
