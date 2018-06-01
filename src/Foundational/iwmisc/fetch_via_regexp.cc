#include <ostream>

#include "fetch_via_regexp.h"

Fetch_via_Regexp::Fetch_via_Regexp()
{
  _n = 0;
  _rx = NULL;
  _matches_found = NULL;

  return;
}

Fetch_via_Regexp::~Fetch_via_Regexp()
{
  if (NULL != _rx)
    delete [] _rx;

  if (NULL != _matches_found)
    delete [] _matches_found;

  return;
}

int
Fetch_via_Regexp::build (Command_Line & cl, const char flag)
{
  assert (0 == _n);

  _n = cl.option_count(flag);

  if (0 == _n)
  {
    cerr << "Fetch_via_Regexp::build:no option '" << flag << "' specified\n";
    return 0;
  }

  _rx = new IW_Regular_Expression[_n];

  for (int i = 0; i < _n; ++i)
  {
    const_IWSubstring s;
    cl.value(flag, s);

    if (! _rx[i].set_pattern(s))
    {
      cerr << "Fetch_via_Regexp::build:invalid regular expression '" << s << "'\n";
      return 0;
    }
  }

  _matches_found = new_int(_n);

  return _n;
}

int
Fetch_via_Regexp::matches (const IWString & s)
{
  for (int i = 0; i < _n; ++i)
  {
    if (_rx[i].matches(s))
    {
      _matches_found[i]++;
      return 1;
    }
  }

  return 0;
}

int
Fetch_via_Regexp::report (std::ostream & output) const
{
  output << "Fetch_via_Regexp::report\n";
  for (int i = 0; i < _n; ++i)
  {
    output << " rx   " << _rx[i].source() << ' ' << _matches_found[i] << '\n';
  }

  return 1;
}

