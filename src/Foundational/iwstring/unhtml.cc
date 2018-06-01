#include "iwstring.h"

/*
  change things like
  foo&nbsp;bar
  to
  foo bar

  We translate
  &lt;   <
  &nbsp; ' '
  &gt;   >
  &amp;  &
  &quot; "
  &apos; -
*/

static int
next_instance(const char * s,
              const int n,
              const char iseek)
{
  for (int i = 0; i < n; ++i)
  {
    if (iseek == s[i])
      return i;
  }

  return -1;
}

int
IWString::unhtml()
{
  int ndx = 0;

  int rc = 0;

  for (int i = 0; i < _number_elements; ++i)
  {
    if ('&' != _things[i])
    {
      if (i > ndx)
        _things[ndx] = _things[i];
      ndx++;
      continue;
    }

    int nlook = _number_elements - i;
    if (nlook > 6)
      nlook = 6;

    int semicolon = next_instance(_things + i, nlook, ';');

//  cerr << "Examining " << nlook << " characters starting at " << i << " semicolon " << semicolon << endl;

    if (semicolon < 0)
    {
      if (i > ndx)
        _things[ndx] = _things[i];
      ndx++;
      continue;
    }

    const_IWSubstring s;
    from_to(i + 1, i + semicolon - 1, s);
//  cerr << "Checking s '" << s << "'\n";

    if ("nbsp" == s)
    {
      _things[ndx] = ' ';
    }
    else if ("lt" == s)
    {
      _things[ndx] = '<';
    }
    else if ("gt" == s)
    {
      _things[ndx] = '>';
    }
    else if ("amp" == s)
    {
      _things[ndx] = '&';
    }
    else if ("quot" == s)
    {
      _things[ndx] = '"';
    }
    else if ("apos" == s)
    {
      _things[ndx] = '\'';
    }
    else
    {
      cerr << "IWString::unhtml:unrecognised directive '" << s << "'\n";
      return 0;
    }

    rc++;
    i += semicolon;
    ndx++;
  }

  _number_elements = ndx;

  return rc;
}
