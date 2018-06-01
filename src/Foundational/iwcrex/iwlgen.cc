#include <stdlib.h>

#include <libgen.h>

#include "iwlgen.h"

IW_lgen::IW_lgen ()
{
  _compiled_regular_expression = NULL;
}

IW_lgen::~IW_lgen ()
{
  if (NULL != _compiled_regular_expression)
    delete _compiled_regular_expression;

  return;
}

int
IW_lgen::ok () const
{
  return 1;
}

void
IW_lgen::set_use_extended (int s)
{
  if (s)
    cerr << "IW_lgen::set_use_extended: extended regular expressions not supported\n";

  return;
}

int
IW_lgen::set_pattern (const char * s)
{
  if (NULL != _compiled_regular_expression)
    delete _compiled_regular_expression;

  _compiled_regular_expression = regcmp (s, NULL);
  if (NULL == _compiled_regular_expression)
  {
    cerr << "IW_lgen::IW_lgen: cannot parse '" << s << "'\n";
    return 0;
  }

  return 1;
}

int
IW_lgen::set_pattern (const char * s, int len)
{
  if (NULL != _compiled_regular_expression)
    delete _compiled_regular_expression;

  IWString tmp;
  tmp.resize (len + 1);
  
  tmp.strncpy (s, len);

  return set_pattern (tmp.null_terminated_chars ());
}

int
IW_lgen::set_pattern (const IWString & s)
{
  if (NULL != _compiled_regular_expression)
    delete _compiled_regular_expression;

  IWString tmp;
  tmp.resize (s.length () + 1);

  tmp = s;

  return set_pattern (tmp.null_terminated_chars ());
}

int
IW_lgen::set_pattern (const const_IWSubstring & s)
{
  if (NULL != _compiled_regular_expression)
    delete _compiled_regular_expression;

  IWString tmp;
  tmp.resize (s.length () + 1);

  tmp = s;

  return set_pattern (tmp.null_terminated_chars ());
}

int
IW_lgen::matches (const char * s)
{
  if (NULL == _compiled_regular_expression)
  {
    cerr << "IW_lgen::matches: regular expression not defined\n";
    return 0;
  }

  return (NULL != regex (_compiled_regular_expression, s));
}

int
IW_lgen::matches (const char * s, int len)
{
  if (NULL == _compiled_regular_expression)
  {
    cerr << "IW_lgen::matches: regular expression not defined\n";
    return 0;
  }

  IWString tmp;
  tmp.resize (len + 1);

  tmp.strncpy (s, len);

  return (NULL != regex (_compiled_regular_expression, tmp.null_terminated_chars ()));
}

int
IW_lgen::matches (const const_IWSubstring & s)
{
  if (NULL == _compiled_regular_expression)
  {
    cerr << "IW_lgen::matches: regular expression not defined\n";
    return 0;
  }

  IWString tmp;
  tmp.resize (s.length () + 1);

  tmp = s;

  return (NULL != regex (_compiled_regular_expression, tmp.null_terminated_chars ()));
}

int
IW_lgen::matches (const IWString & s)
{
  if (NULL == _compiled_regular_expression)
  {
    cerr << "IW_lgen::matches: regular expression not defined\n";
    return 0;
  }

  IWString tmp;
  tmp.resize (s.length () + 1);

  tmp = s;

  return (NULL != regex (_compiled_regular_expression, tmp.null_terminated_chars ()));
}

int
IW_lgen::matches_save_subexpressions (const char * s, int len)
{
  cerr << "IW_lgen::matches_save_subexpressions: sorry not implemented\n";
  abort ();

  return 0;
}
int
IW_lgen::matches_save_subexpressions (const char * s)
{
  cerr << "IW_lgen::matches_save_subexpressions: sorry not implemented\n";
  abort ();

  return 0;
}

int
IW_lgen::matches_save_subexpressions (const IWString & s)
{
  cerr << "IW_lgen::matches_save_subexpressions: sorry not implemented\n";
  abort ();

  return 0;
}

int
IW_lgen::matches_save_subexpressions (const const_IWSubstring & s)
{
  cerr << "IW_lgen::matches_save_subexpressions: sorry not implemented\n";
  abort ();

  return 0;
}

int
IW_lgen::number_subexpressions () const 
{
  cerr << "IW_lgen::matches_save_subexpressions: sorry not implemented\n";
  abort ();

  return 0;
}

int
IW_lgen::dollar (int whichdollar, const_IWSubstring & zresult) const
{
  cerr << "IW_lgen::dollar: sorry not implemented\n";
  abort ();

  return 0;
}

const const_IWSubstring
IW_lgen::dollar (int d) const
{
  cerr << "IW_lgen::dollar: sorry not implemented\n";
  abort ();

  return 0;
}

int
IW_lgen::s (IWString & s, const const_IWSubstring & f)
{
  cerr << "IW_lgen::s: sorry not implemented\n";
  abort ();

  return 0;
}

int
IW_lgen::s (const IWString & s, const const_IWSubstring & f, IWString & zresult)
{
  cerr << "IW_lgen::matches_save_subexpressions: sorry not implemented\n";
  abort ();

  return 0;
}

int
IW_lgen::dollar (int whichdollar, int & zresult)
{
  cerr << "IW_lgen::dollar: sorry, not implemented\n";
  abort ();

  return 0;
}

int
IW_lgen::dollar (int whichdollar, float & zresult)
{
  cerr << "IW_lgen::dollar: sorry, not implemented\n";
  abort ();

  return 0;
}

int
IW_lgen::dollar (int whichdollar, double & zresult)
{
  cerr << "IW_lgen::dollar: sorry, not implemented\n";
  abort ();

  return 0;
}
