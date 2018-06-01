#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <regex.h>

#include "iwlregex.h"

void
IW_lregex::_default_values ()
{
  _regcomp_flags = REG_EXTENDED;
  _reg_eflags = 0;

  _regmatch = NULL;
  _regmatch_valid = 0;

  _compiled = 0;

  return;
}

IW_lregex::IW_lregex ()
{
  _default_values ();
}

IW_lregex::IW_lregex (const char * s)
{
  _default_values ();

  set_pattern (s, strlen (s));
  
  return;
}

IW_lregex::IW_lregex (const const_IWSubstring & s)
{
  _default_values ();

  set_pattern (s.rawchars (), s.nchars ());
  
  return;
}

IW_lregex::IW_lregex (const IWString & s)
{
  _default_values ();

  set_pattern (s.rawchars (), s.nchars ());
  
  return;
}

IW_lregex::~IW_lregex ()
{
  if (NULL != _regmatch)
    delete [] _regmatch;

  if (_compiled)
    regfree (& _compiled_regular_expression);

  return;
}

int
IW_lregex::ok () const
{
  if (0 == _compiled_regular_expression.re_nsub && NULL != _regmatch)
  {
    cerr << "IW_lregex::ok: re_nsub/_regmatch mismatch\n";
    return 0;
  }

  if (_compiled_regular_expression.re_nsub && NULL == _regmatch)
  {
    cerr << "IW_lregex::ok: _regmatch/re_nsub mismatch\n";
    return 0;
  }

  return 1;
}

void
IW_lregex::set_use_extended (int s)
{
  if (s)
    _regcomp_flags |= REG_EXTENDED;
  else
    _regcomp_flags &= (~REG_EXTENDED);

  return;
}

int
IW_lregex::set_pattern (const char * s)
{
  if (_compiled)
  {
    regfree (& _compiled_regular_expression);
    _compiled = 0;
  }

  if (NULL != _regmatch)
  {
    delete [] _regmatch;
    _regmatch = NULL;
  }

  int rc = regcomp (&_compiled_regular_expression, s, _regcomp_flags);
  if (0 != rc)
  {
    cerr << "IW_lregex::set_pattern: cannot compile '" << s << "'\n";
    return 0;
  }

  _compiled = 1;

#define DEBUG_SET_PATTERN
#ifdef DEBUG_SET_PATTERN
  cerr << "Pattern '" << s << "', re_nsub = " << _compiled_regular_expression.re_nsub << endl;
  cerr << "Flags " << _regcomp_flags << endl;
#endif

  if (0 == _compiled_regular_expression.re_nsub)    // no fields, we are done
    return 1;

  _regmatch = new regmatch_t[_compiled_regular_expression.re_nsub + 1];

  assert (NULL != _regmatch);

  return 1;
}

int
IW_lregex::set_pattern (const IWString & s)
{
  return set_pattern (s.rawchars (), s.length ());
}

int
IW_lregex::set_pattern (const const_IWSubstring & s)
{
  return set_pattern (s.rawchars (), s.length ());
}

int
IW_lregex::set_pattern (const char * s, int len)
{
  IWString tmp;
  
  tmp.resize (len + 1);
  
  tmp.strncpy (s, len);

  return set_pattern (tmp.null_terminated_chars ());
}

int
IW_lregex::matches (const IWString & s)
{
  return matches (s.rawchars (), s.length ());
}

int
IW_lregex::matches (const const_IWSubstring & s)
{
  return matches (s.rawchars (), s.length ());
}

int
IW_lregex::matches (const char * s, int len)
{
  if (! _compiled)
    return 0;

  IWString tmp;
  tmp.resize (len + 1);

  tmp.strncpy (s, len);

  return matches (tmp.null_terminated_chars ());
}

int
IW_lregex::matches (const char * s)
{
  assert (ok ());

  if (! _compiled)    // no expression compiled
    return 0;

  _regmatch_valid = 0;    // we are not saving subexpressions

  if (0 == regexec (&_compiled_regular_expression, s, 
                     0, NULL,
                     _reg_eflags))
  {
    return 1;    // match
  }

  return 0;
}

int
IW_lregex::matches_save_subexpressions (const char * s, int len)
{
  IWString tmp;
  tmp.resize (len + 1);

  tmp.strncpy (s, len);

  return _matches_save_subexpressions (tmp.null_terminated_chars (), len);
}

int
IW_lregex::matches_save_subexpressions (const const_IWSubstring & s)
{
  return matches_save_subexpressions (s.rawchars (), s.nchars ());
}

int
IW_lregex::matches_save_subexpressions (const IWString & s)
{
  return matches_save_subexpressions (s.rawchars (), s.nchars ());
}

int
IW_lregex::matches_save_subexpressions (const char * s)
{
  return _matches_save_subexpressions (s, -1);    // length unknown right now
}

int
IW_lregex::_matches_save_subexpressions (const char * s, int len)
{
  assert (ok ());

  if (! _compiled)    // no expression compiled
    return 0;

  if (0 == regexec (&_compiled_regular_expression, s,
                     _compiled_regular_expression.re_nsub + 1, _regmatch,
                     _reg_eflags))
  {
    _regmatch_valid = 1;

    _string_which_was_matched.set (s, len);

#ifdef DEBUG_MATCHES_SAVE_SUBEXPRESSIONS
    for (int i = 0; i <= _compiled_regular_expression.re_nsub; i++)
    {
      cerr << "regmatch[" << i << "] from " << _regmatch[i].rm_so << " to " << _regmatch[i].rm_eo << endl;
    }
#endif

    return 1;    // match
  }

  _regmatch_valid = 0;

  return 0;
}

/*
  After a match, the caller wants to know the components which matched
*/

int
IW_lregex::dollar (int whichdollar, const_IWSubstring & zresult) const
{
  assert (whichdollar >= 0 && whichdollar <= int (_compiled_regular_expression.re_nsub));
  assert (NULL != _regmatch);
  assert ("" != _string_which_was_matched);
  assert (_regmatch_valid);

  const regmatch_t & rm = _regmatch[whichdollar];

#ifdef DEBUG_DOLLAR
  cerr << "Looking for dollar " << whichdollar << " in '" << _string_which_was_matched << "', range " << _regmatch[0].rm_so << " - " << _regmatch[0].rm_eo << endl;
  cerr << "Match starts at " << rm.rm_so << " to " << rm.rm_eo << endl;
#endif

  zresult.set (_string_which_was_matched.rawchars () + rm.rm_so, rm.rm_eo - rm.rm_so);

  return 1;
}

const const_IWSubstring 
IW_lregex::dollar (int whichdollar) const
{
  assert (whichdollar >= 0 && whichdollar <= int (_compiled_regular_expression.re_nsub));
  assert (NULL != _regmatch);
  assert ("" != _string_which_was_matched);
  assert (_regmatch_valid);

  const regmatch_t & rm = _regmatch[whichdollar];

  return const_IWSubstring (_string_which_was_matched.rawchars () + rm.rm_so, rm.rm_eo - rm.rm_so);
}

int
IW_lregex::dollar (int whichdollar, int & zresult)
{
  const_IWSubstring tmp;

  if (! dollar (whichdollar, tmp))
    return 0;

  return tmp.numeric_value (zresult);
}

int
IW_lregex::dollar (int whichdollar, float & zresult)
{
  const_IWSubstring tmp;

  if (! dollar (whichdollar, tmp))
    return 0;

  return tmp.numeric_value (zresult);
}

int
IW_lregex::dollar (int whichdollar, double & zresult)
{
  const_IWSubstring tmp;

  if (! dollar (whichdollar, tmp))
    return 0;

  return tmp.numeric_value (zresult);
}

const regmatch_t *
IW_lregex::regmatch () const
{
  assert (_compiled);    // must be active
  assert (ok ());

  return _regmatch;
}

int
IW_lregex::number_subexpressions () const
{
  assert (ok ());
  assert (_compiled_regular_expression.re_nsub >= 0);

  return _compiled_regular_expression.re_nsub;
}

int
IW_lregex::s (const IWString & s, const const_IWSubstring & zto,
          IWString & zresult)
{
  if (! matches_save_subexpressions (s))
    return 0;

  int nmatch = _compiled_regular_expression.re_nsub;

  zresult.resize (s.nchars () + (nmatch - 1) * zto.nchars ());

  int jptr = 0;       // the char in S being processed

  for (int i = 1; i <= nmatch; i++)
  {
    if (_regmatch[i].rm_so - jptr > 0)
      zresult.strncat (s.rawchars () + jptr, _regmatch[i].rm_so - jptr);

    zresult += zto;
    jptr = _regmatch[i].rm_eo;
  }

// Is there anything beyond the end

  if (jptr < s.nchars ())
    zresult.strncat (s.rawchars () + jptr, s.nchars () - jptr);

  return nmatch;
}

/*
  In place transformation
*/

int
IW_lregex::s (IWString & s, const const_IWSubstring & zto)
{
  IWString tmp;

  int rc = IW_lregex::s (s, zto, tmp);

  if (rc)
    s = tmp;

  return rc;
}

int
iwcrex_s (IWString & s, const const_IWSubstring & regexp, const const_IWSubstring & zto)
{
  IW_lregex rx;

  if (! rx.set_pattern (regexp))
  {
    cerr << "iwcrex_s: cannot parse regexp pattern '" << regexp << "'\n";
    return 0;
  }

  return rx.s (s, zto);
}
