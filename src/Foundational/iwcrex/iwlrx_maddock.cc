#include <stdlib.h>
#include <assert.h>

#include "iwlrx_maddock.h"

void
IW_lrx_Maddock::_default_values ()
{
  _regcomp_flags = REG_NOSUB;
  _reg_eflags = 0;

  _regmatch = NULL;
  _regmatch_valid = 0;

  _compiled = 0;
}

IW_lrx_Maddock::IW_lrx_Maddock ()
{
  _default_values ();
}

IW_lrx_Maddock::IW_lrx_Maddock (const char * s)
{
  _default_values ();

  set_pattern (s, strlen (s));
  
  return;
}

IW_lrx_Maddock::IW_lrx_Maddock (const const_IWSubstring & s)
{
  _default_values ();

  set_pattern (s.rawchars (), s.nchars ());
  
  return;
}

IW_lrx_Maddock::IW_lrx_Maddock (const IWString & s)
{
  _default_values ();

  set_pattern (s.rawchars (), s.nchars ());
  
  return;
}

IW_lrx_Maddock::~IW_lrx_Maddock ()
{
  if (NULL != _regmatch)
    delete [] _regmatch;

  if (_compiled)
    regfree (& _compiled_regular_expression);

  return;
}

int
IW_lrx_Maddock::ok () const
{
  if (0 == _compiled_regular_expression.re_nsub && NULL != _regmatch)
  {
    cerr << "IW_lrx_Maddock::ok: re_nsub/_regmatch mismatch\n";
    return 0;
  }

  if (_compiled_regular_expression.re_nsub && NULL == _regmatch)
  {
    cerr << "IW_lrx_Maddock::ok: _regmatch/re_nsub mismatch\n";
    return 0;
  }

  return 1;
}

void
IW_lrx_Maddock::set_use_extended (int s)
{
  if (s)
    _regcomp_flags |= REG_EXTENDED;
  else
    _regcomp_flags &= (~REG_EXTENDED);
}

//#define DEBUG_SET_PATTERN

int
IW_lrx_Maddock::set_pattern (const char * s, int nchars)
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

  int rc = regncomp (&_compiled_regular_expression, s, nchars, _regcomp_flags);
  if (0 != rc)
  {
    cerr << "IW_lrx_Maddock::set_pattern: cannot compile '";
    cerr.write (s, nchars);
    cerr << "'\n";
    return 0;
  }

#ifdef DEBUG_SET_PATTERN
  cerr << "Compiled '";
  cerr.write (s, nchars);
  cerr << "', re_nsub = " << _compiled_regular_expression.re_nsub << endl;
#endif

  _compiled = 1;

  if (0 == _compiled_regular_expression.re_nsub)    // no fields, we are done
    return 1;

  _regmatch = new regmatch_t[_compiled_regular_expression.re_nsub + 1];

  assert (NULL != _regmatch);

  return 1;
}

int
IW_lrx_Maddock::set_pattern (const IWString & s)
{
  return set_pattern (s.rawchars (), s.length ());
}

int
IW_lrx_Maddock::set_pattern (const const_IWSubstring & s)
{
  return set_pattern (s.rawchars (), s.length ());
}

int
IW_lrx_Maddock::set_pattern (const char * s)
{
  return set_pattern (s, ::strlen (s));
}

int
IW_lrx_Maddock::matches (const char * s)
{
  return matches (s, strlen (s));
}


int
IW_lrx_Maddock::matches (const IWString & s)
{
  return matches (s.rawchars (), s.length ());
}

int
IW_lrx_Maddock::matches (const const_IWSubstring & s)
{
  return matches (s.rawchars (), s.length ());
}

int
IW_lrx_Maddock::matches (const char * s, int nchars)
{
  assert (ok ());

#ifdef DEBUG_LRX_MADDOCK_MATCHES
  cerr << "Looking for match '";
  cerr.write (s, nchars) << "'\n";
#endif

  if (! _compiled)
    return 0;

  _regmatch_valid = 0;    // we are not saving subexpressions

  if (0 == regnexec (&_compiled_regular_expression, s, nchars, 
                     0, NULL,
                     _reg_eflags))
  {
    return 1;    // match
  }

  return 0;
}

int
IW_lrx_Maddock::matches_save_subexpressions (const char * s)
{
  return matches_save_subexpressions (s, ::strlen (s));
}

int
IW_lrx_Maddock::matches_save_subexpressions (const const_IWSubstring & s)
{
  return matches_save_subexpressions (s.rawchars (), s.nchars ());
}

int
IW_lrx_Maddock::matches_save_subexpressions (const IWString & s)
{
  return matches_save_subexpressions (s.rawchars (), s.nchars ());
}

int
IW_lrx_Maddock::matches_save_subexpressions (const char * s, int nchars)
{
  assert (ok ());

  if (! _compiled)
    return 0;

  if (0 == regnexec (&_compiled_regular_expression, s, nchars, 
                     _compiled_regular_expression.re_nsub + 1, _regmatch,
                     _reg_eflags))
  {
    _regmatch_valid = 1;
    _string_which_was_matched.set (s, nchars);
    return 1;    // match
  }

  _regmatch_valid = 0;

  return 0;
}

/*
  After a match, the caller wants to know the components which matched
*/

//#define DEBUG_DOLLAR

int
IW_lrx_Maddock::dollar (int whichdollar, const_IWSubstring & zresult) const
{
  assert (whichdollar >= 0 && whichdollar <= int (_compiled_regular_expression.re_nsub));
  assert (NULL != _regmatch);
  assert ("" != _string_which_was_matched);
  assert (_regmatch_valid);

#ifdef DEBUG_DOLLAR
  cerr << "Retrieving dollar " << whichdollar << endl;
  for (unsigned int i = 0; i <= _compiled_regular_expression.re_nsub; i++)
  {
    const regmatch_t & rm = _regmatch[i];

    cerr << "Dollar " << i << " starts at " << rm.rm_so << " and goes to " << rm.rm_eo << endl;
  }
#endif

  const regmatch_t & rm = _regmatch[whichdollar];

  zresult.set (_string_which_was_matched.rawchars () + rm.rm_so, rm.rm_eo - rm.rm_so);

  return 1;
}

const const_IWSubstring 
IW_lrx_Maddock::dollar (int whichdollar) const
{
  assert (whichdollar >= 0 && whichdollar <= int (_compiled_regular_expression.re_nsub));
  assert (NULL != _regmatch);
  assert ("" != _string_which_was_matched);
  assert (_regmatch_valid);

  const regmatch_t & rm = _regmatch[whichdollar];

  return const_IWSubstring (_string_which_was_matched.rawchars () + rm.rm_so, rm.rm_eo - rm.rm_so);
}

const regmatch_t *
IW_lrx_Maddock::regmatch () const
{
  assert (_compiled);    // must be active
  assert (ok ());

  return _regmatch;
}

int
IW_lrx_Maddock::number_subexpressions () const
{
  assert (ok ());
  assert (_compiled);

  return _compiled_regular_expression.re_nsub;
}

int
IW_lrx_Maddock::s (const IWString & s, const const_IWSubstring & zto,
                   IWString & zresult)
{
  if (! matches_save_subexpressions (s))
    return 0;

  int nmatch = _compiled_regular_expression.re_nsub + 1;

  zresult.resize_keep_storage (s.nchars () + (nmatch - 1) * zto.nchars ());

  int jptr = 0;       // the char in S being processed

  for (int i = 1; i < nmatch; i++)
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
IW_lrx_Maddock::s (IWString & s, const const_IWSubstring & zto)
{
  IWString tmp;

  int rc = IW_lrx_Maddock::s (s, zto, tmp);

  if (rc)
    s = tmp;

  return rc;
}

/*template <class X>
int
IW_lrx_Maddock::qdollar (int whichdollar, X & zresult)
{
  const_IWSubstring tmp;

  if (! dollar (whichdollar, tmp))
    return 0;

  return tmp.numeric_value (zresult);
}*/

int
IW_lrx_Maddock::dollar (int whichdollar, int & zresult)
{
  const_IWSubstring tmp;

  if (! dollar (whichdollar, tmp))
    return 0;

  return tmp.numeric_value (zresult);
}

int
IW_lrx_Maddock::dollar (int whichdollar, float & zresult)
{
  const_IWSubstring tmp;

  if (! dollar (whichdollar, tmp))
    return 0;

  return tmp.numeric_value (zresult);
}

int
IW_lrx_Maddock::dollar (int whichdollar, double & zresult)
{
  const_IWSubstring tmp;

  if (! dollar (whichdollar, tmp))
    return 0;

  return tmp.numeric_value (zresult);
}

int
iwcrex_s (IWString & s, const const_IWSubstring & regexp, const const_IWSubstring & zto)
{
  IW_lrx_Maddock rx;

  if (! rx.set_pattern (regexp))
  {
    cerr << "iwcrex_s: cannot parse regexp pattern '" << regexp << "'\n";
    return 0;
  }

  return rx.s (s, zto);
}
