#ifndef IWCREX_IMPLEMENTATION_H
#define IWCREX_IMPLEMENTATION_H

#include "iwcrex.h"


template <typename T>
void
IW_Regular_Expression_Template<T>::_default_values ()
{
  _matcher = NULL;
  _use_extended = 1;
  _ignore_case = 0;

  return;
}

template <typename T>
IW_Regular_Expression_Template<T>::IW_Regular_Expression_Template ()
{
  _default_values ();

  return;
}

template <typename T>
IW_Regular_Expression_Template<T>::IW_Regular_Expression_Template (const char * s)
{
  _default_values ();

  (void) set_pattern (s);

  return;
}

template <typename T>
IW_Regular_Expression_Template<T>::IW_Regular_Expression_Template (const char * s, int len)
{
  _default_values ();

  const const_IWSubstring tmp;

  (void) set_pattern (s, len);

  return;
}

template <typename T>
IW_Regular_Expression_Template<T>::IW_Regular_Expression_Template (const const_IWSubstring & s)
{
  _default_values ();

  (void) set_pattern (s);

  return;
}

template <typename T>
IW_Regular_Expression_Template<T>::IW_Regular_Expression_Template (const IWString & s)
{
  _default_values ();

  (void) set_pattern (s);

  return;
}

template <typename T>
IW_Regular_Expression_Template<T>::~IW_Regular_Expression_Template ()
{
  if (NULL != _matcher)
    delete _matcher;

  return;
}

template <typename T>
int
IW_Regular_Expression_Template<T>::ok () const
{
  if (0 == _source.length () && NULL == _matcher)
    return 1;

  if (_source.length () > 0 && NULL != _matcher)
    return 1;

  cerr << "Source '" << _source << "' and Matcher " << (NULL != _matcher) << endl;

  return 0;
}

template <typename T>
int
IW_Regular_Expression_Template<T>::debug_print (ostream & os) const
{
  os << "Details on IW_Regular_Expression_Template:";

  if (! ok ())
    os << " OK fails";

  if (NULL == _matcher)
  {
    os << " inactive\n";
    return os.good ();
  }

  os << " source '" << _source << "'\n";

  return os.good ();
}

template <typename T>
int
IW_Regular_Expression_Template<T>::set_pattern (const char * s)
{
  assert (ok ());

//cerr << "Setting patterh to '" << s << "'\n";

  if (NULL == _matcher)
    _matcher = new T;

  _matcher->set_use_extended (_use_extended);
  _matcher->set_ignore_case (_ignore_case);

  _source.resize_keep_storage (0);

//cerr << "Outer wrapper setting pattern to '" << s << "'\n";
  if (! _matcher->set_pattern (s))
  {
    cerr << "IW_Regular_Expression_Template<T>::set_pattern: cannot parse '" << s << "'\n";

    delete _matcher;

    _matcher = NULL;

    return 0;
  }

  _source = s;

  return 1;
}

template <typename T>
int
IW_Regular_Expression_Template<T>::set_pattern (const char * s, int len)
{
  assert (ok ());

  const const_IWSubstring tmp (s, len);

  return set_pattern (tmp);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::set_pattern (const const_IWSubstring & s)
{
  assert (ok ());

  if (NULL == _matcher)
    _matcher = new T;

  _matcher->set_use_extended (_use_extended);

  _source.resize_keep_storage (0);

  if (! _matcher->set_pattern (s))
  {
    cerr << "IW_Regular_Expression_Template<T>::set_pattern: cannot parse '" << s << "'\n";

    delete _matcher;

    _matcher = NULL;

    return 0;
  }

  _source = s;

  return 1;
}

template <typename T>
int
IW_Regular_Expression_Template<T>::set_pattern (const IWString & s)
{
  assert (ok ());

  if (NULL == _matcher)
    _matcher = new T;

  _source.resize_keep_storage (0);

  if (! _matcher->set_pattern (s))
  {
    cerr << "IW_Regular_Expression_Template<T>::set_pattern: cannot parse '" << s << "'\n";

    delete _matcher;

    _matcher = NULL;

    return 0;
  }

  _source = s;

  return 1;
}

template <typename T>
int
IW_Regular_Expression_Template<T>::matches (const char * s)
{
  if (NULL == _matcher)
  {
    cerr << "IW_Regular_Expression_Template::matches: matcher not defined\n";
    return 0;
  }

  return _matcher->matches (s);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::matches (const char * s, int len)
{
  if (NULL == _matcher)
  {
    cerr << "IW_Regular_Expression_Template::matches: matcher not defined\n";
    return 0;
  }

  return _matcher->matches (s, len);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::matches (const const_IWSubstring & s)
{
  if (NULL == _matcher)
  {
    cerr << "IW_Regular_Expression_Template::matches: matcher not defined\n";
    return 0;
  }

  return _matcher->matches (s);
}
template <typename T>
int
IW_Regular_Expression_Template<T>::matches (const IWString & s)
{
  if (NULL == _matcher)
  {
    cerr << "IW_Regular_Expression_Template::matches: matcher not defined\n";
    return 0;
  }

  return _matcher->matches (s);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::matches_save_subexpressions (const char * s)
{
  if (NULL == _matcher)
  {
    cerr << "IW_Regular_Expression_Template<T>::matches_save_subexpressions: matcher not defined\n";
    return 0;
  }

  return _matcher->matches_save_subexpressions (s);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::matches_save_subexpressions (const char * s, int len)
{
  if (NULL == _matcher)
  {
    cerr << "IW_Regular_Expression_Template<T>::matches_save_subexpressions: matcher not defined\n";
    return 0;
  }

  return _matcher->matches_save_subexpressions (s, len);
}


template <typename T>
int
IW_Regular_Expression_Template<T>::matches_save_subexpressions (const const_IWSubstring & s)
{
  if (NULL == _matcher)
  {
    cerr << "IW_Regular_Expression_Template<T>::matches_save_subexpressions: matcher not defined\n";
    return 0;
  }

  return _matcher->matches_save_subexpressions (s);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::matches_save_subexpressions (const IWString & s)
{
  if (NULL == _matcher)
  {
    cerr << "IW_Regular_Expression_Template<T>::matches_save_subexpressions: matcher not defined\n";
    return 0;
  }

  return _matcher->matches_save_subexpressions (s);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::number_subexpressions () const
{
  assert (NULL != _matcher);

  return _matcher->number_subexpressions ();
}

template <typename T>
int
IW_Regular_Expression_Template<T>::dollar (int whichdollar, const_IWSubstring & zresult) const
{
  assert (NULL != _matcher);
  return _matcher->dollar (whichdollar, zresult);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::dollar (int whichdollar, IWString & zresult) const
{
  assert (NULL != _matcher);

  return _matcher->dollar (whichdollar, zresult);
}

template <typename T>
const const_IWSubstring
IW_Regular_Expression_Template<T>::dollar (int d) const
{
  assert (NULL != _matcher);

  return _matcher->dollar (d);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::dollar (int whichdollar, int & zresult)
{
  assert (NULL != _matcher);

  return _matcher->dollar (whichdollar, zresult);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::dollar (int whichdollar, float & zresult)
{
  assert (NULL != _matcher);

  return _matcher->dollar (whichdollar, zresult);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::dollar (int whichdollar, double & zresult)
{
  assert (NULL != _matcher);

  return _matcher->dollar (whichdollar, zresult);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::s (IWString & s, const const_IWSubstring & f)
{
  assert (NULL != _matcher);

  return _matcher->s (s, f);
}

template <typename T>
int
IW_Regular_Expression_Template<T>::s (IWString & s, const const_IWSubstring & f, IWString & zresult)
{
  assert (NULL != _matcher);

  return _matcher->s (s, f, zresult);
}

template<typename T>
void
IW_Regular_Expression_Template<T>::set_use_extended (int s)
{
  _use_extended = s;
}

#endif
