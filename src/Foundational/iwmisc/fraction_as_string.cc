#include <stdlib.h>

#include "iwdigits.h"

Fraction_as_String::Fraction_as_String()
{
  _fraction = nullptr;

  _nbuckets = 0;

  return;
}

Fraction_as_String::~Fraction_as_String()
{
  if (nullptr != _fraction)
    delete [] _fraction;

  return;
}

int
Fraction_as_String::initialise(float minvl, float maxvl, int digits)
{
  if (nullptr != _fraction)
    delete [] _fraction;

  assert(digits > 0);

  _nbuckets = 1;
  for (int i = 0; i < digits; i++)
  {
    _nbuckets = 10 * _nbuckets;
  }

//std::cerr << "Allocating " << _nbuckets << " buckets\n";

  _fraction = new IWString[_nbuckets + 1];

  if (nullptr == _fraction)
  {
    std::cerr << "Fraction_as_String::initialise:cannot allocate " << _nbuckets << " string representations\n";
    return 0;
  }

  _minval = minvl;
  _maxval = maxvl;
  _dx = (_maxval - _minval) / static_cast<double>(_nbuckets);
  _digits = digits;

  return _fill_string_data();
}

int
Fraction_as_String::_fill_string_data()
{
  assert(nullptr != _fraction);

//std::cerr << "Fraction_as_String:_fill_string_data:leading space '" << _leading_space << "'\n";

  for (int i = 0; i < _nbuckets; i++)
  {
    float d = static_cast<float>(_minval + static_cast<double>(i) * _dx);

    IWString & f = _fraction[i];

    f.resize(_digits + 2 + _leading_space.length());

    if (_leading_space.length())
      f << _leading_space;

    f.append_number(d, _digits);

    if (0.0F != d)
    {
      if (! f.contains('.'))
        f << '.';

      while (f.length() < (_digits + 2 + _leading_space.length()))
      {
        f << '0';
      }
    }
//  std::cerr << " i = " << i << " string '" << _fraction[i] << "'\n";
  }

  if (_leading_space.length())
    _fraction[_nbuckets] << _leading_space;

  _fraction[_nbuckets].append_number(_maxval, _digits);

  return 1;
}

void
Fraction_as_String::_append_number_no_string_rep(IWString & s,
                                                 float f) const
{
  if (_leading_space.length())
    s += _leading_space;

  if (_digits > 0)
    s.append_number(f, _digits);
  else
    s.append_number(f);

  return;
}

void
Fraction_as_String::append_number(IWString & s,
                                  float f) const
{
  if (nullptr == _fraction)
  {
    _append_number_no_string_rep(s, f);
    return;
  }

  if (f < _minval || f > _maxval)
  {
    _append_number_no_string_rep(s, f);
    return;
  }

  int i = static_cast<int>( (f - _minval) / _dx + 0.49999);

  assert (i >= 0 && i <= _nbuckets);

  s << _fraction[i];

  return;
}

int
Fraction_as_String::set_include_leading_space(int s)
{
  if (s && ' ' == _leading_space)    // no change
    return 1;

  if (! s && 0 == _leading_space.length())    // no change
    return 1;

  if (s)
    _leading_space = ' ';
  else
    _leading_space.resize(0);

  if (nullptr == _fraction)
    return 1;

  return _fill_string_data();
}


int
Fraction_as_String::set_leading_string(const const_IWSubstring & s)
{
  if (_leading_space == s)   // no change
    return 1;

  _leading_space = s;

  if (nullptr == _fraction)
    return 1;

  for (int i = 0; i < _nbuckets; ++i)
  {
    _fraction[i].resize_keep_storage(0);
  }

  return _fill_string_data();
}

int
Fraction_as_String::set_leading_string(char c) {
  IWString tmp(c);
  return set_leading_string(tmp);
}

const IWString &
Fraction_as_String::string_for_fraction(float f) const
{
  assert(f >= _minval && f <= _maxval);
  assert(nullptr != _fraction);

  int i = static_cast<int>( (f - _minval) / _dx + 0.49999);

  return _fraction[i];
}

int
Fraction_as_String::append_to_each_stored_string(const const_IWSubstring & s)
{
  if (nullptr == _fraction)
  {
    std::cerr << "Fraction_as_String::append_to_each_stored_string:no strings\n";
    return 0;
  }

  for (int i = 0; i <= _nbuckets; i++)
  {
    _fraction[i] << s;
  }

  return _nbuckets;
}
