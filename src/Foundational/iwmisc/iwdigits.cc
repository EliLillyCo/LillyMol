#include <stdlib.h>

#include "misc.h"
#include "iwdigits.h"

void
IWDigits::_default_values()
{
  return;
}

IWDigits::IWDigits()
{
}

IWDigits::IWDigits(int mysize) : iwaray<IWString>(mysize)
{
  _fill_in_the_digits();
}

void
IWDigits::debug_print(std::ostream & output) const
{
  output << "IWDigits::storing " << _number_elements << " precomputed digits\n";
  if (_leading_space.length())
    output << " leading  string '" << _leading_space << "'\n";
  if (_appended_string.length())
    output << " appended string '" << _appended_string << "'\n";

  int istop = _number_elements;
  if (istop > 11)
    istop = 11;

  for (int i = 0; i < istop; ++i)
  {
    output << " i = " << i << " '" << _things[i] << "'\n";
  }

  return;
}

int
IWDigits::set_include_leading_space(int s)
{
  if (s && ' ' == _leading_space)   // no change
    return 1;

  if (! s && 0 == _leading_space.length())   // no change
    return 1;

  if (s)
    _leading_space = ' ';
  else
    _leading_space.resize(0);

  _fill_in_the_digits();

  return 1;
}

int
IWDigits::initialise(int mysize)
{
  if (! iwaray<IWString>::resize(mysize))
  {
    std::cerr << "IWDigits::initialise: cannot size to " << mysize << '\n';
    return 0;
  }

  _fill_in_the_digits();

  return 1;
}

void
IWDigits::_fill_in_the_digits()
{
  for (int i = 0; i < _number_elements; i++)
  {
    IWString & d = _things[i];

    if (_leading_space.length())
      d = _leading_space;

    d.append_number(i);

    if (_appended_string.length())
      d << _appended_string;

//  cerr << " i = " << i << " value '" << d << "'\n";
  }

  return;
}

/*template <typename T>
int
IWDigits::append_number(T & buffer, int zdigit) const
{
  if (zdigit < 0)
  {
    zdigit = - zdigit;     // convert to something >= 0

    if (_leading_space.length() > 0)
    {
      buffer << _leading_space;
      buffer << '-';
      buffer.append_number( zdigit);
      return 1;
    }
    else
      buffer << '-';
  }

  if (zdigit < _number_elements)
  {
    buffer << _things[zdigit];
    return 1;
  }

  if (_leading_space.length())
    buffer << _leading_space << zdigit;
  else
    buffer.append_number(zdigit);

  if (_appended_string.length())
    buffer << _appended_string;

  return 1;
}*/

template <typename T, typename I>
int
IWDigits::append_number(T & buffer, I zdigit) const
{
  if (zdigit < 0)
    return _append_negative_number(buffer, zdigit);

  return _append_number(buffer, zdigit);
}

/*
  We really do not handle negative numbers - that could be changed with a separate array...
*/

template <typename T, typename I>
int
IWDigits::_append_negative_number(T & buffer, I zdigit) const
{
  if (0 == _leading_space.length())
  {
    buffer << '-';
    return _append_number(buffer, - zdigit);
  }

  buffer << _leading_space;

  buffer << zdigit;    // just use default capability of underlying stream

  if (_appended_string.length())
    buffer << _appended_string;

  return 1;
}

template <typename T, typename I>
int
IWDigits::_append_number(T & buffer, I zdigit) const
{
  assert (zdigit >= 0);

  if (zdigit < static_cast<I>(_number_elements))
    buffer << _things[zdigit];
  else
  {
    if (_leading_space.length())
      buffer << _leading_space;
    buffer << zdigit;
    if (_appended_string.length() > 0)
      buffer << _appended_string;
  }

  return 1;
}

template int IWDigits::append_number<std::ostream>(std::ostream &, int) const;
template int IWDigits::append_number(IWString &, int) const;
template int IWDigits::append_number(IWString_and_File_Descriptor &, int) const;
template int IWDigits::append_number(IWString_and_File_Descriptor &, unsigned int) const;
template int IWDigits::_append_number(IWString_and_File_Descriptor &, unsigned int) const;
template int IWDigits::_append_number(IWString &, unsigned int) const;
template int IWDigits::_append_number(std::ostream &, unsigned int) const;
template int IWDigits::_append_negative_number(std::ostream &, int) const;
template int IWDigits::_append_negative_number(IWString &, int) const;
template int IWDigits::_append_number<IWString_and_File_Descriptor, int>(IWString_and_File_Descriptor&, int) const;
template int IWDigits::_append_number<IWString, int>(IWString&, int) const;
template int IWDigits::append_number<IWString, unsigned int>(IWString&, unsigned int) const;
template int IWDigits::append_number<std::ostream, unsigned int>(std::ostream&, unsigned int) const;

const IWString &
IWDigits::string_for_digit(int zdigit) const
{
  assert (zdigit >= 0 && zdigit < _number_elements);

  return _things[zdigit];
}

int
IWDigits::set_leading_string(const const_IWSubstring & s)
{
  if (_leading_space == s)    // no change
    return 1;

  _leading_space = s;

  if (nullptr != _things)
    _fill_in_the_digits();
    
  return 1;
}

int
IWDigits::set_leading_string(char s) {
  IWString tmp(s);
  return set_leading_string(tmp);
}

int
IWDigits::append_to_each_stored_string(const const_IWSubstring & s)
{
  for (int i = 0; i < _number_elements; i++)    // may be empty
  {
    _things[i] << s;
  }

  _appended_string = s;

  return 1;
}

int
IWDigits::append_to_each_stored_string(char s) {
  for (int i = 0; i < _number_elements; i++)    // may be empty
  {
    _things[i] << s;
  }

  _appended_string = s;

  return 1;
}

