#ifndef MASQUERADING_AS_BYTE_H
#define MASQUERADING_AS_BYTE_H

/*
  When storing distance matrices, we can save a lot of space by encoding
  floats as bytes. We need to know the range
*/

#include "iwstring_data_source.h"

#include <iostream>
using std::cerr;
using std::endl;

template <typename T>
class Masquerading_as_Byte
{
  protected:
    double _minval;
    double _maxval;
    double _dx;

//  We can save lots of time by pre-computing the float/double that corresponds
//  to each of the bytes

    T _byte_to_T[256];

//  during output operations, it can be handy to have the text representation of each number

    IWString _string_value[256];

//  byte value 0 and byte value 255 are special, and correspond to missing values

    T _missing0;
    T _missing255;

//  private functions

    int _parse_map_directive (const const_IWSubstring &);

  public:
    Masquerading_as_Byte ();

    int debug_print (std::ostream &) const;
    int do_write (std::ostream &) const;
    int do_write (IWString_and_File_Descriptor &) const;

    int set_range (T mi, T ma);
    void set_min (T mi) { _minval = mi;};
    void set_max (T ma) { _maxval = ma;};

    int in_range (T) const;

    int recognised (const const_IWSubstring &, const const_IWSubstring &, int &);
    int fully_specified ();

    unsigned char convert_to_byte (T) const;
    int convert_to_byte (T, unsigned char &) const;

    T convert_from_byte (unsigned char b) const { return _byte_to_T[b];}

    const IWString & string_value (int i) const { return _string_value[i];}

    int initialise_from_distance_matrix_header_records (const const_IWSubstring & fname);
    int initialise_from_distance_matrix_header_records (iwstring_data_source &);
};

#ifdef MASQUERADING_AS_BYTE_IMPLEMENTATION

#include "misc.h"

#include "Masquerading_as_Byte.h"

template <typename T>
Masquerading_as_Byte<T>::Masquerading_as_Byte ()
{
  _minval = 0.0;
  _maxval = 0.0;
  _dx     = 0.0;

  _missing0   = static_cast<T> (0);
  _missing255 = static_cast<T> (0);

  set_vector (_byte_to_T, 256, static_cast<T> (0.0));
}

template <typename T>
int
Masquerading_as_Byte<T>::debug_print (std::ostream & os) const
{
  T zero = static_cast<T> (0);

  os << "Masquerading_as_Byte<T>::debug_print:";
  if (zero == _minval && zero == _maxval)
  {
    os << " not initialised\n";
    return os.good ();
  }

  os << " minval " << _minval << " maxval " << _maxval << endl;

  os << "Translation table\n";
  for (int i = 0; i < 256; i++)
  {
    os << "  byte " << i << " translated to " << _byte_to_T[i] << endl;
  }

  return os.good ();
}

template <typename T>
int
Masquerading_as_Byte<T>::in_range (T v) const
{
  if (v < _minval)
    return 0;

  if (v > _maxval)
    return 0;

  return 1;
}

template <typename T>
int
Masquerading_as_Byte<T>::convert_to_byte (T v,
                                          unsigned char & b) const
{
  if (v < _minval)
  {
    cerr << "Masquerading_as_Byte::convert_to_byte: value " << v << " out of range, minval = " << _minval << endl;
    return 0;
  }

  if (v > _maxval)
  {
    cerr << "Masquerading_as_Byte::convert_to_byte: value " << v << " out of range, maxval = " << _maxval << endl;
    return 0;
  }

  b = convert_to_byte (v);

  return 1;
}

/*
  We have a value and need to know which of the 254 buckets it should go in
*/

template <typename T>
unsigned char
Masquerading_as_Byte<T>::convert_to_byte (T f) const
{
  assert (0.0 != _dx);
  
  if (f < _minval || f > _maxval)
  {
    cerr << "Masquerading_as_Byte::_convert_to_byte: out of range " << f << " range is " << _minval << " to " << _maxval << endl;
    abort ();
  }

// Kludge around to get the nearest int

  T v = (f - _minval) / _dx;

  int i1 = static_cast<int> (v);
  int i2 = static_cast<int> (v + static_cast<T> (0.05));

  if (i2 > i1)
    return i2 + 1;
  else
    return i1 + 1;
}

template <typename T>
int
Masquerading_as_Byte<T>::do_write (std::ostream & output) const
{
  if (0.0 != _minval || 0.0 != _maxval)
  {
    output << "minval " << _minval << '\n';
    output << "maxval " << _maxval << '\n';
  }

  for (int i = 0; i < 256; i++)
  {
    output << "map " << i << " " << _byte_to_T[i] << '\n';
  }

  return output.good ();
}

template <typename T>
int
Masquerading_as_Byte<T>::do_write (IWString_and_File_Descriptor & output) const
{
  if (0.0 != _minval || 0.0 != _maxval)
  {
    output << "minval " << _minval << '\n';
    output << "maxval " << _maxval << '\n';
  }

  for (int i = 0; i < 256; i++)
  {
    output << "map " << i << " " << _byte_to_T[i] << '\n';
  }

  return output.good ();
}

template <typename T>
int
Masquerading_as_Byte<T>::set_range (T mi, T ma)
{
  assert (mi < ma);

  _minval = mi;
  _maxval = ma;
  _dx     = (_maxval - _minval) / 254.0;

  assert (static_cast<T> (0) != _dx);

  _byte_to_T[0] = _missing0;

  for (int i = 1; i < 255; i++)
  {
    _byte_to_T[i] = _minval + static_cast<T> (i - 1) * _dx;
  }
  _byte_to_T[255] = _missing255;

  for (int i = 0; i < 256; i++)
  {
    IWString & s = _string_value[i];
    if (0 == s.length ())
      s << _byte_to_T[i];
  }

  return 1;
}

template <typename T>
int
Masquerading_as_Byte<T>::_parse_map_directive (const const_IWSubstring & buffer)
{
  const_IWSubstring s, sv;
  if (! buffer.split (s, ' ', sv))
  {
    cerr << "Masquerading_as_Byte::_parse_map_directive: must be at least two tokens in a map directive\n";
    cerr << buffer << endl;
    return 0;
  }

  int ndx;
  if (! s.numeric_value (ndx) || ndx < 0 || ndx > 255)
  {
    cerr << "Masquerading_as_Byte::_parse_map_directive: invalid index '" << buffer << "'\n";
    return 0;
  }

  T v;
  if (! sv.numeric_value (v))
  {
    cerr << "Masquerading_as_Byte::_parse_map_directive: invalid numeric '" << buffer << "'\n";
    return 0;
  }

  _byte_to_T[ndx] = v;
  _string_value[ndx] = sv;

  return 1;
}

template <typename T>
int
Masquerading_as_Byte<T>::recognised (const const_IWSubstring & zdirective,
                                  const const_IWSubstring & zvalue,
                                  int & error_encountered)
{
  error_encountered = 0;

  if ("minval" == zdirective)
  {
    if (! zvalue.numeric_value (_minval))
      error_encountered = 1;
  }
  else if ("maxval" == zdirective)
  {
    if (! zvalue.numeric_value (_maxval))
      error_encountered = 1;
  }
  else if ("map" == zdirective)
  {
    if (! _parse_map_directive (zvalue))
      error_encountered = 1;
  }
  else
    return 0;

  return 1;
}

template <typename T>
int
Masquerading_as_Byte<T>::fully_specified ()
{
  if (0.0 == _minval && 0.0 == _maxval)
  {
    cerr << "Masquerading_as_Byte::fully_specified: zero extremeties\n";
    return 0;
  }

  if (_minval >= _maxval)
  {
    cerr << "Masquerading_as_Byte::fully_specified: invalid extremeties " <<_minval << " to " << _maxval << endl;
    return 0;
  }

  return set_range (_minval, _maxval);
}

template <typename T>
int
Masquerading_as_Byte<T>::initialise_from_distance_matrix_header_records (const const_IWSubstring & fname)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Masquerading_as_Byte::initialise_from_distance_matrix_header_records:cannot open '" << fname << "'\n";
    return 0;
  }

  return initialise_from_distance_matrix_header_records (input);
}

template <typename T>
int
split_into_directive_and_value (const const_IWSubstring & buffer,
                                const_IWSubstring & directive,
                                T & zvalue)
{
  const_IWSubstring token;

  if (! buffer.split (directive, ' ', token))
    return 0;

  return token.numeric_value (zvalue);
}

template <typename T>
int
Masquerading_as_Byte<T>::initialise_from_distance_matrix_header_records (iwstring_data_source & input)
{
  set_vector (_byte_to_T, 256, static_cast<T> (0.0));

  _minval = -9.8;    // just some arbitrary number
  _maxval = -9.8;
  _dx = -9.8;

  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if ('|' == buffer)
      break;

    const_IWSubstring directive, zvalue;

    if (! buffer.split (directive, ' ', zvalue))
      continue;

    int error_encountered;

    if (! recognised (directive, zvalue, error_encountered))
      continue;

    if (error_encountered)
    {
      cerr << "Masquerading_as_Byte::initialise_from_distance_matrix_header_records:invalid input '" << buffer << "'\n";
      return 0;
    }
  }

  if (_minval < 0.0 && _maxval < 0.0)
  {
    cerr << "Masquerading_as_Byte::initialise_from_distance_matrix_header_records:not initialised\n";
    return 0;
  }

  _dx = (_maxval - _minval) / 254.0;

  return 1;
}

#endif
#endif
