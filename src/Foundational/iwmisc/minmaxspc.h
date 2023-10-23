#ifndef FOUNDATIONAL_IWMISC_MINMAXSPC_H
#define FOUNDATIONAL_IWMISC_MINMAXSPC_H

#include <iostream>

#include "Foundational/iwmisc/iwconfig.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/iwstring.h"
#include "set_or_unset.h"
#include "iwarchive.h"

/*
  This specifier consists of ONE OF
    a minimum value
    a maximum value
    a list of allowed values.
*/

template <typename T>
class Min_Max_Specifier : public iwarchive<T>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using iwarchive<T>::_number_elements;
    using iwarchive<T>::_things;
    using iwarchive<T>::_match_any;
#endif
  private:

//  For efficiency we record whether or not we are set

    int _is_set;

    Set_or_Unset<T> _min_val;
    Set_or_Unset<T> _max_val;
#ifdef MINMAXSPC_HAS_NAME
    IWString        _name;
#endif

  public:
    Min_Max_Specifier();
#ifdef MINMAXSPC_HAS_NAME
    Min_Max_Specifier(const char * name);
#endif
    Min_Max_Specifier(const T v);
    ~Min_Max_Specifier();

    int ok() const;
    int debug_print(std::ostream &) const;

    int is_set() const { return _is_set;}

    void reset();

    int initialise(const const_IWSubstring &);

    int add(T);    // we override this so we can maintain _is_set
    int add_if_not_already_present(T);  // must override.

    int set_min(T v);
    int set_max(T v);

    int adjust_to_accommodate(const T v);

    int min(T & v) const { return _min_val.value(v);}
    int max(T & v) const { return _max_val.value(v);}

    int matches(T) const;

    int write_msi(std::ostream &, const char *, int = 0) const;

    int write_compact_representation(std::ostream &) const;  // something like '3,4', '>1', '<6'

    Min_Max_Specifier<T> & operator=(const Min_Max_Specifier<T> &);
    bool operator==(const Min_Max_Specifier<T>& rhs) const;
};

template <typename T>  std::ostream &
      operator << (std::ostream &, const Min_Max_Specifier<T> &);

#if (IW_IMPLEMENTATIONS_EXPOSED) || defined(MINMAXSPC_IMPLEMENTATION)

#include <assert.h>

template <typename T>
Min_Max_Specifier<T>::Min_Max_Specifier()
{
  _is_set = 0;

  return;
}

template <typename T>
Min_Max_Specifier<T>::Min_Max_Specifier(const T v)
{
  add(v);

  _is_set = 1;

  return;
}

#ifdef MINMAXSPC_HAS_NAME
template <typename T>
Min_Max_Specifier<T>::Min_Max_Specifier(const char * name)
{
  _name = name;

  return;
}
#endif

template <typename T>
Min_Max_Specifier<T>::~Min_Max_Specifier()
{
//cerr << "Deleting Min_Max_Specifier with " << _number_elements << " items, _is_set?" << _is_set << " _match_any? " << _match_any << '\n';

  if (!ok()) {
    std::cerr << "Deleting bad Min_Max_Specifier\n";
    std::cerr << _number_elements << " items\n";
    std::cerr << _min_val.is_set() << " _min_val.is_set()\n";
    std::cerr << _max_val.is_set() << " _max_val.is_set()\n";
    std::cerr << _match_any << " _match_any\n";

  }

  assert (ok());

  return;
}
template <typename T>
int
Min_Max_Specifier<T>::ok() const
{
  if (_number_elements)
  {
    if (0 == _is_set)
      return 0;
    if (_min_val.is_set())
      return 0;
    if (_max_val.is_set())
      return 0;
    if (_match_any)
      return 0;
  }

  T min;

  if (_min_val.value(min))
  {
    if (0 == _is_set)
      return 0;
    T tmp;
    if (! _max_val.value(tmp))
      return 1;

    if (tmp < min)
      return 0;

    if (_match_any)
      return 0;

    return 1;
  }

  return iwarchive<T>::ok();
}

template <typename T>
int 
Min_Max_Specifier<T>::debug_print(std::ostream & os) const
{
  os << "Details on Min_Max_Specifier<T>::";
  if (_is_set)
    os << " set.";
  else
    os << " not set.";
  if (_match_any)
    os << " matches all values";
  os << '\n';

  if (! ok())
    os << "Warning, OK fails\n";

  T tmp;
  if (_min_val.value(tmp))
    os << "Minimum value is " << tmp << '\n';
  if (_max_val.value(tmp))
    os << "Maximum value is " << tmp << '\n';

  os << "Contains " << _number_elements << " items:";
  for (int i = 0; i < _number_elements; i++) {
    os << " " << _things[i];
  }
  os << '\n';
  
  return 1;
}

template <typename T>
void
Min_Max_Specifier<T>::reset()
{
  iwarchive<T>::resize_keep_storage(0);
  _min_val.unset();
  _max_val.unset();

  _is_set = 0;

  return;
}

template <typename T>
int
Min_Max_Specifier<T>::add(T extra)
{
  if (_min_val.is_set())
    _min_val.unset();

  if (_max_val.is_set())
    _max_val.unset();

  iwarchive<T>::add(extra);

  _is_set = 1;

  _match_any = 0;

  return 1;
}

template <typename T>
int
Min_Max_Specifier<T>::add_if_not_already_present(T extra)
{
  if (iwarchive<T>::contains(extra)) {
    std::cerr << "Min_Max_Specifier:add_if_not_already_present:already present\n";
    return 0;
  }

  return this->add(extra);
}

/*
  _is_set is true if any component has been specified
*/

/*template <typename T>
int
Min_Max_Specifier<T>::is_set() const
{
//assert (ok ());

  if (_match_any)
    return 0;

  if (_number_elements)
    return 1;

  if (_min_val.is_set())
    return 1;

  if (_max_val.is_set())
    return 1;

  return 0;
}*/

template <typename T>
int
Min_Max_Specifier<T>::set_min(T v)
{
  if (_number_elements)
    this->resize(0);

  _min_val.set(v);

  _match_any = 0;

  _is_set = 1;

  return 1;
}

template <typename T>
int
Min_Max_Specifier<T>::set_max(T v)
{
  if (_number_elements) 
    this->resize(0);

  _max_val.set(v);

  _match_any = 0;

  _is_set = 1;

  return 1;
}

/*
  Remember, that we specify either MIN and/or MAX, or values in the
  archive. So, something which is OK with _min_val, never gets
  checked against the archive.
*/

template <typename T>
int
Min_Max_Specifier<T>::matches(const T v) const
{
//assert (ok());

  if (_match_any)
    return 1;

  T tmp;
  int check_archive;

  if (_min_val.value(tmp))
  {
    if (v < tmp)
      return 0;
    check_archive = 0;
  }
  else
    check_archive = 1;

// At this stage, either _min_val was not set, or v is ok wrt _min_val

  if (_max_val.value(tmp))
    return v <= tmp;
  else if (0 == check_archive)
    return 1;
  else if (0 == _number_elements)
    return 0;
  else
    return iwarchive<T>::matches(v);
}

template <typename T>
int
Min_Max_Specifier<T>::write_msi(std::ostream & os, const char * name,
                                int indentation) const
{
  assert(ok());
  assert(os.good());

  IWString ind;
  if (indentation)
    ind.extend(indentation, ' ');

  T tmp;

  if (_min_val.value(tmp))
    os << ind << "(A I min_" << name << " " << tmp << ")\n";

  if (_max_val.value(tmp))
    os << ind << "(A I max_" << name << " " << tmp << ")\n";

  if (0 == _number_elements)
    return os.good();

  os << ind << "(A I " << name << ' ';
  if (1 == _number_elements)
    os << _things[0];
  else
  {
    os << '(';

    for (int i = 0; i < _number_elements; i++)
    {
      if (i)
        os << ' ';
      os << _things[i];
    }
    os << ')';
  }
  os << ")\n";

  return os.good();
}

/*
  We need to adjust anything set to ensure that V will match.
*/

template <typename T>
int
Min_Max_Specifier<T>::adjust_to_accommodate(const T v)
{
  assert (ok());

  if (matches(v))    // nothing to do!
    return 0;

// If neither min nor max are specified, we will need to add this to the archive.

  int check_archive = 1;

// Does not match, something is wrong.

  T minv;
  T maxv;
  if (_min_val.value(minv) && v < minv)
  {
    minv = v;
    _min_val.set(minv);
    check_archive = 0;
    _match_any = 0;
  }

  if (_max_val.value(maxv) && v > maxv)
  {
    maxv = v;
    _max_val.set(maxv);
    _match_any = 0;
    check_archive = 0;
  }

  if (! check_archive)
    return 1;

  if (! this->contains(v))
    add(v);

  _match_any = 0;

  return 1;
}

template <typename T>
Min_Max_Specifier<T> &
Min_Max_Specifier<T>::operator = (const Min_Max_Specifier<T> & other)
{
  assert(ok());
  assert(other.ok());

  _min_val = other._min_val;
  _max_val = other._max_val;
#ifdef MINMAXSPC_HAS_NAME
  _name = other._name;
#endif
  _is_set = other._is_set;

  iwarchive<T>::operator=(other);

  return *this;
}

template <typename T>
bool
Min_Max_Specifier<T>::operator==(const Min_Max_Specifier<T>& rhs) const {
  if (_min_val != rhs._min_val) {
    return false;
  }
  if (_max_val != rhs._max_val) {
    return false;
  }
  return iwarchive<T>::operator==(rhs);
}

/*
  An initialising string can look like

  1      a digit
  2,3,4  a sequence of digits
  >1     a range
  <2     a range

  TODO:ianwatson branch on type_traits to properly handle
  int and float values.
*/

template <typename T>
int
Min_Max_Specifier<T>::initialise(const const_IWSubstring & buffer)
{
  reset();

  if (buffer.starts_with("<="))
  {
    const_IWSubstring mybuffer = buffer;
    mybuffer.remove_leading_chars(2);

    T v;
    if (! mybuffer.numeric_value(v))
    {
      std::cerr << "Min_Max_Specifier::initialise: invalid numeric '" << buffer << "'\n";
      return 0;
    }

    set_max(v);

    return 1;
  }

  if (buffer.starts_with('<'))
  {
    const_IWSubstring mybuffer = buffer;
    mybuffer.remove_leading_chars(1);

    int m;
    if (! mybuffer.numeric_value(m))
    {
      std::cerr << "Min_Max_Specifier::initialise: invalid numeric '" << buffer << "'\n";
      return 0;
    }

    set_max(m - 1);    // works for ints, but not floats. TODO:ianwatson: std::type_traits

    return 1;
  }

  if (buffer.starts_with(">="))
  {
    const_IWSubstring mybuffer = buffer;
    mybuffer.remove_leading_chars(2);

    T v;
    if (! mybuffer.numeric_value(v))
    {
      std::cerr << "Min_Max_Specifier::initialise: invalid numeric '" << buffer << "'\n";
      return 0;
    }

    set_min(v);

    return 1;
  }

  if (buffer.starts_with('>'))
  {
    const_IWSubstring mybuffer = buffer;
    mybuffer.remove_leading_chars(1);

    int m;
    if (! mybuffer.numeric_value(m))
    {
      std::cerr << "Min_Max_Specifier::initialise: invalid numeric '" << buffer << "'\n";
      return 0;
    }

    set_min(m + 1);      // works with integers but not floats

    return 1;
  }

  int i = 0;
  const_IWSubstring token;

  char separator;
  if (buffer.contains(','))
    separator = ',';
  else
    separator = ' ';

  while (buffer.nextword(token, i, separator))
  {
    if (token.starts_with('='))
      token.remove_leading_chars(1);

    int j = token.index('-');
    if (j > 0)    // looks like range specification
    {
      const_IWSubstring r1, r2;
      if (! token.split(r1, '-', r2) || 0 == r1.length() || 0 == r2.length())
      {
        std::cerr << "Min_Max_Specifier::initialise:invalid range specification '" << token << "'\n";
        return 0;
      }

      T vr1, vr2;
      if (! r1.numeric_value(vr1) || ! r2.numeric_value(vr2) || vr1 > vr2)
      {
        std::cerr << "Min_Max_Specifier::initialise:invalid range specification '" << token << "'\n";
        return 0;
      }
      while (vr1 <= vr2)    // theoretically T could be float or double type...
      {
        iwarchive<T>::add(vr1);
        vr1 += static_cast<T>(1);
      }
    }
    else
    {
      T m;

      if (! token.numeric_value(m))
      {
        std::cerr << "Min_Max_Specifier::initialise: invalid numeric '" << buffer << "'\n";
        return 0;
      }

      iwarchive<T>::add(m);
    }
  }

  _is_set = _number_elements;

  return _is_set;
}

template <typename T>
int
Min_Max_Specifier<T>::write_compact_representation(std::ostream & os) const
{
  T tmp;

  if (_min_val.value(tmp))
  {
    os << ">=" << tmp;
    return os.good();
  }

  if (_max_val.value(tmp))
  {
    os << "<=" << tmp;
    return os.good();
  }

  for (int i = 0; i < _number_elements; i++)
  {
    if (i > 0)
      os << ',';

    os << _things[i];
  }

  return os.good();
}

#endif

#if (IW_IMPLEMENTATIONS_EXPOSED) || defined(MINMAXSPC_OP_IMPLEMENTATION)

// This is a separate include file because I cannot get the templates to instantiate

template <class T>
std::ostream &
operator << (std::ostream & os, const Min_Max_Specifier<T> & qq) 
{
  os << "minmax: ";

  T tmp;
  if (qq.min(tmp))
    os << "min " << tmp;
  if (qq.max(tmp))
    os << " max " << tmp;
  
// This next line commented out because of template instantiation problems.
// The way things work now, this forces instantiation of the iwarchive<T>
// class, which then conflicts with the iwarchive<T> instantiation in its
// file.

//os << " " <<  (iwarchive<T> &) qq;

  if (qq.match_any())
    os << " archive matches any value";
  else
  {
    os << " archive matches these values";
    for (int i = 0; i < qq.number_elements(); i++)
      os << " " << qq[i];
  }

  return os;
}

#endif  // IW_IMPLEMENTATIONS_EXPOSED
#endif  //FOUNDATIONAL_IWMISC_MINMAXSPC_H
