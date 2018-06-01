#ifndef IW_ARCHIVE1_H
#define IW_ARCHIVE1_H

#include <iostream>
using std::cerr;
using std::endl;

class msi_object;
class const_IWSubstring;

#include "iwaray.h"

template <typename T>
class iwarchive : public resizable_array<T>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using resizable_array_base<T>::_number_elements;
    using resizable_array_base<T>::_elements_allocated;
    using resizable_array_base<T>::_things;
#endif

  protected:
    int _match_any;

  public:
    iwarchive ();
    iwarchive (const iwarchive<T> &);

    int ok () const;
    int debug_print (std::ostream & = cerr) const;

    int set_match_any (int i);
    int match_any () const { return _match_any;};

    int add_values_from_msi_object (const msi_object *, const char *);
    int write_msi (std::ostream &, const char *, int) const;

    int add (T);      // we need to override the default add ()

    int matches (T) const;
    int matches (const resizable_array<T> &) const;

    iwarchive<T> & operator= (const iwarchive<T> &);

//  Specifications look like 'x' or 'x,y'

    int specification_from_string (const const_IWSubstring &);
};

template <typename T>  std::ostream &
      operator << (std::ostream &, const iwarchive<T> &);

#if (IW_IMPLEMENTATIONS_EXPOSED) || defined(IWARCHIVE_IMPLEMENTATION)

#include "iwstring.h"

template <typename T>
iwarchive<T>::iwarchive ()
{
  _match_any = 1;
}

template <typename T>
iwarchive<T>::iwarchive (const iwarchive<T> & other)
{
  operator= (other);
}

template <typename T>
int
iwarchive<T>::ok () const
{
  if (_match_any && _number_elements)
    return 0;

  return resizable_array<T>::ok ();
}

template <typename T>
int
iwarchive<T>::debug_print (std::ostream & os) const
{
  os << "iwarchive<T>::debug_print:\n";
  if (! ok ())
    os << "Warning, ok fails\n";

  if (_match_any)
    os << "This archive will match any value\n";
  
  for (int i = 0; i < _number_elements; i++)
  {
    os << "Item " << i << " will match " << _things[i] << endl;
  }

  return 1;
}

template <typename T>
int
iwarchive<T>::set_match_any (int i)
{
  if (0 == i)
    return _match_any = i;

  if (_number_elements)    // should be fatal, ok() will fail
    cerr << "iwarchive<T>::set_match_any: warning, " << _number_elements << " values present\n";

  return _match_any = i;
}

template <typename T>
int
iwarchive<T>::add (T i)
{
  _match_any = 0;
  return resizable_array<T>::add (i);
}

template <typename T>
int
iwarchive<T>::matches (T x) const
{
  if (_match_any)
    return 1;

  for (int i = 0; i < _number_elements; i++)
  {
    if (x == _things[i])
      return 1;
  }

  return 0;
}

template <typename T>
int
iwarchive<T>::matches (const resizable_array<T> & x) const
{
  if (_match_any)
    return 1;

  int xn = x.number_elements ();
  for (int i = 0; i < _number_elements; i++)
  {
    const T tmp = _things[i];
    for (int j = 0; j < xn; j++)
    {
      if (x[j] == tmp)
        return 1;
    }
  }

  return 0;
}

template <typename T>
iwarchive<T> &
iwarchive<T>::operator = (const iwarchive<T> & other)
{
  _match_any = other._match_any;

  this->resize (other._elements_allocated);
  for (int i = 0; i < other._number_elements; i++)
  {
    _things[i] = other._things[i];
  }

  _number_elements = other._number_elements;

  return * this;
}

template <typename T>
int
iwarchive<T>::specification_from_string (const const_IWSubstring & s)
{
  int i = 0;
  const_IWSubstring token;

  while (s.nextword (token, i, ','))
  {
    T v;
    if (! token.numeric_value (v))
    {
      cerr << "iwarchive::specification_from_string:invalid numeric '" << token << "'\n";
      return 0;
    }

    resizable_array<T>::add (v);
  }

  return _number_elements;
}

#endif

#if (IW_IMPLEMENTATIONS_EXPOSED) || defined(IWARCHIVE_MSI_IMPLEMENTATION)

#include "msi_object.h"

template <typename T>
int
iwarchive<T>::add_values_from_msi_object (const msi_object * msi,
                                          const char * attribute_name)
{
  int i = 0;
  const msi_attribute * msi_attribute;
  while (NULL != (msi_attribute = msi->attribute (attribute_name, i++)))
  {
    T x;
    if (! msi_attribute->value (x))
    {
      cerr << "iwarchive<T>::add_values_from_msi_object: invalid value " <<
             (*msi_attribute) << endl;
      return 0;
    }

    add (x);
  }

  if (_number_elements)
    _match_any = 0;

  return _number_elements;
}

template <typename T>
int
iwarchive<T>::write_msi (std::ostream & os,
                         const char * attribute_name,
                         int indentation) const
{
  assert (ok ());
  assert (os.good ());

  IWString ind;
  if (indentation)
    ind.extend (indentation, ' ');

  for (int i = 0; i < _number_elements; i++)
  {
    os << ind << "(A I " <<  attribute_name << " " << _things[i] << ")\n";    // only works with T == int
  }

  return os.good ();
}

#endif

#if (IW_IMPLEMENTATIONS_EXPOSED) || defined(IWARCHIVE_OP_IMPLEMENTATION)

template <class T>
std::ostream &
operator << (std::ostream & os, const iwarchive<T> & qq)
{
  os << "iwarchive: ";
  if (qq.match_any ())
    os << "matches any value";
  else if (0 == qq.number_elements ())
    os << "nothing set";
  else
  {
    os << "matches";
    for (int i = 0; i < qq.number_elements (); i++)
      os << " " << qq[i];
  }

  return os;
}

#endif

#endif
