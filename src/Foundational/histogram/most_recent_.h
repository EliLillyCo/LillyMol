#ifndef MR_IMPLEMENTATION_H
#define MR_IMPLEMENTATION_H

#include <iostream>
#include "most_recent.h"

template <typename T>
void
IWMost_Recent<T>::_default_values()
{
  _things = nullptr;

  _elements_allocated = 0;

  _items_added = 0;

  _next_ptr = 0;

  return;
}

template <typename T>
IWMost_Recent<T>::IWMost_Recent()
{
  _default_values();

  return;
} 

template <typename T>
IWMost_Recent<T>::IWMost_Recent(int initial_size)
{
  _default_values();

  resize(initial_size);

  return;
}

template <typename T>
IWMost_Recent<T>::IWMost_Recent(const IWMost_Recent<T> & rhs)     // copy constructor
{
  _default_values();

  operator = (rhs);

  return;
}

template <typename T>
int
IWMost_Recent<T>::resize(int new_size)
{
  assert (ok ());

  if (0 == _items_added)
  {
    _things = new T[new_size];

    _elements_allocated = new_size;

    return 1;
  }

  assert (NULL == "implement this some time");

  T * new_things = new T[new_size];

  if (_items_added < _elements_allocated)
  {
  }

  delete _things;

  _things = new_things;

  _elements_allocated = new_size;

  return 1;
}

template <typename T>
void
IWMost_Recent<T>::clear()
{
  assert (ok());

  _items_added = 0;
  _next_ptr = 0;

  return;
}

template <typename T>
IWMost_Recent<T> &
IWMost_Recent<T>::operator =(const IWMost_Recent<T> & rhs)
{
  assert (ok());
  assert (rhs.ok());

  clear();

  if (_elements_allocated != rhs._elements_allocated && nullptr != _things)    // we will need to resize ourselves
  {
    delete _things;
    _things = nullptr;
    _elements_allocated = 0;
  }

  if (0 == rhs._elements_allocated)     // rhs is empty
  {
    return *this;
  }

  _things = new T[rhs._elements_allocated];

  _elements_allocated = rhs._elements_allocated;

  for (int i = 0; i < _elements_allocated; i++)
  {
    _things[i] = rhs._things[i];
  }

  _items_added = rhs._items_added;

  _next_ptr = rhs._next_ptr;

  return *this;
}

template <typename T>
int
IWMost_Recent<T>::ok() const
{
  if (0 == _elements_allocated && nullptr == _things)    // probably empty
  {
    return (0 == _items_added && 0 == _next_ptr);
  }

  if (_next_ptr >= _elements_allocated)
    return 0;

  if (0 == _elements_allocated || nullptr == _things)
    return 0;

  return 1;
}

template <typename T>
int
IWMost_Recent<T>::debug_print(ostream & os) const
{
  os << "Most Recent object with " << _elements_allocated << " items allocated, sampled " << _items_added << endl;

  int istop;

  if (_items_added < _elements_allocated)
    istop = _items_added;
  else
    istop = _elements_allocated;

  for (int i = 0; i < istop; i++)
  {
    os << " i = " << i << " value " << _things[i];
    if (i == _next_ptr)
      os << "  <- next";
    os << endl;
  }

  return os.good();
}

template <typename T>
void
IWMost_Recent<T>::extra(T e)
{
  assert (ok ());
  assert (_elements_allocated > 0);

  _items_added++;

  _things[_next_ptr] = e;
  _next_ptr++;
  if (_next_ptr == _elements_allocated)
    _next_ptr = 0;

  return;
}

template <typename T>
int
IWMost_Recent<T>::items_stored() const
{
  if (_items_added <= _elements_allocated)
    return _items_added;

  return _elements_allocated;
}

/*
  One of the operator []'s has been called. What is the real index.
*/

template <typename T>
int
IWMost_Recent<T>::_get_index(int ndx) const
{
  assert (ndx >= 0);

  if (_items_added <= _elements_allocated)
  {
    assert (ndx < _items_added);

    return ndx;
  }

  ndx = _next_ptr + ndx;
  if (ndx >= _elements_allocated)
    ndx -= _elements_allocated;

  return ndx;
}

template <typename T>
typename IWMost_Recent<T>::const_reference
IWMost_Recent<T>::operator [] (int ndx) const
{
  ndx = _get_index(ndx);

  return _things[ndx];
}

template <typename T>
typename IWMost_Recent<T>::reference
IWMost_Recent<T>::operator [] (int ndx)
{
  ndx = _get_index(ndx);

  return _things[ndx];
}

#endif
