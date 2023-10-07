#ifndef IWARAYQ_IMPLEMENTATION_H
#define IWARAYQ_IMPLEMENTATION_H

#include <stdlib.h>
#include <assert.h>

#include "iwaray.h"

template <typename T>
iwaray<T>::iwaray ()
{
  _number_elements = 0;
  _things = nullptr;

  return;
}

template <typename T>
iwaray<T>::iwaray (int qsize)
{
  assert (qsize > 0);

  _things = new T [qsize];

  _number_elements = qsize;
}

template <typename T>
iwaray<T>::~iwaray ()
{
  if (nullptr != _things)
    delete [] _things;

  _number_elements = -1;

  return;
}

template <typename T>
int
iwaray<T>::resize (int qsize)
{
  if (_number_elements)
    delete [] _things;

  _things = new T[qsize];
  _number_elements = qsize;

  return 1;
}

template <typename T>
const T &
iwaray<T>::operator [] (int i) const
{
  assert (i >= 0 && i < _number_elements);

  return _things[i];
}

template <typename T>
T &
iwaray<T>::operator [] (int i)
{
  assert (i >= 0 && i < _number_elements);

  return _things[i];
}

#ifdef IWARAYQ_ASSIGN

template <typename T>
iwaray<T> &
iwaray<T>::operator = (T rhs)
{
  for (int i = 0; i < _number_elements; i++)
  {
    _things[i] = rhs;
  }

  return *this;
}

#endif

#endif
