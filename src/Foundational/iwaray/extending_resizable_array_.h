#ifndef EXTENDING_RESIZABLE_ARRAY_H
#define EXTENDING_RESIZABLE_ARRAY_H

#include "iwaray_.h"

template <typename T>
extending_resizable_array<T>::extending_resizable_array (T ini)
{
  _initialiser = ini;

  return;
}

template <typename T>
T &
extending_resizable_array<T>::operator [] (int i)
{
  assert (i >= 0);

  if (i < _number_elements)
    return _things[i];

  if (i < 10)
    resizable_array<T>::extend (10, _initialiser);
  else if (i < 100)
    resizable_array<T>::extend (100, _initialiser);
  else
    resizable_array<T>::extend (i + 100, _initialiser);

  return _things[i];
}

template <typename T>
int
extending_resizable_array<T>::extend (int new_size)
{
  return resizable_array<T>::extend (new_size, _initialiser);
}

#endif
