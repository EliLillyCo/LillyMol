#ifndef IWARAY_NEXT_IMPLEMENTATION_H
#define IWARAY_NEXT_IMPLEMENTATION_H

#include "iwaray.h"

template <typename T> template <typename F>
int
resizable_array_p<T>::next (F & f,
                            typename resizable_array_base<T>::next_iterator & ndx) const
{
  if (ndx >= _number_elements)
    return 0;

  int i;

  if (ndx < 0)
    i = 0;
  else
    i = ndx + 1;
    
  while (i < _number_elements)
  {
    if (f (*(_things[i])))
    {
      ndx = i;
      return 1;
    }

    i++;
  }

  ndx = _number_elements;

  return 0;
}

template <typename T> template <typename F>
int
resizable_array_p<T>::next (F & f,
                            typename resizable_array_base<T>::next_iterator & ndx)
{
  if (ndx >= _number_elements)
    return 0;

  int i;

  if (ndx < 0)
    i = 0;
  else
    i = ndx + 1;
    
  while (i < _number_elements)
  {
    if (f (*(_things[i])))
    {
      ndx = i;
      return 1;
    }

    i++;
  }

  ndx = _number_elements;

  return 0;
}

template <typename T> template <typename F>
int
resizable_array_p<T>::next (const F & f,
                            typename resizable_array_base<T>::next_iterator & ndx) const
{
  if (ndx >= _number_elements)
    return 0;

  int i;

  if (ndx < 0)
    i = 0;
  else
    i = ndx + 1;
    
  while (i < _number_elements)
  {
    if (f (*(_things[i])))
    {
      ndx = i;
      return 1;
    }

    i++;
  }

  ndx = _number_elements;

  return 0;
}

template <typename T> template <typename F>
int
resizable_array_p<T>::next (const F & f,
                            typename resizable_array_base<T>::next_iterator & ndx)
{
  if (ndx >= _number_elements)
    return 0;

  int i;

  if (ndx < 0)
    i = 0;
  else
    i = ndx + 1;
    
  while (i < _number_elements)
  {
    if (f (*(_things[i])))
    {
      ndx = i;
      return 1;
    }

    i++;
  }

  ndx = _number_elements;

  return 0;
}

template <typename T> template <typename F>
void
resizable_array_p<T>::each (F & f) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    f (*(_things[i]));
  }

  return;
}

template <typename T> template <typename F>
void
resizable_array_p<T>::each (const F & f) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    f (*(_things[i]));
  }

  return;
}
template <typename T> template <typename F>
void
resizable_array_p<T>::each (F & f)
{
  for (int i = 0; i < _number_elements; i++)
  {
    f (*(_things[i]));
  }

  return;
}

template <typename T> template <typename F>
void
resizable_array_p<T>::each (const F & f)
{
  for (int i = 0; i < _number_elements; i++)
  {
    f (*(_things[i]));
  }

  return;
}

template <typename T> template <typename F>
int
resizable_array<T>::next (F & f,
                            typename resizable_array_base<T>::next_iterator & ndx)
{
  if (ndx >= _number_elements)
    return 0;

  int i;

  if (ndx < 0)
    i = 0;
  else
    i = ndx + 1;
    
  while (i < _number_elements)
  {
    if (f (_things[i]))
    {
      ndx = i;
      return 1;
    }

    i++;
  }

  ndx = _number_elements;

  return 0;
}

template <typename T> template <typename F>
int
resizable_array<T>::next (F & f,
                            typename resizable_array_base<T>::next_iterator & ndx) const
{
  if (ndx >= _number_elements)
    return 0;

  int i;

  if (ndx < 0)
    i = 0;
  else
    i = ndx + 1;
    
  while (i < _number_elements)
  {
    if (f (_things[i]))
    {
      ndx = i;
      return 1;
    }

    i++;
  }

  ndx = _number_elements;

  return 0;
}

template <typename T> template <typename F>
int
resizable_array<T>::next (const F & f,
                            typename resizable_array_base<T>::next_iterator & ndx)
{
  if (ndx >= _number_elements)
    return 0;

  int i;

  if (ndx < 0)
    i = 0;
  else
    i = ndx + 1;
    
  while (i < _number_elements)
  {
    if (f (_things[i]))
    {
      ndx = i;
      return 1;
    }

    i++;
  }

  ndx = _number_elements;

  return 0;
}

template <typename T> template <typename F>
int
resizable_array<T>::next (const F & f,
                            typename resizable_array_base<T>::next_iterator & ndx) const
{
  if (ndx >= _number_elements)
    return 0;

  int i;

  if (ndx < 0)
    i = 0;
  else
    i = ndx + 1;
    
  while (i < _number_elements)
  {
    if (f (_things[i]))
    {
      ndx = i;
      return 1;
    }

    i++;
  }

  ndx = _number_elements;

  return 0;
}

template <typename T> template <typename F>
void
resizable_array<T>::each (F & f) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    f (_things[i]);
  }

  return;
}

template <typename T> template <typename F>
void
resizable_array<T>::each (const F & f) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    f (_things[i]);
  }

  return;
}
template <typename T> template <typename F>
void
resizable_array<T>::each (F & f)
{
  for (int i = 0; i < _number_elements; i++)
  {
    f (_things[i]);
  }

  return;
}

template <typename T> template <typename F>
void
resizable_array<T>::each (const F & f)
{
  for (int i = 0; i < _number_elements; i++)
  {
    f (_things[i]);
  }

  return;
}

#endif
