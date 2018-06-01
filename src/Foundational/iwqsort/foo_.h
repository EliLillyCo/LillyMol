#ifndef IWFOO_IMPLEMENTATION_H
#define IWFOO_IMPLEMENTATION_H

#include "foo.h"

template <typename T>
Foo<T>::Foo ()
{
  _value = static_cast<T> (0);

  _notused = static_cast<T> (8);   // just some random value

  _static_variable = static_cast<T> (3);

  for (int i = 0; i < FOO_ARRAY_SIZE; i++)
  {
    _foo_array[i] = static_cast<T> (i);
  }

  return;
}

template <typename T>
int
Foo<T>::iwqsortcompare (const Foo & rhs) const
{
  if (_value < rhs._value)
    return -1;

  if (_value > rhs._value)
    return 1;

  return 0;
}

template <typename T>
int
foo_comparitor (const Foo<T> & f1, const Foo<T> & f2)
{
  return f1.iwqsortcompare (f2);
}

/*template <typename T>
int
Foo<T>::iwqsort_mfn (Foo & rhs, Foo & lhs)
{
  return rhs.iwqsortcompare (lhs);
}*/

#endif
