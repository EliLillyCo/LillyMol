#ifndef IWFACTORIAL_H
#define IWFACTORIAL_H

template <typename T>
class IW_Factorial
{
  private:
    int _size;
    T * _f;

//  private functions

    void _default_values ();

  public:
    IW_Factorial ();
    IW_Factorial (int);

    ~IW_Factorial ();

    int resize (int);

    int max_index () const { return _size + 1;}

    inline T operator [] (int i) const { return _f[i];}
};

#ifdef IWFACTORIAL_IMPLEMENTATION

#include <stdlib.h>
#include <iostream>
using namespace std;
#include <assert.h>

template <typename T>
void
IW_Factorial<T>::_default_values ()
{
  _size = 0;

  _f = NULL;

  return;
}

template <typename T>
IW_Factorial<T>::IW_Factorial ()
{
  _default_values ();

  return;
}

template <typename T>
IW_Factorial<T>::IW_Factorial (int new_size)
{
  _default_values ();

  (void) resize (new_size);

  return;
}

template <typename T>
IW_Factorial<T>::~IW_Factorial ()
{
  assert (_size >= 0);

  if (NULL != _f)
  {
    delete [] _f;
    _f = NULL;
  }

  _size = -1;

  return;
}

template <typename T>
int
IW_Factorial<T>::resize (int new_size)
{
  if (NULL != _f)
    delete [] _f;

  if (0 == new_size)
  {
    _f = NULL;
    _size = 0;

    return 1;
  }

  assert (new_size > 0);

  _f = new T[new_size + 1];

  if (NULL == _f)
  {
    cerr << "IW_Factorial::resize: cannot allocate " << new_size << " items\n";
    return 0;
  }

  _f[0] = 1;

  for (int i = 1; i <= new_size; i++)
  {
    _f[i] = _f[i - 1] * static_cast<T> (i);
  }

  _size = new_size;

  return 1;
}
#endif

#endif
