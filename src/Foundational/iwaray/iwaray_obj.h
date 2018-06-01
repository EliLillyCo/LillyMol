#ifndef IWARAY_OBJ_H
#define IWARAY_OBJ_H

/*
  This file contains templates for functions which only make sense with true
  objects. This is because a resizable_array_p<Foo> has member functions,
  whereas a resizable_array_p<int> does not.

  Include this file for instantiation of the template where needed.

  NO, this doesn't work. The problem is that T is now an Object *
*/

#include "iwaray.h"

template <typename T>
void
resizable_array_p<T>::each (void (T:: * xx) ())
{
  assert (ok ());
  for (int i = 0; i < _number_elements; i++)
  {
    T * tmp = _things[i];
    xx (tmp);
    _things[i]->xx ();
  }

  return;
}

template <typename T>
int
resizable_array_p<T>::each (int (T:: * xx) ())
{
  assert (ok ());

  int rc = 0;
  for (int i = 0; i < _number_elements; i++)
    rc += _things[i]->xx ();

  return rc;
}

#endif
