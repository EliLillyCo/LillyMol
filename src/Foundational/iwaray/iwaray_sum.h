#ifndef IWARAY_SUM_H
#define IWARAY_SUM_H

#include "iwaray.h"
#include "assert.h"

/*
  Special include file for the sum () function. Only relevant to the
  resizable_array_base<int, float, double> class.
*/

template <typename T>
T
resizable_array_base<T>::sum () const
{
  assert (ok ());

  T rc = T (0);
  for (int i = 0; i < _number_elements; i++)
  {
    rc += _things[i];
  }

  return rc;
}

#endif
