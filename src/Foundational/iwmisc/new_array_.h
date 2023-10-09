/*
  Implementation used by the various new_ functions
*/

#include <algorithm>

#ifndef IW_NEW_ARRAY_H
#define IW_NEW_ARRAY_H

#include "misc.h"

template <typename T>
T *
iw_new_array(int size, T initial_value)
{
  T * rc = new T[size];

  std::fill_n(rc, size, initial_value);
//set_vector (rc, size, initial_value);

  return rc;
}

#endif
