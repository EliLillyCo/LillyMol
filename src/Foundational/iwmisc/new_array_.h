/*
  Implementation used by the various new_ functions
*/

#ifndef IW_NEW_ARRAY_H
#define IW_NEW_ARRAY_H

#include "misc.h"

template <typename T>
T *
iw_new_array (int size, T initial_value)
{
  T * rc = new T[size];

  set_vector (rc, size, initial_value);

  return rc;
}

#endif
