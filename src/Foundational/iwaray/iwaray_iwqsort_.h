#ifndef IWARAY_IWQSORT_IMPLEMENTATION_H
#define IWARAY_IWQSORT_IMPLEMENTATION_H

#include "iwaray.h"

#include "iwqsort_fo_.h"

template <typename T> template <typename C>
void
resizable_array_base<T>::iwqsort (C & comparator)
{
  if (_number_elements < 2)
    return;

  ::iwqsort (_things, _number_elements, comparator);

  return;
}

#endif
