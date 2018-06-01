#ifndef IWARAY_OPERATORS_H
#define IWARAY_OPERATORS_H

#ifndef IWARAY_H
#include "iwaray.h"
#endif

// These operator << functions should be in iwaray_.h, but I cannot get them
// instantiated, so there are here.

template <typename T>
std::ostream &
operator << (std::ostream & os, const resizable_array<T> & qq)
{
  int items = qq.number_elements ();
  os << "Array contains " << items << " items:";
  for (int i = 0; i < items; i++)
  {
    os << ' ' << qq[i];
  }

  return os;
}

template <typename T>
std::ostream &
operator << (std::ostream & os, const resizable_array_p<T> & qq)
{
  int items = qq.number_elements ();
  os << "Pointer array contains " << items << " items\n";
  for (int i = 0; i < items; i++)
  {
    const T * t = qq[i];
    os << i << " " << *t << "\n";
  }

  return os;
}

#endif
