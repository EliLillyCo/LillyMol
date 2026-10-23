#include <stdlib.h>
#include <iomanip>

#include "iwstring.h"

template <>
int
append_number (IWString & os,
               int n,
               int width)
{
  if (width <= 0)
  {
    os.append_number(n);
    return 1;
  }

  os.make_room_for_extra_items(width);

  IWString tmp;
  tmp.append_number(n);

  int nspaces = width - tmp.length();

  for (int i = 0; i < nspaces; i++)
  {
    os.add(' ');
  }

  os << tmp;

  return os.length();
}

template <>
int
append_number (std::ostream & os,
               int n,
               int width)
{
  os << std::setw(width) << n;

  return 1;
}
