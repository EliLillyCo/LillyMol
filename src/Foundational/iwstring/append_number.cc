#include <stdlib.h>
#include <iomanip>

#include "iwstring.h"

static int power10 [] = {1,                 // 0
                         10,                // 1
                         100,               // 2
                         1000,              // 3
                         10000,             // 4
                         100000,            // 5
                         1000000,           // 6
                         10000000,          // 7
                         100000000,         // 8
                         1000000000};
 
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

  int chars_for_sign;
  if (n < 0)
    chars_for_sign = 1;
  else
    chars_for_sign = 0;

  if (n >= power10[width + chars_for_sign])  // number cannot fit
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
