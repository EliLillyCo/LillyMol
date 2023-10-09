#include <iostream>

#include "iwbits.h"

int
main ()
{
  for (unsigned int i = 0; i < 255; i++)
  {
    char c[8];

    for (int j = 0; j < 8; j++)
    {
      if (one_bit_32[32 - j - 1] & i)
      {
        c[8 - j - 1] = '1';
      }
      else
      {
        c[8 - j - 1] = '0';
      }
    }

    cout << i << ' ';
    const int * iptr = reinterpret_cast<const int *> (&c);
    cout << (*iptr) << ' ';
    cout.write (c, 8);
    cout << endl;
  }

  return 0;
}
