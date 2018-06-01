/*

*/

#include <stdlib.h>

#include "iwbits.h"

using std::cout;

int
main ()
{
  cout << "static int * first_bit_set = {\n";
  cout <<   "       0,     /* this is wrong, 0 has no bits set\n";
  for (unsigned int i = 0; i < 255; i++)
  {
    unsigned int j = i;

    for (int k = 0; k < 8; k++)
    {
      if (one_bit_8[k] & j)
      {
        cout << "      " << k << ",  /* " << i << " */\n";
        break;
      }
    }
  }

  cout << "};\n";

  return 0;
}
