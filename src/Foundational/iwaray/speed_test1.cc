/*
  Want to know if there is any difference between add (T) and add (const T &)
*/

#include <stdlib.h>

#include "iwaray.h"

int
main (int argc, char ** argv)
{
  resizable_array<int> p;

  for (int i = 0; i < 600000; i++)
  {
    for (int j = 0; j < 1000; j++)
    {
      p.add (j);
    }

    p.resize_keep_storage (0);
  }

  return 0;
}
