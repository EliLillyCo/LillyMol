#include <stdlib.h>

#include "misc.h"

int
Int_Comparator_Larger::operator() (int i1, int i2) const
{
  if (i1 < i2)
    return -1;

  if (i1 > i2)
    return 1;

  return 0;
}

int
Int_Comparator_Smaller::operator() (int i1, int i2) const
{
  if (i1 < i2)
    return 1;

  if (i1 > i2)
    return -1;

  return 0;
}
