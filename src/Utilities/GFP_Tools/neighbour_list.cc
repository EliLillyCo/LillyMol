#include <stdlib.h>

//#define NEIGHBOUR_LIST_IMPLEMENTATION

#include "neighbour_list.h"

int
SID_Comparator::operator() (const Smiles_ID_Dist * sid1, const Smiles_ID_Dist * sid2) const
{
  float d1 = sid1->distance();
  float d2 = sid2->distance();

  if (d1 < d2)
    return -1;

  if (d1 > d2)
    return 1;

  return 0;
}
