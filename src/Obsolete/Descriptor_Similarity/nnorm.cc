#include <stdlib.h>
#include <math.h>

#include "nnorm.h"

double
NNorm::operator() (float f1, float f2) const
{
  return pow(f1 - f2, _n);
}

double
NNorm::final(double rc,
             int n,
             int valid_pairs) const
{
  return rc;
}
