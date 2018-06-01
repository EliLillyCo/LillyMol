#include <stdlib.h>

#include "iwhistogram.h"

Histogram_with_Overlapping_Buckets::Histogram_with_Overlapping_Buckets ()
{
  _delta = 0.0;

  return;
}

int
Histogram_with_Overlapping_Buckets::extra (double e)
{
  int b;
  if (! IWHistogram::which_bucket (e, b))
    return 0;

  _count[b]++;
  _nsamples++;

  if (0.0 == _delta)
    return 1;

// Maybe we need to populate some other buckets as well

  int b2;

  if (b < (_nbuckets - 1) && IWHistogram::which_bucket (e + _delta, b2) && b2 != b)
    _count[b2]++;

  if (b > 0  && IWHistogram::which_bucket (e - _delta, b2) && b2 != b)
    _count[b2]++;

  return 1;
}
