#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;

#include "iwhistogram.h"

Fuzzy_Histogram::Fuzzy_Histogram ()
{
   _bucket_minus = 0;
   _bucket_plus = 0;
}

int
Fuzzy_Histogram::extra (double e)
{
  if (_min == _max || 0.0 == _dx)
  {
    cerr << "IWHistogram::extra: not initialised\n";
    return 0;
  }

  if (e < _min)
  {
    return 0;
  }

  if (e > _max)
  {
    return 0;
  }

  int bucket = int (rint ((e - _min) / _dx));

  _count[bucket]++;
  _nsamples++;

  if (_bucket_minus > 0)
  {
    int istop = bucket - _bucket_minus;
    if (istop < 0)
      istop = 0;

    for (int i = bucket - 1; i >= istop; i--)
    {
      _count[i]++;
    }
  }

  if (_bucket_plus > 0)
  {
    int istop = bucket + _bucket_plus;

    if (istop > _nbuckets)
      istop = _nbuckets;

    for (int i = bucket + 1; i < istop; i++)
    {
      _count[i]++;
    }
  }

  return 1;
}
