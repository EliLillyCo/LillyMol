#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <limits>

#include "iwhistogram.h"

Resizable_Histogram::Resizable_Histogram ()
{
  _hard_min = std::numeric_limits<float>::min();
  _hard_max = std::numeric_limits<float>::max();

  return;
}

int
Resizable_Histogram::ok () const
{
  if (! IWHistogram::ok ())
    return 0;

  return 1;
}

int
Resizable_Histogram::extra (double e)
{
  assert (ok ());

  if (e >= _min && e <= _max)
    return IWHistogram::extra (e);

  if (e < _hard_min)
    return 0;

  if (e > _hard_max)
    return 0;

// We need to resize ourselves

  int new_nbuckets;
  unsigned int * newcount;

  if (e < _min)
  {
    int extra_buckets = int ((_min - e) / _dx) + 1;

    new_nbuckets = _nbuckets + extra_buckets;

#ifdef DEBUG_EXTRA
    cerr << "New value " << e << " need " << extra_buckets << " extra buckets " << _min << ',' << _max << endl;
#endif

    newcount = new unsigned int[new_nbuckets];

    for (int i = 0; i < extra_buckets; i++)
    {
      newcount[i] = 0;
    }
    for (int i = 0; i < _nbuckets; i++)
    {
      newcount[i + extra_buckets] = _count[i];
    }
    _min = _min - extra_buckets * _dx;
  }
  else
  {
    int extra_buckets = int ((e - _max) / _dx) + 1;

    new_nbuckets = _nbuckets + extra_buckets;

#ifdef DEBUG_EXTRA
    cerr << "New value " << e << " need " << extra_buckets << " extra buckets " << _min << ',' << _max << endl;
#endif

    newcount = new unsigned int[new_nbuckets];

    for (int i = 0; i < _nbuckets; i++)
    {
      newcount[i] = _count[i];
    }
    for (int i = _nbuckets; i < new_nbuckets; i++)
    {
      newcount[i] = 0;
    }
    _max = _max + extra_buckets * _dx;
  }

  delete _count;

  _count = newcount;

  _nbuckets = new_nbuckets;

#ifdef DEBUG_EXTRA
  assert (int (rint ((_max - _min) / _dx)) + 1 == _nbuckets);
  cerr << "Modified: " << _min << " to " << _max << ", buckets now " << _nbuckets << endl;
#endif

  if (0 ==  IWHistogram::extra (e))
  {
    cerr << "Yipes, the histogram didn't take " << e << endl;
    return 0;
  }

  return 1;
}

