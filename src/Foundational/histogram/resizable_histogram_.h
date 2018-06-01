#ifndef RSHIST_IMPLEMENT_H
#define RSHIST_IMPLEMENT_H

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <iostream.h>
#include <values.h>

#include "iwhistogram.h"

template <typename T>
Resizable_Histogram<T>::Resizable_Histogram ()
{
  _hard_min = RESIZABLE_HISTOGRAM_HARD_MIN;
  _hard_max = RESIZABLE_HISTOGRAM_HARD_MAX;

  return;
}

template <typename T>
int
Resizable_Histogram<T>::ok () const
{
  if (! IWHistogram<T>::ok ())
    return 0;

  return 1;
}

template <typename T>
int
Resizable_Histogram<T>::extra (T e)
{
  assert (ok ());

  if (e >= _min && e <= _max)
    return IWHistogram<T>::extra (e);

  if (e < _hard_min)
    return 0;

  if (e > _hard_max)
    return 0;

// We need to resize ourselves

  int new_nbuckets;
  int * newcount;

  if (e < _min)
  {
    int extra_buckets = int (double (_min - e) / double (_dx)) + 1;

    new_nbuckets = _nbuckets + extra_buckets;

    cerr << "New value " << e << " need " << extra_buckets << " extra buckets " << _min << ',' << _max << endl;

    newcount = new int[new_nbuckets];

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
    int extra_buckets = int (double (e - _max) / double (_dx)) + 1;

    new_nbuckets = _nbuckets + extra_buckets;

    cerr << "New value " << e << " need " << extra_buckets << " extra buckets " << _min << ',' << _max << endl;

    newcount = new int[new_nbuckets];

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

  return IWHistogram<T>::extra (e);
}

#endif
