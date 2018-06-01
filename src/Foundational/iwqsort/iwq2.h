#ifndef IWQSORT_H
#define IWQSORT_H 1

template <typename T> void iwqsort (T *, int);
template <typename T> void iwqsort (T *, int, int (*) (T &, T &) );
template <typename T> void iwqsort (T *, int, int (*) (const T &, const T &) );

// With a function object comparitor

template <typename T, typename C> void iwqsort (T *, int, C &);

#ifdef RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION

#include "iwaray.h"

template <typename T> template <typename C>
void
resizable_array_base<T>::iwqsort (C & comparator)
{
  if (_number_elements < 2)
    return;

  ::iwqsort (_things, _number_elements, comparator);

  return;
}

#endif

#ifdef IWQSORT_FO_IMPLEMENTATION

#include <assert.h>

//#define DEBUG_IWQSORT

template <typename T>
void 
swap_elements (T & t1, T & t2,
               void * tmp)
{
//cerr << "Swapping " << t1 << " and " << t2 << ", size " << sizeof (T) << endl;
  memcpy (tmp, &t1, sizeof (T));
  memcpy (&t1, &t2, sizeof (T));
  memcpy (&t2, tmp, sizeof (T));

  return;
}

template <typename T, typename C>
void
compare_two_items (T * t,
                   C & comparitor,
                   void * tmp)
{
  int c = comparitor (t[0], t[1]);

  if (c <= 0)
    return;

  swap_elements (t[0], t[1], tmp);

  return;
}

/*
  Scan left until we find something less than pivot
*/

template <typename T, typename C>
void
move_in_from_right (T * t, 
                    int & low,
                    int & high,
                    C & comparitor)
{
  T & pivot = t[low];

  while (high > low)
  {
    int c = comparitor (pivot, t[high]);
#ifdef DEBUG_IWQSORT
    cerr << "Right between " << pivot << " and " << t[high] << " comparison " << c << endl;
#endif

    if (c <= 0)
      high--;
    else if (c > 0)
      break;
  }

  return;
}

/*
  Scan right until we find something greater than pivot
*/

template <typename T, typename C>
void
move_in_from_left (T * t,
                   int & low,
                   int & left,
                   int n,
                   C & comparitor,
                   void * tmp)
{
  T & pivot = t[low];

  while (left < n)
  {
    int c = comparitor (pivot, t[left]);
#ifdef DEBUG_IWQSORT
    cerr << "Left " << left << " between " << pivot << " and " << t[left] << " comparison " << c << endl;
#endif

    if (c > 0)
      left++;
    else if (0 == c)
    {
      low++;
      if (left > low)
        swap_elements (t[low], t[left], tmp);
      left++;
    }
    else
      break;
  }

  return;
}

template <typename T, typename C>
void
iwqsort (T * t, int n,
         C & comparitor,
         void * tmp)
{
#ifdef DEBUG_IWQSORT
  cerr << "On entry n = " << n << endl;
  for (int i = 0; i < n; i++)
  {
    cerr << " i = " << i << " value " << t[i] << endl;
  }
#endif

  if (n < 1)
    return;

  if (2 == n)
  {
    compare_two_items (t, comparitor, tmp);
    return;
  }

  int low = 0;
  int left = 1;
  int right = n - 1;
  while (1)
  {
    move_in_from_left (t, low, left, n, comparitor, tmp);
    move_in_from_right (t, low, right, comparitor);

//  cerr << "Low " << low << " Left " << left << " right " << right << endl;
    if (left < right)
      swap_elements (t[left], t[right], tmp);
    else
      break;
  }

  if (right > low)
  {
    if (low == n - 1)     // all values in this chunk constant
      return;

#ifdef DEBUG_IWQSORT
    cerr << "N = " << n << " before moving, low = " << low << " left " << left << " right " << right << endl;
    for (int i = 0; i < n; i++)
    {
      cerr << " i " << i << " value " << t[i] << endl;
    }
#endif

    for (int i = 0; i < low + 1; i++)
    {
      swap_elements (t[i], t[right - i], tmp);
    }
  }

#ifdef DEBUG_IWQSORT
  cerr << "Starting recursion, low " << low << " left " << left << " right " << right << endl;
  for (int i = 0; i < n; i++)
  {
    cerr << " i " << i << " value " << t[i] << endl;
  }
#endif

  iwqsort (t, right - low, comparitor, tmp);
  iwqsort (t + right + 1, n - right - 1, comparitor, tmp);

  return;
}

template <typename T, typename C>
void
iwqsort (T * t, int n,
         C & comparitor)
{
  assert (n >= 0);

  if (n < 2)
    return;

  unsigned char * tmp = new unsigned char[sizeof (T)];

  iwqsort (t, n, comparitor, tmp);

  delete [] tmp;

  return;
}

#endif

#ifdef IWQSORT_IMPLEMENTATION

#include <assert.h>

#define IWSORT_FN iwqsort
#define IWDEREF(o) o.

//#define DEBUG_IWQSORT

template <typename T>
void 
swap_elements (T * t1, T * t2)
{
  T tmp;
  memcpy (&tmp, t1, sizeof (T));
  memcpy (t1, t2, sizeof (T));
  memcpy (t2, &tmp, sizeof (T));

  return;
}

template <typename T>
void
compare_two_items (T * t,
                   int (*mfn) (T &, T &))
{
  int c = mfn (t[0], t[1]);

  if (c <= 0)
    return;

  swap_elements(t, t+ 1);

  return;
}

/*
  Scan left until we find something less than pivot
*/

template <typename T>
void
move_in_from_right (T * t, 
                    int low,
                    int & high,
                    int (*mfn) (T &, T &))
{
  T & pivot = t[low];

  while (high > low)
  {
    int c = mfn (pivot, t[high]);

//  cerr << "Right between " << IWDEREF(pivot)zvalue () << " and " << IWDEREF(t[high])zvalue () << " comparison " << c << endl;
//  cerr << "Right between " << pivot << " and " << t[high] << " comparison " << c << endl;

    if (c <= 0)
      high--;
    else if (c > 0)
      break;
  }

  return;
}

/*
  Scan right until we find something greater than pivot
*/

template <typename T>
void
move_in_from_left (T * t,
                   int & low,
                   int & left,
                   int n,
                   int (*mfn) (T &, T &))
{
  T & pivot = t[low];

  while (left < n)
  {
    int c = mfn (pivot, t[left]);

//  cerr << "Left " << left << " between " << IWDEREF(pivot)zvalue () << " and " << IWDEREF(t[left])zvalue () << " comparison " << c << endl;

    if (c > 0)
      left++;
    else if (0 == c)
    {
      low++;
      if (left > low)
        swap_elements (t + low, t + left);
      left++;
    }
    else
      break;
  }

  return;
}

template <typename T>
void
IWSORT_FN (T * t, int n,
         int (*mfn) (T &, T &))
{
#ifdef DEBUG_IWQSORT
  cerr << "On entry n = " << n << endl;
  for (int i = 0; i < n; i++)
  {
    cerr << " i = " << i << " value " << t[i] << endl;
  }
#endif
  assert (n >= 0);

  if (n < 2)
    return;

  if (2 == n)
  {
    compare_two_items (t, mfn);
    return;
  }

  int low = 0;
  int left = 1;
  int right = n - 1;
  while (1)
  {
    move_in_from_left (t, low, left, n, mfn);
    move_in_from_right (t, low, right, mfn);

//  cerr << "Low " << low << " Left " << left << " right " << right << endl;
    if (left < right)
      swap_elements (t + left, t + right);
    else
      break;
  }

  if (right > low)
  {
    if (low == n - 1)     // all values in this chunk constant
      return;

#ifdef DEBUG_IWQSORT
    cerr << "N = " << n << " before moving, low = " << low << " left " << left << " right " << right << endl;
    for (int i = 0; i < n; i++)
    {
      cerr << " i " << i << " value " << t[i] << endl;
    }
#endif

    for (int i = 0; i < low + 1; i++)
    {
      swap_elements (t + i, t + right - i);
    }
  }

#ifdef DEBUG_IWQSORT
  cerr << "Starting recursion, low " << low << " left " << left << " right " << right << endl;
  for (int i = 0; i < n; i++)
  {
    cerr << " i " << i << " value " << t[i] << endl;
  }
#endif

  IWSORT_FN (t, right - low, mfn);
  IWSORT_FN (t + right + 1, n - right - 1, mfn);

  return;
}

#endif
#endif
