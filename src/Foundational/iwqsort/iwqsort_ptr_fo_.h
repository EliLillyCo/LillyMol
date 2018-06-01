#ifndef IWQSORT_PFNO_IMPLEMENTATION_H
#define IWQSORT_PFNO_IMPLEMENTATION_H

#include <assert.h>

//#define DEBUG_IWQSORT

template <typename T>
void 
swap_elements (T *& t1, T *& t2)
{
//cerr << "Swapping " << t1 << " and " << t2 << ", size " << sizeof (T) << endl;
  T * tmp = t1;
  tmp = t1;
  t1 = t2;
  t2 = tmp;

  return;
}

template <typename T, typename C>
void
compare_two_items (T ** t,
                   C & comparitor)
{
  int c = comparitor (t[0], t[1]);

  if (c <= 0)
    return;

  swap_elements (t[0], t[1]);

  return;
}

/*
  Scan left until we find something less than pivot
*/

template <typename T, typename C>
void
move_in_from_right (T ** t, 
                    int & low,
                    int & high,
                    C & comparitor)
{
  T * pivot = t[low];

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
move_in_from_left (T ** t,
                   int & low,
                   int & left,
                   int n,
                   C & comparitor)
{
  T * pivot = t[low];

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
        swap_elements (t[low], t[left]);
      left++;
    }
    else
      break;
  }

  return;
}

template <typename T, typename C>
void
iwqsort_ptr_fo (T ** t, int n,
         C & comparitor)
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
    compare_two_items (t, comparitor);
    return;
  }

  int low = 0;
  int left = 1;
  int right = n - 1;
  while (1)
  {
    move_in_from_left (t, low, left, n, comparitor);
    move_in_from_right (t, low, right, comparitor);

//  cerr << "Low " << low << " Left " << left << " right " << right << endl;
    if (left < right)
      swap_elements (t[left], t[right]);
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
      swap_elements (t[i], t[right - i]);
    }
  }

#ifdef DEBUG_IWQSORT
  cerr << "Starting recursion, low " << low << " left " << left << " right " << right << endl;
  for (int i = 0; i < n; i++)
  {
    cerr << " i " << i << " value " << t[i] << endl;
  }
#endif

  iwqsort_ptr_fo (t, right - low, comparitor);
  iwqsort_ptr_fo (t + right + 1, n - right - 1, comparitor);

  return;
}

template <typename T, typename C>
void
iwqsort_ptr (T ** t, int n,
         const C & comparitor)
{
  assert (n >= 0);

  if (n < 2)
    return;

  iwqsort_ptr_fo (t, n, comparitor);

  return;
}

#endif
