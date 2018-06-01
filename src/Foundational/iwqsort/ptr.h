#ifndef IWQSORT_PTR_IMPLEMENTATION_H
#define IWQSORT_PTR_IMPLEMENTATION_H


#include "iwqsort_ptr.h"

//#define DEBUG_IWQSORT

template <typename T>
void 
swap_elements (T * t1, T * t2)
{
  T * tmp;
  memcpy (&tmp, t1, sizeof (T *));
  memcpy (t1, t2, sizeof (T *));
  memcpy (t2, &tmp, sizeof (T *));

  return;
}

template <typename T>
#ifdef IWQS_USER_FUNCTION_SUPPLIED
void
compare_two_items (T ** t,
                   int (*mfn) (T &, T &))
{
  int c = mfn (*(t[0]), *(t[1]));
#else
void
compare_two_items (T ** t)
{
  int c = t[0]->iwqsortcompare (*(t[1]));
#endif

  if (c <= 0)
    return;

  swap_elements (t, t + 1);

  return;
}

/*
  Scan left until we find something less than pivot
*/

template <typename T>
void
#ifdef IWQS_USER_FUNCTION_SUPPLIED
move_in_from_right (T ** t, 
                    int & low,
                    int & high,
                    int (*mfn) (T &, T &))
#else
move_in_from_right (T ** t, 
                    int & low,
                    int & high)
#endif
{
  T & pivot = *(t[low]);

  while (high > low)
  {
#ifdef IWQS_USER_FUNCTION_SUPPLIED
    int c = mfn (pivot, *(t[high]));
#else
    int c = pivot.iwqsortcompare (*(t[high]));
#endif
//  cerr << "Right between " << pivot->zvalue () << " and " << t[high]->zvalue () << " comparison " << c << endl;

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
#ifdef IWQS_USER_FUNCTION_SUPPLIED
move_in_from_left (T ** t,
                   int & low,
                   int & left,
                   int n,
                   int (*mfn) (T &, T &))
#else
move_in_from_left (T ** t,
                   int & low,
                   int & left,
                   int n)
#endif
{
  T & pivot = *(t[low]);

  while (left < n)
  {
#ifdef IWQS_USER_FUNCTION_SUPPLIED
    int c = mfn (pivot, *(t[left]));
#else
    int c = pivot.iwqsortcompare (*(t[left]));
#endif
//  cerr << "Left between " << pivot->zvalue () << " and " << t[left]->zvalue () << " comparison " << c << endl;

    if (c > 0)
      left++;
    else if (0 == c)
    {
      low++;
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
#ifdef IWQS_USER_FUNCTION_SUPPLIED
iwqsort_ptr (T ** t, int n)
         int (*mfn) (T &, T &))
#else
iwqsort_ptr (T ** t, int n)
#endif
{
#ifdef DEBUG_IWQSORT
  cerr << "On entry n = " << n << endl;
  for (int i = 0; i < n; i++)
  {
    cerr << " i = " << i << " value " << t[i]->zvalue () << endl;
  }
#endif
  assert (n >= 0);

  if (n < 2)
    return;

  if (2 == n)
  {
#ifdef IWQS_USER_FUNCTION_SUPPLIED
    compare_two_items (t, mfn);
#else
    compare_two_items (t);
#endif
    return;
  }

  int low = 0;
  int left = 1;
  int right = n - 1;
  while (left < right)
  {
#ifdef IWQS_USER_FUNCTION_SUPPLIED
    move_in_from_left (t, low, left, n, mfn);
    move_in_from_right (t, low, right, mfn);
#else
    move_in_from_left (t, low, left, n);
    move_in_from_right (t, low, right);
#endif

//  cerr << "Left " << left << " right " << right << endl;
    if (left < right)
      swap_elements (t + left, t + right);
    else
      break;
  }

  if (right > 0)
  {
    if (low == right)     // all values in this chunk constant
      return;

#ifdef DEBUG_IWQSORT
    cerr << "N = " << n << " before moving, low = " << low << " left " << left << " right " << right << endl;
    for (int i = 0; i < n; i++)
    {
      cerr << " i " << i << " value " << t[i]->zvalue () << endl;
    }
#endif

    for (int i = 0; i < low + 1; i++)
    {
      swap_elements (t + i, t + right - i);
    }

#ifdef DEBUG_IWQSORT
    cerr << "N = " << n << " after moving, low = " << low << " left " << left << " right " << right << endl;
    for (int i = 0; i < n; i++)
    {
      cerr << " i " << i << " value " << t[i]->zvalue () << endl;
    }
#endif
  }

#ifdef DEBUG_IWQSORT
  cerr << "Starting recursion, low " << low << " left " << left << " right " << right << endl;
  for (int i = 0; i < n; i++)
  {
    cerr << " i " << i << " value " << t[i]->zvalue () << endl;
  }
#endif

#ifdef IWQS_USER_FUNCTION_SUPPLIED
  iwqsort_ptr (t, right - low, mfn);
  iwqsort_ptr (t + right + 1, n - right - 1, mfn);
#else
  iwqsort_ptr (t, right - low);
  iwqsort_ptr (t + right + 1, n - right - 1);
#endif

  return;
}

#endif
