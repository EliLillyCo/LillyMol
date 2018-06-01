#ifndef IWQSORT_MBMFN_H
#define IWQSORT_MBMFN_H 1

/*
  Sorting function where we use the member function

  iwqsortcompare

  Will only work with classes that have that specific member function
*/

#include <assert.h>

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
compare_two_items (T * t)
{
  int c = t[0].iwqsortcompare (t[1]);

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
move_in_from_right (T * t, 
                    int low,
                    int & high)
{
  T & pivot = t[low];

  while (high > low)
  {
    int c = pivot.iwqsortcompare (t[high]);

//  cerr << "Right between " << pivot.zvalue () << " and " << t[high].zvalue () << " comparison " << c << endl;

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
                   int n)
{
  T & pivot = t[low];

  while (left < n)
  {
    int c = pivot.iwqsortcompare (t[left]);

//  cerr << "Left " << left << " between " << pivot.zvalue () << " and " << t[left].zvalue () << " comparison " << c << endl;

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
iwqsort (T * t, int n)
{
#ifdef DEBUG_IWQSORT
  cerr << "On entry n = " << n << endl;
  for (int i = 0; i < n; i++)
  {
    cerr << " i = " << i << " value " << IWDEREF(t[i])zvalue () << endl;
  }
#endif
  assert (n >= 0);

  if (n < 2)
    return;

  if (2 == n)
  {
    compare_two_items (t);

    return;
  }

  int low = 0;
  int left = 1;
  int right = n - 1;
  while (1)
  {
    move_in_from_left (t, low, left, n);
    move_in_from_right (t, low, right);

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
      cerr << " i " << i << " value " << IWDEREF(t[i])zvalue () << endl;
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
    cerr << " i " << i << " value " << t[i].zvalue () << endl;
  }
#endif

  iwqsort (t, right - low);
  iwqsort (t + right + 1, n - right - 1);

  return;
}

#endif
