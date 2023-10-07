#ifndef IWQSORT_H
#define IWQSORT_H 1

#include <string.h>

template <typename T> void iwqsort (T *, int);
template <typename T> void iwqsort (T *, int, int (*) (T &, T &) );
template <typename T> void iwqsort (T *, int, int (*) (const T &, const T &) );

// With a function object comparitor

template <typename T, typename C> void iwqsort (T *, int, C &);

#ifdef RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION

#include "Foundational/iwaray/iwaray.h"

template <typename T> template <typename C>
void
resizable_array_base<T>::iwqsort (C & comparator)
{
  if (_number_elements < 2)
    return;

  ::iwqsort (_things, _number_elements, comparator);

  return;
}

template <typename T> template <typename C>
void
resizable_array_base<T>::iwqsort_lambda(C comparator)
{
  if (_number_elements < 2)
    return;

  ::iwqsort (_things, _number_elements, comparator);

  return;
}

#endif

#ifdef IWQSORT_FO_IMPLEMENTATION

#include <string.h>

#include <assert.h>

//#define DEBUG_IWQSORT

// This does not work. Tried various things to make this adapt
// to the compiler name. Got close, but not quite.
// https://stackoverflow.com/questions/3030099/pragma-in-define-macro
// seemed promising. For now the code below is hard coded for gcc.
#ifdef __clang__
#define COMPILING_ME clang
#else
#define COMPILING_ME gcc
#endif

template <typename T>
void 
swap_elements (T & t1, T & t2,
               void * tmp)
{
//cerr << "Swapping " << t1 << " and " << t2 << ", size " << sizeof (T) << endl;
#ifdef IWQS_USE_OPERATOR
  *tmp = t1;
  t1 = t2;
  t2 = *tmp;
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
  ::memcpy (tmp, &t1, sizeof (T));
  ::memcpy (&t1, &t2, sizeof (T));
  ::memcpy (&t2, tmp, sizeof (T));
#pragma GCC diagnostic pop
#endif

  return;
}
#undef COMPILING_ME

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

    int nmove;
    if (low + low >= right)
      nmove = right - low - 1;
    else
      nmove = low;

    for (int i = 0; i <= nmove ; i++)
    {
//#define CHECK_USELESS_MOVES
#ifdef CHECK_USELESS_MOVES
      if (i >= (right-i))
        cerr << "useless swap, i = " << i << " low " << low << " right " << right << endl;
      else
#endif

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
  ::memcpy (&tmp, t1, sizeof (T));
  ::memcpy (t1, t2, sizeof (T));
  ::memcpy (t2, &tmp, sizeof (T));

  return;
}

template <typename T>
void
#ifdef IWQS_USER_FUNCTION_SUPPLIED
compare_two_items (T * t,
                   int (*mfn) (T &, T &))
{
  int c = mfn (t[0], t[1]);

  if (c <= 0)
    return;

  swap_elements(t, t+ 1);

  return;
}
#elif defined(IWQS_USE_OPERATOR)
compare_two_items (T * t)
{
  if (t[0] <= t[1])
    return;

  swap_elements(t, t + 1);

  return;
}
#else
compare_two_items (T * t)
{
  int c = IWDEREF(t[0])iwqsortcompare (t[1]);

  if (c <= 0)
    return;

  swap_elements (t, t + 1);

  return;
}
#endif

/*
  Scan left until we find something less than pivot
*/

template <typename T>
void
#ifdef IWQS_USER_FUNCTION_SUPPLIED
move_in_from_right (T * t, 
                    int low,
                    int & high,
                    int (*mfn) (T &, T &))
#else
move_in_from_right (T * t, 
                    int low,
                    int & high)
#endif
{
  T & pivot = t[low];

  while (high > low)
  {
#ifdef IWQS_USER_FUNCTION_SUPPLIED
    int c = mfn (pivot, t[high]);
#elif defined(IWQS_USE_OPERATOR)
    int c;
    if (pivot <= t[high])
      c = 0;
    else
      c = 1;
#else
    int c = pivot.iwqsortcompare (t[high]);
#endif
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
#ifdef IWQS_USER_FUNCTION_SUPPLIED
move_in_from_left (T * t,
                   int & low,
                   int & left,
                   int n,
                   int (*mfn) (T &, T &))
#else
move_in_from_left (T * t,
                   int & low,
                   int & left,
                   int n)
#endif
{
  T & pivot = t[low];

  while (left < n)
  {
#ifdef IWQS_USER_FUNCTION_SUPPLIED
    int c = mfn (pivot, t[left]);
#elif defined(IWQS_USE_OPERATOR)
    int c;
    if (pivot == t[left])
      c = 0;
    else if (pivot < t[left])
      c = -1;
    else
      c = 1;
#else
    int c = pivot.iwqsortcompare (t[left]);
#endif
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
#ifdef IWQS_USER_FUNCTION_SUPPLIED
IWSORT_FN (T * t, int n,
         int (*mfn) (T &, T &))
#else
IWSORT_FN (T * t, int n)
#endif
{
#ifdef DEBUG_IWQSORT
  cerr << "On entry n = " << n << endl;
  for (int i = 0; i < n; i++)
  {
#ifdef IWQS_USE_OPERATOR
    cerr << " i = " << i << " value " << t[i] << endl;
#else
    cerr << " i = " << i << " value " << IWDEREF(t[i])zvalue () << endl;
#endif
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
  while (1)
  {
#ifdef IWQS_USER_FUNCTION_SUPPLIED
    move_in_from_left (t, low, left, n, mfn);
    move_in_from_right (t, low, right, mfn);
#else
    move_in_from_left (t, low, left, n);
    move_in_from_right (t, low, right);
#endif

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
#ifdef IWQS_USE_OPERATOR
      cerr << " i " << i << " value " << t[i] << endl;
#else
      cerr << " i " << i << " value " << IWDEREF(t[i])zvalue () << endl;
#endif
    }
#endif

    int nmove;
    if (low + low >= right)
      nmove = right - low - 1;
    else
      nmove = low;


    for (int i = 0; i <= nmove; i++)
    {
#ifdef CHECK_USELESS_MOVES
      if (i >= right - i)
        cerr << "Useless swap, i = " << i << " low " << low << " right " << right << endl;
      else
#endif

      swap_elements (t + i, t + right - i);
    }
  }

#ifdef DEBUG_IWQSORT
  cerr << "Starting recursion, low " << low << " left " << left << " right " << right << endl;
  for (int i = 0; i < n; i++)
  {
#ifdef IWQS_USE_OPERATOR
    cerr << " i " << i << " value " << t[i] << endl;
#else
    cerr << " i " << i << " value " << IWDEREF(t[i])zvalue () << endl;
#endif
  }
#endif

#ifdef IWQS_USER_FUNCTION_SUPPLIED
  IWSORT_FN (t, right - low, mfn);
  IWSORT_FN (t + right + 1, n - right - 1, mfn);
#else
  IWSORT_FN (t, right - low);
  IWSORT_FN (t + right + 1, n - right - 1);
#endif

  return;
}

#endif
#endif
