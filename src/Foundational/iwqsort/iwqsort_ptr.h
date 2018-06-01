#ifndef IWQSORT_PTR_H
#define IWQSORT_PTR_H 1

//template <typename T> void iwqsort (T *, int, int (T::*) (const T *) );
template <typename T> void iwqsort_ptr (T **, int);
template <typename T, typename C> void iwqsort_ptr (T **, int, const C &);

#endif
