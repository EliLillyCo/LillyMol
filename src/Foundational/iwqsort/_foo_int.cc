#include <stdlib.h>

#include "foo_.h"

template class Foo<int>;

//Foo<int>::(Foo<int>::*iwqsort_mfn) (const Foo<int> &) const = NULL;
int (Foo<int>::*iwqsort_mfn) = NULL;

template <> int Foo<int>::_static_variable = 0;

//  static int (Foo<T>::*iwqsort_mfn) (const Foo<T> &) const;

int
foo_comparitor_int (const Foo<int> & f1, const Foo<int> & f2)
{
  return foo_comparitor (f1, f2);
}

