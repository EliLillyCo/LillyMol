#include <stdlib.h>

#include "foo_.h"

template class Foo<float>;

template <> float Foo<float>::_static_variable = 0.0f;

int
foo_comparitor_float (const Foo<float> & f1, const Foo<float> & f2)
{
  return foo_comparitor (f1, f2);
}
