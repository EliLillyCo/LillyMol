#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"

int
main ()
{
  resizable_array<int> foo;
  foo.add (1);
  foo.add (2);

  cerr << "Foo is " << foo << endl;

  return 0;
}
