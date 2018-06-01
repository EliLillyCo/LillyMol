#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"

class Foo
{
  private:
    int _v;

  public:
    Foo(int v) : _v(v) {}
    ~Foo();
};

Foo::~Foo()
{
  cerr << "Deleting foo " << _v << endl;
}

int
main()
{
  resizable_array_p<Foo> foos;

  for (int i = 0; i < 10; i++)
  {
    Foo * f = new Foo(i);
    foos.add(f);
  }

  resizable_array_base<Foo *> & b = foos;
  b.resize_keep_storage(0);
  return 1;
}

