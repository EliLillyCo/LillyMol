#include <stdlib.h>
#include <iostream.h>

class Foo
{
  private:
    int _a1;
    int _a2;
  public:
    Foo ();

    int a1 () const { return _a1;}
    int a2 () const { return _a2;}
};

Foo::Foo ()
{
  _a1 = 1;
  _a2 = 2;
}

static int
test_foo (const Foo & f,
          int (Foo::* fn) () const)
{
  cout << "Value " << (f.*f.fn) () << endl;

  return 1;
}

int
main ()
{
  Foo foo;

  test_foo (foo, &(Foo::a1));
  test_foo (foo, &(Foo::a2));
}
