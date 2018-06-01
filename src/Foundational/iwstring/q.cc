#include <stdlib.h>
#include <iostream.h>
#include "iwaray.h"

class foo
{
  private:
    int _i;
  public:
    foo (int);
    ~foo ();
};

foo::foo (int i)
{
  _i = i;
  return;
}


foo::~foo ()
{
  cout << "Destroying foo " << _i << "\n";
}

template class resizable_array_p<foo>;   // instantiate the template
template class resizable_array_base<foo *>;

int 
main ()
{
  resizable_array_p<foo> mfoos (3);
//resizable_array_p<int> ints (3);

  int i = 0;
  int j = 1;

//ints.add (&i);
//ints.add (&j);

  foo f1 (1);
  foo f2 (2);
  cerr << "foos declared\n";

  mfoos.add (&f1);
  mfoos.add (&f2);
}
