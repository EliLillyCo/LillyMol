/*
  Test the C++ function interface
*/

#include <stdlib.h>
#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define IWARAY_EACH_IMPLEMENTATION

#include "iwaray.h"

template <typename T>
class SillyClass
{
  private:
    T _myvalue;

  public:
    SillyClass(int);

    T myvalue () const { return _myvalue;}
    void set_value(T s) { _myvalue = s;}
};

template <typename T>
SillyClass<T>::SillyClass(int i)
{
  _myvalue = static_cast<T>(i);

  return;
}


template <typename T, typename OP>
int
test_remove_items(int n, OP opq)
{
  resizable_array_p<SillyClass<T>> v;

  for (int i = 0; i < n; ++i)
  {
    v.add(new SillyClass<T>(i));
  }

  int rmd = v.remove_items_fn(opq);

  cerr << "Removed " << rmd << " items\n";

  int rc = 1;    // success

  for (int i = 0; i < v.number_elements(); ++i)
  {
    if (opq(v[i]))
    {
      cerr << "Object " << i << " value " << v[i]->myvalue() << " remains\n";
      rc = 0;
    }
  }

  return rc;
}

template <class T>
int
test_template (int n, T rhs)
{
  resizable_array<T> foo;

  foo.resize (n);

  for (int i = 0; i < n; i++)
  {
    foo.add (T(i));
  }

  foo.each (std::bind2nd (std::plus<T>(), rhs));    // does not work, does not update the value in place, fix sometime

  for (int i = 0; i < n; i++)
  {
    cerr << " i = " << i << " new value " << foo[i] << endl;
  }

  return 1;
}

#ifdef __GNUG__
template int test_template (int, int);
template void resizable_array<int>::each<std::binder2nd<std::plus<int> > >(std::binder2nd<std::plus<int> > const&);
#endif

static int
test_functions ()
{
  (void) test_template<int> (10, 7);

  test_remove_items<int>(100, [] (const SillyClass<int> * s) { if (0 == s->myvalue()%5) return 1; return 0;});
  test_remove_items<float>(20, [] (const SillyClass<float> * s) { return s->myvalue() > 5.0;});

  return 1;
}

int
main ()
{
  test_functions ();

  return 0;
}
