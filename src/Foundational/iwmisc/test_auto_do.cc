/*
  Tester for auto_do
*/

#include <stdlib.h>
#include <iostream>
#include <functional>

#include "iw_auto_do.h"

template <typename T> 
class IWChange
{
  private:
    const T & _inc;
    binary_function<T, T, T> _op;
  public: 
    IWChange (const T & t, binary_function<T, T, T> & c) : _inc (t), _op (c) {}

    void operator () (T & t) { t = _op (t, _inc);}
};

int
main (int argc, char ** argv)
{
  int j = 0;

  typedef IWChange<int> IWadd;

  typedef plus<int> Plus;

//IWadd inc (2);

  IWChange<int> inc (2, Plus );

  for (int i = 0; i < 10; i++)
  {
    IW_Auto_Do<int, IWadd> ado (j, inc);
  }

  cerr << "At end j = " << j << endl;

  return 0;
}
