#include <stdlib.h>
#include <iostream>
using namespace std;

#define SET_OR_UNSET_IMPLEMENTATION
#include "set_or_unset.h"

template class Set_or_Unset<float>;

// Aug 99, cannot get template instantiation to work

template ostream & operator << (ostream &, const Set_or_Unset<float> &);

ostream &
operator << (ostream & os, const Set_or_Unset<float> & qq)
{
  assert (os.good ());

  float tmp;
  if (qq.value (tmp))
    os << "value is " << tmp;
  else
    os << "value not set.";

  return os;
}

