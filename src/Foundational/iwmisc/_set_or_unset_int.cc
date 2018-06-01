#include <stdlib.h>
#include <iostream>
using namespace std;

#define SET_OR_UNSET_IMPLEMENTATION
#include "set_or_unset.h"

template class Set_or_Unset<int>;

template ostream & operator << (ostream &, const Set_or_Unset<int> &);
