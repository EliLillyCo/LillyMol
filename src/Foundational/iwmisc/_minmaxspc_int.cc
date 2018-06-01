#include <iostream>
using namespace std;

#define MINMAXSPC_IMPLEMENTATION
#define MINMAXSPC_OP_IMPLEMENTATION
#include "minmaxspc.h"

template class Min_Max_Specifier<int>;
template ostream & operator << (ostream &, const Min_Max_Specifier<int> &);
