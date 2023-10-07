#include <iostream>

#define MINMAXSPC_IMPLEMENTATION
#define MINMAXSPC_OP_IMPLEMENTATION
#include "minmaxspc.h"

template class Min_Max_Specifier<double>;

template std::ostream & operator << (std::ostream &, const Min_Max_Specifier<double> &);
