#include <iostream>

#define ACCUMULATOR_V2_IMPLEMENTATION
#include "accumulator_v2.h"

template class Accumulator_V2<unsigned int, int>;

template std::ostream & operator << (std::ostream &, const Accumulator_V2<unsigned int, int> &);
