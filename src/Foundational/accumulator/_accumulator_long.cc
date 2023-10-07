#include <iostream>

#define ACCUMULATOR_IMPLEMENTATION
#include "accumulator.h"

template class Accumulator_Int<long>;
template class Accumulator_Base<long, long>;

template std::ostream & operator << (std::ostream &, const Accumulator_Int<long> &);
