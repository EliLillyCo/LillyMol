#include <iostream>

#define ACCUMULATOR_IMPLEMENTATION
#include "accumulator.h"

template class Accumulator_Int<unsigned int>;
template class Accumulator_Base<unsigned int, unsigned int>;

template std::ostream & operator << (std::ostream &, const Accumulator_Int<int> &);
