#include <iostream>

#define ACCUMULATOR_IMPLEMENTATION
#include "accumulator.h"

template class Accumulator<double>;
template class Accumulator_Base<double, KahanSum>;

template std::ostream & operator << (std::ostream &, const Accumulator<double> &);
