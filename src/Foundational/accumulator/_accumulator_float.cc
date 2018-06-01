#include <iostream>
using namespace std;

#define ACCUMULATOR_IMPLEMENTATION
#include "accumulator.h"

template class Accumulator<float>;
template class Accumulator_Base<float, KahanSum>;
template class Accumulator_Base<float, float>;

#ifdef __GNUG__

template ostream & operator << (std::ostream &, const Accumulator<float> &);

#else

template ostream & operator << (ostream &, const Accumulator<float> &);

#endif
