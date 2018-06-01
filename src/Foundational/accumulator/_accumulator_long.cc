#include <iostream>
using namespace std;

#define ACCUMULATOR_IMPLEMENTATION
#include "accumulator.h"

template class Accumulator_Int<long>;
template class Accumulator_Base<long, long>;

template ostream & operator << (ostream &, const Accumulator_Int<long> &);
