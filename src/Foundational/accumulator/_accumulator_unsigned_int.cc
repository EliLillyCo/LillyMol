#include <iostream>
using namespace std;

#define ACCUMULATOR_IMPLEMENTATION
#include "accumulator.h"

template class Accumulator_Int<unsigned int>;
template class Accumulator_Base<unsigned int, unsigned int>;

template ostream & operator << (ostream &, const Accumulator_Int<int> &);
