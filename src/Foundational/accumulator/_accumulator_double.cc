#include <iostream>
using namespace std;

#define ACCUMULATOR_IMPLEMENTATION
#include "accumulator.h"

template class Accumulator<double>;
template class Accumulator_Base<double, KahanSum>;

template ostream & operator << (ostream &, const Accumulator<double> &);
