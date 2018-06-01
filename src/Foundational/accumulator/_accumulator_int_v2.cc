#include <iostream>
using namespace std;

#define ACCUMULATOR_V2_IMPLEMENTATION
#include "accumulator_v2.h"

template class Accumulator_V2<unsigned int, int>;

template ostream & operator << (ostream &, const Accumulator_V2<unsigned int, int> &);
