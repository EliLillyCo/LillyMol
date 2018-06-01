#include <iostream>
using namespace std;

#define ACCUMULATOR_IMPLEMENTATION
#include "accumulator.h"

template class Accumulator_Int<int>;
template class Accumulator_Base<int, int>;
template class Accumulator_Base<long, unsigned long>;
//template Accumulator_Base<int, KahanSum>::Accumulator_Base();
//template unsigned int Accumulator_Base<int, KahanSum>::extra(int);
//template double Accumulator_Base<int, int>::average() const;
//template void Accumulator_Base<int, int>::_default_values();

template std::ostream & operator << (std::ostream &, const Accumulator_Int<int> &);
