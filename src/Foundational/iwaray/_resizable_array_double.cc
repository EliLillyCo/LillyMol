#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwaray_op.h"

template class resizable_array<double>;
template class resizable_array_base<double>;

template std::ostream & operator << (std::ostream &, const resizable_array<double> &);

template int operator== (int, const resizable_array<double> &);
