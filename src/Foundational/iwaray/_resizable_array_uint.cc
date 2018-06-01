#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwaray_op.h"

template class resizable_array<unsigned int>;
template class resizable_array_base<unsigned int>;

template std::ostream & operator << (std::ostream &, const resizable_array<unsigned int> &);

