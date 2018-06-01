#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwaray_op.h"

template class resizable_array<long>;
template class resizable_array_base<long>;

template std::ostream & operator << (std::ostream &, const resizable_array<long> &);

