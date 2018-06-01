#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwaray_op.h"

typedef long long longlong;

template class resizable_array<longlong>;
template class resizable_array_base<longlong>;

template std::ostream & operator << (std::ostream &, const resizable_array<longlong> &);

