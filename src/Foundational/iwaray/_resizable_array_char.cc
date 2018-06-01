#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwaray_op.h"

template class resizable_array<const char *>;
template class resizable_array_base<const char *>;

template std::ostream & operator << (std::ostream &, const resizable_array<char> &);

template class resizable_array<char>;
template class resizable_array_base<char>;

template int operator== (int, const resizable_array<char> &);
