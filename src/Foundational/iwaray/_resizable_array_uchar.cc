#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwaray_op.h"

template class resizable_array<unsigned char>;
template class resizable_array_base<unsigned char>;

//template class ostream & operator << (ostream &, const resizable_array<unsigned char> &);

template int operator== (int, const resizable_array<unsigned char> &);
