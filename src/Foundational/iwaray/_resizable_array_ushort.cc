#include <iostream>


#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwaray_op.h"

template class resizable_array<unsigned short>;
template class resizable_array_base<unsigned short>;

//template ostream & operator << (ostream &, const resizable_array<int> &);

//template int operator== (int, const resizable_array<int> &);


