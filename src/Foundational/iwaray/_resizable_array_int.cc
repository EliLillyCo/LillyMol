#include <iostream>

using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwaray_op.h"

template class resizable_array<int>;
template class resizable_array_base<int>;

template ostream & operator << (ostream &, const resizable_array<int> &);

//template int resizable_array_base<int>::sum () const;

template int operator== (int, const resizable_array<int> &);


