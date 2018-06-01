#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "rmele.h"

template class resizable_array_p<Element_to_Remove>;
template class resizable_array_base<Element_to_Remove *>;
