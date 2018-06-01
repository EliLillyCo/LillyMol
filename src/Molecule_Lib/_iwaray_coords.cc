#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "coordinates.h"

template class resizable_array_p<Coordinates>;
template class resizable_array_base<Coordinates *>;
