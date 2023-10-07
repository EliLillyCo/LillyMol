#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/iwaray/iwaray.h"
#include "coordinates.h"

template class resizable_array_p<Coordinates>;
template class resizable_array_base<Coordinates *>;
