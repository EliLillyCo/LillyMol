#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/iwaray/iwaray.h"
#include "rmele.h"

template class resizable_array_p<Element_to_Remove>;
template class resizable_array_base<Element_to_Remove *>;
