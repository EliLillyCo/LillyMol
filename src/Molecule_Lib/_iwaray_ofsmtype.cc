#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/iwaray/iwaray.h"
#include "ostream_and_type.h"

template class resizable_array_p<ofstream_and_type>;
template class resizable_array_base<ofstream_and_type *>;
