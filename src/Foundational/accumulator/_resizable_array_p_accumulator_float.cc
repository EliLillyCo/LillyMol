#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/iwaray/iwaray.h"
#include "accumulator.h"

template class resizable_array_p<Accumulator<float> >;
template class resizable_array_base<Accumulator<float> * >;
