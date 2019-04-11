#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "leader.h"

template class resizable_array<GFP_L *>;
template class resizable_array_base<GFP_L *>;
