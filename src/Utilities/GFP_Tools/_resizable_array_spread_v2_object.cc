#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "spread_v2.h"

template class resizable_array<Spread_Object *>;
template class resizable_array_base<Spread_Object *>;
