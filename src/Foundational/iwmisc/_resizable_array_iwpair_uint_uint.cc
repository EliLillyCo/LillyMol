#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwpair.h"

typedef IW_Pair<unsigned int, unsigned int> puipui;

template class resizable_array_p<puipui>;
template class resizable_array_base<puipui *>;

