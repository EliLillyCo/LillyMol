#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/iwaray/iwaray.h"
#include "minmaxspc.h"

template class resizable_array_p<Min_Max_Specifier<int> >;
template class resizable_array_base<Min_Max_Specifier<int> * >;

