#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwbits.h"

template class resizable_array_p<IW_Bits_Base>;
template class resizable_array_base<IW_Bits_Base *>;
