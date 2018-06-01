#include <stdlib.h>
#include <assert.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwsed-4.0.5.h"
#include "iwcrex.h"

template class resizable_array_p<IW_Regular_Expression_Template<IW_sed_405_regex> >;
template class resizable_array_base<IW_Regular_Expression_Template<IW_sed_405_regex> * >;
