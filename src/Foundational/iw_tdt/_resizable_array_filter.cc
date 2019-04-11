#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"

#include "iw_tdt_filter.h"

#if defined(__INTEL_COMPILER)
#else
template class resizable_array_p<IW_TDT_Filter_Base>;
template class resizable_array_base<IW_TDT_Filter_Base *>;
#endif
