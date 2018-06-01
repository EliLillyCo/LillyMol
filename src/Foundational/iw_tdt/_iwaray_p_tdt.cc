#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iw_tdt.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class resizable_array_p<IW_TDT>;
template class resizable_array_base<IW_TDT *>;

#else

#endif
