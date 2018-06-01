#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "dy_pool.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class resizable_array_p<Dy_Fingerprint_Pool_Component>;
template class resizable_array_base<Dy_Fingerprint_Pool_Component *>;

#else

static void
unused ()
{
  resizable_array_p<Dy_Fingerprint_Pool_Component>;
}

#endif
