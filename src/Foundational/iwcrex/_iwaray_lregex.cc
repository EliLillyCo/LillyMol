#include <stdlib.h>
#include <assert.h>

#define IWARAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwlregex.h"
#include "iwcrex.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class iwaray<IW_Regular_Expression_Template<IW_lregex> >;

#else

static void
unused ()
{
}

#endif
