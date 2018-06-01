#include <stdlib.h>

#define IWCREX_IMPLEMENTATION
#include "iwcrex.h"
#include "iwlrx_maddock.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class IW_Regular_Expression_Template<IW_lrx_Maddock>;

#else

static void
unused ()
{
  IW_Regular_Expression_Template<IW_lrx_Maddock> foo;
}

#endif
