#include <stdlib.h>

#define IWCREX_IMPLEMENTATION
#include "iwcrex.h"
#include "iwlgen.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class IW_Regular_Expression_Template<IW_lgen>;

#else

static void
unused ()
{
  IW_Regular_Expression_Template<IW_lgen> foo;

}

#endif
