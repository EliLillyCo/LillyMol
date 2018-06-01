#include <stdlib.h>

#define IWCREX_IMPLEMENTATION
#include "iwcrex.h"
#include "iwlrx.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class IW_Regular_Expression_Template<IW_lrx>;

#else

static void
unused ()
{
  IW_Regular_Expression_Template<IW_lrx> foo;

  const_IWSubstring p;
  foo.set_pattern (p);
}

#endif
