#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

#else

static void
unused ()
{
  resizable_array_base<int *> foo;

  foo.item (0);
  foo.add (static_cast<int *> (nullptr));
}

#endif
