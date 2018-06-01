#include "misc.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template int locate_item_in_array (int, int, const int *, int);

#else

static void
unused ()
{
  int ii[2];
  (void) locate_item_in_array (5, 0, ii);
}

#endif
