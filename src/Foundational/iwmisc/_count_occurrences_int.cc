#include "misc.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template int count_occurrences_of_item_in_array (int, int, const int *);

#else

static void
unused ()
{
  int ii[2];
  (void) count_occurrences_of_item_in_array (0, 0, ii);
}

#endif
