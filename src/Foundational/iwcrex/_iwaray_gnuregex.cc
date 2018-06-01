#include <stdlib.h>
#include <assert.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwgnuregex.h"
#include "iwcrex.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class resizable_array_p<IW_Regular_Expression_Template<IW_gnu_regex> >;
template class resizable_array_base<IW_Regular_Expression_Template<IW_gnu_regex> * >;

#else

static void
unused ()
{
  resizable_array_p<IW_Regular_Expression_Template<IW_gnu_regex> > foo;
  foo.add (NULL);

  foo[0];
}

#endif
