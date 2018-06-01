#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "qry_and_demerit.h"
#include "path.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)
template class resizable_array_p<Query_and_Demerit_Value>;
template class resizable_array_base<Query_and_Demerit_Value *>;
#else

static void
unused ()
{
  resizable_array_p<Query_and_Demerit_Value> v;
}

#endif

