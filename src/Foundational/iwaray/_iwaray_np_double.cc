#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwaray_np_.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class resizable_array_np<double>;
template class resizable_array_iterator<double>;

#else

static void
unused ()
{
  resizable_array_np<double> f1;
//resizable_array_iterator<double> f1;
}

#endif
