#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwaray_np_.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class resizable_array_np<float>;
template class resizable_array_iterator<float>;

#else

static void
unused ()
{
  resizable_array_np<float> a1;
//resizable_array_iterator<float> a2;
}

#endif
