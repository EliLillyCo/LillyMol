#include "iwaray.h"
#include "iwaray_np_.h"

#if defined(__GNUG__) || defined (__SUNPRO_CC)

template class resizable_array_np<char>;
template class resizable_array_iterator<char>;

#else

static void
unused ()
{
   resizable_array_np<char> a1;
   resizable_array_iterator<char> a2 (a1);
}

#endif
