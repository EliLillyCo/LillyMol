#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "molecule.h"

//template class resizable_array_p<List_of_Ring_Sizes>;
//template class resizable_array_base<List_of_Ring_Sizes *>;
#if ! defined(NDEBUG)
template List_of_Ring_Sizes * & resizable_array_base<List_of_Ring_Sizes*>::operator[](int) const;
#endif
