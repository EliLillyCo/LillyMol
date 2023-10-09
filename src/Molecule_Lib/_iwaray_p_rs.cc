#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "substructure.h"

template class resizable_array_p<Substructure_Ring_Specification>;
template class resizable_array_base<Substructure_Ring_Specification *>;

template class resizable_array_p<Substructure_Ring_System_Specification>;
template class resizable_array_base<Substructure_Ring_System_Specification *>;
