#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "substructure.h"
#include "path.h"

template class resizable_array_p<Substructure_Query>;
template class resizable_array_base<Substructure_Query *>;

