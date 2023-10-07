#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/iwaray/iwaray.h"
#include "path.h"
#include "substructure.h"

template class resizable_array_p<Single_Substructure_Query>;
template class resizable_array_base<Single_Substructure_Query *>;
