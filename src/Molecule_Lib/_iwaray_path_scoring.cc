#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"

#include "path_scoring.h"

template class resizable_array_p<Path_Scoring>;
template class resizable_array_base<Path_Scoring *>;

template class resizable_array_p<Atomic_Numbers_Encounterd>;
template class resizable_array_base<Atomic_Numbers_Encounterd *>;
