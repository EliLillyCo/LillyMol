#include <stdlib.h>

#define IWFACTORIAL_IMPLEMENTATION
#include "iwfactorial.h"

// Be careful, this typedef is in iwmtypes in the molecule directory

#ifndef _IW_BASIC_MOLECULE_TYPES
typedef unsigned long long iw_uint64_t;
#endif

template class IW_Factorial<iw_uint64_t>;
