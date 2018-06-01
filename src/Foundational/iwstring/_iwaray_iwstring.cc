#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwstring.h"

template class resizable_array<const IWString *>;
template class resizable_array_base<const IWString *>;

template class resizable_array<IWString *>;
template class resizable_array_base<IWString *>;
