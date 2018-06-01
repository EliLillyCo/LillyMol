#include <stdlib.h>
#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "msi_object.h"

template class resizable_array_p<msi_attribute>;
template class resizable_array_base<msi_attribute *>;

template class resizable_array_p<msi_object>;
template class resizable_array_base<msi_object *>;
