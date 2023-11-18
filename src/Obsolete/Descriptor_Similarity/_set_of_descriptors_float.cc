#include <stdlib.h>
#include <assert.h>

#define IWARAY_IMPLEMENTATION
#include "iwaray.h"

#define SET_IF_DESCRIPTORS_IMPLEMENTATION
#include "set_of_descriptors.h"

#include "iwdescriptor.h"

template class Set_of_Descriptors<IWDescriptors<float, float> >;

template class iwaray<IWDescriptors<float, float> >;
