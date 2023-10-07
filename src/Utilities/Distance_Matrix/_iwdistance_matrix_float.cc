#include <stdlib.h>

#define IWQSORT_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

#define IWDISTANCE_MATRIX_IMPLEMENTATION
#include "IWDistanceMatrixBase.h"

template class IWDistanceMatrixBase<float>;
template class ID_Distance<float>;
template void iwqsort (ID_Distance<float> *, int);
