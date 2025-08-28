#include <stdlib.h>

#define IWQSORT_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"
#define IWDISTANCE_MATRIX_IMPLEMENTATION
#include "IWDistanceMatrixBase.h"

template class IWDistanceMatrixBase<unsigned char>;
template class ID_Distance<unsigned char>;
template void iwqsort (ID_Distance<unsigned char> *, int);
