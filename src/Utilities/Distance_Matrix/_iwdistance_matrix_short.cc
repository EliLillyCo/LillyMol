#include <stdlib.h>

#define IWQSORT_IMPLEMENTATION
#include "iwqsort.h"

#define IWDISTANCE_MATRIX_IMPLEMENTATION
#include "IWDistanceMatrixBase.h"

template class IWDistanceMatrixBase<unsigned short>;
template class ID_Distance<unsigned short>;
template void iwqsort (ID_Distance<unsigned short> *, int);
