#include <stdlib.h>

#define IWQSORT_IMPLEMENTATION
#include "iwqsort.h"
#define IWDISTANCE_MATRIX_IMPLEMENTATION
#include "IWDistanceMatrixBase.h"

template class IWDistanceMatrixBase<int>;
template class ID_Distance<int>;
template void iwqsort (ID_Distance<int> *, int);
