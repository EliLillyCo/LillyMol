#include <stdlib.h>

#define IWQSORT_IMPLEMENTATION
#include "iwqsort.h"

#define IWDISTANCE_MATRIX_IMPLEMENTATION
#include "IWDistanceMatrixBase.h"

template class IWDistanceMatrixBase<double>;
template class ID_Distance<double>;
template void iwqsort (ID_Distance<double> *, int);
