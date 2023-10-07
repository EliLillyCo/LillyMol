#include <stdlib.h>
#include <string.h>
#include <iostream>

#define IWQS_USE_OPERATOR
#define IWQSORT_IMPLEMENTATION

#include "iwqsort.h"

template void iwqsort (int *, int);
template void move_in_from_right<int>(int *, int, int &);
template void swap_elements<int>(int *, int *);
template void compare_two_items<int>(int *);
template void move_in_from_left<int>(int*, int&, int&, int);

