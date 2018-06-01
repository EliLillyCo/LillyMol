#include <iostream>

#define SPACE_VECTOR_IMPLEMENTATION
#include "space_vector.h"

template class Space_Vector<float>;
template std::ostream & operator << (std::ostream &, const Space_Vector<float> &);
