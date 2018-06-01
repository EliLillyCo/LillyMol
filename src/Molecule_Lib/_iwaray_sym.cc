#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "molecule.h"
#include "symmetry.h"

template class resizable_array_p<Symmetric_Atoms>;
template class resizable_array_base<Symmetric_Atoms *>;
