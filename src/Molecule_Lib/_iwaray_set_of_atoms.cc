#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "molecule.h"

template class resizable_array_p<Set_of_Atoms>;
template class resizable_array_base<Set_of_Atoms *>;
