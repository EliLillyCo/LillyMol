#include <stdlib.h>

#define COMPILING_STREAMTYPE_TEMPLATE

#include "molecule.h"

#define ISTREAM_AND_TYPE_IMPLEMENTATION

#include "istream_and_type.h"

template class data_source_and_type<Molecule>;
