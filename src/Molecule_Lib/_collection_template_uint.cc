#include <stdlib.h>
#include <iostream>
using namespace std;


#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "collection_template.h"

#ifdef __GNUG__

template class Collection_Template<unsigned int>;

#else

#endif
