#include <iostream>

#define IWARCHIVE_IMPLEMENTATION
#define IWARCHIVE_OP_IMPLEMENTATION

#include "iwarchive.h"

template class iwarchive<int>;

template std::ostream & operator << (std::ostream &, const iwarchive<int> &);
