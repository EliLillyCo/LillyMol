#include <iostream>

#define SET_OR_UNSET_IMPLEMENTATION
#include "set_or_unset.h"

template class Set_or_Unset<int>;

template std::ostream & operator << (std::ostream &, const Set_or_Unset<int> &);
