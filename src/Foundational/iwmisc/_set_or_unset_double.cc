#include <iostream>

#define SET_OR_UNSET_IMPLEMENTATION
#include "set_or_unset.h"

template class Set_or_Unset<double>;

template std::ostream & operator << (std::ostream &, const Set_or_Unset<double> &);
