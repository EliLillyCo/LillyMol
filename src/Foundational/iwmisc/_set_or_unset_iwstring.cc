#include <stdlib.h>
#include <iostream>

#define SET_OR_UNSET_IMPLEMENTATION
#include "set_or_unset.h"

#include "Foundational/iwstring/iwstring.h"

template class Set_or_Unset<IWString>;

template std::ostream & operator << (std::ostream &, const Set_or_Unset<IWString> &);
