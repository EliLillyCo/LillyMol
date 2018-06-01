#include <stdlib.h>

#define SPLIT_DV_IMPLEMENTATION
#include "iwstring.h"

template int IWString::split_into_directive_and_value (const_IWSubstring &, char, int &) const;

template int const_IWSubstring::split_into_directive_and_value (const_IWSubstring &, char, int &) const;
