#include <stdlib.h>

#define SET_DATAITEM_VALUE_IMPLEMENTATION
#include "iw_tdt.h"

template int IW_TDT::set_dataitem_value (const char *, int, const IWString &, int);
template int IW_TDT::set_dataitem_value (const char *, int, const const_IWSubstring &, int);
