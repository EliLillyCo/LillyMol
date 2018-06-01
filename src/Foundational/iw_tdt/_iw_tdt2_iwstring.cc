#include <stdlib.h>

#include "iwstring.h"

#define DATAITEM_VALUE_TEMPLATE_IMPLEMENTATION
#include "iw_tdt.h"

template int IW_TDT2::dataitem_value (const const_IWSubstring &, IWString &, int) const;
