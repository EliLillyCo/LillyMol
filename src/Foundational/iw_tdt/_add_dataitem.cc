#include <stdlib.h>

#define DATAITEM_VALUE_TEMPLATE_IMPLEMENTATION
#define IWTDT_ADD_DATAITEM_IMPLEMENTATION
#include "iw_tdt.h"

template int IW_TDT::add_dataitem (const char * tag,
                      int len_tag,
                      const IWString & s,
                      int where_to_put);
template int IW_TDT::add_dataitem (const char * tag,
                      int len_tag,
                      const const_IWSubstring & s,
                      int where_to_put);
