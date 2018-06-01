#include <stdlib.h>

#define IWTDT_ADD_DATAITEM_IMPLEMENTATION
#include "iw_tdt.h"

template int IW_TDT::add_dataitem (const char * tag,
                      int len_tag,
                      const float & s,
                      int where_to_put);

