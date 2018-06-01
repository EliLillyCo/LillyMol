#include <stdlib.h>
#include <assert.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "iwgrep-2.5.h"
#include "iwcrex.h"

template class resizable_array_p<IW_Regular_Expression_Template<IW_grep_25_regex> >;
template class resizable_array_base<IW_Regular_Expression_Template<IW_grep_25_regex> * >;
