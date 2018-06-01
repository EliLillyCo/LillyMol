#include <stdlib.h>

#define APPEND_INT_FORM_IMPLEMENTATION 1

static const char * digits = "0123456789";

#include "iwstring.h"

template void IWString::_append_int_form (int);
template void IWString::_append_int_form (unsigned int);
template void IWString::_append_int_form (long);
template void IWString::_append_int_form (long long);
template void IWString::_append_int_form (unsigned long long);
template void IWString::_append_int_form (unsigned long int);
