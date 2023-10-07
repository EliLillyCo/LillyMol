#include <stdlib.h>

#include "Foundational/iwstring/iwstring.h"

#define DATAITEM_VALUE_TEMPLATE_IMPLEMENTATION
#define SET_DATAITEM_VALUE_IMPLEMENTATION
#include "iw_tdt.h"

template int IW_TDT::_dataitem_value_string (const char *, int, IWString &,          int) const;
//template int IW_TDT::_dataitem_value_string (const char *, int, const_IWSubstring &, int) const;
template int IW_TDT::_dataitem_value_string<const_IWSubstring>(char const *, int, const_IWSubstring &, int) const;

template int IW_TDT::dataitem (const char * tag, int len_tag, IWString & s, int which_one_to_return) const;
template int IW_TDT::dataitem (const char * tag, int len_tag, const_IWSubstring & s, int which_one_to_return) const;
//template<> int IW_TDT::set_dataitem_value<IWString>(char const*, int, IWString const&, int);
template int IW_TDT::set_dataitem_value(char const*, int, IWString const&, int);
