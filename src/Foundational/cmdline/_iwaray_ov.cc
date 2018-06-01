#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "cmdline.h"

template class resizable_array_p<Option_and_Value>;
template class resizable_array_base<Option_and_Value *>;

template int Option_and_Value::value(int&) const;
template int Option_and_Value::value(unsigned int&) const;
template int Option_and_Value::value(long&) const;
template int Option_and_Value::value(unsigned long&) const;
template int Option_and_Value::value(long long&) const;
template int Option_and_Value::value(unsigned long long&) const;
template int Option_and_Value::value(float&) const;
template int Option_and_Value::value(double&) const;
