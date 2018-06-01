#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "bond.h"

//template class resizable_array_p<Bond>;
//template class resizable_array_base<Bond *>;

template resizable_array_p<Bond>::resizable_array_p();
template resizable_array_p<Bond>::~resizable_array_p();

template resizable_array_base<Bond*>::resizable_array_base();
template resizable_array_base<Bond const*>::resizable_array_base();
template resizable_array_base<Bond*>::~resizable_array_base();

template int resizable_array_p<Bond>::resize_keep_storage(int);
template Bond * resizable_array_base<Bond*>::item(int) const;
template int resizable_array_p<Bond>::remove_item(int);
template Bond * resizable_array_p<Bond>::remove_no_delete(int);
template int resizable_array_base<Bond*>::remove_first(Bond*);
template int resizable_array_base<Bond*>::ok() const;
template void resizable_array_base<Bond*>::sort(int (*)(Bond* const*, Bond* const*));
template int resizable_array_base<Bond*>::add(Bond*);
template int resizable_array_p<Bond>::resize(int);
template int resizable_array_p<Bond>::resize_no_delete(int);
template int resizable_array_base<Bond*>::remove_item(int);
template int resizable_array_base<Bond*>::resize(int);
template Bond * & resizable_array_base<Bond*>::operator[](int) const;
template resizable_array_p<Bond> & resizable_array_p<Bond>::operator=(resizable_array_p<Bond> && rhs);

template resizable_array_base<Bond const*>::~resizable_array_base();
template int resizable_array_base<Bond const*>::resize(int);

#if ! defined(NDEBUG)
template const Bond * & resizable_array_base<Bond const*>::operator[](int) const;
template Bond * resizable_array_base<Bond*>::last_item() const;
#endif
