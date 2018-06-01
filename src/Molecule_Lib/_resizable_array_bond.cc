#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "bond.h"

//template class resizable_array<Bond *>;
//template class resizable_array<const Bond *>;
//template class resizable_array_base<const Bond *>;

template resizable_array<Bond*>::resizable_array();
template resizable_array<Bond*>::~resizable_array();

template int resizable_array_base<Bond const*>::add(Bond const*);
template int resizable_array_base<Bond const*>::insert_before(int, Bond const*);

template resizable_array<Bond const*>::resizable_array();
template resizable_array<Bond const*>::~resizable_array();

template int resizable_array<Bond const*>::add_if_not_already_present(Bond const*);
template int resizable_array_base<Bond*>::insert_before(int, Bond*);
template int resizable_array_base<Bond const*>::swap_elements(int, int);
template int resizable_array<Bond const*>::resize_keep_storage(int);

#if ! defined(NDEBUG)
template Bond * resizable_array_base<Bond*>::last_item() const;
template int resizable_array_base<Bond const*>::ok() const;
#endif

template const Bond * & resizable_array_base<Bond const*>::operator[](int) const;
