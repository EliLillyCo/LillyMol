#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "path.h"

//template class resizable_array<Ring *>;
//template class resizable_array<const Ring *>;
//template class resizable_array_base<const Ring *>;

template resizable_array<Ring const*>::resizable_array();
template resizable_array<Ring const*>::~resizable_array();

template resizable_array<Ring*>::resizable_array();
template resizable_array<Ring*>::~resizable_array();

template int resizable_array<Ring*>::add_if_not_already_present(Ring*);
template int resizable_array_base<Ring const*>::resize(int);
template int resizable_array_base<Ring const*>::add(Ring const*);
template int resizable_array_base<Ring const*>::contains(Ring const*) const;
template resizable_array<Ring*> & resizable_array<Ring*>::operator=(resizable_array<Ring*> const&);
template int resizable_array<Ring const*>::resize_keep_storage(int);

#if ! defined(NDEBUG)
template Ring * & resizable_array_base<Ring*>::operator[](int) const;
#endif
