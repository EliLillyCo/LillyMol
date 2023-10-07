#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "path.h"

//template class resizable_array_p<Ring>;
//template class resizable_array_base<Ring *>;

template resizable_array_p<Ring>::resizable_array_p();
template resizable_array_p<Ring>::~resizable_array_p();
template resizable_array_base<Ring const*>::resizable_array_base();
template resizable_array_base<Ring const*>::~resizable_array_base();
template resizable_array_base<Ring*>::resizable_array_base();
template resizable_array_base<Ring*>::~resizable_array_base();

template int resizable_array_base<Ring*>::contains(Ring*) const;
template Ring * resizable_array_p<Ring>::remove_no_delete(int);
template int resizable_array_p<Ring>::transfer_in(resizable_array_p<Ring>&);
template int resizable_array_base<Ring*>::remove_first(Ring*);
template void resizable_array_base<Ring*>::sort(int (*)(Ring* const*, Ring* const*));
template int resizable_array_p<Ring>::remove_item(int);
template int resizable_array_p<Ring>::transfer_in(resizable_array_p<Ring>&, int);
template int resizable_array_base<Ring*>::add(Ring*);
template int resizable_array_p<Ring>::resize(int);
template int resizable_array_base<Ring*>::remove_item(int);
template int resizable_array_base<Ring*>::resize(int);
template int resizable_array_p<Ring>::resize_no_delete(int);

#if ! defined(NDEBUG)
template Ring * & resizable_array_base<Ring*>::operator[](int) const;
template const Ring * & resizable_array_base<Ring const*>::operator[](int) const;
template int resizable_array_base<Ring*>::ok() const;
template int resizable_array_base<Ring const*>::ok() const;
#endif

