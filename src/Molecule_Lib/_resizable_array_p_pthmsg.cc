#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/iwaray/iwaray.h"

#include "pearlman.h"

//template class resizable_array_p<Path_Message>;
//template class resizable_array_base<Path_Message *>;

template resizable_array_p<Path_Message>::resizable_array_p();
template resizable_array_p<Path_Message>::~resizable_array_p();
template resizable_array_p<Path_Message>::resizable_array_p(int);
template resizable_array_base<Path_Message*>::~resizable_array_base();
template resizable_array_base<Path_Message*>::resizable_array_base();

template int resizable_array_p<Path_Message>::remove_item(int);
template int resizable_array_p<Path_Message>::resize_no_delete(int);
template int resizable_array_p<Path_Message>::resize_keep_storage(int);
template int resizable_array_base<Path_Message*>::add(Path_Message*);
template int resizable_array_p<Path_Message>::resize(int);
template int resizable_array_base<Path_Message*>::resize(int);
template int resizable_array_base<Path_Message*>::remove_item(int);


#if ! defined(NDEBUG)
template Path_Message * & resizable_array_base<Path_Message*>::operator[](int) const;
template const Path_Message * & resizable_array_base<Path_Message const*>::operator[](int) const;
#endif
