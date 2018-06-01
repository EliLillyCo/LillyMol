#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "iwaray.h"
#include "pearlman.h"

//template class resizable_array_p<Beep>;
//template class resizable_array_base<Beep *>;

template resizable_array_p<Beep>::resizable_array_p();
template resizable_array_p<Beep>::resizable_array_p(int);
template resizable_array_p<Beep>::~resizable_array_p();
template resizable_array_base<Beep*>::~resizable_array_base();

template int resizable_array_p<Beep>::remove_item(int);
template int resizable_array_base<Beep*>::add(Beep*);
template void resizable_array_base<Beep*>::sort(int (*)(Beep* const*, Beep* const*));
template int resizable_array_p<Beep>::resize_no_delete(int);
template int resizable_array_p<Beep>::resize(int);
template int resizable_array_base<Beep*>::resize(int);
template resizable_array_base<Beep*>::resizable_array_base();
template int resizable_array_base<Beep*>::remove_item(int);

#if ! defined(NDEBUG)
template Beep * & resizable_array_base<Beep*>::operator[](int) const;
template int resizable_array_base<Beep*>::ok() const;
template int resizable_array_base<Path_Message*>::ok() const;
#endif
