#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "ematch.h"

//template class resizable_array_p<Element_Matcher>;
//template class resizable_array_base<Element_Matcher *>;

template resizable_array_p<Element_Matcher>::resizable_array_p();
template resizable_array_p<Element_Matcher>::resizable_array_p(int);
template resizable_array_p<Element_Matcher>::~resizable_array_p();
template int resizable_array_base<Element_Matcher*>::add(Element_Matcher*);
template int resizable_array_p<Element_Matcher>::resize(int);
template resizable_array_base<Element_Matcher*>::~resizable_array_base();
template int resizable_array_p<Element_Matcher>::resize_no_delete(int);
template int resizable_array_base<Element_Matcher*>::resize(int);
template resizable_array_base<Element_Matcher*>::resizable_array_base();

#if ! defined(NDEBUG)
template int resizable_array_base<Element_Matcher*>::ok() const;
#endif
