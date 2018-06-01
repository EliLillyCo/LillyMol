#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "element.h"

//template class resizable_array_p<Element>;
//template class resizable_array_base<Element *>;

template resizable_array_p<Element>::resizable_array_p();
template resizable_array_p<Element>::resizable_array_p(int);
template resizable_array_p<Element>::~resizable_array_p();
template resizable_array_base<Element const*>::resizable_array_base();
template resizable_array_base<Element const*>::~resizable_array_base();
template resizable_array_base<Element*>::resizable_array_base();
template resizable_array_base<Element*>::~resizable_array_base();

template int resizable_array_p<Element>::resize(int);
template int resizable_array_base<Element*>::add(Element*);
template int resizable_array_base<Element const*>::resize(int);

template int resizable_array_p<Element>::resize_no_delete(int);
template int resizable_array_base<Element*>::resize(int);
template Element * resizable_array_base<Element*>::last_item() const;

#if ! defined(NDEBUG)
template Element * & resizable_array_base<Element*>::operator[](int) const;
template const Element * & resizable_array_base<Element const*>::operator[](int) const;
template int resizable_array_base<Element*>::ok() const;
template const Element * resizable_array_base<Element const*>::item(int) const;
#endif
