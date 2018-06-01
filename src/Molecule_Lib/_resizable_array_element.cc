#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "element.h"

//template class resizable_array<const Element *>;
//template class resizable_array_base<const Element *>;


template resizable_array<Element const*>::resizable_array();
template resizable_array<Element const*>::resizable_array(int);
template resizable_array<Element const*>::~resizable_array();

template int resizable_array_base<Element const*>::add(Element const*);
template void resizable_array_base<Element const*>::operator+=(resizable_array_base<Element const*> const&);
template resizable_array<Element const *> & resizable_array<Element const*>::operator=(resizable_array<Element const*> const&);
template int resizable_array_base<Element const*>::contains(Element const*) const;
template int resizable_array<Element const*>::resize_keep_storage(int);
template int resizable_array<Element const*>::add_if_not_already_present(Element const*);

template int resizable_array_base<Element const*>::ok() const;
template const Element * resizable_array_base<Element const*>::item(int) const;

#if ! defined(NDEBUG)
template const Element * & resizable_array_base<Element const*>::operator[](int) const;
#endif
