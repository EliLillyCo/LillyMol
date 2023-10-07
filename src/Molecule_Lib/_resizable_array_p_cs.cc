#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "chiral_centre.h"

//template class resizable_array_p<Chiral_Centre>;
//template class resizable_array_base<Chiral_Centre *>;


template resizable_array_p<Chiral_Centre>::resizable_array_p();
template resizable_array_p<Chiral_Centre>::resizable_array_p(int);
template resizable_array_p<Chiral_Centre>::~resizable_array_p();
template resizable_array_base<Chiral_Centre*>::resizable_array_base();
template resizable_array_base<Chiral_Centre*>::~resizable_array_base();

template int resizable_array_p<Chiral_Centre>::remove_item(int);
template int resizable_array_base<Chiral_Centre*>::add(Chiral_Centre*);
template int resizable_array_p<Chiral_Centre>::resize_keep_storage(int);
template Chiral_Centre * resizable_array_p<Chiral_Centre>::remove_no_delete(int);
template int resizable_array_p<Chiral_Centre>::resize(int);
template int resizable_array_base<Chiral_Centre*>::resize(int);
template int resizable_array_base<Chiral_Centre*>::remove_item(int);
template resizable_array_p<Chiral_Centre> & resizable_array_p<Chiral_Centre>::operator=(resizable_array_p<Chiral_Centre> && rhs);

template int resizable_array_p<Chiral_Centre>::resize_no_delete(int);

template resizable_array<Chiral_Centre*>::resizable_array();
template resizable_array<Chiral_Centre*>::~resizable_array();
template int resizable_array<Chiral_Centre*>::resize_keep_storage(int);

#if ! defined(NDEBUG)
template Chiral_Centre * & resizable_array_base<Chiral_Centre*>::operator[](int) const;
template Chiral_Centre * resizable_array_base<Chiral_Centre*>::last_item() const;
template int resizable_array_base<Chiral_Centre*>::ok() const;
#endif
