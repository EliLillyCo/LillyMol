#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "atom.h"

//template class resizable_array_p<Atom>;
//template class resizable_array<Atom *>;
//template class resizable_array_base<Atom *>;
//template class resizable_array_p<Connection>;
// template class resizable_array_base<Connection *>;

template resizable_array_p<Atom>::resizable_array_p();
template resizable_array_p<Atom>::~resizable_array_p();
template resizable_array<Atom*>::resizable_array();
template resizable_array<Atom*>::~resizable_array();
template resizable_array_base<Atom*>::~resizable_array_base();
template resizable_array_base<Atom*>::resizable_array_base();

template int resizable_array_base<Atom*>::ok() const;
template Atom * resizable_array_p<Atom>::remove_no_delete(int);
template int resizable_array_p<Atom>::resize(int);
template int resizable_array_p<Atom>::remove_items(int const*);
template int resizable_array_p<Atom>::remove_item(int);
template int resizable_array_base<Atom*>::index(Atom*) const;
template int resizable_array_base<Atom*>::add(Atom*);
template int resizable_array_base<Atom*>::resize(int);
template int resizable_array_p<Atom>::resize_no_delete(int);
template int resizable_array_base<Atom*>::remove_item(int);
template int resizable_array_base<Atom*>::remove_items(int const*);

template resizable_array_p<Connection>::resizable_array_p();
template resizable_array_p<Connection>::~resizable_array_p();
template resizable_array_base<Connection*>::~resizable_array_base();
template resizable_array_base<Connection*>::resizable_array_base();

template int resizable_array_base<Connection*>::add(Connection*);
template int resizable_array_p<Connection>::remove_item(int);
template int resizable_array_p<Connection>::resize(int);
template int resizable_array_base<Connection*>::resize(int);
template int resizable_array_p<Connection>::resize_no_delete(int);
template int resizable_array_base<Connection*>::remove_item(int);
template resizable_array_p<Atom> & resizable_array_p<Atom>::operator=(resizable_array_p<Atom> && rhs);


#if ! defined(NDEBUG)
template Connection * & resizable_array_base<Connection*>::operator[](int) const;
template int resizable_array_base<Connection*>::ok() const;
#endif
