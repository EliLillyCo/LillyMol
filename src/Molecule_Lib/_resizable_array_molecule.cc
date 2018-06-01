#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "molecule.h"

//template class resizable_array_p<Molecule>;
//template class resizable_array_base<Molecule *>;

template resizable_array_p<Molecule>::resizable_array_p();
template resizable_array_p<Molecule>::resizable_array_p(int);
template resizable_array_p<Molecule>::~resizable_array_p();
template resizable_array_base<Molecule*>::~resizable_array_base();
template resizable_array_base<Molecule*>::resizable_array_base();

template int resizable_array_base<Molecule*>::add(Molecule*);
template int resizable_array_p<Molecule>::remove_item(int);
template void resizable_array_base<Molecule*>::sort(int (*)(Molecule* const*, Molecule* const*));
template int resizable_array_base<Molecule*>::resize(int);
template int resizable_array_p<Molecule>::resize(int);
template int resizable_array_p<Molecule>::resize_no_delete(int);
template int resizable_array_p<Molecule>::resize_keep_storage(int);
template int resizable_array_base<Molecule*>::remove_item(int);

#if ! defined(NDEBUG)
template Molecule * & resizable_array_base<Molecule*>::operator[](int) const;
template int resizable_array_base<Molecule*>::ok() const;
#endif
