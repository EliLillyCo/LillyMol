#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "substructure.h"

//template class resizable_array_p<Substructure_Chiral_Centre>;
//template class resizable_array_base<Substructure_Chiral_Centre *>;

template resizable_array_p<Substructure_Chiral_Centre>::resizable_array_p();
template resizable_array_p<Substructure_Chiral_Centre>::resizable_array_p(int);
template resizable_array_p<Substructure_Chiral_Centre>::~resizable_array_p();
template resizable_array_base<Substructure_Chiral_Centre*>::resizable_array_base();
template resizable_array_base<Substructure_Chiral_Centre*>::~resizable_array_base();

template int resizable_array_p<Substructure_Chiral_Centre>::resize(int);
template int resizable_array_base<Substructure_Chiral_Centre*>::resize(int);
template int resizable_array_p<Substructure_Chiral_Centre>::resize_no_delete(int);

template int resizable_array_base<Substructure_Chiral_Centre*>::add(Substructure_Chiral_Centre*);
