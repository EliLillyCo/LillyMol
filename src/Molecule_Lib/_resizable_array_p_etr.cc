#include <iostream>
using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "etrans.h"

//template class resizable_array_p<Element_Transformation>;
//template class resizable_array_base<Element_Transformation *>;

template resizable_array_p<Element_Transformation>::resizable_array_p();
template resizable_array_p<Element_Transformation>::resizable_array_p(int);
template resizable_array_p<Element_Transformation>::~resizable_array_p();
template resizable_array_base<Element_Transformation*>::~resizable_array_base();

template resizable_array_base<Element_Transformation*>::resizable_array_base();

template int resizable_array_p<Element_Transformation>::resize(int);
template int resizable_array_base<Element_Transformation*>::add(Element_Transformation*);
template int resizable_array_p<Element_Transformation>::resize_no_delete(int);
template int resizable_array_base<Element_Transformation*>::resize(int);

//template int resizable_array_p<Element_Transformation>::resize(int);

template int resizable_array_base<Element_Transformation*>::ok() const;
