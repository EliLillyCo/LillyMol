#include <stdlib.h>


#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "mdl_file_data.h"

//template class resizable_array_p<MDL_Atom_Data>;
//template class resizable_array_base<MDL_Atom_Data*>;

template resizable_array_p<MDL_Atom_Data>::resizable_array_p();
template resizable_array_p<MDL_Atom_Data>::resizable_array_p(int);
template resizable_array_p<MDL_Atom_Data>::~resizable_array_p();

template int resizable_array_base<MDL_Atom_Data*>::swap_elements(int, int);
template int resizable_array_p<MDL_Atom_Data>::resize_keep_storage(int);
template int resizable_array_p<MDL_Atom_Data>::remove_item(int);
template int resizable_array_base<MDL_Atom_Data*>::add(MDL_Atom_Data*);

//template class resizable_array_p<MDL_Bond_Data>;
//template class resizable_array_base<MDL_Bond_Data*>;

template resizable_array_p<MDL_Bond_Data>::resizable_array_p();
template resizable_array_p<MDL_Bond_Data>::resizable_array_p(int);
template resizable_array_p<MDL_Bond_Data>::~resizable_array_p();
template resizable_array_base<MDL_Atom_Data*>::resizable_array_base();
template resizable_array_base<MDL_Bond_Data*>::~resizable_array_base();
template resizable_array_base<MDL_Atom_Data*>::~resizable_array_base();
template resizable_array_base<MDL_Bond_Data*>::resizable_array_base();
template int resizable_array_p<MDL_Atom_Data>::resize(int);
template int resizable_array_p<MDL_Bond_Data>::resize_no_delete(int);

template int resizable_array_base<MDL_Bond_Data*>::add(MDL_Bond_Data*);
template int resizable_array_p<MDL_Bond_Data>::remove_item(int);
template int resizable_array_p<MDL_Bond_Data>::resize_keep_storage(int);
template int resizable_array_p<MDL_Bond_Data>::resize(int);
template int resizable_array_base<MDL_Atom_Data*>::remove_item(int);
template int resizable_array_base<MDL_Atom_Data*>::resize(int);
template int resizable_array_p<MDL_Atom_Data>::resize_no_delete(int);
template int resizable_array_base<MDL_Bond_Data*>::remove_item(int);

template int resizable_array_base<MDL_Bond_Data*>::resize(int);
template resizable_array_p<MDL_Atom_Data> & resizable_array_p<MDL_Atom_Data>::operator=(resizable_array_p<MDL_Atom_Data> && rhs);
template resizable_array_p<MDL_Bond_Data> & resizable_array_p<MDL_Bond_Data>::operator=(resizable_array_p<MDL_Bond_Data> && rhs);

#if ! defined(NDEBUG)
template MDL_Atom_Data * & resizable_array_base<MDL_Atom_Data*>::operator[](int) const;
template MDL_Bond_Data * & resizable_array_base<MDL_Bond_Data*>::operator[](int) const;
#endif
