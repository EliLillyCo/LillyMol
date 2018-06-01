#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "molecule.h"

#include "iwstandard.h"

template resizable_array_p<Possible_Lactim_Lactam>::resizable_array_p();
template resizable_array_p<Possible_Lactim_Lactam>::~resizable_array_p();
template resizable_array_p<Possible_Lactim_Lactam>::resizable_array_p(int);
template int resizable_array_base<Possible_Lactim_Lactam*>::add(Possible_Lactim_Lactam*);
template Possible_Lactim_Lactam * & resizable_array_base<Possible_Lactim_Lactam*>::operator[](int) const;
template resizable_array<Possible_Lactim_Lactam*>::resizable_array();
template resizable_array<Possible_Lactim_Lactam*>::~resizable_array();
template int resizable_array_base<Possible_Lactim_Lactam*>::remove_item(int);
template int resizable_array_p<Possible_Lactim_Lactam>::remove_item(int);

#if ! defined(NDEBUG)
#endif
