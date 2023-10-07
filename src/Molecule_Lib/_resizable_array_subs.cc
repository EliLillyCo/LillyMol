#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define EXTENDING_RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/iwaray/iwaray.h"

#include "substructure.h"

//template class resizable_array<Substructure_Bond *>;

template resizable_array<Substructure_Atom*>::resizable_array();
template resizable_array<Substructure_Atom*>::~resizable_array();

template int resizable_array<Substructure_Atom*>::add_if_not_already_present(Substructure_Atom*);
template int resizable_array<Substructure_Atom*>::resize_keep_storage(int);

//template class resizable_array<Substructure_Atom *>;
//template class resizable_array_p<Substructure_Atom>;

template resizable_array_p<Substructure_Atom>::resizable_array_p();
template resizable_array_p<Substructure_Atom>::~resizable_array_p();
template resizable_array_base<Substructure_Atom*>::resizable_array_base();

template int resizable_array_p<Substructure_Atom>::transfer_in(resizable_array_p<Substructure_Atom>&, int);

//template class resizable_array_base<Substructure_Atom *>;

template int resizable_array_base<Substructure_Atom*>::remove_first(Substructure_Atom*);
template Substructure_Atom * resizable_array_base<Substructure_Atom*>::pop();
template int resizable_array_base<Substructure_Atom*>::make_room_for_extra_items(int);
template int resizable_array_base<Substructure_Atom*>::ok() const;
template int resizable_array<Substructure_Atom*>::extend(int, Substructure_Atom*);
template resizable_array_base<Substructure_Atom*>::~resizable_array_base();
template int resizable_array_p<Substructure_Atom>::resize(int);
template int resizable_array_p<Substructure_Atom>::resize_no_delete(int);
//template int resizable_array_p<Substructure_Atom>::remove_no_delete(int);
//template int resizable_array_p<Substructure_Atom>::transfer_in(resizable_array_p<Substructure_Atom>&, int);
template Substructure_Atom * resizable_array_p<Substructure_Atom>::remove_no_delete(int);

//template class extending_resizable_array<Substructure_Atom *>;

template Substructure_Atom * & extending_resizable_array<Substructure_Atom*>::operator[](int);
template extending_resizable_array<Substructure_Atom*>::extending_resizable_array(Substructure_Atom*);

//template class resizable_array_p<Query_Atoms_Matched>;
//template class resizable_array_base<Query_Atoms_Matched *>;

template resizable_array_p<Query_Atoms_Matched>::resizable_array_p();
template resizable_array_p<Query_Atoms_Matched>::~resizable_array_p();
template resizable_array_base<Query_Atoms_Matched*>::resizable_array_base();

template int resizable_array_p<Query_Atoms_Matched>::remove_item(int);
template int resizable_array_p<Query_Atoms_Matched>::resize(int);
template int resizable_array_base<Query_Atoms_Matched*>::add(Query_Atoms_Matched*);
template int resizable_array_p<Query_Atoms_Matched>::resize_keep_storage(int);
template int resizable_array_base<Query_Atoms_Matched*>::resize(int);
template int resizable_array_p<Query_Atoms_Matched>::resize_no_delete(int);
template resizable_array_base<Query_Atoms_Matched*>::~resizable_array_base();
template int resizable_array_base<Query_Atoms_Matched*>::remove_item(int);
template int resizable_array_p<Substructure_Environment>::resize_no_delete(int);

//template class resizable_array_p<Substructure_Bond>;
//template class resizable_array_base<Substructure_Bond *>;

template resizable_array_p<Substructure_Bond>::resizable_array_p();
template resizable_array_p<Substructure_Bond>::~resizable_array_p();
template resizable_array_base<Substructure_Bond*>::resizable_array_base();
template resizable_array_base<Substructure_Bond*>::~resizable_array_base();

template int resizable_array_base<Substructure_Bond*>::add(Substructure_Bond*);
template int resizable_array_base<Substructure_Bond*>::resize(int);
template int resizable_array_p<Substructure_Bond>::resize_no_delete(int);
template int resizable_array_p<Substructure_Bond>::resize(int);

//template class resizable_array<Substructure_Environment *>;
//template class resizable_array_p<Substructure_Environment>;
//template class resizable_array_base<Substructure_Environment *>;

template resizable_array_p<Substructure_Environment>::resizable_array_p();
template resizable_array_p<Substructure_Environment>::~resizable_array_p();

template int resizable_array_base<Substructure_Environment*>::add(Substructure_Environment*);

template resizable_array<Substructure_Environment*>::~resizable_array();
template resizable_array<Substructure_Environment*>::resizable_array();
template int resizable_array_base<Substructure_Environment*>::resize(int);
template resizable_array_base<Substructure_Environment*>::~resizable_array_base();
template int resizable_array_p<Substructure_Environment>::resize(int);
template resizable_array_base<Substructure_Environment*>::resizable_array_base();

//template class resizable_array_p<Substructure_Atom_Specifier>;
//template class resizable_array_base<Substructure_Atom_Specifier *>;

template resizable_array_p<Substructure_Atom_Specifier>::resizable_array_p();
template resizable_array_p<Substructure_Atom_Specifier>::~resizable_array_p();
template resizable_array_base<Substructure_Atom_Specifier*>::resizable_array_base();
template int resizable_array_p<Substructure_Atom_Specifier>::resize(int);
template int resizable_array_p<Substructure_Atom_Specifier>::resize_no_delete(int);


template int resizable_array_base<Substructure_Atom_Specifier*>::add(Substructure_Atom_Specifier*);
template resizable_array_base<Substructure_Atom_Specifier*>::~resizable_array_base();

//template class resizable_array_p<Elements_Needed>;
//template class resizable_array_base<Elements_Needed *>;

template resizable_array_p<Elements_Needed>::resizable_array_p();
template resizable_array_p<Elements_Needed>::~resizable_array_p();
template resizable_array_base<Elements_Needed*>::resizable_array_base();
template resizable_array_base<Elements_Needed*>::~resizable_array_base();

template int resizable_array_base<Elements_Needed*>::add(Elements_Needed*);
template int resizable_array_base<Elements_Needed*>::resize(int);
template int resizable_array_p<Elements_Needed>::resize(int);
template int resizable_array_p<Elements_Needed>::resize_no_delete(int);


//template class resizable_array_p<Link_Atom>;
//template class resizable_array<Link_Atom *>;
//template class resizable_array_base<Link_Atom *>;

//template class resizable_array_p<ISIS_Link_Atom>;
//template class resizable_array_base<ISIS_Link_Atom *>;

template resizable_array_base<Link_Atom*>::resizable_array_base();
template resizable_array_base<Link_Atom*>::~resizable_array_base();
template int resizable_array_p<Link_Atom>::resize_no_delete(int);
template int resizable_array_p<Link_Atom>::resize(int);
template int resizable_array_base<Link_Atom*>::resize(int);

template resizable_array_p<Link_Atom>::resizable_array_p();
template resizable_array_p<Link_Atom>::~resizable_array_p();

template int resizable_array_base<Link_Atom*>::add(Link_Atom*);
template int resizable_array_p<Link_Atom>::transfer_in(resizable_array_p<Link_Atom>&);
template resizable_array<Link_Atom*>::resizable_array();
template resizable_array<Link_Atom*>::~resizable_array();


template resizable_array_p<ISIS_Link_Atom>::resizable_array_p();
template resizable_array_p<ISIS_Link_Atom>::resizable_array_p(int);
template resizable_array_p<ISIS_Link_Atom>::~resizable_array_p();
template resizable_array_base<ISIS_Link_Atom*>::~resizable_array_base();
template resizable_array_base<ISIS_Link_Atom*>::resizable_array_base();

template int resizable_array_p<ISIS_Link_Atom>::resize(int);
template int resizable_array_p<ISIS_Link_Atom>::resize_no_delete(int);
template int resizable_array_base<ISIS_Link_Atom*>::resize(int);
template resizable_array_p<Link_Atom> & resizable_array_p<Link_Atom>::operator=(resizable_array_p<Link_Atom> && rhs);


template int resizable_array_base<Substructure_Atom*>::add(Substructure_Atom*);
template int resizable_array_base<Substructure_Atom*>::resize(int);
template int resizable_array_base<Substructure_Atom*>::remove_two_items(Substructure_Atom*, Substructure_Atom*);
template int resizable_array_base<Substructure_Atom*>::remove_item(int);
//template int resizable_array_p<Substructure_Atom>::remove_no_delete(int);

template int resizable_array_base<Substructure_Atom_Specifier*>::resize(int);


#if ! defined(NDEBUG)
template Substructure_Atom_Specifier * & resizable_array_base<Substructure_Atom_Specifier*>::operator[](int) const;
template Elements_Needed * & resizable_array_base<Elements_Needed*>::operator[](int) const;
template Substructure_Bond * & resizable_array_base<Substructure_Bond*>::operator[](int) const;
template Substructure_Atom * & resizable_array_base<Substructure_Atom*>::operator[](int) const;
template Substructure_Environment * & resizable_array_base<Substructure_Environment*>::operator[](int) const;
template Query_Atoms_Matched * & resizable_array_base<Query_Atoms_Matched*>::operator[](int) const;
template Substructure_Atom * resizable_array_base<Substructure_Atom*>::item(int) const;
template Substructure_Chiral_Centre * & resizable_array_base<Substructure_Chiral_Centre*>::operator[](int) const;
template int resizable_array_base<Substructure_Chiral_Centre*>::ok() const;
template Link_Atom * & resizable_array_base<Link_Atom*>::operator[](int) const;
template int resizable_array_base<Query_Atoms_Matched*>::ok() const;
template int resizable_array_base<Substructure_Bond *>::ok() const;
template int resizable_array_base<Substructure_Environment*>::ok() const;
template int resizable_array_base<Substructure_Atom_Specifier>::ok() const;
template int resizable_array_base<Elements_Needed*>::ok() const;
template int resizable_array_base<Link_Atom*>::ok() const;
template int resizable_array_base<ISIS_Link_Atom*>::ok() const;
template int resizable_array_base<Substructure_Atom_Specifier*>::ok() const;


#endif

