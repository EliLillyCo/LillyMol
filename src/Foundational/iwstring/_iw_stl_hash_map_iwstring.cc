/*
  Template instantiation stuff for hash maps
*/

#include <stdlib.h>

#include "iw_stl_hash_map.h"
#include "iwhash.h"
#include "iwstring.h"

template class IW_STL_Hash_Map<IWString, IWString>;

#ifdef __GNUG__

//template class std::_Hashtable_node<std::pair<IWString const, int> >** std::fill_n<std::_Hashtable_node<std::pair<IWString const, int> >**, unsigned, std::_Hashtable_node<std::pair<IWString const, int> >*>(std::_Hashtable_node<std::pair<IWString const, int> >**, unsigned, std::_Hashtable_node<std::pair<IWString const, int> >* const&);
//template class std::_Hashtable_node<std::pair<IWString const, IWString> >** std::fill_n<std::_Hashtable_node<std::pair<IWString const, IWString> >**, unsigned, std::_Hashtable_node<std::pair<IWString const, IWString> >*>(std::_Hashtable_node<std::pair<IWString const, IWString> >**, unsigned, std::_Hashtable_node<std::pair<IWString const, IWString> >* const&);
//template void std::fill<std::__normal_iterator<std::_Hashtable_node<std::pair<IWString const, IWString> >**, std::vector<std::_Hashtable_node<std::pair<IWString const, IWString> >*, std::allocator<IWString> > >, std::_Hashtable_node<std::pair<IWString const, IWString> >*>(std::__normal_iterator<std::_Hashtable_node<std::pair<IWString const, IWString> >**, std::vector<std::_Hashtable_node<std::pair<IWString const, IWString> >*, std::allocator<IWString> > >, std::__normal_iterator<std::_Hashtable_node<std::pair<IWString const, IWString> >**, std::vector<std::_Hashtable_node<std::pair<IWString const, IWString> >*, std::allocator<IWString> > >, std::_Hashtable_node<std::pair<IWString const, IWString> >* const&);
//template void std::fill<std::__normal_iterator<std::_Hashtable_node<std::pair<IWString const, int> >**, std::vector<std::_Hashtable_node<std::pair<IWString const, int> >*, std::allocator<int> > >, std::_Hashtable_node<std::pair<IWString const, int> >*>(std::__normal_iterator<std::_Hashtable_node<std::pair<IWString const, int> >**, std::vector<std::_Hashtable_node<std::pair<IWString const, int> >*, std::allocator<int> > >, std::__normal_iterator<std::_Hashtable_node<std::pair<IWString const, int> >**, std::vector<std::_Hashtable_node<std::pair<IWString const, int> >*, std::allocator<int> > >, std::_Hashtable_node<std::pair<IWString const, int> >* const&);
//template std::__normal_iterator<std::_Hashtable_node<std::pair<IWString const, IWString> >**, std::vector<std::_Hashtable_node<std::pair<IWString const, IWString> >*, std::allocator<IWString> > > std::fill_n<std::__normal_iterator<std::_Hashtable_node<std::pair<IWString const, IWString> >**, std::vector<std::_Hashtable_node<std::pair<IWString const, IWString> >*, std::allocator<IWString> > >, unsigned, std::_Hashtable_node<std::pair<IWString const, IWString> >*>(std::__normal_iterator<std::_Hashtable_node<std::pair<IWString const, IWString> >**, std::vector<std::_Hashtable_node<std::pair<IWString const, IWString> >*, std::allocator<IWString> > >, unsigned, std::_Hashtable_node<std::pair<IWString const, IWString> >* const&);
//template std::__normal_iterator<std::_Hashtable_node<std::pair<IWString const, int> >**, std::vector<std::_Hashtable_node<std::pair<IWString const, int> >*, std::allocator<int> > > std::fill_n<std::__normal_iterator<std::_Hashtable_node<std::pair<IWString const, int> >**, std::vector<std::_Hashtable_node<std::pair<IWString const, int> >*, std::allocator<int> > >, unsigned, std::_Hashtable_node<std::pair<IWString const, int> >*>(std::__normal_iterator<std::_Hashtable_node<std::pair<IWString const, int> >**, std::vector<std::_Hashtable_node<std::pair<IWString const, int> >*, std::allocator<int> > >, unsigned, std::_Hashtable_node<std::pair<IWString const, int> >* const&);

//template _Hashtable_node<IWString> ** fill_n<_Hashtable_node<IWString> **, unsigned int, _Hashtable_node<IWString> *>(_Hashtable_node<IWString> **, unsigned int, _Hashtable_node<IWString> * const &);
//template void fill<_Hashtable_node<IWString> **, _Hashtable_node<IWString> *>(_Hashtable_node<IWString> **, _Hashtable_node<IWString> **, _Hashtable_node<IWString> * const &);
//template class _Hashtable_const_iterator<IWString, IWString, IWStringHash, _Identity<IWString>, equal_to<IWString>, allocator<IWString> >;
//template class vector<_Hashtable_node<pair<IWString const, IWString> > *, allocator<IWString> >;
//template class hashtable<pair<IWString const, IWString>, IWString, IWStringHash, _Select1st<pair<IWString const, IWString> >, equal_to<IWString>, allocator<IWString> >;
//template _Hashtable_node<pair<IWString const, IWString> > ** fill_n<_Hashtable_node<pair<IWString const, IWString> > **, unsigned int, _Hashtable_node<pair<IWString const, IWString> > *>(_Hashtable_node<pair<IWString const, IWString> > **, unsigned int, _Hashtable_node<pair<IWString const, IWString> > * const &);
//template void fill<_Hashtable_node<pair<IWString const, IWString> > **, _Hashtable_node<pair<IWString const, IWString> > *>(_Hashtable_node<pair<IWString const, IWString> > **, _Hashtable_node<pair<IWString const, IWString> > **, _Hashtable_node<pair<IWString const, IWString> > * const &);

#elif defined(IWSGI)

template class std::hashtable<IWString,IWString,IWStringHash,std::_Identity<IWString>,std::equal_to<IWString>,std::allocator<IWString> >;
template class std::_Hashtable_const_iterator<std::pair<const IWString,int>,IWString,IWStringHash,std::_Select1st<std::pair<const IWString,int> >,std::equal_to<IWString>,std::allocator<int> >;
template class std::hashtable<std::pair<const IWString,int>,IWString,IWStringHash,std::_Select1st<std::pair<const IWString,int> >,std::equal_to<IWString>,std::allocator<int> >;
template class std::vector<std::_Hashtable_node<std::pair<const IWString,int> >*,std::allocator<int> >;

#endif
