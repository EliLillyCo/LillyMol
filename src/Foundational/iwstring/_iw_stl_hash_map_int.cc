/*
  Template instantiation stuff for hash maps
*/

#include <stdlib.h>

#include "iw_stl_hash_map.h"

#ifdef __GNUG__Q
template class vector<_Hashtable_node<pair<IWString const, int> > *, allocator<int> >;
template class  hashtable<pair<IWString const, int>, IWString, IWStringHash, _Select1st<pair<IWString const, int> >, equal_to<IWString>, allocator<int> >;
template _Hashtable_node<pair<IWString const, int> > ** fill_n<_Hashtable_node<pair<IWString const, int> > **, unsigned int, _Hashtable_node<pair<IWString const, int> > *>(_Hashtable_node<pair<IWString const, int> > **, unsigned int, _Hashtable_node<pair<IWString const, int> > * const &);
template class  _Hashtable_const_iterator<pair<IWString const, int>, IWString, IWStringHash, _Select1st<pair<IWString const, int> >, equal_to<IWString>, allocator<int> >;
template void fill<_Hashtable_node<pair<IWString const, int> > **, _Hashtable_node<pair<IWString const, int> > *>(_Hashtable_node<pair<IWString const, int> > **, _Hashtable_node<pair<IWString const, int> > **, _Hashtable_node<pair<IWString const, int> > * const &);
#endif
