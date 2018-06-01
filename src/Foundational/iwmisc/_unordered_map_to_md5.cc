#define IW_IMPLEMENTATIONS_EXPOSED

#include "sparse_fp_creator.h"

template int unordered_map_to_md5<int, IWString_and_File_Descriptor>(std::unordered_map<unsigned int, int, std::hash<unsigned int>, std::equal_to<unsigned int>, std::allocator<std::pair<unsigned int const, int> > > const&, IWString_and_File_Descriptor&);
template int unordered_map_to_md5<int, IWString>(std::unordered_map<unsigned int, int, std::hash<unsigned int>, std::equal_to<unsigned int>, std::allocator<std::pair<unsigned int const, int> > > const&, IWString&);
