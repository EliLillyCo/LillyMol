#ifndef IWDMSUPPORT_H
#define IWDMSUPPORT_H

#include "iw_stl_hash_map.h"

extern void write_index_and_id (int i,
                    const IWString * _id,
                    std::ostream & output);
extern void
write_smiles_and_id (const IWString & id,
                     const IW_STL_Hash_Map_String & id_to_smiles,
                     std::ostream & output);

template <typename T> int parse_directive (const IWString &, T &);
#endif
