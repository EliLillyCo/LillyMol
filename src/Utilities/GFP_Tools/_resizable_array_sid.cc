#include <stdlib.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "smiles_id_dist.h"

template class resizable_array_p<Smiles_ID_Dist>;
template class resizable_array_base<Smiles_ID_Dist *>;
