#ifndef UTILITIES_GFP_TOOLS_NNDATA_H
#define UTILITIES_GFP_TOOLS_NNDATA_H

#include "Foundational/iwstring/iwstring.h"

#include "Utilities/GFP_Tools/nearneighbours.pb.h"

namespace gfp {

// Write `proto` to `output` as text_format.
int WriteNNData(const nnbr::NearNeighbours& proto,
                IWString_and_File_Descriptor& output);

}  // namespace gfp
#endif  // UTILITIES_GFP_TOOLS_NNDATA_H
