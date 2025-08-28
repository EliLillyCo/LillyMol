#ifndef UTILITIES_GFP_TOOLS_SPREAD_WEIGHT_H_
#define UTILITIES_GFP_TOOLS_SPREAD_WEIGHT_H_

#include <optional>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

namespace spread_weights {

// Used with various spread tools to hold information about any
// distance scaling that might be in effect.

class ScalingFactor {
  private:
    // The scaling factor can be specified several ways.
    // As a column in the name.
    int _scaling_factor_column;

    // A file containing a mapping from id to weight.
    int _column_in_weight_file;
    IW_STL_Hash_Map_float _id_to_weight;

    // Can also be a TDT tag;
    IWString _tag;

    // A quick check to as to whether we have been initialised or not.
    bool _active;

    // Private functions
    int ReadWeightRecord(const const_IWSubstring& buffer);
    int ReadWeights(iwstring_data_source& input, int weight_file_has_header);
    int ReadWeights(const IWString& fname, int weight_file_has_header);

    std::optional<float> WeightFromTokenInName(const IWString& name) const;

  public:
    ScalingFactor();

    int Initialise(Command_Line& cl, char flag, bool verbose);

    bool active() const {
      return _active;
    }

    // Once the means of scaling has been established, for any identifier,
    // return the associated weight.
    // Fails if the column is missing, or if the id is not in the hash.
    std::optional<float> Weight(const IWString& name) const;

};  // class ScalingFactor

}  // namespace spread_weights

#endif // UTILITIES_GFP_TOOLS_SPREAD_WEIGHT_H_
