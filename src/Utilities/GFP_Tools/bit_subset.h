#ifndef UTILITIES_GFP_TOOLS_BIG_SUBSET_H
#define UTILITIES_GFP_TOOLS_BIG_SUBSET_H

#include <cstdint>
#include <unordered_set>

#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Utilities/GFP_Tools/gfp_to_svm_lite.pb.h"
#include "Utilities/GFP_Tools/dyfp.h"
#include "Utilities/GFP_Tools/gfp.h"

namespace bit_subset {

// THe type of the bits in a .gfp
using gfp_bit_type = uint32_t;

// The BitSubset class is used to create a specific subset of
// the bits in an IW_General_Fingerprint.
// This is complicated by the fact that the class can be
// instantiated from a proto, before any .gfp files have
// been processed, so what fingerprints are in the input
// is not necessarily known when this is constructed.
class BitSubset {
  private:
    // Makking from tag to subset of bits.
    using Subset = std::unordered_set<uint32_t>;
    IW_STL_Hash_Map<IWString, Subset> _tag_to_bit_subset;

    int _flatten_sparse_counted;

    // Once the .gfp environment is set up, masks for the fixed width
    // fingerprints.
    int _nfixed;
    IWDYFP* _fixed_mask;

    // Once the .gfp environment is set up, pointers to the Subset objects
    // in _bit_subset.
    int _nsparse;
    Subset** _sparse_subset;

  public:
    BitSubset();
    ~BitSubset();

    // Initialise from a GfpBitSubset proto.
    int Build(IWString& fname);
    // Initialise from proto.
    int Build(const GfpBitSubset::GfpBitSubset& proto);

    // private functions.
    // Once the GFP environment is set up, initialise.
    // Fails if any fingerprints in _bit_to_feature are not present.
    int InitialiseGfpKnown(const IW_General_Fingerprint& gfp);

    // Has this been initialized.
    int Active() const { 
      return !_tag_to_bit_subset.empty();
    }

    // Remove bits that are not retained from `gfp`
    // THe only reason it is non const is that on the first call it will
    // initialise internal structures.
    // If that fails, it will return -1.
    // Otherwise returns the number of features now in `gfp`.
    int MakeSubset(IW_General_Fingerprint& gfp);
};

// The BitXref class is used to support conversion of .gfp data
// to .svml files, where each bit must be mapped to a feature number.
// This may be instantiated before the GFP environment is established.
class BitXref {
  private:
    using BitToFeature = std::unordered_map<gfp_bit_type, uint32_t>;

    // Mapping from tag to feature mapping.
    IW_STL_Hash_Map<IWString, BitToFeature> _tag_to_bit_to_feature;

    // Once the gfp environment is set up, we can form an array of
    // pointers to the BitToFeature for each dense fingerprint.
    // Same for sparse fingerprints
    int _nfixed;
    BitToFeature** _fixed;
    int _nsparse;
    BitToFeature** _sparse;

    // when doing tsv output, we need to know the highest feature number.
    uint32_t _highest_feature_number;

    // Read from the proto.
    int _flatten_sparse_counted;

    // private functions.
  public:
    BitXref();
    ~BitXref();

    // Initialisation from a file containing a GfpBitSubset::GfpBitToFeature proto.
    int Build(IWString& fname);
    // Initialisation from proto.
    int Build(const GfpBitSubset::GfpBitToFeature& proto);

    // Once the GFP environment is set up, initialise.
    // Fails if any fingerprints in _bit_to_feature are not present.
    int InitialiseGfpKnown(const IW_General_Fingerprint& gfp);

    // Across all the BitToFeature hashes, the highest value.
    uint32_t HighestFeatureNumber() const;

    // Convert bit numbers in `gfp` to feature numbers in `features`.
    // Returns the number of items set in `features`.
    template <typename T>
    int PopulateFeatureVector(const IW_General_Fingerprint& gfp,
                               std::vector<T>& features) const;

    // Has this been initialized.
    int Active() const { 
      return _tag_to_bit_to_feature.size() > 0;
    }

    // Write all the feature:count pairs to `output`.
    // Note that the first token (class/activity) is not written, and
    // neither is the #... Just the fetures.
    // THe only reason it is non const is that on the first call it will
    // initialise internal structures.
    // If that fails, it will return -1.
    // Otherwise returns the number of features written.
    // Returning 0 just means that none of the bits in `gfp` were written.
    int WriteSvmlFeatures(const IW_General_Fingerprint& gfp, IWString_and_File_Descriptor& output);

    // Write the bits in `gfp` in tabular form to `output`.
    // The number of columns is _highest_feature_number.
    // Note there is no leading `output_separator` written. No newline.
    int WriteDsv(const IW_General_Fingerprint& gfp, char output_separator,
                 IWString_and_File_Descriptor& output);
};

}  // namespace bit_subset

#endif // UTILITIES_GFP_TOOLS_BIG_SUBSET_H
