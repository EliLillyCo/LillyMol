#ifndef MOLECULE_TOOLS_MULTI_ARM_FINGERPRINT_H_
#define MOLECULE_TOOLS_MULTI_ARM_FINGERPRINT_H_

#include <cstdint>

#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/molecule.h"

namespace multi_arm_fingerprint {

struct Options {
  int _min_total_length = 3;
  int _max_total_length = 6;
  // Zero means bounded only by max_total_length and the number of arms.
  int _max_arm_length = 0;
  // Zero enables both three- and four-connected centres.
  int _centre_degree = 0;
  bool _ring_closures = true;
  bool _presence_only = false;
  // Expand all arms at the same radius, retaining shorter terminal paths.
  bool _synchronised = false;
  bool Valid() const;
};

struct Statistics {
  uint64_t _centres = 0;
  uint64_t _arm_paths = 0;
  uint64_t _combinations = 0;
  uint64_t _distinct_bits = 0;
};

// Reusable generator: no command-line, query, preprocessing, or output state.
// Atom typing can be configured here, or supplied by the caller. Molecules may
// acquire aromaticity caches during generation. Separate instances are needed
// when using the internal atom typer concurrently.
class Generator {
 private:
  Options _options;
  Atom_Typing_Specification _atom_typing;

 public:
  Generator();
  // Invalid options leave the existing configuration unchanged.
  bool SetOptions(const Options& options);

  const Options&
  options() const {
    return _options;
  }

  int SetAtomTyping(const const_IWSubstring& specification);

  // A null roots mask means all eligible centres. Otherwise only nonzero entries
  // are roots; arms may traverse all heavy atoms. Output is replaced, not added
  // to. Molecules containing explicit hydrogen atoms are rejected.
  // Returns zero on failure, leaving an empty fingerprint.
  int Fingerprint(Molecule& m, Sparse_Fingerprint_Creator& result,
                  const int* roots = nullptr, Statistics* statistics = nullptr);
  int Fingerprint(Molecule& m, const uint32_t* atom_types,
                  Sparse_Fingerprint_Creator& result, const int* roots = nullptr,
                  Statistics* statistics = nullptr) const;
};

}  // namespace multi_arm_fingerprint
#endif  // MOLECULE_TOOLS_MULTI_ARM_FINGERPRINT_H_
