#ifndef MOLECULE_TOOLS_UNIQUE_MOLECULES_H_
#define MOLECULE_TOOLS_UNIQUE_MOLECULES_H_

// This was initially built as a demo for the python api.
// And while it does work, it is significantly slower than using a set() in
// python! Not sure I understand that...

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include "absl/container/btree_set.h"
#include "absl/container/flat_hash_set.h"

#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/etrans.h"

namespace unique_molecules {

class UniqueMolecules {
  private:
    int _exclude_chiral_info;
    int _exclude_cis_trans_bonding_info;
    int _strip_to_largest_fragment;
    int _ignore_isotopes;
    isotope_t _constant_isotope;

    Mol2Graph _mol2graph;

    Element_Transformations _etrans;

    //resizable_array_p<IWString_STL_Hash_Set> _smiles_hash;
    resizable_array_p<absl::flat_hash_set<std::string>> _smiles_hash;

    int _use_atom_hash;
    resizable_array_p<std::unordered_set<uint64_t>> _atom_hash;

    uint32_t _molecules_processed;
    uint32_t _duplicates_found;

  // Private functions.

    bool InternalIsUnique(Molecule&& m);
    bool InternalIsUniqueInner(Molecule& m);

  public:
    UniqueMolecules();

    // Setters for matching conditions.
    void set_exclude_chiral_info(int s) {
      _exclude_chiral_info = s;
    }
    void set_exclude_cis_trans_bonding_info(int s) {
      _exclude_cis_trans_bonding_info = s;
    }
    void set_strip_to_largest_fragment(int s) {
      _strip_to_largest_fragment = s;
    }
    void set_ignore_isotopes(int s) {
      _ignore_isotopes = s;
    }
    void set_constant_isotope(isotope_t s) {
      _constant_isotope = s;
    }
    // Fully expose these so they can be initialised externally.
    Mol2Graph& graph_specifications() {
      return _mol2graph;
    }
    Element_Transformations& element_transformations() {
      return _etrans;
    }

    // For pre-populating the hashes. Does not update the counters.
    void AddToHashes(const Molecule& m);

    // True if `m` has not been seen before.
    // Note that `m` is not changed, internally a copy is made,
    // and any transformations applied to the copy.
    bool IsUnique(const Molecule& m);

    int Report(std::ostream& output) const;
};

}  // namespace unique_molecules

#endif // MOLECULE_TOOLS_UNIQUE_MOLECULES_H_
