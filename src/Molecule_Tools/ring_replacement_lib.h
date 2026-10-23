#ifndef MOLECULE_TOOLS_RING_REPLACEMENT_LIB_H_
#define MOLECULE_TOOLS_RING_REPLACEMENT_LIB_H_

#include <cstdint>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/replacement_ring.pb.h"
#else
#include "replacement_ring.pb.h"
#endif

namespace ring_replacement {

void Usage(int rc);

class Replacement;

class RingReplacement {
  private:
    int _verbose;

    uint64_t _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    int _write_parent_molecule;

    Atom_Typing_Specification _atype;

    int _preserve_connectivity;

    // One per -R option.
    int _nreplacements;
    resizable_array_p<Replacement>* _rings;

    // Identify the (ring) atoms in the starting molecule that
    // are available for replacement.
    resizable_array_p<Substructure_Query> _queries;

    // Sometimes it might be useful to just skip molecules not matching
    // any of the queries.
    int _ignore_molecules_not_matching_queries;
    // By default, molecules that do not match any queries, and are
    // ignored, are not written.
    int _write_non_matching_molecules;

    resizable_array_p<Substructure_Query> _products_must_have;
    resizable_array_p<Substructure_Query> _products_must_not_have;

    uint64_t _molecules_discarded_by_must_have_queries = 0;
    uint64_t _molecules_discarded_by_must_not_have_queries = 0;

    // As the replacement rings are read, we can filter to a subset of them.
    // From the -F and -D options.
    resizable_array_p<Substructure_Query> _replacement_rings_must_have;
    resizable_array_p<Substructure_Query> _replacement_rings_must_not_have;

    uint64_t _molecules_with_invalid_valences_suppressed;

    int _unique_molecules_only = 0;

    int _min_examples_needed = 0;
    int _replacement_rings_discarded_for_count = 0;

    // Specified via the -f option. The maximum formula difference between
    // what is removed and the replacement ring.
    std::optional<uint32_t> _max_formula_difference;

    IW_STL_Hash_Set _seen;
    int _duplicates_suppressed;

    uint64_t _molecules_formed;

    extending_resizable_array<int> _number_variants;

    // Remove isotopes from product molecules.
    int _remove_isotopes;

    // Translate all isotopes to a single value.
    isotope_t _translate_isotopes;

    int _sort_output_by_precedent = 0;

    Chemical_Standardisation _chemical_standardisation;

    IWString_and_File_Descriptor _stream_for_unchanged_molecules;

  // Private functions.

    int ReadReplacementRings(IWString& fname, int ndx);
    int ReadReplacementRings(iwstring_data_source& input, int ndx);

    int ReadFileOfReplacementRings(IWString& fname, int ndx);
    int ReadFileOfReplacementRings(const IWString& dirname,
                        int ndx, iwstring_data_source& input);
    int ReadSetOfReplacementRings(const IWString& dirname,
                IWString& fname, int ndx);

    int OkSupport(const Replacement& r);
    int IsUnique(Molecule& m);
    int OkProductQueryConstraints(Molecule& m);
    // These are applied to replacement ring candidates.
    int OkWithQueryRequirements(Molecule_to_Match& target);
    int OkWithQueryRequirements(Replacement& r,
                const RplRing::ReplacementRing& proto);

    int Process(resizable_array_p<Molecule>& mols, int ndx, const uint32_t* atypes,
                const int* process_atom);
    int Write(const resizable_array_p<Molecule>& mols, IWString_and_File_Descriptor& output);
    int IdentifyMatchedAtoms(Molecule& m, resizable_array_p<Substructure_Query>& queries,
                             int* process_atom);

  public:
    RingReplacement();
    ~RingReplacement();

    int Initialise(Command_Line& cl);

    // Mostly for python bindings.
    bool set_ring_atom_smarts(const std::string& smarts);
    void set_unique_molecules_only(bool s) {
      _unique_molecules_only = s;
    }
    void clear_unique_molecule_cache() {
      _seen.clear();
    }
    void set_min_support_requirement(uint32_t s);
    void set_max_formula_difference(uint32_t s); 
    // The python binding only supports one set of replacement rings.
    // Because the _rings array is an array of objects and I don't want to resize it.
    uint32_t ReadReplacementRings(const std::string& fname);
    uint32_t number_replacement_rings() const {
      return _rings[0].size();
    }
    void set_remove_isotopes(bool s) {
      _remove_isotopes = s;
    }
    void set_preserve_connectivity(bool s);
    // This is what does the work.
    std::vector<Molecule> Process(Molecule& m);
    // end python related.

    int Preprocess(Molecule& m);

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output);
};

void set_display_valence_error_messages(int s);
void set_warn_loss_of_aromaticity(int s);

}  // namespace ring_replacement

#endif
