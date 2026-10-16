#ifndef MOLECULE_TOOLS_MEDCHEMWIZARD_LIB_H_
#define MOLECULE_TOOLS_MEDCHEMWIZARD_LIB_H_

#include <cstdint>
#include <functional>
#include <limits>
#include <memory>

#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"

namespace medchemwizard {

// Called once for each generated product that passes filtering. The callback can leave
// ownership with the library, or release the unique_ptr into a longer-lived container.
// If recursive generation is enabled, released products must remain alive until
// Process returns; the built-in ProcessToArray and ProcessToStream adapters satisfy this.
using ProductCallback = std::function<int(std::unique_ptr<Molecule>&)>;

struct Options {
  int verbose = 0;
  int reduce_to_largest_fragment = 0;
  int max_depth = 0;
  int max_hits = std::numeric_limits<int>::max();
  int report_multiple_hits_truncated = 1;
  int unique_across_all_molecules = 0;
  int unique_within_molecule = 0;
  int postprocess_molecules_produced = 0;
  int discard_bad_valences = 0;
  int max_atoms = std::numeric_limits<int>::max();
  int min_atoms = 1;
  int remove_all_chiral_centres = 0;
  int append_names = 0;
  int ignore_do_not_change_queries_not_matching = 0;
  IWString sep;
  Chemical_Standardisation chemical_standardisation;
};

struct Stats {
  uint64_t molecules_read = 0;
  uint64_t molecules_produced = 0;
  uint64_t truncated_to_max_hits = 0;
  uint64_t duplicate_molecules_suppressed = 0;
  uint64_t bad_valences_discarded = 0;
  uint64_t discarded_for_too_many_atoms = 0;
  uint64_t discarded_for_too_few_atoms = 0;
  uint64_t molecules_not_matching_do_not_change_queries = 0;
  uint64_t embeddings_rejected_for_changing_protected_atoms = 0;
};

class MedchemWizard {
 private:
  Options _options;
  Stats _stats;
  resizable_array_p<IWReaction> _rxn;
  resizable_array_p<Substructure_Query> _do_not_change_query;
  std::unique_ptr<int[]> _molecules_hitting_reaction;
  IW_STL_Hash_Set _smiles_hash;
  IWString_and_File_Descriptor* _stream_for_bad_valence = nullptr;

  void Preprocess(Molecule& m);
  int IdentifyAtomsThatMustNotChange(Molecule& m, std::unique_ptr<int[]>& do_not_change);
  int MaybeTruncateNhits(int nhits, IWReaction& reaction, Molecule& m);
  int CheckUniqueness(Molecule& m, int query_just_run, IW_STL_Hash_Set& smiles_hash);
  int GenerateProducts(Molecule& m, const IWString& mname, int depth, const int* do_not_change,
                       IW_STL_Hash_Set& smiles_hash, ProductCallback& callback);

 public:
  MedchemWizard();

  Options& options() {
    return _options;
  }
  const Options& options() const {
    return _options;
  }

  void set_bad_valence_stream(IWString_and_File_Descriptor* stream) {
    _stream_for_bad_valence = stream;
  }

  resizable_array_p<Substructure_Query>& do_not_change_queries() {
    return _do_not_change_query;
  }
  const resizable_array_p<Substructure_Query>& do_not_change_queries() const {
    return _do_not_change_query;
  }

  int ReadReactions(const char* fname);
  // Append a file, resolving resources beside it. Failed reads leave the set unchanged.
  int ReadReactionFile(const char* fname, bool textproto);
  int InitialiseFromEnvironment();
  int number_reactions() const {
    return _rxn.number_elements();
  }

  int Process(Molecule& m, ProductCallback callback);
  int ProcessToArray(Molecule& m, resizable_array_p<Molecule>& output);
  int ProcessToStream(Molecule& m, IWString_and_File_Descriptor& output);

  const Stats& stats() const {
    return _stats;
  }
  const int* molecules_hitting_reaction() const {
    return _molecules_hitting_reaction.get();
  }
  const IWReaction* reaction(int ndx) const {
    return _rxn[ndx];
  }
};

}  // namespace medchemwizard

#endif  // MOLECULE_TOOLS_MEDCHEMWIZARD_LIB_H_
