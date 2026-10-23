#include "absl/container/flat_hash_map.h"

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/absl_hash.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/gfp_standard.h"

#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#ifdef BUILD_BAZEL
#include "Retrosynthesis/retrosynthesis.pb.h"
#else
#include "retrosynthesis.pb.h"
#endif

namespace retrosynthesis {

// The overall class which parses command lines may have some structure
// altering directives. These need to be passed down to all the objects
// that get created by it.
class Preprocessor {
  private:
    int _remove_chirality;

    Chemical_Standardisation _chemical_standardisation;
  public:
    Preprocessor();

    // Only the -c and -g options are recognised.
    int Initialise(Command_Line& cl);

    int Process(Molecule& m);
};

// Molecules need to be decomposed.
class MoleculeAndFragmentation : public Molecule {
  private:
    MoleculeAndFragmentation* _parent;

    // THe _found molecules will have the ID of the reagent.
    resizable_array_p<MoleculeAndFragmentation> _found;

    resizable_array_p<MoleculeAndFragmentation> _notfound;

    // We use the user_specified_void_ptr of each atom so that fragments
    // know their origins.
    int* _atom_numbers;

    // When dealing with products, each product must know the reaction uinque id
    // that generated it.

    int _generated_by;

  public:
    MoleculeAndFragmentation();
    ~MoleculeAndFragmentation();

    int set_parent_molecule(const Molecule& m);

    int DebugPrint(std::ostream& output) const;
    int DebugPrint(const IWString& indentation, std::ostream& output);

    bool empty() const;

    resizable_array_p<MoleculeAndFragmentation>& found() {
      return _found;
    }
    resizable_array_p<MoleculeAndFragmentation>& notfound() {
      return _notfound;
    }

    uint32_t number_found() const {
      return _found.size();
    }
    uint32_t number_not_found() const {
      return _notfound.size();
    }

    void set_generated_by(int s) {
      _generated_by = s;
    }

    int generated_by() const {
      return _generated_by;
    }

    // Assume ownership of `m` in the _found array.
    void BeenFound(MoleculeAndFragmentation* m) {
      _found << m;
    }
    void NotFound(MoleculeAndFragmentation* m) {
      _notfound << m;
    }

    // recursively look at the _found arrays and figure out the maximum
    // number of parent atoms covered.
    int BestCoverage() const;
};

std::ostream&
operator<<(std::ostream& output, MoleculeAndFragmentation& mf);

class SetOfFingerprints {
  private:
    IW_General_Fingerprint* _fp;
    int _number_fingerprints;
    isotope_t _isotope;

    // a mapping from the unique smiles of a building block to its name.
    absl::flat_hash_map<IWString, IWString> _usmi_to_id;

    // We keep track of the smallest and largest building block.
    // We can avoid searches when a query molecule is out of range.
    int _min_natoms;
    int _max_natoms;

    IWString _name;

  // Private functions
    int ReadFingerprints(iwstring_data_source& input,
        Preprocessor& preprocess,
        bool proto_has_isotope);
    int ReadFingerprints(const IWString& dirname, const std::string& fname,
                Preprocessor& preprocess, bool proto_has_isotope);

  public:
    SetOfFingerprints();
    ~SetOfFingerprints();

    int Build(const IWString& dirname, const FingerprintData& proto, Preprocessor& preprocess);
    int Build(iwstring_data_source& input, Preprocessor& preprocess);

    void set_name(const IWString& s) {
      _name = s;
    }
    const IWString& name() const {
      return _name;
    }

    isotope_t isotope() const {
      return _isotope;
    }

    // If the m.unique_smiles() is in _usmi_to_id, set m.name();
    int InDatabase(Molecule& m);
};

// A reversed reaction generates 
class Reaction {
  private:
    // Each reaction is assigned a unique id.
    int _unique_id;

    // Perhaps it will be necessary to have multiple reactions to
    // handle a specific set of building blocks - boronic...
    IWReaction* _reaction;
    int _number_reactions;

    // The name of the first _reaction. It is an error to NOT have a name.
    IWString _name;

    // All molecules in the set will have an atom with this isotope.
    isotope_t _isotope;

    // A mapping from isotope number to a set of fingerprints.
    absl::flat_hash_map<int, SetOfFingerprints*> _iso_to_fp;

  // private functions
//  int ReadFingerprints(const IWString& dirname, const std::string& fname, Preprocessor& preprocess);
//  int ReadFingerprints(iwstring_data_source & input, Preprocessor& preprocess);
    int ReadFingerprints(const IWString& dirname, const FingerprintData& proto, Preprocessor& preprocess);

    int ProcessMatch(Molecule& m, MoleculeAndFragmentation& results);
    int InDatabase(Molecule& m);

  public:
    Reaction();
    ~Reaction();

    int Build(const ReactionData& proto, Preprocessor& preprocess, const IWString& dirname);

    const IWString& name() const {
      return _name;
    }

    void set_unique_id(int s) {
      _unique_id = s;
    }
    int unique_id() const {
      return _unique_id;
    }

    int Process(Molecule& m, MoleculeAndFragmentation& result);
};

class Retrosynthesis {
  private:
    Reaction* _reaction;
    int _number_reactions;

    int _max_number_steps;

    int _verbose;

  // Private functions
    int Build(const RetrosynthesisData& proto, Preprocessor& preprocess, const IWString& dirname);

    int ReadFingerprints(IWString& fname);
    int ReadFingerprints(iwstring_data_source& input);

    int Process2(Molecule& m, int recursion, MoleculeAndFragmentation& results);

  public:
    Retrosynthesis();
    ~Retrosynthesis();

    int Initialise(Command_Line& cl, Preprocessor& preprocess);

    int number_reactions() const {
      return _number_reactions;
    }

    int Process(int rxn, Molecule& m, MoleculeAndFragmentation& results);
};  // class Retrosynthesis

}  // namespace retrosynthesis
