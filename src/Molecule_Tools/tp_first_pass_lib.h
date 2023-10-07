// Support functions for the first phase of the Medchem Rules

#include <limits>
#include <unordered_map>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rmele.h"

namespace lilly_medchem_rules {

// Structure to hold either a count or a fraction.
// We do no checking to guard against both being set.
struct CountOrFraction {
  int count = -1;
  float fraction = 0.0f;

  // Construct from something that is either (\d+) or (0\.\d+),
  // setting either count or fraction depending on what matches.
  // Returns non zero if successful.
  int Build(const const_IWSubstring& token);

  // Returns true if `c` is either > count if set, or
  // > fraction if set.
  bool GreaterThan(int c, int tot) const;
};

// Counts of the number of molecules discarded for various reasons
struct MCFirstPassCounter {
  int molecules_processed = 0;
  int molecules_rejected = 0;
  int empty_molecule = 0;
  int molecules_below_atom_count_cutoff = 0;
  int molecules_above_atom_count_cutoff = 0;
  int molecules_with_fragment_too_large = 0;
  int mixtures_rejected = 0;
  int molecules_with_no_interesting_atoms = 0;
  int molecules_with_too_few_interesting_atoms = 0;
  int molecules_containing_non_allowed_atom_types = 0;
  int molecules_containing_non_periodic_table_elements = 0;
  int molecules_containing_colvalent_non_organics = 0;
  int molecules_with_too_few_rings = 0;
  int molecules_with_too_many_rings = 0;
  int molecules_containing_isotopes = 0;
  int molecules_with_isotopes_converted = 0;
  int molecules_with_abnormal_valences = 0;
  int molecules_with_two_connected_hydrogens = 0;
  int molecules_with_bad_max_ring_bond_ratios = 0;
  int molecules_with_bad_ring_system_size = 0;
  int molecules_with_ring_sizes_out_of_range = 0;
};

// A class that does the first phase of processing.
class MCFirstPass {
  private:
    // Normally processing is silent. Incrementing verbose will increase
    // the number of diagnostic messages.
    int _verbose = 0;

    // By default. non periodic table elements are excluded.
    int _exclude_molecules_containing_non_periodic_table_elements = 1;
 
    // It might be OK to have a non periodic table element if it is
    // not connected to anything.
    int _allow_non_periodic_table_elements_if_not_connected = 0;

    // Range of allowed atom counts. Note that this is applied on
    // a per fragment basis.
    int _lower_atom_count_cutoff = 0;
    int _upper_atom_count_cutoff = 0;

    // A quick way to check which elements are OK. Those elements that are
    // permitted will have their value in the array set to 1.
    int _ok_elements[HIGHEST_ATOMIC_NUMBER + 1];

    float _min_fraction_interesting_atoms = 0.0f;

    int _reject_if_any_fragment_larger_than = std::numeric_limits<int>::max();

    int _exclude_molecules_with_no_interesting_atoms = 0;

    int _lower_ring_count_cutoff = 0;
    int _upper_ring_count_cutoff = 0;

    int _exclude_isotopes = 1;

    int _convert_isotopes = 0;

    int _skip_molecules_with_abnormal_valences = 1;

    int _upper_ring_size_cutoff = 0;

    // Constraint on the number of rings in a system.
    int _upper_ring_system_size = std::numeric_limits<int>::max();

    // Initialised out of range.
    float _max_ring_bond_ratio = 2.0f;

    // Remove or translate non printing characters in the name.
    int _remove_non_printing_chars = 0;
    // Non printing characters are converted to this if specified.
    char _translate_non_printing_chars = '\0';

    // Applied near the start of processing.
    Elements_to_Remove _elements_to_remove;
    Element_Transformations _element_transformations;

    // A map from atomic numbers to rejection criteria.
    std::unordered_map<atomic_number_t, CountOrFraction> _element_fraction;

    // Private functions
    // Most of these implement a specific rule.

    // Examine m.name() and depending on the values of
    // params.remove_non_printing_chars and 
    // params.translate_non_printing_chars
    // alter the name.
    // Returns the number of non printing characters initially in m.name().
    int NonPrintingCharactersInName(Molecule & m) const;

    // Return true if any of the atoms in `m` are not consistent with
    // our element specifications.
    bool ExcludeForAtomType(const Molecule& m,
                             MCFirstPassCounter& counter,
                             IWString& rejection_reason) const;

    // Returns true if `m` has below params.min_fraction_interesting_atoms
    // interesting atoms in the largest fragment.
    int ExcludeForTooFewInterestingAtoms(Molecule & m,
                                 MCFirstPassCounter& counter,
                                 IWString& rejection_reason) const;

    // Returns true if `m` has no interesting atoms in the largest
    // fragment.
    int ExcludeForNoInterestingAtoms(Molecule & m,
                             MCFirstPassCounter& counter,
                             IWString& rejection_reason) const;
    // Returns 1 if the fragments in `m` are consistent with
    // our settings. If not, `counter` will be updated and
    // `rejection_reason1` will be set.
    int FragmentsAreOK(Molecule & m,
                        MCFirstPassCounter& counter,
                        IWString& rejection_reason) const;

    int OkLowerUpperAtomCountCutoff(const int natoms,
                                         MCFirstPassCounter& counter,
                                         IWString& rejection_reason) const;

    int ConvertIsotopes(Molecule & m, MCFirstPassCounter& counter) const;

    int RejectForRingSizeCondition(Molecule & m,
                               MCFirstPassCounter& counter,
                               IWString& rejection_reason) const;
    int RejectedForRingSystemSize(Molecule& m,
                MCFirstPassCounter& counter,
                IWString& rejection_reason) const;

    int RejectForRingBondRatio(Molecule & m,
                        MCFirstPassCounter& counter,
                        IWString& rejection_reason) const;

    int RejectedForElementFraction(Molecule& m,
                MCFirstPassCounter& count,
                IWString& rejection_reason) const;

    int RejectedInner(Molecule& m,
                      MCFirstPassCounter& counter,
                      IWString& rejection_reason);
  public:
    MCFirstPass();

    int Build(Command_Line& cl);

    // Examine `m` and return true if it is rejected for any reason.
    // If rejected, the reason will be added to `counter` and 
    // `rejection_reason` will be filled.
    int Rejected(Molecule& m, MCFirstPassCounter& counter, IWString& rejection_reason);

    int Report(const MCFirstPassCounter& counter, std::ostream& output) const;

    // Setters for all the parameters.
    void set_verbose(int v) { _verbose = v;}
    int  verbose() const { return _verbose;}
    void set_exclude_molecules_containing_non_periodic_table_elements(int s) {
      _exclude_molecules_containing_non_periodic_table_elements = s;
    }
    void set_allow_non_periodic_table_elements_if_not_connected(int s) {
      _allow_non_periodic_table_elements_if_not_connected = s;
    }
    void set_lower_atom_count_cutoff(int s) {
      _lower_atom_count_cutoff = s;
    }
    void set_upper_atom_count_cutoff(int s) {
      _upper_atom_count_cutoff = s;
    }
    void set_min_fraction_interesting_atoms(float s) {
      _min_fraction_interesting_atoms = s;
    }
    void set_reject_if_any_fragment_larger_than(int s) {
      _reject_if_any_fragment_larger_than = s;
    }
    void set_exclude_molecules_with_no_interesting_atoms(int s) {
      _exclude_molecules_with_no_interesting_atoms = s;
    }
    void set_lower_ring_count_cutoff(int s) {
      _lower_ring_count_cutoff = s;
    }
    void set_upper_ring_count_cutoff(int s) {
      _upper_ring_count_cutoff = s;
    }
    void set_exclude_isotopes(int s) {
      _exclude_isotopes = s;
    }
    void set_convert_isotopes(int s) {
      _exclude_isotopes = 0;
      _convert_isotopes = s;
    }
    void set_skip_molecules_with_abnormal_valences(int s) {
      _skip_molecules_with_abnormal_valences = s;
    }
    void set_upper_ring_size_cutoff(int s) {
      _upper_ring_size_cutoff = s;
    }
    void set_max_ring_bond_ratio(float s) {
      _max_ring_bond_ratio = s;
    }
    void set_remove_non_printing_chars(int s) {
      _remove_non_printing_chars = s;
    }
    void set_translate_non_printing_chars(char s) {
      _translate_non_printing_chars = s;
    }
    int set_ok_element(const Element * e);
};

// Exposed just for testing.
int InterestingAtoms(const Molecule & m);
}  // namespace lilly_medchem_rules
