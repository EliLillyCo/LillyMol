#ifndef MOLECULE_TOOLS_FILECONV_OPTS_H
#define MOLECULE_TOOLS_FILECONV_OPTS_H

#include <iostream>
#include <limits>
#include <memory>
#include <unordered_map>

#include "Foundational/accumulator/accumulator.h"

#include "Molecule_Lib/allowed_elements.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/numass.h"
#include "Molecule_Lib/rmele.h"

#include "fix_structures.h"
#include "known_fragment_data.h"

namespace fileconv {

// FileconvConfig::Process() uses a FileconvResult struct to
// pass data back to the caller.
struct FileconvResult {
  // The molecule encountered an error, the state is unclear.
  int error = 0;

  // Is the molecule rejected.
  int rejected = 0;
  // If rejected, the reason for that rejection.
  IWString rejection_reason;

  // Whether or not the molecule is changed by Process(m).
  // The number of different changes applied.
  int molecule_changed = 0;

  // Before filters were applied, the state of the molecule.
  // Will always be filled, even if the molecule is rejected.
  // This is to support part of the selimsteg processing, kind
  // of a kludge.
  IWString smiles_before_filters;
};

// Data structure defining fileconv invocations.
// Should be a class, but is a struct to avoid having to implement a large
// number of setter methods.
struct FileconvConfig {
  int audit_input = 0;

  int verbose = 0;

  int debug_print_each_molecule = 0;

  // Geometric values can be printed.
  int print_bond_lengths = 0;
  int print_bond_angles = 0;
  int print_torsions = 0;
  int print_max_atom_separation = 0;

  // If we are printing the max atom separation for each molecule
  // may as well form summary stats across the entire input.
  Accumulator<double> acc_longest_sep;

  // Number of times Process() is called.
  int molecules_processed = 0;

  int molecules_changed = 0;

  // We can append a string to the name of any molecule whose connection table is changed

  IWString molecule_changed_string;

  Chemical_Standardisation chemical_standardisation;

  Charge_Assigner charge_assigner;

  Donor_Acceptor_Assigner donor_acceptor_assigner;

  Structure_Fixing structure_fixing;

  Atom_Typing_Specification atom_typing;

  int fragment_count = 0;
  int reduce_to_largest_fragment = 0;
  int reduce_to_all_largest_fragments = 0;
  int reduce_to_largest_organic_fragment = 0;
  int reduce_to_largest_organic_fragment_carefully = 0;
  bool keep_all_organic_fragments = false;
  int molecules_with_nonorganic_fragments_removed = 0;
  int maximum_fragment_count = 0;
  int min_size_max_fragment_count = 0;
  int molecules_chopped = 0;
  int molecules_with_too_many_components = 0;
  int molecules_with_large_fragments = 0;
  int remove_duplicate_fragments = 0;
  int molecules_with_duplicate_fragments_removed = 0;
  int remove_fragments_this_size_or_smaller = 0;
  int molecules_with_very_small_fragments_removed = 0;
  int remove_largest_fragment = 0;
  int remove_smallest_fragment = 0;
  int keep_smallest_fragments = 0;
  int remove_molecules_with_non_largest_fragment_natoms = -1;
  int strip_to_n_largest_fragments = 0;
  int sort_by_fragment_size = 0;
  IWString tag_for_removed_fragments;

  // Jul 2022. Introduce the idea of a window of atom counts, and keep only
  // fragments that are within range.
  int fragment_selector_lower_atom_count = 0;
  int fragment_selector_upper_atom_count = std::numeric_limits<int>::max();

  // People are registering mixtures. This is very bad for us.
  // If a molecule contains more than 1 fragment with more than
  // a given number of atoms, discard the molecule

  int discard_molecule_if_multiple_fragments_larger_than = 0;

  int molecules_with_multiple_large_fragments = 0;

  // Sometimes if things are ambiguous, just call it a mixture if the two
  // largest fragments differ by a given amount or less

  int mixture_if_largest_frags_differ_by = -1;

  int molecules_declared_mixtures_by_atom_count_difference = 0;

  Known_Fragment_Data known_fragment_data;

  // We can keep track of how many atoms are lost

  Accumulator_Int<int> atoms_lost;
  extending_resizable_array<int> initial_fragment_count;

  // For efficiency we have a single variable which indicates whether or not
  // the fragment code must be called.

  int need_to_call_process_fragments = 0;

  int skip_molecules_with_abnormal_valences = 0;
  int ok_bad_valence_on_isotopically_labelled = 0;
  int molecules_with_abnormal_valences = 0;

  int max_path_length = 0;
  int molecules_with_longest_path_too_long = 0;

  // Will be set if any isotope related directive has been specified.
  int need_to_consider_isotopes = 0;

  int exclude_isotopes = 0;
  int molecules_containing_isotopes = 0;  // only if the isotopic attributes are changed

  int convert_isotopes = 0;

  int convert_all_isotopes_to = 0;

  int remove_isotopic_atoms = 0;
  int isotope_is_atom_type = 0;

  int convert_isotopes_to_atom_map_numbers = 0;
  int convert_atom_map_numbers_to_isotopes = 0;

  resizable_array<int> convert_specific_isotopes;
  resizable_array<int> convert_specific_isotopes_new_isotope;

  resizable_array_p<Substructure_Query> convert_specific_isotopes_query;
  resizable_array<int> convert_specific_isotopes_query_new_isotope;

  // For a given isotope number, a fragment that is added to that
  // atom. Could also do a reaction, but this might be faster.
  // Adds via the first atom in the fragment.
  std::unordered_map<int, Molecule> add_to_isotopic_atom;

  int output_organic_only = 0;
  int non_organic_molecules_found = 0;

  int exclude_non_real_elements = 0;
  int non_real_elements_found = 0;

  int molecules_excluded_for_non_allowed_elements = 0;

  double lower_molecular_weight_cutoff = 0.0;
  double upper_molecular_weight_cutoff = 0.0;

  int molecules_below_molecular_weight_cutoff = 0;
  int molecules_above_molecular_weight_cutoff = 0;

  Accumulator<molecular_weight_t> amw_accumulator;
  Accumulator_Int<int> natoms_accumulator;
  extending_resizable_array<int> atom_count;

  int lower_atom_count_cutoff = 0;
  int upper_atom_count_cutoff = 0;

  /*
    I never implemented these two, they are mostly taken care
    of by atom_count_includes_only_atoms_in_largest_fragment
    If that ever becomes a problem, go ahead and implement them

    static int lower_atom_count_cutoff_applies_to_largest_fragment = 0;
    static int upper_atom_count_cutoff_applies_to_largest_fragment = 0;
  */

  int include_implicit_hydrogens_in_upper_atom_count_comparison = 0;

  int atom_count_includes_only_atoms_in_largest_fragment = 0;

  int molecules_below_atom_count_cutoff = 0;
  int molecules_above_atom_count_cutoff = 0;

  // Only write molecules where the name matches a regular expression.
  std::unique_ptr<RE2> name_rx;
  // Discard molecules whose name matches a regular expression.
  std::unique_ptr<RE2> grep_v_name_rx;
  int molecules_discarded_for_name_mismatch = 0;

  std::unique_ptr<RE2> change_name_rx;
  int change_name_rx_must_match = 0;

  // If set, names are left padded with zero's until they are `zero_pad` wide.
  // And zero_pad_name can be negative, in which case the first `zero_pad_name`
  // characters are removed from the name.
  int zero_pad_name = 0;

  // When reducing fragment count, we can keep the smallest or largest
  // fragment which matches a given query.

  resizable_array_p<Substructure_Query> smallest_fragment_queries;
  resizable_array_p<Substructure_Query> largest_fragment_queries;
  resizable_array_p<Substructure_Query> keep_fragment_queries;
  resizable_array_p<Substructure_Query> remove_fragment_queries;

  int molecules_with_fragments_reduced_by_query = 0;
  int molecules_not_matching_fragment_queries = 0;

  int lower_ring_count_cutoff = 0;
  int molecules_with_too_few_rings = 0;
  int upper_ring_count_cutoff = -1;
  int molecules_with_too_many_rings = 0;

  int min_ring_systems_required = -1;
  int max_ring_systems_allowed = -1;
  int ring_systems_include_spiro = 0;
  int molecules_with_too_few_ring_systems = 0;
  int molecules_with_too_many_ring_systems = 0;

  int min_aliphatic_ring_count = 0;
  int max_aliphatic_ring_count = -1;
  int min_aromatic_ring_count = 0;
  int max_aromatic_ring_count = -1;

  int molecules_with_too_few_aliphatic_rings;
  int molecules_with_too_few_aromatic_rings;
  int molecules_with_too_many_aliphatic_rings;
  int molecules_with_too_many_aromatic_rings;

  int upper_ring_size_cutoff = 0;
  int molecules_with_large_rings = 0;

  int max_rings_in_a_ring_system = -1;
  int molecules_with_ring_systems_too_large = 0;

  // We can translate molecules if we wish

  coord_t dx = 0.0;
  coord_t dy = 0.0;
  coord_t dz = 0.0;

  // We can speed things up by having a single variable that indicates whether a translation
  // is needed

  int translation_specified = 0;
  int translate_to_origin = 0;

  // We can rotate the molecule
  Coordinates rotation;
  angle_t rotation_angle = 0.0f;
  int rotation_specified = 0;

  // What kind of partial charge do we want

  int fileconv_partial_charge_type = 0;

  enum ChargeType {
    kNoChargeCalculation,
    kGasteiger,
    kGasteigerHuckel,
    kAbraham,
  };

  // We can optionally assign each molecule written a number R(number).

  Number_Assigner number_assigner;

  Elements_to_Remove elements_to_remove;

  // We can transform element types.

  Element_Transformations element_transformations;

  // With the -H option, we make implicit hydrogens explicit
  int make_all_implicit_hydrogens_explicit = 0;

  // We can add explicit hydrogen atoms to specified atoms.
  resizable_array_p<Substructure_Query> atoms_for_implicit_hydrogens;

  // A count of the number of molecules to which hydrogens are added.
  int molecules_to_which_hydrogens_were_added = 0;

  // If this contains a non zero atomic number, then those elements will
  // be moved to the end of the connection table. This is useful when
  // writing .sdf files.
  int hydrogens_last = 0;

  // Various things for chirality and stereo centres

  int find_all_chiral_centres = 0;
  int find_all_ring_chiral_centres = 0;
  int invert_all_chiral_centres = 0;
  int reflect_coordinates = 0;
  int chiral_centres_inverted = 0;
  int remove_invalid_chiral_centres = 0;
  int molecules_with_invalid_chiral_centres = 0;
  int molecules_with_chiral_centres = 0;
  int remove_chiral_data_from_all_molecules = 0;
  int remove_non_organic_chirality = 0;
  int max_chiral_centres = 0;
  int molecules_with_too_many_chiral_centres = 0;
  extending_resizable_array<int> chiral_centre_count;

  resizable_array_p<Substructure_Query> remove_chiral_centres_on;

  int chiral_centres_removed_by_rmchiral = 0;
  int molecules_with_chiral_centres_removed_by_rmchiral = 0;

  int remove_directional_bonds_from_input = 0;

  int remove_invalid_directional_bonds_from_input = 0;

  int change_double_bonds_between_permanent_aromatic_to_single = 0;

  // With the -O (organic) switch, we can specify a list of non-organic
  // elements which are in fact OK

  Set_of_Element_Matches ok_non_organics;

  // Dec 2008. I want to be able to reject anything that contains certain elements
  // Enabled by '-O def'

  int filter_for_disallowed_elements = 0;

  Allowed_Elements allowed_elements;

  // we can allow O-B-O type Boron atoms as a special case. We could
  // create a general mechanism for this, but today I do not forsee other
  // cases.
  int ok_boron_special_case;

  int compute_molecular_weight_for_each = 0;

  int compute_molecular_weight_based_on_largest_fragment = 0;

  int appends_to_be_done = 0;

  int do_appends_as_prepends = 0;

  int appended_properties_from_largest_fragment = 0;

  int append_molecular_weight_to_name = 0;

  int append_molecular_weight_ok_isotope_to_name = 0;

  int append_exact_mass_to_name = 0;

  int append_heteratom_count_to_name = 0;

  int append_molecular_formula_to_name = 0;

  int append_isis_molecular_formula_to_name = 0;

  int append_nrings_to_name = 0;

  int append_aromatic_ring_count_to_name = 0;

  int append_natoms_to_name = 0;

  int append_nbonds_to_name = 0;

  int append_net_formal_charge_to_name = 0;

  int append_clnd_count_to_name = 0;

  molecular_weight_t lower_amw_cutoff = -1.0;
  int molecules_below_amw_cutoff = 0;
  molecular_weight_t upper_amw_cutoff = -1.0;
  int molecules_above_amw_cutoff = 0;

  IWString substitute_for_whitespace_in_name;
  int truncate_names_to_first_token = 0;

  char truncate_name_at_first = '\0';
  char truncate_name_at_last = '\0';

  int name_token = -1;

  IWString prepend_to_name;

  resizable_array_p<Substructure_Query> aromatise_these_rings;

  int molecules_changed_by_aromatising_rings = 0;

  int remove_unnecessary_square_brackets = 0;

  int remove_all_possible_square_brackets = 0;

  // The number of molecules found with implicit valences satisfied
  // with implicit Hydrogen atoms
  int molecules_with_alternate_valence_via_implicit_h = 0;

  int reset_atom_map_numbers = 0;

  int molecule_to_fragments = 0;

  // If set, the smiles of the molecule before any filters are applied
  // will be stored in the result.
  int return_smiles_before_filters = 0;

  // May 2022. Any atom(s) can be excised.
  resizable_array_p<Substructure_Query> remove_atom;
  int molecules_with_removed_atoms = 0;

  // Private functions.
  private:
  void DefaultValues();

  // Called by Build, to construct the object.
  int ParseOrganicSpecification(Command_Line& cl, char flag);
  int ParseMolecularWeightSpecifications(Command_Line& cl);
  int ParseRingFiltering(Command_Line& cl);
  int ParseAtomCountDirectives(Command_Line& cl);
  int ParseBadHandlingoptions(Command_Line& cl, char flag);
  int ParseFragmentSizeWindow(const const_IWSubstring& f);
  int GetTranslations(Command_Line& cl, const char cflag);
  int GetChargeCalculation(Command_Line& cl, char flag);
  int ParseImplicitHydrogenDirectives(Command_Line& cl, char flag);
  int GetIsotopeDirectives(Command_Line& cl, char flag);
  int ParseMiscOptions(Command_Line& cl, char flag);
  int GetChiralityInstructions(Command_Line& cl, char flag);
  int GetFragmentSpecifications(Command_Line& cl);
  int GatherAppendSpecifications(Command_Line& cl, char flag);
  int ParseMkFragOptions(Command_Line& cl, char flag);
  int ParseFragmentAddToIsotope(const const_IWSubstring& s);
  int InitialiseRotation(const_IWSubstring buffer);

  // Functions that compute values, print things, change or filter the molecule.
  void DoAppends(Molecule& m, IWString& extra_stuff);
  int PrintTorsion(const Molecule& m, const Bond& b, std::ostream& output);
  int PrintBondAngle(const Molecule& m,
                     const atom_number_t a1,
                     const atom_number_t a2,
                     std::ostream& output);
  int IdentifyMatchedAtomsWithChiralCentres(const Molecule& m,
                                            const Set_of_Atoms& e,
                                            Set_of_Atoms& atoms_with_chiral_centres_to_be_removed);

  int AddFragmentToIsotopicAtoms(Molecule& m);
  int IsotopeRelated(Molecule& m);
  int RemoveChiralCentresOnMatchedAtoms(
      Molecule& m, const resizable_array_p<Substructure_Query>& remove_chiral_centres_on);
  int CountRingSystems(Molecule& m);
  int AromatiseTheseRings(Molecule& m, int* procecss_these_atoms);
  int AromatiseTheseRings(Molecule& m,
                          int* ring_labels,
                          int* procecss_these_atoms,
                          const Set_of_Atoms& e);
  int AromatiseTheseRings(Molecule& m,
                          resizable_array_p<Substructure_Query>& aromatise_these_rings);
  int ExtractSmallerFragmentsIntoNameTag(Molecule& m, const IWString& tag_for_removed_fragments);
  int RemoveByFragment(Molecule& m, const resizable_array<int>& fragments_to_remove);
  int RemoveFragmentsThatMatchQuery(Molecule& m, resizable_array_p<Substructure_Query>& queries);
  int RemoveLargestFragments(Molecule& m, int fragments_to_remove, int larger_or_smaller);
  int TooManyLargeFragments(Molecule& m, int min_size_max_fragment_count, int maximum_fragment_count);

  int RejectedForDisallowedAtoms(Molecule& m);

  void RemoveNonorganicFragments(Molecule& m, const int* frag_membership, int& fragments_removed);
  void RemoveNonorganicFragments(Molecule& m, int& fragments_removed);
  int IdentifyFragmentsToBeKept(Molecule& m, resizable_array_p<Substructure_Query>& q);
  int IdentifyFragmentByQuery(Molecule& m,
                              resizable_array_p<Substructure_Query>& queries,
                              int largest);
  int AtomsInNonLargestFragmentExceed(Molecule& m, int mxfs);
  int TrimToFirstNFragments(Molecule& m, const int* frag_membership, int fragment_count);
  int StripToNLargestFragments(Molecule& m, int n);
  template <typename C> int RemoveFragmentsByAtomCount(Molecule& m,
                        int threshold, C cmp);
  int LooksLikeMixtureByAtomCount(Molecule& m);

  void DoTranslation(Molecule& m) const;

  int ChangeName(Molecule& m, re2::RE2& change_name_rx);
  int NameMatchesNameRx(const IWString& mname);
  void SubstituteForWhitespaceInName(Molecule& m);
  int ZeroPadName(Molecule& m);
  int DebugPrint(Molecule& m, std::ostream& os) const;

  int ApplyAllFiltersInner(Molecule& m, IWString& rejection_reason, int& molecule_changed);
  int ApplyAllFilters(Molecule& m, int molecule_number, int& molecule_changed);

  int ComputeClnd(const Molecule& m);
  void DoAppends(Molecule& m);
  int RemoveFragmentsThisSizeOrSmaller(Molecule& m);
  int ComputePartialCharges(Molecule& m);
  int PrintMaxAtomSeparation(const Molecule& m, std::ostream& output);
  int PrintTorsions(Molecule& m, std::ostream& output);
  int PrintBondAngles(const Molecule& m, std::ostream& output);
  int PrintBondLengths(const Molecule& m, std::ostream& output);
  int InvertAllChiralCentres(Molecule& m);
  int ReflectCoordinates(Molecule& m);
  int FindAllChiralCentres(Molecule& m);
  int RemoveIsotopicAtoms(Molecule& m);
  int ConvertIsotopesToAtomMapNumbers(Molecule& m);
  int ConvertAtomMapNumbersToIsotopes(Molecule& m);
  int AtomTypeToIsotope(Molecule& m);
  int MakeImplicitHydrogensExplicit(Molecule& m);
  int RemoveNonOrganicChirality(Molecule& m);
  void SortByFragmentSize(Molecule& m);
  int ReduceToAllLargestFragments(Molecule& m);
  int RemoveLargestFragment(Molecule& m, std::function<bool(int, int)> less_or_greater);
  int ProcessFragments(Molecule& m);
  int ExcludeForNonOrganicAndNonPeriodicTable(const Molecule& m);
  int ExcludeForNonOrganic(const Molecule& m);
  int ExcludeForNonRealElements(const Molecule& m);
  int ExcludeForAtomTypes(const Molecule& m);
  int ValenceCheckOk(Molecule& m);
  int TooManyLargeFragments(Molecule& m);
  int RemoveAtoms(Molecule& m);

 public:
  // By default, nothing will be set. Useless in default state.
  FileconvConfig();

  // Use default command line options to initialise.
  int Build(Command_Line& cl);

  // Examines `m` and possibly transforms it in place.
  // Information about the processing is returned in the result.
  FileconvResult Process(Molecule& m);

  // Report results of potentially multiple invocations, drawn from
  // cl.size() input files.
  int ReportResults(const Command_Line& cl, std::ostream& output) const;
};

// Display a usage message to std::cerr and exit(rc)
void Usage(int rc);

} // namespace fileconv

#endif  // MOLECULE_TOOLS_FILECONV_OPTS_H
