#include <algorithm>
#include <iostream>
#include <memory>

// Performs reactions based on a reaction file.

#include "re2/re2.h"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/numass.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/rxn_file.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;

const char* prog_name = nullptr;

static int verbose = 0;

static int molecules_processed = 0;

static int products_written = 0;

static Number_Assigner number_assigner;

static Chemical_Standardisation chemical_standardisation;

static int suppress_invalid_valences = 0;

static Molecule_Output_Object stream_for_invalid_valence;

static Molecule_Output_Object stream_for_scaffolds_not_reacting;

static int strip_products_to_largest_fragment = 0;

static int remove_hydrogens_from_product_molecules = 0;

// Remove any Hydrogen that, if removed, would change an
// atom to an OK valence state.
static int remove_hydrogens_causing_valence_errors = 0;

static int avoid_overlapping_scaffold_changes = 0;

static int max_atoms_in_product = 0;

static int products_discarded_for_too_many_atoms = 0;

// Oct 2023. It is common to generate molecules that might
// have fragments that are too large or too small. These can
// be complex to filter out with other tools, so enable someone
// to specify the min and max fragment size in the product.
// Any product violating this is discarded.

static int need_to_check_product_fragment_sizes = 0;
static int min_allowed_fragment_size_in_product = 0;
static int max_allowed_fragment_size_in_product = std::numeric_limits<int>::max();
static int products_discarded_for_violating_fragment_specifications = 0;

// Nov 2023.
// A common operation is to fragment a molecule, and the next step will
// involve splitting the molecule into components.
// Note that each fragment must satisfy min_allowed_fragment_size_in_product and
// max_allowed_fragment_size_in_product.
static int write_multi_fragment_products_as_separate_molecules = 0;

// If set, the 'hits in scaffold' message is issued.
static int display_multiple_scaffold_hits_message = 1;

/*
  There are subtle differences between these two directives.
  make_implicit_hydrogens_explicit is more about the reaction and the scaffold
  whereas
  make_implicit_hydrogens_explicit_on_all_reagents applies just to the
  reagents.
*/

static int make_implicit_hydrogens_explicit = 0;

static int make_implicit_hydrogens_explicit_on_all_reagents = 0;

/*
  We can only react molecules where the name matches a regular expression
*/

static std::unique_ptr<RE2> only_react;

#define SMILES_TAG "$SMI<"

/*
  In order to make processing as simple as possible, the programme
  performs any number of reactions, but in the following manner.

    The -r option specifies the (single) main reaction to perform.
    The -y option specifies any number of secondarY reactions to perform.
           All secondary reactions must be self contained
*/

int number_secondary_reactions = 0;
static IWReaction* secondary_reactions = nullptr;

/*
  The reaction object holds details of element transformations which are
  applied to the matched atoms in the scaffold and sidechain(s). In addition
  we can apply regular element transformations to the resulting molecules
*/

static Element_Transformations etrans;

static int write_molecules_not_reacting = 0;

static extending_resizable_array<int> hit_statistics;

static int function_as_filter = 0;

static IWString identifier_tag("PCN<");

/*
  Sometimes when dealing with "reactions" which just remove things, we still
  want a means of knowing whether or not a molecule was modified. We can
  append a string
*/

static IWString append_to_changed_molecules;

/*
  Found it useful to be able append something only for those cases with multiple scaffold
  matches
*/

static int append_text_only_for_multiple_hits_in_scaffold = 0;

/*
  When we have multiple matches in the scaffold, we produce names like
  S.0, S.1
  We set up the separator, and if it is empty, we don't do that processing
*/

static IWString multiple_scaffold_match_name_separator('.');

static int apply_isotopes_for_atom_numbers = 0;

static int convert_isotopes = 0;

static int keep_atom_numbers_from_products = 0;

/*
  During debugging, it can be helpful to write scaffolds to the output stream
*/

static int write_scaffolds_to_output_stream = 0;
static int write_sidechains_to_output_stream = 0;

static Elements_to_Remove elements_to_remove;

/*
  When suppressing duplicates, we adopt a two-phase approach. First
  compare the smiles, then the unique smiles
*/

static int suppress_duplicate_molecules = 0;

static int duplicate_molecules_suppressed = 0;

static IW_STL_Hash_Set smiles_generated_current_molecule;
static IW_STL_Hash_Set unique_smiles_generated_current_molecule;

static int all_scaffold_possibilities_enumeration = 0;

static int unset_unnecessary_implicit_hydrogens_known_values = 0;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off

  cerr << "  -r <file>      specify single reaction file\n";
//cerr << "  -y <file>      specify secondary reaction(s)\n";
  cerr << "  -D <file>      specify ISIS reaction file\n";
  cerr << "  -D rmfrag      remove small fragments from the results of ISIS reaction file reactions\n";
  cerr << "  -D rmunmp      remove unmapped elements shown on LHS but absent on RHS\n";
  cerr << "  -P <file>      reaction as proto file\n";
  cerr << "  -K <smirks>    reaction as smirks. F:<fname> means smirks is in file\n";
  cerr << "  -z i           ignore molecules not reacting\n";
  cerr << "  -z w           write molecules not reacting\n";
  cerr << "  -Z             ignore sidechains not reacting\n";
  cerr << "  -C <string>    append <string> to the name of all changed molecules\n";
  cerr << "  -C ifmult      only append the -C string in the case of multiple scaffold matches\n";
  cerr << "  -m <number>    the maximum number of scaffold reaction sites to process\n";
  cerr << "  -m do=number   process site number <number> in the scaffold\n";
  cerr << "  -m each        enumerate each scaffold hit separately\n";
  cerr << "  -m RMX         ignore scaffolds that generate multiple substructure hits\n";
  cerr << "  -X <symbol>    extract/remove all atoms of type <symbol>. No bonds changed\n";
  cerr << "  -I             change isotopes to natural form in product molecules\n";
  cerr << "  -M all         generate all regio-isomers from multiple sidechain matches\n";
  cerr << "  -M do=number   process site number <number> in the sidechains\n";
  cerr << "  -M mskip=text  append <text> to names where just one possible sidechain attachment chosen\n";
  cerr << "  -M write=file  write non-reacting sidechains to <file>\n";
  cerr << "  -M RMX         ignore any sidechains with multiple substructure matches\n";
  cerr << "  -V <fname>     product molecules with invalid valences ignored and written to <fname> (NONE means no output)\n";
  cerr << "  -l             strip reagents to largest fragment\n";
  cerr << "  -L             strip products to largest fragment\n";
  cerr << "  -f             function as a TDT filter\n";
  cerr << "  -n <xxx>       number assigner options, enter '-n help' for details\n";
  cerr << "  -W <string>    token put between names of products (default \" + \")\n";
  cerr << "  -u             one embedding per start atom\n";
  cerr << "  -d             suppress duplicate molecules - only checks current molecule\n";
  cerr << "  -k             don't perceive symmetry equivalents in the scaffold\n";
  cerr << "  -J ...         various special purpose options, enter '-J help' for details\n";
  cerr << "  -S <string>    create output files with name stem <string>\n";
  cerr << "  -o <type>      specify output file type(s)\n";
  cerr << "  -E ...         standard element options, enter '-E help' for info\n";
  (void) display_standard_aromaticity_options(cerr);
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static void
reset_duplicate_hash_sets() {
  smiles_generated_current_molecule.clear();
  unique_smiles_generated_current_molecule.clear();

  return;
}

static int
molecule_is_duplicate(Molecule& m) {
  if (smiles_generated_current_molecule.contains(m.smiles())) {
    duplicate_molecules_suppressed++;
    return 1;
  }

  smiles_generated_current_molecule.insert(m.smiles());

  if (unique_smiles_generated_current_molecule.contains(m.unique_smiles())) {
    duplicate_molecules_suppressed++;
    return 1;
  }

  unique_smiles_generated_current_molecule.insert(m.unique_smiles());

  return 0;
}

static int
read_reaction_smiles(iwstring_data_source& input, RXN_File& rxn_smiles) {
  const_IWSubstring buffer;

  if (!input.next_record(buffer)) {
    cerr << "Empty reaction smiles file\n";
    return 0;
  }

  if (!rxn_smiles.build_from_reaction_smiles(buffer, 1)) {
    cerr << "Cannot interpret reaction smiles " << buffer << '\n';
    return 0;
  }

  return 1;
}

static int
read_reaction_smiles(const char* fname, RXN_File& rxn_smiles) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open reaction smiles file '" << fname << "'\n";
    return 0;
  }

  return read_reaction_smiles(input, rxn_smiles);
}

/*
  The reason that there can be only one secondary reaction is because
  of the need to produce the result in RESULT. If there were two
  secondary reactions, we would need to do some copying.
*/

static int
perform_secondary_reactions(Molecule* m, Molecule& result) {
  assert(1 == number_secondary_reactions);

  Molecule_to_Match target(m);

  int rc = 0;
  for (int i = 0; i < number_secondary_reactions; i++) {
    cerr << "Performing secondary reaction " << i << '\n';
    Substructure_Results sresults;

    int nhits = secondary_reactions[i].substructure_search(target, sresults);

    if (0 == nhits) {
      cerr << "No hits to secondary reaction " << i << " ignored\n";
      continue;
    }

    if (!secondary_reactions[i].perform_reaction(m, sresults, result)) {
      cerr << "Yipes, cannot perform secondary reaction " << i << '\n';
      return 0;
    }

    rc++;
  }

  if (rc) {
    result.set_name(m->name());
  }

  return 1;
}

static void
do_append_text(Molecule& m, const IWString& to_append) {
  IWString tmp = m.name();

  if (!tmp.ends_with(' ')) {
    tmp += ' ';
  }

  tmp += to_append;

  m.set_name(tmp);

  return;
}

// Collect the atom numbers of all explicit hydrogens
// attached to `atom` and put in `hydrogens`.
// We need to ensure that the atoms added to `hydrogens`
// are sorted, since they are removed by sequentially
// pop'ing the array.
// Returns true of we find any explicit Hydrogens.
static int
GatherExplicitHydrogens(const Molecule& m, atom_number_t zatom, Set_of_Atoms& hydrogens) {
  hydrogens.resize_keep_storage(0);
  const Atom& a = m.atom(zatom);
  for (const Bond* b : a) {
    atom_number_t h = b->other(zatom);
    if (m.atomic_number(h) == 1) {
      hydrogens << h;
    }
  }

  if (hydrogens.empty()) {
    return 0;
  }

  // No sorting needed.
  if (hydrogens.size() == 1) {
    return 1;
  }

  // Sorting is easy.
  if (hydrogens.size() == 2) {
    if (hydrogens[0] > hydrogens[1]) {
      hydrogens.swap_elements(0, 1);
    }
  }

  hydrogens.iwqsort_lambda([](int v1, int v2) {
    if (v1 < v2) {
      return -1;
    } else {
      return 1;
    }
  });

  return 1;
}

// Examine `m` for atoms with bad valence. If such an
// atom has explicit hydrogens attached, successively
// remove them to see if we end up with an OK valence.
static int
RemoveHydrogensCausingValenceErrors(Molecule& m) {
  const int matoms = m.natoms();

  Set_of_Atoms hydrogens;
  for (int i = 0; i < matoms; ++i) {
    if (m.valence_ok(i)) {
      continue;
    }

    if (!GatherExplicitHydrogens(m, i, hydrogens)) {
      continue;
    }

    // This seems reasonable.
    m.unset_all_implicit_hydrogen_information(i);

    do {
      m.remove_atom(hydrogens.pop());
    } while (!m.valence_ok(i) && !hydrogens.empty());
  }

  // Yes there could be a more informative return value,
  // but not needed now.
  return 1;
}

static int
OkFragmentsInProduct(Molecule& m) {
  const int nfrag = m.number_fragments();
  // It is unclear what should happen in this case.
  if (nfrag == 1) {
    return 1;
  }

  for (int i = 0; i < nfrag; ++i) {
    const int nat = m.atoms_in_fragment(i);
    if (nat < min_allowed_fragment_size_in_product) {
      return 0;
    }
    if (nat > max_allowed_fragment_size_in_product) {
      return 0;
    }
  }

  return 1;
}

// Remove from `fragments` all molecules that fail
// min_allowed_fragment_size_in_product or
// max_allowed_fragment_size_in_product.
static int
RemoveSizeNotOk(resizable_array_p<Molecule>& fragments) {
  return fragments.remove_items_fn([](const Molecule* m) {
    return m->natoms() > max_allowed_fragment_size_in_product ||
           m->natoms() < min_allowed_fragment_size_in_product;
  });
}

static int
do_write(Molecule_and_Embedding* sidechain, Molecule& product, int nhits,
         Molecule_Output_Object& output) {
  if (strip_products_to_largest_fragment) {
    product.reduce_to_largest_fragment();
  }

  if (make_implicit_hydrogens_explicit) {
    product.remove_all_atoms_with_isotope(make_implicit_hydrogens_explicit);
  }

  if (remove_hydrogens_causing_valence_errors) {
    RemoveHydrogensCausingValenceErrors(product);
  }

  if (remove_hydrogens_from_product_molecules) {
    product.remove_all(1);
  }

  if (max_atoms_in_product == 0) {
    // no checking
  }
  else if (write_multi_fragment_products_as_separate_molecules) {
    // will be checked later on a per fragment basis.
  } else if (product.natoms() > max_atoms_in_product) {
    products_discarded_for_too_many_atoms++;

    if (verbose > 1) {
      cerr << "Product '" << product.name() << "' too many atoms " << product.natoms()
           << '\n';
    }

    return 1;
  }

  if (suppress_duplicate_molecules && molecule_is_duplicate(product)) {
    return 1;
  }

  if (write_sidechains_to_output_stream && nullptr != sidechain) {
    output.write(*sidechain);
  }

  if (elements_to_remove.active()) {
    elements_to_remove.process(product);
  }

  if (etrans.number_elements()) {
    etrans.process(product);
  }

  if (convert_isotopes) {
    product.transform_to_non_isotopic_form();
  }

  if (!keep_atom_numbers_from_products) {
    product.reset_all_atom_map_numbers();
    //  product.unset_unnecessary_implicit_hydrogens_known_values();    not sure if this
    //  is needed or not
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(product);
  }

  if (suppress_invalid_valences && !product.valence_ok()) {
    cerr << "Molecule '" << product.name() << "' has an invalid valence\n";

    if (stream_for_invalid_valence.active()) {
      stream_for_invalid_valence.write(product);
    }

    return 1;
  }

  if (unset_unnecessary_implicit_hydrogens_known_values) {
    product.unset_unnecessary_implicit_hydrogens_known_values();
  }

  if (need_to_check_product_fragment_sizes) {
    if (! OkFragmentsInProduct(product)) {
      ++products_discarded_for_violating_fragment_specifications;
      return 1;
    }
  }


  if (number_assigner.active()) {
    number_assigner.process(product);
  }

  if (nhits && append_to_changed_molecules.nchars() &&
      0 == append_text_only_for_multiple_hits_in_scaffold) {
    do_append_text(product, append_to_changed_molecules);
  }

  if (verbose > 1) {
    cerr << "Writing " << product.name() << '\n';
  }

  ++products_written;

  if (write_multi_fragment_products_as_separate_molecules &&
      product.number_fragments() > 1) {
    resizable_array_p<Molecule> fragments;
    product.create_components(fragments);
    RemoveSizeNotOk(fragments);
    for (Molecule* frag : fragments) {
      frag->set_name(product.name());
      output.write(*frag);
    }
    return 1;
  }

  return output.write(product);
}

class Set_of_Atoms_Comparator {
 private:
  const int _matoms;

 public:
  Set_of_Atoms_Comparator(const int s);

  int operator()(const Set_of_Atoms* s1, const Set_of_Atoms* s2) const;
};

Set_of_Atoms_Comparator::Set_of_Atoms_Comparator(const int s) : _matoms(s) {
}

int
Set_of_Atoms_Comparator::operator()(const Set_of_Atoms* e1,
                                    const Set_of_Atoms* e2) const {
  const int n = e1->number_elements();
  assert(e2->number_elements() == n);

  uint64_t s1 = 0;
  uint64_t s2 = 0;

  for (int i = 0; i < n; ++i) {
    s1 = _matoms * s1 + e1->item(i);
    s2 = _matoms * s2 + e2->item(i);
  }

  return s1 < s2;
}

static int
do_all_scaffold_possibilities_enumeration(Molecule& scaffold, IWReaction& reaction,
                                          Substructure_Results& sresults,
                                          Molecule_Output_Object& output) {
  const int nhits = sresults.number_embeddings();
  assert(nhits > 1);

  if (write_scaffolds_to_output_stream) {
    output.write(scaffold);
  }

  const Set_of_Atoms** e = new const Set_of_Atoms*[nhits];
  std::unique_ptr<const Set_of_Atoms*[]> free_e(e);

  Set_of_Atoms_Comparator soac(scaffold.natoms());

  for (int i = 0; i < nhits; ++i) {
    e[i] = sresults.embedding(i);
  }

  std::sort(e, e + nhits, soac);

  do {
    sresults.set_embeddings(e, nhits);
#ifdef ECHO_NEW_PERMUTATION
    cerr << "New permutation";
    for (int i = 0; i < nhits; ++i) {
      cerr << ' ' << *sresults.embedding(i);
    }
    cerr << '\n';
#endif

    Reaction_Iterator iterator;
    for (iterator.initialise(reaction); iterator.active(); iterator++) {
      Molecule result;
      if (!reaction.perform_reaction(&scaffold, sresults, iterator, result)) {
        cerr << "Yipes, cannot react '" << scaffold.name() << " with " << iterator
             << '\n';
        return 0;
      }

      do_write(nullptr, result, nhits, output);
    }
  } while (std::next_permutation(e, e + nhits, soac));

  return 1;
}

/*
  If we are setting the name, set it. Otherwise do nothing
*/

static void
set_molecule_name(Molecule& m, const IWString& name_stem, int& ndx) {
  if (name_stem.empty()) {
    return;
  }

  IWString tmp(name_stem);
  tmp << ndx;
  ndx++;

  m.set_name(tmp);

  return;
}

static int
process_no_reagents_enumerate_scaffold_hits_combinatorial(
    const int depth, Molecule& m, IWString& name_stem, int& ndx, IWReaction& reaction,
    const Substructure_Results& sresults, const int istart,
    Molecule_Output_Object& output) {
  const int nhits = sresults.number_embeddings();

  for (int i = istart; i < nhits; ++i) {
    set_molecule_name(m, name_stem, ndx);

    Molecule result;

    if (!reaction.perform_reaction(&m, sresults.embedding(i), result)) {
      return 0;
    }

    do_write(nullptr, result, 1, output);

    if (i != nhits - 1 && depth > 0) {
      process_no_reagents_enumerate_scaffold_hits_combinatorial(
          depth - 1, result, name_stem, ndx, reaction, sresults, i + 1, output);
    }
  }

  return 1;
}

static int
process_no_reagents_enumerate_scaffold_hits_combinatorial(
    const int depth, Molecule& m, IWReaction& reaction,
    const Substructure_Results& sresults, Molecule_Output_Object& output) {
  if (write_scaffolds_to_output_stream) {
    output.write(m);
  }

  IWString mname;
  if (multiple_scaffold_match_name_separator.length() > 0) {
    mname << m.name() << multiple_scaffold_match_name_separator;
  }

  int ndx = 0;

  const int nhits = sresults.number_embeddings();

  for (int i = 0; i < nhits; ++i) {
    set_molecule_name(m, mname, ndx);

    Molecule result;
    if (!reaction.perform_reaction(&m, sresults.embedding(i), result)) {
      return 0;
    }

    do_write(nullptr, result, 1, output);

    if (i != nhits - 1 && depth > 0) {
      process_no_reagents_enumerate_scaffold_hits_combinatorial(
          depth - 1, result, mname, ndx, reaction, sresults, i + 1, output);
    }
  }

  return 1;
}

/*
  There are multiple hits to the scaffold query, we must enumerate them
*/

static int
process_no_reagents_enumerate_scaffold_hits_individually(
    Molecule& m, IWReaction& reaction, const Substructure_Results& sresults,
    Molecule_Output_Object& output) {
  if (write_scaffolds_to_output_stream) {  // or should this go before each product
    output.write(m);
  }

  IWString mname(m.name());
  if (multiple_scaffold_match_name_separator.length() > 0) {
    mname << multiple_scaffold_match_name_separator;
  }

  const int nhits = sresults.number_embeddings();
  for (int i = 0; i < nhits; i++) {
    if (multiple_scaffold_match_name_separator.length() > 0) {
      IWString tmp(mname);
      tmp << i;
      m.set_name(tmp);
    }

    Molecule result;
    if (!reaction.perform_reaction(&m, sresults.embedding(i), result)) {
      return 0;
    }

    do_write(nullptr, result, 1, output);
  }

  return output.good();
}

/*
  There are no sidechains to be added.
*/

static int
process_no_reagents(Molecule& m, IWReaction& reaction,
                    const Substructure_Results& sresults,
                    Molecule_Output_Object& output) {
  int nhits = sresults.number_embeddings();

  Molecule result;

  int rc;

  const Scaffold_Match_Conditions& smc = reaction.scaffold_match_conditions();
  ;

  if (smc.process_hit_number() >= 0) {
    uint32_t mdo = smc.process_hit_number();

    if (mdo >= sresults.number_embeddings()) {
      cerr << "Request to process embedding " << mdo << " but query produced "
           << sresults.number_embeddings() << " embeddings\n";
      return 0;
    }

    const Set_of_Atoms* e = sresults.embedding(mdo);
    rc = reaction.perform_reaction(&m, e, result);

    nhits = 1;
  } else if (1 == nhits) {
    rc = reaction.perform_reaction(&m, sresults.embedding(0), result);
  } else if (smc.enumerate_scaffold_hits_individually()) {
    return process_no_reagents_enumerate_scaffold_hits_individually(m, reaction, sresults,
                                                                    output);
  } else if (smc.combinatorial_expansion_of_scaffold_hits() >= 0) {
    return process_no_reagents_enumerate_scaffold_hits_combinatorial(
        smc.combinatorial_expansion_of_scaffold_hits(), m, reaction, sresults, output);
  } else if (avoid_overlapping_scaffold_changes) {
    rc = reaction.perform_reaction_recheck_scaffold_reactivity(&m, sresults, result);
  } else {
    rc = reaction.perform_reaction(&m, sresults, result);
  }

  if (0 == rc) {
    return 0;
  }

  if (write_scaffolds_to_output_stream) {
    output.write(m);
  }

  return do_write(nullptr, result, nhits, output);
}

static int
handle_scaffolds_not_reacting(Molecule& m,
                              const Scaffold_Match_Conditions& scaffold_match_conditions,
                              Molecule_Output_Object& output) {
  if (write_molecules_not_reacting) {
    return output.write(m);
  }

  if (stream_for_scaffolds_not_reacting.active()) {
    stream_for_scaffolds_not_reacting.write(m);
  }

  if (scaffold_match_conditions.ignore_not_reacting()) {
    return 1;
  }

  cerr << "No hits to scaffold '" << m.name() << "'\n";

  return 0;
}

static int
do_enumerate_scaffold_hits_combinatorial(const int depth, Molecule& m,
                                         const IWString& mname, int& ndx,
                                         IWReaction& reaction,
                                         const Substructure_Results& sresults,
                                         const int istart,
                                         Molecule_Output_Object& output) {
  const int nhits = sresults.number_embeddings();

  for (int i = istart; i < nhits; ++i) {
    set_molecule_name(m, mname, ndx);

    Reaction_Iterator iterator;
    for (iterator.initialise(reaction); iterator.active(); iterator++) {
      Molecule result;

      if (reaction.perform_reaction(&m, sresults.embedding(i), iterator, result)) {
      } else if (reaction.has_sidechain_isotope_requirement()) {
        // may have failed due to a mismatch on isotope matches.
        // Or it may have actually failed, cannot tell. Not great...
        continue;
      } else {
        cerr << "Reaction with " << i << "'th scaffold embedding and reagent " << iterator
             << " failed\n";
        return 0;
      }

      do_write(reaction.reagent(iterator), result, nhits, output);

      if (i != nhits - 1 && depth > 0) {
        do_enumerate_scaffold_hits_combinatorial(depth - 1, result, mname, ndx, reaction,
                                                 sresults, i + 1, output);
      }
    }
  }

  return 1;
}

static int
trxn(Molecule& m, IWReaction& reaction, Molecule_Output_Object& output) {
  Substructure_Results sresults;

  const int nhits = reaction.determine_matched_atoms(m, sresults);

  hit_statistics[nhits]++;

  const Scaffold_Match_Conditions& smc = reaction.scaffold_match_conditions();

  if (nhits == 0) {
    if (0 == verbose && 0 == smc.ignore_not_reacting()) {
      cerr << molecules_processed << " '" << m.name() << "' only matched "
           << sresults.max_query_atoms_matched_in_search() << " query atoms\n";
    }

    return handle_scaffolds_not_reacting(m, smc, output);
  } else if (verbose > 1 || (nhits > 1 && display_multiple_scaffold_hits_message)) {
    cerr << molecules_processed << " '" << m.name() << "' " << nhits
         << " hits in scaffold\n";
  }

  if (nhits > smc.suppress_if_more_than_this_many_substructure_search_hits()) {
    if (0 == verbose || smc.ignore_not_reacting()) {
      cerr << molecules_processed << " '" << m.name() << "' "
           << " " << nhits << " hits in scaffold, more than threshold "
           << smc.suppress_if_more_than_this_many_substructure_search_hits() << '\n';
    }

    return handle_scaffolds_not_reacting(m, smc, output);
  }

// #define ECHO_SCAFFOLD_EMBEDDINGS
#ifdef ECHO_SCAFFOLD_EMBEDDINGS
  for (int i = 0; i < nhits; i++) {
    const Set_of_Atoms* e = sresults.embedding(i);

    cerr << " scaffold embedding " << (*e) << '\n';
  }
#endif

  // Append multiple matches text here

  if (append_text_only_for_multiple_hits_in_scaffold && nhits > 1) {
    do_append_text(m, append_to_changed_molecules);
  }

  int nr = reaction.number_reagents();

  if (0 == nr) {  // all the sidechains must have their own fragment to add
    return process_no_reagents(m, reaction, sresults, output);
  }

  if (nhits > 1 && all_scaffold_possibilities_enumeration) {
    return do_all_scaffold_possibilities_enumeration(m, reaction, sresults, output);
  }

  int mdo = smc.process_hit_number();

  if (mdo >= nhits) {  // requested scaffold hit not present
    cerr << nhits << " hits to scaffold query, cannot process hit number " << mdo << '\n';

    return handle_scaffolds_not_reacting(m, smc, output);
  }

  // The cases in which we are processing just a single scaffold hit

  if (1 == nhits || mdo >= 0) {
    int scaffold_hit_to_process;
    if (mdo >= 0) {
      scaffold_hit_to_process = mdo;
    } else {
      scaffold_hit_to_process = 0;
    }

    const Set_of_Atoms* e = sresults.embedding(scaffold_hit_to_process);

    if (write_scaffolds_to_output_stream) {
      output.write(m);
    }

    Reaction_Iterator iterator;
    for (iterator.initialise(reaction); iterator.active();
         iterator++) {  // loop over reagents in the reaction
      Molecule result;
      if (reaction.perform_reaction(&m, e, iterator, result)) {
      } else if (reaction.has_sidechain_isotope_requirement()) {
        continue;
      } else {
        cerr << "Yipes, cannot react '" << m.name() << " with " << iterator << '\n';
        return 0;
      }

      do_write(reaction.reagent(iterator), result, nhits, output);
    }

    return output.ok();
  }

  // Now the more difficult case of multiple hits in the scaffold.
  // One possibility is to enumerate all possible scaffold hits

  if (smc.enumerate_scaffold_hits_individually()) {
    if (write_scaffolds_to_output_stream) {
      output.write(m);
    }

    IWString mname(m.name());
    if (multiple_scaffold_match_name_separator.length() > 0) {
      mname << multiple_scaffold_match_name_separator;
    }

    for (int i = 0; i < nhits; i++) {  // loop over scaffold hits
      const Set_of_Atoms* e = sresults.embedding(i);

      if (multiple_scaffold_match_name_separator.length() > 0) {
        IWString tmp(mname);
        tmp << i;
        m.set_name(tmp);
      }

      if (0 == nr) {
        Molecule result;

        if (reaction.perform_reaction(&m, e, result)) {
        } else if (reaction.has_sidechain_isotope_requirement()) {
          continue;
        } else {
          cerr << "Yipes, failed to perform reaction\n";
          return 0;
        }

        do_write(nullptr, result, nhits, output);
      } else {
        Reaction_Iterator iterator;
        for (iterator.initialise(reaction); iterator.active(); iterator++) {
          Molecule result;
          if (reaction.perform_reaction(&m, e, iterator, result)) {
          } else if (reaction.has_sidechain_isotope_requirement()) {
            continue;
          } else {
            cerr << "Reaction with " << i << "'th scaffold embedding and reagent "
                 << iterator << " failed\n";
            return 0;
          }

          do_write(reaction.reagent(iterator), result, nhits, output);
        }
      }
    }

    return output.ok();
  }

  if (smc.combinatorial_expansion_of_scaffold_hits() >= 0) {
    if (write_scaffolds_to_output_stream) {
      output.write(m);
    }

    IWString mname;
    if (multiple_scaffold_match_name_separator.length() > 0) {
      mname << m.name() << multiple_scaffold_match_name_separator;
    }

    int ndx = 0;

    // loop over scaffold hits
    for (int i = 0; i < nhits; i++) {
      set_molecule_name(m, mname, ndx);

      if (0 == nr) {
        Molecule result;

        if (reaction.perform_reaction(&m, sresults.embedding(i), result)) {
        } else if (reaction.has_sidechain_isotope_requirement()) {
          continue;
        } else {
          cerr << "Yipes, failed to perform reaction\n";
          return 0;
        }

        do_write(nullptr, result, nhits, output);
        if (i != nhits - 1 && smc.combinatorial_expansion_of_scaffold_hits() > 0) {
          process_no_reagents_enumerate_scaffold_hits_combinatorial(
              smc.combinatorial_expansion_of_scaffold_hits(), result, mname, ndx,
              reaction, sresults, i + 1, output);
        }
      } else {
        Reaction_Iterator iterator;
        for (iterator.initialise(reaction); iterator.active(); iterator++) {
          Molecule result;
          if (reaction.perform_reaction(&m, sresults.embedding(i), iterator, result)) {
          } else if (reaction.has_sidechain_isotope_requirement()) {
            continue;
          } else {
            cerr << "Reaction with " << i << "'th scaffold embedding and reagent "
                 << iterator << " failed\n";
            return 0;
          }

          do_write(reaction.reagent(iterator), result, nhits, output);

          if (i != nhits - 1 && smc.combinatorial_expansion_of_scaffold_hits() > 0) {
            do_enumerate_scaffold_hits_combinatorial(
                smc.combinatorial_expansion_of_scaffold_hits(), result, mname, ndx,
                reaction, sresults, i + 1, output);
          }
        }
      }
    }
  }

  if (avoid_overlapping_scaffold_changes) {
    Reaction_Iterator iterator;
    for (iterator.initialise(reaction); iterator.active(); iterator++) {
      Molecule result;

      if (reaction.perform_reaction_recheck_scaffold_reactivity(&m, sresults, iterator,
                                                                result)) {
      } else if (reaction.has_sidechain_isotope_requirement()) {
        continue;
      } else {
        cerr << "Reaction filed " << iterator << '\n';
        return 0;
      }

      if (write_scaffolds_to_output_stream) {
        output.write(m);
      }

      do_write(reaction.reagent(iterator), result, nhits, output);
    }

    return output.ok();
  }

  // The default behaviour is to process all reaction sites on the scaffold.
  // The reaction object does that.

  Reaction_Iterator iterator;
  for (iterator.initialise(reaction); iterator.active(); iterator++) {
    Molecule result;
    if (reaction.perform_reaction(&m, sresults, iterator, result)) {
    } else if (reaction.has_sidechain_isotope_requirement()) {
      continue;
    } else {
      cerr << "Yipes, cannot react '" << m.name() << " with " << iterator << '\n';
      return 0;
    }

    if (write_scaffolds_to_output_stream) {
      output.write(m);
    }

    do_write(reaction.reagent(iterator), result, nhits, output);
  }

  return 1;
}

/*
  We have changed the molecule and need to change the $SMI dataitem in
  the TDT and write it
*/

static int
do_write(IW_TDT& tdt, Molecule& result, std::ostream& output) {
  tdt.set_dataitem_value(SMILES_TAG, result.smiles());

  if (append_to_changed_molecules.empty()) {
    output << tdt;

    return output.good();
  }

  if (tdt.index_of_dataitem(identifier_tag) < 0) {
    if (append_to_changed_molecules.empty()) {
      tdt.set_dataitem_value(identifier_tag, result.name());
    } else {
      IWString tmp(result.name());
      tmp << ' ' << append_to_changed_molecules;

      tdt.set_dataitem_value(identifier_tag, tmp);
    }
  } else {
    IWString tmp(result.name());
    if (append_to_changed_molecules.length()) {
      tmp << ' ' << append_to_changed_molecules;
    }

    tdt.add_dataitem(identifier_tag, tmp);
  }

  output << tdt;

  return output.good();
}

static int
handle_scaffolds_not_reacting_filter(
    IW_TDT& tdt, Molecule& m, const Scaffold_Match_Conditions& scaffold_match_conditions,
    std::ostream& output) {
  if (stream_for_scaffolds_not_reacting.active()) {
    stream_for_scaffolds_not_reacting.write(m);
  }

  if (write_molecules_not_reacting) {
    output << tdt;
    return output.good();
  }

  if (scaffold_match_conditions.ignore_not_reacting()) {
    return output.good();
  }

  cerr << "No hits to scaffold query\n";

  return 0;
}

static int
trxn_filter(IW_TDT& tdt, IWReaction& rxn, std::ostream& output) {
  const_IWSubstring smiles;
  if (!tdt.dataitem_value(SMILES_TAG, smiles)) {
    cerr << "Cannot extract smiles from TDT\n";
    cerr << tdt;
    return 0;
  }

  Molecule m;
  if (!m.build_from_smiles(smiles)) {
    cerr << "Cannot parse smiles '" << smiles << "'\n";
    return 0;
  }

  molecules_processed++;

  Substructure_Results sresults;
  int nhits = rxn.determine_matched_atoms(m, sresults);

  hit_statistics[nhits]++;

  if (verbose) {
    const_IWSubstring id;
    (void)tdt.dataitem_value(identifier_tag, id);
    cerr << "Molecule " << molecules_processed;
    if (id.length()) {
      cerr << " '" << id << "'";
    }
    cerr << " found " << nhits << " hits\n";
  }

  const Scaffold_Match_Conditions& smc = rxn.scaffold_match_conditions();

  if (0 == nhits) {
    return handle_scaffolds_not_reacting_filter(tdt, m, smc, output);
  }

  assert(
      0 ==
      rxn.number_reagents());  // all the sidechains must have their own fragment to add

  Molecule result;

  int rc;

  int mdo = smc.process_hit_number();

  if (mdo >= 0) {
    const Set_of_Atoms* e = sresults.embedding(mdo);
    rc = rxn.perform_reaction(&m, e, result);
  } else {
    rc = rxn.perform_reaction(&m, sresults, result);
  }

  if (0 == rc) {
    return 0;
  }

  return do_write(tdt, result, output);
}

static int
trxn_filter(iwstring_data_source& input, IWReaction& rxn, std::ostream& output) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    if (!trxn_filter(tdt, rxn, output)) {
      return 0;
    }
  }

  return 1;
}

/*
  the only kind of reaction we can do when working as a filter is to
  change the molecule. We cannot do any enumeration
*/

static int
trxn_filter(const char* fname, IWReaction& rxn, std::ostream& output) {
  int nr = rxn.number_reagents();

  if (nr > 1) {
    cerr << "Sorry, no ability to enumerate molecules when working as a filter\n";
    return 0;
  }

  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open TDT file '" << fname << "'\n";
    return 0;
  }

  return trxn_filter(input, rxn, output);
}

static int
do_apply_isotopes_for_atom_numbers(Molecule& m) {
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    m.set_isotope(i, i);
  }

  return 1;
}

static int
OkWithOnlyReact(std::unique_ptr<RE2>& only_react, const IWString& name) {
  if (!only_react) {  // Not specified, OK by definition.
    return 1;
  }

  return iwre2::RE2PartialMatch(name, *only_react);
}

static int
trxn(data_source_and_type<Molecule>& input, IWReaction& rxn,
     Molecule_Output_Object& output) {
  Make_Implicit_Hydrogens_Explicit mihe;

  if (make_implicit_hydrogens_explicit) {
    mihe.set_isotope(make_implicit_hydrogens_explicit);
  }

  Molecule* m;

  while (nullptr != (m = input.next_molecule())) {
    molecules_processed++;

    //  if (only_react.active() && ! only_react.matches(m->name()))
    if (!OkWithOnlyReact(only_react, m->name())) {
      handle_scaffolds_not_reacting(*m, rxn.scaffold_match_conditions(), output);
      delete m;
      continue;
    }

    if (suppress_duplicate_molecules) {
      reset_duplicate_hash_sets();
    }

    if (make_implicit_hydrogens_explicit) {
      mihe.reset();
      m->make_implicit_hydrogens_explicit(mihe);
    } else if (make_implicit_hydrogens_explicit_on_all_reagents) {
      m->make_implicit_hydrogens_explicit();
    }

    //  cerr << "Processing scaffold '" << m->smiles() << "'\n";

    if (apply_isotopes_for_atom_numbers) {
      do_apply_isotopes_for_atom_numbers(*m);
    }

    if (number_secondary_reactions) {  // do these first
      Molecule* newmolecule = new Molecule;
      int ps = perform_secondary_reactions(m, *newmolecule);

      delete m;  // finished with this one, the new one is (hopefully) formed

      if (0 == ps) {  // yipes, failed to do the secondary reactions
        delete newmolecule;
        return 0;
      }

      m = newmolecule;  // so the rest of the programme processes the changed molecule
    }

    int rc = trxn(*m, rxn, output);

    delete m;
    if (0 == rc) {
      return 0;
    }
  }

  return 1;
}

static int
trxn(const char* fname, IWReaction& rxn, FileType input_type,
     Molecule_Output_Object& output) {
  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return trxn(input, rxn, output);
}

static void
display_dash_j_qualifiers(std::ostream& os) {
  // clang-format off
  os << R"(
 -J wrscaf      write scaffolds  to the output stream
 -J wrsdch      write sidechains to the output stream
 -J blki        bonds lose Kekule identity
 -J appinfo     append text info from reagents to products
 -J onlyreact=RX only react molecules where the name matches RX
 -J maxat=nn    discard any product molecule having more than <nn> atoms
 -J rmncm       ignore multiple substructure matches involving non changing atoms
 -J rmovm       ignore multiple substructure matches involving     changing atoms
 -J exph        make implicit Hydrogen atoms explicit (changes reaction)
 -J exphR       make implicit Hydrogen atoms explicit on reagents (reaction not changed)
 -J rmph        remove explicit hydrogen atoms from product molecules
 -J rcksm       when multiple scaffold hits present, re-check matches for activity
 -J numok       keep non-unique embeddings - default is unique embeddings only
 -J isonum      isotopically label atoms with their initial atom number
 -J msm=<s>     text designating multiple scaffold matches, use NONE to skip
 -J marvin      the input reaction file (-D) has come from Marvin
 -J keepatmn    retain any atom map numbers in output molecules
 -J larf        in smirks, if an atom is lost, remove the fragment
 -J rmhsqb      remove unnecessary [] in product molecules
 -J rmxhbv      remove explicit hydrogens causing bad valences
 -J coords      include coordinates with smiles output
 -J minpfs=<n>  discard products with a fragment with < minpfs atoms
 -J maxpfs=<n>  discard products with a fragment with > maxpfs atoms
 -J mfpseparate write multi fragment products as separate molecules
 -J nomshmsg    do NOT write 'hits in scaffold' messages for multiple scaffold query hits
 -J noschmsg    do NOT write warning messages about no sidechain substructure matches
)";
  // clang-format on

  return;
}

static int
ReadReaction(IWString& fname, IWReaction& rxn) {
  std::optional<ReactionProto::Reaction> maybe_proto =
      iwmisc::ReadTextProtoCommentsOK<ReactionProto::Reaction>(fname);
  if (!maybe_proto) {
    cerr << "ReadReaction:cannot read reaction proto from '" << fname << "'\n";
    return 0;
  }

  if (!rxn.ConstructFromProto(*maybe_proto, fname)) {
    cerr << "ReadReaction:cannot parse reaction proto\n";
    cerr << maybe_proto->ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}

static int
trxn(int argc, char** argv) {
  Command_Line cl(argc, argv, "r:i:vS:A:E:n:o:z:ZX:t:m:C:M:V:g:y:fIlLW:ukJ:D:P:dK:R:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_elements(cl)) {
    usage(2);
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    usage(5);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose, 'g')) {
      cerr << "Cannot construct chemical standardisation object (-g option)\n";
      usage(7);
    }
  }

  if (cl.option_present('L')) {
    strip_products_to_largest_fragment = 1;
    if (verbose) {
      cerr << "Products will be stripped to their largest fragment\n";
    }
  }

  if (cl.option_present('l')) {
    set_strip_reagents_to_largest_fragment(1);
    if (verbose) {
      cerr << "Reagents will be stripped to their largest fragment\n";
    }
  }

  if (cl.option_present('n')) {
    if (!number_assigner.initialise(cl, 'n', verbose)) {
      cerr << "Cannot process -n option\n";
      usage(51);
    }
  }

  if (cl.option_present('W')) {
    const_IWSubstring w = cl.string_value('W');

    set_component_separator(w);

    if (verbose) {
      cerr << "Component separator set to '" << w << "'\n";
    }
  }

  if (cl.option_present('t')) {
    if (!etrans.construct_from_command_line(cl, verbose, 't')) {
      cerr << "Cannot parse -t options\n";
      usage(3);
    }
  }

  if (cl.option_present('C')) {
    int i = 0;
    const_IWSubstring c;
    while (cl.value('C', c, i++)) {
      if ("ifmult" == c) {
        append_text_only_for_multiple_hits_in_scaffold = 1;
      } else {
        append_to_changed_molecules = c;
      }
    }

    if (append_to_changed_molecules.empty()) {
      cerr << "No text to append specified (-C option)\n";
      usage(6);
    }

    if (0 == verbose) {
      ;
    } else if (append_text_only_for_multiple_hits_in_scaffold) {
      cerr << "Will append " << append_to_changed_molecules
           << " to molecules with multiple scaffold hits\n";
    } else {
      cerr << "Will append '" << append_to_changed_molecules
           << "' to all changed molecules\n";
    }

    //  It is convenient to insert a leading space

    append_to_changed_molecules.shift(1, ' ');
  }

  if (cl.option_present('I')) {
    convert_isotopes = 1;
    if (verbose) {
      cerr << "Isotopes in product molecules will be converted\n";
    }
  }

  if (cl.option_present('V')) {
    const_IWSubstring fname;
    cl.value('V', fname);

    suppress_invalid_valences = 1;

    if (fname == "NONE") {
      if (verbose) {
        cerr << "Products with invalid valences not written\n";
      }
    } else {
      stream_for_invalid_valence.add_output_type(FILE_TYPE_SMI);

      if (!stream_for_invalid_valence.new_stem(fname)) {
        cerr << "Cannot set invalid valence stream stem to '" << fname << "'\n";
        return 62;
      }

      if (verbose) {
        cerr << "Molecules with invalid valences written to '" << fname << ".smi'\n";
      }
    }
  }

  if (cl.option_present('f')) {
    if (cl.option_present('i')) {
      cerr
          << "The -f (function as filter) and -i (input type) options are incompatible\n";
      usage(14);
    }

    function_as_filter = 1;
    if (verbose) {
      cerr << "Will function as a TDT filter\n";
    }
  }

  if (cl.option_present('d')) {
    suppress_duplicate_molecules = 1;

    if (verbose) {
      cerr << "Will suppress duplicate molecules - within current generation\n";
    }
  }

  IWReaction rxn;
  rxn.set_append_names(1);
  Scaffold_Match_Conditions& scaffold_match_conditions = rxn.scaffold_match_conditions();
  Sidechain_Match_Conditions sidechain_match_conditions;

  sidechain_match_conditions.set_verbose(verbose);

  int ignore_multiple_matches_involving_atoms_not_changing = 0;
  int ignore_multiple_matches_involving_changing_atoms = 0;

  int processing_marvin_file = 0;

  if (cl.option_present('J')) {
    int i = 0;
    const_IWSubstring j;
    while (cl.value('J', j, i++)) {
      if ("help" == j) {
        display_dash_j_qualifiers(cerr);
        return 0;
      } else if ("appinfo" == j) {
        set_reaction_transfer_text_info_to_products(1);
        if (verbose) {
          cerr << "Will transfer text info to products\n";
        }
      } else if ("wrscaf" == j) {
        write_scaffolds_to_output_stream = 1;
        if (verbose) {
          cerr << "Will write scaffolds to output stream\n";
        }
      } else if ("wrsdch" == j) {
        write_sidechains_to_output_stream = 1;
        if (verbose) {
          cerr << "Will write sidechains to the output stream\n";
        }
      } else if ("blki" == j) {
        set_aromatic_bonds_lose_kekule_identity(1);
        if (verbose) {
          cerr << "Aromatic bonds lose Kekule identity\n";
        }
      } else if (j.starts_with("onlyreact=")) {
        IWString tmp = j;
        tmp.remove_up_to_first('=');
        if (!iwre2::RE2Reset(only_react, tmp)) {
          cerr << "Cannot initialise only react regular expression from '" << tmp
               << "'\n";
          return 9;
        }

        if (verbose) {
          cerr << "Will only react scaffolds whose names match '" << only_react->pattern()
               << "'\n";
        }
      } else if (j.starts_with("maxat=")) {
        j.remove_leading_chars(6);

        if (!j.numeric_value(max_atoms_in_product) || max_atoms_in_product < 1) {
          cerr << "The maximum atoms in a product value must be a whole number > 0\n";
          usage(3);
        }

        if (verbose) {
          cerr << "Will discard products having more than " << max_atoms_in_product
               << " atoms\n";
        }
      } else if ("rmncm" == j) {
        ignore_multiple_matches_involving_atoms_not_changing = 1;

        if (verbose) {
          cerr << "Will ignore multiple matches when non changing atoms are involved\n";
        }
      } else if ("rmovm" == j) {
        ignore_multiple_matches_involving_changing_atoms = 1;

        if (verbose) {
          cerr << "Will ignore multiple matches when changing atoms are involved\n";
        }
      } else if ("exph" == j) {
        // an unlikely number so we can pull them off at the end
        make_implicit_hydrogens_explicit = 96632;
        rxn.set_make_implicit_hydrogens_explicit(make_implicit_hydrogens_explicit);

        if (verbose) {
          cerr << "Implicit Hydrogens will be made explicit\n";
        }
      } else if ("exphR" == j) {
        make_implicit_hydrogens_explicit_on_all_reagents = 1;

        if (verbose) {
          cerr << "Will convert implicit hydrogens to explicit on all reagents\n";
        }
      } else if ("rmph" == j) {
        remove_hydrogens_from_product_molecules = 1;

        if (verbose) {
          cerr << "Will remove explicit hydrogens from all product molecules\n";
        }
      } else if (j == "rmxhbv") {
        remove_hydrogens_causing_valence_errors = 1;
        if (verbose) {
          cerr << "Will remove explicit hydrogens inducing valence errors\n";
        }
      } else if ("rcksm" == j) {
        avoid_overlapping_scaffold_changes = 1;

        if (verbose) {
          cerr << "Will re-check scaffold matches after each hit\n";
        }
      } else if ("numok" == j) {
        scaffold_match_conditions.set_find_unique_embeddings_only(0);
        sidechain_match_conditions.set_find_unique_embeddings_only(0);
      } else if ("isonum" == j) {
        apply_isotopes_for_atom_numbers = 1;
      } else if (j.starts_with("msm=")) {
        j.remove_leading_chars(4);
        if ("NONE" == j) {
          multiple_scaffold_match_name_separator = "";
          if (verbose) {
            cerr << "No designation of multiple scaffold matches\n";
          }
        } else {
          multiple_scaffold_match_name_separator = j;
          if (verbose) {
            cerr << "Multiple scaffold matches numbered with separator '"
                 << multiple_scaffold_match_name_separator << "'\n";
          }
        }
      } else if ("marvin" == j) {
        processing_marvin_file = 1;
      } else if ("keepatmn" == j) {
        set_include_atom_map_with_smiles(1);
        keep_atom_numbers_from_products = 1;
        if (verbose) {
          cerr << "Will keep atom number information from products\n";
        }
      } else if ("larf" == j) {
        set_smirks_lost_atom_means_remove_frgment(1);
        if (verbose) {
          cerr << "In smirks, lost atoms mean remove fragment\n";
        }
      } else if (j == "rmhsqb") {
        unset_unnecessary_implicit_hydrogens_known_values = 1;
        if (verbose) {
          cerr << "Square brackets removed from products if possible\n";
        }
      } else if (j == "coords") {
        set_append_coordinates_after_each_atom(1);
      } else if (j.starts_with("minpfs=")) {
        j.remove_leading_chars(7);
        if (! j.numeric_value(min_allowed_fragment_size_in_product) || 
            min_allowed_fragment_size_in_product < 1) {
          cerr << "The min fragment size in product molecule directive 'minpfs=' must be a while +ve number\n";
          return 1;
        }
        if (verbose) {
          cerr << "Will discard molecules where any fragment has less than " << min_allowed_fragment_size_in_product << " atoms\n";
        }
        need_to_check_product_fragment_sizes = 1;
      } else if (j.starts_with("maxpfs=")) {
        j.remove_leading_chars(7);
        if (! j.numeric_value(max_allowed_fragment_size_in_product) ||
          max_allowed_fragment_size_in_product < min_allowed_fragment_size_in_product) {
          cerr << "The max fragment size in product molecule directive 'maxpfs=' must be a whole +ve number\n";
          return 1;
        }
        if (verbose) {
          cerr << "Will discard molecules where any fragment has more than " << max_allowed_fragment_size_in_product << " atoms\n";
        }
        need_to_check_product_fragment_sizes = 1;
      } else if (j == "mfpseparate") {
        write_multi_fragment_products_as_separate_molecules = 1;
        if (verbose) {
          cerr << "Will write multi fragment products as separate molecules\n";
        }
      } else if (j == "nomshmsg") {
        display_multiple_scaffold_hits_message = 0;
        if (verbose) {
          cerr << "Will NOT write messagea about multiple scaffold query hits\n";
        }
      } else if (j == "noschmsg") {
        if (verbose) {
          cerr << "Will NOT warn about no sidechain query matches\n";
        }
        scaffold_match_conditions.set_issue_sidechain_no_match_warnings(0);
      } else {
        cerr << "Unrecognised -J qualifier '" << j << "'\n";
        display_dash_j_qualifiers(cerr);
        return 1;
      }
    }
  }

  if (convert_isotopes && apply_isotopes_for_atom_numbers) {
    cerr << "Doesn't make sense to remove isotopes and to isotopically label atoms by "
            "atom number\n";
    usage(5);
  }

  if (1 == cl.option_count('r')) {  // a single reaction
    ;
  } else if (cl.option_present('D')) {  // an ISIS reaction file - perhaps with qualifiers
    ;
  } else if (cl.option_present('K')) {  // smirks
    ;
  } else if (cl.option_present('R')) {  // reaction smiles
    ;
  } else if (cl.option_present('P')) {  // reaction as proto
    ;
  } else {
    cerr << "Must specify reaction file via a single '-r' or '-P' option\n";
    usage(7);
  }

  // Process the -M option before the query is read.
  // I ran into a problem with
  //     (A C smiles "CC")
  //     (A C smarts "C")
  // which failed because there were two matches to the smarts,
  // therefore we handle -M before -r.

  if (cl.option_present('M')) {
    const_IWSubstring m;

    cl.value('M', m);

    if ("all" == m) {
      sidechain_match_conditions.set_make_new_reagent_for_each_hit(1);
      if (verbose) {
        cerr << "All regio-isomers will be generated\n";
      }
    } else if (m.starts_with("do=")) {
      m.remove_leading_chars(3);

      int mm;

      if (!m.numeric_value(mm) || mm < 0) {
        cerr << "The '-M do=' option must be followed by a whole non-negative number\n";
        usage(99);
      }

      sidechain_match_conditions.set_process_hit_number(mm);
    } else if (m.starts_with("mskip=")) {
      const_IWSubstring ms = m;
      ms.remove_leading_chars(6);
      sidechain_match_conditions.set_multiple_match_string(ms);
      if (verbose) {
        cerr << "Where just one of several possible sidechain matches used, will append '"
             << ms << "'\n";
      }
    } else if (m.starts_with("write=")) {
      m.remove_leading_chars(6);
      if (!set_stream_for_sidechains_not_matching_query(m)) {
        cerr << "Cannot initialise stream for non reacting sidechains '" << m << "'\n";
        return 0;
      }

      if (verbose) {
        cerr << "Sidechains not reacting written to '" << m << ".smi\n";
      }

      sidechain_match_conditions.set_ignore_not_reacting(1);
    } else if ("RMX" == m) {
      sidechain_match_conditions.set_ignore_multiple_substucture_matches(1);
    } else {
      cerr << "Unrecognised -M qualifier '" << m << "'\n";
      usage(23);
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  // if (cl.option_present('h'))    // never want this

  rxn.set_query_files_in_current_directory(0);

  if (cl.option_present('r')) {
    const_IWSubstring fname;
    cl.value('r', fname);
    if (!rxn.do_read(fname, sidechain_match_conditions)) {
      cerr << "Cannot read main reaction from '" << fname << "'\n";
      return 43;
    }

    if (verbose > 2) {
      rxn.debug_print(cerr);
    }
  } else if (cl.option_present('D')) {
    set_auto_create_new_elements(1);

    RXN_File ISIS_rxn;

    if (processing_marvin_file) {
      ISIS_rxn.set_interpret_atom_alias_as_smarts(0);
    }

    if (!parse_isis_rxn_file_options(cl, 'D', ISIS_rxn)) {
      cerr << "Cannot parse -D option(s)\n";
      return 5;
    }

    RXN_File_Create_Reaction_Options rxnfcro;

    if (!ISIS_rxn.create_reaction(rxn, rxnfcro)) {
      cerr << "cannot parse ISIS reaction '" << ISIS_rxn.fname() << "'\n";
      return 8;
    }
  } else if (cl.option_present('K')) {
    set_auto_create_new_elements(1);

    const char* k = cl.option_value('K');

    if (!rxn.construct_from_smirks(k)) {
      cerr << "Cannot read smirks from '" << k << "'\n";
      return 2;
    }
  } else if (cl.option_present('R')) {
    const char* fname = cl.option_value('R');

    set_auto_create_new_elements(1);

    RXN_File rxn_smiles;

    if (!read_reaction_smiles(fname, rxn_smiles)) {
      cerr << "Cannot read reaction smiles file '" << fname << "'\n";
      return 1;
    }

    RXN_File_Create_Reaction_Options rxnfcro;

    if (!rxn_smiles.create_reaction(rxn, rxnfcro)) {
      cerr << "Cannot create reaction from reaction smiles " << fname << '\n';
      return 1;
    }
  } else if (cl.option_present('P')) {
    IWString fname = cl.option_value('P');

    set_auto_create_new_elements(1);

    if (!ReadReaction(fname, rxn)) {
      cerr << "Cannot read reaction proto file -P " << fname << '\n';
      return 1;
    }
  }

  if (verbose > 2) {
    rxn.write_msi(cerr);
  }

  if (scaffold_match_conditions.find_unique_embeddings_only()) {
    rxn.set_find_unique_embeddings_only(1);
  }

  if (cl.option_present('u')) {
    rxn.set_find_one_embedding_per_atom(1);
    rxn.set_find_unique_embeddings_only(1);
    sidechain_match_conditions.set_one_embedding_per_start_atom(1);
    if (verbose) {
      cerr << "One embedding per start atom\n";
    }
  }

  if (cl.option_present('k')) {
    rxn.set_do_not_perceive_symmetry_equivalent_matches(1);
    sidechain_match_conditions.set_ignore_symmetry_related_matches(1);

    if (verbose) {
      cerr << "Symmetric equivalents in the scaffold not perceived\n";
    }
  }

  if (ignore_multiple_matches_involving_atoms_not_changing ||
      ignore_multiple_matches_involving_changing_atoms) {
    rxn.setup_to_skip_multiple_embeddings_involving_non_changing_atoms();

    if (ignore_multiple_matches_involving_atoms_not_changing) {
      rxn.set_ignore_multiple_matches_involving_atoms_not_changing(1);
    }
    if (ignore_multiple_matches_involving_changing_atoms) {
      rxn.set_ignore_multiple_matches_involving_changing_atoms(1);
    }
  }

  if (verbose > 2) {
    rxn.write_msi(cerr);
  }

  if (cl.option_present('y')) {
    number_secondary_reactions = cl.option_count('y');

    assert(1 == number_secondary_reactions);  // limitation for now, fix this sometime

    secondary_reactions = new IWReaction[number_secondary_reactions];

    for (int i = 0; i < number_secondary_reactions; i++) {
      const_IWSubstring fname;
      cl.value('y', fname, i);

      if (!secondary_reactions[i].do_read(fname, sidechain_match_conditions)) {
        cerr << "Cannot create secondary reaction " << i << " from '" << fname << "'\n";
        return 8;
      }

      secondary_reactions[i].set_append_names(1);
    }

    if (verbose) {
      cerr << "Created " << number_secondary_reactions << " secondary reactions\n";
    }
  }

  // Apply -m after reading from file so we override any value in the file
  // Note that the variable mdo is global to this file

  if (cl.option_present('m')) {
    int i = 0;
    const_IWSubstring m;
    int mmx = -1;
    int mdo = -1;
    while (cl.value('m', m, i++)) {
      if (m.starts_with("do=")) {
        m += 3;
        if (!m.numeric_value(mdo) || mdo < 0) {
          cerr << "The '-m do=' option must be followed by a whole non-negative number\n";
          usage(29);
        }
      } else if ("each" == m) {
        scaffold_match_conditions.set_enumerate_scaffold_hits_individually(1);
        rxn.set_do_not_perceive_symmetry_equivalent_matches(1);
        if (verbose) {
          cerr << "Will enumerate multiple scaffold hits individually\n";
        }
      } else if (m.starts_with("comb=")) {
        m.remove_leading_chars(5);
        int c;
        if (!m.numeric_value(c) || c < 0) {
          cerr << "The number of scaffold hit combinations to combinatorially enumerate "
                  "(comb=nn) must be a whole number greater than 1\n";
          usage(1);
        }

        scaffold_match_conditions.set_combinatorial_expansion_of_scaffold_hits(c);
      } else if (m.starts_with("comb")) {
        scaffold_match_conditions.set_combinatorial_expansion_of_scaffold_hits(
            std::numeric_limits<int>::max());
        //      rxn.set_do_not_perceive_symmetry_equivalent_matches(1);
        if (verbose) {
          cerr << "Will combinatorially enumerate scaffold hits\n";
        }
      } else if ("RMX" == m) {
        scaffold_match_conditions.set_ignore_multiple_substucture_matches(1);
        if (verbose) {
          cerr << "Will discard scaffolds that show multiple substructure search hits\n";
        }
      }
#ifdef M_ALL_DIRECTIVE
      else if ("all" ==
               m)  // this is not necessary. When there are multiple scafold matches, they
                   // will be discovered during substructure searching
      {
        all_scaffold_possibilities_enumeration = 1;
        if (verbose) {
          cerr << "WIll enumerate all scaffold possibilities\n";
        }
      }
#endif
      else if (m.numeric_value(mmx)) {
        if (!m.numeric_value(mmx) || mmx < 1) {
          cerr
              << "The '-m <number>' option must be followed by a whole positive number\n";
          usage(91);
        }
      } else {
        cerr << "Unrecognised -m qualifier '" << m << "'\n";
        usage(37);
      }
    }

    //  If multiple forms are specified, they must be compatible

    if (scaffold_match_conditions.enumerate_scaffold_hits_individually() && mdo >= 0) {
      cerr << "Enumerating multiple scaffold hits (-m each) is incompatible with '-m "
              "do=nn'\n";

      usage(51);
    }

    if (mmx > 0 && mdo >= 0) {
      if (mdo >= mmx) {
        cerr << "You asked to process scaffold hit number " << mdo << '\n';
        cerr << "You asked to process a maximum of " << mmx << " scaffold hits\n";
        cerr << "Impossible\n";
        usage(28);
      }
    }

    if (mmx > 0) {
      rxn.set_max_matches_to_find(mmx);
      if (verbose) {
        cerr << "Max matches to find is " << m << '\n';
      }
    }

    if (mdo >= 0) {
      scaffold_match_conditions.set_process_hit_number(mdo);
      if (verbose) {
        cerr << "Will process hit number " << mdo << " of scaffold hits\n";
      }
    }
  }

  if (cl.option_present('Z')) {
    sidechain_match_conditions.set_ignore_not_reacting(1);

    if (verbose) {
      cerr << "Reagents not matching the query will be ignored\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (function_as_filter) {
    ;
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  // If a sidechain already has a molecule - it was in the msi file, then it doesn't need
  // a file on the command line

  int sidechains_that_are_just_single_reagents = rxn.number_sidechains_with_reagents();

  int files_on_command_line =
      cl.number_elements() - 1;  // the first one is for the scaffold

  if (files_on_command_line !=
      rxn.number_sidechains() - sidechains_that_are_just_single_reagents) {
    cerr << "Reaction contains " << rxn.number_sidechains() << " sidechains, "
         << sidechains_that_are_just_single_reagents << " have just a single reagent\n";
    cerr << "But " << files_on_command_line << " files specified on the command line\n";
    cerr << "Impossible\n";
    return 9;
  }

  int j = 1;  // index into Command_Line object
  for (int i = 0; i < rxn.number_sidechains(); i++) {
    Sidechain_Reaction_Site* s = rxn.sidechain(i);

    if (s->single_reagent_only()) {
      continue;
    }

    if (make_implicit_hydrogens_explicit_on_all_reagents) {
      s->set_make_implicit_hydrogens_explicit(1);
    }

    if (!s->add_reagents(cl[j], input_type, sidechain_match_conditions)) {
      cerr << "Bad news, sidechain reaction " << i << " could not read reagents from '"
           << cl[j] << "'\n";
      return i + 1;
    }

    j++;
  }

  if (verbose) {
    cerr << "Reaction '" << rxn.comment() << "' read " << rxn.number_reagents()
         << " reagents\n";
    cerr << "will create " << rxn.number_products_per_scaffold_embedding()
         << " products per scaffold embedding\n";

    if (verbose > 1) {
      rxn.debug_print(cerr);
    }
  }

  if (cl.option_present('z')) {
    const_IWSubstring z;
    int i = 0;
    while (cl.value('z', z, i++)) {
      if ('w' == z) {
        write_molecules_not_reacting = 1;
        scaffold_match_conditions.set_ignore_not_reacting(1);
        if (verbose) {
          cerr << "Molecules not reacting will be written to the output stream\n";
        }
      } else if ('i' == z) {
        scaffold_match_conditions.set_ignore_not_reacting(1);
        if (verbose) {
          cerr << "Molecules not reacting will be ignored\n";
        }
      } else if (z.starts_with("write=")) {
        z.remove_leading_chars(6);
        if (cl.option_present('o')) {
          if (!stream_for_scaffolds_not_reacting.determine_output_types(cl)) {
            cerr << "Cannot determine output type(s) for -z write= file\n";
            return 5;
          }
        } else {
          stream_for_scaffolds_not_reacting.add_output_type(FILE_TYPE_SMI);
        }

        if (stream_for_scaffolds_not_reacting.would_overwrite_input_files(cl, z)) {
          cerr << "Cannot overwrite input file(s), '" << z << "'\n";
          return 5;
        }

        if (!stream_for_scaffolds_not_reacting.new_stem(z)) {
          cerr << "Cannot open stream for scaffolds not reacting '" << z << "'\n";
          return 3;
        }

        if (verbose) {
          cerr << "Scaffolds not reacting written to '" << z << "'\n";
        }
      } else {
        cerr << "Unrecognised -z directive '" << z << "'\n";
        usage(19);
      }
    }
  }

  if (cl.option_present('X')) {
    if (!elements_to_remove.construct_from_command_line(cl, verbose, 'X')) {
      cerr << "Cannot discern elements to remove from -X switch\n";
      usage(18);
    }
  }

  Molecule_Output_Object output;

  if (function_as_filter) {  // output object not really used
    ;
  } else if (!cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SMI);
    if (verbose) {
      cerr << "Smiles output by default\n";
    }
  } else if (!output.determine_output_types(cl)) {
    cerr << "Cannot determine output type(s)\n";
    usage(6);
  }

  int rc;
  if (function_as_filter) {
    rc = trxn_filter(cl[0], rxn, std::cout);
  } else {
    IWString output_file_name;
    if (cl.option_present('S')) {
      cl.value('S', output_file_name);
    } else {
      output_file_name = "trxn";
    }

    for (int i = 0; i < cl.number_elements(); i++) {
      if (output.would_use_name(output_file_name, cl[i])) {
        cerr << "Output stem '" << output_file_name << " would overwrite input file '"
             << cl[i] << "'\n";
        return i + 1;
      }
    }

    if (!output.new_stem(output_file_name, 1)) {
      cerr << "Cannot direct output stream(s) to '" << output_file_name << "'\n";
      return 61;
    }

    if (verbose) {
      cerr << "Output written to file(s) with stem '" << output_file_name << "'\n";
    }

    rc = trxn(cl[0], rxn, input_type, output);
  }

  if (verbose) {
    cerr << "Processed " << molecules_processed << " molecules\n";
    for (int i = 0; i < hit_statistics.number_elements(); i++) {
      if (hit_statistics[i]) {
        cerr << hit_statistics[i] << " scaffolds had " << i << " hits\n";
      }
    }

    if (max_atoms_in_product > 0) {
      cerr << products_discarded_for_too_many_atoms
           << " products discarded for having more than " << max_atoms_in_product
           << " atoms\n";
    }

    if (chemical_standardisation.active()) {
      chemical_standardisation.report(cerr);
    }

    if (need_to_check_product_fragment_sizes) {
      cerr << products_discarded_for_violating_fragment_specifications <<
              " product molecules violated for fragment size constraints\n";
    }

    if (suppress_duplicate_molecules) {
      cerr << duplicate_molecules_suppressed << " duplicate products suppressed\n";
    }

    cerr << products_written << " products written\n";
  }

  if (nullptr != secondary_reactions) {
    delete[] secondary_reactions;
  }

  if (rc) {
    return 0;
  } else {
    return 1;
  }
}

#ifndef DLL_FLAG
int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = trxn(argc, argv);

  return rc;
}
#endif
