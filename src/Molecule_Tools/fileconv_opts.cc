/*
  File conversion utility.
*/

#include <assert.h>

#include <iostream>
#include <limits>
#include <memory>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/compile_time.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/allowed_elements.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/chiral_centre.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mdl.h"
#include "Molecule_Lib/numass.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "fileconv_opts.h"
#include "remove_duplicate_fragments.h"

namespace fileconv {

using std::cerr;

FileconvConfig::FileconvConfig() { DefaultValues(); }

void
FileconvConfig::DefaultValues() {
  audit_input = 0;
  verbose = 0;
  print_bond_lengths = 0;
  print_bond_angles = 0;
  print_torsions = 0;
  print_max_atom_separation = 0;
  acc_longest_sep.reset();
  molecules_processed = 0;
  molecules_changed = 0;
  molecule_changed_string.resize_keep_storage(0);

  chemical_standardisation.deactivate();
  charge_assigner.deactivate();
  donor_acceptor_assigner.deactivate();

  fragment_count = 0;
  reduce_to_largest_fragment = 0;
  reduce_to_all_largest_fragments = 0;
  reduce_to_largest_organic_fragment = 0;
  reduce_to_largest_organic_fragment_carefully = 0;
  keep_all_organic_fragments = false;
  molecules_with_nonorganic_fragments_removed = 0;
  maximum_fragment_count = 0;
  min_size_max_fragment_count = 0;
  molecules_chopped = 0;
  molecules_with_too_many_components = 0;
  molecules_with_large_fragments = 0;
  remove_duplicate_fragments = 0;
  molecules_with_duplicate_fragments_removed = 0;
  remove_fragments_this_size_or_smaller = 0;
  remove_largest_fragment = 0;
  keep_smallest_fragments = 0;
  remove_smallest_fragment = 0;
  molecules_with_very_small_fragments_removed = 0;
  remove_molecules_with_non_largest_fragment_natoms = -1;
  strip_to_n_largest_fragments = 0;
  sort_by_fragment_size = 0;
  tag_for_removed_fragments.resize_keep_storage(0);

  known_fragment_data.deactivate();

  atoms_lost.reset();
  initial_fragment_count.resize(0);

  need_to_call_process_fragments = 0;

  skip_molecules_with_abnormal_valences = 0;
  ok_bad_valence_on_isotopically_labelled = 0;
  molecules_with_abnormal_valences = 0;

  max_path_length = 0;
  molecules_with_longest_path_too_long = 0;

  need_to_consider_isotopes = 0;

  exclude_isotopes = 0;
  molecules_containing_isotopes = 0;

  convert_isotopes = 0;
  convert_specific_isotopes.resize(0);
  convert_specific_isotopes_new_isotope.resize(0);
  convert_all_isotopes_to = 0;

  convert_isotopes_to_atom_map_numbers = 0;
  convert_atom_map_numbers_to_isotopes = 0;

  convert_specific_isotopes_query.resize(0);
  convert_specific_isotopes_query_new_isotope.resize(0);

  remove_isotopic_atoms = 0;
  isotope_is_atom_type = 0;

  output_organic_only = 0;
  non_organic_molecules_found = 0;

  exclude_non_real_elements = 0;
  non_real_elements_found = 0;

  allowed_elements.reset_to_defaults();

  molecules_excluded_for_non_allowed_elements = 0;

  lower_molecular_weight_cutoff = 0.0;
  upper_molecular_weight_cutoff = 0.0;

  molecules_below_molecular_weight_cutoff = 0;
  molecules_above_molecular_weight_cutoff = 0;

  amw_accumulator.reset();
  natoms_accumulator.reset();
  atom_count.resize(0);

  lower_atom_count_cutoff = 0;
  upper_atom_count_cutoff = 0;

  // int lower_atom_count_cutoff_applies_to_largest_fragment = 0;
  // int upper_atom_count_cutoff_applies_to_largest_fragment = 0;

  include_implicit_hydrogens_in_upper_atom_count_comparison = 0;

  atom_count_includes_only_atoms_in_largest_fragment = 0;

  molecules_below_atom_count_cutoff = 0;
  molecules_above_atom_count_cutoff = 0;

  name_rx.reset(nullptr);
  grep_v_name_rx.reset(nullptr);
  molecules_discarded_for_name_mismatch = 0;

  smallest_fragment_queries.resize(0);
  largest_fragment_queries.resize(0);
  keep_fragment_queries.resize(0);
  remove_fragment_queries.resize(0);

  molecules_with_fragments_reduced_by_query = 0;
  molecules_not_matching_fragment_queries = 0;

  lower_ring_count_cutoff = 0;
  molecules_with_too_few_rings = 0;
  upper_ring_count_cutoff = -1;
  molecules_with_too_many_rings = 0;

  ring_systems_include_spiro = 0;
  min_ring_systems_required = -1;
  max_ring_systems_allowed = -1;
  molecules_with_too_many_ring_systems = 0;
  molecules_with_too_few_ring_systems = 0;

  min_aliphatic_ring_count = 0;
  max_aliphatic_ring_count = -1;
  min_aromatic_ring_count = 0;
  max_aromatic_ring_count = -1;

  molecules_with_too_few_aliphatic_rings = 0;
  molecules_with_too_few_aromatic_rings = 0;
  molecules_with_too_many_aliphatic_rings = 0;
  molecules_with_too_many_aromatic_rings = 0;

  dx = 0.0;
  dy = 0.0;
  dz = 0.0;

  translation_specified = 0;
  translate_to_origin = 0;

  fileconv_partial_charge_type = kNoChargeCalculation;

  number_assigner.deactivate();

  elements_to_remove.resize(0);

  element_transformations.resize(0);

  make_all_implicit_hydrogens_explicit = 0;

  atoms_for_implicit_hydrogens.resize(0);

  molecules_to_which_hydrogens_were_added = 0;

  hydrogens_last = 0;

  find_all_chiral_centres = 0;
  find_all_ring_chiral_centres = 0;
  invert_all_chiral_centres = 0;
  reflect_coordinates = 0;
  chiral_centres_inverted = 0;
  remove_invalid_chiral_centres = 0;
  molecules_with_invalid_chiral_centres = 0;
  molecules_with_chiral_centres = 0;
  remove_chiral_data_from_all_molecules = 0;
  remove_non_organic_chirality = 0;
  max_chiral_centres = 0;
  molecules_with_too_many_chiral_centres = 0;
  chiral_centre_count.resize(0);

  remove_chiral_centres_on.resize(0);

  chiral_centres_removed_by_rmchiral = 0;
  molecules_with_chiral_centres_removed_by_rmchiral = 0;

  remove_directional_bonds_from_input = 0;

  remove_invalid_directional_bonds_from_input = 0;

  ok_non_organics.resize(0);

  ok_boron_special_case = 0;

  // Not much we can do to reset these. Maybe close?
  // reject_log;
  // rejections_output_object;

  compute_molecular_weight_for_each = 0;

  compute_molecular_weight_based_on_largest_fragment = 0;

  appends_to_be_done = 0;

  do_appends_as_prepends = 0;

  appended_properties_from_largest_fragment = 0;

  append_molecular_weight_to_name = 0;

  append_molecular_weight_ok_isotope_to_name = 0;

  append_exact_mass_to_name = 0;

  append_heteratom_count_to_name = 0;

  append_molecular_formula_to_name = 0;

  append_isis_molecular_formula_to_name = 0;

  append_nrings_to_name = 0;

  append_aromatic_ring_count_to_name = 0;

  append_natoms_to_name = 0;

  append_nbonds_to_name = 0;

  append_net_formal_charge_to_name = 0;

  append_clnd_count_to_name = 0;

  lower_amw_cutoff = -1.0;
  molecules_below_amw_cutoff = 0;
  upper_amw_cutoff = -1.0;
  molecules_above_amw_cutoff = 0;

  substitute_for_whitespace_in_name = "";
  truncate_names_to_first_token = 0;
  truncate_name_at_first = '\0';
  truncate_name_at_last = '\0';

  name_token = -1;

  IWString prepend_to_name = "";

  aromatise_these_rings.resize(0);

  molecules_changed_by_aromatising_rings = 0;

  remove_unnecessary_square_brackets = 0;

  remove_all_possible_square_brackets = 0;

  molecules_with_alternate_valence_via_implicit_h = 0;

  reset_atom_map_numbers = 0;

  molecule_to_fragments = 0;

  return;
}

void
Usage(int rc = 1) {
  // clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "The following options are recognised\n";
  cerr << "  -f <directive> fragment selection, enter '-f help' for details\n";
  cerr << "  -F <number>    exclude molecules having more than <number> fragments\n";
  cerr << "  -O none        \"organic\" molecules only. '-O def' removes all bad elements\n";
  cerr << "  -O <el>        \"organic\" molecules only, but 'el' is OK (repeat for each OK ele)\n";
  cerr << "  -E <symbol>    create element <symbol>, use 'autocreate' for all\n";
  cerr << "  -w <amw, -W <amw>  specify lower(-w) and upper (-W) amw limits\n";
  cerr << "  -W LARGE       compute molecular weight based on largest fragment\n";
  cerr << "  -c <number>    exclude molecules with atom count below <number>\n";
  cerr << "  -C <number>    exclude molecules with atom count above <number>\n";
  cerr << "  -C implicit    when computing atom count for -C, include implicit Hydrogens\n";
  cerr << "  -r <n>         omit molecules with fewer than <n> rings\n";
  cerr << "  -R <n>         omit molecules with more  than <n> rings. Use MRS:n for max rings in a system\n";
  cerr << "                 both -r and -R take S:n for ring systems. Use 'spiro' to span spiro fusions\n";
  cerr << "  -m <ringsize>  discard molecules with any ring larger than <ringsize>\n";
  cerr << "  -X <symbol>    extract/remove all atoms of type <symbol>. No bonds changed\n";
  cerr << "  -Q <type>      compute partial charges - enter '-Q help' for info\n";
  cerr << "  -I ...         isotope options, enter '-I help' for info'\n";
  cerr << "  -s <qual>      chirality options, enter '-s help' for details\n";
  cerr << "  -e             discard molecules with non periodic table elements\n";
  cerr << "  -n ...         number assigner options, enter '-n help' for more info\n";
  cerr << "  -o <type>      specify output file type(s), enter '-o help' for details\n";
  cerr << "  -a             audit input only, no output\n";
  cerr << "  -V             skip any molecule with abnormal valences\n";
  cerr << "  -T <dx,dy,dz>  translate atoms (dx,dy,dz), '-T origin' move to 0,0,0\n";
  cerr << "  -L <fname>     write rejected molecules to <fname>\n";
  cerr << "  -B ...         handling for otherwise fatal problems, enter '-B help'\n";
  cerr << "  -h ...         options relating to hydrogen atoms, '-h help' for info\n";
#ifdef COMPILE_CHANGE_MOLECULE
  cerr << "  -P             invoke compiled-in molecule changer\n";
#endif
  cerr << "  -S <string>    create output files with name stem <string>\n";
  cerr << "  -p <string>    append various things to the name. Enter '-p help' for info\n";
  cerr << "  -J <...>       fix obvious structure errors, enter '-J help' for info\n";
  cerr << "  -Y <...>       miscellaneous options, enter '-Y help' for info\n";
  cerr << "  -t E1=E2       standard element transformation options, enter '-t help'\n";
  cerr << "  -i <type>      specify input file type. Enter '-i help' for details\n";
  (void)display_standard_charge_assigner_options(cerr, 'N');
  cerr << "  -H <..>        donor acceptor assignment, enter '-H help' for info\n";
  (void)display_standard_chemical_standardisation_options(cerr, 'g');
  (void)display_standard_aromaticity_options(cerr);
  (void)display_standard_smiles_options(cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

void
DisplayDashhOptions(std::ostream& output) {
  // clang-format off
  output << " -h <all>       make implicit hydrogens explicit\n";
  output << " -h last=<z>    move atoms with atomic number <z> to end of the connection table\n";
  output << " -h <query>     make implicit hydrogens on atoms matching <query> explicit\n";
  output << "                -h SMARTS:smt, -h PROTO:query.textproto also recognised\n";
  // clang-format on

  ::exit(0);
}

void
DisplayfOptions(int rc) {
  // clang-format off
  cerr << "  -f large       trim to largest fragment (can abbreviate to '-f l')\n";
  cerr << "  -f alarge      determine largest fragment. Keep all fragments with that number of atoms\n";
  cerr << "  -f lo          trim to fragment with most organic atoms\n";
  cerr << "  -f lod         trim to fragment with most organic atoms and desirable features\n";
  cerr << "  -f allo        keep all organic fragments\n";
  cerr << "  -f RMDUP       remove duplicate fragments\n";
  cerr << "  -f rmle=nn     remove fragments with NN or fewer atoms\n";
  cerr << "  -f Q:qfile     keep largest  frag which matches query in <qfile>\n";
  cerr << "  -f q:qfile     keep smallest frag which matches query in <qfile>\n";
  cerr << "  -f Q:F:file    keep largest  frag which matches queries in <file>\n";
  cerr << "  -f q:F:file    keep smallest frag which matches queries in <file>\n";
  cerr << "  -f SMARTS:smt  keep largest  frag which matches smarts\n";
  cerr << "  -f smarts:smt  keep smallest frag which matches smarts\n";
  cerr << "  -f ALL:smt     keep all fragments that match smarts <smarts>\n";
  cerr << "  -f rm:smt      remove all fragments that match <smarts>\n";
  cerr << "  -f saltfile=<file> smiles file of known salts - always removed even\n";
  cerr << "                     if the largest fragment\n";
  cerr << "  -f parentfile=<file> file of known parent molecules - never removed as salts\n";
  cerr << "  -f kmfok       compare known salts and parents by molecular formula only - not unique smiles\n";
  cerr << "  -f kpallsalt   do not change a molecule if every fragment is a known salt\n";
  cerr << "  -f rmxt=<n>    discard molecules with >1 fragment with more than n atoms\n";
  cerr << "  -f rmxt        discard molecules with >1 fragment with more than 16 atoms\n";
  cerr << "  -f sfs         sort fragments by size\n";
  cerr << "  -f dmxt=<d>    discard molecules where largest fragments differ by <d> atoms or fewer\n";
  cerr << "  -f manlf=<d>   discard molecules that have a non-largest fragment with more than <d> atoms\n";
  cerr << "  -f klf=<d>     discard all but the <n> largest fragments\n";
  cerr << "  -f RMF=<tag>   when processing TDT forms, write removed fragments to <tag>\n";
  cerr << "  -f rmlarge     remove the largest fragment (arbitrary if two frags of same size\n";
  cerr << "  -f rmlarge=<n> remove the largest <n> fragments (arbitrary if frags of the same size)\n";
  cerr << "  -f rmsmall     remove the smallst fragment (arbitrary if two frags of same size\n";
  cerr << "  -f rmsmall=<n> remove the smallst <n> fragments (arbitrary if frags of the same size)\n";
  cerr << "  -f keepsmall   remove all but the smallest fragment (arbitrary if frags of the same size)\n";
  cerr << "  -f keepsmall=<n> remove all but the smallest <n> fragments\n";
  cerr << "  -f <number>    remove fragments so that all written molecules\n";
  cerr << "                 have no more than <number> fragments\n";
  // clang-format on


  exit(rc);
}

void
DisplayFOptions(int rc) {
  // clang-format off
  cerr << "  -F <number>      discard molecules with more than <number> fragments\n";
  cerr << "  -F mnsz=<n>      when counting fragments only count those with > <n> atoms\n";
  // clang-format on

  exit(rc);
}

void
DisplayPDirectives(int rc) {
  // clang-format off
  cerr << " Appends various molecular properties to the molecule name\n";
  cerr << "  -p AMW         append molecular weight to the name (isotopic atoms fail molecule)\n";
  cerr << "  -p AMWI        append molecular weight to the name (isotopic atoms handled)\n";
  cerr << "  -p MF          append molecular formula to the name\n";
  cerr << "  -p ISISMF      append ISIS-like molecular formula to the name\n";
  cerr << "  -p NATOMS      append number of atoms to the name\n";
  cerr << "  -p NRINGS      append number of rings to the name\n";
  cerr << "  -p AROMR       append number of aromatic rings to the name\n";
  cerr << "  -p HTROAC      append number of heteroatoms to the name\n";
  cerr << "  -p EXACT       append exact mass to the name\n";
  cerr << "  -p NFC         append net formal charge to the name\n";
  cerr << "  -p CLND        append Chemiluminescent Nitrogen Detection\n";
  cerr << "  -p PREPEND     do a prepend rather than append\n";
  cerr << "  -p LARGE       computed properties derived from the largest fragment\n";
  // clang-format on

  exit(rc);
}

void
DisplayDashYOptions(std::ostream& os, int rc) {
  // clang-format off
  os << "-Y nbvm          No Bad Valence Messages, or strange electron counts\n";
  os << "-Y okbvi         during valence check, ok to have bad valence on isotopes\n";
  os << "-Y appchg=xxxx   append 'xxxx' to the name of molecules that are changed\n";
  os << "-Y pblen         print all bond lengths in the molecules\n";
  os << "-Y pbang         print all bond angles in the molecules\n";
  os << "-Y ptor          print all torsion angles in the molecules\n";
  os << "-Y pmaxd         print the max interatomic distance in each molecule\n";
  os << "-Y dbg           debug print each molecule\n";
  os << "-Y namerx=<rx>   discard molecules unless the molecule name matches <rx>\n";
  os << "-Y grepvname=<rx> discard molecules if the molecule name matches <rx>\n";
  os << "-Y ftn           keep only the first token in molecule names\n";
  os << "-Y nsubw=c       translate all whitespace in molecule names to 'c'\n";
  os << "-Y chname=rx     change name to what is matched by <rx> (optional match)\n";
  os << "-Y CHNAME=rx     change name to what is matched by <rx> (MUST match)\n";
  os << "-Y tfirst=char   truncate name at first <char>\n";
  os << "-Y tlast=char    truncate name at last <char>\n";
  os << "-Y NPREPEND=s    prepend <s> to each name\n";
  os << "-Y ntoken=n      the output name is word 'n' in the input name\n";
  os << "-Y maxplen=<n>   discard molecules with max path length > <n>\n";
  // os << "-Y B4F=<fname>   write frag stripped smiles before filtering to <fname>\n";
  os << "-Y aclf          atom counts are for the largest fragment only\n";
  os << "-Y nhsqb         explicit hydrogen atoms in smiles written without square brackets\n";
  os << "-Y rmsqb         remove unnecessary square bracketed atoms  - Hcount is OK as specified\n";
  os << "-Y FHrmsqb       in atoms like [C] free the H count to what is computed. Square brackets removed\n";
  os << "-Y Xihaltvalence discard molecules like =S- where an implicit hydrogen satisfies an alternate valence\n";
  os << "-Y rmamap        remove atom map numbers\n";
  os << "-Y fixarom=smarts call find_kekule_form on the ring (system) matched by the first atom in <smarts>\n";
  os << "-Y rmatom=smarts remove all atoms that match <smarts>\n";
  os << "-Y zpad=width    left pad the name with 0's to <width>. Negative numbers remove leading chars from name\n";
  os << "-Y rot=x,y,z,angle rotate molecules `angle` degrees around x,y,z\n";
  os << "-Y atype=<atype> atom typing specification - use with '-I atype'\n";
  os << "-Y  help          this message\n";
  // clang-format on

  exit(rc);
}

void
DisplayDashIOptions(char flag, std::ostream& os) {
  // clang-format off
  os << " -" << flag << " 0             discard molecules containing any isotopic atoms\n";
  os << " -" << flag << " change        change any isotopic atoms to normal form\n";
  os << " -" << flag << " change=<n>    change any isotope <n> atoms to normal form\n";
  os << " -" << flag << " change=<i,j>  change any isotope <i> atoms to isotope <j>\n";
  os << " -" << flag << " CHANGE        change to normal form. Free implicit H. Should be default\n";
  os << " -" << flag << " alliso=<i>    change any isotopic atoms to isotope <i>\n";
  os << " -" << flag << " smarts:<i>    change any atoms matching <smarts> to isotope <i>, e.g. '[2C],0' or '[#16],4'\n";
  os << " -" << flag << " remove        remove any isotopically labelled atoms\n";
  os << " -" << flag << " help          this message\n";
  // clang-format on

  exit(1);
}

void
DisplayDashOOptions(char flag, std::ostream& os) {
  // clang-format off
  os << "Specifies elements allowed as ''organic'\n";
  os << "Operates in two different modes\n";
  os << "In the historical mode, fragment stripping is done, then remaining atoms\n";
  os << "checked for organic property\n";
  os << "to enable this use '-" << flag << " none, or -" << flag << " ele' to make elements organic\n";
  os << "The other mode is to first check the entire molecule for disallowed elements\n";
  os << "After fragment stripping is done, then check the largest fragment for non-organic elements\n";
  // clang-format on

  exit(1);
}


// Recursively identify a contiguous group of Nitrogen atoms.
int
identify_ngroup(const Molecule& m, atom_number_t n, int* ngroup) {
  ngroup[n] = n + 1;

  int rc = 1;

  const Atom* a = m.atomi(n);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++) {
    atom_number_t j = a->other(n, i);

    if (ngroup[j])
      continue;

    if (7 != m.atomic_number(j))
      continue;

    rc += identify_ngroup(m, j, ngroup);
  }

  return rc;
}

int
IsAzide(const Molecule& m, atom_number_t n1, int* ngroup) {
  const Atom* an1 = m.atomi(n1);

  atom_number_t n2 = an1->other(n1, 0);

  const Atom* an2 = m.atomi(n2);

  if (7 != an2->atomic_number())
    return 0;

  if (2 != an2->ncon())
    return 0;

  atom_number_t n3 = an2->other(n2, 0);
  if (n3 == n1)
    n3 = an2->other(n2, 1);

  if (7 != m.atomic_number(n3))
    return 0;

  if (2 != m.ncon(n3))
    return 0;

  // We have 3 nitrogens. They must be [N-]=[N+]=N or N#N=N

  const Bond* b12 = an1->item(0);
  const Bond* b23 = an2->item(0);
  if (n1 == b23->other(n2))
    b23 = an2->item(1);

  if (!b23->is_double_bond())
    return 0;

  if (b12->is_triple_bond())
    ;
  else if (b12->is_double_bond() && -1 == an1->formal_charge() && 1 == an2->formal_charge())
    ;
  else
    return 0;

  ngroup[n1] = n1 + 1;
  ngroup[n2] = n1 + 1;
  ngroup[n3] = n1 + 1;

  return 1;
}

int
FileconvConfig::ComputeClnd(const Molecule& m) {
  int rc = 0;

  int matoms = m.natoms();

  int* ngroup = new_int(matoms);
  std::unique_ptr<int[]> free_ngroup(ngroup);

  // First look for [N-]#[N+]=N because they count for 1 each (why I have no idea).

  for (int i = 0; i < matoms; i++) {
    if (ngroup[i])  // already counted as part of another group
      continue;

    if (7 != m.atomic_number(i))
      continue;

    if (1 != m.ncon(i))
      continue;

    if (IsAzide(m, i, ngroup))
      rc += 1;
  }

  for (int i = 0; i < matoms; i++) {
    if (ngroup[i])
      continue;

    if (7 != m.atomic_number(i))
      continue;

    int gsize = identify_ngroup(m, i, ngroup);

    rc += gsize / 2;
    if (1 == gsize % 2)
      rc++;
  }

  return rc;
}

void
FileconvConfig::DoAppends(Molecule& m, IWString& extra_stuff) {
  if (append_molecular_formula_to_name) {
    extra_stuff += " MF = ";
    extra_stuff += m.molecular_formula();
  }

  if (append_isis_molecular_formula_to_name) {
    extra_stuff += " ISISMF = ";
    IWString tmp;
    m.isis_like_molecular_formula_dot_between_fragments(tmp);
    extra_stuff += tmp;
  }

  if (append_natoms_to_name)
    extra_stuff << " NATOMS = " << m.natoms();

  if (append_nbonds_to_name)
    extra_stuff << " NBONDS = " << m.nedges();

  if (append_nrings_to_name)
    extra_stuff << " NRINGS = " << m.nrings();

  if (append_aromatic_ring_count_to_name) {
    m.compute_aromaticity_if_needed();

    int nr = m.nrings();
    int nar = 0;
    for (int i = 0; i < nr; i++) {
      const Ring* ri = m.ringi(i);
      if (ri->is_aromatic())
        nar++;
    }

    extra_stuff << " AROMRING = " << nar;
  }

  if (append_molecular_weight_to_name)
    extra_stuff << " AMW = " << m.molecular_weight();
  else if (append_molecular_weight_ok_isotope_to_name) {
    Molecular_Weight_Control mwc;
    Molecular_Weight_Calculation_Result mwcr;
    mwc.set_ignore_isotopes(0);
    (void)m.molecular_weight(mwc, mwcr);
    extra_stuff << " AMW = " << mwcr.amw();
  }

  if (append_exact_mass_to_name) {
    exact_mass_t x;
    if (!m.exact_mass(x))
      cerr << "Warning, molecule '" << m.name() << "' partial result for exact mass\n";

    extra_stuff << " EXACT_MASS = " << x;
  }

  if (append_heteratom_count_to_name) {
    extra_stuff << " HTROAC = " << (m.natoms() - m.natoms(6) - m.natoms(1));
  }

  if (append_net_formal_charge_to_name)
    extra_stuff << " FORMAL_CHARGE = " << m.formal_charge();

  if (append_clnd_count_to_name)
    extra_stuff << " CLND = " << ComputeClnd(m);

  return;
}

void
FileconvConfig::DoAppends(Molecule& m) {
  const IWString& mname = m.name();

  IWString extra_stuff;
  extra_stuff.resize(mname.length() + 200);

  if (!do_appends_as_prepends)  // we are appending, existing name at the beginning
    extra_stuff << mname;

  if (appended_properties_from_largest_fragment && m.number_fragments() > 1) {
    Molecule tmp(m);
    tmp.reduce_to_largest_fragment();

    DoAppends(tmp, extra_stuff);
  } else
    DoAppends(m, extra_stuff);

  if (do_appends_as_prepends)  // name goes at the end
    extra_stuff << ' ' << mname;

  m.set_name(extra_stuff);

  return;
}

int
FileconvConfig::RemoveFragmentsThisSizeOrSmaller(Molecule& m) {
  int nf = m.number_fragments();

  int matoms = m.natoms();

  Set_of_Atoms atoms_to_be_removed;
  atoms_to_be_removed.resize(matoms);

  int number_fragments_removed = 0;

  for (int i = 0; i < nf; i++) {
    int aif = m.atoms_in_fragment(i);

    if (aif <= remove_fragments_this_size_or_smaller) {
      number_fragments_removed++;

      for (int j = 0; j < matoms; j++) {
        int f = m.fragment_membership(j);
        if (i == f)
          atoms_to_be_removed.add(j);
      }
    }
  }

  if (atoms_to_be_removed.empty())
    return 1;

  if (matoms == atoms_to_be_removed.number_elements()) {
    cerr << m.name() << " all fragments have " << remove_fragments_this_size_or_smaller
         << " or fewer atoms\n";
  } else {
    if (verbose > 1)
      cerr << "Removed " << number_fragments_removed << " fragments and "
           << atoms_to_be_removed.number_elements() << " atoms from '" << m.name() << "'\n";

    molecules_with_very_small_fragments_removed++;

    return m.remove_atoms(atoms_to_be_removed);
  }

  return 1;
}

int
FileconvConfig::ComputePartialCharges(Molecule& m) {

  switch (fileconv_partial_charge_type) {
    case kNoChargeCalculation: {
      return 1;
    }
    case kGasteiger: {
      return m.compute_Gasteiger_partial_charges();
      break;
    }
    case kGasteigerHuckel: {
      return m.compute_Gasteiger_Huckel_partial_charges();
      break;
    }
    case kAbraham: {
      return m.compute_Abraham_partial_charges();
      break;
    }
  }

  return 0;  // should never happen.
}

void
AppendBond(const Bond& b, std::ostream& output) {
  if (b.is_aromatic())
    output << ':';
  else if (b.is_single_bond())
    output << '-';
  else if (b.is_double_bond())
    output << '=';
  else if (b.is_triple_bond())
    output << '#';
  else
    output << '?';

  if (b.is_aromatic())
    ;
  else if (b.nrings())
    output << 'R';

  return;
}

// Report max interatomic separation within `m`.
int
FileconvConfig::PrintMaxAtomSeparation(const Molecule& m, std::ostream& output) {
  const int matoms = m.natoms();
  if (matoms < 2) {
    return 1;
  }
  float max_sep = 0.0;
  atom_number_t max1 = INVALID_ATOM_NUMBER;
  atom_number_t max2 = INVALID_ATOM_NUMBER;
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      const float d = m.distance_between_atoms(i, j);
      if (d > max_sep) {
        max_sep = d;
        max1 = i;
        max2 = j;
      }
    }
  }
  acc_longest_sep.extra(max_sep);
  output << m.name() << " max sep " << max_sep << " btw atoms " << max1 << ' '
         << m.smarts_equivalent_for_atom(max1) << " and" << max2 << ' '
         << m.smarts_equivalent_for_atom(max2) << '\n';
  return output.good();
}

int
FileconvConfig::PrintTorsion(const Molecule& m, const Bond& b, std::ostream& output) {
  const atom_number_t a2 = b.a1();
  const atom_number_t a3 = b.a2();

  const Atom* aa2 = m.atomi(a2);
  const Atom* aa3 = m.atomi(a3);

  const int a2con = aa2->ncon();

  if (1 == a2con)
    return output.good();

  const int a3con = aa3->ncon();

  if (1 == a3con)
    return output.good();

  IWString a2smarts = m.smarts_equivalent_for_atom(a2);
  IWString a3smarts = m.smarts_equivalent_for_atom(a3);

  char sep = ' ';

  for (int i = 0; i < a2con; i++) {
    const Bond* b21 = aa2->item(i);

    const atom_number_t a1 = b21->other(a2);

    if (a1 == a3)
      continue;

    const IWString a1smarts = m.smarts_equivalent_for_atom(a1);

    for (int k = 0; k < a3con; k++) {
      const Bond* b34 = aa3->item(k);

      const atom_number_t a4 = b34->other(a3);

      if (a4 == a2)
        continue;

      const IWString a4smarts = m.smarts_equivalent_for_atom(a4);

      const angle_t tmp = m.dihedral_angle(a1, a2, a3, a4);
      output << m.name() << sep << "dihedral " << a1 << sep << a2 << sep << a3 << sep << a4 << sep
             << tmp << " (" << (tmp * RAD2DEG) << " deg)";
      output << sep << m.atomic_number(a1) << sep << a1smarts << sep;
      AppendBond(*b21, output);
      output << sep << m.atomic_number(a2) << sep << a2smarts << sep;
      AppendBond(b, output);
      output << sep << m.atomic_number(a3) << sep << a3smarts << sep;
      AppendBond(*b34, output);
      output << sep << m.atomic_number(a4) << sep << a4smarts;

      output << "\n";
    }
  }

  return output.good();
}

int
FileconvConfig::PrintTorsions(Molecule& m, std::ostream& output) {
  m.compute_aromaticity_if_needed();

  if (m.highest_coordinate_dimensionality() < 3)
    cerr << "do_print_torsions: WARNING not 3D " << m.name() << '\n';

  const int nb = m.nedges();

  for (int i = 0; i < nb; i++) {
    const Bond* b = m.bondi(i);

    if (!PrintTorsion(m, *b, output))
      return 0;
  }

  return output.good();
}

int
FileconvConfig::PrintBondAngle(const Molecule& m,
                                const atom_number_t a1,
                                const atom_number_t a2,
                                std::ostream& output) {
  const Atom* aa2 = m.atomi(a2);

  const int ncon2 = aa2->ncon();

  for (int i = 0; i < ncon2; i++) {
    const atom_number_t a3 = aa2->other(a2, i);

    if (a3 == a1)
      continue;

    if (a3 < a1)
      continue;

    const angle_t theta = m.bond_angle(a1, a2, a3);
    output << " atoms " << a1 << " " << a2 << " " << a3 << " angle " << theta << " ("
           << (theta * RAD2DEG) << " deg)\n";
  }

  return output.good();
}

int
FileconvConfig::PrintBondAngles(const Molecule& m, std::ostream& output) {
  const int nb = m.nedges();

  for (int i = 0; i < nb; i++) {
    const Bond* b = m.bondi(i);

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    PrintBondAngle(m, a1, a2, output);
    PrintBondAngle(m, a2, a1, output);
  }

  return output.good();
}

int
FileconvConfig::PrintBondLengths(const Molecule& m, std::ostream& output) {
  float max_bond_length = 0.0f;
  atom_number_t max1 = INVALID_ATOM_NUMBER;
  atom_number_t max2 = INVALID_ATOM_NUMBER;
  // In malformed molecules, it can be useful to know
  // the number of possibly questionable bonds.
  constexpr float threshold = 2.0f;
  int number_above_threshold = 0;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; i++) {
    const Atom* ai = m.atomi(i);
    const int icon = ai->ncon();
    for (int j = 0; j < icon; j++) {
      const Bond* b = ai->item(j);

      const atom_number_t k = b->other(i);

      const Atom* ak = m.atomi(k);
      output << "Atom " << i << " (" << ai->atomic_symbol() << ") ";
      if (b->is_aromatic())
        output << "aromatic";
      else if (b->is_single_bond())
        output << "single";
      else if (b->is_double_bond())
        output << "double";
      else if (b->is_triple_bond())
        output << "triple";
      else
        output << "huh";

      const float d = ai->distance(*ak);
      output << " bond to " << k << " (" << ak->atomic_symbol() << ") dist = " << d << '\n';
      if (d > threshold) {
        number_above_threshold++;
      }
      if (d > max_bond_length) {
        max_bond_length = d;
        max1 = i;
        max2 = k;
      }
    }
  }

  output << "max dist " << max_bond_length << " btw " << max1 << ' '
         << m.smarts_equivalent_for_atom(max1) << " and " << max2 << ' '
         << m.smarts_equivalent_for_atom(max2) << '\n';
  if (number_above_threshold) {
    output << number_above_threshold << " bond lengths above " << threshold << '\n';
  }

  return output.good();
}

int
FileconvConfig::InvertAllChiralCentres(Molecule& m) {
  const int nc = m.chiral_centres();
  if (0 == nc)
    return 1;

  int rc = 0;

  for (int i = 0; i < nc; i++) {
    Chiral_Centre* c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);

    c->invert();
    if (0 == rc)
      molecules_with_chiral_centres++;
  }

  chiral_centres_inverted += nc;

  return nc;
}

int
FileconvConfig::ReflectCoordinates(Molecule& m) {
  if (m.highest_coordinate_dimensionality() < 3) {
    cerr << "do_reflect_coordinates:inadequate dimensionality "
         << m.highest_coordinate_dimensionality() << ", '" << m.name() << "'\n";
    return 0;
  }

  // coord_t xmin, xmax, ymin, ymax, zmin, zmax;
  // m.spatial_extremeties (xmin, xmax, ymin, ymax, zmin, zmax);

  const auto matoms = m.natoms();

  // m.translate_atoms(0.0, 0.0, -zmin);    // ensure all on one side of xy plane

  for (auto i = 0; i < matoms; ++i) {
    Coordinates c;
    m.get_coords(i, c);
    m.setxyz(i, -c.x(), -c.y(), -c.z());
  }

  // m.translate_atoms(0.0, 0.0, zmin);

  return InvertAllChiralCentres(m);
}

int
FileconvConfig::FindAllChiralCentres(Molecule& m) {
  int matoms = m.natoms();

  const int* symmetry = m.symmetry_classes();

  int rc = 1;
  for (int i = 0; i < matoms; i++) {
    if (nullptr != m.chiral_centre_at_atom(i))
      continue;

    if (find_all_ring_chiral_centres && m.is_non_ring_atom(i))
      continue;

    Atom* a = const_cast<Atom*>(m.atomi(i));

    int lp;

    if (4 == a->ncon())
      ;
    else if (3 == a->ncon() && 1 == a->implicit_hydrogens())
      ;
    else if (3 == a->ncon() && a->lone_pair_count(lp) && 1 == lp)
      ;
    else
      continue;

    resizable_array<int> symmetries_found;
    symmetries_found.resize(a->ncon());

    int different_symmetries = 1;
    for (int j = 0; j < a->ncon(); j++) {
      atom_number_t k = a->other(i, j);
      if (0 == symmetries_found.add_if_not_already_present(symmetry[k])) {
        different_symmetries = 0;
        break;
      }
    }

    if (!different_symmetries)
      continue;

    if (verbose > 1)
      cerr << "Placing chiral centre at atom " << i << " (" << a->atomic_symbol() << ") "
           << a->ncon() << " connections\n";

    if (nullptr == m.create_chiral_centre(i)) {
      cerr << "Yipes, could not place chiral center atom atom " << i << '\n';
      rc = 0;
    }
  }

  return rc;
}

int
FileconvConfig::RemoveIsotopicAtoms(Molecule& m) {
  int rc = 0;

  for (int i = m.natoms() - 1; i >= 0; --i) {
    if (0 == m.isotope(i))
      continue;

    m.remove_atom(i);
    rc++;
  }

  if (rc)
    molecules_containing_isotopes++;

  return rc;
}

// Determine atom types and set isotopes.
// Note that there is potential for problems here.
// Isotopes are signed integers, but atom types
// can be unsigned. But I am not anxious to change
// the type of isotopes.
int
FileconvConfig::AtomTypeToIsotope(Molecule& m) {
  if (! atom_typing.active()) {
    cerr << "FileconvConfig::AtomTypeToIsotope:atom typing not active\n";
    return 0;
  }

  const int matoms = m.natoms();

  std::unique_ptr<atom_type_t[]> atype = std::make_unique<atom_type_t[]>(matoms);

  atom_typing.assign_atom_types(m, atype.get());

  for (int i = 0; i < matoms; ++i) {
    m.set_isotope(i, atype[i]);
  }

  return 1;
}

int
FileconvConfig::ConvertIsotopesToAtomMapNumbers(Molecule& m) {
  int rc = 0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    const int iso = m.isotope(i);

    if (0 == iso)
      continue;

    m.set_atom_map_number(i, iso);
    rc++;
  }

  return rc;
}

int
FileconvConfig::ConvertAtomMapNumbersToIsotopes(Molecule& m) {
  int rc = 0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    const int amap = m.atom_map_number(i);

    if (0 == amap)
      continue;

    m.set_isotope(i, amap);
    rc++;
  }

  return rc;
}

int
FileconvConfig::MakeImplicitHydrogensExplicit(Molecule& m) {
  int rc = 0;  // the number of atoms to which we add a Hydrogen

  Molecule_to_Match target(&m);

  for (Substructure_Query* q : atoms_for_implicit_hydrogens) {
    Substructure_Results sresults;

    if (q->substructure_search(target, sresults) == 0) {
      continue;
    }

    for (const Set_of_Atoms* e : sresults.embeddings()) {
      const atom_number_t a = e->item(0);  // we process only the first matched atom

      if (0 == m.implicit_hydrogens(a))
        continue;

      rc += m.make_implicit_hydrogens_explicit(a);
    }
  }

  return rc;
}

/*
  Note that the next two functions have been copied verbatim to id_chirality.cc
*/

int
FileconvConfig::IdentifyMatchedAtomsWithChiralCentres(
    const Molecule& m,
    const Set_of_Atoms& e,
    Set_of_Atoms& atoms_with_chiral_centres_to_be_removed) {
  int n = e.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++) {
    atom_number_t j = e[i];

    if (nullptr == m.chiral_centre_at_atom(j))  // atom J does not have a chiral centre
      continue;

    if (atoms_with_chiral_centres_to_be_removed.contains(j))
      continue;

    atoms_with_chiral_centres_to_be_removed.add(j);
    rc++;
  }

  return rc;
}

int
FileconvConfig::RemoveChiralCentresOnMatchedAtoms(
    Molecule& m, const resizable_array_p<Substructure_Query>& remove_chiral_centres_on) {
  int nchiral = m.chiral_centres();

  if (0 == nchiral)
    return 0;

  Molecule_to_Match target(&m);

  int nq = remove_chiral_centres_on.number_elements();

  Set_of_Atoms atoms_with_chiral_centres_to_be_removed;

  for (int i = 0; i < nq; i++) {
    Substructure_Results sresults;

    int nhits = remove_chiral_centres_on[i]->substructure_search(target, sresults);

    for (int j = 0; j < nhits; j++) {
      const Set_of_Atoms* e = sresults.embedding(j);

      IdentifyMatchedAtomsWithChiralCentres(m, *e, atoms_with_chiral_centres_to_be_removed);
    }
  }

  int rc = atoms_with_chiral_centres_to_be_removed.number_elements();

  if (0 == rc)
    return 0;

  molecules_with_chiral_centres_removed_by_rmchiral++;
  chiral_centres_removed_by_rmchiral += rc;

  for (int i = 0; i < rc; i++) {
    atom_number_t j = atoms_with_chiral_centres_to_be_removed[i];

    m.remove_chiral_centre_at_atom(j);
  }

  return rc;
}

int
FileconvConfig::RemoveNonOrganicChirality(Molecule& m) {
  int rc = 0;

  for (int i = m.chiral_centres() - 1; i >= 0; --i) {
    const Chiral_Centre* c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);

    const atom_number_t a = c->a();

    if (m.elementi(a)->organic())
      continue;

    m.remove_chiral_centre_at_atom(a);
    rc++;
  }

  if (rc)
    molecules_changed++;

  return rc;
}

int
FileconvConfig::AromatiseTheseRings(Molecule& m, int* procecss_these_atoms) {
  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i) {
    if (0 == procecss_these_atoms[i])
      continue;

    m.set_implicit_hydrogens_known(i, 0);

    const Atom* a = m.atomi(i);

    const auto acon = a->ncon();

    for (int j = 0; j < acon; ++j) {
      const Bond* b = a->item(j);

      if (b->is_single_bond())
        continue;

      atom_number_t k = a->other(i, j);
      if (procecss_these_atoms[k])
        m.set_bond_type_between_atoms(i, k, SINGLE_BOND);
    }

    if (m.formal_charge(i)) {
      m.set_formal_charge(i, 0);
      m.recompute_implicit_hydrogens(i);
    }
  }

  if (m.find_kekule_form(procecss_these_atoms))
    return 1;

  cerr << "Warning, did not find Kekule form for '" << m.name()
       << "', bonding has been destroyed!\n";
  return 0;
}

int
FileconvConfig::AromatiseTheseRings(Molecule& m,
                                     int* ring_labels,
                                     int* procecss_these_atoms,
                                     const Set_of_Atoms& e) {
  if (!m.is_ring_atom(e[0]))
    return 0;

  const auto matoms = m.natoms();

  const auto s = ring_labels[e[0]];

  for (auto i = 0; i < matoms; ++i) {
    if (s == ring_labels[i])
      procecss_these_atoms[i] = 1;
    else
      procecss_these_atoms[i] = 0;
  }

  return AromatiseTheseRings(m, procecss_these_atoms);
}

/*
  We need two arrays. One that holds the ring infromation for the molecule.
  The other is what gets passed to find_kekule_form.
  For efficiency, we just allocate one array and pass pieces of it to the other functions
*/

int
FileconvConfig::AromatiseTheseRings(Molecule& m,
                                     resizable_array_p<Substructure_Query>& aromatise_these_rings) {
  int rc = 0;

  Molecule_to_Match target(&m);

  const int matoms = m.natoms();

  int* tmp = new_int(matoms + matoms);
  std::unique_ptr<int[]> free_tmp(tmp);

  m.label_atoms_by_ring_system(tmp);

  for (int i = 0; i < aromatise_these_rings.number_elements(); ++i) {
    Substructure_Results sresults;

    int nhits = aromatise_these_rings[i]->substructure_search(target, sresults);

    for (auto j = 0; j < nhits; j++) {
      if (AromatiseTheseRings(m, tmp, tmp + matoms, *(sresults.embedding(i))))
        rc++;
    }
  }

  if (rc)
    molecules_changed_by_aromatising_rings++;

  return rc;
}

int
FileconvConfig::ExtractSmallerFragmentsIntoNameTag(Molecule& m,
                                                    const IWString& tag_for_removed_fragments) {
  const int matoms = m.natoms();

  Set_of_Atoms to_be_removed;
  int notused;

  (void)m.identify_largest_organic_fragment(to_be_removed, notused);

  int* subset = new_int(matoms);
  std::unique_ptr<int[]> free_subset(subset);

  to_be_removed.set_vector(subset, 1);

  Molecule f;
  m.create_subset(f, subset);

  m.remove_atoms(to_be_removed);

  IWString tmp(m.name());

  tmp << ' ' << tag_for_removed_fragments << ':' << f.smiles();

  m.set_name(tmp);

  return 1;
}

void
FileconvConfig::SortByFragmentSize(Molecule& m) {
  int nf = m.number_fragments();

  if (1 == nf)
    return;

  resizable_array<int> fragment_size;

  for (int i = 0; i < nf; i++) {
    int a = m.atoms_in_fragment(i);

    fragment_size.add_if_not_already_present(a);
  }

  int n = fragment_size.number_elements();

  if (1 == n)  // all frags the same size
    return;

  Int_Comparator_Smaller ics;

  fragment_size.iwqsort(ics);

#ifdef DEBUG_DO_SORT_BY_FRAGMENT_SIZE
  cerr << n << " different fragment sizes\n";
  for (int i = 0; i < n; i++) {
    cerr << fragment_size[i] << " atoms\n";
  }
#endif

  // Figure out the order in which existing atoms will be placed in the final molecule

  const int matoms = m.natoms();

  resizable_array<int> new_atom_order(matoms);

  for (int i = 0; i < n; i++) {
    const int a1 = fragment_size[i];

    for (int j = 0; j < matoms; j++) {
      int f = m.fragment_membership(j);
      if (a1 == m.atoms_in_fragment(f))
        new_atom_order.add(j);
    }
  }

#ifdef DEBUG_DO_SORT_BY_FRAGMENT_SIZE
  for (int i = 0; i < matoms; i++) {
    cerr << "i = " << i << " new_atom_order " << new_atom_order[i] << '\n';
  }
#endif

  assert(new_atom_order.number_elements() == matoms);

  // Each atom knows it's starting atom number

  int* tmp = new int[matoms];
  std::unique_ptr<int[]> free_tmp(tmp);

  for (int i = 0; i < matoms; i++) {
    tmp[i] = i;
    m.set_user_specified_atom_void_ptr(i, tmp + i);
  }

  // This is really hairy.... Make a swap, then keep swapping atoms
  // until the destination gets it's final atom in place

  for (int i = 0; i < matoms; i++) {
    int j = new_atom_order[i];

    atom_number_t zatom = i;

    while (1) {
      const int* u = reinterpret_cast<const int*>(m.user_specified_atom_void_ptr(zatom));

      //    cerr << "zatom " << zatom << " compare " << *u << " and " << j << '\n';
      if (*u == j)
        break;

      //    cerr << "i = " << i << " swap " << zatom << " and " << j << '\n';

      m.swap_atoms(zatom, j);

      for (int k = i + 1; k < matoms; k++) {
        if (new_atom_order[k] == *u) {
          zatom = k;
          break;
        }
      }

      if (zatom == j)
        break;
    }
  }

  return;
}

int
FileconvConfig::RemoveByFragment(Molecule& m, const resizable_array<int>& fragments_to_remove) {
  if (fragments_to_remove.empty())
    return 1;

  Set_of_Atoms atoms_to_remove;

  for (int i = 0; i < fragments_to_remove.number_elements(); i++) {
    m.add_atoms_in_fragment(atoms_to_remove, fragments_to_remove[i]);
  }

  return m.remove_atoms(atoms_to_remove);
}

int
FileconvConfig::RemoveFragmentsThatMatchQuery(Molecule& m,
                                               resizable_array_p<Substructure_Query>& queries) {
  int nq = queries.number_elements();

  Molecule_to_Match target(&m);

  resizable_array<int> fragments_to_be_removed;

  for (int i = 0; i < nq; i++) {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    nhits = sresults.number_embeddings();

    if (0 == nhits) {
      cerr << "Fragment selection query produced hits but no embeddings\n";
      cerr << "If compound query 'smarts1&&0smarts2', try reordering '0smarts2&&smarts1'\n";
    }

    for (int j = 0; j < nhits; j++) {
      const Set_of_Atoms* e = sresults.embedding(j);

      atom_number_t a0 = e->item(0);

      int f = m.fragment_membership(a0);

      fragments_to_be_removed.add_if_not_already_present(f);
    }
  }

  if (fragments_to_be_removed.empty())
    return 0;

  if (m.number_fragments() == fragments_to_be_removed.number_elements()) {
    if (verbose > 1)
      cerr << "All fragments match remove queries, ignored\n";
    return 0;
  }

  return RemoveByFragment(m, fragments_to_be_removed);
}

int
FileconvConfig::ReduceToAllLargestFragments(Molecule& m) {
  resizable_array_p<Molecule> frags;

  m.create_components(frags);

  int nf = frags.number_elements();

  int atoms_in_largest_frag = frags[0]->natoms();

  for (int i = 1; i < nf; i++) {
    int m = frags[i]->natoms();

    if (m > atoms_in_largest_frag)
      atoms_in_largest_frag = m;
  }

  resizable_array<int> fragments_to_remove;

  for (int i = 0; i < nf; i++) {
    int m = frags[i]->natoms();

    if (m != atoms_in_largest_frag)
      fragments_to_remove.add(i);
  }

  return RemoveByFragment(m, fragments_to_remove);
}

/*
  remove the remove_largest_fragment largest fragments
  Jun 2022, adapt to do either largest or smallest fragment.
  The only difference is the ordering of the sorted array of
  fragment sizes.
*/

int
FileconvConfig::RemoveLargestFragments(Molecule& m, int fragments_to_remove,
                        int larger_or_smaller) {
  int nf = m.number_fragments();

  if (fragments_to_remove >= nf) {
    return 0;
  }

  int* atoms_in_fragment = new_int(nf);
  std::unique_ptr<int[]> free_atoms_in_fragment(atoms_in_fragment);

  for (int i = 0; i < nf; i++) {
    atoms_in_fragment[i] = m.atoms_in_fragment(i);
  }

  // Sort the atoms in fragment array ascending or descending depending..
  // Also set a comparitor that will be used to compare individual fragment
  // sizes with the threshold for removal.
  std::function<bool(int, int)> size_compare;
  if (larger_or_smaller > 0) {
    Int_Comparator_Larger cmp;
    iwqsort(atoms_in_fragment, nf, cmp);
    size_compare = std::less<int>();

  } else {
    Int_Comparator_Smaller cmp;
    iwqsort(atoms_in_fragment, nf, cmp);
    size_compare = std::greater<int>();
  }

  int threshold_for_removal = atoms_in_fragment[nf - fragments_to_remove];

  // The first fragments_to_remove fragments will be removed

  int* fragment_being_removed = new_int(nf);
  std::unique_ptr<int[]> free_fragment_being_removed(fragment_being_removed);

  for (int i = 0; i < nf; i++) {
    int aif = m.atoms_in_fragment(i);

    if (size_compare(aif, threshold_for_removal)) {
      continue;
    }

    fragment_being_removed[i] = 1;
    fragments_to_remove--;
    if (0 == fragments_to_remove) {
      break;
    }
  }

  assert(0 == fragments_to_remove);

  int matoms = m.natoms();

  int* atoms_to_be_removed = new_int(matoms);
  std::unique_ptr<int[]> free_atoms_to_be_removed(atoms_to_be_removed);

  for (int i = 0; i < matoms; i++) {
    int f = m.fragment_membership(i);

    if (fragment_being_removed[f]) {
      atoms_to_be_removed[i] = 1;
    }
  }

  return m.remove_atoms(atoms_to_be_removed);
}

/*template void iwqsort<int, Int_Comparator_Larger>(int*, int, Int_Comparator_Larger&, void*);
template void iwqsort<int>(int*, int, Int_Comparator_Larger &);
template void swap_elements<int>(int&, int&, void*);
template void move_in_from_left<int, Int_Comparator_Larger>(int*, int&, int&, int,
Int_Comparator_Larger&, void*); template void move_in_from_right<int, Int_Comparator_Larger>(int*,
int&, int&, Int_Comparator_Larger&); template void compare_two_items<int,
Int_Comparator_Larger>(int*, Int_Comparator_Larger&, void*);*/

int
FileconvConfig::TooManyLargeFragments(Molecule& m,
                                       int min_size_max_fragment_count,
                                       int maximum_fragment_count) {
  int nf = m.number_fragments();

  int rc = 0;

  for (int i = 0; i < nf; i++) {
    int aif = m.atoms_in_fragment(i);

    if (aif >= min_size_max_fragment_count)
      rc++;
  }

  if (rc > maximum_fragment_count)
    return 1;

  return 0;
}

// Remove the largest or smallest fragment. Argument `less_or_greater` should
// be one of std::less<int> or std::greater<int> and that governs whether we
// are finding the largest or smallest fragment to remove.
int
FileconvConfig::RemoveLargestFragment(Molecule& m,
                                      std::function<bool(int, int)> less_or_greater) {
  int nf = m.number_fragments();

  int atoms_in_largest_fragment = m.atoms_in_fragment(0);
  int largest_fragment = 0;

  for (int i = 1; i < nf; i++) {
    int aif = m.atoms_in_fragment(i);

    if (less_or_greater(aif, atoms_in_largest_fragment)) {
      atoms_in_largest_fragment = aif;
      largest_fragment = i;
    }
  }

  int matoms = m.natoms();

  int* atoms_to_be_removed = new_int(matoms);
  std::unique_ptr<int[]> free_atoms_to_be_removed(atoms_to_be_removed);

  for (int i = 0; i < matoms; i++) {
    if (m.fragment_membership(i) == largest_fragment) {
      atoms_to_be_removed[i] = 1;
    }
  }

  return m.remove_atoms(atoms_to_be_removed);
}

void
FileconvConfig::RemoveNonorganicFragments(Molecule& m,
                                           const int* frag_membership,
                                           int& fragments_removed) {
  int nf = m.number_fragments();

  resizable_array<int> fragments_to_remove;

  for (int i = 0; i < nf; i++) {
    int aif = m.atoms_in_fragment(i);

    Molecule tmp;
    m.create_subset(tmp, frag_membership, i);

    bool delete_fragment = false;
    for (int j = 0; j < aif; j++) {
      const Element* e = tmp.elementi(j);
      if (e == nullptr) {
        cerr << "Null pointer found for element " << j << " in molecule " << m.name() << '\n';
        cerr << "This should not happen, please contact c3tk" << '\n';
        return;
      }
      if (e->organic()) {
        continue;
      }
      if (!ok_non_organics.matches(e)) {
        delete_fragment = true;
        break;
      }
      if (!e->is_in_periodic_table()) {
        delete_fragment = true;
        break;
      }
    }

    if (delete_fragment) {
      if (verbose)
        cerr << "Removing fragment " << i << " " << tmp.smiles() << '\n';

      fragments_to_remove.add(i);
      fragments_removed++;
    }
  }

  if (nf == fragments_to_remove.number_elements()) {
    m.resize(0);  // remove everything from the molecule
    return;
  }

  RemoveByFragment(m, fragments_to_remove);
  return;
}

void
FileconvConfig::RemoveNonorganicFragments(Molecule& m, int& fragments_removed) {
  fragments_removed = 0;
  int matoms = m.natoms();

  int* frag_membership = new int[matoms];
  std::unique_ptr<int[]> free_frag_membership(frag_membership);
  m.fragment_membership(frag_membership);

  RemoveNonorganicFragments(m, frag_membership, fragments_removed);

  if (verbose > 1)
    cerr << "new number of atoms: " << m.natoms() << '\n';
}

int
at_least_one_of_these_queries_matches(Molecule& m, resizable_array_p<Substructure_Query>& q) {
  Molecule_to_Match target(&m);

  int nq = q.number_elements();

  for (int i = 0; i < nq; i++) {
    if (q[i]->substructure_search(target))
      return 1;
  }

  return 0;
}

int
FileconvConfig::IdentifyFragmentsToBeKept(Molecule& m, resizable_array_p<Substructure_Query>& q) {
  resizable_array_p<Molecule> frags;
  m.create_components(frags);

  resizable_array<int> fragments_to_remove;

  int nf = frags.number_elements();
  for (int i = 0; i < nf; i++) {
    if (!at_least_one_of_these_queries_matches(*(frags[i]), q))
      fragments_to_remove.add(i);
  }

  if (nf == fragments_to_remove.number_elements()) {
    m.resize(0);  // remove everything from the molecule
    return 1;
  }

  return RemoveByFragment(m, fragments_to_remove);
}

int
IdentifyFragmentByQuery(resizable_array_p<Molecule>& components,
                        Substructure_Query* query,
                        int largest) {
  int atoms_in_previous_best_match = largest ? 0 : std::numeric_limits<int>::max();
  int best_match_component = -1;

  for (int i = 0; i < components.number_elements(); i++) {
    if (query->substructure_search(components[i])) {
      int nc = components[i]->natoms();
      if (largest && nc > atoms_in_previous_best_match) {
        atoms_in_previous_best_match = nc;
        best_match_component = i;
      } else if (!largest && nc < atoms_in_previous_best_match) {
        atoms_in_previous_best_match = nc;
        best_match_component = i;
      }
    }
  }

  return best_match_component;
}

int
FileconvConfig::IdentifyFragmentByQuery(Molecule& m,
                                         resizable_array_p<Substructure_Query>& queries,
                                         int largest) {
  resizable_array_p<Molecule> components;

  m.create_components(components);

  for (int i = 0; i < queries.number_elements(); i++) {
    int j = fileconv::IdentifyFragmentByQuery(components, queries[i], largest);
    if (j >= 0) {
      m.delete_all_fragments_except(j);
      if (verbose > 1)
        cerr << largest << " query matches fragment " << j << ", " << m.natoms() << " atoms\n";
      molecules_with_fragments_reduced_by_query++;
      return 1;
    }
  }

  molecules_not_matching_fragment_queries++;

  if (verbose)
    cerr << "None of " << components.number_elements() << " components matched queries '"
         << m.name() << "'\n";

  return 0;
}

int
FileconvConfig::AtomsInNonLargestFragmentExceed(Molecule& m, int mxfs) {
  int nf = m.number_fragments();

  int atoms_in_largest_fragment = m.atoms_in_fragment(0);

  // cerr << "Checking " << nf << " fragments\n";

  for (int i = 1; i < nf; i++) {
    int a = m.atoms_in_fragment(i);

    if (a > atoms_in_largest_fragment)
      atoms_in_largest_fragment = a;
  }

  for (int i = 0; i < nf; i++) {
    int a = m.atoms_in_fragment(i);

    if (a == atoms_in_largest_fragment)
      continue;

    if (a > mxfs)
      return 1;
  }

  return 0;
}

/*
  The molecule has too many fragments. Trim to the first FRAGMENT_COUNT
*/

int
FileconvConfig::TrimToFirstNFragments(Molecule& m,
                                       const int* frag_membership,
                                       int fragment_count) {
  int natoms = m.natoms();
  for (int i = 0, j = 0; j < natoms; j++) {
    if (frag_membership[j] + 1 > fragment_count)
      m.remove_atom(i);
    else
      i++;
  }

  return 1;
}

int
FileconvConfig::StripToNLargestFragments(Molecule& m, int n) {
  resizable_array<int> fragment_sizes;

  int nf = m.number_fragments();

  if (nf <= n)  // nothing to do
    return 1;

  for (int i = 0; i < nf; i++) {
    int a = m.atoms_in_fragment(i);
    fragment_sizes.add(a);
  }

  Int_Comparator_Smaller icl;

  fragment_sizes.iwqsort(icl);

#ifdef DEBUG_DO_STRIP_TO_N_LARGEST_FRAGMENTS
  for (int i = 0; i < fragment_sizes.number_elements(); i++) {
    cerr << "Before erasure " << fragment_sizes[i] << '\n';
  }
#endif

  fragment_sizes.erase(0, n - 1);

#ifdef DEBUG_DO_STRIP_TO_N_LARGEST_FRAGMENTS
  cerr << "Remaining fragment sizes " << fragment_sizes.number_elements() << '\n';
  for (int i = 0; i < fragment_sizes.number_elements(); i++) {
    cerr << "After erasure " << fragment_sizes[i] << '\n';
  }
#endif

  int* remove_fragment = new_int(nf);
  std::unique_ptr<int[]> free_remove_fragment(remove_fragment);

  for (int i = 0; i < nf; i++) {
    int a = m.atoms_in_fragment(i);

    if (!fragment_sizes.remove_first(a))
      continue;

    remove_fragment[i] = 1;
  }

  m.delete_fragments(remove_fragment);

  return 1;
}

// Remove all fragments that have more than, or less than `threshold` atoms.
// cmp must be something like std::less or std::greater.
template <typename C>
int
FileconvConfig::RemoveFragmentsByAtomCount(Molecule& m,
                        int threshold,
                        C cmp) {
  const int nf = m.number_fragments();

  std::unique_ptr<int[]> remove_fragment(new_int(nf));

  for (int i = 0; i < nf; i++) {
    if (cmp(m.atoms_in_fragment(i), threshold)) {
      remove_fragment[i] = 1;
    }
  }

  const int matoms = m.natoms();
  std::unique_ptr<int[]> atoms_to_remove(new_int(matoms));

  for (int i = 0; i < matoms; ++i) {
    if (remove_fragment[m.fragment_membership(i)]) {
      atoms_to_remove[i] = 1;
    }
  }

  return m.remove_atoms(atoms_to_remove.get());
}

/*
  Return 1 if the molecule should be kept, 0 if it should be discarded
*/

int
FileconvConfig::ProcessFragments(Molecule& m) {
  if (0 == maximum_fragment_count)
    ;
  else if (min_size_max_fragment_count) {
    if (TooManyLargeFragments(m, min_size_max_fragment_count, maximum_fragment_count)) {
      if (verbose > 1)
        cerr << "Too many large components " << m.number_fragments() << " max is "
             << maximum_fragment_count << '\n';
      molecules_with_too_many_components++;
      return 0;
    }
  } else if (m.number_fragments() > maximum_fragment_count) {
    if (verbose > 1)
      cerr << "Too many components " << m.number_fragments() << " max is " << maximum_fragment_count
           << '\n';
    molecules_with_too_many_components++;
    return 0;
  }

  if (remove_molecules_with_non_largest_fragment_natoms >= 0 &&
      AtomsInNonLargestFragmentExceed(m, remove_molecules_with_non_largest_fragment_natoms)) {
    if (verbose > 1) {
      cerr << "Contains large fragment(s)\n";
    }
    molecules_with_large_fragments++;
    return 0;
  }

  if (remove_duplicate_fragments) {
    int fragments_removed;

    remove_duplicate_fragments::RemoveDuplicateFragments(m, fragments_removed);

    if (fragments_removed) {
      if (verbose > 1) {
        cerr << fragments_removed << " duplicate fragments removed\n";
      }
      molecules_with_duplicate_fragments_removed++;
    }

    if (1 == m.number_fragments()) {
      return 1;
    }
  }

  if (remove_fragments_this_size_or_smaller) {
    RemoveFragmentsThisSizeOrSmaller(m);
    if (1 == m.number_fragments())
      return 1;
  }

  if (strip_to_n_largest_fragments > 0) {
    StripToNLargestFragments(m, strip_to_n_largest_fragments);
  }

  if (known_fragment_data.active()) {
    known_fragment_data.process(m);
    if (1 == m.number_fragments())
      return 1;
  }

  // If we have some fragments to definitely remove, process them

  if (remove_fragment_queries.number_elements()) {
    RemoveFragmentsThatMatchQuery(m, remove_fragment_queries);
    if (m.number_fragments() == 1) {
      return 1;
    }
  }

  if (reduce_to_largest_fragment) {  // the most common case
    if (verbose <= 1)
      return m.reduce_to_largest_fragment();

    int initial_matoms = m.natoms();

    int rc = m.reduce_to_largest_fragment();
    if (verbose > 1)
      cerr << "Stripped to largest fragment, lost " << (initial_matoms - m.natoms()) << " atoms\n";

    return rc;
  }

  if (reduce_to_all_largest_fragments) {
    return ReduceToAllLargestFragments(m);
  }

  if (remove_largest_fragment > 0 || remove_smallest_fragment > 0) {
    if (1 == remove_largest_fragment) {
      RemoveLargestFragment(m, std::greater<int>());
      if (m.number_fragments() == 1) {
        return 1;
      }
    }

    if (remove_largest_fragment > 1) {
      RemoveLargestFragments(m, remove_largest_fragment, 1);
      if (m.number_fragments() == 1) {
        return 1;
      }
    }

    if (remove_smallest_fragment == 1) {
      return RemoveLargestFragment(m, std::less<int>());
    }

    if (remove_smallest_fragment > 1) {
      RemoveLargestFragments(m, remove_smallest_fragment, -1);
    }

    return 1;
  }

  if (keep_smallest_fragments > 0) {
    const int nf = m.number_fragments();
    if (nf <= keep_smallest_fragments) {
      return 1;
    }
    // If we are keeping all but the largest fragment, remove it.
    if (keep_smallest_fragments + 1 == nf) {
      return RemoveLargestFragment(m, std::greater<int>());
    }
    return RemoveLargestFragments(m, nf - keep_smallest_fragments, 1);
  }

  if (reduce_to_largest_organic_fragment) {
    if (verbose > 1)
      cerr << "Stripped to largest organic fragment\n";

    if (tag_for_removed_fragments.length())
      return ExtractSmallerFragmentsIntoNameTag(m, tag_for_removed_fragments);
    else
      return m.reduce_to_largest_organic_fragment();
  }

  if (reduce_to_largest_organic_fragment_carefully) {
    if (verbose > 1)
      cerr << "Stripped to largest organic fragment with desirable features\n";

    return m.reduce_to_largest_fragment_carefully();
  }

  if (keep_all_organic_fragments) {
    int fragments_removed;

    RemoveNonorganicFragments(m, fragments_removed);

    if (fragments_removed > 0) {
      if (verbose) {
        cerr << fragments_removed << " fragments removed for " << m.name() << '\n';
      }
      molecules_with_nonorganic_fragments_removed++;
      if (m.number_fragments() < 1) {
        return 0;
      }
    }

    return 1;
  }

  if (smallest_fragment_queries.number_elements())
    return IdentifyFragmentByQuery(m, smallest_fragment_queries, 0);

  if (largest_fragment_queries.number_elements())
    return IdentifyFragmentByQuery(m, largest_fragment_queries, 1);

  if (keep_fragment_queries.number_elements()) {
    return IdentifyFragmentsToBeKept(m, keep_fragment_queries);
  }

  if (fragment_selector_lower_atom_count > 0) {
    RemoveFragmentsByAtomCount(m, fragment_selector_lower_atom_count, std::less<int>());
  }

  if (fragment_selector_upper_atom_count < std::numeric_limits<int>::max()) {
    RemoveFragmentsByAtomCount(m, fragment_selector_upper_atom_count, std::greater<int>());
  }

  if (sort_by_fragment_size) {
    SortByFragmentSize(m);
  }

  if (0 == fragment_count)
    return 1;

  int components = m.number_fragments();

  int* frag_membership = new int[m.natoms()];
  std::unique_ptr<int[]> free_frag_membership(frag_membership);
  m.fragment_membership(frag_membership);

  int rc = TrimToFirstNFragments(m, frag_membership, fragment_count);

  if (verbose > 1)
    cerr << "reduced from " << components << " fragments\n";

  return rc;
}

int
MaxRingSize(Molecule& m) {
  int nr = m.nrings();

  if (0 == nr)
    return 0;

  int rc = 0;

  for (int i = nr - 1; i >= 0; i--) {
    const Ring* r = m.ringi(i);

    if (r->number_elements() > rc)
      rc = r->number_elements();
  }

  return rc;
}

/*
  Some tests need a check of each atom.
  For efficiency we have a function which does both ok_non_organics and
  the non periodic table atoms
*/

int
FileconvConfig::ExcludeForNonOrganicAndNonPeriodicTable(const Molecule& m) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; i++) {
    const Element* e = m.elementi(i);
    if (e->organic())
      continue;

    if (!ok_non_organics.matches(e)) {
      non_organic_molecules_found++;
      return 1;
    }

    if (!e->is_in_periodic_table()) {
      non_real_elements_found++;
      return 1;
    }
  }

  return 0;
}

int
FileconvConfig::ExcludeForNonOrganic(const Molecule& m) {
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++) {
    const Element* e = m.elementi(i);
    if (e->organic())
      continue;

    if (!ok_non_organics.matches(e)) {
      non_organic_molecules_found++;
      return 1;
    }
  }

  return 0;
}

int
FileconvConfig::ExcludeForNonRealElements(const Molecule& m) {
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++) {
    const Element* e = m.elementi(i);

    if (!e->is_in_periodic_table()) {
      non_real_elements_found++;
      return 1;
    }
  }

  return 0;
}

/*
  for (hopefull) better efficiency, exclude_for_non_organic is written
  to do both organics and non_periodic table elements
*/

int
FileconvConfig::ExcludeForAtomTypes(const Molecule& m) {
  if (exclude_non_real_elements && output_organic_only)
    return ExcludeForNonOrganic(m);

  if (output_organic_only)
    return ExcludeForNonOrganic(m);

  return ExcludeForNonRealElements(m);
}

int
FileconvConfig::ValenceCheckOk(Molecule& m) {
  if (m.valence_ok())
    return 1;

  if (!ok_bad_valence_on_isotopically_labelled)
    return 0;

  const int matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i) {
    if (m.valence_ok(i))
      continue;

    if (0 == m.isotope(i))  // bad valence, but isotope, therefore OK
      return 0;
  }

  return 1;
}

int
FileconvConfig::LooksLikeMixtureByAtomCount(Molecule& m) {
  int nf = m.number_fragments();

  int atoms_in_largest_fragment = m.atoms_in_fragment(0);
  int atoms_in_second_largest_fragment = m.atoms_in_fragment(1);

  if (atoms_in_largest_fragment < atoms_in_second_largest_fragment)
    std::swap(atoms_in_largest_fragment, atoms_in_second_largest_fragment);

  // cerr << "First two counts " << atoms_in_largest_fragment << " and " <<
  // atoms_in_second_largest_fragment << '\n';

  for (int i = 2; i < nf; i++) {
    int a = m.atoms_in_fragment(i);

    if (a > atoms_in_largest_fragment) {
      atoms_in_second_largest_fragment = atoms_in_largest_fragment;
      atoms_in_largest_fragment = a;
    } else if (a > atoms_in_second_largest_fragment)
      atoms_in_second_largest_fragment = a;
  }

  if (atoms_in_largest_fragment - atoms_in_second_largest_fragment <=
      mixture_if_largest_frags_differ_by)
    return 1;

  return 0;
}

int
FileconvConfig::TooManyLargeFragments(Molecule& m) {
  int rc = 0;

  int nf = m.number_fragments();

  for (int i = 0; i < nf; i++) {
    if (m.atoms_in_fragment(i) > discard_molecule_if_multiple_fragments_larger_than) {
      rc++;

      if (rc > 1)
        return 1;
    }
  }

  return 0;
}

int
FileconvConfig::CountRingSystems(Molecule& m) {
  int* notused = new int[m.natoms()];
  std::unique_ptr<int[]> free_notused(notused);

  if (ring_systems_include_spiro) {
    return m.label_atoms_by_ring_system_including_spiro_fused(notused);
  } else {
    return m.label_atoms_by_ring_system(notused);
  }
}

int
LargestNumberRingsInASystem(Molecule& m) {
  int nr = m.nrings();

  if (0 == nr)
    return 0;

  if (1 == nr)
    return 1;

  if (2 == nr)
    return m.ringi(0)->is_fused() ? 2 : 1;

  // Now things get more complicated

  int* ring_already_done = new_int(nr);
  std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  int rc = 1;

  for (int i = 0; i < nr; i++) {
    if (ring_already_done[i])
      continue;

    const Ring* ri = m.ringi(i);

    if (!ri->is_fused())
      continue;

    int system_size = 1;

    for (int j = i + 1; j < nr; j++) {
      if (ring_already_done[j])
        continue;

      const Ring* rj = m.ringi(j);

      if (rj->fused_system_identifier() != ri->fused_system_identifier())
        continue;

      system_size++;

      ring_already_done[j] = 1;
    }

    if (system_size > rc)
      rc = system_size;
  }

  return rc;
}

int
AtomsInLargestFragment(Molecule& m) {
  int nf = m.number_fragments();

  if (1 == nf)
    return m.natoms();

  int rc = m.atoms_in_fragment(0);

  for (int i = 1; i < nf; i++) {
    int tmp = m.atoms_in_fragment(i);

    if (tmp > rc)
      rc = tmp;
  }

  return rc;
}

// Remove any atoms that match the `remove_atom` queries from `m`.
int
FileconvConfig::RemoveAtoms(Molecule& m) {
  const int matoms = m.natoms();
  std::unique_ptr<int[]> atoms_to_be_removed(new_int(matoms));

  Molecule_to_Match target(&m);
  Substructure_Results sresults;
  int matches_found = 0;
  for (Substructure_Query* q : remove_atom) {
    int nhits = q->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }

    ++matches_found;
    sresults.each_embedding_set_vector(atoms_to_be_removed.get(), 1);
  }

  if (matches_found == 0) {
    return 0;
  }

  m.remove_atoms(atoms_to_be_removed.get());
  ++molecules_with_removed_atoms;
  return 1;
}

// Return true if there is a Boron atom that is NOT attached to at
// least two oxygen atoms
int
ContainsDisallowedBoron(const Molecule& m) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m[i];
    if (a.atomic_number() != 5) {
      continue;
    }

    if (a.ncon() != 3) {  // Probably correct.
      return 1;
    }

    int attached_oxygen_count = 0;
    for (const Bond* b : a) {
      if (! b->is_single_bond()) {
        return 1;
      }

      const atom_number_t o = b->other(i);
      if (m.atomic_number(o) == 8) {
        ++attached_oxygen_count;
      }
    }
    if (attached_oxygen_count < 2) {
      return 1;
    }
  }

  return 0;
}

int
FileconvConfig::RejectedForDisallowedAtoms(Molecule& m) {
  if (ok_boron_special_case && ContainsDisallowedBoron(m)) {
    return 1;
  }

  return allowed_elements.contains_non_allowed_atoms(m);
}

int
FileconvConfig::ApplyAllFiltersInner(Molecule& m,
                                      IWString& rejection_reason,
                                      int& molecule_changed) {
  assert(m.ok());

  int matoms = m.natoms();

  int nf = m.number_fragments();

  initial_fragment_count[nf]++;

  if (audit_input)
    ;
  else if (!filter_for_disallowed_elements)
    ;
  else if (RejectedForDisallowedAtoms(m)) {
    molecules_excluded_for_non_allowed_elements++;
    rejection_reason << "NonAllowedElement";
    return 0;
  }

  // Apr 2013. User wants largest/smallest fragment that contains a query. But if the molecule
  // has just one fragment, we still want that to work.

  if (1 == nf && smallest_fragment_queries.number_elements() &&
      !at_least_one_of_these_queries_matches(m, smallest_fragment_queries)) {
    molecules_not_matching_fragment_queries++;
    rejection_reason << "NoMatchSmlNeededFrag";
    return 0;
  }

  if (1 == nf && largest_fragment_queries.number_elements() &&
      !at_least_one_of_these_queries_matches(m, largest_fragment_queries)) {
    molecules_not_matching_fragment_queries++;
    rejection_reason << "NoMatchLrgNeededFrag";
    return 0;
  }

  if ((nf < 2) && (!keep_all_organic_fragments))
    // cannot reduce the fragment count, except that we may want to remove the
    // only fragment when it is not organic
    ;
  else if (discard_molecule_if_multiple_fragments_larger_than && TooManyLargeFragments(m)) {
    molecules_with_multiple_large_fragments++;
    rejection_reason << "Mixture";
    return 0;
  } else if (mixture_if_largest_frags_differ_by >= 0 && LooksLikeMixtureByAtomCount(m)) {
    molecules_declared_mixtures_by_atom_count_difference++;
    rejection_reason = "Mixture";
    return 0;
  } else if (need_to_call_process_fragments) {
    int initial_matoms = m.natoms();

    if (!ProcessFragments(m)) {
      if (verbose > 1)
        cerr << "Molecule does not match fragment criteria\n";

      rejection_reason = "no match to fragment criteria";
      return 0;
    }

    matoms = m.natoms();

    if (initial_matoms != matoms) {
      molecule_changed++;
      molecules_chopped++;
    }

    atoms_lost.extra(initial_matoms - matoms);
  }

  if (aromatise_these_rings.number_elements()) {
    AromatiseTheseRings(m, aromatise_these_rings);
  }

  if (matoms > 0) {
    molecule_changed += elements_to_remove.process(m);

    molecule_changed += element_transformations.process(m);

    if (structure_fixing.active())
      molecule_changed += structure_fixing.process(m);

    matoms = m.natoms();  // in case the atom count changed
  }

  int implicit_hydrogens;

  if (include_implicit_hydrogens_in_upper_atom_count_comparison)
    implicit_hydrogens = m.implicit_hydrogens();
  else
    implicit_hydrogens = 0;

  int keep = 0;

  if (atom_count_includes_only_atoms_in_largest_fragment)
    matoms = AtomsInLargestFragment(m);

  if (lower_atom_count_cutoff > 0 && matoms < lower_atom_count_cutoff) {
    molecules_below_atom_count_cutoff++;
    if (verbose > 1)
      cerr << "Molecule contains " << m.natoms() << " atoms, which is below cutoff\n";

    rejection_reason = "too few atoms";
  } else if (upper_atom_count_cutoff > 0 &&
             (matoms + implicit_hydrogens) > upper_atom_count_cutoff) {
    molecules_above_atom_count_cutoff++;
    if (verbose > 1)
      cerr << "Molecule contains " << m.natoms() << " atoms, which is above cutoff\n";
    rejection_reason = "too many atoms";
  } else if (lower_molecular_weight_cutoff > 0.0 &&
             m.molecular_weight_ignore_isotopes() < lower_molecular_weight_cutoff) {
    molecules_below_molecular_weight_cutoff++;
    cerr << "Molecular weight is " << m.molecular_weight_ignore_isotopes()
         << " which is below cutoff\n";
    rejection_reason = "AMW too low";
  } else if (upper_molecular_weight_cutoff > 0.0 &&
             m.molecular_weight_ignore_isotopes() > upper_molecular_weight_cutoff) {
    molecules_above_molecular_weight_cutoff++;
    cerr << "Molecular weight is " << m.molecular_weight_ignore_isotopes()
         << " which is above cutoff\n";
    rejection_reason = "AMW too high";
  } else if ((output_organic_only || exclude_non_real_elements) && ExcludeForAtomTypes(m)) {
    if (verbose > 1)
      cerr << "Non organic or non periodic table atoms found\n";
    rejection_reason = "non organic";
  } else if (lower_ring_count_cutoff && m.nrings() < lower_ring_count_cutoff) {
    molecules_with_too_few_rings++;
    if (verbose > 1)
      cerr << "Molecule contains " << m.nrings() << " rings, which is below cutoff\n";
    rejection_reason = "too few rings";
  } else if (upper_ring_count_cutoff >= 0 && m.nrings() > upper_ring_count_cutoff) {
    molecules_with_too_many_rings++;
    if (verbose > 1)
      cerr << "Molecule contains " << m.nrings() << " rings, which is above cutoff\n";
    rejection_reason = "too many rings";
  } else if (max_ring_systems_allowed >= 0 && CountRingSystems(m) > max_ring_systems_allowed) {
    molecules_with_too_many_ring_systems++;
    if (verbose > 1)
      cerr << "Molecule contains too many ring systems\n";
    rejection_reason << "too many ring systems";
  } else if (min_ring_systems_required >= 0 && CountRingSystems(m) < min_ring_systems_required) {
    molecules_with_too_few_ring_systems++;
    if (verbose > 1)
      cerr << "Molecule contains too few ring systems\n";
    rejection_reason << "too few ring systems";
  } else if (upper_ring_size_cutoff > 0 && MaxRingSize(m) > upper_ring_size_cutoff) {
    molecules_with_large_rings++;
    if (verbose > 1)
      cerr << "Molecule contains ring too large\n";
    rejection_reason = "large ring";
  } else if (max_rings_in_a_ring_system > 0 &&
             LargestNumberRingsInASystem(m) > max_rings_in_a_ring_system) {
    molecules_with_ring_systems_too_large++;
    if (verbose > 1)
      cerr << "Molecule has ring system with too many rings\n";
    rejection_reason = "large ring system";
  } else if (min_aliphatic_ring_count > 0 && 
             (m.nrings() - m.aromatic_ring_count()) < min_aliphatic_ring_count) {
    ++molecules_with_too_few_aliphatic_rings;
    if (verbose > 1) {
      cerr << "Too few aliphatic rings\n";
    }
    rejection_reason << "too few aliphatic rings";
  } else if (min_aromatic_ring_count > 0 && 
             m.aromatic_ring_count() < min_aromatic_ring_count) {
    ++molecules_with_too_few_aromatic_rings;
    if (verbose > 1) {
      cerr << "Too few aliphatic rings\n";
    }
    rejection_reason << "too few aromatic rings";
  } else if (max_aliphatic_ring_count >= 0 && 
             (m.nrings() - m.aromatic_ring_count()) > max_aliphatic_ring_count) {
    ++molecules_with_too_many_aliphatic_rings;
    if (verbose > 1) {
      cerr << "Too many aliphatic rings\n";
    }
    rejection_reason << "too many aliphatic rings";
  } else if (max_aromatic_ring_count >= 0 && 
             m.aromatic_ring_count() > max_aromatic_ring_count) {
    ++molecules_with_too_many_aromatic_rings;
    if (verbose > 1) {
      cerr << "Too many aromatic rings " << m.aromatic_ring_count() << '\n';
    }
    rejection_reason << "too many aromatic rings";
  } else if (exclude_isotopes && m.number_isotopic_atoms() > 0) {
    molecules_containing_isotopes++;
    if (verbose > 1)
      cerr << "Molecule contains isotopes\n";
    rejection_reason = "isotopes";
  } else if (skip_molecules_with_abnormal_valences && !ValenceCheckOk(m)) {
    molecules_with_abnormal_valences++;
    if (verbose > 1)
      cerr << "Molecule contains abnormal valence(s)\n";
    rejection_reason = "bad valence";
  } else if (max_path_length > 0 && m.longest_path() > max_path_length) {
    molecules_with_longest_path_too_long++;
    if (verbose > 1)
      cerr << "Longest path too long\n";
    rejection_reason << "long path";
  } else
    keep = 1;  // molecule is good!

  if (!keep)
    return 0;

  // It is an open question whether these atom removal functions
  // should be before the fragment selection code.
  if (!remove_atom.empty()) {
    molecule_changed += RemoveAtoms(m);
  }

  if (remove_isotopic_atoms) {
    molecule_changed += RemoveIsotopicAtoms(m);
  }

  if (isotope_is_atom_type) {
    molecule_changed += AtomTypeToIsotope(m);
  }

  // After we have done all other manipulations with isotopes

  if (convert_isotopes_to_atom_map_numbers)
    ConvertIsotopesToAtomMapNumbers(m);
  else if (convert_atom_map_numbers_to_isotopes)
    ConvertAtomMapNumbersToIsotopes(m);

  matoms = m.natoms();  // recompute, as it may have been changed by transformations

  if (compute_molecular_weight_for_each || lower_amw_cutoff >= 0.0 || upper_amw_cutoff >= 0.0) {
    molecular_weight_t amw = m.molecular_weight();

    amw_accumulator.extra(amw);
    atom_count[matoms]++;

    if (0 == verbose && compute_molecular_weight_for_each)
      cerr << m.name() << ' ';
    if (verbose > 1 || compute_molecular_weight_for_each)
      cerr << "AMW " << amw;

    if (lower_amw_cutoff >= 0.0 && amw < lower_amw_cutoff) {
      if (verbose > 1)
        cerr << " below cutoff\n";
      else if (compute_molecular_weight_for_each)
        cerr << '\n';
      molecules_below_amw_cutoff++;
      return 0;
    }

    if (upper_amw_cutoff > 0.0 && amw > upper_amw_cutoff) {
      if (verbose > 1)
        cerr << " above cutoff\n";
      else if (compute_molecular_weight_for_each)
        cerr << '\n';
      molecules_above_amw_cutoff++;
      return 0;
    }

    if (verbose > 1 || compute_molecular_weight_for_each)
      cerr << '\n';
  }

  if (verbose)
    chiral_centre_count[m.chiral_centres()]++;

  if (max_chiral_centres > 0 && m.chiral_centres() > max_chiral_centres) {
    if (verbose > 1)
      cerr << " too many chiral centres " << max_chiral_centres << '\n';
    molecules_with_too_many_chiral_centres++;
    rejection_reason << "too many chiral centres " << m.chiral_centres();
    return 0;
  }

  return 1;
}

int
do_convert_specific_isotopes_query(Molecule& m,
                                   Molecule_to_Match& target,
                                   Substructure_Query& q,
                                   const int new_isotope) {
  Substructure_Results sresults;

  const int nhits = q.substructure_search(target, sresults);

  if (0 == nhits) {
    return 0;
  }

  int rc = 0;

  for (int i = 0; i < nhits; ++i) {
    const Set_of_Atoms* e = sresults.embedding(i);

    const int n = e->number_elements();

    for (int j = 0; j < n; ++j) {
      const atom_number_t k = e->item(j);

      const int iso = m.isotope(k);

      if (iso == new_isotope)
        continue;

      m.set_isotope(k, new_isotope);
      rc++;
    }
  }

  return rc;
}

int
do_convert_specific_isotopes_query(
    Molecule& m,
    const resizable_array_p<Substructure_Query>& convert_specific_isotopes_query,
    const resizable_array<int>& convert_specific_isotopes_query_new_isotope) {
  assert(convert_specific_isotopes_query.number_elements() ==
         convert_specific_isotopes_query_new_isotope.number_elements());
  int rc = 0;

  Molecule_to_Match target(&m);

  for (int i = 0; i < convert_specific_isotopes_query.number_elements(); ++i) {
    if (do_convert_specific_isotopes_query(m,
                                           target,
                                           *convert_specific_isotopes_query[i],
                                           convert_specific_isotopes_query_new_isotope[i]))
      rc++;
  }

  return rc;
}

int
do_convert_all_isotopes_to(Molecule& m, const int convert_all_isotopes_to)

{
  const int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    if (0 == m.isotope(i))
      continue;

    m.set_isotope(i, convert_all_isotopes_to);
    rc++;
  }

  return rc;
}

int
do_convert_specific_isotopes(Molecule& m,
                             const resizable_array<int>& convert_specific_isotopes,
                             const resizable_array<int>& convert_specific_isotopes_new_isotope) {
  const int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; ++i) {
    const int iso = m.isotope(i);

    if (0 == iso)
      continue;

    const int f = convert_specific_isotopes.index(iso);
    //  cerr << "Do we need to convert isotope " << iso << " f = " << f << '\n';
    if (f < 0)
      continue;

    m.set_isotope(i, convert_specific_isotopes_new_isotope[f]);
    rc++;
  }

  return rc;
}

int
contains_covalently_bonded_non_organic(const Molecule& m) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    const Atom* a = m.atomi(i);

    if (a->element()->organic())
      continue;

    if (1 != a->ncon())  // non organic bonded to something
      return 1;
  }

  return 0;
}

/*
  How many of the leading characters are whitespace
*/

int
count_leading_whitespace(const char* s, int lens) {
  for (int i = 0; i < lens; i++) {
    if (!isspace(s[i]))
      return i;
  }

  return lens;
}

int
do_truncate_name(Molecule& m, char truncate_name_at_first, char truncate_name_at_last) {
  IWString newname = m.name();

  if ('\0' != truncate_name_at_first) {
    int i = newname.index(truncate_name_at_first);
    if (i >= 0)
      newname.iwtruncate(i);
  }

  if ('\0' != truncate_name_at_last) {
    int i = newname.index(truncate_name_at_last);
    if (i >= 0)
      newname.iwtruncate(i);
  }

  if (0 == newname.length()) {
    cerr << "do_truncate_name:cannot empty name '" << m.name() << "'\n";
    return 1;
  }

  if (newname.length() == m.name().length())  // no change
    return 1;

  m.set_name(newname);

  return 1;
}

void
do_truncate_names_to_first_token(Molecule& m) {
  const IWString& mname = m.name();

  int n = mname.length();

  if (0 == n)  // molecule with no name
    return;

  const char* s = mname.rawchars();

  int lws = count_leading_whitespace(s, n);

  if (lws) {
    s += lws;
    n -= lws;

    if (0 == n)  // name was all whitespace
    {
      m.set_name("");
      return;
    }
  }

  for (int i = 0; i < n; i++) {
    if (!isspace(s[i]))
      continue;

    IWString tmp(s, i);
    m.set_name(tmp);
    return;
  }

  return;
}

// Not sure why we make allowance for multiple capture groups
// since we only use the first one.
int
FileconvConfig::ChangeName(Molecule& m, re2::RE2& change_name_rx) {
  const IWString& mname = m.name();

  const int ncaptures = change_name_rx.NumberOfCapturingGroups();
  assert(ncaptures > 0);

  std::vector<std::string> captures(ncaptures);

  RE2::Arg* args = new RE2::Arg[ncaptures];
  for (int i = 0; i < ncaptures; ++i) {
    args[i] = &captures[i];
  }

  re2::StringPiece tmp(mname.data(), mname.length());
  if (RE2::PartialMatchN(tmp, change_name_rx, &args, ncaptures))
    ;
  else if (change_name_rx_must_match) {
    cerr << "change_name:change name rx not matched '" << change_name_rx.pattern()
         << "', molecule name '" << m.name() << "'\n";
    return 0;
  } else  // ok for no matches
    return 1;

  const std::string x = captures[0];
  IWString newname(x);

  if (0 == newname.length()) {
    cerr << "change_name:cannot set name to length zero, rx '" << change_name_rx.pattern()
         << "', mname '" << m.name() << "'\n";
    return 0;
  }

  m.set_name(newname);

  return 1;
}

void
do_extract_token_as_name(Molecule& m, int name_token) {
  IWString token;
  if (!m.name().word(name_token, token)) {
    cerr << "do_extract_token_as_name:cannot extract token " << name_token << " from '" << m.name()
         << "'\n";
    return;
  }
  m.set_name(token);
}

void
FileconvConfig::SubstituteForWhitespaceInName(Molecule& m) {
  assert(1 == substitute_for_whitespace_in_name.length());

  IWString tmp(m.name());

  int n = tmp.length();

  int need_to_reset_name = 0;

  for (int i = 0; i < n; i++) {
    if (!isspace(tmp[i]))
      continue;

    tmp[i] = substitute_for_whitespace_in_name[0];
    need_to_reset_name = 1;
  }

  if (need_to_reset_name)
    m.set_name(tmp);

  return;
}

int
FileconvConfig::ZeroPadName(Molecule& m) {
  constexpr char kZero = '0';

  IWString mname = m.name();

  if (zero_pad_name < 0) {
    if (mname.length() < -zero_pad_name) {
      cerr << "FileconvConfig::ZeroPadName:cannot remove " << zero_pad_name <<
              " leading chars from '" << mname << '\n';
      return 0;
    }
    mname.remove_leading_chars(-zero_pad_name);
    m.set_name(mname);
    return 1;
  }

  // If already that wide, cannot do anything.
  if (mname.length() >= zero_pad_name) {
    return 0;
  }

  mname.shift(zero_pad_name - mname.length(), kZero);
  m.set_name(mname);
  return 1;
}

int
do_remove_unnecessary_square_brackets(Molecule& m) {
  return m.unset_unnecessary_implicit_hydrogens_known_values();
}

int
do_remove_all_possible_square_brackets(Molecule& m) {
  return m.remove_hydrogens_known_flag_to_fix_valence_errors();
}

// Returns false if `m` contains isotopic atoms or
// covalently bonded non organics.
bool
IsOrganic(const Molecule& m) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m.atom(i);
    if (a.isotope() > 0) {  // Isotope.
      return false;
    }
    if (a.element()->organic()) {
      continue;
    }
    if (a.ncon() > 0) {  // Covalently bonded non organic
      return false;
    }
  }

  return true;
}

void
FileconvConfig::DoTranslation(Molecule& m) const {
  if (translate_to_origin) {
    Coordinates centroid;

    static constexpr int kFrag = 0;
    m.centroid(centroid, kFrag);
    m.translate_atoms(-centroid);
    return;
  }

  m.translate_atoms(dx, dy, dz);
}

int
FileconvConfig::NameMatchesNameRx(const IWString& mname) {
  if (grep_v_name_rx) {
    if (iwre2::RE2PartialMatch(mname, *grep_v_name_rx)) {
      return 0;
    }
  }

  if (!name_rx)  // not active, always matches.
    return 1;

  return iwre2::RE2PartialMatch(mname, *name_rx);
}

// Unset isotopes, while clearing any residual implicit hydrogens
// known information.
int
ConvertIsotopes(Molecule& m) {
  int rc = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.isotope(i) == 0) {
      continue;
    }
    m.set_isotope(i, 0);
    m.unset_all_implicit_hydrogen_information(i);
    ++rc;
  }
  return rc;
}
// Return 1 if any change is made.
int
FileconvConfig::AddFragmentToIsotopicAtoms(Molecule& m) {
  int rc = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const isotope_t iso = m.isotope(i);
    if (iso == 0) {
      continue;
    }

    auto iter = add_to_isotopic_atom.find(iso);
    if (iter == add_to_isotopic_atom.end()) {
      continue;
    }

    // Interesting design decision. What to do if someone has entered
    // something like [C] - Carbon with no Hydrogens. Let's be forgiving
    // and see if we can liberate some Hydrogens.
    m.unset_all_implicit_hydrogen_information(i);

    if (m.hcount(i) == 0) {
      continue;
    }

    const int initial_matoms = m.natoms();
    m.add_molecule(&iter->second);
    m.add_bond(i, initial_matoms, SINGLE_BOND);
    rc = 1;
  }

  return rc;
}

// Returns 1 if we change `m` and increment molecules_containing_isotopes
int
FileconvConfig::IsotopeRelated(Molecule& m) {
  int rc = 0;

  if (convert_isotopes == 1) {  // put here so chemical standardisation can get rid of converted Deuterium
    if (m.transform_to_non_isotopic_form()) {
      ++rc;
      molecules_containing_isotopes++;
    }
  } else if (convert_isotopes == 2) {
    if (ConvertIsotopes(m)) {
      ++rc;
      molecules_containing_isotopes++;
    }
  } else if (convert_all_isotopes_to) {
    if (do_convert_all_isotopes_to(m, convert_all_isotopes_to)) {
      ++rc;
      molecules_containing_isotopes++;
    }
  } else {
    if (convert_specific_isotopes.number_elements()) {
      if (do_convert_specific_isotopes(m,
                                       convert_specific_isotopes,
                                       convert_specific_isotopes_new_isotope)) {
        ++rc;
        molecules_containing_isotopes++;
      }
    }

    if (convert_specific_isotopes_query.number_elements()) {
      if (do_convert_specific_isotopes_query(m,
                                             convert_specific_isotopes_query,
                                             convert_specific_isotopes_query_new_isotope)) {
        ++rc;
        molecules_containing_isotopes++;
      }
    }

    if (add_to_isotopic_atom.size()) {
      rc += AddFragmentToIsotopicAtoms(m);
    }
  }

  return rc;
}
FileconvResult
FileconvConfig::Process(Molecule& m) {
  ++molecules_processed;

  FileconvResult result;

  if (!NameMatchesNameRx(m.name())) {
    result.rejected = 1;
    result.rejection_reason = "name_mismatch";
    ++molecules_discarded_for_name_mismatch;
    return result;
  }

  // if (verbose > 1)
  //   cerr << "Molecule " << input.molecules_read() << " finishes at line " << input.lines_read()
  //   << '\n';

  if (print_bond_lengths)
    (void)PrintBondLengths(m, std::cout);

  if (print_bond_angles)
    (void)PrintBondAngles(m, std::cout);

  if (print_torsions)
    (void)PrintTorsions(m, std::cout);

  if (print_max_atom_separation)
    (void)PrintMaxAtomSeparation(m, std::cout);

  if (verbose) {
    int matoms = m.natoms();
    natoms_accumulator(matoms);
    if (verbose > 1)
      atom_count[matoms]++;
  }

  if (audit_input) {
    if (appends_to_be_done) {
      DoAppends(m);
      cerr << molecules_processed << ' ' << m.name() << '\n';
    }

    return result;
  }

  int& molecule_changed = result.molecule_changed;  // Save typing
  molecule_changed = 0;

  if (remove_directional_bonds_from_input) {
    if (m.revert_all_directional_bonds_to_non_directional()) {
      molecule_changed++;
    }
  }

  if (need_to_consider_isotopes) {
    molecule_changed += IsotopeRelated(m);
  }

  if (reset_atom_map_numbers)
    m.reset_all_atom_map_numbers();

  // Jan 2004. Invalid directional bonds on guanidine C(/N)(\N)=N/N=C/c(ccc(c1Cl)O)c1
  // Need do clean those up before trying to standardise a guanidine

  if (remove_invalid_directional_bonds_from_input)
    molecule_changed = m.remove_invalid_directional_bonds();

  if (chemical_standardisation.active())
    molecule_changed += chemical_standardisation.process(m);

  if (change_double_bonds_between_permanent_aromatic_to_single)
    molecule_changed += m.change_double_bonds_between_permanent_aromatic_to_single(0);

  if (remove_unnecessary_square_brackets)
    molecule_changed += do_remove_unnecessary_square_brackets(m);

  if (remove_all_possible_square_brackets)
    molecule_changed += do_remove_all_possible_square_brackets(m);

  // This is a kludge for the selimsteg infrastructure

  if (!return_smiles_before_filters) {  // nothing to do
  } else if (! IsOrganic(m)) {  // not writing this one
  } else { // organic subset, store it.
    result.smiles_before_filters << m.smiles() << ' ' << m.name();
  }

  if (remove_invalid_chiral_centres) {
    int nr = do_remove_invalid_chiral_centres(m);
    if (nr) {
      molecules_with_invalid_chiral_centres++;

      if (verbose > 1)
        cerr << "Removed " << nr << " invalid chiral centres\n";

      molecule_changed = 1;
    }
  } else if (remove_chiral_data_from_all_molecules)
    molecule_changed += m.remove_all_chiral_centres();
  else if (remove_chiral_centres_on.number_elements())
    RemoveChiralCentresOnMatchedAtoms(m, remove_chiral_centres_on);

  if (remove_non_organic_chirality)
    RemoveNonOrganicChirality(m);

  if (!ApplyAllFiltersInner(m, result.rejection_reason, molecule_changed)) {
    result.rejected = 1;
    return result;
  }

  if (truncate_name_at_first || truncate_name_at_last) {
    do_truncate_name(m, truncate_name_at_first, truncate_name_at_last);
  } else if (truncate_names_to_first_token) {
    do_truncate_names_to_first_token(m);
  } else if (substitute_for_whitespace_in_name.length()) {
    SubstituteForWhitespaceInName(m);
  } else if (name_token >= 0) {
    do_extract_token_as_name(m, name_token);
  } else if (change_name_rx.get() == nullptr) {  // not active.
  } else if (ChangeName(m, *change_name_rx)) {  // good, worked.
  } else {
    if (verbose)
      cerr << "Cannot change name according to rx '" << change_name_rx->pattern() << "'\n";
    result.error = 1;
    return result;
  }

  if (zero_pad_name != 0) {
    ZeroPadName(m);
  }

  if (prepend_to_name.length() > 0) {
    IWString tmp(prepend_to_name);
    tmp << m.name();
    m.set_name(tmp);
  }

  if (number_assigner.active())
    number_assigner.process(m);

  if (reflect_coordinates)
    molecule_changed += ReflectCoordinates(m);
  else if (invert_all_chiral_centres)
    molecule_changed += InvertAllChiralCentres(m);
  else if (find_all_chiral_centres || find_all_ring_chiral_centres)
    molecule_changed += FindAllChiralCentres(m);

  if (translation_specified) {
    DoTranslation(m);
    ++molecule_changed;
  }

  if (rotation_specified) {
    m.rotate_atoms(rotation, rotation_angle);
    ++molecule_changed;
  }

  if (charge_assigner.active()) {
    molecule_changed += charge_assigner.process(m);
  }

  if (donor_acceptor_assigner.active()) {
    molecule_changed += donor_acceptor_assigner.process(m);
  }

  if (make_all_implicit_hydrogens_explicit) {
    if (m.make_implicit_hydrogens_explicit()) {
      molecules_to_which_hydrogens_were_added++;
      molecule_changed++;
    }
  } else if (atoms_for_implicit_hydrogens.number_elements()) {
    if (MakeImplicitHydrogensExplicit(m)) {
      molecules_to_which_hydrogens_were_added++;
      molecule_changed++;
    }
  }

  if (hydrogens_last) {
    if (m.MoveToEndOfConnectionTable(hydrogens_last)) {
      ++molecule_changed;
    }
  }

  if (fileconv_partial_charge_type != kNoChargeCalculation)
    ComputePartialCharges(m);

  if (appends_to_be_done)
    DoAppends(m);

  if (molecule_changed) {
    if (molecule_changed_string.length()) {
      IWString tmp(m.name());
      tmp.append_with_spacer(molecule_changed_string);
      m.set_name(tmp);
    }
  }

#ifdef IMPLEMENT_THIS_SOMETIME
  the problem is that this function no longer sees the Molecule_Output_Object
  if (molecule_to_fragments && m.number_fragments() > 1) {
    resizable_array_p<Molecule> components;
    m.create_components(components);

    const IWString name = m.name();
    for (int i = 0; i < components.number_elements(); ++i) {
      IWString myname = name;
      myname << " FRAG:" << i;
      components[i]->set_name(myname);
      output_object.write(*components[i]);
    }
  }
#endif

  return result;
}

bool
ValidAtomicNumber(const atomic_number_t z) {
  return z > 0 && z <= HIGHEST_ATOMIC_NUMBER;
}

const Element*
RecogniseAsElement(const const_IWSubstring& o) {
  const Element* rc = get_element_from_symbol_no_case_conversion(o);

  if (nullptr == rc)
    return nullptr;

  atomic_number_t z = rc->atomic_number();
  if (! ValidAtomicNumber(z)) {
    return nullptr;
  }

  return rc;
}

const Element*
RecogniseAsElementOrAtomicNumber(const const_IWSubstring& o) {
  int z;
  if (!o.numeric_value(z)) {
    return RecogniseAsElement(o);
  }
  if (! ValidAtomicNumber(z)) {
    return nullptr;
  }

  return get_element_from_atomic_number(z);
}

int
FileconvConfig::ParseOrganicSpecification(Command_Line& cl, char flag) {
  int i = 0;
  const_IWSubstring o;
  while (cl.value(flag, o, i++)) {
    if ("DEF" == o || "def" == o) {
      filter_for_disallowed_elements = 1;
    } else if (o.starts_with("allow:")) {
      o.remove_leading_chars(6);

      const Element* e = RecogniseAsElementOrAtomicNumber(o);

      if (nullptr == e) {
        cerr << "Unrecognised element allow:'" << o << "'\n";
        return 0;
      }

      allowed_elements.set_allow(e->atomic_number(), 1);

      if (verbose)
        cerr << "Atomic number " << e->atomic_number() << " allowed\n";
    } else if ("none" == o) {
      output_organic_only = 1;
      continue;
    } else if (o == "nometal") {
      allowed_elements.exclude_metals();
      filter_for_disallowed_elements = 1;
      continue;
    } else if (o == "okOBO") {
      ok_boron_special_case = 1;
      allowed_elements.set_allow(5, 1);
      filter_for_disallowed_elements = 1;
      continue;
    } else if ("help" == o) {
      cerr << "The -O option is troublesome due to it's long history and evolution\n";
      cerr << "Options can also be order dependent\n";
      cerr << " -O def         reject if any of the non-OK non-organics are present\n";
      cerr << " -O allow:El    temporarily allow El as an OK non-organic\n";
      cerr << " -O none        reject if any non-organic atoms present\n";
      cerr << " -O El          element El becomes fully organic for the course of the run\n";
      cerr << " -O nometal     exclude any molecule containing metal atoms\n";
      cerr << " -O okOBO       exclude Boron, except if it is as O-B-O form\n";
      exit(1);
    } else {
      Element_Matcher* e = new Element_Matcher;
      if (!e->construct_from_string(o)) {
        cerr << "Unrecognised -O qualifier '" << o << "'\n";
        delete e;
        return 0;
      }

      ok_non_organics.add(e);

      const Element* ele = e->element();
      if (nullptr != ele)
        allowed_elements.set_allow(ele->atomic_number(), 1);

      if (verbose)
        cerr << "ok_non_organics added '" << o << "'\n";

      output_organic_only = 1;
    }
  }

  if (verbose)
    cerr << "Initialised organic restrictions\n";

  return 1;
}

/*
  The syntax for the -T switch is -T dx,dy,dz
*/

int
FileconvConfig::GetTranslations(Command_Line& cl, const char cflag) {
  if (! cl.option_present(cflag)) {
    return 1;
  }

  const IWString token = cl.option_value(cflag);
  if (token == "origin") {
    translate_to_origin = 1;
    translation_specified = 1;
    return 1;
  }

  int i = 0;
  IWString d;
  if (!token.nextword(d, i, ',')) {
    cerr << "Cannot extract dx from -T switch\n";
    return 0;
  }

  if (!d.numeric_value(dx)) {
    cerr << "The value for dx must be a floating point number\n";
    return 0;
  }

  if (!token.nextword(d, i, ',')) {
    cerr << "Cannot extract dy from -T switch\n";
    return 0;
  }

  if (!d.numeric_value(dy)) {
    cerr << "The value for dy must be a floating point number\n";
    return 0;
  }

  if (!token.nextword(d, i, ',')) {
    cerr << "Cannot extract dz from -T switch\n";
    return 0;
  }

  if (!d.numeric_value(dz)) {
    cerr << "The value for dz must be a floating point number\n";
    return 0;
  }

  if (verbose) {
    cerr << "Molecules translated " << dx << ", " << dy << " ," << dz << '\n';
  }

  translation_specified = 1;

  return 1;
}

// `s` must be of the form <number>-<smiles>
//  isotope number and a molecule.
int
FileconvConfig::ParseFragmentAddToIsotope(const const_IWSubstring& s) {
  const_IWSubstring iso, smiles;
  isotope_t isotope;
  if (! s.split(iso, '-', smiles) || iso.empty() || smiles.empty() ||
      ! iso.numeric_value(isotope)) {
    return 0;
  }

  Molecule m;
  if (! m.build_from_smiles(smiles)) {
    cerr << "FileconvConfig::ParseFragmentAddToIsotope:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  add_to_isotopic_atom[isotope] = m;

  return 1;
}

int
FileconvConfig::GetIsotopeDirectives(Command_Line& cl, char flag) {
  if (!cl.option_present(flag)) {
    return 1;
  }

  IWString tmp;
  for (int i = 0; cl.value(flag, tmp, i); ++i) {
    if ('0' == tmp || tmp.starts_with("excl") || "none" == tmp) {
      exclude_isotopes = 1;
      if (verbose)
        cerr << "Molecules containing isotopes will be excluded\n";
    } else if ('1' == tmp || tmp.starts_with("incl")) {
      exclude_isotopes = 0;
      if (verbose)
        cerr << "No action taken on molecules containing isotopes\n";
    } else if ("convert" == tmp || "change" == tmp) {
      convert_isotopes = 1;
      if (verbose)
        cerr << "All isotopic atoms will be converted to their non isotopic form\n";
    } else if ("CONVERT" == tmp || "CHANGE" == tmp) {
      convert_isotopes = 2;
      if (verbose)
        cerr << "All isotopic atoms will be converted to their non isotopic form, implicit H recomputed\n";
    } else if (tmp.starts_with("add:")) {
      tmp.remove_leading_chars(4);
      if (! ParseFragmentAddToIsotope(tmp)) {
        cerr << "Invalid fragment add to isotope directive '" << tmp << "'\n";
        return 0;
      }
    } else if (tmp.starts_with("convert=") || tmp.starts_with("change=")) {
      if (tmp.starts_with("convert="))
        tmp.remove_leading_chars(8);
      else
        tmp.remove_leading_chars(7);

      int a;      // the starting isotope
      int b = 0;  // the new isotopic value
      if (tmp.contains(',')) {
        const_IWSubstring s1, s2;
        tmp.split(s1, ',', s2);
        if (!s1.numeric_value(a) || a < 0 || !s2.numeric_value(b) || b < 0) {
          cerr << "The isotopic conversion directive must contain a valid set of isotopes '" << tmp
               << "' invalid\n";
          DisplayDashIOptions('I', cerr);
        }
      } else {
        if (!tmp.numeric_value(a) || a < 0) {
          cerr << "The isotopic conversion directive must contain a valid set of isotopes '" << tmp
               << "' invalid\n";
          DisplayDashIOptions('I', cerr);
        }
      }

      if (verbose)
        cerr << "Will convert all isotope " << a << " values to isotope " << b << "\n";

      convert_specific_isotopes.add(a);
      convert_specific_isotopes_new_isotope.add(b);
    } else if (tmp.starts_with("smarts=") ||
               tmp.starts_with("smarts:"))  // must be of the form 'smarts,new'
    {
      tmp.remove_leading_chars(7);

      const int ndx = tmp.rindex(',');

      if (ndx < 0) {
        cerr << "THe smarts for converting isotopes must have be 'smarts,new_isotope', '" << tmp
             << "' invalid\n";
        DisplayDashIOptions('I', cerr);
      }

      const_IWSubstring s;
      tmp.from_to(0, ndx - 1, s);

      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();

      if (!q->create_from_smarts(s)) {
        cerr << "Invalid smarts in isotopic directives '" << tmp << "'\n";
        return 1;
      }

      tmp.from_to(ndx + 1, tmp.length() - 1, s);
      int x;
      if (!s.numeric_value(x) || x < 0) {
        cerr << "Invalid isotope in isotopic directives '" << tmp << "' '" << s << "'\n";
        return 1;
      }

      convert_specific_isotopes_query.add(q.release());
      convert_specific_isotopes_query_new_isotope.add(x);
    } else if (tmp.starts_with("alliso=")) {
      tmp.remove_leading_chars(7);
      if (!tmp.numeric_value(convert_all_isotopes_to) || convert_all_isotopes_to < 0) {
        cerr << "The convert all isotopes 'alliso=' directive must have a valid isotopic value, '"
             << tmp << "' invalid\n";
        DisplayDashIOptions('I', cerr);
      }

      if (verbose)
        cerr << "All existing isotopes converted to " << convert_all_isotopes_to << '\n';
    } else if ("remove" == tmp) {
      remove_isotopic_atoms = 1;
      if (verbose)
        cerr << "All isotopically atoms will be removed\n";
    } else if (tmp == "atype") {
      isotope_is_atom_type = 1;
      if (verbose) {
        cerr << "Atoms will be labelled with their atom type\n";
      }
    } else if ("help" == tmp) {
      DisplayDashIOptions('I', cerr);
    } else {
      cerr << "Unrecognised -I qualifier '" << tmp << "'\n";
      Usage(1);
    }
  }

  need_to_consider_isotopes = 1;

  return 1;
}

int
FileconvConfig::GetChargeCalculation(Command_Line& cl, char flag) {
  if (!cl.option_present(flag)) {
    return 1;
  }

  const_IWSubstring q;
  for (int i = 0; cl.value(flag, q, i); ++i) {
    if ("help" == q) {
      cerr << "The following -Q qualifiers are recognised\n";
      cerr << " -Q gast    Gasteiger partial charges\n";
      cerr << " -Q gh      Gasteiger Huckel partial charges\n";
      cerr << " -Q abraham Abraham partial charges\n";

      return 0;
    } else if (q.starts_with("gast")) {
      fileconv_partial_charge_type = kGasteiger;
      if (verbose)
        cerr << "Will compute Gasteiger partial charges\n";
    } else if ("gh" == q) {
      fileconv_partial_charge_type = kGasteigerHuckel;
      if (verbose)
        cerr << "Will compute Gasteiger Huckel partial charges\n";
    } else if ("abraham" == q) {
      fileconv_partial_charge_type = kAbraham;
      if (verbose)
        cerr << "Will compute Abraham partial charges\n";
    } else {
      cerr << "Unrecognised partial charge specification '" << q << "'\n";
      Usage(12);
    }
  }
  return 1;
}

int
FileconvConfig::ParseMiscOptions(Command_Line& cl, char flag) {
  if (!cl.option_present(flag)) {
    return 1;
  }

  const_IWSubstring y;
  for (int i = 0; cl.value(flag, y, i); ++i) {
    if (y.starts_with("appchg=")) {
      molecule_changed_string = y;
      molecule_changed_string.remove_leading_chars(7);

      if (verbose)
        cerr << "Will append '" << molecule_changed_string
             << "' to molecules whose connection table is changed\n";
    } else if ("pblen" == y) {
      print_bond_lengths = 1;
      if (verbose)
        cerr << "Bond lengths will be printed\n";
    } else if ("pbang" == y) {
      print_bond_angles = 1;
      if (verbose)
        cerr << "Bond angles will be printed\n";
    } else if ("ptor" == y) {
      print_torsions = 1;
      if (verbose)
        cerr << "Torsions be printed\n";
    } else if (y == "pmaxd") {
      print_max_atom_separation = 1;
      if (verbose) {
        cerr << "Max atom separation printed\n";
      }
    } else if (y.starts_with("namerx=")) {
      y.remove_leading_chars(7);
      re2::StringPiece tmp(y.data(), y.length());
      name_rx.reset(new RE2(tmp));
      if (!name_rx->ok()) {
        cerr << "Invalid name regular expression '" << y << "'\n";
        return 0;
      }

      if (verbose)
        cerr << "Will only process molecules whose name matches '" << name_rx->pattern() << "'\n";
    } else if (y.starts_with("grepvname=")) {
      y.remove_leading_chars(10);
      re2::StringPiece tmp(y.data(), y.length());
      grep_v_name_rx.reset(new RE2(tmp));
      if (! grep_v_name_rx->ok()) {
        cerr << "Invalid discard name regular expression '" << y << "'\n";
        return 0;
      }
      if (verbose) {
        cerr << "Will skip molecules whose name matches '" << grep_v_name_rx->pattern() << "'\n";
      }
    } else if ("ftn" == y) {
      truncate_names_to_first_token = 1;

      if (verbose)
        cerr << "Will truncate molecule names to the first token\n";
    } else if (y.starts_with("tfirst=")) {
      y.remove_leading_chars(7);
      truncate_name_at_first = y[0];
    } else if (y.starts_with("tlast=")) {
      y.remove_leading_chars(6);
      truncate_name_at_last = y[0];
    } else if (y.starts_with("chname=") || y.starts_with("CHNAME=")) {
      if ('C' == y[0])  // uppercase form
        change_name_rx_must_match = 1;

      y.remove_leading_chars(7);
      re2::StringPiece tmp(y.data(), y.length());
      change_name_rx.reset(new RE2(tmp));
      if (!change_name_rx->ok()) {
        cerr << "Invalid change name regexp '" << y << "'\n";
        return 3;
      }

      if (0 == change_name_rx->NumberOfCapturingGroups()) {
        cerr << "The chname regular expression must have capturing goup(s) '"
             << change_name_rx->pattern() << "' invalid\n";
        return 3;
      }

      if (verbose)
        cerr << "Will change molecule names to part that matches rexexp '"
             << change_name_rx->pattern() << "'\n";

      if (change_name_rx->NumberOfCapturingGroups() > 10) {
        cerr << "Too many capturing groups in chname=" << y << "\n";
        return 1;
      }
    } else if (y.starts_with("zpad=")) {
      y.remove_leading_chars(5);
      if (! y.numeric_value(zero_pad_name) || zero_pad_name == 0) {
        cerr << "Invalid zero pad directive '" << y << "'\n";
        return 0;
      }
      if (verbose)  {
        cerr << "Will zero pad names to " << zero_pad_name << " width\n";
      }
    } else if (y.starts_with("nsubw=")) {
      y.remove_leading_chars(6);
      substitute_for_whitespace_in_name = y;

      if (1 != substitute_for_whitespace_in_name.length()) {
        cerr << "Sorry, the substitute for whitespace value must be a single character only '" << y
             << "' is invalid\n";
        return 3;
      }

      if (verbose)
        cerr << "Whitespace in names changed to '" << substitute_for_whitespace_in_name << "'\n";
    } else if (y.starts_with("NPREPEND=")) {
      y.remove_leading_chars(9);
      prepend_to_name = y;
      if (verbose)
        cerr << "Will prepend '" << prepend_to_name << "' to all names\n";
    } else if (y.starts_with("ntoken=")) {
      y.remove_leading_chars(7);
      if (!y.numeric_value(name_token) || name_token < 1) {
        cerr << "Invalid name token specification '" << y << "'\n";
        return 1;
      }
      if (verbose) {
        cerr << "Will use the " << name_token << " token of the input name as the output name\n";
      }
      name_token--;
    } else if (y.starts_with("maxplen=")) {
      y.remove_leading_chars(8);
      if (!y.numeric_value(max_path_length) || max_path_length < 1) {
        cerr << "The maxplen= directive must be followed by a +ve whole number\n";
        return 5;
      }

      if (verbose)
        cerr << "Will discard molecules with longest path > " << max_path_length << "\n";
    } else if ("ada" == y) {
      change_double_bonds_between_permanent_aromatic_to_single = 1;
      if (verbose)
        cerr << "Will change double bonds between permanent aromatic atoms to single\n";
    } else if (y.starts_with("B4F=")) {
      return_smiles_before_filters = 1;
      if (verbose) {
        cerr << "Pre filtered smiles will be written\n";
      }
    } else if ("aclf" == y) {
      atom_count_includes_only_atoms_in_largest_fragment = 1;

      if (verbose)
        cerr << "All atom counts will be with regard to larges fragment only\n";
    } else if ("nbvm" == y) {
      set_display_abnormal_valence_messages(0);
      set_display_strange_chemistry_messages(0);
      if (verbose)
        cerr << "Abnormal valence messages suppressed\n";
    } else if ("okbvi" == y) {
      ok_bad_valence_on_isotopically_labelled = 1;
      skip_molecules_with_abnormal_valences = 1;
      if (verbose)
        cerr << "During valence check will ignore bad valences on isotopically labelled atoms\n";
    } else if ("nhsqb" == y) {
      set_explicit_hydrogens_need_square_brackets_in_smiles(0);
      if (verbose)
        cerr << "Explicit hydrogens in smiles written without square brackets\n";
    } else if ("rmsqb" == y) {
      remove_unnecessary_square_brackets = 1;
      if (verbose)
        cerr << "Will remove unnecessary square bracketed elements\n";
    } else if ("FHrmsqb" == y) {
      remove_all_possible_square_brackets = 1;
      if (verbose)
        cerr << "Will assume all H counts are OK and remove smiles square brackets\n";
    } else if ("rmamap" == y) {
      reset_atom_map_numbers = 1;
      if (verbose)
        cerr << "Will not include atom map with smiles\n";
    } else if (y.starts_with("fixarom=")) {
      y.remove_leading_chars(8);

      Substructure_Query* q = new Substructure_Query;
      if (!q->create_from_smarts(y)) {
        cerr << "Invalid smarts specification of atom in ring to be fixed '" << y << "'\n";
        delete q;
        return 1;
      }

      aromatise_these_rings.add(q);
    } else if ("nbrminv" == y) {
      set_invalidate_bond_list_ring_info_during_invalidate_ring_info(0);
    } else if ("iso2amap" == y) {
      convert_isotopes_to_atom_map_numbers = 1;
      if (verbose)
        cerr << "Will convert isotopes to atom map numbers\n";
    } else if ("amap2iso" == y) {
      convert_atom_map_numbers_to_isotopes = 1;
      if (verbose)
        cerr << "Will convert atom map numbers to isotopes\n";
    } else if (y.starts_with("rmatom=")) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      y.remove_leading_chars(7);
      if (! q->create_from_smarts(y)) {
        cerr << "FileconvConfig::invalid smarts for atoms to excise '" << y << "'\n";
        return 0;
      }
      remove_atom << q.release();
      if (verbose) {
        cerr << "Will remove atoms matching '" << y << "'\n";
      }
    } else if (y.starts_with("rot=")) {
      y.remove_leading_chars(4);
      if (! InitialiseRotation(y)) {
        cerr << "FileconvConfig::invalid rotation '" << y << "'\n";
        return 0;
      }
    } else if (y.starts_with("atype=")) {
      y.remove_leading_chars(6);
      if (! atom_typing.build(y)) {
        cerr << "FileconvConfig::cannot initialise atom typing '" << y << "'\n";
        return 0;
      }
    } else if ("help" == y) {
      DisplayDashYOptions(cerr, 2);
    } else {
      cerr << "Unrecognised -Y qualifier '" << y << "'\n";
      return 6;
    }
  }

  return 1;
}

// `buffer` must be of the form x,y,z,angle
// Use these to initialise rotation and rotation_angle.
// Note the argument is passed by value so we have our own copy.
int
FileconvConfig::InitialiseRotation(const_IWSubstring buffer) {
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; buffer.nextword(token, i, ','); ++col) {
    float value;
    if (token.empty() || ! token.numeric_value(value)) {
      cerr << "FileconvConfig::InitialiseRotation:invalid rotation '" << token << "'\n";
      return 0;
    }
    if (col == 0) {
      rotation.set_x(value);
    } else if (col == 1) {
      rotation.set_y(value);
    } else if (col == 2) {
      rotation.set_z(value);
    } else if (col == 3) {
      rotation_angle = value * DEG2RAD;
    } else {
      cerr << "FileconvConfig::InitialiseRotation:bad rotation '" << token << "'\n";
      return 0;
    }
  }

  rotation.normalise();

  rotation_specified = 1;

  return 1;
}

void
DisplayChiralityOptions() {
//cerr << "  -s discard     discard all chiral information in the input\n";
  cerr << "  -s good        ignore erroneous chiral input\n";
  cerr << "  -s find        identify all stereo centres in the molecules (random chirality)\n";        
  cerr << "  -s rfind       identify stereo centres for ring atoms only\n";
  cerr << "  -s invert      invert all stereo centres\n";
  cerr << "  -s 1           include chiral info in output (default)\n";
  cerr << "  -s 0           exclude chiral info from output\n";
  cerr << "  -s remove      remove all chirality upon read\n";
  cerr << "  -s rmchiral=... remove chiral centres on atoms that match query(s)\n";
  cerr << "  -s 'rmchiral=SMARTS:[ND2H]' will remove 2 connected Nitrogen chiral centres\n";
  cerr << "  -s 'rmchiral=SMARTS:[AR3]' will remove wedge bonds in adamantyl systems\n";
  cerr << "  -s rmbad       remove chiral centres that really aren't chiral\n";
  cerr << "  -s xctb        remove all cis-trans bonding information\n";
  cerr << "  -s rmbctb      remove invalid cis-trans bonding information\n";
//cerr << "  -s d3d         discern chirality from the 3D structure\n";
  cerr << "  -s nonocm      no non organic chirality messages\n";
  cerr << "  -s rmnoc       remove chirality from non organic atoms\n";
  cerr << "  -s maxc=<n>    discard molecules with more than <n> chiral centres\n";

  return;
}

int
FileconvConfig::GetChiralityInstructions(Command_Line& cl, char flag) {
  if (!cl.option_present(flag)) {
    return 1;
  }

  MDL_File_Supporting_Material* mdlfos = global_default_MDL_File_Supporting_Material();

  IWString tmp;
  for (int i = 0; cl.value(flag, tmp, i); ++i) {
    if ("help" == tmp) {
      DisplayChiralityOptions();
      exit(1);
    }

    if ('1' == tmp || tmp.starts_with("incl")) {
      set_include_chiral_info_in_smiles(1);
      if (verbose)
        cerr << "Chiral information will be written\n";
    } else if ('0' == tmp || tmp.starts_with("excl")) {
      set_include_chiral_info_in_smiles(0);
      mdlfos->set_include_chiral_info_in_mdl_outputs(0);
      if (verbose)
        cerr << "Chiral information will NOT be written\n";
    } else if ("find" == tmp) {
      find_all_chiral_centres = 1;
      if (verbose)
        cerr << "Chiral centres will be identified\n";
    } else if ("rfind" == tmp) {
      find_all_ring_chiral_centres = 1;
      if (verbose)
        cerr << "Chiral centres in ring systems will be identified\n";
    } else if ("invert" == tmp) {
      invert_all_chiral_centres = 1;
      if (verbose)
        cerr << "Chiral centres will be inverted\n";
    } else if ("reflect" == tmp) {
      reflect_coordinates = 1;
      if (verbose)
        cerr << "Will reflect coordinates and invert chiral centres\n";
    } else if ("discard" == tmp) {
      set_ignore_all_chiral_information_on_input(1);
      if (verbose)
        cerr << "All chiral information in input will be discarded\n";
    } else if ("remove" == tmp) {
      remove_chiral_data_from_all_molecules = 1;
      if (verbose)
        cerr << "Will strip all chiral information from molecules\n";
    } else if ("good" == tmp || "ignore" == tmp) {
      set_ignore_incorrect_chiral_input(1);
      if (verbose)
        cerr << "Incorrect chiral specifications will be ignored\n";
    } else if ("rmbad" == tmp) {
      remove_invalid_chiral_centres = 1;

      if (verbose)
        cerr << "Invalid chiral centres will be removed\n";
    } else if (tmp.starts_with("rmchiral=")) {
      const_IWSubstring foo(tmp);  // process_cmdline_token may change its argument

      foo.remove_leading_chars(9);

      if (!process_cmdline_token(' ', foo, remove_chiral_centres_on, verbose)) {
        cerr << "Cannot read chiral centre removal query(s) '" << tmp << "'\n";
        return 3;
      }
    } else if ("rmctb" == tmp || "xctb" == tmp) {
      remove_directional_bonds_from_input = 1;

      if (verbose)
        cerr << "Directional bonds removed from molecules\n";
    } else if ("rmbctb" == tmp) {
      remove_invalid_directional_bonds_from_input = 1;

      if (verbose)
        cerr << "Invalid directional bonds will be removed\n";
    } else if ("nonocm" == tmp) {
      mdlfos->set_display_non_organic_chirality_messages(0);
      if (verbose)
        cerr << "Will not display messages about non organic chirality\n";
    } else if ("rmnoc" == tmp) {
      remove_non_organic_chirality = 1;

      if (verbose)
        cerr << "Chirality will be removed from non organic atoms\n";
    } else if (tmp.starts_with("maxc=")) {
      tmp.remove_leading_chars(5);
      if (!tmp.numeric_value(max_chiral_centres) || max_chiral_centres < 0) {
        cerr << "The max chiral centre count (-s maxc=) directive must be a whole number\n";
        return 1;
      }

      if (verbose)
        cerr << "Will discard molecules having more than " << max_chiral_centres
             << " chiral centres\n";
    } else {
      cerr << "Unrecognised stereo/chiral specifier (-s) '" << tmp << "'\n";
      Usage(32);
    }
  }

  for (int i = 0; i < remove_chiral_centres_on.number_elements(); i++) {
    remove_chiral_centres_on[i]->set_find_unique_embeddings_only(1);
  }

  return 1;
}

int
FetchFragSmarts(const const_IWSubstring f,
                const char* largest_or_smallest,
                resizable_array_p<Substructure_Query>& queries,
                int verbose) {
  Substructure_Query* tmp = new Substructure_Query(f);
  if (!tmp->create_from_smarts(f)) {
    cerr << "Cannot create query from smarts '" << f << "'\n";
    delete tmp;
    return 0;
  }

  queries.add(tmp);

  if (verbose)
    cerr << "Will keep " << largest_or_smallest << " fragment which matches query '" << f << "'\n";

  return 1;
}

int
FetchFragQueries(const_IWSubstring f,
                 const char* largest_or_smallest,
                 resizable_array_p<Substructure_Query>& queries,
                 int verbose) {
  if (f.starts_with("F:")) {
    f.remove_leading_chars(2);
    if (!queries_from_file(f,
                           queries,
                           1,
                           verbose))  // let's assume queries in same directory as file
    {
      cerr << "Could not process 'F:" << f << "'\n";
      return 0;
    }

    return 1;
  }

  if (f.starts_with("M:")) {
    f.remove_leading_chars(2);

    if (!queries_from_file_of_molecules(f, queries, verbose))
      return 0;

    return queries.number_elements();
  }

  Substructure_Query* tmp = new Substructure_Query(f);
  if (!tmp->read(f)) {
    cerr << "FetchFragQueries: cannot read query from '" << f << "'\n";
    delete tmp;
    return 0;
  }

  queries.add(tmp);
  if (verbose)
    cerr << "Will keep " << largest_or_smallest << " fragment which matches query '" << f << "'\n";

  return 1;
}

int
FileconvConfig::ParseFragmentSizeWindow(const const_IWSubstring& f) {
  if (! f.contains(',')) {
    cerr << "FileconvConfig::ParseFragmentSizeWindow:must be of the form min,max\n";
    return 0;
  }
  const_IWSubstring fmin, fmax;
  if (! f.split(fmin, ',', fmax)) {
    cerr << "FileconvConfig::ParseFragmentSizeWindow:cannot split on , char\n";
    return 0;
  }
  if (fmin.length() > 0) {
    if (! fmin.numeric_value(fragment_selector_lower_atom_count) ||
          fragment_selector_lower_atom_count < 1) {
      cerr << "FileconvConfig:ParseFragmentSizeWindow:invalid lower fragment atom count '" << f << "'\n";
      return 0;
    }
  }

  if (fmax.length() > 0) {
    if (! fmax.numeric_value(fragment_selector_upper_atom_count) ||
          fragment_selector_upper_atom_count < fragment_selector_lower_atom_count) {
      cerr << "FileconvConfig:ParseFragmentSizeWindow:invalid upper fragment atom count '" << f << "'\n";
      return 0;
    }
  }

  if (verbose) {
    cerr << "Fragment count window " << fragment_selector_lower_atom_count << ',' << fragment_selector_upper_atom_count << '\n';
  }

  return 1;
}

int
FileconvConfig::GetFragmentSpecifications(Command_Line& cl) {
  if (cl.option_present('f')) {
    const_IWSubstring f;
    int i = 0;
    while (cl.value('f', f, i++)) {
      if ('l' == f || f.starts_with("large") || f.starts_with("LARGE")) {
        reduce_to_largest_fragment = 1;
        if (verbose)
          cerr << "The largest fragment will be retained\n";
      } else if ("alarge" == f) {
        reduce_to_all_largest_fragments = 1;
        if (verbose)
          cerr << "Will keep all fragments the same size as the largest\n";
      } else if ("rmlarge" == f) {
        remove_largest_fragment = 1;
        if (verbose)
          cerr << "The largest fragment(s) will be removed\n";
      } else if (f.starts_with("rmlarge=")) {
        f.remove_leading_chars(8);
        if (!f.numeric_value(remove_largest_fragment) || remove_largest_fragment < 1) {
          cerr << "The rmlarge= directive must be followed by a whole +ve integer\n";
          return 0;
        }
        if (verbose)
          cerr << "Will remove the " << remove_largest_fragment << " largest fragment(s)\n";
      } else if ("rmsmall" == f) {
        remove_smallest_fragment = 1;
        if (verbose)
          cerr << "The smallest fragment will be removed\n";
      } else if (f.starts_with("rmsmall=")) {
        f.remove_leading_chars(8);
        if (!f.numeric_value(remove_smallest_fragment) || remove_smallest_fragment < 1) {
          cerr << "The rmsmall= directive must be followed by a whole +ve integer\n";
          return 0;
        }
        if (verbose)
          cerr << "Will remove the " << remove_smallest_fragment << " smallest fragment(s)\n";
      } else if ("lo" == f) {
        reduce_to_largest_organic_fragment = 1;

        if (verbose)
          cerr << "Will reduce to the largest organic fragment\n";
      } else if (f == "keepsmall") {
        keep_smallest_fragments = 1;
        if (verbose) {
          cerr << "Will retain on the smallest fragment\n";
        }
      } else if (f.starts_with("keepsmall=")) {
        f.remove_leading_chars(10);
        if (! f.numeric_value(keep_smallest_fragments) || keep_smallest_fragments <= 0) {
          cerr << "The keepsmallest= directive must be followed by a whole +ve integer\n";
          return 0;
        }
        if (verbose) {
          cerr << "Will retain only the " << keep_smallest_fragments << " smallest fragment(s)\n";
        }
      } else if ("lod" == f) {
        reduce_to_largest_organic_fragment_carefully = 1;

        if (verbose)
          cerr << "Will reduce to largest organic fragment using desirability rules\n";
      } else if ("allo" == f) {
        keep_all_organic_fragments = true;

        if (verbose)
          cerr << "Will remove all fragments with non organic atoms\n";
      } else if (f.starts_with("Q:")) {
        f.remove_leading_chars(2);
        if (!FetchFragQueries(f, "largest", largest_fragment_queries, verbose))
          Usage(44);
      } else if (f.starts_with("q:")) {
        f.remove_leading_chars(2);
        if (!FetchFragQueries(f, "smallest", smallest_fragment_queries, verbose))
          Usage(45);
      } else if (f.starts_with("SMARTS:")) {
        f.remove_leading_chars(7);
        if (!FetchFragSmarts(f, "largest", largest_fragment_queries, verbose))
          Usage(46);
      } else if (f.starts_with("smarts:")) {
        f.remove_leading_chars(7);
        if (!FetchFragSmarts(f, "smallest", smallest_fragment_queries, verbose))
          Usage(47);
      } else if (f.starts_with("ALL:")) {
        f.remove_leading_chars(4);
        if (!FetchFragSmarts(f, "all", keep_fragment_queries, verbose))
          Usage(32);
      } else if (f.starts_with("rm:")) {
        f.remove_leading_chars(3);
        if (!FetchFragSmarts(f, "rm", remove_fragment_queries, verbose))
          Usage(33);
      } else if (f.starts_with("rmle=")) {
        f.remove_up_to_first('=');
        if (!f.numeric_value(remove_fragments_this_size_or_smaller) ||
            remove_fragments_this_size_or_smaller < 1) {
          cerr << "The rmle= qualifier must be a whole positive number\n";
          return 11;
        }

        if (verbose)
          cerr << "Will remove all fragments with " << remove_fragments_this_size_or_smaller
               << " or fewer atoms\n";
      } else if (f.numeric_value(fragment_count)) {
        if (fragment_count < 1) {
          cerr << "The -f <number> directive must be followed by a whole positive number\n";
          Usage(8);
        }
        if (verbose)
          cerr << "A maximum of " << fragment_count << " fragments will be written\n";
      } else if ("RMDUP" == f) {
        remove_duplicate_fragments = 1;

        if (verbose)
          cerr << "Will remove duplicate fragments\n";
      } else if (f.starts_with("saltfile=")) {
        f.remove_leading_chars(9);
        if (!known_fragment_data.read_known_salts(f)) {
          cerr << "Cannot read known salts '" << f << "'\n";
          return 9;
        }
      } else if (f.starts_with("parentfile=")) {
        f.remove_leading_chars(11);
        if (!known_fragment_data.read_known_parents(f)) {
          cerr << "Cannot read known parents '" << f << "'\n";
          return 9;
        }
      } else if ("kmfok" == f) {
        known_fragment_data.set_only_check_molecular_formula(1);
      } else if ("kpallsalt" == f) {
        known_fragment_data.set_remove_everything_if_all_fragments_match(0);
        if (verbose)
          cerr << "Will not change molecule consisting of all salts\n";
      } else if (f.starts_with("rmxt=")) {
        f.remove_leading_chars(5);

        if (!f.numeric_value(discard_molecule_if_multiple_fragments_larger_than) ||
            discard_molecule_if_multiple_fragments_larger_than < 1) {
          cerr << "The rmxt qualifier must be followed by a whole +ve number\n";
          return 3;
        }

        if (verbose)
          cerr << "Will discard molecules with multiple fragments with more than "
               << discard_molecule_if_multiple_fragments_larger_than << " atoms\n";
      } else if ("rmxt" == f) {
        discard_molecule_if_multiple_fragments_larger_than = 16;

        if (verbose)
          cerr << "Will discard molecules with multiple fragments with more than "
               << discard_molecule_if_multiple_fragments_larger_than << " atoms\n";
      } else if (f.starts_with("dmxt=")) {
        f.remove_leading_chars(5);
        if (!f.numeric_value(mixture_if_largest_frags_differ_by) ||
            mixture_if_largest_frags_differ_by < 0) {
          cerr << "The delta for mixture (dmxt=) specifier must be a whole non-negative number\n";
          return 3;
        }

        if (verbose)
          cerr << "A molecule where the largest fragments differ by only "
               << mixture_if_largest_frags_differ_by << " atoms will be considered a mixture\n";
      } else if (f.starts_with("manlf=")) {
        f.remove_leading_chars(6);
        if (!f.numeric_value(remove_molecules_with_non_largest_fragment_natoms) ||
            remove_molecules_with_non_largest_fragment_natoms < 0) {
          cerr << "The max atoms in non-largest fragment directive (manlf) must be a whole +ve "
                  "number\n";
          return 2;
        }

        if (verbose)
          cerr << "Will discard molecules where a non-largest fragment has more than "
               << remove_molecules_with_non_largest_fragment_natoms << " atoms\n";
      } else if (f.starts_with("klf=")) {
        f.remove_leading_chars(4);
        if (!f.numeric_value(strip_to_n_largest_fragments) || strip_to_n_largest_fragments < 1) {
          cerr << "The strip to N largest fragments 'klf=' directive must have a whole +ve "
                  "integer\n";
          return 2;
        }

        if (verbose)
          cerr << "Will strip to the " << strip_to_n_largest_fragments << " largest fragments\n";
      } else if ("sfs" == f) {
        sort_by_fragment_size = 1;

        if (verbose)
          cerr << "Will sort fragments by size\n";
      } else if (f.starts_with("RMF=")) {
        f.remove_leading_chars(4);
        tag_for_removed_fragments = f;

        if (0 == tag_for_removed_fragments.length()) {
          cerr << "The tag for removed fragments (RMF=) must be non zero length\n";
          return 2;
        }
      } else if (f.starts_with("WINDOW=")) {
        f.remove_leading_chars(7);
        if (! ParseFragmentSizeWindow(f)) {
          cerr << "Invalid fragment size window '" << f << "'\n";
          return 0;
        }
      } else if ("help" == f) {
        DisplayfOptions(1);
      } else {
        cerr << "Unrecognised -f directive '" << f << "'\n";
        return 0;
      }
    }

    //  We can have only one selection active - except for Known_Fragment_Data which can work in
    //  combination with other fragment selections

    int nsel =
        (fragment_count != 0) + (reduce_to_largest_fragment) + (reduce_to_all_largest_fragments) +
        (reduce_to_largest_organic_fragment) + (reduce_to_largest_organic_fragment_carefully) +
        (keep_all_organic_fragments) + (largest_fragment_queries.number_elements() > 0) +
        (smallest_fragment_queries.number_elements() > 0) +
        (remove_fragment_queries.number_elements() > 0) + (strip_to_n_largest_fragments > 0) +
        (keep_fragment_queries.number_elements() > 0) +
        (discard_molecule_if_multiple_fragments_larger_than > 0) + (known_fragment_data.active() +
        (remove_largest_fragment > 0) + (remove_smallest_fragment > 0) +
        (keep_smallest_fragments > 0) +
        (fragment_selector_lower_atom_count > 0 ||
         fragment_selector_upper_atom_count < std::numeric_limits<int>::max()) +
        remove_duplicate_fragments);

    //  cerr << "Fragment nsel " << nsel << " remove_fragment_queries " <<
    //  remove_fragment_queries.size() << ' ' << fragment_count << ' ' << reduce_to_largest_fragment
    //  << ' ' << reduce_to_all_largest_fragments << ' ' << reduce_to_largest_organic_fragment << '
    //  ' << reduce_to_largest_organic_fragment_carefully << ' ' <<
    //  largest_fragment_queries.number_elements() << ' ' <<
    //  smallest_fragment_queries.number_elements() << ' ' <<
    //  remove_fragment_queries.number_elements() << ' ' << strip_to_n_largest_fragments << ' ' <<
    //  keep_fragment_queries.number_elements() << ' ' << known_fragment_data.active() << '\n';

    if (1 == nsel)  // good, the easy case
      ;
    else if (remove_fragment_queries.number_elements() > 0 &&
             (reduce_to_largest_fragment || reduce_to_all_largest_fragments ||
              reduce_to_largest_organic_fragment || reduce_to_largest_organic_fragment_carefully ||
              keep_all_organic_fragments))
      ;
    else if (0 == nsel && sort_by_fragment_size)
      ;
    else if (discard_molecule_if_multiple_fragments_larger_than > 0 && 2 == nsel)
      ;
    else if (remove_largest_fragment > 0 || remove_smallest_fragment > 0)
      ;
    else if (keep_smallest_fragments > 0)
      ;
    else
      cerr << "You have specified more than one fragment selection criterion " << nsel
           << " beware of problems...\n";

    if (nsel || known_fragment_data.active() ||
        remove_molecules_with_non_largest_fragment_natoms >= 0 || sort_by_fragment_size) {
      need_to_call_process_fragments = 1;
    }

    if (tag_for_removed_fragments.length() && !reduce_to_largest_organic_fragment) {
      cerr << "Sorry, the tag for removed fragments option 'RMF=' only works with the 'lo' "
              "fragment selection qualifier\n";
      return 1;
    }

    //  if (known_fragment_data.active())
    //    known_fragment_data.debug_print (cerr);
  }

  if (cl.option_present('F')) {
    int i = 0;
    const_IWSubstring f;
    while (cl.value('F', f, i++)) {
      if (f.starts_with("mnsz=")) {
        f.remove_leading_chars(5);

        if (!f.numeric_value(min_size_max_fragment_count) || min_size_max_fragment_count < 1) {
          cerr << "The min size for the too many fragments (-F mnsz=) combination must be a whole "
                  "+ve number\n";
          Usage(3);
        }

        if (verbose)
          cerr << "For max fragment count, will only consider fragments with > "
               << min_size_max_fragment_count << " atoms\n";
      } else if ("help" == f) {
        DisplayFOptions(1);
      } else if (f.numeric_value(maximum_fragment_count)) {
        if (maximum_fragment_count < 1) {
          cerr << "The maximum fragment count (-F) must be a whole +ve number\n";
          Usage(3);
        }
        if (verbose)
          cerr << "Molecules containing more than " << maximum_fragment_count
               << " fragments will be skipped\n";
      } else {
        cerr << "Unrecognised -F qualitifier '" << f << "'\n";
        return 3;
      }
    }

    need_to_call_process_fragments = 1;
  }

  return 1;
}

int
FileconvConfig::GatherAppendSpecifications(Command_Line& cl, char flag) {
  if (!cl.option_present(flag)) {
    return 1;
  }

  const_IWSubstring p;
  for (int i = 0; cl.value(flag, p, i); ++i) {
    if ("AMW" == p) {
      append_molecular_weight_to_name = 1;
      if (verbose)
        cerr << "The molecular weight will be appended to the name\n";

      appends_to_be_done = 1;
    } else if ("AMWI" == p) {
      append_molecular_weight_ok_isotope_to_name = 1;
      if (verbose)
        cerr << "The isotopically aware molecular weight will be appended to the name\n";

      appends_to_be_done = 1;
    } else if ("MF" == p) {
      append_molecular_formula_to_name = 1;
      if (verbose)
        cerr << "The molecular formula will be appended to the name\n";

      appends_to_be_done = 1;
    } else if ("ISISMF" == p) {
      append_isis_molecular_formula_to_name = 1;
      if (verbose)
        cerr << "The ISIS-like molecular formula will be appended to the name\n";

      appends_to_be_done = 1;
    } else if ("NATOMS" == p) {
      append_natoms_to_name = 1;
      if (verbose)
        cerr << "The number of atoms will be appended to the name\n";

      appends_to_be_done = 1;
    } else if ("NBONDS" == p) {
      append_nbonds_to_name = 1;
      if (verbose)
        cerr << "The number of bonds will be appended to the name\n";

      appends_to_be_done = 1;
    } else if ("NRINGS" == p) {
      append_nrings_to_name = 1;
      if (verbose)
        cerr << "The number of rings will be appended to the name\n";

      appends_to_be_done = 1;
    } else if ("AROMR" == p) {
      append_aromatic_ring_count_to_name = 1;
      if (verbose)
        cerr << "The number of aromatic rings will be appended to the name\n";

      appends_to_be_done = 1;
    } else if ("EXACT" == p) {
      append_exact_mass_to_name = 1;
      if (verbose)
        cerr << "The exact mass will be appended to the name\n";

      appends_to_be_done = 1;
    } else if ("NFC" == p) {
      append_net_formal_charge_to_name = 1;
      if (verbose)
        cerr << "The net formal charge will be appended to the name\n";

      appends_to_be_done = 1;
    } else if ("CLND" == p) {
      append_clnd_count_to_name = 1;
      if (verbose)
        cerr << "The CLND Nitrogen count will be appended to the name\n";

      appends_to_be_done = 1;
    } else if ("PREPEND" == p) {
      do_appends_as_prepends = 1;
      if (verbose)
        cerr << "appends done as prepends\n";
    } else if (p.starts_with("LARGE")) {
      appended_properties_from_largest_fragment = 1;
      if (verbose)
        cerr << "Appended properties derived from largest fragment\n";
    } else if (p.starts_with("HETERO") || p.starts_with("HTROAC")) {
      append_heteratom_count_to_name = 1;
      if (verbose)
        cerr << "Will append the heteroatom count to the name\n";

      appends_to_be_done = 1;
    } else if ("help" == p) {
      DisplayPDirectives(0);
    } else {
      cerr << "Unrecognised -p qualifier '" << p << "'\n";
      DisplayPDirectives(1);
    }
  }

  return 1;
}

int
FileconvConfig::ParseMolecularWeightSpecifications(Command_Line& cl) {
  if (cl.option_present('w')) {
    IWString mtmp;
    int i = 0;
    while (cl.value('w', mtmp, i++)) {
      if ('*' == mtmp) {
        compute_molecular_weight_for_each = 1;
        if (verbose)
          cerr << "The molecular weight of each molecule will be reported\n";
      } else if ("LARGE" == mtmp) {
        compute_molecular_weight_based_on_largest_fragment = 1;
        if (verbose)
          cerr << "Molecular weights will be for the largest fragment\n";
      } else {
        double tmp;
        if (!mtmp.numeric_value(tmp) || tmp <= 0.0) {
          cerr << "The lower molecular weight cutoff must be a positive number '" << mtmp << "'\n";
          Usage(21);
        }

        lower_amw_cutoff = tmp;
        if (verbose)
          cerr << "Molecules with amw below " << lower_amw_cutoff << " will be excluded\n";
      }
    }

    atom_count.resize(300);
  }

  if (cl.option_present('W')) {
    IWString mtmp;
    int i = 0;
    while (cl.value('W', mtmp, i++)) {
      if ('*' == mtmp) {
        compute_molecular_weight_for_each = 1;
        if (verbose)
          cerr << "The molecular weight of each molecule will be reported\n";
      } else if ("LARGE" == mtmp) {
        compute_molecular_weight_based_on_largest_fragment = 1;
        if (verbose)
          cerr << "Molecular weights will be for the largest fragment\n";
      } else {
        double tmp;
        if (!mtmp.numeric_value(tmp) || tmp <= 0.0) {
          cerr << "The lower molecular weight cutoff must be a positive number '" << mtmp << "'\n";
          Usage(23);
        }
        if (lower_amw_cutoff && tmp < lower_amw_cutoff) {
          cerr << "The upper amw cutoff must be greater than the lower " << tmp << '\n';
          Usage(22);
        }

        upper_amw_cutoff = tmp;
        if (verbose)
          cerr << "Molecules with amw above " << upper_amw_cutoff << " will be excluded\n";
      }
    }
  }

  return 1;
}

int
FileconvConfig::ParseRingFiltering(Command_Line& cl) {
  if (cl.option_present('r')) {
    if (audit_input)
      cerr << "Warning, ring cutoff with -a does not do anything\n";

    const_IWSubstring r;
    for (int i = 0; cl.value('r', r, i); ++i) {
      if (r == "spiro" || r == "SPIRO") {
        ring_systems_include_spiro = 1;
        if (verbose) {
          cerr << "Ring systems span spiro fusions\n";
        }
      }  else if (r.starts_with("S:")) {
        r.remove_leading_chars(2);
        if (! r.numeric_value(min_ring_systems_required) || min_ring_systems_required < 0) {
          cerr << "The min ring systems allowed (-r S:) option must be a whole non negative "
                  "number\n";
          Usage(2);
        }
        if (verbose) {
          cerr << "Molecules with fewer than " << min_ring_systems_required <<
                  " will be discarded\n";
        }
      } else if (r.starts_with("aliph:")) {
        r.remove_leading_chars(6);
        if (! r.numeric_value(min_aliphatic_ring_count) || min_aliphatic_ring_count < 0) {
          cerr << "The min aliphatic ring 'aliph:' directive must be followed by a whole +ve integer\n";
          return 0;
        }
        if (verbose) {
          cerr << "Will reject molecules with fewer than " << min_aliphatic_ring_count << " aliphatic rings\n";
        }
      } else if (r.starts_with("arom:")) {
        r.remove_leading_chars(5);
        if (! r.numeric_value(min_aromatic_ring_count) || min_aromatic_ring_count < 0) {
          cerr << "The min aromatic ring 'arom:' directive must be followed by a whole +ve integer\n";
          return 0;
        }
        if (verbose) {
          cerr << "Will reject molecules with fewer than " << min_aromatic_ring_count << " aromatic rings\n";
        }
      } else if (!r.numeric_value(lower_ring_count_cutoff) || lower_ring_count_cutoff < 1) {
        cerr << "The lower ring count cutoff (-r) option must be a whole +ve number\n";
        Usage(2);
      } else if (verbose)
        cerr << "Molecules containing fewer than " << lower_ring_count_cutoff
             << " rings will be ignored\n";
    }
  }

  if (cl.option_present('R')) {
    if (audit_input)
      cerr << "Warning, ring cutoff with -a does not do anything\n";

    const_IWSubstring r;

    for (int i = 0; cl.value('R', r, i); i++) {
      if (r == "spiro" || r == "SPIRO") {
        ring_systems_include_spiro = 1;
        if (verbose) {
          cerr << "Ring systems span spiro fusions\n";
        }
      }  else if (r.starts_with("S:")) {
        r.remove_leading_chars(2);
        if (!r.numeric_value(max_ring_systems_allowed) || max_ring_systems_allowed < 0) {
          cerr << "The max ring systems allowed (-R S:) option must be a whole non negative "
                  "number\n";
          Usage(2);
        }

        if (verbose)
          cerr << "Will ignore molecules containing more than " << max_ring_systems_allowed
               << " ring systems\n";
      } else if (r.starts_with("MRS:")) {
        r.remove_leading_chars(4);
        if (!r.numeric_value(max_rings_in_a_ring_system) || max_rings_in_a_ring_system < 1) {
          cerr << "The Max Rings in a System (MRS:) directive must be followed by a whole +ve "
                  "number\n";
          return 3;
        }

        if (verbose)
          cerr << "Will discard any molecule with more than " << max_rings_in_a_ring_system
               << " rings in a system\n";
      } else if (r.starts_with("aliph:")) {
        r.remove_leading_chars(6);
        if (! r.numeric_value(max_aliphatic_ring_count) || max_aliphatic_ring_count < 0) {
          cerr << "The max aliphatic ring 'aliph:' directive must be followed by a whole +ve integer\n";
          return 0;
        }
        if (verbose) {
          cerr << "Will reject molecules with more than " << max_aliphatic_ring_count << " aliphatic rings\n";
        }
      } else if (r.starts_with("arom:")) {
        r.remove_leading_chars(5);
        if (! r.numeric_value(max_aromatic_ring_count) || max_aromatic_ring_count < 0) {
          cerr << "The max aromatic ring 'arom:' directive must be followed by a whole +ve integer\n";
          return 0;
        }
        if (verbose) {
          cerr << "Will reject molecules with more than " << max_aromatic_ring_count << " aromatic rings\n";
        }
      } else if (upper_ring_count_cutoff > 0) {
        cerr << "Already specified upper ring count cutoff (-R)\n";
        return 2;
      } else if (!r.numeric_value(upper_ring_count_cutoff) || upper_ring_count_cutoff < 1) {
        cerr << "The upper ring count cutoff (-R) option must be a whole +ve number\n";
        Usage(2);
      } else if (verbose)
        cerr << "Molecules containing more than " << upper_ring_count_cutoff
             << " rings will be ignored\n";
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', upper_ring_size_cutoff) || upper_ring_size_cutoff < 3) {
      cerr << "The upper ring size cutoff value (-m) must be a valid ring size\n";
      Usage(3);
    }

    if (verbose)
      cerr << "Will discard molecules having ring sizes > " << upper_ring_size_cutoff << '\n';
  }

  return 1;
}

int
FileconvConfig::ParseAtomCountDirectives(Command_Line& cl) {
  if (cl.option_present('c')) {
    const_IWSubstring c = cl.string_value('c');

    //  if (c.starts_with("LF:"))
    //  {
    //    c.remove_leading_chars(3);
    //    lower_atom_count_cutoff_applies_to_largest_fragment = 1;
    //  }

    if (!c.numeric_value(lower_atom_count_cutoff)) {
      cerr << "The lower atom count cutoff (-c) must be a whole +ve number\n";
      Usage(2);
    } else if (lower_atom_count_cutoff < 0) {
      cerr << "The lower atom count cutoff (-c) must be a whole +ve number\n";
      Usage(2);
    } else if (verbose)
      cerr << "Will exclude molecules with fewer than " << lower_atom_count_cutoff << " atoms\n";
  }

  if (cl.option_present('C')) {
    int i = 0;
    const_IWSubstring c;
    while (cl.value('C', c, i++)) {
      if ("implicit" == c) {
        include_implicit_hydrogens_in_upper_atom_count_comparison = 1;
        if (verbose)
          cerr << "Will include implicit Hydrogens in the upper atom count computation\n";
        continue;
      }

      //    if (c.starts_with("LF:"))
      //    {
      //      c.remove_leading_chars(3);
      //      upper_atom_count_cutoff_applies_to_largest_fragment = 1;
      //    }

      if (!c.numeric_value(upper_atom_count_cutoff) ||
          upper_atom_count_cutoff < lower_atom_count_cutoff) {
        cerr << "The upper atom count cutoff must be a whole positive number >= "
             << lower_atom_count_cutoff << '\n';
        cerr << "'" << c << "' is invalid\n";
        Usage(5);
      } else if (verbose)
        cerr << "Will exclude molecules with more than " << upper_atom_count_cutoff << " atoms\n";
    }
  }

  return 1;
}

int
FileconvConfig::ParseImplicitHydrogenDirectives(Command_Line& cl, char flag) {
  if (!cl.option_present(flag)) {
    return 1;
  }

  const_IWSubstring h;
  for (int i = 0; cl.value(flag, h, i); ++i) {
    if ('*' == h || "all" == h) {
      make_all_implicit_hydrogens_explicit = 1;
      if (verbose)
        cerr << "All implicit hydrogens make explicit\n";
    } else if (h.starts_with("last=")) {
      h.remove_leading_chars(5);
      if (! h.numeric_value(hydrogens_last)) {
        cerr << "Invalid hydrogens last directove '" << h << "'\n";
        return 0;
      }
      if (verbose) {
        cerr << "Will move elements with atomic number " << hydrogens_last <<
                " to the end of the connection table\n";
      }
    } else if (h == "help") {
      DisplayDashhOptions(cerr);
    } else {
      if (! process_cmdline_token('*', h, atoms_for_implicit_hydrogens, verbose)) {
        cerr << "Cannot process implicit hydrogen query qualifier '" << h << "'\n";
        return 0;
      }
    }
  }

  if (make_all_implicit_hydrogens_explicit && atoms_for_implicit_hydrogens.number_elements()) {
    cerr << "Cannot specify both '-h query' and '-h all' is impossible\n";
    Usage(72);
  }

  return 1;
}

int
FileconvConfig::ParseMkFragOptions(Command_Line& cl, char flag) {
  if (!cl.option_present(flag)) {
    return 1;
  }

  IWString d;
  for (int i = 0; cl.value(flag, d, i); ++i) {
    if (d == "def") {
      molecule_to_fragments = std::numeric_limits<int>::max();
    } else if (d == "help") {
    } else {
      cerr << "Unrecognised -" << flag << " specification '" << d << "'\n";
      return 0;
    }
  }

  return 1;
}

int
FileconvConfig::ReportResults(const Command_Line& cl, std::ostream& output) const {
  if (cl.number_elements() > 1)
    cerr << molecules_processed << " molecules processed, ";
  if (molecules_processed == 0) {
    return 1;
  }

  if (natoms_accumulator.n() > 0) {
    cerr << natoms_accumulator.n() << " molecules had between " << natoms_accumulator.minval()
         << " and " << natoms_accumulator.maxval() << " atoms.";
    if (natoms_accumulator.n() > 1)
      cerr << " Average " << natoms_accumulator.average() << '\n';
    else
      cerr << '\n';
  }

  if (verbose > 1) {
    for (int i = 0; i < atom_count.number_elements(); i++) {
      if (atom_count[i])
        cerr << atom_count[i] << " molecules had " << i << " atoms\n";
    }
  }

  if (audit_input) {
    return 1;
  }

  if (molecules_changed) {
    cerr << ' ' << molecules_changed << " molecules changed";
  }

  for (int i = 0; i < initial_fragment_count.number_elements(); i++) {
    if (initial_fragment_count[i]) {
      cerr << initial_fragment_count[i] << " molecules had " << i << " fragments\n";
    }
  }

  if (fragment_count && molecules_chopped) {
    cerr << molecules_chopped << " molecules chopped to ";
    if (reduce_to_largest_fragment)
      cerr << "largest";
    else
      cerr << fragment_count;
    cerr << " fragments\n";

    cerr << "Molecules lost between " << atoms_lost.minval() << " and " << atoms_lost.maxval()
         << " atoms\n";
  }

  if (output_organic_only && non_organic_molecules_found) {
    cerr << "Skipped " << non_organic_molecules_found
         << " molecules containing non organic atoms\n";
  }
  if (non_real_elements_found) {
    cerr << "Skipped " << non_real_elements_found
         << " molecules containing non-periodic table elements\n";
  }
  if (molecules_excluded_for_non_allowed_elements) {
    cerr << "Skipped " << molecules_excluded_for_non_allowed_elements
         << " molecules with non-allowed elements\n";
  }
  if (molecules_below_molecular_weight_cutoff) {
    cerr << "Skipped " << molecules_below_molecular_weight_cutoff
         << " molecules with molecular wieght below " << lower_molecular_weight_cutoff << '\n';
  }
  if (molecules_above_molecular_weight_cutoff) {
    cerr << "Skipped " << molecules_below_molecular_weight_cutoff
         << " molecules with molecular wieght above " << lower_molecular_weight_cutoff << '\n';
  }
  if (molecules_below_atom_count_cutoff) {
    cerr << "Skipped " << molecules_below_atom_count_cutoff
         << " molecules with atom count less than " << lower_atom_count_cutoff << '\n';
  }
  if (molecules_above_atom_count_cutoff) {
    cerr << "Skipped " << molecules_above_atom_count_cutoff
         << " molecules with atom count greater than " << upper_atom_count_cutoff << '\n';
  }
  if (molecules_with_too_many_components) {
    cerr << "Skipped " << molecules_with_too_many_components << " molecules with more than "
         << maximum_fragment_count << " components\n";
  }
  if (molecules_with_large_fragments) {
    cerr << "Skipped " << molecules_with_large_fragments
         << " with non-largest fragments with more than "
         << remove_molecules_with_non_largest_fragment_natoms << " atoms\n";
  }
  if (molecules_with_fragments_reduced_by_query) {
    cerr << molecules_with_fragments_reduced_by_query
         << " molecules chopped to fragment which matches query\n";
  }
  if (molecules_not_matching_fragment_queries) {
    cerr << "Skipped " << molecules_not_matching_fragment_queries
         << " molecules which did not match fragment query specifications\n";
  }
  if (remove_fragments_this_size_or_smaller) {
    cerr << molecules_with_very_small_fragments_removed << " molecules lost fragments with "
         << remove_fragments_this_size_or_smaller << " or fewer atoms\n";
  }
  if (remove_duplicate_fragments) {
    cerr << molecules_with_duplicate_fragments_removed
         << " molecules had duplicate fragments removed\n";
  }
  if (keep_all_organic_fragments) {
    cerr << molecules_with_nonorganic_fragments_removed
         << " molecules had nonorganic fragments removed\n";
  }
  if (discard_molecule_if_multiple_fragments_larger_than) {
    cerr << molecules_with_multiple_large_fragments
         << " molecules had multiple fragments with at least "
         << discard_molecule_if_multiple_fragments_larger_than << '\n';
  }
  if (mixture_if_largest_frags_differ_by >= 0) {
    cerr << molecules_declared_mixtures_by_atom_count_difference
         << " molecules declared mixtures by atom count difference\n";
  }
  if (molecules_with_too_few_rings) {
    cerr << "Skipped " << molecules_with_too_few_rings << " molecules having fewer than "
         << lower_ring_count_cutoff << " rings\n";
  }
  if (molecules_with_too_many_rings) {
    cerr << "Skipped " << molecules_with_too_many_rings << " molecules having more than "
         << upper_ring_count_cutoff << " rings\n";
  }
  if (molecules_with_large_rings) {
    cerr << "Skipped " << molecules_with_large_rings << " molecules containing rings larger than "
         << upper_ring_size_cutoff << " atoms\n";
  }
  if (molecules_with_ring_systems_too_large) {
    cerr << "Skipped " << molecules_with_ring_systems_too_large
         << " molecules containing ring systems with more than " << max_rings_in_a_ring_system
         << " rings\n";
  }
  if (molecules_with_too_few_ring_systems) {
    cerr << "Skipped " << molecules_with_too_few_ring_systems << " with too few ring systems\n";
  }
  if (molecules_with_too_many_ring_systems) {
    cerr << "Skipped " << molecules_with_too_many_ring_systems << " with too many ring systems\n";
  }
  if (max_path_length) {
    cerr << "Skipped " << molecules_with_longest_path_too_long << " molecules with max path above "
         << max_path_length << '\n';
  }

  if (amw_accumulator.n() > 1) {
    cerr << amw_accumulator.n() << " molecular weights. Between " << amw_accumulator.minval()
         << " and " << amw_accumulator.maxval() << " ave = " << amw_accumulator.average()
         << " var = " << amw_accumulator.variance() << '\n';
  }

  if (molecules_below_amw_cutoff) {
    cerr << "Excluded " << molecules_below_amw_cutoff << " molecules with amw below "
         << lower_amw_cutoff << '\n';
  }
  if (molecules_above_amw_cutoff) {
    cerr << "Excluded " << molecules_above_amw_cutoff << " molecules with amw above "
         << upper_amw_cutoff << '\n';
  }

  if (molecules_with_chiral_centres) {
    cerr << molecules_with_chiral_centres << " molecules had chiral atoms\n";
  }
  if (chiral_centres_inverted) {
    cerr << chiral_centres_inverted << " chiral centres inverted\n";
  }

  if (molecules_with_chiral_centres_removed_by_rmchiral) {
    cerr << "Removed " << chiral_centres_removed_by_rmchiral << " chiral centres from "
         << molecules_with_chiral_centres_removed_by_rmchiral << " molecules per rmchiral\n";
  }

  if (molecules_containing_isotopes) {
    cerr << molecules_containing_isotopes
         << " molecules contained isotopic atoms that were processed\n";
  }

  if (molecules_with_invalid_chiral_centres) {
    cerr << molecules_with_invalid_chiral_centres << " with invalid chiral centres\n";
  }

  elements_to_remove.report(cerr);

  if (element_transformations.number_elements()) {
    element_transformations.debug_print(cerr);
  }

  if (structure_fixing.active()){
    structure_fixing.report(cerr);
  }

  if (molecules_with_abnormal_valences) {
    cerr << molecules_with_abnormal_valences << " molecules containing abnormal valences\n";
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.report(cerr);
  }

  if (charge_assigner.active()) {
    charge_assigner.report(cerr);
  }

  if (make_all_implicit_hydrogens_explicit || atoms_for_implicit_hydrogens.number_elements()) {
    cerr << molecules_to_which_hydrogens_were_added
         << " molecules had implicit hydrogens made explicit\n";
  }

  if (aromatise_these_rings.number_elements()) {
    cerr << molecules_changed_by_aromatising_rings << " molecules made aromatic rings\n";
  }

  if (max_chiral_centres > 0) {
    cerr << molecules_with_too_many_chiral_centres << " molecules had more than "
         << max_chiral_centres << '\n';
    for (int i = 0; i < chiral_centre_count.number_elements(); ++i) {
      if (chiral_centre_count[i] > 0)
        cerr << chiral_centre_count[i] << " molecules had " << i << " chiral centres\n";
    }
  }
  if (! remove_atom.empty()) {
    cerr << "atoms were removed from " << molecules_with_removed_atoms << " molecules\n";
  }

  if (acc_longest_sep.n() > 0) {
    cerr << "longest separation N = " << acc_longest_sep.n() <<
       " btw " << acc_longest_sep.minval() << " and " <<
       acc_longest_sep.maxval() << " mean " << acc_longest_sep.average() << '\n';
  }

  if (min_aliphatic_ring_count > 0) {
    cerr << molecules_with_too_few_aliphatic_rings <<
           " molecules with fewer than " << min_aliphatic_ring_count <<
           " aliphatic rings skipped\n";
  }

  if (min_aromatic_ring_count > 0) {
    cerr << molecules_with_too_few_aromatic_rings <<
           " molecules with fewer than " << min_aromatic_ring_count <<
           " aromatic rings skipped\n";
  }

  if (max_aliphatic_ring_count > 0) {
    cerr << molecules_with_too_many_aliphatic_rings <<
           " molecules with more than " << max_aliphatic_ring_count <<
           " aliphatic rings skipped\n";
  }

  if (max_aromatic_ring_count > 0) {
    cerr << molecules_with_too_many_aromatic_rings <<
           " molecules with more than " << max_aromatic_ring_count <<
           " aromatic rings skipped\n";
  }

  if (name_rx || grep_v_name_rx) {
    cerr << molecules_discarded_for_name_mismatch << " molecules discarded for name regular expression mismatch\n";
  }

  return 1;
}

int
FileconvConfig::Build(Command_Line& cl) {
  verbose = cl.option_count('v');

  if (verbose) {
#if defined(GIT_HASH) && defined(TODAY)
    cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
    cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  }

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    return 0;
  }

  if (!process_elements(cl))
    Usage(2);

  if (!process_standard_smiles_options(cl, verbose)) {
    Usage(4);
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('N')) {
    if (!charge_assigner.construct_from_command_line(cl, verbose, 'N')) {
      cerr << "Cannot determine charge assigner from command line\n";
      Usage(77);
    }

    if (Daylight != global_aromaticity_type()) {
      set_global_aromaticity_type(Daylight);
      if (verbose)
        cerr << "Global aromaticity type set to Dayligth for Charge Assigner";
    }
  }

  if (cl.option_present('H')) {
    if (!donor_acceptor_assigner.construct_from_command_line(cl, 'H', verbose)) {
      cerr << "Cannot initialise donor/acceptor assigner (-H option)\n";
      Usage(6);
    }
  }

  if (cl.option_present('t')) {
    if (!element_transformations.construct_from_command_line(cl, verbose, 't'))
      Usage(8);
  }

  // Processing the -O option is hard because it might be either
  //   A set of elements
  //   A set of directives for allowable elements

  if (cl.option_present('O')) {
    if (!ParseOrganicSpecification(cl, 'O')) {
      DisplayDashOOptions('O', cerr);
      return 3;
    }
  }

  if (cl.option_present('e')) {
    exclude_non_real_elements = 1;

    set_auto_create_new_elements(1);

    if (verbose)
      cerr << "Molecules containing elements not in periodic table will be discarded\n";
  }

  if (cl.option_present('a')) {
    audit_input = 1;
    if (verbose)
      cerr << "No output, just audit input\n";
    if (cl.option_present('o') || cl.option_present('S')) {
      cerr << "The -a and -o/-S options are incompatible\n";
      Usage(5);
    }
  }

  if (!GetChiralityInstructions(cl, 's')) {
    return 1;
  }

  if (!GetChargeCalculation(cl, 'Q')) {
    return 1;
  }

  if (!GetFragmentSpecifications(cl)) {
    return 1;
  }

  if (cl.option_present('J')) {
    if (!structure_fixing.initialise(cl, 'J', verbose > 1)) {
      cerr << "Cannot initialise structure fixing (-J option)\n";
      Usage(4);
    }
  }

  if (cl.option_present('X')) {
    if (!elements_to_remove.construct_from_command_line(cl, verbose, 'X')) {
      cerr << "Cannot discern elements to remove from -X switch\n";
      Usage(18);
    }
  }

  // The -a option is special because it suppresses all filtering and output
  // This checking is very incomplete...

  if (!cl.option_present('a'))
    ;
  else if (cl.option_present('c') || cl.option_present('C') || cl.option_present('w') ||
           cl.option_present('W') || cl.option_present('I') || cl.option_present('s') ||
           cl.option_present('S')) {
    cerr << "Filtering and output options like -c, -C -I -s and -S don't make sense with the -a "
            "option\n";
    return 5;
  }

  if (!ParseMolecularWeightSpecifications(cl)) {
    return 1;
  }

  if (!ParseAtomCountDirectives(cl)) {
    return 1;
  }

  if (!number_assigner.initialise(cl, 'n', verbose)) {
    cerr << "Cannot process -n option\n";
    Usage(51);
  }

  if (!ParseRingFiltering(cl)) {
    return 1;
  }

  if (! ParseImplicitHydrogenDirectives(cl, 'h')) {
    return 1;
  }

  if (cl.option_present('V')) {
    skip_molecules_with_abnormal_valences = 1;
    if (verbose)
      cerr << "Molecules containing abnormal valences will be skipped\n";
  }

  if (!GetIsotopeDirectives(cl, 'I')) {
    return 1;
  }

  if (!GetTranslations(cl, 'T')) {
    return 1;
  }

#ifdef NOT_IMPLEMENTED_YET
  if (! ParseMkFragOptions(cl, 'D')) {
    return 1;
  }
#endif

  if (!ParseMiscOptions(cl, 'Y')) {
    return 1;
  }

  if (!GatherAppendSpecifications(cl, 'p')) {
    return 1;
  }

  return 1;
}

}  // namespace fileconv
