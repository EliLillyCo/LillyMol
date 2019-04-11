/*
  File conversion utility.
*/

#include <iostream>
#include <iostream>
#include <memory>
#include <limits>
#include <assert.h>


using std::cerr;
using std::cout;
using std::endl;

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "accumulator.h"
#include "iwstring.h"
#include "cmdline.h"
#include "iwqsort.h"
#include "misc.h"

#include "molecule.h"
#include "substructure.h"
#include "chiral_centre.h"
#include "smiles.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "charge_assigner.h"
#include "donor_acceptor.h"
#include "target.h"
#include "path.h"
#include "is_actually_chiral.h"
#include "fix_structures.h"
#include "do_remove_duplicate_fragments.h"
#include "istream_and_type.h"
#include "output.h"
#include "known_fragment_data.h"
#include "numass.h"
#include "rmele.h"
#include "etrans.h"
#include "mdl.h"
#include "allowed_elements.h"

static const char *prog_name;

namespace Fileconv 
{
static int first_call = 1;

/*
  The programme can operate in the mode of just scanning its input
  and not writing anything
*/

int audit_input = 0;

static int debug_print_each_molecule = 0;

static int print_bond_lengths = 0;
static int print_bond_angles = 0;
static int print_torsions = 0;

static int molecules_read = 0;
static int molecules_written = 0;

static int molecules_changed = 0;

/*
  We can append a string to the name of any molecule whose connection table is changed
*/

static IWString molecule_changed_string;

static Chemical_Standardisation chemical_standardisation;

static Charge_Assigner charge_assigner;

static Donor_Acceptor_Assigner donor_acceptor_assigner;

static Structure_Fixing structure_fixing;

static int fragment_count = 0;
static int reduce_to_largest_fragment = 0;
static int reduce_to_all_largest_fragments = 0;
static int reduce_to_largest_organic_fragment = 0;
static int reduce_to_largest_organic_fragment_carefully = 0;
static bool keep_all_organic_fragments = false;
static int molecules_with_nonorganic_fragments_removed = 0;
static int maximum_fragment_count = 0;
static int min_size_max_fragment_count = 0;
static int molecules_chopped = 0;
static int molecules_with_too_many_components = 0;
static int molecules_with_large_fragments = 0;
static int remove_duplicate_fragments = 0;
static int molecules_with_duplicate_fragments_removed = 0;
static int remove_fragments_this_size_or_smaller = 0;
static int molecules_with_very_small_fragments_removed = 0;
static int remove_largest_fragment = 0;
static int remove_molecules_with_non_largest_fragment_natoms = -1;
static int strip_to_n_largest_fragments = 0;
static int sort_by_fragment_size = 0;
static IWString tag_for_removed_fragments;

/*
  People are registering mixtures. This is very bad for us.
  If a molecule contains more than 1 fragment with more than
  a given number of atoms, discard the molecule
*/

static int discard_molecule_if_multiple_fragments_larger_than = 0;

static int molecules_with_multiple_large_fragments = 0;

/*
  Sometimes if things are ambiguous, just call it a mixture if the two
  largest fragments differ by a given amount or less
*/

static int mixture_if_largest_frags_differ_by = -1;

static int molecules_declared_mixtures_by_atom_count_difference = 0;

static Known_Fragment_Data known_fragment_data;

/*
  We can keep track of how many atoms are lost
*/

static Accumulator_Int<int> atoms_lost;
static extending_resizable_array<int> initial_fragment_count;

/*
  For efficiency we have a single variable which indicates whether or not
  the fragment code must be called.
*/

static int need_to_call_process_fragments = 0;

static int skip_molecules_with_abnormal_valences = 0;
static int ok_bad_valence_on_isotopically_labelled = 0;
static int molecules_with_abnormal_valences = 0;

static int max_path_length = 0;
static int molecules_with_longest_path_too_long = 0;

static int exclude_isotopes = 0;
static int molecules_containing_isotopes = 0;    // only if the isotopic attributes are changed

static int convert_isotopes = 0;

static int convert_all_isotopes_to = 0;

static int remove_isotopic_atoms = 0;

static int convert_isotopes_to_atom_map_numbers = 0;
static int convert_atom_map_numbers_to_isotopes = 0;

static resizable_array<int> convert_specific_isotopes;
static resizable_array<int> convert_specific_isotopes_new_isotope;

static resizable_array_p<Substructure_Query> convert_specific_isotopes_query;
static resizable_array<int> convert_specific_isotopes_query_new_isotope;

static int output_organic_only = 0;
static int non_organic_molecules_found = 0;

static int exclude_non_real_elements = 0;
static int non_real_elements_found = 0;

static int molecules_excluded_for_non_allowed_elements = 0;

static double lower_molecular_weight_cutoff = 0.0;
static double upper_molecular_weight_cutoff = 0.0;

static int molecules_below_molecular_weight_cutoff = 0;
static int molecules_above_molecular_weight_cutoff = 0;

static Accumulator<molecular_weight_t> amw_accumulator;
static Accumulator_Int<int> natoms_accumulator;
static extending_resizable_array<int> atom_count;

static int lower_atom_count_cutoff = 0;
static int upper_atom_count_cutoff = 0;

/*
  I never implemented these two, they are mostly taken care
  of by atom_count_includes_only_atoms_in_largest_fragment
  If that ever becomes a problem, go ahead and implement them

  static int lower_atom_count_cutoff_applies_to_largest_fragment = 0;
  static int upper_atom_count_cutoff_applies_to_largest_fragment = 0;
*/

static int include_implicit_hydrogens_in_upper_atom_count_comparison = 0;

static int atom_count_includes_only_atoms_in_largest_fragment = 0;

static int molecules_below_atom_count_cutoff = 0;
static int molecules_above_atom_count_cutoff = 0;

static IW_Regular_Expression name_rx;

static IW_Regular_Expression change_name_rx;
static int change_name_rx_must_match = 0;

/*
  With the -B option, one specify the number of connection table errors
  allowed before programme exit.
*/

static int connection_table_errors_allowed = 0;

static IWString connection_table_error_file;

static IWString output_file_stem;

/*
  When reducing fragment count, we can keep the smallest or largest
  fragment which matches a given query.
*/

static resizable_array_p<Substructure_Query> smallest_fragment_queries;
static resizable_array_p<Substructure_Query> largest_fragment_queries;
static resizable_array_p<Substructure_Query> keep_fragment_queries;
static resizable_array_p<Substructure_Query> remove_fragment_queries;

static int molecules_with_fragments_reduced_by_query = 0;
static int molecules_not_matching_fragment_queries = 0;

static int lower_ring_count_cutoff = 0;
static int molecules_with_too_few_rings = 0;
static int upper_ring_count_cutoff = -1;
static int molecules_with_too_many_rings = 0;

static int max_ring_systems_allowed = -1;
static int molecules_with_too_many_ring_systems = 0;

static int upper_ring_size_cutoff = 0;
static int molecules_with_large_rings = 0;

static int max_rings_in_a_ring_system = -1;
static int molecules_with_ring_systems_too_large = 0;

/*
  We can translate molecules if we wish
*/

static coord_t dx = 0.0;
static coord_t dy = 0.0;
static coord_t dz = 0.0;

/*
  We can speed things up by having a single variable that indicates whether a translation
  is needed
*/

static int translation_specified = 0;

/*
  What kind of partial charge do we want
*/


static int fileconv_partial_charge_type = 0;

#define FILECONV_PARTIAL_CHARGE_GASTEIGER 1
#define FILECONV_PARTIAL_CHARGE_GASTEIGER_HUCKEL 2
#define FILECONV_PARTIAL_CHARGE_ABRAHAM 3

/*
  We can optionally assign each molecule written a number R(number).
*/

static Number_Assigner number_assigner;

static Elements_to_Remove elements_to_remove;

/*
  We can transform element types.
*/

static Element_Transformations element_transformations;

// With the -H option, we make implicit hydrogens explicit

static int make_all_implicit_hydrogens_explicit = 0;

static resizable_array_p<Substructure_Query> atoms_for_implicit_hydrogens;

static int molecules_to_which_hydrogens_were_added = 0;

/*
  Various things for chirality and stereo centres
*/

static int find_all_chiral_centres = 0;
static int find_all_ring_chiral_centres = 0;
static int invert_all_chiral_centres = 0;
static int reflect_coordinates = 0;
static int chiral_centres_inverted = 0;
static int remove_invalid_chiral_centres = 0;
static int molecules_with_invalid_chiral_centres = 0;
static int molecules_with_chiral_centres = 0;
static int remove_chiral_data_from_all_molecules = 0;
static int remove_non_organic_chirality = 0;
static int max_chiral_centres = 0;
static int molecules_with_too_many_chiral_centres = 0;
static extending_resizable_array<int> chiral_centre_count;

static resizable_array_p<Substructure_Query> remove_chiral_centres_on;

static int chiral_centres_removed_by_rmchiral = 0;
static int molecules_with_chiral_centres_removed_by_rmchiral = 0;

static int remove_directional_bonds_from_input = 0;

static int remove_invalid_directional_bonds_from_input = 0;

static int change_double_bonds_between_permanent_aromatic_to_single = 0;

/*
  With the -O (organic) switch, we can specify a list of non-organic
  elements which are in fact OK
*/

static Set_of_Element_Matches ok_non_organics;

/*
  Dec 2008. I want to be able to reject anything that contains certain elements
  Enabled by '-O def'
*/

static int filter_for_disallowed_elements = 0;

static Allowed_Elements allowed_elements;

static std::ofstream reject_log;

static Molecule_Output_Object rejections_output_object;

/*
  As part of the smiles infrastructure, I wanted to simultaneously
  produce llyg.smi and llygO.smi. This allows that...
*/

static IWString_and_File_Descriptor stream_for_smiles_before_filters;

static int verbose = 0;
static int compute_molecular_weight_for_each = 0;

static int compute_molecular_weight_based_on_largest_fragment = 0;

static int appends_to_be_done = 0;

static int do_appends_as_prepends = 0;

static int appended_properties_from_largest_fragment = 0;

static int append_molecular_weight_to_name = 0;

static int append_molecular_weight_ok_isotope_to_name = 0;

static int append_exact_mass_to_name = 0;

static int append_heteratom_count_to_name = 0;

static int append_molecular_formula_to_name = 0;

static int append_isis_molecular_formula_to_name = 0;

static int append_nrings_to_name = 0;

static int append_aromatic_ring_count_to_name = 0;

static int append_natoms_to_name = 0;

static int append_nbonds_to_name = 0;

static int append_net_formal_charge_to_name = 0;

static int append_clnd_count_to_name = 0;

static molecular_weight_t lower_amw_cutoff = -1.0;
static int molecules_below_amw_cutoff = 0;
static molecular_weight_t upper_amw_cutoff = -1.0;
static int molecules_above_amw_cutoff = 0;

static IWString substitute_for_whitespace_in_name;
static int truncate_names_to_first_token = 0;

static char truncate_name_at_first = '\0';
static char truncate_name_at_last = '\0';

static IWString prepend_to_name;

static resizable_array_p<Substructure_Query> aromatise_these_rings;

static int molecules_changed_by_aromatising_rings = 0;

static int remove_unnecessary_square_brackets = 0;

static int remove_all_possible_square_brackets = 0;

static int reset_atom_map_numbers = 0;

/*
  If we have multiple invocations, we need to reset some things
*/

static void
reset_file_scope_variables()
{
  audit_input = 0;
  debug_print_each_molecule = 0;
  print_bond_lengths = 0;
  print_bond_angles = 0;
  print_torsions = 0;
  molecules_read = 0;
  molecules_written = 0;
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

  exclude_isotopes = 0;
  molecules_containing_isotopes = 0;

  convert_isotopes = 0;
  convert_specific_isotopes.resize(0);
  convert_specific_isotopes_new_isotope.resize(0);
  convert_all_isotopes_to = 0;

  convert_specific_isotopes_query.resize(0);
  convert_specific_isotopes_query_new_isotope.resize(0);

  remove_isotopic_atoms = 0;

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

//int lower_atom_count_cutoff_applies_to_largest_fragment = 0;
//int upper_atom_count_cutoff_applies_to_largest_fragment = 0;

  include_implicit_hydrogens_in_upper_atom_count_comparison = 0;

  atom_count_includes_only_atoms_in_largest_fragment = 0;

  molecules_below_atom_count_cutoff = 0;
  molecules_above_atom_count_cutoff = 0;

  name_rx.deactivate();

  connection_table_errors_allowed = 0;

  connection_table_error_file.resize(0);

  output_file_stem.resize(0);

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
 
  max_ring_systems_allowed = -1;
  molecules_with_too_many_ring_systems = 0;

  dx = 0.0;
  dy = 0.0;
  dz = 0.0;

  translation_specified = 0;

  fileconv_partial_charge_type = 0;

  number_assigner.deactivate();

  elements_to_remove.resize(0);

  element_transformations.resize(0);

  make_all_implicit_hydrogens_explicit = 0;

  atoms_for_implicit_hydrogens.resize(0);

  molecules_to_which_hydrogens_were_added = 0;

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

//reject_log;

//rejections_output_object;

  verbose = 0;
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

  IWString prepend_to_name = "";

  aromatise_these_rings.resize(0);

  molecules_changed_by_aromatising_rings = 0;

  remove_unnecessary_square_brackets = 0;

  remove_all_possible_square_brackets = 0;

  reset_atom_map_numbers = 0;

  convert_isotopes_to_atom_map_numbers = 0;
  convert_atom_map_numbers_to_isotopes = 0;

  return;
}

static void
usage(int rc = 1)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << "\n";
  cerr << "usage: " << prog_name << " -i <input type> -o <output type> file1 file2...\n";
  cerr << "The following options are recognised\n";
  cerr << "  -f <directive> fragment selection, enter '-f help' for details\n";
  cerr << "  -F <number>    exclude molecules having more than <number> fragments\n";
  cerr << "  -O none        \"organic\" molecules only. '-O def' removes all bad elements\n";
  cerr << "  -O <el>        \"organic\" molecules only, but 'el' is OK (repeat for each OK ele)\n";
  cerr << "  -E <symbol>    create element <symbol>, use 'autocreate' for all\n";
//cerr << "  -E autocreate  automatically create new elements when encountered\n";
//cerr << "  -w '*'         compute the molecular weight for each molecule\n";
//cerr << "  -w <number>    exclude molecules with molecular weight below <number>\n";
//cerr << "  -W <number>    exclude molecules with molecular weight above <number>\n";
  cerr << "  -w <amw, -W <amw>  specify lower(-w) and upper (-W) amw limits\n";
  cerr << "  -W LARGE       compute molecular weight based on largest fragment\n";
  cerr << "  -c <number>    exclude molecules with atom count below <number>\n";
  cerr << "  -C <number>    exclude molecules with atom count above <number>\n";
  cerr << "  -C implicit    when computing atom count for -C, include implicit Hydrogens\n";
  cerr << "  -r <n>         omit molecules with fewer than <n> rings\n";
  cerr << "  -R <n>         omit molecules with more  than <n> rings. Use S:n for ring systems\n";
  cerr << "  -m <ringsize>  discard molecules with rings larger than <ringsize>\n";
  cerr << "  -X <symbol>    extract/remove all atoms of type <symbol>. No bonds changed\n";
//cerr << "  -q <name>      only output molecules which   match the query in file <name>\n";
  cerr << "  -Q <type>      compute partial charges - enter '-Q help' for info\n";
  cerr << "  -I ...         isotope options, enter '-I help' for info'\n";
  cerr << "  -s <qual>      chirality options, enter '-s help' for details\n";
  cerr << "  -e             discard molecules with non periodic table elements\n";
//cerr << "  -G             create one file for each molecule (only one output type)\n";
  cerr << "  -n ...         number assigner options, enter '-n help' for more info\n";
  cerr << "  -o <type>      specify output file type(s), enter '-o help' for details\n";
  cerr << "  -a             audit input only, no output\n";
  cerr << "  -V             skip any molecule with abnormal valences\n";
  cerr << "  -T <dx,dy,dz>  translate molecules (dx, dy, dz)\n";
  cerr << "  -L <fname>     write rejected molecules to <fname>\n";
  cerr << "  -B ...         handling for otherwise fatal problems, enter '-B help'\n";
  cerr << "  -h <all>       make implicit hydrogens explicit\n";
  cerr << "  -h <query>     make implicit hydrogens on atoms matching <query> explicit\n";
#ifdef COMPILE_CHANGE_MOLECULE
  cerr << "  -P             invoke compiled-in molecule changer\n";
#endif
  cerr << "  -S <string>    create output files with name stem <string>\n";
  cerr << "  -p <string>    append various things to the name. Enter '-p help' for info\n";
  cerr << "  -i <type>      specify input file type. Enter '-i help' for details\n";
  (void) display_standard_charge_assigner_options(cerr, 'N');
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -H <..>        donor acceptor assignment, enter '-H help' for info\n";
//(void) display_standard_etrans_options(cerr, 't');
  cerr << "  -t E1=E2       standard element transformation options, enter '-t help'\n";
  (void) display_standard_aromaticity_options(cerr);
  (void) display_standard_smiles_options(cerr);
  cerr << "  -J <...>       fix obvious structure errors, enter '-J help' for info\n";
  cerr << "  -Y <...>       miscellaneous options, enter '-Y help' for info\n";
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static void
display_chirality_options()
{
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

static void
display_f_directives (int rc)
{
  cerr << "  -f <number>    remove fragments so that all written molecules\n";
  cerr << "                 have no more than <number> fragments\n";
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
  cerr << "  -f dmxt=<d>    discard molecules where largest fragments differ by <d> atoms or fewer\n";
  cerr << "  -f manlf=<d>   discard molecules that have a non-largest fragment with more than <d> atoms\n";
  cerr << "  -f klf=<d>     discard all but the <n> largest fragments\n";
  cerr << "  -f RMF=<tag>   when processing TDT forms, write removed fragments to <tag>\n";

  exit(rc);
}

static void
display_dash_F_options (std::ostream & os)
{
  os << "  -F <number>      discard molecules with more than <number> fragments\n";
  os << "  -F mnsz=<n>      when counting fragments only count those with > <n> atoms\n";

  exit(0);
}

static void
display_p_directives (int rc)
{
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

  exit(rc);
}

static void
display_dash_y_options (std::ostream & os,
                        char flag,
                        int rc)
{
  os << " -" << flag << " nbvm          No Bad Valence Messages, or strange electron counts\n";
  os << " -" << flag << " okbvi         during valence check, ok to have bad valence on isotopes\n";
  os << " -" << flag << " appchg=xxxx   append 'xxxx' to the name of molecules that are changed\n";
  os << " -" << flag << " pblen         print all bond lengths in the molecules\n";
  os << " -" << flag << " dbg           debug print each molecule\n";
  os << " -" << flag << " namerx=<rx>   discard molecules unless the molecule name matches <rx>\n";
  os << " -" << flag << " ftn           keep only the first token in molecule names\n";
  os << " -" << flag << " nsubw=c       translate all whitespace in molecule names to 'c'\n";
  os << " -" << flag << " chname=rx     change name to what is matched by <rx> (optional match)\n";
  os << " -" << flag << " CHNAME=rx     change name to what is matched by <rx> (MUST match)\n";
  os << " -" << flag << " tfirst=char   truncate name at first <char>\n";
  os << " -" << flag << " tlast=char    truncate name at last <char>\n";
  os << " -" << flag << " NPREPEND=s    prepend <s> to each name\n";
  os << " -" << flag << " maxplen=<n>   discard molecules with max path length > <n>\n";
//os << " -" << flag << " B4F=<fname>   write frag stripped smiles before filtering to <fname> \n";
  os << " -" << flag << " aclf          atom counts are for the largest fragment only\n";
  os << " -" << flag << " nhsqb         explicit hydrogen atoms in smiles written without square brackets\n";
  os << " -" << flag << " rmsqb         remove unnecessary square bracketed atoms  - Hcount is OK as specified\n";
  os << " -" << flag << " FHrmsqb       in atoms like [C] free the H count to what is computed. Square brackets removed\n";
  os << " -" << flag << " rmamap        remove atom map numbers\n";
  os << " -" << flag << " fixarom=smarts call find_kekule_form on the ring (system) matched by the first atom in <smarts>\n";
  os << " -" << flag << " help          this message\n";

  exit(rc);
}
static void
display_dash_I_options (char flag,
                        std::ostream & os)
{
  os << " -" << flag << " 0             discard molecules containing any isotopic atoms\n";
  os << " -" << flag << " change        change any isotopic atoms to normal form\n";
  os << " -" << flag << " change=<n>    change any isotope <n> atoms to normal form\n";
  os << " -" << flag << " change=<i,j>  change any isotope <i> atoms to isotope <j>\n";
  os << " -" << flag << " alliso=<i>    change any isotopic atoms to isotope <i>\n";
  os << " -" << flag << " smarts,<i>    change any atoms matching <smarts> to isotope <i>, e.g. '[2C],0' or '[#16],4'\n";
  os << " -" << flag << " remove        remove any isotopically labelled atoms\n";
  os << " -" << flag << " help          this message\n";

  exit(1);
}

static void
display_dash_O_options (char flag,
                        std::ostream & os)
{
  os << "Specifies elements allowed as ''organic'\n";
  os << "Operates in two different modes\n";
  os << "In the historical mode, fragment stripping is done, then remaining atoms\n";
  os << "checked for organic property\n";
  os << "to enable this use '-" << flag << " none, or -" << flag << " ele' to make elements organic\n";
  os << "The other mode is to first check the entire molecule for disallowed elements\n";
  os << "After fragment stripping is done, then check the largest fragment for non-organic elements\n";

  exit(1);
}

static void
display_dash_b_options (std::ostream & os)
{
  os << " -B <nn>        ignore as many as <nn> otherwise fatal input errors\n";
  os << " -B log=<fname> echo (unchanged) rejected input records to <fname>\n";

  exit(1);
}

static int
identify_ngroup (const Molecule & m,
                 atom_number_t n,
                 int * ngroup)
{
  ngroup[n] = n + 1;

  int rc = 1;

  const Atom * a = m.atomi(n);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(n, i);

    if (ngroup[j])
      continue;

    if (7 != m.atomic_number(j))
      continue;

    rc += identify_ngroup(m, j, ngroup);
  }

  return rc;
}

static int
is_azide (const Molecule & m,
          atom_number_t n1,
          int * ngroup)
{
  const Atom * an1 = m.atomi(n1);

  atom_number_t n2 = an1->other(n1, 0);

  const Atom * an2 = m.atomi(n2);

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

  const Bond * b12 = an1->item(0);
  const Bond * b23 = an2->item(0);
  if (n1 == b23->other(n2))
    b23 = an2->item(1);

  if (! b23->is_double_bond())
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

static int
compute_clnd (const Molecule & m)
{
  int rc = 0;

  int matoms = m.natoms();

  int * ngroup = new_int(matoms); std::unique_ptr<int[]> free_ngroup(ngroup);
  
// First look for [N-]#[N+]=N because they count for 1 each (why I have no idea).

  for (int i = 0; i < matoms; i++)  
  {
    if (ngroup[i])    // already counted as part of another group
      continue;

    if (7 != m.atomic_number(i))
      continue;

    if (1 != m.ncon(i))
      continue;

    if (is_azide(m, i, ngroup))
      rc += 1;
  }

  for (int i = 0; i < matoms; i++)
  {
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

static void
do_appends (Molecule & m, IWString & extra_stuff)
{
  if (append_molecular_formula_to_name)
  {
    extra_stuff += " MF = ";
    extra_stuff += m.molecular_formula();
  }

  if (append_isis_molecular_formula_to_name)
  {
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

  if (append_aromatic_ring_count_to_name)
  {
    m.compute_aromaticity_if_needed();

    int nr = m.nrings();
    int nar = 0;
    for (int i = 0; i < nr; i++)
    {
      const Ring * ri = m.ringi(i);
      if (ri->is_aromatic())
        nar++;
    }

    extra_stuff << " AROMRING = " << nar;
  }

  if (append_molecular_weight_to_name)
    extra_stuff << " AMW = " << m.molecular_weight();
  else if (append_molecular_weight_ok_isotope_to_name)
  {
    Molecular_Weight_Control mwc;
    Molecular_Weight_Calculation_Result mwcr;
    mwc.set_ignore_isotopes(0);
    (void) m.molecular_weight(mwc, mwcr);
    extra_stuff << " AMW = " << mwcr.amw();
  }

  if (append_exact_mass_to_name)
  {
    exact_mass_t x;
    if (! m.exact_mass(x))
      cerr << "Warning, molecule '" << m.name() << "' partial result for exact mass\n";

    extra_stuff << " EXACT_MASS = " << x;
  }

  if (append_heteratom_count_to_name)
  {
    extra_stuff << " HTROAC = " << (m.natoms() - m.natoms(6) - m.natoms(1));
  }

  if (append_net_formal_charge_to_name)
    extra_stuff << " FORMAL_CHARGE = " << m.formal_charge();

  if (append_clnd_count_to_name)
    extra_stuff << " CLND = " << compute_clnd(m);

  return;
}

static void
do_appends (Molecule & m)
{
  const IWString & mname = m.name();

  IWString extra_stuff;
  extra_stuff.resize(mname.length() + 200);

  if (! do_appends_as_prepends)    // we are appending, existing name at the beginning
    extra_stuff << mname;

  if (appended_properties_from_largest_fragment && m.number_fragments() > 1)
  {
    Molecule tmp;
    tmp.add_molecule(&m);    // make a copy
    tmp.reduce_to_largest_fragment();

    do_appends(tmp, extra_stuff);
  }
  else
    do_appends(m, extra_stuff);

  if (do_appends_as_prepends)     // name goes at the end
    extra_stuff << ' ' << mname;

  m.set_name(extra_stuff);

  return;
}

static int
do_remove_fragments_this_size_or_smaller (Molecule & m)
{
  int nf = m.number_fragments();

  int matoms = m.natoms();

  Set_of_Atoms atoms_to_be_removed;
  atoms_to_be_removed.resize(matoms);

  int number_fragments_removed = 0;

  for (int i = 0; i < nf; i++)
  {
    int aif = m.atoms_in_fragment(i);

    if (aif <= remove_fragments_this_size_or_smaller)
    {
      number_fragments_removed++;

      for (int j = 0; j < matoms; j++)
      {
        int f = m.fragment_membership(j);
        if (i == f)
          atoms_to_be_removed.add(j);
      }
    }
  }

  if (0 == atoms_to_be_removed.number_elements())
    return 1;

  if (matoms == atoms_to_be_removed.number_elements())
  {
    cerr << m.name() << " all fragments have " << remove_fragments_this_size_or_smaller << " or fewer atoms\n";
  }
  else
  {
    if (verbose > 1)
      cerr << "Removed " << number_fragments_removed << " fragments and " << atoms_to_be_removed.number_elements() << " atoms from '" << m.name() << "'\n";

    molecules_with_very_small_fragments_removed++;

    return m.remove_atoms(atoms_to_be_removed);
  }

  return 1;
}

static int
do_compute_partial_charges (Molecule & m)
{
  int rc = 0;

  if (FILECONV_PARTIAL_CHARGE_GASTEIGER == fileconv_partial_charge_type)
    rc = m.compute_Gasteiger_partial_charges();
  else if (FILECONV_PARTIAL_CHARGE_GASTEIGER_HUCKEL == fileconv_partial_charge_type)
    rc = m.compute_Gasteiger_Huckel_partial_charges();
  else if (FILECONV_PARTIAL_CHARGE_ABRAHAM == fileconv_partial_charge_type)
    rc = m.compute_Abraham_partial_charges();
  else
  {
    cerr << "Huh, what kind of partial charges do you want " << fileconv_partial_charge_type << endl;
    return 0;
  }

  return rc;
}

static void
append_bond(const Bond & b,
            std::ostream & output)
{
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

static int
do_print_torsions (const Molecule & m,
                   const Bond & b,
                   std::ostream & output)
{
  const atom_number_t a2 = b.a1();
  const atom_number_t a3 = b.a2();

  const Atom * aa2 = m.atomi(a2);
  const Atom * aa3 = m.atomi(a3);

  const int a2con = aa2->ncon();

  if (1 == a2con)
    return output.good();

  const int a3con = aa3->ncon();

  if (1 == a3con)
    return output.good();

  IWString a2smarts = m.smarts_equivalent_for_atom(a2);
  IWString a3smarts = m.smarts_equivalent_for_atom(a3);

  char sep = ' ';

  for (int i = 0; i < a2con; i++)
  {
    const Bond * b21 = aa2->item(i);

    const atom_number_t a1 = b21->other(a2);

    if (a1 == a3)
      continue;

    const IWString a1smarts = m.smarts_equivalent_for_atom(a1);

    for (int k = 0; k < a3con; k++)
    {
      const Bond * b34 = aa3->item(k);

      const atom_number_t a4 = b34->other(a3);

      if (a4 == a2)
        continue;

      const IWString a4smarts = m.smarts_equivalent_for_atom(a4);

      const angle_t tmp = m.dihedral_angle(a1, a2, a3, a4);
      output << m.name() << sep << "dihedral " << a1 << sep << a2 << sep << a3 << sep << a4 << sep << tmp << " (" << (tmp * RAD2DEG) << " deg)";
      output << sep << m.atomic_number(a1) << sep << a1smarts << sep;
      append_bond(*b21, output);
      output << sep << m.atomic_number(a2) << sep << a2smarts << sep;
      append_bond(b, output);
      output << sep << m.atomic_number(a3) << sep << a3smarts << sep;
      append_bond(*b34, output);
      output << sep << m.atomic_number(a4) << sep << a4smarts;

      output << "\n";
    }
  }


  return output.good();
}

static int
do_print_torsions(Molecule & m,
                  std::ostream & output)
{
  m.compute_aromaticity_if_needed();

  if (m.highest_coordinate_dimensionality() < 3)
    cerr << "do_print_torsions: WARNING not 3D " << m.name() << endl;

  const int nb = m.nedges();

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);

    if (! do_print_torsions(m, *b, output))
      return 0;
  }

  return output.good();
}

static int
do_print_bond_angles (const Molecule & m,
                      const atom_number_t a1,
                      const atom_number_t a2,
                      std::ostream & output)
{
  const Atom * aa2 = m.atomi(a2);

  const int ncon2 = aa2->ncon();

  for (int i = 0; i < ncon2; i++)
  {
    const atom_number_t a3 = aa2->other(a2, i);

    if (a3 == a1)
      continue;

    if (a3 < a1)
      continue;

    const angle_t theta = m.bond_angle(a1, a2, a3);
    output << " atoms " << a1 << " " << a2 << " " << a3 << " angle " << theta << " (" << (theta * RAD2DEG) << " deg)\n";
  }

  return output.good();
}

static int
do_print_bond_angles (const Molecule & m,
                      std::ostream & output)
{
  const int nb = m.nedges();

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    do_print_bond_angles(m, a1, a2, output);
    do_print_bond_angles(m, a2, a1, output);
  }

  return output.good();
}

static int
do_print_bond_lengths (const Molecule & m,
                       std::ostream & output)
{
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = m.atomi(i);
    const int icon = ai->ncon();
    for (int j = 0; j < icon; j++)
    {
      const Bond * b = ai->item(j);

      const atom_number_t k = b->other(i);

      const Atom * ak = m.atomi(k);
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

      output << " bond to " << k << " (" << ak->atomic_symbol() << ") dist = " <<
              ai->distance(*ak) << endl;
    }
  }

  return output.good();
}

static int
do_invert_all_chiral_centres (Molecule & m)
{
  const int nc = m.chiral_centres();
  if (0 == nc)
    return 1;

  int rc = 0;

  for (int i = 0; i < nc; i++)
  {
    Chiral_Centre * c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);

    c->invert();
    if (0 == rc)
      molecules_with_chiral_centres++;
  }

  chiral_centres_inverted += nc;

  return nc;
}

static int
do_reflect_coordinates (Molecule & m)
{
  if (m.highest_coordinate_dimensionality() < 3)
  {
    cerr << "do_reflect_coordinates:inadequate dimensionality " << m.highest_coordinate_dimensionality() << ", '" << m.name() << "'\n";
    return 0;
  }

//coord_t xmin, xmax, ymin, ymax, zmin, zmax;
//m.spatial_extremeties (xmin, xmax, ymin, ymax, zmin, zmax);

  const auto matoms = m.natoms();

//m.translate_atoms(0.0, 0.0, -zmin);    // ensure all on one side of xy plane

  for (auto i = 0; i < matoms; ++i)
  {
    Coordinates c;
    m.get_coords(i, c);
    m.setxyz(i, -c.x(), -c.y(), -c.z());
  }

//m.translate_atoms(0.0, 0.0, zmin);

  return do_invert_all_chiral_centres(m);
}

static int
do_find_all_chiral_centres (Molecule & m)
{
  int matoms = m.natoms();

  const int * symmetry = m.symmetry_classes();

  int rc = 1;
  for (int i = 0; i < matoms; i++)
  {
    if (NULL != m.chiral_centre_at_atom(i))
      continue;

    if (find_all_ring_chiral_centres && m.is_non_ring_atom(i))
      continue;

    Atom * a = const_cast<Atom *>(m.atomi(i));

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
    for (int j = 0; j < a->ncon(); j++)
    {
      atom_number_t k = a->other(i, j);
      if (0 == symmetries_found.add_if_not_already_present(symmetry[k]))
      {
        different_symmetries = 0;
        break;
      }
    }

    if (! different_symmetries)
      continue;

    if (verbose > 1)
      cerr << "Placing chiral centre at atom " << i << " (" << a->atomic_symbol() << ") " <<
              a->ncon() << " connections\n";

    if (NULL == m.create_chiral_centre(i))
    {
      cerr << "Yipes, could not place chiral center atom atom " << i << endl;
      rc = 0;
    }
  }

  return rc;
}

static int
do_remove_isotopes(Molecule & m)
{
  int rc = 0;

  for (int i = m.natoms() - 1; i >= 0; --i)
  {
    if (0 == m.isotope(i))
      continue;

    m.remove_atom(i);
    rc++;
  }

  if (rc)
    molecules_containing_isotopes++;

  return rc;
}

static int
do_convert_isotopes_to_atom_map_numbers(Molecule & m)
{
  int rc = 0;

  const int matoms = 0;

  for (int i = 0; i < matoms; ++i)
  {
    const int iso = m.isotope(i);

    if (0 == iso)
      continue;

    m.set_atom_map_number(i, iso);
    rc++;
  }

  return rc;
}

static int
do_convert_atom_map_numbers_to_isotopes(Molecule & m)
{
  int rc = 0;

  const int matoms = 0;

  for (int i = 0; i < matoms; ++i)
  {
    const int amap = m.atom_map_number(i);

    if (0 == amap)
      continue;

    m.set_isotope(i, amap);
    rc++;
  }

  return rc;
}


#ifdef OLD_ISOTOPE_STUFF

Moved the change isotopes to the beginning of processing to avoid this issue

/*
  An unfortunate problem with the variable molecules_containing_isotopes.
  It only gets incremented for this molecule in the upper level functions
  if we are excluding isotopes.

  Deal with the problem of isotopic hydrogens when chemical standardisation
  present
*/

static int
do_convert_isotopes(Molecule & m)
{
  int need_to_run_chemical_standardisation = 0;

  if (chemical_standardisation.active())
  {
    int matoms = m.natoms();

    for (int i = 0; i < matoms; i++)
    {
      const Atom * a = m.atomi(i);

      if (0 == a->isotope())
        continue;

      if (1 == a->atomic_number())
        need_to_run_chemical_standardisation = 1;
    }
  }

  int rc = m.transform_to_non_isotopic_form();

  if (rc)
  {
    molecules_containing_isotopes++;
    isotopes_converted += rc;
    if (need_to_run_chemical_standardisation)
      chemical_standardisation.process(m);
  }

  return rc;
}

#endif

static int
do_make_implicit_hydrogens_explicit (Molecule & m)
{
  int rc = 0;          // the number of atoms to which we add a Hydrogen

  Molecule_to_Match target(&m);

  int nq = atoms_for_implicit_hydrogens.number_elements();
  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = atoms_for_implicit_hydrogens[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      atom_number_t a = e->item(0);      // we process only the first matched atom

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

static int
identify_matched_atoms_with_chiral_centres (const Molecule & m,
                                const Set_of_Atoms & e,
                                Set_of_Atoms & atoms_with_chiral_centres_to_be_removed)
{
  int n = e.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = e[i];

    if (NULL == m.chiral_centre_at_atom(j))    // atom J does not have a chiral centre
      continue;
      
    if (atoms_with_chiral_centres_to_be_removed.contains(j))
      continue;

    atoms_with_chiral_centres_to_be_removed.add(j);
    rc++;
  }

  return rc;
}

static int
do_remove_chiral_centres_on_matched_atoms (Molecule & m,
                                           const resizable_array_p<Substructure_Query> & remove_chiral_centres_on)
{
  int nchiral = m.chiral_centres();

  if (0 == nchiral)
    return 0;

  Molecule_to_Match target(&m);

  int nq = remove_chiral_centres_on.number_elements();

  Set_of_Atoms atoms_with_chiral_centres_to_be_removed;

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = remove_chiral_centres_on[i]->substructure_search(target, sresults);

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      identify_matched_atoms_with_chiral_centres(m, *e, atoms_with_chiral_centres_to_be_removed);
    }
  }

  int rc = atoms_with_chiral_centres_to_be_removed.number_elements();

  if (0 == rc)
    return 0;

  molecules_with_chiral_centres_removed_by_rmchiral++;
  chiral_centres_removed_by_rmchiral += rc;

  for (int i = 0; i < rc; i++)
  {
    atom_number_t j = atoms_with_chiral_centres_to_be_removed[i];

    m.remove_chiral_centre_at_atom(j);
  }

  return rc;
}

static int
do_remove_non_organic_chirality (Molecule & m)
{
  int rc = 0;

  for (int i = m.chiral_centres() - 1; i >= 0; --i)
  {
    const Chiral_Centre * c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);

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

static int
do_aromatise_these_rings (Molecule & m,
                          int * procecss_these_atoms)
{
  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i)
  {
    if (0 == procecss_these_atoms[i])
      continue;

    m.set_implicit_hydrogens_known(i, 0);

    const Atom * a = m.atomi(i);

    const auto acon = a->ncon();

    for (int j = 0; j < acon; ++j)
    {
      const Bond * b = a->item(j);

      if (b->is_single_bond())
        continue;

      atom_number_t k = a->other(i, j);
      if (procecss_these_atoms[k])
        m.set_bond_type_between_atoms(i, k, SINGLE_BOND);
    }

    if (m.formal_charge(i))
    {
      m.set_formal_charge(i, 0);
      m.recompute_implicit_hydrogens(i);
    }
  }

  if (m.find_kekule_form(procecss_these_atoms))
    return 1;

  cerr << "Warning, did not find Kekule form for '" << m.name() << "', bonding has been destroyed!\n";
  return 0;
}

static int
do_aromatise_these_rings (Molecule & m,
                          int * ring_labels,
                          int * procecss_these_atoms,
                          const Set_of_Atoms & e)
{
  if (! m.is_ring_atom(e[0]))
    return 0;

  const auto matoms = m.natoms();

  const auto s = ring_labels[e[0]];

  for (auto i = 0; i < matoms; ++i)
  {
    if (s == ring_labels[i])
      procecss_these_atoms[i] = 1;
    else
      procecss_these_atoms[i] = 0;
  }

  return do_aromatise_these_rings(m, procecss_these_atoms);
}

/*
  We need two arrays. One that holds the ring infromation for the molecule.
  The other is what gets passed to find_kekule_form.
  For efficiency, we just allocate one array and pass pieces of it to the other functions
*/

static int
do_aromatise_these_rings (Molecule & m,
                          resizable_array_p<Substructure_Query> & aromatise_these_rings)
{
  int rc = 0;

  Molecule_to_Match target(&m);

  const int matoms = m.natoms();

  int * tmp = new_int(matoms + matoms); std::unique_ptr<int[]> free_tmp(tmp);

  m.label_atoms_by_ring_system(tmp);

  for (int i = 0; i < aromatise_these_rings.number_elements(); ++i)
  {
    Substructure_Results sresults;

    int nhits = aromatise_these_rings[i]->substructure_search(target, sresults);

    for (auto j = 0; j < nhits; j++)
    {
      if (do_aromatise_these_rings(m, tmp, tmp + matoms, *(sresults.embedding(i))))
        rc++;
    }
  }

  if (rc)
    molecules_changed_by_aromatising_rings++;

  return rc;
}

static int
extract_smaller_fragments_into_name_tag (Molecule & m,
                                         const IWString & tag_for_removed_fragments)
{
  const int matoms = m.natoms();

  Set_of_Atoms to_be_removed;
  int notused;

  (void) m.identify_largest_organic_fragment(to_be_removed, notused);

  int * subset = new_int(matoms); std::unique_ptr<int[]> free_subset(subset);

  to_be_removed.set_vector(subset, 1);

  Molecule f;
  m.create_subset(f, subset);

  m.remove_atoms(to_be_removed);

  IWString tmp(m.name());

  tmp << ' ' << tag_for_removed_fragments << ':' << f.smiles();

  m.set_name(tmp);

  return 1;
}


static void
do_sort_by_fragment_size (Molecule & m)
{
  int nf = m.number_fragments();

  if (1 == nf)
    return;

  resizable_array<int> fragment_size;

  for (int i = 0; i < nf; i++)
  {
    int a = m.atoms_in_fragment(i);

    fragment_size.add_if_not_already_present(a);
  }

  int n = fragment_size.number_elements();

  if (1 == n)   // all frags the same size
    return;

  Int_Comparator_Smaller ics;

  fragment_size.iwqsort(ics);

#ifdef DEBUG_DO_SORT_BY_FRAGMENT_SIZE
  cerr << n << " different fragment sizes\n";
  for (int i = 0; i < n; i++)
  {
    cerr << fragment_size[i] << " atoms\n";
  }
#endif

// Figure out the order in which existing atoms will be placed in the final molecule

  const int matoms = m.natoms();

  resizable_array<int> new_atom_order(matoms);

  for (int i = 0; i < n; i++)
  {
    const int a1 = fragment_size[i];

    for (int j = 0; j < matoms; j++)
    {
      int f = m.fragment_membership(j);
      if (a1 == m.atoms_in_fragment(f))
        new_atom_order.add(j);
    }
  }

#ifdef DEBUG_DO_SORT_BY_FRAGMENT_SIZE
  for (int i = 0; i < matoms; i++)
  {
    cerr << "i = " << i << " new_atom_order " << new_atom_order[i] << endl;
  }
#endif

  assert(new_atom_order.number_elements() == matoms);

// Each atom knows it's starting atom number

  int * tmp = new int[matoms]; std::unique_ptr<int[]> free_tmp(tmp);

  for (int i = 0; i < matoms; i++)
  {
    tmp[i] = i;
    m.set_user_specified_atom_void_ptr(i, tmp + i);
  }

// This is really hairy.... Make a swap, then keep swapping atoms
// until the destination gets it's final atom in place

  for (int i = 0; i < matoms; i++)
  {
    int j = new_atom_order[i];

    atom_number_t zatom = i;

    while (1)
    {
      const int * u = reinterpret_cast<const int *>(m.user_specified_atom_void_ptr(zatom));

//    cerr << "zatom " << zatom << " compare " << *u << " and " << j << endl;
      if (*u == j)
        break;

//    cerr << "i = " << i << " swap " << zatom << " and " << j << endl;

      m.swap_atoms(zatom, j);

      for (int k = i + 1; k < matoms; k++)
      {
        if (new_atom_order[k] == *u)
        {
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

static int
do_remove_by_fragment (Molecule & m,
                       const resizable_array<int> & fragments_to_remove)
{
  if (0 == fragments_to_remove.number_elements())
    return 1;

  Set_of_Atoms atoms_to_remove;

  for (int i = 0; i < fragments_to_remove.number_elements(); i++)
  {
    m.add_atoms_in_fragment(atoms_to_remove, fragments_to_remove[i]);
  }

  return m.remove_atoms(atoms_to_remove);
}

static int
do_remove_fragments_that_match_query (Molecule & m,
                                      resizable_array_p<Substructure_Query> & queries)
{
  int nq = queries.number_elements();

  Molecule_to_Match target(&m);

  resizable_array<int> fragments_to_be_removed;

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    nhits = sresults.number_embeddings();

    if (0 == nhits)
    {
      cerr << "Fragment selection query produced hits but no embeddings\n";
      cerr << "If compound query 'smarts1&&0smarts2', try reordering '0smarts2&&smarts1'\n";
    }

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      atom_number_t a0 = e->item(0);

      int f = m.fragment_membership(a0);

      fragments_to_be_removed.add_if_not_already_present(f);
    }
  }

  if (0 == fragments_to_be_removed.number_elements())
    return 0;

  if (m.number_fragments() == fragments_to_be_removed.number_elements())
  {
    if (verbose > 1)
      cerr << "All fragments match remove queries, ignored\n";
    return 0;
  }

  return do_remove_by_fragment(m, fragments_to_be_removed);
}

static int
do_reduce_to_all_largest_fragments (Molecule & m)
{
  resizable_array_p<Molecule> frags;

  m.create_components(frags);

  int nf = frags.number_elements();

  int atoms_in_largest_frag = frags[0]->natoms();

  for (int i = 1; i < nf; i++)
  {
    int m = frags[i]->natoms();

    if (m > atoms_in_largest_frag)
      atoms_in_largest_frag = m;
  }
  
  resizable_array<int> fragments_to_remove;

  for (int i = 0; i < nf; i++)
  {
    int m = frags[i]->natoms();

    if (m != atoms_in_largest_frag)
      fragments_to_remove.add(i);
  }

  return do_remove_by_fragment(m, fragments_to_remove);
}

/*
  remove the remove_largest_fragment largest fragments
*/

static int
do_remove_largest_fragments (Molecule & m, int remove_largest_fragment)
{
  int nf = m.number_fragments();

  if (nf >= remove_largest_fragment)
    return 0;

  int * atoms_in_fragment = new_int(nf); std::unique_ptr<int[]> free_atoms_in_fragment(atoms_in_fragment);

  for (int i = 0; i < nf; i++)
  {
    atoms_in_fragment[i] = m.atoms_in_fragment(i);
  }

  Int_Comparator_Larger icl;

  iwqsort(atoms_in_fragment, nf, icl);   // now sorted in increasing order

  int threshold_for_removal = atoms_in_fragment[nf - remove_largest_fragment];

// The first remove_largest_fragment fragments will be removed

  int * fragment_being_removed = new_int(nf); std::unique_ptr<int[]> free_fragment_being_removed(fragment_being_removed);

  for (int i = 0; i < nf; i++)
  {
    int aif = m.atoms_in_fragment(i);

    if (aif < threshold_for_removal)
      continue;

    fragment_being_removed[i] = 1;
    remove_largest_fragment--;
    if (0 == remove_largest_fragment)
      break;
  }

  assert(0 == remove_largest_fragment);

  int matoms = m.natoms();

  int * atoms_to_be_removed = new_int(matoms); std::unique_ptr<int[]> free_atoms_to_be_removed(atoms_to_be_removed);
  
  for (int i = 0; i < matoms; i++)
  {
    int f = m.fragment_membership(i);

    if (fragment_being_removed[f])
      atoms_to_be_removed[i] = 1;
  }

  return m.remove_atoms(atoms_to_be_removed);
}

/*template void iwqsort<int, Int_Comparator_Larger>(int*, int, Int_Comparator_Larger&, void*);
template void iwqsort<int>(int*, int, Int_Comparator_Larger &);
template void swap_elements<int>(int&, int&, void*);
template void move_in_from_left<int, Int_Comparator_Larger>(int*, int&, int&, int, Int_Comparator_Larger&, void*);
template void move_in_from_right<int, Int_Comparator_Larger>(int*, int&, int&, Int_Comparator_Larger&);
template void compare_two_items<int, Int_Comparator_Larger>(int*, Int_Comparator_Larger&, void*);*/

static int
too_many_large_fragments (Molecule & m,
                          int min_size_max_fragment_count, 
                          int maximum_fragment_count)
{
  int nf = m.number_fragments();

  int rc = 0;

  for (int i = 0; i < nf; i++)
  {
    int aif = m.atoms_in_fragment(i);

    if (aif >= min_size_max_fragment_count)
      rc++;
  }

  if (rc > maximum_fragment_count)
    return 1;

  return 0;
}

static int
do_remove_largest_fragment(Molecule & m)
{
  int nf = m.number_fragments();

  if (nf >= remove_largest_fragment)
    return 0;

  int atoms_in_largest_fragment = m.atoms_in_fragment(0);
  int largest_fragment = 0;

  for (int i = 1; i < nf; i++)
  {
    int aif = m.atoms_in_fragment(i);

    if (aif > atoms_in_largest_fragment)
    {
      atoms_in_largest_fragment = aif;
      largest_fragment = i;
    }
  }

  int matoms = m.natoms();

  int * atoms_to_be_removed = new_int(matoms); std::unique_ptr<int[]> free_atoms_to_be_removed(atoms_to_be_removed);

  for (int i = 0; i < matoms; i++)
  {
    if (m.fragment_membership(i) == largest_fragment)
      atoms_to_be_removed[i] = 1;
  }

  return m.remove_atoms(atoms_to_be_removed);;
}

static void
_do_remove_nonorganic_fragments (Molecule & m,
                                const int * frag_membership,
                                int & fragments_removed)
{
  int nf = m.number_fragments ();
  int matoms = m.natoms ();

  resizable_array<int> fragments_to_remove;

  for (int i = 0; i < nf; i++)
  {
    int aif = m.atoms_in_fragment (i);

    Molecule tmp;
    m.create_subset (tmp, frag_membership, i);

    bool delete_fragment = false;
    for (int j = 0; j < aif; j++) {
      const Element * e = tmp.elementi(j);
      if (e == NULL) {
        cerr << "NULL pointer found for element " << j << " in molecule " << m.name() << endl;
        cerr << "This should not happen, please contact c3tk" << endl;
        return;
      }
      if (e->organic()) {
        continue;
      }
      if (! ok_non_organics.matches(e)) {
        delete_fragment = true;
        break;
      }
      if (! e->is_in_periodic_table()) {
        delete_fragment = true;
        break;
      }
    }

    if (delete_fragment) {
      if (verbose)
        cerr << "Removing fragment " << i << " " << tmp.smiles() << endl;

      fragments_to_remove.add(i);
      fragments_removed++;
    }
  }

  if (nf == fragments_to_remove.number_elements())
  {
    m.resize(0);    // remove everything from the molecule
    return;
  }

  do_remove_by_fragment(m, fragments_to_remove);
  return;
}

static void
do_remove_nonorganic_fragments (Molecule & m, int & fragments_removed)
{
  fragments_removed = 0;
  int matoms = m.natoms ();

  int * frag_membership = new int[matoms]; 
  std::unique_ptr<int[]> free_frag_membership(frag_membership);
  m.fragment_membership (frag_membership);

  _do_remove_nonorganic_fragments (m, frag_membership, fragments_removed);

  if (verbose > 1)
    cerr << "new number of atoms: " << m.natoms() << endl;
}

static int
at_least_one_of_these_queries_matches (Molecule & m,
                               resizable_array_p<Substructure_Query> & q)
{
  Molecule_to_Match target(&m);

  int nq = q.number_elements();

  for (int i = 0; i < nq; i++)
  {
    if (q[i]->substructure_search(target))
      return 1;
  }

  return 0;
}

static int
identify_fragments_to_be_kept (Molecule & m,
                               resizable_array_p<Substructure_Query> & q)
{
  resizable_array_p<Molecule> frags;
  m.create_components(frags);

  resizable_array<int> fragments_to_remove;

  int nf = frags.number_elements();
  for (int i = 0; i < nf; i++)
  {
    if (! at_least_one_of_these_queries_matches(*(frags[i]), q))
      fragments_to_remove.add(i);
  }

  if (nf == fragments_to_remove.number_elements())
  {
    m.resize(0);    // remove everything from the molecule
    return 1;
  }

  return do_remove_by_fragment(m, fragments_to_remove);
}

static int
identify_fragment_by_query (resizable_array_p<Molecule> & components,
                            Substructure_Query * query,
                            int largest)
{
  int atoms_in_previous_best_match = largest ? 0 : std::numeric_limits<int>::max();
  int best_match_component = -1;

  for (int i = 0; i < components.number_elements(); i++)
  {
    if (query->substructure_search(components[i]))
    {
      int nc = components[i]->natoms();
      if (largest && nc > atoms_in_previous_best_match)
      {
        atoms_in_previous_best_match = nc;
        best_match_component = i;
      }
      else if (! largest && nc < atoms_in_previous_best_match)
      {
        atoms_in_previous_best_match = nc;
        best_match_component = i;
      }
    }
  }

  return best_match_component;
}
                            
 
static int
identify_fragment_by_query (Molecule & m, 
                            resizable_array_p<Substructure_Query> & queries,
                            int largest)
{
  resizable_array_p<Molecule> components;

  m.create_components(components);

  for (int i = 0; i < queries.number_elements(); i++)
  {
    int j = identify_fragment_by_query(components, queries[i], largest);
    if (j >= 0)
    {
      m.delete_all_fragments_except(j);
      if (verbose > 1)
        cerr << largest << " query matches fragment " << j << ", " << m.natoms() << " atoms\n";
      molecules_with_fragments_reduced_by_query++;
      return 1;
    }
  }

  molecules_not_matching_fragment_queries++;

  if (verbose)
    cerr << "None of " << components.number_elements() << " components matched queries '" << m.name() << "'\n";

  return 0;
}

static int
atoms_in_non_largest_fragment_exceed (Molecule & m,
                                      int mxfs)
{
  int nf = m.number_fragments();

  int atoms_in_largest_fragment = m.atoms_in_fragment(0);

//cerr << "Checking " << nf << " fragments\n";

  for (int i = 1; i < nf; i++)
  {
    int a = m.atoms_in_fragment(i);

    if (a > atoms_in_largest_fragment)
      atoms_in_largest_fragment = a;
  }

  for (int i = 0; i < nf; i++)
  {
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

static int
trim_to_first_n_fragments (Molecule & m,
                           const int * frag_membership,
                           int fragment_count)
{
  int natoms = m.natoms();
  for (int i = 0, j = 0; j < natoms; j++)
  {
    if (frag_membership[j] + 1 > fragment_count)
      m.remove_atom(i);
    else
      i++;
  }


  return 1;
}

static int 
do_strip_to_n_largest_fragments (Molecule & m,
                                 int n)
{
  resizable_array<int> fragment_sizes;

  int nf = m.number_fragments();

  if (nf <= n)   // nothing to do
    return 1;

  for (int i = 0; i < nf; i++)
  {
    int a = m.atoms_in_fragment(i);
    fragment_sizes.add(a);
  }

  Int_Comparator_Smaller icl;

  fragment_sizes.iwqsort(icl);

#ifdef DEBUG_DO_STRIP_TO_N_LARGEST_FRAGMENTS
  for (int i = 0; i < fragment_sizes.number_elements(); i++)
  {
    cerr << "Before erasure " << fragment_sizes[i] << endl;
  }
#endif

  fragment_sizes.erase(0, n - 1);

#ifdef DEBUG_DO_STRIP_TO_N_LARGEST_FRAGMENTS
  cerr << "Remaining fragment sizes " << fragment_sizes.number_elements() << endl;
  for (int i = 0; i < fragment_sizes.number_elements(); i++)
  {
    cerr << "After erasure " << fragment_sizes[i] << endl;
  }
#endif

  int * remove_fragment = new_int(nf); std::unique_ptr<int[]> free_remove_fragment(remove_fragment);

  for (int i = 0; i < nf; i++)
  {
    int a = m.atoms_in_fragment(i);

    if (! fragment_sizes.remove_first(a))
      continue;

    remove_fragment[i] = 1;
  }

  m.delete_fragments(remove_fragment);

  return 1;
}

/*
  Return 1 if the molecule should be kept, 0 if it should be discarded
*/

static int
process_fragments (Molecule & m)
{
  if (0 == maximum_fragment_count)
    ;
  else if (min_size_max_fragment_count)
  {
    if (too_many_large_fragments(m, min_size_max_fragment_count, maximum_fragment_count))
    {
      if (verbose > 1)
        cerr << "Too many large components " << m.number_fragments() << " max is " << maximum_fragment_count << endl;
      molecules_with_too_many_components++;
      return 0;
    }
  }
  else if (m.number_fragments() > maximum_fragment_count)
  {
    if (verbose > 1)
      cerr << "Too many components " << m.number_fragments() << " max is " << maximum_fragment_count << endl;
    molecules_with_too_many_components++;
    return 0;
  }
  
  if (remove_molecules_with_non_largest_fragment_natoms >= 0 && atoms_in_non_largest_fragment_exceed(m, remove_molecules_with_non_largest_fragment_natoms))
  {
    if (verbose > 1) {
      cerr << "Contains large fragment(s)\n";
    }
    molecules_with_large_fragments++;
    return 0;
  }

  if (remove_duplicate_fragments)
  {
    int fragments_removed;

    do_remove_duplicate_fragments(m, fragments_removed);

    if (fragments_removed)
    {
      if (verbose > 1) {
        cerr << fragments_removed << " duplicate fragments removed\n";
      }
      molecules_with_duplicate_fragments_removed++;
    }

    if (1 == m.number_fragments())
      return 1;
  }

  if (remove_fragments_this_size_or_smaller)
  {
    do_remove_fragments_this_size_or_smaller(m);
    if (1 == m.number_fragments())
      return 1;
  }

  if (strip_to_n_largest_fragments > 0)
    do_strip_to_n_largest_fragments(m, strip_to_n_largest_fragments);

  if (known_fragment_data.active())
  {
    known_fragment_data.process(m);
    if (1 == m.number_fragments())
      return 1;
  }

// If we have some fragments to definitely remove, process them

  if (remove_fragment_queries.number_elements())
    do_remove_fragments_that_match_query(m, remove_fragment_queries);

  if (reduce_to_largest_fragment)      // the most common case
  {
    if (verbose <= 1)
      return m.reduce_to_largest_fragment();

    int initial_matoms = m.natoms();

    int rc = m.reduce_to_largest_fragment();
    if (verbose > 1)
      cerr << "Stripped to largest fragment, lost " << (initial_matoms - m.natoms()) << " atoms\n";

    return rc;
  }

  if (reduce_to_all_largest_fragments)
    return do_reduce_to_all_largest_fragments(m);

  if (1 == remove_largest_fragment)
    return do_remove_largest_fragment(m);

  if (remove_largest_fragment > 1)
    return do_remove_largest_fragments(m, remove_largest_fragment);

  if (reduce_to_largest_organic_fragment)
  {
    if (verbose > 1)
      cerr << "Stripped to largest organic fragment\n";

    if (tag_for_removed_fragments.length())
      return extract_smaller_fragments_into_name_tag(m, tag_for_removed_fragments);
    else
      return m.reduce_to_largest_organic_fragment();
  }

  if (reduce_to_largest_organic_fragment_carefully)
  {
    if (verbose > 1)
      cerr << "Stripped to largest organic fragment with desirable features\n";

    return m.reduce_to_largest_fragment_carefully();
  }

  if (keep_all_organic_fragments) {
    int fragments_removed;

    do_remove_nonorganic_fragments(m, fragments_removed);

    if (fragments_removed > 0) {
      if (verbose) {
        cerr << fragments_removed << " fragments removed for " << m.name() << endl;
      }
      molecules_with_nonorganic_fragments_removed++;
      if (m.number_fragments() < 1) {
        return 0;
      }
    }

    return 1;
  }
  if (smallest_fragment_queries.number_elements())
    return identify_fragment_by_query(m, smallest_fragment_queries, 0);

  if (largest_fragment_queries.number_elements())
    return identify_fragment_by_query(m, largest_fragment_queries, 1);

  if (keep_fragment_queries.number_elements())
    return identify_fragments_to_be_kept(m, keep_fragment_queries);

  if (sort_by_fragment_size)
    do_sort_by_fragment_size(m);

  if (0 == fragment_count)
    return 1;

  int components = m.number_fragments();

  int * frag_membership = new int[m.natoms()]; std::unique_ptr<int[]> free_frag_membership(frag_membership);
  m.fragment_membership(frag_membership);

  int rc = trim_to_first_n_fragments(m, frag_membership, fragment_count);

  if (verbose > 1)
    cerr << "reduced from " << components << " fragments\n";

  return rc;
}

static int
max_ring_size (Molecule & m)
{
  int nr = m.nrings();

  if (0 == nr)
    return 0;

  int rc = 0;

  for (int i = nr - 1; i>= 0; i--)
  {
    const Ring * r = m.ringi(i);

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

static int
exclude_for_non_organic_and_non_periodic_table (const Molecule & m)
{
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    const Element * e = m.elementi(i);
    if (e->organic())
      continue;

    if (! ok_non_organics.matches(e))
    {
      non_organic_molecules_found++;
      return 1;
    }

    if (! e->is_in_periodic_table())
    {
      non_real_elements_found++;
      return 1;
    }
  }

  return 0;
}

static int
exclude_for_non_organic (const Molecule & m)
{
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    const Element * e = m.elementi(i);
    if (e->organic())
      continue;

    if (! ok_non_organics.matches (e))
    {
      non_organic_molecules_found++;
      return 1;
    }
  }

  return 0;
}

static int
exclude_for_non_real_elements (const Molecule & m)
{
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    const Element * e = m.elementi(i);

    if (! e->is_in_periodic_table())
    {
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

static int
exclude_for_atom_types (const Molecule & m)
{
  if (exclude_non_real_elements && output_organic_only)
    return exclude_for_non_organic_and_non_periodic_table(m);

  if (output_organic_only)
    return exclude_for_non_organic(m);

  return exclude_for_non_real_elements(m);
}

static int
valence_check_ok (Molecule & m)
{
  if (m.valence_ok())
    return 1;

  if (! ok_bad_valence_on_isotopically_labelled)
    return 0;

  const int matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i)
  {
    if (m.valence_ok(i))
      continue;

    if (0 == m.isotope(i))     // bad valence, but isotope, therefore OK
      return 0;
  }

  return 1;
}

static int
looks_like_mixture_by_atom_count (Molecule & m)
{
  int nf = m.number_fragments();

  int atoms_in_largest_fragment = m.atoms_in_fragment(0);
  int atoms_in_second_largest_fragment = m.atoms_in_fragment(1);

  if (atoms_in_largest_fragment < atoms_in_second_largest_fragment)
    std::swap(atoms_in_largest_fragment, atoms_in_second_largest_fragment);

//cerr << "First two counts " << atoms_in_largest_fragment << " and " << atoms_in_second_largest_fragment << endl;

  for (int i = 2; i < nf; i++)
  {
    int a = m.atoms_in_fragment(i);

    if (a > atoms_in_largest_fragment)
    {
      atoms_in_second_largest_fragment = atoms_in_largest_fragment;
      atoms_in_largest_fragment = a;
    }
    else if (a > atoms_in_second_largest_fragment)
      atoms_in_second_largest_fragment = a;
  }

  if (atoms_in_largest_fragment - atoms_in_second_largest_fragment <= mixture_if_largest_frags_differ_by)
    return 1;

  return 0;
}

static int
too_many_large_fragments (Molecule & m)
{
  int rc = 0;

  int nf = m.number_fragments();

  for (int i = 0; i < nf; i++)
  {
    if (m.atoms_in_fragment(i) > discard_molecule_if_multiple_fragments_larger_than)
    {
      rc++;

      if (rc > 1)
        return 1;
    }
  }

  return 0;
}

static int
count_ring_systems (Molecule & m)
{
  int * notused = new int[m.natoms()]; std::unique_ptr<int[]> free_notused(notused);

  return m.label_atoms_by_ring_system(notused);
}

static int
largest_number_rings_in_a_system (Molecule & m)
{
  int nr = m.nrings();

  if (0 == nr)
    return 0;

  if (1 == nr)
    return 1;

  if (2 == nr)
    return m.ringi(0)->is_fused() ? 2 : 1;

// Now things get more complicated

  int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  int rc = 1;

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])
      continue;

    const Ring * ri = m.ringi(i);

    if (! ri->is_fused())
      continue;

    int system_size = 1;

    for (int j = i + 1; j < nr; j++)
    {
      if (ring_already_done[j])
        continue;

      const Ring * rj = m.ringi(j);

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

static int
atoms_in_largest_fragment (Molecule & m)
{
  int nf = m.number_fragments();

  if (1 == nf)
    return m.natoms();

  int rc = m.atoms_in_fragment(0);

  for (int i = 1; i < nf; i++)
  {
    int tmp = m.atoms_in_fragment(i);

    if (tmp > rc)
      rc = tmp;
  }

  return rc;
}

static int 
_apply_all_filters (Molecule & m,
                    IWString & rejection_reason,
                    int & molecule_changed)
{
  assert(m.ok());

  int matoms = m.natoms();

  int nf = m.number_fragments();

  initial_fragment_count[nf]++;

  if (audit_input)
    ;
  else if (! filter_for_disallowed_elements)
    ;
  else if (allowed_elements.contains_non_allowed_atoms(m))
  {
    molecules_excluded_for_non_allowed_elements++;
    rejection_reason << "NonAllowedElement";
    return 0;
  }

// Apr 2013. User wants largest/smallest fragment that contains a query. But if the molecule
// has just one fragment, we still want that to work.

  if (1 == nf && smallest_fragment_queries.number_elements() && ! at_least_one_of_these_queries_matches(m, smallest_fragment_queries))
  {
    molecules_not_matching_fragment_queries++;
    rejection_reason << "NoMatchSmlNeededFrag";
    return 0;
  }

  if (1 == nf && largest_fragment_queries.number_elements() && ! at_least_one_of_these_queries_matches(m, largest_fragment_queries))
  {
    molecules_not_matching_fragment_queries++;
    rejection_reason << "NoMatchLrgNeededFrag";
    return 0;
  }

  if ((nf < 2) && (! keep_all_organic_fragments))
    // cannot reduce the fragment count, except that we may want to remove the
    // only fragment when it is not organic
    ;
  else if (discard_molecule_if_multiple_fragments_larger_than && 
           too_many_large_fragments(m))
  {
    molecules_with_multiple_large_fragments++;
    rejection_reason << "Mixture";
    return 0;
  }
  else if (mixture_if_largest_frags_differ_by >= 0 &&
           looks_like_mixture_by_atom_count(m))
  {
    molecules_declared_mixtures_by_atom_count_difference++;
    rejection_reason = "Mixture";
    return 0;
  }
  else if (need_to_call_process_fragments)
  {
    int initial_matoms = m.natoms();

    if (! process_fragments(m))
    {
      if (verbose > 1)
        cerr << "Molecule does not match fragment criteria\n";

      rejection_reason = "no match to fragment criteria";
      return 0;
    }

    matoms = m.natoms();

    if (initial_matoms != matoms)
    {
      molecule_changed++;
      molecules_chopped++;
    }

    atoms_lost.extra(initial_matoms - matoms);
  }

  if (aromatise_these_rings.number_elements())
  {
    do_aromatise_these_rings(m, aromatise_these_rings);
  }

  if (matoms)
  {
    molecule_changed += elements_to_remove.process(m);

    molecule_changed += element_transformations.process(m);

    if (structure_fixing.active())
      molecule_changed += structure_fixing.process(m);

    matoms = m.natoms();        // in case the atom count changed
  }

  int implicit_hydrogens;

  if (include_implicit_hydrogens_in_upper_atom_count_comparison)
    implicit_hydrogens = m.implicit_hydrogens();
  else
    implicit_hydrogens = 0;

  int keep = 0;

  if (atom_count_includes_only_atoms_in_largest_fragment)
    matoms = atoms_in_largest_fragment(m);

  if (lower_atom_count_cutoff > 0 && matoms < lower_atom_count_cutoff)
  {
     molecules_below_atom_count_cutoff++;
     if (verbose > 1)
       cerr << "Molecule contains " << m.natoms() << " atoms, which is below cutoff\n";

     rejection_reason = "too few atoms";
  }
  else if (upper_atom_count_cutoff > 0 && (matoms + implicit_hydrogens) > upper_atom_count_cutoff)
  {
     molecules_above_atom_count_cutoff++;
     if (verbose > 1)
       cerr << "Molecule contains " << m.natoms() << " atoms, which is above cutoff\n";
     rejection_reason = "too many atoms";
  }
  else if (lower_molecular_weight_cutoff > 0.0 && 
      m.molecular_weight_ignore_isotopes() < lower_molecular_weight_cutoff)
  {
    molecules_below_molecular_weight_cutoff++;
    cerr << "Molecular weight is " << m.molecular_weight_ignore_isotopes() << " which is below cutoff\n";
    rejection_reason = "AMW too low";
  }
  else if (upper_molecular_weight_cutoff > 0.0 &&
      m.molecular_weight_ignore_isotopes() > upper_molecular_weight_cutoff)
  {
    molecules_above_molecular_weight_cutoff++;
    cerr << "Molecular weight is " << m.molecular_weight_ignore_isotopes() << " which is above cutoff\n";
    rejection_reason = "AMW too high";
  }
  else if ((output_organic_only || exclude_non_real_elements) && exclude_for_atom_types(m))
  {
    if (verbose > 1)
      cerr << "Non organic or non periodic table atoms found\n";
    rejection_reason = "non organic";
  }
  else if (lower_ring_count_cutoff && m.nrings() < lower_ring_count_cutoff)
  {
    molecules_with_too_few_rings++;
    if (verbose > 1)
      cerr << "Molecule contains " << m.nrings() << " rings, which is below cutoff\n";
    rejection_reason = "too few rings";
  }
  else if (upper_ring_count_cutoff >= 0 && m.nrings() > upper_ring_count_cutoff)
  {
    molecules_with_too_many_rings++;
    if (verbose > 1)
      cerr << "Molecule contains " << m.nrings() << " rings, which is above cutoff\n";
    rejection_reason = "too many rings";
  }
  else if (max_ring_systems_allowed >= 0 && count_ring_systems(m) > max_ring_systems_allowed)
  {
    molecules_with_too_many_ring_systems++;
    if (verbose > 1)
      cerr << "Molecule contains too many ring systems\n";
    rejection_reason << "too many ring systems";
  }
  else if (upper_ring_size_cutoff > 0 && max_ring_size(m) > upper_ring_size_cutoff)
  {
    molecules_with_large_rings++;
    if (verbose > 1)
      cerr << "Molecule contains ring too large\n";
    rejection_reason = "large ring";
  }
  else if (max_rings_in_a_ring_system > 0 && largest_number_rings_in_a_system(m) > max_rings_in_a_ring_system)
  {
    molecules_with_ring_systems_too_large++;
    if (verbose > 1)
      cerr << "Molecule has ring system with too many rings\n";
    rejection_reason = "large ring system";
  }
  else if (exclude_isotopes && m.number_isotopic_atoms())
  {
    molecules_containing_isotopes++;
    if (verbose > 1)
      cerr << "Molecule contains isotopes\n";
    rejection_reason = "isotopes";
  }
  else if (skip_molecules_with_abnormal_valences &&
           ! valence_check_ok(m))
  {
    molecules_with_abnormal_valences++;
    if (verbose > 1)
      cerr << "Molecule contains abnormal valence(s)\n";
    rejection_reason = "bad valence";
  }
  else if (max_path_length > 0 && m.longest_path() > max_path_length)
  {
    molecules_with_longest_path_too_long++;
    if (verbose > 1)
      cerr << "Longest path too long\n";
    rejection_reason << "long path";
  }
  else
    keep = 1;     // molecule is good!

  if (! keep)
    return 0;

  if (remove_isotopic_atoms)
    molecule_changed += do_remove_isotopes(m);

// After we have done all other manipulations with isotopes

  if (convert_isotopes_to_atom_map_numbers)
    do_convert_isotopes_to_atom_map_numbers(m);
  else if (convert_atom_map_numbers_to_isotopes)
    do_convert_atom_map_numbers_to_isotopes(m);

  matoms = m.natoms();     // recompute, as it may have been changed by transformations

  if (compute_molecular_weight_for_each || lower_amw_cutoff >= 0.0 || upper_amw_cutoff >= 0.0)
  {
    molecular_weight_t amw = m.molecular_weight();

    amw_accumulator.extra(amw);
    atom_count[matoms]++;

    if (0 == verbose && compute_molecular_weight_for_each)
      cerr << m.molecule_name() << ' ';
    if (verbose > 1 || compute_molecular_weight_for_each)
      cerr << "AMW " << amw;

    if (lower_amw_cutoff >= 0.0 && amw < lower_amw_cutoff)
    {
      if (verbose > 1)
        cerr << " below cutoff\n";
      else if (compute_molecular_weight_for_each)
        cerr << endl;
      molecules_below_amw_cutoff++;
      return 0;
    }

    if (upper_amw_cutoff > 0.0 && amw > upper_amw_cutoff)
    {
      if (verbose > 1)
        cerr << " above cutoff\n";
      else if (compute_molecular_weight_for_each)
        cerr << endl;
      molecules_above_amw_cutoff++;
      return 0;
    }

    if (verbose > 1 || compute_molecular_weight_for_each)
      cerr << endl;
  }

  if (verbose)
    chiral_centre_count[m.chiral_centres()]++;

  if (max_chiral_centres > 0 && m.chiral_centres() > max_chiral_centres)
  {
    if (verbose > 1)
      cerr << " too many chiral centres " << max_chiral_centres << endl;
    molecules_with_too_many_chiral_centres++;
    rejection_reason << "too many chiral centres " << m.chiral_centres();
    return 0;
  }

  return 1;
}

static int
apply_all_filters (Molecule & m,
                   int molecule_number,
                   int & molecule_changed)
{
  IWString rejection_reason;
  int rc = _apply_all_filters(m, rejection_reason, molecule_changed);

  if (rc)      // molecule is OK.
    return rc;

  if (reject_log.rdbuf()->is_open() && reject_log.good())
    reject_log << molecule_number << " '" << m.name() << "' REASON " << rejection_reason << endl;

  if (rejections_output_object.good())
    rejections_output_object.write(m);

  return 0;
}

static int
do_convert_specific_isotopes_query(Molecule & m,
                                   Molecule_to_Match & target,
                                   Substructure_Query & q,
                                   const int new_isotope)
{
  Substructure_Results sresults;

  const int nhits = q.substructure_search(target, sresults);

  if (0 == nhits)
    return 0;

  int rc = 0;

  for (int i = 0; i < nhits; ++i)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    const int n = e->number_elements();

    for (int j = 0; j < n; ++j)
    {
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

static int
do_convert_specific_isotopes_query(Molecule & m,
                                   const resizable_array_p<Substructure_Query> & convert_specific_isotopes_query,
                                   const resizable_array<int> & convert_specific_isotopes_query_new_isotope)
{
  assert (convert_specific_isotopes_query.number_elements() == convert_specific_isotopes_query_new_isotope.number_elements());
  int rc = 0;

  Molecule_to_Match target(&m);

  for (int i = 0; i < convert_specific_isotopes_query.number_elements(); ++i)
  {
    if (do_convert_specific_isotopes_query(m, target, *convert_specific_isotopes_query[i],
                                           convert_specific_isotopes_query_new_isotope[i]))
      rc++;
  }

  return rc;
}

static int
do_convert_all_isotopes_to(Molecule & m,
                           const int convert_all_isotopes_to)
                
{
  const int matoms = 0;

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (0 ==  m.isotope(i))
      continue;

    m.set_isotope(i, convert_all_isotopes_to);
    rc++;
  }

  return rc;
}

static int
do_convert_specific_isotopes(Molecule & m,
                             const resizable_array<int> & convert_specific_isotopes,
                             const resizable_array<int> & convert_specific_isotopes_new_isotope)
{
  const int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    const int iso = m.isotope(i);

    if (0 == iso)
      continue;

    const int f = convert_specific_isotopes.index(iso);
//  cerr << "Do we need to convert isotope " << iso << " f = " << f << endl;
    if (f < 0)
      continue;

    m.set_isotope(i, convert_specific_isotopes_new_isotope[f]);
    rc++;
  }

  return rc;
}

static int
contains_covalently_bonded_non_organic (const Molecule & m)
{
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    if (a->element()->organic())
      continue;

    if (1 != a->ncon())    // non organic bonded to something
      return 1;
  }

  return 0;
}

/*
  How many of the leading characters are whitespace
*/

static int
count_leading_whitespace (const char * s,
                          int lens)
{
  for (int i = 0; i < lens; i++)
  {
    if (! isspace(s[i]))
      return i;
  }

  return lens;
}

static int
do_truncate_name (Molecule & m,
                  char truncate_name_at_first, 
                  char truncate_name_at_last)
{
  IWString newname = m.name();

  if ('\0' != truncate_name_at_first)
  {
    int i = newname.index(truncate_name_at_first);
    if (i >= 0)
      newname.iwtruncate(i);
  }

  if ('\0' != truncate_name_at_last)
  {
    int i = newname.index(truncate_name_at_last);
    if (i >= 0)
      newname.iwtruncate(i);
  }

  if (0 == newname.length())
  {
    cerr << "do_truncate_name:cannot empty name '" << m.name() << "'\n";
    return 1;
  }
  
  if (newname.length() == m.name().length())   // no change
    return 1;

  m.set_name(newname);

  return 1;
}

static void
do_truncate_names_to_first_token(Molecule & m)
{
  const IWString & mname = m.name();

  int n = mname.length();

  if (0 == n)    // molecule with no name
    return;

  const char * s = mname.rawchars();

  int lws = count_leading_whitespace(s, n); 

  if (lws)
  {
    s += lws;
    n -= lws;

    if (0 == n)   // name was all whitespace
    {
      m.set_name("");
      return;
    }
  }

  for (int i = 0; i < n; i++)
  {
    if (! isspace(s[i]))
      continue;

    IWString tmp(s, i);
    m.set_name(tmp);
    return;
  }

  return;
}

static int
change_name (Molecule & m,
             IW_Regular_Expression  & change_name_rx)
{
  const IWString & mname = m.name();

  int nmatches = change_name_rx.matches_save_subexpressions(mname);

  if (nmatches)
    ;
  else if (change_name_rx_must_match)
  {
    cerr << "change_name:change name rx not matched '" << change_name_rx.source() << "', molecule name '" << m.name() << "'\n";
    return 0;
  }
  else              // ok for no matches
    return 1;

  IWString newname;
  change_name_rx.dollar(1, newname);

  if (0 == newname.length())
  {
    cerr << "change_name:cannot set name to length zero, rx '" << change_name_rx.source() << "', mname '" << m.name() << "'\n";
    return 0;
  }

  m.set_name(newname);

  return 1;
}

static void
do_substitute_for_whitespace_in_name(Molecule & m)
{
  assert(1 == substitute_for_whitespace_in_name.length());

  IWString tmp(m.name());

  int n = tmp.length();

  int need_to_reset_name = 0;

  for (int i = 0; i < n; i++)
  {
    if (! isspace(tmp[i]))
      continue;

    tmp[i] = substitute_for_whitespace_in_name[0];
    need_to_reset_name = 1;
  }

  if (need_to_reset_name)
    m.set_name(tmp);

  return;
}

static int
do_remove_unnecessary_square_brackets(Molecule & m)
{
  return m.unset_unnecessary_implicit_hydrogens_known_values();
}

static int
do_remove_all_possible_square_brackets(Molecule & m)
{
  return m.remove_hydrogens_known_flag_to_fix_valence_errors();
}

static int
do_debug_print (Molecule & m,
                std::ostream & os)
{
  if (0 == verbose)
    os << "Molecule " << molecules_read << '\n';
  m.ring_membership();
  m.debug_print(os);

  return os.good();
}

static int
fileconv (Molecule & m,
          Molecule_Output_Object & output_object)
{
  if (name_rx.active() && ! name_rx.matches(m.name()))
    return 1;

// if (verbose > 1)
//   cerr << "Molecule " << input.molecules_read() << " finishes at line " << input.lines_read() << endl;

  if (print_bond_lengths)
    (void) do_print_bond_lengths(m, cout);

  if (print_bond_angles)
    (void) do_print_bond_angles(m, cout);

  if (print_torsions)
    (void) do_print_torsions(m, cout);

  if (verbose)
  {
    int matoms = m.natoms();
    natoms_accumulator(matoms);
    if (verbose > 1)
      atom_count[matoms]++;
  }

  if (audit_input)
  {
    if (appends_to_be_done)
    {
      do_appends(m);
      cerr << molecules_read << ' ' << m.name() << endl;
    }

    return 1;
  }

  int molecule_changed = 0;

  if (convert_isotopes)    // put here so chemical standardisation can get rid of converted Deuterium
  {
    if (m.transform_to_non_isotopic_form())
    {
      molecule_changed++;
      molecules_containing_isotopes++;
    }
  }
  else if (convert_all_isotopes_to)
  {
    if (do_convert_all_isotopes_to(m, convert_all_isotopes_to))
    {
      molecules_changed++;
      molecules_containing_isotopes++;
    }
  }
  else
  {
    if (convert_specific_isotopes.number_elements())
    {
      if (do_convert_specific_isotopes(m, convert_specific_isotopes, convert_specific_isotopes_new_isotope))
      {
        molecules_changed++;
        molecules_containing_isotopes++;
      }
    }

    if (convert_specific_isotopes_query.number_elements())
    {
      if (do_convert_specific_isotopes_query(m, convert_specific_isotopes_query, convert_specific_isotopes_query_new_isotope))
      {
        molecules_changed++;
        molecules_containing_isotopes++;
      }
    }
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

  if (! stream_for_smiles_before_filters.is_open())
    ;
  else if (m.number_isotopic_atoms())   // not in the organic subset
    ;
  else if (contains_covalently_bonded_non_organic(m))   // not in the organic subset
    ;
  else   // organic subset
  {
    stream_for_smiles_before_filters << m.smiles();
    stream_for_smiles_before_filters << ' ' << m.name() << "\n";
    stream_for_smiles_before_filters.write_if_buffer_holds_more_than(32768);
  }

  if (remove_invalid_chiral_centres)
  {
    int nr = do_remove_invalid_chiral_centres(m);
    if (nr)
    {
      molecules_with_invalid_chiral_centres++;

      if (verbose > 1)
        cerr << "Removed " << nr << " invalid chiral centres\n";

      molecule_changed = 1;
    }
  }
  else if (remove_chiral_data_from_all_molecules)
    molecule_changed += m.remove_all_chiral_centres();
  else if (remove_chiral_centres_on.number_elements())
    do_remove_chiral_centres_on_matched_atoms(m, remove_chiral_centres_on);
  
  if (remove_non_organic_chirality)
    do_remove_non_organic_chirality(m);

  if (! apply_all_filters(m, molecules_read, molecule_changed))
    return 1;

  if (truncate_name_at_first || truncate_name_at_last)
    do_truncate_name(m, truncate_name_at_first, truncate_name_at_last);
  else if (truncate_names_to_first_token)
    do_truncate_names_to_first_token(m);
  else if (substitute_for_whitespace_in_name.length())
    do_substitute_for_whitespace_in_name(m);
  else if (! change_name_rx.active())
    ;
  else if (change_name(m, change_name_rx))
    ;
  else
  {
    if (verbose)
      cerr << "Cannot change name according to rx '" << change_name_rx.source() << "'\n";
    return 0;
  }

  if (prepend_to_name.length() > 0)
  {
    IWString tmp(prepend_to_name);
    tmp << m.name();
    m.set_name(tmp);
  }

  if (number_assigner.active())
    number_assigner.process(m);

  if (remove_directional_bonds_from_input)
    m.revert_all_directional_bonds_to_non_directional();

  if (reflect_coordinates)
    molecule_changed += do_reflect_coordinates(m);
  else if (invert_all_chiral_centres)
    molecule_changed += do_invert_all_chiral_centres(m);
  else if (find_all_chiral_centres || find_all_ring_chiral_centres)
    molecule_changed += do_find_all_chiral_centres(m);

  if (translation_specified)
    m.translate_atoms(dx, dy, dz);

  if (charge_assigner.active())
  {
    molecule_changed += charge_assigner.process(m);
  }

  if (donor_acceptor_assigner.active())
  {
    molecule_changed += donor_acceptor_assigner.process(m);
  }

  if (make_all_implicit_hydrogens_explicit)
  {
    if (m.make_implicit_hydrogens_explicit())
    {
      molecules_to_which_hydrogens_were_added++;
      molecule_changed++;
    }
  }
  else if (atoms_for_implicit_hydrogens.number_elements())
  {
    if (do_make_implicit_hydrogens_explicit(m))
    {
      molecules_to_which_hydrogens_were_added++;
      molecule_changed++;
    }
  }

  if (0 != fileconv_partial_charge_type)
    do_compute_partial_charges(m);

  if (appends_to_be_done)
    do_appends(m);

  if (molecule_changed)
  {
    molecules_changed++;
    if (molecule_changed_string.length())
    {
      IWString tmp(m.name());
      tmp.append_with_spacer(molecule_changed_string);
      m.set_name(tmp);
    }
  }

  if (! output_object.write(m))
    return 0;

  if (debug_print_each_molecule > 1)
    m.debug_print(cerr);

  molecules_written++;

  return 1;
}

static int
fileconv (data_source_and_type<Molecule> & input,
          Molecule_Output_Object & output_object)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
//  std::unique_ptr<Molecule> free_m(m);
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (debug_print_each_molecule)
    {
      m->compute_aromaticity_if_needed();
      do_debug_print(*m, cerr);
    }
    
    if (! fileconv(*m, output_object))
      return 0;
  }

//cerr << "At end of loop, " << input.connection_table_errors_encountered() << endl;

  if (input.stopped_because_of_error())
    return 0;

  return 1;
}
/*
*/

static int
fileconv (const char *fname, int input_type,
          Molecule_Output_Object & output)
{
  assert(NULL != fname);

  if (0 == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  if (connection_table_errors_allowed)
  {
    input.set_connection_table_errors_allowed(connection_table_errors_allowed);
    if (connection_table_error_file.length())
      input.set_connection_table_error_file(connection_table_error_file);
  }

  if (audit_input)
    ;
  else if (output.name_token_for_file_name() >= 0)    // output file names generated by the output object
    ;
  else if (output_file_stem.length())
    ;
  else if (output.would_use_name(fname))
  {
    cerr << "fileconv: input '" << fname << "' and output must be distinct\n";
    return 0;
  }

// Set up the output object for this file stem
// If there is a new stem for output files, make sure we get it.
// Make sure we deal properly with the case of multiple input files
// and a single output file (via the -S option)

  static int first_call = 1;

  int rc = 0;
  if (audit_input)
    rc = 1;
  else if (output_file_stem.length())
  {
    if (first_call)
      rc = output.new_stem(output_file_stem);
    else
      rc = output.ok();

    first_call = 0;
  }
  else if (output.name_token_for_file_name() >= 0)
    rc = 1;
  else
    rc = output.new_stem(fname);

  if (0 == rc)
  {
    cerr << "Output object could not open file\n";
    return 0;
  }

  return fileconv(input, output);
}

static int
fileconv_list_of_files (iwstring_data_source & input,
                        int input_type,
                        Molecule_Output_Object & output)
{
  input.set_strip_trailing_blanks(1);
  input.set_skip_blank_lines(1);
  input.set_strip_leading_blanks(1);

  IWString buffer;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;

    if (verbose > 1)
      cerr << "Processing '" << buffer << "'\n";
    if (! fileconv(buffer.null_terminated_chars(), input_type, output))
    {
      cerr << "Fatal error processing file '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
fileconv_list_of_files (const char * fname,
                        int input_type,
                        Molecule_Output_Object & output)
{
  const char * tmp = fname + 2;

  iwstring_data_source input(tmp);

  if (! input.good())
  {
    cerr << "Cannot open list of files to process '" << fname << "'\n";
    return 0;
  }

  return fileconv_list_of_files(input, input_type, output);
}

static const Element *
recognise_as_element (const const_IWSubstring & o)
{
  const Element * rc = get_element_from_symbol_no_case_conversion(o);

  if (NULL == rc)
    return NULL;

  atomic_number_t z = rc->atomic_number();

  if (z <= 0)
    return NULL;

  if (z > HIGHEST_ATOMIC_NUMBER)   // not sure how this could happen
    return NULL;

  return rc;
}

static const Element * 
recognise_as_element_or_atomic_number (const const_IWSubstring & o)
{
  int z;
  if (! o.numeric_value(z))
    return recognise_as_element(o);

  if (z <= 0)
    return NULL;

  if (z > HIGHEST_ATOMIC_NUMBER)
    return NULL;

  return get_element_from_atomic_number(z);
}

static int
handle_dash_O_stuff (Command_Line & cl,
                     char flag)
{
  int i = 0;
  const_IWSubstring o;
  while (cl.value(flag, o, i++))
  {
    if ("DEF" == o || "def" == o)
    {
      filter_for_disallowed_elements = 1;
    }
    else if (o.starts_with("allow:"))
    {
      o.remove_leading_chars(6);

      const Element * e = recognise_as_element_or_atomic_number(o);

      if (NULL == e)
      {
        cerr << "Unrecognised element allow:'" << o << "'\n";
        return 0;
      }

      allowed_elements.set_allow(e->atomic_number(), 1);

      if (verbose)
        cerr << "Atomic number " << e->atomic_number() << " allowed\n";
    }
    else if ("none" == o)
    {
      output_organic_only = 1;

      continue;
    }
    else if ("help" == o)
    {
      cerr << "The -O option is troublesome due to it's long history and evolution\n";
      cerr << "Options can also be order dependent\n";
      cerr << " -O def         reject if any of the non-OK non-organics are present\n";
      cerr << " -O allow:El    temporarily allow El as an OK non-organic\n";
      cerr << " -O none        reject if any non-organic atoms present\n";
      cerr << " -O El          element El becomes fully organic for the course of the run\n";
      exit(1);
    }
    else
    {
      Element_Matcher * e = new Element_Matcher;
      if (! e->construct_from_string(o))
      {
        cerr << "Unrecognised -O qualifier '" << o << "'\n";
        delete e;
        return 0;
      }

      ok_non_organics.add(e);

      const Element * ele = e->element();
      if (NULL != ele)
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

static int
get_translations (Command_Line & cl, int verbose, const char cflag)
{
  IWString token;
  (void) cl.value(cflag, token);

  int i = 0;
  IWString d;
  if (! token.nextword(d, i, ','))
  {
    cerr << "Cannot extract dx from -T switch\n";
    return 0;
  }

  double tmp;
  if (! d.numeric_value(tmp))
  {
    cerr << "The value for dx must be a floating point number\n";
    return 0;
  }

  dx = (coord_t) tmp;

  if (! token.nextword(d, i, ','))
  {
    cerr << "Cannot extract dy from -T switch\n";
    return 0;
  }

  if (! d.numeric_value(tmp))
  {
    cerr << "The value for dy must be a floating point number\n";
    return 0;
  }

  dy = (coord_t) tmp;

  if (! token.nextword(d, i, ','))
  {
    cerr << "Cannot extract dz from -T switch\n";
    return 0;
  }

  if (! d.numeric_value(tmp))
  {
    cerr << "The value for dz must be a floating point number\n";
    return 0;
  }
  dz = (coord_t) tmp;

  if (verbose)
    cerr << "Molecules translated " << dx << ", " << dy << " ," << dz << endl;

  translation_specified = 1;

  return 1;
}

static int
all_files_recognised_by_type (const Command_Line & cl)
{
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! discern_file_type_from_name(cl[i]))
    {
      cerr << "Cannot determine file type from '" << cl[i] << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
fetch_frag_smarts (const const_IWSubstring f, const char * largest_or_smallest,
                    resizable_array_p<Substructure_Query> & queries)
{
  Substructure_Query * tmp = new Substructure_Query(f);
  if (! tmp->create_from_smarts(f))
  {
    cerr << "Cannot create query from smarts '" << f << "'\n";
    delete tmp;
    return 0;
  }

  queries.add(tmp);

  if (verbose)
    cerr << "Will keep " << largest_or_smallest << " fragment which matches query '"
         << f << "'\n";

  return 1;
}

static int
fetch_frag_queries (const_IWSubstring f, const char * largest_or_smallest,
                    resizable_array_p<Substructure_Query> & queries)
{
  if (f.starts_with("F:"))
  {
    f.remove_leading_chars(2);
    if (! queries_from_file(f, queries, 1, verbose))   // let's assume queries in same directory as file
    {
      cerr << "Could not process 'F:" << f << "'\n";
      return 0;
    }

    return 1;
  }

  if (f.starts_with("M:"))
  {
    f.remove_leading_chars(2);

    if (! queries_from_file_of_molecules(f, queries, verbose))
      return 0;

    return queries.number_elements();
  }

  Substructure_Query * tmp = new Substructure_Query(f);
  if (! tmp->read(f))
  {
    cerr << "fetch_frag_queries: cannot read query from '" << f << "'\n";
    delete tmp;
    return 0;
  }

  queries.add(tmp);
  if (verbose)
    cerr << "Will keep " << largest_or_smallest << " fragment which matches query '"
         << f << "'\n";

  return 1;
}

static int
fileconv (int argc, char ** argv)
{
  if (! first_call)
    reset_file_scope_variables();

  first_call = 0;

  Command_Line cl(argc, argv, "PN:T:eB:I:s:g:h:H:t:n:L:S:GA:K:w:W:m:X:c:C:Q:O:E:f:F:vVi:ao:r:R:p:J:Y:");
  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (verbose)
    cerr << __FILE__ << " compiled " << __TIME__ << " " << __DATE__ << endl;

  if (! process_elements(cl))
    usage(2);

  if (! process_standard_smiles_options(cl, verbose))
  {
    usage(4);
  }

  if (! process_standard_aromaticity_options(cl, verbose))
  {
    usage(5);
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      usage(6);
    }
  }

  if (cl.option_present('N'))
  {
    if (! charge_assigner.construct_from_command_line(cl, verbose, 'N'))
    {
      cerr << "Cannot determine charge assigner from command line\n";
      usage(77);
    }

    if (Daylight != global_aromaticity_type())
    {
      set_global_aromaticity_type(Daylight);
      if (verbose)
        cerr << "Global aromaticity type set to Dayligth for Charge Assigner";
    }
  }

  if (cl.option_present('H'))
  {
    if (! donor_acceptor_assigner.construct_from_command_line(cl, 'H', verbose))
    {
      cerr << "Cannot initialise donor/acceptor assigner (-H option)\n";
      usage(6);
    }
  }

// By default we echo any chiral info present in the input

  set_include_chiral_info_in_smiles(1);

  if (cl.option_present('w'))
  {
    IWString mtmp;
    int i = 0;
    while (cl.value('w', mtmp, i++))
    {
      if ('*' == mtmp)
      {
        compute_molecular_weight_for_each = 1;
        if (verbose)
          cerr << "The molecular weight of each molecule will be reported\n";
      }
      else if ("LARGE" == mtmp)
      {
        compute_molecular_weight_based_on_largest_fragment = 1;
        if (verbose)
          cerr << "Molecular weights will be for the largest fragment\n";
      }
      else
      {
        double tmp;
        if (! mtmp.numeric_value(tmp) || tmp <= 0.0)
        {
          cerr << "The lower molecular weight cutoff must be a positive number '" << mtmp << "'\n";
          usage(21);
        }

        lower_amw_cutoff = tmp;
        if (verbose)
          cerr << "Molecules with amw below " << lower_amw_cutoff << " will be excluded\n";
      }
    }

    atom_count.resize(300);
  }

  if (cl.option_present('W'))
  {
    IWString mtmp;
    int i = 0;
    while (cl.value('W', mtmp, i++))
    {
      if ('*' == mtmp)
      {
        compute_molecular_weight_for_each = 1;
        if (verbose)
          cerr << "The molecular weight of each molecule will be reported\n";
      }
      else if ("LARGE" == mtmp)
      {
        compute_molecular_weight_based_on_largest_fragment = 1;
        if (verbose)
          cerr << "Molecular weights will be for the largest fragment\n";
      }
      else
      {
        double tmp;
        if (! mtmp.numeric_value(tmp) || tmp <= 0.0)
        {
          cerr << "The lower molecular weight cutoff must be a positive number '" << mtmp << "'\n";
          usage(23);
        }
        if (lower_amw_cutoff && tmp < lower_amw_cutoff)
        {
          cerr << "The upper amw cutoff must be greater than the lower " << tmp << endl;
          usage(22);
        }

        upper_amw_cutoff = tmp;
        if (verbose)
          cerr << "Molecules with amw above " << upper_amw_cutoff << " will be excluded\n";
      }
    }
  }
  if (cl.option_present('t'))
  {
    if (! element_transformations.construct_from_command_line(cl, verbose, 't'))
      usage(8);
  }

// Processing the -O option is hard because it might be either
//   A set of elements
//   A set of directives for allowable elements

  if (cl.option_present('O'))
  {
    if (! handle_dash_O_stuff(cl, 'O'))
    {
      display_dash_O_options('O', cerr);
      return 3;
    }
  }

  if (cl.option_present('e'))
  {
    exclude_non_real_elements = 1;

    set_auto_create_new_elements(1);

    if (verbose)
      cerr << "Molecules containing elements not in periodic table will be discarded\n";
  }

  int input_type = 0;
  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }

  if (0 != input_type)   // great, explicitly specified
    ;
  else if (1 == cl.number_elements() && 0 == strcmp("-", cl[0]))   // reading from a pipe, assume smiles input
  {
    if (verbose)
      cerr << "Assuming smiles input from pipe read\n";
    input_type = SMI;
  }
  else if (all_files_recognised_by_type(cl))
    ;
  else
  {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  if (cl.option_present('a'))
  {
    audit_input = 1;
    if (verbose)
      cerr << "No output, just audit input\n";
    if (cl.option_present('o') || cl.option_present('S'))
    {
      cerr << "The -a and -o options are incompatible\n";
      usage(5);
    }
  }

  if (cl.option_present('S'))
  {
    if (1 != cl.option_count('S'))
    {
      cerr << "Sorry, only one -S option possible\n";
      usage(3);
    }

    cl.value('S', output_file_stem);
    if (verbose)
      cerr << "New files will be created with stem '" << output_file_stem << "'\n";

  }

  MDL_File_Supporting_Material * mdlfos =  global_default_MDL_File_Supporting_Material();

  if (cl.option_present('s'))
  {
    int i = 0;
    IWString tmp;
    while (cl.value('s', tmp, i++))
    {
      if ("help" == tmp)
      {
        display_chirality_options();
        exit(1);
      }

      if ('1' == tmp || tmp.starts_with("incl"))
      {
        set_include_chiral_info_in_smiles(1);
        if (verbose)
          cerr << "Chiral information will be written\n";
      }
      else if ('0' == tmp || tmp.starts_with("excl"))
      {
        set_include_chiral_info_in_smiles(0);
        mdlfos->set_include_chiral_info_in_mdl_outputs(0);
        if (verbose)
          cerr << "Chiral information will NOT be written\n";
      }
      else if ("find" == tmp)
      {
        find_all_chiral_centres = 1;
        if (verbose)
          cerr << "Chiral centres will be identified\n";
      }
      else if ("rfind" == tmp)
      {
        find_all_ring_chiral_centres = 1;
        if (verbose)
          cerr << "Chiral centres in ring systems will be identified\n";
      }
      else if ("invert" == tmp)
      {
        invert_all_chiral_centres = 1;
        if (verbose)
          cerr << "Chiral centres will be inverted\n";
      }
      else if ("reflect" == tmp)
      {
        reflect_coordinates = 1;
        if (verbose)
          cerr << "Will reflect coordinates and invert chiral centres\n";
      }
      else if ("discard" == tmp)
      {
        set_ignore_all_chiral_information_on_input(1);
        if (verbose)
          cerr << "All chiral information in input will be discarded\n";
      }
      else if ("remove" == tmp)
      {
        remove_chiral_data_from_all_molecules = 1;
        if (verbose)
          cerr << "Will strip all chiral information from molecules\n";
      }
      else if ("good" == tmp || "ignore" == tmp)
      {
        set_ignore_incorrect_chiral_input(1);
        if (verbose)
          cerr << "Incorrect chiral specifications will be ignored\n";
      }
      else if ("rmbad" == tmp)
      {
        remove_invalid_chiral_centres = 1;

        if (verbose)
          cerr << "Invalid chiral centres will be removed\n";
      }
      else if (tmp.starts_with("rmchiral="))
      {
        const_IWSubstring foo(tmp);  // process_cmdline_token may change its argument

        foo.remove_leading_chars(9);

        if (! process_cmdline_token(' ', foo, remove_chiral_centres_on, verbose))
        {
          cerr << "Cannot read chiral centre removal query(s) '" << tmp << "'\n";
          return 3;
        }
      }
      else if ("rmctb" == tmp || "xctb" == tmp)
      {
        remove_directional_bonds_from_input = 1;

        if (verbose)
          cerr << "Directional bonds removed from molecules\n";
      }
      else if ("rmbctb" == tmp)
      {
        remove_invalid_directional_bonds_from_input = 1;

        if (verbose)
          cerr << "Invalid directional bonds will be removed\n";
      }
      else if ("nonocm" == tmp)
      {
        mdlfos->set_display_non_organic_chirality_messages(0);
        if (verbose)
          cerr << "Will not display messages about non organic chirality\n";
      }
      else if ("rmnoc" == tmp)
      {
        remove_non_organic_chirality = 1;

        if (verbose)
          cerr << "Chirality will be removed from non organic atoms\n";
      }
      else if (tmp.starts_with("maxc="))
      {
        tmp.remove_leading_chars(5);
        if (! tmp.numeric_value(max_chiral_centres) || max_chiral_centres < 0)
        {
          cerr << "The max chiral centre count (-s maxc=) directive must be a whole number\n";
          return 1;
        }

        if (verbose)
          cerr << "Will discard molecules having more than " << max_chiral_centres << " chiral centres\n";
      }
      else
      {
        cerr << "Unrecognised stereo/chiral specifier (-s) '" << tmp << "'\n";
        usage(32);
      }
    }

    for (int i = 0; i < remove_chiral_centres_on.number_elements(); i++)
    {
      remove_chiral_centres_on[i]->set_find_unique_embeddings_only(1);
    }
  }

  if (cl.option_present('Q'))
  {
    const_IWSubstring q;
    int i = 0;

    while (cl.value('Q', q, i++))
    {
      if ("help" == q)
      {
        cerr << "The following -Q qualifiers are recognised\n";
        cerr << " -Q gast    Gasteiger partial charges\n";
        cerr << " -Q gh      Gasteiger Huckel partial charges\n";
        cerr << " -Q abraham Abraham partial charges\n";

        return 0;
      }
      else if (q.starts_with("gast"))
      {
        fileconv_partial_charge_type = FILECONV_PARTIAL_CHARGE_GASTEIGER;
        if (verbose)
          cerr << "Will compute Gasteiger partial charges\n";
      }
      else if ("gh" == q)
      {
        fileconv_partial_charge_type = FILECONV_PARTIAL_CHARGE_GASTEIGER_HUCKEL;
        if (verbose)
          cerr << "Will compute Gasteiger Huckel partial charges\n";
      }
      else if ("abraham" == q)
      {
        fileconv_partial_charge_type = FILECONV_PARTIAL_CHARGE_ABRAHAM;
        if (verbose)
          cerr << "Will compute Abraham partial charges\n";
      }
      else
      {
        cerr << "Unrecognised partial charge specification '" << q << "'\n";
        usage(12);
      }
    }
  }

  if (cl.option_present('J'))
  {
    if (! structure_fixing.initialise(cl, 'J', verbose > 1))
    {
      cerr << "Cannot initialise structure fixing (-J option)\n";
      usage(4);
    }
  }

  if (cl.option_present('X'))
  {
    if (! elements_to_remove.construct_from_command_line(cl, verbose, 'X'))
    {
      cerr << "Cannot discern elements to remove from -X switch\n";
      usage(18);
    }
  }

  if (cl.option_present('f'))
  {
    const_IWSubstring f;
    int i = 0;
    while (cl.value('f', f, i++))
    {
      if ('l' == f || f.starts_with("large") || f.starts_with("LARGE"))
      {
        reduce_to_largest_fragment = 1;
        if (verbose)
          cerr << "The largest fragment will be retained\n";
      }
      else if ("alarge" == f)
      {
        reduce_to_all_largest_fragments = 1;
        if (verbose)
          cerr << "Will keep all fragments the same size as the largest\n";
      }
      else if ("rmlarge" == f)
      {
        remove_largest_fragment = 1;
        if (verbose)
          cerr << "The largest fragment(s) will be removed\n";
      }
      else if (f.starts_with("rmlarge="))
      {
        f.remove_leading_chars(8);
        if (! f.numeric_value(remove_largest_fragment) || remove_largest_fragment < 1)
        {
          cerr << "The rmlarge= directive must be followed by a whole +ve integer\n";
          return 5;
        }

        if (verbose)
          cerr << "Will remove the " << remove_largest_fragment << " largest fragment(s)\n";
      }
      else if ("lo" == f)
      {
        reduce_to_largest_organic_fragment = 1;

        if (verbose)
          cerr << "Will reduce to the largest organic fragment\n";
      }
      else if ("lod" == f)
      {
        reduce_to_largest_organic_fragment_carefully = 1;

        if (verbose)
          cerr << "Will reduce to largest organic fragment using desirability rules\n";
      }
      else if ("allo" == f)
      {
        keep_all_organic_fragments = true;

        if (verbose)
          cerr << "Will remove all fragments with non organic atoms\n";
      }
      else if (f.starts_with("Q:"))
      {
        f.remove_leading_chars(2);
        if (! fetch_frag_queries(f, "largest", largest_fragment_queries))
          usage(44);
      }
      else if (f.starts_with("q:"))
      {
        f.remove_leading_chars(2);
        if (! fetch_frag_queries(f, "smallest", smallest_fragment_queries))
          usage(45);
      }
      else if (f.starts_with("SMARTS:"))
      {
        f.remove_leading_chars(7);
        if (! fetch_frag_smarts(f, "largest", largest_fragment_queries))
          usage(46);
      }
      else if (f.starts_with("smarts:"))
      {
        f.remove_leading_chars(7);
        if (! fetch_frag_smarts(f, "smallest", smallest_fragment_queries))
          usage(47);
      }
      else if (f.starts_with("ALL:"))
      {
        f.remove_leading_chars(4);
        if (! fetch_frag_smarts(f, "all", keep_fragment_queries))
          usage(32);
      }
      else if (f.starts_with("rm:"))
      {
        f.remove_leading_chars(3);
        if (! fetch_frag_smarts(f, "rm", remove_fragment_queries))
          usage(33);
      }
      else if (f.starts_with("rmle="))
      {
        f.remove_up_to_first('=');
        if (! f.numeric_value(remove_fragments_this_size_or_smaller) || remove_fragments_this_size_or_smaller < 1)
        {
          cerr << "The rmle= qualifier must be a whole positive number\n";
          return 11;
        }

        if (verbose)
          cerr << "Will remove all fragments with " << remove_fragments_this_size_or_smaller << " or fewer atoms\n";
      }
      else if (f.numeric_value(fragment_count))
      {
        if (fragment_count < 1)
        {
          cerr << "The -f <number> directive must be followed by a whole positive number\n";
          usage(8);
        }
        if (verbose)
          cerr << "A maximum of " << fragment_count << " fragments will be written\n";
      }
      else if ("RMDUP" == f)
      {
        remove_duplicate_fragments = 1;

        if (verbose)
          cerr << "Will remove duplicate fragments\n";
      }
      else if (f.starts_with("saltfile="))
      {
        f.remove_leading_chars(9); 
        if (! known_fragment_data.read_known_salts(f))
        {
          cerr << "Cannot read known salts '" << f << "'\n";
          return 9;
        }
      }
      else if (f.starts_with("parentfile="))
      {
        f.remove_leading_chars(11);
        if (! known_fragment_data.read_known_parents(f))
        {
          cerr << "Cannot read known parents '" << f << "'\n";
          return 9;
        }
      }
      else if ("kmfok" == f)
      {
        known_fragment_data.set_only_check_molecular_formula(1);
      }
      else if ("kpallsalt" == f)
      {
        known_fragment_data.set_remove_everything_if_all_fragments_match(0);
        if (verbose)
          cerr << "Will not change molecule consisting of all salts\n";
      }
      else if (f.starts_with("rmxt="))
      {
        f.remove_leading_chars(5);

        if (! f.numeric_value(discard_molecule_if_multiple_fragments_larger_than) || discard_molecule_if_multiple_fragments_larger_than < 1)
        {
          cerr << "The rmxt qualifier must be followed by a whole +ve number\n";
          return 3;
        }

        if (verbose)
          cerr << "Will discard molecules with multiple fragments with more than " << discard_molecule_if_multiple_fragments_larger_than << " atoms\n";
      }
      else if ("rmxt" == f)
      {
        discard_molecule_if_multiple_fragments_larger_than = 16;

        if (verbose)
          cerr << "Will discard molecules with multiple fragments with more than " << discard_molecule_if_multiple_fragments_larger_than << " atoms\n";
      }
      else if (f.starts_with("dmxt="))
      {
        f.remove_leading_chars(5);
        if (! f.numeric_value(mixture_if_largest_frags_differ_by) || mixture_if_largest_frags_differ_by < 0)
        {
          cerr << "The delta for mixture (dmxt=) specifier must be a whole non-negative number\n";
          return 3;
        }

        if (verbose)
          cerr << "A molecule where the largest fragments differ by only " << mixture_if_largest_frags_differ_by << " atoms will be considered a mixture\n";
      }
      else if (f.starts_with("manlf="))
      {
        f.remove_leading_chars(6);
        if (! f.numeric_value(remove_molecules_with_non_largest_fragment_natoms) || remove_molecules_with_non_largest_fragment_natoms < 0)
        {
          cerr << "The max atoms in non-largest fragment directive (manlf) must be a whole +ve number\n";
          return 2;
        }

        if (verbose)
          cerr << "Will discard molecules where a non-largest fragment has more than " << remove_molecules_with_non_largest_fragment_natoms << " atoms\n";
      }
      else if (f.starts_with("klf="))
      {
        f.remove_leading_chars(4);
        if (! f.numeric_value(strip_to_n_largest_fragments) || strip_to_n_largest_fragments < 1)
        {
          cerr << "The strip to N largest fragments 'klf=' directive must have a whole +ve integer\n";
          return 2;
        }

        if (verbose)
          cerr << "Will strip to the " << strip_to_n_largest_fragments << " largest fragments\n";
      }
      else if ("sfs" == f)
      {
        sort_by_fragment_size = 1;

        if (verbose)
          cerr << "Will sort fragments by size\n";
      }
      else if (f.starts_with("RMF="))
      {
        f.remove_leading_chars(4);
        tag_for_removed_fragments = f;

        if (0 == tag_for_removed_fragments.length())
        {
          cerr << "The tag for removed fragments (RMF=) must be non zero length\n";
          return 2;
        }
      }
      else if ("help" == f)
      {
        display_f_directives(0);
      }
      else
      {
        cerr << "Unrecognised -f directive '" << f << "'\n";
        usage(14);
      }
    }

//  We can have only one selection active - except for Known_Fragment_Data which can work in combination
//  with other fragment selections

    int nsel = (fragment_count != 0) + (reduce_to_largest_fragment) + (reduce_to_all_largest_fragments) +
               (reduce_to_largest_organic_fragment) + (reduce_to_largest_organic_fragment_carefully) +
               (keep_all_organic_fragments) + 
               (largest_fragment_queries.number_elements() > 0) +
               (smallest_fragment_queries.number_elements() > 0) +
               (remove_fragment_queries.number_elements() > 0) +
               (strip_to_n_largest_fragments > 0) +
               (keep_fragment_queries.number_elements() > 0) +
               (discard_molecule_if_multiple_fragments_larger_than > 0) +
               (known_fragment_data.active());

//  cerr << "Fragment nsel " << nsel << " remove_fragment_queries " << remove_fragment_queries.size() << ' ' << fragment_count << ' ' << reduce_to_largest_fragment << ' ' << reduce_to_all_largest_fragments << ' ' << reduce_to_largest_organic_fragment << ' ' << reduce_to_largest_organic_fragment_carefully << ' ' << largest_fragment_queries.number_elements() << ' ' << smallest_fragment_queries.number_elements() << ' ' << remove_fragment_queries.number_elements() << ' ' << strip_to_n_largest_fragments << ' ' << keep_fragment_queries.number_elements() << ' ' << known_fragment_data.active() << endl;

    if (1 == nsel)     // good, the easy case
      ;
    else if (remove_fragment_queries.number_elements() > 0 && (reduce_to_largest_fragment || reduce_to_all_largest_fragments || reduce_to_largest_organic_fragment || reduce_to_largest_organic_fragment_carefully || keep_all_organic_fragments))
      ;
    else if (0 == nsel && sort_by_fragment_size)
      ;
    else if (discard_molecule_if_multiple_fragments_larger_than > 0 && 2 == nsel)
      ;
    else
      cerr << "You have specified more than one fragment selection criterion " << nsel << " beware of problems...\n";

    if (nsel || known_fragment_data.active() || remove_molecules_with_non_largest_fragment_natoms >= 0 || sort_by_fragment_size)
    {
      need_to_call_process_fragments = 1;
    }

    if (tag_for_removed_fragments.length() && ! reduce_to_largest_organic_fragment)
    {
      cerr << "Sorry, the tag for removed fragments option 'RMF=' only works with the 'lo' fragment selection qualifier\n";
      return 1;
    }

//  if (known_fragment_data.active())
//    known_fragment_data.debug_print (cerr);
  }

  if (cl.option_present('F'))
  {
    int i = 0;
    const_IWSubstring f;
    while (cl.value('F', f, i++))
    {
      if (f.starts_with("mnsz="))
      {
        f.remove_leading_chars(5);

        if (! f.numeric_value(min_size_max_fragment_count) || min_size_max_fragment_count < 1)
        {
          cerr << "The min size for the too many fragments (-F mnsz=) combination must be a whole +ve number\n";
          usage(3);
        }

        if (verbose)
          cerr << "For max fragment count, will only consider fragments with > " << min_size_max_fragment_count << " atoms\n";
      }
      else if ("help" == f)
      {
        display_dash_F_options(cerr);
      }
      else if (f.numeric_value(maximum_fragment_count))
      {
        if (maximum_fragment_count < 1)
        {
          cerr << "The maximum fragment count (-F) must be a whole +ve number\n";
          usage(3);
        }
        if (verbose)
          cerr << "Molecules containing more than " << maximum_fragment_count << 
                  " fragments will be skipped\n";
      }
      else
      {
        cerr << "Unrecognised -F qualitifier '" << f << "'\n";
        return 3;
      }
    }

    need_to_call_process_fragments = 1;
  }

// The -a option is special because it suppresses all filtering and output

  if (! cl.option_present('a'))
    ;
  else if (cl.option_present('c') || cl.option_present('C') || cl.option_present('I') || cl.option_present('s') || cl.option_present('S'))
  {
    cerr << "Filtering and output options like -c, -C -I -s and -S don't make sense with the -a option\n";
    return 5;
  }

  if (cl.option_present('c'))
  {
    const_IWSubstring c = cl.string_value('c');

//  if (c.starts_with("LF:"))
//  {
//    c.remove_leading_chars(3);
//    lower_atom_count_cutoff_applies_to_largest_fragment = 1;
//  }

    if (! c.numeric_value(lower_atom_count_cutoff))
    {
      cerr << "The lower atom count cutoff (-c) must be a whole +ve number\n";
      usage(2);
    }
    else if (lower_atom_count_cutoff < 0)
    {
      cerr << "The lower atom count cutoff (-c) must be a whole +ve number\n";
      usage(2);
    }
    else if (verbose)
      cerr << "Will exclude molecules with fewer than " << lower_atom_count_cutoff << " atoms\n";
  }

  if (cl.option_present('C'))
  {
    int i = 0;
    const_IWSubstring c;
    while (cl.value('C', c, i++))
    {
      if ("implicit" == c)
      {
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

      if (! c.numeric_value(upper_atom_count_cutoff) || upper_atom_count_cutoff < lower_atom_count_cutoff)
      {
        cerr << "The upper atom count cutoff must be a whole positive number >= " << lower_atom_count_cutoff << endl;
        cerr << "'" << c << "' is invalid\n";
        usage(5);
      }
      else if (verbose)
        cerr << "Will exclude molecules with more than " << upper_atom_count_cutoff << " atoms\n";
    }
  }

  if (! number_assigner.initialise(cl, 'n', verbose))
  {
    cerr << "Cannot process -n option\n";
    usage(51);
  }

  if (cl.option_present('r'))
  {
    if (audit_input)
      cerr << "Warning, ring cutoff with -a does not do anything\n";

    if (! cl.value('r', lower_ring_count_cutoff) ||
          lower_ring_count_cutoff < 1)
    {
      cerr << "-r option needs a whole number > 0\n";
      usage(52);
    }

    if (verbose)
      cerr << "Molecules containing fewer than " << lower_ring_count_cutoff <<
              " will be ignored\n";
  }

#ifdef OLD_DASH_R_FUNCTION
  if (cl.option_present('R'))
  {
    if (audit_input)
      cerr << "Warning, ring cutoff with -a does not do anything\n";

    if (! cl.value('R', upper_ring_count_cutoff) ||
          upper_ring_count_cutoff < lower_ring_count_cutoff)
    {
      cerr << "-R option needs a whole number > " << lower_ring_count_cutoff << endl;
      usage(53);
    }

    if (verbose)
      cerr << "Molecules containing more than " << upper_ring_count_cutoff <<
              " rings will be ignored\n";
  }
#endif

  if (cl.option_present('R'))
  {
    if (audit_input)
      cerr << "Warning, ring cutoff with -a does not do anything\n";

    const_IWSubstring r;

    for (int i = 0; cl.value('R', r, i); i++)
    {
      if (r.starts_with("S:"))
      {
        r.remove_leading_chars(2);
        if (! r.numeric_value(max_ring_systems_allowed) || max_ring_systems_allowed < 0)
        {
          cerr << "The max ring systems allowed (-R S:) option must be a whole non negative number\n";
          usage(2);
        }

        if (verbose)
          cerr << "Will ignore molecules containing more than " << max_ring_systems_allowed << " ring systems\n";
      }
      else if (r.starts_with("MRS:"))
      {
        r.remove_leading_chars(4);
        if (! r.numeric_value(max_rings_in_a_ring_system) || max_rings_in_a_ring_system < 1)
        {
          cerr << "The Max Rings in a System (MRS:) directive must be followed by a whole +ve number\n";
          return 3;
        }

        if (verbose)
          cerr << "Will discard any molecule with more than " << max_rings_in_a_ring_system << " rings in a system\n";
      }
      else if (upper_ring_count_cutoff > 0)
      {
        cerr << "Already specified upper ring count cutoff (-R)\n";
        return 2;
      }
      else if (! r.numeric_value(upper_ring_count_cutoff) || upper_ring_count_cutoff < 1)
      {
        cerr << "The upper ring count cutoff (-R) option must be a whole +ve number\n";
        usage(2);
      }
      else if (verbose)
      cerr << "Molecules containing more than " << upper_ring_count_cutoff <<
              " rings will be ignored\n";
    }
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', upper_ring_size_cutoff) || upper_ring_size_cutoff < 3)
    {
      cerr << "The upper ring size cutoff value (-m) must be a valid ring size\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will discard molecules having ring sizes > " << upper_ring_size_cutoff << endl;
  }

  if (cl.option_present('h'))
  {
    int i = 0;
    const_IWSubstring h;
    while (cl.value('h', h, i++))
    {
      if ('*' == h || "all" == h)
      {
        make_all_implicit_hydrogens_explicit = 1;
        if (verbose)
          cerr << "All implicit hydrogens make explicit\n";
      }
      else
      {
        Substructure_Query * q = new Substructure_Query;
        if (! q->read(h))
        {
          cerr << "Cannot create -h query from '" << h << "'\n";
          delete q;
          return 32;
        }

        atoms_for_implicit_hydrogens.add(q);
      }
    }

    if (make_all_implicit_hydrogens_explicit && atoms_for_implicit_hydrogens.number_elements())
    {
      cerr << "Cannot specify both '-h query' and '-h all' is impossible\n";
      usage(72);
    }
  }

  if (cl.option_present('V'))
  {
    skip_molecules_with_abnormal_valences = 1;
    if (verbose)
      cerr << "Molecules containing abnormal valences will be skipped\n";
  }

  if (cl.option_present('I'))
  {
    int i = 0;
    IWString tmp;
    while (cl.value('I', tmp, i++))
    {
      if ('0' == tmp || tmp.starts_with("excl") || "none" == tmp)
      {
        exclude_isotopes = 1;
        if (verbose)
          cerr << "Molecules containing isotopes will be excluded\n";
      }
      else if ('1' == tmp || tmp.starts_with("incl"))
      {
        exclude_isotopes = 0;
        if (verbose)
          cerr << "No action taken on molecules containing isotopes\n";
      }
      else if ("convert" == tmp || "change" == tmp)
      {
        convert_isotopes = 1;
        if (verbose)
          cerr << "All isotopic atoms will be converted to their non isotopic form\n";
      }
      else if (tmp.starts_with("convert=") || tmp.starts_with("change="))
      {
        if (tmp.starts_with("convert="))
          tmp.remove_leading_chars(8);
        else
          tmp.remove_leading_chars(7);

        int a;       // the starting isotope
        int b = 0;   // the new isotopic value
        if (tmp.contains(','))
        {
          const_IWSubstring s1, s2;
          tmp.split(s1, ',', s2);
          if (! s1.numeric_value(a) || a < 0 || ! s2.numeric_value(b) || b < 0)
          {
            cerr << "The isotopic conversion directive must contain a valid set of isotopes '" << tmp << "' invalid\n";
            display_dash_I_options('I', cerr);
          }
        }
        else
        {
          if (! tmp.numeric_value(a) || a < 0)
          {
            cerr << "The isotopic conversion directive must contain a valid set of isotopes '" << tmp << "' invalid\n";
            display_dash_I_options('I', cerr);
          }
        }

        if (verbose)
          cerr << "Will convert all isotope " << a << " values to isotope " << b << "\n";

        convert_specific_isotopes.add(a);
        convert_specific_isotopes_new_isotope.add(b);
      }
      else if (tmp.starts_with("smarts=") || tmp.starts_with("smarts:"))     // must be of the form 'smarts,new'
      {
        tmp.remove_leading_chars(7);

        const int ndx = tmp.rindex(',');

        if (ndx < 0)
        {
          cerr << "THe smarts for converting isotopes must have be 'smarts,new_isotope', '" << tmp << "' invalid\n";
          display_dash_I_options('I', cerr);
        }

        const_IWSubstring s;
        tmp.from_to(0, ndx - 1, s);

        Substructure_Query * q = new Substructure_Query();

        if (! q->create_from_smarts(s))
        {
          cerr << "Invalid smarts in isotopic directives '" << tmp << "'\n";
          return 1;
        }

        tmp.from_to(ndx+1, tmp.length() - 1, s);
        int x;
        if (! s.numeric_value(x) || x < 0)
        {
          cerr << "Invalid isotope in isotopic directives '" << tmp << "' '" << s << "'\n";
          delete q;
          return 1;
        }

        convert_specific_isotopes_query.add(q);
        convert_specific_isotopes_query_new_isotope.add(x);
      }
      else if (tmp.starts_with("alliso="))
      {
        tmp.remove_leading_chars(7);
        if (! tmp.numeric_value(convert_all_isotopes_to) || convert_all_isotopes_to < 0)
        {
          cerr << "The convert all isotopes 'alliso=' directive must have a valid isotopic value, '" << tmp << "' invalid\n";
          display_dash_I_options('I', cerr);
        }

        if (verbose)
          cerr << "All existing isotopes converted to " << convert_all_isotopes_to << endl;
      }
      else if ("remove" == tmp)
      {
        remove_isotopic_atoms = 1;
        if (verbose)
          cerr << "All isotopically atoms will be removed\n";
      }
      else if ("help" == tmp)
      {
        display_dash_I_options('I', cerr);
      }
      else
      {
        cerr << "Unrecognised -I qualifier '" << tmp << "'\n";
        usage(55);
      }
    }
  }

  if (cl.option_present('B'))
  {
    int i = 0;
    const_IWSubstring b;
    while (cl.value('B', b, i++))
    {
      if (b.starts_with("log="))
      {
        b.remove_leading_chars(4);
        connection_table_error_file = b;

        if (verbose)
          cerr << "Will log connection table errors to '" << connection_table_error_file << "'\n";
      }
      else if ("help" == b)
      {
        display_dash_b_options(cerr);
      }
      else if (b.numeric_value(connection_table_errors_allowed) && connection_table_errors_allowed >= 0)
      {
        if (verbose)
          cerr << connection_table_errors_allowed << " connection table errors allowed\n";
      }
      else
      {
        cerr << "Invalid -B qualifier '" << b << "'\n";
        return 6;
      }
    }
  }

  if (cl.option_present('T'))
  {
    if (! get_translations(cl, verbose, 'T'))
      usage(87);
  }

  if (cl.option_present('p'))
  {
    int i = 0;
    const_IWSubstring p;
    while (cl.value('p', p, i++))
    {
      if ("AMW" == p)
      {
        append_molecular_weight_to_name = 1;
        if (verbose)
          cerr << "The molecular weight will be appended to the name\n";

        appends_to_be_done = 1;
      }
      else if ("AMWI" == p)
      {
        append_molecular_weight_ok_isotope_to_name = 1;
        if (verbose)
          cerr << "The isotopically aware molecular weight will be appended to the name\n";

        appends_to_be_done = 1;
      }
      else if ("MF" == p)
      {
        append_molecular_formula_to_name = 1;
        if (verbose)
          cerr << "The molecular formula will be appended to the name\n";

        appends_to_be_done = 1;
      }
      else if ("ISISMF" == p)
      {
        append_isis_molecular_formula_to_name = 1;
        if (verbose)
          cerr << "The ISIS-like molecular formula will be appended to the name\n";

        appends_to_be_done = 1;
      }
      else if ("NATOMS" == p)
      {
        append_natoms_to_name = 1;
        if (verbose)
          cerr << "The number of atoms will be appended to the name\n";

        appends_to_be_done = 1;
      }
      else if ("NBONDS" == p)
      {
        append_nbonds_to_name = 1;
        if (verbose)
          cerr << "The number of bonds will be appended to the name\n";

        appends_to_be_done = 1;
      }
      else if ("NRINGS" == p)
      {
        append_nrings_to_name = 1;
        if (verbose)
          cerr << "The number of rings will be appended to the name\n";

        appends_to_be_done = 1;
      }
      else if ("AROMR" == p)
      {
        append_aromatic_ring_count_to_name = 1;
        if (verbose)
          cerr << "The number of aromatic rings will be appended to the name\n";

        appends_to_be_done = 1;
      }
      else if ("EXACT" == p)
      {
        append_exact_mass_to_name = 1;
        if (verbose)
          cerr << "The exact mass will be appended to the name\n";

        appends_to_be_done = 1;
      }
      else if ("NFC" == p)
      {
        append_net_formal_charge_to_name = 1;
        if (verbose)
          cerr << "The net formal charge will be appended to the name\n";

        appends_to_be_done = 1;
      }
      else if ("CLND" == p)
      {
        append_clnd_count_to_name = 1;
        if (verbose)
          cerr << "The CLND Nitrogen count will be appended to the name\n";

        appends_to_be_done = 1;
      }
      else if ("PREPEND" == p)
      {
        do_appends_as_prepends = 1;
        if (verbose)
          cerr << "appends done as prepends\n";
      }
      else if (p.starts_with("LARGE"))
      {
        appended_properties_from_largest_fragment = 1;
        if (verbose)
          cerr << "Appended properties derived from largest fragment\n";
      }
      else if (p.starts_with("HETERO") || p.starts_with("HTROAC"))
      {
        append_heteratom_count_to_name = 1;
        if (verbose)
          cerr << "Will append the heteroatom count to the name\n";

        appends_to_be_done = 1;
      }
      else if ("help" == p)
      {
        display_p_directives(0);
      }
      else
      {
        cerr << "Unrecognised -p qualifier '" << p << "'\n";
        display_p_directives(3);
      }
    }
  }

  if (cl.option_present('Y'))
  {
    int i = 0;
    const_IWSubstring y;
    while (cl.value('Y', y, i++))
    {
      if (y.starts_with("appchg="))
      {
        molecule_changed_string = y;
        molecule_changed_string.remove_leading_chars(7);

        if (verbose)
          cerr << "Will append '" << molecule_changed_string << "' to molecules whose connection table is changed\n";
      }
      else if ("pblen" == y)
      {
        print_bond_lengths = 1;
        if (verbose)
          cerr << "Bond lengths will be printed\n";
      }
      else if ("pbang" == y)
      {
        print_bond_angles = 1;
        if (verbose)
          cerr << "Bond angles will be printed\n";
      }
      else if ("ptor" == y)
      {
        print_torsions = 1;
        if (verbose)
          cerr << "Torsions be printed\n";
      }
      else if ("dbg" == y)
      {
        debug_print_each_molecule++;
      }
      else if (y.starts_with("namerx="))
      {
        y.remove_leading_chars(7);
        if (! name_rx.set_pattern(y))
        {
          cerr << "Invalid name regular expression '" << y << "'\n";
          return 8;
        }

        if (verbose)
          cerr << "Will only process molecules whose name matches '" << name_rx.source() << "'\n";
      }
      else if ("ftn" == y)
      {
        truncate_names_to_first_token = 1;

        if (verbose)
          cerr << "Will truncate molecule names to the first token\n";
      }
      else if (y.starts_with("tfirst="))
      {
        y.remove_leading_chars(7);
        truncate_name_at_first = y[0];
      }
      else if (y.starts_with("tlast="))
      {
        y.remove_leading_chars(6);
        truncate_name_at_last = y[0];
      }
      else if (y.starts_with("chname=") || y.starts_with("CHNAME="))
      {
        if ('C' == y[0])    // uppercase form
          change_name_rx_must_match = 1;

        y.remove_leading_chars(7);
        if (! change_name_rx.set_pattern(y))
        {
          cerr << "Invalid change name regexp '" << y << "'\n";
          return 3;
        }

        if (0 == change_name_rx.number_subexpressions())
        {
          cerr << "The chname regular expression must have subexpression(s) '" << change_name_rx.source() << "' invalid\n";
          return 3;
        }

        if (verbose)
          cerr << "Will change molecule names to part that matches rexexp '" << change_name_rx.source() << "'\n";
      }
      else if (y.starts_with("nsubw="))
      {
        y.remove_leading_chars(6);
        substitute_for_whitespace_in_name = y;

        if (1 != substitute_for_whitespace_in_name.length())
        {
          cerr << "Sorry, the substitute for whitespace value must be a single character only '" << y << "' is invalid\n";
          return 3;
        }

        if (verbose)
          cerr << "Whitespace in names changed to '" << substitute_for_whitespace_in_name << "'\n";
      }
      else if (y.starts_with("NPREPEND="))
      {
        y.remove_leading_chars(9);
        prepend_to_name = y;
        if (verbose)
          cerr << "Will prepend '" << prepend_to_name << "' to all names\n";
      }
      else if (y.starts_with("maxplen="))
      {
        y.remove_leading_chars(7);
        if (! y.numeric_value(max_path_length) || max_path_length < 1)
        {
          cerr << "The maxplen= directive must be followed by a +ve whole number\n";
          return 5;
        }

        if (verbose)
          cerr << "Will discard molecules with longest path > " << max_path_length << "\n";
      }
      else if ("ada" == y)
      {
        change_double_bonds_between_permanent_aromatic_to_single = 1;
        if (verbose)
          cerr << "Will change double bonds between permanent aromatic atoms to single\n";
      }
      else if (y.starts_with("B4F="))
      {
        y.remove_leading_chars(4);
        IWString fname(y);
        if ('-' == y)
          ;
        else if (! fname.ends_with(".smi"))
          fname << ".smi";

        if (! stream_for_smiles_before_filters.open(fname.null_terminated_chars()))
        {
          cerr << "Cannot open stream for pre-filtered smiles '" << fname << "'\n";
          return 3;
        }

        if (verbose)
          cerr << "Pre-filtered smiles written to '" << fname << "'\n";
      }
      else if ("aclf" == y)
      {
        atom_count_includes_only_atoms_in_largest_fragment = 1;

        if (verbose)
          cerr << "All atom counts will be with regard to larges fragment only\n";
      }
      else if ("nbvm" == y)
      {
        set_display_abnormal_valence_messages(0);
        set_display_strange_chemistry_messages(0);
        if (verbose)
          cerr << "Abnormal valence messages suppressed\n";
      }
      else if ("okbvi" == y)
      {
        ok_bad_valence_on_isotopically_labelled = 1;
        skip_molecules_with_abnormal_valences = 1;
        if (verbose)
          cerr << "During valence check will ignore bad valences on isotopically labelled atoms\n";
      }
      else if ("nhsqb" == y)
      {
        set_explicit_hydrogens_need_square_brackets_in_smiles(0);
        if (verbose)
          cerr << "Explicit hydrogens in smiles written without square brackets\n";
      }
      else if ("rmsqb" == y)
      {
        remove_unnecessary_square_brackets = 1;
        if (verbose)
          cerr << "Will remove unnecessary square bracketed elements\n";
      }
      else if ("FHrmsqb" == y)
      {
        remove_all_possible_square_brackets = 1;
        if (verbose)
          cerr << "Will assume all H counts are OK and remove smiles square brackets\n";
      }
      else if ("rmamap" == y)
      {
        reset_atom_map_numbers = 1;
        if (verbose)
          cerr << "Will not include atom map with smiles\n";
      }
      else if (y.starts_with("fixarom="))
      {
        y.remove_leading_chars(8);

        Substructure_Query * q = new Substructure_Query;
        if (! q->create_from_smarts(y))
        {
          cerr << "Invalid smarts specification of atom in ring to be fixed '" << y << "'\n";
          delete q;
          return 1;
        }

        aromatise_these_rings.add(q);
      }
      else if ("nbrminv" == y)
      {
        set_invalidate_bond_list_ring_info_during_invalidate_ring_info(0);
      }
      else if ("iso2amap" == y)
      {
        convert_isotopes_to_atom_map_numbers = 1;
        if (verbose)
          cerr << "Will convert isotopes to atom map numbers\n";
      }
      else if ("amap2iso" == y)
      {
        convert_atom_map_numbers_to_isotopes = 1;
        if (verbose)
          cerr << "Will convert atom map numbers to isotopes\n";
      }
      else if ("help" == y)
      {
        display_dash_y_options(cerr, 'Y', 2);
      }
      else
      {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        return 6;
      }
    }
  }

  Molecule_Output_Object output;

  if (! audit_input)    // then we must get our output types
  {
    if (! cl.option_present('o'))
    {
      output.add_output_type(SMI);

      if (verbose)
        cerr << "Default output type smiles\n";
    }
    else if (! output.determine_output_types(cl))
    {
      cerr << "Cannot determine output types\n";
      usage(8);
    }

    if (cl.option_present('G'))
    {
      output.set_molecules_per_file(1);
      if (verbose)
        cerr << "Each molecule written to its own file (-G option)\n";
    }

    if ("-" == output_file_stem)
      ;
    else if (output.would_overwrite_input_files(cl, output_file_stem))
    {
      cerr << "fileconv:cannot overwrite input file(s)\n";
      return 4;
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << prog_name << ": insufficient arguments " << argc << "\n";
    usage(56);
  }

  if (cl.option_present('L'))
  {
    if (! cl.option_present('o'))
      rejections_output_object.add_output_type(SMI);
    else if (! rejections_output_object.determine_output_types(cl))
    {
      cerr << "Cannot discern output types for rejections file\n";
      return 2;
    }

    IWString reject_log_file_name;
    cl.value('L', reject_log_file_name);

    if (rejections_output_object.would_overwrite_input_files(cl, reject_log_file_name))
    {
      cerr << "Reject file '" << reject_log_file_name << "' cannot overwrite input file(s)\n";
      return 0;
    }

    if (! rejections_output_object.new_stem(reject_log_file_name, 1))
    {
      cerr << "Rejections file cannot use stem '" << reject_log_file_name << "'\n";
      return 3;
    }

    reject_log.open(reject_log_file_name.null_terminated_chars(), std::ios::out);
    if (! reject_log.good())
    {
      cerr << "Cannot open reject file '" << reject_log_file_name << "'\n";
      return 2;
    }
    if (verbose)
      cerr << "Rejected structures will be written to '" << reject_log_file_name << "'\n";
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    const char *fname = cl[i];

    if (verbose)
      cerr << "Processing '" << fname << "'\n";

    int myrc;

    if (::strlen(fname) > 2 && 0 == strncmp("F:", fname, 2))
      myrc = fileconv_list_of_files(fname, input_type, output);
    else
      myrc = fileconv(fname, input_type, output);

    if (! myrc)
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    if (cl.number_elements() > 1)
      cerr << molecules_read << " molecules read, ";
    if (! audit_input)
      cerr << molecules_written << " molecules written";
    if (molecules_changed)
      cerr << ' ' << molecules_changed << " molecules changed";
    cerr << endl;

    for (int i = 0; i < initial_fragment_count.number_elements(); i++)
    {
      if (initial_fragment_count[i])
        cerr << initial_fragment_count[i] << " molecules had " << i << " fragments\n";
    }

    if (fragment_count && molecules_chopped)
    {
      cerr << molecules_chopped << " molecules chopped to ";
      if (reduce_to_largest_fragment)
        cerr << "largest";
      else
        cerr <<  fragment_count;
      cerr << " fragments\n";

      cerr << "Molecules lost between " << atoms_lost.minval() << " and " << atoms_lost.maxval() << " atoms\n";
    }

    if (natoms_accumulator.n() > 0)
    {
      cerr << natoms_accumulator.n() << " molecules had between " << natoms_accumulator.minval() << " and " <<
              natoms_accumulator.maxval() << " atoms.";
      if (natoms_accumulator.n() > 1)
        cerr << " Average " << natoms_accumulator.average() << endl;
      else
        cerr << endl;
    }

    if (verbose > 1)
    {
      for (int i = 0; i < atom_count.number_elements(); i++)
      {
        if (atom_count[i])
          cerr << atom_count[i] << " molecules had " << i << " atoms\n";
      }
    }

    if (output_organic_only && non_organic_molecules_found)
      cerr << "Skipped " << non_organic_molecules_found << 
        " molecules containing non organic atoms\n";
    if (non_real_elements_found)
      cerr << "Skipped " << non_real_elements_found << " molecules containing non-periodic table elements\n";
    if (molecules_excluded_for_non_allowed_elements)
      cerr << "Skipped " << molecules_excluded_for_non_allowed_elements << " molecules with non-allowed elements\n";
    if (molecules_below_molecular_weight_cutoff)
      cerr << "Skipped " << molecules_below_molecular_weight_cutoff <<
              " molecules with molecular wieght below " << lower_molecular_weight_cutoff << endl;
    if (molecules_above_molecular_weight_cutoff)
      cerr << "Skipped " << molecules_below_molecular_weight_cutoff <<
              " molecules with molecular wieght above " << lower_molecular_weight_cutoff << endl;
    if (molecules_below_atom_count_cutoff)
      cerr << "Skipped " << molecules_below_atom_count_cutoff << 
              " molecules with atom count less than " << lower_atom_count_cutoff << endl;
    if (molecules_above_atom_count_cutoff)
      cerr << "Skipped " << molecules_above_atom_count_cutoff << 
              " molecules with atom count greater than " << upper_atom_count_cutoff << endl;
    if (molecules_with_too_many_components)
      cerr << "Skipped " << molecules_with_too_many_components <<
              " molecules with more than " << maximum_fragment_count << " components\n";
    if (molecules_with_large_fragments)
      cerr << "Skipped " << molecules_with_large_fragments << " with non-largest fragments with more than " << remove_molecules_with_non_largest_fragment_natoms << " atoms\n";
    if (molecules_with_fragments_reduced_by_query)
      cerr << molecules_with_fragments_reduced_by_query << 
              " molecules chopped to fragment which matches query\n";
    if (molecules_not_matching_fragment_queries)
      cerr << "Skipped " << molecules_not_matching_fragment_queries << 
              " molecules which did not match fragment query specifications\n";
    if (remove_fragments_this_size_or_smaller)
      cerr << molecules_with_very_small_fragments_removed << " molecules lost fragments with " << remove_fragments_this_size_or_smaller << " or fewer atoms\n";
    if (remove_duplicate_fragments)
      cerr << molecules_with_duplicate_fragments_removed << " molecules had duplicate fragments removed\n";
    if (keep_all_organic_fragments)
      cerr << molecules_with_nonorganic_fragments_removed << " molecules had nonorganic fragments removed\n";
    if (discard_molecule_if_multiple_fragments_larger_than)
      cerr << molecules_with_multiple_large_fragments << " molecules had multiple fragments with at least " << discard_molecule_if_multiple_fragments_larger_than << endl;
    if (mixture_if_largest_frags_differ_by >= 0)
      cerr << molecules_declared_mixtures_by_atom_count_difference << " molecules declared mixtures by atom count difference\n";
    if (molecules_with_too_few_rings)
      cerr << "Skipped " << molecules_with_too_few_rings << 
              " molecules having fewer than " << lower_ring_count_cutoff << " rings\n";
    if (molecules_with_too_many_rings)
      cerr << "Skipped " << molecules_with_too_many_rings << 
              " molecules having more than " << upper_ring_count_cutoff << " rings\n";
    if (molecules_with_large_rings)
      cerr << "Skipped " << molecules_with_large_rings << " molecules containing rings larger than " << upper_ring_size_cutoff << " atoms\n";
    if (molecules_with_ring_systems_too_large)
      cerr << "Skipped " << molecules_with_ring_systems_too_large << " molecules containing ring systems with more than " << max_rings_in_a_ring_system << " rings\n";
    if (molecules_with_too_many_ring_systems)
      cerr << "Skipped " << molecules_with_too_many_ring_systems << " with too many ring systems\n";
    if (max_path_length)
      cerr << "Skipped " << molecules_with_longest_path_too_long << " molecules with max path above " << max_path_length << '\n';

    if (amw_accumulator.n() > 1)
    {
      cerr << amw_accumulator.n() << " molecular weights. Between " <<
              amw_accumulator.minval() << " and " << amw_accumulator.maxval() <<
              " ave = " << amw_accumulator.average() << 
              " var = " << amw_accumulator.variance() << endl;
    }

    if (molecules_below_amw_cutoff)
      cerr << "Excluded " << molecules_below_amw_cutoff << " molecules with amw below " <<
              lower_amw_cutoff << endl;
    if (molecules_above_amw_cutoff)
      cerr << "Excluded " << molecules_above_amw_cutoff << " molecules with amw above " <<
              upper_amw_cutoff << endl;
    
    if (molecules_with_chiral_centres)
      cerr << molecules_with_chiral_centres << " molecules had chiral atoms\n";
    if (chiral_centres_inverted)
      cerr << chiral_centres_inverted << " chiral centres inverted\n";

    if (molecules_with_chiral_centres_removed_by_rmchiral)
      cerr << "Removed " << chiral_centres_removed_by_rmchiral << " chiral centres from " << molecules_with_chiral_centres_removed_by_rmchiral << " molecules per rmchiral\n";

    if (molecules_containing_isotopes)
      cerr << molecules_containing_isotopes << " molecules contained isotopic atoms that were processed\n";

    elements_to_remove.report(cerr);

    if (element_transformations.number_elements())
      element_transformations.debug_print(cerr);

    if (structure_fixing.active())
      structure_fixing.report(cerr);

    if (molecules_with_abnormal_valences)
      cerr << molecules_with_abnormal_valences << " molecules containing abnormal valences\n";

    if (chemical_standardisation.active())
      chemical_standardisation.report(cerr);

    if (charge_assigner.active())
      charge_assigner.report(cerr);

    if (make_all_implicit_hydrogens_explicit || atoms_for_implicit_hydrogens.number_elements())
      cerr << molecules_to_which_hydrogens_were_added << " molecules had implicit hydrogens made explicit\n";

    if (aromatise_these_rings.number_elements())
      cerr << molecules_changed_by_aromatising_rings << " molecules made aromatic rings\n";

    if (max_chiral_centres > 0)
    {
      cerr << molecules_with_too_many_chiral_centres << " molecules had more than " << max_chiral_centres << endl;
      for (int i = 0; i < chiral_centre_count.number_elements(); ++i)
      {
        if (chiral_centre_count[i] > 0)
          cerr << chiral_centre_count[i] << " molecules had " << i << " chiral centres\n";
      }
    }

  }

  return rc;
}

}
//#define DLL_FLAG=1
#ifndef DLL_FLAG

int
main (int argc, char **argv)
{
  prog_name = argv[0];

  int rc = Fileconv::fileconv(argc, argv);

  return rc;
}

#endif

#ifdef HZ_CSHARP

extern "C" int
fileconv_csharp(int argc, char **argv)
{
  return Fileconv::fileconv(argc,argv);
}

#endif

