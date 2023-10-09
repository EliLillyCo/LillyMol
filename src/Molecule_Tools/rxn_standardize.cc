/*
  standardizes a reaction
*/

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mdl_molecule.h"
#include "Molecule_Lib/rxn_file.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

static int reactions_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_product = 0;

static int reduce_to_largest_reactant = 0;

static int eliminate_reagents_not_participating = 0;

static int remove_fragments_not_participating = 0;

static int ignore_bad_reactions = 0;

static int bad_reactions_ignored = 0;

static IWString_and_File_Descriptor stream_for_discarded_reactions;

static int skip_reactions_with_multiple_reactants = 0;

static int skip_reactions_with_no_common_mapped_atoms = 0;

static int reactions_with_multiple_reactants_skipped = 0;

static int skip_reactions_with_multiple_products = 0;

static int reactions_with_multiple_products_skipped = 0;

static int auto_fix_orphans = 0;
static int auto_fix_widows = 0;

static int discard_reactions_with_isotopes = 0;

static int reactions_discarded_for_isotopic_atoms = 0;

static int remove_duplicate_reactants = 0;
static int remove_duplicate_products = 0;
static int remove_duplicate_agents = 0;
static int remove_all_agents = 0;
static int remove_unmapped_components = 0;
static int remove_unchanged_components = 0;

static int use_first_token_of_name = 0;
static int gsub_reaction_names = 0;
static int compress_consecutive_underscores = 0;

// static int take_first_of_multiple_reagents = 0;

static int fix_kekule_problems = 0;
static int reactions_with_kekule_problems_fixed = 0;

static int reactions_with_fragments_not_participating = 0;
static int reactions_with_reagents_not_participating = 0;
static int reactions_with_duplicate_reagents = 0;
static int reactions_with_reagent_count_changed = 0;
static int reactions_with_no_common_mapped_atoms = 0;

static int reactions_with_orphan_atoms = 0;
static int reactions_with_counterions_moved_to_orphan = 0;
static int reactions_with_duplicates_atom_maps_ignored = 0;
static int reactions_with_no_reagents = 0;
static int reactions_with_no_products = 0;

static int reactions_written = 0;

static bool use_all_plusses = false;
static bool use_all_dots = false;

static Accumulator_Int<int> acc_natoms;

static int max_atoms_in_reagent_and_product = 0;

static int reactions_discarded_for_too_many_atoms = 0;

static int reading_multi_reaction_rdfile = 0;

static int do_automatic_atom_mapping = 0;

static int discard_reactions_containing_duplicate_atom_map_numbers = 0;
static int reactions_containing_duplicate_atom_map_numbers = 0;

static int unmap_duplicate_atom_map_numbers = 0;

static int assign_unmapped_atoms = 0;

static int skip_reactions_where_largest_fragment_is_unchanged = 0;

static int skip_reactions_containing_aromatic_bonds = 0;

static int reactions_containing_aromatic_bonds = 0;

static int reactions_where_largest_fragment_does_not_change = 0;

// static int input_is_reaction_smiles = 1;

static void
usage() {
  // clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "  -R            reading multi-reaction RDfile\n";
  cerr << "  -s            reading a reaction smiles (default is -s)\n";
  cerr << "  -c            discard chirality on input\n";
  cerr << "  -K            fix alternating Kekule issues\n";
  cerr << "  -o            create reactant fragments that are orphans\n";
  cerr << "  -d d          separate components with . rather than +\n";
  cerr << "  -d p          separate components with + rather than .\n";
  cerr << "  -C <natoms>   discard any reaction where the largest reagent has more than <natoms> atoms\n";
  cerr << "  -f trunc      truncate reaction names to the first token\n";
  cerr << "  -f gsub       replace unusual characters in reaction names with _\n";
  cerr << "  -f __         if the name is changed, compress consecutive underscores\n";
  cerr << "  -X ...        options related to discarding reactions, enter -X help for help\n";
  cerr << "  -D ...        options related to discarding reaction components, enter -D help for help\n";
  cerr << "  -U <fname>    write discarded reactions to <fname>\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -M            generate atom maps\n";
  cerr << "  -q            quiet; suppress warnings about no mtched atoms\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(1);
}

static void
process_reaction_name(IWString& s) {
  if (!gsub_reaction_names && !use_first_token_of_name) {  // nothing to do
    return;
  }

  // cerr << "Starting - name is : " << s << endl;
  const int n = s.number_elements();

  int contains_whitespace = 0;
  int contains_non_alpha = 0;

  int changes_made = 0;

  for (int i = 0; i < n; ++i) {
    const char c = s[i];
    if (isalnum(c)) {
      continue;
    }

    if ('_' == c || '.' == c || '-' == c) {
      continue;
    }

    if (use_first_token_of_name) {
      s.iwtruncate(i);
      changes_made++;
      // cerr << "Truncated: " << s << endl;
      break;
    } else if (gsub_reaction_names) {
      s[i] = '_';
      changes_made++;
      // cerr << "fixing a char - name is now: " << s << endl;
    } else if (isspace(c)) {
      contains_whitespace++;
    } else {
      contains_non_alpha++;
    }
  }

  // cerr <<"name after fixes:: " << s << endl;
  if (contains_whitespace || contains_non_alpha) {
    cerr << "Warning reaction '" << s << "' whitespace or other characters\n";
  }

  if (!changes_made) {
    return;
  }
  // cerr << "compress_consecutive_underscores: " << compress_consecutive_underscores <<
  // endl;
  if (!compress_consecutive_underscores) {
    return;
  }

  bool previous_was_underscore = false;
  int ndx = 0;
  for (int i = 0; i < n; ++i) {
    // cerr << "name: " << s << endl;
    // cerr << "thisChar: " << s[i] << endl;
    if ('_' == s[i]) {
      if (!previous_was_underscore) {
        s[ndx] = s[i];
        ndx++;
        previous_was_underscore = true;
      }
      // else
      // cerr << "skipping" << endl;
    } else {
      s[ndx] = s[i];
      ndx++;
      previous_was_underscore = false;
    }
  }

  s.iwtruncate(ndx);

  // cerr << "Name set to '" << s << "'\n";

  return;
}

static int
contains_aromatic_bonds(const ISIS_RXN_FILE_Molecule& m) {
  const int nedges = m.nedges();

  for (int i = 0; i < nedges; ++i) {
    const Bond* b = m.bondi(i);
    if (b->is_single_bond()) {
      ;
    } else if (b->is_double_bond()) {
      ;
    } else if (b->is_triple_bond()) {
      ;
    } else if (b->is_aromatic()) {
      return 1;
    }
  }

  return 0;
}

static int
contains_aromatic_bonds(RXN_File& rxn) {
  const int ns = rxn.number_reagents();

  for (int i = 0; i < ns; ++i) {
    if (contains_aromatic_bonds(rxn.reagent(i))) {
      return 1;
    }
  }

  return 0;
}

static int
echo_bad_data(iwstring_data_source& input, const off_t initial_offset,
              IWString_and_File_Descriptor& output) {
  const auto current_offset = input.tellg();
  input.seekg(initial_offset);
  input.echo(output, (current_offset - initial_offset));
  input.seekg(current_offset);
  output.write_if_buffer_holds_more_than(4096);

  return input.tellg() == current_offset;
}

static int
rxn_standardize(RXN_File& rxn, iwstring_data_source& input, const off_t initial_offset,
                IWString_and_File_Descriptor& output) {
  const int initial_nr = rxn.number_reagents();

  if (verbose > 1) {
    cerr << endl
         << "Read '" << rxn.name() << endl
         << "' with " << initial_nr << " reagents\n";
  }

  IWString tmp(rxn.name());
  process_reaction_name(tmp);
  rxn.set_name(tmp);

  if (0 == initial_nr) {
    cerr << "Skipping reaction with no reagents " << rxn.name() << endl;
    reactions_with_no_reagents++;
    if (stream_for_discarded_reactions.is_open()) {
      echo_bad_data(input, initial_nr, stream_for_discarded_reactions);
    }
    return 1;
  }

  if (0 == rxn.number_products()) {
    cerr << "Skipping reaction with no products " << rxn.name() << endl;
    reactions_with_no_products++;
    if (stream_for_discarded_reactions.is_open()) {
      echo_bad_data(input, initial_nr, stream_for_discarded_reactions);
    }
    return 1;
  }

#ifdef SHOW_ATOM_NUMBERED_REACTIONS
  for (int i = 0; i < rxn.number_reagents(); ++i) {
    ISIS_RXN_FILE_Molecule& r = rxn.reagent(i);
    r.set_isotope_to_atom_number_no_perturb_canonical_ordering();
    cerr << r.smiles() << " " << rxn.name() << " reagent " << i << endl;
  }
  for (int i = 0; i < rxn.number_products(); ++i) {
    ISIS_RXN_FILE_Molecule& p = rxn.product(i);
    p.set_isotope_to_atom_number_no_perturb_canonical_ordering();
    cerr << p.smiles() << " " << rxn.name() << " product " << i << endl;
  }

  for (int i = 0; i < rxn.number_reagents(); ++i) {
    ISIS_RXN_FILE_Molecule& r = rxn.reagent(i);
    if (i > 0) {
      cerr << '.';
    }
    cerr << r.smiles();
    r.transform_to_non_isotopic_form();
  }
  cerr << ">>";
  for (int i = 0; i < rxn.number_products(); ++i) {
    ISIS_RXN_FILE_Molecule& p = rxn.product(i);
    if (i > 0) {
      cerr << '.';
    }
    cerr << p.smiles();
    p.transform_to_non_isotopic_form();
  }
  cerr << " " << rxn.name() << endl;
#endif

  // the order of operations is important - For example:  if duplicate
  // reactant are to be removed, then we should do that before checking
  // for duplicated mapping numbers.

  // #define DEBUG_RXN_STANDARDIZE
#ifdef DEBUG_RXN_STANDARDIZE
  cerr << "LINE " << __LINE__ << endl;
  rxn.debug_print(cerr);
#endif

  if (max_atoms_in_reagent_and_product > 0 ||
      verbose)  // do first because it may stop various other failures (hide our
                // problems!)
  {
    const int x = rxn.max_atom_in_any_reagent();
    const int y = rxn.max_atom_in_any_product();
    // cerr << "Reagent atom count: " << x << "    Product atom count: " <<y << endl;
    if (verbose) {
      acc_natoms.extra(x);
    }
    if (max_atoms_in_reagent_and_product > 0 && x > max_atoms_in_reagent_and_product &&
        y > max_atoms_in_reagent_and_product) {
      if (verbose > 1) {
        cerr << rxn.name() << " too many atoms in a reagent and a product " << x << endl;
      }

      reactions_discarded_for_too_many_atoms++;
      if (stream_for_discarded_reactions.is_open()) {
        echo_bad_data(input, initial_offset, stream_for_discarded_reactions);
      }
      return 1;
    }
  }

#ifdef DEBUG_RXN_STANDARDIZE
  cerr << "LINE " << __LINE__ << endl;
  rxn.debug_print(cerr);
#endif

  if (assign_unmapped_atoms) {
    rxn.assign_unmapped_atoms();
  }

  if (unmap_duplicate_atom_map_numbers && rxn.contains_duplicate_atom_map_numbers()) {
    rxn.unmap_duplicate_atom_map_numbers();
  }

#ifdef DEBUG_RXN_STANDARDIZE
  cerr << "LINE " << __LINE__ << endl;
  rxn.debug_print(cerr);
#endif

  if (remove_fragments_not_participating && rxn.remove_fragments_not_participating()) {
    if (verbose > 1) {
      cerr << "Removed fragments not reacting from '" << rxn.name() << ", now "
           << rxn.number_reagents() << " reagents\n";
    }
    reactions_with_fragments_not_participating++;

    if (0 == rxn.number_reagents()) {
      cerr << "After removing fragments not reacting from " << rxn.name()
           << ", no reactants remain\n";
      if (stream_for_discarded_reactions.is_open()) {
        echo_bad_data(input, initial_offset, stream_for_discarded_reactions);
      }
      return 1;
    }
  }

#ifdef DEBUG_RXN_STANDARDIZE
  cerr << "LINE " << __LINE__ << endl;
  rxn.debug_print(cerr);
#endif

  if (eliminate_reagents_not_participating &&
      rxn.eliminate_reagents_not_participating()) {
    if (verbose > 1) {
      cerr << "Removed reagents not reacting from '" << rxn.name() << "', now "
           << rxn.number_reagents() << " reagents\n";
    }

    reactions_with_reagents_not_participating++;

    if (0 == rxn.number_reagents()) {
      cerr << "After removing reactants not reacting from " << rxn.name()
           << ", no reactants remain\n";
      if (stream_for_discarded_reactions.is_open()) {
        echo_bad_data(input, initial_offset, stream_for_discarded_reactions);
      }
      return 1;
    }
  }

#ifdef DEBUG_RXN_STANDARDIZE
  cerr << "LINE " << __LINE__ << endl;
  rxn.debug_print(cerr);
#endif

  //  if (rxn.all_reagents_the_same())
  //  {
  //    if (verbose > 1)
  //      cerr << "Removed duplicate reagents from '" << rxn.name() << "' now " <<
  //      rxn.number_reagents() << " reagents\n";
  //
  //    reactions_with_duplicate_reagents++;
  //  }

#ifdef DEBUG_RXN_STANDARDIZE
  cerr << "LINE " << __LINE__ << endl;
  rxn.debug_print(cerr);
#endif

  if (reduce_to_largest_product && rxn.reduce_to_largest_product()) {
    if (verbose > 1) {
      cerr << "Reduced to the largest product fragment from '" << rxn.name() << "\n";
    }
  }

  if (reduce_to_largest_reactant && rxn.reduce_to_largest_reactant()) {
    if (verbose > 1) {
      cerr << "Reduced to the largest reactant fragment from '" << rxn.name() << "\n";
    }
  }

  if (remove_duplicate_reactants && rxn.remove_duplicate_reagents_ignore_atom_map()) {
    if (verbose > 1) {
      cerr << "Removed duplicate reactants from '" << rxn.name() << "' now "
           << rxn.number_reagents() << " reactants\n";
    }
  }

  if (remove_duplicate_products &&
      rxn.remove_duplicate_products_ignore_atom_map())  // qqqwwwqqq
  {
    if (verbose > 1) {
      cerr << "Removed duplicate products from '" << rxn.name() << "' now "
           << rxn.number_reagents() << " reactants\n";
    }
  }

  if (remove_duplicate_agents && rxn.remove_duplicate_agents()) {
    if (verbose > 1) {
      cerr << "Removed duplicate agents from '" << rxn.name() << "' now "
           << rxn.number_agents() << " agents\n";
    }
  }

  if (remove_all_agents && rxn.remove_all_agents()) {
    if (verbose > 1) {
      cerr << "Removed all agents from '" << rxn.name() << "'\n";
    }
  }

  if (remove_unmapped_components && rxn.remove_unmapped_components())  // qqqwwwqqq
  {
    if (verbose > 1) {
      cerr << "Removed all unmapped components from '" << rxn.name() << "'\n";
    }
  }

  if (remove_unchanged_components && rxn.remove_unchanging_components()) {
    if (verbose > 1) {
      cerr << "Removed all reactants that also appear as products from  '" << rxn.name()
           << "'\n";
    }
  }

#ifdef DEBUG_RXN_STANDARDIZE
  cerr << "LINE " << __LINE__ << endl;
  rxn.debug_print(cerr);
#endif

  if (skip_reactions_with_no_common_mapped_atoms &&
      !rxn.at_least_some_mapped_atoms_common_btw_reagents_and_products()) {
    if (verbose > 1) {
      cerr << "No reagent mapped atoms in products " << rxn.name() << endl;
    }
    reactions_with_no_common_mapped_atoms++;
    return 1;
  }

  if (auto_fix_orphans) {
    rxn.set_auto_fix_orphans(1);  // set this here because we only want to fix the orphans
                                  // after other processing
    if (rxn.check_for_widows_and_orphans()) {
      if (verbose > 1) {
        cerr << rxn.name() << " identified orphan atoms\n";
      }
    }
    rxn.set_auto_fix_orphans(0);  // undo the set to fix orphans
  }

  if (auto_fix_widows) {
    rxn.remove_widow_atoms();
  }

  if (fix_kekule_problems) {
    if (rxn.fix_kekule_differences()) {
      reactions_with_kekule_problems_fixed++;

      if (verbose > 1) {
        cerr << "Fixed Kekule problems in " << rxn.name() << endl;
      }
    }
  }

  // rxn.print_atom_map_into(cerr);

  const int nreagents = rxn.number_reagents();

  if (initial_nr != nreagents) {
    if (verbose > 1) {
      cerr << "Reagent count changed " << initial_nr << " to " << nreagents << " in "
           << rxn.name() << endl;
    }
    reactions_with_reagent_count_changed++;
  }

  if (rxn.contains_duplicate_atom_map_numbers()) {
    reactions_containing_duplicate_atom_map_numbers++;
    if (verbose > 1) {
      cerr << rxn.name() << " contains duplicate atom map numbers\n";
    }
    if (discard_reactions_containing_duplicate_atom_map_numbers) {
      if (stream_for_discarded_reactions.is_open()) {
        echo_bad_data(input, initial_offset, stream_for_discarded_reactions);
      }
      return 1;
    }
  }

  if (skip_reactions_with_multiple_reactants && rxn.number_reagents() > 1) {
    reactions_with_multiple_reactants_skipped++;

    if (stream_for_discarded_reactions.is_open()) {
      echo_bad_data(input, initial_offset, stream_for_discarded_reactions);
    }

    return 1;
  }

  if (skip_reactions_with_multiple_products && rxn.number_products() > 1) {
    reactions_with_multiple_products_skipped++;

    if (stream_for_discarded_reactions.is_open()) {
      echo_bad_data(input, initial_offset, stream_for_discarded_reactions);
    }

    return 1;
  }
  if (skip_reactions_where_largest_fragment_is_unchanged &&
      !rxn.largest_fragment_is_changed()) {
    reactions_where_largest_fragment_does_not_change++;
    if (stream_for_discarded_reactions) {
      echo_bad_data(input, initial_offset, stream_for_discarded_reactions);
    }

    return 1;
  }

  if (skip_reactions_containing_aromatic_bonds && contains_aromatic_bonds(rxn)) {
    reactions_containing_aromatic_bonds++;
    if (stream_for_discarded_reactions) {
      echo_bad_data(input, initial_offset, stream_for_discarded_reactions);
    }

    return 1;
  }

  if (discard_reactions_with_isotopes && rxn.contains_isotopic_reagent_atoms()) {
    if (verbose > 1) {
      cerr << rxn.name() << " contains isotopic atoms\n";
    }
    reactions_discarded_for_isotopic_atoms++;
    if (stream_for_discarded_reactions.is_open()) {
      echo_bad_data(input, initial_offset, stream_for_discarded_reactions);
    }
    return 1;
  }

  if (rxn.contains_orphan_atoms()) {
    if (verbose > 1) {
      cerr << rxn.name() << " contains " << rxn.contains_orphan_atoms()
           << " orphan atoms\n";
    }
    reactions_with_orphan_atoms++;
  }

  Reaction_Smiles_Options opts;

  if (use_all_plusses) {
    opts.set_reagent_product_plus_rather_than_dot(
        true);  // controls char between components
    setUseThisCharAsDotInSmiles(
        '+');  // controls char between fragments of each component
  } else if (use_all_dots) {
    opts.set_reagent_product_plus_rather_than_dot(
        false);  // controls char between components
    setUseThisCharAsDotInSmiles(
        '.');  // controls char between fragments of each component
  } else  // normal, defalt mode - plusses between cdomponents and dots  between fragments
          // of each component
  {
    opts.set_reagent_product_plus_rather_than_dot(
        true);  // controls char between components
    setUseThisCharAsDotInSmiles(
        '.');  // controls char between fragments of each component
  }

  // old way
  // opts.set_reagent_product_plus_rather_than_dot(use_all_plusses);
  // opts.set_orphan_use_all_plusses(orphan_use_all_plusses);

  rxn.write_rxn_smiles(opts, output);

  output.write_if_buffer_holds_more_than(4096);

  reactions_written++;

  return 1;
}

static void
initialise_reaction(RXN_File& rxn) {
  rxn.set_do_automatic_atom_mapping(do_automatic_atom_mapping);
  //  if (auto_fix_orphans)
  //    rxn.set_auto_fix_orphans(1);

  return;
}

// static int
// rxn_standardize(iwstring_data_source & input,
//               IWString_and_File_Descriptor & output)
//{
//   RXN_File rxn;
//
//   initialise_reaction(rxn);
//
//   const auto initial_offset = input.tellg();
//
//   if (rxn.do_read(input))
//     ;
//   else if (input.eof())
//     return 1;
//   else if (ignore_bad_reactions)
//   {
//     bad_reactions_ignored++;
//
//     if (stream_for_discarded_reactions.is_open())
//       echo_bad_data(input, initial_offset, stream_for_discarded_reactions);
//     return 1;
//   }
//   else
//   {
//     cerr << "Fatal error reading reaction, now at line " << input.lines_read() << endl;
//     return 0;
//   }
//
//   reactions_read++;
//
//   return rxn_standardize(rxn, input, initial_offset, output);
// }

static int
rxn_standardize_reaction_smiles(iwstring_data_source& input,
                                IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  auto initial_offset = input.tellg();  // most likely 0

  while (input.next_record(buffer)) {
    RXN_File rxn;

    //  cerr << "Processing " << buffer << endl;

    initialise_reaction(rxn);

    if (!rxn.build_from_reaction_smiles(buffer, 1)) {
      cerr << "Cannot interpret " << buffer << endl;
      if (ignore_bad_reactions) {
        bad_reactions_ignored++;
        continue;
      }

      return 0;
    }

    reactions_read++;

    if (!rxn_standardize(rxn, input, initial_offset, output)) {
      cerr << "Cannot process " << buffer << endl;
      return 0;
    }

    initial_offset = input.tellg();
  }

  return 1;
}

/*

RDfile input looks like


RDFILE 1
$DATM 2016/09/16 21:58
$RFMT $RIREG 1
$RXN
EE3-E15410-047
      -NextMove-0822160938

  1  1  3
$MOL
EE3-E15410-047
     RDKit          2D

 19 19  0  0  0  0  0  0  0  0999 V2000
   12.6282   -6.1179    0.0000 C   0  0  0  0  0  0  0  0  1  0  0  0
   11.7646   -4.6121    0.0000 C   0  0  0  0  0  0  0  0  1  0  0  0

...
$RFMT $RIREG 2

*/

static int
get_reaction(iwstring_data_source& input, off_t& initial_offset, RXN_File& rxn)

{
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (buffer.starts_with("$RFMT $RIREG ")) {
      break;
    }
  }

  // cerr << "Starting reaction " << buffer << endl;

  if (input.eof()) {
    return 0;
  }

  initial_offset = input.tellg();

  return rxn.do_read(input);
}

static int
rxn_standardize_rdfile(iwstring_data_source& input,
                       IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  if (!input.next_record(buffer)) {
    cerr << "rxn_standardize_rdfile:cannot read header\n";
    return 0;
  }

  if ("$RDFILE" == buffer) {
  } else if (buffer.starts_with("$RDFILE")) {
  } else {
    cerr << "rxn_standardize_rdfile:warning header not 'RDFILE', got '" << buffer
         << "'\n";
  }

  while (true) {
    RXN_File rxn;

    initialise_reaction(rxn);

    off_t initial_offset = 0;

    if (get_reaction(input, initial_offset, rxn)) {
      ;
    } else if (input.eof()) {
      return 1;
    } else if (ignore_bad_reactions) {
      bad_reactions_ignored++;
      if (stream_for_discarded_reactions.active()) {
        echo_bad_data(input, initial_offset, stream_for_discarded_reactions);
        return 1;
      }
    } else {
      cerr << "Error reading reaction, line " << input.lines_read() << endl;
      if (stream_for_discarded_reactions.active()) {
        echo_bad_data(input, initial_offset, stream_for_discarded_reactions);
      }

      return 0;
    }

    reactions_read++;

    if (!rxn_standardize(rxn, input, initial_offset, output)) {
      return 0;
    }
  }

  return reactions_read;
}

static void
display_misc_X_options(std::ostream& output) {
  output << " -X rmmr         discard reactions that have multiple reactants\n";
  output << " -X rmmp         discard reactions that have multiple products\n";
  output << " -X rmiso        discard reactions containing isotopic atoms\n";
  output << " -X igbad        discard bad reactions (default is to exit)\n";
  output << " -X nclf         discard reactions where there is no change in the largest "
            "fragment\n";
  output << " -X rmncma       discard reactions that have no common mapped atoms in the "
            "reactants and product\n";
  output << " -X rmab         discard reactions that contain aromatic bonds\n";
  output << " -X rmdmap       discard reactions containing duplicate atom map numbers\n";

  exit(0);
}

static void
display_misc_D_options(std::ostream& output) {
  output << " -D rmdr         remove duplicate reactants\n";
  output << " -D rmdp         remove duplicate products\n";
  output << " -D rmda         remove duplicate agents\n";
  output << " -D rmaa         remove all agents\n";
  output << " -D rmnmap       remove any component with no atom mapping at all\n";
  output << " -D rmdup        remove components that are shown as both reactants and "
            "products\n";
  output << " -D rdlp         reduce product to largest product\n";
  output << " -D rdlr         reduce reactant to largest reactant\n";
  output << " -D rmrnp        remove reactants not participating\n";
  output << " -D rmfnp        remove any reactant fragments not participating\n";
  output << " -D unmap        unmap all atoms with duplicate atom map numbers\n";
  output << " -D fmap         add unique mapping numbers for any atoms that are not "
            "mapped\n";

  exit(0);
}

static int rxn_standardize(const char* fname, IWString_and_File_Descriptor& output);

static int
rxn_standardize_file_containing_reaction_files(iwstring_data_source& input,
                                               IWString_and_File_Descriptor& output) {
  IWString buffer;

  while (input.next_record(buffer)) {
    if (!rxn_standardize(buffer.null_terminated_chars(), output)) {
      return 0;
    }
  }

  return 1;
}

static int
rxn_standardize_file_containing_reaction_files(const char* fname,
                                               IWString_and_File_Descriptor& output) {
  IWString tmp(fname);
  assert(tmp.starts_with("F:"));
  tmp.remove_leading_chars(2);

  iwstring_data_source input(tmp.null_terminated_chars());

  if (!input.good()) {
    cerr << "Cannot open file of reactions '" << tmp << "'\n";
    return 0;
  }

  return rxn_standardize_file_containing_reaction_files(input, output);
}

static int
rxn_standardize(const char* fname, IWString_and_File_Descriptor& output) {
  assert(nullptr != fname);

  const_IWSubstring tmp(fname);
  if (tmp.starts_with("F:")) {
    return rxn_standardize_file_containing_reaction_files(fname, output);
  }

  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (reading_multi_reaction_rdfile) {
    return rxn_standardize_rdfile(input, output);
  }

  // if (input_is_reaction_smiles)
  return rxn_standardize_reaction_smiles(input, output);

  // return rxn_standardize(input, output);
}

static int
rxn_standardize(int argc, char** argv) {
  // Command_Line cl(argc, argv, "vA:E:g:lomU:Ibf:Kd:C:eRX:cD:s");
  Command_Line cl(argc, argv, "scKowd:f:D:F:C:X:U:Rg:E:A:Mqv");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage();
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage();
    }
  } else {
    set_global_aromaticity_type(Daylight);
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      usage();
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage();
    }
  }

  if (cl.option_present('o')) {
    auto_fix_orphans = 1;

    if (verbose) {
      cerr << "Will fix orphan conditions\n";
    }
  }

  if (cl.option_present('w')) {
    auto_fix_widows = 1;

    if (verbose) {
      cerr << "Will fix widow conditions\n";
    }
  }

  if (cl.option_present('f')) {
    const_IWSubstring f;
    for (int i = 0; cl.value('f', f, i); ++i) {
      if ("trunc" == f) {
        use_first_token_of_name = 1;
        if (verbose) {
          cerr << "Will truncate reaction names to the first token\n";
        }
      } else if ("gsub" == f) {
        gsub_reaction_names = 1;
        if (verbose) {
          cerr << "Will replace whitespace with underscores in reaction names\n";
        }
      } else if ("__" == f) {
        compress_consecutive_underscores = 1;
        if (verbose) {
          cerr << "Will compress consecutive underscores in reaction names\n";
        }
      } else {
        cerr << "Unrecognised -f qualifier '" << f << "'\n";
        usage();
      }
    }
  }

  if (cl.option_present('K')) {
    fix_kekule_problems = 1;

    if (verbose) {
      cerr << "Will fix Kekule problems\n";
    }
  }

  if (cl.option_present('d')) {
    const_IWSubstring d;
    for (int i = 0; cl.value('d', d, i); ++i) {
      if ('d' == d) {
        if (use_all_plusses) {
          cerr << "The options -d d and -d p cannot both be specified\n";
          usage();
        }

        use_all_dots = true;

        if (verbose) {
          cerr << "reactants and products separated by . rather than +\n";
        }
      } else if ('p' == d) {
        if (use_all_dots) {
          cerr << "The options -d d and -d p cannot both be specified\n";
          usage();
        }
        use_all_plusses = true;
        if (verbose) {
          cerr << "reactants and products separated by + rather than .\n";
        }
      } else {
        cerr << "Unrecognised -d qualifier '" << d << "'\n";
        usage();
      }
    }
  }

  if (cl.option_present('R')) {
    if (cl.option_present('s')) {
      cerr << "The -R (RD file) and -s (SMILES file) options cannot both be specified'\n";
      return 1;
    }
    reading_multi_reaction_rdfile = 1;
    // input_is_reaction_smiles = 0;
    if (verbose) {
      cerr << "Input from multi-reagent RDfile\n";
    }
  } else if (cl.option_present('s')) {
    // input_is_reaction_smiles = 1;   // this is the default - we really do not need to
    // set it
    if (verbose) {
      cerr << "Input is a reaction smiles\n";
    }
  }

  if (cl.option_present('c')) {
    // for MDL mols:
    set_mdl_molecule_discard_chirality(1);
    // for smiles:
    set_ignore_all_chiral_information_on_input(1);

    if (verbose) {
      cerr << "Will discard all chirality input on input\n";
    }
  }
  if (cl.option_present('M')) {
    do_automatic_atom_mapping = 1;
    if (verbose) {
      cerr << "Will do automatic atom mapping\n";
    }
  }

  if (cl.option_present('q')) {
    set_warn_no_mapped_atoms(0);
    if (verbose) {
      cerr << "Will not warn about no mapped atoms\n";
    }
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', max_atoms_in_reagent_and_product) ||
        max_atoms_in_reagent_and_product < 1) {
      cerr << "The maximum atoms in a reagent must be a whole +ve number\n";
      usage();
    }

    if (verbose) {
      cerr << "Will discard reactions with a reagent with more than "
           << max_atoms_in_reagent_and_product << " atoms\n";
    }
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if ("rmmr" == x) {
        skip_reactions_with_multiple_reactants = 1;
        if (verbose) {
          cerr << "Will discard reactions containing multiple reactants\n";
        }
      }

      else if ("rmncma" == x) {
        skip_reactions_with_no_common_mapped_atoms = 1;
        if (verbose) {
          cerr << "Will discard reactions that have not reactant atoms mapped to product "
                  "atoms\n";
        }
      }

      else if ("rmmp" == x) {
        skip_reactions_with_multiple_products = 1;
        if (verbose) {
          cerr << "Will discard reactions containing multiple products\n";
        }
      }

      else if ("rmiso" == x) {
        discard_reactions_with_isotopes = 1;
        if (verbose) {
          cerr << "Will discard reactions containing isotopic atoms\n";
        }
      } else if ("igbad" == x) {
        ignore_bad_reactions = 1;
        if (verbose) {
          cerr << "Will ignore bad reactions\n";
        }
      } else if ("nclf" == x) {
        skip_reactions_where_largest_fragment_is_unchanged = 1;
        if (verbose) {
          cerr << "Will skip reactions where the largest fragment is unchanged\n";
        }
      } else if ("rmab" == x) {
        skip_reactions_containing_aromatic_bonds = 1;
        if (verbose) {
          cerr << "Will skip reactions containing aromatic bonds\n";
        }
      } else if ("rmdmap" == x) {
        discard_reactions_containing_duplicate_atom_map_numbers = 1;

        if (verbose) {
          cerr << "Will discard reactions containing duplicate atom map numbers\n";
        }
      } else if ("unmap" == x) {
        unmap_duplicate_atom_map_numbers = 1;
        if (verbose) {
          cerr << "Will un-map atoms with duplicate atom map numbers\n";
        }
      }

      else if ("help" == x) {
        display_misc_X_options(cerr);
      } else {
        display_misc_X_options(cerr);
      }
    }
  }

  if (cl.option_present('D')) {
    const_IWSubstring d;
    for (int i = 0; cl.value('D', d, i); ++i) {
      if ("rmdr" == d) {
        remove_duplicate_reactants = 1;
        if (verbose) {
          cerr << "Will remove duplicate reactants\n";
        }
      } else if ("rmdp" == d) {
        remove_duplicate_products = 1;
        if (verbose) {
          cerr << "Will remove duplicate products\n";
        }
      } else if ("rmda" == d) {
        remove_duplicate_agents = 1;
        if (verbose) {
          cerr << "Will remove duplicate agents\n";
        }
      } else if ("rmaa" == d) {
        remove_all_agents = 1;
        if (verbose) {
          cerr << "Will remove all agents\n";
        }
      } else if ("rmnmap" == d) {
        remove_unmapped_components = 1;
        if (verbose) {
          cerr << "Will remove all components with no atom mapping\n";
        }
      } else if ("rmdup" == d) {
        remove_unchanged_components = 1;
        if (verbose) {
          cerr << "Will remove all components that are both reactants and product\n";
        }
      } else if ("rdlr" == d) {
        reduce_to_largest_reactant = 1;

        if (verbose) {
          cerr << "Will reduce reactants to largest fragment\n";
        }
      } else if ("rdlp" == d) {
        reduce_to_largest_product = 1;

        if (verbose) {
          cerr << "Will reduce products to largest fragment\n";
        }
      } else if ("rmrnp" == d) {
        eliminate_reagents_not_participating = 1;

        if (verbose) {
          cerr << "Will remove any reactants that are not participating\n";
        }
      } else if ("rmfnp" == d) {
        remove_fragments_not_participating = 1;

        if (verbose) {
          cerr << "Will remove any reactant fragments that are not participating\n";
        }
      } else if ("help" == d) {
        display_misc_D_options(cerr);
      } else if ("unmap" == d) {
        unmap_duplicate_atom_map_numbers = 1;
        if (verbose) {
          cerr << "Will un-map atoms with duplicate atom map numbers\n";
        }
      } else if ("fmap" == d) {
        assign_unmapped_atoms = 1;
        if (verbose) {
          cerr << "Will assign mapping numbers to any atoms that are not mapped\n";
        }
      } else {
        display_misc_D_options(cerr);
      }
    }
  }

  if (0 == cl.number_elements()) {
    cerr << "Insufficient arguments\n";
    usage();
  }

  if (cl.option_present('U')) {
    const char* u = cl.option_value('U');

    if (!stream_for_discarded_reactions.open(u)) {
      cerr << "Cannot open stream for discarded reactions '" << u << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Discarded reactions written to '" << u << "'\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++) {
    if (!rxn_standardize(cl[i], output)) {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << reactions_read << " reactions\n";
    cerr << bad_reactions_ignored << " bad reactions ignored\n";
    if (reactions_with_no_reagents) {
      cerr << reactions_with_no_reagents << " reactions with no reagents\n";
    }
    if (reactions_with_no_products) {
      cerr << reactions_with_no_products << " reactions with no products\n";
    }
    if (acc_natoms.n() > 0) {
      cerr << "input reactions had btw " << acc_natoms.minval() << " and "
           << acc_natoms.maxval() << " ave " << static_cast<float>(acc_natoms.average())
           << " atoms\n";
    }
    if (skip_reactions_with_multiple_reactants) {
      cerr << reactions_with_multiple_reactants_skipped
           << " reactions with multiple reagents skipped\n";
    }
    if (skip_reactions_with_multiple_products) {
      cerr << reactions_with_multiple_products_skipped
           << " reactions with multiple products skipped\n";
    }
    if (discard_reactions_with_isotopes) {
      cerr << reactions_discarded_for_isotopic_atoms
           << " reactions with isotopes skipped\n";
    }
    if (discard_reactions_containing_duplicate_atom_map_numbers) {
      cerr << reactions_containing_duplicate_atom_map_numbers
           << " reactions contained duplicate atom map numbers\n";
    }
    cerr << reactions_with_fragments_not_participating
         << " reactions with fragments not participating\n";
    cerr << reactions_with_reagents_not_participating
         << " reactions with reagents not participating\n";
    cerr << reactions_with_duplicate_reagents << " reactions with duplicate reagents\n";
    cerr << reactions_with_no_common_mapped_atoms
         << " reactions with no mapped atoms in common btw LHS and RHS\n";
    cerr << reactions_with_reagent_count_changed
         << " reactions where the reagent count was changed\n";
    if (skip_reactions_where_largest_fragment_is_unchanged) {
      cerr << reactions_where_largest_fragment_does_not_change
           << " reactions where the largest fragment did not change\n";
    }
    if (auto_fix_orphans) {
      cerr << reactions_with_orphan_atoms << " reactions with orphan atoms\n";
    }
    if (fix_kekule_problems) {
      cerr << reactions_with_kekule_problems_fixed
           << " reactions with Kekule problems fixed\n";
    }
    cerr << reactions_discarded_for_too_many_atoms
         << " reactions discarded for more than " << max_atoms_in_reagent_and_product
         << " atoms in a reagent and a product\n";
    cerr << reactions_containing_aromatic_bonds
         << " reactions discarded for aromatic bonds\n";
    cerr << reactions_written << " reactions written ("
         << (reactions_read - reactions_written) << " discarded)\n";
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = rxn_standardize(argc, argv);

  return rc;
}
