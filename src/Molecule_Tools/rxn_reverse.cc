// First step in retrosynthetic infrastrure.
// Primary task is to take forward reactions and reverse them.
// Performs various transformation and filtering operations.

#include <iostream>
#include <memory>

#include <stdlib.h>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rxn_file.h"

namespace rxn_reverse
{
using std::cerr;
using std::endl;

int verbose = 0;

const char * prog_name;

struct RxnReverseOptions {
  int input_is_reaction_smiles = 0;

  int reactions_that_cannot_be_read = 0;

  int ignore_bad_reactions = 0;

  int discard_reactions_with_isotopes = 0;
  int reactions_discarded_for_isotopic_atoms = 0;

  // Invalid reaction data written to..
  std::unique_ptr<IWString_and_File_Descriptor> stream_for_bad_reactions;

  int reactions_read = 0;

  // Successfully processed reactions written here.
  std::unique_ptr<IWString_and_File_Descriptor> output;

  // How many are written to 'output'.
  int reactions_written = 0;

  char output_separator = ' ';

  int reactions_discarded = 0;

  int reactions_containing_duplicate_atom_map_numbers = 0;

  int reactions_containing_multi_component_products = 0;

  // A reaction might have one component in the product array, but
  // that single component may have multiple fragments.
  int reactions_containing_multi_fragment_products = 0;
    
  // Once unchanging components and fragments are removed, we
  // might be left with nothing.
  int reactions_discarded_for_nothing_changing = 0;

  int component_grouping_is_plus = 0;

  Reaction_Smiles_Options reaction_smiles_options;
};

int
Report(const RxnReverseOptions& options, std::ostream& output)
{
  output << "Read " << options.reactions_read << " reactions, wrote " << options.reactions_written << "\n";
  output << options.reactions_that_cannot_be_read << " reactions could not be read\n";
  output << options.reactions_discarded_for_isotopic_atoms << " reactions_discarded_for_isotopic_atoms\n";
  output << options.reactions_containing_duplicate_atom_map_numbers << " reactions_containing_duplicate_atom_map_numbers\n";
  output << options.reactions_containing_multi_fragment_products << " reactions_containing_multi_fragment_products\n";
  output << options.reactions_containing_multi_component_products << " reactions_containing_multi_component_products\n";
  output << options.reactions_discarded_for_nothing_changing << " reactions_discarded_for_nothing_changing\n";

  return output.good();
}

// Write the reversed form of 'rxn' to 'output'.
int
WriteReversedReaction(RXN_File& rxn,
                      RxnReverseOptions& options,
                      IWString_and_File_Descriptor& output) {
  rxn.product(0).invalidate_smiles();
  output << rxn.product(0).smiles() << ">>";
  for (int i = 0; i < rxn.number_reagents(); ++i) {
    if (i > 0)
      output << '+';
    rxn.reagent(i).invalidate_smiles();
    output << rxn.reagent(i).smiles();
  }

  output <<  options.output_separator << rxn.name();

  output << "\n";

  options.reactions_written++;

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
rxn_reverse(RXN_File& rxn,
            RxnReverseOptions& options)
{
  rxn.remove_cis_trans_bonding();

  // Some always-do transformations.

  rxn.remove_unmapped_components();
  rxn.eliminate_reagents_not_participating();
  if (rxn.remove_fragments_not_participating ())
    cerr << "remove_fragments_not_participating\n";
  rxn.remove_non_participating_fragments ();
  rxn.remove_unchanging_fragments();
  rxn.remove_unchanging_components();

  if (rxn.number_reagents() == 0 || rxn.number_products() == 0)
  {
    if (verbose > 1)
      cerr << rxn.name() << " reaction eliminated due to unchanging items\n";
    options.reactions_discarded_for_nothing_changing++;
    return 0;
  }

  if (options.discard_reactions_with_isotopes && rxn.contains_isotopic_product_atoms())
  {
    if (verbose > 1)
      cerr << rxn.name() << " contains isotopic atoms\n";
    options.reactions_discarded_for_isotopic_atoms++;
    return 0;
  }

  if (rxn.contains_duplicate_atom_map_numbers())
  {
    if (verbose > 1)
      cerr << rxn.name() << " contains duplicate atom map numbers\n";
    options.reactions_containing_duplicate_atom_map_numbers++;
    return 0;
  }

  if (! rxn.at_least_some_mapped_atoms_common_btw_reagents_and_products()) {
    if (verbose > 1)
      cerr << rxn.name() << " no matched atoms across reaction\n";
    return 0;
  }

  // Reaction is still a forward reaction. If there are multiple products with atom maps
  // it cannot be used as a reverse reaction.

  if (rxn.number_products() > 1) {
    if (verbose > 1)
      cerr << rxn.name() << " multi component product, cannot be reversed\n";
    options.reactions_containing_multi_component_products++;
    return 0;
  }

  if (rxn.product(0).number_fragments() > 1) {
    if (verbose > 1)
      cerr << rxn.name() << " multi fragment product, cannot be reversed\n";
    options.reactions_containing_multi_fragment_products++;
    return 0;
  }

  return WriteReversedReaction(rxn, options, *options.output);
}

int
echo_bad_data(iwstring_data_source & input,
              const off_t initial_offset,
              IWString_and_File_Descriptor & output)
{
  const auto current_offset = input.tellg();
  input.seekg(initial_offset);
//cerr << "Bad reaction begin " << initial_offset << " now " << current_offset << endl;
  input.echo(output, (current_offset - initial_offset));
  input.seekg(current_offset, 0);
  assert (input.tellg() == current_offset);
  output.write_if_buffer_holds_more_than(4096);

  return input.tellg() == current_offset;
}

// Functions that get the next reaction need to
// send various states back to their caller.
enum next_reaction_result {
  // All is good, reaction can be processed.
  valid_reaction,

  // normal eof
  eof,
  
  // No reaction, but ok to fetch the next one.
  no_reaction_ok_to_continue,

  // No reaction, fatal error.
  no_reaction_fatal_error
};

next_reaction_result
next_reaction_smiles(iwstring_data_source & input,
                     RxnReverseOptions& options,
                     RXN_File & rxn)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
    return eof;

  cerr << "SMILES INPUT '" << buffer << "'\n";

  if (rxn.build_from_reaction_smiles(buffer, options.component_grouping_is_plus))
    return valid_reaction;

  cerr << "Cannot interpret " << buffer << endl;
  options.reactions_that_cannot_be_read++;
  if (options.ignore_bad_reactions)
    return no_reaction_ok_to_continue;

  return no_reaction_fatal_error;
}

// Get the next reaction as either reaction smiles or rxn file.

next_reaction_result
get_next_reaction(iwstring_data_source & input,
                  RxnReverseOptions& options,
                  RXN_File & rxn)
{
  if (options.input_is_reaction_smiles)
    return next_reaction_smiles(input, options, rxn);

  const auto initial_offset = input.tellg();

  if (rxn.do_read(input))
    return valid_reaction;
  if (input.eof())
    return eof;
  if (options.ignore_bad_reactions)
  {
    options.reactions_that_cannot_be_read++;

    if (options.stream_for_bad_reactions)
      echo_bad_data(input, initial_offset, *options.stream_for_bad_reactions);
    return no_reaction_ok_to_continue;
  }

  cerr << "Fatal error reading reaction, now at line " << input.lines_read() << endl;
  return no_reaction_fatal_error;
}

int
rxn_reverse(iwstring_data_source& input,
            RxnReverseOptions& options)
{
  input.set_translate_tabs(1);

  while (true)
  {
    RXN_File rxn;
    rxn.set_do_automatic_atom_mapping(0);

    const auto status = get_next_reaction(input, options, rxn);
    switch(status) {
      case valid_reaction:
        break;
      case eof:
        return 1;
      case no_reaction_ok_to_continue:
        continue;
      case no_reaction_fatal_error:
        return 0;
    }

    (void) rxn_reverse(rxn, options);  // No checking return code.
  }

  return 1;
}

int rxn_reverse(const char * fname,
                RxnReverseOptions& options)
{
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "rxn_reverse:cannot open '" << fname << "'\n";
    return 0;
  }

  return rxn_reverse(input, options);
}

void
usage(int rc)
{
  exit(rc);
}

void
display_discard_options(std::ostream & output)
{
  exit(1);
}

int
rxn_reverse(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:D:p");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }
  else 
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }
  else
    set_auto_create_new_elements(1);

  RxnReverseOptions options;

  if (cl.option_present('i')) {
    const_IWSubstring i = cl.string_value('i');
    if (i == "rsmi") {
      options.input_is_reaction_smiles = 1;
    } else if (i == "rxn") {
      options.input_is_reaction_smiles = 0;
    } else {
      cerr << "Unrecognised -i qualifier '" << i << "'\n";
      return 1;
    }
  }

  if (cl.option_present('D'))
  {
    const_IWSubstring d;
    for (int i = 0; cl.value('D', d, i); ++i) {
      if ("" == d)
      {
      }
      else if ("help" == d)
      {
        display_discard_options(cerr);
      }
      else
      {
        cerr << "Unrecognised -D qualifier '" << d << "'\n";
        display_discard_options(cerr);
      }
    }
  }

  if (cl.option_present('p')) {
    options.component_grouping_is_plus = 1;
    if (verbose)
      cerr << "Input assumed to be + separated components\n";
  }

  if (cl.option_present('S')) {
    const char * s = cl.option_value('S');
    options.output = std::make_unique<IWString_and_File_Descriptor>();
    if (! options.output->open(s)) {
      cerr << "Cannot open output file '" << s << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Output to '" << s << "'\n";
  }
  else {  // Use stdout.
    options.output = std::make_unique<IWString_and_File_Descriptor>(1);
  }

  for (auto * file_name : cl) {
    if (! rxn_reverse(file_name, options)) {
      cerr << "Fatal error processing '" << file_name << "'\n";
      return 1;
    }
  }

  if (verbose) {
    Report(options, cerr);
  }

  return 0;
}
}

int
main(int argc, char ** argv)
{
  rxn_reverse::prog_name = argv[0];

  int rc = rxn_reverse::rxn_reverse(argc, argv);

  return rc;
}
