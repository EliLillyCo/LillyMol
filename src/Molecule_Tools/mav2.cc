/*
  Molecular abstractions
*/

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/iwmfingerprint.h"
#define ISTREAM_AND_TYPE_IMPLEMENTATION
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

#include "molecular_abstraction_specifications.h"
#include "molecular_abstraction_functions.h"

using std::cerr;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int write_parent_molecule = 0;

static int reduce_to_largest_fragment = 0;

static int remove_all_chiral_centres = 0;

static int remove_all_cis_trans_bonds = 0;

static int unfix_implicit_hydrogens = 0;

static int number_abstraction_sets = 0;

static Set_of_Molecular_Abstractions * mabs = nullptr;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static int writing_fingerprints = 0;

static int work_as_tdt_filter = 0;

/*
  March 2014. We need to do outputs other than smiles.
  With this kind of output, there is only one kind of
  output, the molecule is written at the end of
  all processing.
*/

static Molecule_Output_Object final_molecule_stream;

static int min_atoms_in_final_molecule = 0;

static int flush_after_each_molecule = 0;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "  -a ...        specify abstraction(s) to be created, enter '-a help for info\n";
  cerr << "  -B <fname>    specify abstraction(s) in a file, same syntax\n";
  cerr << "  -p            write the parent molecule\n";
  cerr << "  -z ...        options for what to do when no changes during a stage\n";
  cerr << "  -C            remove invalid chirality from output molecules\n";
  cerr << "  -h            unfix any explicit implicit hydrogen specifications\n";
  cerr << "  -c            remove all chirality from input molecules\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -t            remove cis-trans bonds from input\n";
  cerr << "  -Y            standard fingerprint options, enter '-Y help' for info\n";
  cerr << "  -f            work as a TDT filter\n";
  cerr << "  -F <stem>     write final molecule to file <stem>\n";
  cerr << "  -o <type>     type for -F file\n";
  cerr << "  -X ...        miscellaneous options, enter '-X help' for info\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -K ...        standard smiles options\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
preprocess(Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (remove_all_chiral_centres)
    m.remove_all_chiral_centres();

  if (remove_all_cis_trans_bonds)
    m.revert_all_directional_bonds_to_non_directional();

  if (unfix_implicit_hydrogens)
  {
    int matoms = m.natoms();
    for (int i = 0; i < matoms; i++)
    {
      m.set_implicit_hydrogens_known(i, 0);
    }
  }

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

static int
mav2(Molecule_With_Info_About_Parent & m,
     IWString_and_File_Descriptor & output)
{
  if (writing_fingerprints)
    ;
  else if (write_parent_molecule)
    output << m.smiles() << ' ' << m.name() << '\n';

  for (int i = 0; i < number_abstraction_sets; i++)
  {
    if (! mabs[i].process(m, output))
    {
      cerr << "Fatal error '" << m.name() << "' set " << i << '\n';
      return 0;
    }
  }

  if (final_molecule_stream.active())
  {
    if (m.natoms() >= min_atoms_in_final_molecule)
      return final_molecule_stream.write(m);
  }

  return 1;
}

static int
write_smiles_and_identifier(Molecule & m,
                            IWString_and_File_Descriptor & output)
{
  output << smiles_tag << m.smiles() << ">\n";
  output << identifier_tag << m.name() << ">\n";

  return 1;
}

static void
MaybeFlush(IWString_and_File_Descriptor& output) {
  if (flush_after_each_molecule) {
    output.flush();
  } else {
    output.write_if_buffer_holds_more_than(8192);
  }
}

static int
mav2 (data_source_and_type<Molecule_With_Info_About_Parent> & input,
                IWString_and_File_Descriptor & output)
{
  Molecule_With_Info_About_Parent * m;
  while (nullptr != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule_With_Info_About_Parent> free_m(m);

    preprocess(*m);

    if (writing_fingerprints)
      write_smiles_and_identifier(*m, output);

    if (! mav2(*m, output))
      return 0;

    if (writing_fingerprints) {
      output << "|\n";
    }

    MaybeFlush(output);
  }

  return 1;
}

static int
mav2_filter(const const_IWSubstring & buffer,
            IWString_and_File_Descriptor & output)
{
  assert (buffer.starts_with(smiles_tag));
  assert (buffer.ends_with('>'));

  const_IWSubstring smiles(buffer);

  smiles.remove_leading_chars(smiles_tag.length());
  smiles.chop();

  Molecule_With_Info_About_Parent m;

  if (! m.build_from_smiles(smiles))
  {
    cerr << "Invalid smiles '" << smiles << "'\n";
    return 0;
  }

  return mav2(m, output);
}

static int
mav2_filter(iwstring_data_source & input,
            IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(32768);

    if (! buffer.starts_with(smiles_tag))
      continue;

    if (! mav2_filter(buffer, output))
      return 0;
  }

  return 1;
}

static int
mav2 (const char * fname, FileType input_type, 
                IWString_and_File_Descriptor & output)
{
  assert (nullptr != fname);

  if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule_With_Info_About_Parent> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return mav2(input, output);
}

/*
  Bad idea, you should be able to have abs1.abs2.abs3... in an input file
  Change this sometime to just get rid of trailing #.... things
*/

#ifdef STRIP_TRAILING_COMMENTS_IF_PRESENT_BROKEN
static void
strip_trailing_comments_if_present (const_IWSubstring & buffer)
{
  int n = buffer.length();

  int paren_level = 0;
  for (int i = 0; i < n; i++)
  {
    char c = buffer[i];

    if ('(' == c)
      paren_level++;
    else if (')' == c)
    {
      paren_level--;

      if (paren_level > 0)
        continue;

      if (i == n - 1)
        return;

      buffer.iwtruncate(i + 1);
      return;
    }
  }

  return;    // should not come to here
}
#endif

static int
read_abstraction_directives_from_file (iwstring_data_source & input,
                                       Set_of_Molecular_Abstractions & mabs)
{
  const_IWSubstring buffer;

  IWString abstractions;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;

//  strip_trailing_comments_if_present(buffer);   implementation is wrong, fix sometime

    if (0 == buffer.length())
      continue;

    abstractions.append_with_spacer(buffer, '.');
  }

  if (0 == abstractions.length())
  {
    cerr << "No directives in file\n";
    return 0;
  }

  Molecular_Abstraction_Directives_Node madn;

  cerr << "Building directive from '" << abstractions << "'\n";
  if (! madn.build(abstractions))
  {
    cerr << "INvalid molecular abstraction directive '" << abstractions << "'\n";
    return 0;
  }

  if (! madn.directive_recognised())
  {
    cerr << "Unrecognised directives '" << abstractions << "'\n";
    return 0;
  }

  if (! mabs.build(madn))
  {
    cerr << "Cannot initialise abstraction from '" << abstractions << "'\n";
    return 0;
  }

  return 1;
}

static int
read_abstraction_directives_from_file (const char * fname,
                                       Set_of_Molecular_Abstractions & mabs)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_abstraction_directives_from_file(input, mabs);
}

static void
display_dash_z_options (std::ostream & os)
{
  os << " -z writeif     only follow WRITE directive if changes at that stage\n";
  os << " -z writec      append number of changes to the tag\n";
  os << " -z 0empty      write an empty molecule if the selection criterion does not match\n";

  exit(1);
}

static void
DisplayDashXOptions(std::ostream& output) {
  output << " -X flush        flush output after each molecule\n";
  ::exit(0);
}

static int
mav2 (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:la:B:pY:fcCthK:z:o:F:m:X:");

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

  set_aromatic_rings_must_contain_unsaturation(0);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('K'))
  {
    if (! process_standard_smiles_options(cl, verbose, 'K'))
    {
      cerr << "Cannot initialise standard smiles options (-K)\n";
      return 5;
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('Y'))
  {
    if (! parse_misc_fingerprint_options(cl, 'Y', verbose))
    {
      cerr << "Cannot initialise fingerprint specifications (-Y)\n";
      return 4;
    }
  }

  if (cl.option_present('c'))
  {
    remove_all_chiral_centres = 1;
    if (verbose)
      cerr << "Chirality will be removed\n";
  }
  else if (cl.option_present('C'))
  {
    set_remove_invalid_chiral_centres_before_writing(1);

    if (verbose)
      cerr << "Invalid chirality removed from all output\n";
  }

  if (cl.option_present('t'))
  {
    remove_all_cis_trans_bonds = 1;
    if (verbose)
      cerr << "Cis trans bonds will be removed\n";
  }

  if (cl.option_present('h'))
  {
    unfix_implicit_hydrogens = 1;
    if (verbose)
      cerr << "Any explicit specification of implicit hydrogens will be discarded\n";
  }

  if (cl.option_present('i') && cl.option_present('f'))
  {
    cerr << "Sorry, the -i and -f options are mutually exclusive\n";
    usage(4);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (cl.option_present('f'))
  {
    work_as_tdt_filter = 1;
  }
  else if (1 == cl.number_elements() && 0 == strncmp("-", cl[0], 1))
    input_type = FILE_TYPE_SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.option_present('p'))
  {
    write_parent_molecule = 1;

    if (verbose)
      cerr << "Will write the parent molecule\n";
  }

  if (cl.option_present('z'))
  {
    int i = 0;
    const_IWSubstring z;
    while(cl.value('z', z, i++))
    {
      if ("writeif" == z)
      {
        set_write_only_if_changes(1);
      }
      else if ("writec" == z)
      {
        set_append_count_to_tag(1);
        if (verbose)
          cerr << "Will append count of changes to the tag\n";
      }
      else if ("0empty" == z)
      {
        set_write_empty_molecule_on_no_match(1);
        if (verbose)
          cerr << "Will write an empty molecule on zero query match\n";
      }
      else if ("help" == z)
        display_dash_z_options(cerr);
      else
      {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        display_dash_z_options(cerr);
      }
    }
  }

  number_abstraction_sets = cl.option_count('a') + cl.option_count('B');

  if (0 == number_abstraction_sets)
  {
    cerr << "Must specify one or more abstractions via the -a option\n";
    usage(4);
  }

// this really isn't tested for more than 1 abstraction set

  mabs = new Set_of_Molecular_Abstractions[number_abstraction_sets]; std::unique_ptr<Set_of_Molecular_Abstractions[]> free_mabs(mabs);

  int ndx = 0;

  for (int i = 0; i < cl.option_count('a'); i++)
  {
    const_IWSubstring a = cl.string_value('a', i);

    if (a == "help") {
      DisplayUsageExamples(cerr);
      return 0;
    }

    Molecular_Abstraction_Directives_Node madn;

//  cerr << "Building directive from '" << a << "'\n";
    if (! madn.build(a))
    {
      cerr << "INvalid molecular abstraction directive '" << a << "'\n";
      return i + 1;
    }

    if (! madn.directive_recognised())
    {
      cerr << "Unrecognised directives '" << a << "'\n";
      return 0;
    }

    if (! mabs[ndx].build(madn))
    {
      cerr << "Cannot initialise abstraction from '" << a << "'\n";
      return 4;
    }

    ndx++;
  }

  for (int i = 0; i < cl.option_count('B'); i++)
  {
    const char * fname = cl.option_value('B');
    if (! read_abstraction_directives_from_file(fname, mabs[ndx]))
    {
      cerr << "Invalid abstraction directive(s) in '" << fname << "'\n";
      return 5;
    }
    ndx++;
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "flush") {
        flush_after_each_molecule = 1;
        if (verbose) {
          cerr << "Will flush output after each molecule\n";
        }
      } else if (x == "smiles") {
        for (int i = 0; i < number_abstraction_sets; ++i) {
          mabs[i].set_write_unique_smiles(0);
        }
      } else if (x == "help") {
        DisplayDashXOptions(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(cerr);
      }
    }
  }

  int write_smiles = 0;
  for (int i = 0; i < number_abstraction_sets; i++)
  {
    if (! mabs[i].what_is_being_written(write_smiles, writing_fingerprints))
    {
      cerr << "Abstraction set " << i << " is invalid\n";
      return i + 1;
    }
  }

  if (cl.empty())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('F'))
  {
    if (cl.option_present('o'))
    {
      if (! final_molecule_stream.determine_output_types(cl, 'o'))
      {
        cerr << "Cannot determine output type(s) for -F file\n";
        return 2;
      }
    }
    else
      final_molecule_stream.add_output_type(FILE_TYPE_SMI);

    const_IWSubstring f = cl.string_value('F');

    if (final_molecule_stream.would_overwrite_input_files(cl, f))
    {
      cerr << "Cannot overwrite input file(s) '" << f << "'\n";
      return 2;
    }

    if (! final_molecule_stream.new_stem(f))
    {
      cerr << "Cannot initialise final molecule stem '" << f << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Final processed file written to '" << f << "'\n";

    if (cl.option_present('m'))
    {
      if (! cl.value('m', min_atoms_in_final_molecule) || min_atoms_in_final_molecule < 1)
      {
        cerr << "The min atoms in final molecule (-m) option must be a whole +ve number\n";
        usage(2);
      }

      if (verbose)
        cerr << "Will only write final molecules if they have at least " << min_atoms_in_final_molecule << " atoms\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc;

  if (work_as_tdt_filter)
  {
    iwstring_data_source input("-");
    rc = mav2_filter(input, output);
  }
  else
  {
    rc = 0;
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! mav2(cl[i], input_type, output))
      {
        rc = i + 1;
        break;
      }
    }
  }

  if (output.size())
    output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";

    for (int i = 0; i < number_abstraction_sets; i++)
    {
      mabs[i].report(cerr);
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = mav2(argc, argv);

  return rc;
}
