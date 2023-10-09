/*
  Writes either the unique smiles (if it can be repinterpreted)
  or a non-aromatic "unique" smiles
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/molecule.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static int molecules_failing_reinterpretation = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int remove_chirality = 0;

static IWString_and_File_Descriptor stream_for_failed_reinterpretation;

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
  cerr << "  -U <fname>    write failed interpretation molecules to <fname>\n";
  cerr << "  -c            remove chirality\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (remove_chirality)
    m.remove_all_chiral_centres();

  return;
}

/*
  Here's where you do whatever you want to do with the molecule
  In this case, we count the number of nitrogen atoms
*/

static int
preferred_smiles (Molecule & m,
                  IWString_and_File_Descriptor & output)
{
  IWString s = m.unique_smiles();

  Molecule tmp;

  if (! tmp.build_from_smiles(s))
  {
    output << m.non_aromatic_unique_smiles() << ' ' << m.name() << '\n';
    molecules_failing_reinterpretation++;
    if (verbose > 1)
      cerr << "Failed aromatic form '" << m.name() << "'\n";

    if (stream_for_failed_reinterpretation.is_open())
    {
      stream_for_failed_reinterpretation << m.non_aromatic_unique_smiles() << ' ' << m.name() << '\n';
      stream_for_failed_reinterpretation.write_if_buffer_holds_more_than(32768);
    }
  }
  else
    output << s << ' ' << m.name() << '\n';

  return 1;
}

static int
preferred_smiles (data_source_and_type<Molecule> & input,
                  IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! preferred_smiles(*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
preferred_smiles (const char * fname, FileType input_type, 
                  IWString_and_File_Descriptor & output)
{
  assert (nullptr != fname);

  if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return preferred_smiles(input, output);
}
static int
preferred_smiles (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lU:c");

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

  set_input_aromatic_structures(1);
  set_display_no_kekule_form_message(0);

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

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('c'))
  {
    remove_chirality = 1;

    if (verbose)
      cerr << "Will remove chirality\n";
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = FILE_TYPE_SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('U'))
  {
    IWString u = cl.string_value('U');

    if (! u.ends_with(".smi"))
      u << ".smi";

    if (! stream_for_failed_reinterpretation.open (u.null_terminated_chars()))
    {
      cerr << "Cannot open stream for failed interpretation '" << u << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Molecules for which unique smiles cannot be reinterpreted written to '" << u << "'\n";
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! preferred_smiles(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
    cerr << "Read " << molecules_read << " molecules, " << molecules_failing_reinterpretation << " failed reinterpretation\n";

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = preferred_smiles (argc, argv);

  return rc;
}
