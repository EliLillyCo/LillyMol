/*
  Merges molecules into one
*/

#include <stdlib.h>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "  -S <fname>    output file name\n";
  cerr << "  -o <type>     output type(s)\n";
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

  return;
}

static int
molecular_merge (data_source_and_type<Molecule> & input,
                 Molecule & all_of_them)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    all_of_them.add_molecule(m);
  }

  return 1;
}

static int
molecular_merge (const char * fname, FileType input_type, 
                 Molecule & all_of_them)
{
  assert (nullptr != fname);

  if (input_type == FILE_TYPE_INVALID)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return molecular_merge(input, all_of_them);
}

static int
molecular_merge (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:lS:o:");

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

  Molecule_Output_Object output;

  if (! cl.option_present('o'))
    output.add_output_type(FILE_TYPE_SMI);
  else if (! output.determine_output_types(cl, 'o'))
  {
    cerr << "Cannot determine output type(s) (-o)\n";
    usage(1);
  }

  if (! cl.option_present('S'))
  {
    cerr << "Must specify output file name via the -S option\n";
    usage(1);
  }
  else
  {
    const_IWSubstring s = cl.string_value('S');

    if (output.would_overwrite_input_files(cl, s))
    {
      cerr << "Cannot overwrite input file(s) '" << s << "'\n";
      return 2;
    }

    if (! output.new_stem(s))
    {
      cerr << "Could not initialise output stem '" << s << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Output file(s) created with stem '" << s << "'\n";
  }

  Molecule all_of_them;

  all_of_them.resize(10000);

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! molecular_merge(cl[i], input_type, all_of_them))
      return i + 1;
  }

  if (verbose)
    cerr << "Read " << molecules_read << " molecules, final molecule contains " << all_of_them.natoms() << " atoms\n";

  if (! output.write(all_of_them))
  {
    cerr << "Yipes, could not write resulting molecule\n";
    return 1;
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = molecular_merge (argc, argv);

  return rc;
}
