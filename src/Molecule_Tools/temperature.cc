#include <memory>
#include <limits>
#include <iostream>

using std::cerr;
using std::endl;

/*
  Creates a temperature fingerprint based on the atom count
  Works from a tdt.
*/

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_with_too_many_atoms = 0;

#include "cmdline.h"
#include "accumulator.h"
#include "iwbits.h"
#include "sparse_fp_creator.h"
#include "misc.h"

#include "molecule.h"
#include "aromatic.h"
#include "istream_and_type.h"
#include "path.h"
#include "iwstandard.h"

#include "mpr.h"

static IWString temperature_tag;
static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static Accumulator_Int<int> atom_statistics;

static Accumulator_Int<int> * accumulators = NULL;

static int allow_multi_fragment_molecules = 0;

static int write_array_of_properties = 0;

static int unsaturation_includes_aromatic = 0;

static int function_as_filter = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static unsigned char abbreviated_fingerprints = 0;

#define ABBREV_NATOMS 1
#define ABBREV_LARGEST_RING_SIZE 2
#define ABBREV_NRINGS 4
#define ABBREV_RING_ATOMS 8
#define ABBREV_AROAMTIC_ATOMS 16
#define ABBREV_FUSED_RING_ATOMS 32
#define ABBREV_HETEROATOM_COUNT 64
#define ABBREV_UNSATURATION_COUNT 128

/*
  We can optionally divide the atom count properties into coarser
  divisions.
*/

static int atom_count_divide = 1;

static int bit_replicates = 8;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __TIME__ << " " << __DATE__ << endl;
  cerr << "Usage" << prog_name << " <options> <input_file>\n";

  cerr << "  -J <dataitem>  specify dataitem tag for properties bits (default " << temperature_tag << ")\n";
  cerr << "  -B ...         abbreviated fingerprint specification, enter '-B help'\n";
  cerr << "  -d <num>       divisor for atom count related properties\n";
  cerr << "  -r <num>       bit replicates for abbreviated properties\n";
  cerr << "  -a             output is a descriptor file\n";
  cerr << "  -m             allow multi fragment molecules\n";
  cerr << "  -u             count aromatic atoms as unsaturated\n";
  cerr << "  -f             work as a filter\n";
  cerr << "  -i <type>      specify input file type (except with -f)\n";
  (void) display_standard_aromaticity_options (cerr);
  (void) display_standard_chemical_standardisation_options (cerr, 'g');
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -E <X>         create element 'X', 'autocreate' for auto creation\n";
  cerr << "  -v             verbose output\n";

  exit (rc);
}

static int
write_descriptor_file_header(IWString_and_File_Descriptor & output)
{
  output << "ID natoms lgrsize nrings ring_atoms aroma fused_atoms htroa unsat\n";

  return 1;
}

static int
do_write_fingerprint (Molecule & m,
                      const unsigned char * properties,
                      IWString_and_File_Descriptor & output)
{
  IW_Bits_Base fp;
  fp.construct_from_array_of_bits (properties, NPROPERTIES * IW_BITS_PER_BYTE);

  IWString tmp;
  fp.daylight_ascii_representation_including_nset_info (tmp);
  output << temperature_tag << tmp << ">\n";

  return 1;
}

static void
set_bit_in_each_creator (Sparse_Fingerprint_Creator & sfc,
                         int n,
                         int b,
                         int c)
{
  for (int i = 0; i < n; i++)
  {
    sfc.hit_bit(i * n + b, c);
  }

  return;
}

static int
do_write_abbreviated_properties (const Molecule & m,
                                 const unsigned char * properties,
                                 IWString_and_File_Descriptor & output)
{
  Sparse_Fingerprint_Creator sfc;

//cerr << "To write " << static_cast<int>(abbreviated_fingerprints) << endl;

  if (ABBREV_NATOMS & abbreviated_fingerprints && properties[0] >= atom_count_divide)
  {
    int c = static_cast<int>(properties[0]) / atom_count_divide;

    set_bit_in_each_creator(sfc, bit_replicates, 0, c);
  }

  if (ABBREV_LARGEST_RING_SIZE & abbreviated_fingerprints && properties[1] > 0)
  {
    int c = static_cast<int>(properties[1]);
    set_bit_in_each_creator(sfc, bit_replicates, 1, c);
  }

  if (ABBREV_NRINGS & abbreviated_fingerprints && properties[2] > 0)
  {
    int c = static_cast<int>(properties[2]);
    set_bit_in_each_creator(sfc, bit_replicates, 2, c);
  }

  if (ABBREV_RING_ATOMS & abbreviated_fingerprints && properties[3] > 0)
  {
    int c = static_cast<int>(properties[3]);
    set_bit_in_each_creator(sfc, bit_replicates, 3, c);
  }

  if (ABBREV_AROAMTIC_ATOMS & abbreviated_fingerprints && properties[4] >= atom_count_divide)
  {
    int c = static_cast<int>(properties[4]) / atom_count_divide;
    set_bit_in_each_creator(sfc, bit_replicates, 4, c);
  }

  if (ABBREV_FUSED_RING_ATOMS & abbreviated_fingerprints && properties[5] > 0)
  {
    int c = static_cast<int>(properties[5]);
    set_bit_in_each_creator(sfc, bit_replicates, 5, c);
  }

  if (ABBREV_HETEROATOM_COUNT & abbreviated_fingerprints && properties[6] >= atom_count_divide)
  {
    int c = static_cast<int>(properties[6] / atom_count_divide);
    set_bit_in_each_creator(sfc, bit_replicates, 6, c);
  }

  if (ABBREV_UNSATURATION_COUNT & abbreviated_fingerprints && properties[7] > 0)
  {
    int c = static_cast<int>(properties[7]);    // don't bother dividing
    set_bit_in_each_creator(sfc, bit_replicates, 7, c);
  }

  IWString tmp;

  sfc.daylight_ascii_form_with_counts_encoded(temperature_tag, tmp);

  output << tmp << "\n";

  return 1;
}

static int
do_write_array_of_properties (const Molecule & m,
                              const unsigned char * properties,
                              IWString_and_File_Descriptor & output)
{
  const IWString & mname = m.name();

  int ispace = mname.index(' ');

  if (ispace < 0)
    output << mname;
  else
  {
    const_IWSubstring tmp(mname);
    tmp.iwtruncate(ispace);
    output << tmp;
  }

  for (int i = 0; i < NPROPERTIES; i++)
  {
    output << ' ' << static_cast<int>(properties[i]);
  }

  output << '\n';

  return 1;
}


static int
preprocess(Molecule & m)
{
  if (allow_multi_fragment_molecules)
    ;
  else if (reduce_to_largest_fragment && m.number_fragments() > 1)
    m.reduce_to_largest_organic_fragment();
  else if (m.number_fragments() > 1)
  {
    cerr << "Error, molecule has " << m.number_fragments() << " components\n";
    return 0;
  }

  if (chemical_standardisation.active())
    chemical_standardisation.process (m);

  return 1;
}

/*
  May 2000. Found differences on different platforms because of differing
  SSSR determinations. Therefore change things to use all the rings - doesn't
  matter, all will be treated the same
*/

static Molecular_Properties_Generator mpg;

static int
temperature (Molecule & m,
             IWString_and_File_Descriptor & output)
{
  unsigned char properties[NPROPERTIES];

  int matoms = m.natoms();

  atom_statistics.extra (matoms);

  if (matoms > std::numeric_limits<unsigned char>::max())
  {
    molecules_with_too_many_atoms++;
    if (verbose)
      cerr << "Molecule with too many atoms " << matoms << endl;
  }

  if (! mpg(m, properties))
  {
    cerr << "molecular_properties_generation failed!\n";
    set_vector(properties, NPROPERTIES, static_cast<unsigned char>(0));
  }

  if (NULL != accumulators)
  {
    for (int i = 0; i < NPROPERTIES; i++)
    {
      accumulators[i].extra (int (properties[i]));
    }
  }

  if (write_array_of_properties)
    return do_write_array_of_properties(m, properties, output);
  else if (0 != abbreviated_fingerprints)
    return do_write_abbreviated_properties(m, properties, output);
  else if (temperature_tag.length() > 0)
    return do_write_fingerprint(m, properties, output);

  cerr << "What kind of output am I supposed to be doing???\n";

  return 0;
}

static int
temperature (iwstring_data_source & input,
             IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (! buffer.starts_with (smiles_tag))
    {
      output << buffer << '\n';

      output.write_if_buffer_holds_more_than(32768);

      continue;
    }

    const_IWSubstring smi = buffer.substr (5);
    smi.chop();

    Molecule m;
    if (! m.build_from_smiles (smi))
    {
      cerr << "Yipes, cannot parse smiles, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }

    preprocess(m);

    output << smiles_tag << m.smiles() << ">\n";

    if (! temperature (m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
temperature (data_source_and_type<Molecule> & input,
             IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    if (! preprocess(*m))
      return 0;

    if (! write_array_of_properties)
    {
      output << smiles_tag << m->smiles() << ">\n";
      output << identifier_tag << m->name() << ">\n";
    }

    int rc = temperature(*m, output);

    if (! write_array_of_properties)
      output << "|\n";

    output.write_if_buffer_holds_more_than(32768);

    if (0 == rc)
      return 0;
  }

  //GH: Fix bug# 19048
  if (input.stopped_because_of_error())
  {
	  // Return fail due to error
	  return 0;
  }
  else
  {
	  // Return success
	  return 1;
  }
  //Original: return 1;
}

static int
temperature (const char * fname,
             int input_type,
             IWString_and_File_Descriptor & output)
{
  if (function_as_filter)
  {
    iwstring_data_source input(fname);
    if (! input.ok())
    {
      cerr << "Cannot open '" << fname << "' for input\n";
      return 0;
    }

    return temperature (input, output);
  }

  data_source_and_type<Molecule> input (input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
	//GH: Fix bug# 19048
    return 0; //Original: return 1;
  }

  if (verbose > 1)
    input.set_verbose (1);

  return temperature (input, output);
}

static void
display_dash_B_qualifiers (std::ostream & os)
{
  os << " -B n     natoms\n";
  os << " -B l     largest ring size\n";
  os << " -B r     nrings\n";
  os << " -B g     ring atoms\n";
  os << " -B a     aromatic atoms\n";
  os << " -B f     fused ring atoms\n";
  os << " -B h     heteroatom count\n";
  os << " -B u     unsaturated atoms\n";

  exit(1);
}

static int
temperature (int argc, char ** argv)
{
  //GH: Fix bug# 19048
  //Added for checking the error code
  int error_code = 1;
  Command_Line cl (argc, argv, "J:vA:mauE:sfi:lg:B:d:r:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (! process_elements (cl, verbose))
  {
    cerr << "Cannot parse -E option\n";
    usage (7);
  }

  if (! process_standard_aromaticity_options (cl))
  {
    cerr << "Cannot parse aromaticity options\n";
    usage (8);
  }

  if (cl.option_present ('g'))
  {
    if (! chemical_standardisation.construct_from_command_line (cl, verbose, 'g'))
    {
      usage (6);
    }
  }

  if (cl.option_present ('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will strip to largest fragments\n";
  }

  if (cl.option_present('B'))
  {
    IWString b = cl.string_value('B');
    b.to_lowercase();

    if ("help" == b)
      display_dash_B_qualifiers(cerr);

    for (int i = 0; i < b.length(); i++)
    {
      char bi = b[i];

      if ('u' == bi)
        abbreviated_fingerprints |= ABBREV_UNSATURATION_COUNT;
      else if ('h' == bi)
        abbreviated_fingerprints |= ABBREV_HETEROATOM_COUNT;
      else if ('f' == bi)
        abbreviated_fingerprints |= ABBREV_FUSED_RING_ATOMS;
      else if ('n' == bi)
        abbreviated_fingerprints |= ABBREV_NATOMS;
      else if ('r' == bi)
        abbreviated_fingerprints |= ABBREV_NRINGS;
      else if ('g' == bi)
        abbreviated_fingerprints |= ABBREV_RING_ATOMS;
      else if ('l' == bi)
        abbreviated_fingerprints |= ABBREV_LARGEST_RING_SIZE;
      else if ('a' == bi)
        abbreviated_fingerprints |= ABBREV_AROAMTIC_ATOMS;
      else
      {
        cerr << "Unrecognised -B qualifier '" << b << "'\n";
        cerr << "                           ";
        for (int j = 0; j < i; j++)
        {
          cerr << ' ';
        }
        cerr << "^\n";
        display_dash_B_qualifiers(cerr);
      }
    }

    if (cl.option_present('d'))
    {
      if (! cl.value('d', atom_count_divide) || atom_count_divide < 1)
      {
        cerr << "The atom count divisor (-d) must be a whole +ve number\n";
        usage(2);
      }

      if (verbose)
        cerr << "Atom count related property divisor set to " << atom_count_divide << endl;
    }

    if (cl.option_present('r'))
    {
      if (! cl.value('r', bit_replicates) || bit_replicates < 1)
      {
        cerr << "The number of bit replicates must be a whole +ve number\n";
        usage(2);
      }

      if (verbose)
        cerr << "Will create " << bit_replicates << " bit replicates\n";
    }
  }

  if (cl.option_present ('J'))
  {
    cl.value ('J', temperature_tag);

    if (verbose)
      cerr << "Temperature fingerprints stored as '" << temperature_tag << "' dataitem\n";

    if (! temperature_tag.ends_with ('<'))
      temperature_tag << '<';
  }

  if (cl.option_present ('m'))
  {
    allow_multi_fragment_molecules = 1;
    if (verbose)
      cerr << "Molecules with multiple fragments can be processed\n";
  }

  if (! cl.option_present('a') && ! cl.option_present('J'))
  {
    cerr << "To get output must specify one of -a or -J\n";
    usage(4);
  }

  if (cl.option_present('a') && cl.option_present('J'))
  {
    cerr << "Only one of -a and -J is allowed\n";
    usage(4);
  }

  if (cl.option_present ('a'))
  {
    write_array_of_properties = 1;
    if (verbose)
      cerr << "Molecular properties written as descriptor file\n";
  }

  if (cl.option_present ('i') && cl.option_present ('f'))
  {
    cerr << "The -i (input type) and -f (filter) options are mutually exclusive\n";
    usage (17);
  }

  if (cl.option_present ('f'))
  {
    function_as_filter = 1;
    if (verbose)
      cerr << "Will work as a filter\n";
  }

  int input_type = 0;

  if (cl.option_present ('f'))
    ;
  else if (cl.option_present ('i'))
  {
    if (! process_input_type (cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (! all_files_recognised_by_suffix (cl))
    return 4;

  if (cl.option_present ('u'))
  {
    mpg.set_unsaturation_includes_aromatic(1);
    if (verbose)
      cerr << "Unsaturated aromatic atoms included as unsaturated\n";
  }

  if (cl.option_present ('s'))
  {
    accumulators = new Accumulator_Int<int>[NPROPERTIES];
    if (verbose)
      cerr << "Will gather statistics on themolecules\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (99);
  }

  IWString_and_File_Descriptor output(1);

  if (write_array_of_properties)
    write_descriptor_file_header(output);

  for (int i = 0; i < cl.number_elements(); i++)
  {

	  //GH: Fix bug# 19048
	  //The return value for temperature function will be checked
	  int ret = temperature (cl[i], input_type, output); //Original:   (void) temperature (cl[i], input_type, output);
	  if(0 == ret)
	  {
	       error_code = 0;
	  }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Produced molecular_property fingerprints for " << atom_statistics.n() << " molecules\n";
    if (molecules_with_too_many_atoms)
      cerr << molecules_with_too_many_atoms << " molecules had more than " << std::numeric_limits<unsigned char>::max() << " atoms\n";

    cerr << "Molecules between " << atom_statistics.minval() << " and " << atom_statistics.maxval() << " atoms\n";
    if (atom_statistics.n() > 1)
      cerr << "ave " << atom_statistics.average() << " variance " << atom_statistics.variance();
    cerr << endl;
  }

  if (NULL != accumulators)
  {
    for (int i = 0; i < NPROPERTIES; i++)
    {
      Accumulator_Int<int> & acc = accumulators[i];
      cerr << "Property " << i << " between " << acc.minval() << " and " << acc.maxval();
      if (acc.n() > 1)
        cerr << " ave " << acc.average() << " variance " << acc.variance();
      cerr << endl;
    }

    delete [] accumulators;
  }

  //GH: Fix bug# 19048
  //Return the error code
  return error_code; //original:  return 0
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = temperature (argc, argv);
  //GH: Fix bug# 19048
  //Flip the error return based on Linux convention: 0 is success, otherwise fails
  rc = (0 == rc)?1: 0;
  return rc;
}
