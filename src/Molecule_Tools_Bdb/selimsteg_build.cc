/*
  Build the selimsteg infrastructure.
  Input is raw structure file
*/

#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/standardise.h"

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static int reduce_to_largest_fragment = 0;

static Chemical_Standardisation chemical_standardisation;

static int discard_molecule_if_multiple_fragments_larger_than = 0;

static int lower_atom_count_cutoff = 0;
static int upper_atom_count_cutoff = 0;

static IWString_and_File_Descriptor organic_stream;

static IWString_and_File_Descriptor stream_for_discards;

static int molecules_discarded = 0;

static IWString_STL_Hash_Set lly_seen;

static int perform_duplicate_id_checking = 0;

static int duplicate_identifiers_found = 0;

static int discarded_for_too_few_atoms = 0;
static int discarded_for_too_many_atoms = 0;
static int discarded_for_isotopes = 0;
static int discarded_for_bad_valence = 0;
static int discarded_for_non_organic = 0;
static int discarded_for_being_mixture = 0;

static IWString_and_File_Descriptor stream_for_unstandardised_smiles;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Does processing for selimsteg infrastructure - input is raw structures\n";
  cerr << "  -O <fname>    file name for organic subset\n";
  cerr << "  -c <natoms>   lower atom count for organic (based on atoms in largest fragment)\n";
  cerr << "  -C <natoms>   upper atom count for organic (based on atoms in largest fragment)\n";
  cerr << "  -B <fname>    file for discarded smiles\n";
  cerr << "  -d            perform duplicate ID checking\n";
  cerr << "  -x <natoms>   discard if multiple fragments with more than <natoms> (mixture)\n";
  cerr << "  -U <fname>    stream for unstandardised smiles\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -l            strip to largest fragment\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
preprocess (Molecule & m)
{
  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

static int
write_discard_stream_if_open (Molecule & m,
                              const char * reason)
{
  molecules_discarded++;

  if (! stream_for_discards.is_open())
    return 1;

  stream_for_discards << m.smiles() << ' ' << m.name() << ' ' << reason << '\n';
  stream_for_discards.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int ok_elements[HIGHEST_ATOMIC_NUMBER + 1];

static int
contains_non_organics (const Molecule & m)
{
  if (m.organic_only())
    return 0;

  for (int i = m.natoms() - 1; i >= 0; i--)
  {
    const Element * e = m.elementi(i);

//  cerr << "Organic? '" << e->symbol() << "' " << e->organic() << endl;

    if (e->organic())
      continue;

    if (! ok_elements[e->atomic_number()])
      return 1;

    if (m.ncon(i) > 0)   // cannot tolerate covalent inorganics
      return 1;
  }

  return 0;    // no covalently bonded non organics
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

/*
  Only return 0 on a catastrophic error
*/

static int
selimsteg_build (Molecule & m,
                IWString_and_File_Descriptor & output)
{
  if (m.number_isotopic_atoms() > 0)   // always bad, test first
  {
    discarded_for_isotopes++;
    return write_discard_stream_if_open(m, "isotope");
  }

  if (too_many_large_fragments(m))   // need to do before fragment stripping
  {
    discarded_for_being_mixture++;
    return write_discard_stream_if_open(m, "mixture");
  }

  if (contains_non_organics(m))
  {
    discarded_for_non_organic++;
    return write_discard_stream_if_open(m, "nonorganic");
  }

  int matoms;

  if (reduce_to_largest_fragment)
  {
    m.reduce_to_largest_fragment_carefully();
    matoms = m.natoms();
  }
  else
    matoms = m.atoms_in_largest_fragment();

  if (lower_atom_count_cutoff > 0 && matoms < lower_atom_count_cutoff)
  {
    discarded_for_too_few_atoms++;
    return write_discard_stream_if_open(m, "too few atoms");
  }

  if (upper_atom_count_cutoff > 0 && matoms > upper_atom_count_cutoff)
  {
    discarded_for_too_many_atoms++;
    return write_discard_stream_if_open(m, "too many atoms");
  }

  if (! m.valence_ok())   // should this be done before fragment stripping??
  {
    discarded_for_bad_valence++;
    return write_discard_stream_if_open(m, "bad valence");
  }

  output << m.smiles() << ' ' << m.name() << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
is_duplicate_id (const IWString & id)
{
  if (lly_seen.contains(id))
  {
    cerr << "Duplicate identifier '" << id << "'\n";
    duplicate_identifiers_found++;
    return 1;
  }

  lly_seen.insert(id);

  return 0;
}

static int
selimsteg_build (data_source_and_type<Molecule> & input,
                 IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    if (! perform_duplicate_id_checking)
      ;
    else if (is_duplicate_id(m->name()))
      continue;

    if (stream_for_unstandardised_smiles.is_open())
    {
      stream_for_unstandardised_smiles << m->smiles() << ' ' << m->name() << '\n';
      stream_for_unstandardised_smiles.write_if_buffer_holds_more_than(32768);
    }

    preprocess(*m);

    output << m->smiles() << ' ' << m->name() << '\n';

    output.write_if_buffer_holds_more_than(32768);

    if (! selimsteg_build(*m, organic_stream))   // note different output stream
      return 0;
  }

  return 1;
}

static int
selimsteg_build (const char * fname, FileType input_type, 
                IWString_and_File_Descriptor & output)
{
  assert (NULL != fname);

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

  return selimsteg_build(input, output);
}

static int
selimsteg_build (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:lO:c:C:x:B:dU:");

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
      cerr << "Will reduce to largest 'organic' fragment\n";
  }

  if (cl.option_present('x'))
  {
    if (! cl.value('x', discard_molecule_if_multiple_fragments_larger_than) || discard_molecule_if_multiple_fragments_larger_than < 1)
    {
      cerr << "The discard molecules if mixtures of size (-x) must have a value +ve whole number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Molecules considered mixtures if >1 fragment with " << discard_molecule_if_multiple_fragments_larger_than << " atoms\n";
  }

  if (cl.option_present('d'))
  {
    perform_duplicate_id_checking = 1;

    if (verbose)
      cerr << "WIll perform duplicate id checking\n";
  }

  int need_to_check_atom_count_cutoffs = 0;

  if (cl.option_present('c'))
  {
    if (! cl.value('c', lower_atom_count_cutoff) || lower_atom_count_cutoff < 1)
    {
      cerr << "The lower atom count limit (-c) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard molecules having fewer than " << lower_atom_count_cutoff << " atoms\n";

    need_to_check_atom_count_cutoffs = 1;
  }

  if (cl.option_present('C'))
  {
    if (! cl.value('C', upper_atom_count_cutoff) || upper_atom_count_cutoff < 1)
    {
      cerr << "The upper atom count limit (-C) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard molecules having more than " << upper_atom_count_cutoff << " atoms\n";

    need_to_check_atom_count_cutoffs = 1;
  }

  if (! need_to_check_atom_count_cutoffs)
    ;
  else if (lower_atom_count_cutoff < upper_atom_count_cutoff)
    ;
  else if (0 == upper_atom_count_cutoff)
    ;
  else
  {
    cerr << "Inconsistent lower (" << lower_atom_count_cutoff << ") and upper (" << upper_atom_count_cutoff << ") atom count limits\n";
    return 3;
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
  else if (1 == cl.number_elements() && 0 == ::strncmp(cl[0], "-", 1))
    input_type = FILE_TYPE_SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (! cl.option_present('O'))
  {
    cerr << "MUst specify file name for organic subset via the -O option\n";
    usage(3);
  }

  if (cl.option_present('O'))
  {
    IWString o = cl.string_value('O');

    if (! o.ends_with(".smi"))
      o << ".smi";

    if (! organic_stream.open(o.null_terminated_chars()))
    {
      cerr << "Cannot open stream for organic subset '" << o << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Organic subset written to '" << o << "'\n";
  }

  if (cl.option_present('B'))
  {
    IWString b = cl.string_value('B');

    if (! b.ends_with(".smi"))
      b << ".smi";

    if (! stream_for_discards.open(b.null_terminated_chars()))
    {
      cerr << "Cannot open stream for discards '" << b << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Discarded molecules written to '" << b << "'\n";
  }

  if (cl.option_present('U'))
  {
    IWString u = cl.string_value('U');

    if (! u.ends_with(".smi"))
      u << ".smi";

    if (! stream_for_unstandardised_smiles.open(u.null_terminated_chars()))
    {
      cerr << "Cannot open stream for unstandardised smiles '" << u << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Unstandardised smiles written to '" << u << "'\n";
  }

// initialise the array of allowable elements

  set_vector(ok_elements, HIGHEST_ATOMIC_NUMBER + 1, 0);

  ok_elements[6] = 1;
  ok_elements[7] = 1;
  ok_elements[8] = 1;
  ok_elements[9] = 1;
  ok_elements[15] = 1;
  ok_elements[16] = 1;
  ok_elements[17] = 1;
  ok_elements[35] = 1;
  ok_elements[53] = 1;
  ok_elements[3]  = 1;     // Li
  ok_elements[11] = 1;    // Na
  ok_elements[12] = 1;    // Mg
  ok_elements[19] = 1;    // K
  ok_elements[20] = 1;    // Ca

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! selimsteg_build(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();
  if (organic_stream.is_open())
    organic_stream.close();
  if (stream_for_discards.is_open())
    stream_for_discards.close();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    cerr << molecules_discarded << " molecules discarded\n";
    cerr << discarded_for_too_few_atoms<< " discarded for too few atoms\n";
    cerr << discarded_for_too_many_atoms<< " discarded for too many atoms\n";
    cerr << discarded_for_isotopes<< " discarded for isotopes\n";
    cerr << discarded_for_bad_valence<< " discarded for bad valence\n";
    cerr << discarded_for_non_organic<< " discarded for non organic\n";
    cerr << discarded_for_being_mixture<< " discarded for being mixtures\n";
    if (perform_duplicate_id_checking)
      cerr << duplicate_identifiers_found << " duplicate identifiers found\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = selimsteg_build (argc, argv);

  return rc;
}
