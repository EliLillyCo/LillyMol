/*
  Implementation of John Rose and Scott Frank complexity of
  chiral centres.

  Natural product related chiral centres are considered benign,
  but those not in that form are considered more difficult.
  Especially difficult are those with adjacent chiral centres
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static resizable_array_p<Substructure_Hit_Statistics> amino_acid;

static extending_resizable_array<int> chiral_centre_count;

static extending_resizable_array<int> amino_acid_chiral_count;
static extending_resizable_array<int> adjacent_chiral_centres;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Counts chiral centres of various forms, from John Rose, Scott Frank\n";
  cerr << "  Amino acid type chiral centres\n";
  cerr << "  Other types\n";
  cerr << "  Adjacent chiral centres\n";
  cerr << "  -q ...        specify amino acid queries via query files\n";
  cerr << "  -s ...        specify amino acid queries via smarts\n";
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

/*
*/

static int
complex_chirality (Molecule & m,
                   IWString_and_File_Descriptor & output)
{
  output << m.name() << ' ';

  int nc = m.chiral_centres();

  chiral_centre_count[nc]++;

  if (0 == nc)
  {
    output << "0 0 0\n";
    return 1;
  }

  int matoms = m.natoms();

// First identify all the chiral centres that are amino acids

  int * is_amino_acid = new_int(matoms); std::unique_ptr<int[]> free_is_amino_acid(is_amino_acid);

  int nq = amino_acid.number_elements();

  int amino_acid_chiral_count_this_molecule = 0;

  Molecule_to_Match target(&m);

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = amino_acid[i]->substructure_search(target, sresults);

    if (verbose > 2 && nhits > 0)
      cerr << m.name() << " hits amino acid query\n";

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      int ne = e->number_elements();

      for (int k = 0; k < ne; k++)
      {
        //int l = e->item(k);

        if (NULL != m.chiral_centre_at_atom(k))
        {
          is_amino_acid[k] = 1;
          amino_acid_chiral_count_this_molecule++;
        }
      }
    }
  }

  amino_acid_chiral_count[amino_acid_chiral_count_this_molecule]++;

  int adjacent_chiral_centres_this_molecule = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (is_amino_acid[i])  // adjacent to these does not count
      continue;

    if (NULL == m.chiral_centre_at_atom(i))
      continue;

    const Atom * ai = m.atomi(i);

    int acon = ai->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = ai->other(i, j);
      if (k < i)
        continue;

      if (is_amino_acid[k])
        continue;

      if (NULL == m.chiral_centre_at_atom(k))   // no adjacent chiral centre
        continue;

      adjacent_chiral_centres_this_molecule++;
    }
  }

  output << nc << ' ' << amino_acid_chiral_count_this_molecule << ' ' << adjacent_chiral_centres_this_molecule << '\n';

  adjacent_chiral_centres[adjacent_chiral_centres_this_molecule]++;

  return output.good ();
}

static int
complex_chirality (data_source_and_type<Molecule> & input,
                IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! complex_chirality(*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
complex_chirality (const char * fname, FileType input_type, 
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

  return complex_chirality(input, output);
}

static int
get_smarts (Command_Line & cl,
            const char flag,
            resizable_array_p<Substructure_Hit_Statistics> & queries)
{
  if (! cl.option_present (flag))
    return 1;

  int i = 0;
  const_IWSubstring s;
  while (cl.value (flag, s, i++))
  {
    Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics;

    if (! q->create_from_smarts (s))
    {
      cerr << "Cannot parse smarts '" << s << "'\n";
      return 31;
    }

    queries.add (q);
  }

  return queries.number_elements ();
}

static int
complex_chirality (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:ls:q:");

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

  if (cl.option_present ('q') && ! process_queries (cl, amino_acid, verbose))
  {
    cerr << "Cannot process queries from -q option(s)\n";
    return 6;
  }

  if (! get_smarts (cl, 's', amino_acid))
  {
    cerr << "Cannot parse -s option\n";
    usage (42);
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

  IWString_and_File_Descriptor output(1);

  output << "ID Nchiral Amino Adjacent\n";

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! complex_chirality(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    for (int i = 0; i < chiral_centre_count.number_elements(); i++)
    {
      if (chiral_centre_count[i])
        cerr << chiral_centre_count[i] << " molecules had " << i << " chiral centres\n";
    }
    for (int i = 0; i < amino_acid_chiral_count.number_elements(); i++)
    {
      if (amino_acid_chiral_count[i])
        cerr << amino_acid_chiral_count[i] << " molecules had " << i << " amino acid chiral centres\n";
    }
    for (int i = 0; i < adjacent_chiral_centres.number_elements(); i++)
    {
      if (adjacent_chiral_centres[i])
        cerr << adjacent_chiral_centres[i] << " molecules had " << i << " adjacent chiral centres\n";
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = complex_chirality (argc, argv);

  return rc;
}
