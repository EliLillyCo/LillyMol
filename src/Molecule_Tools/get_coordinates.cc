/*
  Extract the coordinates of atoms matched by a query
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

const char * prog_name = nullptr;

static int verbose = 0;

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;
using std::endl;

static Elements_to_Remove elements_to_remove;

static Charge_Assigner charge_assigner;

static Chemical_Standardisation chemical_standardisation;

static int molecules_read = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static extending_resizable_array<int> queries_hit;

static int ignore_queries_not_matching = 0;

/*
  Apr 99. John Lampe wanted to be able to get out the atoms and
  the queries which hit them to a special file
*/

static std::ofstream John_Lampe_stream;

static int
set_atom_names (const IWString & query_name,
                const Set_of_Atoms & embedding,
                iwaray<IWString> & atom_label)
{
  for (int i = 0; i < embedding.number_elements (); i++)
  {
    atom_number_t j = embedding[i];

    IWString & aname = atom_label[j];
    if (aname.length ())
      aname += ':';

    aname += query_name;
  }

  return 1;
}


static int
set_atom_names (const IWString & query_name,
                const Substructure_Results & sresults,
                iwaray<IWString> & atom_label)
{
  int nhits = sresults.hits_found ();

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * embedding = sresults.embedding (i);
    (void) set_atom_names (query_name, *embedding, atom_label);
  }

  return 1;
}

static int
do_john_lampe_output (const Molecule & m,
                      const iwaray<IWString> & atom_label,
                      std::ostream & output)
{
  int matoms = m.natoms ();

  output << m.name () << '\n';
  output << matoms << " atoms in whole molecule\n";

  int atoms_written = 0;
  for (int i = 0; i < matoms; i++)
  {
    const IWString & aname = atom_label[i];
    if (0 == aname.length ())
      continue;

    const Atom * a = m.atomi (i);

    output << a->x () << ' ' << a->y () << ' ' << a->z () << ' ' << a->atomic_symbol () << ' ' << aname << '\n';
    atoms_written++;
  }

  if (0 == atoms_written && verbose)
    cerr << "Huh, no atoms hit in molecule " << m.name () << endl;

  output << "|\n";

  return output.good ();
}

static int
get_coordinates (const Molecule & m,
                 const Set_of_Atoms * e,
                 std::ostream & output)
{
  int n = e->number_elements ();
  for (int i = 0; i < n; i++)
  {
    atom_number_t j = e->item (i);

    const Atom * a = m.atomi (j);
    output << j << ": " << a->x () << ' ' << a->y () << ' ' << a->z () << '\n';
  }

  return output.good ();
}

static int
get_coordinates (Molecule & m, int nhits,
                 const Substructure_Results & sresults,
                 std::ostream & output)
{
  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding (i);

    (void) get_coordinates (m, e, output);
  }

  return output.good ();
}

static int
get_coordinates (Molecule & m,
                 std::ostream & output)
{
  int nq = queries.number_elements ();

  iwaray<IWString> atom_label;

  if (John_Lampe_stream.rdbuf ()->is_open ())
    atom_label.resize (m.natoms ());
  else
    output << "Name " << m.name () << '\n';

  int rc = 0;

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int tmp = queries[i]->substructure_search (m, sresults);
    if (verbose)
      cerr << ' ' << tmp << " hits to query " << i << endl;

    if (0 == tmp)
    {
      if (ignore_queries_not_matching)
        continue;

      cerr << "Yipes, zero hits to query " << i << " '" << queries[i]->comment () << "'\n";
      return 0;
    }

    if (John_Lampe_stream.rdbuf ()->is_open ())
      (void) set_atom_names (queries[i]->comment (), sresults, atom_label);
    else
      (void) get_coordinates (m, tmp, sresults, output);

    rc++;     // the number of the queries which hit
  }

  queries_hit[rc]++;

  if (John_Lampe_stream.rdbuf ()->is_open ())
    do_john_lampe_output (m, atom_label, John_Lampe_stream);
  else
    output << "|\n";

  return output.good ();
}

static void
preprocess(Molecule & m)
{
  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (elements_to_remove.number_elements ())
    elements_to_remove.process (m);

  if (charge_assigner.active ())
    charge_assigner.process (m);

  return;
}

static int
get_coordinates (data_source_and_type<Molecule> & input, std::ostream & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule ()))
  {
    molecules_read++;

    preprocess(*m);

    std::unique_ptr<Molecule> free_m (m);

    if (! get_coordinates (*m, output))
      return 0;
  }

  return output.good ();

}
static int
get_coordinates (const char * fname, FileType input_type, std::ostream & output)
{
  data_source_and_type<Molecule> input (input_type, fname);
  if (! input.ok ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_verbose (verbose);

  return get_coordinates (input, output);
}

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

  cerr << "Extract coordinates according to a query\n";
  cerr << "Usage: " << prog_name << " <file>\n";
  cerr << "  -q <query>     specify one or more queries\n";
  cerr << "  -s <smarts>    specify query as smarts\n";
  cerr << "  -z             ignore queries not hitting\n";
  cerr << "  -X <symbol>    extract/remove all atoms of type <symbol>. No bonds changed\n";
  cerr << "  -L <fname>     file for John Lampe output\n";
  display_standard_aromaticity_options (cerr);
  display_standard_charge_assigner_options (cerr, 'N');
  cerr << "  -v             verbose output\n";

  exit (rc);
}

static int
get_coordinates (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vi:A:q:zX:s:L:N:g:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (! process_standard_aromaticity_options (cl, verbose))
    usage (6);

  if (cl.option_present ('g'))
  {
    if (! chemical_standardisation.construct_from_command_line (cl, verbose > 1, 'g'))
    {
      usage (6);
    }
  }

  if (cl.option_present ('N'))
  {
    if (! charge_assigner.construct_from_command_line (cl, verbose, 'N'))
    {
      cerr << "Cannot initialise charge assigner (-N option)\n";
      usage (5);
    }
  }

  if (! cl.option_present ('q') && ! cl.option_present ('s'))
  {
    cerr << "Must specify a substructure query via the -q or -s option\n";
    usage (3);
  }

  queries.resize (cl.option_count ('q') + cl.option_count ('s') + 100);

  if (cl.option_present ('q'))
  {
    if (! process_queries (cl, queries, verbose))
    {
      cerr << prog_name << ": cannot process queries from -q option(s)\n";
      return 6;
    }
  }

  if (cl.option_present ('s'))
  {
    const_IWSubstring smarts;
    int i = 0;
    while (cl.value ('s', smarts, i++))
    {
      Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics;
      if (! q->create_from_smarts (smarts))
      {
        cerr << "Cannot parse smarts '" << smarts << "'\n";
        return 62;
      }

      queries.add (q);
    }
  }

  int nq = queries.number_elements ();
  for (int i = 0; i < nq; i++)
  {
    queries[i]->set_find_unique_embeddings_only (1);
  }

  if (cl.option_present ('X'))
  {
    if (! elements_to_remove.construct_from_command_line (cl, verbose, 'X'))
    {
      cerr << "Cannot discern elements to remove from -X switch\n";
      usage (18);
    }
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.option_present ('L'))
  {
    John_Lampe_stream.open (cl.option_value ('L'));

    if (! John_Lampe_stream.good ())
    {
      cerr << "Cannot open John Lampe stream '" << cl.option_value ('L') << "'\n";
      return 21;
    }

    if (verbose)
      cerr << "Molecules in Lampe format written to '" << cl.option_value ('L') << "'\n";
  }

  if (cl.option_present ('z'))
  {
    ignore_queries_not_matching = 1;
    if (verbose)
      cerr << "Will skip queries not matching\n";
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present ('i'))
  {
    if (! process_input_type (cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }

  if (0 == input_type && ! all_files_recognised_by_suffix (cl))
  {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! get_coordinates (cl[i], input_type, std::cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (0 == verbose)
    return rc;

  cerr << "Processed " << molecules_read << " molecules\n";

  if (0 == molecules_read)
    return rc;

  for (int i = 0; i < queries_hit.number_elements (); i++)
  {
    if (queries_hit[i])
      cerr << queries_hit[i] << " queries matched " << i << " of the queries\n";
  }

  for (int i = 0; i < nq; i++)
    queries[i]->report (cerr, verbose - 1);

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = get_coordinates (argc, argv);

  return rc;
}
