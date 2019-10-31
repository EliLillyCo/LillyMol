/*
  Remove matched atoms from molecules
*/

#include <iostream>
#include <memory>
#include <limits>
using std::cerr;
using std::endl;

#include "cmdline.h"
#include "accumulator.h"
#include "misc.h"

#include "istream_and_type.h"
#include "output.h"
#include "molecule.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "substructure.h"
#include "target.h"
#include "spatially_common_matched_atoms.h"

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static int molecules_matching_one_or_more_queries = 0;

static int molecules_written = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int isotope = 0;

#define ISOTOPE_IS_QUERY_ATOM_NUMBER -3

static resizable_array_p<Substructure_Query> queries;

static int max_embeddings_to_process = std::numeric_limits<int>::max();

static int write_parent_molecule = 0;

static int ignore_molecules_not_matching_any_queries = 0;

static int write_molecules_not_matching_any_queries = 0;

static int process_all_queries = 1;

static int invert_selection = 0;

static Molecule capping_group;

static int discern_embedding_via_spatial_location = 0;


static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Removes matched atoms\n";
  cerr << "  -s <smarts>   query as smarts\n";
  cerr << "  -q <qry>      query specification\n";
  cerr << "  -z i          ignore molecules not matching\n";
  cerr << "  -z f          only process first of multiple matches\n";
  cerr << "  -z w          write  molecules not matching\n";
  cerr << "  -p            write the parent molecule for matching molecules\n";
  cerr << "  -r            find one embedding per root atom match\n";
  cerr << "  -b            stop processing queries once there is a match\n";
  cerr << "  -I <iso>      apply isotopic label to remaining attachment points\n";
  cerr << "  -I q          the isotopic label is the matched query atom number\n";
  cerr << "  -C <smiles>   capping group\n";
  cerr << "  -x            invert selection - matched atoms retained\n";
  cerr << "  -w .          (3D) when multiple embeddings, use the one closest to most common location\n";
  cerr << "  -S <stem>     output file name stem\n";
  cerr << "  -o <type>     output type specification\n";
  cerr << "  -i <type>     input  type specification\n";
  cerr << "  -l            reduce to largest fragment\n";
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
do_add_capping_group (Molecule & m,
                      const int * to_be_removed,
                      const Molecule & capping_group)
{
  int rc = 0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (! to_be_removed[i])
      continue;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

      if (to_be_removed[k])
        continue;

      const int initial_atoms = m.natoms();

      m.add_molecule(&capping_group);
      m.add_bond(k, initial_atoms, SINGLE_BOND);
      rc++;
    }
  }

  return rc;
}

static int
do_apply_isotopic_labels (Molecule & m,
                          int * to_be_removed)
{
  int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (! to_be_removed[i])
      continue;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);
      if (to_be_removed[k])
        continue;

      if (isotope != m.isotope(k))
      {
        m.set_isotope(k, isotope);
        rc++;
      }
    }
  }

  return rc;
}

static int
do_apply_isotopic_labels_qam (Molecule & m,
                              const Set_of_Atoms & e,
                              const int * to_be_removed)
{
  int rc = 0;

  const int n = e.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const atom_number_t j = e[i];

    const Atom * aj = m.atomi(j);

    const int jcon = aj->ncon();

    for (int k = 0; k < jcon; ++k)
    {
      const atom_number_t x = aj->other(j, k);

      if (to_be_removed[x])
        continue;

      if ((i+1) != m.isotope(x))
      {
        m.set_isotope(x, i+1);
        rc++;
      }
    }
  }

  return rc;
}

static void
do_invert_selection (int * to_be_removed,
                     const int n)
{
  for (int i = 0; i < n; ++i)
  {
    if (to_be_removed[i])
      to_be_removed[i] = 0;
    else
      to_be_removed[i] = 1;
  }

  return;
}

static int
remove_matched_atoms (Molecule & m,
                      int * to_be_removed,
                      const Set_of_Atoms & e,
                      Molecule_Output_Object & output)
{
  set_vector(to_be_removed, m.natoms(), 0);
  e.set_vector(to_be_removed, 1);

  if (invert_selection)
    do_invert_selection(to_be_removed, m.natoms());

  if (capping_group.natoms() > 0)
    do_add_capping_group(m, to_be_removed, capping_group);

  if (isotope > 0)
    do_apply_isotopic_labels(m, to_be_removed);
  else if (ISOTOPE_IS_QUERY_ATOM_NUMBER == isotope)
    do_apply_isotopic_labels_qam(m, e, to_be_removed);

  m.remove_atoms(to_be_removed);

  molecules_written++;

  return output.write(m);
}

static int
remove_matched_atoms(resizable_array_p<Molecule> & molecules,
                     Substructure_Results * sresults,
                     Molecule_Output_Object & output)
{
  int max_atoms = molecules[0]->natoms();

  const int n = molecules.number_elements();

  for (int i = 1; i < n; ++i)
  {
    const int x = molecules[i]->natoms();

    if (x > max_atoms)
      max_atoms = x;
  }

  int * to_remove = new int[max_atoms]; std::unique_ptr<int[]> free_to_remove(to_remove);

  for (int i = 0; i < n; ++i)
  {
    Molecule & mi = *molecules[i];

    if (write_parent_molecule)
      output.write(mi);

    const Set_of_Atoms * e = sresults[i].embedding(0);

    remove_matched_atoms(mi, to_remove, *e, output);
  }

  return 1;
}

static int
do_discern_embedding_via_spatial_location(data_source_and_type<Molecule> & input,
                                          Molecule_Output_Object & output)
{
  const int n = input.molecules_remaining();

  resizable_array_p<Molecule> molecules;
  molecules.resize(n);

  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules.add(m);
  }

  assert(n == molecules.number_elements());

  Substructure_Results * sresults = new Substructure_Results[n]; std::unique_ptr<Substructure_Results[]> free_sresults(sresults);

  for (int i = 0; i < n; ++i)
  {
    Molecule & mi = *molecules[i];
    const int nhits = queries[0]->substructure_search(mi, sresults[i]);
    if (0 == nhits)
    {
      if (ignore_molecules_not_matching_any_queries)
        continue;

      if (write_molecules_not_matching_any_queries)
        output.write(mi);

      return 0;
    }
  }

  spatially_common_matched_atoms(molecules, sresults);

  return remove_matched_atoms(molecules, sresults, output);
}

static int
remove_matched_atoms (Molecule & m,
                      int * to_be_removed,
                      const Substructure_Results & sresults,
                      int & parent_written,
                      Molecule_Output_Object & output)
{
  if (! parent_written)
  {
    output.write(m);
    molecules_written++;
    parent_written = 1;
  }

  int nhits = sresults.number_embeddings();

  if (nhits > max_embeddings_to_process)
    nhits = max_embeddings_to_process;

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);
    Molecule tmp(m);
    tmp.set_name(m.name());
    remove_matched_atoms(tmp, to_be_removed, *e, output);
  }

  return 1;
}

static int
remove_matched_atoms (Molecule & m,
                      Molecule_Output_Object & output)
{
  int nq = queries.number_elements();

  Molecule_to_Match target(&m);

  int queries_matching_this_molecule = 0;

  int * tmp = new int[m.natoms()]; std::unique_ptr<int[]> free_tmp(tmp);

  int parent_written;
  if (write_parent_molecule)
    parent_written = 0;
  else
    parent_written = 1;

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    if (! remove_matched_atoms(m, tmp, sresults, parent_written, output))
      return 0;

    queries_matching_this_molecule++;

    if (! process_all_queries)
      break;
  }

  if (queries_matching_this_molecule)
  {
    molecules_matching_one_or_more_queries++;
    return 1;
  }

  if (ignore_molecules_not_matching_any_queries)
    return 1;

  if (write_molecules_not_matching_any_queries)
    return output.write(m);

  cerr << "None of " << nq << " queries matched '" << m.name() << "'\n";

  return 0;
}

static int
remove_matched_atoms(data_source_and_type<Molecule> & input,
                     Molecule_Output_Object & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! remove_matched_atoms(*m, output))
      return 0;
  }

  return 1;
}

static int
remove_matched_atoms(const char * fname, int input_type, 
                     Molecule_Output_Object & output)
{
  assert (NULL != fname);

  if (0 == input_type)
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

  if (discern_embedding_via_spatial_location)
    return do_discern_embedding_via_spatial_location(input, output);

  return remove_matched_atoms(input, output);
}

static int
remove_matched_atoms (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lo:S:z:s:q:I:prbC:xw:");

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

  if (cl.option_present('z'))
  {
    int i = 0;
    const_IWSubstring z;
    while (cl.value('z', z, i++))
    {
      if ('i' == z)
      {
        ignore_molecules_not_matching_any_queries = 1;
      }
      else if ('w' == z)
      {
        write_molecules_not_matching_any_queries = 1;
      }
      else if ('f' == z)
      {
        max_embeddings_to_process = 1;
      }
      else
      {
        cerr << "Unrecognised -z directive '" << z << "'\n";
        usage(4);
      }

    }
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
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.option_present('I'))
  {
    const_IWSubstring i = cl.string_value('I');

    if ("q" == i)
      isotope = ISOTOPE_IS_QUERY_ATOM_NUMBER;
    else if (! i.numeric_value(isotope) || isotope <= 0)
    {
      cerr << "The isotopic label to apply (-I) must be a valid isotopic label\n";
      usage(1);
    }
    else if (verbose)
      cerr << "Will apply isotope " << isotope << " to remaining attachment points\n";
  }

  if (cl.option_present('p'))
  {
    write_parent_molecule = 1;

    if (verbose)
      cerr << "Will write the parent molecule of all matched molecules\n";
  }

  int find_one_embedding_per_atom = 0;
  if (cl.option_present('r'))
  {
    find_one_embedding_per_atom = 1;
    if (verbose)
      cerr << "Only one embedding per root atom will be reported\n";
  }

  if (cl.option_present('b'))
  {
    process_all_queries = 0;
    if (verbose)
      cerr << "Will stop matching queries once there is a match\n";
  }

  if (! cl.option_present('q') && ! cl.option_present('s'))
  {
    cerr << "Must specify a substructure query via the -q or -s option\n";
    usage(3);
  }

  queries.resize(cl.option_count('q') + cl.option_count('s') + 100);

  if (cl.option_present('q'))
  {
    if (! process_queries(cl, queries, verbose))
    {
      cerr << prog_name << ": cannot process queries from -q option(s)\n";
      return 6;
    }
  }

  if (cl.option_present('s'))
  {
    const_IWSubstring smarts;
    int i = 0;
    while (cl.value('s', smarts, i++))
    {
      Substructure_Query * q = new Substructure_Query;
      if (! q->create_from_smarts(smarts))
      {
        cerr << "Cannot parse smarts '" << smarts << "'\n";
        return 62;
      }

      queries.add(q);
    }
  }

  int nq = queries.number_elements();

  if (0 == nq)
  {
    cerr << "No queries, cannot continue\n";
    usage(4);
  }

  for (int i = 0; i < nq; i++)
  {
    queries[i]->set_find_unique_embeddings_only(1);
    queries[i]->set_do_not_perceive_symmetry_equivalent_matches(1);

    if (find_one_embedding_per_atom)
      queries[i]->set_find_one_embedding_per_atom(1);
  }

  if (cl.option_present('C'))
  { 
    const_IWSubstring s = cl.string_value('C');

    if (! capping_group.build_from_smiles(s))
    {
      cerr << "invalid capping group smiles '" << s << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Capping group " << capping_group.smiles() << endl;
  }

  if (cl.option_present('x'))
  {
    invert_selection = 1;

    if (verbose)
      cerr << "Will invert selection - matched atoms retained\n";
  }

  if (cl.option_present('w'))
  {
    if (queries.number_elements() > 1)
    {
      cerr << "Sorry, cannot use the -w option with multiple queries, (see Ian)\n";
      return 1;
    }

    discern_embedding_via_spatial_location = 1;

    if (verbose)
      cerr << "When multiple matches, will remove matched atoms closest to most common location\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  Molecule_Output_Object output;

  if (! cl.option_present('o'))
    output.add_output_type(SMI);
  else if (! output.determine_output_types(cl))
  {
    cerr << "Cannot determine output types\n";
    return 0;
  }

  if (cl.option_present('S'))
  {
    const_IWSubstring s = cl.string_value('S');
    if (output.would_overwrite_input_files(cl, s))
    {
      cerr << "Cannot overwrite input file(s) '" << s << "'\n";
      return 19;
    }

    if (! output.new_stem(s))
    {
      cerr << "Cannot open output file(s) '" << s << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Output written to '" << s << "'\n";
  }
  else
  {
    output.new_stem("-");
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! remove_matched_atoms(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    cerr << molecules_matching_one_or_more_queries << " molecules matched one or more queries\n";
    cerr << "wrote " << molecules_written << " molecules\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = remove_matched_atoms(argc, argv);

  return rc;
}
