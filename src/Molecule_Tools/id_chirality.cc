/*
  Third party acquisitions wants a measure of chirality.
  We append a count of the number of chiral centres to the name
*/

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/chiral_centre.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;
static int verbose = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static IWString tag;

static int create_unique_names = 0;

static int molecules_read = 0;

static extending_resizable_array<int> total_chiral_count;
static extending_resizable_array<int> marked_chiral_centres;
static extending_resizable_array<int> unmarked_chiral_centres;

static int isotope_for_unmarked_centres = 0;

/*
  We can get info about which centres to look for
*/

static resizable_array_p<Substructure_Hit_Statistics> queries;

/*
  Optionally we can suppress processing of certain types
*/

static resizable_array_p<Substructure_Hit_Statistics> ignore_queries;

/*
  We can remove chiral centres if needed
*/

static resizable_array_p<Substructure_Query> remove_chiral_centres_on;

/*
  For efficiency we have some hard coded queries
*/

static int hard_coded_queries = 0;

static int hard_coded_query_cd4 = 0;
static int hard_coded_query_cd3h = 0;

static int do_permutations = 8;

static int only_permute_unlabelled_centres = 1;

static int check_that_marked_chiral_centres_are_actually_chiral = 0;

static int remove_invalid_chiral_centres = 0;

static int invalid_chiral_centres_found = 0;

static int max_allowed_unmarked_centres = -1;

static int max_allowed_chiral_centres = -1;

static Molecule_Output_Object stream_for_rejected_molecules;

static int chiral_centres_removed_by_rmchiral;
static int molecules_with_chiral_centres_removed_by_rmchiral;

static int unique_molecules_only = 1;

static IWString prepend_before_sequence('.');

static IWString_and_File_Descriptor stream_for_chirality_count;

static int perceive_ring_chirality = 1;

static int just_write_chiral_centre_count = 0;

/*
  During processing we need to distinguish between those atoms that
  were already chiral (in the connection table) and those that
  have been detected by us
*/

#define EXCLUDED_FROM_CONSIDERATION -1
#define UNMARKED_CHIRAL_ATOM 1
#define MARKED_CHIRAL_ATOM 2

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
// clang-format off
  cerr << "Enumerates chiral variants, possibly based on one or more queries\n";
  cerr << "  -C <string>    chiral centre count appended as 'string:nn' (default '" << tag << "')\n";
  cerr << "                 appends total chiral count and number unlabelled\n";
  cerr << "  -q <query>     queries for which atoms to process\n";
  cerr << "  -s <smarts>    smarts for which atoms to process\n";
  cerr << "  -h .....       hard coded queries - for speed\n";
  cerr << "  -h CD4         aliphatic 4 connected Carbon\n";
  cerr << "  -h CD3H        aliphatic 3 connected Carbon with an implicit Hydrogen\n";
  cerr << "  -X <query>     queries to NOT process\n";
  cerr << "  -p <isotope>   isotopically label unmarked centres with <isotope>\n";
  cerr << "  -P <number>    create permutations for molecules with <number> or fewer\n";
  cerr << "                 unlabelled centres (default 8)\n";
  cerr << "  -b             permute marked as well as unmarked chiral centres\n";
  cerr << "  -u             generate unique names for enumerated stereo-isomers\n";
  cerr << "  -I check       look for atoms that are marked chiral but are not\n";
  cerr << "  -I remove      remove invalid chiral centres\n";
  cerr << "  -R ...         query specification. Remove chiral centres on matched atoms\n";
  cerr << "  -y             do NOT perceive unmarked chiral centres in rings\n";
  cerr << "  -m <num>       reject molecules with more than <num> unmarked centres\n";
  cerr << "  -n <num>       reject molecules with more than <num> chiral centres\n";
  cerr << "  -B <file>      write rejected molecules to <file>\n";
  cerr << "  -S <stem>      write output to <stem>\n";
  cerr << "  -D <fname>     write chirality count descriptor to <fname>\n";
  cerr << "  -j             only write the chirality count descriptors to the -D file - no enumeration\n";
  cerr << "  -i <type>      input type specification\n";
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -o ...         specify output behaviour\n";
  cerr << "  -Y <string>    prepend <string> before tautomer number\n";
  cerr << "  -E ...         standard element options, enter '-E help' for info\n";
  (void) display_standard_aromaticity_options(cerr);
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -v             verbose output\n";
// clang-format on

  exit(rc);
}

static void
do_append_chiral_count_to_name (Molecule & m,
                                int total_chiral_centres,
                                int unlabelled_chiral_centres_found)
{
  IWString tmp = m.name();

  tmp += tag;

  tmp += total_chiral_centres;

  tmp += ' ';

  tmp += unlabelled_chiral_centres_found;

  m.set_name(tmp);

  return;
}

/*
  First task is to create any chiral centres needed. Note that they are
  created with arbitrary directionality
*/

static int
add_missing_chiral_centres (Molecule & m,
                            const int * is_chiral)
{
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (! is_chiral[i] || EXCLUDED_FROM_CONSIDERATION == is_chiral[i])
      continue;

    if (nullptr != m.chiral_centre_at_atom(i))
      continue;

    Chiral_Centre * c = new Chiral_Centre(i);

    const Atom * ai = m.atomi(i);

    int icon = ai->ncon();

    for (int j = 0; j < icon; j++)
    {
      atom_number_t k = ai->other(i, j);

      if (0 == j)
        c->set_left_down(k);
      else if (1 == j)
        c->set_right_down(k);
      else if (2 == j)
        c->set_top_front(k);
      else if (3 == j)
        c->set_top_back(k);
    }

    if (3 == icon)
      c->set_top_back(kChiralConnectionIsImplicitHydrogen);

    if (! c->complete() || ! m.valid_chiral_centre(c))
    {
      cerr << "id_chirality:yipes, incomplete or invalid chiral centre at atom " << i << " in " << m.name() << "\n";
      c->debug_print(cerr);
      return 0;
    }

    m.add_chiral_centre(c);
    rc++;
  }

  return rc;
}

static int
write_chiral_variant_unique_name_ (Molecule & m, 
                                   int & count,
                                   Molecule_Output_Object & output)
{
  IWString tmp = m.name();

  tmp += prepend_before_sequence;
  tmp << count;

  count++;

  m.set_name(tmp);

  return output.write(m);
}

static int 
write_chiral_variant_unique_name (Molecule & m,
                                  int & count,
                                  Molecule_Output_Object & output)
{
  IWString mname = m.name();

  int rc = write_chiral_variant_unique_name_(m, count, output);

  m.set_name(mname);

  return rc;
}

static int
write_chiral_variant (Molecule & m,
                      int & count,
                      IW_STL_Hash_Set & already_created,
                      Molecule_Output_Object & output)
{
  if (unique_molecules_only)
  {
    if (already_created.contains(m.unique_smiles()))
      return 1;

    already_created.insert(m.unique_smiles());

    m.invalidate_smiles();
  }

  if (create_unique_names)
    return write_chiral_variant_unique_name(m, count, output);

  return output.write(m);
}

static int
create_permutations (Molecule & m,
                     atom_number_t a,
                     const int * is_chiral,
                     int chiral_centres_this_molecule,
                     int depth,
                     int & count,
                     IW_STL_Hash_Set & already_created,
                     Molecule_Output_Object & output)
{
  if (depth == chiral_centres_this_molecule)
    write_chiral_variant(m, count, already_created, output);
  else
  {
    int matoms = m.natoms();

    for (int i = a + 1; i < matoms; i++)
    {
      if (EXCLUDED_FROM_CONSIDERATION == is_chiral[i])
        continue;

      if (is_chiral[i])
        create_permutations(m, i, is_chiral, chiral_centres_this_molecule, depth + 1, count, already_created, output);
    }
  }

  if (only_permute_unlabelled_centres && MARKED_CHIRAL_ATOM == is_chiral[a])
    ;
  else
    m.invert_chirality_on_atom(a);

  if (depth == chiral_centres_this_molecule)
    write_chiral_variant(m, count, already_created, output);
  else
  {
    int matoms = m.natoms();

    for (int i = a + 1; i < matoms; i++)
    {
      if (0 == is_chiral[i] || EXCLUDED_FROM_CONSIDERATION == is_chiral[i])
        continue;

      if (only_permute_unlabelled_centres && UNMARKED_CHIRAL_ATOM != is_chiral[i])
        continue;

      create_permutations(m, i, is_chiral, chiral_centres_this_molecule, depth + 1, count, already_created, output);
    }
  }

  return output.good();
}

static int
create_permutations (Molecule & m,
                     const int * is_chiral,
                     int chiral_centres_this_molecule,
                     int unlabelled_chiral_centres_found,
                     Molecule_Output_Object & output)
{
  int matoms = m.natoms();

#ifdef DEBUG_CREATE_PERMUTATIONS
  for (int i = 0; i < matoms; i++)
  {
    cerr << "Atom " << i << " '" << m.smarts_equivalent_for_atom(i) << "' " << is_chiral[i] << endl;
  }
#endif

  if (only_permute_unlabelled_centres)
    chiral_centres_this_molecule = count_occurrences_of_item_in_array(UNMARKED_CHIRAL_ATOM, matoms, is_chiral);

  IW_STL_Hash_Set already_created;

  int count = 0;     // used for assigning unique numbers to the variants

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (0 == is_chiral[i] || EXCLUDED_FROM_CONSIDERATION == is_chiral[i])
      continue;

    if (only_permute_unlabelled_centres && MARKED_CHIRAL_ATOM == is_chiral[i])
      continue;

//  cerr << "Creating permuations for atom " << i <<endl;

    rc += create_permutations(m, i, is_chiral, chiral_centres_this_molecule, 1, count, already_created, output);
  }

  if (0 == rc)
    return output.write(m);

  return rc;
}

/*
  Common code for the last things done before identifying a chiral centre
*/

static int
process_possibly_chiral_atom (Molecule & m,
                              atom_number_t a,
                              int * is_chiral)
{
  if (! is_actually_chiral(m, a))
    return 0;;

  is_chiral[a] = UNMARKED_CHIRAL_ATOM;

  if (0 != isotope_for_unmarked_centres)
    m.set_isotope(a, isotope_for_unmarked_centres);

  return 1;
}

static int
possibly_chiral (Molecule & m, atom_number_t i)
{
  atomic_number_t z = m.atomic_number(i);

// Only Carbon, Nitrogen or Sulphur can be chiral here

  if (6 == z)
    ;
  else if (7 == z)
    ;
  else if (16 == z)
    ;
  else
    return 0;

  int icon = m.ncon(i);

  if (icon < 2)
    return 0;

  int ih = m.implicit_hydrogens(i);

  if (ih > 1)
    return 0;

//cerr << "Atom " << i << " type " << z << " icon = " << icon << " ih = " << ih << endl;

  if (6 == z)
  {
    if (4 == icon)
      return 1;
    if (3 == icon && 1 == ih)
      return 1;

    return 0;
  }
  else if (7 == z)
  {
    if (4 == icon)
      return 1;

    if (3 == icon && 1 == m.hcount(i) && 1 == m.formal_charge(i))
      return 1;

    return 0;
  }

  int lp;
  if (! m.lone_pair_count(i, lp))
    return 0;

  if (lp > 1)
    return 0;

  if (16 == z)
  {
    if (3 == icon && 4 == m.nbonds(i) && 1 == lp)
      return 1;

    return 0;
  }

  return 1;
}

static int
identify_existing_chiral_centres (Molecule & m,
                                  int * is_chiral)
{
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (nullptr == m.chiral_centre_at_atom(i))
      continue;

    if (EXCLUDED_FROM_CONSIDERATION == is_chiral[i])
    {
      m.remove_chiral_centre_at_atom(i);
      continue;
    }

    if (check_that_marked_chiral_centres_are_actually_chiral && ! is_actually_chiral(m, i))
    {
      invalid_chiral_centres_found++;

      cerr << "Atom " << i << " (" << m.atomic_symbol(i) << ' ' << m.ncon(i) << " connections) is not actually chiral\n";

      if (remove_invalid_chiral_centres)
        m.remove_chiral_centre_at_atom(i);

      continue;
    }

    rc++;
    is_chiral[i] = MARKED_CHIRAL_ATOM;
  }

  marked_chiral_centres[rc]++;

  return rc;
}

static int
do_hard_coded_queries (Molecule & m,
                       int * is_chiral)
{
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (is_chiral[i])
      continue;

    const Atom * ai = m.atomi(i);

    if (6 != ai->atomic_number())
      continue;

    if (hard_coded_query_cd4 && 4 == ai->ncon())
      ;
    else if (hard_coded_query_cd3h && 3 == ai->ncon() && 1 == m.implicit_hydrogens(i))
      ;
    else
      continue;

    if (process_possibly_chiral_atom(m, i, is_chiral))
      rc++;
  }

  return rc;
}

static int
assign_atoms_to_ignore (Molecule & m,
                        int * is_chiral)
{
  int nq = ignore_queries.number_elements();

  int rc = 0;

  Molecule_to_Match target(&m);

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = ignore_queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      e->set_vector(is_chiral, EXCLUDED_FROM_CONSIDERATION);
    }

    rc++;
  }

  return rc;
}

static int
examine_query_hits (Molecule & m,
                    int * is_chiral)
{
  int nq = queries.number_elements();

  int rc = 0;

  Molecule_to_Match target(&m);

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      atom_number_t k = e->item(0);

      if (is_chiral[k])
        continue;

      if (process_possibly_chiral_atom(m, k, is_chiral))
        rc++;
    }
  }

  return rc;
}

static int
id_chirality_by_substructure_searches(Molecule & m,
                                       int * is_chiral)
{
  int unlabelled_chiral_centres_found = 0;      // the number of chiral centres we identify

  int nq = queries.number_elements();

  if (nq)
    unlabelled_chiral_centres_found += examine_query_hits(m, is_chiral);

  if (hard_coded_queries)
    unlabelled_chiral_centres_found += do_hard_coded_queries(m, is_chiral);

  unmarked_chiral_centres[unlabelled_chiral_centres_found]++;

  return unlabelled_chiral_centres_found;
}

static int
id_chirality_by_rules (Molecule & m,
                       int * is_chiral)
{
  int matoms = m.natoms();

  int unlabelled_chiral_centres_found = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (is_chiral[i])
      continue;

    if (! possibly_chiral(m, i))
      continue;

    if (process_possibly_chiral_atom(m, i, is_chiral))
      unlabelled_chiral_centres_found++;
  }

  return unlabelled_chiral_centres_found;
}

/*
  Routines for dealing with chirality and queries lifted directly from fileconv. Should
  be shared, but the static fileconv variables inhibit this.
*/

static int
identify_matched_atoms_with_chiral_centres (const Molecule & m,
                                const Set_of_Atoms & e,
                                Set_of_Atoms & atoms_with_chiral_centres_to_be_removed)
{
  int n = e.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = e[i];

    if (nullptr == m.chiral_centre_at_atom(j))    // atom J does not have a chiral centre
      continue;
      
    if (atoms_with_chiral_centres_to_be_removed.contains(j))
      continue;

    atoms_with_chiral_centres_to_be_removed.add(j);
    rc++;
  }

  return rc;
}


static int
do_remove_chiral_centres_on_matched_atoms (Molecule & m,
                                           const resizable_array_p<Substructure_Query> & remove_chiral_centres_on)
{
  int nchiral = m.chiral_centres();

  if (0 == nchiral)
    return 0;

  Molecule_to_Match target(&m);

  int nq = remove_chiral_centres_on.number_elements();

  Set_of_Atoms atoms_with_chiral_centres_to_be_removed;

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = remove_chiral_centres_on[i]->substructure_search(target, sresults);

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      identify_matched_atoms_with_chiral_centres(m, *e, atoms_with_chiral_centres_to_be_removed);
    }
  }

  int rc = atoms_with_chiral_centres_to_be_removed.number_elements();

  if (0 == rc)
    return 0;

  molecules_with_chiral_centres_removed_by_rmchiral++;
  chiral_centres_removed_by_rmchiral += rc;

  if (verbose > 2)
    cerr << m.name() << " had " << rc << " chiral centres removed by rmchiral\n";

  for (int i = 0; i < rc; i++)
  {
    atom_number_t j = atoms_with_chiral_centres_to_be_removed[i];

    m.remove_chiral_centre_at_atom(j);
  }

  return rc;
}

static void
identify_ring_atoms (Molecule & m,
                     int * is_chiral)
{
  const int * ring_membership = m.ring_membership();

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (ring_membership[i])
      is_chiral[i] = EXCLUDED_FROM_CONSIDERATION;
  }

  return;
}

static int
id_chirality(Molecule & m,
             int * is_chiral,
             Molecule_Output_Object & output)
{
  if (remove_chiral_centres_on.number_elements())
    do_remove_chiral_centres_on_matched_atoms(m, remove_chiral_centres_on);
    
  if (ignore_queries.number_elements())
    assign_atoms_to_ignore(m, is_chiral);

  if (! perceive_ring_chirality)
    identify_ring_atoms(m, is_chiral);

  int chiral_centres_already_present = identify_existing_chiral_centres(m, is_chiral);

  int unlabelled_chiral_centres_found;

  if (queries.number_elements() || hard_coded_queries)
    unlabelled_chiral_centres_found = id_chirality_by_substructure_searches(m, is_chiral);
  else
    unlabelled_chiral_centres_found = id_chirality_by_rules(m, is_chiral);

  int total_chiral_centres = chiral_centres_already_present + unlabelled_chiral_centres_found;

  total_chiral_count[total_chiral_centres]++;

  if (stream_for_chirality_count.is_open())
  {
    static const char sep = ' ';
    append_first_token_of_name(m.name(), stream_for_chirality_count);
    stream_for_chirality_count << sep << total_chiral_centres;
    stream_for_chirality_count << sep << unlabelled_chiral_centres_found;
    stream_for_chirality_count << sep << chiral_centres_already_present;
    stream_for_chirality_count << '\n';
    stream_for_chirality_count.write_if_buffer_holds_more_than(4096);

    if (just_write_chiral_centre_count)
      return 1;
  }

  if (0 == total_chiral_centres)
  {
    if (tag.length())
      do_append_chiral_count_to_name(m, 0, 0);

    return output.write(m);
  }

  if (max_allowed_unmarked_centres >= 0 && unlabelled_chiral_centres_found > max_allowed_unmarked_centres)
  {
    if (stream_for_rejected_molecules.active())
      return stream_for_rejected_molecules.write(m);

    return 1;
  }

  if (max_allowed_chiral_centres >= 0 && total_chiral_centres > max_allowed_chiral_centres)
  {
    if (stream_for_rejected_molecules.active())
      return stream_for_rejected_molecules.write(m);

    return 1;
  }

  if (unlabelled_chiral_centres_found && do_permutations)
    add_missing_chiral_centres(m, is_chiral);

  if (tag.length())
    do_append_chiral_count_to_name(m, total_chiral_centres, unlabelled_chiral_centres_found);

// Make sure we get some kind of output

  if (verbose > 2)
  {
    for (int i = 0; i < m.natoms(); i++)
    {
      cerr << "atom " << i << " chirality " << is_chiral[i] << endl;
    }
  }

  if (do_permutations)
  {
    int permutations_available_this_molecule = unlabelled_chiral_centres_found;
    if (! only_permute_unlabelled_centres)
      permutations_available_this_molecule += chiral_centres_already_present;

    if (verbose > 2)
      cerr << permutations_available_this_molecule << " permutations available in '" << m.name() << "'\n";

    if (permutations_available_this_molecule <= do_permutations)
      return create_permutations(m, is_chiral, total_chiral_centres, unlabelled_chiral_centres_found, output);
  }

  return output.write(m);
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
id_chirality (Molecule & m,
              Molecule_Output_Object & output)
{
  int * tmp = new_int(m.natoms()); std::unique_ptr<int[]> free_tmp(tmp);

  return id_chirality(m, tmp, output);
}

static int
id_chirality (data_source_and_type<Molecule> & input,
              Molecule_Output_Object & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    preprocess(*m);

    if (! id_chirality(*m, output))
      return 0;
  }

  return 1;
}

static int
id_chirality (const char * fname, FileType input_type,
              Molecule_Output_Object & output)
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
    cerr << "cannot open '" << fname << "'\n";
    return 1;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return id_chirality(input, output);
}

static int
id_chirality (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:i:o:C:g:E:S:p:q:s:P:auI:h:m:n:B:X:lbR:Y:D:yj");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose))
  {
    usage(5);
  }

  if (! process_elements(cl, verbose > 1, 'E'))
  {
    usage(6);
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose, 'g'))
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

  if (cl.option_present('p'))
  {
    if (! cl.value('p', isotope_for_unmarked_centres) || isotope_for_unmarked_centres < 1)
    {
      cerr << "The isotopic label for unlabelled centres must be > 0\n";
      usage(6);
    }

    if (verbose)
      cerr << "Unlabelled chiral centres isotopically labelled as " << isotope_for_unmarked_centres << endl;
  }

  if (cl.option_present('P'))
  {
    if (! cl.value('P', do_permutations) || do_permutations < 1)
    {
      cerr << "The permutations (-P option) value must be a whole number > 0\n";
      usage(13);
    }

    if (verbose)
      cerr << "Will generate permutations of up to " << do_permutations << " unmarked chiral centres\n";

    if (cl.option_present('b'))
    {
      only_permute_unlabelled_centres = 0;
      if (verbose)
        cerr << "Will also permute marked centres\n";
    }
  }

  queries.resize(cl.option_count('q') + cl.option_count('s'));

  if (cl.option_present('q'))
  {
    if (! process_queries(cl, queries, verbose > 1, 'q'))
    {
      cerr << "Cannot process queries specified by -q option\n";
      usage(11);
    }
  }

  if (cl.option_present('s'))
  {
    IWString s;
    int i = 0;
    while (cl.value('s', s, i++))
    {
      Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics;
      if (! q->create_from_smarts(s))
      {
        cerr << "INvalid smarts '" << s << "'\n";
        return i;
      }

      queries.add(q);
    }
  }

  if (verbose && queries.number_elements())
    cerr << "defined " << queries.number_elements() << " queries to define chiral centres\n";

  if (cl.option_present('X'))
  {
    if (! process_queries(cl, ignore_queries, verbose > 1, 'X'))
    {
      cerr << "Cannot process exclude atom queries (-X option)\n";
      usage(4);
    }

    if (verbose)
      cerr << "Defined " << ignore_queries.number_elements() << " ignore queries\n";
  }

  if (cl.option_present('R'))
  {
    if (! process_queries(cl, remove_chiral_centres_on, verbose > 1, 'R'))
    {
      cerr << "Cannot read remove chiral atoms query specification (-R option)\n";
      return 5;
    }

    if (verbose)
      cerr << "Created " << remove_chiral_centres_on.number_elements() << " queries for removing chiral centres\n";
  }

  if (cl.option_present('h'))
  {
    int i = 0;
    const_IWSubstring h;
    while (cl.value('h', h, i++))
    {
      if ("all" == h)
      {
        hard_coded_query_cd4 = 1;
        hard_coded_query_cd3h = 1;
        hard_coded_queries = 1;

        if (verbose)
          cerr << "Three and four connected carbons examined\n";
      }
      else if ("CD4" == h)
      {
        hard_coded_query_cd4 = 1;
        hard_coded_queries = 1;

        if (verbose)
          cerr << "Four connected carbon atoms examined\n";
      }
      else if ("CD3H" == h)
      {
        hard_coded_query_cd3h = 1;
        hard_coded_queries = 1;

        if (verbose)
          cerr << "Three connected carbon atoms examined\n";
      }
      else
      {
        cerr << "Unrecognised -h qualifier '" << h << "'\n";
        usage(18);
      }
    }
  }

  if (cl.option_present('C'))
  {
    tag = cl.string_value('C');

    if (verbose)
      cerr << "Append string is '" << tag << "'\n";

// Add a leading space to tag and append a colon

    IWString tmp = ' ';
    tmp += tag;
    tmp += ':';
  
    tag = tmp;
  }

  if (cl.option_present('u'))
  {
    if (! cl.option_present('P'))
    {
      cerr << "For the -u option to make sense you need to generate permutations (-P option)\n";
      usage(17);
    }

    create_unique_names = 1;

    if (verbose)
      cerr << "Will create unique names for each stereo version\n";
  }

  if (cl.option_present('I'))
  {
    const_IWSubstring I;
    int i = 0;
    while (cl.value('I', I, i++))
    {
      if ("check" == I)
      {
        check_that_marked_chiral_centres_are_actually_chiral = 1;
        if (verbose)
          cerr << "Will look for invalid chiral centres\n";
      }
      else if ("remove" == I)
      {
        remove_invalid_chiral_centres = 1;
        check_that_marked_chiral_centres_are_actually_chiral = 1;

        if (verbose)
          cerr << "Will remove invalid chiral centres\n";
      }
      else
      {
        cerr << "Unrecognised -I qualifier '" << I << "'\n";
        usage(11);
      }
    }
  }

  if (cl.option_present('Y'))
  {
    cl.value('Y', prepend_before_sequence);

    if (verbose)
      cerr << "Will insert '" << prepend_before_sequence << "' before sequence number\n";
  }

  if (cl.option_present('y'))
  {
    perceive_ring_chirality = 0;

    if (verbose)
      cerr << "Will NOT perceive chirality of ring atoms\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
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
  else if (1 == cl.number_elements() && 0 == ::strncmp(cl[0], "-", 1))
    input_type = FILE_TYPE_SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;


  Molecule_Output_Object output;

  if (! cl.option_present('o'))
  {
    output.add_output_type(FILE_TYPE_SMI);
  }
  else if (! output.determine_output_types(cl))
  {
    cerr << "Cannot determine output type(s)\n";
    usage(4);
  }

  if (! cl.option_present('S'))
  {
    output.new_stem("-");
  }
  else
  {
    IWString s = cl.string_value('S');

    if (output.would_overwrite_input_files(cl, s))
    {
      cerr << "Cannot overwrite input file(s), stem '" << s << "'\n";
      return 6;
    }

    if (! output.new_stem(s))
    {
      cerr << "Cannot open output stream '" << s << "'\n";
      return 4;
    }
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', max_allowed_unmarked_centres) || max_allowed_unmarked_centres < 0)
    {
      cerr << "The maximum allowed number of unmarked chiral centres (-m option) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard molecules having more than " << max_allowed_unmarked_centres << " unmarked chiral centres\n";
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', max_allowed_chiral_centres) || max_allowed_chiral_centres < 0)
    {
      cerr << "The maximum allowed number of chiral centres (-n option) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard molecules having more than " << max_allowed_chiral_centres << " chiral centres\n";
  }

  if (max_allowed_unmarked_centres >= 0 || max_allowed_chiral_centres >= 0)
  {
    if (cl.option_present('B'))
    {
      if (! cl.option_present('o'))
      {
        stream_for_rejected_molecules.add_output_type(FILE_TYPE_SMI);
      }
      else if (! stream_for_rejected_molecules.determine_output_types(cl))
      {
        cerr << "Cannot determine output type(s) for rejections file\n";
        usage(4);
      }

      IWString b = cl.string_value('B');
      if (output.would_overwrite_input_files(cl, b))
      {
        cerr << "Cannot overwrite input file(s), stem '" << b << "'\n";
        return 6;
      }

      if (! stream_for_rejected_molecules.new_stem(b))
      {
        cerr << "Cannot open rejected molecule stream '" << b << "'\n";
        return 4;
      }

      if (verbose)
        cerr << "Rejected molecules written to '" << b << "'\n";
    }
  }

  if (cl.option_present('D'))
  {
    const char * d = cl.option_value('D');
    if (! stream_for_chirality_count.open(d))
    {
      cerr << "Cannot open stream for chirality descriptor '" << d << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Chirality count descriptor for each molecule written to '" << d << "'\n";

    const char sep = ' ';
    stream_for_chirality_count << "ID" << sep << "poss_chiral" << sep << "unmarked" << sep << "marked\n";

    if (cl.option_present('j'))
    {
      just_write_chiral_centre_count = 1;

      if (verbose)
        cerr << "Will write the chirality count descriptor and not do enumeration\n";
    }
  }

  int rc = 0;

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! id_chirality(cl[i], input_type, output))
    {
      rc = i + 1;
      break;
    }
  }
  
  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules";

    if (check_that_marked_chiral_centres_are_actually_chiral)
      cerr << ", " << invalid_chiral_centres_found << " invalid chiral centres";

    cerr << endl;

    for (int i = 0; i < marked_chiral_centres.number_elements(); i++)
    {
      if (marked_chiral_centres[i])
        cerr << marked_chiral_centres[i] << " molecules had " << i << " explicit chiral centres\n";
    }

    for (int i = 0; i < unmarked_chiral_centres.number_elements(); i++)
    {
      if (unmarked_chiral_centres[i])
        cerr << unmarked_chiral_centres[i] << " molecules had " << i << " unmarked chiral centres\n";
    }

    for (int i = 0; i < total_chiral_count.number_elements(); i++)
    {
      if (0 == i || total_chiral_count[i] > 0)
      {
        float f = static_cast<float>(total_chiral_count[i]) / static_cast<float>(molecules_read);

        cerr << total_chiral_count[i] << " molecules had " << i << " chiral centres " << f << "\n";
      }
    }

    if (remove_chiral_centres_on.number_elements())
    {
      cerr << molecules_with_chiral_centres_removed_by_rmchiral << " molecules with chiral centres removed by -R queries\n";
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = id_chirality(argc, argv);

  return rc;
}
