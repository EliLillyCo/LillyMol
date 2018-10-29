#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "cmdline.h"
#include "misc.h"

#include "aromatic.h"
#include "path.h"
#include "element.h"
#include "target.h"
#include "qry_wstats.h"
#include "istream_and_type.h"
#include "iwstandard.h"
#include "output.h"
#include "rotbond_common.h"

static int verbose = 0;

static int molecules_read = 0;

static int molecules_passing_filters = 0;

static extending_resizable_array<int> nrbonds;  // number of rotatable bonds
static extending_resizable_array<int> lfc;      // longest flexible chain
static extending_resizable_array<int> lird;     // longest inter ring distance

static int min_rotatable_bonds = -1;

static int max_rotatable_bonds = 0;

static int maximum_consecutive_flexible_bonds = 0;

/*
  Single bonds in rings of this size or larger are rotatable
*/

static int single_bonds_in_rings_are_flexible = 0;

/*
  Michal Vieth wanted a list of the rotatable atoms
*/

static std::ofstream stream_for_michal_vieth;

/*
  As we determine a long flexible chain, what kinds of atoms terminate a
  flexible chain. By default, we go through any chain atom, but that
  can be restricted
*/

static int max_ncon_along_flexible_chain = 99;

/*
  Another way of stopping a flexible chain is a double bond branching off
*/

static int stop_flexible_chains_at_unsaturation = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

/*
  We need the first two atoms of the query to define the non-rotatable bond,
  but this leads to very inefficient queries '*-[CD4](-[CH3])(-[CH3])-[CH3]'
  for example.
  We can get dramatically faster queries by leaving out the '*-' and having
  the programme look for an unmatched attachment
*/

static int other_end_is_unmatched_atom = 0;

static Chemical_Standardisation chemical_standardisation;

static int mark_rotatable_bonds = 0;

static Molecule_Output_Object stream_for_molecules_passing;
static Molecule_Output_Object stream_for_molecules_failing;

static int append_results_to_name = 0;

static int append_longest_inter_ring_bond_distance_to_name = 0;

static int maximum_inter_ring_bond_distance = 0;

static int molecules_with_inter_ring_chains_too_long = 0;

static int reduce_to_largest_fragment = 0;

static int use_largest_fragment = 0;

/*
  Sometimes a rotatable bond filter will unfairly penalise small molecules
  that might be OK fragments. We can have an atom count below which things
  just get written
*/

static int always_write_threshold = 0;

static int molecules_passed_by_atom_count = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __TIME__ << " " << __DATE__ << endl;
  cerr << " Does a more careful computation of rotatable bonds in a molecule\n";
  cerr << "  -m <number>    discard molecules with <number> or fewer rotatable bonds\n";
  cerr << "  -M <number>    discard molecules with <number> or more  rotatable bonds\n";
  cerr << "  -f <number>    longest allowable consecutive flexible bonds\n";
  cerr << "  -f <maxc=nn>   maximum connectivity along a long flexible chain\n";
  cerr << "  -f <breaku>    stop any flexible chain growth at any unsaturated atom\n";
  cerr << "  -q <query>     query specification (alternate to -s)\n";
  cerr << "  -s <smarts>    bonds across first two matched atoms are NOT rotatable\n";
  cerr << "  -k             bond is unmatched atom attached to first matched atom\n";
  cerr << "                 can use '-s [CD4](-F)(-F)(-F) -k' rather than\n";
  cerr << "                 can use '-s *-[CD4](-F)(-F)(-F)' which is much slower\n";
  cerr << "  -r <number>    place isotope <number> at the ends of rotatable bonds\n";
  cerr << "  -R <size>      single bonds in rings of size <size> are rotatable\n";
  cerr << "  -a             append rotatable bond count to name\n";
  cerr << "  -c             append longest distance between rings to name\n";
  cerr << "  -C <bonds>     discard molecules with longest distance between rings > <bonds>\n";
  cerr << "  -w <natoms>    write molecules with <natoms> or fewer atoms regardless\n";
  cerr << "  -i <type>      specify input file type. Enter '-i help' for details\n";
  cerr << "  -o <type>      specify output file type(s)\n";
  cerr << "  -S <string>    create output files with name stem <string>\n";
  cerr << "  -B <string>    write discarded molecules to <string>\n";
  cerr << "  -V <fname>     list of rotatable bonds written to <fname>\n";
  cerr << "  -e             compute rotbonds on largest fragment only\n";
  cerr << "  -l             reduce to largest fragment (discard counterions)\n";
  cerr << "  -E <symbol>    create an element with symbol <symbol>\n";
  cerr << "  -E autocreate  automatically create new elements when encountered\n";
  (void) display_standard_aromaticity_options(cerr);
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static void
do_append_results_to_name (Molecule & m,
                           int rotbond,
                           int longest_inter_ring_distance)
{
  IWString tmp = m.name();

  tmp << ' ' << " ROTBOND = " << rotbond;

  if (append_longest_inter_ring_bond_distance_to_name)
    tmp << " LINTRBD " << longest_inter_ring_distance;

  m.set_name(tmp);

  return;
}

static int
bonds_to_nearest_ring (const Molecule & m,
                       atom_number_t previous_atom,
                       atom_number_t current_atom,
                       const Atom * const * atom,
                       const int * fsid)
{
  const Atom * a = atom[current_atom];

  int acon = a->ncon();

  int matoms = m.natoms();

  int rc = matoms;    // we look for the shortest path

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(current_atom, i);

    if (previous_atom == j)
      continue;

    if (fsid[j] > 0)
      return 1;

    int tmp = bonds_to_nearest_ring(m, current_atom, j, atom, fsid);

    if (tmp <= 0)
      continue;

    tmp++;

    if (tmp < rc)
      rc = tmp;
  }

  if (matoms == rc)    // didn't find anything
    return -1;

  return rc;
}

static int
compute_between_ring_descriptors (Molecule & m)
{
  int nr = m.nrings();

  if (nr < 2)    // we deal with between ring information
    return 0;

  int matoms = m.natoms();

  int * fsid = new int[matoms]; std::unique_ptr<int[]> free_fsid(fsid);

  if (1 == m.label_atoms_by_ring_system(fsid))   // more than 1 ring, but 1 ring system
    return 0;

  const Atom ** atom = new const Atom *[matoms]; std::unique_ptr<const Atom *[]> free_atom(atom);
  m.atoms(atom);

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (0 == fsid[i])    // not in a ring
      continue;

    const Atom * ai = atom[i];

    int acon = ai->ncon();
    if (2 == acon)
      continue;

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = ai->item(j);

      if (b->nrings())
        continue;

      atom_number_t k = b->other(i);

      int tmp;
      if (fsid[k])
        tmp = 1;
      else
        tmp = bonds_to_nearest_ring(m, i, k, atom, fsid);

      if (tmp > rc)
        rc = tmp;
    }
  }

  return rc;
}

static int
marked_as_non_rotatable (const Bond * b,
                         const resizable_array_p<Bond> & non_rotatable_bonds)
{
  atom_number_t a1 = b->a1();
  atom_number_t a2 = b->a2();

  int nr = non_rotatable_bonds.number_elements();

  for (int i = 0; i < nr; i++)
  {
    const Bond * nrbi = non_rotatable_bonds[i];

    if (a1 == nrbi->a1() && a2 == nrbi->a2())
      return 1;

    if (a1 == nrbi->a2() && a2 == nrbi->a1())
      return 1;
  }

  return 0;
}
//#define DEBUG_COUNT_FLEXIBLE_BONDS

static int
count_flexible_bonds (Molecule & m,
                      atom_number_t previous_atom,
                      atom_number_t current_atom,
                      int * already_done)
{
#ifdef DEBUG_COUNT_FLEXIBLE_BONDS
  cerr << "count_flexible_bonds from " << previous_atom << " to " << current_atom << endl;
#endif

  if (already_done[current_atom])
    return 1;

  already_done[current_atom] = 1;

  if (m.is_ring_atom(current_atom))
    return 1;

  const Atom * ac = m.atomi(current_atom);

  int acon = ac->ncon();

  if (1 == acon)
    return 1;

  if (acon > max_ncon_along_flexible_chain)
    return 1;

  if (stop_flexible_chains_at_unsaturation && ac->nbonds() > acon)
    return 1;

  int rc = 1;     // bond from previous_atom to current_atom continues the chain

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = ac->item(i);

    if (! b->is_single_bond())
      continue;

    atom_number_t j = b->other(current_atom);

    if (previous_atom == j)
      continue;

//  If we run up against something non-rotatable, we are done

    if (already_done[j] || m.is_ring_atom(j))
    {
      if (1 == rc)
        rc = 2;       // we did get this one flexible bond
      continue;
    }

    int down_here = count_flexible_bonds(m, current_atom, j, already_done) + 1;

    if (down_here > rc)
      rc = down_here;
  }

  return rc;
}

static int
identify_long_flexible_chains (Molecule & m,
                               int * already_done)
{
  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (already_done[i])
      continue;

    const Atom * ai = m.atomi(i);

    if (2 != ai->ncon())
      continue;

    if (2 != ai->nbonds())
      continue;

    if (m.is_ring_atom(i))
      continue;

    int f1 = count_flexible_bonds(m, i, ai->other(i, 0), already_done);
    int f2 = count_flexible_bonds(m, i, ai->other(i, 1), already_done);

    int total_length = f1 + f2 - 1;    //

//  cerr << "From atom " << i << " f1 = " << f1 << " f2 = " << f2 << " total " << total_length << endl;

    if (total_length > rc)
      rc = total_length;
  }

  return rc;
}

/*
  Mark the atoms at the ends of non-rotatable bonds as done
*/

static int
identify_long_flexible_chains (Molecule & m,
                               const resizable_array_p<Bond> & non_rotatable_bonds,
                               int * already_done)
{
  int nb = non_rotatable_bonds.number_elements();

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = non_rotatable_bonds[i];
    already_done[b->a1()] = 1;
    already_done[b->a2()] = 1;
  }

  return identify_long_flexible_chains(m, already_done);
}

static int
identify_long_flexible_chains (Molecule & m,
                               const resizable_array_p<Bond> & non_rotatable_bonds)
{
  int * already_done = new_int(m.natoms()); std::unique_ptr<int[]> free_already_done(already_done);

  return identify_long_flexible_chains(m, non_rotatable_bonds, already_done);
}

/*
  Is the ring bond B in a ring of size at least SINGLE_BONDS_IN_RINGS_ARE_FLEXIBLE
*/

static int
in_ring_of_size (Molecule & m,
                 const Bond * b,
                 int single_bonds_in_rings_are_flexible)
{
  atom_number_t a1 = b->a1();
  atom_number_t a2 = b->a2();

  int nr = m.nrings();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (ri->number_elements() < single_bonds_in_rings_are_flexible)    // we are looking for a ring with at least SINGLE_BONDS_IN_RINGS_ARE_FLEXIBLE atoms
      continue;

    if (ri->is_aromatic())
      continue;

    if (! ri->contains(a1))
      continue;

    if (! ri->contains(a2))
      continue;

    return 1;
  }

  return 0;
}

static int
identify_largest_fragment (Molecule & m)
{
  int nf = m.number_fragments();

  if (1 == nf)
    return 0;

  int atoms_in_largest_fragment = m.atoms_in_fragment(0);
  int rc = 0;

  for (int i = 1; i < nf; i++)
  {
    int j = m.atoms_in_fragment(i);

    if (j > atoms_in_largest_fragment)
    {
      atoms_in_largest_fragment = j;
      rc = i;
    }
  }

  return rc;
}

static atom_number_t
do_find_unmatched_atom (const Molecule & m,
                        const Set_of_Atoms & embedding)
{
  atom_number_t a1 = embedding[0];

  const Atom * a = m.atomi(a1);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(a1, i);

    if (! embedding.contains(j))
      return j;
  }

  cerr << "Yipes, no unmatched atoms bonded to atom " << a1 << endl;


  return INVALID_ATOM_NUMBER;
}

static int
compute_rotatable_bonds (Molecule & m,
                         resizable_array_p<Bond> & non_rotatable_bonds,
                         int * isotope)
{
  int nb = m.nedges();

  non_rotatable_bonds.resize(nb);

  (void) m.ring_membership();    // force sssr determination

  int nq = queries.number_elements();

  if (nq)
  {
    Molecule_to_Match target(&m);

    for (int i = 0; i < nq; i++)
    {
      Substructure_Results sresults;
  
      int nhits = queries[i]->substructure_search(target, sresults);
  
      for (int j = 0; j < nhits; j++)
      {
        const Set_of_Atoms * e = sresults.embedding(j);
  
        assert (e->number_elements() >= 2);
  
        atom_number_t a1 = e->item(0);

        atom_number_t a2;
        if (other_end_is_unmatched_atom)
          a2 = do_find_unmatched_atom(m, *e);
        else
          a2 = e->item(1);
  
        Bond * b = new Bond(a1, a2, SINGLE_BOND);
  
        non_rotatable_bonds.add(b);
      }
    }
  }

  int fragment_number_to_process;

  if (use_largest_fragment)
    fragment_number_to_process = identify_largest_fragment(m);
  else
    fragment_number_to_process = -1;

  int rc = 0;

  if (stream_for_michal_vieth.rdbuf()->is_open())
    stream_for_michal_vieth << m.name() << endl;

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);

    if (! b->is_single_bond())
      continue;

    if (0 == b->nrings())    // chain bonds are always flexible
      ;
    else if (0 == single_bonds_in_rings_are_flexible)    // ring bonds are rigid
      continue;
    else if (b->is_aromatic())       // aromatic bonds are always rigid
      continue;
    else if (! in_ring_of_size(m, b, single_bonds_in_rings_are_flexible))
      continue;

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (-1 == fragment_number_to_process)
      ;
    else if (fragment_number_to_process != m.fragment_membership(a1))
      continue;

    if (1 == m.ncon(a1))
      continue;

    if (1 == m.ncon(a2))
      continue;

    if (triple_bond_at_either_end(m, b))
      continue;

    if (marked_as_non_rotatable(b, non_rotatable_bonds))
      continue;

    if (part_of_otherwise_non_rotabable_entity(m, a1, a2))
      continue;

//  Great, we have a rotatable bond

    if (NULL != isotope)
    {
      isotope[a1] = mark_rotatable_bonds;
      isotope[a2] = mark_rotatable_bonds;
    }

    if (stream_for_michal_vieth.rdbuf()->is_open())
    {
      stream_for_michal_vieth << ' ' << m.atomi(a1)->atomic_symbol() << (a1 + 1);
      stream_for_michal_vieth << ' ' << m.atomi(a2)->atomic_symbol() << (a2 + 1);
      stream_for_michal_vieth << endl;
    }

    rc++;
  }

  return rc;
}

int
compute_rotatable_bonds (Molecule & m,
                         resizable_array_p<Bond> & non_rotabable_bonds)
{
  int * isotope;

  if (mark_rotatable_bonds)
    isotope = new_int(m.natoms());
  else
    isotope = NULL;

  int rc = compute_rotatable_bonds(m, non_rotabable_bonds, isotope);

  if (NULL != isotope)
  {
    m.set_isotopes(isotope);

    delete [] isotope;
  }

  return rc;
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment ();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

static int
rotatable_bonds (Molecule & m)
{
  preprocess(m);

  resizable_array_p<Bond> non_rotabable_bonds;

  int rotbond = compute_rotatable_bonds(m, non_rotabable_bonds);

  int longest_flexible_chain = 0;

  if (maximum_consecutive_flexible_bonds > 1 && m.natoms() > maximum_consecutive_flexible_bonds)
    longest_flexible_chain = identify_long_flexible_chains(m, non_rotabable_bonds);

  int longest_inter_ring_distance = -1;

  if (append_longest_inter_ring_bond_distance_to_name || maximum_inter_ring_bond_distance > 0)
  {
    longest_inter_ring_distance = compute_between_ring_descriptors(m);
    lird[longest_inter_ring_distance]++;
  }

  if (verbose > 1)
  {
    cerr << m.name() << " has " << rotbond << " rotatable bonds";
    if (maximum_consecutive_flexible_bonds > 1)
      cerr << ", longest chain " << longest_flexible_chain;
    cerr << endl;
  }

  if (append_results_to_name)
    do_append_results_to_name(m, rotbond, longest_inter_ring_distance);

  nrbonds[rotbond]++;
  lfc[longest_flexible_chain]++;

  int matoms = m.natoms();

  int pass = 1;

  if (matoms <= always_write_threshold)
  {
    pass = 1;
    molecules_passed_by_atom_count++;
  }
  if (maximum_inter_ring_bond_distance > 0 && longest_inter_ring_distance > maximum_inter_ring_bond_distance)
  {
    molecules_with_inter_ring_chains_too_long++;
    pass = 0;
  }

  if (matoms <= always_write_threshold)
    ;
  else if (rotbond <= min_rotatable_bonds)
    pass = 0;

  if (matoms <= always_write_threshold)
    ;
  else if (max_rotatable_bonds > 0 && rotbond >= max_rotatable_bonds)
    pass = 0;

  if (matoms <= always_write_threshold)
    ;
  else if (pass && maximum_consecutive_flexible_bonds > 0 && longest_flexible_chain >= maximum_consecutive_flexible_bonds)
    pass = 0;

  if (pass)
  {
    molecules_passing_filters++;
    if (stream_for_molecules_passing.active())
      stream_for_molecules_passing.write(m);
  }
  else if (stream_for_molecules_failing.active())
    stream_for_molecules_failing.write(m);

  return 1;
}

static int 
rotatable_bonds (data_source_and_type<Molecule> & input)
{
  Molecule * m;

  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    if (single_bonds_in_rings_are_flexible)
      m->compute_aromaticity();

    if (! rotatable_bonds(*m))
      return 0;
  }

  return 1;
}

static int
rotatable_bonds (const char * fname, int input_type)
{
  data_source_and_type<Molecule> input(input_type, fname);

  if (0 == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return rotatable_bonds(input);
}


static int
rotatable_bonds (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:g:q:s:m:M:i:o:S:B:r:R:alef:kV:cC:w:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "unrecognised_options_encountered\n";
    usage(3);
  }

  verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose > 1, 'A'))
  {
    cerr << "Cannot process standard aromaticity options (-A)\n";
    usage(4);
  }

  if (! process_elements(cl, verbose > 1, 'E'))
  {
    cerr << "Cannot initialise elements (-E option)\n";
    usage(2);
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot initialise chemical standardisation\n";
      usage(6);
    }
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', mark_rotatable_bonds) || mark_rotatable_bonds < 1)
    {
      cerr << "The mark rotatable bonds (-r) option must be followed by a whole positive number\n";
      usage(5);
    }

    if (verbose)
      cerr << "Rotatable bonds labelled with " << mark_rotatable_bonds << endl;
  }

  if (cl.option_present('R'))
  {
    if (! cl.value('R', single_bonds_in_rings_are_flexible) || single_bonds_in_rings_are_flexible < 3)
    {
      cerr << "The smallest ring size for flexible single bonds must be > 2\n";
      usage(5);
    }

    if (verbose)
      cerr << "Single bonds in rings of size " << single_bonds_in_rings_are_flexible << " can be flexible\n";
  }

  if (cl.option_present('q'))
  {
    if (! process_queries(cl, queries, verbose > 1, 'q'))
    {
      cerr << "Cannot read queries (-q option)\n";
      usage(5);
    }
  }

  if (cl.option_present('s'))
  {
    int i = 0;
    const_IWSubstring smarts;
    while (cl.value('s', smarts, i++))
    {
      Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics;
      if (! q->create_from_smarts(smarts))
      {
        delete q;
        return 60 + i;
      }

      if (verbose > 1)
        cerr << "Created query from smarts '" << smarts << "'\n";

      queries.add(q);
    }
  }

  int nq = queries.number_elements();

  if (verbose)
    cerr << "Defined " << nq << " queries specifying non-rotatable bonds\n";

// Should put in a check to ensure that each query has a minimum of 2 atoms

  if (nq)
  {
    for (int i = 0; i < nq;i++)
    {
      queries[i]->set_find_unique_embeddings_only(1);

//    put in check for 2 atoms
    }
  }

  if (cl.option_present('k'))
  {
    other_end_is_unmatched_atom = 1;

    if (verbose)
      cerr << "Non-rotatable bonds are across first matched atom and an unmatched atom\n";
  }

  if (cl.option_present('m'))
  {
    if (! cl.value('m', min_rotatable_bonds) || min_rotatable_bonds < 0)
    {
      cerr << "The minimum number of rotatable bonds must be a non-negative whole number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard molecules having " << min_rotatable_bonds << " or fewer rotatable bonds\n";
  }

  if (cl.option_present('M'))
  {
    if (! cl.value('M', max_rotatable_bonds) || max_rotatable_bonds < 1 || max_rotatable_bonds < min_rotatable_bonds)
    {
      cerr << "The minimum number of rotatable bonds must be a non-negative whole number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard molecules having " << max_rotatable_bonds << " or more rotatable bonds\n";
  }

  if (cl.option_present('f'))
  {
    int i = 0;
    const_IWSubstring f;
    while (cl.value('f', f, i++))
    {
      if (f.starts_with("maxc="))
      {
        f.remove_leading_chars(5);

        if (! f.numeric_value(max_ncon_along_flexible_chain) || max_ncon_along_flexible_chain < 2)
        {
          cerr << "The maximum connectivity of atoms along a flexible chain must be a whole positive number > 1\n";
          usage(18);
        }

        if (verbose)
          cerr << "Long flexible chains stop with atoms of connectivity more than " << max_ncon_along_flexible_chain << endl;
      }
      else if ("breaku" == f)
      {
        stop_flexible_chains_at_unsaturation = 1;

        if (verbose)
          cerr << "Flexible chains stop at any double or triple bond\n";
      }
      else if (! f.numeric_value(maximum_consecutive_flexible_bonds)  || maximum_consecutive_flexible_bonds < 1 || maximum_consecutive_flexible_bonds < min_rotatable_bonds)
      {
        cerr << "The maximum consecutive flexible bonds (-f option) must be a positive whole number\n";
        usage(4);
      }
    }

    if (0 == maximum_consecutive_flexible_bonds)
    {
      cerr << "Must specify the maximum number of consecutive flexible bonds\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will discard molecules having " << maximum_consecutive_flexible_bonds << " or more consecutive rotatable bonds\n";

    if (verbose && 99 != max_ncon_along_flexible_chain)
      cerr << "Long flexible chains stop at atoms with " << max_ncon_along_flexible_chain << " or more connections\n";
  }

  if (cl.option_present('a'))
  {
    append_results_to_name = 1;

    if (verbose)
      cerr << "Rotatable bond count appended to name\n";
  }

  if (cl.option_present('c'))
  {
    append_longest_inter_ring_bond_distance_to_name = 1;

    if (verbose)
      cerr << "Will append longest inter ring bond distance to name\n";
  }

  if (cl.option_present('C'))
  {
    if (! cl.value('C', maximum_inter_ring_bond_distance) || maximum_inter_ring_bond_distance < 1)
    {
      cerr << "The maximum permitted inter ring bond distance (-C) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will discard molecules having inter ring bond distances > " << maximum_inter_ring_bond_distance << " bonds\n";
  }

// If they haven't requested any filtering, then just append the number

  if (min_rotatable_bonds < 0 && 0 == max_rotatable_bonds)
    append_results_to_name = 1;

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('e'))
  {
    use_largest_fragment = 1;

    if (verbose)
      cerr << "Computation will be based on the largest fragment\n";
  }

  if (reduce_to_largest_fragment && use_largest_fragment)
  {
    cerr << "The -l (chop to largest fragment) and -e (use largest fragment) options are mutually exclusive\n";
    usage(3);
  }

  if (cl.option_present('w'))
  {
    if (! cl.value('w', always_write_threshold) || always_write_threshold < 1)
    {
      cerr << "The always write atom count threshold (-w) option must be a valid whole number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will always write molecules with " << always_write_threshold << " atoms or fewer\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  IWString s;

  if (cl.option_present('S'))
  {
    s = cl.string_value('S');

    if (! cl.option_present('o'))
      stream_for_molecules_passing.add_output_type(SMI);
    else if (! stream_for_molecules_passing.determine_output_types(cl))
    {
      cerr << "Cannot determine output type(s) for -S stream\n";
      return 6;
    }

    if (stream_for_molecules_passing.would_overwrite_input_files(cl, s))
    {
      cerr << "The -S file cannot overwrite its input\n";
      return 8;
    }

    if (! stream_for_molecules_passing.new_stem(s))
    {
      cerr << "Cannot set -S option to '" << s << "'\n";
      return 7;
    }

    if (verbose)
      cerr << "Molecules passing constraints written to '" << s << "'\n";
  }

// Is there a file for rejections

  if (cl.option_present('B'))
  {
    IWString b = cl.string_value('B');

    if (b == s)
    {
      cerr << "Passing and failing files must be different '" << b << "' vs '" << s << "'\n";
      usage(18);
    }

    if (! cl.option_present('o'))
      stream_for_molecules_failing.add_output_type(SMI);
    else if (! stream_for_molecules_failing.determine_output_types(cl))
    {
      cerr << "Cannot determine output type(s) for -B stream\n";
      return 6;
    }

    if (stream_for_molecules_failing.would_overwrite_input_files(cl, b))
    {
      cerr << "The -B file cannot overwrite its input\n";
      return 8;
    }

    if (! stream_for_molecules_failing.new_stem(b))
    {
      cerr << "Cannot set -B option to '" << b << "'\n";
      return 7;
    }

    if (verbose)
      cerr << "Molecules failing constraints written to '" << b << "'\n";
  }

  if (! cl.option_present('S') && ! cl.option_present('B'))
  {
    stream_for_molecules_passing.add_output_type(SMI);

    stream_for_molecules_passing.new_stem("-");
  }

  if (cl.option_present('V'))
  {
    IWString v = cl.string_value('V');

    stream_for_michal_vieth.open(v.null_terminated_chars(), std::ios::out);

    if (! stream_for_michal_vieth.good())
    {
      cerr << "Cannot open Vieth file '" << v << "'\n";
      return 10;
    }

    if (verbose)
      cerr << "A list of rotatable bonds written to '" << v << "'\n";
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
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! rotatable_bonds(cl[i], input_type))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules, " << molecules_passing_filters << " passed the constraints\n";

    int sum = 0;
    for (int i = 0; i < nrbonds.number_elements(); i++)
    {
      if (nrbonds[i])
      {
        cerr << nrbonds[i] << " molecules had " << i << " rotatable bonds\n";
        sum += i * nrbonds[i];
      }
    }

    float ave = static_cast<float>(sum) / static_cast<float>(molecules_read);

    cerr << "Ave Rot Bond " << ave << endl;

    if (maximum_consecutive_flexible_bonds > 0)
    {
      for (int i = 0; i < lfc.number_elements(); i++)
      {
        if (lfc[i])
          cerr << lfc[i] << " molecules had flexible chains of length " << i << endl;
      }
    }

    for (int i = 0; i < lird.number_elements(); i++)
    {
      if (lird[i])
        cerr << lird[i] << " molecules had longest inter ring distances of " << i << " bonds\n";
    }

    if (always_write_threshold > 0)
      cerr << molecules_passed_by_atom_count << " molecules passed by having " << always_write_threshold << " or fewer atoms\n";

  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = rotatable_bonds(argc, argv);

  return rc;
}
