/*
  We often want to identify elongated molecules,

C(=O)(C1=CC(=CC(=C1)F)F)N(C)CC1CN(CCC2=CC=C(OC)C=C2)CCC1 PBCHM45199241
C(=O)(NC1=CC=C(C=C1)OC)CN1CCN(CC1)CC1=CC2=C(C=C1)OCO2 PBCHM1312610
C1(=NC2=C(S1)C=CC=C2)C1=CC=C(S1)C(=O)OC(C)C(=O)NC1=CC=CC(=C1)Cl PBCHM4537610
N1(C2=CC=C(NC(=O)CC3=CC=CC(=C3)Cl)C=C2)CCC(N[C@@H](C)CC2=CN=CS2)CC1 PBCHM26327783
*/

#include <stdlib.h>
#include <assert.h>
#include <memory>
#include <limits>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static int molecules_passing = 0;
static int molecules_failing = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 1;   // always

static extending_resizable_array<int> longest_path;

static int minimum_longest_path = 0;

static int molecules_shorter_than_min_longest_path = 0;

/* 
  Sometimes we might target only long molecules that have rings in them
*/

static int min_ring_systems = 0;

static int molecules_with_too_few_ring_systems = 0;

/*
  We probably want to insist on taking para substituted traverses through
  rings and ring systems
*/

static int must_take_longest_path_through_rings = 1;

static int max_deviation_from_longest_path = -1;

static double max_average_deviation_from_longest_path = std::numeric_limits<double>::max();

static int everything_smaller_passes = 0;
static int molecules_passed_for_too_few_atoms = 0;

static int isotopically_label_failing_molecules = 0;

static int pass_if_fused_rings_present = 0;

static int molecules_passed_for_having_fused_ring_system = 0;

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
  cerr << "  -m <nbonds>   minimum longest path, all molecules shorter than this pass\n";
  cerr << "  -r <nrsys>    minimum number of ring systems for a failure\n";
  cerr << "  -d <nbonds>   max allowable average distance from longest path\n";
  cerr << "  -D <nbonds>   max allowable distance from longest path\n";
  cerr << "  -x            allow non-max path length paths through rings\n";
  cerr << "  -w <natoms>   all molecules with fewer than <natoms> pass\n";
  cerr << "  -f            all molecules with fused ring systems pass\n";
  cerr << "  -I            isotopically label rejected molecules (requires -U)\n";
  cerr << "  -P <fname>    write passing molecules to <fname>\n";
  cerr << "  -F <fname>    write failing molecules to <fname>\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
preprocess(Molecule & m)
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
contains_fused_rings(Molecule & m,
                     const int * ring_system_membership)
{
  const int nr = m.nrings();

  if (nr < 2)
    return 0;

  for (int i = 0; i < nr; i++)
  {
    if (m.ringi(i)->is_fused())
      return 1;
  }

// No rings have normal fusion, but what about spiro fusion

  int matoms = m.natoms();

  int * tmp = new_int(matoms); std::unique_ptr<int[]> free_tmp(tmp);

  m.ringi(0)->set_vector(tmp, 1);

  for (int i = 1; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (ri->any_members_set_in_array(tmp))
      return 1;

    ri->set_vector(tmp, 1);
  }

  return 0;
}

//#define DEBUG_IDENTIFY_ALL_ATOMS_ON_LONGEST_PATHS

static int
identify_all_atoms_on_longest_paths (Molecule & m,
                                     atom_number_t a1,
                                     atom_number_t a2,     // our target
                                     int lp,
                                     int * on_longest_path)
{
#ifdef DEBUG_IDENTIFY_ALL_ATOMS_ON_LONGEST_PATHS
  cerr << "identify_all_atoms_on_longest_paths:continue with atom " << a1 << " dist " << lp << " end at " << a2 << endl;
#endif

  on_longest_path[a1] = lp;

  const Atom * a = m.atomi(a1);

  const int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(a1, i);

#ifdef DEBUG_IDENTIFY_ALL_ATOMS_ON_LONGEST_PATHS
    cerr << "identify_all_atoms_on_longest_paths:from " << a1 << " continue to " << j << " olp " << on_longest_path[j] << " d2end " << (j == a2 ? 0 : m.bonds_between(j, a2)) << endl;
#endif

    if (on_longest_path[j])
      continue;

    if (j == a2)         // we are done
      return 1;

    if (lp - 1 == m.bonds_between(j, a2))
      identify_all_atoms_on_longest_paths(m, j, a2, lp - 1, on_longest_path);
  }

  return 1;
}

static int
do_pass (Molecule & m,
         Molecule_Output_Object & output)
{
  molecules_passing++;

  if (output.active())
    return output.write(m);

  return 1;
}

static int
do_fail (Molecule & m,
         const Set_of_Atoms & a1,
         const Set_of_Atoms & a2,
         int lp,
         Molecule_Output_Object & output)
{
  molecules_failing++;

  if (output.active())
  {
    if (isotopically_label_failing_molecules)   // choose a pair at random
    {
      m.set_isotope(a1[0], lp);
      m.set_isotope(a2[0], lp);
    }
    return output.write(m);
  }

  return 1;
}

static int
same_atom_or_bonded(const Molecule & m,
                    const atom_number_t a1,
                    const atom_number_t a2)
{
  if (a1 == a2)
    return 1;

  return m.are_bonded(a1, a2);
}

/*
  Consider
  C1(=NC2=CC=CC=C2O1)NC1=NCN(CN1)CC1=CC2=C(C=C1)OCO2 PBCHM1265241

  Which has four paths of length 15 through it. But these are
  largely the same. We need to identify those cases where
  the ends of the longest paths are bonded to each other.

  Let's see if the identified ends are bonded to each other.
*/

static int
ends_of_longest_paths_joined(Molecule & m,
                             const Set_of_Atoms & equivalent_a1,
                             const Set_of_Atoms & equivalent_a2)
{
  const int n = equivalent_a1.number_elements();

  if (n > 4)   // hard to imagine what is happening there
    return 0;

  assert (n == equivalent_a2.number_elements());

  atom_number_t a1 = equivalent_a1[0];
  atom_number_t a2 = equivalent_a2[0];

  for (int i = 1; i < n; i++)
  {
    atom_number_t e1 = equivalent_a1[i];
    atom_number_t e2 = equivalent_a2[i];

    if (same_atom_or_bonded(m, a1, e1) && same_atom_or_bonded(m, a2, e2))
      ;
    else if (same_atom_or_bonded(m, a1, e2) && same_atom_or_bonded(m, a2, e1))
      ;
    else
      return 0;
  }

  return 1;
}

//#define DEBUG_IS_LONGEST_PATH_THROUGH_RINGS

/*
  As the longest path crosses a ring system, does it take the longest
  path through that ring system.
  
  If we get to a ring atom, jump through the ring system to the
  atom furthest from our entry to that ring system
*/

static int
is_longest_path_through_rings(Molecule & m,
                              const atom_number_t a1,
                              const atom_number_t a2,
                              const int lp,
                              const int * on_longest_path,
                              const int * ring_system_membership)
{
  assert (on_longest_path[a1]);
  assert (on_longest_path[a2]);
  assert (a1 != a2);

  const Atom * ai = m.atomi(a1);

  int acon = ai->ncon();

  atom_number_t next_atom = INVALID_ATOM_NUMBER;

#ifdef DEBUG_IS_LONGEST_PATH_THROUGH_RINGS
  cerr << "LP starts at " << lp << " a1 = " << a1 << " olp[a1] " << on_longest_path[a1] << ", a2 = " << a2 << endl;
#endif

  int next_is_across_ring_bond = 0;   // need this for biphenyl type linkages

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = ai->item(i);

    const atom_number_t j = b->other(a1);

#ifdef DEBUG_IS_LONGEST_PATH_THROUGH_RINGS
    cerr << "Connected to " << j << " olp " << on_longest_path[j] << endl;
#endif

    if (j == a2)   // reached target
      return 1;

    if (0 == on_longest_path[j])
      continue;

#ifdef DEBUG_IS_LONGEST_PATH_THROUGH_RINGS
    cerr << "To atom " << j << " olp " << on_longest_path[j] << endl;
#endif

    if ((lp - 1) != on_longest_path[j])  // not moving towards end
      continue;

    next_atom = j;
    next_is_across_ring_bond = b->nrings();
    break;
  }

  assert (INVALID_ATOM_NUMBER != next_atom);

#ifdef DEBUG_IS_LONGEST_PATH_THROUGH_RINGS
  cerr << "next_atom is " << next_atom << endl;
#endif

  if (! m.is_ring_atom(next_atom))    // chain region, keep looking
    return is_longest_path_through_rings(m, next_atom, a2, lp - 1, on_longest_path, ring_system_membership);

  if (! next_is_across_ring_bond)   // biphenyl type linkage
    return is_longest_path_through_rings(m, next_atom, a2, lp - 1, on_longest_path, ring_system_membership);

// We have passed from outside a ring into a ring.  Need to make sure
// our exit is as far away as possible.

  int longest_within_ring_system_distance = 0;
  atom_number_t atom_at_other_end = INVALID_ATOM_NUMBER;
  int atom_at_other_end_on_longest_path = 0;

#ifdef DEBUG_IS_LONGEST_PATH_THROUGH_RINGS
  cerr << "Looking for furthest ring point from " << a1 << " rsys " << ring_system_membership[a1] << endl;
#endif

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (i == a1)
      continue;

#ifdef DEBUG_IS_LONGEST_PATH_THROUGH_RINGS
    cerr << "Checking atom " << i << " rsm " << ring_system_membership[i] << " d = " << m.bonds_between(a1, i) << endl;
#endif

    if (ring_system_membership[i] != ring_system_membership[a1])
      continue;

    int b = m.bonds_between(a1, i);

    if (b < longest_within_ring_system_distance)
      continue;

    if (b > longest_within_ring_system_distance)   // got new longest separation
    {
      longest_within_ring_system_distance = b;
      atom_at_other_end = i;
      atom_at_other_end_on_longest_path = on_longest_path[i];
    }
    else if (0 == atom_at_other_end_on_longest_path)   // same separation, maybe this one on longest path
    {
      atom_at_other_end = i;
      atom_at_other_end_on_longest_path = on_longest_path[i];
    }
  }

#ifdef DEBUG_IS_LONGEST_PATH_THROUGH_RINGS
  cerr << "From atom " << a1 << " atom at furthest end of ring sys is " << atom_at_other_end << " olp " << atom_at_other_end_on_longest_path << endl;
#endif

  if (0 == atom_at_other_end_on_longest_path)
    return 0;

  if (a2 == atom_at_other_end)
    return 1;

// What we have shown here is that if we start at A1, and follow the
// shortest path out, we are on the longest path. But we don't know
// that we have made the maximum possible traversal of the ring system,
// we might have just gone around an ortho substitution. Make sure
// we are the furthest point away from atom_at_other_end

  int longest_separation_in_ring_system = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (ring_system_membership[i] != ring_system_membership[a1])
      continue;

    for (int j = i + 1; j < matoms; j++)
    {
      if (ring_system_membership[j] != ring_system_membership[a1])
        continue;

      if (2 == m.ncon(j))   // no path out of the ring system there
        continue;

      int d = m.bonds_between(i, j);
      if (d > longest_separation_in_ring_system)
        longest_separation_in_ring_system = d;
    }
  }

#ifdef DEBUG_IS_LONGEST_PATH_THROUGH_RINGS
  cerr << " Bonds btw " << a1 << " and " << atom_at_other_end << " is " << m.bonds_between(a1, atom_at_other_end) << " longest sep " << longest_separation_in_ring_system << endl;
#endif

  if (a2 == atom_at_other_end)
    return 1;

  if (m.bonds_between(a1, atom_at_other_end) < longest_separation_in_ring_system)
    return 0;

// Does the atom at the furthest end of the entry have an exit

  const Atom * aoe = m.atomi(atom_at_other_end);

  acon = aoe->ncon();

  int has_branch_outside_ring_system = 0;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = aoe->other(atom_at_other_end, i);

    if (0 == on_longest_path[j])
      continue;

    if (ring_system_membership[a1] == ring_system_membership[j])
      continue;

    has_branch_outside_ring_system = 1;
    break;
  }

  if (! has_branch_outside_ring_system)
    return 0;

  return is_longest_path_through_rings(m, atom_at_other_end, a2, m.bonds_between(atom_at_other_end, a2), on_longest_path, ring_system_membership);
}

static int
compute_ring_systems(Molecule & m,
                     int * tmp)
{
  const int nr = m.nrings();

  if (0 == nr)
    return 0;

  if (1 == nr)
  {
    const Ring * r0 = m.ringi(0);
    r0->set_vector(tmp, 1);
    return 1;
  }

  return m.label_atoms_by_ring_system_including_spiro_fused(tmp);
}

/*
  Here's where you do whatever you want to do with the molecule
  In this case, we count the number of nitrogen atoms
*/

static int
long_molecules(Molecule & m,
               Molecule_Output_Object & passing_molecules,
               Molecule_Output_Object & failing_molecules)
{
  const int matoms = m.natoms();

  if (matoms < everything_smaller_passes)
  {
    molecules_passed_for_too_few_atoms++;
    return do_pass(m, passing_molecules);
  }

  (void) m.ring_membership();

  int lp = m.longest_path();

  longest_path[lp]++;

  if (lp < minimum_longest_path)
  {
    molecules_shorter_than_min_longest_path++;
    return do_pass(m, passing_molecules);
  }

  int * ring_system_identifier = new_int(matoms); std::unique_ptr<int[]> free_rsi(ring_system_identifier);

  int nrs = compute_ring_systems(m, ring_system_identifier);

  if (min_ring_systems > 0 && nrs < min_ring_systems)
  {
    molecules_with_too_few_ring_systems++;
    return do_pass(m, passing_molecules);
  }

  if (pass_if_fused_rings_present && contains_fused_rings(m, ring_system_identifier))
  {
    molecules_passed_for_having_fused_ring_system++;
    return do_pass(m, passing_molecules);
  }

//iw_write_array(ring_system_identifier, matoms, "RSID", cerr);

// Compute the longest path(s)

  Set_of_Atoms equivalent_a1, equivalent_a2;

  int longest_path = 1;
  int pairs_with_longest_path = 0;

  for (int i = 0; i < matoms; i++)
  {
    for (int j = i + 1; j < matoms; j++)
    {
      int d = m.bonds_between(i, j);

      if (d < longest_path)
        continue;

      if (d > longest_path)
      {
        equivalent_a1.resize_keep_storage(0);
        equivalent_a2.resize_keep_storage(0);
        equivalent_a1.add(i);
        equivalent_a2.add(j);
        longest_path = d;
      }
      else
      {
        equivalent_a1.add(i);
        equivalent_a2.add(j);
      }
    }
  }

#ifdef DEBUG_LONG_MOLECULES
  cerr << pairs_with_longest_path << " pairs of atoms on longest paths, length " << longest_path << "\n";
  cerr << "First pair " << a1 << " to " << a2 << endl;
#endif

  int n = equivalent_a1.number_elements();

  if (n > 1)
  {
    if (! ends_of_longest_paths_joined(m, equivalent_a1, equivalent_a2))
      return do_pass(m, passing_molecules);
  }

  int * on_longest_path = new_int(matoms); std::unique_ptr<int[]> free_on_longest_path(on_longest_path);
  int * tmp = new int[matoms]; std::unique_ptr<int[]> free_tmp(tmp);

  int found_longest_path_through_rings = 0;

  for (int i = 0; i < n; i++)
  {
    atom_number_t a1 = equivalent_a1[i];
    atom_number_t a2 = equivalent_a2[i];

    std::fill_n(tmp, matoms, 0);

    tmp[a2] = 1;
    identify_all_atoms_on_longest_paths(m, a1, a2, lp, tmp);

    if (must_take_longest_path_through_rings)
    {
      if (is_longest_path_through_rings (m, a1, a2, lp, tmp, ring_system_identifier))
        found_longest_path_through_rings = 1;
    }

    for (int j = 0; j < matoms; j++)
    {
      if (tmp[j])
        on_longest_path[j] = 1;
    }
  }

  if (must_take_longest_path_through_rings && ! found_longest_path_through_rings)
    return do_pass(m, passing_molecules);

// Now see how far away atoms that are not on the longest path are from the longest path

  int * shortest_dist_to_longest_path = new_int (matoms, matoms + 1); std::unique_ptr<int[]> free_sdlp(shortest_dist_to_longest_path);

  for (int i = 0; i < matoms; i++)   // loops over atoms off the longest path
  {
    if (on_longest_path[i])
      continue;

    for (int j = 0; j < matoms; j++)
    {
      if (! on_longest_path[j])
        continue;

      int d = m.bonds_between(i, j);

      if (d < shortest_dist_to_longest_path[i])
        shortest_dist_to_longest_path[i] = d;
    }
  }

  Accumulator_Int<int> acc;

  for (int i = 0; i < matoms; i++)
  {
    if (on_longest_path[i])
      continue;

#ifdef DEBUG_LONG_MOLECULES
    cerr << "Atom " << i << " off longest path, dist " << shortest_dist_to_longest_path[i] << endl;
#endif

    acc.extra(shortest_dist_to_longest_path[i]);
  }

  if (0 == acc.n())   // whole molecule is longest path
  {
    if (verbose > 1)
      cerr << m.name() << " entire molecule is longest path\n";

    return do_fail(m, equivalent_a1, equivalent_a2, lp, failing_molecules);
  }

#ifdef DEBUG_LONG_MOLECULES
  cerr << "Compare " << acc.maxval() << " with " << max_deviation_from_longest_path << endl;
#endif

  if (acc.maxval() <= max_deviation_from_longest_path)
  {
    if (verbose > 1)
      cerr << m.name() << " rejected for max deviation " << acc.maxval() << endl;

    return do_fail(m, equivalent_a1, equivalent_a2, lp, failing_molecules);
  }

  if (acc.average() > max_average_deviation_from_longest_path)
  {
    if (verbose > 1)
      cerr << m.name() << " rejected for ave deviation " << static_cast<float>(acc.average()) << endl;
    return do_fail(m, equivalent_a1, equivalent_a2, lp, failing_molecules);
  }

  return do_pass(m, passing_molecules);
}

static int
long_molecules (data_source_and_type<Molecule> & input,
                Molecule_Output_Object & passing_molecules,
                Molecule_Output_Object & failing_molecules)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! long_molecules(*m, passing_molecules, failing_molecules))
      return 0;
  }

  return 1;
}

static int
long_molecules (const char * fname, FileType input_type, 
                Molecule_Output_Object & passing_molecules,
                Molecule_Output_Object & failing_molecules)
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

  return long_molecules(input, passing_molecules, failing_molecules);
}

static int
handle_file_opening (Command_Line & cl,
                     const char flag,
                     Molecule_Output_Object & output,
                     const char * passing_or_failing)
{
  const_IWSubstring s = cl.string_value(flag);

  if (cl.option_present('o'))
  {
    if (! output.determine_output_types(cl, 'o'))
    {
      cerr << "Cannot determine output types (-o)\n";
      return 0;
    }
  }
  else
    output.add_output_type(FILE_TYPE_SMI);

  if (output.would_overwrite_input_files(cl, s))
  {
    cerr << "Cannot overwrite input file(s) '" << s << "'\n";
    return 0;
  }

  if (! output.new_stem(s))
  {
    cerr << "Cannot initialise " << passing_or_failing << " stream '" << s << "'\n";
    return 0;
  }

  if (verbose)
    cerr << passing_or_failing << " molecules written to file stem '" << s << "'\n";

  return 1;
}

static int
long_molecules (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:lm:P:F:xD:d:w:Ifr:");

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

  if (cl.option_present('m'))
  {
    if (! cl.value('m', minimum_longest_path) || minimum_longest_path < 1)
    {
      cerr << "The minimum longest path option (-m) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will pass all molecules with longest path shorter than " << minimum_longest_path << endl;
  }

  if (cl.option_present('x'))
  {
    must_take_longest_path_through_rings = 0;

    if (verbose)
      cerr << "Will allow non maximum path traversals through rings\n";
  }

  if (cl.option_present('D'))
  {
    if (! cl.value('D', max_deviation_from_longest_path) || max_deviation_from_longest_path < 0)
    {
      cerr << "The max deviation from longest path option (-D) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will fail molecules with atoms " << max_deviation_from_longest_path << " bonds from the longest path\n";
  }

  if (cl.option_present('d'))
  {
    if (! cl.value('d', max_average_deviation_from_longest_path) || max_average_deviation_from_longest_path < 0.0)
    {
      cerr << "The max average deviation from longest path (-d) option must be a valid float\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will fail molecules if average deviation from longest path greater than " << max_average_deviation_from_longest_path << endl;
  }

  if (cl.option_present('w'))
  {
    if (! cl.value('w', everything_smaller_passes) || everything_smaller_passes < 2)
    {
      cerr << "The everything smaller than option (-w) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "All molecules with fewer than " << everything_smaller_passes << " will pass\n";
  }

  if (cl.option_present('f'))
  {
    pass_if_fused_rings_present = 1;
    if (verbose)
      cerr << "All molecules with fused rings automatically pass\n";
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', min_ring_systems) || min_ring_systems < 1)
    {
      cerr << "The min ring systems needed (-r) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "To fail, a molecule must have " << min_ring_systems << " or more ring systems\n";
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

  Molecule_Output_Object stream_for_molecules_passing, stream_for_molecules_failing;

  if (cl.option_present('P'))
  {
    if (! handle_file_opening(cl, 'P', stream_for_molecules_passing, "passing"))
      return 3;
  }

  if (cl.option_present('F'))
  {
    if (! handle_file_opening(cl, 'F', stream_for_molecules_failing, "failing"))
      return 3;

    if (cl.option_present('I'))
    {
      isotopically_label_failing_molecules = 1;

      if (verbose)
        cerr << "Will isotopically label failing molecules\n";
    }
  }

// Make sure there is something useful done

  if (! cl.option_present('P') && ! cl.option_present('F') && 0 == verbose)
    verbose = 1;

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! long_molecules(cl[i], input_type, stream_for_molecules_passing, stream_for_molecules_failing))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    for (int i = 0; i < longest_path.number_elements(); i++)
    {
      if (longest_path[i])
        cerr << longest_path[i] << " molecules had a longest path of " << i << endl;
    }

    cerr << molecules_passing << " molecules passed, " << molecules_failing << " failed\n";

    if (molecules_passed_for_too_few_atoms)
      cerr << molecules_passed_for_too_few_atoms << " molecules passed for having fewer than " << everything_smaller_passes << " atoms\n";
    if (molecules_passed_for_having_fused_ring_system)
      cerr << molecules_passed_for_having_fused_ring_system << " passed for having fused ring systems\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = long_molecules (argc, argv);

  return rc;
}
