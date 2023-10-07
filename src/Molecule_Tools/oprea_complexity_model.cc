/*
  Implementation of
  Allu and Oprea, J Chem Inf Model 2005, 45 1237-1243
*/

#include <stdlib.h>
#include <iostream>
#include <memory>

#include "cmdline.h"
#include "accumulator.h"
#include "misc.h"

#include "molecule.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "path.h"
#include "target.h"
#include "qry_wstats.h"
#include "istream_and_type.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

static Accumulator<double> complexity_stats;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "$Id: skeleton.cc,v 1.1.1.1 2005/05/10 17:38:33 rx87851 Exp $\n";
  cerr << "Computes molecular complexity metric of Allu and Oprea\n";
  cerr << " J Chem Inf Model, 2005 45 1237-1243\n";
  cerr << "  -q <qry>      query specification for non-complex moities\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit (rc);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment ();

  if (chemical_standardisation.active ())
    chemical_standardisation.process (m);

  return;
}

static double atom_score[HIGHEST_ATOMIC_NUMBER + 1] = { 0.0,
        0.0,      // H
        0.0,      // He
        0.0,
        0.0,
        0.851,    // B
        1.0,      // C
        1.149,    // N
        1.297,    // O
        1.446,    // F
        0.0,      // Ne
        0.0,
        0.0,
        0.0,
        0.0,
        1.086,    // P
        1.235,    // S
        1.384,    // Cl
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.244,    // Br
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.103,     // I
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0};

static double
compute_per_atom_score(const Molecule & m)
{
  int matoms = m.natoms();

  double rc = 0.0;

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t zi = m.atomic_number(i);

//  cerr << "atomic number " << zi << " score " << atom_score[zi] << endl;

    rc += atom_score[zi];
  }

  return rc;
}

static double single_bond_parameter[(HIGHEST_ATOMIC_NUMBER + 1) * (HIGHEST_ATOMIC_NUMBER + 1)];
static double double_bond_parameter[(HIGHEST_ATOMIC_NUMBER + 1) * (HIGHEST_ATOMIC_NUMBER + 1)];
static double triple_bond_parameter[(HIGHEST_ATOMIC_NUMBER + 1) * (HIGHEST_ATOMIC_NUMBER + 1)];
static double aromatic_bond_parameter[(HIGHEST_ATOMIC_NUMBER + 1) * (HIGHEST_ATOMIC_NUMBER + 1)];

static int
initialise_bond_parameter()
{
  set_vector(single_bond_parameter, (HIGHEST_ATOMIC_NUMBER + 1) * (HIGHEST_ATOMIC_NUMBER + 1), 0.0);
  set_vector(double_bond_parameter, (HIGHEST_ATOMIC_NUMBER + 1) * (HIGHEST_ATOMIC_NUMBER + 1), 0.0);
  set_vector(triple_bond_parameter, (HIGHEST_ATOMIC_NUMBER + 1) * (HIGHEST_ATOMIC_NUMBER + 1), 0.0);
  set_vector(aromatic_bond_parameter, (HIGHEST_ATOMIC_NUMBER + 1) * (HIGHEST_ATOMIC_NUMBER + 1), 0.0);

  single_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 6] = 1.0;
  single_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 7] = 0.857;
  single_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 8] = 0.750;
  single_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 9] = 0.667;
  single_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 15] = 0.400;
  single_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 16] = 0.375;
  single_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 17] = 0.353;
  single_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 35] = 0.171;
  single_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 53] = 0.113;

  single_bond_parameter[7 * (HIGHEST_ATOMIC_NUMBER + 1) + 7] = 0.735;
  single_bond_parameter[7 * (HIGHEST_ATOMIC_NUMBER + 1) + 8] = 0.643;

  single_bond_parameter[8 * (HIGHEST_ATOMIC_NUMBER + 1) + 16] = 0.281;

  double_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 6] = 0.500;
  double_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 7] = 0.429;
  double_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 8] = 0.375;
  double_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 16] = 0.188;

  double_bond_parameter[7 * (HIGHEST_ATOMIC_NUMBER + 1) + 7] = 0.367;
  double_bond_parameter[7 * (HIGHEST_ATOMIC_NUMBER + 1) + 8] = 0.321;

  double_bond_parameter[8 * (HIGHEST_ATOMIC_NUMBER + 1) + 16] = 0.141;

  triple_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 6] = 0.333;
  triple_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 7] = 0.286;

  aromatic_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 6] = 0.667;
  aromatic_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 7] = 0.571;
  aromatic_bond_parameter[6 * (HIGHEST_ATOMIC_NUMBER + 1) + 16] = 0.250;

  aromatic_bond_parameter[7 * (HIGHEST_ATOMIC_NUMBER + 1) + 7] = 0.490;
  aromatic_bond_parameter[7 * (HIGHEST_ATOMIC_NUMBER + 1) + 8] = 0.423;

  return 1;
}

static double 
compute_per_bond_score(const Molecule & m)
{
  int nb = m.nedges();

  double rc = 0.0;

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    atomic_number_t z1 = m.atomic_number(a1);
    atomic_number_t z2 = m.atomic_number(a2);

    if (z1 > z2)
      std::swap(z1, z2);

//#define DEBUG_COMPUTE_PER_BOND_SCORE
#ifdef DEBUG_COMPUTE_PER_BOND_SCORE
    if (b->is_aromatic())
      cerr << ' ' << a1 << " (" << m.atomic_symbol(a1) << ") to " << a2 << " (" << m.atomic_symbol(a2) << ") value " << aromatic_bond_parameter[z1 * (HIGHEST_ATOMIC_NUMBER + 1) + z2] << endl;
    else if (b->is_single_bond())
      cerr << ' ' << a1 << " (" << m.atomic_symbol(a1) << ") to " << a2 << " (" << m.atomic_symbol(a2) << ") value " << single_bond_parameter[z1 * (HIGHEST_ATOMIC_NUMBER + 1) + z2] << endl;
    else if (b->is_double_bond())
      cerr << ' ' << a1 << " (" << m.atomic_symbol(a1) << ") to " << a2 << " (" << m.atomic_symbol(a2) << ") value " << double_bond_parameter[z1 * (HIGHEST_ATOMIC_NUMBER + 1) + z2] << endl;
    else if (b->is_triple_bond())
      cerr << ' ' << a1 << " (" << m.atomic_symbol(a1) << ") to " << a2 << " (" << m.atomic_symbol(a2) << ") value " << triple_bond_parameter[z1 * (HIGHEST_ATOMIC_NUMBER + 1) + z2] << endl;
#endif

    if (b->is_aromatic())
      rc += aromatic_bond_parameter[z1 * (HIGHEST_ATOMIC_NUMBER + 1) + z2];
    else if (b->is_single_bond())
      rc += single_bond_parameter[z1 * (HIGHEST_ATOMIC_NUMBER + 1) + z2];
    else if (b->is_double_bond())
      rc += double_bond_parameter[z1 * (HIGHEST_ATOMIC_NUMBER + 1) + z2];
    else if (b->is_triple_bond())
      rc += triple_bond_parameter[z1 * (HIGHEST_ATOMIC_NUMBER + 1) + z2];
  }

  return rc;
}

static double
compute_chiral_centre_score (const Molecule & m)
{
  int n = m.chiral_centres();

  if (0 == n)
    return 0.0;

  double rc = static_cast<double>(2 * n);

  if (1 == n)
    return rc;

// Look for adjacent chiral centres

#ifdef DEBUG_COMPUTE_CHIRAL_CENTRE_SCORE
  cerr << "MOlecule contains " << n << " chiral centres\n";
#endif

  for (int i = 0; i < n; i++)
  {
    const Chiral_Centre * ci = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);

    for (int j = i + 1; j < n; j++)
    {
      const Chiral_Centre * cj = m.chiral_centre_in_molecule_not_indexed_by_atom_number(j);

#ifdef DEBUG_COMPUTE_CHIRAL_CENTRE_SCORE
      cerr << " i = " << i << " j = " << j << " got " << cj << endl;
#endif

      assert (NULL != cj);

      if (ci->involves(cj->a()))
        rc += 1.0;
    }
  }

  return rc;
}

static double
compute_ring_score(Molecule & m)
{
  int nr = m.nrings();

  if (0 == nr)
    return 0.0;

  double rc = 0.0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    int ring_size = ri->size();

    if (5 == ring_size || 6 == ring_size)
      rc += 1.0;
    else if (ring_size >= 10)
      rc += 1.0;
    else
      rc += 2.0;
  }

  return rc;
}

// Lifted from iwdescr

static int
is_spiro_fused (Molecule & m,
                atom_number_t a)
{
  int nr = m.nrings();
  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi (i);

    if (! ri->contains (a))
      continue;

    if (! ri->is_fused())    // not fused, must be spiro
      return 1;

//  RI is fused. Look for the other ring that also contains A

    for (int j = i + 1; j < nr; j++)
    {
      const Ring * rj = m.ringi (j);

      if (! rj->contains (a))
        continue;

      if (! rj->is_fused())     // got it
        return 1;

//    Both rings are fused. If they are in different ring systems, we have it

      if (ri->fused_system_identifier() != rj->fused_system_identifier())
        return 1;

      return 0;       // we have checked 2 rings and no spiro fusion
    }
  }

  return 0;
}

static double
compute_ring_system_score(Molecule & m)
{
  int nr = m.nrings();

  if (0 == nr)
    return 0;

  if (1 == nr)   // must be a single, non-fused ring
    return 1.0;

// Count the number of rings involved in fusions

  double rc = 0.0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (ri->is_fused())
      rc += 1.0;
  }

  if (0.0 == rc)    // no fused rings found
    return static_cast<double>(nr);

  int matoms = m.natoms();

  int * in_system = new_int(matoms); unique_ptr<int> free_in_system(in_system);

  int number_systems = m.label_atoms_by_ring_system_including_spiro_fused(in_system);

//cerr << "Molecule has " << number_systems << " ring systems\n";

  rc += static_cast<double>(number_systems);

// Count spiro carbons

  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = m.atomi(i);

    if (6 != ai->atomic_number())
      continue;

    if (4 != ai->ncon())
      continue;

    if (2 != m.nrings(i))
      continue;

    if (is_spiro_fused(m, i))
      rc += 3.0;
  }

// Now look for bridged atoms

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (1 == ri->largest_number_of_bonds_shared_with_another_ring())
      continue;

    int nbrs = ri->fused_ring_neighbours();

    for (int j = 0; j < nbrs; j++)
    {
      const Ring * rj = ri->fused_neighbour(j);

      if (rj->ring_number() < ri->ring_number())
        continue;

      int shared_bonds = ri->compute_bonds_shared_with(*rj);

      if (shared_bonds > 1)
        rc += static_cast<double>(4 * (shared_bonds - 1));
    }
  }

  return rc;
}

static double
compute_subtraction_for_presence_of_features_no_overlapping_queries(Molecule & m)
{
  Molecule_to_Match target(&m);

  int * already_matched = new_int(m.natoms()); unique_ptr<int> free_already_matched(already_matched);

  int rc = 0;

  int nq = queries.number_elements();

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    if (verbose > 2)
      cerr << " query '" << queries[i]->comment() << "' " << nhits << " hits to '" << m.name() << "'\n";

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      if (e->all_members_set_in_array(already_matched, 1))
        continue;

      e->set_vector(already_matched, 1);

      rc++;
    }
  }

  return static_cast<double>(2 * rc);
}

static double
compute_subtraction_for_presence_of_features(Molecule & m)
{
  Molecule_to_Match target(&m);

  int rc = 0;

  int nq = queries.number_elements();

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    if (verbose > 2)
      cerr << " query '" << queries[i]->comment() << "' " << nhits << " hits to '" << m.name() << "'\n";

    rc += nhits;
  }

  return static_cast<double>(2 * rc);
}

/*
*/

static int
oprea_complexity_model (Molecule & m,
              ostream & output)
{
  m.compute_aromaticity_if_needed();

  double per_atom_score = compute_per_atom_score(m);

  double per_bond_score = compute_per_bond_score(m);

  double chiral_centre_score = compute_chiral_centre_score(m);

  double ring_score = compute_ring_score(m);

  double ring_system_score = compute_ring_system_score(m);

  double subtraction_for_presence_of_features = compute_subtraction_for_presence_of_features(m);

  double complexity = per_atom_score + per_bond_score + chiral_centre_score +
                      ring_score + ring_system_score - subtraction_for_presence_of_features;

  if (verbose > 2)
    cerr << m.name() << ' ' << per_atom_score << ' ' << per_bond_score << ' ' << chiral_centre_score << ' ' << ring_score << ' ' << ring_system_score << " -" << subtraction_for_presence_of_features << " = " << complexity << endl;

  if (complexity < 0.0)
  {
    cerr << "Warning, negative complexity " << complexity << " for '" << m.name() << "'\n";
    cerr << "Subtraction " << subtraction_for_presence_of_features << '\n';
  }

  output << m.name () << " " << complexity << '\n';

  if (verbose)
    complexity_stats.extra(complexity);

  return output.good ();
}

static int
oprea_complexity_model (data_source_and_type<Molecule> & input,
                     ostream & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule ()))
  {
    molecules_read++;

    unique_ptr<Molecule> free_m (m);

    preprocess (*m);

    if (! oprea_complexity_model (*m, output))
      return 0;
  }

  return output.good ();
}

static int
oprea_complexity_model (const char * fname, int input_type, ostream & output)
{
  assert (NULL != fname);

  if (0 == input_type)
  {
    input_type = discern_file_type_from_name (fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input (input_type, fname);
  if (! input.good ())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose (1);

  return oprea_complexity_model (input, output);
}
static int
oprea_complexity_model (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:lq:");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (cl.option_present ('A'))
  {
    if (! process_standard_aromaticity_options (cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage (5);
    }
  }

  if (cl.option_present ('E'))
  {
    if (! process_elements (cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present ('g'))
  {
    if (! chemical_standardisation.construct_from_command_line (cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage (32);
    }
  }

  if (cl.option_present ('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('q'))
  {
    if (! process_queries(cl, queries, verbose, 'q'))
    {
      cerr << "Cannot process queries (-q)\n";
      usage(4);
    }

    int nq = queries.number_elements();

    if (verbose)
      cerr << "Read " << nq << " non-complex queries\n";

    for (int i = 0; i < nq; i++)
    {
      queries[i]->set_find_unique_embeddings_only(1);
    }
  }

  int input_type = 0;

  if (cl.option_present ('i'))
  {
    if (! process_input_type (cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (! all_files_recognised_by_suffix (cl))
    return 4;

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  initialise_bond_parameter();

  cout << "ID OPMC\n";

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! oprea_complexity_model (cl[i], input_type, cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    if (complexity_stats.n() > 1)
      cerr << "Complexity values between " << complexity_stats.minval() << " and " << complexity_stats.maxval() << " ave " << complexity_stats.average() << endl;

    int nq = queries.number_elements();

    for (int i = 0; i < nq; i++)
    {
      queries[i]->report(cerr, verbose);
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = oprea_complexity_model (argc, argv);

  return rc;
}
