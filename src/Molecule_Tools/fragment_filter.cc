/*
  Implementation of Teddy's rules for fragments
*/

#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/allowed_elements.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rotbond_common.h"
#include "Molecule_Lib/standardise.h"

const char * prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static int molecules_passing = 0;

static int molecules_failing = 0;

static int reject_isotopes = 1;

static int lower_atom_count_cutoff = 0;
static int upper_atom_count_cutoff = 0;

static int rejected_for_too_few_atoms = 0;
static int rejected_for_too_many_atoms = 0;
static int rejected_for_too_few_rings = 0;
static int rejected_for_too_many_rings = 0;
static int rejected_for_too_many_ring_systems = 0;
static int rejected_for_ring_system_too_large = 0;
static int rejected_for_min_heteroatom_fraction = 0;
static int rejected_for_max_heteroatom_fraction = 0;
static int rejected_for_too_many_heteroatoms = 0;
static int rejected_for_no_features = 0;
static int rejected_for_donor_acceptor = 0;
static int rejected_for_non_organic = 0;
static int rejected_for_non_allowed_atom = 0;
static int rejected_for_non_allowed_counterion = 0;
static int rejected_for_bad_valence = 0;
static int rejected_for_isotopes = 0;
static int rejected_for_too_much_spinach = 0;
static int rejected_for_between_ring = 0;
static int rejected_for_too_many_halogens = 0;
static int rejected_for_disubstituted_ring_atoms = 0;
static int rejected_for_too_many_connections_outside_ring_system = 0;
static int rejected_for_no_carbon_with_a_proton = 0;
static int rejected_for_too_many_nitrogens = 0;

static int max_halogens_allowed = 2;

static int must_have_carbon_with_a_proton = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int append_rejection_reason_to_failed_molecule_output = 0;

static int min_rings_needed = 0;

static int max_rings_allowed = 0;

static int max_ring_systems_allowed = 0;

static int max_rings_in_a_ring_system = 0;

static int ring_systems_span_spiro_fusions = 1;

static int max_connections_outside_ring_system = 0;

static float max_nitrogen_fraction = static_cast<float>(0.0);

static double min_heteroatom_fraction = 0.0;
static double max_heteroatom_fraction = 1.0;

static int max_any_heteroatom = 0;
static double max_fraction_any_heteroatom = 0.0;

static int cf3_counts_as_single_atom = 0;

static extending_resizable_array<int> max_atoms;

static int max_donor = 0;

static int max_bonds_between_rings = 5;

static Molecule_Output_Object stream_for_passing_molecules;
static Molecule_Output_Object stream_for_rejected_molecules;

static int append_amw_of_counterions_to_name = 0;

static Allowed_Elements allowed_elements;

static IW_STL_Hash_Map_int dis_allowed_counterions;

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
// clang-format off
  cerr << "  -h <fract>    minimum heteroatom fraction\n";
  cerr << "  -H <fract>    maximum heteroatom fraction\n";
  cerr << "  -X <num/frac> max number/fraction of any heteroatom\n";
  cerr << "  -N <fract>    max fraction of nitrogen atoms\n";
  cerr << "  -M <nr=n>     maximum number of atoms for molecules with <nr> ring\n";
  cerr << "  -R <n>        maximum number of rings allowed\n";
  cerr << "  -r <n>        minimum number of rings required\n";
  cerr << "  -Y <n>        maximum number of ring systems allowed\n";
  cerr << "  -G <n>        maximum number of rings in a ring system\n";
  cerr << "  -D <n>        maximum number of 'donors' (heteroatoms with H)\n";
  cerr << "  -c <n>        lower atom count cutoff - atoms in largest frag only\n";
  cerr << "  -C <n>        upper atom count cutoff - atoms in largest frag only\n";
  cerr << "  -b <n>        maximum number of bonds between ring systems\n";
  cerr << "  -y            break ring systems at spiro fusions\n";
  cerr << "  -F            count CF3 groups as a single atom\n";
  cerr << "  -s <n>        max connections from a ring system\n";
  cerr << "  -p <n>        must have at least <n> carbons with a proton\n";
  cerr << "  -S <stem>     survivors output file name stem\n";
  cerr << "  -B <stem>     failed molecule output file name stem\n";
  cerr << "  -a            append rejection reason to name in rejected file\n";
  cerr << "  -f            append counterion AMW to name\n";
  cerr << "  -J <fname>    list if dis-allowed counterions\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -I keep       keep isotopic molecules\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -o <type>     output specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
// clang-format on

  exit(rc);
}

static void
reset_all_formal_charges (Molecule & m)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    m.set_formal_charge_if_different(i, 0);
  }

  return;
}

static int
read_dis_allowed_counterion_record (const const_IWSubstring & buffer,
                             IW_STL_Hash_Map_int & dis_allowed_counterions)
{
  if (buffer.starts_with('#') || 0 == buffer.length())
    return 1;

  Molecule m;

  if (! m.build_from_smiles(buffer))
  {
    cerr << "Invalid smiles '" << buffer << "'\n";
    return 0;
  }

  reset_all_formal_charges(m);

  const IWString & s = m.unique_smiles();

  dis_allowed_counterions[s] = 0;

//cerr << "Adding '" << s << "' to disallowed\n";

  return 1;
}

static int
read_dis_allowed_counterions(iwstring_data_source & input,
                             IW_STL_Hash_Map_int & dis_allowed_counterions)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (! read_dis_allowed_counterion_record(buffer, dis_allowed_counterions))
    {
      cerr << "Cannot read dis allowed ion '" << buffer << "'\n";
      return 0;
    }
  }

  return dis_allowed_counterions.size();
}

static int
read_dis_allowed_counterions(const char * fname,
                             IW_STL_Hash_Map_int & dis_allowed_counterions)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_dis_allowed_counterions(input, dis_allowed_counterions);
}

//#define DEBUG_COUNT_FEATURES_ALONG_PATH

static int
count_features_along_path(Molecule & m,
                          atom_number_t zatom,
                          atom_number_t destination,
                          int & heteroatoms_found,
                          int & unsaturation_found,
                          int & ch2_found,
                          int & non_rotatable_bonds_found)
{
  const Atom * a = m.atomi(zatom);

  int d = m.bonds_between(zatom, destination);

#ifdef DEBUG_COUNT_FEATURES_ALONG_PATH
  cerr << "count_features_along_path atom " << zatom << " so far " << heteroatoms_found << " htroatom, and " << unsaturation_found << " unsaturation, dest " << destination << ", dist " << d << endl;
#endif

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(zatom);

#ifdef DEBUG_COUNT_FEATURES_ALONG_PATH
    if (j == destination)
      cerr << "Got to destination: heteroatoms_found: " << heteroatoms_found << " unsaturation " << unsaturation_found << endl;
#endif

    if (j == destination)
      ;
    else if (d - 1 != m.bonds_between(j, destination))   // move toward destination
      continue;

    if (b->is_single_bond())
    {
      if (is_non_rotatable_amide(m, zatom, j))
        non_rotatable_bonds_found++;
    }
    else if (b->is_triple_bond())    // don't want inter-ring acetylenes
      unsaturation_found += 99;
    else if (b->is_double_bond())
      non_rotatable_bonds_found++;

    if (j == destination)
    {
      if (unsaturation_found > 2)   // otherwise C1(=CC=CO1)C=NN=CC1=CC=CO1 PBCHM21556 gets through
        return 0;

      if (heteroatoms_found > 0 && unsaturation_found > 0)
        return 1;

      if (1 == ch2_found && (heteroatoms_found > 0 || unsaturation_found > 0))
        return 1;

      if (non_rotatable_bonds_found)
        return 1;

      return 0;
    }

    const Atom * aj = m.atomi(j);

    if (6 != aj->atomic_number())
      heteroatoms_found++;

    if (aj->ncon() < aj->nbonds())
      unsaturation_found++;
    else if (6 == aj->atomic_number())
      ch2_found++;

//  if (heteroatoms_found && unsaturation_found)   leave out, otherwise we don't detect C1(=CC=CO1)C=NN=CC1=CC=CO1 18985
//    return 1;

    return count_features_along_path(m, j, destination, heteroatoms_found, unsaturation_found, ch2_found, non_rotatable_bonds_found);
  }

  cerr << "How did we get here!\n";
  return 0;   // not sure how control ever could come to here
}

/*static int
custom_between_ring_filter (Molecule & m,
                            const Ring & r1,
                            const Ring & r2)
{
  int n1 = r1.number_elements();
  int n2 = r2.number_elements();

  int dmin = m.natoms() + 6;
  atom_number_t min1 = 0;  // Only initialise these to keep the compiler quiet.
  atom_number_t min2 = 0;

  for (int i = 0; i < n1; i++)
  {
    atom_number_t a1 = r1[i];

    for (int j = 0; j < n2; j++)
    {
      atom_number_t a2 = r2[j];

      if (a1 == a2)
        continue;

      int d = m.bonds_between(a1, a2);

      if (d < dmin)
      {
        min1 = a1;
        min2 = a2;
        dmin = d;
      }
    }
  }

  cerr << "Min distance between rings " << dmin << endl;
  if (dmin <= 2)
    return 1;

  if (dmin > 4)
    return 0;

// If 3 bonds, there must be at least a heteroatom or an unsaturation

  int heteroatoms_found = 0;
  int unsaturation_found = 0;

  return count_features_along_path(m, min1, min2, heteroatoms_found, unsaturation_found);
}*/

//#define DEBUG_CUSTOM_BETWEEN_RING_FILTER

static int
custom_between_ring_filter (Molecule & m,
                            const int * ring_sys_id,
                            int nrs)
{
  int nr = m.nrings();

  if (nr < 2)
    return 1;

  int matoms = m.natoms();

#ifdef DEBUG_CUSTOM_BETWEEN_RING_FILTER
  for (int i = 0; i < matoms; i++)
  {
    cerr << " atom " << i << " rsid " << ring_sys_id[i] << endl;
  }

  cerr << "Checking " << nrs << " ring systems\n";
#endif

  for (int rsid1 = 1; rsid1 <= nrs; rsid1++)
  {
    for (int rsid2 = rsid1 + 1; rsid2 <= nrs; rsid2++)
    {
      atom_number_t min1 = 0;  // Initialised only to keep the compiler quiet.
      atom_number_t min2 = 0;
      int dmin = matoms;

      for (int i = 0; i < matoms; i++)
      {
        if (rsid1 != ring_sys_id[i])
          continue;

        for (int j = 0; j < matoms; j++)
        {
          if (rsid2 != ring_sys_id[j])
            continue;

          int d = m.bonds_between(i, j);
          if (d < dmin)
          {
            min1 = i;
            min2 = j;
            dmin = d;
          }
        }
      }

#ifdef DEBUG_CUSTOM_BETWEEN_RING_FILTER
      cerr << "Between systems dmin " << dmin << endl;
#endif

      if (dmin > max_bonds_between_rings)
        return 0;

      if (dmin <= 3)   // even two naked methyl groups are OK
        continue;

      int heteroatoms_found = 0;
      int unsaturation_found = 0;
      int ch2_found = 0;
      int non_rotatable_bonds_found = 0;   // mostly amides
      if (! count_features_along_path(m, min1, min2, heteroatoms_found, unsaturation_found, ch2_found, non_rotatable_bonds_found))
        return 0;
    }
  }

  return 1;
}

static int
various_ring_things(Molecule & m,
                    const int * ring_system_identifier,
                    int ring_systems_present,
                    IWString & rejection_reason)
{
  int matoms = m.natoms();

  int * spinach = new int[matoms]; std::unique_ptr<int[]> free_spinach(spinach);

  m.identify_spinach(spinach);

//cerr << "Processing " << m.smiles() << ' ' << m.name() << ", " << ring_systems_present << " ring_systems_present\n";

  for (int i = 1; i <= ring_systems_present; i++)
  {
    int atoms_in_ring_system = 0;

    int outside_ring_system = 0;

    for (int j = 0; j < matoms; j++)
    {
      if (i != ring_system_identifier[j])
        continue;

      atoms_in_ring_system++;

      const Atom * aj = m.atomi(j);

      int jcon = aj->ncon();

      if (2 == jcon)   // definitely nothing outside
        continue;

      int outside_ring_this_atom = 0;

      for (int k = 0; k < jcon; k++)
      {
        const Bond * b = aj->item(k);

        atom_number_t l = b->other(j);

        if (! spinach[l])    // either between rings or in the same ring
          continue;

        if (b->is_double_bond() && 1 == m.ncon(l))   // part of ring system
          continue;

        outside_ring_this_atom++;
      }

      if (0 == outside_ring_this_atom)
        ;
      else if (1 == outside_ring_this_atom)
        outside_ring_system++;
      else if (4 == aj->ncon() && 1 == m.nrings(j) && m.in_ring_of_given_size(j, 3))   // disubstituted on a 3ring is OK
        ;
      else
      {
//      cerr << "Atom " << j << " has " << outside_ring_this_atom << " connections\n";
        rejected_for_disubstituted_ring_atoms++;
        rejection_reason = " disubstituted";
        return 0;
      }
    }

//  cerr << atoms_in_ring_system << " atoms in ring system " << i << " spinach connections " << outside_ring_system << endl;

    if (outside_ring_system > max_connections_outside_ring_system)
    {
      rejected_for_too_many_connections_outside_ring_system++;
      rejection_reason = "ring_substitution";
      return 0;
    }
  }

  return 1;
}

/*static int
halogen_count (const Molecule & m)
{
  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = m.atomic_number(i);

    if (z < 17)
      ;
    else if (17 == z)
      rc++;
    else if (35 == z)
      rc++;
    else if (53 == z)
      rc++;
  }

  return rc;
}*/

/*
  If the molecule has 3 rings, only allow single atom spinach groups
*/

static int
too_much_spinach (Molecule & m)
{
  int matoms = m.natoms();

  int * spinach = new int[matoms]; std::unique_ptr<int[]> free_spinach(spinach);

  int atoms_in_spinach = m.identify_spinach(spinach);

  if (atoms_in_spinach < 2)
    return 0;

  for (int i = 0; i < matoms; i++)
  {
    if (! spinach[i])
      continue;

    if (1 != m.ncon(i))
      return 1;
  }

  return 0;
}

static void
preprocess (Molecule & m)
{
  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (reduce_to_largest_fragment)
  {
    molecular_weight_t initial_amw = static_cast<molecular_weight_t>(0.0);

    if (append_amw_of_counterions_to_name)
      initial_amw = m.molecular_weight_ignore_isotopes();

    m.reduce_to_largest_fragment_carefully();

    if (append_amw_of_counterions_to_name)
    {
      IWString amw_of_counterion;
      amw_of_counterion << " AMWC=" << (initial_amw - m.molecular_weight_ignore_isotopes());
      m.append_to_name(amw_of_counterion);
    }
  }

  return;
}

static int
is_tetrazole (Molecule & m,
              const Ring & r,
              int * is_nitrogen)
{
  int nitrogen_count = 0;
  int non_nitrogen_count = 0;
  atomic_number_t last_nitrogen = INVALID_ATOM_NUMBER;

  for (int i = 0; i < 5; i++)
  {
    atom_number_t j = r[i];

    if (7 == m.atomic_number(j))
    {
      nitrogen_count++;
      last_nitrogen = j;
    }
    else if (6 != m.atomic_number(j))
    {
      non_nitrogen_count++;
      if (non_nitrogen_count > 1)
        return 0;
    }
  }

  if (4 != nitrogen_count)
    return 0;

  r.set_vector(is_nitrogen, 0);

  is_nitrogen[last_nitrogen] = 1;

  return 1;
}

/*
  Straightforward, except that tetrazoles count for just 1 nitrogen
*/

static float
nitrogen_ratio (Molecule & m)
{
  int matoms = m.natoms();

  int nr = m.nrings();

//float rc = 0.0;

  int * is_nitrogen = new_int(matoms); std::unique_ptr<int[]> free_is_nitrogen(is_nitrogen);

  int n = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (7 == m.atomic_number(i))
      n++;
  }

  if (n < 4 || 0 == nr)    // no tetrazoles here
    return static_cast<float>(n) / static_cast<float>(matoms);

// Maybe we have tetrazoles

  m.compute_aromaticity_if_needed();

  int tetrazoles_found = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (5 != ri->number_elements())
      continue;

    if (ri->is_fused())
      continue;

    if (is_tetrazole(m, *ri, is_nitrogen))
      tetrazoles_found++;
  }

//cerr << m.name() << " contains " << tetrazoles_found << " tetrazoles, n = " << n << endl;

  if (tetrazoles_found)
    n = count_non_zero_occurrences_in_array(is_nitrogen, matoms);

  return static_cast<float>(n) / static_cast<float>(matoms);
}

static int
identify_cf3_groups(const Molecule & m,
                    int * is_cf3)
{
  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = m.atomi(i);

    if (6 != ai->atomic_number())
      continue;

    if (4 != ai->ncon())
      continue;

    Set_of_Atoms f, not_f;

    for (int j = 0; j < 4; j++)
    {
      atom_number_t k = ai->other(i, j);

      if (9 == m.atomic_number(k))
        f.add(k);
      else
        not_f.add(k);
    }

    if (3 == f.number_elements())
    {
      f.set_vector(is_cf3, 1);
      rc++;
    }
  }

  return rc;
}

/*
  Identify ring systems. We merge spiro fused rings.
  We also check for any ring systems that are too large
*/

static int
identify_ring_systems (Molecule & m,
                       int * ring_system_identifier,
                       IWString & rejection_reason)
{
  int nr = m.nrings();

  int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  resizable_array<const Ring *> unfused_possibly_spiro_fused_rings;

  int rc = 0;

  int rings_processed = 0;

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])
      continue;

    const Ring * ri = m.ringi(i);

    rc++;

    ri->set_vector(ring_system_identifier, rc);

    rings_processed++;

    int rings_this_system = 1;

//  Include the most common case of rings in the same fused system

    if (ri->is_fused())
    {
      for (int j = i + 1; j < nr; j++)
      {
        if (ring_already_done[j])
          continue;

        const Ring * rj = m.ringi(j);

        if (rj->fused_system_identifier() == ri->fused_system_identifier())
        {
          rj->set_vector(ring_system_identifier, rc);
          ring_already_done[j] = 1;
          rings_processed++;
          rings_this_system++;
        }
      }
    }

//  Now grow this out via any sprio fusions

    while (ring_systems_span_spiro_fusions && rings_processed < nr)
    {
      int continue_looping = 0;    // if we grow the ring system, need to look again

      for (int j = i + 1; j < nr; j++)
      {
        if (ring_already_done[j])
          continue;

        const Ring * rj = m.ringi(j);

        if (! rj->any_members_set_in_array(ring_system_identifier))
          continue;

        rj->set_vector(ring_system_identifier, rc);
        ring_already_done[j] = 1;
        rings_processed++;
        rings_this_system++;
        continue_looping = 1;
      }

      if (! continue_looping)
        break;
    }

//  cerr << rings_this_system << " rings in the ring system\n";
    if (0 == max_rings_in_a_ring_system)
      ;
    else if (rings_this_system > max_rings_in_a_ring_system)
    {
      rejection_reason << "ring system too large " << rings_this_system;
      rejected_for_ring_system_too_large++;
      return 0;
    }
  }

  return rc;
}

/*
  Always append the reason
*/

static int
fragment_filter(Molecule & m,
                IWString & rejection_reason)
{
  if (! m.valence_ok())
  {
    rejection_reason << "bad valence";
    rejected_for_bad_valence++;
    return 0;
  }

  if (! m.organic_only())
  {
    rejection_reason << "non organic";
    rejected_for_non_organic++;
    return 0;
  }

  if (reject_isotopes && m.number_isotopic_atoms())
  {
    rejection_reason << "isotopes";
    rejected_for_isotopes++;
    return 0;
  }

  int matoms = m.natoms();

  if (lower_atom_count_cutoff > 0 && matoms < lower_atom_count_cutoff)
  {
    rejection_reason << "too few atoms " << matoms;
    rejected_for_too_few_atoms++;
    return 0;
  }

  if (upper_atom_count_cutoff > 0 && matoms > upper_atom_count_cutoff)
  {
    rejection_reason << "too many atoms " << matoms;
    rejected_for_too_many_atoms++;
    return 0;
  }

  int nr = m.nrings();

  if (nr < min_rings_needed)
  {
    rejection_reason << "too few rings " << nr;
    rejected_for_too_few_rings++;
    return 0;
  }

  if (nr > max_rings_allowed)
  {
    rejection_reason << "too many rings " << nr;
    rejected_for_too_many_rings++;
    return 0;
  }

  int * is_cf3 = new_int(matoms); std::unique_ptr<int[]> free_is_cf3(is_cf3);

  int halogen = identify_cf3_groups(m, is_cf3);

  if (halogen > 1)
  {
    rejection_reason << "too many CF3";
    rejected_for_too_many_halogens++;
    return 0;
  }

  int atoms_to_check = matoms;

  if (cf3_counts_as_single_atom)
    atoms_to_check -= (halogen * 3);

  if (0 == max_atoms[nr])
    ;
  else if (atoms_to_check > max_atoms[nr])
  {
    rejection_reason << "too many atoms " << matoms;
    rejected_for_too_many_atoms++;
    return 0;
  }

  int carbon = 0;
  int free_nitrogen = 0;
  int carbon_with_a_proton = 0;

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = m.atomic_number(i);
    if (6 == z)
    {
      carbon++;
      if (m.hcount(i))
        carbon_with_a_proton++;
    }
    else if (7 == z && m.ncon(i) == m.nbonds(i))
      free_nitrogen++;
    else if (is_cf3[i])
      ;
    else if (9 == z)
      halogen++;
    else if (z < 17)    // definitely not a halogen
      ;
    else if (17 == z)
      halogen++;
    else if (35 == z)
      halogen++;
    else if (53 == z)
      halogen++;
  }

  if (must_have_carbon_with_a_proton > 0 && carbon_with_a_proton < must_have_carbon_with_a_proton)
  {
    rejection_reason << "Insufficient CH";
    rejected_for_no_carbon_with_a_proton++;
    return 0;
  }

  int heteroatoms = matoms - carbon;

  double htroaf = static_cast<double>(heteroatoms) / static_cast<double>(matoms);

  if (htroaf > max_heteroatom_fraction)
  {
    rejection_reason << "max heteroatom fraction " << static_cast<float>(htroaf);
    rejected_for_max_heteroatom_fraction++;
    return 0;
  }

// Free nitrogens are considered very desirable, so double their contribution

  htroaf = static_cast<double>(heteroatoms + free_nitrogen) / static_cast<double>(matoms);

//cerr << m.smiles() << ' ' << m.name() << ' ' << static_cast<float>(htroaf) << endl;

  if (htroaf < min_heteroatom_fraction)
  {
    rejection_reason << "min heteroatom fraction " << htroaf;
    rejected_for_min_heteroatom_fraction++;
    return 0;
  }

  if (max_nitrogen_fraction > static_cast<float>(0.0))
  {
    float n = nitrogen_ratio(m);

//  cerr << m.name().word(0) << " nitrogen ratio " << n << ", compare " << max_nitrogen_fraction << endl;

    if (n > max_nitrogen_fraction)
    {
      rejection_reason << "too many nitrogens " << n;
      rejected_for_too_many_nitrogens++;
      return 0;
    }
  }

  if (max_halogens_allowed > 0 && halogen > max_halogens_allowed)
  {
    rejection_reason << "too many halogens " << halogen;
    rejected_for_too_many_halogens++;
    return 0;
  }

  if ((max_any_heteroatom > 0 && heteroatoms >= max_any_heteroatom) || 
      (max_fraction_any_heteroatom > 0.0 && htroaf > max_fraction_any_heteroatom))
  {
    extending_resizable_array<int> element_count;
    for (int i = 0; i < matoms; i++)
    {
      atomic_number_t z = m.atomic_number(i);

      if (6 != z)
        element_count[z]++;
    }

    int n = element_count.number_elements();

    if (max_any_heteroatom > 0)
    {
      for (int i = 7; i < n; i++)
      {
        if (element_count[i] > max_any_heteroatom)
        {
          rejection_reason << "too many heteroatoms " << i << " (" << element_count[i] << ')';
          rejected_for_too_many_heteroatoms++;
          return 0;
        }
      }
    }
    else
    {
      for (int i = 7; i < n; i++)
      {
        int c = element_count[i];
        if (0 == c)
          continue;

        if (static_cast<double>(c) / static_cast<double>(matoms) > max_fraction_any_heteroatom)
        {
          rejection_reason << "too many heteroatoms " << i;
          rejected_for_too_many_heteroatoms++;
          return 0;
        }
      }
    }
  }

  int donors = 0;
  int acceptors = 0;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    if (6 == a->atomic_number())
      continue;

    if (m.hcount(i))
      donors++;
    else
      acceptors++;
  }

  if (0 == donors && 0 == acceptors)
  {
    rejection_reason << "no features";
    rejected_for_no_features++;
    return 0;
  }

  if (max_donor > 0 && donors > max_donor)
  {
    rejection_reason << " hbonds";
    rejected_for_donor_acceptor++;
    return 0;
  }

  if (nr)
  {
    int * ring_system_identifier = new_int(matoms); std::unique_ptr<int[]> free_ring_system_identifier(ring_system_identifier);

    int number_ring_systems = identify_ring_systems(m, ring_system_identifier, rejection_reason);

    if (0 == number_ring_systems)
      return 0;

    if (0 == max_ring_systems_allowed)
      ;
    else if (number_ring_systems > max_ring_systems_allowed)
    {
      rejection_reason << "too many ring systems " << number_ring_systems;
      rejected_for_too_many_ring_systems++;
      return 0;
    }

    if (! various_ring_things(m, ring_system_identifier, number_ring_systems, rejection_reason))
      return 0;

//   Ad hoc rule. If there are 3 rings present, don't allow any large spinach
//   groupings

    if (nr > 2)
    {
      if (too_much_spinach(m))
      {
        rejection_reason << " spinach";
        rejected_for_too_much_spinach++;
        return 0;
      }
    }

    if (number_ring_systems > 1 && max_bonds_between_rings > 0)
    {
      if (! custom_between_ring_filter(m, ring_system_identifier, number_ring_systems))
      {
        rejection_reason << " between ring";
        rejected_for_between_ring++;
        return 0;
      }
    }
  }

  return 1;
}

static int
fragment_filter (Molecule & m,
                 std::ostream & output)
{
  IWString rejection_reason;
  rejection_reason << " BF:";

  if (fragment_filter(m, rejection_reason))
  {
    molecules_passing++;
    if (stream_for_passing_molecules.active())
      stream_for_passing_molecules.write(m);
  }
  else
  {
    molecules_failing++;
    if (stream_for_rejected_molecules.active())
    {
      if (append_rejection_reason_to_failed_molecule_output)
        m.append_to_name(rejection_reason);

      stream_for_rejected_molecules.write(m);
    }
  }

  return 1;
}

static int
contains_disallowed_counterions (Molecule & m,
                                 IW_STL_Hash_Map_int & dis_allowed_counterions)
{
  int nf = m.number_fragments();

  if (1 == nf)
    return 0;

  resizable_array_p<Molecule> c;

  m.create_components(c);

  for (int i = 0; i < nf; i++)
  {
    Molecule * ci = c[i];

    reset_all_formal_charges(*ci);

    const IWString & u = ci->unique_smiles();

//  cerr << "Searching for '" << u << "'\n";

    IW_STL_Hash_Map_int::iterator f = dis_allowed_counterions.find(u);

    if (f == dis_allowed_counterions.end())
      continue;

    (*f).second++;

    return 1;
  }

  return 0;           // no disallowed counterions found
}

static int
fragment_filter (data_source_and_type<Molecule> & input,
                     std::ostream & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    if (allowed_elements.contains_non_allowed_atoms(*m))
    {
      rejected_for_non_allowed_atom++;
      if (stream_for_rejected_molecules.active())
        stream_for_rejected_molecules.write(*m);
      continue;
    }

    if (0 == dis_allowed_counterions.size())
      ;
    else if (contains_disallowed_counterions(*m, dis_allowed_counterions))
    {
      rejected_for_non_allowed_counterion++;
      if (stream_for_rejected_molecules.active())
        stream_for_rejected_molecules.write(*m);
      continue;
    }

    assert (reduce_to_largest_fragment);    // too complicated otherwise

    preprocess(*m);

//  cerr << "After preprocessing '" << m->smiles() << "'\n";

    if (! fragment_filter(*m, output))
      return 0;
  }

  return output.good();
}

static int
fragment_filter (const char * fname, FileType input_type, std::ostream & output)
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

  return fragment_filter(input, output);
}
static int
fragment_filter (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:lS:B:o:r:R:h:H:M:D:Y:aG:c:fI:b:s:p:yFJ:N:X:");

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

  if (cl.option_present('b'))
  {
    if (! cl.value('b', max_bonds_between_rings) || max_bonds_between_rings < 2)
    {
      cerr << "The max bonds between ring systems option (-b) must be a whole +ve integer\n";
      usage(4);
    }

    if (verbose)
      cerr << "Max bonds between ring systems " << max_bonds_between_rings << '\n';
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', max_connections_outside_ring_system) || max_connections_outside_ring_system < 1)
    {
      cerr << "The max spinach connections to a ring system option (-s) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "A max of " << max_connections_outside_ring_system << " spinach connections to any ring system allowed\n";
  }

  if (cl.option_present('p'))
  {
    if (! cl.value('p', must_have_carbon_with_a_proton) || must_have_carbon_with_a_proton < 1)
    {
      cerr << "The min number carbons with protons option (-p) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Each molecule must contain at least " << must_have_carbon_with_a_proton << " carbons with one or more protons\n";
  }

  if (cl.option_present('I'))
  {
    const_IWSubstring s;
    int i = 0;
    while (cl.value('I', s, i++))
    {
      if ("keep" == s)
      {
        reject_isotopes = 0;
        if (verbose)
          cerr << "Will retain molecules with isotopically labelled atoms\n";
      }
      else
      {
        cerr << "Unrecognised -I qualifier '" << s << "'\n";
        return 4;
      }
    }
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', lower_atom_count_cutoff) || lower_atom_count_cutoff < 1)
    {
      cerr << "The lower atom count cutoff must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard molecules with fewer than " << lower_atom_count_cutoff << " atoms in the largest fragment\n";
  }

  if (cl.option_present('D'))
  {
    if (! cl.value('D', max_donor) || max_donor < 0)
    {
      cerr << "The max donor count (-D) must be a whole non-negative number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Molecules more than " << max_donor << " heteroatoms with H's will be rejected\n";
  }

  if (cl.option_present('G'))
  {
    if (! cl.option_present('R'))
    {
      cerr << "Must specify max allowable rings (-R) when using the -G option\n";
      usage(5);
    }

    if (! cl.value('G', max_rings_in_a_ring_system) || max_rings_in_a_ring_system < 1)
    {
      cerr << "The maximum number of rings in a ring system (-G) option must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will reject molecules that contain a ring system with more than " << max_rings_in_a_ring_system << " rings\n";
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', min_rings_needed) || min_rings_needed < 0)
    {
      cerr << "The mimimum number of rings needed (-r) must be a whole non-negative number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Molecules must have a minimum of " << min_rings_needed << " rings\n";
  }

  if (cl.option_present('R'))
  {
    if (! cl.value('R', max_rings_allowed) || max_rings_allowed < min_rings_needed)
    {
      cerr << "The maximum number of rings needed (-R) must be a whole number >= " << min_rings_needed << endl;
      usage(2);
    }

    if (verbose)
      cerr << "Molecules must have no more than " << max_rings_allowed << " rings\n";
  }

  if (cl.option_present('Y'))
  {
    if (! cl.option_present('R'))
    {
      cerr << "Must specify max allowable rings (-R) when using the -Y option\n";
      usage(5);
    }

    if (! cl.value('Y', max_ring_systems_allowed) || max_ring_systems_allowed < 0)
    {
      cerr << "The maximum number of ring systems allowed (-Y) must be a whole number >= " << min_rings_needed << endl;
      usage(2);
    }

    if (verbose)
      cerr << "Molecules must have no more than " << max_ring_systems_allowed << " ring systems\n";
  }

  if (cl.option_present('y'))
  {
    ring_systems_span_spiro_fusions = 0;

    if (verbose)
      cerr << "Ring systems will not span spiro fused systems\n";
  }

  if (cl.option_present('h'))
  {
    if (! cl.value('h', min_heteroatom_fraction) || min_heteroatom_fraction < 0.0 || max_heteroatom_fraction > 1.0)
    {
      cerr << "The minimum heteroatom fraction (-h) must be a valid fraction\n";
      usage(5);
    }

    if (verbose)
      cerr << "Molecules with heteroatom fractions less than " << min_heteroatom_fraction << " will be rejected\n";
  }

  if (cl.option_present('H'))
  {
    if (! cl.value('H', max_heteroatom_fraction) || max_heteroatom_fraction < min_heteroatom_fraction || max_heteroatom_fraction > 1.0)
    {
      cerr << "The maximum heteroatom fraction (-H) must be a valid fraction\n";
      usage(5);
    }

    if (verbose)
      cerr << "Molecules with heteroatom fractions above " << max_heteroatom_fraction << " will be rejected\n";
  }

  if (cl.option_present('X'))
  {
    double x;
    cl.value('X', x);

    if (x <= 0.0)
    {
      cerr << "The max heteroatom fraction/count option (-X) must be a positive number\n";
      usage(2);
    }

    double tmp;

    if (x < 1.0f)
      max_fraction_any_heteroatom = x;
    else if (cl.value('X', tmp))
    {
      max_any_heteroatom = static_cast<int>(tmp);   // should check against anyone entering 3.14....
    }
    else
    {
      cerr << "The max number or fraction any heteratom must be either a valid fraction or a number\n";
      usage(2);
    }
  }

  if (cl.option_present('N'))
  {
    if (! cl.value('N', max_nitrogen_fraction) || max_nitrogen_fraction <= 0.0 || max_nitrogen_fraction > 1.0)
    {
      cerr << "The max nitrogen atom count fraction (-N) must be a valid fraction\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will reject molecules with more than " << max_nitrogen_fraction << " of the atoms nitrogens\n";
  }

  if (cl.option_present('M'))
  {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++))
    {
      const_IWSubstring cnr, cna;
      if (! m.split(cnr, '=', cna))
      {
        cerr << "Max atom directives (-M) must be of the form 'nrings=natoms'\n";
        usage(5);
      }

      int nr;
      if (! cnr.numeric_value(nr) || nr < 0)
      {
        cerr << "Invalid nrings specification '" << cnr << "'\n";
        usage(5);
      }

      int na;
      if (! cna.numeric_value(na) || na < 1)
      {
        cerr << "Invalid natoms specification '" << cna << "'\n";
        usage(5);
      }

      max_atoms[nr] = na;

      if (verbose)
        cerr << "Molecules with " << nr << " rings can have at most " << na << " atoms\n";

      if (na > upper_atom_count_cutoff)
        upper_atom_count_cutoff = na;
    }
  }

  if (cl.option_present('F'))
  {
    cf3_counts_as_single_atom = 1;

    if (verbose)
      cerr << "CF3 groups count as single atom\n";
  }

  if (cl.option_present('a'))
  {
    append_rejection_reason_to_failed_molecule_output = 1;

    if (verbose)
      cerr << "Will append the rejection reason to failed molecule's name\n";
  }

  if (cl.option_present('f'))
  {
    append_amw_of_counterions_to_name = 1;
    if (verbose)
      cerr << "Will append the AMW of the counterion(s) to the molecule name\n";
  }

  if (cl.option_present('J'))
  {
    const char * j = cl.option_value('J');

    int esave = auto_create_new_elements();
    if (! read_dis_allowed_counterions (j, dis_allowed_counterions))
    {
      cerr << "Cannot read disallowed counterions from '" << j << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Read " << dis_allowed_counterions.size() << " disallowed counterions from '" << j << "'\n";

    set_auto_create_new_elements(esave);
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
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

/*if (cl.option_present('X'))
  {
    if (! cl.option_present('o'))
      stream_for_rejected_molecules.add_output_type(FILE_TYPE_SMI);
    else if (! stream_for_rejected_molecules.determine_output_types(cl))
    {
      cerr << "Cannot determine output type(s) for rejection file ('X)\n";
      return 5;
    }

    const_IWSubstring x = cl.string_value('X');

    if (stream_for_rejected_molecules.would_overwrite_input_files(cl, x))
    {
      cerr << "Cannot overwrite input file(s) '" << x << "'\n";
      return 7;
    }

    if (! stream_for_rejected_molecules.new_stem(x))
    {
      cerr << "Cannot open rejection stream '" << x << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Rejected molecules written to '" << x << "'\n";
  }*/

  if (cl.option_present('S'))
  {
    if (! cl.option_present('o'))
      stream_for_passing_molecules.add_output_type(FILE_TYPE_SMI);
    else if (! stream_for_passing_molecules.determine_output_types(cl))
    {
      cerr << "Cannot determine output type(s) for passing file ('S)\n";
      return 5;
    }

    const_IWSubstring s = cl.string_value('S');

    if (stream_for_passing_molecules.would_overwrite_input_files(cl, s))
    {
      cerr << "Cannot overwrite input file(s) '" << s << "'\n";
      return 7;
    }

    if (! stream_for_passing_molecules.new_stem(s))
    {
      cerr << "Cannot open passing stream '" << s << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Passing molecules written to '" << s << "'\n";
  }

  if (cl.option_present('B'))
  {
    if (! cl.option_present('o'))
      stream_for_rejected_molecules.add_output_type(FILE_TYPE_SMI);
    else
      stream_for_rejected_molecules.determine_output_types(cl);

    const_IWSubstring b = cl.string_value('B');

    if (stream_for_rejected_molecules.would_overwrite_input_files(cl, b))
    {
      cerr << "Cannot overwrite input file(s) '" << b << "'\n";
      return 7;
    }

    if (! stream_for_rejected_molecules.new_stem(b))
    {
      cerr << "Cannot open rejected stream '" << b << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Rejected molecules written to '" << b << "'\n";
  }

  if (! stream_for_passing_molecules.active() && ! stream_for_rejected_molecules.active())
  {
    cerr << "NO OUTPUTS SELECTED (-S or -B)\n";
    if (0 == verbose)
      verbose = 1;
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! fragment_filter(cl[i], input_type, std::cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    cerr << molecules_passing << " passed, " << molecules_failing << " failed\n";
    cerr << rejected_for_non_allowed_atom << " contained non allowed element\n";
    cerr << rejected_for_non_allowed_counterion << " contained non allowed counterions\n";
    cerr << rejected_for_non_organic << " were non organic\n";
    cerr << rejected_for_isotopes << " contained isotopes\n";
    cerr << rejected_for_bad_valence << " had bad valences\n";
    cerr << rejected_for_too_few_atoms << " had too few atoms\n";
    cerr << rejected_for_too_many_atoms <<  " had too many atoms\n";
    cerr << rejected_for_too_few_rings <<   " had too few rings\n";
    cerr << rejected_for_too_many_rings <<  " had too many rings\n";
    cerr << rejected_for_too_many_ring_systems << " too many ring systems\n";
    cerr << rejected_for_ring_system_too_large << " ring system too large \n";
    cerr << rejected_for_min_heteroatom_fraction << " too few heteroatoms (fraction)\n";
    cerr << rejected_for_max_heteroatom_fraction << " too many heteroatoms (fraction)\n";
    cerr << rejected_for_too_many_heteroatoms << " too many heteroatoms\n";
    cerr << rejected_for_no_features <<    " had no features\n";
    cerr << rejected_for_donor_acceptor << " had too many heteroatoms with H's\n";
    cerr << rejected_for_too_much_spinach << " 3 ring systems with too much spinach\n";
    cerr << rejected_for_between_ring << " between ring\n";
    cerr << rejected_for_too_many_nitrogens << " too many nitrogens\n";

    for (IW_STL_Hash_Map_int::const_iterator i = dis_allowed_counterions.begin(); i != dis_allowed_counterions.end(); ++i)
    {
      if ((*i).second > 0)
        cerr << " discarded " << (*i).second << " molecules containing '" << (*i).first << "'\n";
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = fragment_filter(argc, argv);

  return rc;
}
