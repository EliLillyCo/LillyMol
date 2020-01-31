/*
  Implementation of extended connectivity
  Takes two .sdf files, that must contain the ligand
  and the corresponding protein

  Forms bonds between the ligand and the protein, then forms
  fingerprints
*/

#include <stdlib.h>

#include <assert.h>
#include <memory>
#include <algorithm>
#include <vector>
#include <limits>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION
#include "iwaray.h"

#include "cmdline.h"
#include "iwqsort.h"
#include "sparse_fp_creator.h"
#include "accumulator.h"
#include "misc.h"
#include "iw_stl_hash_map.h"

#include "molecule.h"
#include "iwstandard.h"
#include "aromatic.h"
#include "istream_and_type.h"
#include "atom_typing.h"
#include "qry_wstats.h"
#include "target.h"
#include "rwsubstructure.h"
#include "output.h"
#include "iwmfingerprint.h"

using std::unique_ptr;

static Chemical_Standardisation chemical_standardisation;

static int verbose=0;

static int molecules_read = 0;

static resizable_array_p<Substructure_Hit_Statistics> protein_atom_query;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString ligand_tag, composite_tag;

static Atom_Typing_Specification ligand_atom_type, protein_atom_type;

static int min_shell_radius = 0;
//static int max_shell_radius = 0;
static int max_shell_radius = std::numeric_limits<int>::max();

static int radius_into_protein = 0;

static int additive = 1;

static Accumulator_Int<int> nbits_acc;

static int whole_molecule_closest_protein_atom_only = 0;
static int each_ligand_atom_closest_protein_atom_only = 0;

static int remove_hydrogen_from_ligand = 0;
static int remove_hydrogen_from_protein = 0;

static int produce_radii_fingerprint = 0;

static int produce_atom_pair_fingerprints = 0;
static int produce_atom_path_fingerprints = 0;

static bond_type_t intermolecular_bond = SINGLE_BOND;

static resizable_array<int> bits_weighted_by_distance_to_interacting_atoms;

/*
  Sometimes we want to determine the structural features that give rise
  to certain bits.
*/

static IW_STL_Hash_Map<unsigned int, unsigned int> bits_to_investigate;
static int looking_for_bit_meanings = 0;
static int bits_found = 0;

static IWString_and_File_Descriptor stream_for_bit_meanings;

static distance_t max_intermolecular_distance = static_cast<distance_t>(0.0);

static Molecule_Output_Object stream_for_joined_molecules;

static extending_resizable_array<int> interactions_per_molecule;

static Accumulator<double> coverage;

static std::ofstream descriptor_file_output;
static char dfile_separator = ' ';

/*
  We might want to only generate shells from atoms that are within a given bond range of
  an interacting atom
*/

static int shells_start_with_atoms_in_range = -1;

static int discard_molecules_with_no_bits_set = 1;
static int molecules_discarded_for_no_bits_set = 0;

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Compute the Extended Connectivity fingerprints for molecules\n";
  cerr << "Usage <options> protein ligand1 ligand2 ligand3..., or use the -l option\n";
  cerr << "  -l <query>     ligand and protein in same file, ligand is what matches <query>\n";
  cerr << "  -T <dist>      make intermolecular connections within distance <dist>\n";
  cerr << "  -b             intermolecular bonds are special coordination type bond (single is default)\n";
  cerr << "  -c <n>         for each atom, only consider the <n> closest atoms in the protein\n";
  cerr << "  -C <n>         across all atoms, only consider the <n> closest atoms in the protein\n";
  cerr << "  -p <query>     only join to atoms that match <query> in the protein\n";
  cerr << "  -r <len>       minimum shell radius\n";
  cerr << "  -R <length>    maximum shell radius\n";
  cerr << "  -e <rad>       only fingerprint ligand atoms within <rad> bonds of an interacting atom\n";
  cerr << "  -h <rad>       only protein atoms within <rad> of an ligand atom are included in shells\n";
  cerr << "  -m             multiplicative formation of bits\n";
//cerr << "  -s             each radius gets its own fingerprint\n";
  cerr << "  -J <tag>       set the tag for the name tag of the fingerprints\n";
  cerr << "  -j <tag>       produce extra fingerprint of just the ligand\n";
//cerr << "  -y ...         intermolecular distance fingerprint specification\n";
//cerr << "  -Y <tag>       tag for radius derived fingerprints\n";
  cerr << "  -d <n>         produce bits corresponding to closest distance\n";
  cerr << "  -k             produce atom pair fingerprints (-R is max separation)\n";
  cerr << "  -a             produce atom path fingerprints (-R is max path length)\n";
  cerr << "  -H ...         remove hydrogen from 'ligand', 'protein' or 'both'\n";
  cerr << "  -P ...         atom typing specification (default atomic numbers)\n";
  cerr << "  -I <fname>     write molecules with inter-molecular bonds to <fname>\n";
  cerr << "  -D <fname>     produce a descriptor file with summary data\n";
  cerr << "  -E ...         standard element options\n";
  cerr << "  -i <type>      input type\n";
  cerr << "  -v             verbose output\n";
  
  exit(rc);
}

static int
read_bits_to_investigate (iwstring_data_source & input,
                          IW_STL_Hash_Map<unsigned int, unsigned int> & bits_to_investigate)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#') || 0 == buffer.length())
      continue;
    buffer.truncate_at_first(' ');

    unsigned int b;
    if (! buffer.numeric_value(b))
    {
      cerr << "Invalid bit number '" << buffer << "', line " << input.lines_read() << endl;
      return 0;
    }

    bits_to_investigate[b] = 0;
  }

  return static_cast<int>(bits_to_investigate.size());
}

static int
read_bits_to_investigate (const char * fname,
                          IW_STL_Hash_Map<unsigned int, unsigned int> & bits_to_investigate)
{
  iwstring_data_source input (fname);
  if (! input.good())
  {
    cerr << "Cannot open bits to be examined file name '" << fname << "'\n";
    return 0;
  }

  return read_bits_to_investigate (input, bits_to_investigate);
}

static int
check_against_list (const IWString & name_of_current_molecule,
                    const IWString & smarts,
                    atom_number_t centre_of_shell,
                    unsigned int sum_so_far,
                    int radius)
{
  IW_STL_Hash_Map<unsigned int, unsigned int>::iterator f = bits_to_investigate.find(sum_so_far);

  if (f == bits_to_investigate.end())
    return 0;

  bits_found++;

  (*f).second++;

  stream_for_bit_meanings << name_of_current_molecule << " bit " << sum_so_far << " atom " << centre_of_shell << " " << smarts << " radius " << radius << '\n';

  stream_for_bit_meanings.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
preprocess_molecule (Molecule & m)
{
  m.remove_all_chiral_centres();
  m.revert_all_directional_bonds_to_non_directional();
  
  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return 1;
}

static int
bond_constant (const Bond & b)
{
  int rc;

  if (b.is_aromatic())
    rc = 17;
  else if (b.is_triple_bond())
    rc = 11;
  else if (b.is_double_bond())
    rc = 7;
  else if (COORDINATION_BOND == b.btype())
    rc = 5;
  else
    rc = 3;

  return rc;
}

static int 
bond_constant (const Bond * bondi,
               int iso1,
               int iso2)
{
  int rc;

  if (bondi->is_aromatic())
    rc = 17;
  else if (bondi->is_triple_bond())
    rc = 11;
  else if (bondi->is_double_bond())
    rc = 7;
  else if (COORDINATION_BOND == bondi->btype())
    rc = 5;
  else
    rc = 3;

  if (iso1 != iso2)
    rc = rc + rc;

  return rc;
}

static void
increment(unsigned int & sum_so_far,
          int bc,
          int atom_constant)
{
  if (additive)
    sum_so_far += bc * atom_constant;
  else
    sum_so_far *= bc * atom_constant;

  return;
}

static void
set_bit (Sparse_Fingerprint_Creator & sfc,
         unsigned int b,
         const int replicates)
{
  if (1 == replicates)
    sfc.hit_bit(b);
  else
  {
    for (auto i = 0; i < replicates; ++i)
    {
      sfc.hit_bit(b + i * 39212);
    }
  }
}

/*
  For the mode in which we look for bits, we need to know the atom
  at the centre of the feature. Ideally, this would be passed as
  an argument, but the argument list is already long enough, so
  we do a very bad thing and just make it a file scope static
  variable
*/

static IWString name_of_current_molecule;
static atom_number_t centre_of_shell = INVALID_ATOM_NUMBER;
static IWString smarts_for_centre_of_shell;

#define PROCESSING_FINISHED 1
#define READY_TO_PROCESS 2
#define NEXT_TIME 3

static int
generate_shells(int matoms,
                int central_atom,
                int radius,
                const Atom * const * a,
                const int * atom_constant,
                int * processing_status,
                unsigned int sum_so_far,
                const int * atom_weight,
                Sparse_Fingerprint_Creator & sfc)
{
  radius++;

  if (additive)
    sum_so_far *= 7879;   // an arbitrary prime number

//#define DEBUG_ECFP_BITS 1
#ifdef DEBUG_ECFP_BITS
  cerr << "On entry sum_so_far " << sum_so_far << " radius " << radius << endl;
#endif

  for (auto i = 0; i < matoms; i++)
  {
    if (READY_TO_PROCESS != processing_status[i])
      continue;

    const Atom * ai = a[i];

    const auto acon = ai->ncon();

    const Bond * const * bonds = ai->rawdata();    // for efficiency

    for (auto j = 0; j < acon; j++)
    {
      const Bond * b = bonds[j];

      const auto k = b->other(i);

      if (PROCESSING_FINISHED == processing_status[k])   // we are extending the shell
      {
        int bc = bond_constant(b, a[k]->isotope(), ai->isotope());
        increment(sum_so_far, bc, atom_constant[i]);
//      cerr << "at level " << radius << " added atom " << i << " ssf " << sum_so_far << " be " << bc << endl;
      }
      else if (READY_TO_PROCESS == processing_status[k])
        ;
      else
        processing_status[k] = NEXT_TIME;
    }
  }

  if (radius > min_shell_radius)
  {
#ifdef DEBUG_ECFP_BITS
    cerr << "Hit bit " << sum_so_far << " at radius " << radius << endl;
#endif
    set_bit(sfc, sum_so_far, atom_weight[central_atom]);

    sfc.hit_bit(sum_so_far);
    if (looking_for_bit_meanings)
      check_against_list(name_of_current_molecule, smarts_for_centre_of_shell, centre_of_shell, sum_so_far, radius);
  }

  if (radius >= max_shell_radius)
    return 1;

  int continue_processing = 0;

  for (auto i = 0; i < matoms; i++)
  {
    if (READY_TO_PROCESS == processing_status[i])
      processing_status[i] = PROCESSING_FINISHED;
    else if (NEXT_TIME == processing_status[i])
    {
      processing_status[i] = READY_TO_PROCESS;
      continue_processing = 1;
    }
  }

  if (! continue_processing)
    return 1;

  return generate_shells(matoms, central_atom, radius, a, atom_constant, processing_status, sum_so_far, atom_weight, sfc);
}


static int
write_fingerprint (Sparse_Fingerprint_Creator & sfc,
                   const IWString & tag,
                   IWString & output)
{
  if (verbose)
    nbits_acc.extra(sfc.nbits());

  IWString tmp;

  sfc.daylight_ascii_form_with_counts_encoded(tag, tmp);
  output << tmp << '\n';

  return output.size();
}

static int
add_radius_fingerprints (const std::vector<distance_t> & radii,
                         Sparse_Fingerprint_Creator & sfc)
{
  distance_t r0 = radii[0];
  
  int ndx = (r0 - 1.0) / 0.1;

  auto bstart = 2161;    // just hoping to find a bit not otherwise set

  for (auto i = 0; i < produce_radii_fingerprint; ++i)
  {
    sfc.hit_bit(bstart + i * 722, ndx);
  }

  return 1;
}

static int
generate_atom_paths (Molecule & m,
                     const int * atype,
                     const int ligand_matoms,
                     const bool only_fingerprint_isotopic_atoms,
                     const int * atom_weight,
                     Sparse_Fingerprint_Creator & sfc)
{
  IWMFingerprint fp;

  fp.construct_fingerprint(m, atype);

  const auto nbits = iwmfingerprint_nbits();

  const int * v = fp.vector();

  for (auto i = 0; i < nbits; ++i)
  {
    if (v[i])
      sfc.hit_bit(i, v[i]);
  }

  return sfc.nbits();
}

static void
form_atom_pair_bonded (int t1, int t2,
                       const Bond & b,
                       Sparse_Fingerprint_Creator & sfc)
{
  if (t1 > t2)
    std::swap(t1, t2);

  auto bc = bond_constant(b);

  const unsigned int x = t1 * 7121 + bc + 212 * t2;

  sfc.hit_bit(x);

  return;
}

static void
form_atom_pair_not_bonded(int t1, int t2, 
                          const int d,
                          Sparse_Fingerprint_Creator & sfc)
{
  if (t1 > t2)
    std::swap(t1, t2);

  unsigned int b = t1 * 312 + d * t2;

  sfc.hit_bit(b);

  return;
}

static int
generate_atom_pairs (Molecule & m,
                     const int * atype,
                     const int ligand_matoms,
                     const bool only_fingerprint_isotopic_atoms,
                     const int * atom_weight,
                     Sparse_Fingerprint_Creator & sfc)
{
  const auto matoms = m.natoms();

  for (auto i = 0; i < ligand_matoms; ++i)
  {
    if (only_fingerprint_isotopic_atoms && 0 == m.isotope(i))
      continue;

    for (auto j = i + 1; j < matoms; ++j)
    {
      int d = m.bonds_between(i, j);
      if (d > max_shell_radius)
        continue;

      if (1 == d)
      {
        const Bond * b = m.bond_between_atoms(i, j);
        form_atom_pair_bonded(atype[i], atype[j], *b, sfc);
      }
      else
        form_atom_pair_not_bonded(atype[i], atype[j], d, sfc);
    }
  }

  return 1;
}

/*
*/

static int
produce_fingerprint (Molecule & m,
                     bool only_fingerprint_isotopic_atoms,
                     const int * atom_constant,
                     const int ligand_matoms,
                     const std::vector<distance_t> & radii,
                     const int * atom_weight,
                     const IWString & tag,
                     IWString & output)
{
  const auto matoms = m.natoms();

  Sparse_Fingerprint_Creator sfc;

  if (produce_atom_pair_fingerprints)
    generate_atom_pairs (m, atom_constant, ligand_matoms, only_fingerprint_isotopic_atoms, atom_weight, sfc);
  else if (produce_atom_path_fingerprints)
    generate_atom_paths (m, atom_constant, ligand_matoms, only_fingerprint_isotopic_atoms, atom_weight, sfc);
  else
  {
    int * processing_status = new int[matoms]; unique_ptr<int> free_processing_status(processing_status);

    Atom ** atoms = new Atom * [matoms]; unique_ptr<Atom *> free_atoms(atoms);

    m.atoms ((const Atom **) atoms);   // disregard of const OK

    for (auto i = 0; i < ligand_matoms; i++)
    {
#ifdef DEBUG_ECFP_BITS
      cerr << "only_fingerprint_isotopic_atoms " << only_fingerprint_isotopic_atoms << " atom " << i << " iso " << m.isotope(i) << endl;
#endif

      if (only_fingerprint_isotopic_atoms && 0 == m.isotope(i))
        continue;
    
      unsigned int e = atom_constant[i];

#ifdef DEBUG_ECFP_BITS
      if (0 == min_shell_radius)
        cerr << "Starting with atom " << i << " bit " << e << endl;
#endif
      if (0 == min_shell_radius)
        set_bit(sfc, e, atom_weight[i]);

      set_vector(processing_status, matoms, 0);

      processing_status[i] = PROCESSING_FINISHED;

      const Atom * ai = atoms[i];

      int acon = ai->ncon();

      for (auto j = 0; j < acon; j++)
      {
        atom_number_t k = ai->other(i, j);

        processing_status[k] = READY_TO_PROCESS;
      }

      generate_shells(matoms, i, 0, atoms, atom_constant, processing_status, e, atom_weight, sfc);
    }
  }

  if (0 == sfc.nbits())
    return 0;

  if (produce_radii_fingerprint && radii.size() > 0)
    add_radius_fingerprints(radii, sfc);

  return write_fingerprint (sfc, tag, output);
}

static void
determine_distance_from_interacting_atoms (Molecule & m,
                                           const resizable_array<int> w, 
                                           int * atom_weight)
{
  const auto matoms = m.natoms();

  for (auto i = 0; i < matoms; ++i)
  {
    if (m.isotope(i))
      atom_weight[i] = w[0];
    else
    {
      int shortest_distance = matoms;

      for (auto j = 0; j < matoms; ++j)
      {
        if (i == j)
          continue;

        auto d = m.bonds_between(i, j);

        if (d < shortest_distance)
          shortest_distance = d;
      }

      if (shortest_distance >= w.number_elements())
        atom_weight[i] = 1;
      else
        atom_weight[i] = w[shortest_distance];
    }
  }

  return;
}

struct Atom_and_Distance
{
  atom_number_t _a1;
  atom_number_t _a2;
  distance_t _dist;
};

class Atom_and_Distance_Comparator
{
  private:
  public:
    int operator() (const Atom_and_Distance & ad1, const Atom_and_Distance & ad2) const
    {
      if (ad1._dist < ad2._dist)
        return -1;
      else if (ad1._dist > ad2._dist)
        return 1;
      else
        return 0;
    }
};

static Atom_and_Distance_Comparator adc;

/*
  There are several different ways in which we can add the inter-molecular bonds:

    All separations that are within the distance threshold
    The N closest bonds, regardless of which atom they involve
    For each ligand atom, the N closest atoms
*/

static void
add_intra_molecular_bonds (Molecule & m,
                           const int ligand_matoms,
                           Atom_and_Distance * ad,
                           std::vector<distance_t> & radii)
{
  const auto matoms = m.natoms();

  Accumulator_Int<int> per_atom_interactions;    // keep track of the attachments to the atoms that make interactions

  if (each_ligand_atom_closest_protein_atom_only)
  {
    for (auto i = 0; i < ligand_matoms; ++i)
    {
      unsigned int ndx = 0;

      for (auto j = ligand_matoms; j < matoms; ++j)
      {
//      cerr << "Atoms " << i << " and " << j << " iso " << m.isotope(j) << " d " << m.distance_between_atoms(i, j) << endl;
        if (0 == m.isotope(j))
          continue;

        distance_t d = m.distance_between_atoms(i, j);
        if (d > max_intermolecular_distance)
          continue;

        ad[ndx]._a1 = i;    // not really needed, always the same in these loops
        ad[ndx]._a2 = j;
        ad[ndx]._dist = d;
        ndx++;
      }

      if (0 == ndx)
        continue;

      per_atom_interactions.extra(ndx);

      if (ndx > 1)
        iwqsort(ad, ndx, adc);

      int jstop = each_ligand_atom_closest_protein_atom_only;
      if (jstop > ndx)
        jstop = ndx;

      for (auto j = 0; j < jstop; ++j)
      {
        const Atom_and_Distance & adj = ad[j];
        m.add_bond(adj._a1, adj._a2, intermolecular_bond);
        m.set_isotope(adj._a1, 1);
        radii.push_back(m.distance_between_atoms(adj._a1, adj._a2));
      }
    }
  }
  else if (whole_molecule_closest_protein_atom_only)
  {
    const auto matoms = m.natoms();    // the size of our AD array

    unsigned int ndx = 0;

    for (auto i = 0; i < ligand_matoms; ++i)
    {
      auto interactions_this_atom = 0;

      for (auto j = ligand_matoms; j < matoms; ++j)
      {
//      cerr << "Atoms " << i << " and " << j << " iso " << m.isotope(j) << " d " << m.distance_between_atoms(i, j) << endl;
        if (0 == m.isotope(j))
          continue;

        distance_t d = m.distance_between_atoms(i, j);
        if (d > max_intermolecular_distance)
          continue;

        ad[ndx]._a1 = i;
        ad[ndx]._a2 = j;
        ad[ndx]._dist = d;
        ndx++;
        interactions_this_atom++;
        if (ndx >= matoms)      // we are about to overflow our array
        {
          cerr << "Too many interactions at distance " << max_intermolecular_distance << " search truncated\n";
          break;
        }
      }

      if (0 == interactions_this_atom)
        continue;

      per_atom_interactions.extra(interactions_this_atom);
    }

    if (ndx > 0)
    {
      if (ndx > 1)
        iwqsort(ad, ndx, adc);

      int jstop = whole_molecule_closest_protein_atom_only;
      if (ndx < whole_molecule_closest_protein_atom_only)
        jstop = ndx;

//    cerr << "Adding " << jstop << " bonds\n";

      for (auto j = 0; j < jstop; ++j)
      {
        m.add_bond(ad[j]._a1, ad[j]._a2, intermolecular_bond);
        m.set_isotope(ad[j]._a1, 1);
        radii.push_back(m.distance_between_atoms(ad[j]._a1, ad[j]._a2));
      }
    }
  }
  else
  {
    for (auto i = 0; i < ligand_matoms; ++i)
    {
      unsigned int ndx = 0;

      auto interactions_this_atom = 0;

      for (auto j = ligand_matoms; j < matoms; ++j)
      {
//      cerr << "Atoms " << i << " and " << j << " iso " << m.isotope(j) << " d " << m.distance_between_atoms(i, j) << endl;
        if (0 == m.isotope(j))
          continue;

        distance_t d = m.distance_between_atoms(i, j);
        if (d > max_intermolecular_distance)
          continue;

        m.add_bond(i, j, intermolecular_bond);
        m.set_isotope(i, 1);
        radii.push_back(d);
        interactions_this_atom++;
      }

      if (0 == interactions_this_atom)
        continue;

      per_atom_interactions.extra(interactions_this_atom);
    }
  }

  if (verbose)
  {
    double r = static_cast<double>(per_atom_interactions.n()) / static_cast<double>(ligand_matoms);
    coverage.extra(r);
    cerr << m.name() << ' ' << per_atom_interactions.n() << " of " << ligand_matoms << " atoms made interactions " << static_cast<float>(r) << ", " << radii.size() << " bonds added\n";
  }

  if (descriptor_file_output.is_open())
  {
    write_space_suppressed_string(m.name(), descriptor_file_output);

    double r = static_cast<double>(per_atom_interactions.n()) / static_cast<double>(ligand_matoms);
    descriptor_file_output << dfile_separator << ad[0]._dist;
    descriptor_file_output << dfile_separator << per_atom_interactions.n();
    descriptor_file_output << dfile_separator << static_cast<float>(r);
    descriptor_file_output << dfile_separator << per_atom_interactions.minval();
    descriptor_file_output << dfile_separator << per_atom_interactions.maxval();
    descriptor_file_output << dfile_separator << static_cast<float>(per_atom_interactions.average());

    descriptor_file_output << "\n";
  }

  return;
}


static int
place_isotopes_on_matched_atoms (Molecule & m,
                        const resizable_array_p<Substructure_Hit_Statistics> & queries,
                        int iso)
{
  Molecule_to_Match target(&m);

  int rc = 0;

  for (auto i = 0; i < queries.number_elements(); ++i)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    for (auto j = 0; j < nhits; ++j)
    {
      const Set_of_Atoms * e = sresults.embedding(j);
      m.set_isotope(*e, iso);
      rc = 1;
    }
  }

  return rc;
}

/*
  All atoms that cannot participate in a fingerprint are removed.
  We use the can_be_part_of_shell to form a list of the atoms to
  be removed, so we need to invert it.
  Note that we do not restore it on exit, it is changed
*/

static int
remove_non_participating_atoms (Molecule & m,
                                int & ligand_matoms,
                                int * atom_constant,
                                int * atom_weight,
                                int * can_be_part_of_shell)
{
  const auto matoms = m.natoms();

  int new_ligand_matoms = -1;

  int ndx = 0;

  for (auto i = 0; i < matoms; ++i)
  {
    if (can_be_part_of_shell[i])
    {
      atom_constant[ndx] = atom_constant[i];
      atom_weight[ndx] = atom_weight[i];
      ndx++;
      can_be_part_of_shell[i] = 0;     // will not be removed
    }
    else
      can_be_part_of_shell[i] = 1;    // will be removed

    if (ligand_matoms - 1 == i)
      new_ligand_matoms = ndx;
  }

  ligand_matoms = new_ligand_matoms;

//return m.remove_atoms(can_be_part_of_shell);
  return m.remove_many_atoms(can_be_part_of_shell);
}

/*
  We start with all atoms in the ligand with their REMOVE value set to zero. We want to follow
  all attached atoms out and mark them as not being removed. 
*/

static void
identify_attached_atoms_recursive (const Molecule & m,
                         atom_number_t zatom,
                         int * remove,
                         int radius)
{
  const Atom * ai = m.atomi(zatom);

  int acon = ai->ncon();

  for (auto i = 0; i < acon; ++i)
  {
    atom_number_t j = ai->other(zatom, i);

    if (0 == remove[i])
      continue;

    remove[j] = 0;

    if (radius > 1)
      identify_attached_atoms_recursive (m, j, remove, radius - 1);
  }
}

/*
  A couple of functions that do basically the same thing - extend the ligand into the protein
*/

static void
identify_attached_atoms (const Molecule & m,
                         const int ligand_matoms,
                         int * remove,
                         int radius)
{
  for (auto i = 0; i < ligand_matoms; ++i)
  {
    const Atom * a = m.atomi(i);

    const auto acon = a->ncon();

    for (auto j = 0; j < acon; ++j)
    {
      const auto k = a->other(i, j);

      if (0 == remove[k])   // alread marked for keeping
        continue;

      remove[k] = 0;

      if (radius > 1)
        identify_attached_atoms_recursive(m, k, remove, radius - 1);
//    cerr << "At radius " << radius << " from atom " << i << " marked atom " << k << " to keep\n";
    }
  }

  return;
}

static int
place_isotopes_on_atoms_within_range_of_interacting_atoms_atom (Molecule & m,
                                                  const atom_number_t zatom,
                                                  int shells_start_with_atoms_in_range)
{
  const Atom * a = m.atomi(zatom);

  const auto acon = a->ncon();

  int rc = 1;

  for (auto i = 0; i < acon; ++i)
  {
    const auto j = a->other(zatom, i);

    if (m.isotope(j))
      continue;

    m.set_isotope(j, 1);
    rc++;

    if (shells_start_with_atoms_in_range > 0)
      rc += place_isotopes_on_atoms_within_range_of_interacting_atoms_atom(m, j, shells_start_with_atoms_in_range - 1);
  }

  return rc;
}

static int
place_isotopes_on_atoms_within_range_of_interacting_atoms (Molecule & m,
                                                  const int ligand_matoms,
                                                  int shells_start_with_atoms_in_range)
{
  auto rc = 0;

  for (auto i = 0; i < ligand_matoms; ++i)
  {
    if (0 == m.isotope(i))
      continue;

    rc += place_isotopes_on_atoms_within_range_of_interacting_atoms_atom(m, i, shells_start_with_atoms_in_range);
  }

  return rc;
}

static int
identify_atoms_within_range_of_interacting_atoms_atom (const Molecule & m,    // /ligand + protein
                                                  const atom_number_t zatom,
                                                  int * can_be_part_of_shell,
                                                  int radius_into_protein)
{
  const Atom * a = m.atomi(zatom);

  const auto acon = a->ncon();

  int rc = 0;

  for (auto i = 0; i < acon; ++i)
  {
    atom_number_t j = a->other(zatom, i);

    if (can_be_part_of_shell[j])
      continue;

    rc++;

//  cerr << "From atom " << zatom << " go to atom " << j << endl;

    can_be_part_of_shell[j] = 1;

    if (radius_into_protein > 0)
      rc += identify_atoms_within_range_of_interacting_atoms_atom(m, j, can_be_part_of_shell, radius_into_protein - 1);
  }

  return rc;
}

static int
identify_atoms_within_range_of_interacting_atoms (const Molecule & m,    // /ligand + protein
                                                  const int ligand_matoms,
                                                  int * can_be_part_of_shell,
                                                  int radius_into_protein)
{
  const auto matoms = m.natoms();

  int rc = 0;

  for (auto i = 0; i < ligand_matoms; ++i)
  {
    if (0 == m.isotope(i))
      continue;

    const Atom * a = m.atomi(i);

    const auto acon = a->ncon();

    for (auto j = 0; j < acon; ++j)
    {
      atom_number_t k = a->other(i, j);

      if (can_be_part_of_shell[k])
        continue;

      rc++;

      can_be_part_of_shell[k] = 1;

//    cerr << "From ligand atom " << i << " go to atom " << k << ", rad " << radius_into_protein << endl;

      if (radius_into_protein > 1)
        rc += identify_atoms_within_range_of_interacting_atoms_atom(m, k, can_be_part_of_shell, radius_into_protein - 1);
    }
  }

  return rc;
}

/*
  We are just about to build a molecule by adding protein to ligand - ligand atoms will
  be first in the new molecule, so put those atom types first in the array
*/

static int
do_atom_typing (Molecule & ligand,
                Molecule & protein,
                int * atom_constant,
                Atom_Typing_Specification & ligand_atom_type,
                Atom_Typing_Specification & protein_atom_type)
{
  ligand.compute_aromaticity_if_needed();
  protein.compute_aromaticity_if_needed();

  if (! ligand_atom_type.assign_atom_types(ligand, atom_constant))   // always do this, load at start of array
    return 0;

  if (protein_atom_type.active())
    return protein_atom_type.assign_atom_types(protein, atom_constant + ligand.natoms());

  return ligand_atom_type.assign_atom_types(protein, atom_constant + ligand.natoms());
}

/*
  We have a lot of variables floating around, and it gets awkward to pass them via
  long argument lists, so put them in a class
*/

class IWECFP_Intermolecular : public Molecule
{
  private:
    int _ligand_matoms;

  public:
    IWECFP_Intermolecular (Molecule & m, int l) : Molecule(m), _ligand_matoms(l) {};
    ~IWECFP_Intermolecular ();
};

/* 
  All calls come to this point
*/

static int
iwecfp_intermolecular (Molecule & protein, 
                       Molecule & ligand,
                       IWString_and_File_Descriptor & output)
{
  if (remove_hydrogen_from_ligand)
    ligand.remove_all(1);
  if (remove_hydrogen_from_protein)
    protein.remove_all(1);

  const IWString ligand_smiles(ligand.smiles());

  auto ligand_matoms = ligand.natoms();   // not const because may get changed when we remove non participating atoms

  if (verbose)
    cerr << "Ligand contains " << ligand_matoms << " atoms, protein " << protein.natoms() << endl;

  int * atom_constant = new_int(ligand_matoms + protein.natoms()); unique_ptr<int> free_atom_constant(atom_constant);

  if (! do_atom_typing (ligand, protein, atom_constant, ligand_atom_type, protein_atom_type))
  {
    cerr << "Cannot assign atom types '" << ligand.name() << "'\n";
    return 0;
  }

  std::vector<distance_t> radii;

  IWString buffer;    // store the results so if there are no bits, we can discard it
  buffer.resize(4096);

  buffer << smiles_tag << ligand.smiles() << ">\n";
  buffer << identifier_tag << ligand.name() << ">\n";

// For the things that can later change, set them to default values

  int * atom_weight = new_int(ligand_matoms + protein.natoms()); unique_ptr<int> free_atom_weight(atom_weight);
  std::fill_n(atom_weight, ligand_matoms, 1);

  if (ligand_tag.length())
    produce_fingerprint(ligand, false, atom_constant, ligand_matoms, radii, atom_weight, ligand_tag, buffer);

  if (protein_atom_query.number_elements())
  {
    if (0 == place_isotopes_on_matched_atoms (protein, protein_atom_query, 1))   // applies isotope 1 to all matched atoms
    {
      cerr << "None of " << protein_atom_query.number_elements() << " protein atom queries matched, cannot continue\n";
      return 0;
    }
  }

  Molecule copy_ligand(ligand);    // we might want to know bonds btw atoms in the ligand, but don't want to generate a distance matrix for the whole complex

  int * distance_matrix = NULL;
  if (produce_atom_pair_fingerprints)
  {
    distance_matrix = new int[ligand_matoms * ligand_matoms];
    ligand.recompute_distance_matrix();
//  std::copy_n (ligand.distance_matrix_warning_may_change(), ligand_matoms * ligand_matoms, distance_matrix);
    copy_vector(distance_matrix, ligand.distance_matrix_warning_may_change(), ligand_matoms * ligand_matoms);
  }

  ligand.add_molecule(&protein);

  int matoms = ligand.natoms();

// If there were no protein queries to apply isotopes, we do it

  if (0 == protein_atom_query.number_elements())
  {
    for (auto i = ligand_matoms; i < matoms; ++i)
    {
      ligand.set_isotope(i, 1);
    }
  }

  if (shells_start_with_atoms_in_range < 0)   // not used, so enable all atoms
  {
    for (auto i = 0; i < ligand_matoms; ++i)
    {
      ligand.set_isotope(i, 1);
    }
  }

// Now that the atom types are known, we can start adding bonds

  Atom_and_Distance * ad = new Atom_and_Distance[ligand.natoms()]; unique_ptr<Atom_and_Distance> free_ad(ad);

  add_intra_molecular_bonds(ligand, ligand_matoms, ad, radii);

  const auto bonds_added = radii.size();

  interactions_per_molecule[bonds_added]++;

  if (verbose > 1)
    cerr << "Added " << bonds_added << " inter molecular bonds to '" << ligand.name() << "'\n";

  if (0 == bonds_added && discard_molecules_with_no_bits_set)
  {
    molecules_discarded_for_no_bits_set++;
    return 1;
  }

  output << buffer;      // now it is OK to do some output

  if (shells_start_with_atoms_in_range >= 0)
    place_isotopes_on_atoms_within_range_of_interacting_atoms (ligand, ligand_matoms, shells_start_with_atoms_in_range);

  if (bits_weighted_by_distance_to_interacting_atoms.number_elements())
  {
    for (auto i = 0; i < ligand_matoms; ++i)
    {
      auto iso = ligand.isotope(i);
      if (iso)
        copy_ligand.set_isotope(i, iso);
    }

    determine_distance_from_interacting_atoms (copy_ligand, bits_weighted_by_distance_to_interacting_atoms, atom_weight);
  }

  int * can_be_part_of_shell = new_int(ligand.natoms(), 1); unique_ptr<int> free_can_be_part_of_shell(can_be_part_of_shell);

  std::fill_n(can_be_part_of_shell + ligand_matoms, matoms - ligand_matoms, 0);
//for (auto i = ligand_matoms; i < ligand.natoms(); ++i)
//{
//  can_be_part_of_shell[i] = 0;
//}

  if (radius_into_protein > 0)
    identify_atoms_within_range_of_interacting_atoms(ligand, ligand_matoms, can_be_part_of_shell, radius_into_protein);

  if (std::numeric_limits<int>::max() != max_shell_radius)
    identify_atoms_within_range_of_interacting_atoms(ligand, ligand_matoms, can_be_part_of_shell, max_shell_radius);

  if (verbose > 2)
    cerr << "Atoms in fp " << std::count(can_be_part_of_shell, can_be_part_of_shell + ligand.natoms(), 1) << " atoms\n";

  remove_non_participating_atoms (ligand, ligand_matoms, atom_constant, atom_weight, can_be_part_of_shell);

  if (verbose > 2)
    cerr << "After remove_non_participating_atoms, have " << ligand.natoms() << " atoms, in ligand " << ligand_matoms << endl;

  produce_fingerprint(ligand, true, atom_constant, ligand_matoms, radii, atom_weight, composite_tag, output);

  output << "|\n";

  output.write_if_buffer_holds_more_than(4096);

  if (! stream_for_joined_molecules.active())
    return 1;

  for (auto i = 0; i < ligand.nedges(); ++i)
  {
    const Bond * b = ligand.bondi(i);

    if (IS_COORDINATION_BOND(b->btype()))
    {
      ligand.set_bond_type_between_atoms(b->a1(), b->a2(), SINGLE_BOND);
    }
  }

  stream_for_joined_molecules.write(ligand);

  return 1;
}

static int
iwecfp_intermolecular (Molecule & m,
                       resizable_array_p<Substructure_Hit_Statistics> & ligand_fragment_query,
                       IWString_and_File_Descriptor & output)
{
  int fragment = -1;

  resizable_array_p<Molecule> c;
  m.create_components(c);

  for (auto i = 0; i < c.number_elements() && fragment < 0; ++i)
  {
    Molecule & ci = *(c[i]);

    if (ci.natoms() < 7 || ci.natoms() > 150)
      continue;

    if (first_query_to_match(ci, ligand_fragment_query) < 0)
      continue;

    fragment = i; 

    if (verbose)
      cerr << "Ligand is fragment with " << ci.natoms() << " atoms\n";

    break;
  }

  if (fragment < 0)
  {
    cerr << "Cannot identify ligand based on ligand fragment queries\n";
    return 0;
  }

  Molecule ligand;
  ligand.add_molecule(c[fragment]);
  ligand.set_name(m.name());

  m.delete_fragment(fragment);
  
  return iwecfp_intermolecular (m, ligand, output);
}

static int
iwecfp_intermolecular (data_source_and_type<Molecule> & input, 
                       resizable_array_p<Substructure_Hit_Statistics> & ligand_fragment_query,
                       IWString_and_File_Descriptor & output)
{
  Molecule * m;

  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    unique_ptr<Molecule> free_m(m);

    if (verbose > 1)
      cerr << molecules_read << " processing '" << m->name() << "'\n";

    preprocess_molecule(*m);

    if (! iwecfp_intermolecular (*m, ligand_fragment_query, output))
    {
      cerr << "Fatal error processing '" << m->name() << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
iwecfp_intermolecular (data_source_and_type<Molecule> & protein_input, 
                       data_source_and_type<Molecule> & ligand_input,
                       IWString_and_File_Descriptor & output)
{
  Molecule * protein = protein_input.next_molecule();

  if (NULL == protein)
  {
    cerr << "Cannot read protein molecule\n";
    return 0;
  }

  unique_ptr<Molecule> free_protein(protein);

  Molecule * ligand;

  while (NULL != (ligand = ligand_input.next_molecule()))
  {
    molecules_read++;

    unique_ptr<Molecule> free_ligand(ligand);

    if (verbose > 1)
      cerr << molecules_read << " processing '" << ligand->name() << "'\n";

    preprocess_molecule(*ligand);

    if (! iwecfp_intermolecular (*protein, *ligand, output))
    {
      cerr << "Fatal error processing '" << ligand->name() << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

/*static int
iwecfp_intermolecular (data_source_and_type<Molecule> & protein_input, 
                       data_source_and_type<Molecule> & ligand_input,
                       IWString_and_File_Descriptor & output)
{
  Molecule * ligand, * protein;

  while (NULL != (ligand = ligand_input.next_molecule()) && NULL != (protein = protein_input.next_molecule()))
  {
    molecules_read++;

    unique_ptr<Molecule> free_ligand(ligand);
    unique_ptr<Molecule> free_protein(protein);

    if (verbose > 1)
      cerr << molecules_read << " processing '" << ligand->name() << "'\n";

    preprocess_molecule(*ligand);

    if (! iwecfp_intermolecular (*protein, *ligand, output))
    {
      cerr << "Fatal error processing '" << ligand->name() << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}*/

static int
iwecfp_intermolecular (const char * fname,
                       int input_type,
                       resizable_array_p<Substructure_Hit_Statistics> & ligand_fragment_query,
                       IWString_and_File_Descriptor & output)
{
  data_source_and_type<Molecule> input (input_type, fname);

  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return iwecfp_intermolecular (input, ligand_fragment_query, output);
}

static int
iwecfp_intermolecular (const char * protein_fname,
                       const char * ligand_fname,
                       int input_type,
                       IWString_and_File_Descriptor & output)
{
  data_source_and_type<Molecule> ligand_input (input_type, ligand_fname);
  data_source_and_type<Molecule> protein_input (input_type, protein_fname);

  if (! ligand_input.ok() || ! protein_input.ok())
  {
    cerr << "Cannot open '" << ligand_fname << "' and/or '" << protein_fname << "'\n";
    return 0;
  }

  return iwecfp_intermolecular (protein_input, ligand_input, output);

}

static void
ensure_ends_with (IWString & s, char e)
{
  if (0 == s.length())
    return;

  if (s.ends_with(e))
    return;

  s << e;

  return;
}

static int
iwecfp_intermolecular (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vE:A:g:J:j:i:r:R:P:mT:p:c:C:I:H:l:y:Y:d:D:e:bh:kaW:z");
  
  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');
  
  if (! process_elements (cl))
  {
    usage(2);
  }

  MDL_File_Supporting_Material * mdlfsm = global_default_MDL_File_Supporting_Material ();
  mdlfsm->set_report_unrecognised_records(0);

  set_write_extra_text_info(1);

  if (! cl.option_present('T'))
  {
    cerr << "Must specify max intermolecular distance via the -T option\n";
    usage(2);
  }

  if (cl.option_present('T'))
  {
    if (! cl.value('T', max_intermolecular_distance) || max_intermolecular_distance < 1.0f)
    {
      cerr << "The maximum intermolecular distance considered (-T) option must be a reasonable inter molecular distance\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will only consider inter-molecular atoms separated by at most " << max_intermolecular_distance << " angstroms\n";
  }

  if (! cl.option_present('R') && ! cl.option_present('h'))
  {
    cerr << "Must specify one of -R and/or -h options in order to limit shell size\n";
    usage(2);
  }

// important that the -r options be processed early

  if (cl.option_present('r'))
  {
    if (! cl.value('r', min_shell_radius) || min_shell_radius < 0)
    {
      cerr << "The min shell radius (-r) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will only fingerprint paths larger than " << min_shell_radius << " bonds\n";
  }

  if (cl.option_present('R'))
  {
    if (! cl.value('R', max_shell_radius) || max_shell_radius < 1)
    {
      cerr << "The max shell radius (-R) must be a whole +ve number\n";
      usage(2);
    }

    if (max_shell_radius < min_shell_radius)
    {
      cerr << "Inconsistent min " << min_shell_radius << " and max " << max_shell_radius << " shell radius\n";
      usage(6);
    }

    if (verbose)
      cerr << "Max radius " << max_shell_radius << endl;
  }

  if (cl.option_present('c') && cl.option_present('C'))
  {
    cerr << "Sorry, cannot use both the -c and -C options\n";
    usage(1);
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', each_ligand_atom_closest_protein_atom_only) || each_ligand_atom_closest_protein_atom_only < 1)
    {
      cerr << "The each atom closest atoms specificatin (-c) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "For each ligand atom, will take only the first " << each_ligand_atom_closest_protein_atom_only << " atoms in the protein\n";
  }

  if (cl.option_present('C'))
  {
    if (! cl.value('C', whole_molecule_closest_protein_atom_only) || whole_molecule_closest_protein_atom_only < 1)
    {
      cerr << "The whole ligand closest atoms specificatin (-C) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Across the whole ligand, will take only the first " << whole_molecule_closest_protein_atom_only << " atoms in the protein\n";
  }

  if (cl.option_present('k'))
  { 
    produce_atom_pair_fingerprints = 1;

    if (verbose)
      cerr << "Will produce atom pair fingerprints\n";
  }
  else if (cl.option_present('a'))
  {
    produce_atom_path_fingerprints = 1;

    if (verbose)
      cerr << "Will produce atom path fingerprints\n";

    set_iwmfingerprint_nbits(200000);

    if (cl.option_present('R'))
      set_max_path_length (max_shell_radius);

    set_display_strange_chemistry_messages(0);
  }

  if (cl.option_present('b'))
  {
    intermolecular_bond = COORDINATION_BOND;

    if (verbose)
      cerr << "Coordination bond type used for intermolecular bonds\n";
  }

  if (cl.option_present('h'))
  {
    if (! cl.value('h', radius_into_protein) || radius_into_protein < 1)
    {
      cerr << "The radius into protein (-h) value must be a whole +ve number\n";
      usage(2);
    }

    if (radius_into_protein > max_shell_radius)
    {
      cerr << "The max radius into protein (-h) value must not be greater than the max shell radius " << max_shell_radius << endl;
      return 2;
    }

    if (verbose)
      cerr << "Will only include protein atoms within " << radius_into_protein << " bonds of the ligand\n";
  }

  if (cl.option_present('W'))
  {
    const_IWSubstring w = cl.string_value('W');

    const_IWSubstring token;

    for (auto i = 0; w.nextword(token, i, ','); )
    {
      int j;

      if (! token.numeric_value(j) || j < 0)
      {
        cerr << "Invalid shell relative weight '" << w << "'\n";
        return 2;
      }

      bits_weighted_by_distance_to_interacting_atoms.add(j);
    }
  }

  if (cl.option_present('p'))
  {
    if (! process_queries(cl, protein_atom_query, verbose, 'p'))
    {
      cerr << "Cannot discern protein atom substructure search specifications (-p)\n";
      usage(2);
    }

    if (verbose)
      cerr << "Read " << protein_atom_query.number_elements() << " protein atom query specifications\n";
  }

  if (cl.option_present('z'))
  {
    discard_molecules_with_no_bits_set = 0;

    if (verbose)
      cerr << "Molecules for which there are no 3d bits set will not be written\n";
  }

  resizable_array_p<Substructure_Hit_Statistics> ligand_fragment_query;

  if (cl.option_present('l'))
  {
    if (! process_queries(cl, ligand_fragment_query, verbose, 'l'))
    {
      cerr << "Cannot discern ligand fragment substructure search specifications (-l)\n";
      usage(2);
    }

    if (verbose)
      cerr << "Read " << ligand_fragment_query.number_elements() << " ligand fragment query specifications\n";
  }

  if (cl.option_present('e'))
  {
    if (! cl.value('e', shells_start_with_atoms_in_range) || shells_start_with_atoms_in_range < 0)
    {
      cerr << "The shells within a given distance of an interacting atom (-e) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Shells will start only with atoms within " << shells_start_with_atoms_in_range << " bonds of a fellow ligand atom that binds\n";
  }

  if (cl.option_present('d'))
  {
    if (! cl.value('d', produce_radii_fingerprint) || produce_radii_fingerprint < 1)
    {
      cerr << "The produce radius fingerprints option (-d) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will produce " << produce_radii_fingerprint << " replicates of radius fingerprints\n";
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line (cl, verbose > 1, 'g'))
    {
      cerr << "Cannot initialise chemical standardisation (-g)\n";
      usage(14);
    }
  }

  if (cl.option_present('m'))
  {
    additive = 0;

    if (verbose)
      cerr << "Fingerprints formed with multiplication operations\n";
  }

  if (cl.option_present('I'))
  {
    const_IWSubstring i = cl.string_value('I');

    stream_for_joined_molecules.add_output_type(SDF);

    if (stream_for_joined_molecules.would_overwrite_input_files(cl, i))
    {
      cerr << "Cannot overwrite input file(s) with '" << i << "'\n";
      return 2;
    }

    if (! stream_for_joined_molecules.new_stem(i))
    {
      cerr << "Cannot open joined molecule file '" << i << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Will write joined molecules to file '" << i << ".smi'\n";
  }

  if (cl.option_present('D'))
  {
    const char * d = cl.option_value('D');

    descriptor_file_output.open(d);

    if (! descriptor_file_output.good())
    {
      cerr << "Cannot open descriptor file string '" << d << "'\n";
      return 2;
    }

    if (verbose)
      cerr << "Data about interactions written to '" << d << "'\n";

    descriptor_file_output << "Name";
    descriptor_file_output << dfile_separator << "mindist";
    descriptor_file_output << dfile_separator << "AtomsInteracting";
    descriptor_file_output << dfile_separator << "FractionAtomsInteracting";
    descriptor_file_output << dfile_separator << "SmallestNumberProteinInteractions";
    descriptor_file_output << dfile_separator << "LargestNumberProteinInteractions";
    descriptor_file_output << dfile_separator << "AverageNumberProteinInteractions";
    descriptor_file_output << '\n';
  }

  if (cl.option_present('H'))
  {
    const_IWSubstring h;

    for (auto i = 0; cl.value('H', h, i); ++i)
    {
      if ("ligand" == h)
        remove_hydrogen_from_ligand = 1;
      else if ("protein" == h)
        remove_hydrogen_from_protein = 1;
      else if ("both" == h)
      {
        remove_hydrogen_from_ligand = 1;
        remove_hydrogen_from_protein = 1;
      }
    }
  }

  if (! cl.option_present('J'))
  {
    cerr << "Must specify fingerprint via the -J option - should start with NC\n";
    usage(1);
  }

  if (! cl.option_present('P'))
  {
    ligand_atom_type.set_user_specified_type(IWATTYPE_USP_Y);
    protein_atom_type.set_user_specified_type(IWATTYPE_USP_Y);
    if (verbose)
      cerr << "Using default atom type Y\n";
  }

// If neither atom type nor tag specified, take a default

  if (cl.option_present('P'))
  {
    const_IWSubstring p;
    for (auto i = 0; cl.value('P', p, i); ++i)
    {
      if (p.starts_with("LIG:"))
      {
        p.remove_leading_chars(4);
        if (! ligand_atom_type.build(p))
        {
          cerr << "Invalid ligand atom type specification '" << p << "'\n";
          return 2;
        }
      }
      else if (p.starts_with("PROT:"))
      {
        p.remove_leading_chars(5);
        if (! protein_atom_type.build(p))
        {
          cerr << "Invalid protein atom type specification '" << p << "'\n";
          return 2;
        }
      }
      else
      {
        if (! ligand_atom_type.build(p))
        {
          cerr << "Cannot discern atom typing specification '" << p << "'\n";
          return 1;
        }
      }
    }

    ligand_atom_type.swap_atomic_number_atom_type_to_atomic_number_prime();
    protein_atom_type.swap_atomic_number_atom_type_to_atomic_number_prime();
  }

  if (cl.option_present('J'))
  {
    cl.value('J', composite_tag);

    ensure_ends_with(composite_tag, '<');

    if (verbose)
      cerr<< "Extended connectivity index written as non-colliding sparse fingerprints, tag '"<< composite_tag <<"'\n";
  }
  
  if (cl.option_present('j'))
  {
    cl.value('j', ligand_tag);

    ensure_ends_with(ligand_tag, '<');

    if (verbose)
      cerr << "Fingerprint of ligand also produced\n";
  }

  if (ligand_tag == composite_tag)
  {
    cerr << "Tag for ligand (-j) and composite (-J) must be different\n";
    usage(2);
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
  else if (! all_files_recognised_by_suffix (cl))
  {
    cerr << "Cannot determine input type(s)\n";
    return 7;
  }

  set_global_aromaticity_type(Daylight);

  if (! cl.option_present('A'))
  {
    set_global_aromaticity_type(Daylight);
  }
  else if (! process_standard_aromaticity_options (cl, verbose))
  {
    cerr << "Cannot process aromaticity options (-A)\n";
    usage(5);
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  if (ligand_fragment_query.number_elements())
  {
    for (auto i = 0; i < cl.number_elements(); ++i)
    {
      if (! iwecfp_intermolecular(cl[i], input_type, ligand_fragment_query, output))
        return i + 1;
    }
  }
  else
  {
    for (auto i = 1; i < cl.number_elements(); ++i)
    {
      if (! iwecfp_intermolecular(cl[0], cl[i], input_type, output))
        return 1;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";

    if (nbits_acc.n() > 0)
      cerr << "Fingerprints had between " << nbits_acc.minval() << " and " << nbits_acc.maxval() << " ave " << static_cast<float>(nbits_acc.average_if_available_minval_if_not()) << " bits set\n";

    if (looking_for_bit_meanings)
    {
      cerr << "Found " << bits_found << " bits in lookup file\n";
      for (IW_STL_Hash_Map<unsigned int, unsigned int>::const_iterator i = bits_to_investigate.begin(); i != bits_to_investigate.end(); ++i)
      {
        cerr << "Found " << (*i).second << " instances of " << (*i).first << '\n';
      }
    }

    cerr << "Ligands had between " << coverage.minval() << " and " << coverage.maxval() << " fraction of atoms making interactions, ave " << static_cast<float>(coverage.average()) << endl;

    for (auto i = 0; i < interactions_per_molecule.number_elements(); ++i)
    {
      if (interactions_per_molecule[i])
        cerr << interactions_per_molecule[i] << " molecules had " << i << " interactions within " << max_intermolecular_distance << endl;
    }

    cerr << molecules_discarded_for_no_bits_set << " fingerprints not produced because of zero bits set\n";
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  int rc = iwecfp_intermolecular (argc, argv);
  return rc;
}
