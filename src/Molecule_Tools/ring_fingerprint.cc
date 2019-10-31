/*
  Compute fingerprints describing the aromatic rings of a molecule
*/

#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "cmdline.h"
#include "accumulator.h"
#include "misc.h"

#include "istream_and_type.h"
#include "molecule.h"
#include "path.h"
#include "aromatic.h"
#include "sparse_fp_creator.h"
#include "iwstandard.h"
#include "smiles.h"

static const char * prog_name = NULL;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

class Ring_Fingerprint
{
  private:
    int _verbose;
    int _molecules_read;
    Chemical_Standardisation _chemical_standardisation;
    int _reduce_to_largest_fragment;
    int _work_as_tdt_filter;
    int _form_bits_for_inter_aromatic_ring_distances;
    int _form_aliphatic_ring_bits;
    IWString _tag;

    int _ntest;
    int _echo_fingerprints;

// private functions

    void _default_values();
    void _usage (int rc);
    void _preprocess (Molecule & m);

    int _perform_tests (Molecule & m);

    int _process (Molecule & m, IWString_and_File_Descriptor & output);
    int _process (data_source_and_type<Molecule> & input,
                             IWString_and_File_Descriptor & output);
    int _process (const char * fname, int input_type,
                             IWString_and_File_Descriptor & output);
    int _write_fingerprint (Molecule & m, const Sparse_Fingerprint_Creator & sfc,
                                      IWString_and_File_Descriptor & output) const;

    int _process_filter (const char * fname, IWString_and_File_Descriptor & output);
    int _process_filter (iwstring_data_source & input, IWString_and_File_Descriptor & output);

    void _form_fingerprint (Molecule & m, Sparse_Fingerprint_Creator & sfc) const;

    int  _fingerprint_aromatic_ring  (Molecule & m, const int * arom, const Ring & r, Sparse_Fingerprint_Creator & sfc) const;
    int  _fingerprint_aliphatic_ring (Molecule & m, const Ring & r, Sparse_Fingerprint_Creator & sfc) const;
    void _discern_exocyclic_bond (Molecule & m, const int * arom, const atom_number_t zatom, int & exocyclic_double_bond, int & substituents, int & biphenyl, int & to_aliphatic_ring) const;
    int  _aromatic_ring_fusion_bits (Molecule & m, const Ring & r, Sparse_Fingerprint_Creator & sfc) const;
    void _ring_fusion_bit (int bstart, const int rs1, const int rs2, Sparse_Fingerprint_Creator & sfc) const;
    void _bits_for_distances_between_aromatic_rings (Molecule & m, Sparse_Fingerprint_Creator & sfc) const;
    void _bits_for_distances_aromatic_to_aliphatic (Molecule & m, const int *, Sparse_Fingerprint_Creator & sfc) const;
    int  _shortest_distance_between (Molecule & m, const Set_of_Atoms & r1, const Set_of_Atoms & r2) const;
    int _do_aliphatic_ring_bits (Molecule & m, const int * arom, Sparse_Fingerprint_Creator & sfc) const;
    void _do_single_aliphatic_ring (Molecule & m, const Ring & r, Sparse_Fingerprint_Creator & sfc) const;

    void _discern_exocyclic_bond_aliph (Molecule & m, const atom_number_t zatom,
                                int & exocyclic_double_bond,
                                int & substituents,
                                int & to_aliphatic_ring,
                                int & to_aromatic_ring,
                                int & within_ring_double_bonds) const;
    void _do_planar_aliphatic_system(Molecule & m, const resizable_array<int> & rings_in_system, Sparse_Fingerprint_Creator & sfc) const;
    void _do_strongly_fused_aliphatic_system(Molecule & m, const resizable_array<int> & rings_in_system, Sparse_Fingerprint_Creator & sfc) const;

  public:
    Ring_Fingerprint();

    int operator() (int argc, char ** argv);
};

Ring_Fingerprint::Ring_Fingerprint ()
{
  _default_values();

  return;
}

void
Ring_Fingerprint::_default_values()
{
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;

  _work_as_tdt_filter = 0;

  _form_bits_for_inter_aromatic_ring_distances = 0;
  _form_aliphatic_ring_bits = 0;

  _tag = "NCRNG<";

  _ntest = 0;
  _echo_fingerprints = 0;

  return;
}

void
Ring_Fingerprint::_usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Computes fingerprints based on ring properties - primarily aromatic rings\n";
  cerr << " -p            form aliphatic ring bits\n";
  cerr << " -d            form inter-ring distance bits\n";
  cerr << " -e            echo the bits in each fingerprint (debugging only)\n";
  cerr << " -J <tag>      fingerprint tag\n";
  cerr << " -l            reduce to largest fragment\n";
  cerr << " -i <type>     input specification\n";
  cerr << " -g ...        chemical standardisation options\n";
  cerr << " -E ...        standard element specifications\n";
  cerr << " -A ...        standard aromaticity specifications\n";
  cerr << " -v            verbose output\n";

  exit(rc);
}

void
Ring_Fingerprint::_preprocess (Molecule & m)
{
  if (_reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  m.remove_all_chiral_centres();

  if (_chemical_standardisation.active())
    _chemical_standardisation.process(m);

  return;
}

static int
write_fingerprint_differences (const Sparse_Fingerprint_Creator & sfc1,
                               const Sparse_Fingerprint_Creator & sfc2)
{ 
  const auto b1 = sfc1.bits_found();
  const auto b2 = sfc2.bits_found();

  int failures = 0;

  for (auto f = b1.cbegin(); f != b1.cend(); ++f)
  {
    const auto ib2 = b2.find(f->first);
    if (ib2 == b2.cend())
      cerr << "Bit " << f->first << " count " << f->second << " not found in second fp\n";
    else if (f->second != ib2->second)
      cerr << "bit " << f->first << " count mismatch " << f->second << " vs " << ib2->second << "\n";
    else
      continue;

    failures++;
  }

  for (auto f = b2.cbegin(); f != b2.cend(); ++f)
  {
    const auto ib1 = b1.find(f->first);

    if (f != b1.cend())
      continue;

    cerr << "bit " << f->first << " count " << f->second << " only found in second fp\n";
    failures++;
  }

  return failures;
}


int
Ring_Fingerprint::_perform_tests (Molecule & m)
{
  if (0 == m.aromatic_ring_count())
    return 1;

  const IWString initial_smiles(m.smiles());

  Sparse_Fingerprint_Creator sfc;

  _form_fingerprint (m, sfc);

  for (auto i = 0; i < _ntest; ++i)
  {
    const auto rsmi = m.random_smiles();

    Molecule mcopy;

    if (! mcopy.build_from_smiles(rsmi))
    {
      cerr << "Yipes, smiles interpretation failed, smiles '" << rsmi << "'\n";
      return 0;
    }

//  cerr << "Begin test " << i << endl;

    Sparse_Fingerprint_Creator mysfc;
    _form_fingerprint(mcopy, mysfc);

    if (sfc == mysfc)
      continue;

    cerr << "Test failure: " << initial_smiles << ' ' << m.name() << endl;
    cerr << "Random smiles " << rsmi << " permutation " << i << endl;
    cerr << "Initial bits\n";
    sfc.debug_print(cerr);
    cerr << "Random variant bits\n";
    mysfc.debug_print(cerr);
    write_fingerprint_differences(sfc, mysfc);

    return 0;
  }

  return 1;
}

void
Ring_Fingerprint::_discern_exocyclic_bond_aliph (Molecule & m,
                                const atom_number_t zatom,
                                int & exocyclic_double_bond,
                                int & substituents,
                                int & to_aliphatic_ring,
                                int & to_aromatic_ring,
                                int & within_ring_double_bonds) const
{
  const auto a= m.atomi(zatom);

  const auto acon = a->ncon();

  for (auto i = 0; i < acon; ++i)
  {
    const auto b = a->item(i);

    if (b->nrings())
    {
      if (1 == m.nrings(zatom) && ! b->is_single_bond())
        within_ring_double_bonds++;
      continue;
    }

    if (b->is_double_bond())
      exocyclic_double_bond++;
    else 
      substituents++;

    const auto j = b->other(zatom);

    if (m.is_aromatic(j))
      to_aromatic_ring++;
    else if (m.nrings(j))
      to_aliphatic_ring++;
  }

  return;
}

void
Ring_Fingerprint::_discern_exocyclic_bond (Molecule & m,
                                           const int * arom,
                                           const atom_number_t zatom,
                                           int & exocyclic_double_bond,
                                           int & substituents,
                                           int & biphenyl,
                                           int & to_aliphatic_ring) const
{
  const auto a = m.atomi(zatom);
  assert (3 == a->ncon());

  for (auto i = 0; i < 3; ++i)
  {
    const auto b = a->item(i);

    if (b->nrings())
      continue;

    if (b->is_double_bond())
    {
      exocyclic_double_bond++;
      return;
    }

    substituents++;

    const auto j = b->other(zatom);

    if (arom[j])
      biphenyl++;
    else if (m.nrings(j))
      to_aliphatic_ring++;

    return;
  }

  return;
}

int
Ring_Fingerprint::_write_fingerprint (Molecule & m,
                                      const Sparse_Fingerprint_Creator & sfc,
                                      IWString_and_File_Descriptor & output) const
{
  if (_echo_fingerprints)
    sfc.debug_print(cerr);

  if (! _work_as_tdt_filter)
  {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(_tag, tmp);

  output << tmp << '\n';

  if (! _work_as_tdt_filter)
    output << "|\n";

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

#define AROMATIC_RING_BIT 1
#define ALIPHATIC_RING_BIT 2
#define AROMATIC_CARBON 3
#define AROMATIC_NITROGEN_H0 4
#define AROMATIC_NITROGEN_H 5
#define AROMATIC_OXYGEN 6
#define AROMATIC_SULPHUR 7
#define AROMATIC_OS 8
#define AROMATIC_MORE_THAN_TWO_RINGS 9
#define AROMATIC_EXOCYCLIC_DOUBLE_BOND 10
#define AROMATIC_CONSECUTIVE_HETEROATOMS 11
#define AROMATIC_ISOLATED_RING 12
#define AROMATIC_NITROGEN_ANYH 13
#define AROMATIC_HETEROATOMS 14
#define AROMATIC_HETEROCYCLE 15
#define AROMATIC_HETEROCYCLE_FUSED_AROMATIC_HETEROCYCLE 16
#define AROMATIC_HETEROCYCLE_FUSED_AROMATIC_CARBOCYCLE 17
#define AROMATIC_CARBOCYCLE_FUSED_AROMATIC_CARBOCYCLE 18
#define AROMATIC_SUBSTITUENTS 19

#define AROMATIC_RING_SIZE_START 20
#define ALIPHATIC_RING_SIZE_START 30
#define AROMATIC_HETEROCYCLE_SIZE_BEGIN 40

#define AROMATIC_SUBSTITUENTS_SIZE 50

#define LARGEST_RING_SIZE 60

#define NON_SSSR_RINGS 65

#define ALIPHATIC_ISOLATED_RING 100

#define BETWEEN_AROMATIC_RINGS_MIN 200
#define BETWEEN_AROMATIC_RINGS_MAX 201
#define BETWEEN_AROMATIC_RINGS_AVE 202
#define BETWEEN_AROMATIC_RINGS_BEGIN 205

#define AROMATIC_AROMATIC_RING_FUSION_START 400
#define AROMATIC_ALIPHATIC_RING_FUSION_COUNT 499
#define AROMATIC_ALIPHATIC_RING_FUSION_START 500
#define ALIPHATIC_ALIPHATIC_RING_FUSION_START 600

#define AROMATIC_TO_ALIPHATIC_MIN 700
#define AROMATIC_TO_ALIPHATIC_MAX 701
#define AROMATIC_TO_ALIPHATIC_AVE 702

#define STRONGLY_FUSED_BIT 990
#define STRONGLY_FUSED_SYSTEM_SIZE 991

#define NRINGS_BIT 1000
#define AROMATIC_RINGS_BIT 1001
#define ALIPHATIC_RINGS_BIT 1002

#define AROMATIC_ATOMS_BIT 1010
#define ALIPHATIC_RING_ATOMS_BIT 1011

#define BIPHENYL_BIT 1020
#define TO_ALIPHATIC_RING_BIT 1021

#define MIXED_AROMATIC_ALIPHATIC_SYSTEM 1030
#define MIXED_AROMATIC_ALIPHATIC_SYSTEM_COUNT 1031

#define FUSED_SYSTEM_SIZE_START 1050

#define ALIPHATIC_HETEROATOM_BIT 1500
#define ALPHATIC_WITHIN_RING_DOUBLE_BOND_BIT 1501
#define ALPHATIC_SUBSTITUENTS_BIT 1502
#define ALIPHATIC_TO_AROMATIC_RING_BIT 1503
#define ALIPHATIC_TO_ALIPHATIC_RING_BIT 1504
#define AROMATIC_ATOMS_IN_ALIPHATIC_RING 1505

#define PLANAR_ALIPHATIC_SYSTEM_SIZE 1520

#define RING_BOND_COUNT_2 1560
#define RING_BOND_COUNT_3 1561

int
Ring_Fingerprint::_fingerprint_aromatic_ring (Molecule & m,
                                              const int * arom,
                                              const Ring & r,
                                              Sparse_Fingerprint_Creator & sfc) const
{
  const auto rsize = r.size();

  int c = 0;
  int nh0 = 0;
  int nh = 0;
  int s = 0;
  int o = 0;

  int heteroatoms = 0;
  int consecutive_heteroatoms = 0;
  int prev_was_carbon = (6 == m.atomic_number(r[rsize-1]));
  int exocyclic_double_bond = 0;
  int more_than_two_rings = 0;
  int substituents = 0;

  int biphenyl = 0;
  int to_aliphatic_ring = 0;

  for (unsigned int i = 0; i < rsize; ++i)
  {
    const auto j = r[i];

    const auto z = m.atomic_number(j);

    if (6 == z)
      c++;
    else if (7 == z)
    {
      if (0 == m.hcount(j))
        nh0++;
      else
        nh++;
    }
    else if (8 == z)
      o++;
    else if (16 == z)
      s++;

    if (6 != z)
      heteroatoms++;

    const auto jcon = m.ncon(j);

    if (3 == jcon)
    {
      _discern_exocyclic_bond(m, arom, j, exocyclic_double_bond, substituents, biphenyl, to_aliphatic_ring);
    }

    if (3 == jcon && m.nrings(j) > 1)
      more_than_two_rings++;

    if (! prev_was_carbon && 6 != z)
      consecutive_heteroatoms++;

    prev_was_carbon = (6 == z);
  }

  sfc.hit_bit(AROMATIC_CARBON, c);
  if (o)
    sfc.hit_bit(AROMATIC_OXYGEN, o);
  if (s)
    sfc.hit_bit(AROMATIC_SULPHUR, s);

  if (o + s)
    sfc.hit_bit(AROMATIC_OS, o + s);

  if (nh0)
    sfc.hit_bit(AROMATIC_NITROGEN_H0, nh0);
  if (nh)
    sfc.hit_bit(AROMATIC_NITROGEN_H, nh);

  if (nh0 + nh)
    sfc.hit_bit(AROMATIC_NITROGEN_ANYH, nh0 + nh);

  if (more_than_two_rings)
    sfc.hit_bit(AROMATIC_MORE_THAN_TWO_RINGS, more_than_two_rings);

  if (exocyclic_double_bond)
    sfc.hit_bit(AROMATIC_EXOCYCLIC_DOUBLE_BOND, exocyclic_double_bond);

  if (substituents)
  {
    sfc.hit_bit(AROMATIC_SUBSTITUENTS, substituents);
    sfc.hit_bit(AROMATIC_SUBSTITUENTS_SIZE, rsize * substituents);
  }

  if (consecutive_heteroatoms)
    sfc.hit_bit(AROMATIC_CONSECUTIVE_HETEROATOMS, consecutive_heteroatoms);

  if (heteroatoms)
  {
    sfc.hit_bit(AROMATIC_HETEROCYCLE);
    sfc.hit_bit(AROMATIC_HETEROATOMS, heteroatoms);
    sfc.hit_bit(AROMATIC_HETEROCYCLE_SIZE_BEGIN + 10 * rsize + heteroatoms, heteroatoms);
  }

  if (biphenyl)
    sfc.hit_bit(BIPHENYL_BIT, biphenyl);

  if (to_aliphatic_ring)
    sfc.hit_bit(TO_ALIPHATIC_RING_BIT, to_aliphatic_ring);

  if (0 == r.fused_ring_neighbours())
    sfc.hit_bit(AROMATIC_ISOLATED_RING);
  else
    _aromatic_ring_fusion_bits (m, r, sfc);

  if (m.non_sssr_rings())
    sfc.hit_bit(NON_SSSR_RINGS, m.non_sssr_rings());

  return 1;
}

int
Ring_Fingerprint::_do_aliphatic_ring_bits (Molecule & m,
                                           const int * arom,
                                           Sparse_Fingerprint_Creator & sfc) const
{
  const auto nr = m.nrings();

  int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  const auto matoms = m.natoms();

  int * fsid = new int[matoms]; std::unique_ptr<int[]> free_fsid(fsid);

  const auto ring_systems = m.label_atoms_by_ring_system_including_spiro_fused(fsid);

  for (auto i = 0; i < nr; ++i)
  {
    if (ring_already_done[i])
      continue;

    const auto ri = m.ringi(i);

    if (ri->is_aromatic())
      continue;

    if (! ri->is_fused())
    {
      _do_single_aliphatic_ring(m, *ri, sfc);
      continue;
    }

    resizable_array<int> rings_in_system;
    rings_in_system.add(i);
    ring_already_done[i] = 1;

    int aromatic_rings_in_system = 0;

    int strongly_fused = ri->strongly_fused_ring_neighbours();

    for (auto j = 0; j < nr; ++j)
    {
      if (ring_already_done[j])
        continue;

      const auto rj = m.ringi(j);

      if (fsid[ri->item(0)] != fsid[rj->item(0)])
        continue;

      ring_already_done[j] = 1;

      if (rj->is_aromatic())
        aromatic_rings_in_system++;

      if (rj->strongly_fused_ring_neighbours())
        strongly_fused++;

      rings_in_system.add(j);
    }

    if (aromatic_rings_in_system)   // should do something better with these, they are common
    {
      sfc.hit_bit(MIXED_AROMATIC_ALIPHATIC_SYSTEM, 1);
      sfc.hit_bit(MIXED_AROMATIC_ALIPHATIC_SYSTEM_COUNT, aromatic_rings_in_system);
    }

    if (0 == strongly_fused)
      _do_planar_aliphatic_system(m, rings_in_system, sfc);
    else
      _do_strongly_fused_aliphatic_system(m, rings_in_system, sfc);
  }

  return 1;
}

void
Ring_Fingerprint::_do_single_aliphatic_ring (Molecule & m,
                                           const Ring & r,
                                           Sparse_Fingerprint_Creator & sfc) const
{
  const auto rsize = r.size();

  int heteroatoms = 0;
  int exocyclic_double_bond = 0;
  int substituents = 0;
  int to_aliphatic_ring = 0;
  int to_aromatic_ring = 0;
  int within_ring_double_bonds = 0;
  int aromatic_atoms = 0;

  for (unsigned int i = 0; i < rsize; ++i)
  {
    const auto j = r[i];

    if (m.is_aromatic(j))
      aromatic_atoms++;

    const auto z = m.atomic_number(j);

    if (6 != z)
      heteroatoms++;

    _discern_exocyclic_bond_aliph(m, j, exocyclic_double_bond, substituents, to_aliphatic_ring, to_aromatic_ring, within_ring_double_bonds);
  }

//cerr << "Ring contains " << heteroatoms << " heteroatoms\n";

  if (heteroatoms && 0 == r.strongly_fused_ring_neighbours())
    sfc.hit_bit(ALIPHATIC_HETEROATOM_BIT, heteroatoms);

  if (within_ring_double_bonds)
    sfc.hit_bit(ALPHATIC_WITHIN_RING_DOUBLE_BOND_BIT, within_ring_double_bonds);
  if (substituents)
    sfc.hit_bit(ALPHATIC_SUBSTITUENTS_BIT, substituents);
  if (to_aromatic_ring)
    sfc.hit_bit(ALIPHATIC_TO_AROMATIC_RING_BIT, to_aromatic_ring);
  if (to_aliphatic_ring)
    sfc.hit_bit(ALIPHATIC_TO_ALIPHATIC_RING_BIT, to_aliphatic_ring);
  if (aromatic_atoms)
    sfc.hit_bit(AROMATIC_ATOMS_IN_ALIPHATIC_RING, aromatic_atoms);

  return;
}

void
Ring_Fingerprint::_do_planar_aliphatic_system(Molecule & m,
                                const resizable_array<int> & rings_in_system, 
                                Sparse_Fingerprint_Creator & sfc) const
{
  const auto nr = rings_in_system.size();

  sfc.hit_bit(PLANAR_ALIPHATIC_SYSTEM_SIZE, nr);

  for (unsigned int i = 0; i < nr; ++i)
  {
    const auto r = m.ringi(rings_in_system[i]);

    if (r->is_aromatic())
      continue;

    _do_single_aliphatic_ring(m, *r, sfc);
  }

  return;
}

// implement this sometime

void
Ring_Fingerprint::_do_strongly_fused_aliphatic_system(Molecule & m,
                                        const resizable_array<int> & rings_in_system,
                                        Sparse_Fingerprint_Creator & sfc) const
{
  const auto nr = rings_in_system.size();

  sfc.hit_bit(STRONGLY_FUSED_SYSTEM_SIZE, nr);

  return;
}

void
Ring_Fingerprint::_ring_fusion_bit (int bstart,
                                    const int rs1,
                                    const int rs2,
                                    Sparse_Fingerprint_Creator & sfc) const
{
  sfc.hit_bit(bstart);    // indicates presence of these kinds of fused rings

  if (rs1 < rs2)
    sfc.hit_bit(bstart + rs1 * 8 + rs2);
  else
    sfc.hit_bit(bstart + rs2 * 8 + rs1);

  return;
}

int
Ring_Fingerprint::_aromatic_ring_fusion_bits (Molecule & m,
                                              const Ring & r,
                                              Sparse_Fingerprint_Creator & sfc) const
{
  const auto rsize = r.size();

  const auto nf = r.fused_ring_neighbours();

  if (0 == nf)
    return 1;

  const auto h1 = m.count_heteroatoms(r);

  int aliphatic_fused_neighbours = 0;

  for (auto i = 0; i < nf; ++i)
  {
    const auto fnbr = r.fused_neighbour(i);

//  if (fnbr->ring_number() < r.ring_number())     breaks things, not sure why
//    continue;

    const auto rs2 = fnbr->size();

    if (fnbr->is_aromatic())
    {
      _ring_fusion_bit(AROMATIC_AROMATIC_RING_FUSION_START, rsize, rs2, sfc);
      const auto h2 = m.count_heteroatoms(*fnbr);
      if (h1 && h2)
        sfc.hit_bit(AROMATIC_HETEROCYCLE_FUSED_AROMATIC_HETEROCYCLE);
      else if (h1 || h2)
        sfc.hit_bit(AROMATIC_HETEROCYCLE_FUSED_AROMATIC_CARBOCYCLE);
      else
        sfc.hit_bit(AROMATIC_CARBOCYCLE_FUSED_AROMATIC_CARBOCYCLE);
    }
    else
    {
      aliphatic_fused_neighbours++;
      sfc.hit_bit(AROMATIC_ALIPHATIC_RING_FUSION_START, rs2);
      _ring_fusion_bit(AROMATIC_ALIPHATIC_RING_FUSION_START, rsize, rs2, sfc);
    }
  }

  sfc.hit_bit(AROMATIC_ALIPHATIC_RING_FUSION_COUNT, aliphatic_fused_neighbours);

  return 1;
}

/*
  Here's where you do whatever you want to do with the molecule
  In this case, we count the number of nitrogen atoms
*/

int
Ring_Fingerprint::_process (Molecule & m,
                            IWString_and_File_Descriptor & output)
{
  Sparse_Fingerprint_Creator sfc;

  _form_fingerprint(m, sfc);

  return _write_fingerprint(m, sfc, output);
}

void
Ring_Fingerprint::_form_fingerprint (Molecule & m,
                                     Sparse_Fingerprint_Creator & sfc) const
{

  const auto nr = m.nrings();

  sfc.hit_bit(NRINGS_BIT, nr + 1);

  if (0 == nr)
    return;

  m.compute_aromaticity_if_needed();

  const auto matoms = m.natoms();

  int * arom = new int[matoms]; std::unique_ptr<int[]> free_arom(arom);
  m.aromaticity(arom);

  int aromatic_atoms = 0;
  int aliphatic_ring_atoms = 0;

  int ring_bond_count_2 = 0;
  int ring_bond_count_3 = 0;

  for (auto i = 0; i < matoms; ++i)
  {
    if (arom[i])
      aromatic_atoms++;
    else if (m.nrings(i))
      aliphatic_ring_atoms++;
    else
      continue;

    const int rbc = m.ring_bond_count(i);
    if (2 == rbc)
      ring_bond_count_2++;
    else
      ring_bond_count_3++;
  }

  sfc.hit_bit(AROMATIC_ATOMS_BIT, aromatic_atoms);
  sfc.hit_bit(ALIPHATIC_RING_ATOMS_BIT, aliphatic_ring_atoms);
  sfc.hit_bit(RING_BOND_COUNT_2, ring_bond_count_2);
  sfc.hit_bit(RING_BOND_COUNT_3, ring_bond_count_3);

  extending_resizable_array<int> aromatic_ring_count, aliphatic_ring_count;

  int number_aromatic_rings = 0;
  int largest_ring_size = 0;

  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = m.ringi(i);

    if (ri->strongly_fused_ring_neighbours())
      sfc.hit_bit(STRONGLY_FUSED_BIT, ri->strongly_fused_ring_neighbours());

    const auto rsize = ri->size();

    if (ri->is_aromatic())
    {
      aromatic_ring_count[rsize]++;
      sfc.hit_bit(AROMATIC_RING_BIT);
      _fingerprint_aromatic_ring (m, arom, *ri, sfc);
      number_aromatic_rings++;
    }
    else if (_form_aliphatic_ring_bits)
    {
      if (0 == ri->fused_ring_neighbours())
        sfc.hit_bit(ALIPHATIC_ISOLATED_RING);

      aliphatic_ring_count[rsize]++;
      sfc.hit_bit(ALIPHATIC_RING_BIT);
    }

    if (rsize > largest_ring_size)
      largest_ring_size = rsize;
  }

  sfc.hit_bit(LARGEST_RING_SIZE, largest_ring_size);

  for (unsigned int i = 4; i < aromatic_ring_count.size(); ++i)
  {
    if (0 == aromatic_ring_count[i])
      continue;

    sfc.hit_bit(AROMATIC_RING_SIZE_START + i, aromatic_ring_count[i]);
  }

  for (unsigned int i = 4; i < aliphatic_ring_count.size(); ++i)
  {
    if (0 == aliphatic_ring_count[i])
      continue;

    sfc.hit_bit(ALIPHATIC_RING_SIZE_START + i, aliphatic_ring_count[i]);
  }

  if (_form_aliphatic_ring_bits && number_aromatic_rings < nr)
    _do_aliphatic_ring_bits (m, arom, sfc);

  if (_form_bits_for_inter_aromatic_ring_distances && number_aromatic_rings && number_aromatic_rings < nr)
    _bits_for_distances_aromatic_to_aliphatic(m, arom, sfc);

  if (nr < 2)
    return;

  if (_form_bits_for_inter_aromatic_ring_distances && number_aromatic_rings > 1)
    _bits_for_distances_between_aromatic_rings(m, sfc);

  int * in_system = new_int(matoms); std::unique_ptr<int[]> free_in_system(in_system);

  for (int i = 0; i < nr; ++i)
  {
    const Ring * ri = m.ringi(i);
    if (!  ri->is_fused())
    {
      sfc.hit_bit(FUSED_SYSTEM_SIZE_START + 1);
      continue;
    }
    
    if (ri->any_members_set_in_array(in_system))
      continue;

    ri->set_vector(in_system, 1);

    int rings_in_system = 1;

    for (int j = i + 1; j < nr; ++j)
    {
      const Ring * rj = m.ringi(j);

      if (ri->fused_system_identifier() != rj->fused_system_identifier())
        continue;

      rj->set_vector(in_system, 1);

      rings_in_system++;
    }

    sfc.hit_bit(FUSED_SYSTEM_SIZE_START + rings_in_system);
  }

  return;
}

void
Ring_Fingerprint::_bits_for_distances_between_aromatic_rings (Molecule & m,
                                Sparse_Fingerprint_Creator & sfc) const
{
  const auto nr = m.nrings();
  assert (nr > 1);

  Accumulator_Int<int> between_rings;

  for (auto i = 0; i < nr; ++i)
  {
    const auto ri = m.ringi(i);

    if (! ri->is_aromatic())
      continue;

    for (auto j = i + 1; j < nr; ++j)
    {
      const auto rj = m.ringi(j);

      if (! rj->is_aromatic())
        continue;

      const int d = (ri->is_fused_to(rj) ? 0 : _shortest_distance_between(m, *ri, *rj));

      sfc.hit_bit(BETWEEN_AROMATIC_RINGS_BEGIN + d);

      between_rings.extra(d);
//    cerr << "BTW ring " << i << " and " << j << " dist " << d << endl;
    }
  }

//cerr << "For " << m.name() << " n = " << between_rings.n() << " min " << between_rings.minval() << " max " << between_rings.maxval() << " ave " << between_rings.average() << endl;

  if (0 == between_rings.n())    // should not happen
    return;

  sfc.hit_bit(BETWEEN_AROMATIC_RINGS_MIN, between_rings.minval());
  sfc.hit_bit(BETWEEN_AROMATIC_RINGS_MAX, between_rings.maxval());
  sfc.hit_bit(BETWEEN_AROMATIC_RINGS_AVE, static_cast<int>(between_rings.average() + 0.49999));

  return;
}

void
Ring_Fingerprint::_bits_for_distances_aromatic_to_aliphatic (Molecule & m,
                                                const int * arom,
                                                Sparse_Fingerprint_Creator & sfc) const
{
  const auto matoms = m.natoms();

  Accumulator_Int<int> arom_2_aliph;
  for (auto i = 0; i < matoms; ++i)
  {
    if (! arom[i])
      continue;

    for (auto j = 0; j < matoms; ++j)
    {
      if (arom[j])
        continue;

      if (0 == m.nrings(j))
        continue;

      const auto d = m.bonds_between(i, j);

      arom_2_aliph.extra(d);
    }
  }

  if (0 == arom_2_aliph.n())    // should not happen
  {
    cerr << "No aromatic - aliphatic distances in '" << m.smiles() << endl;
    return;
  }

  sfc.hit_bit(AROMATIC_TO_ALIPHATIC_MIN, arom_2_aliph.minval());
  sfc.hit_bit(AROMATIC_TO_ALIPHATIC_MAX, arom_2_aliph.maxval());
  sfc.hit_bit(AROMATIC_TO_ALIPHATIC_AVE, static_cast<int>(arom_2_aliph.maxval() + 0.4999));

  return;
}

int
Ring_Fingerprint::_shortest_distance_between (Molecule & m,
                                              const Set_of_Atoms & r1,
                                              const Set_of_Atoms & r2) const
{
  const auto rs1 = r1.size();
  const auto rs2 = r2.size();

  int rc = m.natoms();

  for (unsigned int i = 0; i < rs1; ++i)
  {
    const auto j = r1[i];

    for (unsigned int k = 0; k < rs2; ++k)
    {
      const auto d = m.bonds_between(j, r2[k]);

      if (d < rc)
        rc = d;
    }
  }

  return rc;
}

int
Ring_Fingerprint::_process_filter (iwstring_data_source & input,
                                   IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    output << buffer << '\n';

    if (! buffer.starts_with(smiles_tag))
      continue;

    buffer.remove_leading_chars(smiles_tag.length());
    buffer.chop();

    Molecule m;

    if (! m.build_from_smiles(buffer))
    {
      cerr << "Ring_Fingerprint::_process_filter:cannot interpret smiles '" << buffer << "'\n";
      return 0;
    }

    if (! _process (m, output))
      return 0;
  }

  return 1;
}

int
Ring_Fingerprint::_process_filter (const char * fname,
                                   IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good ())
  {
    cerr << "Ring_Fingerprint::_process_filter:cannot open '" << fname << "'\n";
    return 0;
  }

  return _process_filter (input, output);
}

int
Ring_Fingerprint::_process (data_source_and_type<Molecule> & input,
                            IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    _molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    _preprocess(*m);

    if (_ntest > 0)
    {
      if (! _perform_tests(*m))
      {
        cerr << "Test failure on molecule '" << m->name() << "'\n";
        return 0;
      }
    }
    else if (! _process(*m, output))
      return 0;
  }

  return 1;
}

int
Ring_Fingerprint::_process (const char * fname, int input_type,
                 IWString_and_File_Descriptor & output)
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

  if (_verbose > 1)
    input.set_verbose(1);

  return _process(input, output);
}

int
Ring_Fingerprint::operator() (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:lwt:s:depfJ:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    _usage(1);
  }

  _verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, _verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      _usage(5);
    }
  }
  else
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, _verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      _usage(32);
    }
  }

  if (cl.option_present('l'))
  {
    _reduce_to_largest_fragment = 1;

    if (_verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('J'))
  {
    cl.value('J', _tag);

    if (! _tag.ends_with('<'))
      _tag << '<';

    if (_verbose)
      cerr << "Fingerprint written with tag '" << _tag << "'\n";
  }

  int input_type = 0;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      _usage (6);
    }
  }
  else if (cl.option_present('f'))
    _work_as_tdt_filter = 1;
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.option_present('t'))
  {
    if (! cl.value('t', _ntest) || _ntest < 1)
    {
      cerr << "The number of tests to perform (-t) must be a whole +ve number\n";
      _usage(1);
    }

    if (_verbose)
      cerr << "WIll perform " << _ntest << " tests\n";

    if (cl.option_present('s'))
    {
      random_number_seed_t s;
      if (! cl.value('s', s))
      {
        cerr << "The random smiles seed (-s) must be a valid random number seed\n";
        _usage(1);
      }

      set_smiles_random_number_seed(s);
    }
    else
      set_smiles_random_number_seed_random();
  }

  if (cl.option_present('d'))
  {
    _form_bits_for_inter_aromatic_ring_distances = 1;

    if (_verbose)
      cerr << "Will form inter-ring distance bits\n";
  }

  if (cl.option_present('e'))
  {
    _echo_fingerprints = 1;

    if (_verbose)
      cerr << "Will echo computed fingerprints\n";
  }

  if (cl.option_present('p'))
  {
    _form_aliphatic_ring_bits = 1;

    if (_verbose)
      cerr << "Will produce bits based on aliphatic rings\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    _usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  if (_work_as_tdt_filter)
  {
    if (! _process_filter(cl[0], output))
      rc = 1;
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! _process(cl[i], input_type, output))
      {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (_verbose)
  {
    cerr << "Read " << _molecules_read << " molecules\n";
    if (_ntest > 0 && 0 == rc)
      cerr << "All tests successful\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  Ring_Fingerprint ring_fingerprint;
  
  return ring_fingerprint (argc, argv);
}
